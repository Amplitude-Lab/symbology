#pragma once
#include <atomic>
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#ifdef USE_MIMALLOC
#include <mimalloc.h>
#else
#include <malloc.h>
#endif

// Large SparseRREF arrays can be backed by a private scratch file. Shared
// file pages retain exact bytes when residency is reduced; all small arrays
// and GMP/FLINT limb allocations retain their usual allocator.
namespace hexagon_mapped_buffers {
inline void* native_malloc(size_t n) {
#ifdef USE_MIMALLOC
    return mi_malloc(n);
#else
    return std::malloc(n);
#endif
}
inline void native_free(void* p) {
#ifdef USE_MIMALLOC
    mi_free(p);
#else
    std::free(p);
#endif
}
inline void* native_realloc(void* p,size_t n) {
#ifdef USE_MIMALLOC
    return mi_realloc(p,n);
#else
    return std::realloc(p,n);
#endif
}
inline size_t native_size(void* p) {
#ifdef USE_MIMALLOC
    return mi_usable_size(p);
#else
    return malloc_usable_size(p);
#endif
}
inline void* native_expand(void* p,size_t n) {
#ifdef USE_MIMALLOC
    return mi_expand(p,n);
#else
    return std::realloc(p,n);
#endif
}
class Pool {
    struct Block {size_t offset,capacity,requested;};
    int fd=-1; char* base=nullptr;size_t capacity=0,top=0,soft=0;
    std::mutex mutex,eviction;
    std::unordered_map<void*,Block> live;
    std::multimap<size_t,size_t> free_blocks;
    std::atomic<size_t> published_top{0};
    std::atomic<int64_t> next_check{0};
    size_t evictions=0;
    static size_t align(size_t n){return (n+63)&~size_t(63);}
    void* allocate_locked(size_t n){
        size_t bytes=align(n),offset=0;
        auto old=free_blocks.lower_bound(bytes);
        if(old!=free_blocks.end()){
            offset=old->second;size_t remaining=old->first-bytes;free_blocks.erase(old);
            if(remaining)free_blocks.emplace(remaining,offset+bytes);
        }else{
            if(bytes>capacity-top)throw std::runtime_error("Mapped scratch capacity exhausted");
            offset=top;top+=bytes;published_top.store(top,std::memory_order_release);
        }
        void* p=base+offset;live.emplace(p,Block{offset,bytes,n});return p;
    }
    void free_locked(void* p){
        auto it=live.find(p);if(it==live.end())throw std::runtime_error("Invalid mapped buffer owner");
        free_blocks.emplace(it->second.capacity,it->second.offset);live.erase(it);
    }
public:
    Pool(){
        const char* name=std::getenv("HEXAGON_MAPPED_FILE");if(!name||!*name)return;
        const char* gib=std::getenv("HEXAGON_MAPPED_GIB");
        capacity=(gib?std::stoull(gib):192)*size_t(1024)*1024*1024;
        const char* budget=std::getenv("HEXAGON_MAPPED_SOFT_GIB");
        soft=(budget?std::stoull(budget):16)*size_t(1024)*1024*1024;
        std::filesystem::path path(name);std::filesystem::create_directories(path.parent_path());
        fd=::open(name,O_RDWR|O_CREAT|O_EXCL,0600);
        if(fd<0)throw std::runtime_error("Mapped scratch file must be fresh");
        if(::ftruncate(fd,off_t(capacity)))throw std::runtime_error("Cannot size mapped scratch file");
        void* mapping=::mmap(nullptr,capacity,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
        if(mapping==MAP_FAILED)throw std::runtime_error("Cannot reserve mapped scratch address space");
        base=static_cast<char*>(mapping);
        std::fprintf(stdout,"MAPPED_BUFFERS capacity_bytes=%zu soft_rss_bytes=%zu threshold_bytes=4096\n",capacity,soft);
        std::fflush(stdout);
    }
    ~Pool(){if(base)::munmap(base,capacity);if(fd>=0)::close(fd);}
    bool enabled()const{return base!=nullptr;}
    bool owns(void* p)const{
        auto value=reinterpret_cast<uintptr_t>(p),begin=reinterpret_cast<uintptr_t>(base);
        return base&&value>=begin&&value-begin<capacity;
    }
    void* allocate(size_t n){std::lock_guard lock(mutex);return allocate_locked(n);}
    void release(void* p){std::lock_guard lock(mutex);free_locked(p);}
    void* reallocate(void* p,size_t n){
        if(!n){release(p);return nullptr;}
        std::lock_guard lock(mutex);auto it=live.find(p);
        if(it==live.end())throw std::runtime_error("Invalid mapped realloc owner");
        auto old=it->second;
        if(n<=old.capacity){it->second.requested=n;return p;}
        void* q=allocate_locked(n);std::memcpy(q,p,std::min(n,old.requested));free_locked(p);return q;
    }
    void* expand(void* p,size_t n){
        std::lock_guard lock(mutex);auto it=live.find(p);
        if(it==live.end())throw std::runtime_error("Invalid mapped expand owner");
        if(n>it->second.capacity)return nullptr;
        it->second.requested=n;return p;
    }
    void maybe_evict(bool force=false){
        if(!base)return;
        auto now=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
        auto next=next_check.load(std::memory_order_relaxed);
        if(!force&&(now<next||!next_check.compare_exchange_strong(next,now+1000000000LL)))return;
        std::ifstream status("/proc/self/statm");size_t virtual_pages=0,resident=0;status>>virtual_pages>>resident;
        size_t rss=resident*size_t(::sysconf(_SC_PAGESIZE));if(!status||(!force&&rss<soft))return;
        std::unique_lock lock(eviction,std::try_to_lock);if(!lock.owns_lock())return;
        size_t length=published_top.load(std::memory_order_acquire);if(!length)return;
        std::fprintf(stdout,"MAPPED_EVICT begin rss_bytes=%zu used_address_bytes=%zu\n",rss,length);std::fflush(stdout);
        // MAP_SHARED pages are synchronized before dropping residency. The
        // file, pointers and coefficient bytes remain valid afterward.
        if(::msync(base,length,MS_SYNC))throw std::runtime_error("Mapped scratch synchronization failed");
        if(::madvise(base,length,MADV_DONTNEED))throw std::runtime_error("Mapped scratch residency reduction failed");
        ++evictions;std::fprintf(stdout,"MAPPED_EVICT end count=%zu\n",evictions);std::fflush(stdout);
    }
};
inline Pool& pool(){static Pool instance;return instance;}
inline void* allocate(size_t n){
    auto& p=pool();p.maybe_evict();
    return p.enabled()&&n>=4096?p.allocate(n):native_malloc(n);
}
inline void release(void* pointer){
    auto& p=pool();p.maybe_evict();if(p.owns(pointer))p.release(pointer);else native_free(pointer);
}
inline void* reallocate(void* pointer,size_t n){
    auto& p=pool();p.maybe_evict();if(p.owns(pointer))return p.reallocate(pointer,n);
    if(p.enabled()&&pointer&&n>=4096){
        void* next=p.allocate(n);std::memcpy(next,pointer,std::min(n,native_size(pointer)));native_free(pointer);return next;
    }
    return p.enabled()&&!pointer&&n>=4096?p.allocate(n):native_realloc(pointer,n);
}
inline void* expand(void* pointer,size_t n){
    auto& p=pool();p.maybe_evict();return p.owns(pointer)?p.expand(pointer,n):native_expand(pointer,n);
}
}
