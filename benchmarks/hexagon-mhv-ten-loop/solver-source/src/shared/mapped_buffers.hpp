#pragma once
#include <atomic>
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <cerrno>
#include <thread>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <fcntl.h>
#include <linux/falloc.h>
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
    using SizeIndex=std::multimap<size_t,size_t>;
    struct FreeBlock {size_t bytes;SizeIndex::iterator by_size;};
    int fd=-1; char* base=nullptr;size_t capacity=0,top=0,soft=0,page_size=4096,reclaim_min=64*1024;
    size_t live_capacity=0,live_requested=0,shrink_count=0,reclaim_calls=0,reclaim_requested_bytes=0;
    bool reclaim_enabled=true;
    std::mutex mutex,eviction;
    std::unordered_map<void*,Block> live;
    SizeIndex free_blocks;
    std::map<size_t,FreeBlock> free_offsets;
    std::atomic<size_t> published_top{0};
    std::atomic<int64_t> next_check{0};
    std::atomic<size_t> evictions{0};
    std::atomic<bool> stop_watchdog{false};
    std::thread watchdog;
    static size_t align(size_t n){if(n>SIZE_MAX-63)throw std::runtime_error("Mapped allocation size overflow");return (n+63)&~size_t(63);}
    void erase_free(std::map<size_t,FreeBlock>::iterator it){
        free_blocks.erase(it->second.by_size);free_offsets.erase(it);
    }
    void index_free(size_t offset,size_t bytes){
        auto index=free_blocks.emplace(bytes,offset);
        require_unique_free(offset,bytes,index);
    }
    void require_unique_free(size_t offset,size_t bytes,SizeIndex::iterator index){
        if(!free_offsets.emplace(offset,FreeBlock{bytes,index}).second)
            throw std::runtime_error("Duplicate mapped free block");
    }
    void reclaim_free_locked(size_t offset,size_t length){
        if(!reclaim_enabled||!reclaim_min||length<reclaim_min)return;
        // Both boundaries are page aligned and every byte belongs to free
        // storage. Allocation/reuse is excluded by the allocator mutex.
        int result;
        do{result=::fallocate(fd,FALLOC_FL_PUNCH_HOLE|FALLOC_FL_KEEP_SIZE,off_t(offset),off_t(length));}
        while(result&&errno==EINTR);
        if(result){
            int reason=errno;reclaim_enabled=false;
            std::fprintf(stderr,"MAPPED_FREE_RECLAIM disabled errno=%d reason=%s exact_live_bytes_preserved=1\n",reason,std::strerror(reason));
            std::fflush(stderr);return;
        }
        ++reclaim_calls;reclaim_requested_bytes+=length;
        if(reclaim_calls==1||reclaim_calls%4096==0){
            std::fprintf(stdout,"MAPPED_FREE_RECLAIM calls=%zu requested_bytes=%zu live_capacity_bytes=%zu\n",reclaim_calls,reclaim_requested_bytes,live_capacity);
            std::fflush(stdout);
        }
    }
    void insert_free_locked(size_t offset,size_t bytes){
        const size_t released_start=offset,released_end=offset+bytes;
        auto next=free_offsets.lower_bound(offset);
        if(next!=free_offsets.end()&&offset+bytes>next->first)
            throw std::runtime_error("Overlapping mapped free block");
        if(next!=free_offsets.begin()){
            auto previous=std::prev(next);size_t end=previous->first+previous->second.bytes;
            if(end>offset)throw std::runtime_error("Overlapping mapped free predecessor");
            if(end==offset){offset=previous->first;bytes+=previous->second.bytes;erase_free(previous);}
        }
        next=free_offsets.lower_bound(offset);
        if(next!=free_offsets.end()&&offset+bytes==next->first){bytes+=next->second.bytes;erase_free(next);}
        // Reclaim only pages touching the newly freed bytes, including a
        // boundary page that is now entirely free after neighbour coalescing.
        // Avoid repeatedly punching the full, previously reclaimed neighbour.
        size_t start=std::max((released_start/page_size)*page_size,
                              ((offset+page_size-1)/page_size)*page_size);
        size_t end=std::min(((released_end+page_size-1)/page_size)*page_size,
                            ((offset+bytes)/page_size)*page_size);
        if(end>start)reclaim_free_locked(start,end-start);
        if(offset+bytes==top){top=offset;published_top.store(top,std::memory_order_release);}
        else index_free(offset,bytes);
    }
    void* allocate_locked(size_t n){
        size_t bytes=align(n),offset=0;
        auto old=free_blocks.lower_bound(bytes);
        if(old!=free_blocks.end()){
            offset=old->second;size_t remaining=old->first-bytes;erase_free(free_offsets.find(offset));
            if(remaining)index_free(offset+bytes,remaining);
        }else{
            if(bytes>capacity-top)throw std::runtime_error("Mapped scratch capacity exhausted");
            offset=top;top+=bytes;published_top.store(top,std::memory_order_release);
        }
        void* p=base+offset;
        if(!live.emplace(p,Block{offset,bytes,n}).second)throw std::runtime_error("Mapped live block overlap");
        live_capacity+=bytes;live_requested+=n;return p;
    }
    void free_locked(void* p){
        auto it=live.find(p);if(it==live.end())throw std::runtime_error("Invalid mapped buffer owner");
        auto old=it->second;live.erase(it);live_capacity-=old.capacity;live_requested-=old.requested;
        insert_free_locked(old.offset,old.capacity);
    }
public:
    Pool(){
        const char* name=std::getenv("HEXAGON_MAPPED_FILE");if(!name||!*name)return;
        const char* gib=std::getenv("HEXAGON_MAPPED_GIB");
        capacity=(gib?std::stoull(gib):192)*size_t(1024)*1024*1024;
        const char* budget=std::getenv("HEXAGON_MAPPED_SOFT_GIB");
        soft=(budget?std::stoull(budget):16)*size_t(1024)*1024*1024;
        page_size=size_t(::sysconf(_SC_PAGESIZE));
        if(const char* minimum=std::getenv("HEXAGON_MAPPED_RECLAIM_MIN_KIB"))reclaim_min=std::stoull(minimum)*size_t(1024);
        std::filesystem::path path(name);std::filesystem::create_directories(path.parent_path());
        fd=::open(name,O_RDWR|O_CREAT|O_EXCL,0600);
        if(fd<0)throw std::runtime_error("Mapped scratch file must be fresh");
        if(::ftruncate(fd,off_t(capacity)))throw std::runtime_error("Cannot size mapped scratch file");
        void* mapping=::mmap(nullptr,capacity,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
        if(mapping==MAP_FAILED)throw std::runtime_error("Cannot reserve mapped scratch address space");
        base=static_cast<char*>(mapping);
        std::fprintf(stdout,"MAPPED_BUFFERS capacity_bytes=%zu soft_rss_bytes=%zu threshold_bytes=4096\n",capacity,soft);
        std::fprintf(stdout,"MAPPED_FREE_RECLAIM min_bytes=%zu aligned_live_tails_split=1 free_ranges_coalesced=1\n",reclaim_min);
        std::fflush(stdout);
        watchdog=std::thread([this]{
            while(!stop_watchdog.load(std::memory_order_acquire)){
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
                if(stop_watchdog.load(std::memory_order_acquire))break;
                try{maybe_evict(false,true);}
                catch(const std::exception& e){
                    std::fprintf(stderr,"MAPPED_WATCHDOG_ERROR %s\n",e.what());std::fflush(stderr);std::abort();
                }
            }
        });
        std::fprintf(stdout,"MAPPED_WATCHDOG poll_ms=250 checks_per_second=1\n");std::fflush(stdout);
    }
    ~Pool(){stop_watchdog.store(true,std::memory_order_release);if(watchdog.joinable())watchdog.join();if(base)::munmap(base,capacity);if(fd>=0)::close(fd);}
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
        if(n<=old.capacity){
            size_t bytes=align(n);live_requested=live_requested-old.requested+n;it->second.requested=n;
            if(bytes<old.capacity){
                it->second.capacity=bytes;live_capacity-=old.capacity-bytes;++shrink_count;
                insert_free_locked(old.offset+bytes,old.capacity-bytes);
            }
            return p;
        }
        void* q=allocate_locked(n);std::memcpy(q,p,std::min(n,old.requested));free_locked(p);return q;
    }
    void* expand(void* p,size_t n){
        std::lock_guard lock(mutex);auto it=live.find(p);
        if(it==live.end())throw std::runtime_error("Invalid mapped expand owner");
        if(n>it->second.capacity)return nullptr;
        live_requested=live_requested-it->second.requested+n;it->second.requested=n;return p;
    }
    struct Statistics {size_t live_capacity_bytes,live_requested_bytes,reusable_bytes,used_address_bytes,
                              shrink_calls,reclaim_calls,reclaim_requested_bytes;bool reclaim_enabled;};
    Statistics statistics(){
        std::lock_guard lock(mutex);
        return {live_capacity,live_requested,top-live_capacity,top,shrink_count,reclaim_calls,reclaim_requested_bytes,reclaim_enabled&&bool(reclaim_min)};
    }
    size_t eviction_count()const{return evictions.load(std::memory_order_acquire);}
    void maybe_evict(bool force=false,bool from_watchdog=false){
        if(!base||(!force&&!from_watchdog))return;
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
        size_t count=evictions.fetch_add(1,std::memory_order_release)+1;std::fprintf(stdout,"MAPPED_EVICT end count=%zu\n",count);std::fflush(stdout);
        auto s=statistics();
        std::fprintf(stdout,"MAPPED_ALLOCATION live_capacity_bytes=%zu live_requested_bytes=%zu reusable_bytes=%zu used_address_bytes=%zu shrink_calls=%zu reclaim_calls=%zu reclaim_requested_bytes=%zu\n",
                     s.live_capacity_bytes,s.live_requested_bytes,s.reusable_bytes,s.used_address_bytes,s.shrink_calls,s.reclaim_calls,s.reclaim_requested_bytes);
        std::fflush(stdout);
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
