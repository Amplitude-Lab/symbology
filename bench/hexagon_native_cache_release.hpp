#pragma once
#include <barrier>
#include <atomic>
#include <fstream>
#include <sstream>

// FLINT 3.2.2 keeps reusable multiprecision integers in thread-local caches.
// mi_collect alone cannot release their still-owned limbs. At a joined phase
// boundary, release only caches, retaining every live QQ and Fp object.
// Source: FLINT memory documentation and v3.2.2 fmpz/link/fmpz_single.c.
namespace hexagon_native_memory {
inline size_t resident_bytes(){
    std::ifstream in("/proc/self/status");std::string line;
    while(std::getline(in,line)){std::istringstream row(line);std::string key;size_t kib;
        if(row>>key>>kib&&key=="VmRSS:")return kib*1024;}
    return 0;
}
inline void release_unused_integer_caches(thread_pool* pool){
#if !FLINT_USES_TLS
    throw std::runtime_error("Integer cache release requires FLINT thread-local storage");
#else
    if(!pool)throw std::runtime_error("Integer cache release needs the owned pool");
    if(BS::this_thread::get_pool().has_value())
        throw std::runtime_error("Integer cache release must run outside a worker task");
    pool->wait();
    auto start=std::chrono::steady_clock::now();
    auto before=resident_bytes();
    const auto workers=pool->get_thread_count();
    if(!workers)throw std::runtime_error("Integer cache release needs nonzero workers");
    std::barrier rendezvous(static_cast<std::ptrdiff_t>(workers));
    std::atomic<size_t> completed=0;
    // One task per worker, held at a barrier, guarantees distinct workers;
    // merely submitting N short tasks does not guarantee thread coverage.
    for(size_t i=0;i<workers;++i)pool->detach_task([&]{
        rendezvous.arrive_and_wait();
        flint_cleanup();
#ifdef USE_MIMALLOC
        mi_collect(true);
#endif
        ++completed;
    });
    pool->wait();
    flint_cleanup();
#ifdef USE_MIMALLOC
    mi_collect(true);
#endif
    if(completed!=workers)throw std::runtime_error("Incomplete integer-cache worker coverage");
    std::cout<<"NATIVE_INTEGER_CACHE_RELEASE workers="<<completed<<" main=1 rss_before="<<before
             <<" rss_after="<<resident_bytes()<<" seconds="
             <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
#endif
}
} // namespace hexagon_native_memory
