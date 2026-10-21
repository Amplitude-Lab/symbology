#define USE_MIMALLOC
#ifndef HEXAGON_MAPPED_PROBE_HEADER
#define HEXAGON_MAPPED_PROBE_HEADER "shared/mapped_buffers.hpp"
#endif
#include HEXAGON_MAPPED_PROBE_HEADER
#include <thread>
#include <vector>
#include <iostream>
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
int main()try{
    auto& arena=hexagon_mapped_buffers::pool();
    require(arena.enabled(),"Watchdog probe requires file-backed buffers");
    const size_t count=(size_t(1100)*1024*1024)/sizeof(uint64_t);
    auto* p=static_cast<uint64_t*>(hexagon_mapped_buffers::allocate(count*sizeof(uint64_t)));
    volatile uint64_t* writable=p;
    const size_t before=arena.eviction_count();
    std::vector<std::thread> workers;
    for(size_t t=0;t<4;++t)workers.emplace_back([=]{
        for(uint64_t pass=0;pass<8;++pass)
            for(size_t i=count*t/4;i<count*(t+1)/4;++i)
                writable[i]=uint64_t(i)*UINT64_C(0x9e3779b97f4a7c15)+pass;
    });
    for(auto& worker:workers)worker.join();
    // No SparseRREF allocation/free callbacks occur while waiting. A
    // background eviction must therefore handle the newly resident data.
    auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(30);
    while(arena.eviction_count()<=before&&std::chrono::steady_clock::now()<deadline)
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    require(arena.eviction_count()>before,"Continuous pager did not evict without allocator callbacks");
    for(size_t i=0;i<count;++i)
        require(p[i]==uint64_t(i)*UINT64_C(0x9e3779b97f4a7c15)+7,"Background eviction changed coefficient bytes");
    hexagon_mapped_buffers::release(p);
    std::cout<<"WATCHDOG_PROBE PASS no_allocator_callbacks=1 concurrent_writes=1 bytes="
        <<count*sizeof(uint64_t)<<" evictions="<<arena.eviction_count()<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
