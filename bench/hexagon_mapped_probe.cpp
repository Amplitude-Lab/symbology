#include "irr_physics_pipeline.hpp"
#include <thread>

int main()try {
    Flint::set_memory_functions();
    irr::require(hexagon_mapped_buffers::pool().enabled(),"Set HEXAGON_MAPPED_FILE for the allocator probe");
    auto* a=s_malloc<ulong>(4096);
    for(size_t i=0;i<4096;++i)a[i]=ulong(i*i+17);
    a=s_realloc(a,8192);
    for(size_t i=0;i<4096;++i)irr::require(a[i]==ulong(i*i+17),"Mapped grow lost bytes");
    for(size_t i=4096;i<8192;++i)a[i]=ulong(i*i+17);
    hexagon_mapped_buffers::pool().maybe_evict(true);
    for(size_t i=0;i<8192;++i)irr::require(a[i]==ulong(i*i+17),"Mapped eviction lost bytes");
    a=s_realloc(a,2048);auto* expanded=s_expand(a,4096);
    irr::require(expanded==a,"Mapped in-place expansion");
    for(size_t i=0;i<2048;++i)irr::require(a[i]==ulong(i*i+17),"Mapped shrink lost bytes");
    s_free(a);
    auto* b=s_malloc<ulong>(4);for(size_t i=0;i<4;++i)b[i]=i+83;
    b=s_realloc(b,8192);for(size_t i=0;i<4;++i)irr::require(b[i]==i+83,"Heap to mapped conversion");s_free(b);
    auto* q=s_malloc<rat_t>(512);
    for(size_t i=0;i<512;++i)new(q+i) rat_t(std::string("340282366920938463463374607431768211")+std::to_string(i+1)+"/170141183460469231731687303715884"+std::to_string(i+2));
    q=s_realloc(q,1024);
    for(size_t i=512;i<1024;++i)new(q+i) rat_t(std::string("340282366920938463463374607431768211")+std::to_string(i+1)+"/170141183460469231731687303715884"+std::to_string(i+2));
    hexagon_mapped_buffers::pool().maybe_evict(true);
    for(size_t i=0;i<1024;++i){irr::require(q[i]==rat_t(std::string("340282366920938463463374607431768211")+std::to_string(i+1)+"/170141183460469231731687303715884"+std::to_string(i+2)),"Rational relocation or eviction");q[i].~rat_t();}
    s_free(q);
    std::atomic<bool> valid{true};std::vector<std::thread> workers;
    for(size_t t=0;t<4;++t)workers.emplace_back([&,t]{
        for(size_t pass=0;pass<32;++pass){auto* data=s_malloc<ulong>(8192);
            for(size_t i=0;i<8192;++i)data[i]=(t+1)*1000000+pass*8192+i;
            if(pass%8==0)hexagon_mapped_buffers::pool().maybe_evict(true);
            for(size_t i=0;i<8192;++i)if(data[i]!=(t+1)*1000000+pass*8192+i)valid=false;
            s_free(data);
        }
    });
    for(auto& thread:workers)thread.join();irr::require(valid,"Concurrent mapped buffers");
    std::cout<<"MAPPED_PROBE PASS grow=1 shrink=1 heap_transition=1 eviction=1 rational_arrays=1 concurrent_buffers=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
