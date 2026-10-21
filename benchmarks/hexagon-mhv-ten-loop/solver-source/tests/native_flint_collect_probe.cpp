#ifndef USE_MIMALLOC
#define USE_MIMALLOC 1
#endif
#include "irr_projection.hpp"
#include "hex_hash.hpp"
#include "hexagon_native_cache_release.hpp"
#include <atomic>
#include <fstream>
#include <iostream>
#include <sstream>

// Load one real FEC recurrence, keep its Fp cache alive, and collect unused
// native allocator pages. No sewing, fit, or physical-symbol result is made.
struct Memory {size_t rss=0,anon=0,file=0;};
Memory memory(){
    Memory value;std::ifstream in("/proc/self/status");std::string line;
    while(std::getline(in,line)){
        std::istringstream row(line);std::string name;size_t kib;
        if(row>>name>>kib){if(name=="VmRSS:")value.rss=kib*1024;
            else if(name=="RssAnon:")value.anon=kib*1024;
            else if(name=="RssFile:")value.file=kib*1024;}
    }return value;
}
void snapshot(const char* phase){
    auto value=memory();
    std::cout<<"NATIVE_COLLECT_MEMORY phase="<<phase<<" rss_bytes="<<value.rss
             <<" anonymous_bytes="<<value.anon<<" file_bytes="<<value.file<<std::endl;
}
struct Digest {std::array<std::string,6> blocks;size_t nnz=0;};
Digest digest(const std::array<irr_physics::mat,6>& matrices){
    Digest result;
    for(size_t b=0;b<6;++b){
        hex_sha256_state state;const auto& matrix=matrices[b];
        uint64_t dimensions[2]={matrix.nrow,matrix.ncol};
        state.update(reinterpret_cast<const unsigned char*>(dimensions),sizeof(dimensions));
        for(const auto& row:matrix.rows){
            uint64_t count=row.nnz();result.nnz+=count;
            state.update(reinterpret_cast<const unsigned char*>(&count),sizeof(count));
            if(count){
                state.update(reinterpret_cast<const unsigned char*>(row.indices),count*sizeof(int32_t));
                state.update(reinterpret_cast<const unsigned char*>(row.entries),count*sizeof(ulong));
            }
        }
        result.blocks[b]=state.finish();
    }return result;
}
int main(int argc,char** argv)try{
    irr::require(argc==4||argc==5,"native_collect_probe FEC LEC OUTPUT [--purge-delay-zero]");
    bool immediate=argc==5;
    irr::require(!immediate||std::string(argv[4])=="--purge-delay-zero","Unknown native collect option");
    if(immediate){::setenv("MIMALLOC_PURGE_DELAY","0",1);mi_option_set(mi_option_purge_delay,0);}
    std::cout<<"NATIVE_COLLECT_CONFIG purge_delay_ms="<<mi_option_get(mi_option_purge_delay)
             <<" purge_delay_env="<<(std::getenv("MIMALLOC_PURGE_DELAY")?std::getenv("MIMALLOC_PURGE_DELAY"):"unset")
             <<" mimalloc_version="<<mi_version()<<" current_thread_collect=1 worker_collect_tasks_enabled=1"<<std::endl;
    Flint::set_memory_functions();snapshot("initial");
    auto fec=std::filesystem::absolute(argv[1]),lec=std::filesystem::absolute(argv[2]),out=std::filesystem::absolute(argv[3]);
    irr_physics::context c(fec,lec,out,out,17,3,1,true,1152921504606847009);
    c.cache_limit=size_t(4)*1024*1024*1024;
    auto matrix=c.modular_recurrence(17,true);
    irr::require(c.recurrence_cache.count({true,size_t(17)})==1,"FEC17 residues were not retained in4GiB cache");
    snapshot("recurrence_returned_temporary_QQ_destroyed");
    const char* sentinel_text="1972637129837612938761298376129837612983761298376129837612983761298376129837612983761298376129837612983761/827364918273649182736491827364918273649182736491827364918273649182736491827364918273649182736491827364919";
    rat_t main_sentinel(sentinel_text);
    std::vector<rat_t> worker_sentinels(c.opt->pool.get_thread_count());
    for(size_t i=0;i<worker_sentinels.size();++i)c.opt->pool.detach_task([&,i]{worker_sentinels[i]=rat_t(sentinel_text);});
    c.opt->pool.wait();
    auto before=digest(*matrix);snapshot("hash_before_collect");
    mi_collect(true);snapshot("after_main_thread_collect");
    // These tasks additionally collect whichever worker heaps execute them.
    // Thread-pool scheduling does not promise one task per distinct worker;
    // report attempts rather than claiming complete per-thread coverage.
    std::atomic<size_t> collected{0};size_t tasks=c.opt->pool.get_thread_count();
    c.opt->pool.detach_loop(size_t(0),tasks,[&](size_t){mi_collect(true);++collected;});
    c.opt->pool.wait();snapshot("after_worker_collect_tasks");
    hexagon_native_memory::release_unused_integer_caches(&c.opt->pool);
    snapshot("after_FLINT_integer_cache_release");
    irr::require(main_sentinel.get_str()==sentinel_text,"Live main-thread rational changed after FLINT cache release");
    for(auto& value:worker_sentinels)irr::require(value.get_str()==sentinel_text,"Live worker rational changed after FLINT cache release");
    // Reuse the same FLINT runtime after cleanup, including fresh big-integer
    // allocation and arithmetic, while the old rational sentinels remain live.
    rat_t fresh(sentinel_text);auto twice=fresh+fresh;
    irr::require(twice/fresh==rat_t(2),"FLINT arithmetic failed after cache release");
    auto after=digest(*matrix);snapshot("hash_after_collect");
    irr::require(before.nnz==after.nnz&&before.blocks==after.blocks,"Native collection changed retained Fp recurrence");
    auto cached=c.modular_recurrence(17,true);
    irr::require(cached.get()==matrix.get(),"Native collection lost recurrence-cache ownership");
    for(size_t b=0;b<6;++b)std::cout<<"NATIVE_COLLECT_HASH block="<<irr::specs()[b].name<<" sha256="<<before.blocks[b]<<std::endl;
    std::cout<<"NATIVE_COLLECT_PROBE PASS nnz="<<before.nnz<<" cache_bytes="<<c.cache_bytes
             <<" worker_collect_tasks="<<collected<<" fp_hash_unchanged=1 cache_pointer_unchanged=1 live_QQ_preserved=1 flint_reuse=1"<<std::endl;
}catch(const std::exception& error){std::cerr<<"NATIVE_COLLECT_ERROR "<<error.what()<<std::endl;return 1;}
