#include "irr_recurrence_store.hpp"
int main(int argc,char** argv)try {
    Flint::set_memory_functions();irr::require(argc==2||argc==3,"recurrence_store FEC_ROOT [NEW_TEST_DIRECTORY]");
    rref_option_t opt;opt->pool.reset(2);irr::channels_t cg(&opt->pool);
    irr::recurrence_store store(16ULL<<20,16ULL<<20),off(0,0);field_t p(FIELD_Fp,1000003),q(FIELD_Fp,1000033);
    auto a=store.rational(argv[1],2,true,cg,&opt->pool),b=store.rational(argv[1],2,true,cg,&opt->pool);
    irr::require(a.get()==b.get()&&store.qq_hits==1,"QQ cache miss on identical input");
    auto x=store.modular(argv[1],2,true,cg,p,&opt->pool),y=store.modular(argv[1],2,true,cg,p,&opt->pool);
    irr::require(x.get()==y.get(),"Fp cache missed same prime");
    auto z=store.modular(argv[1],2,true,cg,q,&opt->pool),ref=off.modular(argv[1],2,true,cg,q,&opt->pool);
    for(size_t i=0;i<6;++i)irr::require(nmhv_sew_matrices_equal((*z)[i],(*ref)[i]),"Cache changes residues");
    irr::require(z.get()!=x.get()&&store.qq_hits>=3,"Prime transition copied/reloaded QQ");
    irr::recurrence_store tiny(1,1);auto t=tiny.rational(argv[1],2,true,cg,&opt->pool);
    irr::require(!tiny.qq_bytes&&tiny.qq_bypasses==1,"Oversized QQ was cached");
    irr::recurrence_store lru(10,0);
    auto pinned=std::make_shared<const irr::matrices>();
    lru.retain<irr::matrices>(lru.qq,"pinned",pinned,8,10,lru.qq_bytes);
    lru.retain<irr::matrices>(lru.qq,"other",std::make_shared<const irr::matrices>(),8,10,lru.qq_bytes);
    irr::require(lru.qq.count("pinned")&&!lru.qq.count("other"),"Evicted a live borrower");
    pinned.reset();
    lru.retain<irr::matrices>(lru.qq,"other",std::make_shared<const irr::matrices>(),8,10,lru.qq_bytes);
    irr::require(!lru.qq.count("pinned")&&lru.qq.count("other")&&lru.evictions==1,"LRU did not reclaim released data");
    if(argc==3){
        auto dir=std::filesystem::path(argv[2]);irr::require(std::filesystem::create_directory(dir),"Require new test directory");
        auto file=dir/"identity.txt";{std::ofstream out(file);out<<"old";}
        auto before=store.file_hash(file);{std::ofstream out(file);out<<"changed";}
        irr::require(store.file_hash(file)!=before,"Input mutation retained stale identity");
    }
    store.report();std::cout<<"PASS immutable sharing, cache-off, prime transition, pinned lifetime, tiny budget\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
