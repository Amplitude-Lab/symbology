#include "reduced_d3/irr_projection.hpp"
int main(int argc,char** argv)try{
    irr::require(argc==5,"mrk_rhs_batch_probe FEC LEC DATA DIAGNOSTICS");
    Flint::set_memory_functions();
    irr_physics::context c(argv[1],argv[2],argv[4],argv[3],7,3,1,true,101);
    size_t cases=0;
    for(ulong prime:{ulong(101),ulong(1152921504606847009)})for(uint64_t seed:{uint64_t(0x193a4527),uint64_t(977)}){
        c.F=field_t(FIELD_Fp,prime);c.seed=seed;
        for(size_t degree:{size_t(0),size_t(3),size_t(9)})for(auto [offset,count]:std::vector<std::pair<size_t,size_t>>{{0,1},{2,1},{1,1},{3,1},{7,2},{0,17},{131,1},{0,1},{0,129}}){
            c.samples=count;c.sample_offset=offset;auto value=irr_physics::mrk_rhs_reuse_batch(c,degree);
            irr::require(value.size()==count,"MRK batch size changed");
            irr::require(c.mrk_rhs_cache.size()<=128,"MRK cache exceeds its bound");
            for(size_t i=0;i<count;++i){auto slots=nmhv_physical_sample_mrk_slots(seed,offset+i,9,c.F);
                auto reference=nmhv_physical_mrk_rhs(5,slots,irr_physics::symbols(c),degree,c.F);
                irr::require(value[i]==reference,"MRK batch/cache/reference mismatch");++cases;
            }
        }
    }
    std::cout<<"MRK_BATCH_PROBE PASS cases="<<cases<<" prime_seed_degree_offset_cache_checked=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"MRK_BATCH_PROBE_ERROR "<<e.what()<<std::endl;return 1;}
