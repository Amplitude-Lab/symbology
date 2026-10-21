#include "irr_physics_pipeline.hpp"

int main(int argc,char** argv)try {
    require(argc==7,"ope_rank HEL PARENT EXPECTED OUTPUT FEC LEC (L6 8p4)");
    Flint::set_memory_functions();settings s;s.helicity=argv[1];s.h=8;s.k=4;s.input=argv[2];s.output=argv[4];s.fec=argv[5];s.lec=argv[6];s.data="data/ope_L10_main";
    s.conditions=s.helicity=="MHV"?"ope_small_s":"ope";pipeline run(s);auto op=operations(s).front();run.preflight(op);
    auto r=run.recipe(s.input,op,"rank_test");ulong p=n_nextprime(1ULL<<60,0);
    auto current=read_cache(s.input,p,8,4);run.unpack(current);
    affine_rows original(current.nrow-1,run.c.F),grouped(current.nrow-1,run.c.F);
    fs::create_directories(s.output);std::ofstream report(s.output/"rank.tsv");
    run.batches(r,false,[&](const equations& eq){original.add(eq);ope_rank_trace(grouped,eq,report,p);});
    auto a=original.solve(),b=grouped.solve();
    require(a.origin==b.origin&&a.pivot_columns==b.pivot_columns&&nmhv_sew_matrices_equal(a.directions,b.directions),"Grouped affine result changed");
    sparse_mat<ulong,int32_t> transform(b.directions.nrow+1,current.nrow);
    for(size_t i=0;i<b.directions.nrow;++i)transform[i]=b.directions[i];
    for(size_t j=0;j<b.origin.size();++j)if(b.origin[j])transform[b.directions.nrow].push_back(j,b.origin[j]);
    transform[b.directions.nrow].push_back(current.nrow-1,1);
    auto result=sparse_mat_mul(transform,current,run.c.F,&run.c.opt->pool);nmhv_sew_modular_canonicalize(result,&run.c.opt->pool);
    auto expected=read_cache(argv[3],p,8,4);require(nmhv_sew_matrices_equal(result,expected),"Historical sewing differs");
    run.batches(r,true,[&](const equations& eq){require(affine_satisfies(eq,b,run.c.F),"Held-out OPE failed");});
    std::cout<<"L6 OPE RANK AND HISTORICAL SEWING PASS "<<s.helicity<<std::endl;
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
