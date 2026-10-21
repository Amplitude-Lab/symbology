#include "irr_physics_pipeline.hpp"

void self_test() {
    field_t F(FIELD_Fp,1000003);affine_rows rows(2,F);rows.add({{"a",0,{1,2,0},3}});rows.add({{"b",0,{0,1,0},1}});auto x=rows.solve();require(x.origin==std::vector<ulong>{1,1}&&!x.directions.nrow,"Affine toy");
    affine_rows reversed(2,F);reversed.add({{"b",0,{0,1,0},1},{"a",0,{1,2,0},3}});require(reversed.solve().origin==x.origin,"Order toy");
    affine_rows shifted(1,F);shifted.add({{"origin",0,{2,5},9}});require(shifted.solve().origin==std::vector<ulong>{2},"Origin subtraction toy");
    affine_rows unique(0,F);unique.add({{"check",0,{3},3}});bool failed=false;try{unique.add({{"bad",0,{3},4}});}catch(...){failed=true;}require(failed,"Unique candidate must still be checked");
    // Exercise missing/out-of-order pivots, dependent rows and the parallel
    // backward pass against the original affine solver.
    rref_option_t opt;opt->pool.reset(4);equations test;
    for(size_t i=0;i<260;++i){size_t j=(i*137)%260;std::vector<ulong> r(400);
        r[j]=1;for(size_t k=260;k<400;++k)r[k]=(j*19+k*7)%101;
        ulong rhs=0;for(size_t k=0;k<400;++k)rhs=nmod_add(rhs,nmod_mul(r[k],k+1,F.mod),F.mod);
        test.push_back({"regression",i,std::move(r),rhs});}
    test.push_back(test[10]);affine_rows streamed(400,F,&opt->pool);
    for(size_t i=0;i<test.size();i+=32)streamed.add(equations(test.begin()+i,test.begin()+std::min(i+32,test.size())),false);
    auto actual=streamed.solve();
    auto expected=nmhv_physical_solve_affine_space<int32_t>(test,400,F);
    require(actual.origin==expected.origin&&actual.pivot_columns==expected.pivot_columns&&nmhv_sew_matrices_equal(actual.directions,expected.directions),"Affine regression");
    require(irr_stage::append("collinear","mrk")=="collinear_mrk"&&irr_stage::append("","mrk")=="mrk","Stage paths");std::cout<<"AFFINE AND PATH TESTS PASS\n";
    equations ope{{"ope",0,{0,1,0},3,{4,2,2,0,0}}, {"ope",1,{1,0,0},2,{1,1,1,0,0}}, {"ope",2,{2,0,0},4,{2,2,-2,0,0}}, {"ope",3,{1,1,0},5,{6,2,2,1,0}}};
    affine_rows traced(2,F);std::ostringstream report;ope_rank_trace(traced,ope,report,F.mod.n);
    affine_rows plain(2,F);plain.add(ope);auto a=traced.solve(),b=plain.solve();
    require(a.origin==b.origin&&a.pivot_columns==b.pivot_columns&&nmhv_sew_matrices_equal(a.directions,b.directions),"OPE grouped solve changed result");
    for(auto& row:ope)row.coefficients.insert(row.coefficients.end()-1,2*row.coefficients[0]);
    affine_rows partial(3,F),partial_plain(3,F);std::ostringstream partial_report;ope_rank_trace(partial,ope,partial_report,F.mod.n);partial_plain.add(ope);
    auto pa=partial.solve(),pb=partial_plain.solve();require(pa.directions.nrow==1&&pa.origin==pb.origin&&nmhv_sew_matrices_equal(pa.directions,pb.directions),"OPE nonunique affine basis changed");
    for(auto& row:ope)row.coefficients.erase(row.coefficients.end()-2);
    require(report.str().find("\t4\t1\t3\t1\t2\t0\t1")!=std::string::npos,"OPE first unique power");
    auto bad=ope;bad.back().rhs=6;affine_rows inconsistent(2,F);failed=false;
    try{std::ostringstream out;ope_rank_trace(inconsistent,bad,out,F.mod.n);}catch(...){failed=true;}require(failed,"High-power inconsistency must not be skipped");
    affine_rows zero(0,F);std::ostringstream empty;ope_rank_trace(zero,{{"ope",0,{3},3,{1,1,1,0,0}}},empty,F.mod.n);
    require(empty.str().find("\t-1\t0\t0\t0\t0\t0\t1")!=std::string::npos,"Already unique input");
    std::cout<<"OPE RANK TESTS PASS\n";
}
int main(int argc,char** argv)try {
    Flint::set_memory_functions();if(argc==2&&std::string(argv[1])=="--self-test"){self_test();return 0;}
    require(argc>=7,"irr_conditions_v2 MHV|NMHV H K INPUT OUTPUT_PARENT --conditions collinear+mrk [--plan FILE] --fec-root DIR --lec-root DIR [--data DIR] [--sampling adaptive|fixed] [--projection scalar|tail-exact] [--modular-kernel shoup|plain] [--samples N] [--batch N] [--validation-samples N] [--primes auto|N] [--export-all]");
    settings s;s.helicity=argv[1];s.h=irr_io::number(argv[2]);s.k=irr_io::number(argv[3]);s.input=fs::absolute(argv[4]);s.output=fs::absolute(argv[5]);s.data="data";s.fec="output/FEC";s.lec=fs::path("output/LEC")/s.helicity;
    for(int i=6;i<argc;++i){std::string key=argv[i];if(key=="--ope-rank-only"){s.ope_rank_only=true;continue;}if(key=="--export-all"){s.all=true;continue;}if(key=="--resume-legacy"){s.resume_legacy=true;continue;}require(i+1<argc,"Missing option value");std::string v=argv[++i];
        if(key=="--conditions")s.conditions=v;else if(key=="--plan")s.plan=v;else if(key=="--fec-root")s.fec=v;else if(key=="--lec-root")s.lec=v;else if(key=="--data")s.data=v;
        else if(key=="--samples")s.samples=irr_io::number(v);else if(key=="--batch")s.batch=irr_io::number(v);else if(key=="--validation-samples")s.checks=irr_io::number(v);
        else if(key=="--resume-shifted-cache")s.shifted_cache=fs::absolute(v);
        else if(key=="--sampling")s.sampling=v;else if(key=="--projection")s.projection=v;
        else if(key=="--modular-kernel"){require(v=="plain"||v=="shoup","Choose plain or shoup");s.shoup=v=="shoup";}
        else if(key=="--seed")s.seed=irr_io::number(v);else if(key=="--primes"){s.primes=v=="auto"?0:irr_io::number(v);require(v=="auto"||s.primes,"Prime budget must be positive or auto");}else throw std::runtime_error("Unknown option "+key);
    }
    require((s.helicity=="MHV"||s.helicity=="NMHV")&&s.h&&s.k&&(s.h+s.k)%2==0&&s.batch&&s.checks,"Invalid helicity/weights/batch");
    require(s.conditions.empty()||s.plan.empty(),"Choose --conditions or --plan, not both");
    require(s.sampling=="fixed"||s.sampling=="adaptive","Choose fixed or adaptive sampling");
    require(s.projection=="scalar"||s.projection=="tail-exact","Choose scalar or tail-exact projection");
    require(s.projection!="tail-exact"||s.k<=3,"Tail-exact requires k<=3");
    if(!s.samples&&s.sampling=="fixed"){size_t n=0;if(fs::exists(s.input/"stage.tsv"))s.samples=irr_io::number(irr_stage::read(s.input/"stage.tsv").at("samples"));
        else{if(fs::exists(s.input/"object.tsv"))n=irr_io::number(irr_stage::read(s.input/"object.tsv").at("directions"));else{rref_option_t opt;auto t=irr::read_tensor(s.input/(irr_io::stem(s.h,s.k,0)+".wxf"),&opt->pool);n=t.dim(0);}s.samples=std::max<size_t>(384,2*n+32);}}
    auto ops=operations(s);require(s.shifted_cache.empty()||(ops.size()==1&&ops[0].name=="ope"),"Shifted-cache resume supports a single OPE stage");pipeline run(s);auto parent=s.input;std::string history=irr_stage::history(parent);std::vector<fs::path> stages;
    for(const auto& op:ops){run.preflight(op);history=irr_stage::append(history,op.name);auto root=s.output/history;require(fs::weakly_canonical(root)!=fs::weakly_canonical(parent),"Output equals input");run.bind(root,run.resume_recipe(root,run.recipe(parent,op,history)));stages.push_back(root);parent=root;}
    run.prepare_shifted(stages.back());
    if(s.ope_rank_only){
        require(ops.size()==1&&(ops[0].name=="ope"||ops[0].name=="ope_small_s")&&s.primes>0&&s.shifted_cache.empty(),"Rank-only requires one OPE stage, explicit --primes N and no imported OPE cache");
        ulong p=n_nextprime(1ULL<<60,0);for(size_t i=0;i<s.primes;++i){run.observe_ope(stages.back(),p);p=n_nextprime(p,0);}return 0;
    }
    if(s.all)for(auto root:stages)run.reconstruct(root);else run.reconstruct(stages.back());
    return 0;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}
