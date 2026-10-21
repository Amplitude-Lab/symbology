// Read-only validation of a unique rational symbol at explicitly supplied points.
#include "irr_physics_pipeline.hpp"
int main(int argc,char**argv)try{
    irr::require(argc==13,"physical_holdouts H K QQ FEC LEC DATA OUT PRIME SEED SAMPLES MODE MULTIPLIER");
    Flint::set_memory_functions();settings s;s.helicity="MHV";s.h=std::stoul(argv[1]);s.k=std::stoul(argv[2]);s.input=fs::absolute(argv[3]);s.fec=fs::absolute(argv[4]);s.lec=fs::absolute(argv[5]);s.data=fs::absolute(argv[6]);s.output=fs::absolute(argv[7]);
    ulong p=std::stoull(argv[8]);irr::require(n_is_prime(p),"Prime required");s.seed=std::stoull(argv[9]);s.batch=std::stoul(argv[10]);s.checks=s.batch;std::string mode=argv[11];ulong multiplier=std::stoull(argv[12]);
    irr::require(mode=="all"||mode=="collinear"||mode=="mrk"||mode=="ope_small_s","Mode");
    pipeline run(s);run.c.F=field_t(FIELD_Fp,p);auto packed=run.exact_source(s.input,p);irr::require(packed.nrow==1,"Unique candidate required");
    if(multiplier!=1)for(auto& row:packed.rows)for(size_t z=0;z<row.nnz();++z)row[z]=nmod_mul(row[z],multiplier,run.c.F.mod);
    run.unpack(packed);packed.clear();run.c.seed=s.seed;run.c.sample_offset=1000003;
    size_t total=0,bad=0;
    auto check=[&](const equations& eq,const std::string& name){size_t failed=0,fresh_ope=0;std::map<int,size_t> powers;
        for(const auto& r:eq){irr::require(r.coefficients.size()==1,"One candidate coefficient");failed+=r.coefficients[0]!=r.rhs;
            if(name.starts_with("ope")){++powers[r.ope_coordinate[0]];fresh_ope+=r.ope_coordinate[0]>10;}}
        total+=eq.size();bad+=failed;std::cout<<"PHYSICAL_HOLDOUT {\"condition\":\""<<name<<"\",\"prime\":"<<p<<",\"seed\":"<<s.seed<<",\"rows\":"<<eq.size()<<",\"mismatches\":"<<failed<<",\"S_power_above_10\":"<<fresh_ope<<"}"<<std::endl;
        for(auto [power,count]:powers)std::cout<<"OPE_POWER power="<<power<<" rows="<<count<<std::endl;
    };
    for(const std::string name:{"collinear","mrk","ope_small_s"})if(mode=="all"||mode==name){
        auto op=parse_operation(s,name,name);if(name=="ope_small_s")op.coverage="MHV-small-S-T1T2-all";run.preflight(op);
        if(name=="ope_small_s")for(bool role:{false,true})check(irr_physics::ope_small_s(run.c,role,true),role?"ope_original_check_role":"ope_original_solve_role");
        else check(irr_physics::slots(run.c,name=="mrk"),name);
        run.c.cache_summary();
    }
    std::cout<<"PHYSICAL_SUMMARY {\"rows\":"<<total<<",\"mismatches\":"<<bad<<",\"multiplier\":"<<multiplier<<",\"passed\":"<<(bad?"false":"true")<<"}"<<std::endl;
    return bad?2:0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
