#include "irr_physics_pipeline.hpp"
#include "irr_bulk.hpp"

// Reuse the existing physics solver at each split. Only one stage's large
// residue cache is retained at a time; checkpoints are the inter-stage boundary.
struct modular_pipeline {
    settings base;fs::path output;
    std::shared_ptr<irr::recurrence_store> recurrence;
    std::map<std::pair<size_t,size_t>,std::unique_ptr<pipeline>> workers;
    std::set<std::string> active;
    modular_pipeline(settings value,fs::path root):base(std::move(value)),output(std::move(root)){
        size_t fp=1024;if(auto v=std::getenv("PHYSICS_CACHE_MIB"))fp=std::stoull(v);
        recurrence=std::make_shared<irr::recurrence_store>(base.qq_cache_mib*(1ULL<<20),fp*(1ULL<<20));
    }
    pipeline& worker(size_t h,size_t k){auto key=std::make_pair(h,k);
        if(!workers.count(key)){auto s=base;s.h=h;s.k=k;s.output=output/("SEW_"+std::to_string(h)+"p"+std::to_string(k))/s.helicity;
            auto run=std::make_unique<pipeline>(s);
            run->c.shared_recurrence=recurrence;
            run->parent_evaluator=[this](const fs::path& root,ulong p){return evaluate(root,p);};
            run->special_evaluator=[this](const fs::path& root,ulong p){return evaluate(root,p);};
            workers.emplace(key,std::move(run));}
        return *workers.at(key);
    }
    static void release(pipeline& run){
        for(auto& block:run.c.S)block.clear();run.c.recurrence_cache.clear();run.c.step_plans.clear();
        run.c.cache_reported.clear();run.c.exact_tail.clear();run.c.tail_words.clear();run.c.cache_bytes=0;run.c.cache_prime=0;
    }
    sparse_mat<ulong,int32_t> evaluate(const fs::path& input,ulong p){
        recurrence->select_prime(p);
        auto root=fs::weakly_canonical(input);
        const auto key=fs::weakly_canonical(root).generic_string();require(active.size()<64&&active.insert(key).second,"Cyclic stage dependency");
        struct guard{std::set<std::string>& paths;std::string key;~guard(){paths.erase(key);}} cleanup{active,key};
        if(!fs::exists(root/"stage.tsv")){
            auto& run=worker(base.h,base.k);return run.exact_source(root,p);
        }
        auto r=irr_stage::read(root/"stage.tsv");size_t h=irr_io::number(r.at("h")),k=irr_io::number(r.at("k"));
        require(h+k==base.h+base.k&&r.at("helicity")==base.helicity,"Incompatible stage weight/helicity");
        auto& run=worker(h,k);
        if(r.at("operator")!="bulk"&&r.at("operator")!="shift"){
            auto result=run.evaluate(root,p);release(run);return result;
        }
        run.bind(root,r);
        if(r.at("operator")=="bulk"){
            auto condition=root/r.at("condition"),blocks=root/r.at("condition_blocks");
            require(run.input_hash(condition)==r.at("condition_sha256")&&run.input_hash(blocks)==r.at("blocks_sha256"),"Bulk condition changed");
        }else{
            auto parent=root/r.at("parent");require(run.input_digest(parent)==r.at("parent_identity"),"Shift parent changed");
        }
        if(fs::exists(root/"qq_complete.tsv"))return run.exact_source(root,p);
        if(fs::exists(cache_file(root,p))){std::cout<<"PIPELINE_CACHE_HIT stage="<<root<<" prime="<<p<<std::endl;return read_cache(root,p,h,k);}
        sparse_mat<ulong,int32_t> result;
        if(r.at("operator")=="bulk"){
            irr_bulk::source source(h,k,root,base.fec,base.lec,root/r.at("condition"),root/r.at("condition_blocks"),r.at("prereduce"),run.c.opt,recurrence.get());
            result=source.evaluate(p);
        }else{
            result=evaluate(root/r.at("parent"),p);
            irr_run::phase timer("shift_prime_"+std::to_string(p));
            size_t old_h=irr_io::number(r.at("source_h")),old_k=irr_io::number(r.at("source_k"));
            require(old_h>h&&old_h+old_k==h+k,"Invalid left shift");
            field_t F(FIELD_Fp,p);irr::seam_recoupling_t recouple(run.c.cg,run.c.opt);
            for(size_t current=old_h;current>h;--current)result=irr_shift_cache::shift_one(result,current,h+k-current,base.fec,base.lec,run.c.cg,recouple,F,run.c.opt,recurrence.get());
        }
        save_cache(root,result,p,h,k);
        irr_stage::fixed_text(root/"cache_shape.tsv",std::to_string(result.nrow)+"\t"+std::to_string(result.ncol)+"\n");
        std::cout<<"PIPELINE_STAGE operator="<<r.at("operator")<<" h="<<h<<" k="<<k<<" prime="<<p<<" rows="<<result.nrow<<" nnz="<<result.nnz()<<std::endl;
        return result;
    }
};

int main(int argc,char** argv)try {
    Flint::set_memory_functions();
    require(argc>=5,"irr_pipeline MHV|NMHV H K OUTPUT_ROOT [--input QQ_OR_STAGE] [--plan TSV | --conditions collinear+mrk] [--cache-only --primes N] [--export-all] --fec-root DIR --lec-root DIR --data DIR");
    settings s;s.helicity=argv[1];s.h=irr_io::number(argv[2]);s.k=irr_io::number(argv[3]);
    fs::path output=fs::absolute(argv[4]);s.fec="output/FEC";s.lec=fs::path("output/LEC")/s.helicity;s.data="data";
    fs::path condition,blocks;std::string prereduce=s.k==1?"right":"left";bool cache_only=false;
    for(int i=5;i<argc;++i){std::string key=argv[i];
        if(key=="--cache-only"){cache_only=true;continue;}if(key=="--export-all"){s.all=true;continue;}
        require(i+1<argc,"Missing option value");std::string value=argv[++i];
        if(key=="--input")s.input=fs::absolute(value);else if(key=="--plan")s.plan=fs::absolute(value);else if(key=="--conditions")s.conditions=value;
        else if(key=="--fec-root")s.fec=fs::absolute(value);else if(key=="--lec-root")s.lec=fs::absolute(value);else if(key=="--data")s.data=fs::absolute(value);
        else if(key=="--condition")condition=fs::absolute(value);else if(key=="--condition-blocks")blocks=fs::absolute(value);else if(key=="--prereduce")prereduce=value;
        else if(key=="--samples")s.samples=irr_io::number(value);else if(key=="--validation-samples")s.checks=irr_io::number(value);
        else if(key=="--batch")s.batch=irr_io::number(value);else if(key=="--seed")s.seed=irr_io::number(value);
        else if(key=="--projection")s.projection=value;
        else if(key=="--qq-cache-mib")s.qq_cache_mib=irr_io::number(value);
        else if(key=="--primes"){s.primes=value=="auto"?0:irr_io::number(value);require(value=="auto"||s.primes,"Prime count must be positive or auto");}
        else throw std::runtime_error("Unknown option "+key);
    }
    require((s.helicity=="MHV"||s.helicity=="NMHV")&&s.h&&s.k&&(s.h+s.k)%2==0&&s.batch&&s.checks,"Invalid helicity/weights/batch");
    require(s.plan.empty()||s.conditions.empty(),"Choose plan or conditions");require(!cache_only||s.primes,"Cache-only needs a finite prime count");
    require(s.input.empty()||!s.plan.empty()||!s.conditions.empty(),"Input source requires a condition or shift plan; it is never overwritten");
    require(prereduce=="none"||prereduce=="left"||prereduce=="right","Invalid prereduce direction");
    require(s.projection=="scalar"||s.projection=="tail-exact","Choose scalar or tail-exact projection");
    require(s.qq_cache_mib<=(SIZE_MAX>>20),"QQ cache budget overflows bytes");
    if(condition.empty())condition=s.data/"dlogmat_full_ca_pd3.wxf";if(blocks.empty())blocks=s.data/"dlogmat_full_ca_pd3_blocks.tsv";
    modular_pipeline engine(s,output);auto& initial=engine.worker(s.h,s.k);
    fs::path parent=s.input;size_t h=s.h,k=s.k;std::string history;std::vector<fs::path> stages;
    if(parent.empty()){
        parent=initial.s.output/"bulk";
        fields r={{"schema","irr-bulk-pipeline-v1"},{"operator","bulk"},{"history",""},{"helicity",s.helicity},{"h",std::to_string(h)},{"k",std::to_string(k)},
            {"object_kind","homogeneous"},{"basis_convention","kernel-free-coordinate-v1"},{"prereduce",prereduce},
            {"condition",fs::relative(condition,parent).generic_string()},{"condition_blocks",fs::relative(blocks,parent).generic_string()},
            {"condition_sha256",hex_sha256_file(condition)},{"blocks_sha256",hex_sha256_file(blocks)}};
        initial.bind(parent,r);stages.push_back(parent);
    }else{
        history=irr_stage::history(parent);
        // A supplied source is interpreted at the explicit H,K; do not guess its layout.
        if(fs::exists(parent/"stage.tsv")){auto r=irr_stage::read(parent/"stage.tsv");require(irr_io::number(r.at("h"))==h&&irr_io::number(r.at("k"))==k,"Input split mismatch");}
    }
    std::vector<std::tuple<std::string,std::string,std::string>> plan;
    if(!s.plan.empty()){std::ifstream in(s.plan);require(bool(in),"Cannot read plan");std::string line;
        while(std::getline(in,line)){std::istringstream row(line);std::string label,op,value;if(!(row>>label)||label[0]=='#')continue;require(bool(row>>op>>value),"Plan requires label operator config-or-target-head");plan.emplace_back(label,op,value);}}
    else if(!s.conditions.empty()){std::istringstream in(s.conditions);std::string name;while(std::getline(in,name,'+'))plan.emplace_back(name,name,"-");}
    for(auto [label,op,config]:plan){
        if(op=="shift"){
            size_t target=irr_io::number(config);require(target&&target<h,"Shift target must be smaller and positive");size_t tail=h+k-target;
            auto& run=engine.worker(target,tail);auto root=run.s.output/(history.empty()?"bulk":history);
            auto old=fs::exists(parent/"stage.tsv")?irr_stage::read(parent/"stage.tsv"):fields{};
            std::string kind=old.count("object_kind")?old.at("object_kind"):"affine";
            if(old.empty()&&!fs::exists(parent/"object.tsv"))kind="homogeneous";
            else if(old.empty())kind=irr_stage::read(parent/"object.tsv").at("kind");
            fields r={{"schema","irr-shift-pipeline-v1"},{"operator","shift"},{"history",history},{"helicity",s.helicity},{"h",std::to_string(target)},{"k",std::to_string(tail)},
                {"source_h",std::to_string(h)},{"source_k",std::to_string(k)},{"object_kind",kind},{"parent",fs::relative(parent,root).generic_string()},{"parent_identity",run.input_digest(parent)}};
            run.bind(root,r);parent=root;h=target;k=tail;stages.push_back(root);
        }else{
            auto& run=engine.worker(h,k);
            auto operation=parse_operation(run.s,label,op,config=="-"?fs::path{}:s.plan.parent_path()/config);
            run.preflight(operation);
            history=irr_stage::append(history,op);auto root=run.s.output/history;
            require(fs::weakly_canonical(root)!=fs::weakly_canonical(parent),"Output equals parent");
            run.bind(root,run.recipe(parent,operation,history));parent=root;stages.push_back(root);
        }
    }
    require(!parent.empty()&&(!stages.empty()||s.input.empty()),"Empty input plan; source is read-only");
    if(cache_only){ulong p=n_nextprime(1ULL<<60,0);for(size_t i=0;i<s.primes;++i){auto result=engine.evaluate(parent,p);p=n_nextprime(p,0);}
        std::cout<<"PIPELINE CACHE COMPLETE primes="<<s.primes<<std::endl;}
    else if(s.all){for(const auto& root:stages){auto r=irr_stage::read(root/"stage.tsv");engine.worker(irr_io::number(r.at("h")),irr_io::number(r.at("k"))).reconstruct(root);}}
    else engine.worker(h,k).reconstruct(parent);
    return 0;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}
