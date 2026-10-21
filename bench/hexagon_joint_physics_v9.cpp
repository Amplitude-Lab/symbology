#include "irr_physics_pipeline.hpp"
#include "irr_bulk.hpp"
#include "shared/irrep_seam_rows.hpp"
#include "shared/seam_discovery.hpp"
#include "shared/joint_checkpoint.hpp"

using irr::require;
using MM=sparse_mat<ulong,int32_t>;

MM affine_product(const MM& current,const nmhv_physical_affine_solution_t<int32_t>& child,
                  const field_t& F,thread_pool* pool) {
    MM transform(child.directions.nrow+1,current.nrow);
    for(size_t i=0;i<child.directions.nrow;++i)transform[i]=child.directions[i];
    for(size_t j=0;j<child.origin.size();++j)if(child.origin[j])
        transform[child.directions.nrow].push_back(j,child.origin[j]);
    transform[child.directions.nrow].push_back(current.nrow-1,1);
    auto result=sparse_mat_mul(transform,current,F,pool);
    nmhv_sew_modular_canonicalize(result,pool);
    return result;
}

// Transfer the candidate into the six shared-axis blocks row by row, freeing
// each original row immediately. No full second candidate copy is retained.
void consume_unpack(pipeline& run, MM& packed) {
    irr_run::phase timer("consume_unpack");run.c.rows=packed.nrow;
    std::array<size_t,7> boundary{};
    for(size_t b=0;b<6;++b){size_t cols=run.c.sew.first[b]*run.c.sew.last[b];
        boundary[b+1]=boundary[b]+cols;run.c.S[b]=MM(packed.nrow,cols);}
    require(boundary[6]==packed.ncol,"Consuming unpack column mismatch");
    for(size_t r=0;r<packed.nrow;++r){
        std::array<size_t,6> count{};size_t b=0;
        for(auto [j,x]:packed[r]){while(b<6&&size_t(j)>=boundary[b+1])++b;
            require(b<6&&size_t(j)>=boundary[b],"Consuming unpack index");++count[b];}
        for(size_t z=0;z<6;++z)run.c.S[z][r].reserve(count[z]);b=0;
        for(auto [j,x]:packed[r]){while(size_t(j)>=boundary[b+1])++b;
            run.c.S[b][r].push_back(int32_t(size_t(j)-boundary[b]),x);}
        packed[r].clear();
    }
    packed.clear();
}

// The same affine transform acts on every irrep block. Lift and pack one
// block at a time, releasing the previous candidate block after its product.
MM block_affine_product(pipeline& run,const nmhv_physical_affine_solution_t<int32_t>& child) {
    irr_run::phase timer("block_affine_product");
    MM transform(child.directions.nrow+1,run.c.rows);
    for(size_t i=0;i<child.directions.nrow;++i)transform[i]=child.directions[i];
    for(size_t j=0;j<child.origin.size();++j)if(child.origin[j])
        transform[child.directions.nrow].push_back(j,child.origin[j]);
    transform[child.directions.nrow].push_back(run.c.rows-1,1);
    size_t columns=0;for(const auto& b:run.c.S)columns+=b.ncol;
    MM result(transform.nrow,columns);size_t offset=0;
    for(size_t b=0;b<6;++b){size_t width=run.c.S[b].ncol;
        auto next=sparse_mat_mul(transform,run.c.S[b],run.c.F,&run.c.opt->pool);
        run.c.S[b].clear();
        for(size_t r=0;r<next.nrow;++r){for(auto [j,x]:next[r])
            result[r].push_back(int32_t(offset+size_t(j)),x);next[r].clear();}
        offset+=width;
    }
    nmhv_sew_modular_canonicalize(result,&run.c.opt->pool);return result;
}

MM physical_restrict(pipeline& run,MM current,const std::string& name,bool held_out=false) {
    size_t before=current.nrow-1;consume_unpack(run,current);
    fields r={{"operator",name},{"history","early_"+name},{"seed",std::to_string(run.s.seed)},
        {"sampling","adaptive"},{"projection","scalar"},{"samples","0"},
        {"validation_samples",std::to_string(run.s.checks)},{"coverage","MHV-small-S-T1T2-all"}};
    auto child=run.train(r);
    if(held_out)run.batches(r,true,[&](const equations& eq){require(affine_satisfies(eq,child,run.c.F),"Final physical holdout failed");});
    auto result=block_affine_product(run,child);
    std::cout<<"JOINT_PHYSICS operator="<<name<<" before="<<before
        <<" remaining="<<result.nrow-1<<" nnz="<<result.nnz()<<std::endl;
    return result;
}

void physical_validate(pipeline& run,const MM& current,const std::string& name) {
    run.unpack(current);
    fields r={{"operator",name},{"history","final_"+name},{"seed",std::to_string(run.s.seed)},
        {"sampling","adaptive"},{"projection","scalar"},{"samples","0"},
        {"validation_samples",std::to_string(run.s.checks)},{"coverage","MHV-small-S-T1T2-all"}};
    run.batches(r,true,[&](const equations& eq){
        for(const auto& row:eq){
            require(row.coefficients.back()==row.rhs,"Final physical representative failed: "+name);
            for(size_t j=0;j+1<row.coefficients.size();++j)
                require(!row.coefficients[j],"Final physical direction failed: "+name);
        }
    });
    for(auto& b:run.c.S)b.clear();
    std::cout<<"JOINT_PHYSICAL_CHECK operator="<<name<<" passed=1 remaining="<<current.nrow-1<<std::endl;
}

MM finish_seam(IrrepSeamRows& rows,MM current,const field_t& F,rref_option_t opt) {
    // Each complete sweep contracts every original, unsampled sewing row.
    // Restrict after a block changes rank, then replay from the beginning.
    size_t passes=0,visited=0;
    for(;;){
        ++passes;bool restart=false;auto kt=current.transpose();
        affine_rows reduced(current.nrow-1,F,&opt->pool);
        for(size_t block=0;block<rows.residual_blocks();++block){
            equations eq;
            rows.residual_transpose(kt,current.nrow,block,[&](auto&& row){
                ++visited;if(!row.nnz())return;
                nmhv_physical_equation_t e;e.condition="seam";e.sample=visited;e.rhs=0;
                e.coefficients.resize(current.nrow);for(auto [j,x]:row)e.coefficients[j]=x;
                eq.push_back(std::move(e));
            },&opt->pool);
            if(!eq.empty())reduced.add(eq);
            size_t rank=std::count_if(reduced.pivots.begin(),reduced.pivots.end(),[](auto& p){return bool(p);});
            if(rank && (rank>=64 || block+1==rows.residual_blocks())){
                auto child=reduced.solve();current=affine_product(current,child,F,&opt->pool);
                std::cout<<"JOINT_SEAM pass="<<passes<<" block="<<block<<" rank="<<rank
                    <<" remaining="<<current.nrow-1<<" visited="<<visited<<std::endl;
                restart=true;break;
            }
            if(block%128==0)std::cout<<"JOINT_SEAM_CHECK pass="<<passes<<" block="<<block
                <<"/"<<rows.residual_blocks()<<" remaining="<<current.nrow-1<<std::endl;
        }
        if(!restart){std::cout<<"JOINT_SEAM_COMPLETE full_operator=1 rows="<<rows.nrow()
            <<" passes="<<passes<<" remaining="<<current.nrow-1<<std::endl;return current;}
    }
}

int main(int argc,char**argv)try {
    require(argc==9||argc==10||argc==12,"joint_physics H K FEC LEC DATA OUTPUT TARGET cache|reconstruct [left|right] [--prime-index N]");
    Flint::set_memory_functions();settings s;s.helicity="MHV";s.h=std::stoul(argv[1]);s.k=std::stoul(argv[2]);
    s.fec=fs::absolute(argv[3]);s.lec=fs::absolute(argv[4]);s.data=fs::absolute(argv[5]);
    s.output=fs::absolute(argv[6]);size_t target=std::stoul(argv[7]);std::string action=argv[8];
    require(action=="cache"||action=="reconstruct","Invalid action");
    std::string direction=argc>=10?argv[9]:"right";
    size_t prime_index=1;
    if(argc==12){require(action=="cache"&&std::string(argv[10])=="--prime-index","Prime index is for cache mode");
        size_t used=0;prime_index=std::stoul(argv[11],&used);require(used==std::string(argv[11]).size()&&prime_index>=1&&prime_index<=64,"Prime index must be 1..64");}
    require(direction=="left"||direction=="right","Invalid contraction direction");
    pipeline run(s);auto root=s.output/"joint";
    std::string joint_identity;
    for(const auto& name:{"collinear","mrk","ope_small_s"}) {
        auto op=parse_operation(s,name,name);
        if(std::string(name)=="ope_small_s")op.coverage="MHV-small-S-T1T2-all";
        run.preflight(op);
        for(const auto& file:run.inputs(name))joint_identity+=hex_sha256_file(file)+"\n";
    }
    for(const auto& file:{s.data/"dlogmat_full_ca_pd3.wxf",s.data/"dlogmat_full_ca_pd3_blocks.tsv"})
        joint_identity+=hex_sha256_file(file)+"\n";
    fs::create_directories(s.output);irr_stage::fixed_text(s.output/"joint_inputs.txt",joint_identity);
    fields descriptor={{"schema","irr-physics-v2-candidate"},{"operator","bulk"},
        {"history","joint_seam_collinear_mrk_ope_small_s"},{"helicity","MHV"},
        {"h",std::to_string(s.h)},{"k",std::to_string(s.k)},{"object_kind","affine"},
        {"algorithm","joint-discovery-physical-affine-v9-continuous-pager"},{"prereduce",direction},{"target",std::to_string(target)},
        {"joint_input_sha256",hex_sha256_file(s.output/"joint_inputs.txt")},
        {"working_storage",std::getenv("HEXAGON_MAPPED_FILE")?"file-backed":"heap"},{"column_order","last-copy-major-within-irrep"},{"ope_coverage","MHV-small-S-T1T2-all"},
        {"executable_sha256",hex_sha256_file(fs::canonical("/proc/self/exe"))}};
    run.bind(root,descriptor);
    run.special_evaluator=[&](const fs::path& requested,ulong p){
        require(fs::canonical(requested)==fs::canonical(root),"Unexpected joint source");
        if(fs::exists(cache_file(root,p)))return read_cache(root,p,s.h,s.k);
        run.c.clear_ope_head_cache();run.c.F=field_t(FIELD_Fp,p);auto& opt=run.c.opt;auto* pool=&opt->pool;irr::channels_t cg(pool);
        irr_bulk::check_lec_convention(s.lec,cg,pool);std::ofstream discard("/dev/null");
        auto mf=irr_bulk::read_copies(s.fec,s.h,true),ml=irr_bulk::read_copies(s.lec,s.k,false);
        auto of=irr_bulk::read_copies(s.fec,s.h-1,true),ol=irr_bulk::read_copies(s.lec,s.k-1,false);
        std::unique_ptr<IrrepSeamRows> row_storage;
        {
            auto first=irr_bulk::read_chain(s.fec,s.h,true,mf,discard,pool);
            auto last=irr_bulk::read_chain(s.lec,s.k,false,ml,discard,pool);
            auto W=irr::read_tensor(s.data/"dlogmat_full_ca_pd3.wxf",pool);
            auto t=irr::make_seam_template(cg,W,nmhv_sew_read_irrep_axis(s.data/"dlogmat_full_ca_pd3_blocks.tsv"),opt);
            auto layout=irr::seam_layout(t,of,ol,mf,ml);
            row_storage=std::make_unique<IrrepSeamRows>(cg,t,layout,of,ol,first,last,run.c.F,opt,direction=="left");
        } // Release large temporary QQ recurrences before the streaming solve.
        auto& rows=*row_storage;
        auto checkpoint=[&](const std::string& phase,auto compute){
            auto stage=s.output/phase;auto r=descriptor;r["history"]=phase;
            r["integrability_status"]=phase=="seam_complete"?"complete":"partial-candidate";
            run.bind(stage,r);
            if(fs::exists(hexagon_joint_checkpoint::file(stage,p)))return hexagon_joint_checkpoint::load(stage,p,s.h,s.k);
            auto value=compute();hexagon_joint_checkpoint::save(stage,value,p,s.h,s.k);return value;
        };
        auto current=checkpoint("discovery",[&]{
            auto stack=staged_kernel::stack(rows);staged_kernel::Options options;
            options.target_nullity=target;options.max_residual_nullity=target;
            std::vector<int32_t> variable_order(rows.ncol());size_t offset=0;
            for(size_t b=0;b<6;++b){for(size_t f=0;f<mf[b];++f)for(size_t l=0;l<ml[b];++l)
                variable_order[offset+f*ml[b]+l]=int32_t(offset+l*mf[b]+f);
                offset+=mf[b]*ml[b];}
            require(offset==rows.ncol(),"Copy matrix permutation layout");
            auto value=staged_kernel::discover_candidate(stack,run.c.F,opt,options,variable_order);
            value.rows.emplace_back();++value.nrow;return value; // zero affine origin
        });
        for(const auto& name:{"collinear","mrk","ope_small_s"}) {
            auto parent=std::move(current);
            current=checkpoint("early_"+std::string(name),[&]{return physical_restrict(run,std::move(parent),name);});
        }
        {
            auto parent=std::move(current);
            current=checkpoint("seam_complete",[&]{return finish_seam(rows,std::move(parent),run.c.F,opt);});
        }
        // Final fresh physical checks apply on the fully integrable space.
        for(const auto& name:{"collinear","mrk","ope_small_s"})
            physical_validate(run,current,name);
        save_cache(root,current,p,s.h,s.k);
        std::cout<<"JOINT_PRIME_COMPLETE prime="<<p<<" directions="<<current.nrow-1
            <<" columns="<<current.ncol<<" nnz="<<current.nnz()<<std::endl;
        return current;
    };
    if(action=="cache"){ulong prime=n_nextprime(1ULL<<60,0);for(size_t i=1;i<prime_index;++i)prime=n_nextprime(prime,0);
        std::cout<<"JOINT_SELECTED_PRIME index="<<prime_index<<" prime="<<prime<<std::endl;run.evaluate(root,prime);}
    else run.reconstruct(root);
}catch(const std::exception& e){std::cerr<<"JOINT_ERROR "<<e.what()<<std::endl;return 1;}
