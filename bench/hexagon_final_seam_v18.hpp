#pragma once
#include <set>
#include <chrono>
#include <fcntl.h>
#include <unistd.h>

// A plan contains only factored source-row addresses. Every prime regenerates
// the operator coefficients and affine RHS. It never contains field values.
namespace hexagon_final_seam_v18 {
struct Plan {
    size_t columns=0, source_rows=0, initial_rows=0, final_rows=0;
    ulong pilot=0;
    std::vector<size_t> selected;
};
inline void sync_file(const fs::path& path) {
    int fd=::open(path.c_str(),O_RDONLY);
    require(fd>=0,"V18 cannot open durable plan");
    int status=::fsync(fd);::close(fd);require(status==0,"V18 cannot sync plan");
}
inline void save(const fs::path& path,const std::string& identity,const Plan& p) {
    irr_stage::no_links(path);require(!fs::exists(path),"V18 plan already exists");
    fs::create_directories(path.parent_path());auto temp=fs::path(path.string()+".tmp");
    auto side=fs::path(path.string()+".sha256"),side_temp=fs::path(side.string()+".tmp");
    require(!fs::exists(temp)&&!fs::exists(side)&&!fs::exists(side_temp),"V18 incomplete plan retained; use a fresh output");
    {
        std::ofstream out(temp);
        out<<"hexagon-final-seam-row-plan-v18\n"<<identity<<'\n'
           <<p.columns<<' '<<p.source_rows<<' '<<p.initial_rows<<' '<<p.final_rows<<' '<<p.pilot<<' '<<p.selected.size()<<'\n';
        for(auto id:p.selected)out<<id<<'\n';out.flush();require(bool(out),"V18 plan write failed");
    }
    sync_file(temp);irr_stage::fixed_text(side_temp,hex_sha256_file(temp)+"\n");sync_file(side_temp);
    fs::rename(temp,path);fs::rename(side_temp,side);sync_file(path.parent_path());
}
inline bool load(const fs::path& path,const std::string& identity,Plan& p) {
    irr_stage::no_links(path);auto side=fs::path(path.string()+".sha256");
    if(!fs::exists(path)&&!fs::exists(side))return false;
    require(fs::exists(path)&&fs::exists(side),"V18 incomplete row plan");
    require(irr_stage::text(side)==hex_sha256_file(path)+"\n","V18 row plan checksum mismatch");
    std::ifstream in(path);std::string schema,key;std::getline(in,schema);std::getline(in,key);
    require(schema=="hexagon-final-seam-row-plan-v18"&&key==identity,"V18 row plan identity mismatch");
    size_t count=0;in>>p.columns>>p.source_rows>>p.initial_rows>>p.final_rows>>p.pilot>>count;
    require(bool(in)&&p.columns<=INT32_MAX&&p.source_rows<=INT32_MAX&&p.initial_rows&&p.initial_rows<=INT32_MAX
        &&p.final_rows&&p.final_rows<=p.initial_rows&&count==p.initial_rows-p.final_rows,"V18 row plan shape/rank mismatch");
    p.selected.resize(count);size_t previous=0;
    for(size_t i=0;i<count;++i){in>>p.selected[i];require(bool(in)&&p.selected[i]<p.source_rows
        &&(!i||p.selected[i]>previous),"V18 row plan invalid address");previous=p.selected[i];}
    std::string extra;require(!(in>>extra),"V18 row plan trailing data");return true;
}
inline size_t rank(const affine_rows& reduced) {
    return std::count_if(reduced.pivots.begin(),reduced.pivots.end(),[](const auto& p){return bool(p);});
}
template<class Source> void verify(const Source& rows,const MM& current,thread_pool* pool) {
    auto kt=current.transpose();size_t checked=0;
    for(size_t block=0;block<rows.residual_blocks();++block){
        rows.residual_transpose_direct(kt,current.nrow,block,[&](auto&& row){
            ++checked;require(!row.nnz(),"V18 complete final sewing residual failed");
        },pool);
    }
    require(checked==rows.nrow(),"V18 incomplete final sewing sweep");
    std::cout<<"V18_SEAM_VERIFIED full_operator=1 rows="<<checked<<" directions="<<current.nrow-1<<std::endl;
}
template<class Source> MM finish(Source& rows,MM current,const field_t& F,rref_option_t opt,
                                const fs::path& directory,const std::string& identity) {
    const size_t initial=current.nrow;
    auto path=directory/("candidate_"+std::to_string(initial)+".tsv");
    fs::create_directories(directory);irr_stage::lock guard(directory/".active.lock");
    Plan plan;bool replay=load(path,identity,plan);size_t contractions=0,lifts=0;
    auto started=std::chrono::steady_clock::now();
    if(replay){
        require(plan.columns==current.ncol&&plan.source_rows==rows.nrow()&&plan.initial_rows==initial,
                "V18 row plan operator/candidate mismatch");
        affine_rows reduced(initial-1,F,&opt->pool);auto kt=current.transpose();
        // Keep bounded input batches and the independent source order.
        for(size_t begin=0;begin<plan.selected.size();begin+=128){
            std::vector<size_t> ids(plan.selected.begin()+begin,plan.selected.begin()+std::min(begin+128,plan.selected.size()));
            equations eq;
            rows.selected_residual_transpose(kt,current.nrow,ids,[&](size_t id,auto&& row){
                ++contractions;nmhv_physical_equation_t e;e.condition="V18_replay";e.sample=id;e.rhs=0;
                e.coefficients.resize(current.nrow);for(auto[j,x]:row)e.coefficients[j]=x;eq.push_back(std::move(e));
            },&opt->pool);reduced.add(eq);
        }
        require(rank(reduced)==initial-plan.final_rows,"V18 replay rank changed; no result published");
        if(rank(reduced)){current=affine_product(current,reduced.solve(),F,&opt->pool);++lifts;}
    }else{
        plan.columns=current.ncol;plan.source_rows=rows.nrow();plan.initial_rows=initial;plan.pilot=F.mod.n;
        auto kt=current.transpose();std::unique_ptr<affine_rows> reduced=
            std::make_unique<affine_rows>(current.nrow-1,F,&opt->pool,&plan.selected);
        // All old rows remain satisfied under a subsequent affine restriction.
        // Continue at the next block instead of rebuilding the satisfied prefix.
        for(size_t block=0;block<rows.residual_blocks()&&current.nrow>1;++block){
            equations eq;size_t id=block*128;
            rows.residual_transpose_direct(kt,current.nrow,block,[&](auto&& row){
                size_t address=id++;++contractions;if(!row.nnz())return;
                nmhv_physical_equation_t e;e.condition="V18_learning";e.sample=address;e.rhs=0;
                e.coefficients.resize(current.nrow);for(auto[j,x]:row)e.coefficients[j]=x;eq.push_back(std::move(e));
            },&opt->pool);reduced->add(eq);
            auto gained=rank(*reduced);
            if(gained&&(gained>=64||block+1==rows.residual_blocks())){
                current=affine_product(current,reduced->solve(),F,&opt->pool);++lifts;
                kt=current.transpose();reduced=std::make_unique<affine_rows>(current.nrow-1,F,&opt->pool,&plan.selected);
                std::cout<<"V18_SEAM_RESTRICT block="<<block<<" rank="<<gained<<" remaining="<<current.nrow-1<<std::endl;
            }
        }
        plan.final_rows=current.nrow;std::sort(plan.selected.begin(),plan.selected.end());
        require(std::adjacent_find(plan.selected.begin(),plan.selected.end())==plan.selected.end()
            &&plan.selected.size()==initial-current.nrow,"V18 source witness rank mismatch");
    }
    auto before_verify=std::chrono::steady_clock::now();verify(rows,current,&opt->pool);
    auto finished=std::chrono::steady_clock::now();
    if(!replay)save(path,identity,plan);
    std::cout<<"V18_SEAM_COMPLETE replay="<<replay<<" selected_rows="<<plan.selected.size()
        <<" contractions="<<contractions<<" lifts="<<lifts<<" initial_directions="<<initial-1
        <<" final_directions="<<current.nrow-1<<" solve_seconds="
        <<std::chrono::duration<double>(before_verify-started).count()<<" verify_seconds="
        <<std::chrono::duration<double>(finished-before_verify).count()<<std::endl;
    return current;
}
} // namespace hexagon_final_seam_v18
