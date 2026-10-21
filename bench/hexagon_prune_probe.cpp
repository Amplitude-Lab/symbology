#include "irr_physics_pipeline.hpp"
#include "shared/streaming_echelon.hpp"
#include "shared/streaming_echelon_legacy.hpp"
using irr::require;
using MM=sparse_mat<ulong,int32_t>;
using MV=sparse_vec<ulong,int32_t>;
MV row(std::initializer_list<std::pair<int32_t,long>> entries,const field_t& F){
    MV r;for(auto [c,x]:entries)r.push_back(c,x>=0?ulong(x):nmod_neg(ulong(-x),F.mod));return r;
}
void compare(const MM& a,const MM& b){
    require(a.nrow==b.nrow&&a.ncol==b.ncol,"Prune matrix dimensions differ");
    for(size_t r=0;r<a.nrow;++r){require(a[r].nnz()==b[r].nnz(),"Prune matrix sparsity differs");
        for(size_t j=0;j<a[r].nnz();++j)require(a[r](j)==b[r](j)&&a[r][j]==b[r][j],"Prune coefficient differs");}
}
int main()try{
    Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(2);
    opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    for(ulong prime:{ulong(101),ulong(1152921504606847009)})for(int kind=0;kind<4;++kind){
        field_t F(FIELD_Fp,prime);
        streaming_echelon::Reducer newer(24,F,opt,16,1024);
        legacy_streaming_echelon::Reducer older(24,F,opt,16,1024);
        MM seed(4,24);seed[0]=row({{0,1},{10,2},{11,3}},F);
        seed[1]=row({{1,1},{12,4},{13,5}},F);seed[2]=row({{2,1},{14,6},{15,7}},F);
        seed[3]=kind&1?row({{3,1},{16,-2}},F):row({{3,1},{16,8},{17,9}},F);
        MM old_seed=seed;newer.seed(std::move(seed));older.seed(std::move(old_seed));
        auto check=[&]{
            require(newer.prune_zeros()==older.prune_zeros(),"Zero count differs");
            require(newer.substitution().image==older.substitution().image&&newer.substitution().scale==older.substitution().scale,"Short relation differs");
            for(size_t c=0;c<24;++c)require(newer.is_zero(c)==older.is_zero(c),"Known-zero state differs");
        };check();
        auto add=[&](MV r){MV q=r;newer.consume(std::move(r));older.consume(std::move(q));};
        add(row({{10,1},{18,2},{19,3}},F));add(row({{12,1},{18,4},{19,5}},F));
        newer.flush_pending();older.flush_pending();check();
        if(kind&2){add(row({{18,1},{19,-5}},F));add(row({{11,1},{19,6}},F));
            newer.flush_pending();older.flush_pending();check();}
        if(kind==3){add(row({{19,1}},F));newer.flush_pending();older.flush_pending();check();}
        auto a=newer.finish(),b=older.finish();compare(a.first,b.first);
        require(a.second.size()==b.second.size(),"Pivot batch count differs");
        for(size_t batch=0;batch<a.second.size();++batch){
            require(a.second[batch].size()==b.second[batch].size(),"Pivot count differs");
            for(size_t j=0;j<a.second[batch].size();++j){
                auto [ar,ac]=a.second[batch][j];auto [br,bc]=b.second[batch][j];
                require(ar==br&&ac==bc,"Pivot order differs");
            }
        }
    }
    std::cout<<"PRUNE_PROBE PASS cases=8 dense_noop=1 old_short_relations=1 new_doubletons=1 zero_cascade=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
