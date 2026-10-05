#pragma once
#include "symrep.hpp"
#include <queue>

// Exact finite-field elimination with bounded input batches. Retain only an
// independent triangular row basis; redundant input rows never accumulate in
// the global Schur complement. Every input equation is processed.
namespace streaming_echelon {
using namespace symrep;
using ModVec=sparse_vec<ulong,I>;
using ModMat=sparse_mat<ulong,I>;
struct Substitution {std::vector<I> image;std::vector<ulong> scale;};
// Nested, deterministic equation samples. Full operator verification, not
// the sampling probability, decides acceptance of a discovered kernel.
inline bool sampled_equation(size_t id,unsigned shift,uint64_t salt){
    if(!shift)return true;uint64_t x=uint64_t(id)+salt;x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);x^=x>>31;
    return !(x&((UINT64_C(1)<<shift)-1));
}
class Reducer {
    struct Work {std::vector<ulong> dense;std::vector<unsigned char> flags;std::vector<I> touched,heap;
        explicit Work(size_t n):dense(n),flags(n){};};
    ModMat basis,pending;
    std::vector<I> pivot_columns,pivot_row;
    std::vector<unsigned char> known_zero;
    Substitution substitutions;
    std::vector<Work> work;
    field_t fp;
    rref_option* opt;
    size_t last_pruned_rank=0,zero_count=0;
    size_t row_limit,nnz_limit,pending_nnz=0,input_rows=0,input_nnz=0,batches=0,retained_nnz=0;
    void reduce(ModVec& row,Work& w)const{
        auto enqueue=[&](I c){
            if(!(w.flags[c]&1)){w.flags[c]|=1;w.touched.push_back(c);}
            if(pivot_row[c]>=0&&!(w.flags[c]&2)){
                w.flags[c]|=2;w.heap.push_back(pivot_row[c]);std::push_heap(w.heap.begin(),w.heap.end(),std::greater<I>());
            }
        };
        for(auto [c,x]:row)if(substitutions.image[c]>=0){I j=substitutions.image[c];ulong y=substitutions.scale[c]==1?x:nmod_mul(x,substitutions.scale[c],fp.mod);w.dense[j]=nmod_add(w.dense[j],y,fp.mod);enqueue(j);}
        while(!w.heap.empty()){
            std::pop_heap(w.heap.begin(),w.heap.end(),std::greater<I>());I r=w.heap.back();w.heap.pop_back();
            ulong x=w.dense[pivot_columns[r]];if(!x)continue;
            ulong shoup=n_mulmod_precomp_shoup(x,fp.mod.n);
            for(auto [c,y]:basis[r]){
                ulong value=_nmod_sub(w.dense[c],n_mulmod_shoup(x,y,shoup,fp.mod.n),fp.mod);
                if(value)enqueue(c);w.dense[c]=value;
            }
        }
        row.zero();
        for(I c:w.touched){if(w.dense[c]){require(pivot_row[c]<0,"streaming triangular reduction failed");row.push_back(c,w.dense[c]);}w.dense[c]=0;w.flags[c]=0;}
        w.touched.clear();row.sort_indices();row.compress();
    }
    void flush(){
        if(pending.rows.empty())return;
        pending.nrow=pending.rows.size();
        if(basis.nrow){opt->pool.detach_loop(size_t(0),pending.nrow,[&](size_t i){reduce(pending[i],work[thread_id()]);});opt->pool.wait();}
        std::erase_if(pending.rows,[](const auto& row){return !row.nnz();});pending.nrow=pending.rows.size();
        if(pending.nrow){
            pending.sort_rows_by_nnz();auto pivots=sparse_mat_rref_forward(pending,fp,opt);require(!opt->abort,"streaming modular kernel aborted");
            size_t first=basis.rows.size();
            for(const auto& batch:pivots)for(auto [r,c]:batch){
                require(pivot_row[c]<0,"duplicate streaming pivot");pivot_row[c]=index(basis.rows.size());pivot_columns.push_back(c);retained_nnz+=pending[r].nnz();basis.rows.push_back(std::move(pending[r]));
            }
            basis.nrow=basis.rows.size();
            for(size_t r=first;r<basis.nrow;++r){
                require(basis[r].find(pivot_columns[r])&&*basis[r].find(pivot_columns[r])==1,"streaming pivot is not normalized");
                for(auto c:basis[r].index_span())require(pivot_row[c]<0||pivot_row[c]>=index(r),"streaming basis is not triangular");
            }
        }
        pending.clear();pending=ModMat(0,basis.ncol);pending_nnz=0;++batches;
        if((batches==1||batches%16==0)&&basis.nrow!=last_pruned_rank)prune_zeros();
        if(batches==1||batches%32==0)std::cout<<"streaming_batch="<<batches<<" input_rows="<<input_rows<<" rank="<<basis.nrow<<" retained_nnz="<<retained_nnz<<std::endl;
    }
public:
    Reducer(size_t columns,const field_t& field,rref_option_t options,size_t rows=32768,size_t nonzeros=16*1024*1024)
        :basis(0,columns),pending(0,columns),pivot_row(columns,-1),known_zero(columns,0),fp(field),opt(options),row_limit(rows),nnz_limit(nonzeros){
        require(rows>0&&nonzeros>0,"streaming batch limits must be positive");
        substitutions.image.resize(columns);std::iota(substitutions.image.begin(),substitutions.image.end(),0);substitutions.scale.assign(columns,1);
        for(size_t t=0;t<opt->pool.get_thread_count();++t)work.emplace_back(columns);
    }
    void consume(ModVec&& row){if(!row.nnz())return;++input_rows;input_nnz+=row.nnz();pending_nnz+=row.nnz();pending.rows.push_back(std::move(row));
        if(pending.rows.size()>=row_limit||pending_nnz>=nnz_limit)flush();}
    // A global sparse front gives pivot selection the cheapest equations
    // across the entire operator before more expensive rows are streamed.
    void seed(ModMat&& matrix){
        require(!basis.nrow&&pending.rows.empty()&&matrix.ncol==basis.ncol,"streaming seed requires an empty compatible reducer");
        input_rows=matrix.nrow;input_nnz=matrix.nnz();pending_nnz=input_nnz;pending=std::move(matrix);flush();
    }
    size_t rank()const{return basis.nrow;}
    size_t completed_batches()const{return batches;}
    bool is_zero(size_t column)const{return known_zero[column];}
    const Substitution& substitution()const{return substitutions;}
    size_t prune_zeros(){
        if(last_pruned_rank==basis.nrow)return zero_count;
        // Propagate both singleton zeros and doubleton identifications in
        // reverse triangular order. Each rewritten equation differs by a
        // combination of retained later pivots, preserving the exact row space.
        size_t count=0;retained_nnz=0;
        for(size_t r=basis.nrow;r-->0;){auto& row=basis[r];size_t out=0;bool changed=false;
            I pivot=pivot_columns[r];for(size_t p=0;p<row.nnz();++p){I c=row(p);if(c==pivot){row(out)=c;row[out]=row[p];++out;}
                else {changed|=substitutions.image[c]!=c||substitutions.scale[c]!=1;
                    if(substitutions.image[c]>=0){row(out)=substitutions.image[c];row[out]=substitutions.scale[c]==1?row[p]:nmod_mul(row[p],substitutions.scale[c],fp.mod);++out;}}}
            // Most retained rows do not change. Sorting and reallocating them
            // again costs more than the short-relation propagation itself.
            if(changed){row.resize(out);row.sort_indices();out=0;
            for(size_t first=0;first<row.nnz();){size_t end=first+1;I c=row(first);ulong x=row[first];while(end<row.nnz()&&row(end)==c)x=nmod_add(x,row[end++],fp.mod);if(x){row(out)=c;row[out]=x;++out;}first=end;}
            row.resize(out);row.compress();}require(row.find(pivot)&&*row.find(pivot)==1,"lost streaming pivot during short-relation propagation");
            substitutions.image[pivot]=pivot;substitutions.scale[pivot]=1;
            if(out==1){known_zero[pivot]=1;substitutions.image[pivot]=-1;substitutions.scale[pivot]=0;}
            else if(out==2){size_t p=row(0)==pivot?1:0;substitutions.image[pivot]=row(p);substitutions.scale[pivot]=nmod_neg(row[p],fp.mod);}
            retained_nnz+=out;
        }
        for(auto z:known_zero)count+=z;last_pruned_rank=basis.nrow;zero_count=count;return count;
    }
    void flush_pending(){flush();}
    std::pair<ModMat,std::vector<std::vector<pivot_t<I>>>> finish(){
        flush();std::vector<std::vector<pivot_t<I>>> pivots(2);
        for(size_t r=0;r<basis.nrow;++r)pivots[1].emplace_back(index(r),pivot_columns[r]);
        work.clear();
        std::cout<<"streaming_complete rows="<<input_rows<<" input_nnz="<<input_nnz<<" rank="<<basis.nrow<<" retained_nnz="<<basis.nnz()<<std::endl;
        triangular_solver_2(basis,pivots,fp,opt);require(!opt->abort,"streaming back substitution aborted");
        return {std::move(basis),std::move(pivots)};
    }
};
}
