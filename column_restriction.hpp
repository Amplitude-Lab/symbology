#pragma once
#include "staged_kernel.hpp"

namespace staged_kernel {
// Restrict a replayable operator to selected coordinates, without assembling
// its columns or copying its factors. This alone makes no symmetry claim:
// callers must prove that separately solved sectors form a direct sum.
template<class Source> class ColumnRestriction {
    const Source& source;
    std::vector<I> selected, inverse;
    structured_kernel::Relations extend_relations(structured_kernel::Relations& local) const {
        structured_kernel::Relations full(source.ncol());
        for(size_t c=0;c<inverse.size();++c)if(inverse[c]<0){Vec row;row.push_back(index(c),Q(1));full.consume(row);}
        // Star unions preserve the chosen local root: every nonroot is still
        // a singleton when attached to its representative.
        for(size_t c=0;c<selected.size();++c){auto [r,x]=local.image(c);
            if(x==0){Vec row;row.push_back(selected[c],Q(1));full.consume(row);}
            else if(r!=c){Vec row;row.push_back(selected[c],Q(1));row.push_back(selected[r],-x);full.consume(row);}
        }
        return full;
    }
    template<class T> sparse_mat<T,I> lift(const sparse_mat<T,I>& k) const {
        require(k.ncol==ncol(),"restricted kernel dimensions");
        sparse_mat<T,I> full(k.nrow,source.ncol());
        for(size_t r=0;r<k.nrow;++r){full[r].reserve(k[r].nnz());for(auto [c,x]:k[r])full[r].push_back(selected[c],x);full[r].sort_indices();}
        return full;
    }
public:
    ColumnRestriction(const Source& rows,std::vector<I> coordinates):source(rows),selected(std::move(coordinates)),inverse(rows.ncol(),-1){
        for(size_t c=0;c<selected.size();++c){auto j=selected[c];require(j>=0&&size_t(j)<inverse.size()&&inverse[j]<0,"invalid coordinate restriction");inverse[j]=index(c);}
    }
    size_t ncol()const{return selected.size();}
    size_t nrow()const{return source.nrow();}
    size_t residual_blocks()const{return source.residual_blocks();}
    template<class Emit> void operator()(structured_kernel::Relations& relations,size_t limit,Emit&& emit)const{
        auto full=extend_relations(relations);
        source(full,limit,[&](Vec&& row){for(size_t j=0;j<row.nnz();++j){require(inverse[row(j)]>=0,"restricted row escaped coordinates");row(j)=inverse[row(j)];}normalize(row);emit(std::move(row));});
    }
    template<class Emit> void modular(structured_kernel::Relations& relations,const std::vector<I>& columns,
        const field_t& fp,Emit&& emit,thread_pool* pool=nullptr,size_t max_terms=SIZE_MAX,size_t min_terms=0,
        const streaming_echelon::Substitution* substitution=nullptr,std::vector<unsigned char>* seen=nullptr,unsigned sample_shift=0)const{
        auto full=extend_relations(relations);std::vector<I> mapped(source.ncol(),-1);
        for(size_t c=0;c<selected.size();++c)mapped[selected[c]]=columns[c];
        source.modular(full,mapped,fp,std::forward<Emit>(emit),pool,max_terms,min_terms,substitution,seen,sample_shift);
    }
    template<class Emit> void residual(const MM& k,const field_t& fp,size_t block,Emit&& emit,thread_pool* pool)const{
        auto full=lift(k);source.residual(full,fp,block,std::forward<Emit>(emit),pool);
    }
    int_t integral_residual_bound(const Mat& k)const{
        if constexpr(requires(const int_t& h){source.integral_residual_bound_from_maximum(h);}){
            int_t h=0;for(const auto& row:k.rows)for(auto [c,x]:row){require(x.is_integer(),"restricted certificate needs integral candidates");auto a=x.num().abs();if(a>h)h=std::move(a);}
            return source.integral_residual_bound_from_maximum(h);
        }else return source.integral_residual_bound(lift(k));
    }
    void lift_in_place(Mat& k)const{
        require(k.ncol==ncol(),"restricted kernel dimensions");
        for(auto& row:k.rows){for(size_t j=0;j<row.nnz();++j)row(j)=selected[row(j)];row.sort_indices();}k.ncol=source.ncol();
    }
};
}
