#pragma once
#include "symrep_io.hpp"
namespace restricted_projection {
using namespace symrep;
// Project an already restricted sewn space S[n,f,l]. Eliminate its f axis
// before projecting or row-reducing the entire forward basis. No large tensor
// product is allocated; the intermediates retain only n selected combinations.
inline Tensor sew_one(const Tensor& s,const Tensor& f,const Mat& prefix,const Mat& letter,const Mat& last,thread_pool* pool){
    require(s.rank()==3&&f.rank()==3&&s.dim(1)==f.dim(0)&&s.dim(2)==last.nrow,"restricted sewing dimensions mismatch");
    require(f.dim(1)==prefix.nrow&&f.dim(2)==letter.nrow,"restricted projection dimensions mismatch");
    size_t n=s.dim(0),right=last.ncol,left=letter.ncol,old=f.dim(1);
    Mat coefficients(n*right,f.dim(0));
    for(size_t b=0;b<n;++b)for(size_t p=s.rowptr()[b];p<s.rowptr()[b+1];++p){auto i=s.index(p);
        for(auto [j,x]:last[i[1]])coefficients[b*right+j].push_back(i[0],s.val(p)*x);}
    for(auto& row:coefficients.rows)normalize(row);
    auto expanded=mul(coefficients,flatten(f),pool);coefficients.clear();
    Mat projected(n*right*left,old);
    for(size_t r=0;r<expanded.nrow;++r)for(auto [j,x]:expanded[r]){
        size_t o=j/letter.nrow,l=j%letter.nrow;
        for(auto [k,y]:letter[l])projected[r*left+k].push_back(index(o),x*y);}
    expanded.clear();for(auto& row:projected.rows)normalize(row);
    auto reduced=mul(projected,prefix,pool);projected.clear();Mat result(n,prefix.ncol*left*right);
    for(size_t b=0;b<n;++b)for(size_t r=0;r<right;++r)for(size_t l=0;l<left;++l)
        for(auto [o,x]:reduced[(b*right+r)*left+l])result[b].push_back(index((o*left+l)*right+r),x);
    for(auto& row:result.rows)normalize(row);
    return tensor(result,{n,prefix.ncol,left,right},pool);
}
// Keep a linearly dependent right projection in basis coordinates through
// the expensive forward contraction. Expand its letter coordinates only
// after the forward/prefix indices have been reduced.
inline Tensor sew_one_compact_right(const Tensor& s,const Tensor& f,const Mat& prefix,const Mat& letter,const Mat& last,rref_option_t opt){
    auto frame=row_basis(last);
    if(frame.nrow==last.ncol)return sew_one(s,f,prefix,letter,last,&opt->pool);
    Chart chart(frame,opt);auto coordinates=chart.coordinates(last);
    std::cout<<"RIGHT_PROJECTION_FRAME columns="<<last.ncol<<" rank="<<frame.nrow<<std::endl;
    auto compact=sew_one(s,f,prefix,letter,coordinates,&opt->pool);
    Mat expanded(s.dim(0),prefix.ncol*letter.ncol*last.ncol);
    for(size_t b=0;b<compact.dim(0);++b)for(size_t p=compact.rowptr()[b];p<compact.rowptr()[b+1];++p){auto c=compact.index(p);
        for(auto [j,x]:frame[c[2]])expanded[b].push_back(index((size_t(c[0])*letter.ncol+c[1])*last.ncol+j),compact.val(p)*x);}
    compact.clear();for(auto& row:expanded.rows)normalize(row);
    return tensor(expanded,{s.dim(0),prefix.ncol,letter.ncol,last.ncol},&opt->pool);
}
}
