#pragma once
#include "symrep.hpp"

// Generic triangular affine map x_p + sum_j U_pj x_j = b_p.
// Free coordinates and a single shared affine origin are explicit. No physics,
// expected rank, or prime-specific numerical values appear in its structure.
namespace sparse_affine {
using I=int32_t;
using V=sparse_vec<ulong,I>;
using M=sparse_mat<ulong,I>;
struct Factor {
    field_t F;
    size_t n;
    std::vector<V> upper;
    std::vector<ulong> rhs;
    std::vector<I> free,position;
    std::vector<unsigned char> pivot;
    template<class Rows> Factor(const Rows& rows,const field_t& f):F(f),n(rows.size()),upper(n),rhs(n),position(n,-1),pivot(n){
        for(size_t j=0;j<n;++j)if(rows[j]){
            const auto& row=*rows[j];symrep::require(row.size()==n+1&&row[j]==1,"Affine factor pivot shape/normalization");
            for(size_t k=0;k<j;++k)symrep::require(!row[k],"Affine factor requires triangular equations");
            pivot[j]=1;rhs[j]=row[n];
            for(size_t k=j+1;k<n;++k)if(row[k])upper[j].push_back(I(k),row[k]);
        }else{position[j]=I(free.size());free.push_back(I(j));}
    }
    // Dual application: evaluate an arbitrary operator on this affine map.
    // The last output coordinate is its constant value, NOT an extra direction.
    std::vector<ulong> project(std::vector<ulong> row)const{
        symrep::require(row.size()==n+1,"Affine factor projection width");
        for(size_t j=0;j<n;++j)if(pivot[j]&&row[j]){
            ulong x=row[j];row[n]=nmod_add(row[n],nmod_mul(x,rhs[j],F.mod),F.mod);row[j]=0;
            for(auto[k,v]:upper[j])row[k]=nmod_sub(row[k],nmod_mul(x,v,F.mod),F.mod);
        }
        std::vector<ulong> result;result.reserve(free.size()+1);
        for(auto j:free)result.push_back(row[j]);result.push_back(row[n]);return result;
    }
    std::vector<ulong> lift(const std::vector<ulong>& parameters,bool affine=true)const{
        symrep::require(parameters.size()==free.size(),"Affine factor lift width");
        std::vector<ulong> x(n);for(size_t i=0;i<free.size();++i)x[free[i]]=parameters[i];
        for(size_t j=n;j-->0;)if(pivot[j]){
            ulong value=affine?rhs[j]:0;
            for(auto[k,v]:upper[j])value=nmod_sub(value,nmod_mul(v,x[k],F.mod),F.mod);
            x[j]=value;
        }return x;
    }
    // Compose dependency factors in free-coordinate space, never expand an
    // n-by-n reduced matrix. The output itself is unavoidable at the legacy
    // checkpoint boundary. One bounded accumulator is reused for all rows.
    M materialize()const{
        size_t d=free.size();std::vector<V> maps(n);
        std::vector<ulong> sum(d+1);std::vector<unsigned char> used(d+1);std::vector<I> touched;
        auto add=[&](I j,ulong x){if(!x)return;if(!used[j]){used[j]=1;touched.push_back(j);}sum[j]=nmod_add(sum[j],x,F.mod);};
        for(size_t j=n;j-->0;){
            if(!pivot[j]){maps[j].push_back(position[j],1);continue;}
            add(I(d),rhs[j]);
            for(auto[k,v]:upper[j]){
                ulong scale=nmod_neg(v,F.mod),pre=n_mulmod_precomp_shoup(scale,F.mod.n);
                for(auto[q,x]:maps[k])add(q,n_mulmod_shoup(scale,x,pre,F.mod.n));
            }
            std::sort(touched.begin(),touched.end());
            for(auto q:touched){if(sum[q])maps[j].push_back(q,sum[q]);sum[q]=0;used[q]=0;}touched.clear();
        }
        M out(d+1,n);std::vector<size_t> counts(d+1);
        for(const auto& row:maps)for(auto[q,x]:row)++counts[q];
        for(size_t q=0;q<=d;++q)out[q].reserve(counts[q]);
        for(size_t j=0;j<n;++j){for(auto[q,x]:maps[j])out[q].push_back(I(j),x);maps[j].clear();}
        return out;
    }
};
}
