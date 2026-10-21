#pragma once
#include <vector>
#include <optional>
#include <stdexcept>

namespace exact_operator {
using Dense=std::vector<std::vector<ulong>>;
struct Basis {std::vector<size_t> selected;size_t certified_entries=0;};
// Compress a small operator in a caller-supplied feature space. ALL affine
// RHS/origin columns belong to that space. A zero residual is accepted only
// after independently checking its explicit combination of original rows.
inline Basis compress(const Dense& source,const field_t& F){
    Basis out;if(source.empty())return out;size_t n=source.size(),width=source[0].size();
    for(const auto& r:source){if(r.size()!=width)throw std::runtime_error("Operator compression width");
        for(auto x:r)if(x>=F.mod.n)throw std::runtime_error("Operator compression residue");}
    struct Pivot{std::vector<ulong> row,certificate;};std::vector<std::optional<Pivot>> pivots(width);
    for(size_t i=0;i<n;++i){
        if(source[i].size()!=width)throw std::runtime_error("Operator compression width");
        auto row=source[i];std::vector<ulong> certificate(n);certificate[i]=1;
        for(auto x:row)if(x>=F.mod.n)throw std::runtime_error("Operator compression residue");
        for(size_t j=0;j<width;++j)if(row[j]&&pivots[j]){
            ulong x=row[j];const auto& p=*pivots[j];
            for(size_t k=j;k<width;++k)row[k]=nmod_sub(row[k],nmod_mul(x,p.row[k],F.mod),F.mod);
            for(size_t k=0;k<n;++k)certificate[k]=nmod_sub(certificate[k],nmod_mul(x,p.certificate[k],F.mod),F.mod);
        }
        for(size_t k=0;k<width;++k){ulong x=0;
            for(size_t j=0;j<n;++j)x=nmod_add(x,nmod_mul(certificate[j],source[j][k],F.mod),F.mod);
            if(x!=row[k])throw std::runtime_error("Operator row-span certificate");++out.certified_entries;
        }
        size_t j=0;while(j<width&&!row[j])++j;if(j==width)continue;
        ulong inv=n_invmod(row[j],F.mod.n);
        for(auto& x:row)x=nmod_mul(x,inv,F.mod);for(auto& x:certificate)x=nmod_mul(x,inv,F.mod);
        pivots[j]=Pivot{std::move(row),std::move(certificate)};out.selected.push_back(i);
    }return out;
}
}
