#pragma once

// Preserve the complete affine boundary operator while selecting a basis of
// its suffix rows. RHS values participate in the rank test, so a condition
// supported only on the affine RHS cannot disappear. This uses no amplitude.
namespace irr_physics {
inline std::vector<size_t> independent_collinear_suffix(context& c,const axis& last,
    const equations& raw,size_t samples){
    const size_t count=c.tail_words.size(),width=last.dim+samples;
    require(count&&raw.size()==samples*count,"Exact suffix compression layout");
    const size_t q=size_t(1)<<(2*c.k);const auto& F=c.F;
    std::vector<std::vector<ulong>> source(count,std::vector<ulong>(width));
    for(size_t beta=0;beta<count;++beta){
        for(size_t i=0;i<last.dim;++i)source[beta][i]=c.exact_tail[i*q+c.tail_words[beta]];
        for(size_t s=0;s<samples;++s)source[beta][last.dim+s]=raw[s*count+beta].rhs;
    }
    struct Pivot{std::vector<ulong> row,certificate;};
    std::vector<std::optional<Pivot>> pivots(width);std::vector<size_t> selected;
    size_t audited=0;
    for(size_t beta=0;beta<count;++beta){
        auto r=source[beta];std::vector<ulong> certificate(count);certificate[beta]=1;
        for(size_t j=0;j<width;++j)if(r[j]&&pivots[j]){
            auto x=r[j];const auto& p=*pivots[j];
            for(size_t k=j;k<width;++k)r[k]=nmod_sub(r[k],nmod_mul(x,p.row[k],F.mod),F.mod);
            for(size_t k=0;k<count;++k)certificate[k]=nmod_sub(certificate[k],nmod_mul(x,p.certificate[k],F.mod),F.mod);
        }
        // Independently replay the explicit source combination, including all
        // RHS columns. Dependent rows thus have an exact zero certificate.
        for(size_t k=0;k<width;++k){ulong v=0;
            for(size_t i=0;i<count;++i)v=nmod_add(v,nmod_mul(certificate[i],source[i][k],F.mod),F.mod);
            require(v==r[k],"Exact suffix row-span certificate");++audited;
        }
        size_t j=0;while(j<width&&!r[j])++j;
        if(j==width)continue;
        ulong inv=n_invmod(r[j],F.mod.n);
        for(auto& x:r)x=nmod_mul(x,inv,F.mod);for(auto& x:certificate)x=nmod_mul(x,inv,F.mod);
        pivots[j]=Pivot{std::move(r),std::move(certificate)};selected.push_back(beta);
    }
    size_t rank=selected.size();
    // Retain a zero row when the whole operator is zero, to preserve the
    // existing adaptive stop's zero-rank batch accounting.
    if(selected.empty())selected.push_back(0);
    std::cout<<"COL_TAIL_COMPRESSION raw_words="<<count<<" selected_words="<<selected.size()
        <<" affine_rank="<<rank<<" certificate_entries="<<audited<<" rhs_columns="<<samples
        <<" complete_operator_preserved=1"<<std::endl;
    return selected;
}
}
