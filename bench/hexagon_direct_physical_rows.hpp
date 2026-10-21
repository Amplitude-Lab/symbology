#pragma once
#include "reduced_d3/irr_projection.hpp"

// Dual physical equations in the original invariant coefficient space.
// Every row is A*c=b. There is one shared affine axis, and no coefficient
// basis or identity matrix is constructed. Append -b in an auxiliary column
// only when passing the row to a homogeneous streaming reducer.
namespace irr_physics {
inline size_t invariant_columns(const context& c){
    size_t n=0;for(size_t b=0;b<6;++b)n+=c.sew.first[b]*c.sew.last[b];return n;
}
inline void direct_coordinate_space(context& c){
    c.rows=invariant_columns(c);
    for(auto& block:c.S)block.clear();
}
template<class Image>
void apply_direct_invariant_columns(context& c,equations& rows,size_t stride,Image image_of){
    require(c.rows==invariant_columns(c),"Direct physical coefficient width mismatch");
    require(rows.size()==stride,"Direct physical stride mismatch");
    size_t offset=0;
    for(size_t b=0;b<6;++b){
        size_t width=c.sew.first[b]*c.sew.last[b];
        size_t batch=std::max(size_t(1),c.scratch_limit/(sizeof(ulong)*stride));
        for(size_t first=0;first<width;first+=batch){
            size_t count=std::min(batch,width-first);std::vector<ulong> images(count*stride);
            irr_parallel::rows(count,&c.opt->pool,[&](size_t j){
                image_of(b,first+j,images.data()+j*stride);
                for(size_t sample=0;sample<stride;++sample)
                    rows[sample].coefficients[offset+first+j]=images[j*stride+sample];
            });
        }
        offset+=width;
    }
    require(offset==c.rows,"Direct physical final coefficient width mismatch");
}
inline equations direct_slots(context& c,bool mrk,int metric=3){
    size_t w=c.h+c.k,degree=mrk?w-1:0,d=degree+1,stride=c.samples*d;
    std::vector<std::vector<nmhv_limit_slot_weights_t<ulong>>> weights(w,std::vector<nmhv_limit_slot_weights_t<ulong>>(c.samples));
    for(size_t s=0;s<c.samples;++s){uint64_t state=c.seed^(0xd1b54a32d192ed03ULL*(c.sample_offset+s+1));
        for(size_t j=0;j<w;++j)weights[j][s]=nmhv_limit_source_weights<ulong>(mrk?nmhv_limit_kind_t::mrk:nmhv_limit_kind_t::collinear,j+1,state,c.F);}
    auto evolve=[&](bool first){
        std::vector<ulong> prev=first?std::vector<ulong>(stride):terminal(c,mrk,stride);
        if(first)for(size_t s=0;s<c.samples;++s)prev[s*d]=1;
        size_t maximum=first?c.h:c.k;auto root=first?c.fec:c.lec;
        for(size_t j=1;j<=maximum;++j){auto K=c.modular_recurrence(j,first);axis old(irr::chain_copies(root,j-1,first)),next(irr::chain_copies(root,j,first));size_t pos=first?j-1:w-j;
            auto T=[&](const ulong* src,ulong* dst,size_t letter){
                for(size_t s=0;s<c.samples;++s)for(size_t q=0;q<d;++q){
                    dst[s*d+q]=nmod_mul(src[s*d+q],weights[pos][s].regular[letter],c.F.mod);
                    if(q)dst[s*d+q]=nmod_add(dst[s*d+q],nmod_mul(src[s*d+q-1],weights[pos][s].logarithmic[letter],c.F.mod),c.F.mod);
                }};
            prev=step(c,irr::matrices{},old,next,prev,stride,T,(mrk?"mrk_":"collinear_")+std::string(first?"head_":"tail_")+std::to_string(j),K.get());
        }return prev;
    };
    auto head=evolve(true),tail=evolve(false);axis fa(c.sew.first),la(c.sew.last);auto start=c.now();
    equations rows(c.samples*d);auto rhsdata=mrk?&symbols(c):nullptr;
    c.opt->pool.detach_loop(0,c.samples,[&](size_t s){auto rhs=mrk?nmhv_physical_mrk_rhs(w/2,nmhv_physical_sample_mrk_slots(c.seed,c.sample_offset+s,w-1,c.F),*rhsdata,degree,c.F):std::vector<ulong>{nmhv_physical_collinear_rhs(w/2,nmhv_physical_sample_collinear_slots(c.seed,c.sample_offset+s,w,c.F),c.F)};
        for(size_t q=0;q<d;++q){auto& r=rows[s*d+q];r.condition=mrk?"MRK":"collinear";r.sample=s;r.rhs=c.mhv?nmod_mul(rhs[q],n_invmod(2,c.F.mod.n),c.F.mod):rhs[q];r.coefficients.resize(c.rows);}});c.opt->pool.wait();
    c.record(mrk?"mrk_rhs":"collinear_rhs",rows.size(),1,0,rows.size()*8,start);start=c.now();
    // Form each pair's convolution once, not once for every nonzero solution
    // coefficient. Only one pair image is kept, not head x tail x all slots.
    auto pair_image=[&](size_t a,size_t b,ulong* image,ulong scale){
        for(size_t s=0;s<c.samples;++s)for(size_t q=0;q<d;++q){ulong x=0;for(size_t j=0;j<=q;++j)x=nmod_add(x,nmod_mul(head[a*stride+s*d+j],tail[b*stride+s*d+q-j],c.F.mod),c.F.mod);
            image[s*d+q]=nmod_add(image[s*d+q],nmod_mul(scale,x,c.F.mod),c.F.mod);}
    };
    size_t nnz=0;
    apply_direct_invariant_columns(c,rows,stride,[&](size_t b,size_t j,ulong* image){
        for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu)pair_image(fa.off[b]+mu*fa.m[b]+j/la.m[b],la.off[b]+mu*la.m[b]+j%la.m[b],image,mu?metric:1);});
    for(const auto& row:rows)for(auto x:row.coefficients)nnz+=bool(x);
    c.record(mrk?"mrk_direct":"collinear_direct",rows.size(),c.rows,nnz,(head.size()+tail.size()+rows.size()*c.rows)*8,start);return rows;
}

// Keep the last k boundary slots exact, without constructing a first-side
// intertwiner. The RHS uses the SAME suffix coordinates, including outside
// the LEC image. Validation deliberately remains the scalar slots() path.
inline equations direct_collinear_tail_exact(context& c){
    require(c.k<=3,"Exact boundary tail is limited to k<=3; use scalar projection");
    size_t w=c.h+c.k,n=c.samples,q=size_t(1)<<(2*c.k);axis fa(c.sew.first),la(c.sew.last);
    std::vector<std::vector<nmhv_limit_slot_weights_t<ulong>>> weights(c.h,std::vector<nmhv_limit_slot_weights_t<ulong>>(n));
    for(size_t s=0;s<n;++s){uint64_t state=c.seed^(0xd1b54a32d192ed03ULL*(c.sample_offset+s+1));
        for(size_t j=0;j<c.h;++j)weights[j][s]=nmhv_limit_source_weights<ulong>(nmhv_limit_kind_t::collinear,j+1,state,c.F);}
    std::vector<ulong> head(n,1);
    for(size_t j=1;j<=c.h;++j){auto K=c.modular_recurrence(j,true);axis old(irr::chain_copies(c.fec,j-1,true)),next(irr::chain_copies(c.fec,j,true));
        head=step(c,irr::matrices{},old,next,head,n,[&](const ulong* src,ulong* dst,size_t letter){
            for(size_t s=0;s<n;++s)dst[s]=nmod_mul(src[s],weights[j-1][s].regular[letter],c.F.mod);
        },"collinear_head_"+std::to_string(j),K.get());}
    if(c.exact_tail.empty()){
        c.samples=q;auto tail=terminal(c,false,q);c.samples=n;
        for(size_t j=1;j<=c.k;++j){auto K=c.modular_recurrence(j,false);axis old(irr::chain_copies(c.lec,j-1,false)),next(irr::chain_copies(c.lec,j,false));
            tail=step(c,irr::matrices{},old,next,tail,q,[&](const ulong* src,ulong* dst,size_t letter){
                for(size_t word=0;word<q;++word){auto v=nmhv_limit_scalar<ulong>(nmhv_limit_collinear_map()[letter][(word>>(2*(j-1)))&3],1,c.F);
                    dst[word]=nmod_mul(src[word],v,c.F.mod);}
            },"collinear_exact_tail_"+std::to_string(j),K.get());}
        c.exact_tail=std::move(tail);
        for(size_t word=0;word<q;++word){bool rhs=true;for(size_t j=0;j<c.k;++j)rhs=rhs&&((word>>(2*j))&3)<2;
            bool used=rhs;for(size_t j=0;j<la.dim&&!used;++j)used=c.exact_tail[j*q+word]!=0;
            if(used)c.tail_words.push_back(word);}
        std::cout<<"EXACT_TAIL words="<<c.tail_words.size()<<"/"<<q<<" entries="<<c.exact_tail.size()<<std::endl;
    }
    size_t t=c.tail_words.size(),stride=n*t;equations rows(stride);auto start=c.now();
    irr_parallel::rows(n,&c.opt->pool,[&](size_t s){
        auto slots=nmhv_physical_sample_collinear_slots(c.seed,c.sample_offset+s,w,c.F);
        for(size_t beta=0;beta<t;++beta){size_t word=c.tail_words[beta];bool rhs=true;
            for(size_t j=0;j<c.k;++j){size_t letter=(word>>(2*(c.k-j-1)))&3;slots[c.h+j]={};slots[c.h+j][letter]=1;rhs=rhs&&letter<2;}
            auto& row=rows[s*t+beta];row.condition="collinear";row.sample=c.sample_offset+s;row.coefficients.resize(c.rows);
            row.rhs=rhs?nmhv_physical_collinear_rhs(w/2,slots,c.F):0;
            if(c.mhv)row.rhs=nmod_mul(row.rhs,n_invmod(2,c.F.mod.n),c.F.mod);
        }
    });
    c.record("collinear_rhs",rows.size(),1,0,rows.size()*8,start);start=c.now();
    apply_direct_invariant_columns(c,rows,stride,[&](size_t b,size_t j,ulong* image){
        for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){size_t i=fa.off[b]+mu*fa.m[b]+j/la.m[b],l=la.off[b]+mu*la.m[b]+j%la.m[b];
            for(size_t s=0;s<n;++s)for(size_t beta=0;beta<t;++beta){auto v=c.exact_tail[l*q+c.tail_words[beta]];
                auto x=nmod_mul(head[i*n+s],v,c.F.mod);if(mu)x=nmod_mul(x,3,c.F.mod);
                image[s*t+beta]=nmod_add(image[s*t+beta],x,c.F.mod);}
        }
    });
    size_t nnz=0;for(const auto& row:rows)for(auto x:row.coefficients)nnz+=bool(x);
    c.record("collinear_direct_tail_exact",rows.size(),c.rows,nnz,(head.size()+c.exact_tail.size()+rows.size()*c.rows)*8,start);
    return rows;
}


} // namespace irr_physics
