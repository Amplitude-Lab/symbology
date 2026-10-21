#pragma once
#include "reduced_d3/irr_conversion.hpp"
#include "mhv_ope_germ.hpp"
#include "nmhv_ope_germ.hpp"
#include "reduced_d3/platform_resources.hpp"
#include "irr_runtime.hpp"
#include "irr_stage.hpp"
#include "irr_recurrence_store.hpp"
#include <chrono>
#include <thread>
#include <memory>
#include <set>
#include "shared/hexagon_mrk_rhs_reuse.hpp"
#include "shared/hexagon_native_cache_release.hpp"

namespace irr_physics {
namespace fs=std::filesystem;
using irr::require;
using mat=sparse_mat<ulong,int32_t>;
using equations=std::vector<nmhv_physical_equation_t>;
struct axis {
    irr::copies m,off{};size_t dim=0;
    explicit axis(irr::copies v):m(v){for(size_t b=0;b<6;++b){off[b]=dim;dim+=m[b]*irr::specs()[b].irrep_dimension;}}
};
struct context {
    std::shared_ptr<irr::recurrence_store> shared_recurrence;
    fs::path fec,lec,out,data;size_t h,k,samples;const char variant='C';bool mhv;
    field_t F;rref_option_t opt;irr::channels_t cg;irr::sewing_object_t sew;
    std::array<mat,6> S;size_t rows=0,sample_offset=0;
    bool shoup=true;size_t scratch_limit=64ULL<<20;
    // Optional, bounded cache of first-side OPE jets. It contains no sewing
    // coefficients and is discarded before a new prime's large solve.
    size_t ope_head_cache_limit=0,ope_head_cache_bytes=0;ulong ope_head_cache_prime=0;
    std::map<std::string,std::shared_ptr<const std::vector<ulong>>> ope_head_cache;
    void clear_ope_head_cache(){ope_head_cache.clear();ope_head_cache_bytes=0;ope_head_cache_prime=0;}
    struct step_plan {
        std::vector<std::tuple<size_t,size_t,size_t>> decode;
        std::vector<unsigned char> used;
    };
    std::map<std::pair<std::string,size_t>,step_plan> step_plans;
    ulong plan_prime=0;
    std::vector<ulong> exact_tail;
    std::vector<size_t> tail_words;
    std::optional<nmhv_physical_symbol_table_t> mrk_symbols;
    std::map<std::tuple<ulong,uint64_t,size_t,size_t,size_t>,std::vector<ulong>> mrk_rhs_cache;
    std::set<fs::path> recorded_inputs;
    std::map<std::pair<bool,size_t>,std::shared_ptr<std::array<mat,6>>> recurrence_cache;
    ulong cache_prime=0;size_t cache_bytes=0,cache_limit=1024ULL<<20;
    size_t cache_hits=0,cache_misses=0,cache_bypasses=0;
    std::set<std::pair<bool,size_t>> cache_reported;
    std::ofstream metrics,io;size_t peak_nnz=0;uint64_t seed=0x193a4527;
    std::chrono::steady_clock::time_point epoch=std::chrono::steady_clock::now();
    context(fs::path f,fs::path l,fs::path o,fs::path d,size_t H,size_t K,size_t n,bool m,ulong p)
        :fec(f),lec(l),out(o),data(d),h(H),k(K),samples(n),mhv(m),F(FIELD_Fp,p),cg(&opt->pool) {
        irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
        if(auto value=std::getenv("PHYSICS_CACHE_MIB"))cache_limit=std::stoull(value)*(1ULL<<20);
        if(auto value=std::getenv("OPE_HEAD_CACHE_MIB")){
            auto mib=std::stoull(value);require(mib<=(SIZE_MAX>>20),"OPE head cache budget overflow");ope_head_cache_limit=mib*(1ULL<<20);
        }
        std::cout<<"OPE_HEAD_CACHE budget_bytes="<<ope_head_cache_limit<<std::endl;
        std::cout<<"CACHE budget_bytes="<<cache_limit<<std::endl;
        sew.first=irr::chain_copies(fec,h,true);sew.last=irr::chain_copies(lec,k,false);
        fs::create_directories(out);metrics.open(out/"phases.tsv");io.open(out/"inputs.sha256");
        metrics<<"phase\trows\tcolumns\tnnz\tobject_bytes\tseconds\tpeak_rss\tpeak_private\tstart\tend\n";
    }
    void input(const fs::path& p){auto path=fs::absolute(p).lexically_normal();if(recorded_inputs.insert(path).second){io<<hex_sha256_file(path)<<"  "<<path.generic_string()<<'\n';io.flush();}}
    auto now(){return std::chrono::steady_clock::now();}
    void axpy(ulong* dst,const ulong* src,size_t n,ulong v)const{
        if(!v)return;
        if(v==1){for(size_t i=0;i<n;++i)dst[i]=nmod_add(dst[i],src[i],F.mod);return;}
        if(v==F.mod.n-1){for(size_t i=0;i<n;++i)dst[i]=nmod_sub(dst[i],src[i],F.mod);return;}
        if(shoup){auto pre=n_mulmod_precomp_shoup(v,F.mod.n);
            for(size_t i=0;i<n;++i)if(src[i])dst[i]=nmod_add(dst[i],n_mulmod_shoup(v,src[i],pre,F.mod.n),F.mod);
        }else for(size_t i=0;i<n;++i)if(src[i])dst[i]=nmod_add(dst[i],nmod_mul(v,src[i],F.mod),F.mod);
    }
    void record(const std::string& phase,size_t r,size_t c,size_t nnz,size_t bytes,std::chrono::steady_clock::time_point start){
        auto mem=hex_platform::current_process_memory();peak_nnz=std::max(peak_nnz,nnz);
        auto end=now();metrics<<phase<<'\t'<<r<<'\t'<<c<<'\t'<<nnz<<'\t'<<bytes<<'\t'<<std::chrono::duration<double>(end-start).count()<<'\t'<<mem.peak_working_set_bytes<<'\t'<<mem.peak_private_commit_bytes<<'\t'<<std::chrono::duration<double>(start-epoch).count()<<'\t'<<std::chrono::duration<double>(end-epoch).count()<<'\n';metrics.flush();
        std::cout<<phase<<" rows="<<r<<" cols="<<c<<" nnz="<<nnz<<" seconds="<<std::chrono::duration<double>(end-start).count()
            <<" peak_rss="<<mem.peak_working_set_bytes<<" cache_bytes="<<cache_bytes<<" private="<<mem.private_commit_bytes<<std::endl;
    }
    std::shared_ptr<const irr::matrices> recurrence(size_t w,bool first){
        auto root=first?fec:lec;
        for(size_t b=0;b<6;++b)input(irr::tensor_path(root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf")));
        if(shared_recurrence)return shared_recurrence->rational(root,w,first,cg,&opt->pool);
        return std::make_shared<irr::matrices>(irr::read_recurrence(root,w,first,cg,&opt->pool));
    }
    // The legacy context caches residues only; the pipeline's optional shared
    // store can also retain QQ inputs. Neither store reuses Fp across primes.
    void cache_summary(){
        if(shared_recurrence){shared_recurrence->report();return;}
        std::cout<<"CACHE batch prime="<<F.mod.n<<" hits="<<cache_hits<<" misses="<<cache_misses
            <<" not_stored="<<cache_bypasses<<" entries="<<recurrence_cache.size()
            <<" used_bytes="<<cache_bytes<<" limit_bytes="<<cache_limit
            <<" remaining_bytes="<<cache_limit-cache_bytes<<std::endl;
        cache_hits=cache_misses=cache_bypasses=0;
    }
    std::shared_ptr<const std::array<mat,6>> modular_recurrence(size_t w,bool first){
        if(cache_prime!=F.mod.n){
            recurrence_cache.clear();cache_reported.clear();step_plans.clear();exact_tail.clear();tail_words.clear();cache_bytes=0;cache_prime=F.mod.n;
            std::cout<<"CACHE reset prime="<<cache_prime<<" limit_bytes="<<cache_limit<<std::endl;
        }
        if(shared_recurrence)return shared_recurrence->modular(first?fec:lec,w,first,cg,F,&opt->pool);
        auto key=std::make_pair(first,w);auto it=recurrence_cache.find(key);
        if(it!=recurrence_cache.end()){++cache_hits;return it->second;}
        ++cache_misses;
        auto K=recurrence(w,first);auto value=std::make_shared<std::array<mat,6>>();size_t bytes=0;
        for(size_t b=0;b<6;++b){(*value)[b]=nmhv_sew_modular_reduce((*K)[b],F,&opt->pool);
            bytes+=(*value)[b].alloc()*(sizeof(ulong)+sizeof(int32_t))+(*value)[b].rows.capacity()*sizeof((*value)[b].rows[0]);}
        bool stored=bytes<=cache_limit-cache_bytes;
        if(stored){recurrence_cache.emplace(key,value);cache_bytes+=bytes;}
        else ++cache_bypasses;
        if(cache_reported.insert(key).second)
            std::cout<<"CACHE "<<(stored?"stored":"NOT_STORED")<<" prime="<<cache_prime
                <<" tensor="<<(first?"irrFEC_":"irrLEC_")<<w<<" required_bytes="<<bytes
                <<" used_bytes="<<cache_bytes<<" limit_bytes="<<cache_limit
                <<" remaining_bytes="<<cache_limit-cache_bytes
                <<" reason="<<(stored?"fits":cache_limit?"insufficient_capacity":"disabled")<<std::endl;
        if(bytes>=size_t(64)*1024*1024){
            K.reset(); // Release temporary QQ ownership before freeing only unused integer caches.
            hexagon_native_memory::release_unused_integer_caches(&opt->pool);
        }
        return value;
    }
    void load(const fs::path& root){
        auto t=now();sew.h=h;sew.k=k;sew.first=irr::chain_copies(fec,h,true);sew.last=irr::chain_copies(lec,k,false);
        input(root/"columns.tsv");std::ifstream cols(root/"columns.tsv");std::string line;std::getline(cols,line);size_t base=0;
        for(size_t b=0;b<6;++b){
            std::string block,file;size_t nf,nl,off;cols>>block>>nf>>nl>>off>>file;
            require(bool(cols)&&block==irr::specs()[b].name&&nf==sew.first[b]&&nl==sew.last[b]&&off==base,"Invalid sewing metadata");base+=nf*nl;
            input(irr::tensor_path(root/file));auto x=irr::read_tensor(root/file,&opt->pool);
            if(!b)sew.directions=x.dim(0);
            require(x.rank()==3&&x.dim(0)==sew.directions&&x.dim(1)==nf&&x.dim(2)==nl,"Invalid shared solution axis");
            sew.S[b]=joint_eet_flatten_tensor(x,&opt->pool);
        }
        size_t nnz=0;for(auto& s:sew.S)nnz+=s.nnz();record("load_bulk",sew.directions,base,nnz,0,t);
    }
};

inline mat physical_sewing(context& c,int metric=3){
    axis a(c.sew.first),b(c.sew.last);mat S(c.rows,a.dim*b.dim);
    for(size_t q=0;q<6;++q)for(size_t r=0;r<S.nrow;++r)for(auto [j,v]:c.S[q][r])
        for(size_t mu=0;mu<irr::specs()[q].irrep_dimension;++mu){
            size_t x=a.off[q]+mu*a.m[q]+j/b.m[q],y=b.off[q]+mu*b.m[q]+j%b.m[q];
            S[r].push_back(int32_t(x*b.dim+y),nmod_mul(v,mu?metric:1,c.F.mod));
        }
    for(auto& r:S.rows)r.canonicalize();return S;
}

// T maps one prefix state to its state after appending a letter. For a LEC
// slot evaluation the slot number is reversed; the scalar slot factors commute.
template<class Transfer>
std::vector<ulong> step(context& c,const irr::matrices& K,const axis& old,const axis& next,
    const std::vector<ulong>& previous,size_t stride,Transfer T,const std::string& tag,const std::array<mat,6>* modular=nullptr){
    if(c.plan_prime!=c.F.mod.n){c.step_plans.clear();c.plan_prime=c.F.mod.n;}
    auto start=c.now();std::vector<ulong> current(next.dim*stride),image(c.variant=='C'?0:old.dim*stride);
    size_t nnz=0;
    if(c.variant!='C'){
        for(size_t b=0;b<6;++b){
            auto R=nmhv_sew_modular_reduce(irr::expand(c.cg,b,old.m,K[b]),c.F,&c.opt->pool);nnz+=R.nnz();
            for(size_t letter=0;letter<9;++letter){
                c.opt->pool.detach_loop(0,old.dim,[&](size_t r){T(previous.data()+r*stride,image.data()+r*stride,letter);});c.opt->pool.wait();
                c.opt->pool.detach_loop(0,R.nrow,[&](size_t r){
                    auto dst=current.data()+(next.off[b]+r)*stride;
                    for(auto [j,v]:R[r])if(j%9==letter){auto src=image.data()+(j/9)*stride;
                        for(size_t z=0;z<stride;++z)if(src[z])dst[z]=nmod_add(dst[z],nmod_mul(v,src[z],c.F.mod),c.F.mod);}
                });c.opt->pool.wait();
            }
        }
    }else{
        // Stream a bounded batch of channels. Parallel channel evaluation and
        // destination-row gathering need no per-thread copy of the full germ.
        for(size_t b=0;b<6;++b){
            mat temporary;if(!modular)temporary=nmhv_sew_modular_reduce(K[b],c.F,&c.opt->pool);
            const auto& M=modular?(*modular)[b]:temporary;nnz+=M.nnz();
            auto key=std::make_pair(tag,b);auto it=c.step_plans.find(key);
            if(it==c.step_plans.end()){
                context::step_plan plan;
                for(size_t a=0;a<6;++a)for(size_t copy=0;copy<old.m[a];++copy)for(size_t g:c.cg.channels[b][a])plan.decode.push_back({a,copy,g});
                require(plan.decode.size()==M.ncol,"Recurrence channel layout mismatch");
                plan.used.resize(M.ncol);for(const auto& row:M.rows)for(auto [j,v]:row)plan.used[j]=1;
                it=c.step_plans.emplace(key,std::move(plan)).first;
            }
            const auto& decode=it->second.decode;const auto& used=it->second.used;
            // A full small-germ block avoids rescanning CSR rows for 64-column
            // slices. Large OPE germs still use a bounded scratch allocation.
            size_t batch=std::max(size_t(1),c.scratch_limit/(sizeof(ulong)*stride));
            for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){
                for(size_t first=0;first<M.ncol;first+=batch){
                    size_t count=std::min(batch,M.ncol-first);std::vector<ulong> candidates(count*stride);
                    irr_parallel::rows(count,&c.opt->pool,[&](size_t col){
                        if(!used[first+col])return;
                        auto [a,copy,g]=decode[first+col];auto candidate=candidates.data()+col*stride;std::vector<ulong> term_image(stride);
                        for(auto& term:c.cg.cg[b].vectors[g].component[mu]){
                            size_t r=old.off[a]+(term.coordinate/9-c.cg.axis[a].offset)*old.m[a]+copy;
                            T(previous.data()+r*stride,term_image.data(),term.coordinate%9);ulong v=term.coefficient%c.F.mod;
                            c.axpy(candidate,term_image.data(),stride,v);
                        }
                    });
                    irr_parallel::rows(M.nrow,&c.opt->pool,[&](size_t r){
                        auto dst=current.data()+(next.off[b]+mu*next.m[b]+r)*stride;size_t lo=0,hi=M[r].nnz();
                        while(lo<hi){size_t mid=(lo+hi)/2;if(size_t(M[r](mid))<first)lo=mid+1;else hi=mid;}
                        for(size_t p=lo;p<M[r].nnz()&&size_t(M[r](p))<first+count;++p){
                            auto candidate=candidates.data()+(M[r](p)-first)*stride;auto v=M[r][p];
                            c.axpy(dst,candidate,stride,v);
                        }
                    });
                }
            }
        }
    }
    c.record(tag,next.dim,old.dim*9,nnz,(previous.size()+current.size()+image.size())*sizeof(ulong),start);return current;
}

inline std::vector<ulong> terminal(context& c,bool mrk,size_t stride){
    c.input(c.lec/"LEC_0_components.wxf");auto t=irr::read_tensor(c.lec/"LEC_0_components.wxf",&c.opt->pool);
    auto M=joint_eet_flatten_tensor(t,&c.opt->pool);std::vector<ulong> values(M.nrow*stride);
    std::vector<int> weights=c.mhv?std::vector<int>{1}:(mrk?std::vector<int>{1,0,1,2,-1}:std::vector<int>{1,0,1,0,0});
    for(size_t r=0;r<M.nrow;++r){rat_t x=0;for(auto [j,v]:M[r])x+=v*rat_t(weights.at(j));for(size_t s=0;s<c.samples;++s)values[r*stride+s*(mrk?c.h+c.k:1)]=x%c.F.mod;}return values;
}

inline const nmhv_physical_symbol_table_t& symbols(context& c){
    if(c.mrk_symbols)return *c.mrk_symbols;
    nmhv_physical_symbol_table_t out;
    for(size_t l=1;l<=(c.h+c.k)/2;++l)for(size_t k=0;k<l;++k){
        auto p=c.data/(c.mhv?"mrk_mhv":"mrk")/("bfkl_f_L"+std::to_string(l)+"_k"+std::to_string(k)+".tsv");c.input(p);
        std::ifstream in(p);std::string line;std::getline(in,line);require(line=="word\tcoefficient","MRK header");
        while(std::getline(in,line))if(!line.empty()){auto f=nmhv_limit_split_tabs(line);nmhv_physical_symbol_term_t t;for(char x:f.at(0))t.word.push_back(x-'0');t.coefficient=rat_t(f.at(1));out[{l,k}].push_back(std::move(t));}
    }c.mrk_symbols=std::move(out);return *c.mrk_symbols;
}

#include "shared/hexagon_mrk_rhs_batch.inc"

// Solution rows own the destination. Reuse every sparse coefficient over the
// whole sample block, with no sewing transpose and no concurrent scatter-add.
template<class Image>
void apply_rows(context& c,equations& rows,size_t stride,Image image_of){
    for(size_t b=0;b<6;++b){const auto& S=c.S[b];
        size_t batch=std::max(size_t(1),c.scratch_limit/(sizeof(ulong)*stride));
        for(size_t first=0;first<S.ncol;first+=batch){size_t count=std::min(batch,S.ncol-first);
            std::vector<ulong> images(count*stride);
            irr_parallel::rows(count,&c.opt->pool,[&](size_t j){image_of(b,first+j,images.data()+j*stride);});
            irr_parallel::rows(S.nrow,&c.opt->pool,[&](size_t r){
                std::vector<ulong> sums(stride);size_t lo=0,hi=S[r].nnz();
                if(first)while(lo<hi){size_t mid=(lo+hi)/2;if(size_t(S[r](mid))<first)lo=mid+1;else hi=mid;}
                for(size_t p=lo;p<S[r].nnz()&&size_t(S[r](p))<first+count;++p)
                    c.axpy(sums.data(),images.data()+(S[r](p)-first)*stride,stride,S[r][p]);
                for(size_t s=0;s<stride;++s)rows[s].coefficients[r]=nmod_add(rows[s].coefficients[r],sums[s],c.F.mod);
            });
        }
    }
}

inline equations slots(context& c,bool mrk,int metric=3){
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
    equations rows(c.samples*d);auto rhs_values=mrk?mrk_rhs_reuse_batch(c,degree):std::vector<std::vector<ulong>>{};
    c.opt->pool.detach_loop(0,c.samples,[&](size_t s){auto rhs=mrk?rhs_values.at(s):std::vector<ulong>{nmhv_physical_collinear_rhs(w/2,nmhv_physical_sample_collinear_slots(c.seed,c.sample_offset+s,w,c.F),c.F)};
        for(size_t q=0;q<d;++q){auto& r=rows[s*d+q];r.condition=mrk?"MRK":"collinear";r.sample=s;r.rhs=c.mhv?nmod_mul(rhs[q],n_invmod(2,c.F.mod.n),c.F.mod):rhs[q];r.coefficients.resize(c.rows);}});c.opt->pool.wait();
    c.record(mrk?"mrk_rhs":"collinear_rhs",rows.size(),1,0,rows.size()*8,start);start=c.now();
    // Form each pair's convolution once, not once for every nonzero solution
    // coefficient. Only one pair image is kept, not head x tail x all slots.
    auto pair_image=[&](size_t a,size_t b,ulong* image,ulong scale){
        for(size_t s=0;s<c.samples;++s)for(size_t q=0;q<d;++q){ulong x=0;for(size_t j=0;j<=q;++j)x=nmod_add(x,nmod_mul(head[a*stride+s*d+j],tail[b*stride+s*d+q-j],c.F.mod),c.F.mod);
            image[s*d+q]=nmod_add(image[s*d+q],nmod_mul(scale,x,c.F.mod),c.F.mod);}
    };
    size_t nnz=0;for(auto& S:c.S)nnz+=S.nnz();
    apply_rows(c,rows,stride,[&](size_t b,size_t j,ulong* image){
        for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu)pair_image(fa.off[b]+mu*fa.m[b]+j/la.m[b],la.off[b]+mu*la.m[b]+j%la.m[b],image,mu?metric:1);});
    c.record(mrk?"mrk_apply":"collinear_apply",rows.size(),c.rows,nnz,(head.size()+tail.size()+rows.size()*c.rows)*8,start);return rows;
}

// Keep the last k boundary slots exact, without constructing a first-side
// intertwiner. The RHS uses the SAME suffix coordinates, including outside
// the LEC image. Validation deliberately remains the scalar slots() path.
inline equations collinear_tail_exact(context& c){
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
    apply_rows(c,rows,stride,[&](size_t b,size_t j,ulong* image){
        for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){size_t i=fa.off[b]+mu*fa.m[b]+j/la.m[b],l=la.off[b]+mu*la.m[b]+j%la.m[b];
            for(size_t s=0;s<n;++s)for(size_t beta=0;beta<t;++beta){auto v=c.exact_tail[l*q+c.tail_words[beta]];
                auto x=nmod_mul(head[i*n+s],v,c.F.mod);if(mu)x=nmod_mul(x,3,c.F.mod);
                image[s*t+beta]=nmod_add(image[s*t+beta],x,c.F.mod);}
        }
    });
    size_t nnz=0;for(auto& S:c.S)nnz+=S.nnz();c.record("collinear_apply",rows.size(),c.rows,nnz,(head.size()+c.exact_tail.size()+rows.size()*c.rows)*8,start);
    return rows;
}

// OPE propagation goes from the head through the tail. It is not a product
// of head and tail germs. K is transposed first, then the fixed CG is applied.
template<class Germ,class Walk>
void accumulate_ope_seed(std::vector<ulong>& values,const std::vector<ulong>& head,size_t stride,
    const Germ& germ,thread_pool* pool,const Walk& walk){
    size_t workers=pool?std::min(size_t(pool->get_thread_count()),std::max(size_t(1),stride/256)):1;
    irr_parallel::rows(workers,pool,[&](size_t worker){
        size_t begin=worker*stride/workers,end=(worker+1)*stride/workers;
        walk([&](size_t i,size_t j,ulong v){
            auto* dst=values.data()+j*stride;const auto* src=head.data()+i*stride;
            for(size_t z=begin;z<end;++z)dst[z]=germ.add(dst[z],germ.mul(v,src[z]));
        });
    });
}

template<class Germ>
inline std::vector<ulong> ope_values(context& c,const Germ& germ,const std::string& head_identity={}){
    size_t stride=germ.stride;
    if(c.ope_head_cache_prime!=c.F.mod.n){c.clear_ope_head_cache();c.ope_head_cache_prime=c.F.mod.n;}
    auto key=head_identity+":"+std::to_string(c.h)+":"+std::to_string(stride);
    std::shared_ptr<const std::vector<ulong>> head_storage;
    if(!head_identity.empty())if(auto it=c.ope_head_cache.find(key);it!=c.ope_head_cache.end())head_storage=it->second;
    if(head_storage)std::cout<<"OPE_HEAD_CACHE hit prime="<<c.F.mod.n<<" bytes="<<head_storage->size()*sizeof(ulong)<<std::endl;
    else{
        std::vector<ulong> calculated(stride);calculated[germ.index(0,0,0,0)]=1;
        for(size_t w=1;w<=c.h;++w){auto K=c.modular_recurrence(w,true);axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
            irr::matrices unused;calculated=step(c,unused,old,next,calculated,stride,[&](const ulong* a,ulong* b,size_t l){germ.integrate(a,b,l,w);},"ope_head_"+std::to_string(w),K.get());}
        head_storage=std::make_shared<const std::vector<ulong>>(std::move(calculated));
        size_t bytes=head_storage->size()*sizeof(ulong);
        if(!head_identity.empty()&&bytes<=c.ope_head_cache_limit-c.ope_head_cache_bytes){c.ope_head_cache.emplace(key,head_storage);c.ope_head_cache_bytes+=bytes;
            std::cout<<"OPE_HEAD_CACHE stored prime="<<c.F.mod.n<<" bytes="<<bytes<<" used="<<c.ope_head_cache_bytes<<std::endl;}
    }
    const auto& head=*head_storage;
    std::vector<std::shared_ptr<const std::array<mat,6>>> chain;for(size_t k=1;k<=c.k;++k)chain.push_back(c.modular_recurrence(k,false));
    mat expanded;if(c.variant=='A')expanded=physical_sewing(c);
    axis fa(c.sew.first),la(c.sew.last);size_t components=c.mhv?1:5;
    mat terminal_map;
    if(!c.mhv){auto p=c.lec/"LEC_0_components.wxf";c.input(p);auto x=irr::read_tensor(p,&c.opt->pool);terminal_map=nmhv_sew_modular_reduce(joint_eet_flatten_tensor(x,&c.opt->pool),c.F,&c.opt->pool);}
    std::vector<ulong> result(c.rows*components*stride);auto start=c.now();double seed_seconds=0;
    for(size_t solution=0;solution<c.rows;++solution){
        std::vector<ulong> values(la.dim*stride);
        auto seed_start=c.now();
        accumulate_ope_seed(values,head,stride,germ,&c.opt->pool,[&](auto seed){
            if(c.variant=='A'){for(auto [j,v]:expanded[solution])seed(j/la.dim,j%la.dim,v);}
            else for(size_t b=0;b<6;++b)for(auto [j,v]:c.S[b][solution])for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu)
                seed(fa.off[b]+mu*fa.m[b]+j/la.m[b],la.off[b]+mu*la.m[b]+j%la.m[b],germ.mul(v,mu?3:1));
        });
        seed_seconds+=std::chrono::duration<double>(c.now()-seed_start).count();
        for(size_t k=c.k;k>0;--k){
            axis current(irr::chain_copies(c.lec,k,false)),old(irr::chain_copies(c.lec,k-1,false));
            std::vector<ulong> next(old.dim*stride),image(stride);const auto& K=*chain[k-1];size_t w=c.h+c.k-k+1;
            if(c.variant!='C'){
                auto rational=c.recurrence(k,false);
                for(size_t b=0;b<6;++b){auto R=nmhv_sew_modular_reduce(irr::expand(c.cg,b,old.m,(*rational)[b]),c.F,&c.opt->pool);
                    for(size_t r=0;r<R.nrow;++r)for(size_t letter=0;letter<9;++letter){
                        bool used=false;for(auto [j,v]:R[r])if(j%9==letter){used=true;break;}if(!used)continue;
                        germ.integrate(values.data()+(current.off[b]+r)*stride,image.data(),letter,w);
                        for(auto [j,v]:R[r])if(j%9==letter)for(size_t z=0;z<stride;++z)if(image[z])next[(j/9)*stride+z]=germ.add(next[(j/9)*stride+z],germ.mul(v,image[z]));
                    }
                }
            }else for(size_t b=0;b<6;++b){auto M=K[b];M.transpose_replace();
                std::vector<std::pair<size_t,size_t>> copies;irr::copies base{};size_t columns=0;
                for(size_t a=0;a<6;++a){base[a]=columns;columns+=old.m[a]*c.cg.channels[b][a].size();for(size_t copy=0;copy<old.m[a];++copy)copies.push_back({a,copy});}
                for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){
                    // Distinct old copies own disjoint terminal rows. Cap the
                    // two germ scratch arrays per worker, not the solution count.
                    size_t workers=std::max(size_t(1),size_t(16ULL<<20)/(2*sizeof(ulong)*stride));
                    irr_parallel::rows(copies.size(),&c.opt->pool,[&](size_t owner){
                    auto [a,copy]=copies[owner];size_t col=base[a]+copy*c.cg.channels[b][a].size();
                    std::vector<ulong> candidate(stride),term_image(stride);
                    for(size_t g:c.cg.channels[b][a]){
                        std::fill(candidate.begin(),candidate.end(),0);for(auto [r,v]:M[col])for(size_t z=0;z<stride;++z){auto x=values[(current.off[b]+mu*current.m[b]+r)*stride+z];if(x)candidate[z]=germ.add(candidate[z],germ.mul(v,x));}
                        for(auto& term:c.cg.cg[b].vectors[g].component[mu]){
                            germ.integrate(candidate.data(),term_image.data(),term.coordinate%9,w);
                            size_t r=old.off[a]+(term.coordinate/9-c.cg.axis[a].offset)*old.m[a]+copy;ulong v=term.coefficient%c.F.mod;
                            for(size_t z=0;z<stride;++z)if(term_image[z])next[r*stride+z]=germ.add(next[r*stride+z],germ.mul(v,term_image[z]));
                        }++col;
                    }
                    },workers);
                }
            }
            values=std::move(next);
        }
        if(c.mhv){require(values.size()==stride,"MHV terminal dimension");std::copy(values.begin(),values.end(),result.begin()+solution*stride);}
        else {require(values.size()==terminal_map.nrow*stride&&terminal_map.ncol==components,"NMHV terminal dimension");
            for(size_t r=0;r<terminal_map.nrow;++r)for(auto [j,v]:terminal_map[r])c.axpy(result.data()+(solution*components+j)*stride,values.data()+r*stride,stride,v);}
    }
    size_t nnz=expanded.nnz();if(c.variant!='A')for(auto& x:c.S)nnz+=x.nnz();
    std::cout<<"OPE_SEED rows="<<c.rows<<" stride="<<stride<<" nnz="<<nnz<<" seconds="<<seed_seconds<<std::endl;
    c.record("ope_apply",c.rows,stride,nnz,(head.size()+result.size()+la.dim*stride)*8,start);return result;
}

inline equations nmhv_ope(context& c,bool validation){
    auto root=c.data/"ope_nmhv";auto metadata=root/"projection.tsv";c.input(metadata);
    std::ifstream meta(metadata);std::string line;int mt=0,grade=0;
    while(std::getline(meta,line)){line.erase(std::remove(line.begin(),line.end(),'"'),line.end());auto f=nmhv_limit_split_tabs(line);
        if(f.size()==2&&f[0]=="max_t")mt=std::stoi(f[1]);if(f.size()==2&&f[0]=="grade")grade=std::stoi(f[1]);}
    require(mt==2&&grade>=4,"NMHV OPE requires complete T1/T2 BDS-like data; partial framed W is not supported");
    struct target {int t,a,z,lt,lz;ulong rhs;std::string role;};std::vector<target> targets;
    auto path=root/("E1111_L"+std::to_string((c.h+c.k)/2)+".tsv");c.input(path);
    std::ifstream input(path);std::getline(input,line);require(bool(input),"Missing NMHV E1111 OPE RHS");
    size_t solve_count=0,check_count=0;
    while(std::getline(input,line))if(!line.empty()){
        auto f=nmhv_ope_row(line);
        require(f.size()==9,"Invalid NMHV OPE table");target x{std::stoi(f[1]),std::stoi(f[2]),std::stoi(f[3]),std::stoi(f[4]),std::stoi(f[5]),rat_t(f[6]+"/"+f[7])%c.F.mod,f[8]};
        require(std::stoul(f[0])==(c.h+c.k)/2&&x.t>=1&&x.t<=2&&std::abs(x.a)==x.t&&x.z+x.t<=grade&&x.lt>=0&&x.lz>=0&&size_t(x.lt+x.lz)<=c.h+c.k,"Unsupported NMHV OPE coefficient");
        if(x.role=="check")++check_count;else ++solve_count;
        if((x.role=="check")==validation)targets.push_back(x);
    }
    require(solve_count&&check_count,"NMHV OPE table must contain both solve and check rows");
    std::cout<<"NMHV_OPE_TABLE solve="<<solve_count<<" check="<<check_count<<" validation="<<validation<<std::endl;
    require(!targets.empty(),"Empty NMHV OPE solve/check selection");equations out(targets.size());
    for(size_t i=0;i<out.size();++i){auto x=targets[i];out[i].condition="NMHV OPE";out[i].sample=i;out[i].rhs=x.rhs;out[i].coefficients.resize(c.rows);out[i].ope_coordinate={x.z,x.t,x.a,x.lt,x.lz};}
    std::ofstream labels(c.out/(validation?"nmhv_ope_check_rows.tsv":"nmhv_ope_solve_rows.tsv"));labels<<"twist\thelicity\tpower\tlog_t\tlog_x\trhs_mod_p\n";
    for(auto x:targets)labels<<x.t<<'\t'<<x.a<<'\t'<<x.z<<'\t'<<x.lt<<'\t'<<x.lz<<'\t'<<x.rhs<<'\n';
    for(int f=2;f<=6;++f){
        // Interpolate F^2 times the Laurent polynomial, degree <=4. Do not
        // impose the unknown T^2 F^0 sector or assume equal signed helicities.
        std::vector<ulong> polynomial{1};ulong denominator=1;
        for(int other=2;other<=6;++other)if(other!=f){std::vector<ulong> next(polynomial.size()+1);
            for(size_t j=0;j<polynomial.size();++j){next[j]=nmod_sub(next[j],nmod_mul(other,polynomial[j],c.F.mod),c.F.mod);next[j+1]=nmod_add(next[j+1],polynomial[j],c.F.mod);}polynomial=std::move(next);
            denominator=nmod_mul(denominator,rat_t(f-other)%c.F.mod,c.F.mod);}
        for(auto& x:polynomial)x=nmod_mul(x,nmod_mul(f*f,n_invmod(denominator,c.F.mod.n),c.F.mod),c.F.mod);
        struct prefactor {int component,t,z;ulong value;};std::vector<prefactor> prefactors;
        auto pref=root/("prefactors_F"+std::to_string(f)+".tsv");c.input(pref);std::ifstream pin(pref);int comp,t,z;std::string val;
        while(pin>>comp>>t>>z>>val){val.erase(std::remove(val.begin(),val.end(),'"'),val.end());prefactors.push_back({comp,t,z,rat_t(val)%c.F.mod});}
        require(pin.eof()&&!prefactors.empty(),"Missing (1111) prefactors");
        auto kernel=root/("kernels_F"+std::to_string(f)+".tsv");c.input(kernel);
        for(int cycle=0;cycle<3;++cycle){
            std::cout<<"NMHV_OPE F="<<f<<" cyclic="<<cycle<<" validation="<<validation<<" rows="<<c.rows<<std::endl;
            nmhv_ope_germ_t germ(c.h+c.k,c.F.mod.n,kernel,mt,grade,cycle==1?2:cycle==2?1:0);auto values=ope_values(c,germ);
            for(size_t i=0;i<targets.size();++i){auto x=targets[i];auto weight=polynomial[x.a+2];
                for(auto term:prefactors)if(term.component/2==cycle&&germ.valid(x.t-term.t,x.z-term.z)){
                    size_t component=term.component%2?3:0;auto coefficient=germ.mul(weight,term.value);
                    for(size_t r=0;r<c.rows;++r){auto value=values[(r*5+component)*germ.stride+germ.index(x.t-term.t,x.z-term.z,x.lt,x.lz)];out[i].coefficients[r]=germ.add(out[i].coefficients[r],germ.mul(coefficient,value));}
                }
            }
        }
    }
    return out;
}

// MHV is even under F <-> 1/F. At T<=2 interpolate in q=F+1/F:
// the q and q^2 coefficients are the individual helicity +/-1 and +/-2 RHS.
// The unknown T^2 F^0 coefficient drops out of the quadratic coefficient.
inline equations ope_small_s(context& c,bool validation=false,bool include_twist1=false){
    require(c.mhv,"ope_small_s currently supports MHV only");
    auto root=c.data/"ope_mhv_small_s";c.input(root/"projection.tsv");
    std::map<std::string,std::string> config;
    std::ifstream settings(root/"projection.tsv");std::string setting;
    while(std::getline(settings,setting)){if(!setting.empty()&&setting.back()=='\r')setting.pop_back();setting.erase(std::remove(setting.begin(),setting.end(),'"'),setting.end());auto pair=nmhv_limit_split_tabs(setting);require(pair.size()==2,"Invalid small-S projection setting");config.emplace(pair[0],pair[1]);}
    require(config.at("schema")=="mhv-small-s-v1"&&config.at("endpoint")=="SmallS","Invalid small-S OPE configuration");
    int mt=std::stoi(config.at("max_t")),grade=std::stoi(config.at("grade"));
    require(mt==2&&grade>=10,"Unsupported small-S OPE truncation");
    struct target {int t,z,a,b;ulong rhs;};std::vector<target> targets;
    auto path=root/("calE_L"+std::to_string((c.h+c.k)/2)+".tsv");c.input(path);
    std::ifstream in(path);std::string line;std::getline(in,line);size_t solve_count=0,check_count=0;
    while(std::getline(in,line))if(!line.empty()){
        if(line.back()=='\r')line.pop_back();line.erase(std::remove(line.begin(),line.end(),'"'),line.end());
        auto fields=nmhv_limit_split_tabs(line);require(fields.size()==9,"Invalid MHV small-S RHS row");
        require(fields[8]=="solve"||fields[8]=="check","Invalid MHV small-S RHS role");
        int t=std::stoi(fields[1]),z=std::stoi(fields[3]),a=std::stoi(fields[4]),b=std::stoi(fields[5]);
        require(std::stoul(fields[0])==(c.h+c.k)/2&&(t==1||t==2)&&std::stoi(fields[2])==t&&z>=t&&z+t<=grade&&a>=0&&b>=0&&size_t(a+b)<=c.h+c.k,"Unsupported MHV small-S RHS coefficient");
        bool check=fields[8]=="check"||(t==1&&!include_twist1);if(check)++check_count;else ++solve_count;
        if(check==validation)targets.push_back({t,z,a,b,rat_t(fields[6]+"/"+fields[7])%c.F.mod});
    }
    require(in.eof()&&solve_count&&check_count&&!targets.empty(),"Missing MHV small-S solve/check RHS");
    std::cout<<"MHV_SMALL_S solve="<<solve_count<<" check="<<check_count<<" validation="<<validation<<std::endl;
    equations out(targets.size());for(size_t i=0;i<out.size();++i){auto x=targets[i];out[i].condition="MHV small-S OPE";out[i].sample=i;out[i].rhs=x.rhs;out[i].coefficients.resize(c.rows);out[i].ope_coordinate={x.z,x.t,x.t,x.a,x.b};}
    std::array<ulong,3> q;for(size_t j=0;j<3;++j)q[j]=nmod_add(j+2,n_invmod(j+2,c.F.mod.n),c.F.mod);
    for(size_t j=0;j<3;++j){
        size_t l=(j+1)%3,r=(j+2)%3;auto den=nmod_mul(nmod_sub(q[j],q[l],c.F.mod),nmod_sub(q[j],q[r],c.F.mod),c.F.mod);
        std::array<ulong,3> coeff{0,nmod_neg(nmod_add(q[l],q[r],c.F.mod),c.F.mod),1};
        for(auto& x:coeff)x=nmod_mul(x,n_invmod(den,c.F.mod.n),c.F.mod);
        auto kernel=root/("kernels_F"+std::to_string(j+2)+".tsv");c.input(kernel);
        nmhv_ope_germ_t germ(c.h+c.k,c.F.mod.n,kernel,mt,grade);
        std::cout<<"MHV_SMALL_S F="<<j+2<<" validation="<<validation<<" rows="<<c.rows<<std::endl;
        auto values=ope_values(c,germ,"MHV-SmallS:"+hex_sha256_file(kernel)+":"+std::to_string(mt)+":"+std::to_string(grade));
        for(size_t i=0;i<targets.size();++i){auto x=targets[i];for(size_t row=0;row<c.rows;++row)
            out[i].coefficients[row]=germ.add(out[i].coefficients[row],germ.mul(coeff[x.t],values[row*germ.stride+germ.index(x.t,x.z,x.a,x.b)]));}
    }
    return out;
}

inline equations ope(context& c,bool validation=false){
    if(!c.mhv)return nmhv_ope(c,validation);
    std::array<std::vector<ulong>,2> values;
    for(size_t f=0;f<2;++f){auto p=c.data/"ope"/("kernels_F"+std::to_string(f+2)+".tsv");c.input(p);mhv_ope_germ_t g(c.h+c.k,c.F.mod.n,p);values[f]=ope_values(c,g);}
    mhv_ope_germ_t g(c.h+c.k,c.F.mod.n,c.data/"ope/kernels_F2.tsv");
    auto p=c.data/"ope"/("calE_ope_L"+std::to_string((c.h+c.k)/2)+".tsv");c.input(p);std::ifstream in(p);std::string line;std::getline(in,line);equations out;
    require(bool(in),"Missing BDS-like OPE RHS (framed W is not a replacement)");
    auto rational=[&](std::string x){x.erase(std::remove(x.begin(),x.end(),'"'),x.end());return rat_t(x)%c.F.mod;};
    std::ofstream labels(c.out/"ope_rows.tsv");labels<<line<<'\n';
    while(std::getline(in,line))if(!line.empty()){
        labels<<line<<'\n';auto f=nmhv_limit_split_tabs(line);int t=std::stoi(f[0]),a=std::stoi(f[1]),z=std::stoi(f[2]),b=std::stoi(f[3]);require(g.valid(t,z),"OPE support");
        nmhv_physical_equation_t row;row.condition="OPE";row.sample=out.size();row.rhs=rational(f[6]);row.coefficients.resize(c.rows);
        row.ope_coordinate={z,t,0,a,b};
        for(size_t r=0;r<c.rows;++r)for(size_t j=0;j<2;++j)row.coefficients[r]=g.add(row.coefficients[r],g.mul(rational(f[4+j]),values[j][r*g.stride+g.index(t,z,a,b)]));out.push_back(std::move(row));
    }require(!out.empty(),"Empty OPE RHS");return out;
}
}
