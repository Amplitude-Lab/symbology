#pragma once
#include "symrep_io.hpp"
#include "staged_kernel.hpp"

namespace tensor_kernel {
using namespace symrep;
enum class Strategy { streamed, staged };
enum class BoundaryReduction { reduced, raw };
// The equations for X are sum_(b,l,j) F[b,o,l] C[q,l,j] X[b,j]=0.
// They include extension and sewing against any rank-three right recurrence. Sparse
// factors are retained; their potentially enormous product is replayed by row.
class ConstraintRows {
    size_t previous,older,letters,next;
    Mat local;
    std::vector<Vec> slices;
    bool ordered;
public:
    ConstraintRows(const Tensor& expansion,Mat conditions,size_t next_dimension,bool backward=false,bool order_short_rows=false)
        :previous(expansion.dim(0)),older(expansion.dim(backward?2:1)),letters(expansion.dim(backward?1:2)),next(next_dimension),local(std::move(conditions)),slices(older*letters),ordered(order_short_rows){
        require(local.ncol==letters*next,"local condition dimensions mismatch");
        // Row scaling does not change a homogeneous equation. Clear rational
        // denominators before repeated modular conversion and certification.
        for(auto& row:local.rows)vec_cancel_divisor(row);
        std::vector<size_t> counts(slices.size());
        for(size_t p=0;p<expansion.nnz();++p){auto i=expansion.index(p);++counts[i[backward?1:0]*letters+i[backward?0:1]];}
        for(size_t i=0;i<slices.size();++i)slices[i].reserve(counts[i]);
        for(size_t b=0;b<previous;++b)for(size_t p=expansion.rowptr()[b];p<expansion.rowptr()[b+1];++p){auto i=expansion.index(p);slices[i[backward?1:0]*letters+i[backward?0:1]].push_back(index(b),expansion.val(p));}
        // All letters for one older component share a scale: scaling each
        // letter separately would change the equations.
        for(size_t old=0;old<older;++old){int_t denominator=1;
            for(size_t l=0;l<letters;++l)for(auto [b,x]:slices[old*letters+l])denominator=LCM(denominator,x.den());
            if(denominator!=1)for(size_t l=0;l<letters;++l)for(size_t p=0;p<slices[old*letters+l].nnz();++p)slices[old*letters+l][p]*=Q(denominator);
        }
    }
    size_t ncol()const{return previous*next;}
    size_t nrow()const{return older*local.nrow;}
    size_t older_dimension()const{return older;}
    size_t residual_blocks()const{return (older+63)/64;}
    template<class Emit> void residual(const sparse_mat<ulong,I>& k,const field_t& fp,size_t block,Emit&& emit,thread_pool* pool)const{
        require(block<residual_blocks(),"invalid factored residual block");
        modular_residual(k,fp,std::forward<Emit>(emit),pool,block*64,64);
    }
    int_t integral_residual_bound_from_maximum(const int_t& hk)const{
        int_t hf=0,hc=0;
        auto norm=[](const Vec& row){int_t h=0;for(auto [j,x]:row){require(x.is_integer(),"integral residual bound needs integral factors");h+=x.num().abs();}return h;};
        for(const auto& row:slices){auto h=norm(row);if(h>hf)hf=h;}
        for(const auto& row:local.rows){auto h=norm(row);if(h>hc)hc=h;}
        return hf*hc*hk;
    }
    int_t integral_residual_bound(const Mat& kernel)const{
        int_t hk=0;
        for(const auto& row:kernel.rows)for(auto [j,x]:row){require(x.is_integer(),"integral residual bound needs primitive integral candidates");auto h=x.num().abs();if(h>hk)hk=h;}
        return integral_residual_bound_from_maximum(hk);
    }
    // Apply the complete factored operator to a SMALL candidate kernel.
    // Contract the candidate into the recurrence first. Work is proportional
    // to that contraction, not to the huge expanded equation matrix.
    template<class Emit> void modular_residual(const sparse_mat<ulong,I>& k,const field_t& fp,Emit&& emit,thread_pool* pool,size_t first=0,size_t count=SIZE_MAX)const{
        using MV=sparse_vec<ulong,I>;using MM=sparse_mat<ulong,I>;
        require(k.ncol==ncol()&&first<=older,"factored residual dimensions mismatch");size_t end=first+std::min(count,older-first);
        MM by_basis(previous,k.nrow*next);std::vector<size_t> sizes(previous);
        for(const auto& row:k.rows)for(auto c:row.index_span())++sizes[c/next];for(size_t b=0;b<previous;++b)by_basis[b].reserve(sizes[b]);
        for(size_t s=0;s<k.nrow;++s)for(auto [c,x]:k[s])by_basis[c/next].push_back(index(s*next+c%next),x);
        auto convert=[&](const Vec& row){for(auto [j,x]:row)if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("residual factor denominator");return row%fp.mod;};
        std::vector<MV> cm(local.nrow);for(size_t q=0;q<local.nrow;++q)cm[q]=convert(local[q]);
        for(size_t old=first;old<end;){size_t stop=std::min(end,old+64);MM f((stop-old)*letters,previous);
            for(size_t i=0;i<f.nrow;++i)f[i]=convert(slices[old*letters+i]);auto contracted=sparse_mat_mul(f,by_basis,fp,pool);f.clear();
            std::vector<std::vector<MV>> result(stop-old);
            auto evaluate=[&](size_t o){std::vector<ulong> dense(letters*k.nrow*next);
                for(size_t l=0;l<letters;++l)for(auto [j,x]:contracted[o*letters+l])dense[l*k.nrow*next+j]=x;
                for(const auto& equation:cm){MV row;for(size_t s=0;s<k.nrow;++s){ulong value=0;
                        for(auto [lj,a]:equation)value=nmod_add(value,nmod_mul(a,dense[(lj/next)*k.nrow*next+s*next+lj%next],fp.mod),fp.mod);
                        if(value)row.push_back(index(s),value);}
                    if(row.nnz())result[o].push_back(std::move(row));}
            };
            if(pool){pool->detach_loop(size_t(0),stop-old,evaluate);pool->wait();}else for(size_t o=0;o<stop-old;++o)evaluate(o);
            for(auto& block:result)for(auto& row:block)emit(std::move(row));old=stop;
        }
    }
    void profile(structured_kernel::Relations& relations,std::ostream& out)const{
        std::vector<std::vector<size_t>> live(previous);
        for(size_t b=0;b<previous;++b)for(size_t j=0;j<next;++j)if(!relations.is_zero(b*next+j))live[b].push_back(j);
        std::vector<size_t> counts(letters*next);size_t terms=0,nonempty=0,max_terms=0;
        std::map<size_t,std::pair<size_t,size_t>> histogram;
        for(size_t old=0;old<older;++old){std::fill(counts.begin(),counts.end(),0);
            for(size_t l=0;l<letters;++l)for(auto b:slices[old*letters+l].index_span())for(auto j:live[b])++counts[l*next+j];
            for(const auto& equation:local.rows){size_t n=0;for(auto [lj,a]:equation)n+=counts[lj];
                terms+=n;nonempty+=n!=0;max_terms=std::max(max_terms,n);size_t bin=n?std::bit_width(n):0;++histogram[bin].first;histogram[bin].second+=n;
            }
        }
        out<<"PROFILE rows="<<nrow()<<" columns="<<ncol()<<" core_columns="<<relations.dimension()<<" rows_with_terms="<<nonempty<<" terms_before_merge="<<terms<<" max_row_terms="<<max_terms<<" unmerged_fp_bytes="<<terms*(sizeof(ulong)+sizeof(I))<<std::endl;
        for(auto [bin,c]:histogram)out<<"PROFILE_BIN bits="<<bin<<" rows="<<c.first<<" terms="<<c.second<<std::endl;
    }
    // A finite term limit selects an equation subset for experimental rank
    // discovery. Such a subset must NEVER use the full-operator congruence
    // shortcut in verify(); acceptance requires checking every omitted row.
    // Optional per-prime discovery mask: once an original equation has been
    // emitted, later passes need not generate a substituted copy. Reuse the
    // mask only after this enumeration completes (or do a full independent
    // operator check after an early exit).
    template<class Emit> void modular(structured_kernel::Relations& relations,const std::vector<I>& columns,const field_t& fp,Emit&& emit,thread_pool* pool=nullptr,size_t max_terms=SIZE_MAX,size_t min_terms=0,const streaming_echelon::Substitution* substitution=nullptr,std::vector<unsigned char>* seen=nullptr,unsigned sample_shift=0)const{
        require(!seen||seen->size()==nrow(),"equation discovery mask size mismatch");
        require(sample_shift<64,"invalid equation sampling shift");
        using ModVec=sparse_vec<ulong,I>;std::vector<ModVec> fm(slices.size()),cm(local.nrow);
        auto convert=[&](const Vec& row){for(auto [j,x]:row)if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("factor denominator");return row%fp.mod;};
        for(size_t i=0;i<slices.size();++i)fm[i]=convert(slices[i]);
        for(size_t i=0;i<local.nrow;++i)cm[i]=convert(local[i]);
        std::vector<I> mapped(ncol(),-1);std::vector<ulong> scales(ncol());
        for(size_t c=0;c<ncol();++c){auto [r,x]=relations.image(c);if(x==0)continue;if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("relation denominator");mapped[c]=columns[r];scales[c]=x%fp.mod;
            if(substitution&&mapped[c]>=0){I j=mapped[c];mapped[c]=substitution->image[j];scales[c]=nmod_mul(scales[c],substitution->scale[j],fp.mod);}}
        std::vector<size_t> order;
        if(ordered){
            std::vector<size_t> counts(older*letters*next);
            for(size_t old=0;old<older;++old)for(size_t l=0;l<letters;++l)for(auto [b,x]:fm[old*letters+l])for(size_t j=0;j<next;++j)counts[(old*letters+l)*next+j]+=mapped[b*next+j]>=0;
            std::vector<size_t> cost(nrow());order.resize(nrow());std::iota(order.begin(),order.end(),size_t(0));
            for(size_t old=0;old<older;++old)for(size_t q=0;q<cm.size();++q)for(auto [lj,a]:cm[q])cost[old*cm.size()+q]+=counts[old*letters*next+lj];
            std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return cost[a]<cost[b];});
        }
        auto build=[&](size_t sequence){size_t id=ordered?order[sequence]:sequence,old=id/cm.size();if((seen&&(*seen)[id])||!streaming_echelon::sampled_equation(id,sample_shift,UINT64_C(0x9e3779b97f4a7c15)))return ModVec{};const auto& equation=cm[id%cm.size()];ModVec row;
            // Reject long discovery rows using indices only, before modular
            // products, allocation and sorting. Complete residuals are unchanged.
            if(max_terms!=SIZE_MAX){size_t support=0;
                for(auto [lj,a]:equation)for(auto b:fm[old*letters+lj/next].index_span())
                    if(mapped[size_t(b)*next+lj%next]>=0&&++support>max_terms)return ModVec{};
            }
            for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;
                for(auto [b,x]:fm[old*letters+l]){size_t c=b*next+j;if(mapped[c]<0)continue;
                    ulong value=nmod_mul(nmod_mul(a,x,fp.mod),scales[c],fp.mod);if(value)row.push_back(mapped[c],value);if(row.nnz()>max_terms)return ModVec{};}}
            if(row.nnz()<=min_terms){if(seen&&!row.nnz())(*seen)[id]=1;return ModVec{};}row.sort_indices();size_t out=0;
            for(size_t first=0;first<row.nnz();){size_t end=first+1;I c=row(first);ulong value=row[first];
                while(end<row.nnz()&&row(end)==c)value=nmod_add(value,row[end++],fp.mod);
                if(value){row(out)=c;row[out]=value;++out;}first=end;}
            row.resize(out);row.compress();if(seen)(*seen)[id]=1;return row;
        };
        if(!pool||pool->get_thread_count()==1||nrow()<4096){for(size_t sequence=0;sequence<nrow();++sequence){auto row=build(sequence);if(row.nnz())emit(std::move(row));}return;}
        // Parallelize independent rows within a bounded tile; emission remains
        // deterministic and serial, including callbacks that themselves use
        // the pool for elimination. Bound by raw terms before cancellation.
        for(size_t first=0;first<nrow();){size_t end=first,terms=0;std::vector<size_t> costs;
            while(end<nrow()&&end-first<4096){size_t id=ordered?order[end]:end,old=id/cm.size(),cost=0;
                for(auto [lj,a]:cm[id%cm.size()])cost+=fm[old*letters+lj/next].nnz();
                if(max_terms!=SIZE_MAX)cost=std::min(cost,max_terms+1);
                if(end>first&&cost>(2*1024*1024)-std::min(terms,size_t(2*1024*1024)))break;
                terms+=cost;costs.push_back(cost);++end;
            }
            // Condition rows can differ in cost by orders of magnitude. Equal
            // row-count partitions leave most workers idle on long rows.
            size_t target=std::max(size_t(1),(terms+8*pool->get_thread_count()-1)/(8*pool->get_thread_count()));
            std::vector<std::pair<size_t,size_t>> blocks;size_t start=0,weight=0;
            for(size_t i=0;i<costs.size();++i){weight+=costs[i];if(weight>=target){blocks.emplace_back(start,i+1);start=i+1;weight=0;}}
            if(start<costs.size())blocks.emplace_back(start,costs.size());
            std::vector<ModVec> tile(end-first);pool->detach_loop(size_t(0),blocks.size(),[&](size_t b){for(size_t i=blocks[b].first;i<blocks[b].second;++i)tile[i]=build(first+i);},blocks.size());pool->wait();
            for(auto& row:tile)if(row.nnz())emit(std::move(row));first=end;
        }
    }
    // Deterministic certificate. For integral factors and a primitive integral
    // kernel, each residual has absolute value <= ||C_q||_1 max||F_ol||_1
    // max|K|. Zero residues with accumulated modulus larger than that bound prove zero over Q.
    // Nonintegral or large-height inputs use rational arithmetic instead.
    bool verify(const Mat& kernel,rref_option_t opt,bool candidate=false,ulong congruent_prime=0,const int_t* congruent_modulus=nullptr)const{
        require(kernel.ncol==ncol(),"kernel certificate dimension mismatch");
        bool integral=true;int_t hf=0,hc=0,hk=0;
        auto row_height=[&](const Vec& row){int_t h=0;for(auto [j,x]:row){if(!x.is_integer())integral=false;h+=x.num().abs();}return h;};
        for(const auto& row:slices){auto h=row_height(row);if(h>hf)hf=h;}
        for(const auto& row:local.rows){auto h=row_height(row);if(h>hc)hc=h;}
        for(const auto& row:kernel.rows)for(auto [j,x]:row){if(!x.is_integer())integral=false;auto h=x.num().abs();if(h>hk)hk=h;}
        ulong prime=congruent_prime?congruent_prime:n_nextprime(1ULL<<60,0);int_t bound=hf*hc*hk;std::atomic<bool> valid=true;
        // Only a caller holding a modular kernel certificate may supply this
        // prime. The exact relation map has a surviving representative with
        // scale 1 per component; primitive row normalization therefore does
        // not divide by any accepted prime. Congruence survives lifting and scaling.
        const int_t certificate_modulus=congruent_modulus?*congruent_modulus:int_t(prime);
        if(congruent_prime&&integral&&bound<certificate_modulus){
            std::cout<<"exact_certificate=modular_kernel_plus_height bound_bits="<<bound.bits()<<std::endl;return true;
        }
        auto kt=transpose(kernel);
        auto start=std::chrono::steady_clock::now();
        if(integral&&bound<int_t(prime)){
            field_t fp(FIELD_Fp,prime);std::vector<sparse_vec<ulong,I>> fm(slices.size());
            for(size_t i=0;i<slices.size();++i)fm[i]=slices[i]%fp.mod;
            sparse_mat<ulong,I> cm(local.nrow,local.ncol),km(kt.nrow,kt.ncol);
            for(size_t i=0;i<local.nrow;++i)cm[i]=local[i]%fp.mod;
            for(size_t i=0;i<kt.nrow;++i)km[i]=kt[i]%fp.mod;
            opt->pool.detach_loop(size_t(0),older,[&](size_t old){std::vector<ulong> residual(kernel.nrow);
                for(const auto& equation:cm.rows){std::fill(residual.begin(),residual.end(),0);
                    for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;
                        for(auto [b,x]:fm[old*letters+l]){const auto& kr=km[b*next+j];if(!kr.nnz())continue;ulong ax=nmod_mul(a,x,fp.mod);
                            for(auto [s,k]:kr)residual[s]=nmod_add(residual[s],nmod_mul(ax,k,fp.mod),fp.mod);}}
                    for(auto x:residual)if(x)valid=false;
                }});opt->pool.wait();
            std::cout<<"exact_certificate=bounded_modular bound_bits="<<bound.bits();
        }else{
            opt->pool.detach_loop(size_t(0),older,[&](size_t old){std::vector<Q> residual(kernel.nrow);
                for(const auto& equation:local.rows){std::fill(residual.begin(),residual.end(),Q(0));
                    for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;
                        for(auto [b,x]:slices[old*letters+l])for(auto [s,k]:kt[b*next+j])residual[s]+=a*x*k;}
                    for(const auto& x:residual)if(x!=0)valid=false;
                }});opt->pool.wait();std::cout<<"exact_certificate=rational";
        }
        std::cout<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
        if(candidate)return valid;
        require(valid,"nonzero exact residual in streamed kernel");
        require(row_basis(kernel).nrow==kernel.nrow,"streamed kernel rows are dependent");return true;
    }
    template<class Emit> void operator()(structured_kernel::Relations& relations,size_t limit,Emit&& emit)const{
        for(size_t old=0;old<older;++old)for(const auto& equation:local.rows){
            if(limit!=SIZE_MAX){size_t support=0;bool skip=false;
                for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;
                    for(auto b:slices[old*letters+l].index_span())if(!relations.is_zero(b*next+j)&&++support>limit){skip=true;break;}
                    if(skip)break;}
                if(skip)continue;
            }
            Vec row;
            for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;
                for(auto [b,x]:slices[old*letters+l]){size_t c=b*next+j;if(relations.is_zero(c))continue;
                    auto [root,scale]=relations.image(c);row.push_back(index(root),a*x*scale);}}
            normalize(row);emit(std::move(row));
        }
    }
};
inline Mat right_boundary_conditions(const Tensor& d,const Tensor& boundary,BoundaryReduction reduction=BoundaryReduction::reduced){
    auto start=std::chrono::steady_clock::now();
    require(d.rank()==3&&boundary.rank()==3,"right boundary and condition must have rank three");
    require(d.dim(1)==boundary.dim(1),"boundary alphabet mismatch");
    size_t letters=d.dim(0),q=d.dim(2),next=boundary.dim(0),terminal=boundary.dim(2);
    auto lc=flatten(boundary).transpose();Mat local(q*terminal,letters*next);
    sparse_tensor<Q,I,SPARSE_COO> dc(d);
    for(size_t p=0;p<dc.nnz();++p){auto i=dc.index(p);
        for(size_t t=0;t<terminal;++t)for(auto [j,x]:lc[i[1]*terminal+t])
            local[i[2]*terminal+t].push_back(index(i[0]*next+j),dc.val(p)*x);}
    for(auto& row:local.rows)normalize(row);
    std::erase_if(local.rows,[](const auto& row){return !row.nnz();});local.nrow=local.rows.size();
    auto rows=local.nrow,nonzeros=local.nnz();
    if(reduction==BoundaryReduction::reduced)local=row_basis(local);
    std::cout<<"INTERFACE reduction="<<(reduction==BoundaryReduction::raw?"raw":"reduced")
        <<" input_rows="<<rows<<" input_nnz="<<nonzeros<<" rows="<<local.nrow<<" nnz="<<local.nnz()
        <<" elapsed="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
    return local;
}
inline Tensor sew_one(Tensor&& d,Tensor&& f,Tensor&& l,rref_option_t opt,Strategy strategy=Strategy::streamed,const staged_kernel::Options& options={},BoundaryReduction reduction=BoundaryReduction::reduced){
    size_t previous=f.dim(0),last=l.dim(0);auto local=right_boundary_conditions(d,l,reduction);
    ConstraintRows rows(f,std::move(local),last);f.clear();d.clear();l.clear();
    auto k=strategy==Strategy::staged?staged_kernel::solve(staged_kernel::stack(rows),opt,options):structured_kernel::solve_stream(rows.ncol(),rows,opt,nullptr,[&](const Mat& candidate,ulong prime,const int_t& modulus){return rows.verify(candidate,opt,true,prime,&modulus);});
    for(auto& row:k.rows)vec_cancel_divisor(row);k.sort_rows_by_nnz();
    return tensor(k,{k.nrow,previous,last},&opt->pool);
}
// Reduce the side with fewer basis coefficients. The legacy algorithm
// contracts the left side first; the streamed algorithm contracts the right.
// This dimension-based dispatch applies to arbitrary alphabets and terminals.
inline bool prefer_right_boundary(const Tensor& f,const Tensor& l){return l.dim(0)<=f.dim(0);}
inline Mat extension(const Tensor& seed,const Tensor& d,bool backward,rref_option_t opt,structured_kernel::Statistics* stats=nullptr,Strategy strategy=Strategy::streamed,const staged_kernel::Options& options={}){
    require(d.rank()==3&&d.dim(0)==d.dim(1),"extension requires square pair alphabet");size_t letters=d.dim(0);Mat local(d.dim(2),letters*letters);
    sparse_tensor<Q,I,SPARSE_COO> dc(d);for(size_t p=0;p<dc.nnz();++p){auto i=dc.index(p);local[i[2]].push_back(index(i[backward?1:0]*letters+i[backward?0:1]),dc.val(p));}
    for(auto& row:local.rows)normalize(row);ConstraintRows rows(seed,std::move(local),letters,backward);
    auto k=strategy==Strategy::staged?staged_kernel::solve(staged_kernel::stack(rows),opt,options,stats):structured_kernel::solve_stream(rows.ncol(),rows,opt,stats,[&](const Mat& candidate,ulong prime,const int_t& modulus){return rows.verify(candidate,opt,true,prime,&modulus);});
    if(backward)for(auto& row:k.rows){for(size_t p=0;p<row.nnz();++p){size_t c=row(p);row(p)=index((c%letters)*seed.dim(0)+c/letters);}normalize(row);}
    return k;
}
inline Tensor extend_tensor(const Tensor& seed,const Tensor& d,bool backward,rref_option_t opt,Strategy strategy=Strategy::streamed){
    auto k=extension(seed,d,backward,opt,nullptr,strategy);
    for(auto& row:k.rows)vec_cancel_divisor(row);k.sort_rows_by_nnz();
    return tensor(k,backward?std::vector<size_t>{k.nrow,d.dim(0),seed.dim(0)}:std::vector<size_t>{k.nrow,seed.dim(0),d.dim(0)},&opt->pool);
}
}
