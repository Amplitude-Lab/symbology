#pragma once
// Exact rational finite-group representations and kernels over SparseRREF.
// Row convention: a vector transforms as v -> v R(g), and B R(g) = D(g) B.
#include "bootstrap.hpp"
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/nmod_poly_factor.h>
#include <flint/nmod_poly.h>
#include <map>
#include <set>
#include <sstream>
#include <functional>

namespace symrep {
using Q = rat_t;
using I = int32_t;
using Vec = sparse_vec<Q,I>;
using Mat = sparse_mat<Q,I>;
using Tensor = sparse_tensor<Q,I,SPARSE_CSR>;
inline const field_t QQ(FIELD_QQ);
inline void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error("symrep: " + message);
}
inline I index(size_t n) {
    require(n <= size_t(INT32_MAX), "index exceeds int32 capacity");
    return I(n);
}
inline Mat identity(size_t n) {
    Mat a(n,n); for (size_t i=0;i<n;++i) a[i].push_back(index(i),Q(1)); return a;
}
inline Vec packed(const std::map<I,Q>& entries) {
    Vec v; for (const auto& [j,x]:entries) if(x!=0) v.push_back(j,x); return v;
}
inline void normalize(Vec& v) {
    // SparseRREF products and packed rows are already sorted and merged. Do
    // not rebuild a tree of big rationals for every canonical output row.
    bool canonical=true;
    for(size_t p=0;p<v.nnz();++p)if(v[p]==0||(p&&v(p-1)>=v(p))){canonical=false;break;}
    if(canonical)return;
    // Sort and merge in the existing storage. A tree allocated one node per
    // coordinate while retaining the complete unmerged row at the same time.
    v.sort_indices();size_t out=0;
    for(size_t first=0;first<v.nnz();){size_t last=first+1;Q sum=v[first];I col=v(first);
        while(last<v.nnz()&&v(last)==col)sum+=v[last++];
        if(sum!=0){v(out)=col;v[out]=std::move(sum);++out;}first=last;}
    v.resize(out);
}
inline Mat mul(const Mat& a,const Mat& b,thread_pool* pool=nullptr) {
    require(a.ncol==b.nrow,"matrix product dimension mismatch");
    auto c=sparse_mat_mul(a,b,QQ,pool);
    for(auto& row:c.rows) normalize(row);
    return c;
}
inline Vec mul(const Vec& v,const Mat& a) {
    std::map<I,Q> out;
    for(auto [j,x]:v) {require(j>=0&&size_t(j)<a.nrow,"vector index out of range");
        for(auto [k,y]:a[j]) out[k]+=x*y;}
    return packed(out);
}
inline bool equal(const Mat& a,const Mat& b) {
    if(a.nrow!=b.nrow||a.ncol!=b.ncol) return false;
    for(size_t r=0;r<a.nrow;++r) {
        if(a[r].nnz()!=b[r].nnz()) return false;
        for(size_t j=0;j<a[r].nnz();++j)
            if(a[r](j)!=b[r](j)||a[r][j]!=b[r][j]) return false;
    } return true;
}
inline Mat kron(const Mat& a,const Mat& b) {
    index(a.nrow*b.nrow); index(a.ncol*b.ncol);
    Mat c(a.nrow*b.nrow,a.ncol*b.ncol);
    for(size_t i=0;i<a.nrow;++i) for(size_t j=0;j<b.nrow;++j)
        for(auto [k,x]:a[i]) for(auto [l,y]:b[j])
            c[i*b.nrow+j].push_back(index(size_t(k)*b.ncol+l),x*y);
    return c;
}
inline Mat transpose(const Mat& a) {
    Mat out(a.ncol,a.nrow);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])out[j].push_back(index(i),x);
    return out;
}
inline std::string key(const Mat& a) {
    std::ostringstream s; s<<a.nrow<<','<<a.ncol<<';';
    for(size_t i=0;i<a.nrow;++i) for(auto [j,v]:a[i]) s<<i<<','<<j<<'='<<v<<';';
    return s.str();
}
inline Mat rows(const std::vector<Vec>& v,size_t ncol) {
    Mat a(v.size(),ncol); for(size_t i=0;i<v.size();++i) a[i]=v[i]; return a;
}

// Small exact coordinate charts. Large constraint kernels use SparseRREF below.
struct Echelon {
    std::map<I,Vec> pivots;
    Vec reduce(Vec v) const {
        normalize(v);
        for(const auto& [p,row]:pivots) if(auto c=v.find(p); c&&*c!=0) {
            Q coefficient=*c;
            sparse_vec_sub_mul(v,row,coefficient,QQ);
        }
        return v;
    }
    bool add(Vec v) {
        v=reduce(std::move(v)); if(!v.nnz()) return false;
        sparse_vec_rescale(v,Q(1)/v[0],QQ); pivots.emplace(v(0),std::move(v)); return true;
    }
    Mat basis(size_t ncol) const {
        auto p=pivots;
        for(auto it=p.rbegin();it!=p.rend();++it)
            for(auto jt=p.begin();jt!=p.end()&&jt->first<it->first;++jt)
                if(auto v=jt->second.find(it->first);v&&*v!=0) {
                    Q c=*v; sparse_vec_sub_mul(jt->second,it->second,c,QQ);
                }
        Mat b(p.size(),ncol); size_t r=0; for(const auto& [i,v]:p)b[r++]=v; return b;
    }
};
inline Mat row_basis(const Mat& a) { Echelon e; for(auto& v:a.rows)e.add(v); return e.basis(a.ncol); }
inline Mat chart_inverse(const Mat& a,rref_option_t opt){
    if(a.nrow>64)return sparse_mat_inverse(a,QQ,opt);
    require(a.nrow==a.ncol,"coordinate chart must be square");
    if(!a.nrow)return Mat(0,0);
    fmpq_mat_t input,inverse;fmpq_mat_init(input,a.nrow,a.ncol);fmpq_mat_init(inverse,a.nrow,a.ncol);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])fmpq_set(fmpq_mat_entry(input,i,j),x.data());
    bool success=fmpq_mat_inv(inverse,input);Mat result(a.nrow,a.ncol);
    if(success)for(size_t i=0;i<a.nrow;++i)for(size_t j=0;j<a.ncol;++j)if(!fmpq_is_zero(fmpq_mat_entry(inverse,i,j))){Q x;fmpq_set(x.data(),fmpq_mat_entry(inverse,i,j));result[i].push_back(index(j),x);}
    fmpq_mat_clear(input);fmpq_mat_clear(inverse);require(success,"singular coordinate chart");return result;
}
struct Chart {
    Mat basis, inverse;
    std::vector<I> pivots;
    Chart(const Mat& b,rref_option_t opt):basis(b) {
        Echelon e; auto bt=transpose(b);
        // Kernel/RREF bases often expose one private coordinate per row.
        // Prefer those diagonal charts to an arbitrary early dense minor.
        std::vector<I> private_column(b.nrow,-1);
        for(size_t j=0;j<bt.nrow;++j)if(bt[j].nnz()==1&&private_column[bt[j](0)]<0)private_column[bt[j](0)]=index(j);
        for(auto j:private_column)if(j>=0&&e.add(bt[j]))pivots.push_back(j);
        for(size_t j=0;j<bt.nrow&&pivots.size()<b.nrow;++j) if(e.add(bt[j])) pivots.push_back(index(j));
        require(pivots.size()==b.nrow,"dependent basis in coordinate chart");
        Mat square(b.nrow,b.nrow);std::vector<I> column(b.ncol,-1);
        for(size_t j=0;j<pivots.size();++j)column[pivots[j]]=index(j);
        for(size_t i=0;i<b.nrow;++i)for(auto [j,x]:b[i])if(column[j]>=0)square[i].push_back(column[j],x);
        bool diagonal=true;for(size_t i=0;i<square.nrow;++i)if(square[i].nnz()!=1||square[i](0)!=index(i))diagonal=false;
        if(diagonal){inverse=Mat(b.nrow,b.nrow);for(size_t i=0;i<b.nrow;++i)inverse[i].push_back(index(i),Q(1)/square[i][0]);}
        else if(b.nrow) inverse=chart_inverse(square,opt);
        else inverse=Mat(0,0);
    }
    Mat coordinates(const Mat& a,bool verify=true) const {
        require(a.ncol==basis.ncol,"chart ambient dimension mismatch");
        Mat restricted(a.nrow,pivots.size());
        std::vector<I> column(a.ncol,-1);for(size_t j=0;j<pivots.size();++j)column[pivots[j]]=index(j);
        for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])if(column[j]>=0)restricted[i].push_back(column[j],x);
        auto c=mul(restricted,inverse);
        if(verify)require(equal(mul(c,basis),a),"space is not closed under the supplied symmetry");
        return c;
    }
};
inline Mat kernel(Mat a,rref_option_t opt) {
    if(!a.ncol) return Mat(0,0);
    if(!a.nnz()) return identity(a.ncol);
    std::erase_if(a.rows,[](const Vec& row){return !row.nnz();});a.nrow=a.rows.size();
    a.sort_rows_by_nnz();
    auto piv=sparse_mat_rref_reconstruct(a,opt);
    size_t rank=0; for(const auto& p:piv)rank+=p.size();
    if(rank==a.ncol)return Mat(0,a.ncol);
    // Construct the row kernel directly. The library's column kernel followed
    // by transpose holds two copies of every nonzero beside the rational RREF.
    std::vector<I> free(a.ncol,-1),pivot_row(a.ncol,-1);
    for(const auto& batch:piv)for(auto [r,c]:batch){free[c]=-2;pivot_row[c]=r;}
    size_t nullity=0;for(size_t j=0;j<a.ncol;++j)if(free[j]==-1)free[j]=index(nullity++);
    Mat k(nullity,a.ncol);std::vector<size_t> counts(nullity,1);
    for(const auto& batch:piv)for(auto [r,c]:batch)for(auto [j,x]:a[r])if(j!=c&&x!=0){require(free[j]>=0,"kernel requires reduced pivot columns");++counts[free[j]];}
    for(size_t i=0;i<nullity;++i)k[i].reserve(counts[i]);
    // Visit original columns in order, including each private free entry, so
    // every output row is already canonical without a large permutation.
    for(size_t c=0;c<a.ncol;++c){
        if(pivot_row[c]>=0){for(auto [j,x]:a[pivot_row[c]])if(size_t(j)!=c&&x!=0)k[free[j]].push_back(index(c),x);}
        else k[free[c]].push_back(index(c),Q(-1));
    }
    a.clear();
    for(auto& row:k.rows)normalize(row);
    return k;
}
inline Mat short_kernel_basis(const Mat& basis,rref_option_t opt) {
    if(basis.nrow<2||basis.nrow*basis.ncol>4000000)return basis;
    size_t height=0;for(const auto& row:basis.rows)for(auto [j,x]:row)height=std::max(height,size_t(std::max(x.num().bits(),x.den().bits())));
    if(height<8)return basis;
    // Saturation of a large, tall integer lattice can dwarf the actual kernel
    // solve. Shorten bounded independent row batches; their direct sum has
    // exactly the same rational span, without requiring global saturation.
    constexpr size_t batch=128;
    if(basis.nrow>batch) {
        Mat out=basis;
        for(size_t first=0;first<basis.nrow;first+=batch) {
            Mat part(std::min(batch,basis.nrow-first),basis.ncol);
            for(size_t i=0;i<part.nrow;++i)part[i]=basis[first+i];
            auto reduced=short_kernel_basis(part,opt);
            for(size_t i=0;i<part.nrow;++i)out[first+i]=std::move(reduced[i]);
        }return out;
    }
    // Saturate the row lattice first. Primitive RREF rows alone can generate a
    // large-index sublattice, on which LLL cannot recover short integer vectors.
    Mat integral=basis;
    for(auto& row:integral.rows){int_t den=1,divisor=0;for(auto [j,x]:row)den=LCM(den,x.den());
        for(auto [j,x]:row){int_t v=den;v/=x.den();v*=x.num();divisor=GCD(divisor,v);}sparse_vec_rescale(row,Q(den,divisor),QQ);}
    fmpz_mat_t trans,hnf;fmpz_mat_init(trans,basis.ncol,basis.nrow);fmpz_mat_init(hnf,basis.ncol,basis.nrow);
    for(size_t i=0;i<integral.nrow;++i)for(auto [j,x]:integral[i])fmpz_set(fmpz_mat_entry(trans,j,i),fmpq_numref(x.data()));
    // RREF kernel rows have private free-coordinate columns. Their integral
    // values certify an exponent D with D*Z^h contained in the column lattice.
    // FLINT's modular elementary-divisor HNF avoids severe integer growth in
    // the generic tall-matrix HNF routine on these particular lattices.
    auto integral_t=transpose(integral);std::vector<bool> private_column(basis.nrow,false);size_t found=0;int_t exponent=1;
    for(const auto& row:integral_t.rows)if(row.nnz()==1&&!private_column[row(0)]) {
        private_column[row(0)]=true;++found;exponent=LCM(exponent,row[0].num());
    }
    if(found==basis.nrow) {
        if(exponent<0)exponent=-exponent;
        if(exponent==1)fmpz_mat_one(hnf);
        else {fmpz_mat_set(hnf,trans);fmpz_mat_hnf_modular_eldiv(hnf,exponent.data());}
    }else fmpz_mat_hnf(hnf,trans);
    Mat top(basis.nrow,basis.nrow);
    for(size_t i=0;i<basis.nrow;++i)for(size_t j=0;j<basis.nrow;++j)if(!fmpz_is_zero(fmpz_mat_entry(hnf,i,j))){Q x;fmpq_set_fmpz(x.data(),fmpz_mat_entry(hnf,i,j));top[j].push_back(index(i),x);}
    fmpz_mat_clear(trans);fmpz_mat_clear(hnf);
    // top is lower triangular. Solve top*S=integral directly: forming its
    // inverse first can create enormous rational entries that only cancel
    // after multiplication, even though S itself is integral and small.
    Mat saturated(basis.nrow,basis.ncol);
    for(size_t i=0;i<top.nrow;++i) {
        Vec row=integral[i];Q diagonal=0;
        for(auto [j,x]:top[i]) {
            require(size_t(j)<=i,"nontriangular Hermite factor");
            if(size_t(j)==i)diagonal=x;
            else sparse_vec_sub_mul(row,saturated[j],x,QQ);
        }
        require(diagonal!=0,"singular Hermite factor");sparse_vec_rescale(row,Q(1)/diagonal,QQ);saturated[i]=std::move(row);
    }
    fmpz_mat_t lattice,original,u,check;fmpz_mat_init(lattice,basis.nrow,basis.ncol);fmpz_mat_init(original,basis.nrow,basis.ncol);
    fmpz_mat_init(u,basis.nrow,basis.nrow);fmpz_mat_one(u);fmpz_mat_init(check,basis.nrow,basis.ncol);
    for(size_t i=0;i<saturated.nrow;++i)for(auto [j,x]:saturated[i]){require(x.den()==1,"nonintegral saturated kernel lattice");fmpz_set(fmpz_mat_entry(lattice,i,j),fmpq_numref(x.data()));}
    fmpz_mat_set(original,lattice);fmpz_lll_t context;fmpz_lll_context_init(context,0.75,0.51,Z_BASIS,EXACT);fmpz_lll(lattice,u,context);
    fmpz_mat_mul(check,u,original);int_t det;fmpz_mat_det(det.data(),u);
    bool verified=fmpz_mat_equal(check,lattice)&&(det==1||det==-1);
    Mat out(basis.nrow,basis.ncol);
    for(size_t i=0;i<out.nrow;++i)for(size_t j=0;j<out.ncol;++j)if(!fmpz_is_zero(fmpz_mat_entry(lattice,i,j))){Q x;fmpq_set_fmpz(x.data(),fmpz_mat_entry(lattice,i,j));out[i].push_back(index(j),x);}
    fmpz_mat_clear(lattice);fmpz_mat_clear(original);fmpz_mat_clear(u);fmpz_mat_clear(check);
    require(verified,"LLL failed exact lattice transformation verification");return out;
}

struct Group {
    std::vector<std::string> names;
    std::vector<Mat> generators, elements;
    std::vector<std::pair<size_t,size_t>> words; // element = words.first * generator[second]
    std::vector<size_t> generator_ids,inverses;
    std::vector<std::vector<size_t>> table,classes;
    Group(std::vector<std::string> labels,std::vector<Mat> gens,size_t limit=256)
        :names(std::move(labels)),generators(std::move(gens)) {
        require(!generators.empty()&&names.size()==generators.size(),"provide named generator matrices");
        auto n=generators[0].nrow; require(n>0,"empty defining representation");
        std::set<std::string> unique;
        for(size_t i=0;i<generators.size();++i) {
            require(generators[i].nrow==n&&generators[i].ncol==n,"generators must be square with equal dimensions");
            require(!names[i].empty()&&unique.insert(names[i]).second,"duplicate or empty generator name");
            require(row_basis(generators[i]).nrow==n,"singular generator");
        }
        std::map<std::string,size_t> known;
        elements.push_back(identity(n));words.push_back({0,0});known.emplace(key(elements[0]),0);
        std::vector<std::vector<size_t>> edges;
        for(size_t i=0;i<elements.size();++i) {
            std::vector<size_t> edge;
            for(size_t j=0;j<generators.size();++j) {
                auto a=mul(elements[i],generators[j]);auto k=key(a);auto it=known.find(k);
                if(it==known.end()) {
                    require(elements.size()<limit,"group closure exceeds --max-order; group may be infinite");
                    size_t id=elements.size();known.emplace(std::move(k),id);
                    elements.push_back(std::move(a));words.push_back({i,j});edge.push_back(id);
                } else edge.push_back(it->second);
            } edges.push_back(std::move(edge));
        }
        generator_ids=edges[0];size_t order=elements.size();
        table.assign(order,std::vector<size_t>(order));
        for(size_t g=0;g<order;++g) {
            table[g][0]=g;
            for(size_t h=1;h<order;++h)table[g][h]=edges[table[g][words[h].first]][words[h].second];
        }
        inverses.resize(order);
        for(size_t g=0;g<order;++g) {
            auto it=std::find(table[g].begin(),table[g].end(),0);
            require(it!=table[g].end(),"generator closure is not a group");inverses[g]=it-table[g].begin();
        }
        std::set<size_t> used;
        for(size_t g=0;g<order;++g) if(!used.count(g)) {
            std::set<size_t> c;for(size_t h=0;h<order;++h)c.insert(table[table[h][g]][inverses[h]]);
            classes.emplace_back(c.begin(),c.end());used.insert(c.begin(),c.end());
        }
    }
    size_t order() const {return elements.size();}
    std::vector<Mat> representation(const std::vector<Mat>& gens,bool verify=true) const {
        require(gens.size()==names.size(),"wrong number of representation generators");
        size_t n=gens[0].nrow;
        for(const auto& g:gens)require(g.nrow==n&&g.ncol==n,"representation generator shape mismatch");
        std::vector<Mat> rep;rep.push_back(identity(n));
        for(size_t i=1;i<order();++i)rep.push_back(mul(rep[words[i].first],gens[words[i].second]));
        if(verify)for(size_t i=0;i<order();++i)for(size_t j=0;j<gens.size();++j)
            require(equal(mul(rep[i],gens[j]),rep[table[i][generator_ids[j]]]),"representation violates a group relation");
        return rep;
    }
    std::vector<Q> product(const std::vector<Q>& a,const std::vector<Q>& b) const {
        std::vector<Q> c(order(),Q(0));
        for(size_t g=0;g<order();++g)if(a[g]!=0)
            for(size_t h=0;h<order();++h)if(b[h]!=0)c[table[g][h]]+=a[g]*b[h];
        return c;
    }
    Vec translate(const std::vector<Q>& a,size_t g) const {
        std::map<I,Q> x;for(size_t h=0;h<order();++h)if(a[h]!=0)x[index(table[h][g])]+=a[h];return packed(x);
    }
    Mat evaluate(const std::vector<Q>& a,const std::vector<Mat>& rep) const {
        Mat out(rep[0].nrow,rep[0].ncol);
        for(size_t i=0;i<out.nrow;++i) {
            std::map<I,Q> x;
            for(size_t g=0;g<order();++g)if(a[g]!=0)for(auto [j,y]:rep[g][i])x[j]+=a[g]*y;
            out[i]=packed(x);
        } return out;
    }
};

struct Poly {
    fmpq_poly_t p;
    Poly(){fmpq_poly_init(p);} ~Poly(){fmpq_poly_clear(p);}
    Poly(const Poly& q):Poly(){fmpq_poly_set(p,q.p);}
    Poly& operator=(const Poly& q){fmpq_poly_set(p,q.p);return *this;}
};
inline std::vector<Poly> factors(const Mat& a) {
    fmpq_mat_t m;fmpq_mat_init(m,a.nrow,a.ncol);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])fmpq_set(fmpq_mat_entry(m,i,j),x.data());
    Poly cp;fmpq_mat_charpoly(cp.p,m);fmpq_mat_clear(m);
    fmpz_poly_t numerator;fmpz_poly_init(numerator);fmpq_poly_get_numerator(numerator,cp.p);
    fmpz_poly_factor_t fac;fmpz_poly_factor_init(fac);fmpz_poly_factor(fac,numerator);
    std::vector<Poly> result;
    for(slong i=0;i<fac->num;++i){Poly p;fmpq_poly_set_fmpz_poly(p.p,fac->p+i);fmpq_poly_make_monic(p.p,p.p);result.push_back(p);}
    fmpz_poly_factor_clear(fac);fmpz_poly_clear(numerator);return result;
}
inline std::vector<Q> evaluate(const Group& g,const Poly& p,const std::vector<Q>& z) {
    std::vector<Q> out(g.order(),Q(0));
    for(slong i=fmpq_poly_degree(p.p);i>=0;--i) {
        out=g.product(out,z);Q c;fmpq_poly_get_coeff_fmpq(c.data(),p.p,i);out[0]+=c;
    }return out;
}
struct Model {
    size_t degree=0,division_dimension=0,dimension=0;
    std::vector<Q> central,primitive,character;
    std::vector<size_t> orbit;
    std::vector<Mat> actions; // every group element in the canonical rational irrep
};
inline std::vector<std::vector<Q>> central_idempotents(const Group& g,std::vector<size_t>& degrees) {
    // A separating central element identifies the simple rational factors.
    // Exact centre-dimension checks certify separation; no numerical eigenvectors.
    for(size_t attempt=0;attempt<64;++attempt) {
        std::vector<Q> z(g.order(),Q(0));
        uint64_t rng=0x9e3779b97f4a7c15ULL+attempt;
        for(size_t c=0;c<g.classes.size();++c) {
            rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;
            long weight=attempt==0?long(c+1):long(rng%101)-50;
            for(auto h:g.classes[c])z[h]=Q(weight);
        }
        Mat regular(g.order(),g.order());
        for(size_t a=0;a<g.order();++a){std::map<I,Q> row;
            for(size_t h=0;h<g.order();++h)if(z[h]!=0)row[index(g.table[a][h])]+=z[h];regular[a]=packed(row);}
        auto fs=factors(regular);Poly total;fmpq_poly_one(total.p);
        for(const auto& f:fs)fmpq_poly_mul(total.p,total.p,f.p);
        std::vector<std::vector<Q>> result;degrees.clear();bool separated=true;
        for(const auto& f:fs) {
            Poly quotient,remainder,d,s,t,e;
            fmpq_poly_divrem(quotient.p,remainder.p,total.p,f.p);
            fmpq_poly_xgcd(d.p,s.p,t.p,quotient.p,f.p);
            require(fmpq_poly_is_one(d.p),"central factors are not coprime");
            fmpq_poly_mul(e.p,quotient.p,s.p);fmpq_poly_rem(e.p,e.p,total.p);
            auto coeff=evaluate(g,e,z);Echelon center;
            for(const auto& cl:g.classes) {
                std::vector<Q> sum(g.order(),Q(0));for(auto h:cl)sum[h]=Q(1);
                center.add(g.translate(g.product(coeff,sum),0));
            }
            size_t degree=fmpq_poly_degree(f.p);
            if(center.pivots.size()!=degree){separated=false;break;}
            require(g.product(coeff,coeff)==coeff,"central projector is not idempotent");
            result.push_back(std::move(coeff));degrees.push_back(degree);
        }
        if(separated) {
            std::vector<Q> sum(g.order(),Q(0));for(const auto& p:result)for(size_t i=0;i<g.order();++i)sum[i]+=p[i];
            require(sum[0]==1&&std::all_of(sum.begin()+1,sum.end(),[](const Q& x){return x==0;}),"incomplete rational central decomposition");
            return result;
        }
    }
    throw std::runtime_error("symrep: failed to find a separating rational central element");
}
inline Model model_from_projector(const Group& g,const std::vector<Q>& p,size_t degree,rref_option_t opt) {
    Model m;m.degree=degree;m.primitive=p;Echelon rank;std::vector<Vec> basis;
    for(size_t h=0;h<g.order();++h) {auto v=g.translate(p,h);if(rank.add(v)){m.orbit.push_back(h);basis.push_back(v);}}
    m.dimension=basis.size();if(!m.dimension)return m;
    auto b=rows(basis,g.order());Chart chart(b,opt);
    // Generator intertwiners prove closure for every group word. Reuse them
    // instead of resolving the same small coordinate chart for every element.
    std::vector<Mat> defining;
    for(auto h:g.generator_ids) {
        Mat image(b.nrow,b.ncol);
        for(size_t i=0;i<b.nrow;++i){std::map<I,Q> v;for(auto [j,x]:b[i])v[index(g.table[j][h])]+=x;image[i]=packed(v);}
        defining.push_back(chart.coordinates(image));
    }
    for(size_t h=0;h<g.order();++h) {
        m.actions.push_back(h?mul(m.actions[g.words[h].first],defining[g.words[h].second]):identity(m.dimension));Q tr=0;
        for(size_t i=0;i<b.nrow;++i)if(auto v=m.actions.back()[i].find(index(i)))tr+=*v;
        m.character.push_back(tr);
    }
    Q norm=0,fs=0;for(size_t h=0;h<g.order();++h){norm+=m.character[h]*m.character[g.inverses[h]];fs+=m.character[g.table[h][h]];}
    norm/=Q(g.order());fs/=Q(g.order());
    // Schur-index one, or a certified quaternionic Schur-index two component.
    if(norm==Q(degree))m.division_dimension=degree;
    else if(norm==Q(4*degree)&&fs==Q(-2*long(degree)))m.division_dimension=4*degree;
    return m;
}
inline std::vector<Model> models(const Group& g,rref_option_t opt,size_t subgroup_limit=1024) {
    std::vector<size_t> degrees;auto centers=central_idempotents(g,degrees);
    // Try fixed-space projectors for subgroups, constructed from the generators.
    // This is not a parity/dihedral dispatch. Uncertified models are never used.
    std::vector<std::vector<size_t>> subgroups{{0}};std::set<std::vector<size_t>> seen{{0}};
    for(size_t i=0;i<subgroups.size()&&subgroups.size()<subgroup_limit;++i)
        for(size_t h=1;h<g.order()&&subgroups.size()<subgroup_limit;++h) {
            if(std::binary_search(subgroups[i].begin(),subgroups[i].end(),h))continue;
            std::vector<size_t> generators=subgroups[i];generators.push_back(h);
            std::set<size_t> closure{0};std::vector<size_t> todo{0};
            for(size_t a=0;a<todo.size();++a)for(auto b:generators)if(closure.insert(g.table[todo[a]][b]).second)todo.push_back(g.table[todo[a]][b]);
            std::vector<size_t> hs(closure.begin(),closure.end());if(seen.insert(hs).second)subgroups.push_back(std::move(hs));
        }
    std::vector<Model> out;
    for(size_t s=0;s<centers.size();++s) {
        bool found=false;
        for(const auto& hs:subgroups) {
            std::vector<Q> avg(g.order(),Q(0));for(auto h:hs)avg[h]=Q(1,long(hs.size()));
            auto p=g.product(centers[s],avg);auto m=model_from_projector(g,p,degrees[s],opt);
            if(m.division_dimension) {
                m.central=centers[s];out.push_back(std::move(m));found=true;break;
            }
        }
        require(found,"could not certify a rational irreducible model from subgroup projectors; supply a different defining representation is insufficient (additional model construction is needed for this group)");
    }
    size_t total=0;for(const auto& m:out)total+=m.dimension*m.dimension/m.division_dimension;
    require(total==g.order(),"rational irreducible models do not account for the regular representation");
    return out;
}

struct Layout {
    std::vector<size_t> copies;
    size_t dimension(const std::vector<Model>& ms) const {
        require(copies.size()==ms.size(),"layout/model mismatch");size_t n=0;
        for(size_t i=0;i<ms.size();++i)n+=copies[i]*ms[i].dimension;return n;
    }
    std::vector<size_t> offsets(const std::vector<Model>& ms) const {
        std::vector<size_t> out;size_t n=0;for(size_t i=0;i<ms.size();++i){out.push_back(n);n+=copies[i]*ms[i].dimension;}return out;
    }
};
struct Adapted {Mat basis;Layout layout;};

// Select whole simple modules in a primitive-idempotent image. p Q[G] is the
// fixed model: p*g_i -> v*R(g_i) gives the SAME action for every selected copy.
inline Mat copies_from_corner(const Mat& h,const Mat& q,const Model& m,
    const std::function<Vec(const Vec&,size_t)>& act,
    const std::function<Vec(const Vec&)>& project,rref_option_t opt,
    const std::vector<I>* kernel_coordinates=nullptr) {
    // Over a split rational irrep the corner is one-dimensional. Independent
    // corner vectors already certify independent whole copies: re-eliminating
    // their projected orbits only repeats the kernel rank calculation.
    if(m.division_dimension==1) {
        std::vector<Vec> basis;
        for(const auto& row:h.rows) {
            auto v=mul(row,q);std::vector<Vec> one;int_t den=1,divisor=0;
            for(auto word:m.orbit)one.push_back(act(v,word));
            for(const auto& x:one)for(auto [j,c]:x)den=LCM(den,c.den());
            for(const auto& x:one)for(auto [j,c]:x){int_t value=den;value/=c.den();value*=c.num();divisor=GCD(divisor,value);}
            require(divisor!=0,"zero irreducible orbit");
            for(auto& x:one){sparse_vec_rescale(x,Q(den,divisor),QQ);basis.push_back(std::move(x));}
        }return rows(basis,q.ncol);
    }
    Chart chart(q,opt);Echelon covered;std::vector<Vec> basis;
    // Restriction to independent coordinates is injective on the kernel.
    // Rank tests can therefore use nullity columns instead of all candidate
    // columns. The caller obtains this chart BEFORE optional lattice changes.
    std::vector<I> coordinate(h.ncol,-1);
    if(kernel_coordinates)for(size_t j=0;j<kernel_coordinates->size();++j)coordinate.at(kernel_coordinates->at(j))=index(j);
    auto restrict=[&](const Vec& v){if(!kernel_coordinates)return v;Vec out;
        for(auto [j,x]:v)if(coordinate[j]>=0)out.push_back(coordinate[j],x);normalize(out);return out;};
    for(const auto& row:h.rows) {
        auto residual=covered.reduce(restrict(row));if(!residual.nnz())continue;
        auto v=mul(row,q);std::vector<Vec> one;
        for(auto word:m.orbit)one.push_back(act(v,word));
        // One scale for the entire irreducible copy, never per component.
        int_t den=1,divisor=0;
        for(const auto& x:one)for(auto [j,c]:x)den=LCM(den,c.den());
        for(const auto& x:one)for(auto [j,c]:x){int_t value=den;value/=c.den();value*=c.num();divisor=GCD(divisor,value);}
        require(divisor!=0,"zero irreducible orbit");Q scale(den,divisor);
        for(auto& x:one)sparse_vec_rescale(x,scale,QQ);
        Echelon rank;for(const auto& x:one)rank.add(x);
        require(rank.pivots.size()==m.dimension,"invalid primitive-projector orbit");
        Mat projected(one.size(),q.ncol);for(size_t i=0;i<one.size();++i)projected[i]=project(one[i]);
        auto coords=chart.coordinates(projected);size_t before=covered.pivots.size();
        for(auto& x:coords.rows)covered.add(restrict(x));
        require(covered.pivots.size()-before==m.division_dimension,"irreducible copy overlap or invalid division algebra dimension");
        basis.insert(basis.end(),one.begin(),one.end());
    }
    require(covered.pivots.size()==h.nrow,"corner kernel basis has wrong rank");
    return rows(basis,q.ncol);
}
inline Adapted adapt(const Group& g,const std::vector<Model>& ms,const std::vector<Mat>& gens,rref_option_t opt) {
    auto rep=g.representation(gens);Adapted out;out.layout.copies.resize(ms.size());std::vector<Vec> all;
    for(size_t s=0;s<ms.size();++s) {
        auto p=g.evaluate(ms[s].primitive,rep);auto q=row_basis(p);
        auto b=copies_from_corner(identity(q.nrow),q,ms[s],
            [&](const Vec& v,size_t h){return mul(v,rep[h]);},[&](const Vec& v){return mul(v,p);},opt);
        out.layout.copies[s]=b.nrow/ms[s].dimension;all.insert(all.end(),b.rows.begin(),b.rows.end());
    }
    out.basis=rows(all,gens[0].nrow);
    require(out.basis.nrow==out.basis.ncol,"adaptation did not span the input representation");
    return out;
}
inline Mat action(const Layout& layout,const std::vector<Model>& ms,size_t h) {
    size_t n=layout.dimension(ms);Mat out(n,n);size_t off=0;
    for(size_t s=0;s<ms.size();++s)for(size_t a=0;a<layout.copies[s];++a){
        for(size_t i=0;i<ms[s].dimension;++i)for(auto [j,v]:ms[s].actions[h][i])out[off+i].push_back(index(off+j),v);
        off+=ms[s].dimension;
    }return out;
}
inline std::vector<Mat> generators(const Group& g,const std::vector<Model>& ms,const Layout& layout) {
    std::vector<Mat> out;for(auto h:g.generator_ids)out.push_back(action(layout,ms,h));return out;
}
// Apply an irreducible block action without constructing the repeated matrix.
// Input/output are row coefficients in the logical adapted basis. Each copy
// uses only a small local accumulator, independent of total representation size.
inline Vec apply_irrep_action(const Vec& v,const Layout& layout,const std::vector<Model>& ms,size_t h){
    size_t n=layout.dimension(ms);require(!v.nnz()||(v(0)>=0&&size_t(v(v.nnz()-1))<n),"irrep action vector dimension mismatch");
    Vec out;size_t position=0,offset=0;
    for(size_t s=0;s<ms.size();++s){size_t d=ms[s].dimension,end=offset+layout.copies[s]*d;
        require(h<ms[s].actions.size(),"unknown irrep action element");const auto& small=ms[s].actions[h];std::vector<Q> values(d);
        while(position<v.nnz()&&size_t(v(position))<end){size_t begin=offset+((size_t(v(position))-offset)/d)*d;
            std::fill(values.begin(),values.end(),Q(0));
            while(position<v.nnz()&&size_t(v(position))<begin+d){size_t i=size_t(v(position))-begin;
                for(auto [j,x]:small[i])values[j]+=v[position]*x;++position;}
            for(size_t j=0;j<d;++j)if(values[j]!=0)out.push_back(index(begin+j),values[j]);
        }offset=end;
    }return out;
}
inline void verify_adaptation(const Adapted& a,const Group& g,const std::vector<Model>& ms,const std::vector<Mat>& gens) {
    require(a.basis.nrow==a.layout.dimension(ms),"adapted dimension mismatch");
    for(size_t i=0;i<gens.size();++i)require(equal(mul(a.basis,gens[i]),mul(action(a.layout,ms,g.generator_ids[i]),a.basis)),"adapted basis fails exact intertwiner check");
}
inline std::vector<Mat> induced(const Mat& basis,const std::vector<Mat>& ambient,rref_option_t opt) {
    Chart chart(basis,opt);std::vector<Mat> result;
    for(const auto& a:ambient)result.push_back(chart.coordinates(mul(basis,a)));
    return result;
}
inline void verify_equivariance(const Mat& a,const std::vector<Mat>& domain,const std::vector<Mat>& codomain) {
    require(domain.size()==codomain.size(),"equation generator count mismatch");
    for(size_t i=0;i<domain.size();++i)
        require(equal(mul(a,transpose(domain[i])),mul(transpose(codomain[i]),a)),"constraints do not intertwine the supplied domain/equation actions");
}
inline std::vector<Mat> equation_actions(const Mat& a,const std::vector<Mat>& domain,rref_option_t opt) {
    // Redundant equations must first be replaced by an independent row basis.
    Chart chart(a,opt);std::vector<Mat> out;
    for(const auto& r:domain)out.push_back(transpose(chart.coordinates(mul(a,transpose(r)))));return out;
}

struct KernelResult {Adapted solution;std::vector<size_t> corner_columns,corner_nullities;};
struct KernelTimings {double rref=0,shorten=0,lift=0;};
inline KernelResult kernel_blocks(size_t columns,const Group& g,const std::vector<Model>& ms,const Layout& source,
    const std::function<Mat(size_t,const Mat&)>& reduced_matrix,rref_option_t opt,KernelTimings* timings=nullptr) {
    require(columns==source.dimension(ms),"kernel column/layout mismatch");
    KernelResult out;out.solution.layout.copies.resize(ms.size());auto offsets=source.offsets(ms);std::vector<Vec> all;
    for(size_t s=0;s<ms.size();++s) {
        const auto& m=ms[s];auto p=g.evaluate(m.primitive,m.actions);auto local=row_basis(p);
        require(local.nrow==m.division_dimension,"primitive corner has wrong rational dimension");
        Mat q(source.copies[s]*local.nrow,columns);
        for(size_t c=0;c<source.copies[s];++c)for(size_t j=0;j<local.nrow;++j)
            for(auto [k,v]:local[j])q[c*local.nrow+j].push_back(index(offsets[s]+c*m.dimension+k),v);
        auto reduced=reduced_matrix(s,q);require(reduced.ncol==q.nrow,"reduced block column mismatch");
        auto tick=std::chrono::steady_clock::now();auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();};
        auto h=kernel(std::move(reduced),opt);if(timings)timings->rref+=elapsed();tick=std::chrono::steady_clock::now();
        if(std::getenv("SYMREP_TRACE"))std::cerr<<"corner_kernel sector="<<s<<" nullity="<<h.nrow<<" columns="<<h.ncol<<std::endl;
        // Split irreps need scalar multiplicities only; primitive normalization
        // below suffices. Saturating/LLL-reducing every kernel can cost more
        // than SparseRREF and often makes a sparse recurrence denser.
        std::vector<I> kernel_coordinates;
        if(m.division_dimension>1){
            // A SparseRREF kernel has a private free coordinate for each row.
            // These columns give a certified injective chart without another
            // elimination or a potentially dense inverse.
            auto ht=transpose(h);kernel_coordinates.assign(h.nrow,-1);
            for(size_t j=0;j<ht.nrow;++j)if(ht[j].nnz()==1&&kernel_coordinates[ht[j](0)]<0)kernel_coordinates[ht[j](0)]=index(j);
            require(std::find(kernel_coordinates.begin(),kernel_coordinates.end(),I(-1))==kernel_coordinates.end(),"kernel lacks independent free coordinates");
            h=short_kernel_basis(h,opt);
        }
        if(timings)timings->shorten+=elapsed();tick=std::chrono::steady_clock::now();
        out.corner_columns.push_back(q.nrow);out.corner_nullities.push_back(h.nrow);
        auto apply_local=[&](const Vec& v,const Mat& r) {
            std::map<I,Q> entries;
            for(auto [j,x]:v) {
                require(size_t(j)>=offsets[s]&&size_t(j)<offsets[s]+source.copies[s]*m.dimension,"mixed irrep corner");
                size_t copy=(j-offsets[s])/m.dimension,component=(j-offsets[s])%m.dimension;
                for(auto [k,y]:r[component])entries[index(offsets[s]+copy*m.dimension+k)]+=x*y;
            }return packed(entries);
        };
        auto b=copies_from_corner(h,q,m,[&](const Vec& v,size_t word){return apply_local(v,m.actions[word]);},
            [&](const Vec& v){return apply_local(v,p);},opt,m.division_dimension>1?&kernel_coordinates:nullptr);
        out.solution.layout.copies[s]=b.nrow/m.dimension;all.insert(all.end(),b.rows.begin(),b.rows.end());
        if(timings)timings->lift+=elapsed();
    }
    out.solution.basis=rows(all,columns);
    return out;
}
inline KernelResult kernel_adapted(const Mat& a,const Group& g,const std::vector<Model>& ms,const Layout& source,rref_option_t opt) {
    // Residuals alone cannot detect a non-invariant kernel that misses the
    // chosen primitive corner. Establish invariance before reducing columns.
    require(a.ncol==source.dimension(ms),"kernel column/layout mismatch");
    auto independent=row_basis(a);auto domain=generators(g,ms,source);
    auto codomain=equation_actions(independent,domain,opt);verify_equivariance(independent,domain,codomain);
    auto out=kernel_blocks(a.ncol,g,ms,source,[&](size_t,const Mat& q){return mul(a,transpose(q),&opt->pool);},opt);
    require(!mul(a,out.solution.basis.transpose(),&opt->pool).nnz(),"nonzero exact kernel residual (constraints may break symmetry)");
    return out;
}

// Tensor CG changes are direct sums of small invertible matrices, interleaved
// by irrep. Invert connected support blocks, never a weight-sized dense chart.
inline Mat block_inverse(const Mat& a,rref_option_t opt) {
    require(a.nrow==a.ncol,"non-square tensor-product basis");size_t n=a.nrow;
    std::vector<size_t> parent(2*n);std::iota(parent.begin(),parent.end(),0);
    std::function<size_t(size_t)> root=[&](size_t i){return parent[i]==i?i:parent[i]=root(parent[i]);};
    for(size_t i=0;i<n;++i)for(auto [j,x]:a[i])parent[root(n+j)]=root(i);
    std::map<size_t,std::pair<std::vector<size_t>,std::vector<size_t>>> blocks;
    for(size_t i=0;i<n;++i){blocks[root(i)].first.push_back(i);blocks[root(n+i)].second.push_back(i);}
    Mat out(n,n);std::map<std::string,Mat> cache;
    for(const auto& [id,rc]:blocks){const auto& [rs,cs]=rc;require(rs.size()==cs.size(),"singular support block");
        std::map<size_t,size_t> col;for(size_t j=0;j<cs.size();++j)col[cs[j]]=j;
        Mat local(rs.size(),cs.size());for(size_t i=0;i<rs.size();++i)for(auto [j,x]:a[rs[i]])local[i].push_back(index(col.at(j)),x);
        auto k=key(local);auto it=cache.find(k);if(it==cache.end())it=cache.emplace(k,sparse_mat_inverse(local,QQ,opt)).first;
        for(size_t i=0;i<cs.size();++i)for(auto [j,x]:it->second[i])out[cs[i]].push_back(index(rs[j]),x);
    }
    for(auto& row:out.rows)normalize(row);return out;
}

struct ProductCopy {size_t left,right,left_copy,right_copy,channel;};
using ProductCopies=std::vector<std::vector<ProductCopy>>;
class TensorAdapter {
    const Group& group;const std::vector<Model>& model;rref_option* options;
    std::map<std::pair<size_t,size_t>,Adapted> cache;
public:
    TensorAdapter(const Group& g,const std::vector<Model>& m,rref_option_t opt):group(g),model(m),options(opt){}
    Adapted product(const Layout& left,const Layout& right,ProductCopies* channels=nullptr,bool materialize=true) {
        if(channels)channels->assign(model.size(),{});
        auto lo=left.offsets(model),ro=right.offsets(model);size_t nr=right.dimension(model);
        std::vector<std::vector<Vec>> by_irrep(model.size());Adapted out;out.layout.copies.resize(model.size());
        for(size_t a=0;a<model.size();++a)if(left.copies[a])for(size_t b=0;b<model.size();++b)if(right.copies[b]) {
            auto key=std::make_pair(a,b);auto it=cache.find(key);
            if(it==cache.end()) {
                std::vector<Mat> gens;for(auto h:group.generator_ids)gens.push_back(kron(model[a].actions[h],model[b].actions[h]));
                auto adapted=adapt(group,model,gens,options);verify_adaptation(adapted,group,model,gens);
                it=cache.emplace(key,std::move(adapted)).first;
            }
            const auto& cg=it->second;auto co=cg.layout.offsets(model);
            for(size_t s=0;s<model.size();++s)for(size_t x=0;x<left.copies[a];++x)for(size_t y=0;y<right.copies[b];++y) {
                out.layout.copies[s]+=cg.layout.copies[s];
                if(channels)for(size_t k=0;k<cg.layout.copies[s];++k)(*channels)[s].push_back({a,b,x,y,k});
                if(materialize)for(size_t r=0;r<cg.layout.copies[s]*model[s].dimension;++r) {
                    Vec v;
                    for(auto [j,c]:cg.basis[co[s]+r]) {
                        size_t i1=lo[a]+x*model[a].dimension+j/model[b].dimension;
                        size_t i2=ro[b]+y*model[b].dimension+j%model[b].dimension;
                        v.push_back(index(i1*nr+i2),c);
                    }by_irrep[s].push_back(std::move(v));
                }
            }
        }
        std::vector<Vec> all;for(auto& sector:by_irrep)all.insert(all.end(),sector.begin(),sector.end());
        if(materialize)out.basis=rows(all,left.dimension(model)*nr);
        require(out.layout.dimension(model)==left.dimension(model)*nr,"incomplete tensor-product adaptation");return out;
    }
};
} // namespace symrep
