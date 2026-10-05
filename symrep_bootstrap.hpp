#pragma once
#include "symrep_io.hpp"
#include "tensor_kernel.hpp"
#include <optional>
#include <array>

namespace symrep {
inline Mat condition_rows(const Tensor& d) {
    require(d.rank()==3&&d.dim(0)==d.dim(1),"condition tensor must have shape {letters,letters,equations}");
    Mat a(d.dim(2),d.dim(0)*d.dim(1));sparse_tensor<Q,I,SPARSE_COO> c(d);
    for(size_t p=0;p<c.nnz();++p)a[c.index(p)[2]].push_back(index(c.index(p)[0]*d.dim(1)+c.index(p)[1]),c.val(p));
    for(auto& r:a.rows)normalize(r);return a;
}
inline Tensor condition_tensor(const Mat& a,size_t letters) {
    require(a.ncol==letters*letters,"condition column count mismatch");
    sparse_tensor<Q,I,SPARSE_COO> c(std::vector<size_t>{letters,letters,a.nrow});
    for(size_t q=0;q<a.nrow;++q)for(auto [j,x]:a[q])c.push_back({index(j/letters),index(j%letters),index(q)},x);
    c.sort_indices();return Tensor(std::move(c));
}
inline Mat assemble(const Tensor& expansion,const Tensor& d,bool backward) {
    require(expansion.rank()==3&&d.rank()==3,"extension requires rank-three tensors");
    size_t letters=d.dim(0),previous=expansion.dim(0),older=expansion.dim(backward?2:1),q=d.dim(2);
    require(d.dim(1)==letters&&expansion.dim(backward?1:2)==letters,"letter-axis mismatch");
    std::vector<std::vector<std::pair<std::pair<I,I>,Q>>> lookup(letters);
    sparse_tensor<Q,I,SPARSE_COO> dc(d),ec(expansion);
    for(size_t p=0;p<dc.nnz();++p) {
        auto i=dc.index(p);size_t known=i[backward?1:0],next=i[backward?0:1];
        lookup[known].push_back({{index(next),i[2]},dc.val(p)});
    }
    Mat a(older*q,previous*letters);
    for(size_t p=0;p<ec.nnz();++p) {
        auto i=ec.index(p);size_t b=i[0],o=i[backward?2:1],l=i[backward?1:2];
        for(const auto& [where,value]:lookup[l]) {
            size_t row=backward?where.second*older+o:o*q+where.second;
            size_t col=backward?where.first*previous+b:b*letters+where.first;
            a[row].push_back(index(col),ec.val(p)*value);
        }
    }
    for(auto& r:a.rows)normalize(r);return a;
}
inline Mat assemble_reduced(const Tensor& expansion,const Tensor& d,bool backward,
    const Mat& left,const Mat& right,rref_option_t opt) {
    const size_t letters=d.dim(0),q=d.dim(2),older=expansion.dim(backward?2:1),previous=expansion.dim(0);
    require(left.ncol==older*q&&right.ncol==previous*letters,"reduced assembly chart shape mismatch");
    struct Term{size_t basis,letter;Q value;};std::vector<std::vector<Term>> ex(older);
    std::vector<std::vector<std::pair<I,Q>>> equations(q*letters);
    sparse_tensor<Q,I,SPARSE_COO> ec(expansion),dc(d);
    for(size_t p=0;p<ec.nnz();++p){auto i=ec.index(p);ex[i[backward?2:1]].push_back({size_t(i[0]),size_t(i[backward?1:2]),ec.val(p)});}
    for(size_t p=0;p<dc.nnz();++p){auto i=dc.index(p);equations[i[2]*letters+i[backward?1:0]].push_back({i[backward?0:1],dc.val(p)});}
    auto rt=transpose(right);Mat out(left.nrow,right.nrow);
    opt->pool.detach_loop(size_t(0),left.nrow,[&](size_t r){
        std::map<I,Q> sum;
        for(auto [raw,lvalue]:left[r]) {
            size_t old=backward?raw%older:raw/q,condition=backward?raw/older:raw%q;
            for(const auto& term:ex[old])for(const auto& [letter,value]:equations[condition*letters+term.letter]) {
                size_t col=backward?letter*previous+term.basis:term.basis*letters+letter;
                Q x=lvalue*term.value*value;
                for(auto [j,c]:rt[col])sum[j]+=x*c;
            }
        }out[r]=packed(sum);
    });opt->pool.wait();return out;
}
struct Prepared {
    Layout alphabet,equations,seed,terminal;Tensor condition,expansion;bool backward=false;
    std::optional<Mat> reduced_condition,condition_pair_basis;
};

// An equivariant condition map has one End_G(U)-valued coefficient per pair
// of equivalent copies. All its component equations are generated from these
// coefficients. In row conventions D*CG^T consists of transposed commutants.
class EquivariantConditions {
public:
    struct Entry {size_t equation,copy,left,left_copy,right,right_copy,channel,end;Q value;};
    Adapted pair;ProductCopies channels;Mat coefficients;
    std::vector<std::vector<Mat>> endomorphisms;
    std::vector<Entry> entries;
    EquivariantConditions(const Prepared& p,const Group& g,const std::vector<Model>& ms,TensorAdapter& adapter,rref_option_t opt) {
        pair=adapter.product(p.alphabet,p.alphabet,&channels);auto po=pair.layout.offsets(ms),eo=p.equations.offsets(ms);
        if(p.condition_pair_basis)require(equal(pair.basis,*p.condition_pair_basis),"saved condition coupling convention differs from the current decomposition");
        size_t nr=0,nc=0;for(size_t s=0;s<ms.size();++s){nr+=p.equations.copies[s];nc+=pair.layout.copies[s]*ms[s].division_dimension;}
        coefficients=p.reduced_condition?*p.reduced_condition:Mat(nr,nc);
        require(coefficients.nrow==nr&&coefficients.ncol==nc,"reduced condition dimensions disagree with copy layouts");
        Mat transformed;if(!p.reduced_condition)transformed=mul(condition_rows(p.condition),transpose(pair.basis),&opt->pool);
        Mat rebuilt;if(!p.reduced_condition)rebuilt=Mat(transformed.nrow,transformed.ncol);
        size_t rb=0,cb=0;
        for(size_t s=0;s<ms.size();++s){const auto& m=ms[s];auto local=row_basis(g.evaluate(m.primitive,m.actions));std::vector<Mat> algebra;
            for(size_t e=0;e<local.nrow;++e){Mat map(m.dimension,m.dimension);
                for(size_t i=0;i<m.dimension;++i)map[i]=mul(local[e],m.actions[m.orbit[i]]);
                for(auto h:g.generator_ids)require(equal(mul(map,m.actions[h]),mul(m.actions[h],map)),"condition commutant does not intertwine");
                algebra.push_back(std::move(map));}
            endomorphisms.push_back(std::move(algebra));
            for(size_t c=0;c<p.equations.copies[s];++c){
                if(!p.reduced_condition)for(size_t e=0;e<local.nrow;++e)
                    for(auto [j,x]:transformed[eo[s]+c*m.dimension+local[e](0)])
                        if(size_t(j)>=po[s]&&size_t(j)<po[s]+pair.layout.copies[s]*m.dimension&&(j-po[s])%m.dimension==0)
                            coefficients[rb+c].push_back(index(cb+((j-po[s])/m.dimension)*local.nrow+e),x);
                normalize(coefficients[rb+c]);
                for(auto [j,x]:coefficients[rb+c]){require(size_t(j)>=cb&&size_t(j)<cb+pair.layout.copies[s]*local.nrow,"condition mixes inequivalent irreps");
                    size_t copy=(j-cb)/local.nrow,e=(j-cb)%local.nrow;const auto& ch=channels[s][copy];
                    entries.push_back({s,c,ch.left,ch.left_copy,ch.right,ch.right_copy,ch.channel,e,x});
                    if(!p.reduced_condition)for(size_t i=0;i<m.dimension;++i)for(auto [k,v]:endomorphisms[s][e][i])
                        rebuilt[eo[s]+c*m.dimension+k].push_back(index(po[s]+copy*m.dimension+i),x*v);
                }
            }rb+=p.equations.copies[s];cb+=pair.layout.copies[s]*local.nrow;
        }
        if(!p.reduced_condition){for(auto& row:rebuilt.rows)normalize(row);require(equal(rebuilt,transformed),"condition components are not generated by the reduced symmetry coefficients");}
    }
};

// Generic counterpart of the archive's fixed D3 contraction templates. A map
// from a simple module is determined by the image of its cyclic vector. Its
// allowed images are the primitive corner of older x alphabet. This retains
// ALL End_QG(U) coordinates (also number fields and quaternionic commutants).
// Component contractions below happen once per local Hom basis element; the
// growing recurrence only accumulates coefficients between multiplicity copies.
class RecurrenceTemplates {
    const Group& g;const std::vector<Model>& ms;const Prepared& p;TensorAdapter& adapter;rref_option* opt;
    struct Prototype {Adapted space;ProductCopies channels;std::vector<Mat> corners,left;std::vector<std::map<std::tuple<size_t,size_t,size_t>,size_t>> ids;};
    struct Route {I row,column;Q value;};
    using Routes=std::vector<std::vector<Route>>;
    struct Hom {Mat basis;std::unique_ptr<Chart> chart;std::map<size_t,Routes> routes;};
    std::map<size_t,Prototype> candidate,target;
    std::map<std::pair<size_t,size_t>,Hom> hom;
    std::vector<Mat> local;
    EquivariantConditions conditions;
    std::map<std::pair<size_t,size_t>,std::vector<size_t>> condition_index;
    struct SmallProduct {Adapted space;Mat dual;std::vector<Mat> corner,left;};
    std::map<std::pair<size_t,size_t>,SmallProduct> small_products;
    using SmallKey=std::array<size_t,9>;
    std::map<SmallKey,std::vector<Mat>> small_routes;
    struct Entry{size_t a,old,b,previous,channel;Q value;};
    std::vector<Entry> state;bool have_state=false;
    Layout one(size_t a) const {Layout x;x.copies.assign(ms.size(),0);x.copies[a]=1;return x;}
    Prototype& prototype(size_t a,bool equations) {
        auto& cache=equations?target:candidate;auto it=cache.find(a);if(it!=cache.end())return it->second;
        Prototype out;auto fixed=equations?p.equations:p.alphabet;
        out.space=p.backward?adapter.product(fixed,one(a),&out.channels):adapter.product(one(a),fixed,&out.channels);
        out.ids.resize(ms.size());for(size_t s=0;s<ms.size();++s)for(size_t j=0;j<out.channels[s].size();++j){const auto& c=out.channels[s][j];
            out.ids[s].emplace(p.backward?std::make_tuple(c.left,c.left_copy,c.channel):std::make_tuple(c.right,c.right_copy,c.channel),j);}
        auto offsets=out.space.layout.offsets(ms);Mat inverse;if(equations)inverse=transpose(block_inverse(out.space.basis,opt));
        for(size_t s=0;s<ms.size();++s){Mat q(out.space.layout.copies[s]*local[s].nrow,out.space.basis.nrow);
            for(size_t c=0;c<out.space.layout.copies[s];++c)for(size_t j=0;j<local[s].nrow;++j)
                for(auto [k,v]:local[s][j])q[c*local[s].nrow+j].push_back(index(offsets[s]+c*ms[s].dimension+k),v);
            if(!equations)out.corners.push_back(mul(q,out.space.basis,&opt->pool));
            else {Mat left(q.nrow,inverse.ncol);
                for(size_t c=0;c<out.space.layout.copies[s];++c)for(size_t j=0;j<local[s].nrow;++j)
                    left[c*local[s].nrow+j]=inverse[offsets[s]+c*ms[s].dimension+local[s][j](0)];
                out.left.push_back(std::move(left));}
        }return cache.emplace(a,std::move(out)).first->second;
    }
    Hom& intertwiners(size_t a,size_t b) {
        auto key=std::make_pair(a,b);auto it=hom.find(key);if(it!=hom.end())return it->second;
        Hom h;h.basis=prototype(a,false).corners[b];h.chart=std::make_unique<Chart>(h.basis,opt);
        return hom.emplace(key,std::move(h)).first->second;
    }
    SmallProduct& small_product(size_t a,size_t b) {
        auto key=std::make_pair(a,b);auto it=small_products.find(key);if(it!=small_products.end())return it->second;
        SmallProduct out;out.space=adapter.product(one(a),one(b));out.dual=transpose(block_inverse(out.space.basis,opt));auto offsets=out.space.layout.offsets(ms);
        for(size_t s=0;s<ms.size();++s){Mat q(out.space.layout.copies[s]*local[s].nrow,out.space.basis.nrow),left(q.nrow,q.ncol);
            for(size_t c=0;c<out.space.layout.copies[s];++c)for(size_t j=0;j<local[s].nrow;++j){
                for(auto [k,x]:local[s][j])q[c*local[s].nrow+j].push_back(index(offsets[s]+c*ms[s].dimension+k),x);
                left[c*local[s].nrow+j]=out.dual[offsets[s]+c*ms[s].dimension+local[s][j](0)];}
            out.corner.push_back(mul(q,out.space.basis));out.left.push_back(std::move(left));}
        return small_products.emplace(key,std::move(out)).first->second;
    }
    const std::vector<Mat>& recouple(size_t a,size_t b,size_t u,size_t previous_channel,size_t previous_end,const EquivariantConditions::Entry& c) {
        size_t v=p.backward?c.left:c.right,t=c.equation;SmallKey key{a,b,u,previous_channel,previous_end,v,t,c.channel,c.end};
        auto it=small_routes.find(key);if(it!=small_routes.end())return it->second;
        auto start=std::chrono::steady_clock::now();
        auto& prev=p.backward?small_product(u,a):small_product(a,u);
        const auto& seed=prev.corner[b][previous_channel*local[b].nrow+previous_end];Mat expansion(ms[b].dimension,prev.space.basis.ncol);
        for(size_t i=0;i<ms[b].dimension;++i){auto word=ms[b].orbit[i];auto ambient=p.backward?kron(ms[u].actions[word],ms[a].actions[word]):kron(ms[a].actions[word],ms[u].actions[word]);
            expansion[i]=mul(seed,ambient);}
        auto& pair=small_product(c.left,c.right);auto offsets=pair.space.layout.offsets(ms);Mat chart(ms[t].dimension,pair.space.basis.ncol);
        for(size_t i=0;i<chart.nrow;++i)chart[i]=pair.dual[offsets[t]+c.channel*ms[t].dimension+i];
        auto condition=mul(transpose(conditions.endomorphisms[t][c.end]),chart);
        // Only single irreps occur here: no alphabet-copy or equation-copy
        // axis is expanded. Missing component equations follow from the fixed
        // commutant matrices and the two CG couplings.
        struct DTerm{size_t next,equation;Q value;};std::vector<std::vector<DTerm>> lookup(ms[u].dimension);
        for(size_t q=0;q<condition.nrow;++q)for(auto [j,x]:condition[q]){size_t left=j/ms[c.right].dimension,right=j%ms[c.right].dimension;
            lookup[p.backward?right:left].push_back({p.backward?left:right,q,x});}
        Mat raw(ms[a].dimension*ms[t].dimension,ms[b].dimension*ms[v].dimension);
        for(size_t previous=0;previous<expansion.nrow;++previous)for(auto [j,x]:expansion[previous]){
            size_t old=p.backward?j%ms[a].dimension:j/ms[u].dimension,known=p.backward?j/ms[a].dimension:j%ms[u].dimension;
            for(const auto& d:lookup[known])raw[p.backward?d.equation*ms[a].dimension+old:old*ms[t].dimension+d.equation].push_back(
                index(p.backward?d.next*ms[b].dimension+previous:previous*ms[v].dimension+d.next),x*d.value);}
        for(auto& row:raw.rows)normalize(row);
        auto& source=p.backward?small_product(v,b):small_product(b,v);auto& target=p.backward?small_product(t,a):small_product(a,t);
        std::vector<Mat> out;for(size_t s=0;s<ms.size();++s)out.push_back(mul(target.left[s],mul(raw,transpose(source.corner[s]))));
        recoupling_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return small_routes.emplace(key,std::move(out)).first->second;
    }
    const Routes& routes(size_t a,size_t b,size_t channel) {
        auto& h=intertwiners(a,b);auto it=h.routes.find(channel);if(it!=h.routes.end())return it->second;
        const auto ch=prototype(a,false).channels[b][channel/local[b].nrow];size_t u=p.backward?ch.left:ch.right,copy=p.backward?ch.left_copy:ch.right_copy;
        auto& codomain=prototype(a,true);auto& domain=prototype(b,false);std::vector<Mat> full;
        for(size_t s=0;s<ms.size();++s)full.emplace_back(codomain.left[s].nrow,domain.corners[s].nrow);
        for(auto id:condition_index[{u,copy}]){const auto& c=conditions.entries[id];auto& reduced=recouple(a,b,u,ch.channel,channel%local[b].nrow,c);
            size_t v=p.backward?c.left:c.right,next_copy=p.backward?c.left_copy:c.right_copy;
            for(size_t s=0;s<ms.size();++s){size_t d=local[s].nrow;const auto& r=reduced[s];std::vector<size_t> rm(r.nrow/d),cm(r.ncol/d);
                for(size_t i=0;i<rm.size();++i)rm[i]=codomain.ids[s].at({c.equation,c.copy,i});
                for(size_t i=0;i<cm.size();++i)cm[i]=domain.ids[s].at({v,next_copy,i});
                for(size_t i=0;i<r.nrow;++i)for(auto [j,x]:r[i])full[s][rm[i/d]*d+i%d].push_back(index(cm[j/d]*d+j%d),c.value*x);}}
        Routes result(ms.size());for(size_t s=0;s<ms.size();++s)for(size_t i=0;i<full[s].nrow;++i){normalize(full[s][i]);
            for(auto [j,x]:full[s][i])result[s].push_back({index(i),j,x});}
        ++templates_built;return h.routes.emplace(channel,std::move(result)).first->second;
    }
    // Actual tensor-product order may interleave old copies and fixed-factor
    // copies. Resolve it explicitly, never assume component-major D3 ordering.
    using CopyKey=std::tuple<size_t,size_t,size_t,size_t,size_t>;
    CopyKey copy_key(const ProductCopy& c)const {
        return p.backward?CopyKey{c.right,c.right_copy,c.left,c.left_copy,c.channel}:
                          CopyKey{c.left,c.left_copy,c.right,c.right_copy,c.channel};
    }
    using Maps=std::vector<std::vector<std::vector<std::vector<size_t>>>>;
    Maps maps(const Layout& variable,const ProductCopies& global,bool equations) {
        Maps out(ms.size());for(size_t s=0;s<ms.size();++s){out[s].resize(ms.size());std::map<CopyKey,size_t> ids;
            for(size_t j=0;j<global[s].size();++j)ids.emplace(copy_key(global[s][j]),j);
            for(size_t a=0;a<ms.size();++a)if(variable.copies[a]){auto& channels=prototype(a,equations).channels[s];out[s][a].resize(variable.copies[a]);
                for(size_t c=0;c<variable.copies[a];++c)for(auto ch:channels){auto k=copy_key(ch);std::get<1>(k)=c;out[s][a][c].push_back(ids.at(k));}}
        }return out;
    }
public:
    size_t templates_built=0;double template_s=0,compression_s=0,recoupling_s=0;
    RecurrenceTemplates(const Group& group,const std::vector<Model>& models,const Prepared& prepared,TensorAdapter& cg,rref_option_t options)
        :g(group),ms(models),p(prepared),adapter(cg),opt(options),conditions(prepared,group,models,cg,options){
        for(const auto& m:ms)local.push_back(row_basis(g.evaluate(m.primitive,m.actions)));
        for(size_t i=0;i<conditions.entries.size();++i){const auto& c=conditions.entries[i];condition_index[p.backward?std::make_pair(c.right,c.right_copy):std::make_pair(c.left,c.left_copy)].push_back(i);}
        std::cout<<"condition_component_nnz="<<p.condition.nnz()<<" condition_independent_coefficients="<<conditions.coefficients.nnz()<<std::endl;}
    Mat remember(const Adapted& kernel,const Layout& domain,const ProductCopies& channels) {
        state.clear();auto offsets=domain.offsets(ms),ko=kernel.layout.offsets(ms);size_t nr=0,nc=0;
        for(size_t s=0;s<ms.size();++s){nr+=kernel.layout.copies[s];nc+=domain.copies[s]*local[s].nrow;}
        Mat compact(nr,nc);size_t rowbase=0,colbase=0;
        for(size_t b=0;b<ms.size();++b){std::vector<I> pivots(ms[b].dimension,-1);for(size_t j=0;j<local[b].nrow;++j)pivots[local[b][j](0)]=index(j);
            std::map<size_t,std::map<CopyKey,size_t>> ids;
            for(const auto& ch:channels[b]){auto key=copy_key(ch);size_t a=std::get<0>(key);if(ids.count(a))continue;
                const auto& proto=prototype(a,false).channels[b];for(size_t i=0;i<proto.size();++i)ids[a].emplace(copy_key(proto[i]),i);}
            for(size_t c=0;c<kernel.layout.copies[b];++c)for(auto [raw,value]:kernel.basis[ko[b]+c*ms[b].dimension]){
                require(size_t(raw)>=offsets[b]&&size_t(raw)<offsets[b]+domain.copies[b]*ms[b].dimension,"mixed irrep in recurrence kernel");
                size_t column=(raw-offsets[b])/ms[b].dimension,component=(raw-offsets[b])%ms[b].dimension;I j=pivots[component];if(j<0)continue;
                auto key=copy_key(channels[b][column]);auto a=std::get<0>(key),old=std::get<1>(key);std::get<1>(key)=0;
                state.push_back({a,old,b,c,ids.at(a).at(key)*local[b].nrow+size_t(j),value});
                compact[rowbase+c].push_back(index(colbase+column*local[b].nrow+j),value);
            }rowbase+=kernel.layout.copies[b];colbase+=domain.copies[b]*local[b].nrow;
        }have_state=true;return compact;
    }
    void remember_compact(const Mat& compact,const Layout& solution,const Layout& domain,const ProductCopies& channels) {
        state.clear();size_t rb=0,cb=0;
        require(compact.nrow==std::accumulate(solution.copies.begin(),solution.copies.end(),size_t(0)),"compact recurrence row count mismatch");
        for(size_t b=0;b<ms.size();++b){
            require(local[b].nrow==1,"scalar multiplicities require a split rational model");
            std::vector<std::tuple<size_t,size_t,size_t>> decode;
            std::map<size_t,std::map<CopyKey,size_t>> ids;
            for(const auto& ch:channels[b]){auto key=copy_key(ch);auto a=std::get<0>(key),old=std::get<1>(key);std::get<1>(key)=0;
                if(!ids.count(a)){const auto& proto=prototype(a,false).channels[b];
                    for(size_t i=0;i<proto.size();++i)ids[a].emplace(copy_key(proto[i]),i);}
                decode.emplace_back(a,old,ids.at(a).at(key));}
            for(size_t c=0;c<solution.copies[b];++c)for(auto [j,x]:compact[rb+c]){
                require(size_t(j)>=cb&&size_t(j)<cb+domain.copies[b],"compact recurrence mixes irreps");
                auto [a,old,channel]=decode.at(j-cb);state.push_back({a,old,b,c,channel,x});}
            rb+=solution.copies[b];cb+=domain.copies[b];
        }
        require(compact.ncol==cb,"compact recurrence column count mismatch");have_state=true;
    }
    std::vector<Mat> assemble(const Tensor& seed,const Layout& older,const Layout& previous,
        const ProductCopies& source,const ProductCopies& target_channels,
        const std::function<void(size_t,Mat&&,double)>& consume={}) {
        auto tick=std::chrono::steady_clock::now();auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();};
        std::vector<Entry> seed_entries;const auto& entries=have_state?state:seed_entries;
        const size_t letters=p.alphabet.dimension(ms),old_dimension=older.dimension(ms);auto oo=older.offsets(ms),po=previous.offsets(ms);
        std::vector<std::tuple<size_t,size_t,size_t>> decode(old_dimension);
        for(size_t a=0;a<ms.size();++a)for(size_t c=0;c<older.copies[a];++c)for(size_t j=0;j<ms[a].dimension;++j)decode[oo[a]+c*ms[a].dimension+j]={a,c,j};
        if(!have_state){auto ex=flatten(seed);
        for(size_t b=0;b<ms.size();++b)for(size_t c=0;c<previous.copies[b];++c){std::map<std::pair<size_t,size_t>,Vec> blocks;
            for(auto [raw,value]:ex[po[b]+c*ms[b].dimension]){size_t old=p.backward?raw%old_dimension:raw/letters,letter=p.backward?raw/old_dimension:raw%letters;
                auto [a,copy,component]=decode[old];size_t column=p.backward?letter*ms[a].dimension+component:component*letters+letter;
                blocks[{a,copy}].push_back(index(column),value);}
            for(auto& [key,v]:blocks){auto [a,copy]=key;normalize(v);auto& h=intertwiners(a,b);Mat block(1,h.basis.ncol);block[0]=std::move(v);
                auto coordinates=h.chart->coordinates(block);for(auto [j,value]:coordinates[0])seed_entries.push_back({a,copy,b,c,size_t(j),value});}
        }}
        auto rowmaps=maps(older,target_channels,true),colmaps=maps(previous,source,false);compression_s=elapsed();
        tick=std::chrono::steady_clock::now();
        // Prebuild sequentially: route construction itself uses the worker pool.
        recoupling_s=0;for(const auto& e:entries)routes(e.a,e.b,e.channel);template_s=elapsed();
        std::cout<<"  condition_templates="<<templates_built<<" irrep_recouplings="<<small_routes.size()<<" recoupling_s="<<recoupling_s<<std::endl;
        std::vector<Mat> out(consume?0:ms.size());std::vector<size_t> order(ms.size());std::iota(order.begin(),order.end(),0);
        // Stream a completed sector to the solver before assembling another.
        // Independent sectors share this pool; the consumer uses its own local
        // elimination workspace and must not submit work to this same pool.
        if(consume)std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return
            target_channels[a].size()*source[a].size()>target_channels[b].size()*source[b].size();});
        std::vector<std::future<void>> jobs;
        for(auto s:order)jobs.push_back(opt->pool.submit_task([&,s]{auto start=std::chrono::steady_clock::now();size_t d=local[s].nrow;
            Mat matrix(target_channels[s].size()*d,source[s].size()*d);
            for(const auto& e:entries){const auto& r=hom.at({e.a,e.b}).routes.at(e.channel)[s];auto& rm=rowmaps[s][e.a][e.old];auto& cm=colmaps[s][e.b][e.previous];
                for(const auto& x:r)matrix[rm[x.row/d]*d+x.row%d].push_back(index(cm[x.column/d]*d+x.column%d),e.value*x.value);}
            for(auto& row:matrix.rows)normalize(row);
            if(consume)consume(s,std::move(matrix),std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
            else out[s]=std::move(matrix);
        }));
        // Drain every task before unwinding references captured by the workers.
        std::exception_ptr failure;for(auto& job:jobs)try{job.get();}catch(...){if(!failure)failure=std::current_exception();}
        if(failure)std::rethrow_exception(failure);return out;
    }
};

inline Prepared prepare(const Group& g,const std::vector<Model>& ms,const Tensor& original_condition,
    const Tensor& original_seed,bool backward,rref_option_t opt,const fs::path& root,
    const std::vector<Mat>& terminal_generators={}) {
    require(original_seed.rank()==3&&original_condition.rank()==3,"seed and condition must be rank-three tensors");
    const size_t letters=g.generators[0].nrow,terminal=original_seed.dim(backward?2:1);
    require(original_condition.dim(0)==letters&&original_condition.dim(1)==letters&&original_seed.dim(backward?1:2)==letters,"input alphabet dimensions disagree");
    auto alphabet=adapt(g,ms,g.generators,opt);verify_adaptation(alphabet,g,ms,g.generators);
    auto ag=generators(g,ms,alphabet.layout);auto ai=sparse_mat_inverse(alphabet.basis,QQ,opt);
    auto raw=condition_rows(original_condition);auto raw_basis=row_basis(raw);
    // Derive equation actions while their original row basis is still sparse
    // and in RREF. Changing the domain coordinates leaves these actions equal,
    // but deriving them afterwards requires a much denser coordinate inverse.
    std::vector<Mat> original_pair;for(const auto& r:g.generators)original_pair.push_back(kron(r,r));
    auto eqg=equation_actions(raw_basis,original_pair,opt);auto equations=adapt(g,ms,eqg,opt);verify_adaptation(equations,g,ms,eqg);
    auto eq_inverse=sparse_mat_inverse(equations.basis,QQ,opt);
    auto a=mul(raw_basis,transpose(kron(alphabet.basis,alphabet.basis)),&opt->pool);
    auto da=mul(transpose(eq_inverse),a,&opt->pool);
    std::vector<Mat> pair;for(const auto& r:ag)pair.push_back(kron(r,r));
    verify_equivariance(da,pair,generators(g,ms,equations.layout));
    std::vector<Mat> tg=terminal_generators;
    if(tg.empty()){require(terminal==1,"non-scalar terminal requires --terminal-generator for every group generator");for(size_t j=0;j<g.names.size();++j)tg.push_back(identity(1));}
    require(tg.size()==g.names.size()&&tg[0].nrow==terminal,"terminal representation mismatch");
    auto term=adapt(g,ms,tg,opt);verify_adaptation(term,g,ms,tg);
    auto ti=sparse_mat_inverse(term.basis,QQ,opt);auto newtg=generators(g,ms,term.layout);
    Mat seed=mul(flatten(original_seed),backward?kron(ai,ti):kron(ti,ai),&opt->pool);
    std::vector<Mat> ambient;for(size_t j=0;j<ag.size();++j)ambient.push_back(backward?kron(ag[j],newtg[j]):kron(newtg[j],ag[j]));
    auto seedg=induced(seed,ambient,opt);auto sa=adapt(g,ms,seedg,opt);verify_adaptation(sa,g,ms,seedg);
    auto newseed=mul(sa.basis,seed,&opt->pool);
    Prepared out{alphabet.layout,equations.layout,sa.layout,term.layout,condition_tensor(da,letters),
        tensor(newseed,backward?std::vector<size_t>{newseed.nrow,letters,terminal}:std::vector<size_t>{newseed.nrow,terminal,letters}),backward};
    TensorAdapter adapter(g,ms,opt);EquivariantConditions reduced(out,g,ms,adapter,opt);
    out.reduced_condition=reduced.coefficients;out.condition_pair_basis=reduced.pair.basis;
    save_group(root,g,ms);save_space(root,"alphabet",alphabet,g,ms);save_space(root,"equations",equations,g,ms);
    save_space(root,"seed",sa,g,ms);save_space(root,"terminal",term,g,ms);
    write(root/"alphabet_inverse.wxf",ai);write(root/"condition_independent_rows.wxf",raw_basis);
    write(root/"dlogmat.wxf",out.condition);write(root/"seed.wxf",out.expansion);
    write(root/"dlogmat_multiplicity.wxf",reduced.coefficients);write(root/"dlogmat_pair_basis.wxf",reduced.pair.basis);
    save_layout(root/"dlogmat_pair_copies.tsv",reduced.pair.layout,ms);
    std::ofstream reduced_meta(root/"dlogmat_reduced.tsv");reduced_meta<<"symbology-equivariant-condition-v1\n";
    reduced_meta.close();require(bool(reduced_meta),"cannot write reduced condition metadata");
    std::cout<<"condition_component_nnz="<<out.condition.nnz()<<" condition_independent_coefficients="<<reduced.coefficients.nnz()<<std::endl;
    std::ofstream meta(root/"bootstrap.tsv");meta<<"symbology-symrep-bootstrap-v1 "<<(backward?"backward":"forward")<<'\n';meta.close();require(bool(meta),"cannot write bootstrap metadata");
    seal(root);return out;
}
inline Prepared load_prepared(const fs::path& root,const std::vector<Model>& ms,bool component_data=true) {
    std::ifstream in(root/"bootstrap.tsv");std::string schema,direction;in>>schema>>direction;
    require(bool(in)&&schema=="symbology-symrep-bootstrap-v1"&&(direction=="forward"||direction=="backward"),"invalid prepared bootstrap schema");
    Prepared p;p.backward=direction=="backward";
    p.alphabet=load_layout(root/"alphabet_copies.tsv",ms);p.equations=load_layout(root/"equations_copies.tsv",ms);
    p.seed=load_layout(root/"seed_copies.tsv",ms);p.terminal=load_layout(root/"terminal_copies.tsv",ms);
    if(component_data)p.condition=read(root/"dlogmat.wxf");p.expansion=read(root/"seed.wxf");
    if(component_data&&fs::exists(root/"dlogmat_reduced.tsv")){
        std::ifstream meta(root/"dlogmat_reduced.tsv");std::string format;meta>>format;
        require(bool(meta)&&format=="symbology-equivariant-condition-v1","invalid reduced condition schema");
        p.reduced_condition=read_matrix(root/"dlogmat_multiplicity.wxf");p.condition_pair_basis=read_matrix(root/"dlogmat_pair_basis.wxf");
    }
    require(!component_data||(p.condition.rank()==3&&p.condition.dim(0)==p.alphabet.dimension(ms)&&p.condition.dim(1)==p.alphabet.dimension(ms)&&p.condition.dim(2)==p.equations.dimension(ms)),"prepared condition/layout mismatch");
    require(p.expansion.rank()==3&&p.expansion.dim(0)==p.seed.dimension(ms)&&p.expansion.dim(p.backward?1:2)==p.alphabet.dimension(ms)&&p.expansion.dim(p.backward?2:1)==p.terminal.dimension(ms),"prepared seed/layout mismatch");return p;
}
inline double seconds(std::chrono::steady_clock::time_point start){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
// v1 stores expanded irreducible recurrences; v2 keeps a sparse carrier
// recurrence and its exact irreducible frame as separate matrices.
inline size_t chain_weight(const fs::path& chain,const fs::path& prepared,bool* factorized=nullptr,bool* multiplicity=nullptr) {
    verify_seal(chain);std::ifstream in(chain/"chain.tsv");std::string schema,key,hash,format,extra;size_t weight=0;
    in>>schema>>key>>hash>>weight;
    bool factored=schema=="symbology-symrep-chain-v2"||schema=="symbology-symrep-chain-v3";
    bool compact=schema=="symbology-symrep-chain-v4";
    require(bool(in)&&(schema=="symbology-symrep-chain-v1"||factored||compact)&&key=="prepared_sha256"&&weight>0,"invalid recurrence-chain metadata");
    if(factored){in>>format;require(bool(in)&&format==(schema=="symbology-symrep-chain-v3"?"factorized-orbits":"factorized"),"invalid factorized chain format");}
    if(compact){in>>format;require(bool(in)&&format=="multiplicity","invalid compact multiplicity chain format");}
    require(!(in>>extra),"trailing recurrence-chain metadata");
    require(hash==native_cache::digest(prepared/"checksums.tsv"),"recurrence chain belongs to a different prepared bundle");
    if(factorized)*factorized=factored;if(multiplicity)*multiplicity=compact;return weight;
}
#include "symrep_carrier.hpp"
#include "symrep_multiplicity.hpp"
inline void extend(const Prepared& p,const Group& g,const std::vector<Model>& ms,size_t max_weight,
    const fs::path& root,rref_option_t opt,bool compare,const fs::path& prepared,const fs::path& resume={}) {
    auto seed=p.expansion;auto previous=p.seed,older=p.terminal;TensorAdapter adapter(g,ms,opt);
    RecurrenceTemplates recurrence(g,ms,p,adapter,opt);
    const size_t letters=p.alphabet.dimension(ms);auto letter_g=generators(g,ms,p.alphabet),eq_g=generators(g,ms,p.equations);
    size_t first=2;std::ofstream stats(root/"timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tadapt_s\tkernel_s\tordinary_kernel_s\tverification_s\tcorner_columns\ttemplate_s\tcompression_s\trref_s\tshorten_s\tlift_s\n";
    if(resume.empty()) {
        write(root/"w1.wxf",seed);save_layout(root/"w1_copies.tsv",previous,ms);
    }else {
        bool factorized=false;auto last=chain_weight(resume,prepared,&factorized);require(!factorized,"expand the factorized chain before continuing with the multiplicity backend");require(max_weight>last,"--max-weight must exceed the saved chain weight");first=last+1;
        require(equal(flatten(read(resume/"w1.wxf")),flatten(p.expansion)),"saved seed differs from prepared seed");
        for(const auto& entry:fs::directory_iterator(resume))if(entry.is_regular_file()&&entry.path().filename().string().starts_with("w"))fs::copy_file(entry.path(),root/entry.path().filename());
        std::ifstream oldstats(resume/"timings.tsv");require(bool(oldstats),"missing saved timings");std::string line;std::getline(oldstats,line);
        while(std::getline(oldstats,line)){if(line.empty())continue;auto fields=1+std::count(line.begin(),line.end(),'\t');stats<<line;for(size_t j=fields;j<16;++j)stats<<"\tnan";stats<<'\n';}
        seed=read(resume/("w"+std::to_string(last)+".wxf"));previous=load_layout(resume/("w"+std::to_string(last)+"_copies.tsv"),ms);
        if(last>1)older=load_layout(resume/("w"+std::to_string(last-1)+"_copies.tsv"),ms);
        require(seed.rank()==3&&seed.dim(0)==previous.dimension(ms)&&seed.dim(p.backward?1:2)==letters&&seed.dim(p.backward?2:1)==older.dimension(ms),"saved recurrence dimensions disagree with layouts");
    }
    for(size_t w=first;w<=max_weight;++w) {
        auto start=std::chrono::steady_clock::now();
        ProductCopies source_channels,target_channels;
        auto candidates=p.backward?adapter.product(p.alphabet,previous,&source_channels):adapter.product(previous,p.alphabet,&source_channels);
        auto codomain=p.backward?adapter.product(p.equations,older,&target_channels):adapter.product(older,p.equations,&target_channels);
        double adapt_s=seconds(start),assembly_s=0;size_t observed_nnz=0;
        std::cout<<"weight="<<w<<" phase=reduced_assembly"<<std::endl;
        start=std::chrono::steady_clock::now();
        auto reduced_blocks=recurrence.assemble(seed,older,previous,source_channels,target_channels);assembly_s=seconds(start);
        Mat target_chart;std::vector<size_t> target_offsets;
        if(compare){target_chart=transpose(block_inverse(codomain.basis,opt));target_offsets=codomain.layout.offsets(ms);}
        KernelTimings kernel_times;
        auto solution=kernel_blocks(candidates.basis.nrow,g,ms,candidates.layout,[&](size_t sector,const Mat& q){
            auto reduced=std::move(reduced_blocks[sector]);
            if(compare){const auto& m=ms[sector];auto local=row_basis(g.evaluate(m.primitive,m.actions));
                Mat left(codomain.layout.copies[sector]*local.nrow,target_chart.ncol);
                for(size_t c=0;c<codomain.layout.copies[sector];++c)for(size_t j=0;j<local.nrow;++j)
                    left[c*local.nrow+j]=target_chart[target_offsets[sector]+c*m.dimension+local[j](0)];
                auto direct=assemble_reduced(seed,p.condition,p.backward,left,mul(q,candidates.basis,&opt->pool),opt);
                require(equal(reduced,direct),"cached multiplicity contraction differs from component contraction");}
            observed_nnz+=reduced.nnz();
            std::cout<<"  irrep="<<sector<<" rows="<<reduced.nrow<<" cols="<<reduced.ncol<<" nnz="<<reduced.nnz()<<std::endl;
            return reduced;
        },opt,&kernel_times);
        auto compact=recurrence.remember(solution.solution,candidates.layout,source_channels);
        auto coefficients=mul(solution.solution.basis,candidates.basis,&opt->pool);double kernel_s=seconds(start)-assembly_s,ordinary_s=0;
        start=std::chrono::steady_clock::now();
        auto prev_g=generators(g,ms,previous),old_g=generators(g,ms,older);std::vector<Mat> source,target;
        for(size_t j=0;j<g.names.size();++j){source.push_back(p.backward?kron(letter_g[j],prev_g[j]):kron(prev_g[j],letter_g[j]));target.push_back(p.backward?kron(eq_g[j],old_g[j]):kron(old_g[j],eq_g[j]));}
        verify_adaptation({coefficients,solution.solution.layout},g,ms,source);double verification_s=seconds(start);
        if(compare) {
            std::cout<<"weight="<<w<<" phase=ordinary_comparison"<<std::endl;
            start=std::chrono::steady_clock::now();auto a=assemble(seed,p.condition,p.backward);
            verify_equivariance(a,source,target);verify_adaptation(candidates,g,ms,source);
            require(!mul(a,transpose(coefficients),&opt->pool).nnz(),"extension residual is nonzero");verification_s+=seconds(start);
            start=std::chrono::steady_clock::now();auto baseline=kernel(a,opt);ordinary_s=seconds(start);
            start=std::chrono::steady_clock::now();require(baseline.nrow==coefficients.nrow,"ordinary and symmetry kernel dimensions differ");
            Chart reference(baseline,opt);reference.coordinates(coefficients);verification_s+=seconds(start);
        }
        auto next=tensor(coefficients,p.backward?std::vector<size_t>{coefficients.nrow,letters,previous.dimension(ms)}:std::vector<size_t>{coefficients.nrow,previous.dimension(ms),letters});
        write(root/("w"+std::to_string(w)+".wxf"),next);save_layout(root/("w"+std::to_string(w)+"_copies.tsv"),solution.solution.layout,ms);
        write(root/("w"+std::to_string(w)+"_candidate_basis.wxf"),candidates.basis);
        write(root/("w"+std::to_string(w)+"_kernel.wxf"),solution.solution.basis);
        write(root/("w"+std::to_string(w)+"_multiplicity.wxf"),compact);
        stats<<std::setprecision(9)<<w<<'\t'<<codomain.basis.nrow<<'\t'<<candidates.basis.nrow<<'\t'<<observed_nnz<<'\t'<<coefficients.nrow<<'\t'<<assembly_s<<'\t'<<adapt_s<<'\t'<<kernel_s<<'\t'<<ordinary_s<<'\t'<<verification_s<<'\t';
        for(size_t i=0;i<ms.size();++i){if(i)stats<<',';stats<<solution.corner_columns[i];}stats<<'\t'<<recurrence.template_s<<'\t'<<recurrence.compression_s<<'\t'<<kernel_times.rref<<'\t'<<kernel_times.shorten<<'\t'<<kernel_times.lift<<std::endl;
        std::cout<<"weight="<<w<<" dimension="<<coefficients.nrow<<" symmetry_s="<<assembly_s+adapt_s+kernel_s<<" ordinary_kernel_s="<<ordinary_s<<std::endl;
        older=previous;previous=solution.solution.layout;seed=std::move(next);
    }
    stats.close();require(bool(stats),"cannot write timing report");
    std::ofstream meta(root/"chain.tsv");meta<<"symbology-symrep-chain-v1 prepared_sha256 "<<native_cache::digest(prepared/"checksums.tsv")<<' '<<max_weight<<'\n';meta.close();require(bool(meta),"cannot write chain metadata");seal(root);
}
// Baseline uses exactly the same SparseRREF options and rational arithmetic,
// but retains the original alphabet and the ordinary kernel at every weight.
inline void ordinary_chain(Tensor seed,const Tensor& condition,bool backward,size_t maximum,
    const fs::path& root,rref_option_t opt) {
    require(seed.rank()==3&&condition.rank()==3,"ordinary recurrence needs rank-three tensors");
    write(root/"w1.wxf",seed);std::ofstream stats(root/"timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tkernel_s\tverification_s\n";
    for(size_t w=2;w<=maximum;++w){
        auto start=std::chrono::steady_clock::now();auto a=assemble(seed,condition,backward);double assembly_s=seconds(start);
        start=std::chrono::steady_clock::now();auto k=kernel(a,opt);double kernel_s=seconds(start);
        start=std::chrono::steady_clock::now();require(!mul(a,transpose(k),&opt->pool).nnz(),"ordinary kernel residual");double verification_s=seconds(start);
        auto next=tensor(k,backward?std::vector<size_t>{k.nrow,condition.dim(0),seed.dim(0)}:std::vector<size_t>{k.nrow,seed.dim(0),condition.dim(0)});
        write(root/("w"+std::to_string(w)+".wxf"),next);seed=std::move(next);
        stats<<std::setprecision(9)<<w<<'\t'<<a.nrow<<'\t'<<a.ncol<<'\t'<<a.nnz()<<'\t'<<k.nrow<<'\t'<<assembly_s<<'\t'<<kernel_s<<'\t'<<verification_s<<std::endl;
        std::cout<<"ordinary weight="<<w<<" dimension="<<k.nrow<<" solve_s="<<assembly_s+kernel_s<<std::endl;
    }stats.close();require(bool(stats),"cannot write ordinary timings");seal(root);
}
inline void verify_chain(const fs::path& prepared,const fs::path& chain,const fs::path& reference,
    size_t maximum,const Group& g,const std::vector<Model>& ms,rref_option_t opt,const fs::path& output) {
    bool compact=false;require(maximum<=chain_weight(chain,prepared,nullptr,&compact),"requested verification exceeds saved chain");auto p=load_prepared(prepared,ms);auto letter=read_matrix(prepared/"alphabet_basis.wxf");
    TensorAdapter adapter(g,ms,opt);if(compact)certify_scalar_corners(g,ms);
    auto previous=read_matrix(prepared/"terminal_basis.wxf");
    for(size_t w=1;w<=maximum;++w) {
        auto stem="w"+std::to_string(w);auto actual=flatten(compact?materialize_multiplicity(chain,w,p,ms,adapter,opt):read(chain/(stem+".wxf"))),expected=flatten(read(reference/(stem+".wxf")));
        auto ambient=p.backward?kron(letter,previous):kron(previous,letter);
        auto physical=mul(actual,ambient,&opt->pool);
        require(physical.nrow==expected.nrow,"reference dimension mismatch at "+stem);
        Chart chart(expected,opt);auto change=chart.coordinates(physical);
        require(kernel(change,opt).nrow==0,"singular reference basis change at "+stem);
        write(output/(stem+"_basis_change.wxf"),change);previous=std::move(change);
        std::cout<<"EXACT_REFERENCE_PASS weight="<<w<<" dimension="<<expected.nrow<<std::endl;
    }seal(output);
}

// A symmetry-adapted chain already specifies every generator by its irrep
// layout. Exporting these matrices never needs a carrier frame or LEC expansion.
inline void export_chain_actions(const fs::path& prepared,const fs::path& chain,size_t maximum,
    const Group& g,const std::vector<Model>& ms,rref_option_t opt,const fs::path& output,bool carrier_coordinates=false){
    bool factorized=false;auto last=chain_weight(chain,prepared,&factorized);
    require(maximum<=last,"requested actions exceed saved chain");
    require(!carrier_coordinates||factorized,"carrier coordinates require a factorized chain");
    bool backward=false;std::vector<Mat> previous;
    if(carrier_coordinates){std::ifstream direction(chain/"carrier_direction.tsv");std::string name,extra;direction>>name;
        require(bool(direction)&&(name=="forward"||name=="backward")&&!(direction>>extra),"invalid carrier direction metadata");backward=name=="backward";}
    std::ofstream metadata(output/"actions.tsv"),stats(output/"timings.tsv");
    metadata<<"symbology-symrep-actions-v1\ncoordinates\t"<<(carrier_coordinates?"carrier":"adapted")
        <<"\nprepared_sha256\t"<<native_cache::digest(prepared/"checksums.tsv")
        <<"\nchain_sha256\t"<<native_cache::digest(chain/"checksums.tsv")<<"\nmaximum_weight\t"<<maximum<<'\n';
    for(size_t j=0;j<g.names.size();++j)metadata<<"generator\t"<<j<<'\t'<<g.names[j]<<'\n';
    stats<<"weight\tdimension\tconstruction_s\taction_nnz\n";
    for(size_t w=1;w<=maximum;++w){auto stem=weight_name(w);auto layout=load_layout(chain/(stem+"_copies.tsv"),ms);Mat basis;
        if(carrier_coordinates&&w>1)basis=flatten(read(chain/(stem+"_carrier.wxf")));
        auto start=std::chrono::steady_clock::now();std::vector<Mat> matrices;
        if(!carrier_coordinates)matrices=generators(g,ms,layout);
        else if(w==1)matrices=load_carrier_actions(chain,1,g,layout.dimension(ms),opt);
        else matrices=carrier_actions(basis,previous,g.generators,backward,opt);
        auto elapsed=seconds(start);size_t nnz=0;
        for(size_t j=0;j<matrices.size();++j){require(matrices[j].nrow==layout.dimension(ms)&&matrices[j].ncol==layout.dimension(ms),"action/layout dimension mismatch");
            nnz+=matrices[j].nnz();write(output/(stem+"_g"+std::to_string(j)+".wxf"),matrices[j]);}
        save_layout(output/(stem+"_copies.tsv"),layout,ms);
        stats<<std::setprecision(9)<<w<<'\t'<<layout.dimension(ms)<<'\t'<<elapsed<<'\t'<<nnz<<'\n';
        if(carrier_coordinates)previous=std::move(matrices);
    }
    metadata.close();stats.close();require(bool(metadata)&&bool(stats),"cannot write symmetry actions metadata");seal(output);
}

} // namespace symrep
