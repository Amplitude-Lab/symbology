#pragma once
// An irreducible basis represented as frame * sparse carrier. Recursion uses
// the carrier chart and keeps the exact group action and frame alongside it.
using ModVec=sparse_vec<ulong,I>;
using ModMat=sparse_mat<ulong,I>;
struct BadCarrierPrime:std::runtime_error {using std::runtime_error::runtime_error;};
inline void require_carrier_prime(bool value,const char* reason){if(!value)throw BadCarrierPrime(reason);}
struct FiniteSpan {
    field_t field;std::map<I,ModVec> pivots;
    explicit FiniteSpan(ulong prime):field(FIELD_Fp,prime){}
    ModVec reduce(ModVec v)const{while(v.nnz()){auto it=pivots.find(v(0));if(it==pivots.end())break;
            ulong x=v[0];sparse_vec_sub_mul(v,it->second,x,field);}return v;}
    bool add(ModVec v){v=reduce(std::move(v));if(!v.nnz())return false;
        sparse_vec_rescale(v,nmod_inv(v[0],field.mod),field);I lead=v(0);pivots.emplace(lead,std::move(v));return true;}
};
inline ModVec modular_product(const ModVec& v,const ModMat& m,const field_t& fp){
    std::vector<ulong> accum(m.ncol);for(auto [j,x]:v){ulong pre=n_mulmod_precomp_shoup(x,fp.mod.n);for(auto [k,y]:m[j])accum[k]=_nmod_add(accum[k],n_mulmod_shoup(x,y,pre,fp.mod.n),fp.mod);}
    ModVec result;for(size_t j=0;j<accum.size();++j)if(accum[j])result.push_back(index(j),accum[j]);return result;
}
inline ModMat modular_matrix(const Mat& a,const field_t& fp){ModMat out(a.nrow,a.ncol);
    for(size_t i=0;i<a.nrow;++i){out[i]=a[i]%fp.mod;out[i].canonicalize();}return out;}

// A dense accumulator avoids an ordered-map lookup per scalar multiply
// when a sparse group orbit spreads across many carrier coordinates.
inline Vec carrier_product(const Vec& v,const Mat& a){
    if(v.nnz()==1){auto result=a[v(0)];sparse_vec_rescale(result,v[0],QQ);return result;}
    std::vector<Q> accum(a.ncol);
    for(auto [j,x]:v)for(auto [k,y]:a[j])accum[k]+=x*y;
    Vec out;for(size_t j=0;j<accum.size();++j)if(accum[j]!=0)out.push_back(index(j),accum[j]);return out;
}
// Prove the orbit intertwiner in the group algebra itself. Evaluation in any
// representation preserves this identity, independent of carrier dimension.
// models() already establishes these identities while constructing rho; keep
// this inexpensive explicit check to guard the low-level carrier-frame API.
inline void certify_carrier_models(const Group& g,const std::vector<Model>& ms){
    for(const auto& m:ms){Mat orbit(m.dimension,g.order());
        for(size_t i=0;i<m.dimension;++i)orbit[i]=g.translate(m.primitive,m.orbit[i]);
        for(auto h:g.generator_ids){Mat actual(m.dimension,g.order());
            for(size_t i=0;i<m.dimension;++i)actual[i]=g.translate(m.primitive,g.table[m.orbit[i]][h]);
            require(equal(actual,mul(m.actions[h],orbit)),"carrier orbit fails the universal intertwiner identity");}
    }
}
// Kernel rows have private free-variable columns. Reading those columns is
// enough to recover coordinates; no transpose of the full QQ basis is needed.
struct PrivateCarrierChart {
    std::vector<I> columns;
    std::vector<Q> values;
    bool complete=false;
    explicit PrivateCarrierChart(const Mat& carrier,const field_t* fp=nullptr)
        :columns(carrier.nrow,-1),values(carrier.nrow){
        std::vector<I> owner(carrier.ncol,-1);
        for(size_t i=0;i<carrier.nrow;++i)for(auto [j,x]:carrier[i]){
            if(fp)require_carrier_prime(x.den()%fp->mod!=0,"bad carrier denominator prime");
            if(x!=0)owner[j]=owner[j]==-1?index(i):I(-2);
        }
        for(size_t j=0;j<owner.size();++j)if(owner[j]>=0&&columns[owner[j]]<0){
            auto i=owner[j];columns[i]=index(j);values[i]=*carrier[i].find(index(j));}
        complete=std::all_of(columns.begin(),columns.end(),[](I j){return j>=0;});
    }
};
// Construct only (L tensor R)[:, selected]. Traverse incoming entries of
// selected columns, rather than scanning the nonzeros of the whole product.
// Appending columns in order also makes every output row sorted by construction.
template<class T> inline sparse_mat<T,I> tensor_action_columns(const sparse_mat<T,I>& left,
    const sparse_mat<T,I>& right,const std::vector<I>& columns,const std::vector<T>* scales,const field_t& field){
    sparse_mat<T,I> lt(left.ncol,left.nrow),rt(right.ncol,right.nrow);
    for(size_t i=0;i<left.nrow;++i){I row=index(i);for(auto [j,x]:left[i])lt[j].push_back(row,x);}
    for(size_t i=0;i<right.nrow;++i){I row=index(i);for(auto [j,x]:right[i])rt[j].push_back(row,x);}
    sparse_mat<T,I> out(left.nrow*right.nrow,columns.size());
    for(size_t c=0;c<columns.size();++c){I column=index(c);size_t k=size_t(columns[c])/right.ncol,l=size_t(columns[c])%right.ncol;
        for(auto [i,x]:lt[k])for(auto [j,y]:rt[l]){
            T value=scalar_mul(x,y,field);if(scales)value=scalar_mul(value,(*scales)[c],field);
            if(value!=0)out[size_t(i)*right.nrow+j].push_back(column,value);
        }
    }return out;
}
inline std::vector<Mat> carrier_actions(const Mat& carrier,const std::vector<Mat>& older,
    const std::vector<Mat>& letters,bool backward,rref_option_t opt){
    require(older.size()==letters.size(),"carrier generator count mismatch");
    PrivateCarrierChart chart(carrier);std::vector<I> columns=chart.columns;std::vector<Q> scales;Mat inverse;
    if(chart.complete){for(const auto& x:chart.values)scales.push_back(Q(1)/x);}
    else{Chart general(carrier,opt);columns=std::move(general.pivots);inverse=std::move(general.inverse);}
    std::vector<Mat> result;
    for(size_t h=0;h<letters.size();++h){const auto& left=backward?letters[h]:older[h];const auto& right=backward?older[h]:letters[h];
        require(left.nrow==left.ncol&&right.nrow==right.ncol&&left.nrow*right.nrow==carrier.ncol,"carrier action dimension mismatch");
        auto restriction=tensor_action_columns(left,right,columns,chart.complete?&scales:nullptr,QQ);
        auto a=mul(carrier,restriction,&opt->pool);
        result.push_back(chart.complete?std::move(a):mul(a,inverse,&opt->pool));
    }return result;
}
// Normal recursion needs only modular independence certificates. Exact QQ
// transformations are available separately, in either carrier or irrep coordinates.
inline std::vector<ModMat> carrier_actions_mod(const Mat& carrier,const std::vector<ModMat>& older,
    const std::vector<ModMat>& letters,bool backward,const field_t& fp,rref_option_t opt){
    require(older.size()==letters.size(),"modular carrier generator count mismatch");
    size_t n=carrier.nrow,nc=carrier.ncol;PrivateCarrierChart chart(carrier,&fp);
    auto columns=chart.columns;std::vector<ulong> diagonal;ModMat inverse;auto basis=modular_matrix(carrier,fp);
    if(chart.complete){for(const auto& x:chart.values){ulong value=x%fp.mod;require_carrier_prime(value!=0,"bad carrier chart prime");diagonal.push_back(nmod_inv(value,fp.mod));}}
    else{auto echelon=basis;auto pivots=sparse_mat_rref_forward(echelon,fp,opt);columns.clear();
        for(const auto& batch:pivots)for(auto [r,c]:batch)columns.push_back(c);
        require_carrier_prime(columns.size()==n,"bad carrier rank prime");std::sort(columns.begin(),columns.end());
        std::vector<I> selected(nc,-1);for(size_t i=0;i<n;++i)selected[columns[i]]=index(i);
        ModMat minor(n,n);for(size_t i=0;i<n;++i)for(auto [j,x]:basis[i])if(selected[j]>=0)minor[i].push_back(selected[j],x);
        inverse=sparse_mat_inverse(minor,fp,opt);}
    std::vector<ModMat> result;
    for(size_t h=0;h<letters.size();++h){const auto& left=backward?letters[h]:older[h];const auto& right=backward?older[h]:letters[h];
        require(left.nrow==left.ncol&&right.nrow==right.ncol&&left.nrow*right.nrow==nc,"modular carrier action dimension mismatch");
        auto restriction=tensor_action_columns(left,right,columns,chart.complete?&diagonal:nullptr,fp);
        auto a=sparse_mat_mul(basis,restriction,fp,&opt->pool);
        result.push_back(chart.complete?std::move(a):sparse_mat_mul(a,inverse,fp,&opt->pool));
    }return result;
}
struct CarrierFrameTimings {double modular=0,selection=0,lift=0,check=0;};
inline Adapted carrier_frame(const Group& g,const std::vector<Model>& ms,const std::vector<Mat>& generators,rref_option_t opt,CarrierFrameTimings* timing=nullptr,bool verify_full=true){
    require(generators.size()==g.generators.size(),"carrier generator count mismatch");
    size_t n=generators[0].nrow;
    for(const auto& a:generators)require(a.nrow==n&&a.ncol==n,"carrier generators must be square and equally sized");
    Adapted result;result.layout.copies.resize(ms.size());
    if(!n){result.basis=Mat(0,0);return result;}
    // Only modular projector rows and a small corner chart participate in
    // copy selection. A nonzero modular minor certifies QQ independence.
    // Idempotent trace gives its exact rank (an integer in [0,n], n < prime),
    // so covering every corner and n dimensions also certifies completeness.
    // Reconstruct the selected generators over QQ and generate each full
    // irreducible copy by its certified model orbit. Full matrix checks are
    // available independently; the recursive driver uses the universal proof.
    certify_carrier_models(g,ms);
    std::vector<Mat> local_projectors;
    for(const auto& m:ms)local_projectors.push_back(g.evaluate(m.primitive,m.actions));
    auto started=std::chrono::steady_clock::now();
    ulong prime=1ULL<<60;
    for(size_t attempt=0;attempt<16;++attempt){prime=n_nextprime(prime,0);field_t fp(FIELD_Fp,prime);
        bool bad=false;std::vector<ModMat> mod_generators,mod_actions;
        for(const auto& a:generators){for(const auto& row:a.rows)for(auto [j,x]:row)if(x.den()%fp.mod==0)bad=true;if(bad)break;mod_generators.push_back(modular_matrix(a,fp));}
        for(const auto& a:local_projectors)for(const auto& row:a.rows)for(auto [j,x]:row)if(x.den()%fp.mod==0)bad=true;
        if(bad)continue;
        ModMat unit(n,n);for(size_t i=0;i<n;++i)unit[i].push_back(index(i),ulong(1));mod_actions.push_back(std::move(unit));
        for(size_t h=1;h<g.order();++h)mod_actions.push_back(sparse_mat_mul(mod_actions[g.words[h].first],mod_generators[g.words[h].second],fp,&opt->pool));
        if(timing)timing->modular+=seconds(started);started=std::chrono::steady_clock::now();
        struct Seed{size_t sector,row;};std::vector<Seed> seeds;result.layout.copies.assign(ms.size(),0);
        for(size_t s=0;s<ms.size();++s){ModMat pm(n,n);std::vector<ulong> coefficients;
            for(auto& x:ms[s].primitive){if(x.den()%fp.mod==0){bad=true;break;}coefficients.push_back(x%fp.mod);}if(bad)break;
            for(size_t i=0;i<n;++i){std::vector<ulong> values(n);
                for(size_t h=0;h<g.order();++h)if(coefficients[h])for(auto [j,x]:mod_actions[h][i])values[j]=nmod_add(values[j],nmod_mul(coefficients[h],x,fp.mod),fp.mod);
                for(size_t j=0;j<n;++j)if(values[j])pm[i].push_back(index(j),values[j]);}
            auto echelon=pm;
            auto pivots=sparse_mat_rref_forward(echelon,fp,opt);std::vector<I> columns;
            for(auto& batch:pivots)for(auto [r,c]:batch)columns.push_back(c);std::sort(columns.begin(),columns.end());
            ulong trace=0;for(size_t i=0;i<n;++i)if(auto v=pm[i].find(index(i)))trace=nmod_add(trace,*v,fp.mod);
            if(trace!=columns.size()||columns.size()%ms[s].division_dimension){bad=true;break;}
            std::vector<I> selected(n,-1);for(size_t i=0;i<columns.size();++i)selected[columns[i]]=index(i);
            auto restrict=[&](const ModVec& v){ModVec out;for(auto [j,x]:v)if(selected[j]>=0)out.push_back(selected[j],x);return out;};
            std::vector<ModMat> moves;for(auto h:ms[s].orbit){ModMat move(n,columns.size());for(size_t i=0;i<n;++i)move[i]=restrict(mod_actions[h][i]);moves.push_back(std::move(move));}
            auto local=modular_matrix(local_projectors[s],fp);FiniteSpan covered(prime);
            std::vector<size_t> order(n);std::iota(order.begin(),order.end(),0);
            std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return pm[a].nnz()<pm[b].nnz();});
            for(auto i:order){auto& v=pm[i];if(!v.nnz()||!covered.reduce(restrict(v)).nnz())continue;size_t before=covered.pivots.size();
                ModMat orbit(ms[s].dimension,columns.size());for(size_t j=0;j<moves.size();++j)orbit[j]=modular_product(v,moves[j],fp);
                for(const auto& row:local.rows)covered.add(modular_product(row,orbit,fp));
                if(covered.pivots.size()-before!=ms[s].division_dimension){bad=true;break;}
                seeds.push_back({s,i});++result.layout.copies[s];if(covered.pivots.size()==columns.size())break;}
            if(covered.pivots.size()!=columns.size())bad=true;if(bad)break;}
        if(bad||result.layout.dimension(ms)!=n)continue;
        if(timing)timing->selection+=seconds(started);started=std::chrono::steady_clock::now();
        // Batch just the selected unit-vector orbits, so SparseRREF can use
        // the configured worker pool without forming full QQ group matrices.
        std::map<size_t,size_t> selected_rows;
        for(auto seed:seeds)if(!selected_rows.count(seed.row)){size_t pos=selected_rows.size();selected_rows.emplace(seed.row,pos);}
        std::vector<Mat> unit_orbits;Mat units(selected_rows.size(),n);
        for(auto [row,pos]:selected_rows)units[pos].push_back(index(row),Q(1));unit_orbits.push_back(std::move(units));
        for(size_t h=1;h<g.order();++h)unit_orbits.push_back(mul(unit_orbits[g.words[h].first],generators[g.words[h].second],&opt->pool));
        std::vector<Vec> basis;
        for(auto seed:seeds){std::vector<Vec> copy;int_t den=1,divisor=0;
            std::vector<Q> projected(n);size_t pos=selected_rows.at(seed.row);
            for(size_t h=0;h<g.order();++h)if(ms[seed.sector].primitive[h]!=0)
                for(auto [j,x]:unit_orbits[h][pos])projected[j]+=ms[seed.sector].primitive[h]*x;
            Vec v;for(size_t j=0;j<n;++j)if(projected[j]!=0)v.push_back(index(j),projected[j]);
            std::vector<Vec> orbit(g.order());std::vector<bool> ready(g.order(),false);orbit[0]=v;ready[0]=true;
            std::function<const Vec&(size_t)> image=[&](size_t h)->const Vec&{if(!ready[h]){orbit[h]=carrier_product(image(g.words[h].first),generators[g.words[h].second]);ready[h]=true;}return orbit[h];};
            for(auto h:ms[seed.sector].orbit){copy.push_back(image(h));for(auto [j,x]:copy.back())den=LCM(den,x.den());}
            for(const auto& row:copy)for(auto [j,x]:row){int_t z=den;z/=x.den();z*=x.num();divisor=GCD(divisor,z);}
            require(divisor!=0,"zero carrier-frame generator");for(auto& row:copy){sparse_vec_rescale(row,Q(den,divisor),QQ);basis.push_back(std::move(row));}}
        result.basis=rows(basis,n);
        if(timing)timing->lift+=seconds(started);started=std::chrono::steady_clock::now();
        if(verify_full)verify_adaptation(result,g,ms,generators);if(timing)timing->check+=seconds(started);return result;
    }throw std::runtime_error("symrep: no good modular chart for the carrier symmetry frame");
}


#include "symrep_orbits.hpp"
inline std::string weight_name(size_t w){return "w"+std::to_string(w);}
inline std::vector<Mat> load_carrier_actions(const fs::path& root,size_t w,const Group& g,size_t n,rref_option* opt=nullptr){
    if(!fs::exists(root/(weight_name(w)+"_carrier_g0.wxf"))){
        require(w>1&&opt,"carrier actions require recurrence reconstruction");
        std::ifstream direction(root/"carrier_direction.tsv");std::string name,extra;direction>>name;
        require(bool(direction)&&(name=="forward"||name=="backward")&&!(direction>>extra),"invalid carrier direction metadata");
        bool backward=name=="backward";auto carrier=read(root/(weight_name(w)+"_carrier.wxf"));
        require(carrier.rank()==3&&carrier.dim(0)==n&&carrier.dim(backward?1:2)==g.generators[0].nrow,"carrier action recurrence dimensions mismatch");
        auto older=load_carrier_actions(root,w-1,g,carrier.dim(backward?2:1),opt);
        return carrier_actions(flatten(carrier),older,g.generators,backward,opt);
    }
    std::vector<Mat> actions;
    for(size_t j=0;j<g.generators.size();++j){auto a=read_matrix(root/(weight_name(w)+"_carrier_g"+std::to_string(j)+".wxf"));
        require(a.nrow==n&&a.ncol==n,"saved carrier generator dimension mismatch");actions.push_back(std::move(a));}
    return actions;
}
inline Adapted load_carrier_frame(const fs::path& root,size_t w,const Group& g,const std::vector<Model>& ms,
    const std::vector<Mat>& actions,rref_option* opt=nullptr){
    auto layout=load_layout(root/(weight_name(w)+"_copies.tsv"),ms);
    if(fs::exists(root/(weight_name(w)+"_orbits.tsv"))){
        require(opt,"orbit materialization needs solver options");
        auto frame=materialize_orbits(load_orbits(root/(weight_name(w)+"_orbits.tsv"),layout,ms),g,ms,actions,opt);
        verify_adaptation(frame,g,ms,actions);return frame;
    }
    Adapted frame{read_matrix(root/(weight_name(w)+"_frame.wxf")),load_layout(root/(weight_name(w)+"_copies.tsv"),ms)};
    require(frame.basis.nrow==actions[0].nrow&&frame.basis.ncol==actions[0].nrow,"saved carrier frame is not square");
    verify_adaptation(frame,g,ms,actions);return frame;
}
inline void save_carrier(const fs::path& root,size_t w,const Tensor& carrier,const Adapted& frame,
    const std::vector<Mat>& actions,const std::vector<Model>& ms){
    auto stem=weight_name(w);write(root/(stem+"_carrier.wxf"),carrier);write(root/(stem+"_frame.wxf"),frame.basis);
    save_layout(root/(stem+"_copies.tsv"),frame.layout,ms);
    for(size_t j=0;j<actions.size();++j)write(root/(stem+"_carrier_g"+std::to_string(j)+".wxf"),actions[j]);
}
inline void check_carrier_action(const Mat& carrier,const std::vector<Mat>& actions,
    const std::vector<Mat>& older,const Group& g,bool backward,rref_option_t opt){
    for(size_t j=0;j<actions.size();++j){auto ambient=backward?kron(g.generators[j],older[j]):kron(older[j],g.generators[j]);
        require(equal(mul(actions[j],carrier,&opt->pool),mul(carrier,ambient,&opt->pool)),"carrier action fails exact ambient intertwiner check");}
}
inline void save_chain_meta(const fs::path& root,const fs::path& prepared,size_t maximum,bool factorized){
    std::ofstream meta(root/"chain.tsv");meta<<(factorized?"symbology-symrep-chain-v3":"symbology-symrep-chain-v1")
        <<" prepared_sha256 "<<native_cache::digest(prepared/"checksums.tsv")<<' '<<maximum<<(factorized?" factorized-orbits":"")<<'\n';
    meta.close();require(bool(meta),"cannot write recurrence metadata");seal(root);
}
// Each equation row is independent during contraction. Assemble rows in
// parallel and merge sorted entries directly, avoiding a tree allocation for
// every nonzero of the large native-coordinate condition matrix.
inline Mat assemble_carrier(const Tensor& seed,const Tensor& condition,bool backward,rref_option_t opt){
    require(seed.rank()==3&&condition.rank()==3,"carrier recurrence requires rank-three tensors");
    size_t letters=condition.dim(0),old=seed.dim(backward?2:1),previous=seed.dim(0),equations=condition.dim(2);
    require(condition.dim(1)==letters&&seed.dim(backward?1:2)==letters,"carrier condition letter-axis mismatch");
    auto conditions=condition_rows(condition);std::vector<Vec> by_letter(letters*old);sparse_tensor<Q,I,SPARSE_COO> entries(seed);
    for(size_t p=0;p<entries.nnz();++p){auto i=entries.index(p);by_letter[i[backward?1:2]*old+i[backward?2:1]].push_back(i[0],entries.val(p));}
    Mat result(old*equations,previous*letters);
    opt->pool.detach_loop(size_t(0),result.nrow,[&](size_t r){size_t q=backward?r/old:r%equations,o=backward?r%old:r/equations;auto& row=result[r];
        size_t capacity=0;for(auto [pair,x]:conditions[q]){size_t known=backward?pair%letters:pair/letters;capacity+=by_letter[known*old+o].nnz();}row.reserve(capacity);
        for(auto [pair,x]:conditions[q]){size_t next=backward?pair/letters:pair%letters,known=backward?pair%letters:pair/letters;
            for(auto [b,y]:by_letter[known*old+o])row.push_back(index(backward?next*previous+b:b*letters+next),x*y);}
        row.sort_indices();size_t used=0;
        for(size_t k=0;k<row.nnz();++k){I col=row(k);Q value=row[k];if(used&&row(used-1)==col)row[used-1]+=value;
            else{row(used)=col;row[used]=value;++used;}}
        row.resize(used);row.canonicalize();
    });opt->pool.wait();return result;
}
// H_w maps the sparse carrier C_w to an irreducible basis. The logical basis
// is always organized; C_w, rather than expanded H_w C_w, enters the next
// contraction. This backend uses full sparse elimination, not corner solves.
inline void extend_carrier(const Prepared& p,const Group& g,const std::vector<Model>& ms,size_t maximum,
    const fs::path& root,rref_option_t opt,bool compare,const fs::path& prepared,const fs::path& resume={},bool streamed_kernel=false,tensor_kernel::Strategy strategy=tensor_kernel::Strategy::streamed){
    auto letters=read_matrix(prepared/"alphabet_basis.wxf"),terminal=read_matrix(prepared/"terminal_basis.wxf"),seed_frame=read_matrix(prepared/"seed_basis.wxf");
    auto terminal_inverse=sparse_mat_inverse(terminal,QQ,opt);std::vector<Mat> terminal_actions;
    for(auto h:g.generator_ids)terminal_actions.push_back(mul(mul(terminal_inverse,action(p.terminal,ms,h)),terminal));
    auto physical_change=p.backward?kron(letters,terminal):kron(terminal,letters);
    auto seed=mul(mul(sparse_mat_inverse(seed_frame,QQ,opt),flatten(p.expansion)),physical_change);
    bool compact=bool(fixed_corner_plans(g,ms));
    auto seed_tensor=tensor(seed,p.expansion.dims());auto actions=carrier_actions(seed,terminal_actions,g.generators,p.backward,opt);
    verify_adaptation({seed_frame,p.seed},g,ms,actions);
    if(compare)check_carrier_action(seed,actions,terminal_actions,g,p.backward,opt);
    const auto seed_actions=actions;
    ulong prime=1ULL<<60;field_t fp(FIELD_Fp,2);
    std::vector<ModMat> mod_actions,mod_letters;
    // Retry only modular degeneracies. Rebuild inherited actions from saved
    // exact carriers; the QQ kernel never needs to be solved again.
    auto rebuild_modular=[&](size_t last){
        for(size_t attempt=0;attempt<32;++attempt){do{prime=n_nextprime(prime,0);}while(g.order()>1&&prime%g.order()!=1);fp=field_t(FIELD_Fp,prime);bool good=true;
            for(const auto* set:std::array<const std::vector<Mat>*,2>{&seed_actions,&g.generators})
                for(const auto& a:*set)for(const auto& row:a.rows)for(auto [j,x]:row)if(x.den()%fp.mod==0)good=false;
            if(!good)continue;
            try{mod_actions.clear();mod_letters.clear();
                for(const auto& a:seed_actions)mod_actions.push_back(modular_matrix(a,fp));
                for(const auto& a:g.generators)mod_letters.push_back(modular_matrix(a,fp));
                for(size_t w=2;w<=last;++w)mod_actions=carrier_actions_mod(flatten(read(root/(weight_name(w)+"_carrier.wxf"))),mod_actions,mod_letters,p.backward,fp,opt);
                return;
            }catch(const BadCarrierPrime&){continue;}
        }throw std::runtime_error("symrep: no good carrier action prime after 32 attempts");
    };
    if(compact)rebuild_modular(1);
    auto condition=condition_tensor(read_matrix(prepared/"condition_independent_rows.wxf"),letters.ncol);
    size_t first=2;std::ofstream stats(root/"timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tkernel_s\taction_s\tframe_s\tverification_s\tcarrier_nnz\tstored_frame_nnz\torbit_seeds\n";
    std::ofstream direction(root/"carrier_direction.tsv");direction<<(p.backward?"backward":"forward")<<'\n';direction.close();require(bool(direction),"cannot write carrier direction");
    if(resume.empty())save_carrier(root,1,seed_tensor,{seed_frame,p.seed},actions,ms);
    else{
        bool factorized=false;auto last=chain_weight(resume,prepared,&factorized);
        require(factorized,"factorized continuation requires a factorized chain");require(maximum>last,"--max-weight must exceed the saved chain weight");first=last+1;
        auto saved=read(resume/"w1_carrier.wxf");require(saved.dims()==seed_tensor.dims()&&equal(flatten(saved),seed),"saved carrier seed differs from prepared seed");
        require(equal(read_matrix(resume/"w1_frame.wxf"),seed_frame),"saved seed frame differs from prepared seed");
        for(const auto& entry:fs::directory_iterator(resume))if(entry.is_regular_file()&&entry.path().filename().string().starts_with("w"))fs::copy_file(entry.path(),root/entry.path().filename());
        std::ifstream oldstats(resume/"timings.tsv");require(bool(oldstats),"missing saved timings");std::string line;std::getline(oldstats,line);bool old_columns=line.ends_with("\tframe_nnz");
        while(std::getline(oldstats,line))if(!line.empty())stats<<line<<(old_columns?"\t0":"")<<'\n';
        seed_tensor=read(resume/(weight_name(last)+"_carrier.wxf"));
        auto older=last==1?p.terminal:load_layout(resume/(weight_name(last-1)+"_copies.tsv"),ms);
        auto previous=load_layout(resume/(weight_name(last)+"_copies.tsv"),ms);
        require(seed_tensor.rank()==3&&seed_tensor.dim(0)==previous.dimension(ms)&&seed_tensor.dim(p.backward?1:2)==letters.ncol&&seed_tensor.dim(p.backward?2:1)==older.dimension(ms),"saved carrier dimensions disagree with layouts");
        if(compact&&last>1){prime=1ULL<<60;rebuild_modular(last);}
        if(!compact||compare)actions=load_carrier_actions(resume,last,g,seed_tensor.dim(0),opt);
        if(fs::exists(resume/(weight_name(last)+"_orbits.tsv")))load_orbits(resume/(weight_name(last)+"_orbits.tsv"),previous,ms);
        else if(!compact||compare)load_carrier_frame(resume,last,g,ms,actions,opt);
    }
    for(size_t w=first;w<=maximum;++w){auto start=std::chrono::steady_clock::now();Mat a;
        if(!streamed_kernel||compare)a=assemble_carrier(seed_tensor,condition,p.backward,opt);double assembly_s=seconds(start);
        if(compare)require(equal(a,assemble(seed_tensor,condition,p.backward)),"parallel carrier assembly differs from complete contraction");
        auto nr=a.nrow,nc=a.ncol,nnz=a.nnz();
        if(streamed_kernel&&!compare){nr=seed_tensor.dim(p.backward?2:1)*condition.dim(2);nc=seed_tensor.dim(0)*letters.ncol;}
        Mat check;if(compare)check=a;
        start=std::chrono::steady_clock::now();auto next=streamed_kernel?tensor_kernel::extension(seed_tensor,condition,p.backward,opt,nullptr,strategy):kernel(std::move(a),opt);
        opt->pool.detach_loop(size_t(0),next.nrow,[&](size_t i){vec_cancel_divisor(next[i]);});opt->pool.wait();next.sort_rows_by_nnz();double kernel_s=seconds(start);
        const size_t dimension=next.nrow,carrier_nnz=next.nnz();bool saved_final=false;
        auto carrier_path=root/(weight_name(w)+"_carrier.wxf");
        std::vector<Mat> next_actions;std::vector<ModMat> next_mod;
        CarrierFrameTimings phases;std::optional<OrbitFrame> orbits;Adapted frame;double action_s=0,frame_s=0;
        if(compact){
            for(size_t attempt=0;attempt<32;++attempt){start=std::chrono::steady_clock::now();bool selecting=false;
                try{
                    // Only a bad-prime retry needs the final QQ carrier again.
                    if(saved_final)next=flatten(read(carrier_path));
                    next_mod=carrier_actions_mod(next,mod_actions,mod_letters,p.backward,fp,opt);action_s+=seconds(start);
                    if(w==maximum&&!compare){
                        if(!saved_final){MatrixTensorView view(next,p.backward?std::vector<size_t>{dimension,letters.ncol,seed_tensor.dim(0)}:std::vector<size_t>{dimension,seed_tensor.dim(0),letters.ncol});
                            write_tensor_view(carrier_path,view);saved_final=true;}
                        next.clear();
                    }
                    start=std::chrono::steady_clock::now();selecting=true;orbits=select_fixed_orbits_mod(g,ms,next_mod,fp,opt,&phases);frame_s+=seconds(start);break;
                }catch(const BadCarrierPrime&){(selecting?frame_s:action_s)+=seconds(start);start=std::chrono::steady_clock::now();rebuild_modular(w-1);action_s+=seconds(start);}
            }
            require(bool(orbits),"no good orbit selection prime after 32 attempts");frame.layout=orbits->layout;
        }
        if(!compact||compare){start=std::chrono::steady_clock::now();next_actions=carrier_actions(next,actions,g.generators,p.backward,opt);action_s+=seconds(start);}
        if(!compact){start=std::chrono::steady_clock::now();frame=carrier_frame(g,ms,next_actions,opt,&phases,compare);frame_s=seconds(start);}
        double verification_s=0;
        if(compare){start=std::chrono::steady_clock::now();
            if(orbits){frame=materialize_orbits(*orbits,g,ms,next_actions,opt);verify_adaptation(frame,g,ms,next_actions);}
            require(!mul(check,transpose(next),&opt->pool).nnz(),"carrier kernel residual is nonzero");
            check_carrier_action(next,next_actions,actions,g,p.backward,opt);
            require(kernel(frame.basis,opt).nrow==0,"singular carrier frame");verification_s=seconds(start);}
        if(!saved_final){auto t=tensor(next,p.backward?std::vector<size_t>{dimension,letters.ncol,seed_tensor.dim(0)}:std::vector<size_t>{dimension,seed_tensor.dim(0),letters.ncol},&opt->pool);
            write(carrier_path,t);if(w<maximum)seed_tensor=std::move(t);}
        save_layout(root/(weight_name(w)+"_copies.tsv"),frame.layout,ms);
        if(orbits)save_orbits(root/(weight_name(w)+"_orbits.tsv"),*orbits);else write(root/(weight_name(w)+"_frame.wxf"),frame.basis);
        stats<<std::setprecision(9)<<w<<'\t'<<nr<<'\t'<<nc<<'\t'<<(streamed_kernel&&!compare?std::string("nan"):std::to_string(nnz))<<'\t'<<dimension<<'\t'<<assembly_s<<'\t'<<kernel_s<<'\t'<<action_s<<'\t'<<frame_s<<'\t'<<verification_s<<'\t'<<carrier_nnz<<'\t'<<(orbits?0:frame.basis.nnz())<<'\t'<<(orbits?std::accumulate(orbits->layout.copies.begin(),orbits->layout.copies.end(),size_t(0)):0)<<std::endl;
        std::cout<<"factorized weight="<<w<<" dimension="<<dimension<<" assembly_s="<<assembly_s<<" kernel_s="<<kernel_s<<" action_s="<<action_s<<" frame_s="<<frame_s<<" modular_s="<<phases.modular<<" selection_s="<<phases.selection<<" lift_s="<<phases.lift<<" check_s="<<phases.check<<std::endl;
        actions=std::move(next_actions);if(compact)mod_actions=std::move(next_mod);
    }
    stats.close();require(bool(stats),"cannot write factorized recurrence timings");save_chain_meta(root,prepared,maximum,true);
}
// Verify the full recursively represented symbol spaces against an independent
// ordinary/archive chain. Keep carrier-to-reference and irrep-to-reference
// changes separate: only the former belongs in the next carrier contraction.
inline void verify_carrier_chain(const fs::path& prepared,const fs::path& chain,const fs::path& reference,
    size_t maximum,const Group& g,const std::vector<Model>& ms,rref_option_t opt,const fs::path& output){
    bool factorized=false;require(maximum<=chain_weight(chain,prepared,&factorized)&&factorized,"expected a factorized chain covering the requested weights");
    auto p=load_prepared(prepared,ms);auto letter=identity(p.alphabet.dimension(ms));auto previous=identity(p.terminal.dimension(ms));
    for(size_t w=1;w<=maximum;++w){auto stem=weight_name(w);auto t=read(chain/(stem+"_carrier.wxf"));
        require(t.rank()==3&&t.dim(p.backward?1:2)==letter.nrow&&t.dim(p.backward?2:1)==previous.nrow,"factorized reference recurrence shape mismatch");
        auto physical=mul(flatten(t),p.backward?kron(letter,previous):kron(previous,letter),&opt->pool);
        auto expected=flatten(read(reference/(stem+".wxf")));require(physical.nrow==expected.nrow,"reference dimension mismatch at "+stem);
        auto change=Chart(expected,opt).coordinates(physical);require(kernel(change,opt).nrow==0,"singular carrier reference basis change at "+stem);
        auto actions=load_carrier_actions(chain,w,g,physical.nrow,opt);auto frame=load_carrier_frame(chain,w,g,ms,actions,opt);
        require(kernel(frame.basis,opt).nrow==0,"singular saved irreducible frame at "+stem);
        write(output/(stem+"_basis_change.wxf"),mul(frame.basis,change,&opt->pool));previous=std::move(change);
        std::cout<<"EXACT_REFERENCE_PASS weight="<<w<<" dimension="<<expected.nrow<<std::endl;
    }seal(output);
}
inline void expand_carrier_chain(const fs::path& prepared,const fs::path& chain,size_t maximum,const Group& g,
    const std::vector<Model>& ms,rref_option_t opt,const fs::path& output){
    bool factorized=false;require(maximum<=chain_weight(chain,prepared,&factorized)&&factorized,"expand requires a factorized chain covering the requested weights");
    auto p=load_prepared(prepared,ms);auto letter_inverse=read_matrix(prepared/"alphabet_inverse.wxf");
    auto previous_inverse=sparse_mat_inverse(read_matrix(prepared/"terminal_basis.wxf"),QQ,opt);
    std::ofstream stats(output/"timings.tsv"),export_stats(output/"expansion_timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tadapt_s\tkernel_s\tordinary_kernel_s\tverification_s\tcorner_columns\ttemplate_s\tcompression_s\trref_s\tshorten_s\tlift_s\n";
    export_stats<<"weight\tdimension\texpansion_s\texpanded_nnz\n";
    for(size_t w=1;w<=maximum;++w){auto start=std::chrono::steady_clock::now();auto stem=weight_name(w);auto t=read(chain/(stem+"_carrier.wxf"));
        require(t.rank()==3&&t.dim(p.backward?1:2)==letter_inverse.nrow&&t.dim(p.backward?2:1)==previous_inverse.nrow,"factorized expansion shape mismatch");
        auto actions=load_carrier_actions(chain,w,g,t.dim(0),opt);auto frame=load_carrier_frame(chain,w,g,ms,actions,opt);
        auto coefficients=mul(mul(frame.basis,flatten(t),&opt->pool),p.backward?kron(letter_inverse,previous_inverse):kron(previous_inverse,letter_inverse),&opt->pool);
        if(w==1)require(equal(coefficients,flatten(p.expansion)),"expanded seed differs from prepared seed");
        write(output/(stem+".wxf"),tensor(coefficients,t.dims()));save_layout(output/(stem+"_copies.tsv"),frame.layout,ms);
        if(w<maximum)previous_inverse=sparse_mat_inverse(frame.basis,QQ,opt);
        if(w>1){stats<<w<<"\tnan\t"<<coefficients.ncol<<"\tnan\t"<<coefficients.nrow;for(size_t j=5;j<16;++j)stats<<"\tnan";stats<<'\n';}
        export_stats<<w<<'\t'<<coefficients.nrow<<'\t'<<seconds(start)<<'\t'<<coefficients.nnz()<<'\n';
        std::cout<<"EXPANDED weight="<<w<<" dimension="<<coefficients.nrow<<" nnz="<<coefficients.nnz()<<std::endl;
    }
    stats.close();export_stats.close();require(bool(stats)&&bool(export_stats),"cannot write expansion timings");save_chain_meta(output,prepared,maximum,false);
}
