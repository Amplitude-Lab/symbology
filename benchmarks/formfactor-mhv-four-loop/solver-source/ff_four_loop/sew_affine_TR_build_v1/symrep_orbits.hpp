#pragma once
struct OrbitFrame {
    Layout layout;
    std::vector<std::vector<size_t>> seeds;
    size_t dimension=0;
    std::vector<size_t> pre_elements;
    std::vector<long> pre_shifts;
};
inline ulong trace_product(const ModMat& a,const ModMat& bt,const field_t& fp){
    ulong tr=0;
    for(size_t i=0;i<a.nrow;++i){size_t x=0,y=0;const auto& left=a[i];const auto& right=bt[i];
        while(x<left.nnz()&&y<right.nnz()){
            if(left(x)<right(y))++x;else if(left(x)>right(y))++y;
            else{tr=nmod_add(tr,nmod_mul(left[x],right[y],fp.mod),fp.mod);++x;++y;}
        }
    }return tr;
}
inline ModMat mod_transpose(const ModMat& a){ModMat out(a.ncol,a.nrow);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])out[j].push_back(index(i),x);return out;}
struct CharacterChart {
    std::vector<size_t> elements;
    Mat inverse;
    CharacterChart(const Group& g,const std::vector<Model>& ms,rref_option_t opt){
        Echelon e;std::vector<Vec> rows;
        for(size_t h=0;h<g.order();++h){Vec row;
            for(size_t s=0;s<ms.size();++s)if(ms[s].character[h]!=0)row.push_back(index(s),ms[s].character[h]);
            if(e.add(row)){elements.push_back(h);rows.push_back(row);}
            if(rows.size()==ms.size())break;
        }
        require(rows.size()==ms.size(),"rational characters are dependent");inverse=chart_inverse(symrep::rows(rows,ms.size()),opt);
    }
    std::vector<size_t> multiplicities(const Group& g,const std::vector<Model>& ms,const std::vector<ModMat>& generators,
        const field_t& fp,rref_option_t opt)const{
        size_t n=generators[0].nrow;std::vector<ModMat> transposed;
        for(const auto& x:generators)transposed.push_back(mod_transpose(x));
        std::map<size_t,ModMat> cache;ModMat unit(n,n);for(size_t i=0;i<n;++i)unit[i].push_back(index(i),ulong(1));cache.emplace(0,std::move(unit));
        std::function<const ModMat&(size_t)> action=[&](size_t h)->const ModMat&{
            for(size_t j=0;j<g.generator_ids.size();++j)if(h==g.generator_ids[j])return generators[j];
            if(auto it=cache.find(h);it!=cache.end())return it->second;
            auto [parent,gen]=g.words[h];return cache.emplace(h,sparse_mat_mul(action(parent),generators[gen],fp,&opt->pool)).first->second;
        };
        std::vector<ulong> traces;
        for(auto h:elements){if(!h){traces.push_back(n);continue;}
            auto [parent,gen]=g.words[h];traces.push_back(trace_product(action(parent),transposed[gen],fp));}
        std::vector<size_t> result;size_t accounted=0;
        for(size_t s=0;s<ms.size();++s){ulong m=0;for(auto [j,x]:inverse[s]){require_carrier_prime(x.den()%fp.mod!=0,"bad character-chart prime");m=nmod_add(m,nmod_mul(x%fp.mod,traces[j],fp.mod),fp.mod);}
            require_carrier_prime(m<=n/ms[s].dimension,"invalid modular character multiplicity");result.push_back(m);accounted+=m*ms[s].dimension;}
        require_carrier_prime(accounted==n,"character multiplicities do not cover the carrier");return result;
    }
};
struct DenseFiniteSpan {
    field_t field;std::vector<ModVec> pivots;std::vector<ulong> values;size_t rank=0;
    DenseFiniteSpan(size_t n,ulong prime):field(FIELD_Fp,prime),pivots(n),values(n){}
    bool insert(const ModVec& v,bool keep=true){std::fill(values.begin(),values.end(),0);for(auto [j,x]:v)values[j]=x;
        for(size_t j=0;j<values.size();++j)if(values[j]){
            if(!pivots[j].nnz()){
                if(keep){ModVec row;ulong inverse=nmod_inv(values[j],field.mod),pre=n_mulmod_precomp_shoup(inverse,field.mod.n);
                    for(size_t k=j;k<values.size();++k)if(values[k])row.push_back(index(k),n_mulmod_shoup(inverse,values[k],pre,field.mod.n));
                    pivots[j]=std::move(row);++rank;}return true;
            }
            ulong factor=values[j],pre=n_mulmod_precomp_shoup(factor,field.mod.n);
            for(auto [k,x]:pivots[j])values[k]=_nmod_sub(values[k],n_mulmod_shoup(factor,x,pre,field.mod.n),field.mod);
        }return false;
    }
};
struct FixedCornerPlan {
    std::vector<std::pair<size_t,long>> fixed;
    size_t variable=0;
    long shift=0;
    bool scalar=false;
};
inline std::optional<std::vector<FixedCornerPlan>> fixed_corner_plans(const Group& g,const std::vector<Model>& ms){
    std::vector<FixedCornerPlan> plans;
    std::vector<Q> one(g.order());one[0]=1;
    for(size_t s=0;s<ms.size();++s){const auto& m=ms[s];FixedCornerPlan plan;
        if(m.dimension==1){plan.scalar=true;
            for(size_t j=0;j<g.generators.size();++j){auto x=m.character[g.generator_ids[j]];require(x==1||x==-1,"nonrational scalar character");plan.fixed.push_back({j,x==1?1:-1});}
            plans.push_back(plan);continue;
        }
        if(m.division_dimension!=m.degree)return {};
        auto q=one;
        for(size_t j=0;j<g.generators.size();++j){size_t h=g.generator_ids[j];if(g.table[h][h])continue;
            bool commute=true;for(auto [k,value]:plan.fixed)if(g.table[h][g.generator_ids[k]]!=g.table[g.generator_ids[k]][h])commute=false;
            if(!commute)continue;
            auto translated=g.translate(m.primitive,h);long sign=0;
            if(key(rows({translated},g.order()))==key(rows({g.translate(m.primitive,0)},g.order())))sign=1;
            else{auto neg=translated;sparse_vec_rescale(neg,Q(-1),QQ);if(key(rows({neg},g.order()))==key(rows({g.translate(m.primitive,0)},g.order())))sign=-1;}
            if(sign){std::vector<Q> factor(g.order());factor[0]=Q(1,2);factor[h]+=Q(sign,2);q=g.product(q,factor);plan.fixed.push_back({j,sign});}
        }
        if(g.product(q,m.primitive)!=m.primitive||g.product(m.primitive,q)!=m.primitive||row_basis(g.evaluate(q,m.actions)).nrow!=m.division_dimension)return {};
        bool found=false;
        for(size_t j=0;j<g.generators.size()&&!found;++j){std::vector<Q> u(g.order());u[g.generator_ids[j]]=1;
            auto middle=g.product(g.product(q,u),q),corner=g.product(g.product(m.primitive,u),m.primitive);
            Echelon degree;auto power=m.primitive;
            for(size_t k=0;k<m.division_dimension;++k){degree.add(g.translate(power,0));power=g.product(power,corner);}
            if(degree.pivots.size()!=m.division_dimension)continue;
            for(long shift:{0L,1L,-1L}){auto filter=middle;for(size_t h=0;h<g.order();++h)filter[h]-=Q(shift)*q[h];bool valid=true;
                for(size_t t=0;t<ms.size();++t){size_t rank=row_basis(g.evaluate(filter,ms[t].actions)).nrow;
                    if(rank!=(t==s?m.division_dimension:0)){valid=false;break;}}
                if(valid){plan.variable=j;plan.shift=shift;found=true;break;}}
        }
        if(!found)return {};plans.push_back(plan);
    }return plans;
}
inline ModMat fixed_equations(const std::vector<ModMat>& gs,const FixedCornerPlan& plan,const field_t& fp){
    size_t n=gs[0].nrow;ModMat out(n*plan.fixed.size(),n);size_t off=0;
    for(auto [j,sign]:plan.fixed){for(size_t i=0;i<n;++i){out[off+i]=gs[j][i];ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_add_mul(out[off+i],unit,sign==1?fp.mod.n-1:ulong(1),fp);}off+=n;}return out;
}
// A simple residue-field embedding of End_G(U) certifies independence of
// whole QQ copies: a determinant nonzero in one embedding is nonzero in the
// number field and hence invertible there. Other Galois components need not
// undergo a second rank selection. The cyclic corner operator was certified
// to generate the full degree-e field by fixed_corner_plans().
inline std::optional<ulong> corner_embedding(const Group& g,const Model& m,size_t variable,
    const field_t& fp,rref_option_t opt){
    if(m.division_dimension!=m.degree)return {};
    auto p=g.evaluate(m.primitive,m.actions),basis=row_basis(p);
    auto a=Chart(basis,opt).coordinates(mul(mul(basis,m.actions[g.generator_ids[variable]]),p));
    fmpq_mat_t matrix;fmpq_mat_init(matrix,a.nrow,a.ncol);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:a[i])fmpq_set(fmpq_mat_entry(matrix,i,j),x.data());
    Poly characteristic;fmpq_mat_charpoly(characteristic.p,matrix);fmpq_mat_clear(matrix);
    std::vector<ulong> coefficients;
    for(size_t j=0;j<=a.nrow;++j){Q x;fmpq_poly_get_coeff_fmpq(x.data(),characteristic.p,j);
        require_carrier_prime(x.den()%fp.mod!=0,"bad corner-polynomial denominator prime");coefficients.push_back(x%fp.mod);}
    nmod_poly_t polynomial,derivative;nmod_poly_init(polynomial,fp.mod.n);nmod_poly_init(derivative,fp.mod.n);
    for(size_t j=0;j<coefficients.size();++j)nmod_poly_set_coeff_ui(polynomial,j,coefficients[j]);
    nmod_poly_derivative(derivative,polynomial);nmod_poly_factor_t roots;nmod_poly_factor_init(roots);nmod_poly_roots(roots,polynomial,1);
    std::optional<ulong> result;
    for(slong i=0;i<roots->num;++i){auto f=roots->p+i;if(nmod_poly_degree(f)!=1)continue;
        ulong root=nmod_neg(nmod_div(nmod_poly_get_coeff_ui(f,0),nmod_poly_get_coeff_ui(f,1),fp.mod),fp.mod);
        if(nmod_poly_evaluate_nmod(derivative,root)&&(!result||root<*result))result=root;}
    nmod_poly_factor_clear(roots);nmod_poly_clear(polynomial);nmod_poly_clear(derivative);return result;
}
inline OrbitFrame select_fixed_orbits_mod(const Group& g,const std::vector<Model>& ms,const std::vector<ModMat>& gs,
    const field_t& fp,rref_option_t opt,CarrierFrameTimings* timing=nullptr){
    auto plans=fixed_corner_plans(g,ms);require(bool(plans),"no fixed-corner plan");size_t n=gs[0].nrow;
    OrbitFrame answer;answer.dimension=n;answer.seeds.resize(ms.size());answer.layout.copies.resize(ms.size());answer.pre_elements.resize(ms.size());answer.pre_shifts.resize(ms.size());if(!n)return answer;
    auto start=std::chrono::steady_clock::now();CharacterChart chart(g,ms,opt);ulong prime=fp.mod.n;
    auto wanted=chart.multiplicities(g,ms,gs,fp,opt);
    if(timing)timing->modular=seconds(start);start=std::chrono::steady_clock::now();
    // All scalar sectors can share the fixed space of one group element
    // when its 1-eigenspace contains exactly the one-dimensional irreps.
    std::optional<size_t> scalar_element;
    for(size_t h=1;h<g.order();++h){bool valid=true;
        for(const auto& m:ms){auto a=m.actions[h];for(size_t i=0;i<a.nrow;++i){std::map<I,Q> row;for(auto [j,x]:a[i])row[j]=x;row[index(i)]-=1;a[i]=packed(row);}
            size_t rank=row_basis(a).nrow;if(rank!=(m.dimension==1?0:m.dimension)){valid=false;break;}}
        if(valid){scalar_element=h;break;}
    }
    size_t scalar_total=0;for(size_t s=0;s<ms.size();++s)if(ms[s].dimension==1)scalar_total+=wanted[s];
    auto scalars=[&](rref_option* opt){
    if(scalar_element&&scalar_total){std::map<size_t,ModMat> cache;ModMat unit(n,n);for(size_t i=0;i<n;++i)unit[i].push_back(index(i),ulong(1));cache.emplace(0,unit);
        std::function<const ModMat&(size_t)> image=[&](size_t h)->const ModMat&{for(size_t j=0;j<gs.size();++j)if(h==g.generator_ids[j])return gs[j];if(auto it=cache.find(h);it!=cache.end())return it->second;auto [parent,j]=g.words[h];return cache.emplace(h,sparse_mat_mul(image(parent),gs[j],fp,&opt->pool)).first->second;};
        auto equations=image(*scalar_element);for(size_t i=0;i<n;++i)sparse_vec_sub_mul(equations[i],unit[i],ulong(1),fp);
        equations.sort_rows_by_nnz();auto piv=sparse_mat_rref(equations,fp,opt);std::vector<bool> used(n,false);
        for(auto& batch:piv)for(auto [r,c]:batch)used[c]=true;std::vector<size_t> free;
        for(size_t i=0;i<n;++i)if(!used[i])free.push_back(i);
        size_t expected=0;for(size_t s=0;s<ms.size();++s)if(ms[s].dimension==1)expected+=wanted[s];require_carrier_prime(free.size()==expected,"wrong total scalar eigenspace dimension");
        auto inclusion=sparse_mat_rref_kernel(equations,piv,fp,opt);equations.clear();for(auto& row:inclusion.rows)sparse_vec_rescale(row,fp.mod.n-1,fp);
        std::vector<ModMat> small;
        for(size_t j=0;j<gs.size();++j){std::optional<Q> value;bool scalar=true;
            for(const auto& m:ms)if(m.dimension==1){auto x=m.character[g.generator_ids[j]];if(value&&*value!=x)scalar=false;value=x;}
            if(scalar&&value){ModMat diagonal(free.size(),free.size());for(size_t i=0;i<free.size();++i)diagonal[i].push_back(index(i),*value%fp.mod);small.push_back(std::move(diagonal));continue;}
            ModMat left(free.size(),n);for(size_t i=0;i<free.size();++i)left[i]=gs[j][free[i]];
            small.push_back(sparse_mat_mul(left,inclusion,fp,&opt->pool));}
        for(size_t s=0;s<ms.size();++s)if(ms[s].dimension==1&&wanted[s]){auto eq=fixed_equations(small,plans->at(s),fp);eq.sort_rows_by_nnz();auto pv=sparse_mat_rref_forward(eq,fp,opt);
            std::vector<bool> occupied(free.size(),false);for(auto& batch:pv)for(auto [r,c]:batch)occupied[c]=true;
            for(size_t i=0;i<free.size();++i)if(!occupied[i])answer.seeds[s].push_back(free[i]);
            require_carrier_prime(answer.seeds[s].size()==wanted[s],"wrong scalar-sector dimension");answer.layout.copies[s]=wanted[s];}

    }
    };
    std::vector<size_t> higher;
    for(size_t s=0;s<ms.size();++s)if(!plans->at(s).scalar&&wanted[s])higher.push_back(s);
    auto nonscalars=[&](rref_option* opt){
    if(!higher.empty()){
        size_t variable=plans->at(higher[0]).variable;std::optional<ulong> root;

        // A plain generator eigenspace may serve several rational sectors.
        // Search its eigenvalues, then verify the dimensions on every model.
        const auto& ma=ms[higher[0]].actions[g.generator_ids[variable]];
        fmpq_mat_t matrix;fmpq_mat_init(matrix,ma.nrow,ma.ncol);
        for(size_t i=0;i<ma.nrow;++i)for(auto [j,x]:ma[i])fmpq_set(fmpq_mat_entry(matrix,i,j),x.data());
        Poly poly;fmpq_mat_charpoly(poly.p,matrix);fmpq_mat_clear(matrix);
        std::vector<ulong> coefficients;
        for(size_t j=0;j<=ma.nrow;++j){Q x;fmpq_poly_get_coeff_fmpq(x.data(),poly.p,j);require_carrier_prime(x.den()%fp.mod!=0,"bad spectral-polynomial denominator prime");coefficients.push_back(x%fp.mod);}
        nmod_poly_t modular;nmod_poly_init(modular,fp.mod.n);
        for(size_t j=0;j<coefficients.size();++j)nmod_poly_set_coeff_ui(modular,j,coefficients[j]);
        nmod_poly_factor_t factors;nmod_poly_factor_init(factors);nmod_poly_roots(factors,modular,1);std::vector<ulong> roots;
        for(slong i=0;i<factors->num;++i){auto f=factors->p+i;if(nmod_poly_degree(f)==1)roots.push_back(nmod_neg(nmod_div(nmod_poly_get_coeff_ui(f,0),nmod_poly_get_coeff_ui(f,1),fp.mod),fp.mod));}
        nmod_poly_factor_clear(factors);nmod_poly_clear(modular);std::sort(roots.begin(),roots.end());
        std::vector<std::pair<size_t,long>> project;
        for(auto [j,sign]:plans->at(higher[0]).fixed)if(g.table[g.generator_ids[j]][g.generator_ids[variable]]!=g.table[g.generator_ids[variable]][g.generator_ids[j]])project.push_back({j,sign});
        bool compatible=true;
        for(auto s:higher){std::vector<std::pair<size_t,long>> local;for(auto [j,sign]:plans->at(s).fixed)if(g.table[g.generator_ids[j]][g.generator_ids[variable]]!=g.table[g.generator_ids[variable]][g.generator_ids[j]])local.push_back({j,sign});
            if(plans->at(s).variable!=variable||local!=project)compatible=false;}
        if(compatible)for(auto candidate:roots){bool valid=true;
            for(size_t t=0;t<ms.size();++t){if(!wanted[t])continue;size_t d=ms[t].dimension;
                auto equation=modular_matrix(ms[t].actions[g.generator_ids[variable]],fp);DenseFiniteSpan span(d,prime);
                for(size_t i=0;i<d;++i){ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_sub_mul(equation[i],unit,candidate,fp);span.insert(equation[i]);}
                bool target=std::find(higher.begin(),higher.end(),t)!=higher.end();
                if(span.rank!=d-(target?1:0)){valid=false;break;}
                if(target){auto p=modular_matrix(g.evaluate(ms[t].primitive,ms[t].actions),fp);for(const auto& row:p.rows)span.insert(row);if(span.rank!=d){valid=false;break;}}
            }
            if(valid){root=candidate;break;}
        }
        if(root){auto equation=gs[variable];for(size_t i=0;i<n;++i){ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_sub_mul(equation[i],unit,*root,fp);}
            equation.sort_rows_by_nnz();auto piv=sparse_mat_rref(equation,fp,opt);std::vector<bool> used(n,false);
            for(const auto& batch:piv)for(auto [r,c]:batch)used[c]=true;std::vector<size_t> free;
            for(size_t i=0;i<n;++i)if(!used[i])free.push_back(i);size_t total=0;for(auto s:higher)total+=wanted[s];
            require_carrier_prime(free.size()==total,"wrong shared spectral dimension");
            auto inclusion=sparse_mat_rref_kernel(equation,piv,fp,opt);equation.clear();for(auto& row:inclusion.rows)sparse_vec_rescale(row,fp.mod.n-1,fp);
            std::vector<ModMat> small(gs.size(),ModMat(total,total));std::set<size_t> needed;
            for(auto s:higher)for(auto [j,sign]:plans->at(s).fixed)if(g.table[g.generator_ids[j]][g.generator_ids[variable]]==g.table[g.generator_ids[variable]][g.generator_ids[j]])needed.insert(j);
            for(auto j:needed){ModMat left(total,n);for(size_t i=0;i<total;++i)left[i]=gs[j][free[i]];small[j]=sparse_mat_mul(left,inclusion,fp,&opt->pool);}
            std::vector<size_t> candidates=free;
            auto projected=[&](const std::vector<size_t>& selected){ModMat left(selected.size(),n);for(size_t i=0;i<selected.size();++i)left[i].push_back(index(selected[i]),ulong(1));
                for(auto [j,sign]:project){auto next=sparse_mat_mul(left,gs[j],fp,&opt->pool);
                    for(size_t i=0;i<left.nrow;++i){if(sign<0)sparse_vec_rescale(next[i],fp.mod.n-1,fp);sparse_vec_add_mul(next[i],left[i],ulong(1),fp);}left=std::move(next);}
                return sparse_mat_mul(left,inclusion,fp,&opt->pool);};
            auto pairing=projected(candidates);
            auto rank_of=[&](ModMat matrix){auto pivots=sparse_mat_rref_forward(matrix,fp,opt);size_t count=0;for(const auto& batch:pivots)count+=batch.size();return count;};
            size_t count=rank_of(pairing),pre_element=0;long pre_shift=0;
            // If F U F=U^-1, the rows e_i(U-a)q paired with the U-eigenline
            // equal (lambda^-1-a)Y+(lambda-lambda^-1)I. This matrix pencil
            // is regular when lambda != lambda^-1. It repairs a singular
            // free-coordinate pairing without expanding to all carrier rows.
            if(count<total&&project.size()==1){size_t h=g.generator_ids[project[0].first],u=g.generator_ids[variable];
                ulong inverse=nmod_inv(*root,fp.mod);
                if(*root!=inverse&&g.table[u][h]==g.table[h][g.inverses[u]]){
                    for(long shift:{0L,1L,-1L,2L,-2L,3L,-3L}){auto candidate=pairing;ulong a=Q(shift)%fp.mod;
                        for(size_t i=0;i<total;++i){sparse_vec_rescale(candidate[i],nmod_sub(inverse,a,fp.mod),fp);ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_add_mul(candidate[i],unit,nmod_sub(*root,inverse,fp.mod),fp);}
                        if(rank_of(candidate)==total){pairing=std::move(candidate);count=total;pre_element=u;pre_shift=shift;break;}
                    }
                }
            }
            if(count<total){candidates.resize(n);std::iota(candidates.begin(),candidates.end(),0);pairing=projected(candidates);}
            for(auto s:higher){auto plan=plans->at(s);std::erase_if(plan.fixed,[&](auto x){return !needed.count(x.first);});
                auto eq=fixed_equations(small,plan,fp);eq.sort_rows_by_nnz();auto pv=sparse_mat_rref(eq,fp,opt);size_t rank=0;for(const auto& batch:pv)rank+=batch.size();
                require_carrier_prime(total-rank==wanted[s],"wrong small spectral sector dimension");
                auto corner=sparse_mat_rref_kernel(eq,pv,fp,opt);auto rows=mod_transpose(sparse_mat_mul(pairing,corner,fp,&opt->pool));
                auto pivots=sparse_mat_rref_forward(rows,fp,opt);for(const auto& batch:pivots)for(auto [r,c]:batch)answer.seeds[s].push_back(candidates[c]);
                require_carrier_prime(answer.seeds[s].size()==wanted[s],"shared spectral pairing lost rank");
                std::sort(answer.seeds[s].begin(),answer.seeds[s].end());answer.layout.copies[s]=wanted[s];answer.pre_elements[s]=pre_element;answer.pre_shifts[s]=pre_shift;}
        }
    }
    };
    size_t available=opt->pool.get_thread_count();
    if(n>=128&&available>=2&&scalar_element&&scalar_total&&!higher.empty()){
        auto work=[&](bool scalar,size_t workers){rref_option_t local;local->pool.reset(workers);local->verbose=false;
            local->method=opt->method;local->col_weight=opt->col_weight;local->eliminate_one_nnz=opt->eliminate_one_nnz;
            if(scalar)scalars(local);else nonscalars(local);};
        auto scalar_job=std::async(std::launch::async,work,true,available/2);
        auto other_job=std::async(std::launch::async,work,false,available-available/2);
        scalar_job.get();other_job.get();
    }else{scalars(opt);nonscalars(opt);}
    auto sector=[&](size_t s,rref_option* work){const auto& plan=plans->at(s);
        auto equations=fixed_equations(gs,plan,fp);equations.sort_rows_by_nnz();
        auto piv=plan.scalar?sparse_mat_rref_forward(equations,fp,work):sparse_mat_rref(equations,fp,work);
        std::vector<bool> used(n,false);for(const auto& batch:piv)for(auto [r,c]:batch)used[c]=true;
        std::vector<size_t> free;for(size_t i=0;i<n;++i)if(!used[i])free.push_back(i);
        if(plan.scalar){require_carrier_prime(free.size()==wanted[s],"wrong scalar eigenspace dimension");answer.seeds[s]=free;answer.layout.copies[s]=free.size();
            return;}
        auto inclusion=sparse_mat_rref_kernel(equations,piv,fp,work);equations.clear();for(auto& row:inclusion.rows)sparse_vec_rescale(row,fp.mod.n-1,fp);require(inclusion.nrow==n&&inclusion.ncol==free.size(),"wrong fixed-space inclusion shape");
        auto v=sparse_mat_mul(gs[plan.variable],inclusion,fp,&work->pool);std::vector<std::pair<size_t,long>> filters;
        for(auto [j,sign]:plan.fixed)if(g.table[g.generator_ids[j]][g.generator_ids[plan.variable]]!=g.table[g.generator_ids[plan.variable]][g.generator_ids[j]])filters.push_back({j,sign});
        ModMat rep(free.size(),free.size());
        for(size_t k=0;k<filters.size();++k){auto [j,sign]=filters[k];bool last=k+1==filters.size();
            auto left=last?ModMat(free.size(),n):gs[j];if(last)for(size_t i=0;i<free.size();++i)left[i]=gs[j][free[i]];
            auto product=sparse_mat_mul(left,v,fp,&work->pool);
            for(size_t i=0;i<product.nrow;++i){if(sign<0)sparse_vec_rescale(product[i],fp.mod.n-1,fp);
                sparse_vec_add_mul(product[i],v[last?free[i]:i],ulong(1),fp);sparse_vec_rescale(product[i],nmod_inv(2,fp.mod),fp);}
            if(last)rep=std::move(product);else v=std::move(product);
        }
        if(filters.empty())for(size_t i=0;i<free.size();++i)rep[i]=v[free[i]];
        if(auto root=corner_embedding(g,ms[s],plan.variable,fp,work)){
            auto eigen=rep;for(size_t i=0;i<free.size();++i){ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_sub_mul(eigen[i],unit,*root,fp);}
            eigen.sort_rows_by_nnz();auto pivots=sparse_mat_rref_forward(eigen,fp,work);std::vector<bool> occupied(free.size(),false);
            for(const auto& batch:pivots)for(auto [r,c]:batch)occupied[c]=true;
            for(size_t i=0;i<free.size();++i)if(!occupied[i])answer.seeds[s].push_back(free[i]);
            require_carrier_prime(answer.seeds[s].size()==wanted[s],"wrong endomorphism-embedding dimension");answer.layout.copies[s]=wanted[s];return;
        }
        auto filter=rep;for(size_t i=0;i<free.size();++i){if(plan.shift){ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_add_mul(filter[i],unit,plan.shift>0?fp.mod.n-1:ulong(1),fp);}}
        DenseFiniteSpan span(free.size(),prime);std::vector<size_t> order(free.size());std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return filter[a].nnz()<filter[b].nnz();});
        for(auto i:order){if(!span.insert(filter[i],false))continue;size_t before=span.rank;auto vector=filter[i];
            for(size_t k=0;k<ms[s].division_dimension;++k){span.insert(vector);if(k+1<ms[s].division_dimension)vector=modular_product(vector,rep,fp);}
            require_carrier_prime(span.rank-before==ms[s].division_dimension,"nonfree fixed-corner generator");
            answer.seeds[s].push_back(free[i]);++answer.layout.copies[s];if(answer.layout.copies[s]==wanted[s])break;
        }
        require_carrier_prime(answer.layout.copies[s]==wanted[s],"incomplete fixed-corner frame");

     };
    std::vector<size_t> pending;for(size_t s=0;s<ms.size();++s)if(wanted[s]&&answer.layout.copies[s]!=wanted[s])pending.push_back(s);
    size_t workers=std::min(pending.size(),size_t(opt->pool.get_thread_count()));
    if(workers>1&&n>=128){std::vector<std::future<void>> jobs;
        for(size_t worker=0;worker<workers;++worker)jobs.push_back(std::async(std::launch::async,[&,worker]{
            rref_option_t local;local->pool.reset(std::max(size_t(1),size_t(opt->pool.get_thread_count())/workers));
            local->verbose=false;local->method=opt->method;local->col_weight=opt->col_weight;local->eliminate_one_nnz=opt->eliminate_one_nnz;
            for(size_t k=worker;k<pending.size();k+=workers)sector(pending[k],local);
        }));
        for(auto& job:jobs)job.get();
    }else for(auto s:pending)sector(s,opt);
    if(timing)timing->selection=seconds(start);return answer;
}

inline OrbitFrame select_fixed_orbits(const Group& g,const std::vector<Model>& ms,const std::vector<Mat>& generators,
    rref_option_t opt,CarrierFrameTimings* timing=nullptr){
    ulong prime=1ULL<<60;
    for(size_t attempt=0;attempt<32;++attempt){do{prime=n_nextprime(prime,0);}while(g.order()>1&&prime%g.order()!=1);field_t fp(FIELD_Fp,prime);bool good=true;
        for(const auto& a:generators)for(const auto& row:a.rows)for(auto [j,x]:row)if(x.den()%fp.mod==0)good=false;
        if(!good)continue;
        std::vector<ModMat> gs;for(const auto& a:generators)gs.push_back(modular_matrix(a,fp));
        try{return select_fixed_orbits_mod(g,ms,gs,fp,opt,timing);}catch(const BadCarrierPrime&){continue;}
    }throw std::runtime_error("symrep: no good denominator prime for fixed-corner selection");
}

inline Adapted materialize_orbits(const OrbitFrame& frame,const Group& g,const std::vector<Model>& ms,
    const std::vector<Mat>& generators,rref_option_t opt){
    const size_t n=frame.dimension;std::map<size_t,size_t> selected;
    for(const auto& sector:frame.seeds)for(auto row:sector)if(!selected.count(row)){size_t pos=selected.size();selected.emplace(row,pos);}
    std::vector<Mat> units;Mat initial(selected.size(),n);
    for(auto [row,pos]:selected)initial[pos].push_back(index(row),Q(1));units.push_back(std::move(initial));
    for(size_t h=1;h<g.order();++h)units.push_back(mul(units[g.words[h].first],generators[g.words[h].second],&opt->pool));
    std::vector<Vec> basis;
    for(size_t s=0;s<ms.size();++s)for(auto seed:frame.seeds[s]){
        std::vector<Q> values(n);size_t pos=selected.at(seed);
        auto coefficients=ms[s].primitive;
        if(!frame.pre_elements.empty()&&(frame.pre_elements[s]||frame.pre_shifts[s])){std::vector<Q> pre(g.order());pre[frame.pre_elements[s]]=1;pre[0]-=Q(frame.pre_shifts[s]);coefficients=g.product(pre,coefficients);}
        for(size_t h=0;h<g.order();++h)if(coefficients[h]!=0)for(auto [j,x]:units[h][pos])values[j]+=coefficients[h]*x;
        Vec first;for(size_t j=0;j<n;++j)if(values[j]!=0)first.push_back(index(j),values[j]);
        std::vector<Vec> orbit(g.order());std::vector<bool> ready(g.order(),false);orbit[0]=std::move(first);ready[0]=true;
        std::function<const Vec&(size_t)> image=[&](size_t h)->const Vec&{if(!ready[h]){orbit[h]=carrier_product(image(g.words[h].first),generators[g.words[h].second]);ready[h]=true;}return orbit[h];};
        std::vector<Vec> copy;int_t den=1,divisor=0;
        for(auto h:ms[s].orbit){copy.push_back(image(h));for(auto [j,x]:copy.back())den=LCM(den,x.den());}
        for(const auto& row:copy)for(auto [j,x]:row){int_t z=den;z/=x.den();z*=x.num();divisor=GCD(divisor,z);}
        require(divisor!=0,"zero symmetry-orbit seed");
        for(auto& row:copy){sparse_vec_rescale(row,Q(den,divisor),QQ);basis.push_back(std::move(row));}
    }
    return {rows(basis,n),frame.layout};
}
inline void save_orbits(const fs::path& path,const OrbitFrame& frame){
    std::ofstream out(path);out<<"symbology-irrep-orbits-v2 "<<frame.dimension<<" primitive-integer-copy\n";
    out<<"irrep\tcopy\tcarrier_row\tpre_element\tpre_shift\n";
    for(size_t s=0;s<frame.seeds.size();++s)for(size_t c=0;c<frame.seeds[s].size();++c)
        out<<s<<'\t'<<c<<'\t'<<frame.seeds[s][c]<<'\t'<<(frame.pre_elements.empty()?0:frame.pre_elements[s])<<'\t'<<(frame.pre_shifts.empty()?0:frame.pre_shifts[s])<<'\n';
    out.close();require(bool(out),"cannot write symmetry-orbit frame");
}
inline OrbitFrame load_orbits(const fs::path& path,const Layout& layout,const std::vector<Model>& ms){
    std::ifstream in(path);std::string format,normalization,header;size_t n=0;
    in>>format>>n>>normalization;bool extended=format=="symbology-irrep-orbits-v2";
    require(bool(in)&&(extended||format=="symbology-irrep-orbits-v1")&&normalization=="primitive-integer-copy","invalid symmetry-orbit frame schema");
    require(n==layout.dimension(ms),"symmetry-orbit dimension disagrees with layout");std::getline(in,header);std::getline(in,header);
    require(header==(extended?"irrep\tcopy\tcarrier_row\tpre_element\tpre_shift":"irrep\tcopy\tcarrier_row"),"invalid symmetry-orbit columns");
    OrbitFrame result{layout,{},n};result.seeds.resize(ms.size());result.pre_elements.resize(ms.size());result.pre_shifts.resize(ms.size());
    for(size_t s=0;s<ms.size();++s){std::set<size_t> seen;
        for(size_t c=0;c<layout.copies[s];++c){size_t sector,copy,row,element=0;long shift=0;in>>sector>>copy>>row;if(extended)in>>element>>shift;
            require(bool(in)&&sector==s&&copy==c&&row<n&&element<ms[s].actions.size()&&seen.insert(row).second,"invalid symmetry-orbit record");
            if(c)require(element==result.pre_elements[s]&&shift==result.pre_shifts[s],"inconsistent symmetry-orbit prefilter");
            result.pre_elements[s]=element;result.pre_shifts[s]=shift;result.seeds[s].push_back(row);}}
    require(!(in>>header),"trailing symmetry-orbit records");return result;
}
