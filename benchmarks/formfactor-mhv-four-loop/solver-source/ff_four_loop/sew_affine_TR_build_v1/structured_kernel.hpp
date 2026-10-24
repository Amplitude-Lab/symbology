#pragma once
// Exact homogeneous presolve. Degree-one equations fix a variable to zero;
// degree-two equations identify two variables up to a rational scale. These
// operations preserve the complete kernel, including inconsistent cycles.
#include "symrep.hpp"
#include "certified_kernel.hpp"

namespace structured_kernel {
using namespace symrep;
struct Statistics {
    size_t input_rows=0,input_cols=0,input_nnz=0,core_cols=0,core_nnz=0,passes=0;
    double presolve_s=0,solve_s=0,lift_s=0;
};
class Relations {
    std::vector<size_t> parent,size;
    std::vector<Q> scale;
    std::vector<bool> zero;
public:
    explicit Relations(size_t n):parent(n),size(n,1),scale(n,Q(1)),zero(n,false){std::iota(parent.begin(),parent.end(),0);}
    // x_i = scale_i * x_root. Union by size bounds recursion depth.
    size_t root(size_t i){if(parent[i]!=i&&parent[parent[i]]!=parent[i]){auto p=parent[i];parent[i]=root(p);scale[i]*=scale[p];}return parent[i];}
    bool consume(const Vec& row){
        if(row.nnz()==1){auto r=root(row(0));bool changed=!zero[r];zero[r]=true;return changed;}
        if(row.nnz()!=2)return false;
        size_t a=row(0),b=row(1),ra=root(a),rb=root(b);
        if(zero[ra]||zero[rb]){bool changed=!zero[ra]||!zero[rb];zero[ra]=zero[rb]=true;return changed;}
        Q av=row[0]*scale[a],bv=row[1]*scale[b];
        if(ra==rb){if(av+bv==0)return false;bool changed=!zero[ra];zero[ra]=true;return changed;}
        if(size[ra]>size[rb]){std::swap(ra,rb);std::swap(av,bv);}
        parent[ra]=rb;scale[ra]=-bv/av;size[rb]+=size[ra];zero[rb]=zero[ra]||zero[rb];return true;
    }
    void compress(Vec& row){
        size_t n=0;
        for(size_t p=0;p<row.nnz();++p){size_t c=row(p),r=root(c);if(zero[r])continue;
            Q value=row[p]*scale[c];row(n)=index(r);row[n]=std::move(value);++n;}
        row.resize(n);normalize(row);
    }
    std::pair<size_t,Q> image(size_t c){auto r=root(c);return {r,zero[r]?Q(0):scale[c]};}
    bool is_zero(size_t c){return zero[root(c)];}
    size_t dimension(){size_t n=0;for(size_t i=0;i<parent.size();++i)if(root(i)==i&&!zero[i])++n;return n;}
};
// A replayable row source may skip rows whose support bound exceeds limit.
// For the final pass limit=SIZE_MAX, it must emit EVERY transformed row.
// generate(relations, limit, emit) applies the current coordinate map before
// allocation; this lets factored operators avoid building discarded entries.
template<class Generate,class Verify=std::nullptr_t>
inline Mat solve_stream(size_t ncols,Generate&& generate,rref_option_t opt,Statistics* stats=nullptr,Verify verify=nullptr,size_t batch_rows=0){
    Statistics s;s.input_cols=ncols;Relations relations(ncols);
    auto start=std::chrono::steady_clock::now();
    for(size_t pass=0;pass<8;++pass){size_t changed=0;
        generate(relations,size_t(2),[&](Vec&& row){normalize(row);if(row.nnz()<=2)changed+=relations.consume(row);});
        ++s.passes;std::cout<<"presolve_pass="<<s.passes<<" relations="<<changed<<" remaining_columns="<<relations.dimension()<<std::endl;
        if(!changed)break;
    }
    std::vector<I> column(ncols,-1),image(ncols,-1);std::vector<Q> coefficients(ncols);size_t n=0;
    for(size_t c=0;c<ncols;++c){auto [r,x]=relations.image(c);coefficients[c]=x;if(x!=0&&column[r]<0)column[r]=index(n++);}
    for(size_t c=0;c<ncols;++c){auto [r,x]=relations.image(c);if(x!=0)image[c]=column[r];}
    std::vector<std::vector<I>> originals(n);
    for(size_t c=0;c<ncols;++c)if(image[c]>=0)originals[image[c]].push_back(index(c));
    auto lift=[&](const Mat& core){Mat result(core.nrow,ncols);
        opt->pool.detach_loop(size_t(0),core.nrow,[&](size_t i){size_t count=0;
            for(auto c:core[i].index_span())count+=originals[c].size();result[i].reserve(count);
            for(auto [c,x]:core[i])for(auto original:originals[c])result[i].push_back(original,x*coefficients[original]);
            // Relation components are disjoint, so no merge is required.
            result[i].sort_indices();
        });opt->pool.wait();
        return result;};
    auto generate_reduced=[&](auto emit){s.core_nnz=0;generate(relations,SIZE_MAX,[&](Vec&& row){normalize(row);if(!row.nnz())return;
        for(size_t j=0;j<row.nnz();++j)row(j)=column[row(j)];normalize(row);row.compress();s.core_nnz+=row.nnz();emit(std::move(row));});};
    s.core_cols=n;Mat core,verified_lift;bool have_verified_lift=false;
    if constexpr(std::is_same_v<Verify,std::nullptr_t>){Mat a(0,n);generate_reduced([&](Vec&& row){a.rows.push_back(std::move(row));});a.nrow=a.rows.size();
        s.presolve_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
        std::cout<<"core_rows="<<a.nrow<<" core_cols="<<n<<" core_nnz="<<s.core_nnz<<std::endl;
        core=symrep::kernel(std::move(a),opt);
    }else{
        s.presolve_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
        auto test=[&](const Mat& candidate,ulong prime,const int_t& modulus){auto expanded=lift(candidate);
            opt->pool.detach_loop(size_t(0),expanded.nrow,[&](size_t i){vec_cancel_divisor(expanded[i]);});opt->pool.wait();
            bool valid;
            if constexpr(std::is_invocable_r_v<bool,Verify,const Mat&,ulong,const int_t&>)valid=verify(expanded,prime,modulus);
            else if constexpr(std::is_invocable_r_v<bool,Verify,const Mat&,ulong>)valid=verify(expanded,prime);
            else valid=verify(expanded);
            if(valid){verified_lift=std::move(expanded);have_verified_lift=true;}return valid;};
        auto build_mod=[&](const field_t& fp,auto emit){
            if constexpr(requires{generate.modular(relations,column,fp,emit,&opt->pool);}){
                s.core_nnz=0;generate.modular(relations,column,fp,[&](auto&& row){s.core_nnz+=row.nnz();emit(std::move(row));},&opt->pool);
            }else if constexpr(requires{generate.modular(relations,column,fp,emit);}){
                s.core_nnz=0;generate.modular(relations,column,fp,[&](auto&& row){s.core_nnz+=row.nnz();emit(std::move(row));});
            }
            else generate_reduced([&](Vec&& row){for(auto [j,x]:row)if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("input denominator");emit(row%fp.mod);});
        };
        core=certified_kernel::solve(n,generate_reduced,test,opt,build_mod,batch_rows);
    }
    s.solve_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
    auto result=have_verified_lift?std::move(verified_lift):lift(core);
    s.lift_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();if(stats)*stats=s;return result;
}
inline Mat solve(Mat a,rref_option_t opt,Statistics* stats=nullptr){
    Statistics s;s.input_rows=a.nrow;s.input_cols=a.ncol;s.input_nnz=a.nnz();
    auto start=std::chrono::steady_clock::now();
    auto elapsed=[&](){auto end=std::chrono::steady_clock::now();double t=std::chrono::duration<double>(end-start).count();start=end;return t;};
    // A cheap structural dispatch avoids scanning/sorting all coefficients
    // and lifting a second kernel when very few eliminations are available.
    // This is a cost heuristic, not a mathematical requirement; a missed
    // cascade is still solved exactly by the unchanged backend.
    size_t short_rows=0;for(const auto& row:a.rows)short_rows+=row.nnz()>0&&row.nnz()<=2;
    if(short_rows<std::max(size_t(1),a.ncol/100)){
        s.core_cols=a.ncol;s.core_nnz=a.nnz();s.presolve_s=elapsed();
        auto result=symrep::kernel(std::move(a),opt);s.solve_s=elapsed();if(stats)*stats=s;return result;
    }
    Relations relations(a.ncol);
    // Never require a fixed point for correctness: remaining relations can be
    // solved by the core. Bound scans for long, adversarial dependency chains.
    for(size_t pass=0;pass<8;++pass){
        size_t changed=0;
        for(auto& row:a.rows){if(pass)relations.compress(row);else normalize(row);
            if(row.nnz()<=2){changed+=relations.consume(row);row.clear();}}
        ++s.passes;
        if(!changed&&pass==0){s.core_cols=a.ncol;s.core_nnz=a.nnz();s.presolve_s=elapsed();
            auto result=symrep::kernel(std::move(a),opt);s.solve_s=elapsed();if(stats)*stats=s;return result;}
        if(!changed)break;
    }
    std::vector<I> column(a.ncol,-1);std::vector<Q> coefficients(a.ncol);
    size_t n=0;
    for(size_t i=0;i<a.ncol;++i){auto [root,value]=relations.image(i);coefficients[i]=value;
        if(value!=0&&column[root]<0)column[root]=index(n++);}
    std::vector<I> image(a.ncol,-1);
    for(size_t i=0;i<a.ncol;++i){auto [root,value]=relations.image(i);if(value!=0)image[i]=column[root];}
    for(auto& row:a.rows){relations.compress(row);for(size_t p=0;p<row.nnz();++p)row(p)=column[row(p)];normalize(row);row.compress();}
    a.ncol=n;s.core_cols=n;s.core_nnz=a.nnz();s.presolve_s=elapsed();
    auto core=symrep::kernel(std::move(a),opt);s.solve_s=elapsed();
    auto transposed=core.transpose();core.clear();Mat result(transposed.ncol,s.input_cols);
    std::vector<size_t> counts(result.nrow);
    for(size_t c=0;c<s.input_cols;++c)if(image[c]>=0)for(auto [i,x]:transposed[image[c]])++counts[i];
    for(size_t i=0;i<result.nrow;++i)result[i].reserve(counts[i]);
    for(size_t c=0;c<s.input_cols;++c)if(image[c]>=0)for(auto [i,x]:transposed[image[c]])result[i].push_back(index(c),x*coefficients[c]);
    s.lift_s=elapsed();if(stats)*stats=s;return result;
}
}
