#include "shared/streaming_echelon.hpp"
#include "shared/sparse_affine_factor.hpp"
#include <random>
#include <sstream>
using M=sparse_mat<ulong,int32_t>;
using V=sparse_vec<ulong,int32_t>;
using Dense=std::vector<std::vector<ulong>>;
void check(bool x,const char* msg){if(!x)throw std::runtime_error(msg);}
size_t gauss(Dense& a,size_t cols,const field_t& F,bool back=true){
    size_t r=0;for(size_t c=0;c<cols&&r<a.size();++c){size_t s=r;while(s<a.size()&&!a[s][c])++s;if(s==a.size())continue;std::swap(a[r],a[s]);
        ulong inv=n_invmod(a[r][c],F.mod.n);for(auto& v:a[r])v=nmod_mul(v,inv,F.mod);
        for(size_t i=back?0:r+1;i<a.size();++i)if(i!=r&&a[i][c]){ulong x=a[i][c];for(size_t j=c;j<a[i].size();++j)a[i][j]=nmod_sub(a[i][j],nmod_mul(x,a[r][j],F.mod),F.mod);}++r;
    }return r;
}
M sparse(const Dense& a,size_t cols){M m(a.size(),cols);for(size_t r=0;r<a.size();++r)for(size_t c=0;c<cols;++c)if(a[r][c])m[r].push_back(c,a[r][c]);return m;}
Dense dense(const M& a){Dense d(a.nrow,std::vector<ulong>(a.ncol));for(size_t r=0;r<a.nrow;++r)for(auto[c,x]:a[r])d[r][c]=x;return d;}
Dense reduce(Dense input,size_t n,const field_t& F,rref_option_t opt,structural_replay::Trace* trace){
    streaming_echelon::Reducer reducer(n,F,opt,11,1000,trace);
    for(auto& row:sparse(input,n).rows)reducer.consume(std::move(row));
    auto result=reducer.finish();auto a=dense(result.first);gauss(a,n,F);std::erase_if(a,[](auto& r){return std::none_of(r.begin(),r.end(),[](ulong x){return x;});});return a;
}
int main()try{
    Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(4);opt->method=1;opt->verbose=false;
    std::mt19937 rng(624792);size_t cases=0,replays=0,fallbacks=0,factor_checks=0;
    for(size_t n:{size_t(1),size_t(7),size_t(31),size_t(80)})for(size_t density:{size_t(10),size_t(50),size_t(100)})for(size_t variant=0;variant<4;++variant){
        std::vector<std::vector<int>> integers(2*n+3,std::vector<int>(n+1));
        for(auto& row:integers)for(size_t c=0;c<n+1;++c)if(rng()%100<density)row[c]=int(rng()%19)-9;
        if(variant==0)for(auto& row:integers)for(size_t c=n/2;c<n+1;++c)row[c]=0;
        if(variant==1)for(size_t r=0;r<integers.size();++r)integers[r][n]=integers[r][0];
        if(variant==2)for(auto& row:integers)for(size_t c=0;c<n+1;++c)row[c]=c%3==0?int(c+1):0;
        structural_replay::Trace learned;
        for(ulong p:{1152921504606847009UL,1152921504606847067UL}){
            field_t F(FIELD_Fp,p);Dense a(integers.size(),std::vector<ulong>(n+1));
            for(size_t i=0;i<a.size();++i)for(size_t j=0;j<=n;++j)a[i][j]=integers[i][j]<0?p-ulong(-integers[i][j]):ulong(integers[i][j]);
            Dense reference=a;size_t rank=gauss(reference,n+1,F);reference.resize(rank);
            structural_replay::Trace trace;
            if(p==1152921504606847009UL){auto got=reduce(a,n+1,F,opt,&learned);check(got==reference,"learn row space differs from independent dense Gauss");}
            else{std::stringstream stream;learned.write(stream,"synthetic-family");trace=structural_replay::Trace::read(stream,"synthetic-family");
                try{check(reduce(a,n+1,F,opt,&trace)==reference,"replay row space differs");++replays;}
                catch(const structural_replay::Mismatch&){check(reduce(a,n+1,F,opt,nullptr)==reference,"fallback row space differs");++fallbacks;}}
            Dense triangular=a;size_t variable_rank=gauss(triangular,n,F,false);bool consistent=true;
            for(size_t i=variable_rank;i<triangular.size();++i)consistent&=!triangular[i][n];
            if(consistent){
                std::vector<std::optional<std::vector<ulong>>> rows(n);
                for(size_t i=0;i<variable_rank;++i){size_t c=0;while(c<n&&!triangular[i][c])++c;rows[c]=triangular[i];}
                sparse_affine::Factor factor(rows,F);auto map=factor.materialize();auto values=dense(map);
                Dense expected=a;gauss(expected,n,F);std::vector<int> pivot(n,-1);for(size_t i=0;i<variable_rank;++i){size_t c=0;while(!expected[i][c])++c;pivot[c]=i;}
                Dense expected_map(n-variable_rank+1,std::vector<ulong>(n));size_t q=0;
                for(size_t c=0;c<n;++c)if(pivot[c]<0){expected_map[q][c]=1;for(size_t j=0;j<n;++j)if(pivot[j]>=0)expected_map[q][j]=nmod_neg(expected[pivot[j]][c],F.mod);++q;}
                for(size_t j=0;j<n;++j)if(pivot[j]>=0)expected_map.back()[j]=expected[pivot[j]][n];
                check(values==expected_map,"factored affine map differs from dense Gauss");
                std::vector<ulong> params(factor.free.size());for(auto& x:params)x=rng()%p;
                auto point=factor.lift(params);for(auto& row:a){ulong z=0;for(size_t j=0;j<n;++j)z=nmod_add(z,nmod_mul(row[j],point[j],F.mod),F.mod);check(z==row[n],"factored lift violates original equation");}
                std::vector<ulong> op(n+1);for(auto& x:op)x=rng()%p;auto projected=factor.project(op);ulong z=projected.back(),v=op.back();
                for(size_t j=0;j<params.size();++j)z=nmod_add(z,nmod_mul(params[j],projected[j],F.mod),F.mod);
                for(size_t j=0;j<n;++j)v=nmod_add(v,nmod_mul(op[j],point[j],F.mod),F.mod);check(z==v,"factored operator projection mismatch");++factor_checks;
            }++cases;
        }
        std::stringstream bad;learned.write(bad,"a");bool caught=false;try{auto t=structural_replay::Trace::read(bad,"b");}catch(const structural_replay::Mismatch&){caught=true;}check(caught,"identity guard missing");
    }
    // Same support, exceptional prime loses rank, then another field gains a
    // previously discarded equation. Both must fall back without publishing.
    for(bool gain:{false,true}){
        structural_replay::Trace t;field_t F(FIELD_Fp,gain?101:103),G(FIELD_Fp,gain?103:101);Dense a={{1,1,0},{1,102%F.mod.n,0}};reduce(a,3,F,opt,&t);
        std::stringstream s;t.write(s,"exception");auto replay=structural_replay::Trace::read(s,"exception");Dense b={{1,1,0},{1,102%G.mod.n,0}};
        bool caught=false;try{reduce(b,3,G,opt,&replay);}catch(const structural_replay::Mismatch&){caught=true;}
        check(caught,"exceptional rank guard missing");++fallbacks;
    }
    {
        field_t F(FIELD_Fp,101);std::vector<std::optional<std::vector<ulong>>> rows(2);
        rows[0]=std::vector<ulong>{1,2,3};sparse_affine::Factor factor(rows,F);size_t guards=0;
        try{factor.lift({});}catch(const std::runtime_error&){++guards;}
        try{factor.project({1});}catch(const std::runtime_error&){++guards;}
        rows[0]=std::vector<ulong>{2,2,3};try{sparse_affine::Factor bad(rows,F);}catch(const std::runtime_error&){++guards;}
        rows[0].reset();rows[1]=std::vector<ulong>{1,1,3};try{sparse_affine::Factor bad(rows,F);}catch(const std::runtime_error&){++guards;}
        check(guards==4,"affine map shape/triangular guards missing");
    }
    std::cout<<"GENERIC_CONTROL passed=1 field_cases="<<cases<<" successful_replays="<<replays<<" guarded_fallbacks="<<fallbacks<<" factor_cases="<<factor_checks<<std::endl;
}catch(const std::exception& e){std::cerr<<"CONTROL_ERROR "<<e.what()<<std::endl;return 1;}
