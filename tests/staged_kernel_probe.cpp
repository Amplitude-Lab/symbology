#include "staged_matrix.hpp"
#include "product_invariance.hpp"
#include "symrep_bootstrap.hpp"
#include <random>
using namespace symrep;

void check(const Mat& a, const Mat& k) {
    require(k.ncol==a.ncol && !mul(a,transpose(k)).nnz(),"staged kernel exact residual");
    require(row_basis(k).nrow==k.nrow && k.nrow+row_basis(a).nrow==a.ncol,"staged kernel completeness");
}
int main(int argc,char** argv){try{
    rref_option_t opt;opt->pool.reset(3);std::mt19937 rng(20261004);size_t cases=0;
    for(size_t trial=0;trial<100;++trial){
        size_t n=trial%19,m=trial%23;Mat a(m,n);
        for(size_t r=0;r<m;++r)for(size_t c=0;c<n;++c)if(rng()%3==0){Q x(long(rng()%13)-6,long(1+rng()%5));if(x!=0)a[r].push_back(index(c),x);}
        staged_kernel::Options options;options.seed_terms=trial%2?1:4;options.target_nullity=trial%3?2:64;
        options.seed_nonzeros=trial%4?8:1024;options.batch_rows=3;options.batch_nonzeros=16;
        options.complete_through_terms=trial%2?4:0;options.remember_equations=trial%3;options.sample_equations=trial%2;
        auto k=staged_kernel::solve(a,opt,options);check(a,k);++cases;
        // Composition includes a redundant copy; no application-specific solver.
        Mat b=a;for(auto& row:b.rows)for(size_t j=0;j<row.nnz();++j)row[j]*=Q(3,7);
        staged_kernel::MatrixRows first(a),second(b);
        auto joined=staged_kernel::solve(staged_kernel::stack(first,second),opt,options);
        check(a,joined);++cases;
        std::vector<Q> rhs(m);for(size_t r=0;r<m;++r)for(auto [c,x]:a[r])rhs[r]+=x*Q(long(c)+1);
        auto affine=staged_kernel::solve_affine(a,rhs,opt,options);
        require(affine.consistent&&affine.particular.nrow==1,"consistent affine system rejected");
        auto image=mul(a,transpose(affine.particular));
        for(size_t r=0;r<m;++r){Q value=0;if(auto p=image[r].find(0))value=*p;require(value==rhs[r],"affine particular solution mismatch");}
        check(a,affine.directions);++cases;
    }
    // A long dependent plateau precedes a necessary late equation. Early
    // switching must still intersect that equation, never accept the sample.
    Mat plateau(45,12);
    for(size_t r=0;r<44;++r){plateau[r].push_back(0,Q(1));plateau[r].push_back(1,Q(2));plateau[r].push_back(2,Q(3));}
    for(size_t c=3;c<12;++c)plateau[44].push_back(index(c),Q(long(c)));
    staged_kernel::Options adaptive;adaptive.seed_terms=1;adaptive.target_nullity=1;adaptive.batch_rows=1;
    adaptive.stagnation_batches=2;adaptive.stagnation_rank_gain=0;adaptive.max_residual_nullity=32;
    auto pk=staged_kernel::solve(plateau,opt,adaptive);check(plateau,pk);++cases;
    // Contradictions, large nullspaces and difficult primes.
    Mat inconsistent(2,1);inconsistent[0].push_back(0,Q(1));inconsistent[1].push_back(0,Q(1));
    require(!staged_kernel::solve_affine(inconsistent,{Q(1),Q(2)},opt).consistent,"inconsistent affine system accepted");
    Mat wide(3,200);wide[0].push_back(0,Q(1));wide[0].push_back(1,Q(2));wide[0].push_back(2,Q(3));
    auto wk=staged_kernel::solve(wide,opt);check(wide,wk);++cases;
    ulong prime=n_nextprime(1ULL<<60,0);Mat hard(2,5);
    hard[0].push_back(0,Q(1));hard[0].push_back(1,Q(int_t(2).pow(90ul)+int_t(17)));hard[0].push_back(2,Q(3));
    hard[1].push_back(1,Q(1)/Q(prime));hard[1].push_back(3,Q(2));hard[1].push_back(4,Q(1));
    auto hk=staged_kernel::solve(hard,opt);check(hard,hk);++cases;

    // All cuts of weight six for a three-letter commuting alphabet. Expand
    // independently into original words, so agreement is basis-independent.
    const size_t letters=3,weight=6;Mat rules(3,9);size_t q=0;
    for(size_t i=0;i<letters;++i)for(size_t j=i+1;j<letters;++j){rules[q].push_back(index(i*letters+j),Q(1));rules[q++].push_back(index(j*letters+i),Q(-1));}
    auto d=condition_tensor(rules,letters);
    std::vector<Tensor> f(weight),l(weight);std::vector<Mat> fw(weight),lw(weight);
    f[1]=tensor(identity(letters),{letters,1,letters});l[1]=tensor(identity(letters),{letters,letters,1});fw[1]=lw[1]=identity(letters);
    for(size_t w=2;w<weight;++w){
        f[w]=tensor_kernel::extend_tensor(f[w-1],d,false,opt,tensor_kernel::Strategy::staged);
        l[w]=tensor_kernel::extend_tensor(l[w-1],d,true,opt,tensor_kernel::Strategy::staged);
        fw[w]=mul(flatten(f[w]),kron(fw[w-1],identity(letters)));
        lw[w]=mul(flatten(l[w]),kron(identity(letters),lw[w-1]));
    }
    Mat reference;
    for(size_t right=1;right<weight;++right){size_t left=weight-right;
        auto dc=d,fc=f[left],lc=l[right];auto k=tensor_kernel::sew_one(std::move(dc),std::move(fc),std::move(lc),opt,tensor_kernel::Strategy::staged,{},right%2?tensor_kernel::BoundaryReduction::raw:tensor_kernel::BoundaryReduction::reduced);
        auto words=row_basis(mul(flatten(k),kron(fw[left],lw[right])));
        require(words.nrow==28,"arbitrary-cut symmetric tensor dimension");
        if(reference.nrow)require(equal(reference,words),"different sewing cuts changed the complete word space");else reference=words;
    }
    if(argc==2){
        fs::path out=argv[1];fs::create_directories(out);
        write(out/"condition.wxf",d);write(out/"FEC_2.wxf",f[2]);write(out/"LEC_2.wxf",l[2]);
        Mat cycle(3,3);for(size_t i=0;i<3;++i)cycle[i].push_back(index((i+1)%3),Q(1));
        auto fa=Chart(fw[2],opt).coordinates(mul(fw[2],kron(cycle,cycle)));
        auto la=Chart(lw[2],opt).coordinates(mul(lw[2],kron(cycle,cycle)));
        write(out/"left_action.wxf",fa);write(out/"right_action.wxf",la);
        auto dc=d,fc=f[2],lc=l[2];auto sewn=flatten(tensor_kernel::sew_one(std::move(dc),std::move(fc),std::move(lc),opt));
        auto transformed=mul(sewn,kron(fa,la));
        for(size_t r=0;r<sewn.nrow;++r){for(auto [c,x]:sewn[r])transformed[r].push_back(c,-x);normalize(transformed[r]);}
        auto invariant=mul(kernel(transpose(transformed),opt),sewn);
        require(invariant.nrow==5,"cyclic invariant dimension");
        write(out/"reference.wxf",tensor(invariant,{invariant.nrow,f[2].dim(0),l[2].dim(0)}));
    }
    std::cout<<"PASS "<<cases<<" exact generic kernel/composition/affine cases, inconsistency, bad primes and all five weight-six cuts"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
