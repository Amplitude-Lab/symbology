#include "structured_kernel.hpp"
#include <random>
using namespace symrep;
int main(){try{
    rref_option_t opt;opt->pool.reset(2);std::mt19937 rng(73425);size_t cases=0;
    auto check=[&](Mat a,bool test_legacy_backend=true){
        size_t nullity=a.ncol-row_basis(a).nrow;
        if(test_legacy_backend){auto k=structured_kernel::solve(a,opt);
            require(!mul(a,k.transpose()).nnz(),"presolve exact residual failed");
            require(row_basis(k).nrow==k.nrow,"presolve dependent output");
            require(k.nrow==nullity,"presolve incomplete kernel");}
        auto generate=[&](auto emit){for(const auto& row:a.rows)emit(Vec(row));};
        auto modular=certified_kernel::solve(a.ncol,generate,[&](const Mat& candidate){return !mul(a,transpose(candidate)).nnz();},opt);
        require(modular.nrow==nullity,"certified modular rank mismatch");++cases;
        auto batched=certified_kernel::solve(a.ncol,generate,[&](const Mat& candidate){return !mul(a,transpose(candidate)).nnz();},opt,nullptr,4);
        require(batched.nrow==nullity&&row_basis(batched).nrow==nullity,"batched modular rank mismatch");
        auto stream=[&](structured_kernel::Relations& relations,size_t limit,auto emit){for(const auto& original:a.rows){Vec row=original;relations.compress(row);if(row.nnz()<=limit)emit(std::move(row));}};
        auto streamed=structured_kernel::solve_stream(a.ncol,stream,opt,nullptr,[&](const Mat& candidate){return !mul(a,transpose(candidate)).nnz();});
        require(streamed.nrow==nullity,"streamed kernel incomplete");
    };
    check(Mat(0,0));check(Mat(0,7));check(Mat(8,0));check(Mat(8,9));
    for(size_t c=0;c<150;++c){size_t nr=rng()%65,nc=rng()%40;Mat a(nr,nc);
        for(size_t i=0;i<nr;++i){size_t terms=c%3==0?2:1+rng()%8;
            for(size_t j=0;j<terms&&nc;++j){long num=long(rng()%13)-6;Q x=Q(num)/Q(1+rng()%5);if(x!=0)a[i].push_back(index(rng()%nc),x);}normalize(a[i]);}
        check(std::move(a));
    }
    // Non-unit scales, a contradictory homogeneous cycle, and variables
    // absent from every equation. Cycle contradiction means zeros, not error.
    Mat a(4,6);a[0].push_back(0,Q(2));a[0].push_back(1,Q(-3));
    a[1].push_back(1,Q(5));a[1].push_back(2,Q(-7));
    a[2].push_back(0,Q(1));a[2].push_back(2,Q(1));
    a[3].push_back(2,Q(1));a[3].push_back(3,Q(2));a[3].push_back(4,Q(3));check(a);
    ulong p=n_nextprime(1ULL<<60,0);Mat bad(2,2);bad[0].push_back(0,Q(1));bad[0].push_back(1,Q(1));bad[1].push_back(0,Q(1));bad[1].push_back(1,Q(p)+1);check(bad);
    Mat missing(1,2);missing[0].push_back(0,Q(1));missing[0].push_back(1,Q(p));check(missing);
    Mat denominator(1,3);denominator[0].push_back(0,Q(1)/Q(p));denominator[0].push_back(1,Q(1));denominator[0].push_back(2,Q(2));check(denominator,false);
    Mat large(1,2);large[0].push_back(0,Q(1));large[0].push_back(1,Q(int_t(2).pow(180ul)+int_t(37)));check(large);
    Mat mixed(3,6);for(size_t r=0;r<3;++r)mixed[r].push_back(index(r),Q(1));
    mixed[0].push_back(3,Q(1)/Q(2));mixed[1].push_back(4,Q(int_t(2).pow(80ul)+int_t(23)));mixed[2].push_back(5,Q(int_t(2).pow(180ul)+int_t(37)));check(mixed);
    Mat small(1,3);small[0].push_back(0,Q(1));small[0].push_back(1,Q(1)/Q(2));small[0].push_back(2,Q(3)/Q(7));size_t rejected=0;
    auto retry=certified_kernel::solve(small.ncol,[&](auto emit){emit(Vec(small[0]));},[&](const Mat& k){require(!mul(small,transpose(k)).nnz(),"decoded rational CRT rows changed their kernel");return ++rejected>=3;},opt);
    require(rejected==3&&retry.nrow==2,"rational CRT rows lost state after failed certification");
    std::cout<<"PASS "<<cases<<" independent exact rank/residual/kernel-rank cases, including unlucky primes, vanishing supports, prime denominators and multi-prime CRT"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
