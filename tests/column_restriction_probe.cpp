#include "staged_matrix.hpp"
#include "column_restriction.hpp"
#include <random>
using namespace symrep;
int main(){try{
    std::mt19937 rng(20261004);rref_option_t opt;opt->pool.reset(3);
    for(size_t trial=0;trial<120;++trial){
        size_t n=trial%21,m=trial%27;Mat a(m,n);
        for(size_t r=0;r<m;++r)for(size_t c=0;c<n;++c)if(rng()%3==0){Q x(long(rng()%11)-5,long(1+rng()%5));if(x!=0)a[r].push_back(index(c),x);}
        std::vector<I> select;for(size_t c=0;c<n;++c)if(rng()%3)select.push_back(index(c));std::shuffle(select.begin(),select.end(),rng);
        Mat reduced(m,select.size());for(size_t r=0;r<m;++r)for(size_t c=0;c<select.size();++c)if(auto x=a[r].find(select[c]))reduced[r].push_back(index(c),*x);
        staged_kernel::MatrixRows rows(a);staged_kernel::ColumnRestriction restricted(rows,select);
        staged_kernel::Options options;options.seed_terms=trial%2?1:8;options.target_nullity=trial%3?1:32;options.seed_nonzeros=16;options.batch_rows=3;options.batch_nonzeros=32;
        options.complete_through_terms=trial%2?4:0;options.sample_equations=trial%3;
        auto k=staged_kernel::solve(staged_kernel::stack(restricted),opt,options);
        require(!mul(reduced,transpose(k)).nnz(),"restricted residual");
        require(row_basis(k).nrow==k.nrow&&row_basis(reduced).nrow+k.nrow==select.size(),"restricted completeness");
        restricted.lift_in_place(k);require(k.ncol==a.ncol&&!mul(a,transpose(k)).nnz(),"restricted lifted residual");
    }
    std::cout<<"PASS 120 exact coordinate-restriction kernels, reordered/empty selections, short relations and staged paths"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
