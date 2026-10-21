#include "irr_physics_pipeline.hpp"
#include <random>
#include <ctime>
int main()try{
    Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(16);
    std::mt19937 rng(716288);field_t F(FIELD_Fp,1152921504606847009UL);
    for(auto [n,rank,degree]:std::vector<std::tuple<size_t,size_t,size_t>>{{800,720,3},{800,400,8},{500,400,500},{800,40,800},{800,800,3}}){
        affine_rows source(n,F,&opt->pool);
        for(size_t j=0;j<rank;++j){std::vector<ulong> row(n+1);row[j]=1;row[n]=rng()%101;
            if(degree>=n){for(size_t k=j+1;k<n;++k)row[k]=rng()%11;}
            else if(j+1<n)for(size_t z=0;z<degree;++z)row[j+1+rng()%(n-j-1)]=1+rng()%10;
            source.pivots[j]=std::move(row);
        }
        auto old=source,newer=source;setenv("HEXAGON_AFFINE_FACTORS","0",1);auto start=std::clock();auto a=old.solve();double legacy=double(std::clock()-start)/CLOCKS_PER_SEC;
        setenv("HEXAGON_AFFINE_FACTORS","1",1);start=std::clock();auto b=newer.solve();double factored=double(std::clock()-start)/CLOCKS_PER_SEC;
        require(a.origin==b.origin&&a.pivot_columns==b.pivot_columns&&a.directions.nrow==b.directions.nrow,"Synthetic affine result shape/origin mismatch");
        for(size_t r=0;r<a.directions.nrow;++r){const auto& x=a.directions[r];const auto& y=b.directions[r];require(x.nnz()==y.nnz(),"Synthetic affine support mismatch");for(size_t i=0;i<x.nnz();++i)require(x(i)==y(i)&&x[i]==y[i],"Synthetic affine coefficient mismatch");}
        // Check the original triangular equations, including their affine RHS.
        for(size_t j=0;j<rank;++j){ulong x=0;for(size_t k=j;k<n;++k)x=nmod_add(x,nmod_mul((*source.pivots[j])[k],b.origin[k],F.mod),F.mod);require(x==source.pivots[j]->back(),"Synthetic affine origin residual");
            for(const auto& row:b.directions.rows){x=0;for(auto[k,v]:row)x=nmod_add(x,nmod_mul((*source.pivots[j])[k],v,F.mod),F.mod);require(!x,"Synthetic affine direction residual");}}
        std::cout<<"AFFINE_BENCH n="<<n<<" rank="<<rank<<" degree="<<degree<<" legacy_cpu="<<legacy<<" factored_cpu="<<factored<<" coefficient_equality=1 original_equations=1"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<"AFFINE_BENCH_ERROR "<<e.what()<<std::endl;return 1;}
