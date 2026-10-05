#include "staged_matrix.hpp"
#include <random>
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==4,"staged_kernel_benchmark assembled|staged small-kernel|large-kernel THREADS");
    bool staged=std::string(argv[1])=="staged",large=std::string(argv[2])=="large-kernel";
    require(staged||std::string(argv[1])=="assembled","unknown benchmark method");
    require(large||std::string(argv[2])=="small-kernel","unknown benchmark case");
    rref_option_t opt;opt->pool.reset(std::stoul(argv[3]));std::mt19937 rng(20261004);
    size_t columns=512,rank=large?80:508;Mat core(rank,columns);
    for(size_t r=0;r<rank;++r){core[r].push_back(index(r),Q(1));for(size_t j=rank;j<std::min(columns,rank+4);++j)core[r].push_back(index(j),Q(long(rng()%9)+1));}
    Mat a=core;
    for(size_t r=0;r<12000;++r){Vec row;for(size_t j=0;j<80;++j){auto& base=core[rng()%rank];Q sign(rng()%2?1:-1);for(auto [c,x]:base)row.push_back(c,sign*x);}normalize(row);a.rows.push_back(std::move(row));}
    a.nrow=a.rows.size();staged_kernel::MatrixRows source(a);auto operators=staged_kernel::stack(source);
    auto start=std::chrono::steady_clock::now();
    auto k=staged?staged_kernel::solve(operators,opt):structured_kernel::solve_stream(a.ncol,source,opt,nullptr,[&](const Mat& candidate,ulong prime,const int_t& modulus){return staged_kernel::verify(candidate,operators,opt,prime,modulus);});
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    require(k.nrow==columns-rank&&row_basis(k).nrow==k.nrow&&!mul(core,transpose(k)).nnz(),"benchmark exact-space mismatch");
    std::cout<<"BENCH method="<<argv[1]<<" case="<<argv[2]<<" columns="<<columns<<" rows="<<a.nrow<<" input_nnz="<<a.nnz()<<" nullity="<<k.nrow<<" kernel_nnz="<<k.nnz()<<" solve_seconds="<<elapsed<<" exact_space=pass"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
