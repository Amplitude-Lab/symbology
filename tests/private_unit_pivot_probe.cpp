#include "certified_kernel.hpp"
using namespace symrep;
int main(){try{
    using MM=certified_kernel::ModMat;std::vector<std::vector<pivot_t<I>>> pivots;
    MM a(2,5);a[0].push_back(0,ulong(7));a[0].push_back(2,ulong(1));a[1].push_back(0,ulong(5));a[1].push_back(4,ulong(1));
    require(certified_kernel::private_unit_pivots(a,pivots),"private units not recognized");
    auto [row0,column0]=pivots[0][0];auto [row1,column1]=pivots[0][1];require(column0==2&&column1==4,"incorrect private pivots");
    a[1].push_back(2,ulong(1));a[1].sort_indices();require(!certified_kernel::private_unit_pivots(a,pivots),"shared column accepted as private");
    MM empty(0,4);require(certified_kernel::private_unit_pivots(empty,pivots)&&pivots[0].empty(),"empty identity block");
    // A large pivot ratio can multiply numerator heights. Keeping the known
    // unit block reconstructs this basis within three 61-bit primes.
    Q m(int_t(2).pow(70ul)+int_t(3)),n(int_t(2).pow(71ul)+int_t(7)),d(int_t(2).pow(69ul)+int_t(11));
    Mat equation(1,3);equation[0].push_back(0,Q(1)/m);equation[0].push_back(1,n/d);equation[0].push_back(2,Q(1));
    size_t passes=0;rref_option_t opt;opt->pool.reset(2);
    auto generate=[&](auto emit){++passes;emit(Vec(equation[0]));};
    auto verify=[&](const Mat& k){return !mul(equation,transpose(k)).nnz();};
    auto k=certified_kernel::solve(3,generate,verify,opt);
    require(k.nrow==2&&row_basis(k).nrow==2&&verify(k)&&passes<=3,"unit pivot reconstruction regression");
    require(k[0].find(2)&&*k[0].find(2)==Q(1)/m&&k[1].find(2)&&*k[1].find(2)==n/d,"unit coordinate chart was changed");
    size_t repivoted_passes=0;
    auto repivot=[&](const field_t& fp,auto emit){++repivoted_passes;MM reduced(1,3);reduced[0]=equation[0]%fp.mod;
        sparse_mat_rref(reduced,fp,opt);emit(std::move(reduced[0]));};
    auto repivoted=certified_kernel::solve(3,generate,verify,opt,repivot);
    require(repivoted.nrow==2&&row_basis(repivoted).nrow==2&&verify(repivoted)&&repivoted_passes>passes,"large-height control did not exercise avoidable repivoting");
    std::cout<<"PASS private-unit rank certificate, shared-column rejection, empty matrix and large-height reconstruction; primes="<<passes<<" repivoted_primes="<<repivoted_passes<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
