#include "tensor_kernel.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==6||argc==7,"early_sew_check DATA CHAIN OLD_SEW NEW_SEW LEFT_WEIGHT [NEXT_RIGHT_BASIS]");rref_option_t opt;opt->pool.reset(4);
    fs::path data=argv[1],chain=argv[2];size_t weight=std::stoul(argv[5]);auto old=read(argv[3]),fresh=read(argv[4]);
    auto f=read(chain/("FEC_"+std::to_string(weight+1)+".wxf"));auto last=read(data/"LEC_1.wxf"),d=read(data/"dlogmat_E6.wxf");
    auto last2=argc==7?read(argv[6]):tensor_kernel::extend_tensor(last,d,true,opt);size_t previous_right=last2.dim(2);auto right=flatten(last2);
    require(fresh.dim(1)==f.dim(1)&&fresh.dim(2)==last2.dim(0)&&old.dim(1)==f.dim(0)&&old.dim(2)==previous_right,"sewing comparison dimensions");
    Mat newspace(fresh.dim(0),f.dim(1)*right.ncol),oldspace(old.dim(0),newspace.ncol);
    for(size_t b=0;b<fresh.dim(0);++b)for(size_t p=fresh.rowptr()[b];p<fresh.rowptr()[b+1];++p){auto i=fresh.index(p);for(auto [j,x]:right[i[1]])newspace[b].push_back(index(size_t(i[0])*right.ncol+j),fresh.val(p)*x);}
    for(size_t b=0;b<old.dim(0);++b)for(size_t p=old.rowptr()[b];p<old.rowptr()[b+1];++p){auto i=old.index(p);for(size_t q=f.rowptr()[i[0]];q<f.rowptr()[i[0]+1];++q){auto j=f.index(q);oldspace[b].push_back(index((size_t(j[0])*d.dim(0)+j[1])*previous_right+i[1]),old.val(p)*f.val(q));}}
    for(auto& row:newspace.rows)normalize(row);for(auto& row:oldspace.rows)normalize(row);
    require(equal(row_basis(newspace),row_basis(oldspace)),"early sewing differs from independently sewn invariant space");
    std::cout<<"PASS exact complete invariant spaces across "<<weight+1<<"+m and "<<weight<<"+(m+1) splits; dimension="<<fresh.dim(0)<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
