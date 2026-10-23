#include "boundary_tensor_train.hpp"
using namespace boundary_tt;
int main(int argc,char**argv){try{require(argc==2,"usage boundary_H6_train TRAIN_DIRECTORY");fs::path p=argv[1];rref_option_t opt;opt->pool.reset(1);auto e=load(p/"factor_0"),r=load(p/"factor_1");auto e2=shuffle(e,e,opt),e3=shuffle(e2,e,opt),er=shuffle(e,r,opt);auto h=sum({{std::move(er),Q(1)},{std::move(e3),Q(1,6)}},opt);save(p/"H6",h);std::cout<<"H6 ranks=";for(auto d:h.rank)std::cout<<d<<',';std::cout<<std::endl;}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
