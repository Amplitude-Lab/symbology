#include "reduced_d3/irr_recurrence.hpp"
#include <iostream>
int main(int argc,char** argv){try{
	if(argc!=3)throw std::runtime_error("Expected two tensor paths");
	Flint::set_memory_functions();
	auto a=irr::read_tensor(argv[1],nullptr),b=irr::read_tensor(argv[2],nullptr);
	if(a.dims()!=b.dims()||a.nnz()!=b.nnz())throw std::runtime_error("Tensor shape/nnz differs");
	for(size_t i=0;i<a.dim(0)+1;++i)if(a.rowptr()[i]!=b.rowptr()[i])throw std::runtime_error("Row pointers differ");
	for(size_t p=0;p<a.nnz();++p){
		if(a.val(p)!=b.val(p))throw std::runtime_error("Tensor value differs");
		for(size_t j=0;j+1<a.rank();++j)if(a.index(p)[j]!=b.index(p)[j])throw std::runtime_error("Tensor index differs");
	}std::cout<<"PASS exact tensor "<<argv[2]<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
