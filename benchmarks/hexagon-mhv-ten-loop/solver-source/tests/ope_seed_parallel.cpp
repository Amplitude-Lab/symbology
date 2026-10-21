#include "irr_projection.hpp"
#include <atomic>

struct test_germ_t {
	ulong add(ulong a,ulong b)const{return (a+b)%1000003;}
	ulong mul(ulong a,ulong b)const{return (a*b)%1000003;}
};

int main(){try {
	constexpr size_t stride=4096,first=3,last=2;
	std::vector<ulong> head(first*stride),expected(last*stride),actual(last*stride);
	for(size_t i=0;i<head.size();++i)head[i]=(i*113+17)%1000003;
	const std::vector<std::tuple<size_t,size_t,ulong>> terms={
		{0,0,3},{1,1,1000002},{2,0,17},{0,0,29},{1,0,0},{2,1,41}
	};
	test_germ_t germ;
	for(auto [i,j,v]:terms)for(size_t z=0;z<stride;++z)
		expected[j*stride+z]=germ.add(expected[j*stride+z],germ.mul(v,head[i*stride+z]));
	Flint::set_memory_functions();rref_option_t opt;
	for(size_t threads:{size_t(1),size_t(4)}){
		opt->pool.reset(threads);std::fill(actual.begin(),actual.end(),0);
		std::atomic<size_t> active=0,maximum=0;
		auto walk=[&](auto emit){
			auto count=++active,peak=maximum.load();
			while(peak<count&&!maximum.compare_exchange_weak(peak,count)){}
			if(threads>1)std::this_thread::sleep_for(std::chrono::milliseconds(10));
			for(auto [i,j,v]:terms)emit(i,j,v);
			--active;
		};
		irr_physics::accumulate_ope_seed(actual,head,stride,germ,&opt->pool,walk);
		irr::require(actual==expected,"OPE seed changed finite-field coefficients");
		if(threads>1)irr::require(maximum>1,"OPE seed did not use multiple workers");
	}
	std::cout<<"OPE SEED PARALLEL PASS"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
