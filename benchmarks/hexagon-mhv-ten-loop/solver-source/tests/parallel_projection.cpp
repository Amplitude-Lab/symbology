#include "irr_projection.hpp"

// Emit exact finite-field images for old/new implementation comparison. The
// optional solution cap bounds OPE test cost; it does not alter production APIs.
int main(int argc,char** argv){try {
	irr::require(argc==10,"projection HELICITY H K FEC LEC SEW OUTPUT collinear|mrk|ope ROWS");
	Flint::set_memory_functions();std::filesystem::path output=argv[7];
	irr::require(!std::filesystem::exists(output),"Refusing existing projection output");
	irr_physics::context c(argv[4],argv[5],output,"data",std::stoul(argv[2]),std::stoul(argv[3]),16,std::string(argv[1])=="MHV",1000003);
	c.load(argv[6]);c.rows=c.sew.directions;size_t cap=std::stoul(argv[9]);if(cap)c.rows=std::min(c.rows,cap);
	for(size_t b=0;b<6;++b)c.S[b]=nmhv_sew_modular_reduce(c.sew.S[b].take({0,c.rows}),c.F,&c.opt->pool);
	std::ofstream result(output/"image.tsv");std::string mode=argv[8];
	if(mode=="ope"){
		for(size_t f:{2,3}){mhv_ope_germ_t germ(c.h+c.k,c.F.mod.n,c.data/"ope"/("kernels_F"+std::to_string(f)+".tsv"));
			auto values=irr_physics::ope_values(c,germ);for(auto x:values)result<<x<<'\n';
		}
	}else {
		auto rows=irr_physics::slots(c,mode=="mrk");
		for(const auto& row:rows){result<<row.rhs;for(auto x:row.coefficients)result<<'\t'<<x;result<<'\n';}
	}
	std::cout<<"PROJECTION PASS"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
