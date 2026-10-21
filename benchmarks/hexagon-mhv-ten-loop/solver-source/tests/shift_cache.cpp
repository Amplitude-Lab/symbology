#include "irr_shift_cache.hpp"
namespace fs=std::filesystem;
using irr::require;
template<class T> sparse_mat<T,int32_t> pack(const std::array<sparse_mat<T,int32_t>,6>& blocks){
	size_t cols=0;for(const auto& b:blocks)cols+=b.ncol;sparse_mat<T,int32_t> out(blocks[0].nrow,cols);size_t off=0;
	for(const auto& b:blocks){for(size_t i=0;i<b.nrow;++i)for(auto [j,x]:b[i])out[i].push_back(off+j,x);off+=b.ncol;}return out;
}
int main(int argc,char** argv){try{
	require(argc==5,"test_shift_cache FEC LEC UNIQUE_L6 OUTPUT");fs::path fec=argv[1],lec=argv[2],source=argv[3],root=argv[4];
	require(!fs::exists(root),"Use a fresh test directory");fs::create_directories(root);
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=false;opt->pool.reset(4);
	irr::channels_t cg(&opt->pool);irr::seam_recoupling_t recouple(cg,opt);
	// Unit tests of private-column and general modular charts, including bad primes.
	field_t F(FIELD_Fp,101);sparse_mat<ulong,int32_t> chart(2,3);
	chart[0].push_back(0,1);chart[0].push_back(1,1);chart[1].push_back(0,1);chart[1].push_back(1,2);
	irr_shift_cache::right_inverse(chart,F,opt);chart[1]=chart[0];bool caught=false;
	try{irr_shift_cache::right_inverse(chart,F,opt);}catch(const irr_shift_cache::bad_prime&){caught=true;}require(caught,"Rank loss not detected");
	caught=false;try{irr_shift_cache::residue(rat_t(1)/rat_t(101),F);}catch(const irr_shift_cache::bad_prime&){caught=true;}require(caught,"Bad denominator not detected");
	auto unique=irr_io::read_object(source,10,2,fec,lec,&opt->pool);
	for(std::string kind:{"unique","homogeneous","affine"}){
		// Synthetic coordinate-order fixtures from a genuine L6 object. These
		// rows test linearity/origin placement, not independent physical solutions.
		auto obj=unique;obj.kind=kind;obj.directions=kind=="unique"?0:2;
		for(size_t b=0;b<6;++b){obj.S[b]=sparse_mat<rat_t,int32_t>(obj.rows(),unique.S[b].ncol);
			for(size_t i=0;i<obj.rows();++i)for(auto [j,x]:unique.S[b][0])obj.S[b][i].push_back(j,x*rat_t(i+1));}
		auto expected=irr::shift_one(obj,true,fec,lec,cg,recouple,opt);auto expected2=irr::shift_one(expected,true,fec,lec,cg,recouple,opt);
		auto dir=root/kind;fs::create_directories(dir/"cache");irr_io::write_reduced(dir,obj,&opt->pool);
		irr_stage::write(dir/"basis_identity.tsv",irr_stage::basis_identity(fec,lec,10,2));irr_stage::fixed_text(dir/"input_identity.txt","synthetic L6 shift fixture "+kind+"\n");
		fs::create_directories(root/(kind+"_expected"));irr_io::write_reduced(root/(kind+"_expected"),expected2,&opt->pool);
		ulong prime=n_nextprime(1ULL<<60,0);auto matrix=pack(obj.S);
		for(size_t n=0;n<4;++n){field_t field(FIELD_Fp,prime);auto mod=irr_shift_cache::reduce(matrix,field);
			auto shifted=irr_shift_cache::shift_one(mod,10,2,fec,lec,cg,recouple,field,opt);
			require(nmhv_sew_matrices_equal(shifted,irr_shift_cache::reduce(pack(expected.S),field)),"One-step QQ/Fp disagreement");
			shifted=irr_shift_cache::shift_one(shifted,9,3,fec,lec,cg,recouple,field,opt);
			require(nmhv_sew_matrices_equal(shifted,irr_shift_cache::reduce(pack(expected2.S),field)),"Two-step QQ/Fp disagreement");
			nmhv_sew_write_modular_basis(dir/"cache"/("p"+std::to_string(prime)+".bin"),mod,prime,6,10,2,irr_shift_cache::fingerprint(dir));prime=n_nextprime(prime,0);
		}std::cout<<"PASS QQ/Fp L6 two shifts, four primes, kind="<<kind<<std::endl;
	}
	std::cout<<"PASS modular chart and denominator guards"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
