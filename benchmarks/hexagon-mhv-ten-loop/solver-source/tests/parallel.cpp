#include "irr_conversion.hpp"

int main(int argc,char** argv){try {
	if(argc!=2)throw std::runtime_error("test_parallel requires a fresh cache directory");
	std::filesystem::path root=argv[1];
	if(std::filesystem::exists(root))throw std::runtime_error("Refusing existing test directory");
	std::filesystem::create_directories(root);
	Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(4);
	bool caught=false;try{irr_parallel::rows(100,&opt->pool,[](size_t r){if(r==17)throw std::runtime_error("worker");});}catch(const std::runtime_error&){caught=true;}
	irr::require(caught,"Detached worker exception lost");
	sparse_mat<rat_t,int32_t> expected(137,13);
	int_t large("1267650600228229401496703205376");
	ulong prime=n_nextprime(1ULL<<60,0);
	for(size_t r=0;r<expected.nrow;++r)for(size_t c=0;c<expected.ncol;++c)if(r%11&&(r+c)%3)
		expected[r].push_back(int32_t(c),c==1?rat_t(prime):rat_t(large+int_t(r+c))/rat_t(c+1));
	k_crt::Pivots pivots;
	std::vector<std::pair<ulong,std::filesystem::path>> files;
	for(size_t p=0;p<5;++p,prime=n_nextprime(prime,0)){
		field_t F(FIELD_Fp,prime);sparse_mat<ulong,int32_t> modular(expected.nrow,expected.ncol);
		for(size_t r=0;r<expected.nrow;++r){modular[r]=expected[r]%F.mod;modular[r].canonicalize();}
		auto path=root/("prime_"+std::to_string(p)+".bin");k_crt::save(path,modular,prime,pivots);files.push_back({prime,path});
	}
	sparse_mat<rat_t,int32_t> actual;
	auto input_height=sparse_mat_height(expected);
	auto short_files=files;short_files.resize(1);
	irr::require(!k_crt::reconstruct(short_files,actual,137,13,pivots,input_height,&opt->pool),"Insufficient modulus accepted");
	for(size_t threads:{1,4}){
		opt->pool.reset(threads);
		irr::require(k_crt::reconstruct(files,actual,137,13,pivots,input_height,&opt->pool),"CRT failed");
		irr::require(nmhv_sew_matrices_equal(actual,expected),"CRT changed exact entries");
	}
	{std::ofstream append(files.back().second,std::ios::binary|std::ios::app);append.put('x');}
	caught=false;try{k_crt::reconstruct(files,actual,137,13,pivots,input_height,&opt->pool);}catch(const std::runtime_error&){caught=true;}
	irr::require(caught,"Corrupt cache accepted");
	std::cout<<"PARALLEL PASS: exceptions, changing modular support, insufficient modulus, exact CRT, corruption"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
