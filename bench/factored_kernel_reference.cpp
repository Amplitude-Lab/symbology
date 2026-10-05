#include "tensor_kernel.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==7,"factored_kernel_reference forward|backward CONDITION SEED MAX OUTPUT THREADS");
    bool backward=std::string(argv[1])=="backward";auto condition=read(argv[2]),seed=read(argv[3]);
    size_t maximum=std::stoul(argv[4]);fs::path output=argv[5];fresh(output);
    rref_option_t opt;opt->pool.reset(std::stoul(argv[6]));write(output/"w1.wxf",seed);
    std::ofstream stats(output/"timings.tsv");stats<<"weight\tdimension\tsolve_s\n";
    for(size_t w=2;w<=maximum;++w){auto start=std::chrono::steady_clock::now();
        seed=tensor_kernel::extend_tensor(seed,condition,backward,opt);
        stats<<w<<'\t'<<seed.dim(0)<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
        write(output/("w"+std::to_string(w)+".wxf"),seed);
    }stats.close();seal(output);
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
