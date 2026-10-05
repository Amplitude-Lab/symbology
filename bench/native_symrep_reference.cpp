// Fair baseline: call the existing bootstrap.hpp implementation with the same
// SparseRREF options/compiler/worker count as symrep. No symmetry conversion.
#include "symrep_io.hpp"
using namespace symrep;
int main(int argc,char** argv){try {
    require(argc==7,"native_symrep_reference forward|backward CONDITION SEED MAX OUTPUT THREADS");
    bool backward=std::string(argv[1])=="backward";auto condition=read(argv[2]),seed=read(argv[3]);
    size_t maximum=std::stoul(argv[4]);fs::path output=argv[5];fresh(output);
    rref_option_t opt;opt->pool.reset(std::stoul(argv[6]));opt->verbose=false;
    write(output/"w1.wxf",seed);std::ofstream stats(output/"timings.tsv");stats<<"weight\tdimension\tsolve_s\n";
    for(size_t w=2;w<=maximum;++w){auto start=std::chrono::steady_clock::now();auto d=condition;
        seed=backward?extend_backward(std::move(d),std::move(seed),QQ,opt):extend_forward(std::move(d),std::move(seed),QQ,opt);
        stats<<w<<'\t'<<seed.dim(0)<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
        write(output/("w"+std::to_string(w)+".wxf"),seed);
    }stats.close();seal(output);
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
