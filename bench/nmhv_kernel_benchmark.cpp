// Both kernels in one executable, with identical I/O, workers and dependencies.
#include "tensor_kernel.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==8,"nmhv_kernel_benchmark extend|sew original|streamed|staged CONDITION SEED LAST_OR_WEIGHT OUTPUT THREADS");
    rref_option_t opt;opt->pool.reset(std::stoul(argv[7]));bool staged=std::string(argv[2])=="staged",streamed=staged||std::string(argv[2])=="streamed";
    require(streamed||std::string(argv[2])=="original","invalid strategy");auto d=read(argv[3]),seed=read(argv[4]);fs::path out=argv[6];
    if(std::string(argv[1])=="extend"){
        fresh(out);write(out/"w1.wxf",seed);std::ofstream stats(out/"timings.tsv");stats<<"weight\tdimension\tsolve_s\n";
        for(size_t w=2;w<=std::stoul(argv[5]);++w){auto start=std::chrono::steady_clock::now();
            if(streamed)seed=tensor_kernel::extend_tensor(seed,d,true,opt,staged?tensor_kernel::Strategy::staged:tensor_kernel::Strategy::streamed);
            else {auto copy=d;seed=extend_backward(std::move(copy),std::move(seed),QQ,opt);}
            stats<<w<<'\t'<<seed.dim(0)<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
            write(out/("w"+std::to_string(w)+".wxf"),seed);
        }stats.close();seal(out);
    }else{
        require(std::string(argv[1])=="sew","invalid operation");auto last=read(argv[5]);
        auto s=streamed?tensor_kernel::sew_one(std::move(d),std::move(seed),std::move(last),opt,staged?tensor_kernel::Strategy::staged:tensor_kernel::Strategy::streamed):sew_first_last(std::move(d),std::move(seed),std::move(last),QQ,opt);
        std::cout<<"SEW dimension="<<s.dim(0)<<std::endl;write(out,s);
    }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
