#include "symrep_bootstrap.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==4||argc==5,"frame_probe PREPARED CHAIN WEIGHT [PIVOT_METHOD]");rref_option_t opt;opt->pool.reset(2);opt->verbose=false;
    if(argc==5)opt->method=std::stoi(argv[4]);
    auto g=load_group(argv[1]);auto ms=models(g,opt);size_t w=std::stoul(argv[3]);auto expected=load_layout(fs::path(argv[2])/(weight_name(w)+"_copies.tsv"),ms);
    auto gs=load_carrier_actions(argv[2],w,g,expected.dimension(ms),opt);CarrierFrameTimings timing;
    auto start=std::chrono::steady_clock::now();auto result=select_fixed_orbits(g,ms,gs,opt,&timing);
    require(result.layout.copies==expected.copies,"wrong orbit multiplicities");
    std::cout<<"total="<<seconds(start)<<" modular="<<timing.modular<<" selection="<<timing.selection<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
