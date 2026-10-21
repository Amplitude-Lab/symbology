#include "hexagon_mrk_rhs_reuse.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc,char** argv)try{
    if(argc!=7)throw std::runtime_error("mrk_rhs_probe DATA LOOP PRIME SAMPLE MAX_DEGREE old|new|compare");
    Flint::set_memory_functions();const size_t L=std::stoul(argv[2]),sample=std::stoul(argv[4]),degree=std::stoul(argv[5]);
    const ulong prime=std::stoull(argv[3]);field_t F(FIELD_Fp,prime);std::string mode=argv[6];
    nmhv_physical_symbol_table_t data;
    for(size_t l=1;l<=L;++l)for(size_t d=0;d<l;++d){
        auto path=std::filesystem::path(argv[1])/"mrk_mhv"/("bfkl_f_L"+std::to_string(l)+"_k"+std::to_string(d)+".tsv");
        std::ifstream in(path);std::string line;std::getline(in,line);
        if(line!="word\tcoefficient")throw std::runtime_error("Bad MRK input header");
        auto& terms=data[{l,d}];
        while(std::getline(in,line))if(!line.empty()){
            auto a=nmhv_limit_split_tabs(line);nmhv_physical_symbol_term_t t;t.coefficient=rat_t(a.at(1));
            for(char c:a.at(0))t.word.push_back(c-'0');terms.push_back(std::move(t));
        }
    }
    auto slots=nmhv_physical_sample_mrk_slots(0x193a4527,sample,2*L-1,F);
    std::vector<ulong> old_value,new_value;
    if(mode=="old"||mode=="compare"){
        auto start=std::chrono::steady_clock::now();old_value=nmhv_physical_mrk_rhs(L,slots,data,degree,F);
        std::cout<<"MRK_OLD loop="<<L<<" prime="<<prime<<" sample="<<sample<<" degree="<<degree<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
    }
    if(mode=="new"||mode=="compare"){
        auto start=std::chrono::steady_clock::now();hexagon_mrk_reuse::evaluator e(L,slots,data,degree,F);new_value=e.evaluate();
        std::cout<<"MRK_NEW loop="<<L<<" prime="<<prime<<" sample="<<sample<<" degree="<<degree<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" memo_bytes="<<e.memo_bytes()<<std::endl;
    }
    if(mode=="compare"&&old_value!=new_value)throw std::runtime_error("MRK optimized/reference coefficient mismatch");
    const auto& value=new_value.empty()?old_value:new_value;if(value.empty())throw std::runtime_error("Unknown mode");
    std::cout<<"MRK_COEFFICIENTS";for(auto x:value)std::cout<<' '<<x;std::cout<<std::endl;
    std::cout<<"MRK_PROBE PASS mode="<<mode<<" coefficients="<<value.size()<<std::endl;
}catch(const std::exception& error){std::cerr<<"MRK_PROBE_ERROR "<<error.what()<<std::endl;return 1;}
