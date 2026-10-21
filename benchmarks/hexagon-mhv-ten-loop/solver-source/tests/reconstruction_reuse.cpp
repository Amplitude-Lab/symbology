#include "nmhv_sew_modular.hpp"
#include "irr_stage.hpp"
#include <map>

int main(int argc,char** argv){try {
    if(argc!=2&&argc!=3)throw std::runtime_error("reconstruction_reuse DIRECTORY [--resume]");
    std::filesystem::path root=argv[1];
    bool resume=argc==3&&std::string(argv[2])=="--resume";
    if(std::filesystem::exists(root)!=resume)throw std::runtime_error("Directory freshness/resume mismatch");
    Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(2);opt->verbose=false;
    int_t denominator=1;for(size_t i=0;i<80;++i)denominator*=2;denominator+=7;
    std::map<ulong,size_t> calls;
    auto evaluate=[&](ulong p,auto&,auto&){
        ++calls[p];nmhv_sew_prime_result_t<int32_t> out;
        out.basis=sparse_mat<ulong,int32_t>(1,2);out.basis[0].push_back(0,1);
        field_t F(FIELD_Fp,p);out.basis[0].push_back(1,n_invmod(denominator%F.mod,p));out.dimensions={2,1};out.seam_nnz={2};return out;
    };
    for(size_t count=1;count<=3;++count){
        bool complete=false;
        try{auto r=nmhv_sew_reconstruct_with_evaluator<int32_t>(2,3,1,count,root,123,opt,evaluate,true);
            if(r.basis.nrow!=1||r.basis[0][1]!=rat_t(1)/rat_t(denominator))throw std::runtime_error("Wrong reconstructed matrix");complete=true;
        }catch(const std::runtime_error& e){if(std::string(e.what())!="Rational reconstruction failed validation; increase --modular-primes.")throw;}
        if(complete!=(count==3))throw std::runtime_error("Expected two failed reconstructions then success");
    }
    size_t total=0;for(auto [p,n]:calls){std::cout<<"prime="<<p<<" evaluations="<<n<<'\n';total+=n;}
    if(total!=(resume?0:5))throw std::runtime_error("Repeated prime evaluation");
    bool rejected=false;try{nmhv_sew_reconstruct_with_evaluator<int32_t>(2,3,1,3,root,124,opt,evaluate,true);}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Changed cache identity accepted");
    {
        irr_stage::lock first(root/"test.lock",true),second(root/"test.lock",true);
        rejected=false;try{irr_stage::lock writer(root/"test.lock");}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("Writer ignored shared source locks");
    }
    {
        irr_stage::lock writer(root/"test.lock");
        rejected=false;try{irr_stage::lock reader(root/"test.lock",true);}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("Reader ignored writer lock");
    }
    std::cout<<"PASS: each prime evaluated once\n";
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
