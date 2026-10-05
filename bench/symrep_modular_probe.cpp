// Controlled diagnostics: modular prime size and sparse/dense eigenspace RREF.
// This does not change the production backend or its certificate prime.
#include "symrep_bootstrap.hpp"
#include <flint/nmod_mat.h>
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==7,"modular_probe PREPARED CHAIN WEIGHT PRIME_BITS orbit|sparse|dense THREADS");
    rref_option_t opt;opt->pool.reset(std::stoul(argv[6]));
    auto g=load_group(argv[1]);auto ms=load_models(argv[1],g,opt);auto p=load_prepared(argv[1],ms);
    fs::path chain=argv[2];size_t weight=std::stoul(argv[3]),bits=std::stoul(argv[4]);
    require(bits>=16&&bits<=60,"prime bits must be between 16 and 60");
    ulong prime=ulong(1)<<bits;do{prime=n_nextprime(prime,0);}while(prime%g.order()!=1);
    field_t fp(FIELD_Fp,prime);std::vector<ModMat> actions,letters;
    for(size_t j=0;j<g.generators.size();++j){
        actions.push_back(modular_matrix(read_matrix(chain/("w1_carrier_g"+std::to_string(j)+".wxf")),fp));
        letters.push_back(modular_matrix(g.generators[j],fp));}
    for(size_t w=2;w<=weight;++w)actions=carrier_actions_mod(flatten(read(chain/(weight_name(w)+"_carrier.wxf"))),actions,letters,p.backward,fp,opt);
    auto expected=load_layout(chain/(weight_name(weight)+"_copies.tsv"),ms);std::string mode=argv[5];
    auto tick=std::chrono::steady_clock::now();
    if(mode=="orbit"){
        CarrierFrameTimings phases;auto result=select_fixed_orbits_mod(g,ms,actions,fp,opt,&phases);
        require(result.layout.copies==expected.copies,"modular orbit multiplicities differ");
        std::cout<<"RESULT mode=orbit bits="<<bits<<" seconds="<<seconds(tick)<<" trace_s="<<phases.modular<<" selection_s="<<phases.selection<<std::endl;
    }else{
        require(mode=="sparse"||mode=="dense","unknown diagnostic mode");
        auto equation=actions[0];size_t n=equation.nrow;
        for(size_t i=0;i<n;++i){ModVec unit;unit.push_back(index(i),ulong(1));sparse_vec_sub_mul(equation[i],unit,ulong(1),fp);}
        size_t expected_rank=0;
        for(size_t s=0;s<ms.size();++s){auto a=ms[s].actions[g.generator_ids[0]];
            for(size_t i=0;i<a.nrow;++i){Vec unit;unit.push_back(index(i),Q(1));sparse_vec_sub_mul(a[i],unit,Q(1),QQ);}
            expected_rank+=row_basis(a).nrow*expected.copies[s];}
        size_t rank=0,input_nnz=equation.nnz();tick=std::chrono::steady_clock::now();
        if(mode=="sparse"){equation.sort_rows_by_nnz();auto pivots=sparse_mat_rref(equation,fp,opt);for(auto& batch:pivots)rank+=batch.size();}
        else{nmod_mat_t dense;nmod_mat_init(dense,n,n,prime);
            for(size_t i=0;i<n;++i)for(auto [j,x]:equation[i])nmod_mat_entry(dense,i,j)=x;
            rank=nmod_mat_rref(dense);nmod_mat_clear(dense);}
        require(rank==expected_rank,"eigenspace rank differs from rational model prediction");
        std::cout<<"RESULT mode="<<mode<<" bits="<<bits<<" size="<<n<<" input_nnz="<<input_nnz<<" rank="<<rank<<" seconds="<<seconds(tick)<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
