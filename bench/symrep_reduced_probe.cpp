// Dump reduced recurrence blocks or compare SparseRREF strategies on identical inputs.
#include "symrep_bootstrap.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc>=2,"expected dump PREPARED CHAIN WEIGHT OUTPUT, or solve BLOCKS METHOD THREADS");
    rref_option_t opt;
    if(std::string(argv[1])=="dump"){
        require(argc==6,"dump PREPARED CHAIN WEIGHT OUTPUT");opt->pool.reset(2);fs::path prepared=argv[2],chain=argv[3],out=argv[5];size_t w=std::stoul(argv[4]);
        require(w>1&&w<=chain_weight(chain,prepared)+1,"invalid weight");fresh(out);auto g=load_group(prepared);auto ms=load_models(prepared,g,opt);auto p=load_prepared(prepared,ms);
        TensorAdapter adapter(g,ms,opt);RecurrenceTemplates recurrence(g,ms,p,adapter,opt);
        auto previous=load_layout(chain/(weight_name(w-1)+"_copies.tsv"),ms),older=w==2?p.terminal:load_layout(chain/(weight_name(w-2)+"_copies.tsv"),ms);
        ProductCopies source,target;auto domain=p.backward?adapter.product(p.alphabet,previous,&source,false):adapter.product(previous,p.alphabet,&source,false);
        auto codomain=p.backward?adapter.product(p.equations,older,&target,false):adapter.product(older,p.equations,&target,false);
        auto seed=materialize_multiplicity(chain,w-1,p,ms,adapter,opt);auto blocks=recurrence.assemble(seed,older,previous,source,target);
        for(size_t s=0;s<blocks.size();++s)write(out/("block"+std::to_string(s)+".wxf"),blocks[s]);
    }else{
        require(std::string(argv[1])=="solve"&&(argc==5||argc==6),"solve BLOCKS METHOD THREADS [index|degree|reverse-degree|row-cost]");opt->method=std::stoi(argv[3]);opt->pool.reset(std::stoul(argv[4]));
        std::string policy=argc==6?argv[5]:"index";require(policy=="index"||policy=="degree"||policy=="reverse-degree"||policy=="row-cost","unknown column ordering");
        double total=0;size_t nnz=0;std::vector<fs::path> files;for(auto& e:fs::directory_iterator(argv[2]))if(e.path().extension()==".wxf")files.push_back(e.path());std::sort(files.begin(),files.end());
        for(const auto& file:files){auto a=read_matrix(file);std::vector<int64_t> costs(a.ncol);
            for(const auto& row:a.rows)for(auto [j,x]:row)costs[j]+=policy=="row-cost"?row.nnz():1;
            if(policy=="reverse-degree"&&!costs.empty()){auto maximum=*std::max_element(costs.begin(),costs.end());for(auto& c:costs)c=maximum-c;}
            if(policy!="index")opt->col_weight=[&](int64_t j){return costs[j];};
            auto start=std::chrono::steady_clock::now();auto k=kernel(a,opt);double elapsed=seconds(start);total+=elapsed;nnz+=k.nnz();
            require(!mul(a,transpose(k),&opt->pool).nnz(),"nonzero residual");std::cout<<file.filename()<<" seconds="<<elapsed<<" nullity="<<k.nrow<<" kernel_nnz="<<k.nnz()<<'\n';}
        std::cout<<"TOTAL seconds="<<total<<" kernel_nnz="<<nnz<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
