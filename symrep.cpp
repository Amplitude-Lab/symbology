#include "symrep_bootstrap.hpp"
#include <iostream>
using namespace symrep;

struct Args {
    std::map<std::string,std::string> values;
    std::vector<std::string> names,terminal_names;
    std::vector<Mat> gens,terminal_gens;
    bool compare=false;
    Args(int argc,char** argv) {
        for(int i=2;i<argc;++i){std::string k=argv[i];if(k=="--compare"){compare=true;continue;}
            require(i+1<argc,"missing value for "+k);std::string v=argv[++i];
            if(k=="--generator"||k=="--terminal-generator") {
                auto at=v.find('=');require(at!=std::string::npos&&at>0,"generator syntax is NAME=FILE.wxf");auto name=v.substr(0,at);
                require(name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos,"invalid generator name");
                (k=="--generator"?names:terminal_names).push_back(name);
                (k=="--generator"?gens:terminal_gens).push_back(read_matrix(v.substr(at+1)));
            }else{
                require(k=="--output"||k=="--input"||k=="--matrix"||k=="--condition"||k=="--seed"||k=="--direction"||k=="--threads"||k=="--max-order"||k=="--max-weight"||k=="--chain"||k=="--reference"||k=="--kernel-strategy"||k=="--backend"||k=="--coordinates","unknown option "+k);
                require(values.emplace(k,v).second,"duplicate option "+k);
            }
        }
    }
    std::string get(const std::string& k) const {auto it=values.find(k);require(it!=values.end(),"required option "+k);return it->second;}
    size_t number(const std::string& k,size_t fallback) const {auto it=values.find(k);if(it==values.end())return fallback;auto n=strict_integer(it->second);require(n>0,"expected positive "+k);return n;}
};
int main(int argc,char** argv){try {
    if(argc<2||std::string(argv[1])=="--help") {
        std::cout<<"Exact QQ symmetry kernel solver (row actions v -> v R).\n"
          "symrep kernel --matrix A.wxf --generator NAME=R.wxf [...] --output DIR [--threads N]\n"
          "symrep kernel --matrix A.wxf --input PREVIOUS_KERNEL --output DIR\n"
          "symrep adapt --generator NAME=R.wxf [...] --output DIR\n"
          "symrep prepare --generator NAME=R.wxf [...] --condition D.wxf --seed FEC1.wxf --direction forward|backward --output DIR\n"
          "symrep extend --input PREPARED --max-weight N --output DIR [--compare] [--chain PREVIOUS] [--backend auto|multiplicity|factorized] [--kernel-strategy original|streamed|staged]\n"
          "symrep expand --input PREPARED --chain FACTORIZED --output DIR [--max-weight N]\n"
          "symrep actions --input PREPARED --chain DIR --output DIR [--max-weight N] [--coordinates adapted|carrier]\n"
          "symrep ordinary --condition D.wxf --seed FEC1.wxf --direction forward|backward --max-weight N --output DIR\n"
          "symrep verify --input PREPARED --chain DIR --reference DIR --max-weight N --output DIR\n"
          "Default factorized kernel: streamed; split multiplicity backend is retained.\nAll commands: --threads N, --max-order N. prepare: optional --terminal-generator NAME=R.wxf for each generator.\n"
          "Outputs are QQ irreducible copies, including Galois-conjugate components together.\n";return 0;
    }
    std::string command=argv[1];require(command=="kernel"||command=="adapt"||command=="prepare"||command=="extend"||command=="verify"||command=="ordinary"||command=="expand"||command=="actions","unknown command");
    Args args(argc,argv);rref_option_t opt;opt->pool.reset(args.number("--threads",2));opt->verbose=false;
    const std::map<std::string,std::set<std::string>> allowed{
        {"adapt",{}},{"kernel",{"--matrix","--input"}},{"prepare",{"--condition","--seed","--direction"}},
        {"extend",{"--input","--max-weight","--chain","--backend","--kernel-strategy"}},{"expand",{"--input","--chain","--max-weight"}},{"verify",{"--input","--chain","--reference","--max-weight"}},
        {"ordinary",{"--condition","--seed","--direction","--max-weight"}},
        {"actions",{"--input","--chain","--max-weight","--coordinates"}}};
    for(const auto& [key,value]:args.values)require(key=="--output"||key=="--threads"||key=="--max-order"||allowed.at(command).count(key),key+" is not valid for "+command);
    require(!args.compare||command=="extend","--compare is only valid for extend");
    require(args.terminal_gens.empty()||command=="prepare","terminal generators are only valid for prepare");
    require(args.gens.empty()||command=="adapt"||command=="kernel"||command=="prepare","generator arguments are not valid for "+command);
    auto start=std::chrono::steady_clock::now();size_t limit=args.number("--max-order",256);
    if(command=="ordinary") {
        auto direction=args.get("--direction");require(direction=="forward"||direction=="backward","direction must be forward or backward");
        fs::path output=args.get("--output");fresh(output);
        ordinary_chain(read(args.get("--seed")),read(args.get("--condition")),direction=="backward",args.number("--max-weight",2),output,opt);
        std::cout<<"complete total_s="<<seconds(start)<<std::endl;return 0;
    }
    bool saved=command=="extend"||command=="verify"||command=="expand"||command=="actions"||(command=="kernel"&&args.values.count("--input"));
    require(!saved||args.gens.empty(),"--input cannot be combined with --generator");
    auto group=saved?load_group(args.get("--input"),limit):Group(args.names,args.gens,limit);auto ms=saved?load_models(args.get("--input"),group,opt):models(group,opt);
    std::cout<<"group_order="<<group.order()<<" rational_irreps="<<ms.size()<<" group_setup_s="<<seconds(start)<<std::endl;
    fs::path output=args.get("--output");fresh(output);
    if(command=="adapt") {
        auto a=adapt(group,ms,args.gens,opt);verify_adaptation(a,group,ms,args.gens);save_group(output,group,ms);save_space(output,"space",a,group,ms);seal(output);
    }else if(command=="prepare") {
        auto direction=args.get("--direction");require(direction=="forward"||direction=="backward","direction must be forward or backward");
        require(args.terminal_names.empty()||args.terminal_names==args.names,"terminal generator names/order must match alphabet generators");
        prepare(group,ms,read(args.get("--condition")),read(args.get("--seed")),direction=="backward",opt,output,args.terminal_gens);
    }else if(command=="extend") {
        auto prepared=fs::path(args.get("--input"));
        auto resume=args.values.count("--chain")?fs::path(args.get("--chain")):fs::path{};
        auto backend=args.values.count("--backend")?args.get("--backend"):"auto";
        require(backend=="auto"||backend=="multiplicity"||backend=="factorized","backend must be auto, multiplicity, or factorized");
        bool factorized=backend=="factorized";
        if(backend=="auto"){
            if(!resume.empty())chain_weight(resume,prepared,&factorized);
            else factorized=std::any_of(ms.begin(),ms.end(),[](const Model& m){return m.division_dimension>1;});
        }
        auto p=load_prepared(prepared,ms,!factorized);
        std::cout<<"backend="<<(factorized?"factorized":"multiplicity")<<std::endl;
        auto strategy=args.values.count("--kernel-strategy")?args.get("--kernel-strategy"):(factorized?"streamed":"original");
        require(strategy=="original"||strategy=="streamed"||strategy=="staged","kernel strategy must be original, streamed, or staged");
        require(strategy=="original"||factorized,"streamed kernels require the factorized backend");
        if(factorized)extend_carrier(p,group,ms,args.number("--max-weight",2),output,opt,args.compare,prepared,resume,strategy!="original",strategy=="staged"?tensor_kernel::Strategy::staged:tensor_kernel::Strategy::streamed);
        else if(split_models(ms))extend_multiplicity(p,group,ms,args.number("--max-weight",2),output,opt,args.compare,prepared,resume);
        else extend(p,group,ms,args.number("--max-weight",2),output,opt,args.compare,prepared,resume);
    }else if(command=="verify") {
        bool factorized=false;chain_weight(args.get("--chain"),args.get("--input"),&factorized);
        if(factorized)verify_carrier_chain(args.get("--input"),args.get("--chain"),args.get("--reference"),args.number("--max-weight",2),group,ms,opt,output);
        else verify_chain(args.get("--input"),args.get("--chain"),args.get("--reference"),args.number("--max-weight",2),group,ms,opt,output);
    }else if(command=="expand") {
        bool compact=false;auto maximum=chain_weight(args.get("--chain"),args.get("--input"),nullptr,&compact);
        if(compact)expand_multiplicity_chain(args.get("--input"),args.get("--chain"),args.number("--max-weight",maximum),group,ms,opt,output);
        else expand_carrier_chain(args.get("--input"),args.get("--chain"),args.number("--max-weight",maximum),group,ms,opt,output);
    }else if(command=="actions") {
        auto maximum=chain_weight(args.get("--chain"),args.get("--input"));
        auto coordinates=args.values.count("--coordinates")?args.get("--coordinates"):"adapted";
        require(coordinates=="adapted"||coordinates=="carrier","coordinates must be adapted or carrier");
        export_chain_actions(args.get("--input"),args.get("--chain"),args.number("--max-weight",maximum),group,ms,opt,output,coordinates=="carrier");
    }else {
        auto a=read_matrix(args.get("--matrix"));Adapted source;std::vector<Mat> domain;
        if(saved) {
            source.layout=load_layout(fs::path(args.get("--input"))/"kernel_copies.tsv",ms);
            source.basis=identity(source.layout.dimension(ms));domain=generators(group,ms,source.layout);
        }else {source=adapt(group,ms,args.gens,opt);domain=args.gens;}
        require(a.ncol==source.basis.ncol,"matrix columns must match generator dimension");
        verify_adaptation(source,group,ms,domain);
        auto ac=mul(a,transpose(source.basis),&opt->pool);auto result=kernel_adapted(ac,group,ms,source.layout,opt);
        auto basis=mul(result.solution.basis,source.basis,&opt->pool);
        verify_adaptation({basis,result.solution.layout},group,ms,domain);
        require(!mul(a,transpose(basis),&opt->pool).nnz(),"kernel residual");
        save_group(output,group,ms);save_space(output,"source",source,group,ms);
        save_space(output,"kernel",{basis,result.solution.layout},group,ms);
        write(output/"kernel_adapted_coordinates.wxf",result.solution.basis);seal(output);
        std::cout<<"kernel_dimension="<<basis.nrow<<std::endl;
    }
    std::cout<<"complete total_s="<<seconds(start)<<" output="<<output<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
