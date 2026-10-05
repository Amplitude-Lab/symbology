#include "tensor_kernel.hpp"
#include "recursive_word_chart.hpp"
#include "early_sew.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc>=6&&argc<=9,"early_sew_probe actions|profile|count|rank|sew|refine|refine-cache|refine-sample|refine-hybrid DATA CHAIN MAX OUTPUT [MAX_TERMS [RIGHT_WEIGHT [raw|reduced]]]");
    fs::path data=argv[2],chain=argv[3],out=argv[5];fs::create_directories(out);
    size_t maximum=std::stoul(argv[4]),right_weight=argc>=8?std::stoul(argv[7]):2;require(right_weight>0,"right weight must be positive");rref_option_t opt;opt->pool.reset(8);
    std::vector<Mat> letters,previous;std::vector<std::string> names{"cyc","flip","parity"};
    for(auto name:names){letters.push_back(read_matrix(data/(name+"repmat.wxf")));previous.push_back(identity(1));}
    std::string mode=argv[1];
    if(mode!="actions"){
        auto start=std::chrono::steady_clock::now();auto d=read(data/"dlogmat_E6.wxf"),first=read(data/"FEC_1.wxf"),last=read(data/"LEC_1.wxf");
        // These checks establish that the complete recursive spaces are
        // preserved, so selected private columns recover their exact actions.
        auto rules=condition_rows(d);Chart rc(row_basis(rules),opt),fc(flatten(first),opt),lc(flatten(last),opt);
        for(const auto& g:letters){rc.coordinates(mul(rules,transpose(kron(g,g)),&opt->pool));fc.coordinates(mul(flatten(first),g));lc.coordinates(mul(flatten(last),g));}
        auto right=carrier_actions(flatten(last),previous,letters,true,opt);
        auto last2=std::move(last);
        for(size_t weight=2;weight<=right_weight;++weight){
            last2=tensor_kernel::extend_tensor(last2,d,true,opt);
            right=carrier_actions(flatten(last2),right,letters,true,opt);
        }
        auto right_name="LEC_"+std::to_string(right_weight)+".wxf";
        if(!fs::exists(out/right_name))write(out/right_name,last2);
        std::vector<Mat> left;
        for(auto name:names)left.push_back(read_matrix(out/"actions"/("w"+std::to_string(maximum)+"_"+name+".wxf")));
        auto f=read(chain/("FEC_"+std::to_string(maximum)+".wxf"));size_t nf=f.dim(0),nl=last2.dim(0);
        require(argc<9||std::string(argv[8])=="raw"||std::string(argv[8])=="reduced","invalid interface policy");
        auto local=tensor_kernel::right_boundary_conditions(d,last2,argc==9&&std::string(argv[8])=="raw"?tensor_kernel::BoundaryReduction::raw:tensor_kernel::BoundaryReduction::reduced);d.clear();
        tensor_kernel::ConstraintRows equations(f,std::move(local),nl);f.clear();
        early_sew::ProductInvarianceRows invariants(std::move(left),std::move(right));early_sew::Rows source{equations,invariants,argc>=7?std::stoul(argv[6]):SIZE_MAX};
        std::cout<<"EARLY_SOURCE columns="<<source.ncol()<<" setup_seconds="<<seconds(start)<<std::endl;
        if(mode=="profile"||mode=="count"||mode=="rank"){
            structured_kernel::Relations relations(source.ncol());for(size_t pass=0;pass<8;++pass){size_t changed=0;source(relations,2,[&](Vec&& row){if(row.nnz()<=2)changed+=relations.consume(row);});
                std::cout<<"EARLY_PRESOLVE pass="<<pass+1<<" changed="<<changed<<" columns="<<relations.dimension()<<" elapsed="<<seconds(start)<<std::endl;if(!changed)break;}
            equations.profile(relations,std::cout);
            if(mode=="count"||mode=="rank"){std::vector<I> columns(source.ncol(),-1);size_t n=0;for(size_t c=0;c<source.ncol();++c){auto [r,x]=relations.image(c);if(x!=0&&columns[r]<0)columns[r]=index(n++);}
                field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));size_t nr=0,nnz=0;sparse_mat<ulong,I> a(0,n);
                source.modular(relations,columns,fp,[&](auto&& row){++nr;nnz+=row.nnz();if(mode=="rank")a.rows.push_back(std::move(row));},&opt->pool);
                std::cout<<"EARLY_COUNT rows="<<nr<<" cols="<<n<<" nnz="<<nnz<<" payload_bytes="<<nnz*12<<" elapsed="<<seconds(start)<<std::endl;
                if(mode=="rank"){a.nrow=a.rows.size();a.sort_rows_by_nnz();auto pivots=sparse_mat_rref_forward(a,fp,opt);size_t rank=0;for(auto& p:pivots)rank+=p.size();std::cout<<"EARLY_RANK rank="<<rank<<" nullity="<<n-rank<<" elapsed="<<seconds(start)<<std::endl;}}
            return 0;
        }
        require(mode=="sew"||mode=="refine"||mode=="refine-cache"||mode=="refine-sample"||mode=="refine-hybrid","unknown early sewing mode");
        auto verify=[&](const Mat& k,ulong prime,const int_t& modulus){
            auto candidate_path=out/("CANDIDATE_"+std::to_string(maximum)+"p"+std::to_string(right_weight)+"_"+std::to_string(prime)+".wxf");
            write_tensor_view(candidate_path,MatrixTensorView(k,{k.nrow,nf,nl}));std::ofstream metadata(candidate_path.string()+".modulus");metadata<<prime<<'\n'<<modulus.get_str()<<'\n';metadata.close();require(bool(metadata),"cannot write candidate congruence metadata");
            std::cout<<"UNCERTIFIED_CANDIDATE "<<candidate_path<<std::endl;
            return early_sew::verify_with_extra_primes(k,source,opt,prime,modulus);};Mat result;
        if(mode=="sew"){require(source.max_terms==SIZE_MAX,"selected rows cannot use a full-equation certificate");result=structured_kernel::solve_stream(source.ncol(),source,opt,nullptr,verify);}
        else{early_sew::RefinedRows refined{source,opt,source.max_terms==SIZE_MAX?64:source.max_terms,128,mode=="refine-cache",mode=="refine-sample"||mode=="refine-hybrid",mode=="refine-hybrid"?1024u:0u};result=structured_kernel::solve_stream(source.ncol(),refined,opt,nullptr,verify);}
        write_tensor_view(out/("EARLY_"+std::to_string(maximum)+"p"+std::to_string(right_weight)+".wxf"),MatrixTensorView(result,{result.nrow,nf,nl}));
        std::cout<<"EARLY_RESULT dimension="<<result.nrow<<" nnz="<<result.nnz()<<" elapsed="<<seconds(start)<<std::endl;return 0;
    }
    for(size_t w=1;w<=maximum;++w){auto path=w==1?data/"FEC_1.wxf":chain/("FEC_"+std::to_string(w)+".wxf");
        bool cached=true;for(auto name:names)cached &= fs::exists(out/("w"+std::to_string(w)+"_"+name+".wxf"));
        if(cached){previous.clear();for(auto name:names)previous.push_back(read_matrix(out/("w"+std::to_string(w)+"_"+name+".wxf")));continue;}
        auto start=std::chrono::steady_clock::now();auto t=read(path);auto basis=flatten(t);t.clear();
        auto next=carrier_actions(basis,previous,letters,false,opt);basis.clear();
        for(size_t j=0;j<names.size();++j){size_t monomial=0;for(const auto& row:next[j].rows)monomial+=row.nnz()==1;
            std::cout<<"ACTION weight="<<w<<" name="<<names[j]<<" dimension="<<next[j].nrow<<" nnz="<<next[j].nnz()<<" monomial_rows="<<monomial<<" elapsed="<<seconds(start)<<std::endl;
            write_tensor_view(out/("w"+std::to_string(w)+"_"+names[j]+".wxf"),MatrixTensorView(next[j],{next[j].nrow,next[j].ncol}));}
        previous=std::move(next);
    }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
