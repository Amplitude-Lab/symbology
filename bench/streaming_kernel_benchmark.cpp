#include "tensor_kernel.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc>=9&&argc<=11,"streaming_kernel_benchmark forward|backward|sew|forward-staged|backward-staged|sew-staged|profile|count|count-sew CONDITION SEED LAST_OR_DASH OUTPUT BATCH_ROWS METHOD THREADS [natural|short [SEED_TERMS]]");
    rref_option_t opt;opt->pool.reset(std::stoul(argv[8]));opt->method=std::stoi(argv[7]);
    size_t batch=std::stoul(argv[6]);std::string mode=argv[1];bool staged=mode.ends_with("-staged");
    if(staged)mode.resize(mode.size()-7);
    bool backward=mode=="backward",sew=mode=="sew"||mode=="count-sew";
    bool count=mode=="count"||mode=="count-sew",profile=mode=="profile"||count;require(backward||sew||profile||mode=="forward","unknown operation");
    require(!staged||(!profile&&batch==0),"staged mode requires a solve and BATCH_ROWS=0; its own policy bounds batches");
    bool order_short=batch!=0;
    if(argc>=10){std::string order=argv[9];require(order=="natural"||order=="short","equation order must be natural or short");order_short=order=="short";}
    staged_kernel::Options options;
    if(argc>=11){require(staged,"SEED_TERMS requires staged mode");options.seed_terms=std::stoul(argv[10]);require(options.seed_terms>0,"SEED_TERMS must be positive");}
    std::cout<<"POLICY staged="<<staged<<" order="<<(order_short?"short":"natural")<<" seed_terms="<<options.seed_terms<<std::endl;
    auto d=read(argv[2]),seed=read(argv[3]);size_t previous=seed.dim(0),next=d.dim(0),letters=d.dim(0);Mat local;
    if(sew){auto l=read(argv[4]);next=l.dim(0);local=tensor_kernel::right_boundary_conditions(d,l);}
    else{local=Mat(d.dim(2),letters*letters);sparse_tensor<Q,I,SPARSE_COO> dc(d);
        for(size_t p=0;p<dc.nnz();++p){auto i=dc.index(p);local[i[2]].push_back(index(i[backward?1:0]*letters+i[backward?0:1]),dc.val(p));}
        for(auto& row:local.rows)normalize(row);}
    tensor_kernel::ConstraintRows rows(seed,std::move(local),next,backward,order_short);seed.clear();d.clear();structured_kernel::Statistics stats;
    if(profile){structured_kernel::Relations relations(rows.ncol());
        for(size_t pass=0;pass<8;++pass){size_t changed=0;rows(relations,size_t(2),[&](Vec&& row){normalize(row);if(row.nnz()<=2)changed+=relations.consume(row);});
            std::cout<<"profile_presolve_pass="<<pass+1<<" relations="<<changed<<" remaining_columns="<<relations.dimension()<<std::endl;if(!changed)break;}
        rows.profile(relations,std::cout);
        if(count){std::vector<I> columns(rows.ncol(),-1);size_t n=0;
            for(size_t c=0;c<rows.ncol();++c){auto [r,x]=relations.image(c);if(x!=0&&columns[r]<0)columns[r]=index(n++);}
            field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));size_t nr=0,nnz=0;
            rows.modular(relations,columns,fp,[&](auto&& row){++nr;nnz+=row.nnz();},&opt->pool);
            std::cout<<"COUNT prime="<<fp.mod.n<<" rows="<<nr<<" cols="<<n<<" nnz="<<nnz<<" payload_bytes="<<nnz*(sizeof(ulong)+sizeof(I))<<std::endl;
        }return 0;}
    auto k=staged?staged_kernel::solve(staged_kernel::stack(rows),opt,options,&stats)
        :structured_kernel::solve_stream(rows.ncol(),rows,opt,&stats,[&](const Mat& candidate,ulong prime,const int_t& modulus){return rows.verify(candidate,opt,true,prime,&modulus);},batch);
    if(backward)for(auto& row:k.rows){for(size_t p=0;p<row.nnz();++p){size_t c=row(p);row(p)=index((c%letters)*previous+c/letters);}normalize(row);}
    for(auto& row:k.rows)vec_cancel_divisor(row);k.sort_rows_by_nnz();
    std::cout<<"KERNEL dimension="<<k.nrow<<" nnz="<<k.nnz()<<" core_cols="<<stats.core_cols<<" core_nnz="<<stats.core_nnz<<" presolve_s="<<stats.presolve_s<<" solve_s="<<stats.solve_s<<std::endl;
    write_tensor_view(argv[5],MatrixTensorView(k,backward?std::vector<size_t>{k.nrow,letters,previous}:std::vector<size_t>{k.nrow,previous,next}));
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
