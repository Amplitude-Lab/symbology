#include "symrep_bootstrap.hpp"
#include "column_restriction.hpp"
#include "partitioned_tensor_view.hpp"
using namespace symrep;
// Full extension under rational one-dimensional characters. All sectors are
// retained; the original factored operator certifies each restricted kernel.
int main(int argc,char** argv){try{
    require(argc>=7,"split_character_extension PREPARED CHAIN_OR_DASH OUTPUT THREADS staged|assembled|staged-raw|assembled-raw SEED_TERMS [compare] [front_nonzeros=N] [method=N]");
    fs::path prepared=argv[1],resume=argv[2],root=argv[3];fresh(root);
    rref_option_t opt;opt->pool.reset(std::stoul(argv[4]));
    std::string strategy=argv[5];bool raw=strategy.ends_with("-raw");if(raw)strategy.resize(strategy.size()-4);
    require(strategy=="staged"||strategy=="assembled","unknown strategy");
    bool compare=false;staged_kernel::Options options;options.seed_terms=std::stoul(argv[6]);require(options.seed_terms,"positive initial support required");
    std::set<std::string> seen;
    for(int i=7;i<argc;++i){std::string arg=argv[i];auto at=arg.find('=');auto key=arg.substr(0,at);require(seen.insert(key).second,"duplicate option");
        if(arg=="compare")compare=true;
        else if(key=="front_nonzeros"&&at!=std::string::npos){auto n=strict_integer(arg.substr(at+1));require(n>0,"positive front size required");options.seed_nonzeros=size_t(n);}
        else if(key=="method"&&at!=std::string::npos){auto n=strict_integer(arg.substr(at+1));require(n>=0&&n<=2,"method must be 0, 1 or 2");opt->method=int(n);}
        else require(false,"unknown option "+arg);
    }
    std::cout<<"POLICY seed_terms="<<options.seed_terms<<" front_nonzeros="<<options.seed_nonzeros<<" method="<<opt->method<<std::endl;
    auto group=load_group(prepared,256);auto ms=load_models(prepared,group,opt);auto p=load_prepared(prepared,ms);
    require(!p.backward,"pilot currently supports forward extension");
    for(const auto& m:ms)require(m.dimension==1&&m.division_dimension==1,"character extension needs one-dimensional rational models");
    std::vector<std::vector<size_t>> product(ms.size(),std::vector<size_t>(ms.size(),SIZE_MAX));
    for(size_t a=0;a<ms.size();++a)for(size_t b=0;b<ms.size();++b)for(size_t c=0;c<ms.size();++c){bool match=true;
        for(auto id:group.generator_ids)match=match&&equal(mul(ms[a].actions[id],ms[b].actions[id]),ms[c].actions[id]);
        if(match){require(product[a][b]==SIZE_MAX,"ambiguous character product");product[a][b]=c;}}
    for(auto& row:product)for(auto c:row)require(c<ms.size(),"character products not closed");
    auto labels=[&](const Layout& layout){std::vector<size_t> out;for(size_t s=0;s<ms.size();++s)out.insert(out.end(),layout.copies[s],s);return out;};
    size_t weight=2;auto previous=p.seed,older=p.terminal;Tensor seed=p.expansion;TensorAdapter adapter(group,ms,opt);
    if(resume!="-"){
        bool factorized=false,compact=false;size_t last=chain_weight(resume,prepared,&factorized,&compact);require(!factorized,"character extension needs adapted chain");weight=last+1;
        for(size_t w=1;w<=last;++w){auto stem=weight_name(w);auto t=compact?materialize_multiplicity(resume,w,p,ms,adapter,opt):read(resume/(stem+".wxf"));
            if(w==1)require(equal(flatten(t),flatten(p.expansion)),"saved seed differs from prepared seed");
            write(root/(stem+".wxf"),t);fs::copy_file(resume/(stem+"_copies.tsv"),root/(stem+"_copies.tsv"));if(w==last)seed=std::move(t);}
        previous=load_layout(resume/(weight_name(last)+"_copies.tsv"),ms);if(last>1)older=load_layout(resume/(weight_name(last-1)+"_copies.tsv"),ms);
    }else{write(root/"w1.wxf",seed);save_layout(root/"w1_copies.tsv",previous,ms);}
    const auto pl=labels(previous),ol=labels(older),al=labels(p.alphabet),ql=labels(p.equations);
    require(seed.dim(0)==pl.size()&&seed.dim(1)==ol.size()&&seed.dim(2)==al.size(),"recurrence layout dimensions");
    for(size_t b=0;b<seed.dim(0);++b)for(size_t j=seed.rowptr()[b];j<seed.rowptr()[b+1];++j){auto i=seed.index(j);require(product[ol[i[0]]][al[i[1]]]==pl[b],"recurrence mixes characters");}
    auto local=condition_rows(p.condition);
    for(size_t q=0;q<local.nrow;++q)for(auto [j,x]:local[q])require(product[al[j/al.size()]][al[j%al.size()]]==ql[q],"condition mixes characters");
    if(raw){
        auto letters=read_matrix(prepared/"alphabet_basis.wxf");
        auto original=mul(read_matrix(prepared/"condition_independent_rows.wxf"),transpose(kron(letters,letters)),&opt->pool);
        auto eq=read_matrix(prepared/"equations_basis.wxf");
        require(equal(original,mul(transpose(eq),local,&opt->pool)),"original equation basis does not match prepared map");
        std::cout<<"CONDITION_BASIS adapted_nnz="<<local.nnz()<<" original_nnz="<<original.nnz()<<std::endl;
        // The checked invertible equation basis change preserves the kernel.
        // Equations may now mix characters; column restriction performs their
        // character projections without requiring a dense adapted row basis.
        local=std::move(original);
    }
    auto owned_rows=std::make_unique<tensor_kernel::ConstraintRows>(seed,std::move(local),al.size(),false,strategy=="staged");seed.clear();p.condition.clear();
    auto& rows=*owned_rows;Mat complete(0,rows.ncol());std::vector<PartitionedTensorView::Part> parts;
    size_t total_rows=0,total_nonzeros=0;Layout solution;solution.copies.resize(ms.size());std::ofstream stats(root/"timings.tsv");
    stats<<"weight\tsector\tcolumns\tdimension\tnnz\tseconds\n";
    for(size_t s=0;s<ms.size();++s){std::vector<I> columns;
        for(size_t b=0;b<pl.size();++b)for(size_t l=0;l<al.size();++l)if(product[pl[b]][al[l]]==s)columns.push_back(index(b*al.size()+l));
        staged_kernel::ColumnRestriction restricted(rows,std::move(columns));auto tick=std::chrono::steady_clock::now();
        std::cout<<"SECTOR_BEGIN weight="<<weight<<" sector="<<s<<" columns="<<restricted.ncol()<<std::endl;
        auto k=strategy=="staged"?staged_kernel::solve(staged_kernel::stack(restricted),opt,options):
            structured_kernel::solve_stream(restricted.ncol(),restricted,opt,nullptr,[&](const Mat& candidate,ulong prime,const int_t& modulus){return staged_kernel::verify(candidate,restricted,opt,prime,modulus);});
        solution.copies[s]=k.nrow;restricted.lift_in_place(k);
        std::cout<<"SECTOR_COMPLETE weight="<<weight<<" sector="<<s<<" dimension="<<k.nrow<<" nnz="<<k.nnz()<<" seconds="<<seconds(tick)<<std::endl;
        stats<<weight<<'\t'<<s<<'\t'<<restricted.ncol()<<'\t'<<k.nrow<<'\t'<<k.nnz()<<'\t'<<seconds(tick)<<std::endl;
        auto checkpoint=root/(weight_name(weight)+"_sector"+std::to_string(s)+".wxf");
        write_tensor_view(checkpoint,MatrixTensorView(k,{k.nrow,pl.size(),al.size()}));
        std::ofstream receipt(root/(weight_name(weight)+"_sector"+std::to_string(s)+".tsv"));
        receipt<<"certified_sector "<<s<<" dimension "<<k.nrow<<" nnz "<<k.nnz()<<" sha256 "<<native_cache::digest(checkpoint)<<'\n';
        receipt.close();require(bool(receipt),"sector receipt failed");
        PartitionedTensorView::Part part{checkpoint,{}};part.row_nonzeros.reserve(k.nrow);
        for(const auto& row:k.rows)part.row_nonzeros.push_back(row.nnz());parts.push_back(std::move(part));
        total_rows+=k.nrow;total_nonzeros+=k.nnz();
        if(compare){for(auto& row:k.rows)complete.rows.push_back(std::move(row));complete.nrow=complete.rows.size();}
    }
    if(compare){auto ordinary=structured_kernel::solve_stream(rows.ncol(),rows,opt,nullptr,[&](const Mat& k,ulong prime,const int_t& modulus){return rows.verify(k,opt,true,prime,&modulus);});
        // Each sector already has an exact residual certificate and private
        // free columns. Disjoint supports make their union independent, so
        // the independently computed full nullity proves equal row spaces.
        require(ordinary.nrow==complete.nrow,"full-space dimension mismatch");std::cout<<"EXACT_FULL_SPACE_PASS"<<std::endl;}
    complete.clear();owned_rows.reset();
    PartitionedTensorView output(std::move(parts),{pl.size(),al.size()});
    write_tensor_view(root/(weight_name(weight)+".wxf"),output);
    save_layout(root/(weight_name(weight)+"_copies.tsv"),solution,ms);stats.close();require(bool(stats),"timing output failed");
    std::ofstream meta(root/"chain.tsv");meta<<"symbology-symrep-chain-v1 prepared_sha256 "<<native_cache::digest(prepared/"checksums.tsv")<<' '<<weight<<'\n';meta.close();require(bool(meta),"chain output failed");seal(root);
    std::cout<<"FULL_CHARACTER_BASIS_COMPLETE weight="<<weight<<" dimension="<<total_rows<<" nnz="<<total_nonzeros<<" sectors="<<ms.size()<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
