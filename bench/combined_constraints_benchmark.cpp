// Compare solving orders on identical, saved FEC/LEC factors and actions.
// The shared solver contains no polygon, loop-count or expected-nullity input.
#include "symrep_bootstrap.hpp"
#include "product_invariance.hpp"
using namespace symrep;

struct Input {
    Tensor condition, left, right;
    std::vector<Mat> left_actions, right_actions;
};

Input load_input(const fs::path& root) {
    verify_seal(root);
    Input in{read(root/"condition.wxf"),read(root/"left.wxf"),read(root/"right.wxf"),{}, {}};
    std::ifstream meta(root/"input.tsv");std::string schema,extra;size_t count;
    meta>>schema>>count;
    require(bool(meta)&&schema=="combined-constraints-v1"&&count>0&&count<1024,"invalid benchmark input schema");
    require(!(meta>>extra),"trailing input metadata");
    for(size_t g=0;g<count;++g){
        in.left_actions.push_back(read_matrix(root/("left_g"+std::to_string(g)+".wxf")));
        in.right_actions.push_back(read_matrix(root/("right_g"+std::to_string(g)+".wxf")));
        auto& a=in.left_actions.back();auto& b=in.right_actions.back();
        require(a.nrow==in.left.dim(0)&&a.ncol==a.nrow&&b.nrow==in.right.dim(0)&&b.ncol==b.nrow,"action dimensions differ from saved bases");
    }
    return in;
}

void save_input(const fs::path& root,const Input& in) {
    write(root/"condition.wxf",in.condition);write(root/"left.wxf",in.left);write(root/"right.wxf",in.right);
    for(size_t g=0;g<in.left_actions.size();++g){
        write(root/("left_g"+std::to_string(g)+".wxf"),in.left_actions[g]);
        write(root/("right_g"+std::to_string(g)+".wxf"),in.right_actions[g]);
    }
    std::ofstream meta(root/"input.tsv");meta<<"combined-constraints-v1 "<<in.left_actions.size()<<'\n';meta.close();
    require(bool(meta),"cannot write benchmark input metadata");seal(root);
}

Input prepare(Tensor d,Tensor f,Tensor l,const std::vector<Mat>& letters,size_t left_weight,size_t right_weight,rref_option_t opt) {
    require(left_weight>0&&right_weight>0&&!letters.empty(),"positive weights and generators required");
    auto rules=condition_rows(d);Chart rc(row_basis(rules),opt),fc(flatten(f),opt),lc(flatten(l),opt);
    std::vector<Mat> terminal(letters.size(),identity(1));
    for(const auto& g:letters){
        require(g.nrow==d.dim(0)&&g.ncol==g.nrow,"alphabet action dimensions");
        rc.coordinates(mul(rules,transpose(kron(g,g)),&opt->pool));
        fc.coordinates(mul(flatten(f),g));lc.coordinates(mul(flatten(l),g));
    }
    auto fa=carrier_actions(flatten(f),terminal,letters,false,opt);
    auto la=carrier_actions(flatten(l),terminal,letters,true,opt);
    for(size_t w=2;w<=left_weight;++w){f=tensor_kernel::extend_tensor(f,d,false,opt);fa=carrier_actions(flatten(f),fa,letters,false,opt);}
    for(size_t w=2;w<=right_weight;++w){l=tensor_kernel::extend_tensor(l,d,true,opt);la=carrier_actions(flatten(l),la,letters,true,opt);}
    return {std::move(d),std::move(f),std::move(l),std::move(fa),std::move(la)};
}

// Apply K (A tensor B) - K by two factor contractions, without a Kronecker matrix.
// Return its transpose: equations in the coefficients of the kernel rows.
Mat restricted_symmetry(const Mat& k,const Mat& a,const Mat& b,thread_pool* pool) {
    Mat by_left(a.nrow,k.nrow*b.nrow);
    for(size_t s=0;s<k.nrow;++s)for(auto [c,x]:k[s])by_left[size_t(c)/b.nrow].push_back(index(s*b.nrow+size_t(c)%b.nrow),x);
    auto left=mul(transpose(a),by_left,pool);by_left.clear();Mat by_right(k.nrow*a.nrow,b.nrow);
    for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:left[i])by_right[(size_t(j)/b.nrow)*a.nrow+i].push_back(index(size_t(j)%b.nrow),x);
    left.clear();auto transformed=mul(by_right,b,pool);by_right.clear();Mat residual(k.ncol,k.nrow);
    for(size_t s=0;s<k.nrow;++s){
        for(size_t i=0;i<a.nrow;++i)for(auto [j,x]:transformed[s*a.nrow+i])residual[i*b.nrow+j].push_back(index(s),x);
        for(auto [c,x]:k[s])residual[c].push_back(index(s),-x);
    }
    for(auto& row:residual.rows)normalize(row);return residual;
}

template<class Stack> struct CompleteRows {
    const Stack& source;
    template<class Emit> void operator()(structured_kernel::Relations& r,size_t limit,Emit&& emit)const{source(r,limit,emit);}
    template<class Emit> void modular(structured_kernel::Relations& r,const std::vector<I>& columns,const field_t& fp,Emit&& emit,thread_pool* pool)const{
        auto seen=source.masks(false);source.discover(r,columns,fp,staged_kernel::Selection{},seen,emit,pool);
    }
};

template<class Stack> Mat global_solve(const Stack& source,rref_option_t opt,structured_kernel::Statistics* stats=nullptr) {
    CompleteRows<Stack> rows{source};
    return structured_kernel::solve_stream(source.ncol(),rows,opt,stats,[&](const Mat& k,ulong p,const int_t& m){return staged_kernel::verify(k,source,opt,p,m);});
}

int main(int argc,char** argv){try{
    require(argc>=2,"combined_constraints_benchmark: prepare|toy|run|check (see docs/combined-constraints-benchmark.md)");
    std::string mode=argv[1];rref_option_t opt;
    if(mode=="prepare"){
        require(argc>=11,"prepare CONDITION FIRST_SEED LAST_SEED LEFT_WEIGHT RIGHT_WEIGHT scalar OUTPUT THREADS GENERATOR...");
        require(std::string(argv[7])=="scalar","preparation supports scalar terminal seeds only");
        opt->pool.reset(std::stoul(argv[9]));fs::path out=argv[8];fresh(out);std::vector<Mat> gs;
        for(int i=10;i<argc;++i)gs.push_back(read_matrix(argv[i]));
        save_input(out,prepare(read(argv[2]),read(argv[3]),read(argv[4]),gs,std::stoul(argv[5]),std::stoul(argv[6]),opt));
        std::cout<<"PREPARE_PASS"<<std::endl;
    }else if(mode=="toy"){
        require(argc==7,"toy LEFT_WEIGHT RIGHT_WEIGHT permutation|rational OUTPUT THREADS");
        opt->pool.reset(std::stoul(argv[6]));fs::path out=argv[5];fresh(out);
        Mat rules(3,9);size_t q=0;for(size_t i=0;i<3;++i)for(size_t j=i+1;j<3;++j){rules[q].push_back(index(3*i+j),Q(1));rules[q++].push_back(index(3*j+i),Q(-1));}
        Mat cycle(3,3),flip(3,3);for(size_t i=0;i<3;++i){cycle[i].push_back(index((i+1)%3),Q(1));flip[i].push_back(index(2-i),Q(1));}
        auto in=prepare(condition_tensor(rules,3),tensor(identity(3),{3,1,3}),tensor(identity(3),{3,3,1}),{cycle,flip},std::stoul(argv[2]),std::stoul(argv[3]),opt);
        std::string coordinates=argv[4];require(coordinates=="permutation"||coordinates=="rational","unknown toy coordinates");
        if(coordinates=="rational"){
            auto change=[&](Tensor& t,std::vector<Mat>& actions){auto h=identity(t.dim(0));if(t.dim(0)>1)h[0].push_back(1,Q(2,3));auto hi=chart_inverse(h,opt);
                t=tensor(mul(h,flatten(t)),t.dims());for(auto& a:actions)a=mul(mul(h,a),hi);};
            change(in.left,in.left_actions);change(in.right,in.right_actions);
        }
        save_input(out,in);std::cout<<"PREPARE_PASS"<<std::endl;
    }else if(mode=="run"){
        require(argc==7,"run integrability-first|joint-global|joint-staged INPUT OUTPUT THREADS raw|reduced");
        std::string method=argv[2],reduction=argv[6];require(method=="integrability-first"||method=="joint-global"||method=="joint-staged","unknown solving order");
        require(reduction=="raw"||reduction=="reduced","unknown interface reduction");opt->pool.reset(std::stoul(argv[5]));auto in=load_input(argv[3]);
        size_t nl=in.left.dim(0),nr=in.right.dim(0);auto tick=std::chrono::steady_clock::now();
        auto local=tensor_kernel::right_boundary_conditions(in.condition,in.right,reduction=="raw"?tensor_kernel::BoundaryReduction::raw:tensor_kernel::BoundaryReduction::reduced);
        tensor_kernel::ConstraintRows equations(in.left,std::move(local),nr);
        in.condition.clear();in.left.clear();in.right.clear();Mat k;size_t intermediate=0;structured_kernel::Statistics stats;
        if(method=="integrability-first"){
            k=global_solve(staged_kernel::stack(equations),opt,&stats);intermediate=k.nrow;
            std::cout<<"INTEGRABILITY_KERNEL dimension="<<k.nrow<<" nnz="<<k.nnz()<<std::endl;
            for(size_t g=0;g<in.left_actions.size();++g){auto r=restricted_symmetry(k,in.left_actions[g],in.right_actions[g],&opt->pool);auto small=kernel(std::move(r),opt);k=mul(small,k,&opt->pool);}
        }else {
            staged_kernel::ProductInvarianceRows invariants(in.left_actions,in.right_actions);auto together=staged_kernel::stack(invariants,equations);
            if(method=="joint-global")k=global_solve(together,opt,&stats);
            else k=staged_kernel::solve(together,opt,{},&stats);
        }
        double solve_seconds=seconds(tick);for(auto& row:k.rows)vec_cancel_divisor(row);
        write_tensor_view(argv[4],MatrixTensorView(k,{k.nrow,nl,nr}));
        std::cout<<"BENCH_RESULT {\"method\":\""<<method<<"\",\"left_dimension\":"<<nl<<",\"right_dimension\":"<<nr<<",\"initial_columns\":"<<nl*nr<<",\"presolved_columns\":"<<stats.core_cols
                 <<",\"integrability_dimension\":"<<intermediate<<",\"dimension\":"<<k.nrow<<",\"kernel_nnz\":"<<k.nnz()<<",\"solve_seconds\":"<<solve_seconds<<"}"<<std::endl;
    }else if(mode=="check"){
        require(argc>=5,"check INPUT THREADS RESULT...");opt->pool.reset(std::stoul(argv[3]));auto in=load_input(argv[2]);
        size_t nl=in.left.dim(0),nr=in.right.dim(0);
        // Independent original rational sewing assembler, followed by explicit
        // product actions. This is deliberately outside every timed process.
        auto reference=flatten(sew_first_last(std::move(in.condition),std::move(in.left),std::move(in.right),QQ,opt));
        size_t unrestricted=reference.nrow;
        for(size_t g=0;g<in.left_actions.size();++g){auto action=kron(in.left_actions[g],in.right_actions[g]);auto changed=mul(reference,action,&opt->pool);
            for(size_t r=0;r<reference.nrow;++r){for(auto [c,x]:reference[r])changed[r].push_back(c,-x);normalize(changed[r]);}
            auto small=kernel(transpose(changed),opt);reference=mul(small,reference,&opt->pool);
        }
        reference=row_basis(reference);
        for(int i=4;i<argc;++i){auto t=read(argv[i]);require(t.rank()==3&&t.dim(1)==nl&&t.dim(2)==nr,"result coordinate dimensions differ");auto k=flatten(t);
            require(k.nrow==reference.nrow&&equal(row_basis(k),reference),"independent exact full-space comparison failed");}
        std::cout<<"CHECK_RESULT {\"exact_space\":\"pass\",\"integrability_dimension\":"<<unrestricted<<",\"dimension\":"<<reference.nrow<<",\"outputs\":"<<argc-4<<"}"<<std::endl;
    }else require(false,"unknown command");
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
