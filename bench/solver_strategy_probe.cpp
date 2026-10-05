#include "symrep_bootstrap.hpp"
#include "structured_kernel.hpp"
#include "tensor_kernel.hpp"
#include <random>
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc>=3,"strategy_probe MODE ...");rref_option_t opt;opt->pool.reset(8);
    std::string mode=argv[1];
    if(mode=="stream-sew"){
        require(argc==6,"strategy_probe stream-sew CONDITION FEC LEC OUTPUT");
        auto d=read(argv[2]),f=read(argv[3]),l=read(argv[4]);size_t previous=f.dim(0),last=l.dim(0);
        auto local=tensor_kernel::right_boundary_conditions(d,l);
        std::cout<<"LOCAL rows="<<local.nrow<<" columns="<<local.ncol<<" nnz="<<local.nnz()<<std::endl;
        tensor_kernel::ConstraintRows rows(f,std::move(local),last);f.clear();
        structured_kernel::Statistics stats;auto k=structured_kernel::solve_stream(rows.ncol(),rows,opt,&stats);
        std::cout<<"KERNEL rows="<<k.nrow<<" nnz="<<k.nnz()<<" core_cols="<<stats.core_cols<<" core_nnz="<<stats.core_nnz<<" presolve_s="<<stats.presolve_s<<" solve_s="<<stats.solve_s<<std::endl;
        for(auto& row:k.rows)vec_cancel_divisor(row);rows.verify(k,opt);
        write_tensor_view(argv[5],MatrixTensorView(k,{k.nrow,previous,last}));
    }else if(mode=="extension"){
        require(argc==7,"strategy_probe extension CONDITION SEED OUTPUT forward|backward normal|presolve");
        auto a=assemble(read(argv[3]),read(argv[2]),std::string(argv[5])=="backward");
        structured_kernel::Statistics stats;std::cout<<"MATRIX rows="<<a.nrow<<" cols="<<a.ncol<<" nnz="<<a.nnz()<<std::endl;
        auto k=std::string(argv[6])=="presolve"?structured_kernel::solve(std::move(a),opt,&stats):kernel(std::move(a),opt);
        std::cout<<"KERNEL rows="<<k.nrow<<" nnz="<<k.nnz()<<" core_cols="<<stats.core_cols<<" core_nnz="<<stats.core_nnz<<" presolve_s="<<stats.presolve_s<<" solve_s="<<stats.solve_s<<std::endl;
        write(argv[4],k);
    }else if(mode=="sew-one"){
        require(argc==7,"strategy_probe sew-one CONDITION FEC LEC OUTPUT normal|presolve");
        auto d=read(argv[2]),f=read(argv[3]),l=read(argv[4]);
        require(l.dim(2)==1,"sew-one requires a terminal one-dimensional right boundary");
        const size_t letters=d.dim(0),q=d.dim(2),last=l.dim(0),older=f.dim(1);
        auto lc=flatten(l).transpose();Mat local(q,letters*last);
        sparse_tensor<Q,I,SPARSE_COO> dc(d);
        for(size_t p=0;p<dc.nnz();++p){auto i=dc.index(p);for(auto [j,x]:lc[i[1]])local[i[2]].push_back(index(i[0]*last+j),dc.val(p)*x);}
        for(auto& row:local.rows)normalize(row);
        local=row_basis(local);size_t nr=local.nrow;
        std::cout<<"LOCAL rows="<<nr<<" columns="<<local.ncol<<" nnz="<<local.nnz()<<std::endl;
        std::vector<std::vector<std::tuple<size_t,size_t,Q>>> routes(letters);
        for(size_t i=0;i<nr;++i)for(auto [j,x]:local[i])routes[j/last].emplace_back(i,j%last,x);
        Mat a(older*nr,f.dim(0)*last);
        for(size_t b=0;b<f.dim(0);++b)for(size_t p=f.rowptr()[b];p<f.rowptr()[b+1];++p){auto i=f.index(p);
            for(const auto& [r,j,x]:routes[i[1]])a[i[0]*nr+r].push_back(index(b*last+j),f.val(p)*x);}
        for(auto& row:a.rows)normalize(row);
        std::map<size_t,size_t> histogram;for(const auto& row:a.rows)++histogram[row.nnz()];
        std::cout<<"MATRIX rows="<<a.nrow<<" cols="<<a.ncol<<" nnz="<<a.nnz()<<" one="<<histogram[1]<<" two="<<histogram[2]<<std::endl;
        structured_kernel::Statistics stats;
        auto k=std::string(argv[6])=="presolve"?structured_kernel::solve(std::move(a),opt,&stats):kernel(std::move(a),opt);
        std::cout<<"KERNEL rows="<<k.nrow<<" nnz="<<k.nnz()<<" core_cols="<<stats.core_cols<<" core_nnz="<<stats.core_nnz<<" presolve_s="<<stats.presolve_s<<" solve_s="<<stats.solve_s<<std::endl;
        write_tensor_view(argv[5],MatrixTensorView(k,{k.nrow,f.dim(0),last}));
    }else throw std::runtime_error("unknown mode");
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
