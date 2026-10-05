// General matrix workloads: no polygon dimensions or group-specific choices.
#include "structured_kernel.hpp"
#include <random>
#include <sys/resource.h>
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==4,"kernel_benchmark CASE baseline|structured THREADS");std::string name=argv[1];
    rref_option_t opt;opt->pool.reset(std::stoul(argv[3]));std::mt19937 rng(271828);
    Mat a;size_t nullity=0;
    if(name=="chain"||name=="graph"){
        size_t n=12000;a=Mat(name=="chain"?n-1:2*n-1,n);nullity=1;
        for(size_t i=0;i<n-1;++i){a[i].push_back(index(i),Q(1));a[i].push_back(index(i+1),Q(-1));}
        if(name=="graph")for(size_t i=n-1;i<a.nrow;++i){size_t u=rng()%n,v=rng()%n;a[i].push_back(index(u),Q(1));a[i].push_back(index(v),Q(-1));normalize(a[i]);}
    }else if(name=="zero-propagation"){
        size_t n=12000;a=Mat(n,n+8);nullity=8;a[0].push_back(0,Q(1));
        for(size_t i=1;i<n;++i){a[i].push_back(index(i-1),Q(2));a[i].push_back(index(i),Q(3));}
    }else if(name=="rational"||name=="wide"||name=="dense-tail"||name=="blocks"){
        size_t r=name=="wide"?600:2000;nullity=name=="wide"?2400:name=="dense-tail"?100:40;
        a=Mat(2*r,r+nullity);Mat base(r,r+nullity);
        for(size_t i=0;i<r;++i){base[i].push_back(index(i),Q(1));size_t terms=name=="dense-tail"?70:5;
            for(size_t j=0;j<terms;++j){size_t col=r+rng()%nullity;long val=long(rng()%7)-3;if(val)base[i].push_back(index(col),Q(val)/Q(name=="rational"?1+rng()%7:1));}normalize(base[i]);a[i]=base[i];}
        for(size_t i=0;i<r;++i){a[r+i]=base[i];sparse_vec_sub_mul(a[r+i],base[name=="blocks"?(i/20*20+(i+1)%20):(i+1)%r],Q(-1),QQ);}
    }else throw std::runtime_error("unknown case");
    size_t rows=a.nrow,cols=a.ncol,nnz=a.nnz();auto original=a;
    auto start=std::chrono::steady_clock::now();structured_kernel::Statistics stats;
    auto k=std::string(argv[2])=="structured"?structured_kernel::solve(std::move(a),opt,&stats):kernel(std::move(a),opt);
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    require(k.nrow==nullity,"unexpected nullity");require(!mul(original,k.transpose(),&opt->pool).nnz(),"nonzero exact residual");
    // This test certificate is independent of elimination: the original
    // contains an identity pivot block, or a spanning tree, so its rank is known.
    std::vector<I> owner(k.ncol,-1);std::vector<bool> private_coordinate(k.nrow);
    for(size_t i=0;i<k.nrow;++i)for(auto [j,x]:k[i])owner[j]=owner[j]==-1?index(i):I(-2);
    for(auto i:owner)if(i>=0)private_coordinate[i]=true;
    require(std::all_of(private_coordinate.begin(),private_coordinate.end(),[](bool x){return x;}),"kernel lacks independence certificate");
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"RESULT case="<<name<<" method="<<argv[2]<<" rows="<<rows<<" cols="<<cols<<" nnz="<<nnz<<" nullity="<<k.nrow<<" kernel_nnz="<<k.nnz()<<" seconds="<<seconds<<" rss_kib="<<usage.ru_maxrss<<" core_cols="<<stats.core_cols<<" core_nnz="<<stats.core_nnz<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
