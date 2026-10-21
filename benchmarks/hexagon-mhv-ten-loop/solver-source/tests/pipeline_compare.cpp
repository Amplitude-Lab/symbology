#include "irr_physics_pipeline.hpp"

template<typename T>
void compare(sparse_mat<T,int32_t> a,sparse_mat<T,int32_t> b,const field_t& F,rref_option_t& opt,const T& minus_one){
    require(a.nrow==b.nrow&&a.ncol==b.ncol,"Affine dimensions differ");size_t n=a.nrow-1;
    auto rank=[&](sparse_mat<T,int32_t> m){auto pivots=sparse_mat_rref(m,F,opt);return nmhv_sew_modular_rank(pivots);};
    sparse_mat<T,int32_t> directions(2*n,a.ncol);
    for(size_t i=0;i<n;++i){directions[i]=a[i];directions[n+i]=b[i];}
    require(rank(directions)==n,"Direction spaces differ or are not independent");
    sparse_mat<T,int32_t> difference(1,a.ncol);difference[0]=a[n];sparse_vec_add_mul(difference[0],b[n],minus_one,F);difference[0].canonicalize();
    directions.append(difference);require(rank(std::move(directions))==n,"Origins differ outside the direction space");
    std::cout<<"PASS affine spaces directions="<<n<<" columns="<<a.ncol<<std::endl;
}

int main(int argc,char** argv)try {
    require(argc==8,"pipeline_compare A B H K FEC LEC PRIME|QQ");Flint::set_memory_functions();
    size_t h=irr_io::number(argv[3]),k=irr_io::number(argv[4]);
    rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    if(std::string(argv[7])=="QQ"){
        auto load=[&](const fs::path& root){auto obj=irr_io::read_object(root,h,k,argv[5],argv[6],&opt->pool);
            sparse_mat<rat_t,int32_t> out(obj.directions+1,irr_shift_cache::width(obj.first,obj.last));size_t off=0;
            for(const auto& block:obj.S){for(size_t i=0;i<block.nrow;++i)for(auto [j,x]:block[i])out[i].push_back(off+j,x);off+=block.ncol;}return out;};
        compare(load(argv[1]),load(argv[2]),field_t(FIELD_QQ),opt,rat_t(-1));return 0;
    }
    ulong p=std::stoull(argv[7]);field_t F(FIELD_Fp,p);
    auto load=[&](const fs::path& root){
        if(fs::exists(cache_file(root,p)))return read_cache(root,p,h,k);
        auto obj=irr_io::read_object(root,h,k,argv[5],argv[6],&opt->pool);
        sparse_mat<ulong,int32_t> out(obj.directions+1,irr_shift_cache::width(obj.first,obj.last));size_t off=0;
        for(const auto& block:obj.S){auto m=irr_shift_cache::reduce(block,F);for(size_t i=0;i<m.nrow;++i)for(auto [j,x]:m[i])out[i].push_back(off+j,x);off+=m.ncol;}return out;
    };
    compare(load(argv[1]),load(argv[2]),F,opt,nmod_neg(1,F.mod));
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
