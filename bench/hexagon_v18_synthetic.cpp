#define main frozen_v17_main
#include "joint_physics_v17.cpp"
#undef main
#include "shared/hexagon_final_seam_v18.hpp"

struct ToyRows {
    MM matrix;field_t F;
    size_t nrow()const{return matrix.nrow;}size_t ncol()const{return matrix.ncol;}
    size_t residual_blocks()const{return (nrow()+127)/128;}
    sparse_vec<ulong,int32_t> residual_direct(size_t id,const MM& kt,size_t count)const{
        std::vector<ulong> v(count);for(auto[c,x]:matrix[id])for(auto[j,y]:kt[c])v[j]=nmod_add(v[j],nmod_mul(x,y,F.mod),F.mod);
        sparse_vec<ulong,int32_t> row;for(size_t j=0;j<count;++j)if(v[j])row.push_back(j,v[j]);return row;
    }
    template<class Emit>void selected_residual_transpose(const MM& kt,size_t count,const std::vector<size_t>& ids,Emit&& emit,thread_pool*)const{
        for(auto id:ids)emit(id,residual_direct(id,kt,count));
    }
    template<class Emit>void residual_transpose_direct(const MM& kt,size_t count,size_t b,Emit&& emit,thread_pool*)const{
        for(size_t id=b*128;id<std::min(nrow(),(b+1)*128);++id)emit(residual_direct(id,kt,count));
    }
};
int main(int argc,char**argv)try{
    require(argc==2,"v18_synthetic OUTPUT");fs::path root=argv[1];require(!fs::exists(root),"Use fresh synthetic output");
    fs::create_directories(root);rref_option_t opt;opt->pool.reset(2);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    size_t cases=0;
    for(size_t n:{2,9,70,145})for(ulong prime:{1000003UL,1000033UL}){
        field_t F(FIELD_Fp,prime);ToyRows rows{MM(3*(n-1),n),F};
        // Nonzero affine point x_j=(j+2)/(j+1), x_0=1, including redundant rows.
        for(size_t j=1;j<n;++j)for(size_t copy=0;copy<3;++copy){
            auto& row=rows.matrix[3*(j-1)+copy];row.push_back(0,nmod_neg(nmod_div(j+2,j+1,F.mod),F.mod));row.push_back(j,1);}
        MM initial(n,n);for(size_t j=1;j<n;++j)initial[j-1].push_back(j,1);initial[n-1].push_back(0,1);
        auto result=hexagon_final_seam_v18::finish(rows,std::move(initial),F,opt,root/("n"+std::to_string(n)),"toy_"+std::to_string(n));
        require(result.nrow==1&&result[0].nnz()==n,"Toy unique shape");
        for(size_t j=0;j<n;++j){auto value=result[0].find(j);require(value&&*value==(j?nmod_div(j+2,j+1,F.mod):1),"Toy exact rational reference");}
        ++cases;
    }
    size_t rejected=0;
    {
        field_t F(FIELD_Fp,1000033);ToyRows lost{MM(24,9),F};MM initial(9,9);
        for(size_t j=1;j<9;++j)initial[j-1].push_back(j,1);initial[8].push_back(0,1);
        try{hexagon_final_seam_v18::finish(lost,std::move(initial),F,opt,root/"n9","toy_9");}
        catch(const std::exception&e){require(std::string(e.what()).find("replay rank changed")!=std::string::npos,"Wrong lost-rank failure");++rejected;}
    }
    for(ulong prime:{1000003UL,1000033UL}){
        field_t F(FIELD_Fp,prime);ToyRows rows{MM(5,9),F};
        for(size_t j=1;j<=5;++j){rows.matrix[j-1].push_back(0,nmod_neg(nmod_div(j+2,j+1,F.mod),F.mod));rows.matrix[j-1].push_back(j,1);}
        MM initial(9,9);for(size_t j=1;j<9;++j)initial[j-1].push_back(j,1);initial[8].push_back(0,1);
        auto result=hexagon_final_seam_v18::finish(rows,std::move(initial),F,opt,root/"affine_family","family");
        require(result.nrow==4,"Lost affine directions");
        for(size_t j=0;j<=5;++j){auto x=result[3].find(j);require(x&&*x==(j?nmod_div(j+2,j+1,F.mod):1),"Affine rational origin");}
        for(size_t r=0;r<3;++r)require(result[r].nnz()==1&&result[r](0)==r+6&&result[r][0]==1,"Affine free direction");++cases;
    }
    {
        field_t F(FIELD_Fp,1000003);ToyRows impossible{MM(1,1),F};impossible.matrix[0].push_back(0,1);MM point(1,1);point[0].push_back(0,1);
        try{hexagon_final_seam_v18::finish(impossible,std::move(point),F,opt,root/"inconsistent","impossible");}
        catch(const std::exception&){++rejected;}
    }
    auto file=root/"n9/candidate_9.tsv";hexagon_final_seam_v18::Plan plan;
    try{hexagon_final_seam_v18::load(file,"wrong_identity",plan);}catch(const std::exception&){++rejected;}
    {std::ofstream out(file,std::ios::app);out<<"corrupt\n";}
    try{hexagon_final_seam_v18::load(file,"toy_9",plan);}catch(const std::exception&){++rejected;}
    require(rejected==4,"Plan/origin/rank guards");
    std::cout<<"V18_SYNTHETIC_PASS cases="<<cases<<" exact_rational_origins=1 plan_guards="<<rejected<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
