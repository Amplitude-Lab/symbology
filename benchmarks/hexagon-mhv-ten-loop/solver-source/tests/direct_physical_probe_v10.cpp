#include "shared/direct_physical_rows.hpp"
using irr::require;
using PM=irr_physics::mat;

int main(int argc,char**argv)try {
    require(argc==8,"direct_probe H K FEC LEC DATA DIAGNOSTICS prime");
    Flint::set_memory_functions();
    irr_physics::context c(argv[3],argv[4],argv[6],argv[5],std::stoul(argv[1]),std::stoul(argv[2]),2,true,std::stoull(argv[7]));
    // Nonzero offset checks that the direct path keeps the same sample stream.
    c.sample_offset=5;irr_physics::direct_coordinate_space(c);
    size_t n=irr_physics::invariant_columns(c);constexpr size_t affine_rows=4;
    std::array<PM,6> candidate;
    size_t candidate_nnz=0,origin_nnz=0;
    for(size_t b=0;b<6;++b){size_t width=c.sew.first[b]*c.sew.last[b];candidate[b]=PM(affine_rows,width);
        for(size_t r=0;r<affine_rows;++r)for(size_t j=0;j<width;++j){
            // Reproducible sparse directions and a nonzero final-row origin.
            // Include column zero in every nonempty sector, including odd/E.
            if(j&&((j*37+r*13+b*7)%17))continue;
            ulong x=(j*101+r*59+b*23+7)%c.F.mod.n;if(!x)x=1;
            if((j+r+b)&1)x=nmod_neg(x,c.F.mod);
            candidate[b][r].push_back(j,x);++candidate_nnz;if(r+1==affine_rows)++origin_nnz;
        }}
    require(origin_nnz,"Probe needs a nonzero affine origin");
    auto identity=[&]{c.rows=n+1;size_t off=0;for(size_t b=0;b<6;++b){size_t width=c.sew.first[b]*c.sew.last[b];
        c.S[b]=PM(n+1,width);for(size_t j=0;j<width;++j)c.S[b][off+j].push_back(j,1);off+=width;}};
    for(auto name:{"scalar_collinear","scalar_mrk","tail_exact"}){
        bool mrk=std::string(name)=="scalar_mrk",tail=std::string(name)=="tail_exact";
        irr_physics::direct_coordinate_space(c);
        auto dual=tail?irr_physics::direct_collinear_tail_exact(c):irr_physics::direct_slots(c,mrk);
        for(const auto& row:dual)require(row.coefficients.size()==n,"Direct row unexpectedly contains an origin column");
        // Exhaustive small-control comparison with the old identity projection.
        // Identity is confined to this probe; the production helper never builds it.
        identity();auto full=tail?irr_physics::collinear_tail_exact(c):irr_physics::slots(c,mrk);
        require(dual.size()==full.size(),"Direct/identity row count mismatch");
        for(size_t i=0;i<dual.size();++i){
            require(full[i].coefficients.size()==n+1&&!full[i].coefficients.back(),"Identity probe origin mismatch");
            require(std::equal(dual[i].coefficients.begin(),dual[i].coefficients.end(),full[i].coefficients.begin())&&dual[i].rhs==full[i].rhs,"Direct projection differs from identity projection");}
        full.clear();c.rows=affine_rows;c.S=candidate;
        auto primal=tail?irr_physics::collinear_tail_exact(c):irr_physics::slots(c,mrk);
        require(dual.size()==primal.size(),"Direct/candidate row count mismatch");
        size_t nnz=0,min_nnz=SIZE_MAX,max_nnz=0,nonzero_rhs=0,nonzero_origin_images=0;
        for(size_t i=0;i<dual.size();++i){
            require(primal[i].coefficients.size()==affine_rows&&dual[i].rhs==primal[i].rhs,"Candidate RHS/affine width mismatch");
            std::array<ulong,affine_rows> expected{};size_t offset=0;
            for(size_t b=0;b<6;++b){for(size_t r=0;r<affine_rows;++r)for(auto [j,x]:candidate[b][r])
                expected[r]=nmod_add(expected[r],nmod_mul(dual[i].coefficients[offset+j],x,c.F.mod),c.F.mod);
                offset+=candidate[b].ncol;}
            require(std::equal(expected.begin(),expected.end(),primal[i].coefficients.begin()),"Direct rows times affine candidate differ from old projection");
            // The origin is A*o; the affine solve RHS is b-A*o.
            require(nmod_sub(dual[i].rhs,expected.back(),c.F.mod)==nmod_sub(primal[i].rhs,primal[i].coefficients.back(),c.F.mod),"Affine translated RHS mismatch");
            size_t count=0;for(auto x:dual[i].coefficients)count+=bool(x);nnz+=count;min_nnz=std::min(min_nnz,count);max_nnz=std::max(max_nnz,count);
            nonzero_rhs+=bool(dual[i].rhs);nonzero_origin_images+=bool(expected.back());
        }
        std::cout<<"DIRECT_PROBE condition="<<name<<" passed=1 identity_comparison=1 affine_candidate_comparison=1 rows="<<dual.size()<<" columns="<<n
            <<" nnz="<<nnz<<" min_row_nnz="<<min_nnz<<" max_row_nnz="<<max_nnz<<" nonzero_rhs="<<nonzero_rhs
            <<" candidate_rows="<<affine_rows<<" candidate_nnz="<<candidate_nnz<<" origin_nnz="<<origin_nnz
            <<" nonzero_origin_images="<<nonzero_origin_images<<" sample_offset="<<c.sample_offset<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<"DIRECT_PROBE_ERROR "<<e.what()<<std::endl;return 1;}
