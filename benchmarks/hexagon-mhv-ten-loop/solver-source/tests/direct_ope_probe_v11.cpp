#include "shared/direct_ope_rows.hpp"
using irr::require;
using PM=irr_physics::mat;

static ulong random_residue(uint64_t& state,ulong p){state^=state<<13;state^=state>>7;state^=state<<17;return state%p;}
static ulong dot(const std::vector<ulong>& a,const std::vector<ulong>& b,const field_t& F){
    require(a.size()==b.size(),"Probe dot dimensions");ulong value=0;
    for(size_t i=0;i<a.size();++i)value=nmod_add(value,nmod_mul(a[i],b[i],F.mod),F.mod);return value;
}

int main(int argc,char**argv)try {
    require(argc==8,"direct_ope_probe H K FEC LEC DATA DIAGNOSTICS prime");Flint::set_memory_functions();
    irr_physics::context c(argv[3],argv[4],argv[6],argv[5],std::stoul(argv[1]),std::stoul(argv[2]),1,true,std::stoull(argv[7]));
    size_t weight=c.h+c.k;require(weight>=4,"Low-grade OPE probe needs weight at least 4");uint64_t state=0x65b94de3712fca08ULL;
    size_t duality_cases=0;
    for(size_t f=2;f<=4;++f){auto kernel=c.data/"ope_mhv_small_s"/("kernels_F"+std::to_string(f)+".tsv");
        nmhv_ope_germ_t germ(weight,c.F.mod.n,kernel,2,4);irr_physics::require_ope_grade_filtration(germ);
        std::vector<ulong> previous(germ.stride),cotangent(germ.stride),image(germ.stride),adjoint(germ.stride);
        // Arbitrary vectors include unreachable log slots, testing the exact
        // implemented linear map rather than only valid physical input jets.
        for(size_t w=1;w<=weight;++w)for(size_t letter=0;letter<9;++letter){
            for(auto& x:previous)x=random_residue(state,c.F.mod.n);for(auto& x:cotangent)x=random_residue(state,c.F.mod.n);
            germ.integrate(previous.data(),image.data(),letter,w);
            irr_physics::integrate_ope_adjoint(germ,cotangent.data(),adjoint.data(),letter,w);
            require(dot(image,cotangent,c.F)==dot(previous,adjoint,c.F),"OPE integrate primal/adjoint duality failed");++duality_cases;
        }
    }
    std::cout<<"DIRECT_OPE_ADJOINT_PROBE passed=1 cases="<<duality_cases<<" kernels=3 max_t=2 grade=4 max_weight="<<weight<<std::endl;
    constexpr size_t affine_rows=4;std::array<PM,6> candidate;size_t columns=0,origin_nnz=0,candidate_nnz=0;
    for(size_t b=0;b<6;++b){size_t width=c.sew.first[b]*c.sew.last[b];columns+=width;candidate[b]=PM(affine_rows,width);
        for(size_t r=0;r<affine_rows;++r)for(size_t j=0;j<width;++j){if(j&&((j*31+r*11+b*7)%13))continue;
            ulong value=(j*79+r*37+b*23+9)%c.F.mod.n;if(!value)value=1;if((r+j+b)&1)value=nmod_neg(value,c.F.mod);
            candidate[b][r].push_back(j,value);++candidate_nnz;if(r+1==affine_rows)++origin_nnz;
        }}
    require(origin_nnz,"OPE probe needs a nonzero affine origin");
    auto direct=irr_physics::direct_ope_small_s_lowgrade(c,true);
    c.rows=affine_rows;c.S=candidate;auto legacy=irr_physics::ope_small_s(c,false,true);
    std::map<std::array<int,5>,size_t> by_coordinate;
    for(size_t i=0;i<legacy.size();++i)require(by_coordinate.emplace(legacy[i].ope_coordinate,i).second,"Duplicate legacy OPE coordinate");
    size_t matched=0,nonzero_rhs=0,nonzero_origin_images=0,direct_nnz=0;
    for(const auto& row:direct){require(row.coefficients.size()==columns,"Direct OPE row contains an origin column");
        auto it=by_coordinate.find(row.ope_coordinate);require(it!=by_coordinate.end(),"Low-grade direct row absent from existing OPE solve rows");
        const auto& old=legacy[it->second];require(row.rhs==old.rhs&&old.coefficients.size()==affine_rows,"Direct OPE RHS/legacy affine dimensions differ");
        std::array<ulong,affine_rows> expected{};size_t offset=0;
        for(size_t b=0;b<6;++b){for(size_t r=0;r<affine_rows;++r)for(auto [j,x]:candidate[b][r])
            expected[r]=nmod_add(expected[r],nmod_mul(row.coefficients[offset+j],x,c.F.mod),c.F.mod);offset+=candidate[b].ncol;}
        require(std::equal(expected.begin(),expected.end(),old.coefficients.begin()),"Direct low-grade OPE times affine candidate differs from original grade-12 projection");
        require(nmod_sub(row.rhs,expected.back(),c.F.mod)==nmod_sub(old.rhs,old.coefficients.back(),c.F.mod),"Direct OPE translated affine RHS differs");
        ++matched;nonzero_rhs+=bool(row.rhs);nonzero_origin_images+=bool(expected.back());for(auto x:row.coefficients)direct_nnz+=bool(x);
    }
    require(matched==direct.size()&&matched,"No direct/legacy OPE comparisons");
    std::cout<<"DIRECT_OPE_PROBE passed=1 rows="<<matched<<" original_columns="<<columns<<" direct_nnz="<<direct_nnz
        <<" candidate_rows="<<affine_rows<<" candidate_nnz="<<candidate_nnz<<" origin_nnz="<<origin_nnz
        <<" nonzero_origin_images="<<nonzero_origin_images<<" nonzero_rhs="<<nonzero_rhs
        <<" rhs_source=existing_solve_rows full_grade12_projection_comparison=1 heldouts_used_for_fit=0"<<std::endl;
}catch(const std::exception& e){std::cerr<<"DIRECT_OPE_PROBE_ERROR "<<e.what()<<std::endl;return 1;}
