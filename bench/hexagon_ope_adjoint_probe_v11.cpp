#include "hexagon_ope_adjoint_v11.hpp"
using irr::require;
ulong dot(const std::vector<ulong>& a,const std::vector<ulong>& b,const field_t& F){
    require(a.size()==b.size(),"OPE probe dot shape");ulong x=0;
    for(size_t i=0;i<a.size();++i)x=nmod_add(x,nmod_mul(a[i],b[i],F.mod),F.mod);return x;
}
int main(int argc,char**argv)try{
    require(argc==8,"ope_adjoint_probe H K FEC LEC DATA DIAGNOSTICS prime");Flint::set_memory_functions();
    size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]),weight=h+k;ulong prime=std::stoull(argv[7]);
    auto kernel=std::filesystem::path(argv[5])/"ope_mhv_small_s/kernels_F2.tsv";
    for(ulong p:{ulong(101),prime}){nmhv_ope_germ_t germ(weight,p,kernel,2,4);size_t cases=0;
        for(size_t w:{size_t(1),std::min(size_t(3),weight),weight})for(size_t letter=0;letter<9;++letter){
            std::vector<ulong> previous(germ.stride),lambda(germ.stride),value(germ.stride),adjoint(germ.stride);
            for(int t=0;t<=2;++t)for(int z=-t;z<=4-t;++z)
                for(size_t a=0;a<=weight;++a)for(size_t b=0;b<=weight;++b){
                    size_t j=germ.index(t,z,a,b);
                    // Arbitrary dense source and output covector also exercise
                    // zeroed coordinates beyond each step's log-degree cutoff.
                    previous[j]=(j*101+letter*13+w*5+7)%p;
                    lambda[j]=(j*53+letter*29+w*17+11)%p;
                }
            germ.integrate(previous.data(),value.data(),letter,w);
            hexagon_ope_adjoint::integrate(germ,lambda.data(),adjoint.data(),letter,w);
            require(dot(value,lambda,germ.F)==dot(previous,adjoint,germ.F),"Ordered germ transpose duality failed");++cases;
        }
        std::cout<<"OPE_ADJOINT_PROBE prime="<<p<<" transpose_cases="<<cases<<" passed=1"<<std::endl;
    }
    irr_physics::context c(argv[3],argv[4],argv[6],argv[5],h,k,1,true,prime);
    auto dual=hexagon_ope_adjoint::selected_rows(c,4,32,4);
    const size_t n=irr_physics::invariant_columns(c),candidate_rows=4;
    std::array<irr_physics::mat,6> candidate;
    for(size_t b=0;b<6;++b){size_t width=c.sew.first[b]*c.sew.last[b];candidate[b]=irr_physics::mat(candidate_rows,width);
        for(size_t r=0;r<candidate_rows;++r)for(size_t j=0;j<width;++j)
            if(!j||!((j*37+r*13+b*7)%17))candidate[b][r].push_back(j,(j*101+r*59+b*23+7)%prime);
    }
    c.rows=candidate_rows;c.S=candidate;
    auto primal=irr_physics::ope_small_s(c,false,true);
    size_t checked=0;
    for(const auto& row:dual){require(row.coefficients.size()==n,"Direct OPE row width");
        auto found=std::find_if(primal.begin(),primal.end(),[&](const auto& old){return old.ope_coordinate==row.ope_coordinate;});
        require(found!=primal.end()&&found->rhs==row.rhs,"Direct OPE target/RHS differs from full-grade original projection");
        std::array<ulong,candidate_rows> expected{};size_t offset=0;
        for(size_t b=0;b<6;++b){for(size_t r=0;r<candidate_rows;++r)for(auto [j,x]:candidate[b][r])
            expected[r]=nmod_add(expected[r],nmod_mul(row.coefficients[offset+j],x,c.F.mod),c.F.mod);
            offset+=candidate[b].ncol;}
        require(std::equal(expected.begin(),expected.end(),found->coefficients.begin()),"Adjoint OPE row times affine candidate differs from full-grade projection");
        require(nmod_sub(row.rhs,expected.back(),c.F.mod)==nmod_sub(found->rhs,found->coefficients.back(),c.F.mod),"Direct OPE translated RHS differs");++checked;
    }
    std::cout<<"OPE_ADJOINT_PROBE affine_candidate_rows="<<candidate_rows<<" targets="<<checked<<" columns="<<n
             <<" lower_grade_matches_original=1 nonzero_origin_checked=1 passed=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"OPE_ADJOINT_PROBE_ERROR "<<e.what()<<std::endl;return 1;}
