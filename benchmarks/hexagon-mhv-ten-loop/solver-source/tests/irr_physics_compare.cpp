#include "irr_object.hpp"
#include "irr_runtime.hpp"
int main(int argc,char** argv)try{
    irr::require(argc==8,"compare H K OBJECT FEC LEC REFERENCE_WXF exact|space");Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->verbose=false;irr_run::threads(opt);
    auto obj=irr_io::read_object(argv[3],std::stoul(argv[1]),std::stoul(argv[2]),argv[4],argv[5],&opt->pool);
    auto tensor=irr::read_tensor(argv[6],&opt->pool);auto reference=joint_eet_flatten_tensor(tensor,&opt->pool);
    sparse_mat<rat_t,int32_t> packed(obj.rows(),reference.ncol);size_t off=0;
    for(auto& m:obj.S){for(size_t i=0;i<m.nrow;++i){size_t row=i==obj.directions?0:i+1;for(auto [j,x]:m[i])packed[row].push_back(off+j,x);}off+=m.ncol;}
    irr::require(off==reference.ncol&&packed.nrow==reference.nrow,"Shape mismatch");const field_t F(FIELD_QQ);
    if(std::string(argv[7])=="space"){
        // Homogenize affine spaces: origin gets a leading 1, directions a 0.
        auto canonical=[&](sparse_mat<rat_t,int32_t> m){sparse_mat<rat_t,int32_t> a(m.nrow,m.ncol+1);a[0].push_back(0,rat_t(1));for(size_t i=0;i<m.nrow;++i)for(auto [j,x]:m[i])a[i].push_back(j+1,x);
            auto piv=sparse_mat_rref(a,F,opt);std::vector<std::pair<int32_t,int32_t>> order;for(auto& g:piv)for(auto p:g)order.emplace_back(p.c,p.r);std::sort(order.begin(),order.end());
            sparse_mat<rat_t,int32_t> out(order.size(),a.ncol);for(size_t i=0;i<order.size();++i)out[i]=a[order[i].second];return out;};
        packed=canonical(std::move(packed));reference=canonical(std::move(reference));
    }
    irr::require(nmhv_sew_matrices_equal(packed,reference),"Exact comparison failed");std::cout<<"PHYSICS COMPARISON PASS "<<argv[3]<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
