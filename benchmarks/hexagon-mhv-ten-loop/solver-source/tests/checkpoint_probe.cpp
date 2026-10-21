#include "shared/joint_checkpoint.hpp"
#include "shared/seam_discovery.hpp"

int main(int argc,char** argv)try {
    irr::require(argc==2,"checkpoint_probe FRESH_OUTPUT");
    fs::path root=argv[1]; irr::require(!fs::exists(root),"Use a fresh checkpoint probe directory");
    fs::create_directories(root); irr_stage::fixed_text(root/"input_identity.txt","synthetic-sparse-affine-checkpoint-probe-v1\n");
    ulong p=n_nextprime(1ULL<<60,0); field_t F(FIELD_Fp,p);
    sparse_mat<ulong,int32_t> rref(3,6);
    rref[0].push_back(0,3);rref[0].push_back(2,1);rref[0].push_back(4,7);
    rref[1].push_back(1,11);rref[1].push_back(4,9);rref[1].push_back(5,1);
    std::vector<std::vector<pivot_t<int32_t>>> pivots{{{0,2},{1,5}}};
    auto expected=staged_kernel::modular_kernel(rref,pivots,F);
    auto consumed=staged_kernel::consume_modular_kernel(rref,pivots,F);
    irr::require(nmhv_sew_matrices_equal(expected,consumed)&&rref.nrow==0,"Consuming kernel preserves coordinates");
    sparse_mat<ulong,int32_t> empty(0,3);
    expected=staged_kernel::modular_kernel(empty,{},F);consumed=staged_kernel::consume_modular_kernel(empty,{},F);
    irr::require(nmhv_sew_matrices_equal(expected,consumed),"Unconstrained consuming kernel");
    sparse_mat<ulong,int32_t> a(4,5003);
    for(int32_t c=0;c<4097;++c)a[0].push_back(c,ulong(c)+1);
    a[0].push_back(5002,p-1);
    a[1].push_back(3,17); a[2].push_back(2,42); a[2].push_back(9,731);
    // The empty last row is a valid zero affine origin.
    hexagon_joint_checkpoint::save(root,a,p,8,4);
    auto b=hexagon_joint_checkpoint::load(root,p,8,4);
    irr::require(nmhv_sew_matrices_equal(a,b),"Sparse affine zero-origin roundtrip");
    bool identity_rejected=false;
    try { auto bad=hexagon_joint_checkpoint::load(root,p,10,4); }
    catch(const std::exception&) { identity_rejected=true; }
    irr::require(identity_rejected,"Sparse checkpoint must reject a different split");
    auto second=root/"nonzero_origin";fs::create_directories(second);
    irr_stage::fixed_text(second/"input_identity.txt","synthetic-sparse-affine-nonzero-origin-probe-v1\n");
    a[3].push_back(7,120); hexagon_joint_checkpoint::save(second,a,p,8,4);
    b=hexagon_joint_checkpoint::load(second,p,8,4);
    irr::require(nmhv_sew_matrices_equal(a,b),"Sparse affine nonzero-origin roundtrip");
    auto path=hexagon_joint_checkpoint::file(second,p);
    std::fstream damaged(path,std::ios::binary|std::ios::in|std::ios::out);
    damaged.seekp(sizeof(hexagon_joint_checkpoint::Header)+sizeof(uint64_t));
    char byte=127;damaged.write(&byte,1);damaged.close();
    bool damage_rejected=false;
    try { auto bad=hexagon_joint_checkpoint::load(second,p,8,4); }
    catch(const std::exception&) { damage_rejected=true; }
    irr::require(damage_rejected,"Sparse checkpoint must reject changed bytes");
    std::cout<<"CHECKPOINT_PROBE PASS zero_origin=1 nonzero_origin=1 identity_rejection=1 damage_rejection=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
