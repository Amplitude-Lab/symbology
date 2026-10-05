// Compile with the unpacked reference package's src and src/reduced_d3 includes.
// Exports its direct-multiplicity chain to ordinary recurrence tensors for an
// independent, exact change-of-basis comparison with symrep.
#include "irr_recurrence.hpp"
#include "joint_eet_lec.hpp"
#include <iostream>
namespace fs=std::filesystem;
using Q=rat_t;using I=int32_t;
using Mat=sparse_mat<Q,I>;
using Tensor=sparse_tensor<Q,I,SPARSE_CSR>;
void write(const fs::path& p,const Tensor& t) {
    if(fs::exists(p))throw std::runtime_error("Refusing to overwrite "+p.string());
    auto bytes=sparse_tensor_write_wxf(t);std::ofstream out(p,std::ios::binary);
    out.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());out.close();
    if(!out)throw std::runtime_error("Write failed");
}
int main(int argc,char**argv){try {
    if(argc!=6)throw std::runtime_error("hexagon_symrep_reference FEC|LEC CHAIN MAX_WEIGHT DATA OUTPUT");
    bool backward=std::string(argv[1])=="LEC";fs::path root=argv[2],data=argv[4],out=argv[5];size_t max=std::stoul(argv[3]);
    if(fs::exists(out))throw std::runtime_error("Output already exists");fs::create_directories(out);
    rref_option_t opt;opt->pool.reset(2);auto* pool=&opt->pool;irr::channels_t cg(pool);irr::copies previous{1,0,0,0,0,0};
    // NMHV starts from a five-dimensional R-invariant terminal representation.
    // Read its multiplicities rather than assuming a scalar MHV terminal.
    if(backward&&fs::exists(root/"LEC_0_blocks.tsv")) {
        auto blocks=nmhv_sew_read_irrep_axis(root/"LEC_0_blocks.tsv");
        for(size_t b=0;b<6;++b)previous[b]=blocks[b].multiplicity;
    }
    for(size_t w=1;w<=max;++w) {
        std::ostringstream dn;dn<<'w'<<std::setw(2)<<std::setfill('0')<<w;irr::matrices K;irr::copies current{};size_t n=0,old=0;
        for(size_t b=0;b<6;++b){auto t=irr::read_tensor(root/dn.str()/((backward?"irrLEC_":"irrFEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"),pool);
            sparse_tensor<Q,I,SPARSE_COO> c(std::move(t));K[b]=c.to_sparse_mat(pool);current[b]=K[b].nrow;n+=current[b]*irr::specs()[b].irrep_dimension;old+=previous[b]*irr::specs()[b].irrep_dimension;}
        Mat all(n,old*9);size_t off=0;
        for(size_t b=0;b<6;++b){auto ex=irr::expand(cg,b,previous,K[b]);for(size_t r=0;r<ex.nrow;++r)all[off++]=ex[r];}
        sparse_tensor<Q,I,SPARSE_COO> c(all,pool);c.reshape({n,old,9});if(backward)c.transpose_replace({0,2,1});
        if(backward&&w==1&&fs::exists(root/"LEC_0_components.wxf")) {
            c.flatten({{0},{1,2}});auto rows=c.to_sparse_mat(pool);
            auto components=sparse_mat_read_wxf<Q,I>(root/"LEC_0_components.wxf",field_t(FIELD_QQ));
            rows=nmhv_sew_transform_product_basis(rows,joint_eet_identity<Q,I>(9),components,field_t(FIELD_QQ),pool);
            c=sparse_tensor<Q,I,SPARSE_COO>(rows,pool);c.reshape({n,9,old});
        }
        write(out/("w"+std::to_string(w)+".wxf"),Tensor(std::move(c),pool));previous=current;
    }
    auto cyclic=sparse_mat_read_wxf<Q,I>(data/"cycrepmat.wxf",field_t(FIELD_QQ));
    auto square=sparse_mat_mul(cyclic,cyclic,field_t(FIELD_QQ));auto parity=sparse_mat_mul(square,cyclic,field_t(FIELD_QQ));
    write(out/"parity.wxf",Tensor(parity,pool));
    if(backward&&fs::exists(root/"LEC_0_components.wxf")) {
        auto blocks=nmhv_sew_read_irrep_axis(root/"LEC_0_blocks.tsv");
        auto components=sparse_mat_read_wxf<Q,I>(root/"LEC_0_components.wxf",field_t(FIELD_QQ));
        auto inverse=sparse_mat_inverse(components,field_t(FIELD_QQ),opt);
        for(const auto& name:{"Cyclic","Flip","Parity"}) {
            Mat action(components.nrow,components.nrow);size_t offset=0;
            for(size_t b=0;b<6;++b){size_t n=blocks[b].multiplicity*irr::specs()[b].irrep_dimension;
                auto block=irr::apply_action(b,name,joint_eet_identity<Q,I>(n));
                for(size_t r=0;r<n;++r)for(auto [j,x]:block[r])action[offset+r].push_back(I(offset+j),x);offset+=n;}
            action=sparse_mat_mul(sparse_mat_mul(inverse,action,field_t(FIELD_QQ),pool),components,field_t(FIELD_QQ),pool);
            write(out/(std::string(name)+"_terminal.wxf"),Tensor(action,pool));
        }
    }
    std::cout<<"REFERENCE EXPORTED "<<out<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
