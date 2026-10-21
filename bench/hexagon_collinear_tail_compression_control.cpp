#include "irr_physics_pipeline.hpp"
int main(int argc,char**argv)try{
    require(argc==5,"suffix_compression_control FEC LEC DATA FRESH_OUT");
    Flint::set_memory_functions();fs::path out=argv[4];require(!fs::exists(out),"Fresh output");
    for(ulong p:{1152921504606847009UL,1152921504606847067UL}){
        irr_physics::context c(argv[1],argv[2],out/std::to_string(p),argv[3],15,3,2,true,p);
        c.k=2;c.tail_words={0,1,2,3,4,5};irr_physics::axis last({3,0,0,0,0,0});
        c.exact_tail.assign(last.dim*16,0);
        // Two coefficient rows, two dependencies, an RHS-only constraint,
        // and an entirely zero row. Rank in the affine operator is three.
        std::array<std::array<ulong,3>,6> coefficients{{{1,0,0},{0,1,0},{1,1,0},{2,0,0},{0,0,0},{0,0,0}}};
        std::array<std::array<ulong,2>,6> rhs{{{3,4},{5,6},{8,10},{6,8},{7,0},{0,0}}};
        irr_physics::equations rows(12);
        for(size_t beta=0;beta<6;++beta){for(size_t j=0;j<3;++j)c.exact_tail[j*16+beta]=coefficients[beta][j];
            for(size_t s=0;s<2;++s)rows[s*6+beta].rhs=rhs[beta][s];}
        auto selected=irr_physics::independent_collinear_suffix(c,last,rows,2);
        require(selected==std::vector<size_t>{0,1,4},"RHS-only constraint must remain");
        rows[1*6+2].rhs=11;
        selected=irr_physics::independent_collinear_suffix(c,last,rows,2);
        require(selected==std::vector<size_t>{0,1,2,4},"Changed affine RHS rank must remain");
        for(auto& r:rows)r.rhs=0;std::fill(c.exact_tail.begin(),c.exact_tail.end(),0);
        selected=irr_physics::independent_collinear_suffix(c,last,rows,2);
        require(selected==std::vector<size_t>{0},"Zero-operator adaptive stop accounting");
        std::cout<<"SUFFIX_COMPRESSION_CONTROL_PASS prime="<<p<<" cases=3"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
