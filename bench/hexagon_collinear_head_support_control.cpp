#include "irr_physics_pipeline.hpp"
int main(int argc,char**argv)try{
    require(argc==5,"head_support_control FEC LEC DATA FRESH_OUT");
    Flint::set_memory_functions();fs::path out=argv[4];require(!fs::exists(out),"Fresh output");
    for(ulong p:{1152921504606847009UL,1152921504606847067UL}){
        irr_physics::context c(argv[1],argv[2],out/std::to_string(p),argv[3],15,3,2,true,p);
        c.k=2;c.tail_words={0,1,2,3};irr_physics::axis last({2,0,0,0,0,0});
        c.exact_tail.assign(32,0);
        std::array<std::array<ulong,2>,4> coefficients{{{1,0},{0,1},{1,1},{0,0}}};
        std::array<std::array<ulong,2>,4> rhs{{{3,4},{6,7},{9,11},{7,0}}};
        irr_physics::equations rows(8);
        for(size_t beta=0;beta<4;++beta){for(size_t j=0;j<2;++j)c.exact_tail[j*16+beta]=coefficients[beta][j];
            for(size_t s=0;s<2;++s)rows[s*4+beta].rhs=rhs[beta][s];}
        std::vector<unsigned char> active{1,0};
        auto selected=irr_physics::independent_collinear_suffix(c,last,rows,2,&active);
        require(selected==std::vector<size_t>{0,1,3},"Inactive coefficient with independent RHS must remain");
        std::vector<unsigned char> none{0,0};
        selected=irr_physics::independent_collinear_suffix(c,last,rows,2,&none);
        require(selected==std::vector<size_t>{0,1},"RHS-only operator rank must remain");
        for(auto& r:rows)r.rhs=0;
        selected=irr_physics::independent_collinear_suffix(c,last,rows,2,&active);
        require(selected==std::vector<size_t>{0},"Provably dead coefficient column can be omitted");
        std::cout<<"HEAD_SUPPORT_CONTROL_PASS prime="<<p<<" cases=3"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
