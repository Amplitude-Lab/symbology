#include "shared/direct_ope_rows.hpp"
#define nmhv_ope_germ_t nmhv_square_germ_t
#define nmhv_ope_row nmhv_square_ope_row
#include "nmhv_ope_germ_square_reference.hpp"
#undef nmhv_ope_row
#undef nmhv_ope_germ_t
using irr::require;
static ulong next_value(uint64_t& x,ulong prime){x^=x<<13;x^=x>>7;x^=x<<17;return x%prime;}
static ulong dot(const std::vector<ulong>& a,const std::vector<ulong>& b,const field_t& F){
    ulong out=0;require(a.size()==b.size(),"Dot width");
    for(size_t i=0;i<a.size();++i)out=nmod_add(out,nmod_mul(a[i],b[i],F.mod),F.mod);
    return out;
}
int main(int argc,char** argv)try{
    require(argc==2,"ope_triangular_probe OPE_DATA_ROOT");Flint::set_memory_functions();
    size_t cases=0;uint64_t random=0x92874bed1563;
    for(ulong prime:{ulong(101),ulong(1152921504606847009)})for(int grade:{4,12})for(int f=2;f<=4;++f){
        auto path=std::filesystem::path(argv[1])/("kernels_F"+std::to_string(f)+".tsv");
        nmhv_ope_germ_t packed(20,prime,path,2,grade);
        nmhv_square_germ_t square(20,prime,path,2,grade);
        require(packed.stride==size_t(3*(grade+1)*231),"Packed logarithmic triangle size");
        std::vector<ulong> input(packed.stride),output(packed.stride),dual(packed.stride),adjoint(packed.stride);
        std::vector<ulong> reference_in(square.stride),reference_out(square.stride);
        for(size_t w=1;w<=20;++w)for(size_t letter=0;letter<9;++letter){
            for(auto& v:reference_in)v=next_value(random,prime);
            for(int t=0;t<=2;++t)for(int z=-t;z<=grade-t;++z)
                for(size_t a=0;a<=20;++a)for(size_t b=0;b<=20-a;++b)
                    input[packed.index(t,z,a,b)]=reference_in[square.index(t,z,a,b)];
            square.integrate(reference_in.data(),reference_out.data(),letter,w);
            packed.integrate(input.data(),output.data(),letter,w);
            for(int t=0;t<=2;++t)for(int z=-t;z<=grade-t;++z)
                for(size_t a=0;a<=20;++a)for(size_t b=0;b<=20;++b){
                    if(a+b<=20)require(output[packed.index(t,z,a,b)]==reference_out[square.index(t,z,a,b)],"Triangular integration differs from square reference");
                    else require(!reference_out[square.index(t,z,a,b)],"Discarded square slot is nonzero");
                }
            for(auto& v:dual)v=next_value(random,prime);
            irr_physics::integrate_ope_adjoint(packed,dual.data(),adjoint.data(),letter,w);
            require(dot(output,dual,packed.F)==dot(input,adjoint,packed.F),"Packed OPE adjoint is not the transpose");
            ++cases;
        }
    }
    std::cout<<"OPE_TRIANGULAR_PROBE PASS cases="<<cases<<" square_reference=1 every_letter=1 weights=1..20 primes=2 grades=4,12 discarded_slots_zero=1 adjoint_duality=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"OPE_TRIANGULAR_ERROR "<<e.what()<<std::endl;return 1;}
