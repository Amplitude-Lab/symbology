#include "shared/direct_ope_rows.hpp"
#define nmhv_ope_germ_t nmhv_square_germ_t
#define nmhv_ope_row nmhv_square_ope_row
#include "nmhv_ope_germ_square_reference.hpp"
#undef nmhv_ope_row
#undef nmhv_ope_germ_t
using irr::require;
static ulong random_value(uint64_t& x,ulong p){x^=x<<13;x^=x>>7;x^=x<<17;return x%p;}
int main(int argc,char** argv)try{
    require(argc==5,"ope_weight_control OPE_DATA FEC LEC OUTPUT");Flint::set_memory_functions();
    uint64_t rng=0x857aba90231;size_t cases=0,guards=0;
    for(ulong prime:{ulong(101),ulong(1152921504606847009)})for(int grade:{4,12})for(int f=2;f<=4;++f){
        auto file=std::filesystem::path(argv[1])/("kernels_F"+std::to_string(f)+".tsv");
        nmhv_ope_germ_t full(20,prime,file,2,grade);nmhv_square_germ_t square(20,prime,file,2,grade);
        for(size_t w=1;w<=20;++w){auto prev=full.with_log_weight(w-1),cur=full.with_log_weight(w);
            for(size_t letter=0;letter<9;++letter){
                std::vector<ulong> in(prev.stride),out(cur.stride),refin(square.stride),refout(square.stride);
                for(auto& v:in)v=random_value(rng,prime);
                for(int t=0;t<=2;++t)for(int z=-t;z<=grade-t;++z)
                    for(size_t a=0;a<w;++a)for(size_t b=0;b<w-a;++b)
                        refin[square.index(t,z,a,b)]=in[prev.index(t,z,a,b)];
                square.integrate(refin.data(),refout.data(),letter,w);
                cur.integrate_from(prev,in.data(),out.data(),letter,w);
                for(int t=0;t<=2;++t)for(int z=-t;z<=grade-t;++z)
                    for(size_t a=0;a<=20;++a)for(size_t b=0;b<=20;++b)
                        require((a+b<=w?out[cur.index(t,z,a,b)]:0)==refout[square.index(t,z,a,b)],"Weight-dependent integration differs from square reference");
                ++cases;
            }
        }
        try{full.with_log_weight(21);throw std::logic_error("Missing layout guard");}
        catch(const std::runtime_error&){++guards;}
    }
    // Complete recurrence, CG gathering, and in-place terminal restoration.
    for(ulong prime:{ulong(1152921504606847009),ulong(1152921504606847067)}){
        auto root=std::filesystem::path(argv[4])/std::to_string(prime);
        irr_physics::context c(argv[2],argv[3],root,std::filesystem::path(argv[1]).parent_path(),8,3,1,true,prime);
        nmhv_ope_germ_t full(20,prime,std::filesystem::path(argv[1])/"kernels_F2.tsv",2,4);
        auto packed=irr_physics::compact_ope_head(c,full);
        std::vector<ulong> expected(full.stride);expected[full.index(0,0,0,0)]=1;
        for(size_t w=1;w<=c.h;++w){auto K=c.modular_recurrence(w,true);
            irr_physics::axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
            expected=irr_physics::step(c,irr::matrices{},old,next,expected,full.stride,
                [&](const ulong* a,ulong* b,size_t l){full.integrate(a,b,l,w);},"control_"+std::to_string(w),K.get());
        }
        require(packed==expected,"Complete compact OPE head differs from full layout");
        std::cout<<"OPE_WEIGHT_HEAD_CONTROL prime="<<prime<<" coefficients="<<packed.size()<<" equal=1"<<std::endl;
    }
    std::cout<<"OPE_WEIGHT_CONTROL PASS cases="<<cases<<" guards="<<guards<<" complete_heads=2 square_reference=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"OPE_WEIGHT_ERROR "<<e.what()<<std::endl;return 1;}
