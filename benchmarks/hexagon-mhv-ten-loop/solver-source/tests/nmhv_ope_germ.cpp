#include "nmhv_ope_germ.hpp"
#include <iostream>
int main(int argc,char** argv)try{
    if(argc!=2)throw std::runtime_error("test_nmhv_ope_germ DATA");
    Flint::set_memory_functions();std::filesystem::path data=argv[1];size_t checks=0;
    for(int f=2;f<=3;++f){
        mhv_ope_germ_t old(6,1000003,data/"ope"/("kernels_F"+std::to_string(f)+".tsv"));
        nmhv_ope_germ_t now(6,1000003,data/"ope_nmhv"/("kernels_F"+std::to_string(f)+".tsv"),2,18);
        // Compare the old ordered integration on its full common truncation.
        for(size_t first=0;first<9;++first)for(size_t second=0;second<9;++second){
            std::vector<ulong> a(old.stride),b(now.stride),x(old.stride),y(now.stride);
            a[old.index(0,0,0,0)]=b[now.index(0,0,0,0)]=1;
            for(size_t w=1;w<=6;++w){size_t letter=(w%2)?first:second;old.integrate(a.data(),x.data(),letter,w);now.integrate(b.data(),y.data(),letter,w);a.swap(x);b.swap(y);
                for(int t=0;t<=2;++t)for(int z=-t;z<=6-t;++z)for(size_t lt=0;lt<=w;++lt)for(size_t lz=0;lz<=w-lt;++lz){
                    if(a[old.index(t,z,lt,lz)]!=b[now.index(t,z,lt,lz)])throw std::runtime_error("Ordered germ regression");++checks;}
            }
        }
    }
    std::cout<<"ORDERED GERM PASS coefficients="<<checks<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
