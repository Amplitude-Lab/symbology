#include "packed_shuffle.hpp"
#include "tensor_shuffle.h"
#include <random>
using namespace symrep;
using COO=packed_shuffle::COO;
int main(){try{
    std::mt19937 random(51003);
    for(size_t trial=0;trial<64;++trial){size_t letters=2+trial%3,wa=1+trial%3,wb=1+(trial/3)%3;
        COO a(std::vector<size_t>(wa,letters)),b(std::vector<size_t>(wb,letters));
        auto fill=[&](COO& t){for(size_t i=0;i<8;++i){std::vector<I> word(t.rank());for(auto& c:word)c=index(random()%letters);t.push_back(word,Q(long(random()%11)-5,long(1+random()%4)));}t.canonicalize();t.sort_indices();};
        fill(a);fill(b);if(!a.nnz()||!b.nnz())continue;
        auto expected=tensor_shuffle_product_parallel(a,b,QQ,nullptr);Q scale(long(trial%7)-3,5);
        for(size_t p=0;p<expected.nnz();++p)expected.val(p)*=scale;expected.canonicalize();
        packed_shuffle::Accumulator sum(wa+wb,letters);sum.add(a,b,scale+Q(2));sum.add(b,a,Q(-2));auto actual=sum.finish();
        require(equal(flatten(Tensor(actual)),flatten(Tensor(expected))),"packed weighted shuffle differs from independent word-vector implementation");
    }
    COO empty({2}),one({2});one.push_back(std::vector<I>{0},Q(1));packed_shuffle::Accumulator zero(2,2);zero.add(empty,one);require(!zero.finish().nnz(),"empty shuffle was not zero");
    bool overflow=false;try{packed_shuffle::Accumulator too_big(65,2);}catch(const std::exception&){overflow=true;}require(overflow,"packed word overflow was not rejected");
    std::cout<<"PASS 64 independent rational weighted shuffle cases, cancellation, empty factors and key overflow"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
