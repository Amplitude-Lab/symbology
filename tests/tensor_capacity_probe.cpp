#include "bootstrap.hpp"
using Q=rat_t;using I=int32_t;
int main(){
    using COO=sparse_tensor<Q,I,SPARSE_COO>;using CSR=sparse_tensor<Q,I,SPARSE_CSR>;
    // Reproduces the four-loop crash: cancellation leaves far fewer live
    // entries than the COO allocation, then CSR is copied for later use.
    for(size_t rank:{size_t(2),size_t(8)})for(size_t live:{size_t(0),size_t(3),size_t(31)}){
        COO a(std::vector<size_t>(rank,40),20000);
        for(size_t i=0;i<live;++i){std::vector<I> word(rank,0);word[0]=i;word.back()=i;Q large(int_t(2).pow(140ul)+int_t(i));a.push_back(word,large);}
        CSR b(std::move(a));if(b.alloc()!=live)return 1;CSR copied=b,assigned;assigned=b;
        for(size_t i=0;i<live;++i)if(copied.val(i)!=b.val(i)||assigned.val(i)!=b.val(i))return 2;
        if(live){b.reserve(live+20);if(b.val(0)!=copied.val(0))return 3;COO again(std::move(b));CSR roundtrip(std::move(again));if(roundtrip.nnz()!=live)return 4;}
    }
    Q target(int_t(2).pow(180ul));int_t z=int_t(2).pow(180ul);
    for(size_t i=0;i<10000;++i){Q next(int_t(2).pow(180ul)+int_t(i));target=std::move(next);int_t n=int_t(2).pow(180ul)+int_t(i);z=std::move(n);}
    if(target.num()!=z)return 5;
    std::cout<<"PASS tensor capacity after cancellation, copy/assignment/growth, and large scalar moves"<<std::endl;
}
