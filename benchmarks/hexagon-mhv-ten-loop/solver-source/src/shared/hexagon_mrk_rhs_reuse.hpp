#pragma once
#include "joint_eet_mrk.hpp"
#include <limits>

// Evaluate the same MRK shuffle expression with common subexpressions shared.
// No physical coefficient, sample, degree or normalization is changed.
namespace hexagon_mrk_reuse {
using polynomial = std::vector<ulong>;
struct term { size_t word; ulong coefficient; };
struct symbol { size_t weight; std::vector<term> terms; };

class evaluator {
    const std::vector<std::array<ulong,3>>& slots;
    const field_t& field;
    size_t loop,degree,n,states;
    std::vector<std::vector<size_t>> offsets;
    std::vector<ulong> arena;
    std::vector<std::array<ulong,2>> tau;
    std::vector<std::array<ulong,3>> pairs;
    std::map<std::pair<size_t,size_t>,symbol> symbols;
    ulong add(ulong a,ulong b)const{return nmod_add(a,b,field.mod);}
    ulong mul(ulong a,ulong b)const{return nmod_mul(a,b,field.mod);}

    // H(p,M) = (Q shuffle ... shuffle Q)/p! shuffle T^shuffle d
    // on the specified slots, with d=|M|-2p. Choose the least slot:
    // d choices make it a T; otherwise its Q partner determines one pair.
    size_t common(size_t p,size_t mask){
        auto& cached=offsets[p][mask];
        if(cached!=SIZE_MAX)return cached;
        size_t r=std::popcount(mask);
        if(r<2*p)throw std::runtime_error("Invalid MRK common factor weight");
        const size_t logs=r-2*p;
        polynomial value(std::min(degree,r)+1);
        if(!mask)value[0]=1;
        else {
            const size_t i=std::countr_zero(mask),rest=mask^(size_t(1)<<i);
            if(logs){
                size_t child=common(p,rest),last=std::min(degree,r-1);
                for(size_t q=0;q<=last;++q){
                    auto x=mul(logs,arena[child+q]);
                    value[q]=add(value[q],mul(tau[i][0],x));
                    if(q<degree)value[q+1]=add(value[q+1],mul(tau[i][1],x));
                }
            }
            if(p)for(size_t bits=rest;bits;bits&=bits-1){
                size_t j=std::countr_zero(bits),child=common(p-1,rest^(size_t(1)<<j));
                const auto& pair=pairs[i*n+j];size_t last=std::min(degree,r-2);
                for(size_t q=0;q<=last;++q){auto x=arena[child+q];if(!x)continue;
                    for(size_t k=0;k<3&&q+k<=degree;++k)
                        value[q+k]=add(value[q+k],mul(pair[k],x));
                }
            }
        }
        cached=arena.size();arena.insert(arena.end(),value.begin(),value.end());
        return cached;
    }
    ulong evaluate_symbol(const symbol& sym,size_t mask){
        std::vector<size_t> indices;indices.reserve(sym.weight);
        for(size_t bits=mask;bits;bits&=bits-1)indices.push_back(std::countr_zero(bits));
        if(indices.size()!=sym.weight)throw std::runtime_error("MRK symbol mask width");
        std::vector<ulong> products;
        if(sym.terms.size()>=(size_t(1)<<sym.weight)/4){
            products.resize(size_t(1)<<sym.weight);products[0]=1;size_t count=1;
            for(auto i:indices){for(size_t j=count;j-->0;){auto x=products[j];
                products[2*j]=mul(x,slots[i][1]);products[2*j+1]=mul(x,slots[i][2]);}
                count*=2;
            }
        }
        ulong result=0;
        for(const auto& t:sym.terms){ulong x=t.coefficient;
            if(!products.empty())x=mul(x,products[t.word]);
            else for(size_t j=0;j<sym.weight;++j)x=mul(x,slots[indices[j]][1+((t.word>>(sym.weight-1-j))&1)]);
            result=add(result,x);
        }
        return result;
    }
    template<class Function> void subsets(size_t available,size_t count,size_t chosen,Function& function){
        if(!count){function(chosen);return;}
        size_t left=std::popcount(available);
        if(left<count)return;
        if(left==count){function(chosen|available);return;}
        while(left>=count){size_t bit=available&(~available+1);available^=bit;--left;
            subsets(available,count-1,chosen|bit,function);
        }
    }
public:
    evaluator(size_t L,const std::vector<std::array<ulong,3>>& input,
              const nmhv_physical_symbol_table_t& data,size_t D,const field_t& F)
        :slots(input),field(F),loop(L),degree(D),n(input.size()),states(size_t(1)<<n),
         offsets(L,std::vector<size_t>(states,SIZE_MAX)),tau(n),pairs(n*n){
        if(!L||n!=2*L-1||n>=8*sizeof(size_t))throw std::runtime_error("Invalid MRK RHS input");
        const auto inv2=n_invmod(2,F.mod.n),minus2=nmod_neg(2,F.mod);
        for(size_t i=0;i<n;++i)tau[i]={nmod_sub(mul(inv2,slots[i][1]),slots[i][2],F.mod),slots[i][0]};
        for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j){auto& q=pairs[i*n+j];
            q[0]=nmod_sub(mul(minus2,mul(tau[i][0],tau[j][0])),mul(inv2,mul(slots[i][1],slots[j][1])),F.mod);
            q[1]=mul(minus2,add(mul(tau[i][0],tau[j][1]),mul(tau[i][1],tau[j][0])));
            q[2]=mul(minus2,mul(tau[i][1],tau[j][1]));
        }
        for(size_t l=1;l<=L;++l)for(size_t d=0;d<l;++d){
            auto found=data.find({l,d});if(found==data.end())throw std::runtime_error("Missing MRK symbol");
            auto& s=symbols[{l,d}];s.weight=2*l-1-d;
            for(const auto& t:found->second){
                if(t.word.size()!=s.weight)throw std::runtime_error("Wrong MRK word length");
                size_t word=0;for(auto letter:t.word){if(letter>1)throw std::runtime_error("Nonbinary MRK word");word=2*word+letter;}
                auto x=t.coefficient%F.mod;if(x)s.terms.push_back({word,x});
            }
        }
    }
    polynomial evaluate(){
        polynomial result(degree+1);size_t full=states-1;
        const auto quarter=n_invmod(4,field.mod.n);
        for(size_t l=1;l<=loop;++l){size_t p=loop-l;
            auto scale=nmod_neg(mul(quarter,nmhv_physical_detail::power(2,l,field)),field.mod);
            for(size_t d=0;d<l;++d){const auto& sym=symbols.at({l,d});if(sym.terms.empty())continue;
                auto accumulate=[&](size_t chosen){
                    auto x=mul(scale,evaluate_symbol(sym,chosen));if(!x)return;
                    auto remaining=full^chosen,offset=common(p,remaining);
                    size_t last=std::min(degree,size_t(std::popcount(remaining)));
                    for(size_t q=0;q<=last;++q)result[q]=add(result[q],mul(x,arena[offset+q]));
                };
                subsets(full,sym.weight,0,accumulate);
            }
        }
        // Derivative + phase/2 equals the same tau polynomial for the last slot.
        // Q^(L-1)/(L-1)! is already normalized in common().
        for(size_t i=0;i<n;++i){size_t remaining=full^(size_t(1)<<i),offset=common(loop-1,remaining);
            for(size_t q=0;q<=std::min(degree,n-1);++q){auto x=arena[offset+q];
                result[q]=add(result[q],mul(tau[i][0],x));
                if(q<degree)result[q+1]=add(result[q+1],mul(tau[i][1],x));
            }
        }
        for(auto& x:result)x=mul(2,x);
        return result;
    }
    size_t memo_bytes()const{return loop*states*sizeof(size_t)+arena.capacity()*sizeof(ulong);}
};
inline std::vector<ulong> rhs(size_t loop,const std::vector<std::array<ulong,3>>& slots,
    const nmhv_physical_symbol_table_t& symbols,size_t degree,const field_t& field){
    return evaluator(loop,slots,symbols,degree,field).evaluate();
}
} // namespace hexagon_mrk_reuse
