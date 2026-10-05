#pragma once
#include "symrep_io.hpp"
#include <unordered_map>

// Exact word shuffles with one integer key per word. Mixed-radix packing is
// injective only while alphabet^weight fits in uint64_t; overflow is rejected.
// Several weighted products can share the accumulator, avoiding separate
// materialized terms and repeated tensor additions.
namespace packed_shuffle {
using namespace symrep;
using COO=sparse_tensor<Q,I,SPARSE_COO>;
class Accumulator {
    size_t weight,alphabet;
    std::vector<uint64_t> powers;
    std::unordered_map<uint64_t,Q> sums;
public:
    Accumulator(size_t w,size_t a):weight(w),alphabet(a),powers(w){
        require(w>0&&a>0,"packed shuffle needs positive word dimensions");
        uint64_t p=1;for(size_t i=w;i-->0;){powers[i]=p;require(p<=UINT64_MAX/a,"packed shuffle word key overflow");p*=a;}
    }
    size_t stored_words()const{return sums.size();}
    void add(const COO& first,const COO& second,const Q& scale=Q(1)){
        require(first.rank()+second.rank()==weight,"packed shuffle weight mismatch");
        for(auto d:first.dims())require(d<=alphabet,"packed shuffle alphabet mismatch");
        for(auto d:second.dims())require(d<=alphabet,"packed shuffle alphabet mismatch");
        if(scale==0||!first.nnz()||!second.nnz())return;
        const COO& a=first.nnz()<=second.nnz()?first:second;
        const COO& b=first.nnz()<=second.nnz()?second:first;
        std::vector<std::vector<size_t>> ap,bp;std::vector<unsigned char> mask(weight);
        std::fill(mask.begin(),mask.begin()+a.rank(),1);
        do{std::vector<size_t> x,y;for(size_t i=0;i<weight;++i)(mask[i]?x:y).push_back(i);ap.push_back(std::move(x));bp.push_back(std::move(y));}while(std::prev_permutation(mask.begin(),mask.end()));
        const size_t patterns=ap.size();
        // Cache only a bounded block of one factor, regardless of input size.
        const size_t block=std::max(size_t(1),size_t(4*1024*1024)/patterns);
        std::vector<uint64_t> right(patterns);
        for(size_t begin=0;begin<a.nnz();begin+=block){size_t end=std::min(a.nnz(),begin+block);
            std::vector<uint64_t> left((end-begin)*patterns);std::vector<Q> values(end-begin);
            for(size_t i=begin;i<end;++i){auto word=a.index_vector(i);values[i-begin]=a.val(i)*scale;
                for(size_t s=0;s<patterns;++s)for(size_t j=0;j<word.size();++j)left[(i-begin)*patterns+s]+=uint64_t(word[j])*powers[ap[s][j]];}
            for(size_t j=0;j<b.nnz();++j){auto word=b.index_vector(j);std::fill(right.begin(),right.end(),0);
                for(size_t s=0;s<patterns;++s)for(size_t p=0;p<word.size();++p)right[s]+=uint64_t(word[p])*powers[bp[s][p]];
                for(size_t i=begin;i<end;++i){Q value=values[i-begin]*b.val(j);if(value==0)continue;
                    for(size_t s=0;s<patterns;++s)sums[left[(i-begin)*patterns+s]+right[s]]+=value;}
            }
        }
    }
    COO finish(){
        std::erase_if(sums,[](const auto& item){return item.second==0;});
        COO result(std::vector<size_t>(weight,alphabet));result.reserve(sums.size());std::vector<I> word(weight);
        for(auto it=sums.begin();it!=sums.end();){auto key=it->first;for(size_t i=weight;i-->0;){word[i]=index(key%alphabet);key/=alphabet;}
            result.push_back(word,it->second);it=sums.erase(it);}
        sums.rehash(0);result.sort_indices();return result;
    }
};
}
