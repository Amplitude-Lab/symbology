#pragma once
#include "tensor_kernel.hpp"
#include "product_invariance.hpp"

namespace early_sew {
using namespace symrep;
using ProductInvarianceRows = staged_kernel::ProductInvarianceRows;
struct Rows {
    tensor_kernel::ConstraintRows& equations;
    ProductInvarianceRows& invariants;
    size_t max_terms=SIZE_MAX;
    size_t min_terms=0;
    const streaming_echelon::Substitution* substitution=nullptr;
    std::vector<unsigned char>* equation_seen=nullptr;
    std::vector<unsigned char>* invariant_seen=nullptr;
    unsigned sample_shift=0;
    size_t ncol()const{return equations.ncol();}
    template<class Emit> void operator()(structured_kernel::Relations& relations,size_t limit,Emit&& emit)const{
        equations(relations,limit,emit);invariants(relations,limit,emit);
    }
    template<class Emit> void modular(structured_kernel::Relations& relations,const std::vector<I>& columns,const field_t& fp,Emit&& emit,thread_pool* pool=nullptr)const{
        invariants.modular(relations,columns,fp,emit,pool,max_terms,min_terms,substitution,invariant_seen,sample_shift);equations.modular(relations,columns,fp,emit,pool,max_terms,min_terms,substitution,equation_seen,sample_shift);
    }
};

// Compatibility adapter for the experimental driver's discovery modes.
// All elimination, intersection and certification now live in staged_kernel.
struct RefinedRows {
    Rows source;
    rref_option* opt;
    size_t seed_terms=64,target_nullity=128;
    bool remember_equations=false,sample_equations=false;
    size_t sample_after_terms=0;
    template<class Emit> void operator()(structured_kernel::Relations& r,size_t limit,Emit&& emit)const{source(r,limit,emit);}
    template<class Emit> void modular(structured_kernel::Relations& r,const std::vector<I>& columns,const field_t& fp,Emit&& emit,thread_pool* pool)const{
        auto operators=staged_kernel::stack(source.invariants,source.equations);
        staged_kernel::Options options;
        options.seed_terms=seed_terms;options.target_nullity=target_nullity;
        options.remember_equations=remember_equations;options.sample_equations=sample_equations;
        options.complete_through_terms=sample_after_terms;
        staged_kernel::RefinedRows refined{operators,opt,options};
        refined.modular(r,columns,fp,std::forward<Emit>(emit),pool);
    }
};
inline bool verify_with_extra_primes(const Mat& candidate,const Rows& source,rref_option_t opt,ulong prime,const int_t& modulus){
    return staged_kernel::verify(candidate,staged_kernel::stack(source.invariants,source.equations),opt,prime,modulus);
}
}
