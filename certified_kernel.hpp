#pragma once
// Replayable modular kernel reconstruction. Retains a finite-field equation
// matrix, not a second rational copy. A supplied EXACT residual verifier plus
// the modular rank and private free columns certify a complete Q-kernel.
#include "symrep.hpp"
#include "streaming_echelon.hpp"
namespace certified_kernel {
using namespace symrep;
using ModMat=sparse_mat<ulong,I>;
using IntMat=sparse_mat<int_t,I>;
struct BadPrime:std::runtime_error{using std::runtime_error::runtime_error;};
// A private unit column for every row is already an identity pivot block.
// Re-pivoting such a matrix can replace its unit pivots by arbitrary large
// coefficients and make rational reconstruction needlessly expensive.
// Detect the certificate directly; ordinary dependent equation matrices fall
// back to SparseRREF unchanged. No rank guess or supplied pivot metadata is used.
inline bool private_unit_pivots(const ModMat& a,std::vector<std::vector<pivot_t<I>>>& pivots){
    if(a.nrow>a.ncol)return false;
    std::vector<I> owner(a.ncol,-1),chosen(a.nrow,-1);
    std::vector<unsigned char> unit(a.ncol);
    for(size_t r=0;r<a.nrow;++r)for(auto [c,x]:a[r]){
        if(owner[c]==-1){owner[c]=index(r);unit[c]=x==1;}else{owner[c]=-2;unit[c]=0;}
    }
    for(size_t c=0;c<a.ncol;++c)if(owner[c]>=0&&unit[c]&&chosen[owner[c]]<0)chosen[owner[c]]=index(c);
    if(std::find(chosen.begin(),chosen.end(),I(-1))!=chosen.end())return false;
    pivots.assign(1,{});pivots[0].reserve(a.nrow);
    for(size_t r=0;r<a.nrow;++r)pivots[0].emplace_back(index(r),chosen[r]);
    return true;
}
template<class Generate,class Verify,class GenerateMod=std::nullptr_t>
Mat solve(size_t columns,Generate&& generate,Verify&& verify,rref_option_t opt,GenerateMod generate_mod=nullptr,size_t batch_rows=0){
    ulong prime=1ULL<<60;size_t best_rank=0;std::vector<I> pattern;IntMat residues;Mat rational_rows;std::vector<unsigned char> encoded;int_t modulus=1;
    for(size_t attempt=0;attempt<16;++attempt){
        prime=n_nextprime(prime,0);field_t fp(FIELD_Fp,prime);ModMat a(0,columns);bool bad_denominator=false;
        auto begin=std::chrono::steady_clock::now();size_t input_nnz=0;
        std::unique_ptr<streaming_echelon::Reducer> reducer;
        if(batch_rows)reducer=std::make_unique<streaming_echelon::Reducer>(columns,fp,opt,batch_rows);
        auto accept=[&](sparse_vec<ulong,I>&& row){if(!row.nnz())return;input_nnz+=row.nnz();if(reducer)reducer->consume(std::move(row));else a.rows.push_back(std::move(row));};
        if constexpr(std::is_same_v<GenerateMod,std::nullptr_t>){
            generate([&](Vec&& row){if(!row.nnz())return;sparse_vec<ulong,I> converted;converted.reserve(row.nnz());
                for(auto [j,x]:row){if(x.den()%fp.mod==0){bad_denominator=true;continue;}auto v=x%fp.mod;if(v)converted.push_back(j,v);}
                accept(std::move(converted));});
        }else{try{generate_mod(fp,accept);}catch(const BadPrime&){bad_denominator=true;}}
        a.nrow=a.rows.size();
        if(bad_denominator)continue;
        std::vector<std::vector<pivot_t<I>>> pivots;
        if(reducer){auto reduced=reducer->finish();a=std::move(reduced.first);pivots=std::move(reduced.second);reducer.reset();}
        else {std::cout<<"modular_assembly rows="<<a.nrow<<" cols="<<columns<<" nnz="<<input_nnz<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<std::endl;
            if(private_unit_pivots(a,pivots))std::cout<<"modular_private_unit_pivots="<<a.nrow<<std::endl;
            else {a.sort_rows_by_nnz();pivots=sparse_mat_rref(a,fp,opt);}}
        require(!opt->abort,"modular kernel aborted");
        std::cout<<"modular_elimination_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<std::endl;
        std::vector<I> pivot_row(columns,-1),free_index(columns,-1),free_columns;size_t rank=0;
        for(const auto& batch:pivots)for(auto [r,c]:batch){pivot_row[c]=r;++rank;}
        for(size_t c=0;c<columns;++c)if(pivot_row[c]<0){free_index[c]=index(free_columns.size());free_columns.push_back(index(c));}
        std::cout<<"modular_attempt="<<attempt+1<<" rows="<<a.nrow<<" cols="<<columns<<" input_nnz="<<input_nnz<<" rank="<<rank<<std::endl;
        if(rank<best_rank)continue;
        ModMat k(free_columns.size(),columns);
        std::vector<size_t> counts(k.nrow,1);
        for(size_t c=0;c<columns;++c)if(pivot_row[c]>=0)for(auto [j,x]:a[pivot_row[c]])if(size_t(j)!=c&&x){require(free_index[j]>=0,"modular kernel requires back substitution");++counts[free_index[j]];}
        for(size_t r=0;r<k.nrow;++r)k[r].reserve(counts[r]);
        for(size_t c=0;c<columns;++c){if(pivot_row[c]<0)k[free_index[c]].push_back(index(c),prime-1);
            else for(auto [j,x]:a[pivot_row[c]])if(size_t(j)!=c&&x){require(free_index[j]>=0,"modular kernel requires back substitution");k[free_index[j]].push_back(index(c),x);}}
        a.clear();
        if(modulus==1||rank>best_rank||free_columns!=pattern){modulus=1;residues=IntMat(k.nrow,columns);rational_rows=Mat(k.nrow,columns);encoded.assign(k.nrow,0);pattern=free_columns;best_rank=rank;}
        // All coefficients share M and p. Precompute M^{-1} mod p once and
        // update equal supports in place; a second complete CRT matrix is
        // otherwise especially costly when residues have become big integers.
        const bool initial=modulus==1;const ulong inverse=initial?1:n_invmod(modulus%fp.mod,prime);const int_t new_modulus=modulus*prime;
        auto combine=[&](int_t& x,ulong y){if(initial){x=y;return;}
            ulong correction=nmod_mul(nmod_sub(y,x%fp.mod,fp.mod),inverse,fp.mod);
            fmpz_addmul_ui(x.data(),modulus.data(),correction);
        };
        std::atomic<bool> reconstructed=true;
        opt->pool.detach_loop(size_t(0),k.nrow,[&](size_t r){auto& old=residues[r];const auto& now=k[r];
            // Successful rational rows are a compact encoding of their CRT
            // residues. Decode only this row if another prime is necessary.
            if(encoded[r]){old.reserve(rational_rows[r].nnz());
                for(auto [c,x]:rational_rows[r]){int_t denominator_inverse,value;
                    if(x.is_integer())fmpz_mod(value.data(),x.num_data(),modulus.data());
                    else{require(fmpz_invmod(denominator_inverse.data(),x.den_data(),modulus.data()),"reconstructed denominator is not invertible modulo CRT modulus");
                        fmpz_mul(value.data(),x.num_data(),denominator_inverse.data());fmpz_mod(value.data(),value.data(),modulus.data());}
                    if(value!=0)old.push_back(c,std::move(value));}
                rational_rows[r].clear();encoded[r]=0;
            }
            if(old.nnz()==now.nnz()&&std::equal(old.index_span().begin(),old.index_span().end(),now.index_span().begin())){
                for(size_t j=0;j<old.nnz();++j)combine(old[j],now[j]);
            }else{
                size_t i=0,j=0;sparse_vec<int_t,I> joined;joined.reserve(old.nnz()+now.nnz());
                while(i<old.nnz()||j<now.nnz()){
                    I c=std::min(i<old.nnz()?old(i):INT32_MAX,j<now.nnz()?now(j):INT32_MAX);int_t x=0;ulong y=0;
                    if(i<old.nnz()&&old(i)==c)x=std::move(old[i++]);if(j<now.nnz()&&now(j)==c)y=now[j++];
                    combine(x,y);if(x!=0)joined.push_back(c,std::move(x));}
                old=std::move(joined);
            }
            k[r].clear();
            // One-prime residues still fit one machine word. Once any row
            // needs another prime, avoid futile reconstruction of later rows;
            // retaining these integer residues is smaller than retaining Q.
            if(initial&&!reconstructed.load(std::memory_order_relaxed))return;
            Vec candidate;candidate.reserve(old.nnz());bool success=true;
            for(auto [j,x]:old){Q value;if(!rational_reconstruct(value,x,new_modulus)){success=false;break;}if(value!=0)candidate.push_back(j,std::move(value));}
            if(success){rational_rows[r]=std::move(candidate);encoded[r]=1;old.clear();}else reconstructed=false;
        });opt->pool.wait();
        modulus=new_modulus;k.clear();
        std::cout<<"rational_reconstruction modulus_bits="<<modulus.bits()<<" reconstructed_rows="<<std::count(encoded.begin(),encoded.end(),1)
                 <<" total_rows="<<rational_rows.nrow<<" rational_nnz="<<rational_rows.nnz()<<" residue_nnz="<<residues.nnz()<<std::endl;
        if(!reconstructed)continue;
        auto& result=rational_rows;
        // The CRT/reconstruction retains an identity on these private columns.
        // Check it explicitly before accepting the residual certificate.
        for(size_t i=0;i<result.nrow;++i)for(auto [c,x]:result[i])if(free_index[c]>=0)
            require(free_index[c]==index(i)&&x==-1,"lost private free coordinate during reconstruction");
        for(size_t i=0;i<result.nrow;++i)require(result[i].find(pattern[i])&&*result[i].find(pattern[i])==-1,"missing private free coordinate");
        // Reconstruction is congruent to the computed modular kernel. A
        // factored source can use this prime with a rigorous residual-height
        // bound, avoiding a second large matrix/kernel multiplication.
        bool valid;
        if constexpr(std::is_invocable_r_v<bool,Verify,const Mat&,ulong,const int_t&>)valid=verify(result,prime,modulus);
        else if constexpr(std::is_invocable_r_v<bool,Verify,const Mat&,ulong>)valid=verify(result,prime);
        else valid=verify(result);
        if(valid){std::cout<<"reconstruction_modulus_bits="<<modulus.bits()<<std::endl;return std::move(result);}
    }
    throw std::runtime_error("Exact kernel reconstruction did not certify within 16 good-prime attempts");
}
}
