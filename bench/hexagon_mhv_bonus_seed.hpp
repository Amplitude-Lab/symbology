#pragma once
#include "hexagon_mhv_double_bonus.hpp"
#include "shared/streaming_echelon.hpp"

// Optional empirical ansatz: impose the published double-final bonus
// relations as short homogeneous rows in the ORIGINAL invariant coordinates.
// All six irrep sectors retain the same component convention and global
// affine normalization. The ordinary LEC inputs and conversion chain remain
// unchanged. This restriction must be named in the candidate descriptor.
namespace mhv_double_bonus {
inline void seed_original_coordinates(irr_physics::context& c,
    streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate){
    irr::require(c.mhv&&c.k>=2&&c.k<=4,"MHV bonus seed split");
    auto start=std::chrono::steady_clock::now();
    auto words=expand_lec_words(c.lec,c.k,c.cg,c.opt);
    auto B=alphabet_matrix();calibrate_single_final_entry(words[1],B);
    auto ordinary=coefficient_dual_transform(published_qbar_calibration(),B);
    auto calibration=constraint_images(words[2],irr::chain_copies(c.lec,2,false),2,ordinary);
    for(const auto& block:calibration)require(!block.nnz(),"Bonus seed alphabet calibration");
    auto constraints=coefficient_dual_transform(full_d3_orbit(published_bonus_relations()),B);
    auto restriction=restrict_tail(words[c.k],c.sew.last,c.k,constraints,c.opt);
    if(c.k==3)require(dimension(restriction.reduced)==65,"Bonus seed expected 65 triple tails");
    size_t offset=0,submitted=0,nnz=0,rank_before=reducer.rank(),conditional_columns=0;
    for(size_t b=0;b<6;++b){
        auto equations=nmhv_sew_modular_reduce(restriction.images[b].transpose(),c.F,&c.opt->pool);
        equations.clear_zero_row();equations.sort_rows_by_nnz();
        auto pivots=sparse_mat_rref_forward(equations,c.F,c.opt);
        std::vector<size_t> independent;
        for(const auto& batch:pivots)for(auto [r,col]:batch)independent.push_back(r);
        require(independent.size()==restriction.original[b]-restriction.reduced[b],"Bonus copy rank differs at this prime");
        for(size_t f=0;f<c.sew.first[b];++f)for(auto r:independent){
            streaming_echelon::ModVec row;
            for(auto [l,value]:equations[r])row.push_back(coordinate.at(offset+f*c.sew.last[b]+l),value);
            row.sort_indices();nnz+=row.nnz();++submitted;reducer.consume(std::move(row));
        }
        offset+=c.sew.first[b]*c.sew.last[b];
        conditional_columns+=c.sew.first[b]*restriction.reduced[b];
    }
    require(offset==coordinate.size(),"Bonus seed original-coordinate layout");
    reducer.flush_pending();reducer.prune_zeros();
    require(!reducer.is_zero(offset),"Homogeneous bonus rows changed affine normalization");
    std::cout<<"CONDITIONAL_MHV_DOUBLE_BONUS equations="<<submitted<<" input_nnz="<<nnz
        <<" rank_gain="<<reducer.rank()-rank_before<<" rank="<<reducer.rank()
        <<" original_columns="<<offset<<" conditional_tail_columns="<<conditional_columns
        <<" full_D3_orbit=1 shared_E_copy_map=1 empirical_ten_loop_assumption=1 seconds="
        <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
}
} // namespace mhv_double_bonus
