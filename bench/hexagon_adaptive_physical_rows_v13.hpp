#pragma once
#include "shared/direct_physical_rows.hpp"
#include "shared/direct_ope_rows.hpp"
#include "shared/staged_kernel.hpp"

// A bounded direct-physics schedule for the augmented sewing reducer.
// Mathematical equations and final verification are unchanged; this helper
// changes only when more original-coordinate rows are generated and admitted.
namespace hexagon_adaptive_physics {
struct Policy {
    size_t light_rounds=4;
    size_t activate_nullity=50000;
    size_t target_augmented_nullity=9001;
    size_t samples_per_batch=4;
    size_t mrk_samples_per_batch=16;
    size_t ope_target_cap=456;
    size_t collinear_batches_per_cycle=8;
    size_t max_collinear_samples=4096;
    size_t max_mrk_samples=512;
    size_t stagnant_batches=2;
    size_t cycles_before_ope=2;
    double callback_seconds=120;
    size_t maximum_cycles_per_callback=1;
};
class DirectSchedule {
    irr_physics::context& c;
    uint64_t seed;
    Policy policy;
    size_t columns,collinear_offset=0,mrk_offset=0;
    size_t collinear_stagnant=0,mrk_stagnant=0;
    size_t cycles=0;
    bool collinear_done=false,mrk_done=false,activated=false,ope_done=false;
    size_t nullity(const streaming_echelon::Reducer& reducer)const{
        return columns+1-reducer.rank();
    }
    size_t admit(streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate,
                 irr_physics::equations equations,const char* condition,size_t round,
                 size_t sample_offset,size_t samples){
        size_t before=reducer.rank(),nnz=0,nonzero_rhs=0,smallest=SIZE_MAX,largest=0;
        std::vector<staged_kernel::MV> rows;rows.reserve(equations.size());
        for(auto& equation:equations){
            irr::require(equation.coefficients.size()==columns,"Adaptive direct physical width");
            staged_kernel::MV row;
            for(size_t j=0;j<columns;++j)if(equation.coefficients[j])
                row.push_back(coordinate[j],equation.coefficients[j]);
            if(equation.rhs){row.push_back(int32_t(columns),nmod_neg(equation.rhs,c.F.mod));++nonzero_rhs;}
            staged_kernel::normalize_mod(row,c.F);
            nnz+=row.nnz();smallest=std::min(smallest,row.nnz());largest=std::max(largest,row.nnz());
            if(row.nnz())rows.push_back(std::move(row));
        }
        std::sort(rows.begin(),rows.end(),[](const auto& a,const auto& b){return a.nnz()<b.nnz();});
        for(auto& row:rows)reducer.consume(std::move(row));reducer.flush_pending();
        size_t gain=reducer.rank()-before;
        std::cout<<"ADAPTIVE_PHYSICS round="<<round<<" operator="<<condition
                 <<" samples="<<sample_offset<<".."<<sample_offset+samples
                 <<" equations="<<equations.size()<<" nnz="<<nnz
                 <<" min_row_nnz="<<(smallest==SIZE_MAX?0:smallest)<<" max_row_nnz="<<largest
                 <<" nonzero_rhs="<<nonzero_rhs<<" rank_before="<<before
                 <<" rank_after="<<reducer.rank()<<" gain="<<gain
                 <<" augmented_nullity="<<nullity(reducer)<<std::endl;
        return gain;
    }
    void batch(streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate,
               bool mrk,size_t round,size_t count){
        auto started=std::chrono::steady_clock::now();
        auto& offset=mrk?mrk_offset:collinear_offset;
        auto& stagnant=mrk?mrk_stagnant:collinear_stagnant;
        auto& done=mrk?mrk_done:collinear_done;
        size_t maximum=mrk?policy.max_mrk_samples:policy.max_collinear_samples;
        if(done||offset>=maximum){done=true;return;}
        count=std::min(count,maximum-offset);
        irr_physics::direct_coordinate_space(c);c.seed=seed;c.samples=count;c.sample_offset=offset;
        auto equations=mrk?irr_physics::direct_slots(c,true):irr_physics::direct_collinear_tail_exact(c);
        size_t gain=admit(reducer,coordinate,std::move(equations),mrk?"mrk":"collinear_tail_exact",round,offset,count);
        #ifdef USE_MIMALLOC
        mi_collect(true);
#endif
        offset+=count;stagnant=gain?0:stagnant+1;
        auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::cout<<"ADAPTIVE_PHYSICS_COST operator="<<(mrk?"mrk":"collinear_tail_exact")
                 <<" seconds="<<seconds<<" rank_gain="<<gain<<" gain_per_second="<<gain/std::max(seconds,1e-9)<<std::endl;
        done=offset>=maximum||stagnant>=policy.stagnant_batches;
        // t is one global normalization coordinate, never an irrep copy.
        reducer.prune_zeros();irr::require(!reducer.is_zero(columns),"Adaptive physical RHS is inconsistent");
    }
public:
    DirectSchedule(irr_physics::context& context,uint64_t sample_seed,Policy value={})
        :c(context),seed(sample_seed),policy(value),columns(irr_physics::invariant_columns(context)){
        irr::require(c.k<=3,"Adaptive exact-tail physical schedule needs k<=3");
        irr::require(policy.samples_per_batch&&policy.mrk_samples_per_batch&&policy.collinear_batches_per_cycle&&policy.stagnant_batches&&policy.cycles_before_ope&&policy.callback_seconds>0&&policy.maximum_cycles_per_callback,"Invalid adaptive physical batch policy");
        irr::require(policy.target_augmented_nullity&&policy.target_augmented_nullity<=columns+1,"Invalid adaptive physical target");
    }
    void operator()(streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate,size_t round){
        irr::require(coordinate.size()==columns,"Adaptive physical coordinate permutation width");
        if(round<policy.light_rounds){
            // The same disjoint deterministic streams continue into adaptive
            // batches. There is no replay of samples 0..light_rounds-1.
            batch(reducer,coordinate,false,round,1);batch(reducer,coordinate,true,round,1);return;
        }
        if(nullity(reducer)>policy.activate_nullity)return;
        auto started=std::chrono::steady_clock::now();
        auto within_budget=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<policy.callback_seconds;};
        if(!activated){activated=true;
            std::cout<<"ADAPTIVE_PHYSICS_ACTIVATE round="<<round<<" augmented_nullity="<<nullity(reducer)
                     <<" target="<<policy.target_augmented_nullity<<" batch_samples="<<policy.samples_per_batch
                     <<" collinear_max_samples="<<policy.max_collinear_samples
                     <<" mrk_max_samples="<<policy.max_mrk_samples<<std::endl;
        }
        size_t callback_cycles=0;
        while(callback_cycles<policy.maximum_cycles_per_callback&&nullity(reducer)>policy.target_augmented_nullity&&(!collinear_done||!mrk_done)&&within_budget()){
            for(size_t b=0;b<policy.collinear_batches_per_cycle&&!collinear_done
                             &&nullity(reducer)>policy.target_augmented_nullity&&within_budget();++b)
                batch(reducer,coordinate,false,round,policy.samples_per_batch);
            if(!mrk_done&&nullity(reducer)>policy.target_augmented_nullity&&within_budget())
                batch(reducer,coordinate,true,round,policy.mrk_samples_per_batch);
            ++cycles;++callback_cycles;
            if(!ope_done&&cycles>=policy.cycles_before_ope)break;
        }
        if(nullity(reducer)>policy.target_augmented_nullity&&!ope_done&&policy.ope_target_cap&&cycles>=policy.cycles_before_ope&&within_budget()){
            // Generate leading known OPE rows once. Do not retain multi-GB
            // head caches alongside the large triangular sewing basis.
            auto cache_budget=c.ope_head_cache_limit;c.clear_ope_head_cache();c.ope_head_cache_limit=0;
            auto equations=irr_physics::direct_ope_small_s_lowgrade(c,true,0,policy.ope_target_cap);
            c.clear_ope_head_cache();c.ope_head_cache_limit=cache_budget;
            admit(reducer,coordinate,std::move(equations),"ope_lowgrade4",round,0,0);ope_done=true;
#ifdef USE_MIMALLOC
            mi_collect(true);
#endif
            reducer.prune_zeros();irr::require(!reducer.is_zero(columns),"Early OPE RHS is inconsistent");
        }
        std::cout<<"ADAPTIVE_PHYSICS_PAUSE round="<<round<<" augmented_nullity="<<nullity(reducer)
                 <<" collinear_samples="<<collinear_offset<<" mrk_samples="<<mrk_offset
                 <<" collinear_stagnant_batches="<<collinear_stagnant<<" mrk_stagnant_batches="<<mrk_stagnant
                 <<" collinear_done="<<collinear_done<<" mrk_done="<<mrk_done<<std::endl;
        // Stagnation and sample caps are discovery heuristics. They do not
        // certify a full functional physical operator. Final physical checks
        // and the complete original sewing sweep remain mandatory.
    }
};
} // namespace hexagon_adaptive_physics
