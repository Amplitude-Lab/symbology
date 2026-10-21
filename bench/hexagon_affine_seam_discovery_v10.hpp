#pragma once
#include "shared/seam_discovery.hpp"

namespace staged_kernel {
// Intersect a homogeneous augmented kernel with t=1. All coefficient
// coordinates share the same t. The origin is normalized once, placed last,
// and removed from every direction before the auxiliary coordinate is dropped.
inline MM affine_slice_augmented_kernel(MM kernel,const field_t& fp,
                                        const std::vector<I>& variable_order={}){
    require(kernel.ncol>0&&kernel.nrow>0,"Early physical equations have no affine solution");
    size_t columns=kernel.ncol-1,origin_row=SIZE_MAX;
    for(size_t r=0;r<kernel.nrow;++r)if(auto t=kernel[r].find(index(columns));t&&*t){origin_row=r;break;}
    require(origin_row!=SIZE_MAX,"Early physical equations force the affine normalization to zero");
    MV origin;ulong tvalue=*kernel[origin_row].find(index(columns)),inverse=n_invmod(tvalue,fp.mod.n);
    for(auto [j,x]:kernel[origin_row])if(size_t(j)<columns)
        origin.push_back(j,nmod_mul(x,inverse,fp.mod));
    normalize_mod(origin,fp);kernel[origin_row].clear();
    MM current(kernel.nrow,columns);size_t direction=0;
    for(size_t r=0;r<kernel.nrow;++r)if(r!=origin_row){
        ulong value=0;if(auto t=kernel[r].find(index(columns));t)value=*t;
        for(auto [j,x]:kernel[r])if(size_t(j)<columns)current[direction].push_back(j,x);
        if(value)for(auto [j,x]:origin)
            current[direction].push_back(j,nmod_neg(nmod_mul(value,x,fp.mod),fp.mod));
        normalize_mod(current[direction],fp);kernel[r].clear();++direction;
    }
    kernel.clear();require(direction+1==current.nrow,"Affine direction count mismatch");
    current[direction]=std::move(origin);
    if(!variable_order.empty()){
        require(variable_order.size()==columns,"Affine inverse permutation dimension");
        std::vector<I> inverse_order(columns);
        for(size_t c=0;c<columns;++c)inverse_order[variable_order[c]]=index(c);
        for(auto& row:current.rows){for(auto& c:row.index_span())c=inverse_order[c];row.sort_indices();}
    }
    return current;
}

// Early physical rows and sewing rows are admitted to one augmented reducer.
// The source sees only its original columns; the additional final column t
// is reserved for physical RHS values. physical(reducer, columns, round)
// emits [A,-b] rows in the same coefficient permutation and may use bounded
// batches. It is called after the cheapest seam seed and before every next
// seam refinement. Source masks apply only to sewing, never to physical rows.
template<class Source,class Physical>
MM discover_affine_candidate(const Source& source,const field_t& fp,
                             rref_option_t& options_rref,Options options,
                             Physical&& physical,const std::vector<I>& variable_order={}){
    auto opt=options_rref;auto pool=&opt->pool;
    const size_t columns_count=source.ncol(),augmented_count=columns_count+1;
    structured_kernel::Relations relations(columns_count);
    std::vector<I> columns(columns_count);std::iota(columns.begin(),columns.end(),0);
    if(!variable_order.empty()){
        require(variable_order.size()==columns.size(),"Affine variable permutation dimension");
        columns=variable_order;std::vector<unsigned char> used(columns.size());
        for(auto c:columns){require(c>=0&&size_t(c)<columns.size()&&!used[c],"Invalid affine variable permutation");used[c]=1;}
    }
    require(options.seed_terms&&options.seed_nonzeros,"Affine staged front limits must be positive");
    require(augmented_count<=size_t(std::numeric_limits<I>::max()),"Augmented affine coordinate range");
    auto start=std::chrono::steady_clock::now();
    auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    auto seen=source.masks(options.remember_equations||options.sample_equations);
    Selection selected{options.seed_terms};
    streaming_echelon::Reducer reducer(augmented_count,fp,opt,options.batch_rows,options.batch_nonzeros);
    MM a(0,augmented_count);size_t seed_nnz=0,seed_rows=0;bool seeded=false;
    source.discover(relations,columns,fp,selected,seen,[&](MV&& row){
        ++seed_rows;seed_nnz+=row.nnz();
        if(seeded)reducer.consume(std::move(row));
        else{a.rows.push_back(std::move(row));if(seed_nnz>=options.seed_nonzeros){
            a.nrow=a.rows.size();reducer.seed(std::move(a));seeded=true;}}
    },pool);
    if(!seeded){a.nrow=a.rows.size();reducer.seed(std::move(a));}
    else reducer.flush_pending();
    std::cout<<"REFINE_SEED rows="<<seed_rows<<" cols="<<columns_count<<" augmented_cols="<<augmented_count
             <<" nnz="<<seed_nnz<<" bounded="<<seeded<<" elapsed="<<elapsed()<<std::endl;
    size_t physical_round=0;
    auto inject_physical=[&]{
        size_t before=reducer.rank(),round=physical_round++;
        physical(reducer,columns,round);reducer.flush_pending();reducer.prune_zeros();
        require(!reducer.is_zero(columns_count),"Early physical equations are inconsistent with affine normalization");
        std::cout<<"JOINT_EARLY_PHYSICS round="<<round<<" before_rank="<<before<<" rank="<<reducer.rank()
                 <<" nullity="<<augmented_count-reducer.rank()<<" elapsed="<<elapsed()<<std::endl;
    };
    inject_physical();
    struct Enough{};size_t threshold=options.seed_terms;
    bool complete_discovery=threshold==SIZE_MAX,plateau=false,first_refinement=true;
    while(augmented_count-reducer.rank()>options.target_nullity&&!complete_discovery&&!plateau){
        if(!first_refinement)inject_physical();
        first_refinement=false;
        if(augmented_count-reducer.rank()<=options.target_nullity)break;
        size_t zeros=reducer.prune_zeros();auto active_columns=columns;
        for(auto& c:active_columns)if(c>=0&&reducer.is_zero(c))c=-1;
        selected.substitution=&reducer.substitution();
        if(threshold==SIZE_MAX)selected.sample_shift=selected.sample_shift>2?selected.sample_shift-2:0;
        else{threshold=threshold>=16384?SIZE_MAX:threshold*4;
            selected.sample_shift=options.sample_equations&&threshold>options.complete_through_terms?4:0;}
        selected.max_terms=threshold;
        std::cout<<"REFINE_ZERO_SUBSTITUTION columns="<<zeros<<" limit="<<threshold
                 <<" sample_denominator="<<(UINT64_C(1)<<selected.sample_shift)<<std::endl;
        bool stopped=false;size_t last_batch=reducer.completed_batches(),last_rank=reducer.rank(),stagnant=0;
        try{
            source.discover(relations,active_columns,fp,selected,seen,[&](MV&& row){
                reducer.consume(std::move(row));
                if(augmented_count-reducer.rank()<=options.target_nullity)throw Enough{};
                if(reducer.completed_batches()!=last_batch){
                    size_t gain=reducer.rank()-last_rank;
                    stagnant=gain<=options.stagnation_rank_gain?stagnant+1:0;
                    last_batch=reducer.completed_batches();last_rank=reducer.rank();
                    if(options.stagnation_batches&&stagnant>=options.stagnation_batches
                       &&augmented_count-reducer.rank()<=options.max_residual_nullity){
                        plateau=true;std::cout<<"REFINE_STAGNATION nullity="<<augmented_count-reducer.rank()
                            <<" batches="<<stagnant<<" elapsed="<<elapsed()<<std::endl;throw Enough{};
                    }
                }
            },pool);
        }catch(const Enough&){stopped=true;}
        reducer.flush_pending();
        complete_discovery=!stopped&&threshold==SIZE_MAX&&!selected.sample_shift;
        std::cout<<"REFINE_DISCOVERY limit="<<threshold<<" rank="<<reducer.rank()
                 <<" nullity="<<augmented_count-reducer.rank()<<" elapsed="<<elapsed()<<std::endl;
    }
    auto reduced=reducer.finish();
    auto kernel=consume_modular_kernel(reduced.first,reduced.second,fp);
    auto current=affine_slice_augmented_kernel(std::move(kernel),fp,variable_order);
    std::cout<<"JOINT_DISCOVERY directions="<<current.nrow-1<<" complete="<<complete_discovery
             <<" affine_origin=1 nnz="<<current.nnz()<<std::endl;
    return current;
}
} // namespace staged_kernel
