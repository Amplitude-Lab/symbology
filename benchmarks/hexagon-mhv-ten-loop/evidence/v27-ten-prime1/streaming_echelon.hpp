#pragma once
#include "symrep.hpp"
#include <queue>
#include "structural_replay.hpp"
#include "wide_modular_accumulator.hpp"

// Exact finite-field elimination with bounded input batches. Retain only an
// independent triangular row basis; redundant input rows never accumulate in
// the global Schur complement. Every input equation is processed.
namespace streaming_echelon {
using namespace symrep;
using symrep::require;
using ModVec=sparse_vec<ulong,I>;
using ModMat=sparse_mat<ulong,I>;
struct Substitution {std::vector<I> image;std::vector<ulong> scale;};
// Nested, deterministic equation samples. Full operator verification, not
// the sampling probability, decides acceptance of a discovered kernel.
inline bool sampled_equation(size_t id,unsigned shift,uint64_t salt){
    if(!shift)return true;uint64_t x=uint64_t(id)+salt;x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);x^=x>>31;
    return !(x&((UINT64_C(1)<<shift)-1));
}
class Reducer {
    struct Work {std::vector<ulong> dense;std::vector<unsigned char> flags;std::vector<I> touched,heap;
        size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;uint64_t coefficient_updates=0,queue_pushes=0;
        std::vector<__uint128_t> wide;uint64_t wide_updates=0,wide_overflows=0;
        workload::BitQueue bits;
        explicit Work(size_t n):dense(n),flags(n),bits(n){};};
    ModMat basis,pending;
    structural_replay::Trace* trace=nullptr;
    double reduce_cpu=0,forward_cpu=0,verify_cpu=0,back_cpu=0;
    size_t independent_inputs=0,dependent_inputs=0;
    static double cpu_now(){return double(std::clock())/CLOCKS_PER_SEC;}
    std::vector<I> pivot_columns,pivot_row;
    std::vector<unsigned char> known_zero;
    Substitution substitutions;
    std::vector<Work> work;
    field_t fp;
    rref_option* opt;
    size_t last_pruned_rank=0,zero_count=0;
    size_t row_limit,nnz_limit,pending_nnz=0,input_rows=0,input_nnz=0,batches=0,retained_nnz=0;
    size_t pivot_capacity_before=0,pivot_capacity_after=0;
    void reduce_wide(ModVec& row,Work& w,structural_replay::Route* learn,const structural_replay::Route* route,size_t edge_limit)const{
        using Wide=workload::WideModularAccumulator::Wide;
        workload::WideModularAccumulator accumulator(fp.mod.n);
        if(w.wide.empty())w.wide.resize(basis.ncol);
        bool learning=bool(learn),routed=route&&route->complete;if(learn)learn->complete=true;
        auto touch=[&](I c){if(!(w.flags[c]&1)){w.flags[c]|=1;w.touched.push_back(c);}};
        auto enqueue=[&](I c){touch(c);if(pivot_row[c]>=0&&!(w.flags[c]&2)){w.flags[c]|=2;++w.queue_pushes;w.bits.push(pivot_row[c]);}};
        for(auto[c,x]:row)if(substitutions.image[c]>=0){I j=substitutions.image[c];
            ulong y=substitutions.scale[c]==1?x:nmod_mul(x,substitutions.scale[c],fp.mod);
            accumulator.add(w.wide[j],y);if(routed)touch(j);else enqueue(j);
        }
        auto eliminate=[&](I r,bool queued){
            I pivot=pivot_columns[r];ulong x=accumulator.reduce(w.wide[pivot]);w.wide[pivot]=0;
            if(!x)return;x=fp.mod.n-x;
            for(auto[c,y]:basis[r]){if(accumulator.add(w.wide[c],Wide(x)*y))++w.wide_overflows;if(queued)enqueue(c);else touch(c);}
            // Normalized pivot coefficient is one. Its own update is exactly
            // cancelled by the original coefficient modulo the field.
            w.wide[pivot]=0;w.coefficient_updates+=basis[r].nnz();w.wide_updates+=basis[r].nnz();
        };
        if(routed){++w.route_hits;route->each([&](I r){if(size_t(r)<basis.nrow){++w.route_steps;eliminate(r,false);}});
            size_t touched=w.touched.size();for(size_t i=0;i<touched;++i){I c=w.touched[i];if(w.wide[c]&&pivot_row[c]>=0)enqueue(c);}
            if(!w.bits.empty())++w.route_fallbacks;
        }
        while(!w.bits.empty()){
            I r=I(w.bits.pop());++w.heap_steps;
            if(learning)learning=learn->append(r,edge_limit,basis[r].nnz());
            eliminate(r,true);
        }
        if(learn){learn->pivots.shrink_to_fit();learn->packed.shrink_to_fit();}
        row.zero();for(I c:w.touched){ulong x=w.wide[c]?accumulator.reduce(w.wide[c]):0;
            if(x){require(pivot_row[c]<0,"Wide triangular reduction failed");row.push_back(c,x);}w.wide[c]=0;w.flags[c]=0;
        }
        w.touched.clear();row.sort_indices();row.compress();
    }
    void reduce(ModVec& row,Work& w,structural_replay::Route* learn=nullptr,
                const structural_replay::Route* route=nullptr,size_t edge_limit=0)const{
        bool routed=route&&route->complete;
        // Delay modular reduction only where overflow is infrequent. Route
        // replay still scans ALL residual pivots and falls back exactly.
        if(workload::kernels()&&fp.mod.n<=(uint64_t(1)<<61)){reduce_wide(row,w,learn,route,edge_limit);return;}
        bool learning=bool(learn);if(learn){learn->complete=true;learn->pivots.clear();}
        auto touch=[&](I c){if(!(w.flags[c]&1)){w.flags[c]|=1;w.touched.push_back(c);}};
        auto enqueue=[&](I c){touch(c);if(pivot_row[c]>=0&&!(w.flags[c]&2)){
            w.flags[c]|=2;++w.queue_pushes;if(workload::kernels())w.bits.push(pivot_row[c]);else{w.heap.push_back(pivot_row[c]);std::push_heap(w.heap.begin(),w.heap.end(),std::greater<I>());}}};
        for(auto[c,x]:row)if(substitutions.image[c]>=0){I j=substitutions.image[c];ulong y=substitutions.scale[c]==1?x:nmod_mul(x,substitutions.scale[c],fp.mod);
            w.dense[j]=nmod_add(w.dense[j],y,fp.mod);if(routed)touch(j);else enqueue(j);}
        if(routed){++w.route_hits;
            route->each([&](I r){if(size_t(r)>=basis.nrow)return;++w.route_steps;
                ulong x=w.dense[pivot_columns[r]];if(!x)return;ulong shoup=n_mulmod_precomp_shoup(x,fp.mod.n);
                for(auto[c,y]:basis[r]){ulong value=_nmod_sub(w.dense[c],n_mulmod_shoup(x,y,shoup,fp.mod.n),fp.mod);if(value)touch(c);w.dense[c]=value;}w.coefficient_updates+=basis[r].nnz();});
            // Routing is only a performance hint: coefficients left on ANY
            // pivot are eliminated by the ordinary heap. This is an exact
            // local fallback even when a cancellation changes at a new prime.
            size_t touched=w.touched.size();for(size_t i=0;i<touched;++i){I c=w.touched[i];if(w.dense[c]&&pivot_row[c]>=0)enqueue(c);}
            if(workload::kernels()?!w.bits.empty():!w.heap.empty())++w.route_fallbacks;
        }
        while(workload::kernels()?!w.bits.empty():!w.heap.empty()){
            I r;if(workload::kernels())r=I(w.bits.pop());else{std::pop_heap(w.heap.begin(),w.heap.end(),std::greater<I>());r=w.heap.back();w.heap.pop_back();}++w.heap_steps;
            if(learning)learning=learn->append(r,edge_limit,basis[r].nnz());
            ulong x=w.dense[pivot_columns[r]];if(!x)continue;ulong shoup=n_mulmod_precomp_shoup(x,fp.mod.n);
            for(auto[c,y]:basis[r]){ulong value=_nmod_sub(w.dense[c],n_mulmod_shoup(x,y,shoup,fp.mod.n),fp.mod);if(value)enqueue(c);w.dense[c]=value;}w.coefficient_updates+=basis[r].nnz();
        }
        if(learn){learn->pivots.shrink_to_fit();learn->packed.shrink_to_fit();}
        row.zero();for(I c:w.touched){if(w.dense[c]){require(pivot_row[c]<0,"streaming triangular reduction failed");row.push_back(c,w.dense[c]);}w.dense[c]=0;w.flags[c]=0;}
        w.touched.clear();row.sort_indices();row.compress();
    }
    void flush(){
        if(pending.rows.empty())return;
        pending.nrow=pending.rows.size();
        structural_replay::Batch learned;
        const structural_replay::Batch* replay=nullptr;
        ModMat discarded(0,basis.ncol);
        if(trace){
            learned.inputs=pending.nrow;learned.support=structural_replay::support_hash(pending);
            if(trace->replay){
                structural_replay::guard(trace->cursor<trace->batches.size(),"Replay extra batch");
                replay=&trace->batches[trace->cursor];
                structural_replay::guard(replay->inputs==learned.inputs&&replay->support==learned.support,"Replay source support/order changed");
            }
        }
        double reduce_started=cpu_now();
        size_t edge_limit=0;
        if(trace&&!trace->replay&&basis.nrow){
            size_t available=workload::routes()?(8ULL<<20):trace->max_route_edges*sizeof(I)-trace->route_bytes;
            size_t headers=pending.nrow*sizeof(structural_replay::Route);
            if(available>headers)edge_limit=std::min(size_t(16384),(available-headers)/(workload::routes()?1:sizeof(I))/pending.nrow);
            if(edge_limit)learned.routes.resize(pending.nrow);
        }
        if(basis.nrow){auto action=[&](size_t i){
            reduce(pending[i],work[workload::slot()],learned.routes.empty()?nullptr:&learned.routes[i],
                replay&&i<replay->routes.size()?&replay->routes[i]:nullptr,edge_limit);
        };
        if(workload::kernels())workload::balanced(pending.nrow,&opt->pool,[&](size_t i){
            uint64_t n=pending[i].nnz();for(auto[c,x]:pending[i]){I r=pivot_row[c];if(r>=0)n+=basis[r].nnz();}return n;
        },action);else{opt->pool.detach_loop(size_t(0),pending.nrow,action);opt->pool.wait();}}
        if(trace&&!trace->replay){
            if(workload::routes())trace->retain_routes(learned);
            else{trace->route_bytes+=learned.routes.size()*sizeof(structural_replay::Route);for(const auto& r:learned.routes){trace->route_edges+=r.pivots.size();trace->route_bytes+=r.bytes();}}
        }

        reduce_cpu+=cpu_now()-reduce_started;
        std::vector<size_t> source_order;
        if(replay){
            ModMat selected(replay->selected.size(),basis.ncol);std::vector<unsigned char> used(pending.nrow);
            for(size_t i=0;i<selected.nrow;++i){auto r=replay->selected[i];used[r]=1;selected[i]=std::move(pending[r]);}
            for(size_t r=0;r<pending.nrow;++r)if(!used[r]&&pending[r].nnz())discarded.rows.push_back(std::move(pending[r]));
            discarded.nrow=discarded.rows.size();pending=std::move(selected);
        }else{
            for(size_t r=0;r<pending.nrow;++r)if(pending[r].nnz())source_order.push_back(r);
            std::stable_sort(source_order.begin(),source_order.end(),[&](size_t x,size_t y){return pending[x].nnz()<pending[y].nnz();});
            ModMat sorted(source_order.size(),basis.ncol);
            for(size_t i=0;i<source_order.size();++i)sorted[i]=std::move(pending[source_order[i]]);pending=std::move(sorted);
        }
        if(pending.nrow){
            double forward_started=cpu_now();
            auto pivots=replay?structural_replay::forward(pending,*replay,fp,opt):sparse_mat_rref_forward(pending,fp,opt);
            forward_cpu+=cpu_now()-forward_started;
            require(!opt->abort,"streaming modular kernel aborted");
            if(trace&&!replay){learned.groups=pivots;
                for(auto& group:learned.groups)for(auto& p:group){learned.selected.push_back(source_order[p.r]);p.r=index(learned.selected.size()-1);}}

            size_t first=basis.rows.size();
            for(const auto& batch:pivots)for(auto [r,c]:batch){
                // Forward elimination can leave geometric-growth capacity in
                // a retained pivot. Shrink before transferring ownership. The
                // V11 mapped allocator returns the unused tail to its arena.
                size_t count=pending[r].nnz();pivot_capacity_before+=pending[r].alloc();
                if(pending[r].alloc()!=count)pending[r].compress();
                require(pending[r].nnz()<=count,"Pivot compression added stored coefficients");
                pivot_capacity_after+=pending[r].alloc();
                if(replay)structural_replay::guard(pivot_row[c]<0,"Replay duplicate pivot");
                require(pivot_row[c]<0,"duplicate streaming pivot");pivot_row[c]=index(basis.rows.size());pivot_columns.push_back(c);retained_nnz+=pending[r].nnz();basis.rows.push_back(std::move(pending[r]));
            }
            basis.nrow=basis.rows.size();
            independent_inputs+=basis.nrow-first;
            for(size_t r=first;r<basis.nrow;++r){
                if(replay){
                    auto x=basis[r].find(pivot_columns[r]);structural_replay::guard(x&&*x==1,"Replay pivot normalization changed");
                    for(auto c:basis[r].index_span())structural_replay::guard(pivot_row[c]<0||pivot_row[c]>=index(r),"Replay triangular dependency changed");
                }
                require(basis[r].find(pivot_columns[r])&&*basis[r].find(pivot_columns[r])==1,"streaming pivot is not normalized");
                for(auto c:basis[r].index_span())require(pivot_row[c]<0||pivot_row[c]>=index(r),"streaming basis is not triangular");
            }
        }
        if(replay){
            double verify_started=cpu_now();dependent_inputs+=discarded.nrow;
            opt->pool.detach_loop(size_t(0),discarded.nrow,[&](size_t i){reduce(discarded[i],work[workload::slot()]);});opt->pool.wait();
            for(const auto& row:discarded.rows)structural_replay::guard(!row.nnz(),"Replay discarded equation gained rank or changed affine RHS");
            verify_cpu+=cpu_now()-verify_started;
        }
        if(trace){if(!trace->replay)trace->batches.push_back(std::move(learned));++trace->cursor;}
        pending.clear();pending=ModMat(0,basis.ncol);pending_nnz=0;++batches;
        if((batches==1||batches%16==0)&&basis.nrow!=last_pruned_rank)prune_zeros();
        if(batches==1||batches%32==0)std::cout<<"streaming_batch="<<batches<<" input_rows="<<input_rows<<" rank="<<basis.nrow<<" retained_nnz="<<retained_nnz
            <<" retained_alloc="<<basis.alloc()<<" retained_array_bytes="<<basis.alloc()*(sizeof(ulong)+sizeof(I))
            <<" pivot_capacity_freed_slots="<<pivot_capacity_before-pivot_capacity_after<<std::endl;
    }
public:
    Reducer(size_t columns,const field_t& field,rref_option_t options,size_t rows=32768,size_t nonzeros=16*1024*1024,structural_replay::Trace* plan=nullptr)
        :basis(0,columns),pending(0,columns),trace(plan),pivot_row(columns,-1),known_zero(columns,0),fp(field),opt(options),row_limit(rows),nnz_limit(nonzeros){
        if(trace){if(trace->replay)structural_replay::guard(trace->columns==columns,"Replay column count changed");else trace->columns=columns;}
        require(rows>0&&nonzeros>0,"streaming batch limits must be positive");
        substitutions.image.resize(columns);std::iota(substitutions.image.begin(),substitutions.image.end(),0);substitutions.scale.assign(columns,1);
        for(size_t t=0;t<opt->pool.get_thread_count();++t)work.emplace_back(columns);
    }
    void consume(ModVec&& row){if(!row.nnz())return;++input_rows;input_nnz+=row.nnz();pending_nnz+=row.nnz();pending.rows.push_back(std::move(row));
        if(trace&&trace->replay){
            structural_replay::guard(trace->cursor<trace->batches.size(),"Replay extra source rows");
            if(pending.rows.size()>=trace->batches[trace->cursor].inputs)flush();
        }else if(pending.rows.size()>=row_limit||pending_nnz>=nnz_limit)flush();}
    // A global sparse front gives pivot selection the cheapest equations
    // across the entire operator before more expensive rows are streamed.
    void seed(ModMat&& matrix){
        require(!basis.nrow&&pending.rows.empty()&&matrix.ncol==basis.ncol,"streaming seed requires an empty compatible reducer");
        input_rows=matrix.nrow;input_nnz=matrix.nnz();pending_nnz=input_nnz;pending=std::move(matrix);flush();
    }
    size_t input_nonzeros()const{return input_nnz;}
    uint64_t arithmetic_work()const{uint64_t n=0;for(auto& w:work)n+=w.coefficient_updates;return n;}
    workload::GainWork competing_work;
    size_t retained_nonzeros()const{return retained_nnz;}
    size_t rank()const{return basis.nrow;}
    size_t completed_batches()const{return batches;}
    bool is_zero(size_t column)const{return known_zero[column];}
    const Substitution& substitution()const{return substitutions;}
    size_t prune_zeros(){
        if(last_pruned_rank==basis.nrow)return zero_count;
        // Relations only change when a newly retained row can introduce a
        // singleton or doubleton. New rows contain no earlier pivot columns;
        // old substitutions therefore leave every new row of length >2
        // unchanged. Without a new short row, the old full sweep is a no-op.
        bool new_short=false;
        for(size_t r=last_pruned_rank;r<basis.nrow;++r)if(basis[r].nnz()<=2){new_short=true;break;}
        if(!new_short){last_pruned_rank=basis.nrow;return zero_count;}
        // Propagate both singleton zeros and doubleton identifications in
        // reverse triangular order. Each rewritten equation differs by a
        // combination of retained later pivots, preserving the exact row space.
        size_t count=0;retained_nnz=0;
        for(size_t r=basis.nrow;r-->0;){auto& row=basis[r];size_t out=0;bool changed=false;
            I pivot=pivot_columns[r];for(size_t p=0;p<row.nnz();++p){I c=row(p);
                if(c==pivot){if(out!=p){row(out)=c;row[out]=row[p];}++out;}
                else {changed|=substitutions.image[c]!=c||substitutions.scale[c]!=1;
                    if(substitutions.image[c]>=0){
                        if(out!=p||substitutions.image[c]!=c)row(out)=substitutions.image[c];
                        if(out!=p||substitutions.scale[c]!=1)
                            row[out]=substitutions.scale[c]==1?row[p]:nmod_mul(row[p],substitutions.scale[c],fp.mod);
                        ++out;}}}
            // Most retained rows do not change. Sorting and reallocating them
            // again costs more than the short-relation propagation itself.
            if(changed){row.resize(out);row.sort_indices();out=0;
            for(size_t first=0;first<row.nnz();){size_t end=first+1;I c=row(first);ulong x=row[first];while(end<row.nnz()&&row(end)==c)x=nmod_add(x,row[end++],fp.mod);if(x){row(out)=c;row[out]=x;++out;}first=end;}
            row.resize(out);row.compress();}require(row.find(pivot)&&*row.find(pivot)==1,"lost streaming pivot during short-relation propagation");
            substitutions.image[pivot]=pivot;substitutions.scale[pivot]=1;
            if(out==1){known_zero[pivot]=1;substitutions.image[pivot]=-1;substitutions.scale[pivot]=0;}
            else if(out==2){size_t p=row(0)==pivot?1:0;substitutions.image[pivot]=row(p);substitutions.scale[pivot]=nmod_neg(row[p],fp.mod);}
            retained_nnz+=out;
        }
        for(auto z:known_zero)count+=z;last_pruned_rank=basis.nrow;zero_count=count;return count;
    }
    void flush_pending(){flush();}
    std::pair<ModMat,std::vector<std::vector<pivot_t<I>>>> finish(){
        flush();if(trace)trace->complete();std::vector<std::vector<pivot_t<I>>> pivots(2);
        for(size_t r=0;r<basis.nrow;++r)pivots[1].emplace_back(index(r),pivot_columns[r]);
        size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;
        for(const auto& w:work){heap_steps+=w.heap_steps;route_steps+=w.route_steps;route_hits+=w.route_hits;route_fallbacks+=w.route_fallbacks;}
        std::cout<<"DEPENDENCY_PROFILE heap_steps="<<heap_steps<<" route_steps="<<route_steps<<" route_hits="<<route_hits
            <<" local_fallbacks="<<route_fallbacks<<" stored_edges="<<(trace?trace->route_edges:0)
            <<" stored_route_bytes="<<(trace?trace->route_bytes:0)<<" route_payload_cap_bytes="<<(trace?trace->max_route_edges*sizeof(I):0)<<std::endl;
        uint64_t updates=0,pushes=0;for(auto& w:work){updates+=w.coefficient_updates;pushes+=w.queue_pushes;}
        std::cout<<"KERNEL_WORK coefficient_updates="<<updates<<" queue_pushes="<<pushes<<" bitqueue="<<workload::kernels()<<std::endl;
        uint64_t wide_updates=0,overflows=0;for(auto& w:work){wide_updates+=w.wide_updates;overflows+=w.wide_overflows;}
        std::cout<<"WIDE_KERNEL updates="<<wide_updates<<" exact_overflow_corrections="<<overflows<<std::endl;
        work.clear();
        std::cout<<"streaming_complete rows="<<input_rows<<" input_nnz="<<input_nnz<<" rank="<<basis.nrow<<" retained_nnz="<<basis.nnz()
                 <<" retained_alloc="<<basis.alloc()<<" retained_array_bytes="<<basis.alloc()*(sizeof(ulong)+sizeof(I))
                 <<" pivot_capacity_freed_slots="<<pivot_capacity_before-pivot_capacity_after<<std::endl;
        double back_started=cpu_now();triangular_solver_2(basis,pivots,fp,opt);back_cpu+=cpu_now()-back_started;require(!opt->abort,"streaming back substitution aborted");
        std::cout<<"STRUCTURAL_PROFILE replay="<<bool(trace&&trace->replay)<<" reduce_cpu="<<reduce_cpu
            <<" forward_cpu="<<forward_cpu<<" discarded_verify_cpu="<<verify_cpu<<" back_cpu="<<back_cpu
            <<" independent_inputs="<<independent_inputs<<" verified_dependent_inputs="<<dependent_inputs<<std::endl;
        return {std::move(basis),std::move(pivots)};
    }
};
}
