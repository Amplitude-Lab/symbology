#pragma once
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace workload {
inline bool option(const char* name) {
    const char* p=std::getenv(name);
    if(!p)return true;
    if(std::string(p)!="0"&&std::string(p)!="1")throw std::runtime_error(std::string("Invalid switch: ")+name);
    return *p=='1';
}
inline bool kernels(){static const bool x=option("SOLVER_WORK_KERNELS");return x;}
inline bool schedule(){static const bool x=option("SOLVER_WORK_SCHEDULE");return x;}
inline bool compression(){static const bool x=option("SOLVER_EXACT_COMPRESSION");return x;}
inline bool routes(){static const bool x=option("SOLVER_PACKED_ROUTES");return x;}
inline size_t slot(){return BS::this_thread::get_index().value_or(0);}
template<class Scalar>bool prefer_dense(const Scalar* values,size_t n){
    if(n<256)return false;
    size_t nonzero=std::count_if(values,values+n,[](Scalar x){return x!=0;});
    return nonzero>=n-n/8;
}

// Exact ordered integer set. Upper levels mark nonempty words below them.
// No matrix values or ordering assumptions enter the queue implementation.
class BitQueue {
    std::vector<std::vector<uint64_t>> level;
public:
    explicit BitQueue(size_t n){do{n=(n+63)/64;level.emplace_back(std::max(size_t(1),n));}while(n>1);}
    bool empty()const{return !level.back()[0];}
    void push(size_t x){
        for(auto& a:level){size_t j=x/64;uint64_t old=a[j];a[j]|=uint64_t(1)<<(x%64);if(old)break;x=j;}
    }
    size_t pop(){
        if(empty())throw std::runtime_error("Empty bit queue");
        size_t x=0;for(size_t k=level.size();k-->0;)x=64*x+std::countr_zero(level[k][x]);
        size_t result=x;
        for(auto& a:level){size_t j=x/64;a[j]&=~(uint64_t(1)<<(x%64));if(a[j])break;x=j;}
        return result;
    }
};

// Contiguous task ranges with approximately equal predicted arithmetic.
// Workers dynamically claim ranges; output ownership and arithmetic order
// within each row remain fixed. The pool must not already be in use.
template<class Estimate,class Function>
void balanced(size_t count,thread_pool* pool,Estimate estimate,Function function,size_t cap=SIZE_MAX){
    if(!count)return;
    size_t workers=pool?std::min({count,size_t(pool->get_thread_count()),cap}):1;
    std::vector<uint64_t> cost(count);uint64_t total=0;
    for(size_t i=0;i<count;++i){cost[i]=std::max(uint64_t(1),uint64_t(estimate(i)));total+=cost[i];}
    workers=std::min(workers,std::max(size_t(1),size_t(total/16384)));
    if(workers<=1){for(size_t i=0;i<count;++i)function(i);return;}
    uint64_t target=std::max(uint64_t(16384),total/(workers*8));
    std::vector<std::pair<size_t,size_t>> ranges;size_t first=0;uint64_t sum=0;
    for(size_t i=0;i<count;++i){sum+=cost[i];if(sum>=target){ranges.emplace_back(first,i+1);first=i+1;sum=0;}}
    if(first<count)ranges.emplace_back(first,count);
    std::atomic<size_t> next=0;std::atomic<bool> failed=false;std::vector<std::exception_ptr> errors(workers);
    pool->detach_loop(size_t(0),workers,[&](size_t worker){try{
        while(!failed){size_t q=next.fetch_add(1,std::memory_order_relaxed);if(q>=ranges.size())break;
            for(size_t i=ranges[q].first;i<ranges[q].second;++i)function(i);}
    }catch(...){errors[worker]=std::current_exception();failed=true;}});pool->wait();
    for(auto e:errors)if(e)std::rethrow_exception(e);
}

// Measurements here are integer work counters, never clock deadlines.
struct GainWork {
    uint64_t gain=1,work=1;size_t visits=0;
    void observe(uint64_t g,uint64_t w){gain=g;work=std::max(uint64_t(1),w);++visits;}
};
inline bool competitive(const GainWork& candidate,const GainWork& alternative){
    if(!candidate.visits||!alternative.visits)return true;
    if(!alternative.gain)return true;
    return __uint128_t(candidate.gain)*alternative.work>=__uint128_t(alternative.gain)*candidate.work;
}
struct CycleProgress {
    std::vector<bool> seen;size_t completed=0;
    explicit CycleProgress(size_t n):seen(n){}
    void admit(size_t task,const std::vector<bool>& finished){
        if(task>=seen.size()||finished.size()!=seen.size())throw std::runtime_error("Cycle task shape");
        seen[task]=true;bool done=true;for(size_t i=0;i<seen.size();++i)done=done&&(seen[i]||finished[i]);
        if(done){++completed;seen=finished;}
    }
};
inline size_t choose(const std::vector<GainWork>& costs,const std::vector<bool>& enabled,size_t previous,size_t streak){
    size_t best=SIZE_MAX;
    for(size_t i=0;i<costs.size();++i)if(enabled[i]){
        if(!costs[i].visits)return i;
        if(streak>=4&&i==previous&&std::count(enabled.begin(),enabled.end(),true)>1)continue;
        if(best==SIZE_MAX||(__uint128_t(costs[i].gain+1)*costs[best].work>
                           __uint128_t(costs[best].gain+1)*costs[i].work))best=i;
    }
    return best;
}
}
