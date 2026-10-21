#include "shared/streaming_echelon.hpp"
#include "shared/exact_operator_compression.hpp"
#include <random>
#include <set>
#include <sstream>
using Dense=exact_operator::Dense;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
size_t gauss(Dense& a,const field_t& F){
    if(a.empty())return 0;size_t rank=0;
    for(size_t j=0;j<a[0].size()&&rank<a.size();++j){size_t r=rank;while(r<a.size()&&!a[r][j])++r;if(r==a.size())continue;
        std::swap(a[r],a[rank]);ulong inv=n_invmod(a[rank][j],F.mod.n);for(auto& x:a[rank])x=nmod_mul(x,inv,F.mod);
        for(size_t i=0;i<a.size();++i)if(i!=rank&&a[i][j]){ulong x=a[i][j];for(size_t k=j;k<a[i].size();++k)a[i][k]=nmod_sub(a[i][k],nmod_mul(x,a[rank][k],F.mod),F.mod);}++rank;
    }a.resize(rank);return rank;
}
int main()try{
    Flint::set_memory_functions();std::mt19937 rng(7102026);size_t queue_ops=0,compression_cases=0,schedule_cases=0;
    {
        size_t overflows=0;
        for(ulong p:{101UL,103UL,1152921504606847009UL,1152921504606847067UL,18446744073709551557UL}){
            field_t F(FIELD_Fp,p);workload::WideModularAccumulator kernel(p);__uint128_t value=0;ulong reference=0;
            for(size_t i=0;i<10000;++i){ulong a=p-1-(rng()%97),b=p-1-(rng()%97);
                overflows+=kernel.add(value,__uint128_t(a)*b);reference=nmod_add(reference,nmod_mul(a,b,F.mod),F.mod);
                if(i%37==0)check(kernel.reduce(value)==reference,"Delayed modular accumulation differs from canonical arithmetic");
            }check(kernel.reduce(value)==reference,"Delayed accumulation final residue");
        }check(overflows>0,"Accumulator overflow path untested");
        std::cout<<"WIDE_CONTROL overflow_cases="<<overflows<<" canonical_reference_equal=1"<<std::endl;
    }
    for(size_t n:{size_t(1),size_t(63),size_t(64),size_t(65),size_t(4097),size_t(105899),size_t(262145)}){
        workload::BitQueue q(n);std::set<size_t> ref;
        for(size_t i=0;i<100000;++i){if(ref.empty()||rng()%3){size_t j=rng()%n;q.push(j);ref.insert(j);}else{check(q.pop()==*ref.begin(),"Bit queue order");ref.erase(ref.begin());}check(q.empty()==ref.empty(),"Bit queue empty");++queue_ops;}
        while(!ref.empty()){check(q.pop()==*ref.begin(),"Bit queue drain");ref.erase(ref.begin());++queue_ops;}check(q.empty(),"Bit queue retained state");
    }
    rref_option_t opt;opt->pool.reset(4);
    {
        std::vector<ulong> periodic(512);for(size_t i=0;i<periodic.size();i+=16)periodic[i]=1;
        check(!workload::prefer_dense(periodic.data(),periodic.size()),"Density sampling alias");
        std::fill(periodic.begin(),periodic.end(),1);check(workload::prefer_dense(periodic.data(),periodic.size()),"Dense dispatch");
        workload::CycleProgress cycle(2);
        for(size_t callback=0;callback<3;++callback)cycle.admit(0,{false,false});
        check(cycle.completed==0,"Premature cycle");cycle.admit(1,{false,false});check(cycle.completed==1,"Cycle lost across budgets");
        cycle.admit(1,{true,false});check(cycle.completed==2,"Finished task blocks cycle");
        workload::GainWork task,seam;task.observe(7,1000);seam.observe(30,1000);
        check(!workload::competitive(task,seam),"Unproductive early work admitted");
        task.observe(70,1000);check(workload::competitive(task,seam),"Productive work deferred");
        seam.observe(0,1);check(workload::competitive(task,seam),"Stagnant alternative received artificial rank gain");
    }
    for(size_t n:{size_t(0),size_t(1),size_t(37),size_t(1000)}){
        std::vector<size_t> x(n);workload::balanced(n,&opt->pool,[&](size_t i){return i%17?1:100000;},[&](size_t i){++x[i];});
        for(auto v:x)check(v==1,"Balanced ownership");
    }
    bool caught=false;try{workload::balanced(1000,&opt->pool,[](size_t){return 10000;},[](size_t i){if(i==51)throw std::runtime_error("expected");});}catch(const std::runtime_error&){caught=true;}check(caught,"Worker exception swallowed");
    for(ulong p:{101UL,103UL,1152921504606847009UL,1152921504606847067UL})for(size_t n:{size_t(1),size_t(7),size_t(24)})for(size_t cols:{size_t(1),size_t(3),size_t(31)})for(size_t family=0;family<4;++family){
        field_t F(FIELD_Fp,p);Dense a(n,std::vector<ulong>(cols+3));
        for(auto& row:a)for(auto& v:row)v=rng()%p;
        if(family==0)for(auto& row:a)std::fill(row.begin(),row.end(),0);
        if(family==1)for(size_t r=1;r<n;++r)for(size_t c=0;c<cols+3;++c)a[r][c]=nmod_mul(r%p,a[0][c],F.mod);
        if(family==2){for(auto& row:a)std::fill(row.begin(),row.end(),0);a.back().back()=1;}
        if(family==3)for(size_t r=n/2;r<n;++r)a[r]=a[r/2];
        auto basis=exact_operator::compress(a,F);Dense reference=a;gauss(reference,F);Dense selected;for(auto r:basis.selected)selected.push_back(a[r]);gauss(selected,F);check(selected==reference,"Compressed augmented row space");
        // An arbitrary expensive linear projection and three affine columns.
        Dense projected(n,std::vector<ulong>(12)),map(cols+3,std::vector<ulong>(12));for(auto& row:map)for(auto& x:row)x=rng()%p;
        for(size_t r=0;r<n;++r)for(size_t j=0;j<cols+3;++j)for(size_t k=0;k<12;++k)projected[r][k]=nmod_add(projected[r][k],nmod_mul(a[r][j],map[j][k],F.mod),F.mod);
        Dense after;for(auto r:basis.selected)after.push_back(projected[r]);gauss(after,F);gauss(projected,F);check(after==projected,"Preprojection compression changed operator");++compression_cases;
    }
    {
        field_t F(FIELD_Fp,101),G(FIELD_Fp,103);
        check(exact_operator::compress({{1,1},{1,1}},F).selected.size()==1,"Exceptional compression rank loss");
        check(exact_operator::compress({{1,1},{1,102}},G).selected.size()==2,"Exceptional compression rank gain");
        bool bad=false;try{exact_operator::compress({{1,2},{1}},F);}catch(const std::runtime_error&){bad=true;}check(bad,"Ragged compression input");
    }
    for(size_t pattern=0;pattern<64;++pattern){
        auto run=[&](size_t fake_clock){std::vector<size_t> order;std::vector<workload::GainWork> cost(3);size_t previous=SIZE_MAX,streak=0;
            for(size_t t=0;t<48;++t){std::vector<bool> enabled={true,true,bool(pattern&1)};auto chosen=workload::choose(cost,enabled,previous,streak);check(chosen<3&&enabled[chosen],"Ineligible scheduling");order.push_back(chosen);
                cost[chosen].observe((pattern+chosen*7+t)%11,1+(pattern*97+chosen*31+t*13)%1000);
                streak=previous==chosen?streak+1:1;previous=chosen;fake_clock+=rng()%10000;
            }return order;};
        check(run(0)==run(100000),"Schedule depends on time");++schedule_cases;
    }
    {
        structural_replay::Trace t;t.columns=100000;t.max_route_edges=2048;
        for(size_t batch=0;batch<12;++batch){structural_replay::Batch b;b.inputs=20;b.routes.resize(20);
            for(size_t r=0;r<20;++r){auto& route=b.routes[r];route.complete=true;for(int j=0;j<100;++j)route.append(j*(r+1),128,(batch+1)*100);route.packed.shrink_to_fit();}
            t.retain_routes(b);t.batches.push_back(std::move(b));check(t.route_bytes<=8192,"Route budget exceeded");
        }
        check(t.batches.front().routes.empty()&&!t.batches.back().routes.empty(),"Cost-based route eviction");
        t.decision(1);t.decision(3);
        std::stringstream stream;t.write(stream,"packed");auto copy=structural_replay::Trace::read(stream,"packed");check(copy.route_bytes==t.route_bytes&&copy.route_edges==t.route_edges,"Packed roundtrip accounting");
        check(copy.decision(0)==1&&copy.peek_decision()==3&&copy.decision(3)==3,"Recorded task/stop replay");
        copy.cursor=copy.batches.size();copy.complete();
        for(size_t i=0;i<t.batches.size();++i)for(size_t r=0;r<t.batches[i].routes.size();++r){std::vector<int32_t> a,b;t.batches[i].routes[r].each([&](auto j){a.push_back(j);});copy.batches[i].routes[r].each([&](auto j){b.push_back(j);});check(a==b,"Packed route roundtrip");}
        auto broken=t;for(auto& b:broken.batches)if(!b.routes.empty()){b.routes[0].packed.back()|=128;break;}
        std::stringstream bad;broken.write(bad,"packed");bool caught=false;try{structural_replay::Trace::read(bad,"packed");}catch(const structural_replay::Mismatch&){caught=true;}check(caught,"Packed truncation guard");
    }
    {
        structural_replay::Trace t;t.columns=10000;t.max_route_edges=512;
        auto block=[](size_t bytes,uint64_t benefit){structural_replay::Batch b;b.inputs=1;b.routes.resize(1);auto& r=b.routes[0];r.complete=true;r.packed.assign(bytes,1);r.count=bytes;r.benefit=benefit;return b;};
        auto low=block(600,6000);t.retain_routes(low);t.batches.push_back(std::move(low));
        auto high=block(600,30000);t.retain_routes(high);t.batches.push_back(std::move(high));
        size_t before=t.route_bytes;auto rejected=block(1500,32000);t.retain_routes(rejected);
        check(rejected.routes.empty()&&t.route_bytes==before&&!t.batches[0].routes.empty()&&!t.batches[1].routes.empty(),"Rejected route block evicted useful cache");
    }
    std::cout<<"ADDITIONAL_CONTROL passed=1 queue_operations="<<queue_ops<<" compression_cases="<<compression_cases<<" deterministic_schedules="<<schedule_cases<<" packed_budget_and_damage_controls=1 exact_density_and_persistent_cycle_controls=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"CONTROL_ERROR "<<e.what()<<std::endl;return 1;}
