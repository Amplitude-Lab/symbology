#!/usr/bin/env python3
"""Reproducible isolated candidate; never edit the V26 package."""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[2]; S=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--compile-only',action='store_true');a=p.parse_args();O=a.output.resolve()
parent=R/'output_hexagon_v23_20261009/package-v26'
def edit(name,old,new):
    q=O/name;s=q.read_text();assert old in s,(name,old[:100]);q.write_text(s.replace(old,new))
if not a.compile_only:
    assert not O.exists();shutil.copytree(parent,O,ignore=shutil.ignore_patterns('bin','.git','output','log'));(O/'bin').mkdir()
    for name in ['workload_kernels.hpp','exact_operator_compression.hpp','wide_modular_accumulator.hpp']:
        shutil.copy2(S/name,O/'src/shared'/name)
    q=O/'src/shared/structural_replay.hpp';s=q.read_text().replace('#include <stdexcept>','#include <stdexcept>\n#include "workload_kernels.hpp"')
    start=s.index('struct Route ');end=s.index('struct Batch ',start)
    s=s[:start]+'''struct Route {
    bool complete=false;std::vector<I> pivots;std::vector<unsigned char> packed;
    size_t count=0;uint64_t benefit=0;I last=-1;
    size_t bytes()const{return pivots.size()*sizeof(I)+packed.size();}
    bool append(I r,size_t limit,uint64_t work){
        if(!workload::routes()){
            if(pivots.size()>=limit){complete=false;pivots.clear();return false;}
            pivots.push_back(r);return true;
        }
        uint32_t delta=uint32_t(r-last),x=delta;size_t n=1;while(x>=128){x>>=7;++n;}
        if(packed.size()+n>limit)return false;
        while(delta>=128){packed.push_back((delta&127)|128);delta>>=7;}packed.push_back(delta);
        last=r;++count;benefit+=work;return true;
    }
    template<class Function>void each(Function f)const{
        if(packed.empty()){for(I x:pivots)f(x);return;}
        I previous=-1;uint32_t x=0;unsigned shift=0;
        for(auto b:packed){x|=uint32_t(b&127)<<shift;if(b&128)shift+=7;
            else{previous+=I(x);f(previous);x=0;shift=0;}}
    }
};
''' +s[end:]
    s=s.replace('struct Trace {','''struct Trace {
    std::vector<size_t> physical_decisions;size_t decision_cursor=0;
    size_t peek_decision()const{guard(replay&&decision_cursor<physical_decisions.size(),"Replay missing physical stop");return physical_decisions[decision_cursor];}
    size_t decision(size_t proposed){
        if(!replay){physical_decisions.push_back(proposed);return proposed;}
        guard(decision_cursor<physical_decisions.size(),"Replay extra physical decision");return physical_decisions[decision_cursor++];
    }
    static size_t batch_bytes(const Batch& b){size_t n=b.routes.size()*sizeof(Route);for(auto& r:b.routes)n+=r.bytes();return n;}
    static uint64_t batch_benefit(const Batch& b){uint64_t n=0;for(auto& r:b.routes)n+=r.benefit;return n;}
    void retain_routes(Batch& next){
        size_t bytes=batch_bytes(next);if(!bytes)return;
        uint64_t benefit=batch_benefit(next);size_t budget=max_route_edges*sizeof(I);
        auto drop=[](Batch& b){std::vector<Route>().swap(b.routes);};
        if(bytes>budget||!benefit){drop(next);return;}
        // Plan evictions before mutating the cache. A rejected large block
        // must not displace useful smaller blocks along the way.
        struct Victim{size_t index,bytes;uint64_t benefit;};std::vector<Victim> victims;
        if(route_bytes+bytes>budget){
            for(size_t i=0;i<batches.size();++i)if(!batches[i].routes.empty()){
                size_t old_bytes=batch_bytes(batches[i]);uint64_t old_benefit=batch_benefit(batches[i]);
                if(__uint128_t(benefit)*old_bytes>__uint128_t(old_benefit)*bytes)victims.push_back({i,old_bytes,old_benefit});
            }
            std::stable_sort(victims.begin(),victims.end(),[](auto& a,auto& b){return __uint128_t(a.benefit)*b.bytes<__uint128_t(b.benefit)*a.bytes;});
            size_t available=budget-route_bytes,needed=0;while(available<bytes&&needed<victims.size())available+=victims[needed++].bytes;
            if(available<bytes){drop(next);return;}victims.resize(needed);
        }
        for(auto v:victims){route_bytes-=v.bytes;for(auto& r:batches[v.index].routes)route_edges-=r.count?r.count:r.pivots.size();drop(batches[v.index]);
        }
        route_bytes+=bytes;for(auto& r:next.routes)route_edges+=r.count?r.count:r.pivots.size();
    }
''')
    s=s.replace('guard(!replay||cursor==batches.size(),"Replay batch count changed");','guard(!replay||(cursor==batches.size()&&decision_cursor==physical_decisions.size()),"Replay batch/physical decision count changed");')
    s=s.replace('sparse-structural-replay-v2','sparse-structural-replay-v3-packed-work')
    s=s.replace("<<max_route_edges<<'\\n';", "<<max_route_edges<<' '<<physical_decisions.size()<<'\\n';for(auto x:physical_decisions)out<<x<<' ';out<<'\\n';")
    s=s.replace("for(const auto& route:b.routes){out<<route.complete<<' '<<route.pivots.size();for(auto x:route.pivots)out<<' '<<x;out<<'\\n';}",'''for(const auto& route:b.routes){out<<route.complete<<' '<<route.pivots.size();for(auto x:route.pivots)out<<' '<<x;
                out<<' '<<route.count<<' '<<route.benefit<<' '<<route.packed.size()<<' ';
                static const char hex[]="0123456789abcdef";if(route.packed.empty())out<<'-';else for(auto x:route.packed)out<<hex[x>>4]<<hex[x&15];out<<'\\n';}''')
    s=s.replace('t.max_route_edges;','t.max_route_edges;size_t nd;in>>nd;guard(bool(in)&&nd<1000000,"Replay decision count");t.physical_decisions.resize(nd);for(auto& x:t.physical_decisions){in>>x;guard(bool(in),"Replay decision data");}',1)
    s=s.replace('t.route_edges+=n;t.route_bytes+=n*sizeof(I);','''t.route_edges+=n;t.route_bytes+=n*sizeof(I);
                size_t bytes;std::string hex;in>>route.count>>route.benefit>>bytes>>hex;
                guard(bool(in)&&bytes<=t.max_route_edges*sizeof(I)-t.route_bytes&&(valid||(bytes==0&&route.count==0))&&(!bytes||!n),"Replay packed budget/mode");
                guard(bytes?hex.size()==2*bytes:hex=="-","Replay packed shape");
                auto digit=[](char x)->unsigned{guard((x>='0'&&x<='9')||(x>='a'&&x<='f'),"Replay hex digit");return x<='9'?x-'0':x-'a'+10;};
                route.packed.resize(bytes);for(size_t j=0;j<bytes;++j)route.packed[j]=(digit(hex[2*j])<<4)|digit(hex[2*j+1]);
                uint32_t delta=0;unsigned shift=0;uint64_t index=0;size_t decoded=0;
                for(auto b:route.packed){guard(shift<35&&(shift!=28||(b&127)<8),"Replay delta overflow");delta|=uint32_t(b&127)<<shift;
                    if(b&128)shift+=7;else{guard(delta>0,"Replay nonpositive delta");index+=delta;guard(index<=t.columns,"Replay packed index");++decoded;delta=0;shift=0;}}
                guard(!shift&&decoded==route.count,"Replay packed truncation/count");t.route_edges+=route.count;t.route_bytes+=bytes;''')
    q.write_text(s)
    q=O/'src/shared/streaming_echelon.hpp';s=q.read_text()
    s=s.replace('#include "structural_replay.hpp"','#include "structural_replay.hpp"\n#include "wide_modular_accumulator.hpp"')
    s=s.replace('size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;','size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;uint64_t coefficient_updates=0,queue_pushes=0;',1)
    s=s.replace('explicit Work(size_t n):dense(n),flags(n){};', 'workload::BitQueue bits;\n        explicit Work(size_t n):dense(n),flags(n),bits(n){};')
    s=s.replace('workload::BitQueue bits;', 'std::vector<__uint128_t> wide;uint64_t wide_updates=0,wide_overflows=0;\n        workload::BitQueue bits;')
    point=s.index('    void reduce(ModVec& row,Work& w,')
    s=s[:point]+'''    void reduce_wide(ModVec& row,Work& w,structural_replay::Route* learn,const structural_replay::Route* route,size_t edge_limit)const{
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
''' +s[point:]
    s=s.replace('bool routed=route&&route->complete;', '''bool routed=route&&route->complete;
        // Delay modular reduction only where overflow is infrequent. Route
        // replay still scans ALL residual pivots and falls back exactly.
        if(workload::kernels()&&fp.mod.n<=(uint64_t(1)<<61)){reduce_wide(row,w,learn,route,edge_limit);return;}''')
    s=s.replace('if(learn){learn->complete=true;learn->pivots.clear();}', 'bool learning=bool(learn);if(learn){learn->complete=true;learn->pivots.clear();}')
    s=s.replace('w.flags[c]|=2;w.heap.push_back(pivot_row[c]);std::push_heap(w.heap.begin(),w.heap.end(),std::greater<I>());', 'w.flags[c]|=2;++w.queue_pushes;if(workload::kernels())w.bits.push(pivot_row[c]);else{w.heap.push_back(pivot_row[c]);std::push_heap(w.heap.begin(),w.heap.end(),std::greater<I>());}')
    s=s.replace('for(I r:route->pivots){if(size_t(r)>=basis.nrow)continue;', 'route->each([&](I r){if(size_t(r)>=basis.nrow)return;')
    s=s.replace('ulong x=w.dense[pivot_columns[r]];if(!x)continue;', 'ulong x=w.dense[pivot_columns[r]];if(!x)return;',1)
    s=s.replace('if(value)touch(c);w.dense[c]=value;}}', 'if(value)touch(c);w.dense[c]=value;}w.coefficient_updates+=basis[r].nnz();});')
    s=s.replace('if(!w.heap.empty())++w.route_fallbacks;', 'if(workload::kernels()?!w.bits.empty():!w.heap.empty())++w.route_fallbacks;')
    s=s.replace('while(!w.heap.empty()){\n            std::pop_heap(w.heap.begin(),w.heap.end(),std::greater<I>());I r=w.heap.back();w.heap.pop_back();++w.heap_steps;\n            if(learn&&learn->complete){if(learn->pivots.size()<edge_limit)learn->pivots.push_back(r);else{learn->complete=false;learn->pivots.clear();}}', '''while(workload::kernels()?!w.bits.empty():!w.heap.empty()){
            I r;if(workload::kernels())r=I(w.bits.pop());else{std::pop_heap(w.heap.begin(),w.heap.end(),std::greater<I>());r=w.heap.back();w.heap.pop_back();}++w.heap_steps;
            if(learning)learning=learn->append(r,edge_limit,basis[r].nnz());''')
    s=s.replace('if(value)enqueue(c);w.dense[c]=value;}','if(value)enqueue(c);w.dense[c]=value;}w.coefficient_updates+=basis[r].nnz();')
    s=s.replace('if(learn)learn->pivots.shrink_to_fit();','if(learn){learn->pivots.shrink_to_fit();learn->packed.shrink_to_fit();}')
    start=s.index('        size_t edge_limit=0;');end=s.index('\n        reduce_cpu+=',start)
    s=s[:start]+'''        size_t edge_limit=0;
        if(trace&&!trace->replay&&basis.nrow){
            size_t available=workload::routes()?(8ULL<<20):trace->max_route_edges*sizeof(I)-trace->route_bytes;
            size_t headers=pending.nrow*sizeof(structural_replay::Route);
            if(available>headers)edge_limit=std::min(size_t(16384),(available-headers)/(workload::routes()?1:sizeof(I))/pending.nrow);
            if(edge_limit)learned.routes.resize(pending.nrow);
        }
        if(basis.nrow){auto action=[&](size_t i){
            reduce(pending[i],work[thread_id()],learned.routes.empty()?nullptr:&learned.routes[i],
                replay&&i<replay->routes.size()?&replay->routes[i]:nullptr,edge_limit);
        };
        if(workload::kernels())workload::balanced(pending.nrow,&opt->pool,[&](size_t i){
            uint64_t n=pending[i].nnz();for(auto[c,x]:pending[i]){I r=pivot_row[c];if(r>=0)n+=basis[r].nnz();}return n;
        },action);else{opt->pool.detach_loop(size_t(0),pending.nrow,action);opt->pool.wait();}}
        if(trace&&!trace->replay){
            if(workload::routes())trace->retain_routes(learned);
            else{trace->route_bytes+=learned.routes.size()*sizeof(structural_replay::Route);for(const auto& r:learned.routes){trace->route_edges+=r.pivots.size();trace->route_bytes+=r.bytes();}}
        }
''' +s[end:]
    s=s.replace('    size_t rank()const', '''    size_t input_nonzeros()const{return input_nnz;}
    uint64_t arithmetic_work()const{uint64_t n=0;for(auto& w:work)n+=w.coefficient_updates;return n;}
    workload::GainWork competing_work;
    size_t retained_nonzeros()const{return retained_nnz;}
    size_t rank()const''')
    s=s.replace('        work.clear();','''        uint64_t updates=0,pushes=0;for(auto& w:work){updates+=w.coefficient_updates;pushes+=w.queue_pushes;}
        std::cout<<"KERNEL_WORK coefficient_updates="<<updates<<" queue_pushes="<<pushes<<" bitqueue="<<workload::kernels()<<std::endl;
        uint64_t wide_updates=0,overflows=0;for(auto& w:work){wide_updates+=w.wide_updates;overflows+=w.wide_overflows;}
        std::cout<<"WIDE_KERNEL updates="<<wide_updates<<" exact_overflow_corrections="<<overflows<<std::endl;
        work.clear();''')
    s=s.replace('work[thread_id()]','work[workload::slot()]')
    q.write_text(s)
    # General operator compression and workload-aware projections.
    exec((S/'projection_edits.py').read_text())
    exec((S/'schedule_edits.py').read_text())
    s=(parent/'tests/joint_physics_v26.cpp').read_text().replace('v26-64mib-dependency-replay-affine-factors-reclaimed','v27-workload-compression-packed-replay')
    s=s.replace('    run.bind(root,descriptor);','''    for(const char* name:{"SOLVER_WORK_KERNELS","SOLVER_WORK_SCHEDULE","SOLVER_EXACT_COMPRESSION","SOLVER_PACKED_ROUTES"})descriptor[name]=workload::option(name)?"1":"0";
    descriptor["early_physics"]=workload::schedule()?"4 light rounds; deterministic nonzero/work budgets; measured gain per work compared with seam elimination including fill; persistent cycles; recorded task/stop replay":"4 light rounds; original persistent COL/MRK cursor and 120-second deadlines";
    descriptor["operator_compression"]=workload::compression()?"generic exact augmented feature-space row-span certificates; collinear and polynomial convolution training adapters":"V26 collinear suffix compression only";
    descriptor["kernel_policy"]=workload::kernels()?"exact delayed 128-bit accumulation for fields <=61 bits; complete residual route fallback; hierarchical bit queue; weighted tasks; exact vector density; bounded scratch":"legacy heap and task kernels";
    descriptor["route_policy"]=workload::routes()?"delta varints; 64 MiB retained plus 8 MiB learning workspace; block benefit-per-byte admission":"legacy 64 MiB prefix admission";
    run.bind(root,descriptor);''',1)
    s=s.replace('schedule(run.c,run.s.seed,policy);','schedule(run.c,run.s.seed,policy,trace);')
    (O/'tests/joint_physics_v27.cpp').write_text(s)
    s=(parent/'tests/v24_generic_control.cpp').read_text().replace('route.pivots.clear();','{route.pivots.clear();route.packed.clear();route.count=0;}')
    (O/'tests/v27_generic_control.cpp').write_text(s)
    shutil.copy2(S/'additional_controls.cpp',O/'tests/v27_additional_controls.cpp')
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++';commands=[]
for target in ['v27_generic_control','v27_additional_controls','joint_physics_v27']:
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF','-I'+str(runtime/'include'),'tests/'+target+'.cpp','-o','bin/'+target,'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread'];commands.append(cmd);subprocess.run(cmd,cwd=O,check=True)
hashes={str(q.relative_to(O)):hashlib.sha256(q.read_bytes()).hexdigest() for q in O.rglob('*') if q.is_file() and any(x in q.parts for x in ['src','tests','bin'])}
(O/'v27-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=commands,hashes=hashes),indent=2)+'\n');print('BUILT V27',flush=True)
