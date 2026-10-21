#!/usr/bin/env python3
"""Bounded dependency-route replay on the generic V23 discovery reducer."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
parent=R/'output_hexagon_v23_20261009/package-v23b';assert not out.exists();shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
p=out/'src/shared/structural_replay.hpp';s=p.read_text()
s=s.replace('struct Batch {','struct Route {bool complete=false;std::vector<I> pivots;};\nstruct Batch {std::vector<Route> routes;')
s=s.replace('    size_t columns=0,cursor=0;bool replay=false;','    size_t columns=0,cursor=0;bool replay=false;\n    size_t route_edges=0,route_bytes=0,max_route_edges=(256ULL<<20)/sizeof(I);')
s=s.replace('sparse-structural-replay-v1','sparse-structural-replay-v2')
s=s.replace("<<columns<<' '<<batches.size()<<'\\n';", "<<columns<<' '<<batches.size()<<' '<<max_route_edges<<'\\n';")
s=s.replace("for(const auto& g:b.groups){out<<g.size();for(auto[r,c]:g)out<<' '<<r<<' '<<c;out<<'\\n';}}", """for(const auto& g:b.groups){out<<g.size();for(auto[r,c]:g)out<<' '<<r<<' '<<c;out<<'\\n';}
            out<<b.routes.size()<<'\\n';for(const auto& route:b.routes){out<<route.complete<<' '<<route.pivots.size();for(auto x:route.pivots)out<<' '<<x;out<<'\\n';}}
""")
s=s.replace('in>>t.columns>>count;', 'in>>t.columns>>count>>t.max_route_edges;')
s=s.replace('count<1000000,','count<1000000&&t.max_route_edges<=(256ULL<<20)/sizeof(I),')
s=s.replace('guard(total==nr,"Replay missing pivot");}', '''guard(total==nr,"Replay missing pivot");
            size_t routes;in>>routes;guard(bool(in)&&routes<=b.inputs&&routes<=(t.max_route_edges*sizeof(I)-t.route_bytes)/sizeof(Route),"Replay route count/budget");b.routes.resize(routes);t.route_bytes+=routes*sizeof(Route);
            for(auto& route:b.routes){size_t n;unsigned valid;in>>valid>>n;
                guard(bool(in)&&valid<=1&&n<=(t.max_route_edges*sizeof(I)-t.route_bytes)/sizeof(I)&&(valid||n==0),"Replay route length");
                route.complete=valid;route.pivots.resize(n);I previous=-1;
                for(auto& x:route.pivots){in>>x;guard(bool(in)&&x>=0&&size_t(x)<t.columns&&x>previous,"Replay route index/order");previous=x;}
                t.route_edges+=n;t.route_bytes+=n*sizeof(I);
            }}''')
p.write_text(s)
p=out/'src/shared/streaming_echelon.hpp';s=p.read_text()
s=s.replace('explicit Work(size_t n):dense(n),flags(n){};', 'size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;\n        explicit Work(size_t n):dense(n),flags(n){};')
start=s.index('    void reduce(ModVec& row,Work& w)const{');end=s.index('    void flush(){',start)
s=s[:start]+'''    void reduce(ModVec& row,Work& w,structural_replay::Route* learn=nullptr,
                const structural_replay::Route* route=nullptr,size_t edge_limit=0)const{
        bool routed=route&&route->complete;
        if(learn){learn->complete=true;learn->pivots.clear();}
        auto touch=[&](I c){if(!(w.flags[c]&1)){w.flags[c]|=1;w.touched.push_back(c);}};
        auto enqueue=[&](I c){touch(c);if(pivot_row[c]>=0&&!(w.flags[c]&2)){
            w.flags[c]|=2;w.heap.push_back(pivot_row[c]);std::push_heap(w.heap.begin(),w.heap.end(),std::greater<I>());}};
        for(auto[c,x]:row)if(substitutions.image[c]>=0){I j=substitutions.image[c];ulong y=substitutions.scale[c]==1?x:nmod_mul(x,substitutions.scale[c],fp.mod);
            w.dense[j]=nmod_add(w.dense[j],y,fp.mod);if(routed)touch(j);else enqueue(j);}
        if(routed){++w.route_hits;
            for(I r:route->pivots){if(size_t(r)>=basis.nrow)continue;++w.route_steps;
                ulong x=w.dense[pivot_columns[r]];if(!x)continue;ulong shoup=n_mulmod_precomp_shoup(x,fp.mod.n);
                for(auto[c,y]:basis[r]){ulong value=_nmod_sub(w.dense[c],n_mulmod_shoup(x,y,shoup,fp.mod.n),fp.mod);if(value)touch(c);w.dense[c]=value;}}
            // Routing is only a performance hint: coefficients left on ANY
            // pivot are eliminated by the ordinary heap. This is an exact
            // local fallback even when a cancellation changes at a new prime.
            size_t touched=w.touched.size();for(size_t i=0;i<touched;++i){I c=w.touched[i];if(w.dense[c]&&pivot_row[c]>=0)enqueue(c);}
            if(!w.heap.empty())++w.route_fallbacks;
        }
        while(!w.heap.empty()){
            std::pop_heap(w.heap.begin(),w.heap.end(),std::greater<I>());I r=w.heap.back();w.heap.pop_back();++w.heap_steps;
            if(learn&&learn->complete){if(learn->pivots.size()<edge_limit)learn->pivots.push_back(r);else{learn->complete=false;learn->pivots.clear();}}
            ulong x=w.dense[pivot_columns[r]];if(!x)continue;ulong shoup=n_mulmod_precomp_shoup(x,fp.mod.n);
            for(auto[c,y]:basis[r]){ulong value=_nmod_sub(w.dense[c],n_mulmod_shoup(x,y,shoup,fp.mod.n),fp.mod);if(value)enqueue(c);w.dense[c]=value;}
        }
        if(learn)learn->pivots.shrink_to_fit();
        row.zero();for(I c:w.touched){if(w.dense[c]){require(pivot_row[c]<0,"streaming triangular reduction failed");row.push_back(c,w.dense[c]);}w.dense[c]=0;w.flags[c]=0;}
        w.touched.clear();row.sort_indices();row.compress();
    }
''' + s[end:]
old='''        if(basis.nrow){opt->pool.detach_loop(size_t(0),pending.nrow,[&](size_t i){reduce(pending[i],work[thread_id()]);});opt->pool.wait();}'''
new='''        size_t edge_limit=0;
        if(trace&&!trace->replay&&basis.nrow){
            size_t available=trace->max_route_edges*sizeof(I)-trace->route_bytes;
            size_t headers=pending.nrow*sizeof(structural_replay::Route);
            if(available>headers)edge_limit=std::min(size_t(16384),(available-headers)/sizeof(I)/pending.nrow);
            if(edge_limit){learned.routes.resize(pending.nrow);trace->route_bytes+=headers;}
        }
        if(basis.nrow){opt->pool.detach_loop(size_t(0),pending.nrow,[&](size_t i){
            reduce(pending[i],work[thread_id()],learned.routes.empty()?nullptr:&learned.routes[i],
                replay&&i<replay->routes.size()?&replay->routes[i]:nullptr,edge_limit);
        });opt->pool.wait();}
        if(trace&&!trace->replay)for(const auto& route:learned.routes){trace->route_edges+=route.pivots.size();trace->route_bytes+=route.pivots.size()*sizeof(I);}
'''
assert old in s;s=s.replace(old,new)
s=s.replace('        work.clear();','''        size_t heap_steps=0,route_steps=0,route_hits=0,route_fallbacks=0;
        for(const auto& w:work){heap_steps+=w.heap_steps;route_steps+=w.route_steps;route_hits+=w.route_hits;route_fallbacks+=w.route_fallbacks;}
        std::cout<<"DEPENDENCY_PROFILE heap_steps="<<heap_steps<<" route_steps="<<route_steps<<" route_hits="<<route_hits
            <<" local_fallbacks="<<route_fallbacks<<" stored_edges="<<(trace?trace->route_edges:0)
            <<" stored_route_bytes="<<(trace?trace->route_bytes:0)<<" route_payload_cap_bytes="<<(trace?trace->max_route_edges*sizeof(I):0)<<std::endl;
        work.clear();''')
p.write_text(s)
p=out/'tests/joint_physics_v24.cpp';s=(parent/'tests/joint_physics_v23.cpp').read_text().replace('v23b-guarded-profiled-structural-replay-factors','v24-bounded-dependency-structural-replay-factors')
s=s.replace('            try{\n                if(existing){','''            if(!existing&&(fs::exists(path.string()+".tmp")||fs::exists(path.string()+".sha256"))){
                std::cout<<"DISCOVERY_STRUCTURE_FALLBACK reason=incomplete_plan original_plan_retained=1"<<std::endl;return compute(nullptr);
            }
            try{
                if(existing){
                    structural_replay::guard(fs::exists(path.string()+".sha256"),"Replay checksum sidecar missing");''')
p.write_text(s)
control=(R/'bench/hexagon_v23_generic_control.cpp').read_text()
control=control.replace('++replays;}','''++replays;
                    std::stringstream damaged;learned.write(damaged,"synthetic-family");
                    auto incomplete=structural_replay::Trace::read(damaged,"synthetic-family");
                    for(auto& batch:incomplete.batches)for(auto& route:batch.routes)route.pivots.clear();
                    check(reduce(a,n+1,F,opt,&incomplete)==reference,"local route fallback row space differs");
                }''')
(out/'tests/v24_generic_control.cpp').write_text(control)
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++';commands=[]
shutil.copy2(R/'bench/hexagon_v24_affine_benchmark.cpp',out/'tests/v24_affine_benchmark.cpp')
for target in ['v24_generic_control','v24_affine_benchmark','joint_physics_v24']:
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF','-I'+str(runtime/'include'),'tests/'+target+'.cpp','-o','bin/'+target,'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread'];commands.append(cmd);subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*') if q.is_file() and any(x in q.parts for x in ['src','tests','bin'])}
(out/'v24-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=commands,hashes=hashes),indent=2)+'\n');print('BUILT V24',flush=True)
