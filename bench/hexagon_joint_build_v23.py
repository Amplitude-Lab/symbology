#!/usr/bin/env python3
"""Isolated generic structural replay and triangular affine factor candidate."""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
out=a.output.resolve();parent=R/'output_hexagon_v19_20261009/package-v22'
assert not out.exists();shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
for name in ['structural_replay.hpp','sparse_affine_factor.hpp']:
    shutil.copy2(R/'bench'/name,out/'src/shared'/name)
p=out/'src/shared/streaming_echelon.hpp';s=p.read_text().replace('#include <queue>','#include <queue>\n#include "structural_replay.hpp"')
s=s.replace('    ModMat basis,pending;','    ModMat basis,pending;\n    structural_replay::Trace* trace=nullptr;')
s=s.replace('        pending.nrow=pending.rows.size();\n        if(basis.nrow)', '''        pending.nrow=pending.rows.size();
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
        if(basis.nrow)''')
old='''        std::erase_if(pending.rows,[](const auto& row){return !row.nnz();});pending.nrow=pending.rows.size();
        if(pending.nrow){
            pending.sort_rows_by_nnz();auto pivots=sparse_mat_rref_forward(pending,fp,opt);require(!opt->abort,"streaming modular kernel aborted");'''
new='''        std::vector<size_t> source_order;
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
            auto pivots=replay?structural_replay::forward(pending,*replay,fp,opt):sparse_mat_rref_forward(pending,fp,opt);
            require(!opt->abort,"streaming modular kernel aborted");
            if(trace&&!replay){learned.groups=pivots;
                for(auto& group:learned.groups)for(auto& p:group){learned.selected.push_back(source_order[p.r]);p.r=index(learned.selected.size()-1);}}
'''
assert old in s;s=s.replace(old,new)
old='''        pending.clear();pending=ModMat(0,basis.ncol);pending_nnz=0;++batches;'''
new='''        if(replay){
            opt->pool.detach_loop(size_t(0),discarded.nrow,[&](size_t i){reduce(discarded[i],work[thread_id()]);});opt->pool.wait();
            for(const auto& row:discarded.rows)structural_replay::guard(!row.nnz(),"Replay discarded equation gained rank or changed affine RHS");
        }
        if(trace){if(!trace->replay)trace->batches.push_back(std::move(learned));++trace->cursor;}
        pending.clear();pending=ModMat(0,basis.ncol);pending_nnz=0;++batches;'''
assert old in s;s=s.replace(old,new)
s=s.replace('size_t nonzeros=16*1024*1024)','size_t nonzeros=16*1024*1024,structural_replay::Trace* plan=nullptr)')
s=s.replace(':basis(0,columns),pending(0,columns),pivot_row',':basis(0,columns),pending(0,columns),trace(plan),pivot_row')
s=s.replace('        require(rows>0&&nonzeros>0,', '        if(trace){if(trace->replay)structural_replay::guard(trace->columns==columns,"Replay column count changed");else trace->columns=columns;}\n        require(rows>0&&nonzeros>0,')
s=s.replace('if(pending.rows.size()>=row_limit||pending_nnz>=nnz_limit)flush();', '''if(trace&&trace->replay){
            structural_replay::guard(trace->cursor<trace->batches.size(),"Replay extra source rows");
            if(pending.rows.size()>=trace->batches[trace->cursor].inputs)flush();
        }else if(pending.rows.size()>=row_limit||pending_nnz>=nnz_limit)flush();''')
s=s.replace('        flush();std::vector<std::vector<pivot_t<I>>> pivots(2);','        flush();if(trace)trace->complete();std::vector<std::vector<pivot_t<I>>> pivots(2);')
# Measure the expensive elimination components, including replay verification.
s=s.replace('    structural_replay::Trace* trace=nullptr;','''    structural_replay::Trace* trace=nullptr;
    double reduce_cpu=0,forward_cpu=0,verify_cpu=0,back_cpu=0;
    size_t independent_inputs=0,dependent_inputs=0;
    static double cpu_now(){return double(std::clock())/CLOCKS_PER_SEC;}''')
s=s.replace('        if(basis.nrow){opt->pool.detach_loop', '        double reduce_started=cpu_now();\n        if(basis.nrow){opt->pool.detach_loop')
s=s.replace('        std::vector<size_t> source_order;', '        reduce_cpu+=cpu_now()-reduce_started;\n        std::vector<size_t> source_order;')
s=s.replace('            auto pivots=replay?', '            double forward_started=cpu_now();\n            auto pivots=replay?')
s=s.replace('            require(!opt->abort,"streaming modular kernel aborted");', '            forward_cpu+=cpu_now()-forward_started;\n            require(!opt->abort,"streaming modular kernel aborted");')
s=s.replace('                require(pivot_row[c]<0,"duplicate streaming pivot");', '                if(replay)structural_replay::guard(pivot_row[c]<0,"Replay duplicate pivot");\n                require(pivot_row[c]<0,"duplicate streaming pivot");')
s=s.replace('            for(size_t r=first;r<basis.nrow;++r){', '''            independent_inputs+=basis.nrow-first;
            for(size_t r=first;r<basis.nrow;++r){
                if(replay){
                    auto x=basis[r].find(pivot_columns[r]);structural_replay::guard(x&&*x==1,"Replay pivot normalization changed");
                    for(auto c:basis[r].index_span())structural_replay::guard(pivot_row[c]<0||pivot_row[c]>=index(r),"Replay triangular dependency changed");
                }''')
s=s.replace('        if(replay){\n            opt->pool.detach_loop', '        if(replay){\n            double verify_started=cpu_now();dependent_inputs+=discarded.nrow;\n            opt->pool.detach_loop')
s=s.replace('        if(trace){if(!trace->replay)', '        if(trace){if(!trace->replay)')
s=s.replace('            for(const auto& row:discarded.rows)structural_replay::guard(!row.nnz(),"Replay discarded equation gained rank or changed affine RHS");', '            for(const auto& row:discarded.rows)structural_replay::guard(!row.nnz(),"Replay discarded equation gained rank or changed affine RHS");\n            verify_cpu+=cpu_now()-verify_started;')
s=s.replace('        triangular_solver_2(basis,pivots,fp,opt);', '        double back_started=cpu_now();triangular_solver_2(basis,pivots,fp,opt);back_cpu+=cpu_now()-back_started;')
s=s.replace('        return {std::move(basis),std::move(pivots)};', '''        std::cout<<"STRUCTURAL_PROFILE replay="<<bool(trace&&trace->replay)<<" reduce_cpu="<<reduce_cpu
            <<" forward_cpu="<<forward_cpu<<" discarded_verify_cpu="<<verify_cpu<<" back_cpu="<<back_cpu
            <<" independent_inputs="<<independent_inputs<<" verified_dependent_inputs="<<dependent_inputs<<std::endl;
        return {std::move(basis),std::move(pivots)};''')
p.write_text(s)
p=out/'src/shared/affine_seam_discovery.hpp';s=p.read_text()
s=s.replace('double physical_interval_seconds=120){','double physical_interval_seconds=120,structural_replay::Trace* trace=nullptr){')
s=s.replace('options.batch_rows,options.batch_nonzeros);','options.batch_rows,options.batch_nonzeros,trace);');p.write_text(s)
p=out/'src/reduced_d3/irr_physics_pipeline.hpp';s=p.read_text().replace('#include <memory>','#include <memory>\n#include "shared/sparse_affine_factor.hpp"')
needle='''    auto solve() {
        // Back-substitute'''
insert='''    auto solve_factored() {
        sparse_affine::Factor factor(pivots,F);auto map=factor.materialize();
        nmhv_physical_affine_solution_t<int32_t> result;result.origin.resize(n);
        for(auto[j,x]:map[map.nrow-1])result.origin[j]=x;
        for(size_t j=0;j<n;++j)if(pivots[j])result.pivot_columns.push_back(j);
        result.rank=result.pivot_columns.size();map.rows.pop_back();--map.nrow;result.directions=std::move(map);
        return result;
    }
    auto solve() {
        if(auto value=std::getenv("HEXAGON_AFFINE_FACTORS");value&&std::string(value)=="1")return solve_factored();
        // Back-substitute'''
assert needle in s;p.write_text(s.replace(needle,insert))
p=out/'tests/joint_physics_v23.cpp';s=(parent/'tests/joint_physics_v22.cpp').read_text().replace('joint-discovery-physical-affine-v22-weight-dependent-ope-head','joint-discovery-physical-affine-v23b-guarded-profiled-structural-replay-factors')
s=s.replace('    run.bind(root,descriptor);','''    auto switch_value=[](const char* name){auto value=std::getenv(name);std::string x=value?value:"0";require(x=="0"||x=="1","Invalid optimization switch");return x;};
    descriptor["discovery_replay"]=switch_value("HEXAGON_DISCOVERY_REPLAY");
    descriptor["affine_factors"]=switch_value("HEXAGON_AFFINE_FACTORS");
    run.bind(root,descriptor);''',1)
start=s.index('            hexagon_adaptive_physics::DirectSchedule schedule(');end=s.index('\n        });',start)
s=s[:start]+'''            auto compute=[&](structural_replay::Trace* trace){
                hexagon_adaptive_physics::DirectSchedule schedule(run.c,run.s.seed,policy);
                auto early=[&](streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate,size_t round){
                    if(bonus&&round==0)mhv_double_bonus::seed_original_coordinates(run.c,reducer,coordinate);
                    schedule(reducer,coordinate,round);
                };
                return staged_kernel::discover_affine_candidate(stack,run.c.F,opt,options,early,variable_order,120,trace);
            };
            if(descriptor.at("discovery_replay")=="0")return compute(nullptr);
            auto directory=s.output/"discovery_plans";fs::create_directories(directory);irr_stage::lock guard(directory/".active.lock");
            auto path=directory/"structure.tsv";auto identity=hex_sha256_file(root/"input_identity.txt");
            bool existing=fs::exists(path);structural_replay::Trace trace;
            try{
                if(existing){
                    structural_replay::guard(irr_stage::text(path.string()+".sha256")==hex_sha256_file(path)+"\\n","Replay checksum mismatch");
                    std::ifstream in(path);trace=structural_replay::Trace::read(in,identity);
                }
                auto value=compute(&trace);
                if(!existing){auto temp=fs::path(path.string()+".tmp");require(!fs::exists(temp),"Incomplete discovery plan retained");
                    {std::ofstream out(temp);trace.write(out,identity);out.flush();require(bool(out),"Discovery plan write");}
                    hexagon_final_seam_v18::sync_file(temp);auto digest=hex_sha256_file(temp);
                    fs::rename(temp,path);irr_stage::fixed_text(path.string()+".sha256",digest+"\\n");
                    hexagon_final_seam_v18::sync_file(path.string()+".sha256");hexagon_final_seam_v18::sync_file(directory);
                }
                std::cout<<"DISCOVERY_STRUCTURE replay="<<existing<<" batches="<<trace.batches.size()<<" discarded_equations_verified=1"<<std::endl;
                return value;
            }catch(const structural_replay::Mismatch& e){
                std::cout<<"DISCOVERY_STRUCTURE_FALLBACK reason="<<e.what()<<" original_plan_retained=1"<<std::endl;
                return compute(nullptr);
            }
''' + s[end:]
p.write_text(s)
if a.prepare_only:print(out);raise SystemExit()
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF','-I'+str(runtime/'include'),'tests/joint_physics_v23.cpp','-o','bin/joint_physics_v23','-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*') if q.is_file() and any(x in q.parts for x in ['src','tests','bin'])}
(out/'v23-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=[cmd],hashes=hashes),indent=2)+'\n')
print('BUILT V23',flush=True)
