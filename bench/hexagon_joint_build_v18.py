#!/usr/bin/env python3
"""Build an isolated V18 row-replay variant and a matched V17 control."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
out = a.output.resolve()
assert not out.exists(), 'Preserve existing packages; choose a fresh output'
source = ROOT/'output_hexagon_17p3_20261007/joint-v17-schedule-package'
shutil.copytree(source, out, ignore=shutil.ignore_patterns('bin','.git','output','log'))
(out/'bin').mkdir()
physics = out/'src/reduced_d3/irr_physics_pipeline.hpp'
code = physics.read_text()
old = 'field_t F;size_t n;std::vector<std::optional<std::vector<ulong>>> pivots;thread_pool* pool;'
assert code.count(old)==1
code = code.replace(old, old+'std::vector<size_t>* witnesses;')
old = 'affine_rows(size_t count,const field_t& field,thread_pool* p=nullptr):F(field),n(count),pivots(n),pool(p){}'
new = 'affine_rows(size_t count,const field_t& field,thread_pool* p=nullptr,std::vector<size_t>* w=nullptr):F(field),n(count),pivots(n),pool(p),witnesses(w){}'
assert code.count(old)==1
code = code.replace(old,new)
old = 'pivots[j]=std::move(r);'
assert code.count(old)==1
code = code.replace(old,old+'if(witnesses)witnesses->push_back(eq.sample);')
physics.write_text(code)
rows = out/'src/shared/irrep_seam_rows.hpp'
code = rows.read_text()
marker = '    size_t residual_blocks()const{return (rows+127)/128;}'
assert code.count(marker)==1
extension = r'''
    // Exact contraction in the existing factor order, avoiding construction,
    // sorting and merging of a full source equation for every checked row.
    MV residual_direct(size_t id,const MM& kt,size_t count)const{
        symrep::require(kt.nrow==columns&&kt.ncol==count&&id<rows,"V18 direct residual layout");
        const auto& p=part(id);size_t r=(id-p.offset)/p.copies,copy=(id-p.offset)%p.copies;
        std::vector<ulong> values(count);
        for(auto[c,x]:p.local[r]){const auto& v=p.columns[c];
            for(auto[j,y]:first[v.block][open_indices[v.block][v.channel][copy]]){
                size_t col=colbase[v.block]+(contract_left?v.copy*last_copies[v.block]+size_t(j):size_t(j)*last_copies[v.block]+v.copy);
                ulong xy=nmod_mul(x,y,F.mod);
                for(auto[s,z]:kt[col])values[s]=nmod_add(values[s],nmod_mul(xy,z,F.mod),F.mod);
            }
        }
        MV result;for(size_t s=0;s<count;++s)if(values[s])result.push_back(int32_t(s),values[s]);return result;
    }
    template<class Emit> void selected_residual_transpose(const MM& kt,size_t count,
            const std::vector<size_t>& ids,Emit&& emit,thread_pool* pool)const{
        std::vector<MV> tile(ids.size());auto fill=[&](size_t i){tile[i]=residual_direct(ids[i],kt,count);};
        if(pool){pool->detach_loop(size_t(0),tile.size(),fill);pool->wait();}
        else for(size_t i=0;i<tile.size();++i)fill(i);
        for(size_t i=0;i<ids.size();++i)emit(ids[i],std::move(tile[i]));
    }
    template<class Emit> void residual_transpose_direct(const MM& kt,size_t count,
            size_t block,Emit&& emit,thread_pool* pool)const{
        size_t begin=block*128,end=std::min(rows,begin+128);symrep::require(begin<rows,"V18 residual block address");
        std::vector<size_t> ids(end-begin);std::iota(ids.begin(),ids.end(),begin);
        selected_residual_transpose(kt,count,ids,[&](size_t,auto&& row){emit(std::move(row));},pool);
    }
'''
rows.write_text(code.replace(marker,extension+'\n'+marker))
shutil.copy2(ROOT/'bench/hexagon_final_seam_v18.hpp',out/'src/shared/hexagon_final_seam_v18.hpp')
driver = (out/'tests/joint_physics_v17.cpp').read_text()
old = 'joint-discovery-physical-affine-v17-resumable-physics-cycle'
assert driver.count(old)==1
driver = driver.replace(old,'joint-discovery-physical-affine-v18-final-row-replay')
marker = 'int main(int argc,char**argv)try {'
assert driver.count(marker)==1
driver = driver.replace(marker,'#include "shared/hexagon_final_seam_v18.hpp"\n\n'+marker)
old = 'return finish_seam(rows,std::move(parent),run.c.F,opt);'
assert driver.count(old)==1
driver = driver.replace(old,'return hexagon_final_seam_v18::finish(rows,std::move(parent),run.c.F,opt,s.output/"final_seam_plans",hex_sha256_file(root/"input_identity.txt"));')
(out/'tests/joint_physics_v18.cpp').write_text(driver)
for name in ('hexagon_v18_seam_control.cpp','hexagon_v18_synthetic.cpp'):
    shutil.copy2(ROOT/'bench'/name,out/'tests'/name)
runtime = ROOT/'output_hexagon_runtime_recovery_20261008/env'
compiler = ROOT/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
commands=[]
for name in ('hexagon_v18_synthetic','hexagon_v18_seam_control','joint_physics_v17','joint_physics_v18'):
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
         '-I'+str(runtime/'include'),'tests/'+name+'.cpp','-o','bin/'+name,
         '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),
         '-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
    commands.append(cmd)
    (out/'build-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    subprocess.run(cmd,cwd=out,check=True)
    print(name,'built',flush=True)
hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in out.rglob('*') if p.is_file() and (p.parts[-2]=='bin' or 'src' in p.parts or 'tests' in p.parts)}
(out/'v18-build-receipt.json').write_text(json.dumps({'parent':str(source),'changes':'V18 final-seam source-row replay, forward continuation and direct factored residuals; earlier V17 discovery/physics unchanged','commands':commands,'hashes':hashes},indent=2)+'\n')
