#!/usr/bin/env python3
"""Release dead temporary allocations at joined solver phase boundaries."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
parent=R/'output_hexagon_v23_20261009/package-v24';assert not out.exists();shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
q=out/'src/shared/temporary_allocator_release.hpp';q.write_text('''#pragma once
#if defined(__GLIBC__)
#include <malloc.h>
#endif
// Call only after temporary owners have been destroyed and worker tasks joined.
// Live input matrices, recurrence caches and affine coefficients remain owned.
namespace temporary_allocator {
inline bool release_unused_pages(){
#if defined(__GLIBC__)
    return malloc_trim(0)!=0;
#else
    return false;
#endif
}
}
''')
s=(parent/'tests/joint_physics_v24.cpp').read_text().replace('#include "shared/hexagon_mhv_bonus_seed.hpp"','#include "shared/hexagon_mhv_bonus_seed.hpp"\n#include "shared/temporary_allocator_release.hpp"')
s=s.replace('v24-bounded-dependency-structural-replay-factors','v25-bounded-dependency-replay-affine-factors-reclaimed')
old='irr_run::phase phase_timer("compute_"+phase);auto value=compute();hexagon_joint_checkpoint::save(stage,value,p,s.h,s.k);return value;'
new='''irr_run::phase phase_timer("compute_"+phase);auto value=compute();
            {
                // compute() has returned: its temporary owners (including any
                // Trace) are dead; value owns the complete affine map.
                auto before=hexagon_native_memory::resident_bytes();auto start=std::chrono::steady_clock::now();
                bool released=temporary_allocator::release_unused_pages();
                std::cout<<"TEMPORARY_ALLOCATOR_RELEASE phase="<<phase<<" released="<<released<<" rss_before="<<before
                    <<" rss_after="<<hexagon_native_memory::resident_bytes()<<" seconds="
                    <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
            }
            hexagon_joint_checkpoint::save(stage,value,p,s.h,s.k);return value;'''
assert old in s;s=s.replace(old,new);(out/'tests/joint_physics_v25.cpp').write_text(s)
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF','-I'+str(runtime/'include'),'tests/joint_physics_v25.cpp','-o','bin/joint_physics_v25','-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread'];subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*') if q.is_file() and any(x in q.parts for x in ['src','tests','bin'])}
(out/'v25-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=[cmd],hashes=hashes),indent=2)+'\n');print('BUILT V25',flush=True)
