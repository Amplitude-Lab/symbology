#!/usr/bin/env python3
"""Use one conservative dependency-route budget for every workload."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
parent=R/'output_hexagon_v23_20261009/package-v25';assert not out.exists();shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
q=out/'src/shared/structural_replay.hpp';s=q.read_text();assert s.count('256ULL<<20')==2;q.write_text(s.replace('256ULL<<20','64ULL<<20'))
s=(parent/'tests/joint_physics_v25.cpp').read_text().replace('v25-bounded-dependency-replay-affine-factors-reclaimed','v26-64mib-dependency-replay-affine-factors-reclaimed')
s=s.replace('    run.bind(root,descriptor);','    descriptor["discovery_route_budget_bytes"]="67108864";\n    run.bind(root,descriptor);',1)
(out/'tests/joint_physics_v26.cpp').write_text(s)
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++';commands=[]
for target in ['v24_generic_control','joint_physics_v26']:
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF','-I'+str(runtime/'include'),'tests/'+target+'.cpp','-o','bin/'+target,'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread'];commands.append(cmd);subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*') if q.is_file() and any(x in q.parts for x in ['src','tests','bin'])}
(out/'v26-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=commands,hashes=hashes),indent=2)+'\n');print('BUILT V26',flush=True)
