#!/usr/bin/env python3
"""Build isolated exact-tail COL training and instrumented scalar controls."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

R = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
out = a.output.resolve()
assert not out.exists(), 'Choose a fresh package; preserve completed work'
parent = R/'output_hexagon_v18_20261009/package'
shutil.copytree(parent, out, ignore=shutil.ignore_patterns('bin', '.git', 'output', 'log'))
(out/'bin').mkdir()
# Process CPU is the sum over all its threads, including system CPU.
profile = out/'src/reduced_d3/irr_profile.hpp'
code = profile.read_text().replace('#include <iostream>', '#include <iostream>\n#include <sys/resource.h>')
code = code.replace('namespace irr_run {', '''namespace irr_run {
inline double cpu_seconds(){rusage r{};if(getrusage(RUSAGE_SELF,&r))throw std::runtime_error("getrusage failed");
    return r.ru_utime.tv_sec+r.ru_utime.tv_usec*1e-6+r.ru_stime.tv_sec+r.ru_stime.tv_usec*1e-6;}
''')
code = code.replace('std::string name;', 'std::string name; double cpu_start=cpu_seconds();')
code = code.replace('<<" rss_bytes="', '<<" cpu_seconds="<<cpu_seconds()-cpu_start<<" rss_bytes="')
profile.write_text(code)
commands=[]
runtime=R/'output_hexagon_runtime_recovery_20261008/env'
compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
for version in (17,18,19):
    source = 18 if version==19 else version
    driver=(parent/f'tests/joint_physics_v{source}.cpp').read_text()
    if version==19:
        old='{\"sampling\",\"adaptive\"},{\"projection\",\"scalar\"},{\"samples\",\"0\"}'
        assert driver.count(old)==2
        # Only training uses the complete boundary suffix. Independent final
        # validation continues to use the old scalar slot projection.
        driver=driver.replace(old, '{\"sampling\",\"adaptive\"},{\"projection\",name==\"collinear\"?\"tail-exact\":\"scalar\"},{\"samples\",\"0\"}',1)
        driver=driver.replace('joint-discovery-physical-affine-v18-final-row-replay',
                              'joint-discovery-physical-affine-v19-exact-collinear-tail')
        driver=driver.replace('{\"ope_storage\",', '{\"collinear_training\",\"exact-all-tail-slots; scalar independent validation\"},{\"ope_storage\",')
    old='auto value=compute();hexagon_joint_checkpoint::save'
    assert driver.count(old)==1
    driver=driver.replace(old,'irr_run::phase phase_timer("compute_"+phase);auto value=compute();hexagon_joint_checkpoint::save')
    old='physical_validate(run,current,name);'
    assert driver.count(old)==1
    driver=driver.replace(old,'{irr_run::phase phase_timer("check_"+std::string(name));physical_validate(run,current,name);}')
    (out/f'tests/joint_physics_v{version}.cpp').write_text(driver)
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
         '-I'+str(runtime/'include'),f'tests/joint_physics_v{version}.cpp','-o',f'bin/joint_physics_v{version}',
         '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),
         '-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
    commands.append(cmd)
    subprocess.run(cmd,cwd=out,check=True)
    print('BUILT',version,flush=True)
# A converter comparison compiled against the same dependencies.
cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
     '-I'+str(runtime/'include'),'tests/pipeline_compare.cpp','-o','bin/pipeline_compare',
     '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),
     '-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
commands.append(cmd);subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest()
        for q in out.rglob('*') if q.is_file() and ('src' in q.parts or 'tests' in q.parts or 'bin' in q.parts)}
(out/'v19-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=commands,hashes=hashes,
    change='Exact complete collinear suffix for post-discovery training. Original scalar independent validation, all RHS and other phases preserved. CPU instrumentation on all matched controls.'),indent=2)+'\n')
