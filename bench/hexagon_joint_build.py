#!/usr/bin/env python3
"""Build the isolated early-physics hexagon experiment from a v2 package."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--package',required=True,type=Path)
p.add_argument('--output',required=True,type=Path)
p.add_argument('--cxx',default='g++')
p.add_argument('--runtime-prefix',required=True,type=Path)
a=p.parse_args();out=a.output.resolve();prefix=a.runtime_prefix.resolve()
if out.exists():p.error('Use a fresh package directory; existing data is retained.')
shutil.copytree(a.package.resolve(),out,ignore=shutil.ignore_patterns('output','log','bin','.git'))
subprocess.run(['patch','-p1','--input',str(ROOT/'patches/hexagon-v2-staged-sewing.patch')],cwd=out,check=True)
subprocess.run(['patch','-p1','--input',str(ROOT/'patches/hexagon-v2-ope-head-cache.patch')],cwd=out,check=True)
shared=out/'src/shared';shared.mkdir()
for name in ['symrep.hpp','structured_kernel.hpp','certified_kernel.hpp','streaming_echelon.hpp','staged_kernel.hpp']:
    text=(ROOT/name).read_text()
    if name=='symrep.hpp':text=text.replace('#include "bootstrap.hpp"','#include "../bootstrap.hpp"')
    else:text=text.replace('using namespace symrep;','using namespace symrep;\nusing symrep::require;')
    (shared/name).write_text(text)
for name in ['irrep_seam_rows.hpp','seam_discovery.hpp']:
    shutil.copy2(ROOT/'bench'/('hexagon_'+name),shared/name)
shutil.copy2(ROOT/'bench/hexagon_joint_physics.cpp',out/'tests/joint_physics.cpp')
shutil.copy2(ROOT/'bench/hexagon_joint_certificate.cpp',out/'tests/joint_certificate.cpp')
shutil.copy2(ROOT/'bench/hexagon_origin_check.cpp',out/'tests/origin_check.cpp')
(out/'bin').mkdir()
for target in ['joint_physics','joint_certificate','origin_check']:
    subprocess.run([a.cxx,'-O3','-std=c++20','-Isrc','-Isrc/reduced_d3','-Isrc/SparseRREF',
        '-I'+str(prefix/'include'),'tests/'+target+'.cpp','-o','bin/'+target,
        '-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),
        '-lflint','-lmpfr','-lgmp','-lmimalloc','-ltbb','-pthread'],cwd=out,check=True)
print(out/'bin/joint_physics')
