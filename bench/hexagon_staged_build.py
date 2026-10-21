#!/usr/bin/env python3
"""Build a frozen hexagon-v2 experiment with the shared staged sewing solver.

The original package and its outputs are never changed. This adapter preserves
its SparseRREF version, CG convention, physical seed and six-block file format.
"""
import argparse,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--package',required=True,type=Path)
p.add_argument('--output',required=True,type=Path)
p.add_argument('--cxx',default='g++')
p.add_argument('--runtime-prefix',type=Path)
a=p.parse_args();package=a.package.resolve();out=a.output.resolve()
if out.exists():p.error('Use a fresh output directory; existing data is preserved.')
shutil.copytree(package,out,ignore=shutil.ignore_patterns('output','log','bin','.git'))
subprocess.run(['patch','-p1','--input',str(ROOT/'patches/hexagon-v2-staged-sewing.patch')],cwd=out,check=True)
shared=out/'src/shared';shared.mkdir()
for name in ['symrep.hpp','structured_kernel.hpp','certified_kernel.hpp','streaming_echelon.hpp','staged_kernel.hpp']:
    shutil.copy2(ROOT/name,shared/name)
    wrapped=shared/name
    text=wrapped.read_text()
    if name!='symrep.hpp':wrapped.write_text(text.replace('using namespace symrep;','using namespace symrep;\nusing symrep::require;'))
q=shared/'symrep.hpp';q.write_text(q.read_text().replace('#include "bootstrap.hpp"','#include "../bootstrap.hpp"'))
for name in ['irrep_seam_rows.hpp']:
    shutil.copy2(ROOT/'bench'/('hexagon_'+name),shared/name)
shutil.copy2(ROOT/'bench/hexagon_seam_certificate.cpp',out/'tests/seam_certificate.cpp')
shutil.copy2(ROOT/'bench/hexagon_line_check.cpp',out/'tests/line_symbol_check.cpp')
shutil.copy2(ROOT/'bench/hexagon_seam_prime.cpp',out/'tests/seam_prime_benchmark.cpp')
shutil.copy2(ROOT/'bench/hexagon_physical_certificate.cpp',out/'tests/physical_seam_certificate.cpp')
(out/'bin').mkdir()
flags=['-O3','-std=c++20','-Isrc','-Isrc/reduced_d3','-Isrc/SparseRREF']
libs=['-lflint','-lmpfr','-lgmp','-lmimalloc','-ltbb','-pthread']
if a.runtime_prefix:
    prefix=a.runtime_prefix.resolve();flags+=['-I'+str(prefix/'include')];libs=['-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),*libs]
for source,target in [('src/reduced_d3/irr_sew.cpp','irr_sew_staged'),('tests/seam_certificate.cpp','seam_certificate'),('tests/line_symbol_check.cpp','line_symbol_check'),('tests/seam_prime_benchmark.cpp','seam_prime_benchmark'),('tests/physical_seam_certificate.cpp','physical_seam_certificate')]:
    subprocess.run([a.cxx,*flags,source,'-o','bin/'+target,*libs],cwd=out,check=True)
print(out/'bin/irr_sew_staged')
