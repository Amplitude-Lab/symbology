#!/usr/bin/env python3
"""Build exact triangular logarithmic OPE storage from V14; preserve square reference."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--package', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
out = a.output.resolve()
if out.exists():
    p.error('Output must be fresh.')
shutil.copytree(a.package.resolve(), out, ignore=shutil.ignore_patterns('bin', '.git', 'output', 'log'))
(out / 'bin').mkdir()
germ = out / 'src/nmhv_ope_germ.hpp'
source = germ.read_text()
(out / 'src/nmhv_ope_germ_square_reference.hpp').write_text(source)
source = source.replace('stride((mt+1)*(gd+1)*logs*logs)', 'stride((mt+1)*(gd+1)*logs*(logs+1)/2)')
source = source.replace('return ((t*(grade+1)+z+t)*logs+a)*logs+b;',
    'return size_t(t*(grade+1)+z+t)*(logs*(logs+1)/2)+a*logs-a*(a-1)/2+b;')
source = source.replace('t&&size_t(a+1)<logs', 't&&size_t(a+b)<w')
source = source.replace('!t&&size_t(b+1)<logs', '!t&&size_t(a+b)<w')
assert 'logs*logs)' not in source
germ.write_text(source)
adjoint = out / 'src/shared/direct_ope_rows.hpp'
source = adjoint.read_text().replace('t&&a+1<germ.logs', 't&&a+b<w').replace('!t&&b+1<germ.logs', '!t&&a+b<w')
adjoint.write_text(source)
driver = (out / 'tests/joint_physics_v14.cpp').read_text().replace(
    'joint-discovery-physical-affine-v14-bounded-interleaving-cache-release',
    'joint-discovery-physical-affine-v15-bounded-triangular-ope')
driver = driver.replace('{"ope_coverage",', '{"ope_storage","triangular_log_powers_a_plus_b_le_weight"},{"ope_coverage",')
(out / 'tests/joint_physics_v15.cpp').write_text(driver)
shutil.copy2(ROOT / 'bench/hexagon_ope_triangular_probe.cpp', out / 'tests/ope_triangular_probe.cpp')
prefix = Path('/tmp/symbology-audit-toolchain/env')
for name in ['ope_triangular_probe', 'joint_physics_v15']:
    cmd = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
           '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
           '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
           '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
           '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(cmd, flush=True)
    subprocess.run(cmd, cwd=out, check=True)
