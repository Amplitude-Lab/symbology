#!/usr/bin/env python3
"""Build bounded physics interleaving, including callbacks within long seam sweeps."""
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
for source, target in [('hexagon_adaptive_physical_rows_v13.hpp', 'adaptive_physical_rows.hpp'),
                       ('hexagon_affine_seam_discovery_v13.hpp', 'affine_seam_discovery.hpp')]:
    shutil.copy2(ROOT / 'bench' / source, out / 'src/shared' / target)
driver = (out / 'tests/joint_physics_v12.cpp').read_text()
driver = driver.replace('joint-discovery-physical-affine-v12-optional-double-bonus',
    'joint-discovery-physical-affine-v13-bounded-interleaving')
driver = driver.replace('4 light rounds; bounded adaptive collinear/MRK and leading OPE at min(50000,n/5) nullity',
    '4 light rounds; at most one COL/MRK cycle per callback; OPE after two cycles; callbacks between refinements and after 120 seconds within a sweep')
(out / 'tests/joint_physics_v13.cpp').write_text(driver)
shutil.copy2(ROOT / 'bench/hexagon_affine_seam_probe_v13.cpp', out / 'tests/affine_probe_v13.cpp')
prefix = Path('/tmp/symbology-audit-toolchain/env')
for name in ['affine_probe_v13', 'joint_physics_v13']:
    cmd = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
           '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
           '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
           '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
           '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(cmd, flush=True)
    subprocess.run(cmd, cwd=out, check=True)
