#!/usr/bin/env python3
"""Build a separate physics-schedule correction; never modify an active package."""
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
shutil.copy2(ROOT / 'bench/hexagon_physics_work_cursor.hpp', out / 'src/shared/hexagon_physics_work_cursor.hpp')
shutil.copy2(ROOT / 'bench/hexagon_adaptive_physical_rows_v17.hpp', out / 'src/shared/adaptive_physical_rows.hpp')
driver = (out / 'tests/joint_physics_v16.cpp').read_text()
old_algorithm = 'joint-discovery-physical-affine-v16-shared-mrk-shuffles'
new_algorithm = 'joint-discovery-physical-affine-v17-resumable-physics-cycle'
assert driver.count(old_algorithm) == 1
driver = driver.replace(old_algorithm, new_algorithm)
old_schedule = '4 light rounds; at most one COL/MRK cycle per callback; OPE after two cycles; callbacks between refinements and after 120 seconds within a sweep'
new_schedule = '4 light rounds; persistent COL/MRK cycle cursor across 120-second callbacks; due OPE precedes the next cycle'
assert driver.count(old_schedule) == 1
driver = driver.replace(old_schedule, new_schedule)
(out / 'tests/joint_physics_v17.cpp').write_text(driver)
prefix = Path('/tmp/symbology-audit-toolchain/env')
cmd = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
       '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
       '-I' + str(prefix / 'include'), 'tests/joint_physics_v17.cpp', '-o', 'bin/joint_physics_v17',
       '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
       '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
print(cmd, flush=True)
subprocess.run(cmd, cwd=out, check=True)
