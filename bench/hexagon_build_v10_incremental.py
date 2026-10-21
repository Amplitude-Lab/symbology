#!/usr/bin/env python3
"""Freeze and compile V10 against an already verified V8 helper package."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--package', required=True, type=Path)
p.add_argument('--cxx', required=True)
p.add_argument('--runtime-prefix', required=True, type=Path)
a = p.parse_args()
package, prefix = a.package.resolve(), a.runtime_prefix.resolve()
assert (package / 'src/shared/mapped_buffers.hpp').is_file()
for source, target in [
    ('hexagon_direct_physical_rows.hpp', 'direct_physical_rows.hpp'),
    ('hexagon_affine_seam_discovery_v10.hpp', 'affine_seam_discovery.hpp'),
]:
    text = (ROOT / 'bench' / source).read_text().replace(
        '#include "reduced_d3/irr_projection.hpp"',
        '#include "../reduced_d3/irr_projection.hpp"')
    (package / 'src/shared' / target).write_text(text)
for source, target in [
    ('hexagon_joint_physics_v10.cpp', 'joint_physics_v10'),
    ('hexagon_direct_physical_probe.cpp', 'direct_physical_probe_v10'),
    ('hexagon_affine_seam_probe.cpp', 'affine_seam_probe_v10'),
]:
    text = (ROOT / 'bench' / source).read_text().replace(
        '"hexagon_direct_physical_rows.hpp"', '"shared/direct_physical_rows.hpp"').replace(
        '"hexagon_affine_seam_discovery_v10.hpp"', '"shared/affine_seam_discovery.hpp"')
    file = package / 'tests' / (target + '.cpp')
    file.write_text(text)
    command = [a.cxx, '-O3', '-std=c++20', '-Isrc', '-Isrc/reduced_d3',
               '-Isrc/SparseRREF', '-I' + str(prefix / 'include'), str(file),
               '-o', str(package / 'bin' / target), '-L' + str(prefix / 'lib'),
               '-Wl,-rpath,' + str(prefix / 'lib'), '-lflint', '-lmpfr',
               '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(command, flush=True)
    subprocess.run(command, cwd=package, check=True)
