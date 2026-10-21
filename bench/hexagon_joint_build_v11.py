#!/usr/bin/env python3
"""Build a fresh V11 package; leave the active V10 package unchanged."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--package', required=True, type=Path)
p.add_argument('--output', required=True, type=Path)
p.add_argument('--cxx', required=True)
p.add_argument('--runtime-prefix', required=True, type=Path)
a = p.parse_args()
out, prefix = a.output.resolve(), a.runtime_prefix.resolve()
if out.exists():
    p.error('Use a fresh output package; previous packages are preserved.')
shutil.copytree(a.package.resolve(), out, ignore=shutil.ignore_patterns('bin', '.git', 'output', 'log'))
(out / 'bin').mkdir()
shared = out / 'src/shared'
for source, target in [
    ('hexagon_streaming_echelon_v11.hpp', 'streaming_echelon.hpp'),
    ('hexagon_mapped_buffers_v11.hpp', 'mapped_buffers.hpp'),
    ('hexagon_adaptive_physical_rows_v11.hpp', 'adaptive_physical_rows.hpp'),
    ('hexagon_direct_ope_rows.hpp', 'direct_ope_rows.hpp'),
    ('hexagon_streaming_echelon_v8.hpp', 'streaming_echelon_v8_comparison.hpp'),
]:
    text = (ROOT / 'bench' / source).read_text()
    if source.startswith('hexagon_streaming_echelon'):
        text = text.replace('using namespace symrep;', 'using namespace symrep;\nusing symrep::require;')
    text = text.replace('#include "reduced_d3/irr_projection.hpp"', '#include "../reduced_d3/irr_projection.hpp"')
    (shared / target).write_text(text)

driver = (ROOT / 'bench/hexagon_joint_physics_v10.cpp').read_text()
driver = driver.replace('#include "shared/affine_seam_discovery.hpp"',
    '#include "shared/affine_seam_discovery.hpp"\n#include "shared/adaptive_physical_rows.hpp"')
driver = driver.replace('joint-discovery-physical-affine-v10-direct-physical', 'joint-discovery-physical-affine-v11-adaptive-compact')
driver = driver.replace('4 rounds: 1 exact-tail collinear sample + 1 MRK sample per round',
    '4 light rounds; bounded adaptive collinear/MRK and leading OPE at min(50000,n/5) nullity')
driver = driver.replace('options.max_residual_nullity=target+1;', 'options.max_residual_nullity=std::max(target+1,size_t(30001));')
start = driver.index('            auto early=[&](')
end = driver.index('            return staged_kernel::discover_affine_candidate(', start)
driver = driver[:start] + '''            hexagon_adaptive_physics::Policy policy;
            policy.target_augmented_nullity=target+1;
            policy.activate_nullity=std::min(size_t(50000),rows.ncol()/5);
            hexagon_adaptive_physics::DirectSchedule early(run.c,run.s.seed,policy);
''' + driver[end:]
driver = driver.replace('        auto& rows=*row_storage;', '''#ifdef USE_MIMALLOC
        mi_collect(true); // Only unused native allocator pages may be returned.
#endif
        auto& rows=*row_storage;''')
(out / 'tests/joint_physics_v11.cpp').write_text(driver)

targets = [('joint_physics_v11', None),
           ('direct_ope_probe_v11', 'hexagon_direct_ope_probe.cpp'),
           ('capacity_probe_v11', 'hexagon_capacity_v11_probe.cpp'),
           ('mapped_probe_v11', 'hexagon_mapped_probe.cpp'),
           ('affine_probe_v11', 'hexagon_affine_seam_probe.cpp')]
for name, source in targets:
    if source:
        text = (ROOT / 'bench' / source).read_text().replace('"hexagon_direct_ope_rows.hpp"', '"shared/direct_ope_rows.hpp"')
        text = text.replace('"hexagon_streaming_echelon_v8.hpp"', '"shared/streaming_echelon_v8_comparison.hpp"')
        text = text.replace('"hexagon_streaming_echelon_v11.hpp"', '"shared/streaming_echelon.hpp"')
        text = text.replace('"hexagon_affine_seam_discovery_v10.hpp"', '"shared/affine_seam_discovery.hpp"')
        if name == 'mapped_probe_v11':
            text = text.replace('irr::require(expanded==a,"Mapped in-place expansion");',
                'if(!expanded)expanded=s_realloc(a,4096);a=expanded;irr::require(a,"Mapped expansion failed");')
        (out / 'tests' / (name + '.cpp')).write_text(text)
    command = [a.cxx, '-O3', '-std=c++20', '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
               '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
               '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
               '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(command, flush=True)
    subprocess.run(command, cwd=out, check=True)
