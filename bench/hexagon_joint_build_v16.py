#!/usr/bin/env python3
"""Build exact MRK common-factor reuse and bounded RHS prefetch from V15."""
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
for name in ['hexagon_mrk_rhs_reuse.hpp','hexagon_mrk_rhs_batch.inc']:
    shutil.copy2(ROOT / 'bench' / name, out / 'src/shared' / name)
projection = out / 'src/reduced_d3/irr_projection.hpp'
code = projection.read_text().replace('#include <set>', '#include <set>\n#include "shared/hexagon_mrk_rhs_reuse.hpp"')
code = code.replace('    std::optional<nmhv_physical_symbol_table_t> mrk_symbols;',
    '    std::optional<nmhv_physical_symbol_table_t> mrk_symbols;\n    std::map<std::tuple<ulong,uint64_t,size_t,size_t,size_t>,std::vector<ulong>> mrk_rhs_cache;')
marker = '// Solution rows own the destination.'
assert code.count(marker)==1
code = code.replace(marker, '#include "shared/hexagon_mrk_rhs_batch.inc"\n\n'+marker)
projection.write_text(code)
for path in [projection, out / 'src/shared/direct_physical_rows.hpp']:
    code = path.read_text()
    before = 'equations rows(c.samples*d);auto rhsdata=mrk?&symbols(c):nullptr;'
    after = 'equations rows(c.samples*d);auto rhs_values=mrk?mrk_rhs_reuse_batch(c,degree):std::vector<std::vector<ulong>>{};'
    assert code.count(before)==1
    code = code.replace(before, after)
    before = 'nmhv_physical_mrk_rhs(w/2,nmhv_physical_sample_mrk_slots(c.seed,c.sample_offset+s,w-1,c.F),*rhsdata,degree,c.F)'
    assert code.count(before)==1
    path.write_text(code.replace(before, 'rhs_values.at(s)'))
driver = (out / 'tests/joint_physics_v15.cpp').read_text().replace(
    'joint-discovery-physical-affine-v15-bounded-triangular-ope',
    'joint-discovery-physical-affine-v16-shared-mrk-shuffles')
driver = driver.replace('{"ope_storage",',
    '{"mrk_rhs_engine","shared_one_loop_log_shuffle_v1"},{"mrk_rhs_prefetch","4_for_single_sample;bounded_cache_128"},{"ope_storage",')
(out / 'tests/joint_physics_v16.cpp').write_text(driver)
shutil.copy2(ROOT / 'bench/hexagon_mrk_rhs_batch_probe.cpp', out / 'tests/mrk_rhs_batch_probe.cpp')
prefix = Path('/tmp/symbology-audit-toolchain/env')
for name in ['mrk_rhs_batch_probe', 'joint_physics_v16']:
    command = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
               '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
               '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
               '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
               '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(command, flush=True)
    subprocess.run(command, cwd=out, check=True)
