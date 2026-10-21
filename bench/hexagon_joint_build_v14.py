#!/usr/bin/env python3
"""Preserve V13 mathematics and release unused FLINT integer caches after large QQ inputs."""
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
shutil.copy2(ROOT / 'bench/hexagon_native_cache_release.hpp', out / 'src/shared/hexagon_native_cache_release.hpp')
projection = out / 'src/reduced_d3/irr_projection.hpp'
code = projection.read_text().replace('#include <set>', '#include <set>\n#include "shared/hexagon_native_cache_release.hpp"')
old = '        return value;\n    }\n    void load('
new = '''        if(bytes>=size_t(64)*1024*1024){
            K.reset(); // Release temporary QQ ownership before freeing only unused integer caches.
            hexagon_native_memory::release_unused_integer_caches(&opt->pool);
        }
        return value;
    }
    void load('''
assert code.count(old) == 1
projection.write_text(code.replace(old, new))
driver = (out / 'tests/joint_physics_v13.cpp').read_text()
driver = driver.replace('joint-discovery-physical-affine-v13-bounded-interleaving',
    'joint-discovery-physical-affine-v14-bounded-interleaving-cache-release')
driver = driver.replace('''#ifdef USE_MIMALLOC
        mi_collect(true); // Only unused native allocator pages may be returned.
#endif''', '''        hexagon_native_memory::release_unused_integer_caches(&opt->pool);''')
assert 'release_unused_integer_caches(&opt->pool)' in driver
(out / 'tests/joint_physics_v14.cpp').write_text(driver)
shutil.copy2(ROOT / 'bench/hexagon_native_flint_collect_probe.cpp', out / 'tests/native_flint_collect_probe.cpp')
prefix = Path('/tmp/symbology-audit-toolchain/env')
for name in ['native_flint_collect_probe', 'joint_physics_v14']:
    cmd = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
           '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
           '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
           '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
           '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(cmd, flush=True)
    subprocess.run(cmd, cwd=out, check=True)
