#!/usr/bin/env python3
"""Build the optional empirical bonus-row experiment from a validated V11 package."""
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
    p.error('Output must be fresh; preserve previous experiments.')
shutil.copytree(a.package.resolve(), out, ignore=shutil.ignore_patterns('bin', '.git', 'output', 'log'))
(out / 'bin').mkdir()
for name in ['hexagon_mhv_double_bonus.hpp', 'hexagon_mhv_bonus_seed.hpp']:
    shutil.copy2(ROOT / 'bench' / name, out / 'src/shared' / name)
driver = (out / 'tests/joint_physics_v11.cpp').read_text()
driver = driver.replace('#include "shared/adaptive_physical_rows.hpp"',
    '#include "shared/adaptive_physical_rows.hpp"\n#include "shared/hexagon_mhv_bonus_seed.hpp"')
driver = driver.replace('    pipeline run(s);auto root=s.output/"joint";', '''    bool bonus=false;
    if(auto value=std::getenv("HEXAGON_MHV_DOUBLE_BONUS")){
        require(std::string(value)=="0"||std::string(value)=="1","Invalid bonus assumption switch");
        bonus=std::string(value)=="1";
    }
    pipeline run(s);auto root=s.output/"joint";''')
driver = driver.replace('joint-discovery-physical-affine-v11-adaptive-compact',
    'joint-discovery-physical-affine-v12-optional-double-bonus')
driver = driver.replace('{"early_physics",', '''{"empirical_assumptions",bonus?"MHV double-final bonus relations, Dixon-Liu 2308.08199 eq76-79, full D3 orbit; extrapolated to L10":"none"},
        {"early_physics",''')
driver = driver.replace('hexagon_adaptive_physics::DirectSchedule early(run.c,run.s.seed,policy);', '''hexagon_adaptive_physics::DirectSchedule schedule(run.c,run.s.seed,policy);
            auto early=[&](streaming_echelon::Reducer& reducer,const std::vector<int32_t>& coordinate,size_t round){
                if(bonus&&round==0)mhv_double_bonus::seed_original_coordinates(run.c,reducer,coordinate);
                schedule(reducer,coordinate,round);
            };''')
(out / 'tests/joint_physics_v12.cpp').write_text(driver)
for name in ['mhv_double_bonus_audit', 'native_collect_probe']:
    shutil.copy2(ROOT / 'bench' / ('hexagon_' + name + '.cpp'), out / 'tests' / (name + '.cpp'))
prefix = Path('/tmp/symbology-audit-toolchain/env')
for name in ['mhv_double_bonus_audit', 'native_collect_probe', 'joint_physics_v12']:
    cmd = [str(prefix / 'bin/x86_64-conda-linux-gnu-c++'), '-O3', '-std=c++20',
           '-Isrc', '-Isrc/shared', '-Isrc/reduced_d3', '-Isrc/SparseRREF',
           '-I' + str(prefix / 'include'), 'tests/' + name + '.cpp', '-o', 'bin/' + name,
           '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
           '-lflint', '-lmpfr', '-lgmp', '-lmimalloc', '-ltbb', '-pthread']
    print(cmd, flush=True)
    subprocess.run(cmd, cwd=out, check=True)
