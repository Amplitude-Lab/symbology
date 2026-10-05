#!/usr/bin/env python3
"""Measure exact recursive carrier and direct adapted transformations separately."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--probe', type=Path, required=True)
p.add_argument('--prepared-root', type=Path, required=True)
p.add_argument('--chains', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--repeats', type=int, default=5)
a = p.parse_args()
assert a.repeats > 0
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=False)
shutil.copy2(a.probe, out / 'actions-probe')
probe = out / 'actions-probe'
env = {k: v for k, v in os.environ.items() if not k.startswith('SYMREP_')}
report = dict(repeats=a.repeats, threads=2, cases=[],
              probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest(),
              notes=['Fresh sequential processes; rotated method order; reused prepared inputs and solved bases.',
                     'Construction phases exclude input reading, output, applications and group setup.',
                     'Peak RSS covers the whole process, including I/O and application tests.',
                     'Process wall time includes 12 application batches per weight; it is not construction-only time.',
                     'Application phases report 11-sample medians after one warmup.',
                     'Equal coefficient arrays in different bases need not represent the same physical vector.'])


def save():
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


def matrices(path):
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in path.glob('*.wxf')}


for polygon, kind, maximum in [('heptagon', 'fec', 5), ('heptagon', 'lec', 4),
                               ('hexagon', 'fec', 8), ('hexagon', 'lec', 6)]:
    name = f'{polygon}-{kind}-w{maximum}'
    prepared = a.prepared_root.resolve() / f'{polygon}-{kind}-prepared'
    current = a.chains.resolve() / f'{name}-symmetry-0'
    # The carrier coordinate option is specific to factorized heptagon chains.
    # Ordinary chains provide the recursive raw-coordinate comparison for both polygons.
    methods = ['ordinary', 'adapted'] + (['legacy', 'carrier'] if polygon == 'heptagon' else [])
    case = dict(name=name, runs=[])
    report['cases'].append(case)
    for repeat in range(a.repeats):
        run = {}
        for method in methods[repeat % len(methods):] + methods[:repeat % len(methods)]:
            chain = a.chains.resolve() / f'{name}-native-0' if method == 'ordinary' else current
            label = f'{name}-{method}-{repeat}'
            dest, usage = out / label, out / (label + '.usage')
            command = [str(probe), str(prepared), str(chain), str(maximum),
                       'carrier' if method == 'ordinary' else method, str(dest), '2', 'measure']
            start = time.perf_counter()
            with (out / (label + '.log')).open('w') as log:
                subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(usage), *command],
                               env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
            wall = time.perf_counter() - start
            with (dest / 'timings.tsv').open() as stream:
                phases = list(csv.DictReader(stream, delimiter='\t'))
            run[method] = dict(wall_s=wall, peak_rss_kib=int(usage.read_text()), phases=phases,
                               construction_sum_s=sum(float(x['construction_s']) for x in phases),
                               matrix_bytes=sum(f.stat().st_size for f in dest.glob('*.wxf')))
            if repeat:
                assert matrices(dest) == matrices(out / f'{name}-{method}-0')
            print(label, run[method]['construction_sum_s'], run[method]['peak_rss_kib'], flush=True)
        if polygon == 'heptagon':
            assert matrices(out / f'{name}-legacy-{repeat}') == matrices(out / f'{name}-carrier-{repeat}')
        case['runs'].append(run)
        save()
    case['summary'] = {}
    for method in methods:
        rows = [r[method] for r in case['runs']]
        case['summary'][method] = {key: statistics.median(x[key] for x in rows)
                                  for key in ['wall_s', 'peak_rss_kib', 'construction_sum_s', 'matrix_bytes']}
        case['summary'][method]['final_weight'] = {
            key: statistics.median(float(x['phases'][-1][key]) for x in rows)
            for key in ['dimension', 'construction_s', 'apply_s', 'implicit_apply_s', 'action_nnz']}
    case['repeated_matrices_identical'] = True
    if polygon == 'heptagon':
        case['legacy_and_current_exact_matrices_identical'] = True
    save()
print('PASS recursive transformation measurements and exact matrix identity checks')
