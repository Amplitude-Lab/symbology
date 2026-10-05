#!/usr/bin/env python3
"""Sequential before/after audit of equivariant condition contraction.

Uses fresh prepared data, three independent processes per method by default,
the unchanged native bootstrap, and the unpacked hexagon reference executables.
Validates every saved recurrence array against the preceding implementation,
then checks exact recursive changes of basis to independent references.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import shutil
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--before', type=pathlib.Path, required=True)
p.add_argument('--native', type=pathlib.Path, required=True)
p.add_argument('--archive', type=pathlib.Path, required=True)
p.add_argument('--references', type=pathlib.Path, required=True,
               help='prior audit containing hexagon-reference-fec and -lec seeds')
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repeats', type=int, default=3)
args = p.parse_args()
assert args.repeats > 0
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
archive, references = args.archive.resolve(), args.references.resolve()
snapshot = out / 'source'
snapshot.mkdir()
for src in [ROOT / 'symrep', ROOT / 'symrep.cpp', *ROOT.glob('symrep*.hpp'),
            pathlib.Path(__file__), ROOT / 'bench/native_symrep_reference.cpp']:
    shutil.copy2(src, snapshot / src.name)
program = snapshot / 'symrep'
env = dict(os.environ, THREADS='2')
report = dict(threads=2, repeats=args.repeats, cases=[], processes=[],
              source_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                             for f in snapshot.iterdir()},
              notes=['Sequential fresh processes; OS file caches not flushed.',
                     'Preparation and reference verification excluded from chain wall time.',
                     'Both symmetry versions use the same freshly prepared bundle.',
                     'All chain times include startup, mandatory checks, expansion and WXF exports.',
                     'Archive uses its own SparseRREF checkout.'])


def save():
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


def run(label, command, cwd=ROOT):
    usage = out / (label + '.usage.json')
    start = time.perf_counter()
    with (out / (label + '.log')).open('w') as stream:
        subprocess.run(['/usr/bin/time', '-f', '{"peak_rss_kib":%M}', '-o', str(usage),
                        *map(str, command)], cwd=cwd, env=env,
                       stdout=stream, stderr=subprocess.STDOUT, check=True)
    record = dict(label=label, wall_s=time.perf_counter() - start,
                  **json.loads(usage.read_text()))
    report['processes'].append(record)
    print(record, flush=True)
    save()
    return record


for polygon, kind, maximum in [('hexagon', 'fec', 8), ('hexagon', 'lec', 6),
                                ('heptagon', 'fec', 4), ('heptagon', 'lec', 3)]:
    name = f'{polygon}-{kind}'
    direction = 'forward' if kind == 'fec' else 'backward'
    data = archive / 'data' if polygon == 'hexagon' else ROOT / 'data'
    seed = (references / f'hexagon-reference-{kind}/w1.wxf' if polygon == 'hexagon'
            else data / f'{kind.upper()}_1.wxf')
    condition = data / (('dlogmat_ca_1term_pd3.wxf' if kind == 'fec' else
                         'dlogmat_full_ca_pd3.wxf') if polygon == 'hexagon' else 'dlogmat_E6.wxf')
    parity = references / 'hexagon-reference-fec/parity.wxf' if polygon == 'hexagon' else data / 'parityrepmat.wxf'
    gens = [('cyclic', data / 'cycrepmat.wxf'), ('flip', data / 'fliprepmat.wxf'), ('parity', parity)]
    prepared = out / f'{name}-prepared'
    case = dict(name=name, maximum=maximum, runs=[])
    report['cases'].append(case)
    case['preparation'] = run(f'{name}-prepare', [program, 'prepare',
        *[v for key, path in gens for v in ('--generator', f'{key}={path}')],
        '--condition', condition, '--seed', seed, '--direction', direction,
        '--output', prepared, '--threads', 2])
    for repeat in range(args.repeats):
        row = {}
        methods = ['native', 'before', 'after'] + (['archive'] if polygon == 'hexagon' else [])
        for method in methods:
            label = f'{name}-{method}-{repeat}'
            dest = out / label
            if method == 'native':
                command = [args.native.resolve(), direction, condition, seed, maximum, dest, 2]
            elif method == 'archive':
                command = [archive / f'bin/irr_{kind}', *(['MHV'] if kind == 'lec' else []), maximum, dest]
            else:
                command = [args.before.resolve() if method == 'before' else program,
                           'extend', '--input', prepared, '--max-weight', maximum,
                           '--output', dest, '--threads', 2]
            if method == 'after':
                command += ['--backend', 'multiplicity']
            row[method] = run(label, command, archive if method == 'archive' else ROOT)
            if method != 'archive':
                with (dest / 'timings.tsv').open() as stream:
                    row[method]['phases'] = list(csv.DictReader(stream, delimiter='\t'))
        for f in (out / f'{name}-before-{repeat}').glob('w*'):
            assert f.read_bytes() == (out / f'{name}-after-{repeat}' / f.name).read_bytes(), f
        row['all_recurrence_arrays_byte_identical'] = True
        case['runs'].append(row)
        save()
    case['median_wall_s'] = {m: statistics.median(r[m]['wall_s'] for r in case['runs']) for m in methods}
    refs = [('native', out / f'{name}-native-0')]
    if polygon == 'hexagon':
        expanded = out / f'{name}-archive-expanded'
        run(f'{name}-export', [ROOT / 'bench/hexagon_symrep_reference', kind.upper(),
                              out / f'{name}-archive-0', maximum, archive / 'data', expanded])
        refs.append(('archive', expanded))
    for label, ref in refs:
        run(f'{name}-{label}-verify', [program, 'verify', '--input', prepared,
            '--chain', out / f'{name}-after-0', '--reference', ref, '--max-weight', maximum,
            '--output', out / f'{name}-{label}-verified', '--threads', 2])
    case['exact_references_passed'] = [label for label, _ in refs]
    save()
print('ALL EXACT CONDITION BENCHMARK CHECKS PASSED', flush=True)
