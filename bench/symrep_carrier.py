#!/usr/bin/env python3
"""Sequential native/multiplicity/factorized audit, with exact chain validation.

Preparation, factorized-to-expanded export and verification have separate
measurements. Larger heptagon weights test whether the improvement persists.
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
p.add_argument('--references', type=pathlib.Path, required=True)
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
            pathlib.Path(__file__), ROOT / 'bench/native_symrep_reference.cpp',
            ROOT / 'tests/symrep_probe.cpp', ROOT / 'tests/symrep_regression.py']:
    shutil.copy2(src, snapshot / src.name)
program = snapshot / 'symrep'
env = {k: v for k, v in os.environ.items() if not k.startswith('SYMREP_')}
env['THREADS'] = '2'
report = dict(threads=2, repeats=args.repeats, cases=[], processes=[],
              source_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in snapshot.iterdir()},
              executable_sha256={str(f.resolve()): hashlib.sha256(f.read_bytes()).hexdigest()
                                 for f in [args.before, args.native]},
              notes=['Sequential fresh processes; OS file caches not flushed.',
                     'Every method includes startup and its own output representation.',
                     'Factorized output is C plus H, where B=H*C is the irreducible basis.',
                     'The factorized backend uses full sparse elimination internally.',
                     'Preparation, expansion and optional verification are timed separately.',
                     'auto selects factorized for nonsplit rational irreps; multiplicity otherwise.',
                     'Archive has its own SparseRREF and writes compact multiplicities.'])


def save():
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


def run(label, command, cwd=ROOT):
    usage = out / (label + '.usage.json')
    start = time.perf_counter()
    with (out / (label + '.log')).open('w') as stream:
        subprocess.run(['/usr/bin/time', '-f', '{"peak_rss_kib":%M}', '-o', str(usage),
                        *map(str, command)], cwd=cwd, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    record = dict(label=label, wall_s=time.perf_counter() - start, **json.loads(usage.read_text()))
    report['processes'].append(record)
    print(record, flush=True)
    save()
    return record


prepared_cases = {}
for polygon, kind, maximum, extended in [
    ('heptagon', 'fec', 4, False), ('heptagon', 'lec', 3, False),
    ('hexagon', 'fec', 8, False), ('hexagon', 'lec', 6, False),
    ('heptagon', 'fec', 5, True), ('heptagon', 'lec', 4, True),
]:
    name = f'{polygon}-{kind}' + (f'-w{maximum}' if extended else '')
    direction = 'forward' if kind == 'fec' else 'backward'
    data = archive / 'data' if polygon == 'hexagon' else ROOT / 'data'
    seed = (references / f'hexagon-reference-{kind}/w1.wxf' if polygon == 'hexagon'
            else data / f'{kind.upper()}_1.wxf')
    condition = data / (('dlogmat_ca_1term_pd3.wxf' if kind == 'fec' else
                         'dlogmat_full_ca_pd3.wxf') if polygon == 'hexagon' else 'dlogmat_E6.wxf')
    parity = references / 'hexagon-reference-fec/parity.wxf' if polygon == 'hexagon' else data / 'parityrepmat.wxf'
    gens = [('cyclic', data / 'cycrepmat.wxf'), ('flip', data / 'fliprepmat.wxf'), ('parity', parity)]
    case = dict(name=name, maximum=maximum, runs=[])
    report['cases'].append(case)
    key = (polygon, kind)
    if key not in prepared_cases:
        prepared = out / f'{name}-prepared'
        case['preparation'] = run(f'{name}-prepare', [program, 'prepare',
            *[v for key, path in gens for v in ('--generator', f'{key}={path}')],
            '--condition', condition, '--seed', seed, '--direction', direction,
            '--output', prepared, '--threads', 2])
        prepared_cases[key] = prepared
    prepared = prepared_cases[key]
    for repeat in range(args.repeats):
        row = {}
        methods = ['native'] + ([] if extended else ['before']) + ['after'] + (['archive'] if polygon == 'hexagon' else [])
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
            row[method] = run(label, command, archive if method == 'archive' else ROOT)
            if method != 'archive':
                with (dest / 'timings.tsv').open() as stream:
                    row[method]['phases'] = list(csv.DictReader(stream, delimiter='\t'))
        case['runs'].append(row)
        save()
        run(f'{name}-verify-{repeat}', [program, 'verify', '--input', prepared,
            '--chain', out / f'{name}-after-{repeat}', '--reference', out / f'{name}-native-{repeat}',
            '--max-weight', maximum, '--output', out / f'{name}-verified-{repeat}', '--threads', 2])
    if extended:
        case['full_comparison'] = run(f'{name}-full-comparison', [program, 'extend', '--input', prepared,
            '--max-weight', maximum, '--output', out / f'{name}-full-comparison', '--threads', 2, '--compare'])
    case['median_wall_s'] = {m: statistics.median(r[m]['wall_s'] for r in case['runs']) for m in methods}
    case['native_exact_reference_all_runs'] = True
    if polygon == 'hexagon':
        reference = out / f'{name}-archive-expanded'
        run(f'{name}-archive-export', [ROOT / 'bench/hexagon_symrep_reference', kind.upper(),
            out / f'{name}-archive-0', maximum, archive / 'data', reference])
        run(f'{name}-archive-verify', [program, 'verify', '--input', prepared,
            '--chain', out / f'{name}-after-0', '--reference', reference, '--max-weight', maximum,
            '--output', out / f'{name}-archive-verified', '--threads', 2])
        case['archive_exact_reference'] = True
    elif not extended:
        expanded = out / f'{name}-expanded'
        case['expansion'] = run(f'{name}-expand', [program, 'expand', '--input', prepared,
            '--chain', out / f'{name}-after-0', '--output', expanded, '--threads', 2])
        run(f'{name}-expanded-verify', [program, 'verify', '--input', prepared,
            '--chain', expanded, '--reference', out / f'{name}-native-0', '--max-weight', maximum,
            '--output', out / f'{name}-expanded-verified', '--threads', 2])
        case['expanded_exact_reference'] = True
    save()
print('PASS all sequential timing and exact reference checks', flush=True)
