#!/usr/bin/env python3
"""Sequential native/compact-orbit/archive timings with exact reference checks."""
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
p.add_argument('--native', type=pathlib.Path, required=True)
p.add_argument('--prepared-root', type=pathlib.Path, required=True)
p.add_argument('--archive', type=pathlib.Path, required=True)
p.add_argument('--archive-bin', type=pathlib.Path, help='Separate archive build directory')
p.add_argument('--references', type=pathlib.Path, required=True)
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repeats', type=int, default=5)
p.add_argument('--previous', type=pathlib.Path, help='Previous symmetry executable for controlled before/after measurements')
p.add_argument('--validated-bases', type=pathlib.Path,
               help='Reuse exact reference proofs only when prepared checksums and all saved basis bytes match a completed audit')
args = p.parse_args()
assert args.repeats > 0
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
archive, references = args.archive.resolve(), args.references.resolve()
archive_bin = args.archive_bin.resolve() if args.archive_bin else archive / 'bin'
snapshot = out / 'source'
snapshot.mkdir()
for src in [ROOT / 'symrep', ROOT / 'symrep.cpp', *ROOT.glob('symrep*.hpp'),
            pathlib.Path(__file__), ROOT / 'bench/native_symrep_reference.cpp',
            ROOT / 'tests/symrep_probe.cpp', ROOT / 'tests/symrep_regression.py']:
    shutil.copy2(src, snapshot / src.name)
for name in ('SparseRREF', 'patches'):
    dest = snapshot / name
    dest.mkdir()
    for src in (ROOT / name).iterdir():
        if src.is_file() and src.suffix in ('.h', '.hpp', '.patch'):
            shutil.copy2(src, dest / src.name)
program = snapshot / 'symrep'
previous = None
if args.previous:
    previous = snapshot / 'previous-symrep'
    shutil.copy2(args.previous.resolve(), previous)
env = {k: v for k, v in os.environ.items() if not k.startswith('SYMREP_')}
env['THREADS'] = '2'
report = dict(threads=2, repeats=args.repeats, cases=[], processes=[],
    source_sha256={str(f.relative_to(snapshot)): hashlib.sha256(f.read_bytes()).hexdigest() for f in snapshot.rglob('*') if f.is_file()},
    native_sha256=hashlib.sha256(args.native.read_bytes()).hexdigest(),
    notes=['Sequential fresh processes, alternating native/symmetry order; OS caches not flushed.',
           'Times include startup, solving, output and completion checksums.',
           'Prepared bundles are reused. Expanded export and optional verification are separate.',
           'Split models solve and store scalar multiplicities; nonsplit models use certified rational orbit recipes.',
           'The nonsplit factorized backend still uses a full sparse QQ constraint kernel.',
           'Repeated outputs must be byte-identical to the independently verified first run.',
           'Archive executables retain their own libraries, build and compact output format.'])
validated = None
if args.validated_bases:
    validated = args.validated_bases.resolve()
    verified_cases = {c['name']: c for c in json.loads((validated / 'results.json').read_text())['cases']}


def save():
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


def run(label, command, cwd=ROOT):
    usage = out / (label + '.usage.json')
    start = time.perf_counter()
    with (out / (label + '.log')).open('w') as stream:
        subprocess.run(['/usr/bin/time', '-f', '{"peak_rss_kib":%M}', '-o', str(usage),
                        *map(str, command)], cwd=cwd, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    result = dict(label=label, wall_s=time.perf_counter() - start, **json.loads(usage.read_text()))
    report['processes'].append(result)
    print(result, flush=True)
    save()
    return result


def basis_files(directory):
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
            for f in directory.glob('w*') if f.is_file()}


for polygon, kind, maximum in [('heptagon', 'fec', 4), ('heptagon', 'lec', 3),
                               ('heptagon', 'fec', 5), ('heptagon', 'lec', 4),
                               ('hexagon', 'fec', 8), ('hexagon', 'lec', 6)]:
    name = f'{polygon}-{kind}-w{maximum}'
    direction = 'forward' if kind == 'fec' else 'backward'
    data = archive / 'data' if polygon == 'hexagon' else ROOT / 'data'
    seed = (references / f'hexagon-reference-{kind}/w1.wxf' if polygon == 'hexagon'
            else data / f'{kind.upper()}_1.wxf')
    condition = data / (('dlogmat_ca_1term_pd3.wxf' if kind == 'fec' else
                         'dlogmat_full_ca_pd3.wxf') if polygon == 'hexagon' else 'dlogmat_E6.wxf')
    prepared = args.prepared_root.resolve() / f'{polygon}-{kind}-prepared'
    case = dict(name=name, maximum=maximum, runs=[],
                prepared_sha256=hashlib.sha256((prepared / 'checksums.tsv').read_bytes()).hexdigest())
    report['cases'].append(case)
    for repeat in range(args.repeats):
        row = {}
        methods = (['native', 'symmetry'] if repeat % 2 == 0 else ['symmetry', 'native'])
        if polygon == 'hexagon':
            methods += ['archive']
        if previous:
            methods += ['previous']
        # Rotate all methods, including archive/previous, across fresh runs.
        methods = methods[repeat % len(methods):] + methods[:repeat % len(methods)]
        for method in methods:
            label = f'{name}-{method}-{repeat}'
            dest = out / label
            if method == 'native':
                command = [args.native.resolve(), direction, condition, seed, maximum, dest, 2]
            elif method == 'archive':
                command = [archive_bin / f'irr_{kind}', *(['MHV'] if kind == 'lec' else []), maximum, dest]
            else:
                command = [previous if method == 'previous' else program, 'extend', '--input', prepared, '--max-weight', maximum,
                           '--output', dest, '--threads', 2]
            row[method] = run(label, command, archive if method == 'archive' else ROOT)
            row[method]['output_bytes'] = sum(f.stat().st_size for f in dest.rglob('*') if f.is_file())
            if method != 'archive':
                with (dest / 'timings.tsv').open() as stream:
                    row[method]['phases'] = list(csv.DictReader(stream, delimiter='\t'))
                if repeat:
                    assert basis_files(dest) == basis_files(out / f'{name}-{method}-0'), label
        case['runs'].append(row)
        save()
    case['median_wall_s'] = {m: statistics.median(r[m]['wall_s'] for r in case['runs']) for m in methods}
    case['median_peak_rss_kib'] = {m: statistics.median(r[m]['peak_rss_kib'] for r in case['runs']) for m in methods}
    case['median_output_bytes'] = {m: statistics.median(r[m]['output_bytes'] for r in case['runs']) for m in methods}
    if validated:
        prior = verified_cases[name]
        assert prior['prepared_sha256'] == case['prepared_sha256']
        assert prior['native_exact_all_weights']
        for method in ('native', 'symmetry'):
            assert basis_files(out / f'{name}-{method}-0') == basis_files(validated / f'{name}-{method}-0'), (name, method)
        case['native_exact_all_weights'] = True
        case['repeated_bases_byte_identical'] = True
        case['exact_basis_certificate'] = str(validated / 'results.json')
        case['certificate_sha256'] = hashlib.sha256((validated / 'results.json').read_bytes()).hexdigest()
        if polygon == 'hexagon':
            assert prior['archive_exact_all_weights']
            case['archive_exact_all_weights'] = True
        save()
        continue
    run(f'{name}-verify', [program, 'verify', '--input', prepared,
        '--chain', out / f'{name}-symmetry-0', '--reference', out / f'{name}-native-0',
        '--max-weight', maximum, '--output', out / f'{name}-verified', '--threads', 2])
    case['native_exact_all_weights'] = True
    case['repeated_bases_byte_identical'] = True
    if polygon == 'heptagon' and ((kind == 'fec' and maximum == 5) or (kind == 'lec' and maximum == 4)):
        run(f'{name}-full-comparison', [program, 'extend', '--input', prepared,
            '--max-weight', maximum, '--output', out / f'{name}-full-comparison', '--threads', 2, '--compare'])
        case['full_matrix_checks'] = True
    if polygon == 'hexagon':
        reference = out / f'{name}-archive-expanded'
        run(f'{name}-archive-export', [ROOT / 'bench/hexagon_symrep_reference', kind.upper(),
            out / f'{name}-archive-0', maximum, archive / 'data', reference])
        run(f'{name}-archive-verify', [program, 'verify', '--input', prepared,
            '--chain', out / f'{name}-symmetry-0', '--reference', reference, '--max-weight', maximum,
            '--output', out / f'{name}-archive-verified', '--threads', 2])
        case['archive_exact_all_weights'] = True
    save()
print('PASS sequential benchmarks, repeated basis identity and exact reference checks', flush=True)
