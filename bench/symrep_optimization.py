#!/usr/bin/env python3
"""Four-way hexagon audit: native bootstrap, old/new general solver, archive.

Run sequentially; all chains are fresh. Exact recursive changes of basis are
checked separately from timed solves. Build native_symrep_reference.cpp with
the same flags/libraries as symrep and pass its executable via --native.
"""
import argparse, csv, hashlib, json, os, pathlib, shutil, statistics, subprocess, time

ROOT = pathlib.Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--archive', type=pathlib.Path, required=True)
p.add_argument('--before', type=pathlib.Path, required=True)
p.add_argument('--native', type=pathlib.Path, required=True)
p.add_argument('--prepared', type=pathlib.Path, required=True, help='previous validated benchmark directory')
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--fec-weight', type=int, default=8)
p.add_argument('--lec-weight', type=int, default=5)
args = p.parse_args()
out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
archive = args.archive.resolve(); env = dict(os.environ, THREADS='2')
records = []; summary = {'processes': records, 'cases': [], 'threads': 2, 'repeats': args.repeats}
snapshot = out / 'source'; snapshot.mkdir()
for src in [ROOT / 'symrep', *ROOT.glob('symrep*.hpp'), ROOT / 'symrep.cpp', pathlib.Path(__file__), ROOT / 'bench/native_symrep_reference.cpp']:
    shutil.copy2(src, snapshot / src.name)
summary['sha256'] = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in snapshot.iterdir()}
program = snapshot / 'symrep'

def run(label, command, cwd=ROOT):
    usage = out / (label + '.usage.json'); start = time.perf_counter()
    with (out / (label + '.log')).open('w') as stream:
        subprocess.run(['/usr/bin/time', '-f', '{"peak_rss_kib":%M}', '-o', str(usage), *map(str, command)],
                       cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
    rec = dict(label=label, wall_s=time.perf_counter()-start, **json.loads(usage.read_text()))
    records.append(rec); print(rec, flush=True)
    (out / 'results.json').write_text(json.dumps(summary, indent=2)+'\n')
    return rec

def phases(directory):
    with (directory / 'timings.tsv').open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

for kind, maximum in [('fec', args.fec_weight), ('lec', args.lec_weight)]:
    prepared = args.prepared.resolve() / f'hexagon-{kind}-0-prepared'
    seed = args.prepared.resolve() / f'hexagon-reference-{kind}/w1.wxf'
    direction = 'forward' if kind == 'fec' else 'backward'
    condition = archive / 'data' / ('dlogmat_ca_1term_pd3.wxf' if kind == 'fec' else 'dlogmat_full_ca_pd3.wxf')
    case = dict(kind=kind, maximum=maximum, runs=[]); summary['cases'].append(case)
    for repeat in range(args.repeats):
        row = {}
        for method in ['native', 'before', 'after', 'archive']:
            label = f'{kind}-{method}-{repeat}'; dest = out / label
            if method == 'native': command = [args.native.resolve(), direction, condition, seed, maximum, dest, 2]
            elif method == 'archive': command = [archive / 'bin' / f'irr_{kind}', *(['MHV'] if kind == 'lec' else []), maximum, dest]
            else: command = [args.before.resolve() if method == 'before' else program, 'extend', '--input', prepared, '--max-weight', maximum, '--output', dest, '--threads', 2]
            row[method] = run(label, command, archive if method == 'archive' else ROOT)
            if method != 'archive': row[method]['phases'] = phases(dest)
        case['runs'].append(row)
        reference = out / f'{kind}-archive-expanded-{repeat}'
        run(f'{kind}-export-{repeat}', [ROOT/'bench/hexagon_symrep_reference', kind.upper(), out/f'{kind}-archive-{repeat}', maximum, archive/'data', reference])
        for version in ['before', 'after']:
            for ref_name, ref in [('native',out/f'{kind}-native-{repeat}'), ('archive',reference)]:
                label = f'{kind}-{version}-{ref_name}-verify-{repeat}'
                run(label, [program, 'verify', '--input', prepared, '--chain', out/f'{kind}-{version}-{repeat}', '--reference', ref, '--max-weight', maximum, '--output', out/label, '--threads', 2])
    case['median_wall_s'] = {m: statistics.median(r[m]['wall_s'] for r in case['runs']) for m in ['native', 'before', 'after', 'archive']}
    (out/'results.json').write_text(json.dumps(summary, indent=2)+'\n')
print('ALL EXACT REFERENCES PASSED', flush=True)
