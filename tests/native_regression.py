"""Small exact fixtures plus reproducible E6 baselines, without private projects."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--probe', type=Path, default=ROOT / 'tests/numerical_probe')
args = parser.parse_args()
probe = args.probe.resolve()

def call(argv, expected=0, timeout=120):
    p = subprocess.run(list(map(str, argv)), capture_output=True, text=True, timeout=timeout)
    assert p.returncode == expected, (argv, p.returncode, p.stdout[-2000:], p.stderr[-2000:])
    return p.stdout

print(call([probe]).splitlines()[-1])
with tempfile.TemporaryDirectory(prefix='symbology-regression ') as tmp:
    work = Path(tmp); f = work / 'fixtures'; out = work / 'output'
    call([probe, 'gen', f])
    for solver in ('incremental', 'sampled'):
        common = [ROOT/'bootstrap', '--solve-collinear', '--solver', solver, '--output-dir', out]
        for name, pair, code in [
            ('imported', ['--pair-cond', f/'inconsistent.wxf'], 1),
            ('pair', ['--pair', f/'pair-seed.wxf', f/'pair-rhs.wxf', 'identity'], 1),
            ('consistent', ['--pair-cond', f/'consistent.wxf'], 0),
        ]:
            stem = f'{solver}_{name}'
            call([*common, *pair, '--export-conditions', '--out-stem', stem], code)
            exported = out/'collinear'/f'cond_{stem}.wxf'
            dumped = call([probe, 'dump', exported])
            if code:
                assert 'DIMS 2 3' in dumped and 'ENTRY 1 2 = 1' in dumped, dumped
            call([*common, '--pair-cond', exported, '--out-stem', stem+'_again'], code)
            if not code:
                assert (out/'collinear'/f'sol_{stem}.wxf').read_bytes() == (out/'collinear'/f'sol_{stem}_again.wxf').read_bytes()
            else:
                assert not (out/'collinear'/f'sol_{stem}_again.wxf').exists()
    print('PASS: 6 condition export/reimport tests (both solvers; imported and pair contradictions)')
    for n in ('2junk', '1.5', '-1', '999999999999999999999999999'):
        call([ROOT/'tensor_ops', 'power', f/'swap.wxf', n, out/'bad.wxf'], 1)
    for weight in ('abc', '1/0', '2/00', '1/2junk'):
        call([ROOT/'tensor_add', f/'swap.wxf', f/'swap.wxf', weight, '1', out/'bad.wxf'], 1)
        call([ROOT/'tensor_ops', 'shuf', f/'swap.wxf', f/'swap.wxf', weight, out/'bad.wxf'], 1)
    assert not (out/'bad.wxf').exists()
    call([ROOT/'tensor_add', f/'swap.wxf', f/'swap.wxf', '1/1000000000000000000000000000000', '-1/1000000000000000000000000000000', out/'zero.wxf'])
    print('PASS: strict numeric rejection and arbitrary-precision exact cancellation')
    base = json.loads((ROOT/'audits/baseline-2026-09-09/manifest-computerhs.json').read_text())['crc32']
    # Resume from existing FEC_2/FEC_3 is the regression: do not prebuild FEC_4/5.
    for target in ('SEW_3p1', 'SEW_5p1'):
        call([ROOT/'compute_rhs', '--target', target, '--data-dir', ROOT/'data', '--output-dir', out,
              '--letter-projection', out/'collinear/colprojdiv_w1.wxf'])
    for rel, expected in base.items():
        actual = f'{zlib.crc32((work/rel).read_bytes()):08x}'
        assert actual == expected, (rel, actual, expected)
    print(f'PASS: all {len(base)} two-/three-loop reference tensors match after resume')
