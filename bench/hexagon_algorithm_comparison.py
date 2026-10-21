#!/usr/bin/env python3
"""Compare the downloaded homogeneous sewing routes with frozen V17.

Use fresh solve roots, shared read-only recurrences and bounded serial jobs.
The existing exact lower-loop fixtures are validation inputs only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--package', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--threads', type=int, default=4)
p.add_argument('--resume', action='store_true')
p.add_argument('--eight-reference', type=Path)
a = p.parse_args()
pkg = a.package.resolve()
out = a.output.resolve()
out.mkdir(exist_ok=a.resume, parents=True)
fec = ROOT / 'output_hexagon_16p4_20261006/FEC'
lec = ROOT / 'output_hexagon_16p4_20261006/LEC/MHV'
v17 = ROOT / 'output_hexagon_17p3_20261007/joint-v17-schedule-package/bin/joint_physics_v17'
prime = 1152921504606847009
env = dict(os.environ, THREADS=str(a.threads), PHYSICS_CACHE_MIB='256',
           OPE_HEAD_CACHE_MIB='512', HEXAGON_MHV_DOUBLE_BONUS='0',
           LD_LIBRARY_PATH=str(ROOT/'output_hexagon_runtime_recovery_20261008/env/lib'))
for key in ('HEXAGON_MAPPED_FILE', 'HEXAGON_MAPPED_GIB', 'HEXAGON_MAPPED_SOFT_GIB'):
    env.pop(key, None)
records = []
meta = {'schema': 'hexagon-algorithm-comparison-v1', 'started': time.time(),
        'threads': a.threads, 'rss_guard_gib': 3.5, 'virtual_guard_gib': 32,
        'timing_scope': 'Shared workstation; other chat physical checks have priority. Correctness and feasibility, not isolated speedup.',
        'bonus': False, 'records': records,
        'binary_sha256': {str(x): hashlib.sha256(x.read_bytes()).hexdigest()
                          for x in (pkg/'bin/irr_pipeline', v17, pkg/'bin/pipeline_compare')}}
if a.resume:
    old = json.loads((out/'results.json').read_text())
    assert old['binary_sha256'] == meta['binary_sha256']
    meta = old
    records = meta['records']
    meta['resumed'] = time.time()

def save():
    tmp = out / 'results.tmp'
    tmp.write_text(json.dumps(meta, indent=2) + '\n')
    tmp.replace(out / 'results.json')

def run(label, cmd, override=None, timeout=900):
    previous = [x for x in records if x['label']==label]
    if any(x['state']=='passed' for x in previous):
        return
    artifact = label if not previous else label+'-attempt'+str(len(previous)+1)
    assert not (out/(artifact+'.json')).exists()
    record = {'label': label, 'command': list(map(str,cmd)), 'state': 'running'}
    records.append(record)
    meta['state'] = 'running_' + label
    save()
    command = [sys.executable, str(ROOT/'bench/measure_process.py'),
               '--report', str(out/(artifact+'.json')), '--log', str(out/(artifact+'.log')),
               '--status', str(out/(artifact+'.status.json')), '--timeout', str(timeout),
               '--virtual-gib', '32', '--rss-gib', '3.5', '--min-available-gib', '4',
               '--', *map(str,cmd)]
    proc = subprocess.run(command, env=dict(env, **(override or {})), cwd=pkg,
                          stdout=subprocess.DEVNULL)
    result = json.loads((out/(artifact+'.json')).read_text())
    record['artifact_label'] = artifact
    record.update(state='passed' if result['exit_code']==0 else 'failed', measurement=result)
    save()
    if result['exit_code'] or proc.returncode:
        raise RuntimeError(label + ' failed; inspect retained log')
    print(label, round(result['wall_seconds'],3), result['peak_summed_rss_bytes'], flush=True)

plan = out/'physics.tsv'
ope = out/'ope.tsv'
ope.write_text('coverage\tMHV-small-S-T1T2-all\n')
plan.write_text('collinear\tcollinear\t-\nmrk\tmrk\t-\nope\tope_small_s\tope.tsv\n')
try:
    for scheduling in ('parallel','legacy'):
        for test in ('test_sew_solver','test_sew_replay'):
            run(scheduling+'-'+test, [pkg/'bin'/test, out/(scheduling+'-'+test)],
                {'IRR_SEW_SCHEDULING': scheduling})
    # Exact historical coefficient fixtures independently check the imported chain.
    blocks = ('pe_A1','pe_A2','pe_E','po_A1','po_A2','po_E')
    for weight in (3,4):
        for block in blocks:
            name = f'irrFEC_{weight}_{block}.wxf'
            cmd = [pkg/'bin/irr_tensor_compare', pkg/f'tests/reference/FEC/w{weight:02d}'/name,
                   fec/f'w{weight:02d}'/name]
            subprocess.run(list(map(str,cmd)), env=env, check=True, stdout=subprocess.DEVNULL)
    meta['independent_exact_fixture_comparisons'] = 12
    save()
    strategies = {'full': [], 'stream': ['--sew-solver','stream','--sew-reuse','plan'],
                  'G': ['--sew-solver','staged','--sew-reuse','plan','--sew-partials','disk'],
                  'S': ['--sew-solver','staged','--sew-reuse','stages','--sew-partials','disk']}
    for h in (5,9):
        roots = {}
        for mode, flags in strategies.items():
            root = out/f'L{(h+3)//2}-{mode}'
            roots[mode] = root/f'SEW_{h}p3/MHV/collinear_mrk_ope_small_s'
            common = [pkg/'bin/irr_pipeline','MHV',h,3,root,'--fec-root',fec,
                      '--lec-root',lec,'--data',pkg/'data','--plan',plan,
                      '--prereduce','right',*flags]
            if flags:
                common += ['--sew-batch-rows','256']
            if mode in ('G','S'):
                common += ['--sew-target-fraction','0.05']
            for count in (1,2):
                run(f'L{(h+3)//2}-{mode}-prime{count}',
                    [*common,'--cache-only','--primes',count])
        native = out/f'L{(h+3)//2}-V17'
        for count in (1,2):
            run(f'L{(h+3)//2}-V17-prime{count}',
                [v17,h,3,fec,lec,pkg/'data',native,64,'cache','right','--prime-index',count])
        for mode in strategies:
            for pp in (prime,1152921504606847067):
                run(f'L{(h+3)//2}-{mode}-compare-{pp}',
                    [pkg/'bin/pipeline_compare', roots[mode],native/'joint',h,3,fec,lec,pp])
        # Export complete rational L4/L6 with two additional validation primes.
        for mode in ('full','S'):
            root = out/f'L{(h+3)//2}-{mode}'
            flags = strategies[mode]
            common = [pkg/'bin/irr_pipeline','MHV',h,3,root,'--fec-root',fec,
                      '--lec-root',lec,'--data',pkg/'data','--plan',plan,
                      '--prereduce','right',*flags]
            if flags:
                common += ['--sew-batch-rows','256','--sew-target-fraction','0.05']
            run(f'L{(h+3)//2}-{mode}-QQ',common)
        run(f'L{(h+3)//2}-V17-QQ', [v17,h,3,fec,lec,pkg/'data',native,64,'reconstruct','right'])
        for mode in ('full','S'):
            run(f'L{(h+3)//2}-{mode}-compare-QQ',
                [pkg/'bin/pipeline_compare',roots[mode],native/'joint',h,3,fec,lec,'QQ'])
    if a.eight_reference:
        h = 13
        root = out/'L8-S'
        run('L8-S-prime1', [pkg/'bin/irr_pipeline','MHV',h,3,root,
            '--fec-root',fec,'--lec-root',lec,'--data',pkg/'data','--plan',plan,
            '--prereduce','right','--sew-solver','staged','--sew-reuse','stages',
            '--sew-partials','disk','--sew-batch-rows','256',
            '--sew-target-fraction','0.05','--cache-only','--primes','1'])
        run('L8-S-compare-reference', [pkg/'bin/pipeline_compare',
            root/'SEW_13p3/MHV/collinear_mrk_ope_small_s',
            a.eight_reference.resolve(),h,3,fec,lec,prime])
        meta['eight_loop_reference'] = str(a.eight_reference.resolve())
    meta['state'] = 'complete'
    meta.pop('error', None)
except BaseException as exc:
    meta['state'] = 'failed'
    meta['error'] = repr(exc)
    raise
finally:
    meta['finished'] = time.time()
    save()
