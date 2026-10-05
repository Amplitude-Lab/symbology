#!/usr/bin/env python3
"""Fresh, sequential NMHV LEC benchmarks; original and streamed in one binary.

Build bench/nmhv_kernel_benchmark and bench/nmhv_kernel_check first. Generate
the heptagon seed with `nmhv_kernel_check fixture DIR`. The hexagon seed is
w1.wxf exported by hexagon_symrep_reference from the archive's NMHV chain.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--hexagon', type=Path, required=True)
p.add_argument('--hexagon-seed', type=Path, required=True)
p.add_argument('--heptagon-seed', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--program', type=Path, default=ROOT/'bench/nmhv_kernel_benchmark')
p.add_argument('--before-program', type=Path,
               help='Compare two revisions of the streamed kernel instead of original versus streamed')
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--threads', type=int, default=8)
p.add_argument('--timeout', type=float, default=300)
a = p.parse_args()
assert a.repeats > 0 and a.threads > 0
out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
program = a.program.resolve()
before_program = a.before_program.resolve() if a.before_program else None
methods = ('before', 'after') if before_program else ('original', 'streamed')
cases = [('hexagon', a.hexagon.resolve()/'data/dlogmat_full_ca_pd3.wxf', a.hexagon_seed.resolve(), 5),
         ('heptagon', ROOT/'data/dlogmat_E6.wxf', a.heptagon_seed.resolve(), 4)]
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
report = {'threads': a.threads, 'executable_sha256': sha(program),
          'inputs': {str(x): sha(x) for _, d, s, _ in cases for x in (d, s)}, 'runs': []}
if before_program:
    report['before_executable_sha256'] = sha(before_program)

def limits():
    resource.setrlimit(resource.RLIMIT_AS, (28*1024**3, 28*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

for repeat in range(a.repeats):
    for polygon, condition, seed, weight in cases:
        for method in (methods if repeat % 2 == 0 else methods[::-1]):
            label = f'{polygon}-{method}-{repeat}'
            executable = before_program if method == 'before' else program
            strategy = 'streamed' if before_program else method
            command = [executable, 'extend', strategy, condition, seed, weight, out/label, a.threads]
            start = time.perf_counter()
            with (out/f'{label}.log').open('w') as log:
                with subprocess.Popen(['/usr/bin/time', '-f', '%M', '-o', str(out/f'{label}.rss'),
                                       *map(str, command)], stdout=log, stderr=subprocess.STDOUT,
                                      preexec_fn=limits, start_new_session=True, cwd=ROOT) as child:
                    try:
                        code = child.wait(timeout=a.timeout)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                        raise
                    if code:
                        raise subprocess.CalledProcessError(code, command)
            record = dict(case=polygon, method=method, repeat=repeat,
                          wall_s=time.perf_counter()-start,
                          peak_rss_kib=int((out/f'{label}.rss').read_text()))
            report['runs'].append(record)
            (out/'results.json').write_text(json.dumps(report, indent=2)+'\n')
            print(record, flush=True)
for polygon, *_ in cases:
    for method in methods:
        runs = [r for r in report['runs'] if r['case'] == polygon and r['method'] == method]
        print(polygon, method, {key: statistics.median(r[key] for r in runs)
                              for key in ('wall_s', 'peak_rss_kib')})
