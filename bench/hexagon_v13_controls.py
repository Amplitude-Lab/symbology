#!/usr/bin/env python3
"""Sequential, bounded known-symbol controls before any V13 ten-loop solve."""
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
OUT = BASE / 'v13-controls'
PACKAGE = BASE / 'joint-v13-interleaved-package'
INPUT = ROOT / 'output_hexagon_16p4_20261006'
FEC, LEC, DATA = INPUT / 'FEC', INPUT / 'LEC/MHV', INPUT / 'data'
OLD_UNIT = 'symbology-hexagon-ten-v10-20261007.service'
CGROUP = Path('/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice') / OLD_UNIT
active = None
cancelled = False

def stop(signum, frame):
    global cancelled
    cancelled = True
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)

def atomic(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)

def main():
    global active
    if OUT.exists():
        raise RuntimeError('Use a fresh controls directory; do not overwrite completed or interrupted work.')
    for path in [OUT, OUT / 'log', OUT / 'audit', OUT / 'scratch']:
        path.mkdir()
    state = subprocess.check_output(['systemctl', '--user', 'show', OLD_UNIT,
                                    '--property=FreezerState', '--value'], text=True).strip()
    if state != 'frozen' or int((CGROUP / 'memory.current').read_text()) > 16_000_000_000:
        raise RuntimeError('Expected preserved baseline frozen with <=16GB charged memory.')
    for path in [BASE / 'audit/v13-build.json', BASE / 'audit/bonus-double-audit-repair.json']:
        if json.loads(path.read_text())['exit_code'] != 0:
            raise RuntimeError('Required build or conditional-tail audit has not passed.')
    result = {'status': 'starting', 'started_unix_seconds': time.time(), 'steps': [],
              'ten_loop_native_primes': 0, 'baseline_kept_frozen': True,
              'conditional_assumptions': 'Dixon-Liu 2308.08199 eq76-79, full D3 orbit; extrapolated to L10'}

    def measured(label, args, timeout=900, mapped=False, soft=12):
        global active
        if cancelled:
            raise RuntimeError('Controls cancelled')
        cmd = [sys.executable, str(ROOT / 'bench/measure_process.py'),
               '--report', str(OUT / 'audit' / (label + '.json')),
               '--log', str(OUT / 'log' / (label + '.log')),
               '--status', str(OUT / 'audit' / (label + '.status.json')),
               '--timeout', str(timeout), '--virtual-gib', '96', '--rss-gib', '30',
               '--min-available-gib', '3', '--', 'env', 'THREADS=8',
               'PHYSICS_CACHE_MIB=4096', 'OPE_HEAD_CACHE_MIB=2048',
               'MIMALLOC_PURGE_DELAY=0', 'HEXAGON_MHV_DOUBLE_BONUS=1']
        if mapped:
            cmd += ['HEXAGON_MAPPED_FILE=' + str(OUT / 'scratch' / (label + '.mmap')),
                    'HEXAGON_MAPPED_GIB=48', 'HEXAGON_MAPPED_SOFT_GIB=' + str(soft)]
        cmd += list(map(str, args))
        entry = {'name': label, 'state': 'running', 'command': cmd}
        result['steps'].append(entry)
        result['status'] = 'running_' + label
        atomic(OUT / 'progress.json', result)
        with (OUT / 'log' / (label + '.manager.log')).open('w') as stream:
            active = subprocess.Popen(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            code = active.wait()
            active = None
        report = json.loads((OUT / 'audit' / (label + '.json')).read_text())
        entry.update(state='passed' if not code and not report['exit_code'] else 'failed', report=report)
        atomic(OUT / 'progress.json', result)
        if code or report['exit_code'] or cancelled:
            raise RuntimeError('Control failed: ' + label)

    def compare(label, actual, expected):
        left, right = actual.read_bytes()[72:], expected.read_bytes()[72:]
        if left != right:
            raise RuntimeError(label + ' differs from independently completed symbol')
        result[label] = {'every_coefficient_equal': True, 'payload_sha256': hashlib.sha256(left).hexdigest(),
                         'actual': str(actual), 'reference': str(expected)}
        atomic(OUT / 'progress.json', result)

    try:
        measured('affine', [PACKAGE / 'bin/affine_probe_v13'])
        measured('eight', [PACKAGE / 'bin/joint_physics_v13', 13, 3, FEC, LEC, DATA,
                          OUT / 'eight', 2048, 'cache', 'right', '--prime-index', 1],
                 timeout=1800, mapped=True, soft=1)
        compare('eight_comparison', OUT / 'eight/joint/cache/p1152921504606847009.bin',
                BASE / 'eight-13p3-control/joint/cache/p1152921504606847009.bin')
        measured('native-collect', [BASE / 'joint-v12-bonus-package/bin/native_collect_probe',
                                   FEC, LEC, OUT / 'native-collect', '--purge-delay-zero'], timeout=900, mapped=True, soft=12)
        measured('nine', [PACKAGE / 'bin/joint_physics_v13', 15, 3, FEC, LEC, DATA,
                         OUT / 'nine', 6000, 'cache', 'right', '--prime-index', 1],
                 timeout=7200, mapped=True, soft=12)
        measured('nine-shift', [INPUT / 'joint-v8-package/bin/irr_convert', 'shift-cache',
                               'MHV', 15, 3, OUT / 'nine/joint', OUT / 'nine-shifted-14p4',
                               '--head', 14, '--fec-root', FEC, '--lec-root', LEC, '--kind', 'unique'])
        compare('nine_comparison', OUT / 'nine-shifted-14p4/cache/p1152921504606847009.bin',
                ROOT / 'output_hexagon_14p4_20261006/cached-nine/joint/cache/p1152921504606847009.bin')
        result['status'] = 'controls_passed;ten_loop_selection_pending'
    except Exception as error:
        result.update(status='control_failed', error=str(error))
    result['finished_unix_seconds'] = time.time()
    atomic(OUT / 'progress.json', result)
    print(json.dumps(result), flush=True)
    return 0 if result['status'].startswith('controls_passed') else 1

if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
