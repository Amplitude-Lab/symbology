#!/usr/bin/env python3
"""Bounded V11 controls, preserving the requested baseline freezer state."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]
BASE = REPO / 'output_hexagon_17p3_20261007'
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, default=BASE / 'v11-validation')
p.add_argument('--package', type=Path, default=BASE / 'joint-v11-package')
p.add_argument('--keep-baseline-frozen', action='store_true')
a = p.parse_args()
OUT = a.output.resolve()
PACKAGE = a.package.resolve()
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

def run(argv):
    return subprocess.run(argv, text=True, capture_output=True, check=True)

def state():
    return run(['systemctl', '--user', 'show', OLD_UNIT, '--property=FreezerState', '--value']).stdout.strip()

def main():
    global active
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ['audit', 'log', 'scratch']:
        (OUT / name).mkdir(exist_ok=True)
    lock = (OUT / 'validation.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    result = {'status': 'starting', 'started_unix_seconds': time.time(),
              'preserved_unit': OLD_UNIT, 'steps': [], 'v11_production_started': False}
    atomic(OUT / 'progress.json', result)
    main_pid = run(['systemctl', '--user', 'show', OLD_UNIT, '--property=MainPID', '--value']).stdout.strip()
    command = Path('/proc', main_pid, 'cmdline').read_bytes()
    expected_state = 'frozen' if a.keep_baseline_frozen else 'running'
    if b'hexagon_ten_loop_runner_v10.py' not in command or state() != expected_state:
        raise RuntimeError('Expected this task\'s V10 supervisor with state ' + expected_state)
    result['baseline_should_remain_frozen'] = a.keep_baseline_frozen
    frozen = False

    def measured(label, arguments, timeout=900, mapped=False):
        global active
        if cancelled:
            raise RuntimeError('Validation cancelled')
        command = [sys.executable, str(REPO / 'bench/measure_process.py'),
                   '--report', str(OUT / 'audit' / (label + '.json')),
                   '--log', str(OUT / 'log' / (label + '.log')),
                   '--status', str(OUT / 'audit' / (label + '.status.json')),
                   '--timeout', str(timeout), '--virtual-gib', '96',
                   '--rss-gib', '24', '--min-available-gib', '3', '--',
                   'env', 'THREADS=8', 'PHYSICS_CACHE_MIB=4096',
                   'OPE_HEAD_CACHE_MIB=2048', 'MIMALLOC_PURGE_DELAY=0']
        if mapped:
            command += ['HEXAGON_MAPPED_FILE=' + str(OUT / 'scratch' / (label + '.mmap')),
                        'HEXAGON_MAPPED_GIB=32', 'HEXAGON_MAPPED_SOFT_GIB=1']
        command += list(map(str, arguments))
        entry = {'name': label, 'command': command, 'state': 'running'}
        result['steps'].append(entry)
        result['status'] = 'running_' + label
        atomic(OUT / 'progress.json', result)
        with (OUT / 'log' / (label + '.manager.log')).open('w') as stream:
            active = subprocess.Popen(command, cwd=REPO, stdout=stream, stderr=subprocess.STDOUT)
            code = active.wait()
            active = None
        report = json.loads((OUT / 'audit' / (label + '.json')).read_text())
        entry.update(state='passed' if code == 0 and report['exit_code'] == 0 else 'failed', report=report)
        atomic(OUT / 'progress.json', result)
        if code or report['exit_code'] or cancelled:
            raise RuntimeError('Validation failed: ' + label)

    try:
        run(['systemctl', '--user', 'freeze', OLD_UNIT])
        frozen = True
        if state() != 'frozen':
            raise RuntimeError('V10 service did not freeze')
        result['frozen_unix_seconds'] = time.time()
        result['memory_before_reclaim'] = int((CGROUP / 'memory.current').read_text())
        # Reclaim only file pages from the frozen task. Mapped coefficient
        # bytes remain in its scratch file; no live allocation is modified.
        for amount in [34_359_738_368, 8_589_934_592]:
            if int((CGROUP / 'memory.current').read_text()) <= 16_000_000_000:
                break
            try:
                (CGROUP / 'memory.reclaim').write_text(str(amount))
            except OSError as error:
                result.setdefault('reclaim_diagnostics', []).append(str(error))
        result['memory_after_reclaim'] = int((CGROUP / 'memory.current').read_text())
        if result['memory_after_reclaim'] > 16_000_000_000:
            raise RuntimeError('Insufficient reclaimed cache to fit validation within total50GB budget')
        atomic(OUT / 'progress.json', result)
        measured('build', [sys.executable, REPO / 'bench/hexagon_joint_build_v11.py',
            '--package', REPO / 'output_hexagon_16p4_20261006/joint-v8-package',
            '--output', PACKAGE, '--cxx', '/tmp/symbology-audit-toolchain/env/bin/x86_64-conda-linux-gnu-c++',
            '--runtime-prefix', '/tmp/symbology-audit-toolchain/env'], timeout=900)
        measured('capacity', [PACKAGE / 'bin/capacity_probe_v11'], mapped=True)
        measured('mapped', [PACKAGE / 'bin/mapped_probe_v11'], mapped=True)
        measured('affine', [PACKAGE / 'bin/affine_probe_v11'])
        fec = REPO / 'output_hexagon_16p4_20261006/FEC'
        lec = REPO / 'output_hexagon_16p4_20261006/LEC/MHV'
        data = REPO / 'output_hexagon_16p4_20261006/data'
        measured('direct-ope', [PACKAGE / 'bin/direct_ope_probe_v11', 5, 3,
            fec, lec, data, OUT / 'direct-ope-diagnostics', 1152921504606847009], mapped=True)
        measured('eight', [PACKAGE / 'bin/joint_physics_v11', 13, 3,
            fec, lec, data, OUT / 'eight', 2048, 'cache', 'right', '--prime-index', 1],
            timeout=1800, mapped=True)
        source = OUT / 'eight/joint/cache/p1152921504606847009.bin'
        reference = BASE / 'eight-v10-13p3-control/joint/cache/p1152921504606847009.bin'
        left, right = source.read_bytes()[72:], reference.read_bytes()[72:]
        if left != right:
            raise RuntimeError('V11 eight-loop coefficient payload differs')
        result['eight_coefficient_payload_equal'] = True
        result['eight_payload_sha256'] = hashlib.sha256(left).hexdigest()
        result['status'] = 'controls_passed;production_comparison_pending'
    except Exception as error:
        result.update(status='validation_failed', error=str(error))
    finally:
        if frozen and not a.keep_baseline_frozen:
            thaw = subprocess.run(['systemctl', '--user', 'thaw', OLD_UNIT], text=True, capture_output=True)
            result['thaw_exit_code'] = thaw.returncode
            result['thaw_stderr'] = thaw.stderr
            result['thawed_unix_seconds'] = time.time()
            result['post_validation_freezer_state'] = state()
        elif frozen:
            result['post_validation_freezer_state'] = state()
        result['finished_unix_seconds'] = time.time()
        atomic(OUT / 'progress.json', result)
        print(json.dumps(result), flush=True)
    return 0 if result['status'].startswith('controls_passed') else 1

if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
