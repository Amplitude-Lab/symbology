#!/usr/bin/env python3
"""Validate a separate scheduling repair while giving the live pilot priority."""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
INPUT = ROOT / 'output_hexagon_16p4_20261006'
OUT = BASE / 'v17-schedule-controls'
PACKAGE = BASE / 'joint-v17-schedule-package'
PILOT = 'symbology-hexagon-ten-v16-pilot-20261007.service'
OWN = 'symbology-hexagon-v17-controls-20261007.service'
active = None
cancelled = False


def stop(signum, frame):
    global cancelled
    cancelled = True
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)


def memory_path(unit):
    value = subprocess.check_output(
        ['systemctl', '--user', 'show', unit, '--property=ControlGroup', '--value'], text=True).strip()
    return Path('/sys/fs/cgroup' + value) / 'memory.current' if value else None


def main():
    global active
    if (OUT / 'progress.json').exists() or PACKAGE.exists():
        raise RuntimeError('Validation requires fresh outputs; retain prior evidence')
    pilot_memory, own_memory = memory_path(PILOT), memory_path(OWN)
    result = {'status': 'starting', 'started': time.time(), 'steps': [],
              'timing_scope': 'Correctness controls overlap the productive ten-loop pilot; not an isolated speed comparison.',
              'control_memory_max_bytes': 3_000_000_000,
              'combined_sampled_stop_bytes': 46_000_000_000,
              'production_selection': 'not selected; active V16 package and runner remain unchanged'}

    def save():
        p = OUT / 'progress.tmp'
        p.write_text(json.dumps(result, indent=2) + '\n')
        p.replace(OUT / 'progress.json')

    def aggregate_memory():
        return sum(int(p.read_text()) for p in (pilot_memory, own_memory) if p and p.exists())

    def measured(label, command, timeout=900):
        global active
        if cancelled or aggregate_memory() >= 43_000_000_000:
            raise RuntimeError('Control deferred to preserve production memory headroom')
        cmd = [sys.executable, str(ROOT / 'bench/measure_process.py'),
               '--report', str(OUT / (label + '.json')), '--log', str(OUT / (label + '.log')),
               '--status', str(OUT / (label + '.status.json')), '--timeout', str(timeout),
               '--virtual-gib', '96', '--rss-gib', '2.7', '--min-available-gib', '5',
               '--', *map(str, command)]
        step = {'name': label, 'state': 'running', 'command': cmd}
        result['steps'].append(step)
        result['status'] = 'running_' + label
        save()
        stopped_for_headroom = False
        with (OUT / (label + '.manager.log')).open('w') as stream:
            active = subprocess.Popen(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            while active.poll() is None:
                combined = aggregate_memory()
                result['peak_combined_sampled_bytes'] = max(result.get('peak_combined_sampled_bytes', 0), combined)
                if combined >= 46_000_000_000:
                    stopped_for_headroom = True
                    active.send_signal(signal.SIGTERM)
                    break
                time.sleep(.1)
            code = active.wait()
            active = None
        report = json.loads((OUT / (label + '.json')).read_text())
        step.update(state='passed' if not code and not report['exit_code'] else 'failed',
                    stopped_for_production_headroom=stopped_for_headroom, report=report)
        save()
        if code or report['exit_code'] or cancelled or stopped_for_headroom:
            raise RuntimeError('Control failed or deferred: ' + label)

    try:
        measured('cursor', [OUT / 'cursor_probe'])
        measured('build', [sys.executable, ROOT / 'bench/hexagon_joint_build_v17.py',
                          '--package', BASE / 'joint-v16-mrk-reuse-package', '--output', PACKAGE])
        label = 'eight-t4'
        measured(label, ['env', 'THREADS=4', 'PHYSICS_CACHE_MIB=4096', 'OPE_HEAD_CACHE_MIB=2048',
                         'MIMALLOC_PURGE_DELAY=0', 'HEXAGON_MHV_DOUBLE_BONUS=1',
                         'HEXAGON_MAPPED_FILE=' + str(OUT / 'eight.mmap'),
                         'HEXAGON_MAPPED_GIB=32', 'HEXAGON_MAPPED_SOFT_GIB=1',
                         PACKAGE / 'bin/joint_physics_v17', 13, 3,
                         INPUT / 'FEC', INPUT / 'LEC/MHV', INPUT / 'data', OUT / label,
                         2048, 'cache', 'right', '--prime-index', 1])
        reference = (BASE / 'eight-13p3-control/joint/cache/p1152921504606847009.bin').read_bytes()[72:]
        candidate = (OUT / label / 'joint/cache/p1152921504606847009.bin').read_bytes()[72:]
        if candidate != reference:
            raise RuntimeError('Known eight-loop coefficient mismatch')
        log = (OUT / (label + '.log')).read_text()
        if 'operator=ope_lowgrade4' not in log or 'JOINT_PRIME_COMPLETE' not in log:
            raise RuntimeError('Required early OPE or final completion evidence is absent')
        result.update(status='cursor_and_eight_controls_passed', finished=time.time(),
                      every_eight_loop_coefficient_equal=True,
                      payload_sha256=hashlib.sha256(candidate).hexdigest(),
                      next_validation='Review full-scale pilot and nine-loop regression before production selection')
    except Exception as error:
        result.update(status='failed_or_deferred', error=str(error), finished=time.time())
    save()
    return 0 if result['status'] == 'cursor_and_eight_controls_passed' else 1


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
