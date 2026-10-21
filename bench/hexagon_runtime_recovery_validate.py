#!/usr/bin/env python3
"""Validate the recovered runtime on known results before resuming V17."""
import fcntl
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import time

from hexagon_ten_loop_runner_v17 import sha, atomic_json

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
RECOVERY = ROOT / 'output_hexagon_runtime_recovery_20261008'
CONTROL = BASE / 'ten-v17-runner'
ENV = dict(os.environ, LD_LIBRARY_PATH=str(RECOVERY / 'env/lib'))
PROGRESS = RECOVERY / 'validation-progress.json'
READY = RECOVERY / 'validation-passed.json'


def verify_ready():
    ready = json.loads(READY.read_text())
    assert ready['status'] == 'controls_passed'
    assert ready['validator_sha256'] == sha(Path(__file__))
    for filename, digest in ready['binaries_and_libraries'].items():
        assert sha(Path(filename)) == digest, filename
    for item in ready['comparisons']:
        assert item['every_coefficient_equal']
        assert sha(Path(item['actual'])) == item['actual_sha256']
        assert sha(Path(item['reference'])) == item['reference_sha256']


def main():
    if sys.argv[1:] == ['--verify-ready']:
        verify_ready()
        return
    assert not sys.argv[1:]
    assert not READY.exists()
    locks = []
    for path in [CONTROL / 'runner.lock', BASE / 'ten-v16-runner/runner.lock']:
        lock = path.open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        locks.append(lock)
    binaries = [BASE / 'joint-v17-schedule-package/bin/joint_physics_v17',
                ROOT / 'output_hexagon_16p4_20261006/joint-v8-package/bin/irr_convert',
                ROOT / 'output_hexagon_16p4_20261006/joint-v8-package/bin/origin_check',
                ROOT / 'output_hexagon_certificate_provenance_20261007/bin/hexagon_shifted_joint_certificate']
    identities = {}
    for binary in binaries:
        result = subprocess.run(['ldd', str(binary)], env=ENV, text=True, capture_output=True)
        assert result.returncode == 0 and 'not found' not in result.stdout + result.stderr
        identities[str(binary)] = sha(binary)
        for line in result.stdout.splitlines():
            match = re.search(r'(?:=>\s*)?(/[^\s]+)\s+\(', line)
            if match:
                library = Path(match[1]).resolve()
                identities[str(library)] = sha(library)
    state = {'schema': 'hexagon-recovered-runtime-controls-v1', 'status': 'running',
             'started': time.time(), 'validator_sha256': sha(Path(__file__)),
             'binaries_and_libraries': identities, 'steps': [], 'comparisons': []}
    atomic_json(PROGRESS, state)

    def measured(name, command):
        report = RECOVERY / (name + '.json')
        log = RECOVERY / (name + '.log')
        status = RECOVERY / (name + '.status.json')
        assert not report.exists() and not log.exists()
        state['current_step'] = name
        atomic_json(PROGRESS, state)
        args = ['/usr/bin/python3', str(ROOT / 'bench/measure_process.py'),
                '--report', str(report), '--log', str(log), '--status', str(status),
                '--timeout', '5400', '--virtual-gib', '256', '--rss-gib', '40',
                '--min-available-gib', '3', '--'] + command
        with (RECOVERY / (name + '.manager.log')).open('x') as stream:
            result = subprocess.run(args, env=ENV, stdout=stream, stderr=subprocess.STDOUT)
        measurement = json.loads(report.read_text())
        state['steps'].append({'name': name, 'measurement': measurement,
                               'report': str(report), 'report_sha256': sha(report),
                               'log': str(log), 'log_sha256': sha(log)})
        atomic_json(PROGRESS, state)
        assert result.returncode == measurement['exit_code'] == 0
        assert not any(measurement[k] for k in ('timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal'))

    def compare(actual, reference):
        assert sha(actual) == actual.with_suffix('.bin.sha256').read_text().strip()
        assert sha(reference) == reference.with_suffix('.bin.sha256').read_text().strip()
        a, b = actual.read_bytes(), reference.read_bytes()
        ha, hb = struct.unpack('<8s8Q', a[:72]), struct.unpack('<8s8Q', b[:72])
        assert ha[:8] == hb[:8] and ha[0] == b'NMHVSEW\0' and ha[5] == 1
        assert len(a) == len(b) == 72 + 8 * ha[5] * ha[6]
        assert a[72:] == b[72:]
        state['comparisons'].append({'actual': str(actual), 'actual_sha256': sha(actual),
                                    'reference': str(reference), 'reference_sha256': sha(reference),
                                    'columns': ha[6], 'every_coefficient_equal': True})
        atomic_json(PROGRESS, state)

    try:
        eight = json.loads((BASE / 'v17-schedule-controls/progress.json').read_text())
        eight_command = next(s['report']['command'] for s in eight['steps'] if s['name'] == 'eight-t4')
        eight_command = [s.replace(str(BASE / 'v17-schedule-controls/eight-t4'), str(RECOVERY / 'eight'))
                          .replace(str(BASE / 'v17-schedule-controls/eight.mmap'), str(RECOVERY / 'eight.mmap'))
                          .replace('THREADS=4', 'THREADS=16') for s in eight_command]
        measured('eight-native', eight_command)
        prime = 'p1152921504606847009.bin'
        compare(RECOVERY / 'eight/joint/cache' / prime, BASE / 'v17-schedule-controls/eight-t4/joint/cache' / prime)
        nine = json.loads((BASE / 'v17-nine-control/progress.json').read_text())
        for step in ['nine-native', 'nine-shift']:
            command = next(s['report']['command'] for s in nine['steps'] if s['name'] == step)
            command = [s.replace(str(BASE / 'v17-nine-control/nine'), str(RECOVERY / 'nine')) for s in command]
            measured(step, command)
            child = 'nine/joint/cache' if step == 'nine-native' else 'nine-shifted-14p4/cache'
            compare(RECOVERY / child / prime, BASE / 'v17-nine-control' / child / prime)
        assert len(state['comparisons']) == 3
        state.update(status='controls_passed', finished=time.time(),
                     scope='Same frozen executables; recovered FLINT 3.2.2 and separately pinned compatible runtime libraries. Every known eight-loop native and nine-loop native/shifted coefficient matches. No ten-loop certification follows.')
        atomic_json(PROGRESS, state)
        atomic_json(READY, state)
        verify_ready()
    except Exception as error:
        state.update(status='failed', error=str(error), finished=time.time())
        atomic_json(PROGRESS, state)
        raise
    finally:
        for lock in locks:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()
    subprocess.run(['systemctl', '--user', 'start', 'symbology-hexagon-ten-v17-20261007.service'], check=True)


if __name__ == '__main__':
    main()
