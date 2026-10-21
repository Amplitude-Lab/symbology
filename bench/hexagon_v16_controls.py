#!/usr/bin/env python3
"""Build and compare V16 sequentially while preserving the frozen ten-loop pilot."""
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
INPUT = ROOT / 'output_hexagon_16p4_20261006'
OUT = BASE / 'v16-controls'
PACKAGE = BASE / 'joint-v16-mrk-reuse-package'
UNIT = 'symbology-hexagon-ten-v15-20261007.service'
active = None
cancelled = False

def stop(signum, frame):
    global cancelled
    cancelled = True
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)

def main():
    global active
    proof = BASE / 'mrk-rhs-reuse/control/progress.json'
    if json.loads(proof.read_text())['status'] != 'all_coefficients_equal':
        raise RuntimeError('Exact MRK controls must finish first')
    if OUT.exists() or PACKAGE.exists():
        raise RuntimeError('V16 controls require fresh directories')
    OUT.mkdir()
    result = {'status': 'starting', 'started': time.time(), 'steps': [],
              'rhs_controls_sha256': hashlib.sha256(proof.read_bytes()).hexdigest(),
              'timing_scope': 'sequential; pilot frozen; compilation completed before benchmarks'}
    def save():
        p = OUT / 'progress.tmp'
        p.write_text(json.dumps(result, indent=2)+'\n')
        p.replace(OUT / 'progress.json')
    def measured(label, command, timeout=1800):
        global active
        if cancelled:
            raise RuntimeError('V16 controls cancelled')
        cmd = [sys.executable, str(ROOT / 'bench/measure_process.py'), '--report', str(OUT / (label+'.json')),
               '--log', str(OUT / (label+'.log')), '--status', str(OUT / (label+'.status.json')),
               '--timeout', str(timeout), '--virtual-gib', '96', '--rss-gib', '7',
               '--min-available-gib', '3', '--', *map(str, command)]
        step = {'name': label, 'command': cmd, 'state': 'running'}
        result['steps'].append(step)
        result['status'] = 'running_' + label
        save()
        with (OUT / (label+'.manager.log')).open('w') as stream:
            active = subprocess.Popen(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            code = active.wait()
            active = None
        report = json.loads((OUT / (label+'.json')).read_text())
        step.update(state='passed' if not code and not report['exit_code'] else 'failed', report=report)
        save()
        if code or report['exit_code'] or cancelled:
            raise RuntimeError('V16 control failed: '+label)
    frozen = False
    try:
        subprocess.run(['systemctl', '--user', 'freeze', UNIT], check=True)
        frozen = True
        result['freeze_started'] = time.time()
        measured('build', [sys.executable, ROOT / 'bench/hexagon_joint_build_v16.py', '--package',
                           BASE / 'joint-v15-triangular-package', '--output', PACKAGE])
        measured('batch-probe', ['env', 'THREADS=8', PACKAGE / 'bin/mrk_rhs_batch_probe',
                                  INPUT / 'FEC', INPUT / 'LEC/MHV', INPUT / 'data', OUT / 'batch-probe-data'])
        reference = (BASE / 'eight-13p3-control/joint/cache/p1152921504606847009.bin').read_bytes()[72:]
        for label, package, binary, threads in [
                ('v15-eight-t8', BASE / 'joint-v15-triangular-package', 'joint_physics_v15', 8),
                ('v16-eight-t8', PACKAGE, 'joint_physics_v16', 8),
                ('v16-eight-t16', PACKAGE, 'joint_physics_v16', 16)]:
            measured(label, ['env', f'THREADS={threads}', 'PHYSICS_CACHE_MIB=4096', 'OPE_HEAD_CACHE_MIB=2048',
                             'MIMALLOC_PURGE_DELAY=0', 'HEXAGON_MHV_DOUBLE_BONUS=1',
                             f'HEXAGON_MAPPED_FILE={OUT / (label+".mmap")}', 'HEXAGON_MAPPED_GIB=32',
                             'HEXAGON_MAPPED_SOFT_GIB=1', package / 'bin' / binary, 13, 3,
                             INPUT / 'FEC', INPUT / 'LEC/MHV', INPUT / 'data', OUT / label,
                             2048, 'cache', 'right', '--prime-index', 1])
            payload = (OUT / label / 'joint/cache/p1152921504606847009.bin').read_bytes()[72:]
            if payload != reference:
                raise RuntimeError('Known eight-loop coefficient mismatch: '+label)
            result[label] = {'every_coefficient_equal': True, 'payload_sha256': hashlib.sha256(payload).hexdigest()}
            save()
        result.update(status='controls_passed', finished=time.time())
    except Exception as error:
        result.update(status='failed', error=str(error), finished=time.time())
    finally:
        if frozen:
            thaw = subprocess.run(['systemctl', '--user', 'thaw', UNIT], capture_output=True, text=True)
            result.update(thaw_exit_code=thaw.returncode, freeze_finished=time.time())
            if thaw.returncode:
                result.update(status='failed', thaw_error=thaw.stderr)
        save()
    return 0 if result['status']=='controls_passed' else 1

if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
