#!/usr/bin/env python3
"""Isolated exact MRK controls; freeze/thaw the numerical pilot around timing."""
import hashlib
import json
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
OUT = BASE / 'mrk-rhs-reuse/control'
BINARY = BASE / 'mrk-rhs-reuse/bin/mrk_rhs_probe'
DATA = ROOT / 'output_hexagon_16p4_20261006/data'
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
    if OUT.exists():
        raise RuntimeError('Refusing to replace an existing MRK control')
    OUT.mkdir(parents=True)
    result = {'status': 'starting', 'started': time.time(), 'cases': [],
              'binary_sha256': hashlib.sha256(BINARY.read_bytes()).hexdigest(),
              'source_sha256': hashlib.sha256((ROOT / 'bench/hexagon_mrk_rhs_reuse.hpp').read_bytes()).hexdigest(),
              'pilot_unit': UNIT, 'timing_scope': 'sequential controls with pilot cgroup frozen'}
    def save():
        temp = OUT / 'progress.tmp'
        temp.write_text(json.dumps(result, indent=2) + '\n')
        temp.replace(OUT / 'progress.json')
    save()
    frozen = False
    try:
        state = subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'ActiveState', '--value'], text=True).strip()
        if state != 'active':
            raise RuntimeError('Expected an active pilot to freeze for isolated timing')
        subprocess.run(['systemctl', '--user', 'freeze', UNIT], check=True)
        frozen = True
        result['freeze_started'] = time.time()
        prime1, prime2 = 1152921504606847009, 1152921504606847067
        cases = [(L, p, s, 2*L-1) for L in range(1, 7) for p in (101, prime1) for s in (0, 7)]
        cases += [(6, 101, 3, 0), (6, prime1, 3, 3), (7, prime1, 3, 13),
                  (8, 101, 0, 15), (8, prime1, 7, 15), (9, prime1, 3, 17),
                  (10, prime1, 0, 19), (10, prime2, 7, 19)]
        for index, (L, prime, sample, degree) in enumerate(cases):
            if cancelled:
                raise RuntimeError('MRK validation cancelled')
            label = f'{index:02d}-L{L}-p{prime}-s{sample}-d{degree}'
            report, log = OUT / (label + '.json'), OUT / (label + '.log')
            command = [sys.executable, str(ROOT / 'bench/measure_process.py'), '--report', str(report),
                       '--log', str(log), '--timeout', '900', '--virtual-gib', '16', '--rss-gib', '6',
                       '--min-available-gib', '3', '--', str(BINARY), str(DATA), str(L), str(prime),
                       str(sample), str(degree), 'compare']
            result.update(status='running', active_case=label)
            save()
            with (OUT / (label + '.manager.log')).open('w') as stream:
                active = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
                code = active.wait()
                active = None
            measurement = json.loads(report.read_text())
            text = log.read_text()
            if code or measurement['exit_code'] or 'MRK_PROBE PASS mode=compare' not in text:
                raise RuntimeError('MRK comparison failed: ' + label)
            times = dict(re.findall(r'MRK_(OLD|NEW).*?seconds=([0-9.e+\-]+)', text))
            result['cases'].append({'loop': L, 'prime': prime, 'sample': sample, 'degree': degree,
                                    'every_coefficient_equal': True, 'seconds': times,
                                    'peak_rss_bytes': measurement['peak_summed_rss_bytes'], 'log': str(log)})
            save()
        result.update(status='all_coefficients_equal', finished=time.time())
    except Exception as error:
        result.update(status='failed', error=str(error), finished=time.time())
    finally:
        if frozen:
            thaw = subprocess.run(['systemctl', '--user', 'thaw', UNIT], capture_output=True, text=True)
            result.update(thaw_exit_code=thaw.returncode, freeze_finished=time.time())
            if thaw.returncode:
                result.update(status='failed', thaw_error=thaw.stderr)
        save()
    return 0 if result['status'] == 'all_coefficients_equal' else 1

if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
