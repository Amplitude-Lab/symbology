#!/usr/bin/env python3
"""Queue the nine-loop scheduler regression after a verified V16 pilot exit."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
INPUT = ROOT / 'output_hexagon_16p4_20261006'
CONTROL = BASE / 'ten-v16-runner'
OUT = BASE / 'v17-nine-control'
PACKAGE = BASE / 'joint-v17-schedule-package'
PILOT = 'symbology-hexagon-ten-v16-pilot-20261007.service'
PRIME = 1152921504606847009
active = None
cancelled = False


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic(path, value):
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def stop(signum, frame):
    global cancelled
    cancelled = True
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)


def service():
    query = subprocess.run(['systemctl', '--user', 'show', PILOT,
        '--property=LoadState,ActiveState,SubState,MainPID,Result,ExecMainStatus'],
        text=True, capture_output=True)
    state = dict(line.split('=', 1) for line in query.stdout.splitlines())
    if query.returncode and state.get('LoadState') != 'not-found':
        raise RuntimeError('Cannot inspect pilot service: ' + query.stderr)
    return state


def compare_cache(candidate, reference):
    if sha(candidate) != Path(str(candidate) + '.sha256').read_text().strip():
        raise RuntimeError('Control cache checksum mismatch')
    actual, expected = candidate.read_bytes(), reference.read_bytes()
    # The identity fingerprint differs for separate validated algorithms.
    if actual[:64] != expected[:64] or actual[72:] != expected[72:]:
        raise RuntimeError('Control cache header or coefficient mismatch')


def preflight():
    evidence = BASE / 'v17-schedule-controls/progress.json'
    if json.loads(evidence.read_text())['status'] != 'cursor_and_eight_controls_passed':
        raise RuntimeError('Scheduler and complete eight-loop controls have not passed')
    manifest = json.loads((CONTROL / 'manifest.json').read_text())
    pinned = {
        ROOT / 'bench/hexagon_ten_loop_runner_v16.py': manifest['runner_sha256'],
        ROOT / 'bench/measure_process.py': manifest['measure_sha256'],
        Path(manifest['driver']): manifest['driver_sha256'],
        Path(manifest['converter']): manifest['converter_sha256'],
    }
    for path, digest in pinned.items():
        if sha(path) != digest:
            raise RuntimeError('Pinned production identity changed: ' + str(path))
    references = [BASE / 'v13-controls/nine/joint/cache' / f'p{PRIME}.bin',
                  ROOT / 'output_hexagon_14p4_20261006/cached-nine/joint/cache' / f'p{PRIME}.bin']
    for path in references:
        if sha(path) != Path(str(path) + '.sha256').read_text().strip():
            raise RuntimeError('Known nine-loop cache checksum mismatch')
    return {'v16_pinned_files': {str(k): v for k, v in pinned.items()},
            'v17_driver_sha256': sha(PACKAGE / 'bin/joint_physics_v17'),
            'eight_control_sha256': sha(evidence),
            'nine_reference_sha256': {str(p): sha(p) for p in references}}, references


def main():
    global active
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--preflight-only', action='store_true')
    options = args.parse_args()
    identity, references = preflight()
    if options.preflight_only:
        print(json.dumps({'preflight_passed': True, **identity}, indent=2))
        return 0
    if (OUT / 'progress.json').exists():
        raise RuntimeError('Prior control exists; preserve and inspect it before another attempt')
    OUT.mkdir(exist_ok=True)
    result = {'status': 'waiting_for_pilot_terminal_result', 'started': time.time(),
              'identity': identity, 'steps': [], 'wait_deadline': time.time() + 7200,
              'production_selection': 'V17 not selected; nine-loop regression only',
              'memory_limit_bytes': 20_000_000_000}

    def save():
        atomic(OUT / 'progress.json', result)

    def measured(label, command):
        global active
        if cancelled:
            raise RuntimeError('Control cancelled')
        cmd = [sys.executable, str(ROOT / 'bench/measure_process.py'),
               '--report', str(OUT / (label + '.json')), '--log', str(OUT / (label + '.log')),
               '--status', str(OUT / (label + '.status.json')), '--timeout', '5400',
               '--virtual-gib', '96', '--rss-gib', str(19_000_000_000 / 1024**3),
               '--min-available-gib', '5', '--', *map(str, command)]
        step = {'name': label, 'command': cmd, 'state': 'running'}
        result['steps'].append(step)
        result['status'] = 'running_' + label
        save()
        with (OUT / (label + '.manager.log')).open('w') as stream:
            active = subprocess.Popen(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            code = active.wait()
            active = None
        report = json.loads((OUT / (label + '.json')).read_text())
        step.update(state='passed' if not code and not report['exit_code'] else 'failed', report=report)
        save()
        if code or report['exit_code'] or cancelled:
            raise RuntimeError('Nine-loop control failed: ' + label)

    lock = None
    save()
    try:
        while True:
            if cancelled or time.time() >= result['wait_deadline']:
                raise RuntimeError('Cancelled or pilot wait deadline reached; inspect before resuming')
            state = service()
            if state['ActiveState'] not in ('active', 'activating', 'deactivating'):
                break
            time.sleep(5)
        # Share the production supervisor lock for serialization, without
        # modifying its contents or bypassing any stage producer lock.
        lock = (CONTROL / 'runner.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if preflight()[0] != identity:
            raise RuntimeError('Validation inputs changed during the wait')
        progress = json.loads((CONTROL / 'progress.json').read_text())
        if progress['state'] != 'prime_budget_exhausted_cache_only':
            raise RuntimeError('Pilot did not exit in the expected verified cache-only state; review required')
        if state.get('LoadState') != 'not-found' and int(state['ExecMainStatus']) != 2:
            raise RuntimeError('Unexpected pilot service exit code')
        spec = importlib.util.spec_from_file_location('checked_v16_runner', ROOT / 'bench/hexagon_ten_loop_runner_v16.py')
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        ledger = json.loads((CONTROL / 'ledger.json').read_text())
        entry = ledger['primes']['1']
        native = runner.check_cache(BASE / 'ten-v16-main/joint', PRIME, 17, 3, 275246)
        shifted = runner.check_cache(BASE / 'ten-v16-shifted-16p4', PRIME, 16, 4, 518048)
        if native != entry['native'] or shifted != entry['conversion']:
            raise RuntimeError('Pilot committed caches changed')
        if (BASE / 'ten-v16-shifted-16p4' / f'source_p{PRIME}.sha256').read_text().strip() != native['sha256']:
            raise RuntimeError('Pilot conversion source receipt mismatch')
        attempt = ledger['attempts'][entry['native_attempt']]
        report = json.loads(Path(attempt['report']).read_text())
        log = Path(attempt['log']).read_text()
        if report['exit_code'] or any(report.get(k) for k in ['timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal']):
            raise RuntimeError('Pilot native execution report failed')
        for condition in ['collinear', 'mrk', 'ope_small_s']:
            if f'JOINT_PHYSICAL_CHECK operator={condition} passed=1 remaining=0' not in log:
                raise RuntimeError('Pilot independent physical check missing: ' + condition)
        if 'JOINT_SEAM_COMPLETE full_operator=1' not in log:
            raise RuntimeError('Pilot original seam completion missing')
        result['pilot_review'] = {'service': state, 'native': native, 'shifted': shifted,
            'native_report': report, 'original_seam_and_fresh_physical_checks_passed': True,
            'scope': 'One conditional finite-field result; no rational ten-loop result yet.'}
        save()
        measured('nine-native', ['env', 'THREADS=16', 'PHYSICS_CACHE_MIB=4096', 'OPE_HEAD_CACHE_MIB=2048',
            'MIMALLOC_PURGE_DELAY=0', 'HEXAGON_MHV_DOUBLE_BONUS=1',
            'HEXAGON_MAPPED_FILE=' + str(OUT / 'nine.mmap'), 'HEXAGON_MAPPED_GIB=48',
            'HEXAGON_MAPPED_SOFT_GIB=12', PACKAGE / 'bin/joint_physics_v17', 15, 3,
            INPUT / 'FEC', INPUT / 'LEC/MHV', INPUT / 'data', OUT / 'nine', 6000, 'cache', 'right', '--prime-index', 1])
        native_path = OUT / 'nine/joint/cache' / f'p{PRIME}.bin'
        compare_cache(native_path, references[0])
        measured('nine-shift', ['env', 'THREADS=16', INPUT / 'joint-v8-package/bin/irr_convert',
            'shift-cache', 'MHV', 15, 3, OUT / 'nine/joint', OUT / 'nine-shifted-14p4',
            '--head', 14, '--fec-root', INPUT / 'FEC', '--lec-root', INPUT / 'LEC/MHV', '--kind', 'unique'])
        shifted_path = OUT / 'nine-shifted-14p4/cache' / f'p{PRIME}.bin'
        compare_cache(shifted_path, references[1])
        result.update(status='nine_loop_regression_passed', finished=time.time(),
            every_native_coefficient_equal=True, every_shifted_coefficient_equal=True,
            native_payload_sha256=hashlib.sha256(native_path.read_bytes()[72:]).hexdigest(),
            shifted_payload_sha256=hashlib.sha256(shifted_path.read_bytes()[72:]).hexdigest(),
            next_action='Review production selection and an explicit cache-provenance plan. Preserve the completed V16 prime.')
    except Exception as error:
        result.update(status='stopped_for_review', error=str(error), finished=time.time())
    finally:
        save()
        if lock is not None:
            lock.close()
    return 0 if result['status'] == 'nine_loop_regression_passed' else 1


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    sys.exit(main())
