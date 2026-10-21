#!/usr/bin/env python3
"""Fault-injection controls for restart recovery, including frozen V27D at L8."""
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from hexagon_boot_recovery import atomic_json, flush_completed, recover_files, sha

REPO = Path(__file__).resolve().parents[1]


def small_controls():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / 'owned'
        archive = Path(tmp) / 'retained'
        cache = root / 'stage/cache'
        cache.mkdir(parents=True)
        good = cache / 'p101.csr'
        good.write_bytes(b'completed coefficients')
        good.with_name(good.name + '.sha256').write_text(sha(good) + '\n')
        (cache / 'p103.csr').write_bytes(b'incomplete transaction')
        (cache / 'p107.csr.tmp').write_bytes(b'incomplete write')
        (cache / 'p109.csr.sha256').write_text('0' * 64 + '\n')
        retained = recover_files(root, archive)
        assert len(retained) == 3 and good.exists()
        assert all(Path(p).exists() for p in retained)
        records = flush_completed([root], Path(tmp) / 'durable.json')
        assert records[str(good)]['sha256'] == sha(good)
        good.write_bytes(b'corrupted completed coefficients')
        try:
            recover_files(root, archive)
        except RuntimeError as error:
            assert 'Checksum mismatch' in str(error) and good.exists()
        else:
            raise AssertionError('corruption was not rejected')
        lockpath = Path(tmp) / 'lock'
        with lockpath.open('a') as first, lockpath.open('a') as second:
            fcntl.flock(first, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                fcntl.flock(second, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                pass
            else:
                raise AssertionError('duplicate owner accepted')
        with lockpath.open('a') as new:
            fcntl.flock(new, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return {'missing_commit_retained': True, 'temporary_retained': True,
                'orphan_retained': True, 'corruption_rejected_without_deletion': True,
                'durability_receipt_verified': True, 'exclusive_lock_and_stale_path': True}


def integration():
    base = REPO / 'output_hexagon_ten_v27d_20261009/recovery-control'
    base.mkdir(parents=True, exist_ok=False)
    out = base / 'L8'
    profile = json.loads((REPO / 'output_hexagon_v27_20261009/current-nine-loop-profile.json').read_text())
    env = {k: v for k, v in os.environ.items() if not k.startswith(('HEXAGON_', 'SOLVER_'))}
    env.update({k: str(v) for k, v in profile['environment'].items()})
    for key in profile['unset_environment']:
        env.pop(key, None)
    shared = REPO / 'output_hexagon_16p4_20261006'
    cmd = [profile['binary'], '13', '3', str(shared / 'FEC'), str(shared / 'LEC/MHV'),
           str(shared / 'data'), str(out), '2048', 'cache', 'right', '--prime-index', '1']
    target = out / 'discovery/cache/p1152921504606847009.csr.sha256'
    with (base / 'interrupted.log').open('w') as stream:
        child = subprocess.Popen(cmd, env=env, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        start = time.monotonic()
        while child.poll() is None and not target.exists() and time.monotonic() - start < 120:
            time.sleep(.01)
        if child.poll() is not None or not target.exists():
            raise RuntimeError('Did not reach interruptible completed checkpoint')
        os.killpg(child.pid, signal.SIGKILL)
        child.wait()
    first_digest = sha(Path(str(target)[:-7]))
    retained = recover_files(out, base / 'retained')
    flush_completed([out], base / 'durable.json')
    with (base / 'resumed.log').open('w') as stream:
        child = subprocess.run(cmd, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=120)
    if child.returncode:
        raise RuntimeError('Resumed lower-loop solve failed; inspect resumed.log')
    actual = out / 'joint/cache/p1152921504606847009.bin'
    reference = REPO / 'output_hexagon_v27_20261009/suite-d/L8-V27/V27/joint/cache/p1152921504606847009.bin'
    text = (base / 'resumed.log').read_text()
    assert 'PROFILE phase=sparse_checkpoint_read' in text
    assert sha(Path(str(target)[:-7])) == first_digest
    assert actual.read_bytes() == reference.read_bytes()
    for op in ('collinear', 'mrk', 'ope_small_s'):
        assert f'JOINT_PHYSICAL_CHECK operator={op} passed=1 remaining=0' in text
    result = {'small_controls': small_controls(), 'frozen_driver_sha256': sha(profile['binary']),
              'actual_process_killed_after_discovery_checkpoint': True,
              'completed_checkpoint_unchanged_after_resume': True,
              'resumed_every_native_coefficient_equal': True, 'all_fresh_physical_checks_passed': True,
              'retained_files': retained, 'log': str(base / 'resumed.log'),
              'scope': 'SIGKILL/relaunch fault injection; machine reboot not performed'}
    atomic_json(base / 'passed.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    if '--small-only' in sys.argv:
        print(json.dumps(small_controls(), indent=2))
    else:
        integration()
