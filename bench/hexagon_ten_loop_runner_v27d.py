#!/usr/bin/env python3
"""Fresh V27D/16 ten-loop solve, reboot recovery, reconstruction and certificates.

Run in the installed user service. Deterministic failures stop for review;
shutdown interruptions resume at verified stage boundaries. No old prime is
imported. Reference coefficients are compared only after a fresh solve.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import signal
import struct
import subprocess
import sys
import time
import uuid
from hexagon_boot_recovery import atomic_json, flush_completed, process_cpu, quarantine, recover_files, sha, sync

REPO = Path(__file__).resolve().parents[1]
EXPERIMENT = REPO / 'output_hexagon_ten_v27d_20261009'
CONTROL = EXPERIMENT / 'runner'
NATIVE = EXPERIMENT / 'native-17p3'
SHIFTED = EXPERIMENT / 'shifted-16p4'
EXTERNAL = EXPERIMENT / 'asymmetric-qq-16p4'
PROFILE = REPO / 'output_hexagon_v27_20261009/current-nine-loop-profile.json'
SHARED = REPO / 'output_hexagon_16p4_20261006'
FEC, LEC, DATA = SHARED / 'FEC', SHARED / 'LEC/MHV', SHARED / 'data'
REFERENCE = REPO / 'output_hexagon_17p3_20261007/ten-v17-main/joint'
CONVERTER = SHARED / 'joint-v8-package/bin/irr_convert'
CERTIFICATE = REPO / 'output_hexagon_certificate_provenance_20261007/bin/hexagon_shifted_joint_certificate'
ORIGIN = SHARED / 'joint-v8-package/bin/origin_check'
CASES = REPO / 'audits/hexagon-ten-loop-2026-10-06/origin-cases.tsv'
MEASURE = REPO / 'bench/measure_process.py'
HEADER = struct.Struct('<8s8Q')
ASSUMPTIONS = 'MHV double-final bonus relations, Dixon-Liu 2308.08199 eq76-79, full D3 orbit; extrapolated to L10'
active = None
cancel_signal = 0


def tsv(path):
    return dict(line.split('\t', 1) for line in path.read_text().splitlines())


def check_cache(root, prime, h, k, columns):
    path = root / 'cache' / f'p{prime}.bin'
    if path.is_symlink() or not path.is_file():
        raise RuntimeError('Missing ordinary cache: ' + str(path))
    digest = sha(path)
    if path.with_name(path.name + '.sha256').read_text().strip() != digest:
        raise RuntimeError('Cache checksum mismatch: ' + str(path))
    with path.open('rb') as stream:
        header = HEADER.unpack(stream.read(HEADER.size))
    expected = (b'NMHVSEW\0', 1, 10, h, k, 1, columns, prime,
                int(sha(root / 'input_identity.txt')[:16], 16))
    if header != expected or path.stat().st_size != HEADER.size + 8 * columns:
        raise RuntimeError('Cache identity, dimensions or uniqueness mismatch: ' + str(path))
    return {'path': str(path), 'sha256': digest, 'prime': prime, 'rows': 1,
            'directions': 0, 'columns': columns}


def check_qq(root=None):
    root = root or (EXTERNAL if EXTERNAL.exists() else SHIFTED)
    receipt = root / 'qq_complete.tsv'
    if not receipt.exists():
        return False
    fields = tsv(receipt)
    if fields.get('schema') != 'irr-physics-qq-v1' or fields.get('identity') != sha(root / 'input_identity.txt'):
        raise RuntimeError('QQ receipt identity mismatch')
    if not {'object.tsv', 'columns.tsv', 'manifest.json', 'stage.tsv'} <= fields.keys():
        raise RuntimeError('Incomplete QQ receipt')
    for name, digest in fields.items():
        if name in ('schema', 'identity'):
            continue
        if Path(name).name != name or name in ('.', '..'):
            raise RuntimeError('Unsafe QQ receipt name')
        path = root / name
        if path.is_symlink() or sha(path) != digest:
            raise RuntimeError('QQ checksum mismatch: ' + name)
    manifest, obj = json.loads((root / 'manifest.json').read_text()), tsv(root / 'object.tsv')
    if (manifest.get('complete') is not True or manifest.get('validation_primes') != 2 or
            [obj.get(k) for k in ('kind', 'directions', 'h', 'k')] != ['unique', '0', '16', '4'] or
            sum(n.endswith(('_origin.wxf', '_origin.wxf.gz')) for n in fields) != 6):
        raise RuntimeError('QQ object/completion scope mismatch')
    return True


def native_receipts(text, index, checkpoint_resumed=False):
    selected = re.findall(r'JOINT_SELECTED_PRIME index=(\d+) prime=(\d+)', text)
    if len(selected) != 1 or int(selected[0][0]) != index:
        raise RuntimeError('Native prime-selection receipt missing')
    prime = int(selected[0][1])
    for name in ('collinear', 'mrk', 'ope_small_s'):
        if f'JOINT_PHYSICAL_CHECK operator={name} passed=1 remaining=0' not in text:
            raise RuntimeError('Fresh native check missing: ' + name)
    # A resumed seam checkpoint is verified by hash, metadata and full sparse
    # load. Otherwise the just-computed full seam must have its own receipt.
    if not checkpoint_resumed and 'V18_SEAM_VERIFIED full_operator=1' not in text:
        raise RuntimeError('Full native seam verification missing')
    if not re.search(rf'JOINT_PRIME_COMPLETE prime={prime} directions=0 columns=275246\b', text):
        raise RuntimeError('Native completion receipt missing')
    return prime


def stop(signum, _frame):
    global cancel_signal
    cancel_signal = signum
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)


def preflight():
    profile = json.loads(PROFILE.read_text())
    if profile['selected'] != 'V27D' or profile['threads'] != 16:
        raise RuntimeError('Selected profile is not V27D/16')
    package, driver = Path(profile['package']), Path(profile['binary'])
    if sha(driver) != profile['binary_sha256']:
        raise RuntimeError('Selected solver changed')
    build = json.loads(Path(profile['build_receipt']).read_text())
    for name, digest in build['hashes'].items():
        if sha(package / name) != digest:
            raise RuntimeError('Frozen package changed: ' + name)
    suite = REPO / 'output_hexagon_v27_20261009/suite-d/progress.json'
    if json.loads(suite.read_text())['state'] != 'passed':
        raise RuntimeError('V27D controls must pass first')
    runtime = REPO / 'output_hexagon_runtime_recovery_20261008/validation-passed.json'
    for name, digest in json.loads(runtime.read_text())['binaries_and_libraries'].items():
        if sha(name) != digest:
            raise RuntimeError('Validated runtime changed: ' + name)
    cases = [line.split() for line in CASES.read_text().splitlines() if line.split()[:1] == ['10']]
    if len(cases) != 44 or any(len(c) != 3 for c in cases):
        raise RuntimeError('Expected 44 reserved origin cases')
    files = [Path(__file__), Path(__file__).with_name('hexagon_boot_recovery.py'), MEASURE,
             PROFILE, Path(profile['build_receipt']), suite, runtime, CONVERTER, CERTIFICATE,
             ORIGIN, CASES, REFERENCE / 'basis_identity.tsv', REFERENCE.parent / 'joint_inputs.txt']
    files += [REPO / 'bench' / name for name in ('hexagon_crt_six_prime_candidate_probe.py',
        'hexagon_crt_candidate_fresh_check.py', 'hexagon_crt_asymmetric_probe.py',
        'hexagon_crt_progress_audit.py', 'hexagon_v27d_qq_export.py', 'hexagon_candidate_export_controls.wls')]
    files += [REPO / 'output_hexagon_candidate_export_20261009' / name for name in
              ('candidate_qq_export', 'validation-passed.json')]
    recovery = EXPERIMENT / 'recovery-control/passed.json'
    proof = json.loads(recovery.read_text())
    if (proof['frozen_driver_sha256'] != sha(driver) or
            not proof['resumed_every_native_coefficient_equal'] or
            not proof['all_fresh_physical_checks_passed']):
        raise RuntimeError('Checkpoint restart control did not pass')
    files.append(recovery)
    identity = {'schema': 'hexagon-ten-v27d-reboot-v1', 'hashes': {str(p): sha(p) for p in files},
                'driver': str(driver), 'driver_sha256': sha(driver),
                'profile_environment': profile['environment'], 'unset': profile['unset_environment'],
                'native': str(NATIVE), 'shifted': str(SHIFTED), 'fec': str(FEC), 'lec': str(LEC),
                'data': str(DATA), 'h': 17, 'k': 3, 'target': 9000, 'prereduce': 'right',
                'memory_max_bytes': 50_000_000_000, 'memory_high_bytes': 48_000_000_000,
                'max_primes': 12, 'timeout_per_phase_seconds': 172800,
                'reconstruction': 'freeze six-prime asymmetric candidate before primes 7/8; two fresh checks, independent WXF decoder, exact seam and origin checks; balanced fallback if asymmetric recovery incomplete',
                'empirical_assumptions': ASSUMPTIONS, 'imports_old_coefficients': False,
                'qualification_baseline': 'V16; separate nine-loop qualification remains unverified'}
    return identity


def main():
    global active
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    identity = preflight()
    if args.preflight_only:
        print(json.dumps({'preflight_passed': True, 'configuration': identity}, indent=2))
        return 0
    cgroup_relative = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup') / cgroup_relative.lstrip('/')
    limit = (cgroup / 'memory.max').read_text().strip()
    if limit == 'max' or int(limit) > 50_000_000_000:
        raise RuntimeError('Start through the installed service with its 50 GB cgroup limit')
    for name in ('log', 'audit', 'interrupted'):
        (CONTROL / name).mkdir(parents=True, exist_ok=True)
    lock = (CONTROL / 'runner.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (CONTROL / 'attention-required.json').exists():
        raise RuntimeError('Previous failure requires review; see attention-required.json')
    manifest_path = CONTROL / 'manifest.json'
    if manifest_path.exists():
        if json.loads(manifest_path.read_text()) != identity:
            raise RuntimeError('Execution identity changed; separate reviewed experiment required')
    else:
        if NATIVE.exists() or SHIFTED.exists():
            raise RuntimeError('Fresh execution must not import an existing output chain')
        atomic_json(manifest_path, identity)
    ledger_path = CONTROL / 'ledger.json'
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else {'primes': {}, 'attempts': []}
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    env = {k: v for k, v in os.environ.items() if not k.startswith(('HEXAGON_', 'SOLVER_'))}
    env.update({k: str(v) for k, v in identity['profile_environment'].items()})
    for key in identity['unset']:
        env.pop(key, None)

    def commit():
        atomic_json(ledger_path, ledger)

    def progress(state, **fields):
        completed = [a['measurement'] for a in ledger['attempts'] if 'measurement' in a]
        cpu = sum(m['cpu_user_seconds'] + m['cpu_system_seconds'] for m in completed)
        partial = sum(a.get('last_cpu_seconds', 0) for a in ledger['attempts'] if 'measurement' not in a)
        value = {'schema': 'hexagon-ten-progress-v27d', 'state': state, 'boot_id': boot,
                 'runner_pid': os.getpid(), 'updated_unix_seconds': time.time(),
                 'native_primes_committed': len(ledger['primes']), 'max_primes': 12,
                 'accounted_cpu_hours': (cpu + partial) / 3600,
                 'interrupted_cpu_seconds_lower_bound': partial,
                 'final_symbol_certified': False, **fields}
        atomic_json(CONTROL / 'progress.json', value)
        print(json.dumps(value), flush=True)

    def recover():
        # Own each output lock during recovery, then release before the solver
        # takes the same locks. The outer controller lock remains held.
        locks = []
        try:
            for root in (NATIVE, SHIFTED):
                for path in sorted(root.rglob('.active.lock')):
                    handle = path.open('a')
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    locks.append(handle)
            retained = []
            for root in (NATIVE, SHIFTED):
                retained += recover_files(root, CONTROL / 'interrupted' / root.name)
            atomic_json(CONTROL / 'audit' / f'recovery-{uuid.uuid4().hex}.json', {
                'boot_id': boot, 'time': time.time(), 'retained_incomplete_files': retained})
        finally:
            for handle in locks:
                handle.close()
        for attempt in ledger['attempts']:
            if attempt['state'] == 'running':
                report = Path(attempt['report'])
                if report.exists():
                    try:
                        attempt['measurement'] = json.loads(report.read_text())
                    except json.JSONDecodeError:
                        pass
                attempt['state'] = 'interrupted_before_ledger_commit'
        commit()

    def successful(attempt):
        m = attempt.get('measurement', {})
        return m.get('exit_code') == 0 and not any(m.get(k) for k in (
            'cancel_signal', 'timeout', 'rss_limit_exceeded', 'memory_pressure'))

    def measured(phase, index, arguments, extra=None):
        global active
        if cancel_signal:
            raise InterruptedError('Shutdown requested')
        token = f'p{index:02d}-{phase}-{uuid.uuid4().hex[:12]}'
        report, log = CONTROL / 'audit' / (token + '.json'), CONTROL / 'log' / (token + '.log')
        status = CONTROL / 'audit' / (token + '.status.json')
        command = [sys.executable, str(MEASURE), '--report', str(report), '--log', str(log),
                   '--status', str(status), '--timeout', '172800', '--virtual-gib', '256',
                   '--rss-gib', str(50_000_000_000 / 1024**3), '--min-available-gib', '3', '--', *arguments]
        attempt = {'phase': phase, 'prime_index': index, 'command': command, 'report': str(report),
                   'log': str(log), 'status': str(status), 'boot_id': boot, 'started': time.time(),
                   'state': 'running', **(extra or {})}
        ledger['attempts'].append(attempt)
        commit()
        with (CONTROL / 'log' / (token + '.manager.log')).open('w') as console:
            active = subprocess.Popen(command, cwd=REPO, env=env, stdout=console, stderr=subprocess.STDOUT)
            progress('running_' + phase, prime_index=index, manager_pid=active.pid,
                     solver_log=str(log), resource_status=str(status))
            try:
                while active.poll() is None:
                    if cancel_signal:
                        active.send_signal(signal.SIGTERM)
                    try:
                        active.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        pass
                    attempt['last_cpu_seconds'] = max(attempt.get('last_cpu_seconds', 0), process_cpu(active.pid))
                    flush_completed((NATIVE, SHIFTED), CONTROL / 'durable-checkpoints.json')
                    commit()
                code = active.wait()
            finally:
                if active.poll() is None:
                    active.send_signal(signal.SIGTERM)
                    active.wait(timeout=20)
                active = None
        attempt['finished'] = time.time()
        attempt['manager_exit_code'] = code
        if report.exists():
            attempt['measurement'] = json.loads(report.read_text())
            sync(report)
        sync(log)
        if cancel_signal:
            attempt['state'] = 'interrupted'
            commit()
            raise InterruptedError('Shutdown requested; verified checkpoints retained')
        if code or not successful(attempt):
            attempt['state'] = 'failed'
            commit()
            raise RuntimeError(f'{phase} failed; review {report} and {log}; no unchanged retry')
        attempt['state'] = 'process_succeeded'
        commit()
        return attempt

    def validate_native(attempt, index):
        text = Path(attempt['log']).read_text()
        selected = re.findall(r'JOINT_SELECTED_PRIME index=\d+ prime=(\d+)', text)
        resumed = len(selected) == 1 and int(selected[0]) in attempt.get('seam_checkpoint_primes_at_start', [])
        prime = native_receipts(text, index, resumed)
        cache = check_cache(NATIVE / 'joint', prime, 17, 3, 275246)
        stage = tsv(NATIVE / 'joint/stage.tsv')
        expected = {'executable_sha256': identity['driver_sha256'], 'working_storage': 'heap',
                    'target': '9000', 'empirical_assumptions': ASSUMPTIONS,
                    'algorithm': 'joint-discovery-physical-affine-v27-workload-compression-packed-replay'}
        if any(stage.get(k) != v for k, v in expected.items()):
            raise RuntimeError('Native solver identity differs from V27D profile')
        if sha(NATIVE / 'joint_inputs.txt') != sha(REFERENCE.parent / 'joint_inputs.txt') or sha(NATIVE / 'joint/basis_identity.tsv') != sha(REFERENCE / 'basis_identity.tsv'):
            raise RuntimeError('Native input/basis identity differs from preserved reference')
        reference = REFERENCE / 'cache' / f'p{prime}.bin'
        comparison = {'available': reference.exists()}
        if reference.exists():
            check_cache(REFERENCE, prime, 17, 3, 275246)
            equal = Path(cache['path']).read_bytes()[HEADER.size:] == reference.read_bytes()[HEADER.size:]
            comparison.update({'every_coefficient_equal': equal, 'reference': str(reference), 'sha256': sha(reference)})
            if not equal:
                raise RuntimeError('Fresh V27D coefficients differ from preserved ten-loop reference')
        attempt['state'] = 'verified_native'
        return {'native': cache, 'attempt': ledger['attempts'].index(attempt),
                'producer_algorithm': 'V27D', 'reference_comparison': comparison}

    def asymmetric(index):
        """Reuse the independently tested reconstruction strategy, not old coefficients."""
        frozen = CONTROL / 'asymmetric-candidate.json'
        if index == 6 and not frozen.exists() and not (CONTROL / 'asymmetric-incomplete.json').exists():
            # Uncommitted candidate files from shutdown remain in the audit
            # directory; a new unique filename never overwrites them.
            draft = CONTROL / 'audit' / ('candidate-' + uuid.uuid4().hex + '.json')
            attempt = measured('asymmetric-freeze', 6, [sys.executable,
                str(REPO / 'bench/hexagon_crt_six_prime_candidate_probe.py'), str(CONTROL), str(draft)])
            candidate = json.loads(draft.read_text())
            if candidate['candidate_complete']:
                atomic_json(frozen, candidate)
                attempt['state'] = 'candidate_frozen_before_new_check_primes'
                commit()
            else:
                atomic_json(CONTROL / 'asymmetric-incomplete.json', candidate)
        if index == 8 and frozen.exists() and not EXTERNAL.exists():
            fresh = CONTROL / 'audit' / ('fresh-checks-' + uuid.uuid4().hex + '.json')
            measured('asymmetric-fresh-checks', 8, [sys.executable,
                str(REPO / 'bench/hexagon_crt_candidate_fresh_check.py'), str(CONTROL), str(frozen), str(fresh), '7', '8'])
            if not json.loads(fresh.read_text())['both_fresh_checks_pass']:
                raise RuntimeError('Frozen rational candidate fails its two fresh checks')
            measured('asymmetric-export', 8, [sys.executable, str(REPO / 'bench/hexagon_v27d_qq_export.py'),
                str(CONTROL), str(SHIFTED), str(EXTERNAL), str(frozen), str(fresh)])
        return EXTERNAL.exists() and check_qq(EXTERNAL)

    def certify():
        if not check_qq():
            raise RuntimeError('QQ object incomplete')
        qqroot = EXTERNAL if EXTERNAL.exists() else SHIFTED
        seam = measured('exact-target-seam', 0, [str(CERTIFICATE), '16', '4', str(qqroot), str(FEC), str(LEC), str(DATA)])
        rows = re.findall(r'^SHIFTED_JOINT_CERTIFICATE (.+)$', Path(seam['log']).read_text(), re.M)
        if len(rows) != 1:
            raise RuntimeError('Exact certificate missing')
        certificate = json.loads(rows[0])
        if (certificate.get('exact_QQ_original_target_seam') is not True or
            certificate.get('shift_source_hash_receipts_verified') is not True or
            [certificate.get(k) for k in ('source_h', 'source_k', 'target_h', 'target_k', 'directions', 'columns')] != [17, 3, 16, 4, 0, 518048]):
            raise RuntimeError('Exact certificate scope mismatch')
        origin = measured('exact-origin-holdouts', 0, [str(ORIGIN), '16', '4', str(qqroot), str(FEC), str(LEC), str(CASES)])
        text = Path(origin['log']).read_text()
        if re.findall(r'^ORIGIN_CHECK PASS words=(\d+)$', text, re.M) != ['44'] or len(re.findall(r'^ORIGIN_WORD loop=10 ', text, re.M)) != 44:
            raise RuntimeError('44 exact origin holdouts did not pass')
        final = {'schema': 'hexagon-ten-v27d-certified-v1', 'qq_root': str(qqroot),
                 'qq_receipt_sha256': sha(qqroot / 'qq_complete.tsv'), 'native_prime_producers': 'all V27D; no imported coefficients',
                 'empirical_assumptions': ASSUMPTIONS, 'unrestricted_uniqueness_claimed': False,
                 'exact_target_seam': certificate, 'exact_origin_holdouts': 44,
                 'physical_scope': 'full supplied COL/MRK/OPE checks each native prime; finite kinematics and OPE data',
                 'seam_log': seam['log'], 'origin_log': origin['log'], 'created': time.time()}
        atomic_json(CONTROL / 'final-verification.json', final)
        progress('complete_conditional_candidate_certified', final_symbol_certified=True, certificate=final)
        return 0

    try:
        recover()
        if (CONTROL / 'final-verification.json').exists():
            final = json.loads((CONTROL / 'final-verification.json').read_text())
            if not check_qq() or sha(Path(final['qq_root']) / 'qq_complete.tsv') != final['qq_receipt_sha256']:
                raise RuntimeError('Completed result changed')
            progress('complete_conditional_candidate_certified', final_symbol_certified=True, certificate=final)
            return 0
        if check_qq():
            return certify()
        for index in range(1, 13):
            if cancel_signal:
                raise InterruptedError('Shutdown requested')
            entry = ledger['primes'].get(str(index))
            if entry:
                if check_cache(NATIVE / 'joint', entry['native']['prime'], 17, 3, 275246) != entry['native']:
                    raise RuntimeError('Committed prime changed')
            else:
                prior = next((a for a in reversed(ledger['attempts']) if a['phase'] == 'native' and a['prime_index'] == index and successful(a)), None)
                if prior is None:
                    # A final cache without completed process accounting must
                    # not bypass fresh checks on recovery. Retain it, then let
                    # the solver reload the complete seam checkpoint instead.
                    committed = {Path(e['native']['path']).name for e in ledger['primes'].values()}
                    for path in (NATIVE / 'joint/cache').glob('*.bin'):
                        if path.name not in committed:
                            quarantine(path, NATIVE, CONTROL / 'interrupted/native', 'uncommitted final cache; rerun final checks')
                            side = path.with_name(path.name + '.sha256')
                            if side.exists():
                                quarantine(side, NATIVE, CONTROL / 'interrupted/native', 'uncommitted final cache checksum')
                    # Only a completed, hash-verified seam for the selected
                    # prime permits omission of a new seam-computation marker.
                    candidates = []
                    for p in (NATIVE / 'seam_complete/cache').glob('*.csr'):
                        if p.with_name(p.name + '.sha256').is_file():
                            candidates.append(int(p.stem[1:]))
                    prior = measured('native', index, [identity['driver'], '17', '3', str(FEC), str(LEC), str(DATA),
                        str(NATIVE), '9000', 'cache', 'right', '--prime-index', str(index)],
                        {'seam_checkpoint_primes_at_start': candidates})
                    selected = re.findall(r'JOINT_SELECTED_PRIME index=\d+ prime=(\d+)', Path(prior['log']).read_text())
                    prior['seam_checkpoint_resumed'] = len(selected) == 1 and int(selected[0]) in candidates
                entry = validate_native(prior, index)
                ledger['primes'][str(index)] = entry
                commit()
                progress('native_prime_committed', prime_index=index, prime=entry['native']['prime'])
            if entry.get('conversion'):
                if check_cache(SHIFTED, entry['native']['prime'], 16, 4, 518048) != entry['conversion']:
                    raise RuntimeError('Committed shifted prime changed')
                if asymmetric(index):
                    return certify()
                continue
            converted = measured('convert', index, [str(CONVERTER), 'shift-cache', 'MHV', '17', '3',
                str(NATIVE / 'joint'), str(SHIFTED), '--head', '16', '--fec-root', str(FEC), '--lec-root', str(LEC), '--kind', 'unique', '--reconstruct'])
            prime = entry['native']['prime']
            entry['conversion'] = check_cache(SHIFTED, prime, 16, 4, 518048)
            if (SHIFTED / f'source_p{prime}.sha256').read_text().strip() != entry['native']['sha256']:
                raise RuntimeError('Shift source hash mismatch')
            qq = check_qq(SHIFTED)
            if not qq and ('CACHE_ONLY completed_primes=' not in Path(converted['log']).read_text() or
                           json.loads((SHIFTED / 'manifest.json').read_text()).get('complete') is not False):
                raise RuntimeError('Neither cache-only nor QQ conversion complete')
            converted['state'] = 'verified_conversion'
            commit()
            if qq or asymmetric(index):
                return certify()
            progress('converted_prime_committed_cache_only', prime_index=index, prime=prime)
        progress('prime_budget_exhausted_cache_only')
        return 2
    except InterruptedError as error:
        progress('interrupted_resumable', error=str(error), cancel_signal=cancel_signal)
        return 0
    except Exception as error:
        atomic_json(CONTROL / 'attention-required.json', {'error': str(error), 'time': time.time(), 'boot_id': boot})
        progress('attention_required', error=str(error))
        return 78


if __name__ == '__main__':
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)
    try:
        sys.exit(main())
    except Exception as error:
        print('RUNNER_ERROR ' + str(error), file=sys.stderr, flush=True)
        sys.exit(78)
