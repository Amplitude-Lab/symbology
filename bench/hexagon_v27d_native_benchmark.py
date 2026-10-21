#!/usr/bin/env python3
"""Reconstruct and certify the fresh V27D native reference without new solves.

Run once in an owned cgroup with the existing 50 GB memory ceiling. Every
numerical stage is measured separately. Existing producers remain immutable.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

R = Path(__file__).resolve().parents[1]
F = R / 'output_hexagon_ten_v27d_20261009'
C = F / 'runner'
O = F / 'native-benchmark-20261010'
Q = O / 'qq-17p3'
P = R / 'output_hexagon_16p4_20261006'
V = R / 'output_hexagon_validation_20261009'
EXPORT = R / 'output_hexagon_candidate_export_20261009/candidate_qq_export'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1048576), b''):
            h.update(data)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temp.replace(path)


def cgroup_stats():
    return {k: int(v) for k, v in
            (line.split() for line in (CG / 'cpu.stat').read_text().splitlines())}


def measured(name, command, timeout=2400):
    state['stage'] = name
    save(O / 'progress.json', state)
    before = cgroup_stats()
    print('START', name, flush=True)
    command = list(map(str, command))
    report = O / (name + '.resources.json')
    log = O / (name + '.log')
    result = subprocess.run([
        'python3', str(R / 'bench/measure_process.py'),
        '--report', str(report), '--log', str(log),
        '--status', str(O / 'resource-status.json'), '--timeout', str(timeout),
        '--virtual-gib', '128', '--rss-gib', '46.566128730773926',
        '--min-available-gib', '3', '--', *command], cwd=R, env=env)
    after = cgroup_stats()
    stats = load(report)
    entry = {'name': name, 'report': str(report), 'report_sha256': sha(report),
             'log': str(log), 'log_sha256': sha(log), 'measurement': stats,
             'cgroup_cpu_seconds_including_monitor':
             (after['usage_usec'] - before['usage_usec']) / 1e6}
    state['stages'].append(entry)
    save(O / 'progress.json', state)
    assert result.returncode == stats['exit_code'] == 0, name
    assert not any(stats[k] for k in
                   ['timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal']), name
    print('PASS', name, flush=True)
    return log.read_text()


def main():
    assert not O.exists(), 'Preserve any previous attempt; use a reviewed new output path'
    assert subprocess.check_output([
        'systemctl', '--user', 'show', 'symbology-hexagon-ten-v27d-20261009.service',
        '-p', 'MainPID', '--value'], text=True).strip() == '0'
    assert load(C / 'progress.json')['final_symbol_certified']
    # cgroup files report page-rounded values, while systemd retains decimal
    # requested limits. Require a rounded-down cap without expanding it.
    page_size = os.sysconf('SC_PAGE_SIZE')
    for name, requested in [('memory.max', 50000000000), ('memory.high', 48000000000)]:
        actual = int((CG / name).read_text())
        assert requested - page_size < actual <= requested, (name, actual)
    assert (CG / 'memory.swap.max').read_text().strip() == '0'
    frozen = load(C / 'manifest.json')
    inputs = dict(frozen['hashes'])
    inputs[frozen['driver']] = frozen['driver_sha256']
    for path, digest in inputs.items():
        assert sha(path) == digest, path
    build = load(V / 'validation-build-binding.json')
    for path in [V / 'hexagon_joint_certificate', R / 'bench/hexagon_joint_certificate.cpp']:
        assert sha(path) == build['artifacts'][str(path)]['sha256']
    for path in [__file__, C / 'manifest.json', C / 'ledger.json',
                 C / 'audit/final-independent-review.json',
                 R / 'bench/hexagon_result_reconstruct.py',
                 V / 'hexagon_joint_certificate', V / 'hexagon_joint_certificate.build.json',
                 V / 'validation-build-binding.json', R / 'bench/hexagon_joint_certificate.cpp',
                 C / 'asymmetric-candidate.json']:
        inputs[str(path)] = sha(path)
    source = F / 'native-17p3/joint'
    for name in ['stage.tsv', 'basis_identity.tsv', 'input_identity.txt']:
        inputs[str(source / name)] = sha(source / name)
    candidate = load(C / 'asymmetric-candidate.json')
    assert sha(candidate['candidate_path']) == candidate['candidate_sha256']
    inputs[candidate['candidate_path']] = candidate['candidate_sha256']
    O.mkdir()
    save(O / 'plan.json', {
        'schema': 'v27d-native-reference-plan-v1', 'created': time.time(),
        'source_native_root': str(source), 'output': str(Q), 'input_hashes': inputs,
        'reconstruction': {'CRT_primes': 7, 'numerator_bound_bits': 96,
                           'check_prime_index': 8,
                           'check_scope': 'Already-computed prime excluded from CRT; retrospective holdout'},
        'solver': 'V27D/16', 'memory_high_bytes': 48000000000,
        'memory_max_bytes': 50000000000, 'swap_max_bytes': 0,
        'producer_modified': False,
        'qualification_baseline': 'V16; >15% CPU saving and <=5% RSS increase',
        'selection_scope': 'Primary ten-loop reference representation, not a new optimization qualification'})
    measured('reconstruction', ['python3', R / 'bench/hexagon_result_reconstruct.py',
        '--control', C, '--kind', 'native', '--output', O / 'reconstruction.json',
        '--crt-primes', 7, '--numerator-bits', 96, '--check-primes', 1])
    rec = load(O / 'reconstruction.json')
    assert rec['passed'] and rec['holdout_mismatches'] == [0]
    vector = Path(rec['vector'])
    assert sha(vector) == rec['vector_sha256']
    measured('export', [EXPORT, 17, 3, P / 'FEC', P / 'LEC/MHV', vector, Q])
    for name in ['input_identity.txt', 'basis_identity.tsv', 'stage.tsv']:
        shutil.copy2(source / name, Q / name)
    (Q / 'cache').mkdir()
    for cache in rec['cache_receipts']:
        path = Path(cache['path'])
        assert sha(path) == cache['sha256']
        for src in [path, Path(str(path) + '.sha256')]:
            inputs[str(src)] = sha(src)
            shutil.copy2(src, Q / 'cache' / src.name)
    (Q / 'reconstruction_primes.tsv').write_text(''.join(
        f"p{x['prime']}\t{x['role']}\n" for x in rec['cache_receipts']))
    (Q / 'summary.tsv').write_text(
        'crt_primes\t7\nretrospective_holdout_primes\t1\nmodulus_bits\t421\ndirections\t0\n')
    save(Q / 'manifest.json', {
        'schema': 'fresh-v27d-native-QQ-export-v1', 'created': time.time(),
        'source_native_root': str(source), 'reconstruction_receipt': str(O / 'reconstruction.json'),
        'reconstruction_receipt_sha256': sha(O / 'reconstruction.json'),
        'vector_sha256': sha(vector), 'input_hashes': inputs,
        'old_coefficients_imported': False, 'producer_modified': False,
        'physics_certificate': False,
        'assumption': frozen['empirical_assumptions']})
    measured('seal', [EXPORT, '--seal', Q])
    decoded = measured('independent-native-decode', ['wolframscript', '-file',
        R / 'bench/hexagon_candidate_export_controls.wls', 'verify', Q, vector])
    assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=275246 nnz=246941' in decoded
    seamlog = measured('exact-native-seam', [V / 'hexagon_joint_certificate',
        17, 3, Q, P / 'FEC', P / 'LEC/MHV', P / 'data'])
    seam = json.loads(next(line[len('JOINT_CERTIFICATE '):]
        for line in seamlog.splitlines() if line.startswith('JOINT_CERTIFICATE ')))
    assert seam['exact_QQ_original_seam'] and seam['directions'] == 0
    assert seam['original_rows'] == 2053259 and seam['columns'] == 275246
    assert seam['modulus_bits'] > seam['bound_bits']
    originlog = measured('exact-native-origin', [P / 'joint-v8-package/bin/origin_check',
        17, 3, Q, P / 'FEC', P / 'LEC/MHV',
        R / 'audits/hexagon-ten-loop-2026-10-06/origin-cases.tsv'])
    assert 'ORIGIN_CHECK PASS words=44' in originlog
    shifted = O / 'exact-shift-16p4'
    measured('exact-shift', [P / 'joint-v8-package/bin/irr_convert', 'shift', 'MHV',
        17, 3, Q, shifted, '--head', 16, '--fec-root', P / 'FEC', '--lec-root', P / 'LEC/MHV'])
    equality = measured('independent-split-equality', ['wolframscript', '-file',
        R / 'bench/hexagon_candidate_export_controls.wls', 'verify', shifted, candidate['candidate_path']])
    assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=518048 nnz=491197' in equality
    historical = {}
    for file in Q.glob('*.wxf'):
        assert sha(file) == sha(V / 'native-qq-17p3' / file.name)
        historical[file.name] = sha(file)
    assert len(historical) == 6
    for path, digest in inputs.items():
        assert sha(path) == digest, 'Input changed: ' + path
    state.update(state='passed', finished=time.time(), exact_native_seam=seam,
        independent_WXF_decode=True, exact_origin_holdouts=44,
        exact_equality_with_fresh_certified_16p4=True,
        exact_historical_native_WXF_hashes=historical,
        reconstruction=rec, qq_root=str(Q), qq_seal_sha256=sha(Q / 'qq_complete.tsv'),
        producer_modified=False, input_hashes=inputs,
        assumptions=frozen['empirical_assumptions'], unrestricted_uniqueness_claimed=False,
        cgroup_memory_peak_bytes=int((CG / 'memory.peak').read_text()),
        cgroup_memory_events=(CG / 'memory.events').read_text())
    save(O / 'verification.json', state)
    save(O / 'progress.json', state)
    print('NATIVE_BENCHMARK_VALIDATION_PASS', flush=True)


if __name__ == '__main__':
    CG = Path('/sys/fs/cgroup') / next(line.split(':', 2)[2].lstrip('/')
        for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    env = os.environ.copy()
    env.update(THREADS='16', LD_LIBRARY_PATH=str(R / 'output_hexagon_runtime_recovery_20261008/env/lib'))
    state = {'schema': 'v27d-native-benchmark-validation-v1', 'state': 'running',
             'started': time.time(), 'stages': []}
    try:
        main()
    except BaseException as error:
        state.update(state='failed', error=repr(error))
        if O.exists():
            save(O / 'attention-required.json', state)
            save(O / 'progress.json', state)
        raise
