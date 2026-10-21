#!/usr/bin/env python3
"""Audit and select the completed native V27D reference, preserving history."""
from array import array
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import time

from hexagon_v27d_native_benchmark import R, F, C, O, Q, sha, load, save


def main():
    selection = R / 'benchmarks/hexagon-mhv-ten-loop/current-benchmark.json'
    assert not selection.exists(), 'Preserve an existing selection before a reviewed update'
    assert not (O / 'benchmark.json').exists()
    verified = load(O / 'verification.json')
    assert verified['state'] == 'passed'
    assert not (O / 'attention-required.json').exists()
    unit = 'symbology-hexagon-v27d-native-qq-20261010-r1.service'
    properties = subprocess.check_output(['systemctl', '--user', 'show', unit,
        '-p', 'LoadState', '-p', 'ActiveState', '-p', 'SubState', '-p', 'Result', '-p', 'MainPID',
        '-p', 'ExecMainStatus', '-p', 'CPUUsageNSec', '-p', 'MemoryHigh',
        '-p', 'MemoryMax', '-p', 'MemorySwapMax'], text=True)
    service = dict(line.split('=', 1) for line in properties.splitlines())
    assert service['ActiveState'] == 'inactive' and service['Result'] == 'success'
    assert service['MainPID'] == '0' and service['ExecMainStatus'] == '0'
    plan = load(O / 'plan.json')
    assert plan['memory_max_bytes'] == 50000000000 and plan['swap_max_bytes'] == 0
    if service['LoadState'] != 'not-found':
        assert service['MemoryMax'] == '50000000000'
    service['run_memory_cap_source'] = ('Frozen plan and successful runtime cgroup preflight assertions; '
        'a collected transient unit reports default limits after exit, not its former runtime limits.')
    # Some systemd versions clear the live counter when the cgroup vanishes;
    # its structured completion journal retains the exact nanosecond total.
    journal = [json.loads(line) for line in subprocess.check_output([
        'journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    completions = [entry for entry in journal if entry.get('USER_UNIT') == unit
                   and 'CPU_USAGE_NSEC' in entry]
    assert len(completions) == 1
    service['completed_CPUUsageNSec'] = int(completions[0]['CPU_USAGE_NSEC'])
    service['CPU_source'] = 'structured systemd completion journal CPU_USAGE_NSEC'
    save(O / 'service-journal.json', journal)
    save(O / 'service-completion.json', service)
    for path, digest in verified['input_hashes'].items():
        assert sha(path) == digest, path
    assert sha(Q / 'qq_complete.tsv') == verified['qq_seal_sha256']
    historical_root = R / 'output_hexagon_validation_20261009/native-qq-17p3'
    assert len(verified['exact_historical_native_WXF_hashes']) == 6
    for name, digest in verified['exact_historical_native_WXF_hashes'].items():
        assert sha(Q / name) == digest == sha(historical_root / name)
    for line in (Q / 'qq_complete.tsv').read_text().splitlines():
        name, digest = line.split('\t')
        if name == 'schema':
            assert digest == 'irr-physics-qq-v1'
        else:
            file = Q / ('input_identity.txt' if name == 'identity' else name)
            assert sha(file) == digest, name
    for line in (Q / 'cache_files.tsv').read_text().splitlines():
        name, digest = line.split('\t')
        assert sha(Q / 'cache' / name) == digest
    rec = verified['reconstruction']
    assert rec == load(O / 'reconstruction.json') and rec['passed']
    assert sha(rec['vector']) == rec['vector_sha256']
    ledger = load(C / 'ledger.json')
    residues = []
    for index, cache in enumerate(rec['cache_receipts'], 1):
        source = ledger['primes'][str(index)]
        assert source['producer_algorithm'] == 'V27D'
        assert source['native']['path'] == cache['path']
        assert sha(cache['path']) == cache['sha256'] == source['native']['sha256']
        raw = Path(cache['path']).read_bytes()
        assert raw == (Q / 'cache' / Path(cache['path']).name).read_bytes()
        header = struct.unpack_from('<8s8Q', raw)
        prime = cache['prime']
        assert header == (b'NMHVSEW\0', 1, 10, 17, 3, 1, 275246, prime,
                          int(sha(Q / 'input_identity.txt')[:16], 16))
        values = array('Q')
        values.frombytes(raw[72:])
        if sys.byteorder != 'little':
            values.byteswap()
        assert len(values) == 275246 and all(x < prime for x in values)
        residues.append((prime, values))
    rows = nnz = numerator_bits = denominator_bits = 0
    for line in Path(rec['vector']).read_text().splitlines():
        if line.startswith('#'):
            continue
        i, n, d = map(int, line.split())
        assert i == rows and d > 0 and math.gcd(n, d) == 1
        assert all(d % prime and (n - d * values[i]) % prime == 0
                   for prime, values in residues)
        rows += 1
        nnz += n != 0
        numerator_bits = max(numerator_bits, abs(n).bit_length())
        denominator_bits = max(denominator_bits, d.bit_length())
    assert (rows, nnz, numerator_bits, denominator_bits) == (275246, 246941, 45, 315)
    stages = {}
    for stage in verified['stages']:
        report = load(stage['report'])
        assert report == stage['measurement'] and report['exit_code'] == 0
        assert sha(stage['report']) == stage['report_sha256']
        assert sha(stage['log']) == stage['log_sha256']
        stages[stage['name']] = {
            'reported_process_CPU_hours': (report['cpu_user_seconds'] + report['cpu_system_seconds']) / 3600,
            'inclusive_cgroup_CPU_hours': stage['cgroup_cpu_seconds_including_monitor'] / 3600,
            'wall_seconds': report['wall_seconds'],
            'peak_sampled_RSS_bytes': report['peak_summed_rss_bytes'],
            'report': stage['report'], 'report_sha256': stage['report_sha256']}
    seam_lines = (O / 'exact-native-seam.log').read_text().splitlines()
    seam = json.loads(next(x[len('JOINT_CERTIFICATE '):] for x in seam_lines
                           if x.startswith('JOINT_CERTIFICATE ')))
    assert seam == verified['exact_native_seam'] and seam['exact_QQ_original_seam']
    assert seam['original_rows'] == 2053259 and seam['directions'] == 0
    assert seam['native_congruence_primes'] == 8 and seam['extra_residual_primes'] == 4
    assert seam['modulus_bits'] > seam['bound_bits']
    assert 'ORIGIN_CHECK PASS words=44' in (O / 'exact-native-origin.log').read_text()
    assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=275246 nnz=246941' in (O / 'independent-native-decode.log').read_text()
    assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=518048 nnz=491197' in (O / 'independent-split-equality.log').read_text()
    previous = load(C / 'audit/final-independent-review.json')
    assert previous['passed']
    native_stages = [a['measurement'] for a in ledger['attempts'] if a['phase'] == 'native']
    assert len(native_stages) == 8
    native_cpu = sum(m['cpu_user_seconds'] + m['cpu_system_seconds'] for m in native_stages) / 3600
    assert math.isclose(native_cpu, previous['phase_totals']['native']['cpu_hours'], rel_tol=1e-12)
    frozen = load(C / 'manifest.json')
    audit = {
        'schema': 'hexagon-v27d-native-benchmark-selection-v1', 'created': time.time(),
        'status': 'selected_verified_working_benchmark', 'version': 'V27D', 'workers': 16,
        'primary_split': [17, 3], 'primary_QQ_root': str(Q),
        'primary_vector': rec['vector'], 'primary_vector_sha256': rec['vector_sha256'],
        'primary_QQ_seal_sha256': sha(Q / 'qq_complete.tsv'),
        'secondary_split': [16, 4], 'secondary_QQ_root': previous['qq_root'],
        'solver': frozen['driver'], 'solver_sha256': frozen['driver_sha256'],
        'profile_environment': frozen['profile_environment'], 'unset': frozen['unset'],
        'frozen_run_manifest': str(C / 'manifest.json'),
        'frozen_run_manifest_sha256': sha(C / 'manifest.json'),
        'columns': rows, 'nonzero': nnz, 'max_numerator_bits': numerator_bits,
        'max_denominator_bits': denominator_bits,
        'fresh_V27D_prime_solves': 8, 'CRT_primes': 7, 'retrospective_check_primes': 1,
        'holdout_scope': 'Prime8 was already computed; excluded from native CRT. No claim that it was newly generated after this reconstruction.',
        'old_coefficients_imported': False,
        'validation': {'original_native_seam': seam, 'exact_origin_holdouts': 44,
            'independent_native_WXF_decode': True, 'exact_equality_with_certified_fresh_16p4': True,
            'every_native_rational_congruence_rechecked': rows * 8,
            'historical_native_WXF_hashes': verified['exact_historical_native_WXF_hashes'],
            'receipt': str(O / 'verification.json'), 'receipt_sha256': sha(O / 'verification.json'),
            'previous_fresh_campaign_review': str(C / 'audit/final-independent-review.json'),
            'previous_fresh_campaign_review_sha256': sha(C / 'audit/final-independent-review.json')},
        'performance': {'eight_native_solves_process_CPU_hours': native_cpu,
            'native_solve_sampled_peak_RSS_bytes': previous['peak_sampled_process_RSS_bytes'],
            'native_solve_kernel_peak_RSS_bytes': previous['native_kernel_peak_RSS_bytes'],
            'additional_reconstruction_and_validation_stages': stages,
            'additional_service_CPU_hours': service['completed_CPUUsageNSec'] / 3.6e12,
            'additional_service_cgroup_memory_peak_bytes': verified['cgroup_memory_peak_bytes'],
            'additional_sampled_peak_RSS_bytes': max(x['peak_sampled_RSS_bytes'] for x in stages.values()),
            'memory_limit_bytes': 50000000000,
            'scope': 'Eight sequential fresh native solves with prepared bases, first learning and seven replay primes. Native reconstruction/validation is measured afterward. Exact cross-split conversion is an additional acceptance audit. No new cold end-to-end timing or causal speedup experiment.',
            'CPU_accounting': 'Process user+system totals are reported per stage. Cgroup CPU includes measurement overhead and Wolfram descendants not fully reaped by the wrapper. Do not silently combine these accounting scopes.'},
        'future_comparison_protocol': {
            'metric': 'total process user+system CPU-hours, including all descendants; also record inclusive service CPU and peak RSS',
            'workload': 'Same native17+3 bases, equations, RHS, eight primes and exact validation scope; learning and replay costs separate',
            'representation_changes': 'Measure all reconstruction and conversion costs if choosing another internal split; compare exported results to this exact native17+3 reference',
            'qualification_baseline': 'V16', 'CPU_saving_required_strictly_greater_than': 0.15,
            'RSS_ratio_max': 1.05, 'nine_loop_V16_qualification': 'unverified',
            'ten_loop_first_prime_V16_comparison': 'CPU improvement observed; RSS tolerance failed',
            'working_reference_selection_is_not_new_qualification': True},
        'assumptions': frozen['empirical_assumptions'], 'unrestricted_uniqueness_claimed': False,
        'frozen_producer_modified': False, 'historical_benchmark_manifest_preserved': True,
        'audit_source_sha256': sha(__file__)}
    save(O / 'benchmark.json', audit)
    save(selection, {'schema': 'hexagon-current-ten-loop-benchmark-v1',
        'selected': 'V27D/16, native17+3', 'primary_QQ_root': str(Q),
        'benchmark': str(O / 'benchmark.json'), 'benchmark_sha256': sha(O / 'benchmark.json'),
        'qualification_baseline': 'V16', 'previous_manifest_is_historical': True})
    lessons_path = R / 'docs/hexagon_optimization_lessons.json'
    lessons = load(lessons_path)
    save(O / 'lessons-before-selection.json', lessons)
    lessons['updated'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    lessons['scope'] = ('Fresh V27D/16 ten-loop computation and exact native17+3 and shifted16+4 reconstructions are complete. '
        'Native17+3 is the selected primary benchmark; both representations and all historical receipts are preserved. '
        'Fixed V16 >15% CPU / <=5% RSS qualification is unchanged: direct nine-loop comparison remains unverified, '
        'and the historical ten-loop first-prime comparison fails the RSS tolerance. No universal optimality claim.')
    lessons['current_ten_loop_benchmark'] = str(selection)
    lessons['lessons'].append({
        'id': 'fresh-v27d-native17p3-primary-benchmark-20261010',
        'finding': 'Direct reconstruction from seven fresh V27D native prime caches gives all275246 rational coordinates,246941 nonzero,max45 numerator bits and315 denominator bits. Prime8 excluded from CRT agrees. Independent WXF decode,2053259 exact original native sewing rows,44 origin checks and exact conversion equality with the certified fresh16+4 result pass. Native17+3 is selected as primary benchmark by user request.',
        'practice': 'Preserve native and shifted references. Choose any reconstruction split by measured complete cost and coefficient heights. Native7+1 uses an already-computed retrospective check; do not relabel it as a new postfreeze prime. Working benchmark selection does not reset the V16 qualification baseline.',
        'audit': str(O / 'benchmark.json'), 'audit_sha256': sha(O / 'benchmark.json'),
        'performance': audit['performance'], 'producer_modified': False,
        'preflight_correction': 'The first wrapper launch failed before any numerical work because an exact decimal cgroup-limit assertion ignored page rounding. Accept only rounded-down values within one page; caps remain48GB/50GB and swap0. Failure receipt is preserved under runner/audit/native-benchmark-preflight-20261010.json.',
        'resource_accounting_lesson': 'Collected transient systemd units return unavailable CPU and default unlimited memory properties after completion. Persist runtime counters or use the structured completion journal CPU_USAGE_NSEC; do not misread post-collection defaults as the limits used by the run. A postprocessing assertion was corrected accordingly; no numerical work was rerun.'})
    save(lessons_path, lessons)
    print(json.dumps({'selected': str(selection), 'performance': audit['performance'],
                      'validation': audit['validation']}, indent=2))


if __name__ == '__main__':
    main()
