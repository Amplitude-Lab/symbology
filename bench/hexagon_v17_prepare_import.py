#!/usr/bin/env python3
"""Create a new, auditable V17 chain without changing the frozen V16 result."""
import fcntl
import importlib.util
import json
from pathlib import Path
import sys
import time
import hexagon_v17_cache_migration as m

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
OLD = BASE / 'ten-v16-runner'
CONTROL = BASE / 'ten-v17-runner'
NATIVE = BASE / 'ten-v17-main/joint'
SHIFTED = BASE / 'ten-v17-shifted-16p4'
OLD_NATIVE = BASE / 'ten-v16-main/joint'
OLD_SHIFTED = BASE / 'ten-v16-shifted-16p4'
PACKAGE = BASE / 'joint-v17-schedule-package'
DRIVER = PACKAGE / 'bin/joint_physics_v17'
PRIME = 1152921504606847009


def main():
    CONTROL.mkdir(exist_ok=True)
    locks = [(p / 'runner.lock').open('a') for p in [OLD, CONTROL]]
    for lock in locks:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    m.require(not (CONTROL / 'v16-prime-import.json').exists(), 'Import already committed; verify and reuse it')
    manifest = json.loads((OLD / 'manifest.json').read_text())
    fixed = {ROOT / 'bench/hexagon_ten_loop_runner_v16.py': manifest['runner_sha256'],
             ROOT / 'bench/measure_process.py': manifest['measure_sha256'],
             Path(manifest['driver']): manifest['driver_sha256'],
             Path(manifest['converter']): manifest['converter_sha256']}
    for path, digest in fixed.items():
        m.require(m.sha(path) == digest, 'Frozen V16 executable changed: ' + str(path))
    m.require(json.loads((OLD / 'progress.json').read_text())['state'] == 'prime_budget_exhausted_cache_only',
              'V16 pilot is not in its expected terminal state')
    eight = BASE / 'v17-schedule-controls/progress.json'
    nine = BASE / 'v17-nine-control/progress.json'
    controls = BASE / 'v17-migration-controls/progress.json'
    e, n, c = [json.loads(p.read_text()) for p in [eight, nine, controls]]
    m.require(e['status'] == 'cursor_and_eight_controls_passed' and e['every_eight_loop_coefficient_equal'], 'Eight-loop validation missing')
    m.require(n['status'] == 'nine_loop_regression_passed' and n['every_native_coefficient_equal'] and n['every_shifted_coefficient_equal'], 'Nine-loop validation missing')
    m.require(c['status'] == 'migration_controls_passed' and c['migration_module_sha256'] == m.sha(Path(m.__file__)) and
              c['driver_sha256'] == m.sha(DRIVER) == n['identity']['v17_driver_sha256'], 'Migration validation identity mismatch')
    old_package = BASE / 'joint-v16-mrk-reuse-package'
    differences = []
    for path in sorted((old_package / 'src').rglob('*')):
        if path.is_file():
            relative = path.relative_to(old_package)
            if m.sha(path) != m.sha(PACKAGE / relative):
                differences.append(str(relative))
    m.require(differences == ['src/shared/adaptive_physical_rows.hpp'], 'Unexpected equation-generation source change')
    m.require(m.sha(PACKAGE / 'src/shared/adaptive_physical_rows.hpp') == m.sha(ROOT / 'bench/hexagon_adaptive_physical_rows_v17.hpp'), 'Scheduler source mismatch')
    m.require(m.sha(PACKAGE / 'src/shared/hexagon_physics_work_cursor.hpp') == m.sha(ROOT / 'bench/hexagon_physics_work_cursor.hpp'), 'Cursor source mismatch')
    old_driver = (old_package / 'tests/joint_physics_v16.cpp').read_text()
    old_schedule = m.fields(OLD_NATIVE / 'stage.tsv')['early_physics']
    expected_driver = old_driver.replace(m.OLD_ALGORITHM, m.NEW_ALGORITHM).replace(old_schedule, m.NEW_SCHEDULE)
    m.require(expected_driver == (PACKAGE / 'tests/joint_physics_v17.cpp').read_text(), 'Unexpected driver changes')
    spec = importlib.util.spec_from_file_location('original_runner', ROOT / 'bench/hexagon_ten_loop_runner_v16.py')
    original = importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    ledger = json.loads((OLD / 'ledger.json').read_text())
    m.require(set(ledger['primes']) == {'1'}, 'Expected exactly the one completed V16 prime')
    entry = ledger['primes']['1']
    m.require(original.check_cache(OLD_NATIVE, PRIME, 17, 3, 275246) == entry['native'], 'Original native receipt mismatch')
    m.require(original.check_cache(OLD_SHIFTED, PRIME, 16, 4, 518048) == entry['conversion'], 'Original shifted receipt mismatch')
    attempt = ledger['attempts'][entry['native_attempt']]
    report = json.loads(Path(attempt['report']).read_text()); log = Path(attempt['log']).read_text()
    m.require(report['exit_code'] == 0 and report['peak_summed_rss_bytes'] < 50_000_000_000 and not any(report.get(k) for k in
        ['timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal']), 'Original execution failed or exceeded budget')
    for condition in ['collinear', 'mrk', 'ope_small_s']:
        m.require(f'JOINT_PHYSICAL_CHECK operator={condition} passed=1 remaining=0' in log, 'Original physical check missing')
    m.require('JOINT_SEAM_COMPLETE full_operator=1' in log and f'JOINT_PRIME_COMPLETE prime={PRIME} directions=0 columns=275246' in log,
              'Original seam or native completion missing')
    native = m.import_native(OLD_NATIVE, NATIVE, DRIVER, PRIME)
    shifted = m.import_shifted(OLD_SHIFTED, SHIFTED, OLD_NATIVE, NATIVE, PRIME)
    pinned_paths = [*fixed.keys(), OLD / 'manifest.json', OLD / 'ledger.json', Path(attempt['report']), Path(attempt['log']),
        eight, nine, controls, Path(m.__file__), Path(__file__), DRIVER, BASE / 'v17-schedule-controls/source-diff.json',
        OLD_NATIVE.parent / 'joint_inputs.txt', NATIVE.parent / 'joint_inputs.txt']
    for directory in [OLD_NATIVE, OLD_SHIFTED, NATIVE, SHIFTED]:
        pinned_paths.extend(directory / name for name in ['input_identity.txt', 'basis_identity.tsv', 'stage.tsv'])
    pinned_paths.extend([OLD_SHIFTED / f'source_p{PRIME}.sha256', SHIFTED / f'source_p{PRIME}.sha256'])
    receipt = {'schema': 'hexagon-explicit-cache-migration-v1', 'created': time.time(), 'prime': PRIME,
        'original_producer_algorithm': m.OLD_ALGORITHM, 'destination_pipeline_algorithm': m.NEW_ALGORITHM,
        'native_prime_recomputed': False, 'shift_map_recomputed': False, 'caches': [native, shifted],
        'original_runner_manifest': str(OLD / 'manifest.json'), 'original_ledger': str(OLD / 'ledger.json'),
        'changed_descriptor_keys': sorted(m.ALLOWED), 'unchanged': 'all coefficient payloads, geometric/physical inputs, D3 conventions, affine origin, basis identities',
        'meaning': 'The stage descriptor selects V17 for future native evaluations. Prime 1 was produced by V16 and is explicitly imported, not relabeled as a V17 calculation.',
        'validation': {'scheduler_and_eight': str(eight), 'nine': str(nine), 'migration': str(controls)},
        'pinned_files': {str(p): m.sha(m.ordinary(p)) for p in pinned_paths}}
    m.install(CONTROL / 'v16-prime-import.json', json.dumps(receipt, indent=2) + '\n')
    m.verify_receipt(CONTROL / 'v16-prime-import.json')
    new_entry = {'native': original.check_cache(NATIVE, PRIME, 17, 3, 275246),
        'conversion': original.check_cache(SHIFTED, PRIME, 16, 4, 518048),
        'native_identity_sha256': m.sha(NATIVE / 'input_identity.txt'),
        'native_basis_identity_sha256': m.sha(NATIVE / 'basis_identity.tsv'),
        'shifted_identity_sha256': m.sha(SHIFTED / 'input_identity.txt'),
        'producer_algorithm': 'V16; explicit validated import', 'import_receipt': str(CONTROL / 'v16-prime-import.json'),
        'original_native_attempt': attempt}
    m.install(CONTROL / 'ledger.json', json.dumps({'primes': {'1': new_entry}, 'attempts': []}, indent=2) + '\n')
    print(json.dumps({'import_verified': True, 'native_recomputed': False, 'shifted_recomputed': False,
                      'native_payload_sha256': native['payload_sha256'], 'shifted_payload_sha256': shifted['payload_sha256']}))


if __name__ == '__main__':
    main()
