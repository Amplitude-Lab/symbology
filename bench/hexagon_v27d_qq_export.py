#!/usr/bin/env python3
"""Transactional export of a newly frozen, six-plus-two checked V27D vector."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid
from hexagon_boot_recovery import atomic_json, sha, sync

REPO = Path(__file__).resolve().parents[1]


def main():
    control, shifted, target, frozen, fresh = map(Path, sys.argv[1:])
    if target.exists():
        raise RuntimeError('Refuse to replace a rational output')
    candidate, checks = json.loads(frozen.read_text()), json.loads(fresh.read_text())
    vector = Path(candidate['candidate_path'])
    assert candidate['candidate_complete'] and candidate['prime7_not_committed_at_freeze']
    assert candidate['candidate_rows'] == 518048 and sha(vector) == candidate['candidate_sha256']
    assert checks['both_fresh_checks_pass'] and checks['candidate_sha256'] == sha(vector)
    assert [c['prime_index'] for c in checks['checks']] == [7, 8]
    export_root = REPO / 'output_hexagon_candidate_export_20261009'
    validation_path = export_root / 'validation-passed.json'
    validation = json.loads(validation_path.read_text())
    exporter = export_root / 'candidate_qq_export'
    assert validation['passed'] and sha(exporter) == validation['binary_sha256']
    assert sha(REPO / 'bench/hexagon_candidate_qq_export.cpp') == validation['source_sha256']
    for row in validation['known_result_controls']:
        for name, digest in row['files'].items():
            assert sha(Path(row['reference']) / name) == digest == sha(Path(row['export']) / name)
    ledger = json.loads((control / 'ledger.json').read_text())
    pinned = {str(p): sha(p) for p in (frozen, fresh, vector, validation_path, exporter, Path(__file__))}
    primes, modulus = [], 1
    for index in range(1, 9):
        entry = ledger['primes'][str(index)]
        assert entry['producer_algorithm'] == 'V27D'
        for key in ('native', 'conversion'):
            cache = Path(entry[key]['path'])
            assert sha(cache) == entry[key]['sha256'] == Path(str(cache) + '.sha256').read_text().strip()
            pinned[str(cache)] = sha(cache)
        prime = entry['native']['prime']
        assert prime == entry['conversion']['prime']
        assert (shifted / f'source_p{prime}.sha256').read_text().strip() == entry['native']['sha256']
        if index <= 6:
            assert candidate['crt_primes'][index - 1] == prime
            assert candidate['cache_receipts'][index - 1]['sha256'] == entry['conversion']['sha256']
            modulus *= prime
        else:
            check = checks['checks'][index - 7]
            assert check['prime'] == prime and check['sha256'] == entry['conversion']['sha256']
            assert check['mismatches'] == check['noninvertible'] == 0
        primes.append(entry)
    temporary = target.with_name(target.name + '.incomplete-' + uuid.uuid4().hex)
    shared = REPO / 'output_hexagon_16p4_20261006'
    subprocess.run([str(exporter), '16', '4', str(shared / 'FEC'), str(shared / 'LEC/MHV'),
                    str(vector), str(temporary)], check=True)
    for name in ('input_identity.txt', 'basis_identity.tsv', 'stage.tsv'):
        shutil.copy2(shifted / name, temporary / name)
    (temporary / 'cache').mkdir()
    roles = {}
    for index, entry in enumerate(primes, 1):
        cache = Path(entry['conversion']['path'])
        for path in (cache, Path(str(cache) + '.sha256')):
            shutil.copy2(path, temporary / 'cache' / path.name)
        prime = entry['native']['prime']
        name = f'source_p{prime}.sha256'
        shutil.copy2(shifted / name, temporary / name)
        roles[f'p{prime}'] = 'CRT' if index <= 6 else 'validation'
    (temporary / 'reconstruction_primes.tsv').write_text(''.join(k + '\t' + v + '\n' for k, v in sorted(roles.items())))
    (temporary / 'summary.tsv').write_text(f'crt_primes\t6\nvalidation_primes\t2\nmodulus_bits\t{modulus.bit_length()}\ndirections\t0\n')
    atomic_json(temporary / 'manifest.json', {
        'schema': 'hexagon-external-asymmetric-QQ-export-v1', 'complete': True,
        'crt_primes': 6, 'validation_primes': 2, 'candidate_sha256': sha(vector),
        'source_shifted_root': str(shifted), 'source_native_root': str(Path(primes[0]['native']['path']).parent.parent),
        'candidate_receipt': str(frozen), 'fresh_check_receipt': str(fresh),
        'exporter_binary_sha256': sha(exporter), 'input_hashes': pinned,
        'all_native_producers': 'V27D', 'physics_certified': False,
        'unrestricted_uniqueness_claimed': False, 'created': time.time()})
    subprocess.run([str(exporter), '--seal', str(temporary)], check=True)
    # Independent decoder verifies every rational in the six written WXF blocks.
    subprocess.run(['/usr/bin/wolframscript', '-file', str(REPO / 'bench/hexagon_candidate_export_controls.wls'),
                    'verify', str(temporary), str(vector)], check=True)
    for path, digest in pinned.items():
        assert sha(path) == digest, 'Input changed during export: ' + path
    for path in temporary.rglob('*'):
        if path.is_file():
            sync(path)
    sync(temporary / 'cache')
    sync(temporary)
    temporary.replace(target)
    sync(target.parent)
    print('V27D_QQ_EXPORT_COMMITTED columns=518048 crt_primes=6 validation_primes=2', flush=True)


if __name__ == '__main__':
    main()
