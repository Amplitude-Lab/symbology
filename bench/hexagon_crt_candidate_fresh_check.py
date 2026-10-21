#!/usr/bin/env python3
"""Check a frozen speculative vector against subsequently committed primes.

This writes only a diagnostic receipt, never a production cache or QQ export.
"""
from array import array
import json
import math
import os
from pathlib import Path
import resource
import struct
import sys
import time

from hexagon_crt_progress_audit import sha


def main():
    control, candidate_receipt, output = map(Path, sys.argv[1:4])
    indices = list(map(int, sys.argv[4:]))
    assert indices and len(indices) == len(set(indices))
    os.nice(10)
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    started = time.monotonic()
    candidate = json.loads(candidate_receipt.read_text())
    assert candidate['schema'] == 'hexagon-six-crt-speculative-candidate-v1'
    assert candidate['candidate_complete'] and candidate['prime7_not_committed_at_freeze']
    assert candidate['candidate_rows'] == candidate['columns'] == 518048
    assert not candidate['failed_columns']
    assert sha(Path(__file__).with_name('hexagon_crt_six_prime_candidate_probe.py')) == candidate['source_sha256']
    for name, digest in candidate['helper_sha256'].items():
        assert Path(name).name == name
        assert sha(Path(__file__).with_name(name)) == digest
    for entry in candidate['cache_receipts']:
        assert sha(Path(entry['path'])) == entry['sha256']
    vector = Path(candidate['candidate_path'])
    assert vector.is_file() and not vector.is_symlink()
    assert sha(vector) == candidate['candidate_sha256']
    assert vector.stat().st_size == candidate['candidate_bytes']
    ledger = json.loads((control / 'ledger.json').read_text())
    checks = []
    arrays = []
    for index in indices:
        assert index in candidate['future_fresh_check_prime_indices']
        entry = ledger['primes'][str(index)]['conversion']
        path = Path(entry['path'])
        assert path.is_file() and not path.is_symlink()
        assert path.stat().st_mtime > candidate['recorded_unix_seconds']
        assert sha(path) == entry['sha256'] == path.with_suffix('.bin.sha256').read_text().strip()
        with path.open('rb') as stream:
            header = struct.unpack('<8s8Q', stream.read(72))
            assert header[:7] == (b'NMHVSEW\0', 1, 10, 16, 4, 1, 518048)
            prime = header[7]
            assert prime == entry['prime'] and prime not in candidate['crt_primes']
            assert header[8] == int(sha(path.parent.parent / 'input_identity.txt')[:16], 16)
            assert path.stat().st_size == 72 + 8 * 518048
            values = array('Q')
            values.fromfile(stream, 518048)
        if sys.byteorder != 'little':
            values.byteswap()
        assert all(value < prime for value in values)
        arrays.append(values)
        checks.append({'prime_index': index, 'prime': prime, 'path': str(path),
                       'sha256': entry['sha256'], 'mismatches': 0, 'noninvertible': 0,
                       'first_failing_columns': []})
    count = 0
    with vector.open() as stream:
        assert next(stream).startswith('# speculative')
        for line in stream:
            column, numerator, denominator = map(int, line.split())
            assert column == count and denominator > 0 and math.gcd(numerator, denominator) == 1
            for check, values in zip(checks, arrays):
                prime = check['prime']
                noninvertible = denominator % prime == 0
                mismatch = (numerator - values[column] * denominator) % prime != 0
                check['noninvertible'] += noninvertible
                check['mismatches'] += mismatch
                if (noninvertible or mismatch) and len(check['first_failing_columns']) < 12:
                    check['first_failing_columns'].append(column)
            count += 1
    assert count == 518048
    all_pass = all(c['mismatches'] == c['noninvertible'] == 0 for c in checks)
    report = {'schema': 'hexagon-frozen-candidate-fresh-check-v1',
              'source_sha256': sha(Path(__file__)), 'recorded_unix_seconds': time.time(),
              'candidate_receipt': str(candidate_receipt.resolve()),
              'candidate_receipt_sha256': sha(candidate_receipt),
              'candidate_sha256': candidate['candidate_sha256'], 'columns': count,
              'checks': checks, 'all_requested_checks_pass': all_pass,
              'both_fresh_checks_pass': all_pass and sorted(indices) == [7, 8],
              'production_modified': False, 'certified': False,
              'scope': 'Fresh modular vector checks only. Independently validated QQ export, exact original 16+4 sewing and 44 exact origin holdouts remain required. Bonus assumptions remain conditional.',
              'wall_seconds': time.monotonic() - started}
    with output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
