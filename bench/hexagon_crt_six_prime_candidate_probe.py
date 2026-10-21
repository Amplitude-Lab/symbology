#!/usr/bin/env python3
"""Freeze a speculative six-CRT coefficient vector for future fresh-prime checks; no production export."""
from array import array
import json
import math
import os
from pathlib import Path
import resource
import struct
import sys
import time
from hexagon_crt_asymmetric_probe import reconstruct, controls
from hexagon_crt_progress_audit import sha


def main():
    control, output = map(Path, sys.argv[1:])
    os.nice(10)
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    started = time.monotonic()
    checked = controls()
    ledger_bytes = (control / 'ledger.json').read_bytes()
    ledger = json.loads(ledger_bytes)
    assert '7' not in ledger['primes'], 'Candidate must be frozen before prime 7 is committed'
    caches, primes, receipts = [], [], []
    for index in ('1', '2', '3', '4', '5', '6'):
        entry = ledger['primes'][index]['conversion']
        path = Path(entry['path'])
        assert path.is_file() and not path.is_symlink()
        digest = sha(path)
        assert digest == entry['sha256'] == path.with_suffix('.bin.sha256').read_text().strip()
        with path.open('rb') as stream:
            header = struct.unpack('<8s8Q', stream.read(72))
            assert header[:7] == (b'NMHVSEW\0', 1, 10, 16, 4, 1, 518048)
            assert header[7] == entry['prime']
            assert header[8] == int(sha(path.parent.parent / 'input_identity.txt')[:16], 16)
            assert path.stat().st_size == 72 + 8 * header[6]
            values = array('Q')
            values.fromfile(stream, header[6])
        if sys.byteorder != 'little':
            values.byteswap()
        assert all(value < header[7] for value in values)
        caches.append(values)
        primes.append(header[7])
        receipts.append({'path': str(path), 'sha256': digest})
    residues = list(caches[0])
    modulus = primes[0]
    for prime, values in zip(primes[1:6], caches[1:6]):
        inverse = pow(modulus, -1, prime)
        residues = [value + modulus * ((other - value) * inverse % prime)
                    for value, other in zip(residues, values)]
        modulus *= prime
    numerator_bound = (1 << 48) - 1
    denominator_bound = (modulus - 1) // (2 * numerator_bound)
    candidate_path = output.with_suffix('.candidate.tsv')
    failed = []
    height_n = height_d = 0
    nonzero = 0
    with candidate_path.open('x') as stream:
        stream.write('# speculative shifted 16+4 vector; zero-based column, numerator, denominator; NOT CERTIFIED\n')
        for column, residue in enumerate(residues):
            candidate = reconstruct(residue, modulus, numerator_bound, denominator_bound)
            if candidate is None:
                failed.append(column)
                continue
            numerator, denominator = candidate
            assert denominator > 0 and math.gcd(numerator, denominator) == 1
            assert abs(numerator) <= numerator_bound and denominator <= denominator_bound
            assert all(denominator % p and (numerator - cache[column] * denominator) % p == 0
                       for p, cache in zip(primes, caches))
            height_n = max(height_n, abs(numerator).bit_length())
            height_d = max(height_d, denominator.bit_length())
            nonzero += numerator != 0
            stream.write(f'{column}\t{numerator}\t{denominator}\n')
        stream.flush()
        os.fsync(stream.fileno())
    # Independent text round trip and modular consistency against all six input caches.
    count = 0
    with candidate_path.open() as stream:
        assert next(stream).startswith('# speculative')
        previous = -1
        for line in stream:
            column, numerator, denominator = map(int, line.split())
            assert column > previous and 0 <= column < 518048
            previous = column
            assert denominator > 0 and math.gcd(numerator, denominator) == 1
            assert all(denominator % p and (numerator - cache[column] * denominator) % p == 0
                       for p, cache in zip(primes, caches))
            count += 1
    assert count + len(failed) == 518048
    assert '7' not in json.loads((control / 'ledger.json').read_text())['primes']
    report = {'schema': 'hexagon-six-crt-speculative-candidate-v1',
              'source_sha256': sha(Path(__file__)),
              'helper_sha256': {name: sha(Path(__file__).with_name(name)) for name in
                                ('hexagon_crt_asymmetric_probe.py', 'hexagon_crt_progress_audit.py')},
              'controls_passed': checked, 'cache_receipts': receipts,
              'crt_primes': primes, 'modulus_bits': modulus.bit_length(),
              'numerator_bound_bits': 48, 'denominator_bound_bits': denominator_bound.bit_length(),
              'columns': 518048, 'candidate_rows': count, 'failed_columns': failed,
              'nonzero_coefficients': nonzero, 'max_numerator_bits': height_n, 'max_denominator_bits': height_d,
              'candidate_path': str(candidate_path.resolve()), 'candidate_sha256': sha(candidate_path),
              'candidate_bytes': candidate_path.stat().st_size,
              'all_six_input_residue_checks_and_text_roundtrip': True,
              'candidate_complete': len(failed) == 0, 'fresh_check_primes_passed': [],
              'future_fresh_check_prime_indices': [7, 8], 'prime7_not_committed_at_freeze': True,
              'scope': 'Speculative coefficient vector with fixed 48-bit numerator bound, frozen before primes 7 and 8. Six input residues are fitting data, not independent checks. No QQ export, physical certificate, unrestricted uniqueness or certified result. Requires future checks, independently validated export, original exact 16+4 sewing and 44 exact origin holdouts.',
              'recorded_unix_seconds': time.time(), 'wall_seconds': time.monotonic() - started,
              'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
              'production_modified': False, 'virtual_limit_bytes': 512 << 20, 'cpu_limit_seconds': 120}
    with output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k: report[k] for k in ('candidate_complete', 'candidate_rows', 'failed_columns',
          'max_numerator_bits', 'max_denominator_bits', 'wall_seconds', 'candidate_sha256')}, indent=2))


if __name__ == '__main__':
    main()
