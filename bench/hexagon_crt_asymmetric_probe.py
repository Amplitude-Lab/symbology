#!/usr/bin/env python3
"""Bounded read-only probe of asymmetric rational reconstruction; no export."""
import argparse
from array import array
from math import gcd
import json
import os
from pathlib import Path
import resource
import struct
import sys
import time
from hexagon_crt_progress_audit import sha


def reconstruct(residue, modulus, numerator_bound, denominator_bound):
    assert 2 * numerator_bound * denominator_bound < modulus
    a, b, u, v = modulus, residue, 0, 1
    while b > numerator_bound:
        q, remainder = divmod(a, b)
        a, b, u, v = b, remainder, v, u - q * v
    if v < 0:
        b, v = -b, -v
    if not (0 < v <= denominator_bound and abs(b) <= numerator_bound and gcd(b, v) == 1):
        return None
    assert (b - residue * v) % modulus == 0
    return b, v


def controls():
    count = 0
    for modulus in (31, 101, 211):
        for numerator_bound in (1, 2, 3, 5, 7, 11, 17, 29):
            denominator_bound = (modulus - 1) // (2 * numerator_bound)
            if denominator_bound == 0:
                continue
            expected = {}
            for d in range(1, denominator_bound + 1):
                for n in range(-numerator_bound, numerator_bound + 1):
                    if gcd(n, d) != 1:
                        continue
                    residue = n * pow(d, -1, modulus) % modulus
                    assert residue not in expected
                    expected[residue] = n, d
            for residue in range(modulus):
                assert reconstruct(residue, modulus, numerator_bound, denominator_bound) == expected.get(residue)
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('control', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    os.nice(10)
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    start = time.monotonic()
    checked = controls()
    ledger = json.loads((args.control / 'ledger.json').read_text())
    caches, primes, receipts = [], [], []
    for index in ('1', '2', '3', '4'):
        entry = ledger['primes'][index]['conversion']
        path = Path(entry['path'])
        digest = sha(path)
        assert digest == entry['sha256'] == path.with_suffix('.bin.sha256').read_text().strip()
        with path.open('rb') as stream:
            h = struct.unpack('<8s8Q', stream.read(72))
            assert h[:7] == (b'NMHVSEW\0', 1, 10, 16, 4, 1, 518048)
            assert h[7] == entry['prime'] and h[8] == int(sha(path.parent.parent / 'input_identity.txt')[:16], 16)
            assert path.stat().st_size == 72 + 8 * h[6]
            values = array('Q')
            values.fromfile(stream, h[6])
        if sys.byteorder != 'little':
            values.byteswap()
        assert all(v < h[7] for v in values)
        caches.append(values)
        primes.append(h[7])
        receipts.append({'path': str(path), 'sha256': digest})
    # Fixed, uniformly spaced sample; selection is independent of coefficient values.
    columns = sorted({i * (518048 - 1) // (8192 - 1) for i in range(8192)})
    inverse = pow(primes[0], -1, primes[1])
    modulus = primes[0] * primes[1]
    residues = {c: caches[0][c] + primes[0] * ((caches[1][c] - caches[0][c]) * inverse % primes[1]) for c in columns}
    cases = []
    for bits in (24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90):
        numerator_bound = (1 << bits) - 1
        denominator_bound = (modulus - 1) // (2 * numerator_bound)
        row = {'numerator_bound_bits': bits, 'denominator_bound_bits': denominator_bound.bit_length(),
               'sample_columns': len(columns), 'no_bounded_candidate': 0, 'validation_mismatch': 0,
               'passed_both_checks': 0, 'passed_nonzero': 0}
        for column, residue in residues.items():
            q = reconstruct(residue, modulus, numerator_bound, denominator_bound)
            if q is None:
                row['no_bounded_candidate'] += 1
            elif any(q[1] % p == 0 or (q[0] - cache[column] * q[1]) % p for p, cache in zip(primes[2:], caches[2:])):
                row['validation_mismatch'] += 1
            else:
                row['passed_both_checks'] += 1
                row['passed_nonzero'] += q[0] != 0
        cases.append(row)
    report = {'schema': 'hexagon-asymmetric-reconstruction-probe-v1', 'source_sha256': sha(Path(__file__)),
              'helper_sha256': sha(Path(__file__).with_name('hexagon_crt_progress_audit.py')),
              'controls_passed': checked, 'cache_receipts': receipts, 'crt_primes': primes[:2],
              'check_primes': primes[2:], 'cases': cases,
              'scope': 'A fixed 8192-column sample compares bounds under 2*N*D<M. Existing check primes guide the diagnostic; this is not a new independent full-symbol validation, proof of physical correctness, or production selection.',
              'wall_seconds': time.monotonic() - start, 'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
              'production_modified': False, 'virtual_limit_bytes': 512 << 20, 'cpu_limit_seconds': 120}
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
