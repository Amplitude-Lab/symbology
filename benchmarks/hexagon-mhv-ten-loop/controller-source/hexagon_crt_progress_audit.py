#!/usr/bin/env python3
"""Read-only, independent balanced rational-reconstruction diagnostics.

This does not publish a symbol or replace the producer's reconstruction and
certification. It reads committed shifted caches named in a runner ledger.
"""
import argparse
from array import array
import hashlib
import json
from math import gcd, isqrt
import os
from pathlib import Path
import resource
import struct
import sys
import time


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def reconstruct(residue, modulus, bound):
    a, b, u, v = modulus, residue, 0, 1
    while b > bound:
        q, remainder = divmod(a, b)
        a, b, u, v = b, remainder, v, u - q * v
    if v < 0:
        b, v = -b, -v
    if not (0 < v <= bound and abs(b) <= bound and gcd(b, v) == 1):
        return None
    if (b - residue * v) % modulus:
        raise AssertionError('Reconstruction congruence failed')
    return b, v


def controls():
    count = 0
    for modulus in (31, 101, 1009, 1000003, 1152921504606847009):
        bound = isqrt((modulus - 1) // 2)
        limit = min(bound, 40)
        for denominator in range(1, limit + 1):
            for numerator in range(-limit, limit + 1):
                if gcd(numerator, denominator) != 1:
                    continue
                residue = numerator * pow(denominator, -1, modulus) % modulus
                assert reconstruct(residue, modulus, bound) == (numerator, denominator)
                count += 1
    # Exhaustively compare every residue with all admissible reduced fractions.
    for modulus in (31, 101):
        bound = isqrt((modulus - 1) // 2)
        expected = {}
        for denominator in range(1, bound + 1):
            for numerator in range(-bound, bound + 1):
                if gcd(numerator, denominator) != 1:
                    continue
                residue = numerator * pow(denominator, -1, modulus) % modulus
                assert residue not in expected
                expected[residue] = numerator, denominator
        for residue in range(modulus):
            assert reconstruct(residue, modulus, bound) == expected.get(residue)
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
    started = time.monotonic()
    count = controls()
    ledger_path = args.control / 'ledger.json'
    ledger_bytes = ledger_path.read_bytes()
    ledger = json.loads(ledger_bytes)
    committed = sorted(ledger['primes'].items(), key=lambda item: int(item[0]))
    caches, primes, receipts = [], [], []
    shape = None
    for index, entry in committed:
        cache = entry['conversion']
        path = Path(cache['path'])
        assert path.is_file() and not path.is_symlink()
        digest = sha(path)
        assert digest == cache['sha256'] == path.with_suffix('.bin.sha256').read_text().strip()
        identity = sha(path.parent.parent / 'input_identity.txt')
        with path.open('rb') as stream:
            magic, version, loop, head, tail, rows, cols, prime, fingerprint = struct.unpack('<8s8Q', stream.read(72))
            assert (magic, version, loop, head, tail, rows, fingerprint) == (b'NMHVSEW\0', 1, 10, 16, 4, 1, int(identity[:16], 16))
            assert cols == cache['columns'] and prime == cache['prime']
            assert path.stat().st_size == 72 + 8 * rows * cols
            values = array('Q')
            values.fromfile(stream, rows * cols)
        if sys.byteorder != 'little':
            values.byteswap()
        assert all(value < prime for value in values)
        current_shape = rows, cols, identity
        if shape is not None:
            assert shape == current_shape
        shape = current_shape
        caches.append(values)
        primes.append(prime)
        receipts.append({'index': int(index), 'path': str(path), 'sha256': digest})
    assert len(caches) >= 3
    columns = shape[1]
    residues = list(caches[0])
    modulus = primes[0]
    cases = []
    for crt_count in range(1, len(caches) - 1):
        if crt_count > 1:
            prime = primes[crt_count - 1]
            inverse = pow(modulus, -1, prime)
            residues = [value + modulus * ((other - value) * inverse % prime)
                        for value, other in zip(residues, caches[crt_count - 1])]
            modulus *= prime
        bound = isqrt((modulus - 1) // 2)
        checks = list(zip(primes[crt_count:crt_count + 2], caches[crt_count:crt_count + 2]))
        result = {'crt_primes': primes[:crt_count], 'validation_primes': [p for p, _ in checks],
                  'modulus_bits': modulus.bit_length(), 'balanced_bound_bits': bound.bit_length(),
                  'columns': columns, 'no_bounded_candidate': 0, 'validation_mismatch': 0,
                  'validation_denominator_noninvertible': 0, 'passed_both_fresh_checks': 0,
                  'max_passed_numerator_bits': 0, 'max_passed_denominator_bits': 0,
                  'first_failure_columns_zero_based': []}
        for column, residue in enumerate(residues):
            candidate = reconstruct(residue, modulus, bound)
            if candidate is None:
                category = 'no_bounded_candidate'
            else:
                numerator, denominator = candidate
                if any(denominator % prime == 0 for prime, _ in checks):
                    category = 'validation_denominator_noninvertible'
                elif any((numerator - values[column] * denominator) % prime for prime, values in checks):
                    category = 'validation_mismatch'
                else:
                    category = 'passed_both_fresh_checks'
                    result['max_passed_numerator_bits'] = max(result['max_passed_numerator_bits'], abs(numerator).bit_length())
                    result['max_passed_denominator_bits'] = max(result['max_passed_denominator_bits'], denominator.bit_length())
            result[category] += 1
            if category != 'passed_both_fresh_checks' and len(result['first_failure_columns_zero_based']) < 8:
                result['first_failure_columns_zero_based'].append({'column': column, 'category': category})
        cases.append(result)
        print(json.dumps(result), flush=True)
    report = {'schema': 'hexagon-independent-crt-progress-v1', 'scope': 'Independent bounded Euclidean reconstruction diagnostic over every shifted coefficient. This is not the production reconstruction receipt or a rational symbol certificate.',
              'source_sha256': sha(Path(__file__)), 'ledger_snapshot_sha256': hashlib.sha256(ledger_bytes).hexdigest(),
              'synthetic_and_exhaustive_controls_passed': count, 'cache_receipts': receipts,
              'cases': cases, 'wall_seconds': time.monotonic() - started,
              'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
              'resource_limits': {'virtual_bytes': 512 << 20, 'cpu_seconds': 120, 'nice_increment': 10},
              'production_modified': False}
    assert not args.output.exists(), 'Preserve existing diagnostic evidence'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'receipt': str(args.output), 'wall_seconds': report['wall_seconds'], 'peak_rss_bytes': report['peak_rss_bytes']}))


if __name__ == '__main__':
    main()
