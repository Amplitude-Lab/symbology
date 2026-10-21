#!/usr/bin/env python3
"""Read-only four-CRT varying-bound diagnostic; never exports a symbol."""
from array import array
import json
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
    for prime, values in zip(primes[1:4], caches[1:4]):
        inverse = pow(modulus, -1, prime)
        residues = [value + modulus * ((other - value) * inverse % prime)
                    for value, other in zip(residues, values)]
        modulus *= prime
    columns = sorted({i * (518048 - 1) // (8192 - 1) for i in range(8192)})

    def evaluate(bits, selected):
        numerator_bound = (1 << bits) - 1
        denominator_bound = (modulus - 1) // (2 * numerator_bound)
        row = {'numerator_bound_bits': bits, 'denominator_bound_bits': denominator_bound.bit_length(),
               'columns': len(selected), 'no_bounded_candidate': 0, 'validation_mismatch': 0,
               'passed_both_checks': 0, 'passed_nonzero': 0,
               'max_passed_numerator_bits': 0, 'max_passed_denominator_bits': 0}
        for column in selected:
            candidate = reconstruct(residues[column], modulus, numerator_bound, denominator_bound)
            if candidate is None:
                row['no_bounded_candidate'] += 1
            elif any(candidate[1] % p == 0 or (candidate[0] - cache[column] * candidate[1]) % p
                     for p, cache in zip(primes[4:], caches[4:])):
                row['validation_mismatch'] += 1
            else:
                numerator, denominator = candidate
                row['passed_both_checks'] += 1
                row['passed_nonzero'] += numerator != 0
                row['max_passed_numerator_bits'] = max(row['max_passed_numerator_bits'], abs(numerator).bit_length())
                row['max_passed_denominator_bits'] = max(row['max_passed_denominator_bits'], denominator.bit_length())
        return row

    ordered_bits = [36, 48, 24, 60, 12] + [bits for bits in range(6, 241, 6) if bits not in (36, 48, 24, 60, 12)]
    pending = set(range(518048))
    accepted = {}
    union_cases = []
    for bits in ordered_bits:
        nb = (1 << bits) - 1
        db = (modulus - 1) // (2 * nb)
        if db < 1:
            continue
        added = []
        for column in pending:
            candidate = reconstruct(residues[column], modulus, nb, db)
            if candidate is not None and all(candidate[1] % p and (candidate[0] - cache[column] * candidate[1]) % p == 0
                                               for p, cache in zip(primes[4:], caches[4:])):
                accepted[column] = candidate
                added.append(column)
        pending.difference_update(added)
        union_cases.append({'numerator_bound_bits': bits, 'denominator_bound_bits': db.bit_length(),
                            'newly_passed_both_checks': len(added), 'remaining': len(pending)})
        if not pending:
            break
    full = {'columns': 518048, 'passed_both_checks': len(accepted), 'remaining': len(pending),
            'passed_nonzero': sum(n != 0 for n, d in accepted.values()),
            'max_passed_numerator_bits': max(abs(n).bit_length() for n, d in accepted.values()),
            'max_passed_denominator_bits': max(d.bit_length() for n, d in accepted.values())}
    report = {'schema': 'hexagon-four-crt-varying-bound-diagnostic-v1',
              'source_sha256': sha(Path(__file__)),
              'helper_sha256': {name: sha(Path(__file__).with_name(name)) for name in
                                ('hexagon_crt_asymmetric_probe.py', 'hexagon_crt_progress_audit.py')},
              'controls_passed': checked, 'cache_receipts': receipts,
              'crt_primes': primes[:4], 'check_primes': primes[4:], 'modulus_bits': modulus.bit_length(),
              'bound_cases': union_cases, 'full_union': full,
              'scope': 'Existing check primes guide the bound comparison and full diagnostic. This is not a new independent validation, production selection, exported symbol or exact certificate.',
              'wall_seconds': time.monotonic() - started,
              'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
              'production_modified': False, 'virtual_limit_bytes': 512 << 20, 'cpu_limit_seconds': 120}
    with output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'receipt': str(output), 'full': full,
                      'wall_seconds': report['wall_seconds']}, indent=2))


if __name__ == '__main__':
    main()
