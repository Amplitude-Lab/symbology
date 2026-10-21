#!/usr/bin/env python3
"""Explicit V16 -> V17 cache import, preserving all coefficients and originals.

Only the validated scheduling policy, executable identity and cache fingerprint
change. This does not claim that the imported prime was computed by V17.
Callers must hold both the old and new supervisor locks while installing it.
"""
import hashlib
import json
import os
from pathlib import Path
import struct

HEADER = struct.Struct('<8s8Q')
OLD_ALGORITHM = 'joint-discovery-physical-affine-v16-shared-mrk-shuffles'
NEW_ALGORITHM = 'joint-discovery-physical-affine-v17-resumable-physics-cycle'
NEW_SCHEDULE = '4 light rounds; persistent COL/MRK cycle cursor across 120-second callbacks; due OPE precedes the next cycle'
ALLOWED = {'algorithm', 'early_physics', 'executable_sha256'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def ordinary(path):
    path = Path(path)
    require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents]),
            'Expected ordinary file without symlink ancestors: ' + str(path))
    return path


def fields(path):
    rows = [line.split('\t', 1) for line in ordinary(path).read_text().splitlines()]
    require(all(len(row) == 2 for row in rows), 'Invalid TSV')
    require(len(dict(rows)) == len(rows), 'Duplicate TSV key')
    return dict(rows)


def tsv(value):
    return ''.join(k + '\t' + v + '\n' for k, v in sorted(value.items()))


def install(path, data):
    path = Path(path)
    if isinstance(data, str):
        data = data.encode()
    require(not any(p.is_symlink() for p in [path, *path.parents]), 'Symlink destination')
    if path.exists():
        require(ordinary(path).read_bytes() == data, 'Conflicting import destination: ' + str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    with tmp.open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    tmp.replace(path)


def cache(path, identity_path, h, k, prime):
    data = ordinary(path).read_bytes()
    require(ordinary(str(path) + '.sha256').read_text().strip() == sha(path), 'Cache checksum mismatch')
    require(len(data) >= HEADER.size, 'Truncated cache')
    header = HEADER.unpack_from(data)
    magic, version, loop, head, tail, rows, columns, p, fingerprint = header
    require((magic, version, loop, head, tail, rows, p, fingerprint) ==
            (b'NMHVSEW\0', 1, (h+k)//2, h, k, 1, prime, int(sha(identity_path)[:16], 16)),
            'Unique cache header/identity mismatch')
    require(columns > 0 and len(data) == HEADER.size + 8 * columns, 'Cache length mismatch')
    require(all(x[0] < prime for x in struct.iter_unpack('<Q', data[72:])), 'Noncanonical residue')
    return data


def rebind_cache(source, target, source_identity, target_identity, h, k, prime):
    data = cache(source, source_identity, h, k, prime)
    rebound = data[:64] + struct.pack('<Q', int(sha(target_identity)[:16], 16)) + data[72:]
    install(target, rebound)
    install(str(target) + '.sha256', sha(target) + '\n')
    require(cache(target, target_identity, h, k, prime) == rebound, 'Import readback mismatch')
    return {'source': str(source), 'destination': str(target), 'source_sha256': sha(source),
            'destination_sha256': sha(target), 'payload_sha256': hashlib.sha256(data[72:]).hexdigest(),
            'header_bytes_changed': [64, 72], 'every_coefficient_unchanged': True}


def native_identity(stage, old_identity):
    old_lines = old_identity.splitlines(keepends=True)
    require(old_lines and old_lines[0].startswith('irr-bulk-v1 irr-physics-stream-v1 '), 'Unknown native identity schema')
    first = 'irr-bulk-v1 irr-physics-stream-v1' + ''.join(' ' + k + '=' + v for k, v in sorted(stage.items())) + '\n'
    return first + ''.join(old_lines[1:])


def import_native(source, destination, driver, prime):
    source, destination, driver = map(Path, [source, destination, driver])
    stage = fields(source / 'stage.tsv')
    require(stage['algorithm'] == OLD_ALGORITHM, 'Only validated V16 imports are supported')
    require(stage['object_kind'] == 'affine' and stage['helicity'] == 'MHV', 'Unsupported object convention')
    old_identity = ordinary(source / 'input_identity.txt').read_text()
    require(native_identity(stage, old_identity) == old_identity, 'Source descriptor/identity disagreement')
    new_stage = dict(stage, algorithm=NEW_ALGORITHM, early_physics=NEW_SCHEDULE, executable_sha256=sha(ordinary(driver)))
    require({k for k in stage if stage[k] != new_stage[k]} == ALLOWED, 'Unexpected descriptor difference')
    joint_inputs = ordinary(source.parent / 'joint_inputs.txt')
    require(sha(joint_inputs) == stage['joint_input_sha256'], 'Physical input identity mismatch')
    install(destination.parent / 'joint_inputs.txt', joint_inputs.read_bytes())
    install(destination / 'input_identity.txt', native_identity(new_stage, old_identity))
    install(destination / 'stage.tsv', tsv(new_stage))
    install(destination / 'basis_identity.tsv', ordinary(source / 'basis_identity.tsv').read_bytes())
    install(destination / 'manifest.json', json.dumps({'schema': 'hexagon-irr-physics-v1', 'complete': False,
        'status': 'cache-only', 'history': stage['history']}, separators=(',', ':')) + '\n')
    return rebind_cache(source / 'cache' / f'p{prime}.bin', destination / 'cache' / f'p{prime}.bin',
        source / 'input_identity.txt', destination / 'input_identity.txt', int(stage['h']), int(stage['k']), prime)


def import_shifted(source, destination, old_native, new_native, prime):
    source, destination, old_native, new_native = map(Path, [source, destination, old_native, new_native])
    stage = fields(source / 'stage.tsv')
    require(stage['schema'] == 'irr-shift-cache-v1' and stage['operator'] == 'shift-cache', 'Unsupported shift schema')
    require(stage['kind'] == 'unique' and stage['rows'] == '1' and Path(stage['parent']) == old_native,
            'Shift parent or origin convention mismatch')
    identity = ordinary(source / 'input_identity.txt').read_text().splitlines(keepends=True)
    require(identity[0] == 'irr-shift-cache-v1\n' and identity[1] == sha(old_native / 'input_identity.txt') + '\n',
            'Shift source identity mismatch')
    require((old_native / 'basis_identity.tsv').read_bytes() == (new_native / 'basis_identity.tsv').read_bytes(),
            'Recurrence basis changed')
    require(ordinary(source / f'source_p{prime}.sha256').read_text() == sha(old_native / 'cache' / f'p{prime}.bin') + '\n',
            'Original shift provenance mismatch')
    identity[1] = sha(new_native / 'input_identity.txt') + '\n'
    install(destination / 'input_identity.txt', ''.join(identity))
    install(destination / 'basis_identity.tsv', ordinary(source / 'basis_identity.tsv').read_bytes())
    install(destination / 'stage.tsv', tsv(dict(stage, parent=str(new_native))))
    manifest = json.loads(ordinary(source / 'manifest.json').read_text())
    require(manifest.get('complete') is False, 'Only cache-only shift imports are supported')
    install(destination / 'manifest.json', ordinary(source / 'manifest.json').read_bytes())
    install(destination / f'source_p{prime}.sha256', sha(new_native / 'cache' / f'p{prime}.bin') + '\n')
    return rebind_cache(source / 'cache' / f'p{prime}.bin', destination / 'cache' / f'p{prime}.bin',
        source / 'input_identity.txt', destination / 'input_identity.txt', int(stage['h']), int(stage['k']), prime)


def verify_receipt(path):
    proof = json.loads(ordinary(path).read_text())
    require(proof['schema'] == 'hexagon-explicit-cache-migration-v1', 'Migration schema mismatch')
    for name, digest in proof['pinned_files'].items():
        require(sha(ordinary(name)) == digest, 'Migration artifact changed: ' + name)
    for entry in proof['caches']:
        source, destination = Path(entry['source']), Path(entry['destination'])
        actual, imported = ordinary(source).read_bytes(), ordinary(destination).read_bytes()
        require(sha(source) == entry['source_sha256'] and sha(destination) == entry['destination_sha256'], 'Migration cache changed')
        require(actual[:64] == imported[:64] and actual[72:] == imported[72:], 'Migration changed coefficients or dimensions')
        require(hashlib.sha256(actual[72:]).hexdigest() == entry['payload_sha256'], 'Migration payload changed')
        for file in [source, destination]:
            require(ordinary(str(file) + '.sha256').read_text().strip() == sha(file), 'Migration checksum changed')
    require(proof['native_prime_recomputed'] is False and proof['original_producer_algorithm'] == OLD_ALGORITHM,
            'Migration provenance incomplete')
    return proof
