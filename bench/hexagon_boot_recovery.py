"""Durable, conservative recovery helpers for an exclusively owned experiment."""
import hashlib
import json
import os
from pathlib import Path
import re
import time
import uuid


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def sync(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
    with temp.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temp.replace(path)
    sync(path.parent)


def quarantine(path, owned_root, archive, reason):
    """Retain uncommitted bytes; never erase or repair a checksum."""
    path, owned_root, archive = map(Path, (path, owned_root, archive))
    relative = path.relative_to(owned_root)
    if path.is_symlink() or path.resolve() != owned_root.resolve() / relative:
        raise RuntimeError('Refusing linked recovery path: ' + str(path))
    destination = archive / uuid.uuid4().hex / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    path.replace(destination)
    sync(path.parent)
    sync(destination.parent)
    atomic_json(destination.with_name(destination.name + '.recovery.json'), {
        'original': str(path), 'retained': str(destination), 'reason': reason,
        'sha256': sha(destination), 'time': time.time()})
    return str(destination)


def recover_files(root, archive):
    """Call only with no producer running and its stage locks held.

    Missing/truncated checksum markers are incomplete transactions. A valid
    marker with different data is corruption and stops execution for review.
    The solver checks identities, headers and every sparse residue on load.
    """
    root = Path(root)
    retained = []
    if not root.exists():
        return retained
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise RuntimeError('Linked calculation output: ' + str(path))
        if not path.is_file():
            continue
        if path.name.endswith('.tmp'):
            retained.append(quarantine(path, root, archive, 'interrupted temporary write'))
            continue
        if path.suffix not in ('.bin', '.csr') or path.parent.name != 'cache':
            continue
        side = path.with_name(path.name + '.sha256')
        digest = side.read_text().strip() if side.exists() else ''
        if not re.fullmatch('[0-9a-f]{64}', digest):
            retained.append(quarantine(path, root, archive, 'checksum commit incomplete'))
            if side.exists():
                retained.append(quarantine(side, root, archive, 'checksum commit incomplete'))
        elif sha(path) != digest:
            raise RuntimeError('Checksum mismatch; preserved for review: ' + str(path))
    # A checksum with no data must not silently certify a subsequently regenerated file.
    for side in sorted(root.rglob('*.sha256')):
        if side.parent.name == 'cache' and side.name.endswith(('.bin.sha256', '.csr.sha256')):
            if not Path(str(side)[:-7]).exists():
                retained.append(quarantine(side, root, archive, 'orphan checksum'))
    return retained


def flush_completed(roots, receipt_path):
    """Fsync newly completed cache transactions without changing solver files.

    The producer has already hashed/read back these files. Hash verification
    is repeated at recovery and by the solver; the receipt is durability
    bookkeeping, never a replacement for those checks.
    """
    receipt_path = Path(receipt_path)
    records = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    changed = False
    for root in map(Path, roots):
        for side in sorted(root.rglob('*.sha256')):
            path = Path(str(side)[:-7])
            if not path.is_file() or side.is_symlink() or path.is_symlink():
                continue
            digest = side.read_text().strip()
            if not re.fullmatch('[0-9a-f]{64}', digest):
                continue
            key = str(path)
            if key in records:
                if records[key]['sha256'] != digest:
                    raise RuntimeError('Completed checkpoint checksum changed: ' + key)
                continue
            sync(path)
            sync(side)
            for meta in ('input_identity.txt', 'basis_identity.tsv', 'stage.tsv'):
                if (path.parent.parent / meta).is_file():
                    sync(path.parent.parent / meta)
            for directory in (path.parent, path.parent.parent, root):
                sync(directory)
            records[key] = {'sha256': digest, 'bytes': path.stat().st_size,
                            'fsynced_unix_seconds': time.time()}
            changed = True
    if changed:
        atomic_json(receipt_path, records)
    return records


def process_cpu(pid):
    """Live user+system seconds for the solver and live descendants (lower bound)."""
    total, pending = 0, [int(pid)]
    while pending:
        current = pending.pop()
        try:
            fields = Path(f'/proc/{current}/stat').read_text().rsplit(') ', 1)[1].split()
            total += int(fields[11]) + int(fields[12])
            pending.extend(map(int, Path(f'/proc/{current}/task/{current}/children').read_text().split()))
        except (OSError, ValueError):
            pass
    return total / os.sysconf('SC_CLK_TCK')
