"""Remove only the reviewed abandoned Git objects and expired mmap arenas.

This is a dated, one-use cleanup receipt generator, not a general prune tool.
Numerical checkpoint, coefficient, input, source, and provenance files are excluded.
"""
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import time

ROOT = Path('/home/ana/Documents/symbology')
AUDIT = Path(__file__).resolve().parent


def save(name, data):
    (AUDIT / name).write_text(json.dumps(data, indent=2) + '\n')


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True, check=True)


def identity(path):
    s = path.lstat()
    return dict(inode=s.st_ino, device=s.st_dev, bytes=s.st_size,
                allocated=s.st_blocks * 512, mtime_ns=s.st_mtime_ns,
                links=s.st_nlink, mode=s.st_mode)


def active_uses(targets):
    hits = []
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            comm = (proc / 'comm').read_text().strip()
            if comm == 'git' or comm.startswith('git-'):
                hits.append(dict(pid=proc.name, reason='Git process active'))
            for fd in (proc / 'fd').iterdir():
                try:
                    p = os.readlink(fd)
                    if p in targets:
                        hits.append(dict(pid=proc.name, path=p, reason='open file'))
                except OSError:
                    pass
            for line in (proc / 'maps').read_text().splitlines():
                fields = line.split(maxsplit=5)
                if len(fields) == 6 and fields[5] in targets:
                    hits.append(dict(pid=proc.name, path=fields[5], reason='memory mapping'))
        except (OSError, PermissionError):
            pass
    return hits


def git_state():
    status = git('status', '--short', '--untracked-files=no').stdout
    return dict(HEAD=git('rev-parse', 'HEAD').stdout.strip(),
                refs=git('show-ref').stdout, tracked_status=status,
                index_sha256=hashlib.sha256((ROOT / '.git/index').read_bytes()).hexdigest())


assert not (AUDIT / 'deleted.jsonl').exists(), 'One-use cleanup already started'
inventory = json.loads((AUDIT / 'large_files_before.json').read_text())
garbage = json.loads((AUDIT / 'git-garbage-candidates.json').read_text())
scratch = [x for x in inventory if x['path'].endswith('.mmap')]
assert len(garbage) == 29 and len(scratch) == 8
allowed_scratch_roots = {
    'output_hexagon_17p3_20261007',
    'output_hexagon_runtime_recovery_20261008',
    'output_hexagon_validation_20261009',
}
current_garbage = git('count-objects', '-vH').stderr
current_garbage_paths = {line.split('garbage found:', 1)[1].strip()
                         for line in current_garbage.splitlines() if 'garbage found:' in line}
assert current_garbage_paths == {x['path'] for x in garbage}
tracked = set(git('ls-files', '-z').stdout.split('\0'))
records = []
for category, candidates in [('abandoned_git_temporary_file', garbage),
                             ('expired_solver_scratch_arena', scratch)]:
    for candidate in candidates:
        rel = candidate['path']
        p = ROOT / rel
        s = identity(p)
        assert rel not in tracked
        assert stat.S_ISREG(s['mode']) and s['links'] == 1
        assert s['bytes'] == candidate['bytes'] and s['allocated'] == candidate['allocated']
        assert abs(s['mtime_ns'] / 1e9 - candidate['mtime']) < 1e-5
        assert time.time() - s['mtime_ns'] / 1e9 > 86400
        assert p.resolve() == p
        if category == 'abandoned_git_temporary_file':
            assert rel.startswith('.git/objects/') and p.name.startswith(('tmp_pack_', 'tmp_obj_'))
            reason = 'Git explicitly reports garbage; no Git process, lock, open descriptor or mapping; connectivity checked before removal.'
        else:
            assert Path(rel).parts[0] in allowed_scratch_roots
            reason = ('Expired process-local memory arena from completed/interrupted historical run. '
                      'bench/hexagon_mapped_buffers.hpp requires O_CREAT|O_EXCL and cannot resume this file; '
                      'bench/hexagon_ten_loop_runner_v17.py explicitly unlinks scratch. '
                      'Durable checkpoint/cache files and all run receipts remain.')
        records.append(dict(path=rel, category=category, reason=reason, **s))

targets = {str(ROOT / x['path']) for x in records}
assert not active_uses(targets), 'Candidate is in use or Git is running'
assert not list((ROOT / '.git').rglob('*.lock')), 'Git lock present'
fsck = git('fsck', '--connectivity-only', '--no-dangling')
save('git-connectivity-before.json', dict(exit_code=fsck.returncode, stdout=fsck.stdout, stderr=fsck.stderr))
before = git_state()
protected = {x['path']: identity(ROOT / x['path']) for x in inventory
             if str(ROOT / x['path']) not in targets and (ROOT / x['path']).exists()}
documents = {}
for name in ('cluster_symbols.tex', 'cluster_symbols.pdf',
             'cluster_symbols_all_loop.tex', 'cluster_symbols_all_loop.pdf'):
    p = ROOT / 'research/4pformfactor/tex' / name
    if p.exists():
        documents[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
save('protected-large-file-stats.json', protected)
save('deletion-manifest.json', dict(files=records,
     evidence=['git-count-before.txt', 'pre_cleanup_state.json', 'git-connectivity-before.json',
               'bench/hexagon_mapped_buffers.hpp:98', 'bench/hexagon_ten_loop_runner_v17.py:261'],
     hash_policy='Metadata recorded; no content hashing of disposable scratch arenas or incomplete Git temporaries. Sparse arenas have enormous unallocated regions.'))
assert not active_uses(targets), 'Candidate became active'
v = os.statvfs(ROOT)
free_before = v.f_bavail * v.f_frsize
started = time.time()
with (AUDIT / 'deleted.jsonl').open('x') as journal:
    for record in records:
        p = ROOT / record['path']
        assert identity(p) == {k: record[k] for k in identity(p)}, 'File changed before deletion'
        p.unlink()
        journal.write(json.dumps(dict(**record, removed_unix_seconds=time.time())) + '\n')
        journal.flush()
        os.fsync(journal.fileno())

after = git_state()
fsck = git('fsck', '--connectivity-only', '--no-dangling')
save('git-connectivity-after.json', dict(exit_code=fsck.returncode, stdout=fsck.stdout, stderr=fsck.stderr))
count = git('count-objects', '-vH')
(AUDIT / 'git-count-after.txt').write_text(count.stdout + count.stderr)
changed_protected = [p for p, s in protected.items() if not (ROOT / p).exists() or identity(ROOT / p) != s]
changed_documents = [p for p, h in documents.items() if hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h]
v = os.statvfs(ROOT)
by_category = {}
for record in records:
    c = by_category.setdefault(record['category'], dict(files=0, allocated_bytes_removed=0))
    c['files'] += 1
    c['allocated_bytes_removed'] += record['allocated']
receipt = dict(status='completed', files_removed=len(records),
               allocated_bytes_removed=sum(x['allocated'] for x in records),
               logical_bytes_removed=sum(x['bytes'] for x in records),
               observed_free_bytes_change=v.f_bavail * v.f_frsize - free_before,
               categories=by_category, wall_seconds=time.time() - started,
               git_before=before, git_after=after, git_state_unchanged=(before == after),
               git_connectivity_passed=(fsck.returncode == 0),
               remaining_git_garbage_report=count.stderr,
               protected_large_files_checked=len(protected),
               changed_protected_large_files=changed_protected,
               document_sha256=documents, changed_documents=changed_documents,
               preserved='All durable numerical checkpoints, coefficients, final outputs, sources and provenance; entire current V27D run.')
save('cleanup_receipt.json', receipt)
assert before == after and not changed_protected and not changed_documents
assert not count.stderr
print(json.dumps({k: v for k, v in receipt.items() if k not in ('git_before', 'git_after', 'document_sha256')}, indent=2))
