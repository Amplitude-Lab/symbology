"""Content identities and isolated, per-attempt output publication.

The compiler lists file/directory inputs. Directory-driven native commands
need a conservative snapshot; ordinary tensor operations copy only their
explicit inputs. No hard links: a child must never modify a live input.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


class project_lease:
    """Shared admission lock for local and exported runners on one machine."""
    def __init__(self, root):
        self.file = None
        (root / 'runs').mkdir(parents=True, exist_ok=True)
        f = (root / 'runs' / '.execution.lock').open('a+b')
        try:
            if os.name == 'posix':
                import fcntl
                fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            else:
                import msvcrt
                if not f.tell():
                    f.write(b'0'); f.flush()
                f.seek(0)
                msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as exc:
            f.close()
            raise ValueError('Another local or exported run is using this project; wait for it to finish.') from exc
        self.file = f

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def __enter__(self): return self
    def __exit__(self, *exc): self.close()
    def __del__(self): self.close()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def input_manifest(argv, cwd, outputs):
    root = Path(cwd)
    out = {str((root / o).resolve()) for o in outputs}
    files, dirs = set(), set()
    if '-script' in argv:
        dirs.update(str(root / d) for d in ('data', 'output'))
    for i, arg in enumerate(argv):
        p = Path(str(arg))
        p = p if p.is_absolute() else root / p
        try:
            resolved = str(p.resolve())
            is_file = p.is_file()
        except (OSError, ValueError):
            # Literal coefficients/program arguments can exceed filename limits.
            # They still participate in the command fingerprint.
            continue
        if resolved in out:
            continue
        if i and argv[i-1] in ('--data-dir', '--output-dir'):
            dirs.add(str(p.resolve()))
        elif i == 0 or p.suffix in ('.wxf', '.wl') or is_file:
            files.add(str(p.resolve()))
    return sorted(files), sorted(dirs)


def input_files(step, root):
    files, dirs = input_manifest(step.get('argv', []), step.get('cwd') or root, step.get('outputs', []))
    files = set(step.get('inputs', files))
    dirs = step.get('input_dirs', dirs)
    outputs = [(root / o).resolve() for o in step.get('outputs', [])]
    outputs += [Path(str(o) + '.result.json') for o in outputs]
    for directory in dirs:
        p = Path(directory)
        if p.is_dir():
            files.update(str(f) for f in p.rglob('*') if f.is_file()
                         and not f.name.endswith(('.sig', '.meaning.json', '.tmp')))
        else:
            files.add(str(p))
    return sorted(Path(f) for f in files if not any(Path(f).resolve() == o or o in Path(f).resolve().parents for o in outputs))


def fingerprint(step, root):
    h = hashlib.sha256(b'symbology-cache-v3\0')
    h.update('\0'.join(step.get('argv', [])).encode())
    for p in input_files(step, root):
        h.update(str(p).encode())
        h.update((digest(p) if p.is_file() else 'missing').encode())
    return h.hexdigest()


def step_cached(step, root):
    outputs = step.get('outputs') or []
    if not outputs:
        return False
    fp = fingerprint(step, root)
    for o in outputs:
        p = root / o
        try:
            record = json.loads(Path(str(p) + '.sig').read_text())
            if not p.is_file() or not p.stat().st_size or record != {'input': fp, 'output': digest(p)}:
                return False
        except (OSError, ValueError):
            return False
    return True


def record_sigs(step, root, fp=None):
    fp = fp or fingerprint(step, root)
    for o in step.get('outputs') or []:
        p = root / o
        if p.is_file():
            sig = Path(str(p) + '.sig')
            tmp = Path(str(sig) + '.tmp')
            tmp.write_text(json.dumps({'input': fp, 'output': digest(p)}))
            os.replace(tmp, sig)


def run_isolated(run, step, root, launch):
    outputs = step.get('outputs', [])
    if any(Path(o).is_absolute() or '..' in Path(o).parts for o in outputs):
        raise ValueError('Step outputs must stay inside the project.')
    before = fingerprint(step, root)
    with tempfile.TemporaryDirectory(prefix=f'.attempt-{run.run_id}-', dir=root / 'runs') as tmp:
        stage = Path(tmp)
        def staged(p):
            return stage / p.relative_to(root)
        originals = {}
        for p in input_files(step, root):
            if not p.is_file() or not p.is_relative_to(root):
                continue
            dest = staged(p)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest)
            originals[dest.relative_to(stage)] = digest(p)
        for directory in ('data', 'output', 'wolfram_gen'):
            (stage / directory).mkdir(exist_ok=True)
        for o in outputs:
            p = stage / o
            p.parent.mkdir(parents=True, exist_ok=True)
        # Generated Wolfram programs embed absolute project paths.
        for p in stage.rglob('*.wl'):
            text = p.read_text()
            for old, new in step.get('path_replacements', []):
                text = text.replace(old, new)
            p.write_text(text.replace(str(root), str(stage)))
        argv = [a.replace(str(root) + os.sep, str(stage) + os.sep) if i else a
                for i, a in enumerate(step['argv'])]
        argv = [str(stage) if a == str(root) else a for a in argv]
        attempt = {**step, 'argv': argv, 'cwd': str(stage)}
        rc = launch(run, attempt)
        if 'result' in attempt:
            step['result'] = attempt['result']
        if rc or run.cancel_requested:
            return rc or 130
        missing = [o for o in outputs if not (stage / o).exists()
                   or ((stage / o).is_file() and (stage / o).stat().st_size == 0)]
        if missing:
            raise ValueError('Step exited 0 but its declared outputs are missing or empty: ' + ', '.join(missing))
        if fingerprint(step, root) != before:
            raise ValueError('Inputs changed while this step was running; outputs were not published. Run again.')
        publish = []
        for folder in ('data', 'output'):
            for p in (stage / folder).rglob('*'):
                if p.is_file():
                    rel = p.relative_to(stage)
                    if p.suffix == '.wl' and rel in originals:
                        continue  # relocated input programs are never published
                    if originals.get(rel) != digest(p):
                        publish.append((p, root / rel))
        # Directory-driven commands can invalidate old derived outputs (for
        # example a collinear projection becoming empty). Publish those
        # removals too; data inputs and relocated programs remain protected.
        for rel in originals:
            if rel.parts[0] == 'output' and not (stage / rel).exists() and rel.suffix != '.wl':
                publish.append((None, root / rel))
        # Each file is atomically replaced on the same filesystem. Retain
        # rollback copies until the entire validated set has been published.
        backups, committed = [], []
        with run.cond:
            if run.cancel_requested:
                return 130
            try:
                for i, (src, dst) in enumerate(publish):
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    if dst.exists():
                        backup = stage / f'.rollback-{i}'
                        shutil.copy2(dst, backup)
                        backups.append((backup, dst))
                    if src is None:
                        dst.unlink(missing_ok=True)
                    else:
                        os.replace(src, dst)
                    committed.append(dst)
                for o in outputs:
                    if (stage / o).is_dir():
                        (root / o).mkdir(parents=True, exist_ok=True)
            except Exception:
                for dst in committed:
                    dst.unlink(missing_ok=True)
                for backup, dst in backups:
                    os.replace(backup, dst)
                raise
        return 0
