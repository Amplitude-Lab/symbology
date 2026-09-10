from __future__ import annotations

import json
import re
import shutil
import threading
import functools
import inspect
import hashlib
from datetime import datetime, timezone
from pathlib import Path

from .config import PROJECTS_DIR
from .validation import portable_stem

_lock = threading.RLock()
_project_locks: dict[str, threading.RLock] = {}


def transaction_for(pid: str):
    """Serialize a complete read/modify/write operation in this server process."""
    with _lock:
        lock = _project_locks.setdefault(pid, threading.RLock())
    def decorate(fn):
        @functools.wraps(fn)
        def wrapped(*args, **kwargs):
            with lock:
                return fn(*args, **kwargs)
        wrapped.__signature__ = inspect.signature(fn, eval_str=True)
        return wrapped
    return decorate


def transaction(fn):
    signature = inspect.signature(fn, eval_str=True)
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        pid = signature.bind(*args, **kwargs).arguments.get("pid", "__create__")
        return transaction_for(pid)(fn)(*args, **kwargs)
    wrapped.__signature__ = signature
    return wrapped


def alphabet_version(alpha: dict) -> str:
    fields = ("letters", "variables", "expressions", "expr_loader", "roots")
    return hashlib.sha256(json.dumps({k: alpha.get(k) for k in fields},
                                    sort_keys=True).encode()).hexdigest()

_PID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _valid_pid(pid: str) -> bool:
    return bool(_PID_RE.match(pid or ""))


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _slugify(name: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9_-]+", "-", name.strip()).strip("-").lower()
    return portable_stem(slug or "project", "project")


def project_dir(pid: str) -> Path:
    return PROJECTS_DIR / pid


def _project_file(pid: str) -> Path:
    return project_dir(pid) / "project.json"


def _new_id(base: str) -> str:
    slug = _slugify(base)
    pid, i = slug, 2
    while project_dir(pid).exists():
        pid = f"{slug}-{i}"
        i += 1
    return pid


def init_dirs() -> None:
    PROJECTS_DIR.mkdir(parents=True, exist_ok=True)


def create_project(name: str) -> dict:
    with _lock:
        pid = _new_id(name)
        pdir = project_dir(pid)
        (pdir / "data").mkdir(parents=True)
        (pdir / "output").mkdir(parents=True)
        (pdir / "runs").mkdir(parents=True)
        (pdir / "wolfram_gen").mkdir(parents=True)
        proj = {
            "id": pid,
            "name": name,
            "created_at": _now(),
            "alphabets": [],
            "flows": [],
        }
        save_project(proj)
        return proj


def load_project(pid: str) -> dict | None:
    with _lock:
        if not _valid_pid(pid):
            return None
        f = _project_file(pid)
        if not f.exists():
            return None
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError as exc:
            # A raw 500 with a stack trace helps nobody; name the broken file.
            raise ValueError(
                f"project file is corrupted ({f}: {exc}) — restore it from a "
                "backup or delete the project directory to start over"
            ) from exc


def save_project(proj: dict) -> None:
    with _lock:
        f = _project_file(proj["id"])
        tmp = f.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(proj, indent=2))
        tmp.replace(f)


def delete_project(pid: str) -> bool:
    with _lock:
        if not _valid_pid(pid):
            return False
        pdir = project_dir(pid)
        if not pdir.exists():
            return False
        shutil.rmtree(pdir)
        return True


def list_projects() -> list:
    with _lock:
        out = []
        if not PROJECTS_DIR.exists():
            return out
        for d in sorted(PROJECTS_DIR.iterdir()):
            f = d / "project.json"
            if not f.exists():
                continue
            try:
                p = json.loads(f.read_text())
            except Exception:
                continue
            out.append({
                "id": p["id"],
                "name": p["name"],
                "created_at": p.get("created_at"),
                "n_alphabets": len(p.get("alphabets", [])),
                "n_flows": len(p.get("flows", [])),
            })
        return out


def find_alphabet(proj: dict, aid: str) -> dict | None:
    for a in proj.get("alphabets", []):
        if a["id"] == aid:
            return a
    return None


def find_flow(proj: dict, fid: str) -> dict | None:
    for fl in proj.get("flows", []):
        if fl["id"] == fid:
            return fl
    return None


def find_property(alphabet: dict, prop_id: str) -> dict | None:
    for p in alphabet.get("properties", []):
        if p["id"] == prop_id:
            return p
    return None
