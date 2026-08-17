from __future__ import annotations

import os
import shutil
from pathlib import Path

SERVER_DIR = Path(__file__).resolve().parents[1]
FRONTEND_DIR = SERVER_DIR.parent
REPO_ROOT = FRONTEND_DIR.parent
PROJECTS_DIR = FRONTEND_DIR / "projects"
WEB_DIST = FRONTEND_DIR / "web" / "dist"

HOST = "127.0.0.1"
PORT = 8321

_WOLFRAM_CANDIDATES = [
    "wolframscript",
    "/Applications/Mathematica.app/Contents/MacOS/wolframscript",
    "/Applications/Wolfram.app/Contents/MacOS/wolframscript",
    "/usr/local/bin/wolframscript",
    "/opt/Wolfram/WolframEngine/Executables/wolframscript",
    "C:\\Program Files\\Wolfram Research\\Mathematica\\13.0\\wolframscript.exe",
    "C:\\Program Files\\Wolfram Research\\Mathematica\\14.0\\wolframscript.exe",
]


def find_wolframscript() -> str | None:
    env = os.environ.get("WOLFRAMSCRIPT")
    if env and Path(env).exists():
        return env
    for cand in _WOLFRAM_CANDIDATES:
        found = shutil.which(cand) if os.sep not in cand else (cand if Path(cand).exists() else None)
        if found:
            return found
    return None


def find_bootstrap() -> str | None:
    for name in ("bootstrap", "bootstrap.exe"):
        p = REPO_ROOT / name
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    return None


def find_compute_rhs() -> str | None:
    for name in ("compute_rhs", "compute_rhs.exe"):
        p = REPO_ROOT / name
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    return None


def find_tensor_add() -> str | None:
    for name in ("tensor_add", "tensor_add.exe"):
        p = REPO_ROOT / name
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    return None


def find_tensor_ops() -> str | None:
    for name in ("tensor_ops", "tensor_ops.exe"):
        p = REPO_ROOT / name
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    return None


def env_status() -> dict:
    ws = find_wolframscript()
    bs = find_bootstrap()
    return {
        "repo_root": str(REPO_ROOT),
        "wolframscript": {"found": ws is not None, "path": ws},
        "bootstrap": {"found": bs is not None, "path": bs},
        "compute_rhs": {"found": find_compute_rhs() is not None, "path": find_compute_rhs()},
        "tensor_add": {"found": find_tensor_add() is not None, "path": find_tensor_add()},
        "tensor_ops": {"found": find_tensor_ops() is not None, "path": find_tensor_ops()},
        "projects_dir": str(PROJECTS_DIR),
    }
