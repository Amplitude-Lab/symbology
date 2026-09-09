#!/usr/bin/env python3
"""Regression gate: re-run a recorded flow in a hermetic sandbox and compare CRC32s.

Usage:  python3 scripts/regression_check.py [manifest.json]

The manifest (audits/baseline-*/manifest.json) records the exact argv of every
step of a previously verified run, the seed inputs to copy, the expected CRC32
of every produced output, and strings the solver output must contain (the
known-correct solution vector). The script:

  1. creates a fresh sandbox dir (mktemp -d),
  2. copies the recorded seed inputs from the live project into it,
  3. runs each step sequentially (pre-creating output parent dirs the way the
     front-end engine does), failing fast on a nonzero return code,
  4. verifies the solution-marker strings and every recorded CRC32,
  5. removes the sandbox (keep it with KEEP_SANDBOX=1).

Exit code 0 only when every output is bit-identical. This is the gate for any
change to the C++ calculation core: `make regression` after `make`.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = ROOT / "audits" / "baseline-2026-09-09" / "manifest.json"


def die(msg):
    print(f"REGRESSION FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    manifest_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MANIFEST
    manifest = json.loads(manifest_path.read_text())
    # Some baselines read the repo's own data/ + output/ chain (repo_root_inputs)
    # instead of a front-end project directory.
    proj = ROOT if manifest.get("repo_root_inputs") else ROOT / "front-end" / "projects" / manifest["project"]

    for binary in sorted({s["argv"][0] for s in manifest["steps"]}):
        p = Path(binary.replace("@ROOT@", str(ROOT)))
        if not (p.exists() and os.access(p, os.X_OK)):
            die(f"binary {p} missing or not executable — run `make` first")

    sandbox = Path(tempfile.mkdtemp(prefix="symbology-regression."))
    try:
        # Seed inputs: directories are copied wholesale, files individually.
        for rel in manifest["seed_inputs"]:
            src = proj / rel
            dst = sandbox / rel
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            elif src.is_file():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
            else:
                die(f"seed input {rel} missing from {proj}")

        root_s, sandbox_s = str(ROOT), str(sandbox)
        stdout_all = []
        for step in manifest["steps"]:
            argv = [a.replace("@ROOT@", root_s).replace("@PROJ@", sandbox_s) for a in step["argv"]]
            print(f"  [{step['id']}] {step['label']}")
            # The front-end engine pre-creates output parent dirs before each
            # step; replicate that so a missing dir is never what we test.
            for a in argv[1:]:
                if a.endswith(".wxf") and a.startswith(str(sandbox)):
                    Path(a).parent.mkdir(parents=True, exist_ok=True)
            r = subprocess.run(argv, capture_output=True, text=True)
            stdout_all.append((step["id"], r.stdout))
            if r.returncode != 0:
                print(r.stdout[-3000:], r.stderr[-3000:], sep="\n")
                die(f"step {step['id']} exited {r.returncode}")

        check = manifest.get("solution_check") or {}
        for step_id, out in stdout_all:
            if step_id != check.get("step"):
                continue
            for needle in check.get("must_contain", []):
                if needle not in out:
                    die(f"solution marker '{needle}' not printed by {step_id}")

        bad = 0
        for rel, want in sorted(manifest["crc32"].items()):
            f = sandbox / rel
            if not f.is_file():
                print(f"  MISSING  {rel}")
                bad += 1
                continue
            got = f"{zlib.crc32(f.read_bytes()):08x}"
            mark = "ok      " if got == want else "CHANGED "
            print(f"  {mark} {rel}  {got} (want {want})")
            bad += got != want
        if bad:
            die(f"{bad} output(s) differ from the baseline")
        print("REGRESSION PASS: all outputs bit-identical to baseline")
    finally:
        if os.environ.get("KEEP_SANDBOX"):
            print(f"sandbox kept: {sandbox}")
        else:
            shutil.rmtree(sandbox, ignore_errors=True)


if __name__ == "__main__":
    main()
