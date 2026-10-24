#!/usr/bin/env python3
"""Resolve the hexagon ten-loop benchmark inputs.

Three sources, in order of preference for ``auto``:

local  The archived production packages on the workstation
       (``output_hexagon_16p4_20261006``, ``output_hexagon_17p3_20261007``,
       ``output_hexagon_runtime_recovery_20261008``,
       ``output_hexagon_validation_20261009``).
bench  The split xz parts kept locally under
       ``benchmarks/hexagon-mhv-ten-loop/inputs/``,
       materialized once into ``inputs/.materialized/`` and verified against
       ``SHA256SUMS``.  These large inputs are deliberately not pushed; a fresh
       clone needs the parts installed separately or the original local archives.

The materialized layout mirrors the production packages so that path-building
code can treat both identically.
"""
from __future__ import annotations

import hashlib
import lzma
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BENCH_INPUTS = REPO / "benchmarks/hexagon-mhv-ten-loop/inputs"

# component name -> (archive stem, production-relative location)
COMPONENTS = {
    "fec": "fec.tar.xz",
    "lec": "lec.tar.xz",
    "data": "data.tar.xz",
    "early_ope": "early_ope.tar.xz",
    "joint": "joint.tar.xz",
    "solver": "solver.tar.xz",
    "runtime": "runtime.tar.xz",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _parts(stem: str) -> list[Path]:
    parts = sorted(BENCH_INPUTS.glob(stem + ".part.*"))
    if not parts:
        raise SystemExit(f"no local input parts matching {stem}.part.* under {BENCH_INPUTS}; install the split inputs separately")
    return parts


def materialize(force: bool = False) -> Path:
    """Materialize every local split component into inputs/.materialized."""
    cache = BENCH_INPUTS / ".materialized"
    sums = BENCH_INPUTS / "SHA256SUMS"
    if not sums.is_file():
        raise SystemExit(f"missing {sums}")
    expected = {}
    for line in sums.read_text().splitlines():
        digest, name = line.split(None, 1)
        expected[name.strip()] = digest
    cache.mkdir(parents=True, exist_ok=True)
    marker = cache / ".verified"
    if marker.is_file() and not force:
        return cache
    for name in list(COMPONENTS) + ["hexagon_seam_replay"]:
        if name == "hexagon_seam_replay":
            src = BENCH_INPUTS / "bin/hexagon_seam_replay"
            (cache / "bin").mkdir(exist_ok=True)
            if sha256(src) != expected["bin/hexagon_seam_replay"]:
                raise SystemExit(f"sha256 mismatch for committed {src}")
            target = cache / "bin/hexagon_seam_replay"
            target.write_bytes(src.read_bytes())
            target.chmod(0o755)
            continue
        stem = COMPONENTS[name]
        parts = _parts(stem)
        for p in parts:
            if p.name not in expected or sha256(p) != expected[p.name]:
                raise SystemExit(f"sha256 mismatch for committed part {p}")
        out = cache / (stem + ".tar")
        if not out.exists():
            print(f"materializing {stem} from {len(parts)} parts ...", file=sys.stderr)
            with out.open("wb") as f:
                for p in parts:
                    with p.open("rb") as part:
                        f.write(part.read())
        tar = out.with_suffix("")  # drop .xz suffix artifact naming
        import subprocess
        # Decompress xz stream then untar into the component root.
        comp_root = cache / name
        if not comp_root.exists():
            comp_root.mkdir(parents=True)
            with lzma.open(out, "rb") as src, subprocess.Popen(
                    ["tar", "-C", str(comp_root), "-xf", "-"], stdin=subprocess.PIPE) as tar_p:
                while True:
                    block = src.read(1 << 22)
                    if not block:
                        break
                    tar_p.stdin.write(block)
                tar_p.stdin.close()
                if tar_p.wait():
                    raise SystemExit(f"untar failed for {stem}")
        out.unlink(missing_ok=True)
    marker.write_text("verified against SHA256SUMS\n")
    return cache


def roots(mode: str = "auto") -> dict:
    """Return the input roots used by the benchmark runners.

    Keys: fec, lec, data, solver, reference, library, seam_binary,
    local_archives (bool).
    """
    local = {
        "fec": REPO / "output_hexagon_16p4_20261006/FEC",
        "lec": REPO / "output_hexagon_16p4_20261006/LEC/MHV",
        "data": REPO / "output_hexagon_16p4_20261006/data",
        "solver": REPO / "output_hexagon_17p3_20261007/joint-v17-schedule-package/bin/joint_physics_v17",
        "reference": REPO / "output_hexagon_17p3_20261007/ten-v17-main/joint",
        "early_ope": REPO / "output_hexagon_17p3_20261007/ten-v17-main/early_ope_small_s",
        "library": REPO / "output_hexagon_runtime_recovery_20261008/env/lib",
        "seam_binary": REPO / "output_hexagon_validation_20261009/hexagon_seam_replay",
    }
    if mode == "local":
        missing = [str(p) for p in local.values() if not p.exists()]
        if missing:
            raise SystemExit("local archives missing: " + ", ".join(missing))
        return {**local, "local_archives": True}
    if mode in ("bench", "auto"):
        if mode == "auto" and all(p.exists() for p in local.values()):
            return {**local, "local_archives": True}
        cache = materialize()
        # The archives were created with the production directory name as the
        # top-level tar entry (fec.tar.xz -> FEC/, solver.tar.xz ->
        # joint-v17-schedule-package/, runtime.tar.xz -> lib/, ...).
        bench = {
            "fec": cache / "fec/FEC",
            "lec": cache / "lec/LEC/MHV",
            "data": cache / "data/data",
            "solver": cache / "solver/joint-v17-schedule-package/bin/joint_physics_v17",
            "reference": cache / "joint/joint",
            "early_ope": cache / "early_ope/early_ope_small_s",
            "library": cache / "runtime/lib",
            "seam_binary": cache / "bin/hexagon_seam_replay",
        }
        missing = [str(p) for p in bench.values() if not p.exists()]
        if missing:
            raise SystemExit("materialized inputs incomplete: " + ", ".join(missing))
        return {**bench, "local_archives": False}
    raise SystemExit(f"unknown input source mode: {mode}")
