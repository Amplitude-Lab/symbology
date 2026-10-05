#!/usr/bin/env python3
"""Reproduce the MHV heptagon bootstrap from repository seeds, without caches.

Linux runner: every numerical stage shares one wall-time/RAM budget. The small
committed solution is a comparison target, never an input to the solve.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "benchmarks/heptagon-mhv-five-loop"
PROGRAMS = ("bootstrap", "compute_rhs", "bench/early_sew_probe",
            "bench/early_mhv_bootstrap", "bench/early_mhv_words")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def check(run, phase, loop):
    """Bind saved coordinates before using/comparing the reference tensor."""
    reference = json.loads((run / "reference/manifest.json").read_text())
    records = []

    def require_hash(relative, expected):
        path = run / relative
        actual = sha(path)
        if actual != expected:
            raise RuntimeError(f"Reference mismatch for {relative}: {actual} != {expected}. "
                               "Do not interpret the archived coefficients in a different basis.")
        records.append(dict(path=relative, sha256=actual, bytes=path.stat().st_size))

    for relative, expected in reference["seed_sha256"].items():
        require_hash(relative, expected)
    if phase in ("bases", "final"):
        for weight in range(2, 2 * loop - 1):
            relative = f"chain/FEC_{weight}.wxf"
            require_hash(relative, reference["basis_sha256"][relative])
        for lower_loop in range(2, loop):
            for name in ("E", "R"):
                relative = f"chain/{lower_loop}loop/{name}{lower_loop}.wxf"
                require_hash(relative, reference["lower_sha256"][relative])
    if phase == "final":
        require_hash("early/LEC_2.wxf", reference["lec2_sha256"])
        if loop == 5:
            expected = reference["solution"]["sha256"]
            require_hash("reference/hepMHV_5L_recursive.wxf", expected)
            require_hash("amplitude/hepMHV_5L_recursive.wxf", expected)
        else:
            require_hash("amplitude/E4_candidate.wxf", reference["lower_sha256"]["chain/4loop/E4.wxf"])
    save(run / (phase + "-binding.json"), dict(phase=phase, loop=loop, state="pass", files=records))
    print(f"PASS {phase} reference binding", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run", help="fresh computation from the tracked seeds")
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--loops", type=int, choices=(4, 5), default=5,
                            help="4 is the smaller end-to-end control; default 5")
    run_parser.add_argument("--hours", type=float, default=24)
    run_parser.add_argument("--ram-gib", type=int, default=50)
    run_parser.add_argument("--min-available-gib", type=float, default=4)
    check_parser = sub.add_parser("check", help="verify a run's recorded coordinate/reference hashes")
    check_parser.add_argument("--output", type=Path, required=True)
    check_parser.add_argument("--loops", type=int, choices=(4, 5), required=True)
    check_parser.add_argument("--phase", choices=("seeds", "bases", "final"), required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    if args.command == "check":
        check(root, args.phase, args.loops)
        return 0
    if args.hours <= 0 or args.ram_gib <= 0 or args.min_available_gib < 0:
        parser.error("positive time/RAM limits and a nonnegative available-memory floor are required")
    if not sys.platform.startswith("linux"):
        parser.error("the measured process controller requires Linux; use Linux or WSL")
    for name in PROGRAMS:
        if not os.access(ROOT / name, os.X_OK):
            parser.error(f"missing executable {name}; run make five-loop-tools first")
    reference = json.loads((RELEASE / "manifest.json").read_text())
    for relative, expected in reference["seed_sha256"].items():
        if sha(ROOT / relative) != expected:
            parser.error(f"the benchmark seed differs from the recorded input: {relative}")
    if sha(RELEASE / "reference/hepMHV_5L_recursive.wxf") != reference["solution"]["sha256"]:
        parser.error("committed reference solution checksum mismatch")

    root.mkdir(parents=True, exist_ok=False)
    for directory in ("bin", "data", "chain", "early/actions", "amplitude", "reference"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    for name in PROGRAMS:
        shutil.copy2(ROOT / name, root / "bin" / Path(name).name)
    for name in ("run_pipeline.py", "measure_process.py", "five_loop_heptagon.py"):
        shutil.copy2(ROOT / "bench" / name, root / "bin" / name)
    for relative in reference["seed_sha256"]:
        shutil.copy2(ROOT / relative, root / relative)
    shutil.copy2(RELEASE / "manifest.json", root / "reference/manifest.json")
    shutil.copy2(RELEASE / "reference/hepMHV_5L_recursive.wxf", root / "reference/hepMHV_5L_recursive.wxf")

    def git_output(*arguments):
        result = subprocess.run(["git", *arguments], cwd=ROOT, text=True, capture_output=True)
        return result.stdout.strip() if not result.returncode else None

    source_files = [*ROOT.glob("*.hpp"), *ROOT.glob("*.h"),
                    *ROOT.glob("SparseRREF/*.h"), *ROOT.glob("SparseRREF/*.hpp"),
                    *ROOT.glob("patches/*.patch"), ROOT / "Makefile",
                    ROOT / "bootstrap.cpp", ROOT / "compute_rhs.cpp",
                    *[ROOT / (name + ".cpp") for name in PROGRAMS if name.startswith("bench/")],
                    ROOT / "bench/five_loop_heptagon.py", ROOT / "bench/run_pipeline.py",
                    ROOT / "bench/measure_process.py", ROOT / "scripts/setup-sparserref.sh"]
    save(root / "manifest.json", dict(
        schema="heptagon-mhv-cold-benchmark-v1", loops=args.loops, threads=8,
        git_head=git_output("rev-parse", "HEAD"), git_status=git_output("status", "--short"),
        started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_files},
        frozen_sha256={str(p.relative_to(root)):sha(p) for p in root.rglob("*") if p.is_file()},
        meminfo_at_start=Path("/proc/meminfo").read_text(),
        cpuinfo=Path("/proc/cpuinfo").read_text(),
        input_policy="Only tracked seed data; all bases, actions and lower-loop amplitudes are regenerated. "
                     "The archived solution is used only for the final comparison."))

    binary, data, chain = root / "bin", root / "data", root / "chain"
    early, amplitude = root / "early", root / "amplitude"
    loop = args.loops
    forward = 2 * loop - 2
    stages = []

    def stage(name, command, outputs=(), markers=()):
        stages.append(dict(name=name, command=[str(x) for x in command],
                           required_outputs=[str(x) for x in outputs], success_markers=list(markers)))

    def binding(phase):
        stage(phase + "-binding", [sys.executable, binary / "five_loop_heptagon.py", "check",
                                  "--output", root, "--loops", loop, "--phase", phase],
              [root / (phase + "-binding.json")], [f"PASS {phase} reference binding"])

    binding("seeds")
    lower_outputs = [chain / f"{l}loop" / f"{name}{l}.wxf" for l in range(2, loop) for name in ("E", "R")]
    stage("lower-loops", [binary / "compute_rhs", "--target", f"SEW_{forward-1}p1",
                          "--data-dir", data, "--output-dir", chain,
                          "--letter-projection", chain / "collinear/colprojdiv_w1.wxf",
                          "--threads", 8, "--kernel-strategy", "streamed",
                          "--sew-strategy", "auto", "--projection-strategy", "restricted"], lower_outputs)
    stage("forward-basis", [binary / "bootstrap", "--extend", "-c", data / "dlogmat_E6.wxf",
                            "-f", chain / f"FEC_{forward-1}.wxf", "-o", chain / f"FEC_{forward}.wxf",
                            "--threads", 8, "--kernel-strategy", "streamed"], [chain / f"FEC_{forward}.wxf"])
    binding("bases")
    stage("actions", [binary / "early_sew_probe", "actions", data, chain, forward, early / "actions"],
          [early / "actions" / f"w{forward}_{g}.wxf" for g in ("cyc", "flip", "parity")])
    stage("boundary", [binary / "early_mhv_bootstrap", "boundary", data, chain, early, chain, amplitude, loop],
          [amplitude / f"boundary_{loop}L.wxf"], ["PACKED_BOUNDARY_RESULT"])
    stage("kernel", [binary / "early_sew_probe", "refine-hybrid", data, chain, forward, early],
          [early / f"EARLY_{forward}p2.wxf", early / "LEC_2.wxf"],
          ["EARLY_RESULT dimension=", "exact_certificate=complete_factored_residues_plus_height"])
    stage("coefficients", [binary / "early_mhv_bootstrap", "coefficients", data, chain, early, chain, amplitude, loop],
          [amplitude / "candidate_coefficients.wxf", amplitude / "candidate_recursive.wxf"],
          ["unique=1", "CANDIDATE saved"])
    stage("verify", [binary / "early_mhv_bootstrap", "verify", data, chain, early, chain, amplitude, loop],
          [amplitude / f"hepMHV_{loop}L_recursive.wxf", amplitude / f"E{loop}_candidate.wxf"],
          ["PASS complete collinear divergent part", f"CERTIFIED unique {loop}-loop MHV symbol"])
    stage("words", [binary / "early_mhv_words", data, chain, early / "LEC_2.wxf",
                    amplitude / f"hepMHV_{loop}L_recursive.wxf", loop],
          markers=["PASS four independent published original-alphabet word coefficients"])
    binding("final")
    config = dict(run_directory=str(root), working_directory=str(binary),
                  measure_program=str(binary / "measure_process.py"),
                  deadline_unix_seconds=time.time() + args.hours * 3600,
                  virtual_gib=args.ram_gib, rss_gib=args.ram_gib,
                  min_available_gib=args.min_available_gib, stages=stages)
    save(root / "config.json", config)
    print(f"Cold {loop}-loop benchmark: {root}\nProgress: {root / 'status.json'}", flush=True)
    # Replace this process so cancellation reaches the shared-budget controller.
    os.execv(sys.executable, [sys.executable, str(binary / "run_pipeline.py"), str(root / "config.json")])


if __name__ == "__main__":
    raise SystemExit(main())
