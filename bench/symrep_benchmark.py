#!/usr/bin/env python3
"""Reproducible, sequential exact benchmarks against ordinary and archive chains.

Requires built symrep, the unpacked hexagon irr_fec/irr_lec executables, and
hexagon_symrep_reference compiled against the archive headers (see docs/symrep.md).
All outputs are new directories; no input data or reference program is modified.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import platform
import shutil
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--hexagon", type=pathlib.Path, required=True)
parser.add_argument("--exporter", type=pathlib.Path, default=ROOT / "bench/hexagon_symrep_reference")
parser.add_argument("--solver", type=pathlib.Path, default=ROOT / "symrep")
parser.add_argument("--output", type=pathlib.Path, required=True)
parser.add_argument("--repeats", type=int, default=3)
parser.add_argument("--threads", type=int, default=2)
parser.add_argument("--heptagon-fec-weight", type=int, default=4)
parser.add_argument("--heptagon-lec-weight", type=int, default=3)
args = parser.parse_args()
assert args.repeats > 0 and args.threads > 0
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
archive, exporter = args.hexagon.resolve(), args.exporter.resolve()
program = args.solver.resolve()
snapshot = out / "benchmark-source"
snapshot.mkdir()
for source in [ROOT / "symrep.cpp", *ROOT.glob("symrep*.hpp"), pathlib.Path(__file__).resolve()]:
    shutil.copy2(source, snapshot / source.name)
shutil.copy2(program, snapshot / "symrep")
program = snapshot / "symrep"
env = dict(os.environ, THREADS=str(args.threads))
records = []


def measured(label, command, cwd=ROOT):
    log, usage = out / f"{label}.log", out / f"{label}.usage.json"
    start = time.perf_counter()
    with log.open("w") as stream:
        subprocess.run(["/usr/bin/time", "-f", '{"peak_rss_kib":%M}', "-o", str(usage), *map(str, command)], cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
    record = dict(label=label, wall_s=time.perf_counter() - start, **json.loads(usage.read_text()))
    records.append(record)
    print(f"{label}: {record['wall_s']:.6f} s, {record['peak_rss_kib']} KiB", flush=True)
    return record


def solver(label, *command):
    return measured(label, [program, *command, "--threads", args.threads])


def timings(directory):
    with (directory / "timings.tsv").open() as stream:
        return [{key: (value if key == "corner_columns" else float(value)) for key, value in row.items()} for row in csv.DictReader(stream, delimiter="\t")]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


report = {
    "schema": "symrep-benchmark-v1",
    "machine": dict(platform=platform.platform(), processor=platform.processor(), logical_cpus=os.cpu_count()),
    "threads": args.threads, "repeats": args.repeats,
    "source_base": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "source_sha256": {p.name: sha(p) for p in [ROOT / "symrep.cpp", *ROOT.glob("symrep*.hpp")]},
    "executable_sha256": sha(program),
    "notes": ["Fresh output directories and sequential processes; OS file caches are not flushed.",
              "QQ throughout. Symmetry times include automatic CG construction, assembly and kernel basis reduction.",
              "Preparation and verification are reported separately. Process wall times include startup, WXF I/O and mandatory checks.",
              "Ordinary uses the same executable, SparseRREF checkout, compiler and thread count, retaining the original input basis.",
              "Archive is an independent specialized solver with its own bundled SparseRREF; its wall time is contextual, not an isolated decomposition comparison.",
              "No --compare in timed symmetry runs; exact reference changes of basis verify every weight afterward."],
    "cases": [], "processes": records,
}
cpuinfo = pathlib.Path("/proc/cpuinfo")
if cpuinfo.exists():
    report["machine"]["cpu_model"] = next((l.split(":", 1)[1].strip() for l in cpuinfo.read_text().splitlines() if l.startswith("model name")), "unknown")

# Run the supplied solver itself and expand its multiplicity tensors for exact comparison.
for kind, maximum in [("fec", 4), ("lec", 3)]:
    for repeat in range(args.repeats):
        target = out / f"archive-{kind}-{repeat}"
        command = [archive / "bin" / f"irr_{kind}"] + (["MHV"] if kind == "lec" else []) + [maximum, target]
        measured(f"archive-{kind}-{repeat}", command, archive)
    measured(f"archive-export-{kind}", [exporter, kind.upper(), out / f"archive-{kind}-0", maximum, archive / "data", out / f"hexagon-reference-{kind}"])

heptagon_generators = [("cyclic", ROOT / "data/cycrepmat.wxf"), ("flip", ROOT / "data/fliprepmat.wxf"), ("parity", ROOT / "data/parityrepmat.wxf")]
hexagon_generators = [("cyclic", archive / "data/cycrepmat.wxf"), ("flip", archive / "data/fliprepmat.wxf"), ("parity", out / "hexagon-reference-fec/parity.wxf")]
cases = [
    ("hexagon-fec", "forward", 4, archive / "data/dlogmat_ca_1term_pd3.wxf", out / "hexagon-reference-fec/w1.wxf", hexagon_generators),
    ("hexagon-lec", "backward", 3, archive / "data/dlogmat_full_ca_pd3.wxf", out / "hexagon-reference-lec/w1.wxf", hexagon_generators),
    ("heptagon-fec", "forward", args.heptagon_fec_weight, ROOT / "data/dlogmat_E6.wxf", ROOT / "data/FEC_1.wxf", heptagon_generators),
    ("heptagon-lec", "backward", args.heptagon_lec_weight, ROOT / "data/dlogmat_E6.wxf", ROOT / "data/LEC_1.wxf", heptagon_generators),
]
for name, direction, maximum, condition, seed, generators in cases:
    generator_options = [item for key, path in generators for item in ("--generator", f"{key}={path}")]
    item = dict(name=name, maximum_weight=maximum, inputs={str(p.name): sha(p) for p in [condition, seed, *[p for _, p in generators]]}, runs=[])
    report["cases"].append(item)
    for repeat in range(args.repeats):
        prefix = f"{name}-{repeat}"
        prepared, adapted, ordinary = [out / f"{prefix}-{s}" for s in ("prepared", "adapted", "ordinary")]
        preparation = solver(f"{prefix}-prepare", "prepare", *generator_options, "--condition", condition, "--seed", seed, "--direction", direction, "--output", prepared)
        baseline = solver(f"{prefix}-ordinary", "ordinary", "--condition", condition, "--seed", seed, "--direction", direction, "--max-weight", maximum, "--output", ordinary)
        symmetry = solver(f"{prefix}-adapted", "extend", "--input", prepared, "--max-weight", maximum, "--output", adapted)
        verified = solver(f"{prefix}-verify", "verify", "--input", prepared, "--chain", adapted, "--reference", ordinary, "--max-weight", maximum, "--output", out / f"{prefix}-verified")
        if name.startswith("hexagon"):
            solver(f"{prefix}-archive-verify", "verify", "--input", prepared, "--chain", adapted, "--reference", out / name.replace("hexagon-", "hexagon-reference-"), "--max-weight", maximum, "--output", out / f"{prefix}-archive-verified")
        item["runs"].append(dict(preparation=preparation, ordinary=baseline, symmetry=symmetry, verification=verified, ordinary_weights=timings(ordinary), symmetry_weights=timings(adapted)))
        (out / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    item["median_wall_s"] = {key: statistics.median(r[key]["wall_s"] for r in item["runs"]) for key in ("preparation", "ordinary", "symmetry", "verification")}
    print(name, item["median_wall_s"], flush=True)
    (out / "results.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"Exact verification complete. Results: {out / 'results.json'}")
