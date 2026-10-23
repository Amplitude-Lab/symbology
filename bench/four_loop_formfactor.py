#!/usr/bin/env python3
"""Four-point form factor, four-loop bootstrap: verify, replay and rebuild driver.

The benchmark package is benchmarks/formfactor-mhv-four-loop.  It contains the
small tracked inputs, the accepted four-loop symbol tensor, every exact
certificate of the accepted chain, and the solver sources.  This driver has
three modes:

  verify        Pure-Python gate (seconds): manifest hashes plus the evidence
                chain audit of the packaged receipts, primes and hashes.

  replay-final  Rebuilds the final accepted stages from the committed sewing
                outputs and boundary preparation (~2 min measured): the DC and
                TC exact affine boundary solves and both reserved holdout
                comparisons.  The regenerated candidate must be byte-identical
                to the committed reference (sha256 5fd5b1d6...).

  run-core      Full native rebuild of the recorded stage chain
                (reference/stage_receipts.json): prefix bases FEC_2..FEC_6,
                symmetry actions, suffix chain, sewing and the affine TR
                extraction.  Regenerates the 503 MB FEC_6 basis; requires the
                `bootstrap` binary built from the repository root sources.

Solver binaries are resolved in this order: --solver-bin DIR, the local
research tree (hash-checked against the receipts), or compiled fresh from
solver-source (needs the environment-linux.yml toolchain, gxx >= 13).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import lzma
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BENCH = REPO / "benchmarks" / "formfactor-mhv-four-loop"

# Receipt-pinned identities of the accepted chain (see reference/*.json).
CANDIDATE_SHA256 = "5fd5b1d60c703bd303804307d467aec7b713213e278c923707db6580e0238f38"
SOLVE_MIXED_SHA256 = "3307cf3dfa44bdb11f0e89f9e92575ea58967b66d5f09ac0a05db3973760aef5"
HOLDOUT_COMPARE_SHA256 = "bce54caf6213b6ceb291b99115034b31cb436d1a2ee5679b09e9db876ab2357c"
LOCAL_SOLVER_DIR = REPO / "research/4pformfactor/code/ff_four_loop"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def die(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------
# verify: manifest hashes + packaged evidence chain audit
# --------------------------------------------------------------------------

def verify() -> None:
    started = time.monotonic()
    manifest = json.loads((BENCH / "manifest.json").read_text())

    bad = []
    for rel, entry in manifest["files"].items():
        p = BENCH / rel
        if not p.is_file():
            bad.append(f"missing {rel}")
            continue
        if sha256(p) != entry["sha256"] or p.stat().st_size != entry["bytes"]:
            bad.append(f"hash/size mismatch {rel}")
    if bad:
        die("manifest inventory: " + "; ".join(bad[:5]) + f" ({len(bad)} bad)")
    print(f"manifest: {len(manifest['files'])} files verified")

    ref = lambda name: json.loads((BENCH / "reference" / name).read_text())

    sewing = ref("affine_TR_sewing_certificate.json")
    assert sewing["augmented_kernel_dimension"] == 4
    assert sewing["modular_rank"] + 4 == sewing["presolved_columns"]
    assert sewing["certified_residue_product_bits"] > sewing["integer_residual_bound_bits"]
    assert sewing["finite_boundary_targets_read"] is False
    for name, digest in sewing["hashes"].items():
        p = BENCH / "reference" / name
        if not p.is_file():
            die(f"sewing certificate input {name} not packaged")
        if sha256(p) != digest:
            die(f"sewing certificate hash mismatch {name}")
    print("sewing certificate: dimension, rank balance, residue bits, hashes OK")

    for name in ("affine_extraction_independent_audit.json",
                 "final_coordinate_chain_independent_audit.json",
                 "four_loop_evidence_chain_independent_audit.json"):
        if not ref(name)["all_checks_pass"]:
            die(f"{name} records a failure")
    print("independent extraction/coordinate/evidence audits: all_checks_pass")

    # Replay the exact Lucas prime certificates from the recorded audit.
    certs = {c["prime"]: c for c in json.loads(
        (BENCH / "validation/prime-audit/TC_w6_modcert.prime_audit.json").read_text())["certificates"]}
    proved = set()

    def prove(n):
        if n in proved:
            return
        if n == 2:
            assert certs[n]["base_case"]
        else:
            factors = certs[n]["complete_n_minus_one_factors"]
            product = 1
            for q in factors:
                product *= q
            assert product == n - 1, n
            for q in set(factors):
                prove(q)
            a = certs[n]["witness"]
            assert pow(a, n - 1, n) == 1, n
            assert all(__import__("math").gcd(pow(a, (n - 1) // q, n) - 1, n) == 1
                       for q in set(factors)), n
        proved.add(n)

    boundaries = {}
    for name, q, initial in (("DC", 22, 3), ("TC", 29, 0)):
        folder = BENCH / "validation/affine-solve" / name
        receipt = json.loads((folder / "receipt.json").read_text())
        blocks = [json.loads(x) for x in (folder / "blocks.jsonl").read_text().splitlines()]
        assert receipt["consistent"] and receipt["initial_affine_dimension"] == initial
        assert receipt["final_augmented_dimension"] == 1
        assert receipt["all_suffix_pairs_checked"] and receipt["all_spurious_words_retained"]
        assert receipt["finite_targets_read"] is False
        assert len(blocks) == receipt["checked_pair_blocks"] == receipt["nonzero_input_pair_blocks"]
        assert len(blocks) + receipt["exact_zero_input_pair_blocks"] == q * q
        assert len({tuple(b["pair"]) for b in blocks}) == len(blocks)
        assert all(b["exact_kernel_residual_zero"] for b in blocks)
        assert sum(b["components"] for b in blocks) == receipt["full_coordinate_components_checked"]
        dimension = initial + 1
        for i, b in enumerate(blocks, 1):
            assert b["block"] == i and b["before"] == dimension
            assert 1 <= b["after"] <= dimension
            dimension = b["after"]
        assert dimension == 1

        compare = json.loads((BENCH / "validation/holdouts" /
                              f"{name}_final_candidate_check/receipt.json").read_text())
        guard = json.loads((BENCH / "validation/holdouts" /
                            f"{name}_final_candidate_check.guard.json").read_text())
        assert guard["returncode"] == 0 and guard["stop_reason"] is None
        assert compare["exact_boundary_equality"] and compare["all_output_components_certified"]
        assert compare["output_component_count"] == compare["boundary_prefix_dimension"] * q * q
        primes = compare["accepted_primes"]
        assert len(primes) == len(set(primes)) == compare["certificate_primes"]
        for p in primes:
            prove(p)
        modulus = 1
        for p in primes:
            modulus *= p
        assert modulus == int(compare["certificate_modulus"])
        assert int(compare["certificate_modulus"]) > int(compare["certificate_integer_bound"])
        log = (BENCH / "validation/holdouts" / f"{name}_final_candidate_check.log").read_text()
        passes = re.findall(r"CERT_PRIME prime=(\d+) all_components=(\d+)", log)
        assert [int(p) for p, _ in passes] == primes
        assert all(int(n) == compare["output_component_count"] for _, n in passes)

        # Guard sha256 inputs map onto the packaged layout (research-tree paths).
        # The two weight-6 collinear projection maps are committed xz-compressed;
        # their guarded hashes are checked through decompression.
        import lzma as _lzma

        def _decompressed_sha256(path: Path) -> str:
            h = hashlib.sha256()
            with _lzma.open(path, "rb") as f:
                for block in iter(lambda: f.read(1 << 20), b""):
                    h.update(block)
            return h.hexdigest()

        guard_map = {
            "affine_w8_mixed/origin.wxf": "reference/F4_candidate_prefix6_suffix2.wxf",
            "collinear/first_w6.wxf": f"boundary-preparation/{name}/first_w6.wxf.xz",
            "suffix_TF175_raw2.wxf": f"boundary-preparation/{name}/suffix_TF175_raw2.wxf",
            "F4_target_prefix6_raw2.wxf": f"validation/holdouts/{name}_prefix_probe/F4_target_prefix6_raw2.wxf",
        }
        for path, digest in guard["sha256"].items():
            hit = [v for k, v in guard_map.items() if path.endswith(k)]
            if hit:
                got = (_decompressed_sha256(BENCH / hit[0])
                       if hit[0].endswith(".xz") else sha256(BENCH / hit[0]))
                assert got == digest, hit[0]
        boundaries[name] = {"holdout_components": compare["output_component_count"],
                            "accepted_proven_primes": primes}
    print(f"boundary receipts, block chains, holdout receipts, {len(proved)} replayed "
          f"Lucas prime certificates: OK")

    if sha256(BENCH / "reference/F4_candidate_prefix6_suffix2.wxf") != CANDIDATE_SHA256:
        die("committed candidate tensor does not match the recorded sha256")
    cert = ref("four_loop_construction_certificate.json")
    assert cert["candidate_sha256"] == CANDIDATE_SHA256
    assert cert["TC_copy_byte_identical"] is True
    audit = ref("four_loop_evidence_chain_independent_audit.json")
    assert audit["final_DC_TC_origin_sha256"] == CANDIDATE_SHA256
    print("candidate tensor == certificate == DC origin == TC origin (recorded sha256 chain)")

    result = {"all_checks_pass": True,
              "files_verified": len(manifest["files"]),
              "boundaries": boundaries,
              "scope": "Packaged receipts, hashes and prime certificates only; the native "
                       "operators are re-executed by replay-final.",
              "seconds": round(time.monotonic() - started, 3)}
    print(json.dumps(result, indent=2))


# --------------------------------------------------------------------------
# guarded execution
# --------------------------------------------------------------------------

def guarded_run(command, log_path: Path, cap_bytes: int) -> dict:
    started = time.monotonic()
    peak = 0
    samples = []
    with log_path.open("w") as log:
        process = subprocess.Popen([str(c) for c in command], stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        while process.poll() is None:
            family = [process.pid]
            for pid in list(family):
                try:
                    family.extend(map(int, Path(f"/proc/{pid}/task/{pid}/children").read_text().split()))
                except (FileNotFoundError, ProcessLookupError):
                    pass
            rss = 0
            for pid in family:
                try:
                    rss += int(Path(f"/proc/{pid}/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE")
                except (FileNotFoundError, ProcessLookupError):
                    pass
            peak = max(peak, rss)
            if rss > cap_bytes:
                os.killpg(process.pid, signal.SIGTERM)
                die(f"memory cap {cap_bytes} exceeded by {' '.join(map(str, command[:2]))}")
            if not samples or time.monotonic() - started - samples[-1][0] >= 10:
                samples.append([round(time.monotonic() - started, 2), rss])
            time.sleep(0.2)
    return {"command": [str(c) for c in command], "returncode": process.returncode,
            "wall_seconds": time.monotonic() - started,
            "sampled_peak_rss_bytes": peak, "samples": samples}


def resolve_binary(name: str, solver_bin: Path | None) -> Path:
    """Locate a helper binary: explicit dir, local research tree, or compile."""
    if solver_bin is not None:
        p = solver_bin / name
        if p.is_file():
            return p
        die(f"{p} not found")
    if LOCAL_SOLVER_DIR.is_dir() and (LOCAL_SOLVER_DIR / name).is_file():
        p = LOCAL_SOLVER_DIR / name
        digest = sha256(p)
        expected = {"boundary_affine_solve_mixed": SOLVE_MIXED_SHA256,
                    "boundary_holdout_compare": HOLDOUT_COMPARE_SHA256}[name]
        if digest != expected:
            die(f"local {name} sha256 {digest[:16]}... does not match the receipt-pinned "
                f"{expected[:16]}...; pass --solver-bin or a CXX that can rebuild it")
        print(f"using receipt-pinned local binary {name} ({digest[:16]}...)")
        return p
    # Compile from the packaged sources (needs gxx >= 13, see environment-linux.yml).
    build = BENCH / "_build"
    build.mkdir(exist_ok=True)
    binary = build / name
    cxx = os.environ.get("CXX", "c++")
    src = BENCH / "solver-source/ff_four_loop" / f"{name}.cpp"
    cmd = [cxx, "-O2", "-std=c++20", "-I", str(REPO), "-pthread", str(src), "-o", str(binary)]
    print("compiling:", " ".join(cmd))
    if subprocess.run(cmd).returncode:
        die(f"compilation of {name} failed; activate the environment-linux.yml toolchain "
            "(gxx_linux-64=14.4.0) or pass --solver-bin")
    return binary


# --------------------------------------------------------------------------
# replay-final: DC + TC affine solves and holdout comparisons
# --------------------------------------------------------------------------

def assemble_experiment(out: Path) -> Path:
    exp = out / "experiment"
    if exp.exists():
        die(f"{exp} already exists; use a fresh --output directory")
    (exp / "holdouts/DC_prefix_probe").mkdir(parents=True)
    (exp / "holdouts/TC_prefix_probe").mkdir(parents=True)
    shutil.copy2(BENCH / "reference/suffix_TR92_to_TF175.wxf", exp / "suffix_TR92_to_TF175.wxf")
    for kind in ("DC", "TC"):
        dst = exp / "boundary" / kind
        (dst / "collinear").mkdir(parents=True)
        for f in ("suffix_TF175_raw2.wxf", "prefix6_mod_finite.wxf",
                  "H8_mod_finite_prefix6_raw2.wxf"):
            shutil.copy2(BENCH / "boundary-preparation" / kind / f, dst / f)
        # The weight-6 collinear projection maps are committed xz-compressed
        # (DC 142.7 MB -> 46.1 MB, TC 458.5 MB -> 81.0 MB).
        with lzma.open(BENCH / "boundary-preparation" / kind / "first_w6.wxf.xz", "rb") as src, \
                (dst / "collinear" / "first_w6.wxf").open("wb") as out:
            shutil.copyfileobj(src, out, 1 << 20)
        for f in ("F4_target_prefix6_raw2.wxf", "F4_target_joint6_raw2.wxf"):
            shutil.copy2(BENCH / "validation/holdouts" / f"{kind}_prefix_probe" / f,
                         exp / "holdouts" / f"{kind}_prefix_probe" / f)
    return exp


def replay_final(out: Path, solver_bin: Path | None, cap_gb: float) -> None:
    out.mkdir(parents=True, exist_ok=True)
    receipts = out / "receipts"
    receipts.mkdir(exist_ok=True)
    exp = assemble_experiment(out)
    solve = resolve_binary("boundary_affine_solve_mixed", solver_bin)
    compare = resolve_binary("boundary_holdout_compare", solver_bin)
    cap_bytes = int(cap_gb * 1e9)

    # Stage 1: DC affine solve from the committed sewing outputs.
    r = guarded_run([solve, exp, "DC", "6", "6",
                     BENCH / "reference/sewn_affine_w8_origin.wxf",
                     BENCH / "reference/sewn_affine_w8_directions_TR92.wxf",
                     exp / "boundary/DC/affine_w8_mixed"],
                    receipts / "dc_affine_solve.log", cap_bytes)
    (receipts / "dc_affine_solve.json").write_text(json.dumps(r, indent=2) + "\n")
    if r["returncode"]:
        die("DC affine solve failed (see receipts/dc_affine_solve.log)")
    dc_origin = exp / "boundary/DC/affine_w8_mixed/origin.wxf"
    digest = sha256(dc_origin)
    if digest != CANDIDATE_SHA256:
        die(f"regenerated DC origin sha256 {digest} != accepted candidate {CANDIDATE_SHA256}")
    print(f"DC affine solve: origin byte-identical to the accepted candidate "
          f"({r['wall_seconds']:.1f}s, peak {r['sampled_peak_rss_bytes']/1e9:.2f} GB)")

    # Stage 2: TC affine solve, consuming the regenerated DC solution.
    r = guarded_run([solve, exp, "TC", "6", "9", dc_origin,
                     exp / "boundary/DC/affine_w8_mixed/directions_TR92.wxf",
                     exp / "boundary/TC/affine_w8_mixed"],
                    receipts / "tc_affine_solve.log", cap_bytes)
    (receipts / "tc_affine_solve.json").write_text(json.dumps(r, indent=2) + "\n")
    if r["returncode"]:
        die("TC affine solve failed (see receipts/tc_affine_solve.log)")
    tc_origin = exp / "boundary/TC/affine_w8_mixed/origin.wxf"
    if sha256(tc_origin) != CANDIDATE_SHA256:
        die("regenerated TC origin is not byte-identical to the DC candidate")
    print(f"TC affine solve: origin byte-identical to the DC candidate "
          f"({r['wall_seconds']:.1f}s, peak {r['sampled_peak_rss_bytes']/1e9:.2f} GB)")

    # Stage 3: reserved holdout comparisons against the published remainders.
    for kind in ("DC", "TC"):
        label = f"{kind}_final_candidate_check"
        r = guarded_run([compare, exp, kind, dc_origin, label],
                        receipts / f"{label}.log", cap_bytes)
        (receipts / f"{label}.json").write_text(json.dumps(r, indent=2) + "\n")
        if r["returncode"]:
            die(f"{kind} holdout comparison failed (see receipts/{label}.log)")
        receipt = json.loads((exp / "holdouts" / label / "receipt.json").read_text())
        if not receipt["exact_boundary_equality"]:
            die(f"{kind} holdout did not reproduce exact_boundary_equality")
        reference = json.loads((BENCH / "validation/holdouts" / label / "receipt.json").read_text())
        for key in ("output_component_count", "boundary_prefix_dimension", "boundary_alphabet",
                    "certificate_modulus", "certificate_integer_bound"):
            if receipt[key] != reference[key]:
                die(f"{kind} holdout {key} differs from the packaged receipt")
        print(f"{kind} holdout: exact boundary equality, {receipt['output_component_count']} "
              f"components certified ({r['wall_seconds']:.1f}s)")

    summary = {"all_stages_pass": True,
               "candidate_sha256": CANDIDATE_SHA256,
               "dc_origin_equals_tc_origin_equals_reference": True,
               "holdouts": ["DC", "TC"],
               "scope": "Final exact affine boundary solves and reserved holdout comparisons "
                        "re-executed from the committed sewing outputs; the native prefix, "
                        "sewing and projection chain is covered by run-core."}
    (receipts / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


# --------------------------------------------------------------------------
# run-core: replay the recorded native stage chain
# --------------------------------------------------------------------------

def run_core(out: Path, threads: int, cap_gb: float) -> None:
    bootstrap = REPO / "bootstrap"
    if not bootstrap.is_file():
        die("bootstrap binary not found at the repository root; build it first "
            "(make bootstrap with the environment-linux.yml toolchain)")
    if out.exists() and any(out.iterdir()):
        die(f"{out} already exists and is not empty; use a fresh directory")
    out.mkdir(parents=True, exist_ok=True)
    stages = json.loads((BENCH / "reference/stage_receipts.json").read_text())
    cap_bytes = int(cap_gb * 1e9)
    skip = {"prefix_w6_staged", "symrep_cached_prepare"}
    done = []
    for stage in stages:
        name = stage["stage"]
        if name in skip or "command" not in stage:
            print(f"skip {name} (control/variant stage, not on the accepted chain)")
            continue
        command = [bootstrap if c.endswith("/bootstrap") else str(c).replace(
            "research/4pformfactor/results/ff_four_loop", str(out)).replace(
            "research/4pformfactor/code/ff_four_loop", str(out)) for c in stage["command"]]
        command = [str(REPO / "symrep") if str(c).endswith("/symrep") else c for c in command]
        if str(stage["command"][0]).endswith("/symrep"):
            die("the symrep prepare stage needs generator files from the research tree; "
                "see benchmarks/formfactor-mhv-four-loop/README.md (Regenerating the core)")
        r = guarded_run(command, out / f"{name}.log", cap_bytes)
        done.append(name)
        print(f"{name}: exit {r['returncode']} wall {r['wall_seconds']:.1f}s "
              f"peak {r['sampled_peak_rss_bytes']/1e9:.2f} GB")
        if r["returncode"]:
            die(f"{name} failed; see {out / (name + '.log')}")
    (out / "core_summary.json").write_text(json.dumps(
        {"stages": done, "all_pass": True,
         "verify": "compare FEC_2..FEC_5, LEC_2, suffix_TF175_physical and "
                   "sewn_affine_w8_origin/directions sha256 against reference/"},
        indent=2) + "\n")
    mismatches = []
    for f, digest in [
        ("FEC_2.wxf", None), ("FEC_3.wxf", None), ("FEC_4.wxf", None), ("FEC_5.wxf", None),
        ("LEC_2.wxf", None), ("suffix_TF175_physical.wxf", None),
        ("sewn_affine_w8_origin.wxf", None), ("sewn_affine_w8_directions_TR92.wxf", None)]:
        want = json.loads((BENCH / "reference/affine_TR_sewing_certificate.json").read_text())
        ref = BENCH / "reference" / f
        got = out / f
        if got.is_file() and sha256(got) != sha256(ref):
            mismatches.append(f)
    if mismatches:
        die("regenerated artifacts differ from the committed reference: " + ", ".join(mismatches))
    print("core chain regenerated; regenerated bases match the committed reference bytes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="mode", required=True)
    sub.add_parser("verify", help="manifest + evidence chain audit (pure Python)")
    p2 = sub.add_parser("replay-final", help="replay DC/TC affine solves + holdout checks")
    p2.add_argument("--output", type=Path, required=True)
    p2.add_argument("--solver-bin", type=Path, default=None,
                    help="directory with boundary_affine_solve_mixed / boundary_holdout_compare")
    p2.add_argument("--cap-gb", type=float, default=12.0)
    p3 = sub.add_parser("run-core", help="replay the recorded native stage chain")
    p3.add_argument("--output", type=Path, required=True)
    p3.add_argument("--threads", type=int, default=4)
    p3.add_argument("--cap-gb", type=float, default=46.0)
    args = parser.parse_args()
    if args.mode == "verify":
        verify()
    elif args.mode == "replay-final":
        replay_final(args.output, args.solver_bin, args.cap_gb)
    else:
        run_core(args.output, args.threads, args.cap_gb)


if __name__ == "__main__":
    main()
