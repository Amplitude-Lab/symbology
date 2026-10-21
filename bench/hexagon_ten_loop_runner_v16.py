#!/usr/bin/env python3
"""Sequential V16 conditional, interleaved-physics 17+3 caches and 16+4 reconstruction.

Run this supervisor in an independently managed systemd user service. A rational
reconstruction is followed by an independent exact sewing certificate and
44 reserved origin checks. Unrestricted uniqueness is not asserted.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import struct
import subprocess
import sys
import time
import uuid

REPO = Path(__file__).resolve().parents[1]
EXPERIMENT = REPO / "output_hexagon_17p3_20261007"
SHARED = REPO / "output_hexagon_16p4_20261006"
NATIVE = EXPERIMENT / "ten-v16-main"
SHIFTED = EXPERIMENT / "ten-v16-shifted-16p4"
CONTROL = EXPERIMENT / "ten-v16-runner"
DRIVER = EXPERIMENT / "joint-v16-mrk-reuse-package/bin/joint_physics_v16"
CONVERTER = SHARED / "joint-v8-package/bin/irr_convert"
FEC, LEC, DATA = SHARED / "FEC", SHARED / "LEC/MHV", SHARED / "data"
MEASURE = REPO / "bench/measure_process.py"
CERTIFICATE = REPO / "output_hexagon_certificate_provenance_20261007/bin/hexagon_shifted_joint_certificate"
ORIGIN_CHECK = SHARED / "joint-v8-package/bin/origin_check"
ORIGIN_CASES = REPO / "audits/hexagon-ten-loop-2026-10-06/origin-cases.tsv"
RSS_GIB = 50_000_000_000 / 1024**3
PHYSICAL_STRATEGY = "4 light rounds; at most one COL/MRK cycle per callback; OPE after two cycles; callbacks between refinements and after 120 seconds within a sweep"
ASSUMPTIONS = "MHV double-final bonus relations, Dixon-Liu 2308.08199 eq76-79, full D3 orbit; extrapolated to L10"
HEADER = struct.Struct("<8s8Q")
active = None
cancel_signal = 0


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def tsv(path):
    return dict(line.split("\t", 1) for line in path.read_text().splitlines())


def check_cache(root, prime, h, k, columns):
    path = root / "cache" / f"p{prime}.bin"
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Missing ordinary cache: {path}")
    digest = sha(path)
    if path.with_name(path.name + ".sha256").read_text().strip() != digest:
        raise RuntimeError(f"Cache checksum mismatch: {path}")
    with path.open("rb") as stream:
        magic, version, loop, head, tail, rows, cols, p, identity = HEADER.unpack(stream.read(HEADER.size))
    expected_id = int(sha(root / "input_identity.txt")[:16], 16)
    if (magic, version, loop, head, tail, cols, p, identity) != (
            b"NMHVSEW\0", 1, 10, h, k, columns, prime, expected_id):
        raise RuntimeError(f"Cache header/identity mismatch: {path}")
    if path.stat().st_size != HEADER.size + 8 * rows * cols:
        raise RuntimeError(f"Cache length mismatch: {path}")
    if rows != 1:
        raise RuntimeError(f"Native space still has {rows - 1} directions; more verified physical constraints are required")
    if root == NATIVE / "joint":
        stage = tsv(root / "stage.tsv")
        if (stage.get("algorithm"), stage.get("early_physics"), stage.get("target")) != (
                "joint-discovery-physical-affine-v16-shared-mrk-shuffles", PHYSICAL_STRATEGY, "9000"):
            raise RuntimeError("Native V16 physical-strategy identity mismatch")
    if root == NATIVE / "joint" and stage.get("empirical_assumptions") != ASSUMPTIONS:
        raise RuntimeError("Native conditional-assumption identity mismatch")
    if root == NATIVE / "joint" and stage.get("ope_storage") != "triangular_log_powers_a_plus_b_le_weight":
        raise RuntimeError("Native OPE storage identity mismatch")
    if root == NATIVE / "joint" and (stage.get("mrk_rhs_engine"), stage.get("mrk_rhs_prefetch")) != (
            "shared_one_loop_log_shuffle_v1", "4_for_single_sample;bounded_cache_128"):
        raise RuntimeError("Native MRK generator identity mismatch")
    return {"path": str(path), "sha256": digest, "prime": prime,
            "rows": rows, "directions": rows - 1, "columns": cols}


def check_qq():
    receipt = SHIFTED / "qq_complete.tsv"
    if not receipt.exists():
        return False
    fields = tsv(receipt)
    if fields.get("schema") != "irr-physics-qq-v1" or fields.get("identity") != sha(SHIFTED / "input_identity.txt"):
        raise RuntimeError("Shifted QQ receipt identity mismatch")
    required = {"object.tsv", "columns.tsv", "manifest.json", "stage.tsv"}
    if not required <= fields.keys():
        raise RuntimeError("Shifted QQ receipt is incomplete")
    for name, digest in fields.items():
        if name in ("schema", "identity"):
            continue
        if Path(name).name != name or name in (".", ".."):
            raise RuntimeError("Unsafe QQ receipt filename")
        path = SHIFTED / name
        if path.is_symlink() or sha(path) != digest:
            raise RuntimeError(f"Shifted QQ checksum mismatch: {name}")
    manifest = json.loads((SHIFTED / "manifest.json").read_text())
    obj = tsv(SHIFTED / "object.tsv")
    if (manifest.get("complete") is not True or manifest.get("validation_primes") != 2 or
            (obj.get("kind"), obj.get("directions"), obj.get("h"), obj.get("k")) != ("unique", "0", "16", "4")):
        raise RuntimeError("Shifted QQ completion metadata mismatch")
    if sum(name.endswith("_origin.wxf") or name.endswith("_origin.wxf.gz") for name in fields) != 6:
        raise RuntimeError("QQ receipt must cover all six origin blocks")
    return True


def request_stop(signum, _frame):
    global cancel_signal
    cancel_signal = signum
    if active is not None and active.poll() is None:
        active.send_signal(signal.SIGTERM)


def main():
    global active
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-primes", type=int, default=12)
    parser.add_argument("--preflight-only", action="store_true", help="Verify completed controls and required artifacts without starting a calculation")
    args = parser.parse_args()
    if not 1 <= args.max_primes <= 64:
        parser.error("--max-primes must be 1..64")
    controls = json.loads((EXPERIMENT / "v13-controls/progress.json").read_text())
    if not controls.get("status", "").startswith("controls_passed") or not all(
            controls.get(name, {}).get("every_coefficient_equal") for name in ("eight_comparison", "nine_comparison")):
        raise RuntimeError("V13 known eight- and nine-loop controls must pass before ten-loop execution")
    if tsv(EXPERIMENT / "bonus-double-audit-repair/audit_complete.tsv").get("complete") != "true":
        raise RuntimeError("Conditional tail audit must be complete")
    for name in ("v14-native-probe.json", "v14-eight.json"):
        proof = json.loads((EXPERIMENT / "audit" / name).read_text())
        if proof.get("exit_code") != 0:
            raise RuntimeError("V14 engineering and known-symbol controls must pass first")
    for name in ("v14-eight-comparison.json", "v14-native-comparison.json"):
        proof = json.loads((EXPERIMENT / "audit" / name).read_text())
        if proof.get("exact_payload_equal") is not True:
            raise RuntimeError("V14 exact coefficient comparisons must pass first")
    for name in ("v15-build.json", "v15-triangular-probe.json", "v15-eight.json"):
        if json.loads((EXPERIMENT / "audit" / name).read_text()).get("exit_code") != 0:
            raise RuntimeError("V15 build, triangular primitive and known-symbol controls must pass first")
    if json.loads((EXPERIMENT / "audit/v15-eight-comparison.json").read_text()).get("exact_payload_equal") is not True:
        raise RuntimeError("V15 exact known-symbol comparison must pass first")
    rhs_controls = json.loads((EXPERIMENT / "mrk-rhs-reuse/control/progress.json").read_text())
    integration = json.loads((EXPERIMENT / "v16-controls/progress.json").read_text())
    if rhs_controls.get("status") != "all_coefficients_equal" or len(rhs_controls.get("cases", [])) != 32:
        raise RuntimeError("Exact MRK generator controls must pass first")
    if integration.get("status") != "controls_passed" or not all(integration.get(name, {}).get("every_coefficient_equal")
            for name in ("v15-eight-t8", "v16-eight-t8", "v16-eight-t16")):
        raise RuntimeError("V16 integrated known-symbol controls must pass first")
    cases = [line.split() for line in ORIGIN_CASES.read_text().splitlines() if line.strip() and line.split()[0] == "10"]
    if len(cases) != 44 or any(len(row) != 3 or row[0] != "10" for row in cases):
        raise RuntimeError("Expected the 44 reserved ten-loop origin cases")
    CONTROL.mkdir(parents=True, exist_ok=True)
    for name in ("audit", "log", "scratch"):
        (CONTROL / name).mkdir(exist_ok=True)
    lock = (CONTROL / "runner.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    identity = {
        "schema": "hexagon-ten-sequential-v16", "repo": str(REPO),
        "runner_sha256": sha(Path(__file__)), "measure_sha256": sha(MEASURE),
        "driver": str(DRIVER), "driver_sha256": sha(DRIVER),
        "converter": str(CONVERTER), "converter_sha256": sha(CONVERTER),
        "native": str(NATIVE), "shifted": str(SHIFTED),
        "fec": str(FEC), "lec": str(LEC), "data": str(DATA),
        "h": 17, "k": 3, "target": 9000, "prereduce": "right",
        "threads": 16, "physics_cache_mib": 4096, "ope_head_cache_mib": 8192,
        "empirical_assumptions": ASSUMPTIONS,
        "known_symbol_controls_sha256": sha(EXPERIMENT / "v13-controls/progress.json"),
        "v14_eight_comparison_sha256": sha(EXPERIMENT / "audit/v14-eight-comparison.json"),
        "v14_native_comparison_sha256": sha(EXPERIMENT / "audit/v14-native-comparison.json"),
        "v15_eight_comparison_sha256": sha(EXPERIMENT / "audit/v15-eight-comparison.json"),
        "v15_triangular_probe_sha256": sha(EXPERIMENT / "audit/v15-triangular-probe.log"),
        "ope_storage": "triangular_log_powers_a_plus_b_le_weight",
        "mrk_rhs_engine": "shared_one_loop_log_shuffle_v1",
        "mrk_rhs_prefetch": "4_for_single_sample;bounded_cache_128",
        "v16_integration_control_sha256": sha(EXPERIMENT / "v16-controls/progress.json"),
        "mrk_exact_control_sha256": sha(EXPERIMENT / "mrk-rhs-reuse/control/progress.json"),
        "certificate_binary_sha256": sha(CERTIFICATE),
        "origin_binary_sha256": sha(ORIGIN_CHECK),
        "origin_cases_sha256": sha(ORIGIN_CASES),
        "native_cache_release": "thread-local FLINT caches, joined pool; live QQ/Fp preserved",
        "conditional_tail_audit_sha256": sha(EXPERIMENT / "bonus-double-audit-repair/audit_complete.tsv"),
        "rss_limit_bytes": 50_000_000_000, "virtual_limit_gib": 256,
        "min_available_gib": 3, "mapped_capacity_gib": 192,
        "mapped_soft_rss_gib": 40, "timeout_seconds": 172800,
        "physical_strategy": PHYSICAL_STRATEGY,
        "physical_strategy_placement": "after cheap seam seed, between refinement stages and within long seam sweeps",
        "completion_scope": "conditional QQ candidate; exact original 16+4 seam, native-prime physical checks and 44 exact origin holdouts; unrestricted uniqueness not claimed",
    }
    if args.preflight_only:
        print(json.dumps({"preflight_passed": True, "configuration": identity}, indent=2))
        return 0
    identity_path = CONTROL / "manifest.json"
    if identity_path.exists():
        if json.loads(identity_path.read_text()) != identity:
            raise RuntimeError("Runner identity changed; use a reviewed separate experiment")
    else:
        if (NATIVE / "joint/cache").exists() or SHIFTED.exists():
            raise RuntimeError("Unmanaged existing result directories; refusing to adopt them")
        atomic_json(identity_path, identity)
    ledger_path = CONTROL / "ledger.json"
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else {"primes": {}, "attempts": []}
    progress_path = CONTROL / "progress.json"

    def progress(state, **fields):
        value = {"schema": "hexagon-ten-progress-v16", "state": state,
                 "runner_pid": os.getpid(), "updated_unix_seconds": time.time(),
                 "native_primes_committed": len(ledger["primes"]),
                 "max_primes": args.max_primes, "final_symbol_certified": False, **fields}
        atomic_json(progress_path, value)
        print(json.dumps(value), flush=True)

    def commit():
        atomic_json(ledger_path, ledger)

    def cleanup_scratch(attempt):
        path = Path(attempt["scratch"])
        if path.exists():
            if path.is_symlink() or path.parent != CONTROL / "scratch" or path.suffix != ".mmap":
                raise RuntimeError("Unsafe temporary scratch cleanup path")
            stat = path.stat()
            attempt["scratch_cleanup"] = {"logical_bytes": stat.st_size,
                                           "allocated_bytes": stat.st_blocks * 512,
                                           "removed_unix_seconds": time.time()}
            path.unlink()
            commit()

    def measured(phase, index, arguments):
        global active
        if cancel_signal:
            raise RuntimeError("Runner cancelled before starting next process")
        token = f"pindex{index:02d}-{phase}-{uuid.uuid4().hex[:12]}"
        scratch = CONTROL / "scratch" / (token + ".mmap")
        env_args = ["env", "THREADS=16", "PHYSICS_CACHE_MIB=4096", "OPE_HEAD_CACHE_MIB=8192",
                    "MIMALLOC_PURGE_DELAY=0", "HEXAGON_MHV_DOUBLE_BONUS=1",
                    f"HEXAGON_MAPPED_FILE={scratch}", "HEXAGON_MAPPED_GIB=192", "HEXAGON_MAPPED_SOFT_GIB=40"]
        report = CONTROL / "audit" / (token + ".json")
        log = CONTROL / "log" / (token + ".log")
        status = CONTROL / "audit" / (token + ".status.json")
        command = [sys.executable, str(MEASURE), "--report", str(report), "--log", str(log),
                   "--status", str(status), "--timeout", "172800", "--virtual-gib", "256",
                   "--rss-gib", str(RSS_GIB), "--min-available-gib", "3", "--", *env_args, *arguments]
        attempt = {"phase": phase, "prime_index": index, "command": command,
                   "scratch": str(scratch), "report": str(report), "log": str(log),
                   "status": str(status), "started_unix_seconds": time.time(), "state": "running"}
        ledger["attempts"].append(attempt)
        commit()
        with (CONTROL / "log" / (token + ".manager.log")).open("w") as console:
            active = subprocess.Popen(command, cwd=REPO, stdout=console, stderr=subprocess.STDOUT)
            if cancel_signal:
                active.send_signal(signal.SIGTERM)
            progress("running_" + phase, prime_index=index, manager_pid=active.pid,
                     command=command, resource_status=str(status), solver_log=str(log))
            code = active.wait()
            active = None
        attempt["finished_unix_seconds"] = time.time()
        attempt["manager_exit_code"] = code
        if report.exists():
            attempt["measurement"] = json.loads(report.read_text())
        measurement = attempt.get("measurement", {})
        if (code or measurement.get("exit_code") != 0 or any(measurement.get(key) for key in
                ("timeout", "memory_pressure", "rss_limit_exceeded", "cancel_signal")) or cancel_signal):
            attempt["state"] = "failed"
            commit()
            raise RuntimeError(f"{phase} stopped at prime index {index}; see {log} and {report}. No automatic retry")
        attempt["state"] = "process_succeeded_pending_output_validation"
        commit()
        return attempt

    def certify():
        if not check_qq():
            raise RuntimeError("Cannot certify an incomplete rational reconstruction")
        attempt = measured("exact-target-seam", 0,
                           [str(CERTIFICATE), "16", "4", str(SHIFTED), str(FEC), str(LEC), str(DATA)])
        records = re.findall(r"^SHIFTED_JOINT_CERTIFICATE (.+)$", Path(attempt["log"]).read_text(), re.M)
        if len(records) != 1:
            raise RuntimeError("Independent target-seam certificate receipt missing")
        certificate = json.loads(records[0])
        if (certificate.get("exact_QQ_original_target_seam") is not True or
                certificate.get("shift_source_hash_receipts_verified") is not True or
                certificate.get("directions") != 0 or certificate.get("columns") != 518048 or
                [certificate.get(k) for k in ("source_h", "source_k", "target_h", "target_k")] != [17, 3, 16, 4]):
            raise RuntimeError("Independent target-seam certificate scope mismatch")
        attempt["state"] = "verified_exact_original_target_seam"
        attempt["certificate"] = certificate
        commit()
        cleanup_scratch(attempt)
        origin = measured("exact-origin-holdouts", 0,
                          [str(ORIGIN_CHECK), "16", "4", str(SHIFTED), str(FEC), str(LEC), str(ORIGIN_CASES)])
        origin_text = Path(origin["log"]).read_text()
        if (re.findall(r"^ORIGIN_CHECK PASS words=(\d+)$", origin_text, re.M) != ["44"] or
                len(re.findall(r"^ORIGIN_WORD loop=10 ", origin_text, re.M)) != 44):
            raise RuntimeError("Exact 44-case origin holdout receipt missing")
        origin["state"] = "verified_exact_origin_holdouts"
        commit()
        cleanup_scratch(origin)
        final = {"schema": "hexagon-ten-conditional-certificate-v1", "qq_root": str(SHIFTED),
                 "qq_receipt_sha256": sha(SHIFTED / "qq_complete.tsv"),
                 "empirical_assumptions": ASSUMPTIONS, "unrestricted_uniqueness_claimed": False,
                 "exact_original_target_seam": certificate, "exact_origin_holdouts": 44,
                 "native_prime_physical_validation": "full COL, MRK and supplied OPE fit and fresh check rows on each committed prime",
                 "physical_scope": "finite supplied kinematic/OPE checks plus exact origin holdouts; no claim of all physical identities over QQ",
                 "certificate_log": attempt["log"], "origin_log": origin["log"],
                 "created_unix_seconds": time.time()}
        atomic_json(CONTROL / "final-verification.json", final)
        progress("conditional_candidate_certified_report_update_pending", qq_root=str(SHIFTED),
                 final_symbol_certified=True, certificate_scope=final)
        return 0

    try:
        if check_qq():
            return certify()
        for index in range(1, args.max_primes + 1):
            if cancel_signal:
                raise RuntimeError("Runner cancelled")
            entry = ledger["primes"].get(str(index))
            if entry:
                native = check_cache(NATIVE / "joint", entry["native"]["prime"], 17, 3, 275246)
                if native != entry["native"]:
                    raise RuntimeError("Committed native prime changed")
            else:
                # Recover a fully validated process that ended just before its
                # output was committed. Never rerun that completed native solve.
                attempt = None
                for prior in reversed(ledger["attempts"]):
                    if prior["phase"] != "native" or prior["prime_index"] != index:
                        continue
                    report = Path(prior["report"])
                    if not report.exists():
                        continue
                    result = json.loads(report.read_text())
                    if result.get("exit_code") == 0 and not any(result.get(key) for key in
                            ("timeout", "memory_pressure", "rss_limit_exceeded", "cancel_signal")):
                        attempt = prior
                        break
                if attempt is None:
                    attempt = measured("native", index, [str(DRIVER), "17", "3", str(FEC), str(LEC), str(DATA),
                        str(NATIVE), "9000", "cache", "right", "--prime-index", str(index)])
                text = Path(attempt["log"]).read_text()
                selected = re.findall(r"JOINT_SELECTED_PRIME index=(\d+) prime=(\d+)", text)
                if len(selected) != 1 or int(selected[0][0]) != index:
                    raise RuntimeError("Native prime-selection receipt missing or inconsistent")
                prime = int(selected[0][1])
                if not re.search(rf"JOINT_PRIME_COMPLETE prime={prime} directions=0 columns=275246\b", text):
                    raise RuntimeError("Native completion receipt with zero directions is absent")
                native = check_cache(NATIVE / "joint", prime, 17, 3, 275246)
                entry = {"native": native, "native_attempt": ledger["attempts"].index(attempt),
                         "native_identity_sha256": sha(NATIVE / "joint/input_identity.txt"),
                         "native_basis_identity_sha256": sha(NATIVE / "joint/basis_identity.tsv")}
                ledger["primes"][str(index)] = entry
                attempt["state"] = "verified_native_unique"
                commit()
                cleanup_scratch(attempt)
                progress("native_prime_committed", prime_index=index, prime=prime)
            prime = native["prime"]
            if entry.get("conversion"):
                converted = check_cache(SHIFTED, prime, 16, 4, 518048)
                if converted != entry["conversion"] or (SHIFTED / f"source_p{prime}.sha256").read_text().strip() != native["sha256"]:
                    raise RuntimeError("Committed conversion/provenance changed")
                continue
            attempt = measured("convert", index, [str(CONVERTER), "shift-cache", "MHV", "17", "3",
                str(NATIVE / "joint"), str(SHIFTED), "--head", "16", "--fec-root", str(FEC),
                "--lec-root", str(LEC), "--kind", "unique", "--reconstruct"])
            converted = check_cache(SHIFTED, prime, 16, 4, 518048)
            if (SHIFTED / f"source_p{prime}.sha256").read_text().strip() != native["sha256"]:
                raise RuntimeError("Shifted source-prime receipt mismatch")
            manifest = json.loads((SHIFTED / "manifest.json").read_text())
            qq = check_qq()
            if not qq and manifest.get("complete") is not False:
                raise RuntimeError("Conversion is neither valid cache-only nor verified QQ")
            text = Path(attempt["log"]).read_text()
            if not qq and "CACHE_ONLY completed_primes=" not in text:
                raise RuntimeError("Conversion cache-only completion receipt missing")
            entry["conversion"] = converted
            entry["conversion_attempt"] = len(ledger["attempts"]) - 1
            entry["shifted_identity_sha256"] = sha(SHIFTED / "input_identity.txt")
            attempt["state"] = "verified_conversion_qq" if qq else "verified_conversion_cache_only"
            commit()
            cleanup_scratch(attempt)
            if qq:
                progress("qq_reconstructed_awaiting_independent_certification", prime_index=index,
                         qq_root=str(SHIFTED), reconstruction_manifest=manifest)
                return certify()
            progress("converted_prime_committed_cache_only", prime_index=index, prime=prime)
        progress("prime_budget_exhausted_cache_only", shifted_root=str(SHIFTED))
        return 2
    except Exception as error:
        progress("interrupted" if cancel_signal else "failed", error=str(error), cancel_signal=cancel_signal)
        return 1


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    sys.exit(main())
