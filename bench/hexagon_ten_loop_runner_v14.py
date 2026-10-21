#!/usr/bin/env python3
"""Sequential V14 conditional, interleaved-physics 17+3 caches and 16+4 reconstruction.

Run this supervisor in an independently managed systemd user service. A rational
reconstruction is an intermediate result awaiting independent certification.
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
NATIVE = EXPERIMENT / "ten-v14-main"
SHIFTED = EXPERIMENT / "ten-v14-shifted-16p4"
CONTROL = EXPERIMENT / "ten-v14-runner"
DRIVER = EXPERIMENT / "joint-v14-cache-release-package/bin/joint_physics_v14"
CONVERTER = SHARED / "joint-v8-package/bin/irr_convert"
FEC, LEC, DATA = SHARED / "FEC", SHARED / "LEC/MHV", SHARED / "data"
MEASURE = REPO / "bench/measure_process.py"
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
                "joint-discovery-physical-affine-v14-bounded-interleaving-cache-release", PHYSICAL_STRATEGY, "9000"):
            raise RuntimeError("Native V14 physical-strategy identity mismatch")
    if root == NATIVE / "joint" and stage.get("empirical_assumptions") != ASSUMPTIONS:
        raise RuntimeError("Native conditional-assumption identity mismatch")
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
    CONTROL.mkdir(parents=True, exist_ok=True)
    for name in ("audit", "log", "scratch"):
        (CONTROL / name).mkdir(exist_ok=True)
    lock = (CONTROL / "runner.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    identity = {
        "schema": "hexagon-ten-sequential-v14", "repo": str(REPO),
        "runner_sha256": sha(Path(__file__)), "measure_sha256": sha(MEASURE),
        "driver": str(DRIVER), "driver_sha256": sha(DRIVER),
        "converter": str(CONVERTER), "converter_sha256": sha(CONVERTER),
        "native": str(NATIVE), "shifted": str(SHIFTED),
        "fec": str(FEC), "lec": str(LEC), "data": str(DATA),
        "h": 17, "k": 3, "target": 9000, "prereduce": "right",
        "threads": 16, "physics_cache_mib": 4096, "ope_head_cache_mib": 16384,
        "empirical_assumptions": ASSUMPTIONS,
        "known_symbol_controls_sha256": sha(EXPERIMENT / "v13-controls/progress.json"),
        "v14_eight_comparison_sha256": sha(EXPERIMENT / "audit/v14-eight-comparison.json"),
        "v14_native_comparison_sha256": sha(EXPERIMENT / "audit/v14-native-comparison.json"),
        "native_cache_release": "thread-local FLINT caches, joined pool; live QQ/Fp preserved",
        "conditional_tail_audit_sha256": sha(EXPERIMENT / "bonus-double-audit-repair/audit_complete.tsv"),
        "rss_limit_bytes": 50_000_000_000, "virtual_limit_gib": 256,
        "min_available_gib": 3, "mapped_capacity_gib": 192,
        "mapped_soft_rss_gib": 40, "timeout_seconds": 172800,
        "physical_strategy": PHYSICAL_STRATEGY,
        "physical_strategy_placement": "after cheap seam seed, between refinement stages and within long seam sweeps",
        "completion_scope": "reconstructed QQ awaiting independent target seam, origin and physical-scope verification",
    }
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
        value = {"schema": "hexagon-ten-progress-v14", "state": state,
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
        env_args = ["env", "THREADS=16", "PHYSICS_CACHE_MIB=4096", "OPE_HEAD_CACHE_MIB=16384",
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

    try:
        if check_qq():
            progress("qq_reconstructed_awaiting_independent_certification", qq_root=str(SHIFTED))
            return 0
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
                return 0
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
