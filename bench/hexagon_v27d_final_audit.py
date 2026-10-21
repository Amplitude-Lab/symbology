"""Read-only review of the completed frozen V27D ten-loop experiment.

Rechecks receipts, every cached coefficient, rational congruences and exact
reference tensor equality. Does not run or change a numerical producer.
"""
from array import array
from collections import defaultdict
from pathlib import Path
import json
import math
import re
import struct
import time

from hexagon_boot_recovery import atomic_json, sha

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "output_hexagon_ten_v27d_20261009"
CONTROL = RUN / "runner"
OUTPUT = CONTROL / "audit/final-independent-review.json"
assert not OUTPUT.exists(), "Preserve the existing review"
bindings = {}


def bind(path):
    path = Path(path)
    assert path.is_file() and not path.is_symlink(), path
    bindings[str(path)] = sha(path)
    return path


def read_json(path):
    return json.loads(bind(path).read_text())


def fields(path):
    return dict(line.split("\t", 1) for line in bind(path).read_text().splitlines())


manifest = read_json(CONTROL / "manifest.json")
ledger = read_json(CONTROL / "ledger.json")
progress = read_json(CONTROL / "progress.json")
final = read_json(CONTROL / "final-verification.json")
frozen = read_json(CONTROL / "asymmetric-candidate.json")
exit_receipt = read_json(CONTROL / "last-service-exit.json")
assert progress["state"] == "complete_conditional_candidate_certified"
assert progress["final_symbol_certified"] and progress["certificate"] == final
assert exit_receipt["SERVICE_RESULT"] == "success" and exit_receipt["EXIT_STATUS"] == "0"
assert not (CONTROL / "attention-required.json").exists()
assert list(ledger["primes"]) == list(map(str, range(1, 9)))
for path, digest in manifest["hashes"].items():
    assert sha(bind(path)) == digest, path
assert sha(bind(manifest["driver"])) == manifest["driver_sha256"] == (
    "493bbc1a5712a27060e3f3641b75d1820ff031fb4658d022e7c4df733f2d1091")

qq = Path(final["qq_root"])
export = read_json(qq / "manifest.json")
assert export["all_native_producers"] == "V27D"
for path, digest in export["input_hashes"].items():
    assert sha(bind(path)) == digest, path
assert sha(bind(qq / "qq_complete.tsv")) == final["qq_receipt_sha256"]
for name, digest in fields(qq / "qq_complete.tsv").items():
    if name == "schema":
        assert digest == "irr-physics-qq-v1"
    elif name == "identity":
        assert sha(bind(qq / "input_identity.txt")) == digest
    else:
        assert sha(bind(qq / name)) == digest, name
for name, digest in fields(qq / "cache_files.tsv").items():
    assert sha(bind(qq / "cache" / name)) == digest, name
roles = fields(qq / "reconstruction_primes.tsv")
assert list(roles.values()).count("CRT") == 6
assert list(roles.values()).count("validation") == 2

H = struct.Struct("<8s8Q")


def cache(path, h, k, columns, prime, expected=None):
    path = bind(path)
    data = path.read_bytes()
    header = H.unpack_from(data)
    assert header[:8] == (b"NMHVSEW\0", 1, 10, h, k, 1, columns, prime)
    assert header[8] == int(sha(bind(path.parent.parent / "input_identity.txt"))[:16], 16)
    assert len(data) == 72 + 8*columns
    assert sha(path) == bind(path.with_name(path.name + ".sha256")).read_text().strip()
    if expected:
        assert sha(path) == expected
    return data


prime_reviews = []
shifted_values = []
for index, entry in ledger["primes"].items():
    assert entry["producer_algorithm"] == "V27D"
    prime = entry["native"]["prime"]
    assert prime == entry["conversion"]["prime"]
    matches = {}
    for kind, h, k, columns, reference in [
        ("native", 17, 3, 275246, Path(entry["reference_comparison"]["reference"])),
        ("conversion", 16, 4, 518048,
         ROOT / f"output_hexagon_17p3_20261007/ten-v17-shifted-16p4/cache/p{prime}.bin")]:
        item = entry[kind]
        actual = cache(item["path"], h, k, columns, prime, item["sha256"])
        expected = cache(reference, h, k, columns, prime)
        assert actual[72:] == expected[72:]
        matches[kind] = {"columns": columns, "all_coefficients_equal": True,
                         "path": item["path"], "reference": str(reference)}
        if kind == "conversion":
            assert bind(qq / f"source_p{prime}.sha256").read_text().strip() == entry["native"]["sha256"]
            assert (qq / "cache" / f"p{prime}.bin").read_bytes() == actual
            values = array("Q")
            values.frombytes(actual[72:])
            assert all(value < prime for value in values)
            shifted_values.append((prime, values))
    prime_reviews.append({"prime_index": int(index), "prime": prime, "comparisons": matches})

vector = bind(frozen["candidate_path"])
assert sha(vector) == frozen["candidate_sha256"] == export["candidate_sha256"]
assert frozen["candidate_complete"] and not frozen["failed_columns"]
rows = nonzero = numerator_bits = denominator_bits = 0
with vector.open() as stream:
    assert next(stream).startswith("# speculative")
    for line in stream:
        index, n, d = map(int, line.split())
        assert index == rows and d > 0 and math.gcd(n, d) == 1
        assert all(d % p and (n - d*values[index]) % p == 0 for p, values in shifted_values)
        rows += 1
        nonzero += n != 0
        numerator_bits = max(numerator_bits, abs(n).bit_length())
        denominator_bits = max(denominator_bits, d.bit_length())
assert (rows, nonzero, numerator_bits, denominator_bits) == (518048, 491197, 44, 231)
checks = read_json(export["fresh_check_receipt"])
assert checks["both_fresh_checks_pass"] and checks["candidate_sha256"] == sha(vector)
assert checks["candidate_receipt_sha256"] == sha(CONTROL / "asymmetric-candidate.json")
assert [row["prime_index"] for row in checks["checks"]] == [7, 8]
for row in checks["checks"]:
    assert row["mismatches"] == row["noninvertible"] == 0
    assert row["prime"] not in frozen["crt_primes"]
    assert sha(bind(row["path"])) == row["sha256"]

phase_totals = defaultdict(lambda: {"cpu_hours": 0, "wall_seconds": 0,
                                   "peak_sampled_rss_bytes": 0})
native_profiles = {}
kernel_peak = 0
attempts = ledger["attempts"]
for attempt in attempts:
    measurement = read_json(attempt["report"])
    assert measurement == attempt["measurement"] and measurement["exit_code"] == 0
    assert not any(measurement.get(k) for k in ("timeout", "memory_pressure", "rss_limit_exceeded", "cancel_signal"))
    text = bind(attempt["log"]).read_text()
    total = phase_totals[attempt["phase"]]
    total["cpu_hours"] += (measurement["cpu_user_seconds"] + measurement["cpu_system_seconds"])/3600
    total["wall_seconds"] += measurement["wall_seconds"]
    total["peak_sampled_rss_bytes"] = max(total["peak_sampled_rss_bytes"], measurement["peak_summed_rss_bytes"])
    if attempt["phase"] == "native":
        assert "V18_SEAM_VERIFIED full_operator=1 rows=989800 directions=0" in text
        for family in ("collinear", "mrk", "ope_small_s"):
            assert f"JOINT_PHYSICAL_CHECK operator={family} passed=1 remaining=0" in text
        native_profiles[str(attempt["prime_index"])] = {
            match[0]: float(match[1]) for match in re.findall(
                r"^PROFILE phase=(compute_\S+) seconds=\S+ cpu_seconds=(\S+)", text, re.M)}
        kernel_peak = max(kernel_peak, *(int(v) for v in re.findall(r"process_peak_bytes=(\d+)", text)))
    elif attempt["phase"] == "asymmetric-export":
        assert "INDEPENDENT_WOLFRAM_VERIFY_PASS columns=518048 nnz=491197" in text
        assert "V27D_QQ_EXPORT_COMMITTED columns=518048 crt_primes=6 validation_primes=2" in text
    elif attempt["phase"] == "exact-target-seam":
        certificate = json.loads(re.findall(r"^SHIFTED_JOINT_CERTIFICATE (.+)$", text, re.M)[0])
        assert certificate == final["exact_target_seam"]
        assert certificate["exact_QQ_original_target_seam"] and certificate["shift_source_hash_receipts_verified"]
        assert [certificate[k] for k in ("source_h", "source_k", "target_h", "target_k", "columns", "directions", "original_rows")] == [17,3,16,4,518048,0,4196956]
        assert certificate["modulus_bits"] > certificate["bound_bits"]
        prime_lines = re.findall(r"^ORIGINAL_TARGET_SEAM_PRIME prime=(\d+) rows=4196956 affine_rows=1 derived_cache_congruence=([01]) native_target_solve=0 passed=1$", text, re.M)
        assert len(prime_lines) == len(set(p for p, _ in prime_lines)) == 9
        assert sum(int(derived) for _, derived in prime_lines) == 8
    elif attempt["phase"] == "exact-origin-holdouts":
        words = re.findall(r"^ORIGIN_WORD loop=10 word=(\d+) expected=(-?\d+) actual=(-?\d+)$", text, re.M)
        assert len(words) == len(set(w for w, _, _ in words)) == 44
        assert all(expected == actual for _, expected, actual in words)
        assert "ORIGIN_CHECK PASS words=44" in text

freeze = next(a for a in attempts if a["phase"] == "asymmetric-freeze")
for index in (7,8):
    native = next(a for a in attempts if a["phase"] == "native" and a["prime_index"] == index)
    assert frozen["recorded_unix_seconds"] <= freeze["finished"] < native["started"]
assert max(a["finished"] for a in attempts if a["phase"] in ("native", "convert")) < next(a["started"] for a in attempts if a["phase"] == "asymmetric-fresh-checks")

old_root = ROOT / "output_hexagon_17p3_20261007/ten-asymmetric-qq-16p4"
exact_matches = []
assert (qq / "columns.tsv").read_bytes() == (old_root / "columns.tsv").read_bytes()
assert (qq / "basis_identity.tsv").read_bytes() == (old_root / "basis_identity.tsv").read_bytes()
for path in sorted(qq.glob("*.wxf")):
    old = bind(old_root / path.name)
    assert bind(path).read_bytes() == old.read_bytes()
    exact_matches.append({"filename": path.name, "sha256": sha(path), "all_tensor_bytes_equal": True})
assert len(exact_matches) == 6

cpu = sum(total["cpu_hours"] for total in phase_totals.values())
assert abs(cpu - progress["accounted_cpu_hours"]) < 1e-9
events = dict((k,int(v)) for k,v in (line.split() for line in exit_receipt["cgroup_counters"]["memory.events"].splitlines()))
assert events["max"] == events["oom"] == events["oom_kill"] == 0
service_cpu = dict((k,int(v)) for k,v in (line.split() for line in exit_receipt["cgroup_counters"]["cpu.stat"].splitlines()))
service_cpu_hours = (service_cpu["user_usec"] + service_cpu["system_usec"])/3.6e9
report = {
    "schema": "hexagon-v27d-final-independent-review-v1", "created": time.time(), "passed": True,
    "qq_root": str(qq), "source_split": [17,3], "reconstructed_split": [16,4],
    "fresh_native_primes": 8, "reconstruction_primes": 6, "reserved_check_primes": 2,
    "native_coefficients_per_prime": 275246, "exact_rational_coordinates": rows,
    "nonzero_exact_rational_coordinates": nonzero, "max_numerator_bits": numerator_bits,
    "max_denominator_bits": denominator_bits,
    "all_native_and_shifted_reference_comparisons": prime_reviews,
    "every_rational_congruence_rechecked": 8*rows,
    "reserved_primes_frozen_before_computation": True,
    "independent_Wolfram_decode_verified": True,
    "six_exact_tensor_files_identical_to_previous_result": exact_matches,
    "exact_seam_certificate": final["exact_target_seam"], "exact_origin_holdouts": 44,
    "each_native_prime_full_seam_rows": 989800, "supplied_COL_MRK_OPE_passed_every_prime": True,
    "phase_totals": dict(phase_totals), "process_user_plus_system_cpu_hours": cpu,
    "whole_service_user_plus_system_cpu_hours": service_cpu_hours,
    "CPU_accounting_scope": "Measured pipeline commands versus inclusive service-cgroup total. The small difference includes controller/measurement overhead and any descendants not reaped by command wrappers; do not assign it to a specific numerical phase.",
    "wall_hours": (attempts[-1]["finished"]-attempts[0]["started"])/3600,
    "peak_sampled_process_RSS_bytes": max(t["peak_sampled_rss_bytes"] for t in phase_totals.values()),
    "native_kernel_peak_RSS_bytes": kernel_peak,
    "service_cgroup_peak_bytes": int(exit_receipt["cgroup_counters"]["memory.peak"]),
    "memory_events": events, "configured_memory_max_bytes": 50000000000, "swap_bytes": 0,
    "completed_native_phase_cpu_seconds": native_profiles,
    "validation_scope": "Exact original16+4 sewing residuals (9-prime integer bound certificate), 44 exact reserved origin words, two reserved reconstruction primes and supplied per-prime physical checks. No unrestricted uniqueness or all-orders proof of bonus relations.",
    "assumptions_recorded_by_frozen_solver": final["empirical_assumptions"],
    "qualification": "V16 remains fixed baseline (>15% CPU saving, <=5% RSS increase). The historical ten-loop first-prime comparison saves37.3%CPU but exceeds RSS tolerance; direct nine-loop V16 qualification remains unverified. Correct completion does not establish performance optimality.",
    "optimization_priorities": ["dominant discovery coefficient-update arithmetic and locality", "exact-capacity pivot-only transpose/backsolve storage", "bounded allocator lifetime and reclaim policy", "adaptive reconstruction representation selected by measured heights and conversion cost"],
    "all_inputs_checkpoints_and_logs_preserved": True, "numerical_producer_modified": False,
    "audit_method": "Recomputed CPU totals, modular congruences and every saved native/shifted reference payload; independently rebound certificates to inputs, output seals and logs; exact old/new WXF tensor equality. Did not rerun the completed bootstrap or physical certificate programs.",
    "input_hashes": bindings, "audit_source_sha256": sha(Path(__file__))}
atomic_json(OUTPUT, report)
print(json.dumps({k:report[k] for k in ["passed", "process_user_plus_system_cpu_hours", "whole_service_user_plus_system_cpu_hours", "wall_hours", "peak_sampled_process_RSS_bytes", "native_kernel_peak_RSS_bytes", "service_cgroup_peak_bytes", "every_rational_congruence_rechecked"]}, indent=2))
print("FINAL_REVIEW", OUTPUT)
