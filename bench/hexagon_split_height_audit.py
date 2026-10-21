"""Measure canonical rational heights in fixed hexagon split conventions.

Only reads existing exact coefficient vectors; never fits or reconstructs one.
"""
import csv
import hashlib
import json
import math
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output_hexagon_split_heights_20261010"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_vector(path):
    result = []
    with path.open() as stream:
        for fields in csv.reader(stream, delimiter="\t"):
            if not fields or fields[0].startswith("#"):
                continue
            index, num, den = map(int, fields)
            assert index == len(result) and den > 0
            assert math.gcd(num, den) == 1
            result.append((num, den))
    assert result
    return result


def measure(job):
    path = Path(job["vector"])
    vector = read_vector(path)
    common = math.lcm(*(den for _, den in vector))
    max_num = max(abs(num) for num, _ in vector)
    max_den = max(den for _, den in vector)
    report = dict(job, coefficients=len(vector),
                  nonzero=sum(num != 0 for num, _ in vector),
                  max_numerator_bits=max_num.bit_length(),
                  max_denominator_bits=max_den.bit_length(),
                  max_denominator_decimal=str(max_den),
                  common_denominator_bits=common.bit_length(),
                  common_denominator_decimal=str(common),
                  maximum_individual_2ND_bits=max((2*abs(n)*d).bit_length() for n,d in vector),
                  global_box_2ND_bits=(2*max_num*max_den).bit_length(),
                  vector_sha256=sha(path))
    root = Path(job["root"])
    report["object_hashes"] = {
        str(p.relative_to(root)): sha(p) for p in sorted(root.iterdir())
        if p.is_file() and (p.name in ["columns.tsv", "object.tsv", "basis_identity.tsv",
                                     "inputs.sha256", "qq_complete.tsv"]
                            or p.suffix == ".wxf" or p.name.endswith(".wxf.shape"))
    }
    if job["split"][1] == 3 and job["loop"] in (8, 9):
        cache_root = ROOT / "output_hexagon_17p3_20261007"
        cache = cache_root / ("v17-schedule-controls/eight-t4/joint/cache/p1152921504606847009.bin"
                             if job["loop"] == 8 else
                             "v17-nine-control/nine/joint/cache/p1152921504606847009.bin")
        data = cache.read_bytes()
        magic, version, loop, h, k, rows, cols, prime, identity = struct.unpack_from("<8s8Q", data)
        assert magic == b"NMHVSEW\0" and version == 1
        assert (loop,h,k,rows,cols) == (job["loop"], *job["split"], 1, len(vector))
        assert len(data) == 72 + 8*len(vector)
        mismatches = sum(n % prime * pow(d, -1, prime) % prime != value[0]
                         for (n,d), value in zip(vector, struct.iter_unpack("<Q", data[72:])))
        assert mismatches == 0
        report["native_cache_check"] = dict(path=str(cache), sha256=sha(cache),
                                            prime=prime, mismatches=mismatches,
                                            coefficients_checked=len(vector))
    return report


jobs = json.loads((OUT / "decode-jobs.json").read_text())
for loop, name, root in [
    (8,"eight",ROOT/"output_hexagon_12p4_20261006/modular_direct/SEW_12p4/MHV/collinear_mrk_ope_small_s"),
    (9,"nine",ROOT/"output_hexagon_14p4_20261006/cached-nine/joint")]:
    jobs.append(dict(loop=loop, split=[2*loop-4,4], root=str(root),
                     vector=str(ROOT/f"output_hexagon_candidate_export_20261009/{name}-reference.canonical.tsv")))
jobs.extend([
    dict(loop=10, split=[17,3], root=str(ROOT/"output_hexagon_validation_20261009/native-qq-17p3"),
         vector=str(ROOT/"output_hexagon_validation_20261009/native-reconstruction-7p96.tsv")),
    dict(loop=10, split=[16,4], root=str(ROOT/"output_hexagon_17p3_20261007/ten-asymmetric-qq-16p4"),
         vector=str(ROOT/"output_hexagon_17p3_20261007/audit/crt-six-prime-speculative-candidate-20261008T1532Z.candidate.tsv"))])
records = [measure(job) for job in sorted(jobs, key=lambda x:(x["loop"],x["split"][1]))]
pairs = []
for loop in sorted(set(job["loop"] for job in jobs)):
    native, shifted = [row for row in records if row["loop"] == loop]
    d3, d4 = (int(row["common_denominator_decimal"]) for row in (native, shifted))
    pairs.append(dict(loop=loop, shifted_common_denominator_divides_native=(d3%d4==0),
                      native_common_denominator_divides_shifted=(d4%d3==0),
                      shifted_has_smaller_maximum_denominator=(int(shifted["max_denominator_decimal"])<int(native["max_denominator_decimal"])),
                      shifted_has_more_coordinates=shifted["coefficients"]>native["coefficients"]))
conversion_receipts = []
for name in ["L8-right-shift", "L9-native", "L4-shifted", "L6-shifted"]:
    resource = OUT/f"{name}.resources.json"
    log = OUT/f"{name}.log"
    measurement = json.loads(resource.read_text())
    assert measurement["exit_code"] == 0
    assert "SHIFT exact_factorization=PASS" in log.read_text()
    conversion_receipts.append(dict(path=str(resource),sha256=sha(resource),
                                    log=str(log),log_sha256=sha(log),
                                    measurement=measurement))
converter = ROOT/"output_hexagon_16p4_20261006/joint-v8-package/bin/irr_convert"
receipt=dict(schema="hexagon-exact-split-height-audit-v1", records=records, comparisons=pairs,
             exact_conversion_receipts=conversion_receipts,
             converter=dict(path=str(converter),sha256=sha(converter)),
             independent_decoder_log=dict(path=str(OUT/"decode.log"),sha256=sha(OUT/"decode.log")),
             scope="Fixed existing parity x D3, physical-copy-primitive-v1 conventions. L4/L6 exact left conversion; L8/L9 exact right conversion with all-coordinate native modular comparison; L10 historical certified QQ vectors. This does not certify the currently running fresh ten-loop calculation.",
             interpretation="Finite-loop empirical comparison, not an all-loop theorem or basis-independent rule.",
             scripts={str(p):sha(p) for p in [Path(__file__),ROOT/"bench/hexagon_split_height_decode.wls"]})
(OUT/"height-comparison.json").write_text(json.dumps(receipt,indent=2)+"\n")
for row in records:
    print(row["loop"],row["split"],"count",row["coefficients"],"nonzero",row["nonzero"],
          "Nbits",row["max_numerator_bits"],"Dbits",row["max_denominator_bits"],
          "LCMbits",row["common_denominator_bits"],"2ND",row["maximum_individual_2ND_bits"])
print(json.dumps(pairs,indent=2))
