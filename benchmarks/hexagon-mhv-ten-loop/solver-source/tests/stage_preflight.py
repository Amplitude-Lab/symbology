#!/usr/bin/env python3
"""CLI rejection tests use copies; no historical RHS or output is modified."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess

parser = argparse.ArgumentParser()
for key in ("exe", "fec", "lec", "data", "output"):
    parser.add_argument("--" + key, required=True, type=Path)
a = parser.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
for case in ("missing", "grade", "role", "holdout", "duplicate", "receipt", "hash"):
    root = a.output / case
    data = root / "data" / "ope_mhv_small_s"
    shutil.copytree(str(a.data / "ope_L10_main" / "ope_mhv_small_s"), str(data))
    # Legacy <=9-loop views need no receipt. Rename our copied receipt so each
    # malformed semantic case reaches the table checker, not just a hash error.
    receipt = data / "coverage_identity.tsv"
    if case not in ("receipt", "hash"):
        receipt.rename(data / "coverage_identity.saved")
    table = data / "calE_L2.tsv"
    if case == "receipt":
        receipt.write_text(receipt.read_text().replace("complete", "partial"))
    elif case == "missing":
        table.rename(data / "calE_L2.saved")
    elif case == "grade":
        text = (data / "projection.tsv").read_text().replace("grade\t12", "grade\t4").replace('"grade"\t12', '"grade"\t4')
        (data / "projection.tsv").write_text(text)
    else:
        with table.open(newline="") as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        if case == "role":
            rows[1][-1] = "unknown"
        elif case == "hash":
            rows[1][6] = str(int(rows[1][6]) + 1)
        elif case == "duplicate":
            rows.append(rows[1])
        else:
            rows = [rows[0]] + [r for r in rows[1:] if r[-1] == "solve" and r[1] == "2"]
        with table.open("w", newline="") as stream:
            csv.writer(stream, delimiter="\t", lineterminator="\r\n").writerows(rows)
    (root / "condition.tsv").write_text("data\tdata\n")
    (root / "plan.tsv").write_text("ope ope_small_s condition.tsv\n")
    result = subprocess.run([str(a.exe.resolve()), "MHV", "3", "1", str(root / "output"), "--plan", str(root / "plan.tsv"), "--fec-root", str(a.fec.resolve()), "--lec-root", str(a.lec.resolve()), "--data", str(a.data.resolve()), "--cache-only", "--primes", "1"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=90)
    (root / "result.log").write_text(result.stdout)
    assert result.returncode != 0 and "BULK_MATRIX" not in result.stdout, result.stdout
    expected = {"missing": "Missing stage input", "grade": "exceeds kernel/log coverage", "role": "OPE role", "holdout": "nonempty solve and held-out", "duplicate": "Duplicate OPE coefficient", "receipt": "Incomplete OPE coverage", "hash": "OPE coverage file changed"}[case]
    assert expected in result.stdout, result.stdout
    print("PASS preflight", case)
