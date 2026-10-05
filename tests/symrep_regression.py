#!/usr/bin/env python3
"""Public CLI fixtures: exact recursion, continuation, terminal actions and failures."""
import csv
import argparse
import pathlib
import shutil
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--solver", type=pathlib.Path, default=ROOT / "symrep")
parser.add_argument("--probe", type=pathlib.Path, default=ROOT / "tests/symrep_probe")
args = parser.parse_args()


def run(*args, fail=None):
    result = subprocess.run([str(SOLVER), *map(str, args)], text=True, capture_output=True)
    if fail is None:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0 and fail in result.stderr, result.stdout + result.stderr
    return result.stdout


def dimension(path):
    with path.open() as stream:
        return sum(int(r["copies"]) * int(r["dimension"]) for r in csv.DictReader(stream, delimiter="\t"))


SOLVER = args.solver.resolve()
with tempfile.TemporaryDirectory(prefix="symrep-regression-") as name:
    work = pathlib.Path(name)
    fixture = work / "fixtures"
    subprocess.run([str(args.probe.resolve()), "fixtures", str(fixture)], check=True)
    generators = ["--generator", f"cyclic={fixture / 'cycle.wxf'}", "--generator", f"flip={fixture / 'flip.wxf'}"]
    for matrix, expected in [("sum", 2), ("identity", 0), ("empty", 3)]:
        output = work / matrix
        run("kernel", *generators, "--matrix", fixture / f"{matrix}.wxf", "--output", output)
        assert dimension(output / "kernel_copies.tsv") == expected
        run("kernel", *generators, "--matrix", fixture / f"{matrix}.wxf", "--output", output, fail="already exists")
    run("kernel", *generators, "--matrix", fixture / "bad.wxf", "--output", work / "bad", fail="not closed")
    run("kernel", "--input", work / "sum", "--matrix", fixture / "empty2.wxf", "--output", work / "recursive-free")
    assert dimension(work / "recursive-free/kernel_copies.tsv") == 2
    run("kernel", "--input", work / "recursive-free", "--matrix", fixture / "identity2.wxf", "--output", work / "recursive-zero")
    assert dimension(work / "recursive-zero/kernel_copies.tsv") == 0
    run("kernel", "--input", work / "recursive-zero", "--matrix", fixture / "empty0.wxf", "--output", work / "recursive-empty")
    assert dimension(work / "recursive-empty/kernel_copies.tsv") == 0
    run("adapt", *generators, "--threads", "0", "--output", work / "invalid", fail="positive")
    run("adapt", *generators, "--max-order", "2", "--output", work / "limit", fail="exceeds")

    for direction, seed in [("forward", "fec"), ("backward", "lec")]:
        prepared, partial, continued, ordinary = [work / f"{direction}-{tag}" for tag in ("prepared", "partial", "continued", "ordinary")]
        run("prepare", *generators, "--condition", fixture / "condition.wxf", "--seed", fixture / f"{seed}.wxf", "--direction", direction, "--output", prepared)
        assert (prepared / "dlogmat_multiplicity.wxf").exists()
        assert (prepared / "dlogmat_pair_basis.wxf").exists()
        assert (prepared / "dlogmat_reduced.tsv").read_text() == "symbology-equivariant-condition-v1\n"
        run("extend", "--input", prepared, "--max-weight", 2, "--output", partial, "--compare")
        run("extend", "--input", prepared, "--chain", partial, "--max-weight", 3, "--output", continued, "--compare")
        assert [dimension(continued / f"w{w}_copies.tsv") for w in range(1, 4)] == [3, 6, 10]
        assert "multiplicity" in (continued / "chain.tsv").read_text()
        assert not (continued / "w3.wxf").exists()
        assert (partial / "w2_multiplicity.wxf").read_bytes() == (continued / "w2_multiplicity.wxf").read_bytes()
        run("ordinary", "--condition", fixture / "condition.wxf", "--seed", fixture / f"{seed}.wxf", "--direction", direction, "--max-weight", 3, "--output", ordinary)
        streamed = work / f"{direction}-streamed"
        run("extend", "--input", prepared, "--backend", "factorized", "--kernel-strategy", "streamed", "--max-weight", 3, "--compare", "--output", streamed)
        run("verify", "--input", prepared, "--chain", streamed, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{direction}-streamed-verified")
        staged = work / f"{direction}-staged"
        run("extend", "--input", prepared, "--backend", "factorized", "--kernel-strategy", "staged", "--max-weight", 3, "--compare", "--output", staged)
        run("verify", "--input", prepared, "--chain", staged, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{direction}-staged-verified")
        run("verify", "--input", prepared, "--chain", continued, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{direction}-verified")
        expanded = work / f"{direction}-multiplicity-expanded"
        run("expand", "--input", prepared, "--chain", continued, "--output", expanded)
        run("verify", "--input", prepared, "--chain", expanded, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{direction}-multiplicity-expanded-verified")
        run("extend", "--input", prepared, "--chain", expanded, "--max-weight", 4, "--compare", "--output", work / f"{direction}-multiplicity-expanded-continued")
        # Compare and normal recursion must save the identical complete recipe.
        for workers in (1, 3):
            direct = work / f"{direction}-multiplicity-direct-{workers}"
            run("extend", "--input", prepared, "--max-weight", 3, "--threads", workers, "--output", direct)
            for file in direct.glob("w*"):
                assert file.read_bytes() == (continued / file.name).read_bytes(), file.name
        resumed = work / f"{direction}-multiplicity-direct-continued"
        run("extend", "--input", prepared, "--chain", direct, "--max-weight", 4, "--output", resumed)
        assert (resumed / "w4_multiplicity.wxf").read_bytes() == (work / f"{direction}-multiplicity-expanded-continued/w4_multiplicity.wxf").read_bytes()
        with (resumed / "timings.tsv").open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                assert None not in row and float(row["pipeline_s"]) >= float(row["planning_s"])
        run("extend", "--input", prepared, "--max-weight", 1, "--output", work / f"{direction}-multiplicity-seed")
        run("extend", "--input", prepared, "--chain", work / f"{direction}-multiplicity-seed", "--max-weight", 3, "--compare", "--output", work / f"{direction}-multiplicity-seed-continued")
    run("extend", "--input", work / "backward-prepared", "--chain", work / "forward-partial", "--max-weight", 3, "--output", work / "mismatch", fail="different prepared bundle")

    prepared = work / "zero-prepared"
    run("prepare", *generators, "--condition", fixture / "full_condition.wxf", "--seed", fixture / "fec.wxf", "--direction", "forward", "--output", prepared)
    run("extend", "--input", prepared, "--max-weight", 3, "--compare", "--output", work / "zero-chain")
    assert dimension(work / "zero-chain/w3_copies.tsv") == 0
    unconstrained = work / "free-prepared"
    run("prepare", *generators, "--condition", fixture / "no_condition.wxf", "--seed", fixture / "fec.wxf", "--direction", "forward", "--output", unconstrained)
    run("extend", "--input", unconstrained, "--max-weight", 3, "--compare", "--output", work / "free-chain")
    assert dimension(work / "free-chain/w3_copies.tsv") == 27

    terminals = ["--terminal-generator", f"cyclic={fixture / 'cycle.wxf'}", "--terminal-generator", f"flip={fixture / 'flip.wxf'}"]
    prepared = work / "vector-prepared"
    run("prepare", *generators, *terminals, "--condition", fixture / "condition.wxf", "--seed", fixture / "vector_seed.wxf", "--direction", "forward", "--output", prepared)
    run("extend", "--input", prepared, "--max-weight", 3, "--compare", "--output", work / "vector-chain")
    assert dimension(work / "vector-chain/w3_copies.tsv") == 30
    backward_vector = work / "backward-vector-prepared"
    run("prepare", *generators, *terminals, "--condition", fixture / "condition.wxf", "--seed", fixture / "vector_seed.wxf", "--direction", "backward", "--output", backward_vector)
    run("extend", "--input", backward_vector, "--max-weight", 3, "--compare", "--output", work / "backward-vector-chain")
    assert dimension(work / "backward-vector-chain/w3_copies.tsv") == 30

    # The factorized logical irrep basis must survive both recursive directions,
    # saved continuation, non-scalar terminal actions and zero/full kernels.
    for case, direction, seed, condition in [
        ("forward", "forward", "fec", "condition"),
        ("backward", "backward", "lec", "condition"),
        ("vector", "forward", "vector_seed", "condition"),
        ("backward-vector", "backward", "vector_seed", "condition"),
        ("zero", "forward", "fec", "full_condition"),
        ("free", "forward", "fec", "no_condition"),
    ]:
        prepared = work / f"{case}-prepared"
        partial, continued, ordinary, expanded = [work / f"{case}-factorized-{tag}" for tag in ("partial", "continued", "ordinary", "expanded")]
        run("extend", "--input", prepared, "--backend", "factorized", "--max-weight", 2, "--compare", "--output", partial)
        run("extend", "--input", prepared, "--chain", partial, "--max-weight", 3, "--compare", "--output", continued)
        assert "factorized" in (continued / "chain.tsv").read_text()
        # Normal solving releases the final QQ carrier before orbit selection;
        # its saved result must still match the full check and support resume.
        direct_partial, direct_continued = [work / f"{case}-factorized-direct-{tag}" for tag in ("partial", "continued")]
        run("extend", "--input", prepared, "--backend", "factorized", "--max-weight", 2, "--output", direct_partial)
        run("extend", "--input", prepared, "--chain", direct_partial, "--max-weight", 3, "--output", direct_continued)
        for file in direct_continued.glob("w*"):
            assert file.read_bytes() == (continued / file.name).read_bytes(), file.name
        for f in partial.glob("w*"):
            assert f.read_bytes() == (continued / f.name).read_bytes()
        run("ordinary", "--condition", fixture / f"{condition}.wxf", "--seed", fixture / f"{seed}.wxf", "--direction", direction, "--max-weight", 3, "--output", ordinary)
        run("verify", "--input", prepared, "--chain", continued, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{case}-factorized-verified")
        run("expand", "--input", prepared, "--chain", continued, "--output", expanded)
        run("verify", "--input", prepared, "--chain", expanded, "--reference", ordinary, "--max-weight", 3, "--output", work / f"{case}-expanded-verified")
        assert (expanded / "w1.wxf").read_bytes() == (prepared / "seed.wxf").read_bytes()
        run("extend", "--input", prepared, "--chain", expanded, "--max-weight", 4, "--compare", "--output", work / f"{case}-expanded-continued")
    run("extend", "--input", work / "forward-prepared", "--chain", work / "forward-factorized-partial", "--backend", "multiplicity", "--max-weight", 3, "--output", work / "wrong-backend", fail="expand the factorized chain")
    run("extend", "--input", work / "forward-prepared", "--backend", "invalid", "--output", work / "bad-backend", fail="backend must be")
    run("expand", "--input", work / "forward-prepared", "--chain", work / "forward-multiplicity-expanded", "--output", work / "wrong-expansion", fail="expand requires")

    # Public nonsplit example: C3's two-dimensional QQ irrep has End dimension
    # two, so auto must use the factorized format even without private heptagon data.
    cyclic = work / "cyclic-prepared"
    run("prepare", *generators[:2], "--condition", fixture / "condition.wxf", "--seed", fixture / "fec.wxf", "--direction", "forward", "--output", cyclic)
    run("extend", "--input", cyclic, "--max-weight", 1, "--output", work / "cyclic-seed")
    assert "factorized" in (work / "cyclic-seed/chain.tsv").read_text()
    run("extend", "--input", cyclic, "--chain", work / "cyclic-seed", "--max-weight", 3, "--compare", "--output", work / "cyclic-chain")
    run("extend", "--input", cyclic, "--max-weight", 3, "--output", work / "cyclic-direct")
    for file in (work / "cyclic-direct").glob("w*"):
        assert file.read_bytes() == (work / "cyclic-chain" / file.name).read_bytes(), file.name
    run("verify", "--input", cyclic, "--chain", work / "cyclic-chain", "--reference", work / "forward-ordinary", "--max-weight", 3, "--output", work / "cyclic-verified")

    # Check action export in both coordinate systems against the actual recursive
    # tensors/frames, including non-split irreps, vector terminals and zero kernels.
    for case, chain_name in [("forward", "forward-continued"), ("backward-vector", "backward-vector-chain"),
                             ("zero", "zero-chain"), ("forward", "forward-factorized-direct-continued"),
                             ("backward", "backward-factorized-direct-continued"), ("cyclic", "cyclic-direct")]:
        source = work / chain_name
        for coordinates in (["adapted", "carrier"] if "factorized" in (source / "chain.tsv").read_text() else ["adapted"]):
            target = work / f"{chain_name}-actions-{coordinates}"
            run("actions", "--input", work / f"{case}-prepared", "--chain", source, "--coordinates", coordinates, "--output", target)
            subprocess.run([str(args.probe.resolve()), "verify-actions", str(work / f"{case}-prepared"), str(source), str(target), coordinates], check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    run("actions", "--input", cyclic, "--chain", work / "cyclic-direct", "--max-weight", 4, "--output", work / "actions-too-far", fail="exceed")
    run("actions", "--input", cyclic, "--chain", work / "cyclic-direct", "--coordinates", "unknown", "--output", work / "actions-bad-coordinates", fail="coordinates must be")
    run("actions", "--input", work / "forward-prepared", "--chain", work / "forward-continued", "--coordinates", "carrier", "--output", work / "actions-no-carrier", fail="factorized chain")

    tampered = work / "tampered"
    shutil.copytree(prepared, tampered)
    with (tampered / "seed.wxf").open("ab") as stream:
        stream.write(b"changed")
    run("extend", "--input", tampered, "--output", work / "rejected", fail="changed bundle file")
    truncated = work / "truncated-seal"
    shutil.copytree(prepared, truncated)
    seal = truncated / "checksums.tsv"
    seal.write_text(seal.read_text().splitlines()[0] + "\n")
    run("extend", "--input", truncated, "--output", work / "rejected-seal", fail="inventory is incomplete")

print("PASS symrep CLI: kernels, FEC/LEC, continuation, vector terminal, empty spaces, factorized frames/export, recursive action exports, bundle integrity")
