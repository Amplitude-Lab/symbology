# Reproduce the four-loop bootstrap of the four-point form factor

This benchmark rebuilds the weight-eight (four-loop) four-point MHV form-factor
symbol in the finite normalization. It sews the FEC prefix recursion against
the two-letter suffix space on the retained 88-letter alphabet, imposes
integrability, first-entry, cyclic/flip/Galois symmetry and the proposed
adjacency pair spaces, and fixes the remaining three purely logarithmic
ambiguities by exact double-collinear normal regularity and
descendant-alphabet support. The triple-collinear image is then checked on the
complete boundary. Both finite four-loop descendants are **reserved holdouts**:
they are compared only after the candidate is fixed, against the published
four-loop three-point form-factor and six-point amplitude remainders.

The accepted candidate
[reference/F4_candidate_prefix6_suffix2.wxf](reference/F4_candidate_prefix6_suffix2.wxf)
is committed: **22,457,737 bytes**, shape `1 × 35598 × 175`, 988,559 nonzero
coefficients, SHA256
`5fd5b1d60c703bd303804307d467aec7b713213e278c923707db6580e0238f38`.
It is byte-identical to the DC affine-solve origin and to the TC copy
(`TC_copy_byte_identical` in
[reference/four_loop_construction_certificate.json](reference/four_loop_construction_certificate.json)).
The decoding bases FEC_1..FEC_5, LEC_2 and `suffix_TF175_physical.wxf` are
committed in `reference/`; the 503 MB FEC_6 basis is regenerated, not
uploaded. No Wolfram installation is needed for `verify` or `replay-final`.

**Scope.** Unique tensor within the explicit proposed-adjacency and continued
twelve-dimensional literature-remainder (TR12) hypotheses, at symbol level in
the finite normalization. External Steinmann equations were not imposed and
finite four-loop targets did not enter the solving. This does not claim
uniqueness in the larger extended-Steinmann space, any all-loop statement, or
equality to a separately known full four-loop parent symbol. A published
four-loop MHV form-factor symbol was not available for comparison; the two
published finite descendants pass exactly.

## Package layout

| Directory | Contents |
| --- | --- |
| `data/` | Small tracked inputs: interior constraint tensor, adjacency pair spaces, seed FEC_1/LEC_1, symmetry matrices, alphabet88 |
| `reference/` | Accepted candidate, decoding bases, sewing outputs (affine origin + TR92 directions), all consolidated certificates, recorded stage receipts |
| `boundary-preparation/DC/`, `boundary-preparation/TC/` | Boundary inputs of the final solves; the two weight-6 collinear projection maps are committed xz-compressed (142.7 MB → 46.1 MB and 458.5 MB → 81.0 MB) |
| `validation/` | DC/TC affine-solve receipts and block chains, holdout receipts/guards/logs and published targets, ambiguity certificates, replayed Lucas prime certificates, H8 boundary factor exports |
| `solver-source/` | Helper `.cpp`/`.py`/`.wls` sources, drivers, audits and H8 factor build. The accepted solver is frozen with all its headers in `ff_four_loop/sew_affine_TR_build_v1/`, exactly as recovered from the certificate-pinned research-tree source; fresh builds use this copy. The later `build_affine_tr_sew.py` guarded successor and `build_masked_sew.py` experiment are retained only as provenance |

## Verify the packaged result (seconds, pure Python)

```bash
python3 bench/four_loop_formfactor.py verify
```

Checks every file against [manifest.json](manifest.json), replays the
evidence-chain audit of the packaged receipts (sewing certificate dimension
and residue-bit balance, affine-solve block chains, holdout receipts, the
complete `CERT_PRIME` logs), re-proves the two accepted primes with exact
Lucas certificates, and verifies the guard hashes of the holdout inputs
(through decompression for the two xz maps).

## Replay the final accepted stages (~2 min measured)

Rebuilds the DC and TC exact affine boundary solves and both reserved holdout
comparisons from sewing outputs and committed boundary preparation. With no
`--sewing-dir`, it uses the committed reference outputs:

```bash
python3 bench/four_loop_formfactor.py replay-final --output output/formfactor-4l-replay
```

To consume a fresh `run-core` result instead, append
`--sewing-dir output/formfactor-4l-core`.

The regenerated DC origin must be **byte-identical** to the committed
candidate, the TC solve (which consumes the regenerated DC solution) must
reproduce it again, and both holdout comparisons must return
`exact_boundary_equality` with the packaged component counts and certificate
moduli. Solver binaries are resolved from `--solver-bin DIR`, or from the
local research tree after checking them against the receipt-pinned SHA256
values, or compiled fresh from `solver-source/` (needs the toolchain below).
Every stage runs under a sampled memory cap (default 12 GB) and writes
receipts to the output directory.

## Accepted replay: the only normal core path

Run the following from the repository root with fresh output directories:

```bash
python3 bench/four_loop_formfactor.py run-core \
  --output output/formfactor-4l-core
python3 bench/four_loop_formfactor.py replay-final \
  --sewing-dir output/formfactor-4l-core \
  --output output/formfactor-4l-core-final
```

`run-core` uses an explicit allowlist, not “every historical receipt except a
few skips.” A newly encountered receipt is never executed by default. The
accepted sequence is:

1. FEC_3→FEC_6 with `bootstrap --extend`, starting from committed FEC_2;
2. LEC_2→LEC_4;
3. `galois`, `galois_suffix`, `actions_structural` for prefix actions, and the
   low-cost `actions` helper only for suffix actions through weight 2;
4. the original certificate-pinned `sew_affine_TR_build_v1` solver;
5. `affine_TR_extract OUT 6 SEW_6p2_affine_TR_matrix.wxf
   affine_tail_w8_invariant_TF175_origin.wxf tr`;
6. DC then TC exact affine boundary solves and the two reserved holdouts.

Every accepted output with a committed copy is compared byte-for-byte:
FEC_3/4/5, LEC_2, the affine-TR matrix, affine origin and TR92 directions. DC
and TC must both regenerate the candidate hash shown above. New runs write
`core_stage_receipts.jsonl` with each command, wall time and sampled peak RSS.

### Excluded paths (removed from normal replay, retained as provenance)

| Historical route | Status and why it is not replayed |
| --- | --- |
| `prefix_w6_staged` | Stopped prefix control; the completed streamed FEC_6 baseline is used. |
| `cyclic_LEC_2` | Failed with exit 1; its output has no consumer. |
| All-row `actions ... prefix 6` | Interrupted expensive control, about 2,767 s per generator in a completed comparison. Replaced by the certified `actions_structural` construction; outputs were independently byte-identical. The suffix-only `actions` call is still accepted. |
| `sew_w8_6p2` | Terminated with exit 143; no accepted output. |
| `sew_w8_masked_6p2` | Terminated with exit 143. The unfinished necessary-TF discovery was explicitly superseded by validated affine-TR coordinates. The separate `sew_masked` executable is never run. The internal `mode=masked` label inside `sew_affine_TR` is expected and is not this abandoned stage. |
| Symrep/isotypic preparation | No downstream consumer; archived under `experimental/symrep-isotypic-sewing/`. |

These routes are removed from the executable replay path, not deleted from the
historical record. Their receipts and sources remain available for audit.

### Measured accepted replay and parameter counts

The table reports the 2026-10-24 accepted replay. “Fixed” means independent
coordinates eliminated by that stage, not merely nonzero output entries.
Basis stages use the exact condition-matrix rank; action and validation stages
add no elimination. Memory is sampled process-tree RSS in decimal GB.

| Stage | Coordinates before → after | Parameters fixed | Wall seconds | Peak RSS, GB |
| --- | ---: | ---: | ---: | ---: |
| FEC_3 | 949 → 288 | 661 | 0.20 | <0.01 |
| FEC_4 | 6,919 → 1,548 | 5,371 | 0.20 | 0.01 |
| FEC_5 | 47,460 → 7,701 | 39,759 | 2.60 | 0.19 |
| FEC_6 | 276,267 → 35,598 | 240,669 | 100.10 | 6.27 |
| LEC_2 | 681 → 175 | 506 | 0.20 | <0.01 |
| LEC_3 | 6,974 → 1,506 | 5,468 | 0.20 | <0.01 |
| LEC_4 | 64,575 → 10,578 | 53,997 | 26.10 | 1.83 |
| `galois` | induced maps on FEC_1…6 | 0 | 6.40 | 4.211* |
| `galois_suffix` | induced maps on LEC_1…2 | 0 | 0.20 | 0.000217* |
| `actions_structural` prefix | induced maps at 35,598 prefix coordinates | 0 | 14.80 | 4.248* |
| `actions` suffix | induced maps at 175 suffix coordinates | 0 | 0.20 | 0.00322* |
| Affine-TR sewing | 3,275,017 augmented → 4 homogeneous (3 affine) | 3,275,013 | 4,564.90 | 10.09 |
| Affine extraction | 4 augmented → one origin + 3 directions | 0 | 3.60 | 0.45 |
| DC boundary | affine dimension 3 → 0 | 3 | 26.46 | 1.56 |
| TC boundary after DC | affine dimension 0 → 0 | 0 | 52.73 | 3.85 |
| DC reserved holdout | validation only; 7,901,300 components | 0 | 7.02 | 1.57 |
| TC reserved holdout | validation only; 8,697,622 components | 0 | 10.83 | 3.86 |

*Preparation peaks come from the immediate same-path preparation control
`output/formfactor-4l-prep-replay-20261024/prep_receipts.json`; every output
was byte-identical to the main core run. Main core walls are from
`output/formfactor-4l-accepted-core`; final-stage walls/peaks are from
`output/formfactor-4l-accepted-final/receipts/`.

The complete core took **4,719.870 s** and peaked at **10.09 GB**. The final
replay took **102.152 s** and peaked at **3.86 GB**. Total measured accepted
core-plus-final time was **4,822.022 s (80 min 22 s)**, ending with zero free
affine parameters and candidate SHA256
`5fd5b1d60c703bd303804307d467aec7b713213e278c923707db6580e0238f38`.

The 503 MB FEC_6 basis is regenerated, not committed. The `bootstrap` binary
is built with `make bootstrap`; the recovered FLINT runtime under
`output_hexagon_runtime_recovery_20261008/env/lib` is used automatically when
present. Boundary training and Wolfram-side H8 exports are not repeated by
this packaged replay: their committed maps and factors are explicit inputs.
The replay therefore covers the accepted packaged pipeline, not raw
re-derivation of every earlier research-stage constraint.

## Environment and machine parameters

The measured accepted replay used the following complete environment. Earlier
receipts used a 44 GB own-process guard; the 2026-10-24 replay used the runner's
46 GB core cap and 12 GB final-stage cap.

| Parameter | Value |
| --- | --- |
| Host | ThinkStation P3 Ultra SFF G2, NVMe storage |
| CPU | Intel Core Ultra 7 265; 20 cores, no SMT; 800–6,600 MHz reported range |
| Memory / swap | 66,681,028,608 bytes RAM (62.10 GiB); no swap |
| OS / kernel | Ubuntu 22.04.5 LTS; Linux 6.8.0-138-generic, x86_64 |
| Python | 3.10.12; standard-library runner |
| Compiler | conda-forge GCC/G++ 14.4.0 (`gxx_linux-64=14.4.0`) |
| Numerical libraries | FLINT 3.2.2, GMP 6.3.0, MPFR 4.2.2, TBB 2023.1.0, mimalloc 3.4.1 |
| Thread settings | `OMP_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1` |
| Memory guards | 46 GB decimal core cap; 12 GB final cap; 0.2 s RSS sampling |
| Runtime | `output_hexagon_runtime_recovery_20261008/env/lib` on `LD_LIBRARY_PATH` when present |
| Accepted sewing identity | Original binary SHA256 `2cfe28203e1d22c74c8e3b58cea50182e84903496fb9d02c279006c6960df644`; frozen source `solver-source/ff_four_loop/sew_affine_TR_build_v1/` |
| Fresh-build check | The frozen v1 source compiled and linked independently with GCC 14.4.0; no research-tree source access |

The repository [environment-linux.yml](../../environment-linux.yml) pins the
supported build environment. The C++ helpers compile from `solver-source/`;
receipt-pinned local binaries are checked before use. The optional historical
H8 factor exports used WolframScript 1.11.0 (Wolfram 14), but Wolfram is not
needed by any normal replay mode. No external dataset or research-tree input
is needed: every packaged input is committed and hash-bound in
[manifest.json](manifest.json).

## Toolchain

Beyond the environment above: Python 3 standard library for `verify`, and
`xz` (or Python's `lzma`, used by the runner) to decompress the two
projection maps. A licensed Wolfram kernel (`wolframscript`) is needed only
to regenerate the H8 boundary factor exports from
`solver-source/ff_four_loop_H8_*.wls`; the committed factor files let every
runner mode skip that step. The upstream research tree that produced the
result is kept local and is not part of this package; every input it
contributed is hashed in the manifest.

## What is regenerated versus committed

| Object | Policy |
| --- | --- |
| FEC_6 forward basis (503 MB) | regenerated by `run-core` |
| DC/TC weight-6 collinear projection maps | committed xz-compressed, decompressed by the runner |
| Affine-TR matrix, origin and TR92 directions | regenerated and byte-compared by `run-core`; committed copies retained |
| Candidate tensor, FEC_2..FEC_5, LEC_2, suffix bases | committed, hash-verified |
| All receipts, certificates, prime proofs, holdout targets | committed |
| Solver binaries | compiled from committed sources (or checked against receipt-pinned hashes) |

See [VALIDATION.md](VALIDATION.md) for the full validation record and the
measured resource envelope, and
[manifest.json](manifest.json) for the complete SHA256 inventory.
