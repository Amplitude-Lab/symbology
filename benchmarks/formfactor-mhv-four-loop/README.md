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
| `solver-source/` | All helper `.cpp`/`.py`/`.wls` sources of the accepted chain (39 C++ programs, drivers, audits, H8 factor build). `sew_masked` and `sew_affine_TR` have no static source: the packaged `build_masked_sew.py` / `build_affine_tr_sew.py` generate them by patching the repository headers, reproducing the recorded sources byte-for-byte (verified at packaging time); the runner invokes them automatically. The remaining `.py`/`.wls` drivers are provenance copies of the research-tree originals |

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
comparisons from the committed sewing outputs and boundary preparation:

```bash
python3 bench/four_loop_formfactor.py replay-final --output output/formfactor-4l-replay
```

The regenerated DC origin must be **byte-identical** to the committed
candidate, the TC solve (which consumes the regenerated DC solution) must
reproduce it again, and both holdout comparisons must return
`exact_boundary_equality` with the packaged component counts and certificate
moduli. Solver binaries are resolved from `--solver-bin DIR`, or from the
local research tree after checking them against the receipt-pinned SHA256
values, or compiled fresh from `solver-source/` (needs the toolchain below).
Every stage runs under a sampled memory cap (default 12 GB) and writes
receipts to the output directory.

## Regenerate the full core chain

`python3 bench/four_loop_formfactor.py run-core --output DIR` replays the
recorded native stage commands of `reference/stage_receipts.json`:

1. prefix bases FEC_3→FEC_6 (`bootstrap --extend`, starting from the
   committed hash-pinned FEC_2 basis),
2. `symrep` representation preparation from the committed generators,
3. the suffix chain LEC_2→LEC_4 and the cyclic LEC_2 projection,
4. the symmetry-action files (`galois`, `galois_suffix`, `actions` helpers
   from `solver-source`, deterministic from `data/`),
5. the weight-8 sewing: invariant, masked and affine-TR stages.

Every regenerable base with a committed copy — FEC_3, FEC_4, FEC_5, LEC_2 and
the affine-TR sewing matrix — is then compared byte-for-byte against
`reference/`. The 503 MB FEC_6 basis is regenerated, not committed. The
`bootstrap` and `symrep` binaries are built from the repository root sources
(`make bootstrap symrep` with the environment-linux.yml toolchain); a locally
prebuilt binary linking the recovered FLINT runtime under
`output_hexagon_runtime_recovery_20261008/env/lib` is picked up automatically
when present. Recorded cost of the accepted chain: **8,157 s (about 2.3 h)**
in total, dominated by the three sewing stages (1,069 / 1,910 / 4,445 s),
with sampled peaks at most 11.2 GB on adopted stages (the 11.9 GB spike in
the receipts is the staged control that was **not** adopted). One receipt
caveat: the 371.6 s `symrep` entry predates the final condition tensor (see
VALIDATION.md, "symrep timing anomaly"); with the committed tensor that
stage takes hours, and since its output is read by no downstream stage it is
auxiliary to the compared artifacts.

The boundary training that produced the two committed collinear projection
maps (`bootstrap --project --symmetry collinear --target FEC_6`, driver
`solver-source/ff_four_loop_boundary_run.py`) and the Wolfram-side H8 factor
exports (`solver-source/ff_four_loop_H8_*.wls`, a licensed `wolframscript`,
measured 4.9 s / 0.3 GB) can be regenerated the same way; the committed maps
let `replay-final` skip them.

## Environment

All recorded measurements ran on one workstation: Intel Core Ultra 7 265
(20 cores, no SMT), 62 GiB RAM, Ubuntu 22.04.5 LTS, kernel 6.8.0-138-generic.
The native stages use four threads (`OMP_NUM_THREADS=4` in the recorded
receipts) under the recorded 44 GB own-process cap; adopted stages peaked at
≤ 11.2 GB sampled RSS. The optional Wolfram-side H8 factor exports were
measured with WolframScript 1.11.0 (Wolfram 14) at 4.9 s / 0.3 GB. The
`bootstrap` and `symrep` binaries are built with the repository
[environment-linux.yml](../../environment-linux.yml) toolchain
(gxx_linux-64 14.4.0, libflint 3.2.2, gmp, tbb); the C++ helpers compile
from `solver-source/` with the same toolchain. No external dataset or
research-tree access is needed: every input is committed and hash-bound in
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
| Candidate tensor, FEC_2..FEC_5, LEC_2, suffix bases, sewing outputs | committed, hash-verified |
| All receipts, certificates, prime proofs, holdout targets | committed |
| Solver binaries | compiled from committed sources (or checked against receipt-pinned hashes) |

See [VALIDATION.md](VALIDATION.md) for the full validation record and the
measured resource envelope, and
[manifest.json](manifest.json) for the complete SHA256 inventory.
