# Validation record — four-point form factor four-loop bootstrap

This file records what is verified, how, and with which measured resources.
Numbers are taken from the packaged receipts; the runner re-derives the
final stages from the committed inputs.

## Result being validated

`reference/F4_candidate_prefix6_suffix2.wxf`, SHA256
`5fd5b1d60c703bd303804307d467aec7b713213e278c923707db6580e0238f38`,
22,457,737 bytes, shape `1 × 35598 × 175`, 988,559 nonzero coefficients,
decoding `S4 = Σ C[1,i,a] P6_i ⊗ Q2_a` with the committed FEC_1..FEC_6 prefix
recursion and `suffix_TF175_physical` (see
`reference/candidate_manifest.json` for every basis hash).

## Chain of evidence (all committed)

1. **Sewing** (`reference/affine_TR_sewing_certificate.json`): complete
   augmented kernel dimension 4; `modular_rank + 4 == presolved_columns`;
   certified residue product 121 bits > 91-bit integer bound; finite
   boundary targets never read. All seven hashed inputs verify.
2. **Ambiguity structure** (`validation/ambiguity/`): the three remaining
   directions are exactly `span_Q{a⁸+b⁸, a⁶b²+a²b⁶, a⁴b⁴}` with
   `a = log(u1/u3)`, `b = log(u2/u4)`, factorized and certified.
3. **DC boundary solve** (`validation/affine-solve/DC/`): affine dimension
   3 → 0 through 141 nonzero pair blocks (343 verified zero-input blocks,
   484 = 22² total), 2,273,673 full-coordinate components, all exact
   kernels verified, coordinate changes transport the affine RHS.
   Recorded 25.8 s wall.
4. **TC boundary solve** (consuming the DC solution): affine dimension 0
   retained; 203 nonzero blocks, 638 zero-input blocks (841 = 29²),
   2,089,544 components. Recorded 52.2 s wall, 3.8 GB sampled peak RSS,
   byte-identical origin copy.
5. **Reserved holdouts** (`validation/holdouts/`): the DC and TC finite
   four-loop images equal the published four-loop three-point form-factor
   and six-point amplitude remainders exactly — 7,901,300 and 8,697,622
   components certified over two proven primes with product exceeding the
   integer bounds (10,051,647,050,603,419,729,920 and
   11,397,249,185,853,671,345,356,80). The targets did not enter the
   solving (`finite_targets_read` false in every solve receipt).
6. **Independent audits** (`reference/`): coordinate extraction, final
   coordinate transport, and the consolidated evidence chain all record
   `all_checks_pass`, with the boundary-block receipts, complete prime
   logs, heights and hashes checked by independent Python rational replay.

The two accepted primes are re-proved from the recorded Lucas order
certificates (every product, recursive hypothesis, modular power and gcd
re-checked; factorizations are not trusted).

## Gates provided by the runner

| Gate | Cost | What it proves |
| --- | --- | --- |
| `verify` | seconds | every packaged file matches the manifest; all receipts, block chains, prime proofs and guard hashes are internally consistent; the candidate equals the certificate identity |
| `replay-final` | ~2 min, ≤ 12 GB cap | the DC and TC exact affine solves re-executed from selected committed or regenerated sewing outputs regenerate the candidate **byte-identically**; both holdout comparisons re-pass with the recorded component counts and moduli |
| `run-core` | 4,564 s recorded (~76 min) plus action preparation/extraction; 10.2 GB sampled sewing peak | the accepted native stages rebuild FEC_3, FEC_4, FEC_5, LEC_2, the affine-TR matrix, affine origin and TR92 directions bit-for-bit; terminated plain/masked sewing experiments are not run |

Measured envelope of the accepted chain (from the packaged receipts):
prefix and suffix stages 118.7 s total on four threads; affine-TR sewing
4,445.3 s at 10.2 GB sampled peak; extraction 3.5 s CPU at 442,056 KiB maximum
RSS; DC solve 25.8 s; TC solve 52.2 s at 3.8 GB; holdout comparisons 6.0 s and
6.1 s at ≤ 1.6 GB own RSS. Accepted structural prefix-action construction used
13.30 s user CPU, 3.08 s system CPU and 4,132,904 KiB maximum RSS; the
interrupted all-row control is not part of the accepted path. The heavier
boundary-training projections and H8 Wolfram exports are bypassed by
the committed maps; their receipts and sources are included.

## Not claimed

- No uniqueness in the larger extended-Steinmann space (external Steinmann
  equations were not imposed).
- No all-loop statement; the 96-letter enlargement agrees through weight 8
  but all-loop 96→88 exclusion is open.
- No comparison against a published full four-loop parent symbol (none was
  available); the two finite descendants are the published cross-checks.
- `run-core` replays the adopted affine-TR baseline; the staged prefix control
  and terminated plain/masked sewing attempts are not part of the accepted chain.
- No function-level (beyond-symbol) or full-kinematic-parent claim.

## Accepted legacy source recovery (2026-10-24)

The accepted sewing certificate pins `sew_affine_TR_build_v1` and executable
SHA256 `2cfe28203e1d22c74c8e3b58cea50182e84903496fb9d02c279006c6960df644`.
The complete v1 source, headers, generator and provenance manifest were found
in the research tree and are now frozen under
`solver-source/ff_four_loop/sew_affine_TR_build_v1/`; fresh builds use that
copy, not current repository headers. The later `build_affine_tr_sew.py`
guarded successor and `build_masked_sew.py` remain as provenance only. The
original 2026-10-23 generator self-check concerned those later build products,
not the certificate-pinned v1 binary.

## Sewing-stage input closure (2026-10-23)

A procedure-fidelity audit (prompted by re-examination of the recorded
commands) confirmed that the abandoned `sew_masked` reads
`boundary_local_rref.wxf`, while accepted `sew_affine_TR` reads
`boundary_local_TR_rref.wxf`, `affine_TR_rhs_w8_residual.wxf` and the seven
`actions/suffix_TR_*_w2.wxf` variants from the experiment root, none of
which are products of the recorded stage commands. The accepted inputs and
extraction inputs are committed and seeded by `run-core`; their producers
(`local_TR_gate.cpp`, `suffix_TR_embedding.cpp`) are committed in
`solver-source` for regeneration. The masked-only local RREF remains committed
as experimental provenance but is not needed by the accepted run.

## Archived: isotypic sewing experiment (2026-10-23)

The auxiliary `symrep prepare` stages were removed from the recorded
stage list: the chain uses no irreducible-representation decomposition
anywhere (the FEC/LEC bases are plain condition-tensor kernel recursions,
and the sewing imposes the seven symmetry generators as explicit
invariance rows). The killed/failed receipts, the timing anomaly history
(the 371.6 s receipt predates the adjacency condition merge) and the
experiment sources are archived under
`experimental/symrep-isotypic-sewing/` as a documented alternative
approach to optimizing the sewing.

## Sew-stage trajectories are deterministic; nondeterminism claim retracted (2026-10-23)

This file earlier recorded "solver path nondeterminism observed with
byte-identical inputs" for `sew_w8_6p2`. That claim is retracted: it was
written from a single rerun, and every same-day rerun — six launches in
total, the first being the rerun the retracted claim was based on; one
(`ff-exp-taskset`) ran with restricted CPU placement and ~4x slower
timings — reproduced **bit-identical trajectories**. Every rerun of the
invariant sewing prints the same numbers — pass 1 `relations=4879894
remaining_columns=1349756`, and two independent runs agree through pass 8
(`relations=4 remaining_columns=1255224`, `REFINE_SEED rows=3473367
cols=1255224 nnz=58523063`) with only timing strings
(`elapsed=`/`setup_seconds=`) differing. The masked
sewing reproduces the same way: the 2026-10-23 `run-core` rerun repeated
the original recording's eight passes number-for-number
(`4883699/1345951 … 57/1250698`) and its streaming batches, ranks and nnz.
No run-to-run divergence has been observed in any log.

The recording itself was also misread. The receipts show `exit_code 143`
for both `sew_w8_6p2` (wall 1,068.96 s) and `sew_w8_masked_6p2`
(1,909.70 s): both recorded runs were terminated by signal 15, and
`SEW_6p2_invariant.wxf` was never written — the invariant log stops at
`sew_strategy=staged product_generators=7` with no presolve output. The
earlier statement that "the original completed in 1,069 s" is wrong:
these walls are kill times, not completions, and are not cost estimates.
Only `sew_w8_affine_TR_6p2` completed (exit 0, 4,445.31 s).

With determinism established, the same-day eliminations stand: inputs
(byte-identical interior/FEC_6/LEC_2/actions), argv, library runtime,
environment, scheduling/serialization and timing-based switching were all
excluded, so the recording was not produced by the current binary and
inputs under any of the exercised conditions (concurrent load, solo
relaunch, restricted-core placement). The standing (unproven)
hypothesis is that the binary that produced the Oct-10 receipts differed
from the on-disk `bootstrap` (untracked; receipts pin no binary hash).
`run-core` skips the stage on the no-consumers ground (99538b1b); the
byte-comparison gate remains the arbiter. This retraction also removes
the alternative explanation previously offered for the auxiliary symrep
timing anomaly (archived with that experiment under
`experimental/symrep-isotypic-sewing/`), leaving the pre-merge file
chronology as its sole recorded attribution.

## Accepted-path cold replay (2026-10-24)

A fresh `run-core` used the explicit accepted-stage allowlist and frozen v1
solver, with no plain/masked sewing or all-row prefix-action control. Structural
prefix actions completed in 14.8 s; affine-TR sewing completed in 4,564.9 s at
10.09 GB sampled peak RSS; extraction completed in 3.6 s. The regenerated
FEC_3, FEC_4, FEC_5, LEC_2, affine-TR matrix, affine origin and TR92 directions
are all byte-identical to `reference/`. See
`output/formfactor-4l-accepted-core/core_summary.json`. The complete core took
4,719.870 s. Exact per-stage ranks, fixed-parameter counts, wall times and peak
RSS are tabulated in [README.md](README.md); preparation-stage memory was
independently measured in `output/formfactor-4l-prep-replay-20261024/`.

`replay-final --sewing-dir output/formfactor-4l-accepted-core` then consumed
these regenerated outputs. DC (26.5 s) and TC (52.7 s) both regenerated the
accepted candidate byte-identically, and the reserved DC/TC holdouts certified
exact equality for 7,901,300 and 8,697,622 components. See
`output/formfactor-4l-accepted-final/receipts/summary.json`. The final replay
took 102.152 s, giving 4,822.022 s for accepted core plus final. This is a cold
replay of the packaged accepted path, not regeneration of the committed
boundary-preparation maps or earlier research-stage constraint inputs. The
interrupted all-row retry is preserved separately at
`output/formfactor-4l-allrow-control-incomplete`.
