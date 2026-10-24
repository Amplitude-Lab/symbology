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
| `replay-final` | ~2 min, ≤ 12 GB cap | the DC and TC exact affine solves re-executed from committed inputs regenerate the candidate **byte-identically**; both holdout comparisons re-pass with the recorded component counts and moduli |
| `run-core` | 8,157 s recorded (~2.3 h), ≤ 11.2 GB peak, 503 MB FEC_6 regenerated | the recorded native stage commands rebuild every committed base bit-for-bit (FEC_3, FEC_4, FEC_5, LEC_2, affine-TR sewing matrix) |

Measured envelope of the accepted chain (from the packaged receipts):
prefix chain 94 s total on four threads; DC solve 25.8 s; TC solve 52.2 s at
3.8 GB; holdout comparisons 6.0 s and 6.1 s at ≤ 1.6 GB own RSS. The heavier
boundary-training projections and H8 Wolfram exports (minutes at ≤ 12 GB) are
bypassed by the committed maps; their receipts and sources are included.

## Not claimed

- No uniqueness in the larger extended-Steinmann space (external Steinmann
  equations were not imposed).
- No all-loop statement; the 96-letter enlargement agrees through weight 8
  but all-loop 96→88 exclusion is open.
- No comparison against a published full four-loop parent symbol (none was
  available); the two finite descendants are the published cross-checks.
- `run-core` replays the adopted streamed baseline; the staged control that
  measured 11.9 GB was stopped and is not part of the accepted chain.
- No function-level (beyond-symbol) or full-kinematic-parent claim.

## Packaging self-checks (2026-10-23)

The packaged `build_masked_sew.py` and `build_affine_tr_sew.py` were made
location-independent and re-executed from `solver-source/`: the regenerated
`sew_masked.cpp`, `sew_affine_TR.cpp` and patched `tensor_kernel.hpp` are
byte-identical to the recorded research-tree build products. The runner
invokes these generators automatically when no prebuilt sew helper is
available, so `run-core` needs only the repository, the
`environment-linux.yml` toolchain and `make bootstrap`.

## Sewing-stage input closure (2026-10-23)

A procedure-fidelity audit (prompted by re-examination of the recorded
commands) confirmed that `sew_masked` additionally reads
`boundary_local_rref.wxf` and `sew_affine_TR` reads
`boundary_local_TR_rref.wxf`, `affine_TR_rhs_w8_residual.wxf` and the seven
`actions/suffix_TR_*_w2.wxf` variants from the experiment root, none of
which are products of the recorded stage commands. These files (2.1 MB +
0.8 MB + 55 MB + 19 KB) are now committed as staged inputs and seeded by
`run-core`; their producers (`local_TR_gate.cpp`, `suffix_TR_embedding.cpp`)
are committed in `solver-source` for regeneration.

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
