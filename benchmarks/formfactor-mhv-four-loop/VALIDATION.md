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
`environment-linux.yml` toolchain and `make bootstrap symrep`.

## Sewing-stage input closure (2026-10-23)

A procedure-fidelity audit (prompted by re-examination of the recorded
commands) confirmed that `sew_masked` additionally reads
`boundary_local_rref.wxf` and `sew_affine_TR` reads
`boundary_local_TR_rref.wxf`, `affine_TR_rhs_w8_residual.wxf` and the seven
`actions/suffix_TR_*_w2.wxf` variants from the experiment root, none of
which are products of the recorded stage commands. These files (2.1 MB +
0.8 MB + 55 MB + 19 KB) are now committed as staged inputs and seeded by
`run-core`; their producers (`local_TR_gate.cpp`, `suffix_TR_embedding.cpp`)
are committed in `solver-source` for regeneration. Also verified: the two
`symrep prepare` variants write to different directories
(`symrep_prepared`/`symrep_cached_prepared`), neither is read by any
downstream stage (the sewing consumes only `data/`, the bases and the
action files), and skipping the cached variant and the stopped staged
control matches the accepted chain.

## symrep timing anomaly: root cause (2026-10-23)

The recorded `symrep_prepare` receipt (371.6 s wall, 370.97 s CPU, flat
20 MB RSS) is not reproducible with the committed inputs: a verbatim replay
of the exact recorded command, with hash-verified identical binary, data,
runtime libraries and environment, exceeded the recorded total CPU within
~11 minutes and both observed instances compute steadily at one full core
(the probe's core verified at 5.1 GHz boost, exonerating frequency
throttling; the runtime libraries pass their own recovery-validation gate).

The file chronology of the recording day explains it. `interior.wxf` and
`FEC_2.wxf` were regenerated at Oct 10 19:46:23.85 (the adjacency condition
merge; `adjacency_symmetric.json` written at 19:46:01, superseding the
05:39-era condition tensor). The recorded symrep run started ~19:46:19 —
four seconds **before** the rewrite — so it computed against the older,
smaller pre-merge tensor; every downstream stage (prefix bases 19:46:52→
19:48:25, suffix chain 20:03, sewing 20:16→22:19) ran against the final
tensor. `stage_receipts.json` is an append-only log and carries no
timestamps, so the 371.6 s entry is a receipt of the earlier input
generation. The `symrep_prepared` output is read by no downstream stage
(the sewing consumes only `data/`, the bases and the action files), so the
divergence has no effect on artifact identity — which the byte-identical
regeneration of FEC_3/FEC_4/FEC_5 and the exact-size FEC_6 independently
confirm. Consequence for `run-core`: with the committed final tensor the
symrep stage legitimately takes hours (measured: 3+ h CPU, still running at
the time of this entry); its recorded 371.6 s should not be used as a cost
estimate, and the stage is auxiliary to the compared artifacts.

## Solver path nondeterminism observed with byte-identical inputs (2026-10-23)

A rerun of the recorded `sew_w8_6p2` stage with hash-verified identical
binary, `interior.wxf`, FEC_6 (sha d116600a... — the candidate-manifest
pin), LEC_2 and all ten action files (each compared byte-for-byte against
the originals) entered a multi-pass presolve trajectory the recorded run
never printed: the original completed in 1,069 s after the interface
reduction with no presolve output, while the rerun was still reducing
~1.25 M remaining columns after six presolve passes and 2+ h. The presolve
skips rows against a relation state that evolves during parallel emission,
so scheduling order can change the trajectory; the exact kernel, when
reached, is unchanged. Consequences recorded honestly: (1) the receipt
timings are single samples of a path-nondeterministic solver, not strict
upper bounds; (2) run-to-run wall time can vary by orders of magnitude
under identical inputs; (3) this finding also weakens the earlier
"pre-merge fossil" attribution for the symrep anomaly to one of two
consistent explanations (input-generation difference, as the file
chronology suggests, or the same path nondeterminism). A rerun taking the
recorded fast path would reproduce the original timings; the final
byte-comparison gate remains the arbiter of correctness either way.
