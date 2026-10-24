# Validation record — four-loop four-point form factor

This file records certificates, accepted-route provenance, excluded historical
routes and independent checks. Operational instructions and measured
performance are in [README.md](README.md).

## Result and mathematical scope

`reference/F4_candidate_prefix6_suffix2.wxf` has SHA256
`5fd5b1d60c703bd303804307d467aec7b713213e278c923707db6580e0238f38`,
22,457,737 bytes, shape `1 × 35598 × 175`, and 988,559 nonzero coefficients.
It decodes as `S4 = Σ C[1,i,a] P6_i ⊗ Q2_a` using the committed prefix
recursion and `suffix_TF175_physical`; all basis hashes are in
`reference/candidate_manifest.json`.

Uniqueness is conditional on the proposed-adjacency and continued TR12
hypotheses at symbol level. External Steinmann equations were not imposed.
There is no all-loop claim, larger-ES uniqueness claim, full-parent comparison,
or function-level statement. Finite four-loop targets were reserved holdouts,
not solving inputs.

## Chain of evidence

1. **Affine-TR sewing** — `reference/affine_TR_sewing_certificate.json`:
   augmented kernel dimension 4; `modular_rank + 4 == presolved_columns`;
   residue product 121 bits above the 91-bit integer residual bound; finite
   boundary targets never read; all seven hashed inputs verify.
2. **Ambiguity structure** — `validation/ambiguity/`: the three directions are
   `span_Q{a⁸+b⁸, a⁶b²+a²b⁶, a⁴b⁴}` with `a=log(u1/u3)` and
   `b=log(u2/u4)`, factorized and certified.
3. **DC affine solve** — `validation/affine-solve/DC/`: affine dimension 3→0
   through 141 nonzero pair blocks and 343 zero-input blocks, covering all
   484 pairs and 2,273,673 full-coordinate components. Every exact kernel and
   affine-RHS coordinate transport is checked.
4. **TC affine solve after DC** — `validation/affine-solve/TC/`: affine
   dimension 0 retained through 203 nonzero and 638 zero-input blocks,
   covering all 841 pairs and 2,089,544 components; the origin is a
   byte-identical copy of the DC candidate.
5. **Reserved holdouts** — `validation/holdouts/`: DC and TC finite images
   equal the published four-loop three-point form-factor and six-point
   amplitude remainders exactly. They certify 7,901,300 and 8,697,622
   components over two independently proven primes with products above the
   integer bounds. Every solve receipt records `finite_targets_read=false`.
6. **Independent audits** — `reference/`: extraction, final coordinate
   transport and consolidated evidence-chain audits all record
   `all_checks_pass`. The runner also replays every Lucas prime certificate
   rather than trusting the stated factorizations.

## Accepted replay and exclusions

The normal runner has an explicit accepted-stage allowlist. Unknown receipts
are never executed by default. The accepted route is FEC/LEC bases, Galois
helpers, `actions_structural` prefix actions, suffix `actions`, original
`sew_affine_TR_build_v1`, `affine_TR_extract`, DC, TC, then reserved holdouts.

The following are not accepted stages and cannot be selected by normal replay.
Their files remain only as historical or experimental evidence.

| Excluded route | Verified status |
| --- | --- |
| `prefix_w6_staged` | Stopped control; completed streamed FEC_6 is the accepted baseline. |
| `cyclic_LEC_2` | Exit 1; output has no downstream consumer. |
| All-row `actions ... prefix 6` | Interrupted expensive control. The accepted structural construction is independently byte-identical to the completed all-row outputs. Suffix-only `actions` remains accepted. |
| `sew_w8_6p2` | Exit 143; no accepted output. Its wall time is a kill time, not a completed-solve cost. |
| `sew_w8_masked_6p2` | Exit 143; unfinished necessary-TF discovery superseded by validated affine-TR coordinates. The separate `sew_masked` executable is never used. The internal `mode=masked` label inside `sew_affine_TR` is expected. |
| Symrep/isotypic preparation | No downstream consumer; archived under `experimental/symrep-isotypic-sewing/`. |

No nondeterminism, completed plain/masked sewing, or accepted isotypic stage is
claimed. Terminated attempts are excluded from accepted-chain cost.

## Accepted legacy source and inputs

The accepted sewing certificate pins `sew_affine_TR_build_v1` and executable
SHA256 `2cfe28203e1d22c74c8e3b58cea50182e84903496fb9d02c279006c6960df644`.
The complete original source, headers, generator and provenance manifest are
frozen under `solver-source/ff_four_loop/sew_affine_TR_build_v1/`. Fresh builds
use this copy, not current repository headers or the later guarded successor.
The frozen source compiled independently with GCC 14.4.0.

Accepted sewing reads committed `boundary_local_TR_rref.wxf`,
`affine_TR_rhs_w8_residual.wxf` and seven `actions/suffix_TR_*_w2.wxf` files.
Extraction additionally reads `affine_tail_w8_invariant_TF175_origin.wxf` and
`suffix_TR92_to_TF175.wxf`. Their producers are committed in `solver-source`.
The masked-only `boundary_local_rref.wxf` is retained as experimental evidence
and is not seeded into the accepted run.

The two DC/TC collinear projection maps are committed compressed inputs. The
replay does not claim to regenerate them or earlier research-stage constraint
preparations. No external dataset or research-tree input is needed by the
packaged commands.

## Replayed gates and measured acceptance

| Gate | Result |
| --- | --- |
| `verify` | All 207 packaged files, receipts, hashes, block chains, prime certificates and guards pass. |
| `run-core` | All seven comparable accepted artifacts regenerate byte-identically: FEC_3/4/5, LEC_2, affine-TR matrix, affine origin and TR92 directions. |
| `replay-final` | DC and TC regenerate the candidate byte-identically and both reserved holdouts pass with the packaged component counts and moduli. |

The fresh 2026-10-24 core completed in 4,719.870 s at 10.09 GB sampled peak
RSS. Structural prefix actions took 14.8 s, affine-TR sewing 4,564.9 s, and
extraction 3.6 s. `core_summary.json` records `all_pass=true`; new runs also
write `core_stage_receipts.jsonl`. Exact per-stage ranks, eliminated-parameter
counts, wall times and memory are tabulated in [README.md](README.md).

The final replay consumed those regenerated sewing outputs. DC took 26.46 s,
TC 52.73 s, and the complete final replay 102.152 s. Both reserved holdouts
passed. Accepted core plus final took 4,822.022 s. Receipts are preserved at:

- `output/formfactor-4l-accepted-core/core_summary.json`;
- `output/formfactor-4l-accepted-final/receipts/summary.json`;
- `output/formfactor-4l-prep-replay-20261024/prep_receipts.json`.

The interrupted all-row retry is preserved separately at
`output/formfactor-4l-allrow-control-incomplete`; it is not an accepted result
or timing measurement.
