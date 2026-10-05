# Five-loop MHV heptagon: one-day workstation trial

The user authorized a fresh attempt with one day of wall time and 50 GiB
of RAM. The run started on 4 October 2026 at 00:28:32 PDT; its fixed
deadline was 5 October 2026 at 00:28:32 PDT. **The five-loop MHV heptagon
symbol completed and passed every pipeline check** on 4 October at 02:38:28
PDT, after 7,796.38 seconds (2 h 9 min 56 s), with 19.40 GiB sampled peak
RAM. The follow-up heartbeat has been paused. This is a recursive symbol
result; function-level constants are outside this calculation.

The accepted result is
`output_fiveloop_day_20261004/amplitude5/hepMHV_5L_recursive.wxf`.
It requires the saved FEC1..8 and LEC1..2 bases and uses the **8+2** sewing
cut. It is not a standalone expansion in all original-alphabet words.

The small accepted tensor is now tracked in the
[published reproduction benchmark](../../benchmarks/heptagon-mhv-five-loop/README.md).
That package rebuilds every omitted prerequisite from tracked seeds and has
passed a complete cold reproduction with the current shared kernel library;
see its [validation record](../../benchmarks/heptagon-mhv-five-loop/VALIDATION.md).

| Completed stage | Wall seconds | Sampled peak RSS |
| --- | ---: | ---: |
| Exact invariant kernel | 7,341.009 | 19.399 GiB |
| Unique physical coefficients | 71.505 | 12.045 GiB |
| Complete collinear verification | 366.057 | 13.816 GiB |
| Independent published words | 17.192 | 11.899 GiB |
| Controller total, including output hashes | 7,796.382 | 19.399 GiB |

All stages returned zero with no time, RAM or host-memory guard triggered.
The exact invariant kernel has dimension four. Complete collinear
verification checked all **55,728,602 divergent words** in the computed
88,620,356-term collinear expression against the independently computed
lower-loop boundary. The four independent original-alphabet coefficients
are **−26,880, −26,880, +26,880, −1,680**, agreeing with equation (12) and
the following paragraph of
[He et al.](https://arxiv.org/html/2511.09669v2). These values were checked
after construction and were not fitting inputs. This does not claim a
coefficient-by-coefficient comparison with the entire published symbol.

The accepted recursive file has 6,825,548 bytes and SHA256
`a06256013cbb03a7a1b67a18f930e02639cf2e777d48196ae36ba00355661691`.
The verified collinear expression remains named `E5_candidate.wxf` to
preserve the driver's recorded output path; its certification is recorded
in `verify.log` and `status.json`. The final audit evidence is in
[`evidence.json`](evidence.json), including hashes and stage measurements.

The 2 h 10 min / 19.40 GiB result is the **new solve and verification run
from cached prerequisites**, not a cold rebuild of all lower weights.
It reuses FEC8 (previously 629.34 s / 32.67 GiB), the recursive actions
(70.18 s / 12.13 GiB), the independent boundary (492.64 s / 14.07 GiB), and
accepted lower-loop data. The prefix projection was rebuilt inside this
run and is included in the verification time. Thus the largest measured
peak across this run and its separately measured new prerequisites is
32.67 GiB, incurred when building FEC8; no fresh end-to-end cold benchmark
is claimed. All measurements are under CPU contention on this workstation.

## What produced the optimization

This experiment combines a change of mathematical formulation with general
sparse-solver engineering. It does not establish that an irreducible
representation basis alone produced the observed improvements.

1. Sew FEC8 with LEC2 directly, imposing local integrability and physical
   invariance together. This avoids first constructing unrestricted FEC9.
   The previous trial measured 82.55 GiB for that unrestricted equation
   payload alone. Changing the sewing cut does not require an irrep basis.
2. Apply the supplied cyclic, flip and parity actions to the product space
   early. The group is D7 × C2. Nontrivial representations on the two sides
   are retained because they can pair to produce invariant tensors. This
   is the symmetry-specific reduction, but the current five-loop driver
   uses generator constraints rather than an entirely decomposed irrep
   multiplicity solver.
3. Keep the large equation operators factored, propagate exact short
   relations, and build a sparse triangular system from selected equations.
   Avoiding the complete equation matrix is the principal storage change.
   Omitted discovery equations are mandatory in the subsequent complete
   factored residual intersection; sampling is not a correctness shortcut.
4. Use finite-field elimination and rational reconstruction with exact
   certification. Extra verification primes operate on the small candidate
   without repeating rank elimination solely to enlarge a certificate bound.
   Packed word accumulation and a compact projected right basis reduce
   separate boundary/projection costs. These improvements are independent
   of an irrep decomposition.

The preceding four-loop measurements reduced direct unrestricted 6+2
sewing from 13.89 GiB to about 0.50 GiB with the combined method. Those
measurements include multiple changes and CPU contention, so they are not
an ablation attributing the gain to symmetry. The completed five-loop
measurement above establishes feasibility, but not an isolated speed ratio
against an ordinary solver that did not complete this workload. See the
[preceding detailed audit](../early-sewing-2026-10-03/REPORT.md).

## Run and acceptance conditions

All paths below are relative to the repository. The run directory is
`output_fiveloop_day_20261004/`.

* `config.json` records the fixed shared deadline, limits and complete argv
  for each stage. `controller-launch.json` records the launch identity and
  configuration hash. The initial controller PID is 1052842; the initial
  solver PID is 1052844. Always check process command lines before acting
  on a recorded PID.
* `manifest.json` records the branch (`symrep-solver-dev`), base commit,
  executable/source snapshots and 63 input or cached-result hashes checked
  against the previous audit. Existing FEC8, generator actions, the complete
  five-loop boundary and the collinear prefix through FEC7 are reused.
  Collinear cache files that may be rewritten were copied into this run.
* `bin/` holds fixed executable and controller snapshots. The mathematical
  algorithm is the already validated `refine-hybrid` implementation. This
  attempt extends its budget; it introduces no new mathematical solver.
* `status.json` identifies the current stage. Each stage has a `.log`, a
  live JSON updated every ten seconds, and a final measurement JSON. The
  RSS measurement samples the sum over the process tree every 100 ms;
  shared pages can be counted more than once. User/system CPU times are
  recorded separately from wall time.
* `monitor.json` identifies the half-hour thread heartbeat. It reports
  meaningful progress or problems and stops monitoring after completion or
  expiration of the original budget. The controller advances stages even
  between heartbeat checks.

The sequential stages share a single 86,400-second allowance:

1. Exact kernel: require `early/EARLY_8p2.wxf`, its saved LEC2 basis and the
   complete factored-residue/height certificate in `kernel.log`.
2. Coefficients: determine a unique candidate using small projections.
   Files named `candidate_*` remain uncertified physical outputs.
3. Verification: test the complete collinear divergent part against the
   independent lower-loop boundary before writing
   `amplitude5/hepMHV_5L_recursive.wxf`.
4. Independent words: check the four published original-alphabet word
   coefficients. These values are validation data, not fitting inputs.

The final recursive symbol, if accepted, uses the **8+2** cut and its saved
FEC1..8/LEC1..2 bases. It is a symbol-level result, not a determination of
function-level constants.

Each stage has a hard 50 GiB per-process address-space cap, a sampled
50 GiB aggregate process-tree RSS stop, and a host guard that stops this
run if available RAM drops below 4 GiB. The RSS guard is sampled, not a
kernel-enforced aggregate memory limit; address space can also be exhausted
before RSS reaches 50 GiB. The machine has about 62 GiB physical RAM and no
swap. An unrelated CPU-intensive user job remains untouched, so timings
are contended. The shared deadline terminates the owned stage process group
and prevents later stages from starting. It is not 24 hours per stage.

## Validation before launch

`python3 tests/process_budget_probe.py` passed all seven tests: live status,
timeout, descendant memory accounting and termination, sequential accepted
outputs, rejecting missing certificates, a shared deadline, and cancellation
propagation. Evidence: `output_fiveloop_day_20261004/budget-tests.log`.

The fixed executable snapshots then ran the full four-loop pipeline:

| Stage | Wall seconds | Sampled peak RSS |
| --- | ---: | ---: |
| Exact hybrid kernel | 17.493 | 0.500 GiB |
| Unique coefficients | 0.305 | 96.3 MiB |
| Complete collinear verification | 2.211 | 208.9 MiB |
| Independent published words | 0.202 | 95.7 MiB |

The controller completed in 20.335 seconds. The kernel has dimension three
and passed exact certification. Full collinear verification and all four
published word checks passed. The complete E4 output is byte-identical to
`output_nmhv_kernel_20261002/default-fourloop/4loop/E4.wxf`. Evidence and
output hashes are in `output_fiveloop_day_20261004/control4/`.

The five-loop process was verified running after launch. Its early observed
14.24 GiB peak is only a partial-run measurement, not its eventual peak.
Consult the live state and completed stage reports for newer evidence.

## First extended-run milestone

At approximately 01:33 PDT on 4 October, the run had passed the previous
forty-minute trial's stopping point and entered discovery using unrestricted
row lengths with a 1/16 selection fraction. The last inspected rank was
456,346 of 478,037 (21,691 unresolved directions), compared with 377,212
at the end of the complete 1,024-term band. Sampled peak RSS was 17.20 GiB
after about 64 minutes. This is intermediate finite-field discovery, not a
certified kernel. The active process identities and original limits were
verified; the computation was left running unchanged. Detailed observations
are appended to `heartbeat-observations.jsonl` in the run directory.

## Exact kernel and unique candidate completed

At approximately 02:35 PDT on 4 October, the controller had completed the
kernel and coefficient stages and automatically started full verification.
The completed kernel stage took 7,341.009 seconds (2 h 2 min 21 s), with
19.399 GiB sampled peak RSS. Its exact invariant dimension is **four**, with
684,334 nonzero rational coefficients in the FEC8 × LEC2 representation.
The saved kernel is `early/EARLY_8p2.wxf` (4,565,763 bytes), SHA256
`86d9731ff62aa6b39495fb3a6b3b53549442475e7551177cce0abbe4e47a32e0`.

Discovery stopped at 128 remaining directions after consuming 8,262,891
equations with 6,813,101,616 input nonzeros across the streamed equations.
These nonzeros were not stored simultaneously. The triangular system retained
211,386,644 nonzeros. Complete factored residual intersection then reduced
128 → 5 → 4. The final modular rank is 478,033 of 478,037 reduced columns.
Rational reconstruction succeeded at the first 61-bit prime. Because the
rigorous residual bound required 76 bits, a second good prime checked the
small rational candidate; the combined verified modulus has 121 bits.
`exact_certificate=complete_factored_residues_plus_height` was recorded.

Coefficient discovery took 71.505 seconds and peaked at 12.045 GiB. The
two-letter projection did not determine a unique solution; the three-letter
projection did. A candidate recursive symbol was saved, with SHA256
`a06256013cbb03a7a1b67a18f930e02639cf2e777d48196ae36ba00355661691`.
These projected equations alone did not certify the physical amplitude.
The subsequent complete collinear divergent-part verification and the
independent published word checks both passed, as recorded above.

The copied prefix projection cache was rebuilt automatically by the public
projection command in the new output directory (about 146 seconds), rather
than silently accepting receipts associated with the old location. This
cost is included in the new full-verification stage. Old accepted files
were not overwritten.
