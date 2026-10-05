# Early invariant sewing: five-loop trial

Bounded trial on `symrep-solver-dev`, 3 October 2026. This audit concerns
MHV symbols, not function-level constants. **The five-loop bootstrap did not
complete.** No certified weight-ten kernel or five-loop amplitude was saved.
The completed four-loop controls and independent five-loop boundary are
useful results, but are not a five-loop solution.

| Stage | Wall time | Peak RSS | Outcome |
| --- | ---: | ---: | --- |
| Recursive FEC8 generator actions | 70.18 s | 12.13 GiB | Complete |
| Full five-loop collinear boundary | 492.64 s | 14.07 GiB | Complete |
| Collinear prefix projection through FEC7 | 164.85 s | 5.64 GiB | Complete |
| Latest five-loop hybrid kernel trial | 2,400.59 s | 14.24 GiB | Time limit; no certified kernel |

These are individual stage measurements under CPU contention. Earlier
exploratory runs also failed to produce a five-loop kernel; their costs are
recorded below, not included in the latest trial's forty-minute allowance.
The prerequisite FEC8 basis was reused from the preceding trial, where its
construction took 629.34 s and peaked at 32.67 GiB; see the
[earlier feasibility audit](../fiveloop-2026-10-03/REPORT.md).
The successful stages fit this machine. The peak memory and time required
for a **completed** five-loop solve remain unestablished.

## Formulation

Use FEC8 × LEC2 for weight ten, avoiding unrestricted FEC9. The general
product-invariance source takes the supplied left and right generator
matrices and imposes `X (A_g tensor B_g) = X`. It keeps nontrivial
representations on both sides: their pairings can produce invariants.
The heptagon trial uses cyclic, flip and parity, generating D7 × C2.

Exact one- and two-term presolve is followed by a finite-field sparse front.
The front initially selects equations with at most 64 terms. An incremental
triangular basis processes successively larger equations. Singleton zeros
and doubleton relations found during elimination are propagated back into
equation generation. Discovery does not assume an expected physical nullity.

Once the candidate kernel is small, the complete original factored operators
are applied to it. New residual conditions restrict this kernel; all earlier
conditions remain satisfied under restriction. This checks every original
condition, including every row omitted during discovery. The resulting
finite-field row space feeds the existing CRT/rational reconstruction and
rigorous residual-height certificate. Modular rank and independent private
coordinates certify completeness over Q. Equation selection alone is never
accepted as a certificate.

The `refine-sample` experiment additionally remembers consumed original
equations and selects a nested, deterministic fraction of the remaining
rows: initially 1/16, then 1/4 and finally all rows if necessary. The cheap
64-term front remains complete. Sampling is only a discovery heuristic;
the complete factored intersection described above is unchanged. On the
four-loop control this processes 770,161 rows rather than 1,524,437 and
retains 296,130 triangular nonzeros rather than 543,916. These structural
reductions do not depend on the contended wall-clock measurements.

That improvement did **not** scale directly to weight ten. Sampling the
cheap bands skipped useful sparse pivots: at rank 367,601 it retained
364,232,619 nonzeros, versus roughly 40 million at rank 377,251 with complete
short bands. The sampled five-loop experiment was stopped after 773.87 s,
at 15.93 GiB peak RSS, without a kernel. This is a failed exploratory run,
not evidence that sampling is universally faster. The `refine-hybrid`
experiment instead completes bands through 1,024 terms before sampling
larger equations. It caches consumed rows to avoid the earlier repeated
reductions. Its four-loop exact space comparison passes (19.274 s,
0.505 GiB). The five-loop hybrid trial reached its 2,400-second limit
without a rational candidate. It peaked at 14.24 GiB, below its 32 GiB
address-space limit, and its host-memory guard did not fire.
The last completed band had rank 377,212 of 478,037, leaving 100,825
directions. It retained 41,807,117 nonzeros, after processing over 7.37
million discovery rows. The last two logged blocks (262,144 input rows)
added no rank. It had entered the 4,096-term band when stopped; its final
partial-band rank was not logged. Evidence: `early8-hybrid.{json,log}`.

A separate four-loop control enlarged the initial globally eliminated front
from 64 to 512 terms. Its exact space comparison also passes, but cost rises
to 28.290 s and 2.254 GiB, compared with 14.662 s and 0.511 GiB for the small
sampled front. It is not adopted as a default, and that control alone does
not establish its performance at five loops. Evidence: `early6-front512.*`
and `fourloop-front512-compare.log`.

Reconstructed rational candidates are now checkpointed before certification.
If their rigorous residual bound exceeds the CRT reconstruction modulus,
the verifier checks the small candidate at extra good primes, without
repeating rank elimination. Acceptance requires the product of all verified
moduli to exceed the integer residual bound. A nonzero residual rejects the
candidate and resumes ordinary CRT reconstruction. Tests include a false
candidate congruent to a kernel at the first prime, and a symmetry action
whose denominator makes the next prime unusable. Files named `CANDIDATE_*`
are explicitly uncertified, including their saved congruence metadata.

The factored check contracts the candidate into the recurrence before the
local conditions. It does not form the large equation matrix. Symmetry is
also checked as a product of the two small factors, never as a full Kronecker
matrix. Large coefficient denominators are cleared per output equation for
the symmetry height bound, rather than by a global least common multiple.

## Four-loop control

All measurements below are contended: an unrelated, approximately 16-core
workload and other user work remained active. Some development checks and
low-priority compiler jobs also ran alongside exploratory five-loop trials.
These are feasibility measurements, not isolated-machine benchmarks.

* Direct unrestricted 6+2 sewing: 115.92 s, 13.89 GiB sampled peak RSS,
  four-dimensional space. Its equation input has 669,677,041 nonzeros.
* Early symmetry leaves 48,566 columns after exact presolve, versus 89,588
  for the unrestricted 6+2 system. Its full equation matrix would still have
  500,348,634 nonzeros.
* The initial sparse front contains 14,645,170 nonzeros and has rank 39,152.
* Staged solving plus complete factored verification: 25.41 s, 0.506 GiB,
  complete invariant dimension three.
* Adding propagation of finite-field two-variable relations: 19.686 s,
  0.502 GiB.
  This test ran at reduced CPU priority while a five-loop experiment ran;
  its timing is not a controlled speed ratio.
* Sampling with the complete residual check: 17.068 s and 0.506 GiB.
  A fresh run of the final auxiliary-prime/checkpoint version took 14.662 s
  and 0.511 GiB. Both reproduce the exact accepted three-dimensional space.
  The latter's physical coefficient solve took 0.307 s and 67.2 MiB; full
  collinear verification took 2.417 s and 179.9 MiB. Its complete E4 is
  byte-identical to the accepted result, and all four independent published
  word coefficients pass again.
* Both final spaces agree exactly with the accepted 7+1 invariant space
  after expanding the single changed sewing interface.
* The new collinear coefficient stage finds a unique four-loop candidate.
  Full verification checks all 533,360 divergent words in the collinear
  expression. The complete 985,265-word collinear symbol is byte-identical
  to the accepted four-loop E4, including its finite part.
* Coefficient discovery took 0.404 s and 79.1 MiB sampled peak RSS. Full
  collinear verification took 2.519 s and 214.1 MiB, including its child
  projection command.
* An independent original-alphabet word checker reproduces the four
  published four-loop coefficients 960, 960, -960 and 120. The formulas are
  [equation (12) and the following paragraph of He et al.](https://arxiv.org/html/2511.09669v2).
  They are validation data, never inputs to the coefficient solve.

Evidence: `output_earlysew_20261003/`, especially `sew6p2.*`,
`early6-final.*`, `early6-pairs.*`, `fourloop-pairs-compare.log`, and
`amplitude4-{coefficients,verify}.*`.

## Recursive actions and five-loop inputs

The previously certified FEC8 has dimension 17,471. Its cyclic, flip and
parity matrices were built recursively from private carrier coordinates in
70.18 s, with 12.13 GiB sampled peak RSS. The action nonzero counts are
42,636,958, 52,037,405 and 32,937,307 respectively. All original local
conditions and both boundaries are checked for closure under the supplied
generators before sewing.

The weight-ten 8+2 coefficient space has 2,061,578 variables. Exact local and
symmetry presolve leaves 478,037. The first sparse front has 1,083,241 rows
and 18,270,204 nonzeros, giving modular rank 213,487. No full FEC9 tensor or
weight-ten word expansion is used in this stage.

Two earlier five-loop variants were deliberately stopped to apply validated
improvements. The first ran 1,329.8 s and peaked at 14.24 GiB before stopping;
it produced no certified weight-ten kernel. The zero-only variant was also
stopped. Their JSON reports record nonzero exit status and are not successes.
The current variant adds doubleton substitution. Its result is pending.

A separate discovery-mask experiment avoids regenerating original equations
already consumed in earlier bands of the same prime. Its four-loop control
passes the complete exact space comparison, collinear check, byte-identical
E4 comparison and published word checks. It reduces discovery rows from
1,524,437 to 1,288,533, but its contended time of 19.98 s does not establish
a speed improvement over 19.69 s, and its returned basis is denser. It is
therefore optional (`refine-cache`), rather than the default `refine` route
used by the current five-loop trial. Both routes perform the same final
complete factored check.

## Boundary memory bottleneck

The independent five-loop boundary calculation initially failed with SIGSEGV
under a 20 GiB address-space limit, after 174.92 s and 14.33 GiB sampled peak
RSS. It had not completed its first shuffle term. The existing shuffle
implementation allocates a vector for each word in a large hash table and
then materializes an additional tensor. A failed allocation at the imposed
limit is consistent with this failure; no stack trace was captured.

`packed_shuffle.hpp` replaces those vector keys with checked, injective
64-bit word keys and accumulates the entire weighted boundary in one table.
It caches bounded blocks of partial shuffle keys, preserves all interleaving
multiplicities, and rejects key overflow. Output extraction releases hash
nodes while building the tensor. This is a general exact rational word
operation, not a table of heptagon coefficients. It is currently used by the
experimental validation driver, leaving the public shuffle route unchanged.

Sixty-four independent rational weighted-product tests compare it with the
existing word-vector implementation, including cancellations, repeated
letters, empty factors and key overflow. The complete four-loop boundary is
byte-identical (1.507 s, 116.4 MiB sampled peak RSS). The **complete five-loop
boundary successfully finished** in 492.64 s with 14.07 GiB sampled peak RSS,
under the same 20 GiB address-space limit. Its saved file is
`output_earlysew_20261003/amplitude5/boundary_5L.wxf`: 112,344,203 nonzeros and
1,265,737,523 file bytes. This boundary comes from
the accepted lower-loop E/R symbols and does not depend on a five-loop
candidate. Having the boundary is not yet a five-loop amplitude result.

The independent collinear projection through FEC7 also completed: 164.85 s
and 5.64 GiB sampled peak RSS, including the lower projection chain. The
weight-seven image has dimension 3,279, with 8,429,985 coefficients in its
coordinate map and 11,393,050 in its recursive basis. It is saved under
`amplitude5/projection/collinear/` for the final candidate check.

The projected LEC2 has rank 49 in its 121 letter-pair coordinates. The
experimental driver now keeps a rank-49 right coordinate frame through the
large forward contraction and expands the two right letters afterward.
`restricted_projection::sew_one_compact_right` implements this general
change of contraction order. Its coordinate factorization is checked exactly.
Independent five-factor tests cover zero, deficient and full right rank.
The complete four-loop physical check again passes and produces byte-identical
E4. That control took 2.711 s and 205.1 MiB, versus 2.519 s and 214.1 MiB
previously; these contended measurements establish no timing improvement.
The structural saving is the reduced number of right coordinates in the
expensive forward multiplication. The ordinary public projection path is
unchanged.

The kernel run initially had a 90-minute wall allowance. A separate monitor
permits at most three hours without changing its 42 GiB address-space limit.
If that extension is used, the waiting original measurement process is
paused while the solver continues; the replacement monitor records process
RSS/high-water RSS and enforces the new wall deadline. The original waiting
parent resumes to collect the real child exit status. Use
`early8-pairs-extended.json` for the complete peak and time; the original
report's RSS sampling excludes any paused interval. Successful completion
past the old deadline and termination at the replacement deadline were both
checked on separate owned test processes before using this monitor.

The original pair-propagation solver was later paused while the sampled and
hybrid alternatives were tested, then deliberately stopped to release RAM.
Its last reported rank was 474,680 of 478,037; no rational kernel was saved.
Peak RSS was 18.71 GiB. The recorded 9,606.89 seconds includes the paused
interval and must not be presented as uninterrupted solve time. The resumed
original sampler re-evaluated its superseded 90-minute deadline and marked
`timeout=true`, although the solver was deliberately terminated before the
replacement three-hour limit. The raw report is preserved; the explicit
interpretation is in `early8-pairs-extended-audited.json` and the termination
record is in `early8-pairs-stop-note.json`.

## Collinear validation

`bench/early_mhv_bootstrap.cpp` uses small alphabet projections to identify
candidate coefficients. These are necessary conditions from the lower-loop
shuffle recursion; their full-rank system establishes at most one candidate.
It then projects the physical candidate to the full collinear alphabet and
checks **every word containing any divergent letter** against the original
shuffle boundary. Only after that check is the final recursive amplitude
written. Discovery coefficients alone do not certify an amplitude.

A proposed shortcut based on commuting divergent letters failed an exact
row-space membership prerequisite on these data. It was rejected and is not
used by the verification path. No conclusion about the physical amplitude
was drawn from that failed prerequisite.

## Tests and status

Sixteen exact small product-space tests compare the new kernel with an
independently assembled rational matrix, including nontrivial S3 irrep
pairings, rational changes of basis, empty kernels and staged restriction.
Factored residuals are compared entry by entry with the explicit operator.
The rebuilt structured-kernel suite passes 160 independent exact cases,
including unlucky primes, prime denominators, vanishing coefficient supports,
multi-prime reconstruction, and the experimental batched backend. The rebuilt
factored-kernel suite passes 102 rational extension/sewing cases against
independently assembled equations. A fresh four-loop boundary calculation
also matches the accepted boundary byte for byte.
The existing default solver paths remain available; this experimental
MHV sewing route has not been made the default.

Hardware performance counters were unavailable (`perf_event_paranoid=4`);
no system permission or other user's process was changed.

## What the unsuccessful trial establishes

Materializing unrestricted FEC9 remains excluded by its measured 82.55 GiB
equation payload alone. Early sewing avoids that allocation, and the full
five-loop collinear boundary now completes within ordinary workstation RAM.
However, these results do not establish that the complete bootstrap fits
the same peak, or that its solve time is moderate. No five-loop amplitude
or published five-loop word-coefficient match is claimed.

The measured obstacle is the balance between dependent-equation work and
elimination fill. Sampling cheap equations too early lost useful sparse
pivots and increased retained storage almost ninefold at a comparable rank.
Completing cheap bands preserves sparsity, but later spends minutes on
blocks that add few or no independent constraints. Caching row identities
removes exact replays but cannot remove general linear dependencies.

The next algorithmic experiment should keep the cheap pivot stage, then
adapt the number of selected equations to the remaining dimension and
spread selections across all symmetry generators and local-condition
blocks. Its complete factored residual check must remain mandatory. Saving
the triangular finite-field state at band boundaries would also avoid
repeating the expensive short-equation stage when changing later discovery
policies. These are identified next steps, **not implemented improvements
or promised five-loop timings**. Increasing the initial batch alone was
tested on the control and did not improve its time or memory.

## Reproduction

Build `bench/early_sew_probe`, `bench/early_sew_check`,
`bench/early_mhv_bootstrap`, and `bench/early_mhv_words` with the repository
Makefile. `make check-early-sew` runs the small exact product-space checks.
The numbered experimental binaries preserve the precise code used by each
timing run; use the final canonical binary for a fresh reproduction.

The trial's `chain/` directory points to the previously accepted FEC2 through
FEC7 and to the newly certified FEC8. `actions/` holds recursively computed
generator matrices. Each sewing output directory has an `actions` symlink
to that shared action directory. With fresh output paths:

```bash
./bench/early_sew_probe actions data FEC_CHAIN 8 ACTIONS
mkdir -p EARLY
ln -s /absolute/path/to/ACTIONS EARLY/actions
./bench/early_sew_probe refine-hybrid data FEC_CHAIN 8 EARLY
./bench/early_mhv_bootstrap boundary data FEC_CHAIN EARLY LOWER_OUTPUT AMPLITUDE 5
./bench/early_mhv_bootstrap coefficients data FEC_CHAIN EARLY LOWER_OUTPUT AMPLITUDE 5
./bench/early_mhv_bootstrap verify data FEC_CHAIN EARLY LOWER_OUTPUT AMPLITUDE 5
./bench/early_mhv_words data FEC_CHAIN EARLY/LEC_2.wxf AMPLITUDE/hepMHV_5L_recursive.wxf 5
```

`LOWER_OUTPUT` supplies the accepted `2loop/`, `3loop/`, and `4loop/` E/R
collinear tensors. The independent boundary can be computed before sewing
finishes. Only successful full verification writes the final recursive
amplitude. Discovery files named `candidate_*` are not certified outputs.
The new recursive tensor uses the **8+2** cut and must be interpreted with
the saved FEC1..8 and LEC1..2 bases. It is not a 9+1 tensor and cannot be
substituted into that interface without a change of coordinates.
