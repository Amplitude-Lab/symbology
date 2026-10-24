# Validation record and changelog

Validation scope and measured optimization history, newest entries last. The
[README](README.md) carries the committed result and reproduction instructions
only; solver-optimization campaigns and their interim statuses live here.

The native and shifted results use the same fixed recurrence bases and CG
normalization. All statements below concern ordinary symbols. Conditional
uniqueness assumes the extrapolated double-final-entry bonus relations.

## Exact checks already completed

| Check | Coverage | Wall seconds | Peak RSS, GB |
| --- | --- | ---: | ---: |
| Native → shifted rational conversion | Exact factorization, all 518,048 shifted coordinates | 170.324 | 14.285 |
| Native original seam | 2,053,259 equations; 700-bit bound, 721-bit modulus | 116.652 | 14.334 |
| Native origin | 44 independently predicted coefficients | 145.819 | 18.896 |
| Shifted original seam | 4,196,956 equations; 514-bit bound, 541-bit modulus | 67.667 | 4.398 |
| Shifted origin | The same 44 predictions in the other representation | 44.998 | 5.624 |

The native seam certificate checks all eight saved native congruences and
four additional residual-only primes. The latter do not rerun a native solve.
The shifted certificate uses eight derived congruences and one additional
residual-only prime. In each case the denominator-cleared integer residual
bound is less than the product modulus, proving exact vanishing.

The rational exporter has independent eight-/nine-loop controls and WXF
readback. Wolfram's decoder additionally compared all native coefficients
with the reconstructed TSV and all converted coefficients with the accepted
shifted TSV. These tests detect output/coordinate errors independently of the
exporter's own sparse writer.

## New physical checker controls

Before ten-loop checks, the new checker passed:

| Saved reference | Collinear rows | MRK rows | Extended OPE coefficients |
| --- | ---: | ---: | ---: |
| Eight loops | 8 | 128 | 2,100 |
| Nine loops | 8 | 144 | 2,618 |

A deliberately doubled eight-loop symbol failed all 64 MRK equations in a
separate negative control. Successful controls and this rejection are both
required. The new checker reuses the existing projection engines; new points
and unused OPE data do not constitute an independent implementation of every
physical formula.

The extended ten-loop OPE table contains 3,192 coefficients. Its overlap of
2,280 coefficients with the original table agrees exactly. The remaining
912 coefficients are at T1 S11/S13 and T2 S12/S14, with all supplied log
powers. They were absent from the original fit and its holdouts. The input
plan pins two new primes, seed 202610090731 and sample offset 1000003 before
the ten-loop checks. Physical checks remain finite-field tests; they are not
a rational identity certificate for all collinear, MRK or OPE data. The OPE
scope is maximal helicity at twists one and two; helicity-zero T2 and higher
twists are not covered. The original fit included solve-role T1 rows as well
as T2 rows, so all T1 data must not be mislabeled as unused holdouts.

Both new primes passed all requested rows: **9,072 finite-field checks** in
total, including the same **912 distinct unused OPE coefficients at each
prime**. There were no mismatches. These are validation equations, not 9,072
new independent constraints used in fitting.

| New prime | Collinear seconds / peak GB | MRK seconds / peak GB | Extended OPE seconds / peak GB |
| --- | ---: | ---: | ---: |
| 1 | 88.187 / 16.962 | 275.950 / 16.912 | 1157.125 / 18.908 |
| 2 | 93.250 / 16.438 | 288.675 / 16.905 | 1083.189 / 18.995 |

All six ten-loop validation jobs took 2986.376 s in sum, with an
18.995 GB sampled peak. Known-reference controls and exact certificates
are accounted for separately. See `evidence/physical-review.json`.

## Adaptive rational reconstruction

The uniform 48-bit numerator bound failed on 62 native coefficients when
using six CRT primes: 31 had no candidate and 31 failed the two saved modular
checks. The verified native denominators reach 315 bits; the fixed trial's
denominator bound was only 312 bits.

A retrospective minimum-product Euclidean-convergent experiment recovers
all 275,246 native coefficients using six primes. Both saved check caches
agree, and all 265,143 coefficients across the eight-/nine-loop exports are
also reproduced. Maximum individual 2|n|d needs 329 bits, while the global
box formed from separate numerator and denominator maxima needs 361 bits.
The measured repetition, including controls, took **9.754 s and 0.183 GB**.

This is a useful candidate-selection heuristic, not a proof of general
reconstruction uniqueness. Freeze the method before future validation primes.
The current native export retains its explicit seven-prime reconstruction
provenance; the experiment does not retroactively create fresh check primes.

## Further solver optimization

Early physical equations can reduce dimension while increasing elimination
fill. In the recorded V17 pilot comparison, discovery dimension fell from
22,489 to 22,032 while nonzeros rose from 384,425,474 to 400,158,698. Later
collinear matrices became much sparser, offsetting some of that cost. Overall
wall time improved by only 0.66%, at different primes and background loads;
this is not a controlled speedup result.

The V17 prime-2 intervals between checkpoints were 3,982.214 s for collinear,
124.299 s for MRK, 934.674 s for OPE and 1,251.481 s for final sewing. They
include checkpoint work. Discovery and collinear fitting remain the main
whole-run targets. New higher-S rows should be tested for rank gain, resulting
nonzeros and total remaining solve time before being admitted earlier. A
useful split is to test S11/S12 for fitting while retaining S13/S14 as checks.
This scheduling experiment has not yet been measured.

The completed same-checkpoint comparison reproduces all 275,246 native
coefficients with both backends. No equation or finishing algorithm changed.

| Backend | Final-seam seconds | Whole process seconds | CPU seconds | Peak process RSS, GB |
| --- | ---: | ---: | ---: | ---: |
| Mapped | 1276.920 | 1349.263 | 2781.930 | 13.959 |
| Heap | 222.071 | 289.422 | 2268.137 | 14.410 |

The finishing phase is **5.750× faster**, with 82.61% less phase wall time.
Including input loading, the process is 4.662× faster. The sampled peak grows
by 0.451 GB. All runs stay within the 50 GB cap; live cgroup snapshots show
zero high/max/OOM events and zero memory pressure at the observed times.
The transient services were collected before the post-exit `systemctl show`
queries, so the default infinity limits in those files are not measurements
of the running jobs. Preserved live snapshots verify their actual caps.

The mapped log records at least 8,310,784 free-range reclamation calls. The
allocator performs hole punching while holding its mutex. System CPU time
falls from 877.627 s with mapped storage to 82.527 s with heap storage, while
whole-process CPU utilization rises from 2.06 to 7.84 cores on average. These
observations support allocation and reclamation overhead as a substantial
cost in this phase. They do not identify every remaining bottleneck.

Use heap storage for a separately launched finishing replay of this tested
checkpoint. A future full solver should benchmark a controlled switch after
large spaces shrink, preserving ownership of live buffers or restarting from
a verified checkpoint. This does not justify replacing mapped storage in
large discovery phases. Both measurements used the same O2 binary, prime,
checkpoint and 16-thread setting, ran serially in mapped-then-heap order,
and used warm OS/input caches on a shared workstation. There is one trial
per backend; this is a phase measurement, not a whole-prime speedup or a
global optimum. See `evidence/replay-comparison.json`.


The remaining algorithmic comparisons, in order of the observed expensive
phases, are:

1. **Reuse a verified elimination order across primes.** Discovery and the
   collinear fit dominate. Reusing which original rows and pivots were useful
   may avoid repeated discovery, while recomputing their values at each new
   prime. Every pivot must be checked for invertibility, with fallback on an
   exceptional prime. Reusing solution coefficients as fitting equations would
   be circular and is not proposed. The V18 follow-up below measures final-seam replay at eight and nine loops;
   V26 below now tests discovery replay at eight and nine loops. Its ten-loop
   cost and fully implicit cross-stage physical projection remain unmeasured.
2. **Benchmark earlier higher-S OPE rows by total cost.** The new unused data
   provide a concrete test set. Rank gain alone is inadequate: record output
   nonzeros, projection cost, subsequent seam time and peak memory. Keep
   separate physical data reserved for validation after changing the fit.
3. **Use weight-dependent log truncation in OPE heads.** The existing storage
   uses the total weight-20 triangular log bound. A weight-j head can potentially
   use a+b<=j before later letters raise the weight. At j=17 the triangular
   slot counts are 171 versus 231. This is a possible storage reduction, not
   a ten-loop speedup. The V22 follow-up implements this layout and validates it
   against square operators, complete FEC8 heads and full nine-loop results.
   Ten-loop performance remains unmeasured.
4. **Separate large discovery and small finishing allocations.** The completed
   replay demonstrates a benefit for the tested finishing checkpoint. Any future implementation
   must preserve ownership of live mapped buffers when changing allocation
   policy. Low instantaneous RSS does not authorize exceeding the 50 GB cap.

These are concrete remaining comparisons. The solver is feasible and has
validated optimizations, but calling it fully optimized would be unsupported.


A separate package comparison is archived in
`output_hexagon_algorithm_comparison_20261008/COMPARISON.txt`. Its L4/L6
rational coefficients and all 38,989 L8 modular coefficients agree with the
comparison references. That work led to the separately validated V18 structural replay
implementation summarized below. The allocator measurements here retain frozen
V17 mathematics; V18 claims use their own matched controls.


A further finishing-stage candidate is recorded in
`output_hexagon_validation_20261009/rank-saturation-observation.json`.
The current restriction trigger waits for 64 pivots or the end of the sweep.
When fewer than 64 parameters remain, a `rank == remaining_dimension` trigger
could restrict immediately once a unique representative is determined. The
subsequent complete zero-residual sweep must still run. The current logs do
not record where rank saturation occurs, so no time saving is established.
This is an unselected proposal for controlled testing, not a change to V17.

## V18 follow-up: measured eight- and nine-loop controls

A separate source package now implements direct contraction of seam factors
into candidate coordinates, a forward learning scan that avoids revisiting
satisfied prefixes, and reuse of selected source-row addresses across primes.
All modular coefficients and right-hand sides are recomputed. The selected
rows determine one affine lift, followed by a complete seam check. Discovery,
early physical scheduling, CG conventions and the bonus ansatz remain those
of V17. In V18, structural replay is limited to final sewing.

The saved 19-stage validation and 208-file build receipt were independently
reviewed here. Every first-prime native coefficient was compared again with
the original saved payload. The logs confirm all current physical holdouts
and both primes' shifted comparisons against the saved reference values.
Correction: the original comparator could select existing finite-field caches,
so those calls did not establish that the rational files were reduced anew.
The V22 follow-up below verifies forced rational-reference reduction at nine
loops using a snapshot without caches. Phase controls also retain original-implementation
checks of all 133,214/372,388 factored sewing rows. These finite-field tests
are distinct from the ten-loop exact original-seam certificates.

| Complete native first prime | V17 seconds / peak GB | V18 seconds / peak GB | Time reduction |
| --- | ---: | ---: | ---: |
| Eight loops, 2 GiB OPE cache each | 60.861 / 1.871 | 58.907 / 1.874 | 3.21% |
| Nine loops, 2 GiB OPE cache each | 511.502 / 5.957 | 495.895 / 5.950 | 3.05% |
| Nine loops, V18 enlarged to 3 GiB | 511.502 / 5.957 | 467.490 / 5.984 | 8.60% |

The last row changes both algorithm and cache configuration. At nine loops,
each terminal OPE head occupies 876,573,360 bytes; three exceed 2 GiB but fit
within 3 GiB. The larger budget reduces terminal generations from four to
three, saving about 29.9 seconds of logged head-generation work. It does not
establish a universal memory guarantee or a peak-memory reduction.

On the identical nine-loop finishing checkpoint, V17 takes 22.075 seconds,
V18 learning 6.989 seconds, and V18 replay 3.391 seconds. The actual later
prime reuses 380 selected rows and performs one lift: its small solve takes
0.040 seconds and its complete final verification 3.337 seconds. At eight
loops the checkpoint already has zero directions, so its finishing gain
comes from verification rather than row learning.

These controls use matched O3 builds, native heap, 16 threads and warm bases,
with single complete-prime measurements. They ran after the ten-loop
allocator comparison, under a 20 GB hard limit with no swap and zero recorded
limit/OOM events. Their compiler setting differs from the O2 ten-loop replay;
the speedup factors must not be multiplied or treated as a measured ten-loop
V18 gain. A ten-loop pilot and timing replication remain outstanding.
The V26 follow-up below subsequently tests discovery replay at lower loops.

Evidence is archived in `evidence/v18-lower-loop`, including the independent
receipt review. The implementation remains in
`output_hexagon_v18_20261009/package` at the repository root. The accepted
17+3/16+4 ten-loop references and their producer identities are unchanged.

## V22 follow-up: complete nine-loop CPU and memory measurements

The historical V19–V22 study selected **V22 with 16 threads for the measured
nine-loop workload**. Its then-applied gate required more than 15% fewer total
user-plus-system CPU-hours than original V17, with at most 5% higher sampled
peak RSS. The memory comparison uses the lower of the two V17 peaks. We
recomputed the means and gate, verified 14 run receipts and 834 source/build
file bindings, and repeated every native payload comparison.
These historical V17 comparisons do not establish the current V16-baseline
qualification described in the selection correction below.

| Configuration | First-prime trials | Mean CPU-hours | Mean wall seconds | All-trial peak RSS range, GB |
| --- | ---: | ---: | ---: | ---: |
| V17 / 16 threads | 2 | 1.481645 | 508.864 | 5.903–6.011 |
| V20 / 16 threads | 2 | 1.049393 | 384.887 | 5.877–5.935 |
| V21 / 16 threads | 3 | 1.039490 | 380.377 | 5.849–5.939 |
| V21 / 8 threads | 1 | 1.178127 | 684.994 | 5.323 |
| V22 / 16 threads | 2 | 0.918766 | 347.391 | 5.840–5.890 |

V22 saves **37.99% CPU time versus V17** and **11.61% versus V21/16**.
Second-prime runs test correctness and replay; they are not repetitions of
first-prime timing. RSS ranges include both primes and all repetitions.
The eight-thread alternative loses to V21/16 on CPU and wall time. Its
wall-time-driven early-physics scheduling also changes the discovered space,
so this is a full-configuration comparison, not isolated hardware scaling.

The main algorithmic change compresses the affine collinear suffix operator
*with all sampled RHS columns*. Selected original rows come with explicit
row-span certificates. The 36 suffix rows first reduce to 14; exact zero-head
support, recalculated separately for every sample batch and field, reduces
them to nine. RHS-only equations remain. This removes redundant large sparse
products while preserving the complete constraint system and final checks.

V22 additionally stores logarithm slots only for a+b<=w at intermediate head
weight w, retaining the full T/S jet and restoring the existing external
layout before projection/cache use. At weight 15, 136 slots replace 190.
Controls cover 2,160 square-reference operator cases and all 1,178,100 FEC8
head coefficients at each of two primes. Mean nine-loop OPE CPU time falls
from 1,313.38 s in V21 to 1,000.55 s in V22. Collinear time also improves;
the shared recurrence's separate input/output strides changed, and no
experiment isolates that compiler/layout contribution. The entire speedup
must not be attributed solely to OPE storage.

All 105,898 native coefficients and current COL/MRK/OPE and complete
372,388-row factored seam checks pass. V20 shifted outputs match all 195,155
coordinates from actual rational-reference reductions at both primes.
Those reference snapshots contain no cache directory and match the original
rational files by hash. V22 native payload equality with V20 transfers this
check. The earlier V18 wording about direct QQ reduction was too strong;
its saved evidence is retained, with this correction. The ten-loop independent
Wolfram comparison and exact sewing certificates are unaffected.

The trials use prepared bases, heap allocation, 16 threads except the stated
control, 4 GiB physical cache and 3 GiB OPE cache (original V17 uses 2 GiB).
They were serialized under a 20 GB hard cgroup cap, no swap, and zero recorded
limit/OOM events. There is no confidence interval, cold-basis timing or V22
ten-loop measurement. The production limit remains 50 GB. Discovery still
uses about 1,033 CPU seconds per nine-loop run and is the next useful target
for guarded row/pivot replay or precomposition of the double-bonus tail basis.

The separate mathematical audit adds a physical cross-check: the complete
21-term [Eq. (92) and its cyclic images](https://arxiv.org/html/2308.08199)
pass on the completed ten-loop native symbol at two primes, after fitting.
They add rank three beyond the local double-bonus source (65 to 62 local
triple directions), but that source does not imply them. Exact local witnesses
are not physical-amplitude counterexamples. These rows were not used in V22;
there is no new ten-loop rational triple certificate or all-loop triple proof.
The earlier truncated transcription was invalid and has been superseded by
the complete v3 checks. An unrelated quadruple calibration, 169 versus 166,
remains unresolved and is not used to constrain or judge the amplitude.

The reviewed receipts are in `evidence/v22-nine-loop`; the pinned nine-loop
profile and reproduction command are in its `accepted-nine-loop-profile.json`
and `REPORT.txt`. No new ten-loop run or producer change follows from this
review. The native 17+3 and equivalent 16+4 references remain unchanged.

## V26 follow-up: generic discovery replay and triangular affine factors

The new study extends structural replay into discovery and stores triangular
affine equations before expanding a full coordinate map. We independently
reviewed all **21 native run receipts**, the selected build's **221 files**,
and every recorded native coefficient comparison. CPU means and RSS gates
were recomputed, with 1,270 input hash bindings recorded in
`evidence/v26-generic/v26-generic-review.json`. The nine-loop reference chain
retains the earlier forced rational-reference reduction audit.

Discovery plans store source positions, pivot groups and bounded dependency
routes. Each field regenerates every matrix coefficient and RHS and checks
all discarded equations against its new basis. Missing routes use exact local
elimination; changed source identity, row streams, pivots or rank cause a
fresh discovery solve. Prior solution values never become constraints.
The affine factor stores equations x_j = b_j - sum_{k>j} U_jk x_k, composing
free-coordinate maps with a shared constant coordinate. This avoids an
intermediate expanded reduced matrix. The current driver still materializes
complete affine maps at its existing checkpoint boundaries.

| Complete native solve | Learning CPU-h | Replay CPU-h | Mean CPU-h | Mean wall seconds | Conservative peak GB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Eight-loop V26 / 16 workers | 0.097601 | 0.098072 | 0.097837 | 47.654 | 1.883 |
| Nine-loop fresh V22 / 16 workers | 0.915311 | 0.907920 | 0.911615 | 342.590 | 5.949 |
| Nine-loop V26 / 16 workers | 0.880155 | 0.864388 | 0.872271 | 324.935 | 5.943 |

Each mean covers **one learning prime and one replay prime**, not repeated
timings of the same prime. V26 saves 41.13% CPU time against the historical
nine-loop V17 mean (1.481645 CPU-h, two first-prime repetitions), and **4.32%**
against the fresh matched V22 pair. The eight-loop saving versus its original
V17 control is 23.34%. These comparisons include learning costs. Original
V17 uses a 2 GiB OPE cache; the optimized profiles use 3 GiB, so the historical
comparison measures the whole configuration, not only the new kernels.

The new RSS metric is the larger of 100 ms sampled RSS and the process kernel
high-water marks captured in phase logs. At nine loops it is 5,942,951,936
bytes, 0.672% above the conservative historical V17 sampled reference of
5,903,253,504 bytes. At eight loops it is 1,883,369,472 bytes, 0.648% above
its historical reference. Earlier tables retain their explicitly sampled
metric. All selected runs meet the greater-than-15% CPU saving and
at-most-5% RSS growth gate then applied against original V17. These measured
V17 results cannot be relabeled as V16 qualification.

Profiling attributed 825 CPU seconds to incoming-row reduction versus 27 to
new pivot search, motivating dependency-route reuse. A 256 MiB route budget
was correct but failed the eight-loop memory gate: kernel peak RSS was
2.186 GB, or 2.060 GB after allocator reclamation. The 100 ms sampler missed
the latter spike, reporting only 1.884 GB. V26 therefore uses **64 MiB for
every input**, reclaiming unused allocator pages only after temporary owners
are destroyed at joined phase boundaries. V24/V25 have slightly lower
nine-loop CPU means but fail this cross-instance memory requirement; those
small time differences do not establish a performance ordering.

Generic controls pass 96 field cases, 48 structural replays and 72 affine-map
comparisons with independent dense elimination, including deficient and
inconsistent systems and two exceptional-prime rank cases. Deliberately
empty routes pass 48 fallback controls. Five other matrix families agree
with both the previous implementation and every original equation. Both
selected fields reproduce all 38,989 eight-loop and 105,898 nine-loop native
coordinates; complete seam sweeps and fresh COL/MRK/OPE checks also pass.
No new bonus relation was imposed.

The trials were serialized with 16 workers, native heap, 4 GiB physical cache,
3 GiB OPE cache, an 18 GB pressure threshold, a 20 GB hard cap and no swap.
Timings include loading, discovery, physical fitting, checkpoint and plan I/O,
and final checks; they exclude compilation, cold input generation and a new
rational reconstruction. The whole study consumed 10.688 CPU-hours across
21 fresh native runs. Controls establish correctness and these measured
workloads, not universal speedup, confidence intervals or global optimality.

This study motivated compressed dependency routes, selection by work saved
per byte, deterministic work-based physical scheduling, and carrying
implicit affine factors through physical projections. V27 below tests the
first three directions; complete implicit projection remains unmeasured. A changed
time-based schedule already falls back safely but can erase replay savings.
V26 ten-loop performance remains unmeasured; its lower-loop savings are not
applied to the recorded ten-loop cost. The existing 50 GB production limit,
accepted ten-loop coefficient files and producer identities are preserved.
The pinned profile, source receipt, failed alternatives and all run logs are
archived in `evidence/v26-generic`; full source packages and structural plans
remain in `output_hexagon_v23_20261009` at the repository root.

## V27 follow-up: measurements and corrected working selection

**V27D with 16 workers is the current working benchmark.** The qualification
threshold remains strictly more than 15% CPU saving against **V16**, with at
most 5% RSS growth. The earlier recommendation to retain V26 incorrectly
treated the requested V26 comparison as a new qualification threshold.
It was not a correctness rejection. Nine-loop qualification against V16
remains unverified: available historical qualification measurements used V17.
Selecting a working benchmark does not assert that this missing comparison
has passed.

Independent review checks all 27 complete native runs, 909 source/build file
bindings across four preserved packages, timing intervals, every coefficient
payload, and the recalculated measurements and ablations. Its 2,352 original
hash bindings remain historical evidence. The separate selection amendment
supersedes the earlier recommendation without changing any measurement or
physical result. The superseded execution gate is retained byte for byte;
its relocation and the amended current gate are explicitly verified.

| Complete native workload | Mean CPU-hours/prime | Mean wall seconds/prime | Peak RSS, GB |
| --- | ---: | ---: | ---: |
| Nine-loop fresh V26 / 16 workers | 0.874501 | 325.347 | 5.912 |
| Nine-loop V27D / 16 workers | 0.867992 | 297.189 | 6.010 |
| Nine-loop V27D / 8 workers | 0.838348 | 468.654 | 5.374 |
| Eight-loop fresh V26 / 16 workers | 0.099578 | 48.714 | 1.902 |
| Eight-loop V27D / 16 workers | 0.063614 | 26.606 | 1.784 |

Each row averages one learning and one replay prime, not statistical timing
repetitions. RSS uses the larger of 100 ms sampled RSS and kernel process
high-water marks. V27D/16 saves 0.74% CPU time and 8.65% wall time, and raises
peak RSS by 1.65% against fresh V26/16. The sub-percent CPU difference is not
an established statistical speedup. V27D/8 saves 4.13% CPU against that 16-worker baseline,
but takes longer. A separate matched eight-worker, first-prime comparison
gives V26 0.869639 versus V27D 0.845487 CPU-hours, a 2.78% saving; it is not
the two-prime comparison. The 36.12% eight-loop saving against V26 does not
establish nine-loop qualification against V16. A historical eight-loop V16
control gives 49.32% CPU saving and 1.44% RSS growth for V27D, but compares
one original mapped/cache-configured V16 prime with a new heap learning/replay
pair. It is neither a matched nine-loop experiment nor a substitute for one.

Four reusable implementation directions were tested. Workload kernels use
ordered bit queues, balanced row ranges, reusable worker scratch and exact
density-based sparse/dense dispatch. Delayed 128-bit accumulation corrects
overflow exactly; fields larger than 61 bits retain canonical reduction.
The scheduling policy compares physical rank gain per arithmetic work with
seam work, including fill, and records decisions for replay. A reusable
augmented row-span compressor certifies discarded equations with all RHS
and origin columns included. Finally, delta-varint dependency routes are
selected by estimated updates saved per byte, retaining a 64 MiB accounting
budget and an 8 MiB logical learning budget. Incomplete routes always finish
with exact residual elimination. These changes introduce no new physical
equations or bonus assumptions.

Local nine-loop first-prime ablations save 3.25% CPU for the kernel bundle,
1.43% for the compression extension and 0.64% for packed routes. The new
schedule costs 0.13% more CPU in that comparison. Over learning plus replay,
packed routes save 1.65% versus legacy routes. These interacting, single
comparisons are not additive gains or confidence intervals; disabling packed
routes restores the legacy cache, not the absence of structural replay.

The unsuccessful alternatives explain the limits. Uniform density sampling
aliased polynomial layouts, requiring exact density inspection. The initial
work-budget-only schedule admitted expensive MRK/OPE constraints too early;
V27B used 1.103377 CPU-hours/prime, substantially worse than V26. Comparing
physical gain/work against seam work recovered most of this regression.
The MRK feature-space operator stayed full rank, so that compression adapter
saved no rows. All completed variants remained numerically correct and their
separate binaries and measurements are preserved.

The physical-constraint requirement distinguishes **availability from
execution**. Register every required equation or lazy block generator from
the start alongside integrability and symmetry constraints. The scheduler
may assemble and eliminate costly blocks later, when the affine space is
smaller or their expected benefit improves. For x = o + Ny, a deferred block
Ax = b becomes ANy = b - Ao; changes of coordinates must update both its
linear map and its common origin. Pending blocks remain required, and
completion accounting and final checks must prevent their accidental omission
or reuse of stale projected rows. V27B's regression concerns execution order;
it does not justify withholding physical constraints from the initial problem.
The completed V27 tests do not establish a generic pending-block scheduler.

More CPU utilization is not automatically a CPU-hour improvement. In the
first-prime discovery phase, V26/16 used 930.056 CPU seconds in 126.644 wall
seconds; V27D/16 used 1,112.650 CPU seconds in 104.756 wall seconds. V27D/8
used 868.333 CPU seconds in 140.247 wall seconds. Eight workers also lower
MRK CPU cost but increase collinear and OPE CPU cost. **Phase-dependent worker
selection is an untested next comparison.** Hardware counters were unavailable;
these phase measurements do not identify cache-miss or bandwidth rates.
Profile records are nested and must not all be summed as disjoint costs.

Every run matches all 38,989 eight-loop or 105,898 nine-loop native reference
coefficients and passes complete sewing and fresh COL/MRK/OPE checks. The
prior forced-QQ-only reference chain is preserved. New and legacy kernel
switches pass 96 field cases, 48 structural replays, two exceptional-prime
fallbacks and 72 affine-map cases. Further controls include 144 augmented
compression cases, 64 deterministic schedules, 756,718 queue operations,
packed-route budget/damage tests, exact-density and persistent-cycle tests,
and 5,077 exercised overflow corrections against canonical arithmetic.

The 27 completed native trials consumed 16.098 CPU-hours, serialized with
warm prepared inputs, 4 GiB physical cache, 3 GiB OPE cache and heap storage.
Each used a 20 GB hard cgroup cap, 18 GB pressure threshold and no swap; all
recorded cgroup events were zero. Timings include loading, checkpoints,
plan I/O and final checks, excluding compilation, cold input generation and
new rational reconstruction. Full source and plans remain in
`output_hexagon_v27_20261009`; reviewed receipts are in
`evidence/v27-experimental`, with current selection in
`evidence/v27-selection-amendment/current-nine-loop-profile.json`. The original
`REPORT.txt`, `RESULTS.json`, completion and receipt review retain the old
decision as history; `selection-amendment.json` supersedes it. Legacy `passed`
fields in the copied comparison data refer to the superseded V26 rule, not
current V16 qualification. No ten-loop timing, producer or accepted result
changed at that selection stage. The subsequent first-prime measurement is
recorded below. The 50 GB limit remains, and global optimality is not established.

## First fresh V27D ten-loop prime: faster computation, higher peak memory

Prime 1152921504606847009 has completed a fresh native V27D/16 solve and
17+3 to 16+4 conversion. Independent receipt review rechecks **all 275,246
native and 518,048 shifted coefficients** against the preserved modular
references. The complete **989,800-row factored seam** and fresh supplied
collinear/MRK/OPE checks pass over this field. This is a new modular control,
distinct from the original unreduced-seam certificates over QQ for the
already certified result. The new run's rational reconstruction and exact
certificates are still pending. At this first-prime snapshot, the second
prime was running; its completed replay measurement follows below.

| Complete native learning prime | CPU-hours | Wall seconds | Sampled peak RSS, GB |
| --- | ---: | ---: | ---: |
| Historical V16 / 16 workers | 37.285539 | 12,711.891 | 34.628 |
| Fresh V27D / 16 workers | 23.379927 | 6,599.910 | 46.188 |

CPU falls **37.295%** and wall time **48.081%**, while sampled peak RSS grows
**33.383%**. Thus the greater-than-15% CPU target passes but the at-most-5%
RSS tolerance fails; the **combined ten-loop criterion fails** for this
configuration. The absolute 50 GB cap is a separate condition. The V27D
kernel high-water mark is slightly larger, 46.208 GB; the table compares
sampled RSS on both sides. The service-wide recorded memory peak was
47.774 GB, with zero high/max/OOM events at the comparison snapshot.

V16 used mapped storage and an 8 GiB OPE-head cache; V27D uses heap storage
and a 3 GiB OPE-head cache. Both use 16 workers and prepared basis inputs.
This historical configuration comparison is neither isolated algorithm
attribution, statistical repetition, nor the missing nine-loop V16
qualification. Replay and complete reconstruction costs cannot yet be
extrapolated from it. Conversion adds 65.740 wall seconds and 0.031945
CPU-hours, with 10.785 GB sampled peak RSS; it is excluded from the native
comparison above.

| Native phase | CPU-hours | Wall seconds |
| --- | ---: | ---: |
| Discovery | 13.395889 | 3,691.310 |
| Collinear restriction | 4.871917 | 1,266.440 |
| MRK restriction | 0.479636 | 125.562 |
| Low-S OPE restriction | 2.555983 | 809.382 |
| Complete final sewing | 0.092571 | 42.516 |

Discovery and collinear restriction account for approximately 57.3% and
20.8% of whole-native CPU cost. These top-level phase measurements do not
cover every loading/checkpoint/check cost; nested profile records must not
be added again. The final MRK restriction had zero rank gain in this run.
That observation alone permits neither omission of required equations nor
removal of independent MRK checks; skipping redundant elimination would
need reusable exact span/equation certificates.

The transient memory peak makes back substitution a concrete next target.
The frozen `streaming_echelon.hpp` calls `triangular_solver_2` after retaining
1,649,199,276 entries, or 19.790 GB of row payload. The transpose builder
in `SparseRREF/sparse_mat.h` appends incidence for every stored column using
mutex-protected `push_back`. At four bytes per row index, the full transpose
requires at least 6.597 GB of indices before capacity slack and metadata.
This is a storage lower bound, not the amount a proposed change would save.
The inspected recursive and scalar backsolve consumers access pivot columns.

A **pivot-column-only transpose with counted exact capacity** is therefore
a generic candidate: prove that all consumers need only these columns,
count incidence before allocation, and validate the fill phase and indexing
on finite-field, rational, affine, deficient-rank and exceptional-prime
controls. Actual pivot-column incidence and allocator retention must be
measured before predicting the memory benefit. A 20-second sample during
backsolve used approximately 15.06 active-core equivalents, so this evidence
does not indicate low task parallelism. Live transpose storage, later affine
materialization and retained allocations are not yet isolated.

No solver, controller, input, checkpoint or completed coefficient was
changed for this review, and no competing heavy benchmark was launched.
The new run continues under 48 GB MemoryHigh, 50 GB MemoryMax and no swap.
All original physical results and conditional bonus assumptions remain.
Archived evidence is in `evidence/v27-ten-prime1`; the active run is under
`output_hexagon_ten_v27d_20261009` at the repository root. The proposed
transpose change is unimplemented and unmeasured.

## Second fresh V27D prime: replay helps the seam, discovery still dominates

The replay solve at prime 1152921504606847067 passes the recorded complete
989,800-row factored seam and fresh supplied collinear/MRK/OPE checks.
Independent receipt review again compares all **275,246 native** and
**518,048 shifted** coefficients with the preserved same-prime reference.
At this snapshot a third V27D prime had been launched; no prime-3 receipt is
archived in this package. Fresh rational reconstruction, new-prime validation
and final exact certificates for this V27D chain remain pending.

| Native V27D/16 solve | CPU-hours | Wall seconds | Sampled peak RSS, GB |
| --- | ---: | ---: | ---: |
| Prime 1, learning | 23.379927 | 6,599.910 | 46.188 |
| Prime 2, replay | 23.329189 | 6,540.289 | 47.601 |

The observed whole-native CPU reduction is **0.2170%**, the wall reduction
0.9034%, and sampled peak RSS growth 3.059%. These are two different primes,
without repeated timings or a same-prime control disabling replay. They do
not isolate a causal replay speedup. Both native solves and their conversions
together consume **46.770291 CPU-hours**; prime-2 conversion adds 64.429 wall
seconds and 0.029231 CPU-hours. This is not the completed reconstruction cost.

| Phase | Learning CPU-hours | Replay CPU-hours | Learning wall seconds | Replay wall seconds |
| --- | ---: | ---: | ---: | ---: |
| Discovery | 13.395889 | 13.385500 | 3,691.310 | 3,655.980 |
| Complete final sewing | 0.092571 | 0.046964 | 42.516 | 21.163 |

Final-seam replay reduces contractions from 989,800 to 2,877 and affine
lifts from 34 to one, then retains the complete 989,800-row zero-residual
verification. Its CPU cost falls **49.27%**, but that difference represents
only **0.1951% of the first prime's whole-native CPU cost**. Verification
still takes about 20.6 wall seconds, while the replay solve takes 0.558 s.
This substantial phase improvement is too small to change the dominant cost.

Discovery CPU falls only **0.0776%**. The arithmetic counter records
12,329,356,644,947 coefficient updates when learning and 12,332,936,520,495
when replaying: both are about **12.33 trillion**. The retained routes execute
64,870,291 dependency steps alongside 3,231,942,664 ordinary queue steps,
about **1.968%** of their combined count. This fraction measures step coverage,
not saved CPU or arithmetic. The log reports 3,739 local residual completions
among 4,402 route uses. Source inspection confirms that routed prefixes are
followed by exact elimination of all remaining pivots. These completions are
not global replay failure; the discovery profile explicitly records replay.

At the comparison snapshot, the service had 4,481 MemoryHigh events, a
48.019 GB cgroup peak and zero MemoryMax/OOM events. Counts are cumulative
service events, not seconds spent stalled or a per-phase attribution.
The process sampled RSS maximum was 47.601 GB (kernel high-water mark
47.601 GB to this precision). The 48 GB soft threshold, 50 GB hard cap and
zero swap allowance remain unchanged. This adds evidence of memory pressure;
it does not repair the relative RSS qualification failure from prime 1.

For historical context, the same-prime original producer was **V17**:
38.571675 CPU-hours and 38.581 GB sampled peak RSS. V27D uses 39.52% less
CPU but 23.38% more sampled RSS in that configuration comparison. V17 cannot
be relabeled V16; the fixed V16 qualification policy and unverified direct
nine-loop comparison remain unchanged.

The measured priorities are discovery update arithmetic and locality,
backsolve transpose storage, and allocation lifetimes. Update counts alone
do not diagnose arithmetic versus bandwidth limits. Bounded route admission
should be judged by net queue/decode/residual work avoided, rather than
coverage or stored edges alone. The pivot-column-only exact-capacity
transpose remains a candidate requiring independent generic validation.
No larger route cache, producer change or competing heavy benchmark was
introduced. All required physical equations, affine RHS/origin and final
checks remain. Receipts and the independent review are in
`evidence/v27-ten-prime2`; the ongoing calculation remains under its frozen
V27D controller and existing limits.

## Input staging audit and gate (2026-10-23)

All reproduction inputs (FEC1–17 bases, LEC bases, physical data, OPE
checkpoint inputs, reference caches, the frozen V17 solver package, its
runtime libraries and the recorded seam-replay binary; 10 GB raw) were
packaged as split xz archives under `inputs/` — 48 parts, 4.1 GB, every
part pinned in `inputs/SHA256SUMS`. They are **kept local and not pushed**
(too large for the repository); the runners resolve them via
`bench/hexagon_bench_inputs.py`: `--source bench` forces the split copies,
`--source local` the original archives, `auto` (default) prefers the
originals. Historical absolute paths in `evidence/`, `controller-source/`
and receipt records are provenance of the original measurements and are not
read by the runners; no executable document depends on a local path.

Gate executed locally with `--source bench` (machine frequency-throttled,
concurrent single-core job): materialization of all components from the
split parts plus the heap seam replay completed with `REPLAY_PASS`,
`all_coefficients_equal=true` on all 275,246 columns (phase 243.3 s; the
same replay measured 222.1 s in the recorded serial pair and 323.6 s under
load with the original archives).

## Committed cold-regeneration path for the warm bases (2026-10-23)

The base-generation programs (`irr_fec`, `irr_lec`, `irr_sew`, `irr_convert`,
`irr_conditions`, `irr_pipeline`, `irr_cache`), the package `Makefile` and
`run_workflow.sh` stage entries, and the 60 MB dlogmat/cyclic-matrix seed
data (`solver-source/data/`, per-file checksums) are now committed, so the
FEC/LEC chains are regenerable from the repository alone. The mechanism is
not speculative: the audit receipt `audits/hexagon-ten-loop-2026-10-06/
fec16.json` records a full `irr_fec 16` rerun from the seed data (422 s
wall, 1,169 CPU-s, 13.7 GB peak, 8 threads). A complete cold FEC1–17 chain
rebuild has not been timed end-to-end. Of the package's build-time
`SOURCE_SHA256.tsv` pins, 406 of the 466 pinned files are present in the
packaged tree; 401 of those verify against their build-time hashes and
five carry post-freeze edits (`src/nmhv_ope_germ.hpp`,
`src/reduced_d3/irr_bulk.hpp`, `src/reduced_d3/irr_projection.hpp`,
`src/reduced_d3/irr_sew.cpp`, `src/SparseRREF/sparse_rref.h`), each
identical between the committed snapshot and the current local package
and matching the package manifest's committed-state hashes; the
remaining 60 pins reference solver data files not shipped in the
package. The OPE small-s
checkpoint and per-prime joint caches remain production outputs of the
solver phases; the committed `reference/` symbols are the certified
endpoints.
