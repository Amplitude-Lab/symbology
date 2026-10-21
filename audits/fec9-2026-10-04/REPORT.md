# Full FEC9 construction attempt — 2026-10-04

The requested target is the **complete weight-nine forward basis**, retaining
all symmetry sectors, as a prerequisite for comparing a weight-twelve `9+3`
cut with `8+4`. This is not a regrouping of the known five-loop amplitude or
an invariant subset of FEC9.

**Current state:** attempt5 completed and certified the full parity-adapted
FEC8 basis, but its FEC9 stage stopped on the host-memory floor after
18,696.76 seconds (5.19 hours), with 37.08 GiB sampled peak RSS. No FEC9
sector was certified. The complete FEC8 output and all logs are preserved.
The general simultaneous-constraint benchmark and the five-loop publication
are now complete. The latter passed a cold reproduction in 8,280.83 seconds
and 31.48 GiB and is published at `7f18e01` on `symrep-solver-dev`. After
bounded exact controls, attempt6 used compact integer factors and larger
pivot batches, reusing the certified parity FEC8. It stopped at 04:59:51 PDT
on 5 October on the host-memory floor after 16,336.40 seconds, with 33.67 GiB
sampled peak RSS. No FEC9 sector was certified. No solver is now active; the
original 09:22:13 PDT deadline and resource limits remain unchanged.

## Initial attempt — stopped by the host-memory guard

A frozen eight-worker staged-kernel run started on 4 October 2026 at
09:22:13 PDT, with a fixed deadline of 5 October at 09:22:13 PDT. Its output
directory is `output_fec9_20261004/`. It stopped at 10:19:35 PDT after
3,442.36 seconds, with 22.018 GiB sampled peak RSS. The 4 GiB host-available
memory floor triggered; neither the 50 GiB RSS ceiling nor the time limit
was reached. Its exit status is -15 (the guard's termination signal).
**No FEC9 basis or complete kernel certificate was produced.** The initial
`status.json`, `fec9-live.json`, `fec9.json` and `fec9.log` are preserved.
The `full-fec9-basis-run` thread heartbeat checks every 30 minutes, reports
meaningful changes, and stops following this attempt after completion or the
original deadline. The resource/deadline guards run independently between
those checks.

The controller accepts only a successful process exit, an exact complete
kernel certificate in the log, and the completed `FEC_9.wxf` file. It records
the output hash after success. The limits are 50 GiB per-process address
space, 50 GiB sampled process-tree RSS, and a 4 GiB host-available-memory
floor. There is one shared one-day deadline. At launch another user job held
roughly 23 GiB and the host had about 31 GiB available; the host guard can
therefore stop this attempt substantially below its nominal 50 GiB ceiling.
That unrelated job is left untouched.

`manifest.json` binds the condition tensor and FEC8 hash to the accepted
five-loop inputs, records source snapshots and binary hashes, and identifies
the algorithm. The fixed executables and controller are in `bin/`; source
snapshots are in `sources/`. No expected FEC9 dimension is supplied.

## Method and controls

`bench/streaming_kernel_benchmark` now exposes `forward-staged`,
`backward-staged` and `sew-staged` through the existing shared library. Like
the original benchmark path, it releases the loaded recurrence after creating
the factored operator, avoiding retention of an extra complete FEC8 tensor
through elimination and reconstruction. It uses the shared default discovery
policy, including bounded fronts, complete cheap bands, deterministic expensive
samples, and complete-pass fallback for large kernels.

Before launch, the frozen driver passed exact space comparisons against its
ordinary finite-field path for forward extension (911-dimensional output),
backward extension, and sewing. All seven process-budget/controller tests
passed. Commands, logs, elapsed times and exit codes are in `controls.json`.
No calculation-core or SparseRREF algorithm was changed for this initial run.

The last completed discovery band (16,384-term limit, 1/16 sampling) had rank
151,391 of 266,452, leaving 115,061 candidate directions. The latest logged
triangular count was already 1,005,540,275 nonzeros at rank 140,308; that is a
lower observation, not the peak retained count. Thus avoiding assembly alone
did not eliminate the fill bottleneck. At the follow-up inspection, the other
user job held roughly 30 GiB and host available memory was about 24 GiB after
our process had exited. The failed trial does not establish that the same
algorithm would exhaust 50 GiB on an otherwise free host.

Further trials retain the original deadline. The benchmark driver now exposes
the existing short-equation ordering and staged initial-support bound as
explicit optional arguments. Exact small-basis comparisons and bounded FEC7
controls are used to choose a retry policy before another expensive FEC9 run.
These are generic ordering controls, not polygon-specific constraints or an
invariant-space restriction. Results are retained in
`output_fec9_20261004/policy-controls/`.

The prior FEC9 count found 733,782 initial variables, 266,452 after exact short
relations, and 7,386,376,131 nonzeros in a fully assembled modular input. Its
82.55 GiB payload is a measurement of that representation, not a lower bound
on staged-solver memory. The present attempt tests the optimized representation
directly. Retained elimination fill and complete rational output size remain
to be measured; the small final kernel shortcut from the invariant five-loop
solve is not assumed applicable to the full basis.

For prior evidence, see the [FEC9 count and FEC8 construction audit](../fiveloop-2026-10-03/REPORT.md),
the [completed five-loop solve](../fiveloop-day-2026-10-04/REPORT.md), and the
[shared staged-library tests](../staged-kernel-2026-10-04/REPORT.md).

## First extended-run milestone

At about 09:58 PDT on 4 October (36 minutes after launch), the original
controller and solver were verified active with unchanged limits. The complete
1,024-term discovery band had finished at modular rank 110,584 of 266,452
columns, leaving 155,868 candidate directions. This is an intermediate upper
bound on kernel dimension, not the final FEC9 dimension. The run had moved to
the 4,096-term band with a 1/16 discovery sample.

The last logged triangular front retained 530,885,683 nonzeros. Sampled peak
RSS was about 17.6 GiB, and the host still had about 10 GiB available. The run
remained productive and was left unchanged. Its fill is already greater than
the final triangular fill of the invariant five-loop run, which underscores
why equal-order input sizes do not determine peak memory. No rational FEC9
basis or exact final certificate exists at this milestone. Timestamped state
and process observations are retained in `observations.jsonl`.

## Ordering controls and second attempt

After the first memory stop, six policies were tested in fresh processes on
full FEC5 and full FEC7 using the same input bases and eight workers. Every
FEC5 variant was compared by exact full row space with the ordinary result.
All FEC7 variants produced the complete 6,826-dimensional kernel with an exact
certificate. These are single policy pilots under concurrent user load, not
repeated timing medians.

| Full FEC7 policy | Wall s | Peak MiB | Output nonzeros |
| --- | ---: | ---: | ---: |
| Ordinary assembled finite-field core | 15.491 | 1,855.18 | 20,515,434 |
| Staged, natural order, initial bound 64 | 61.509 | 2,915.97 | 54,655,594 |
| Staged, support order, initial bound 64 | 23.937 | 1,583.36 | 30,687,758 |
| Staged, support order, initial bound 256 | 22.012 | 1,333.78 | 26,756,402 |
| Staged, support order, initial bound 1024 | 20.792 | 1,402.88 | 24,196,855 |
| Staged, natural order, initial bound 256 | 43.781 | 1,847.19 | 34,429,370 |

The 1024/support policy was selected for the retry because it reduced retained
triangular nonzeros from 21,667,987 to 6,765,818 and returned a smaller rational
basis. Its small peak-RSS disadvantage to the 256/support control was weighed
against the lower fill and output size. The ordinary assembled control remains
faster at FEC7; it is not used for FEC9 because its already counted input alone
exceeds the memory ceiling. No universal default change follows from these
pilots.

The initial global front still has the same 64-million-nonzero budget. If it
fills, remaining initial-band equations are streamed. All equations and kernel
directions remain required. Only driver policy arguments changed; the shared
kernel and arithmetic implementation are unchanged from the first frozen run.
`policy-controls/results.json` contains commands, bounds, process measurements,
output sizes and certificate markers, and all raw logs remain under the run
family. Attempt2 has a separate configuration, source/executable manifest,
launch identity and controller. Its deadline is inherited, not restarted.

At about 11:10 PDT, attempt2 had run for 32 minutes and reached a logged
modular rank of 153,803, beyond the first attempt's last recorded 151,391.
Its last logged triangular front had 950,310,093 nonzeros; peak RSS was about
20.9 GiB. It had completed the 4,096- and 16,384-term sampled bands and was
processing unrestricted row lengths at a 1/16 discovery fraction. The host
had about 20 GiB available at this check, so the free-memory guard was not
close to its threshold. Process identities, the unchanged configuration hash
and the original deadline were checked. The productive run was left running.
These remain partial modular discovery measurements, not a completed basis
or a prediction of final memory use.


## Second attempt — stopped by the host-memory guard

Attempt2 stopped at 11:39:51 PDT on 4 October after 3,736.57 seconds
(62.28 minutes), at 28.40 GiB sampled peak RSS. The host's 4 GiB available
memory floor triggered again; the 50 GiB process limit and original deadline
were not reached. The last logged rank was 175,544 of 266,452 and retained
triangular nonzeros had reached 1,421,619,420. Those are intermediate
observations, not a final dimension or peak fill count. Neither a FEC9 output
nor an exact complete-kernel certificate was produced. The controller and
solver have exited. `attempt2-result.json` preserves the result and log tail.

Support ordering improved the pilot and progressed farther on FEC9, but did
not remove the large-kernel fill problem. A further trial needs a substantive
change, rather than another identical run under the same shared-host limits.
The next candidate is a decomposition under commuting involutions, retaining
every character: smaller independent linear systems may reduce both fill and
peak RAM. This requires rebuilding and preserving the corresponding complete
recursive chain and validating it against ordinary full spaces.


## Bounded subgroup investigation

A generic replayable column-restriction adapter now transfers exact short
relations to the original coordinates, emits only the selected columns, and
retains the complete operator certificate. It passed 120 independent exact
rational-matrix tests, including reordered and empty selections. This adapter
does not by itself prove that separate kernels exhaust a full problem. The
experimental character-extension driver checks the seed/condition grading,
retains every character and certifies every sector independently. It preserves
the adapted chain, writes bounded-memory WXF with exact byte readback and
records separately hashed sector checkpoints. No default solver is changed.

For the heptagon, the flip/parity subgroup has four rational one-dimensional
characters. Exact ordinary-space comparisons passed through weight five with
the existing compact backend, and through weight six with the new streamed
column restriction. The old compact backend took 248.07 seconds and 5.54 GiB
for weights six and seven together, producing 24,329,265 FEC7 nonzeros. Its
weight-seven blocks alone assembled 361,615,877 rational entries.

The new factored-column driver produced all 6,826 FEC7 directions in 144.99
seconds at 2.566 GiB, with 25,971,935 output nonzeros, starting from its
validated weight-six chain. Every sector needed two reconstruction primes;
coefficient heights grew relative to the ordinary basis. These pilots are not
equivalent starting-coordinate benchmarks, nor evidence that full FEC9 fits.
A parity-only split is also being measured because it preserves a sparser
recurrence and smaller coefficients at weight six. All logs and resources are
in `output_fec9_20261004/subgroup-controls/`. The original deadline persists.


The parity-only input chain passed exact ordinary comparisons through weight
five, and the column-restricted full weight-six solve passed an independent
full-operator kernel comparison (2,555 directions). At weight seven:

| Policy / starting basis | Wall s | Peak GiB | FEC7 nonzeros |
| --- | ---: | ---: | ---: |
| Ordinary, original cached FEC6 | 15.491 | 1.812 | 20,515,434 |
| Four characters, adapted FEC6, assembled modular blocks | 144.988 | 2.566 | 25,971,935 |
| Parity, adapted FEC6, assembled modular blocks | 32.086 | 1.845 | 17,536,214 |
| Parity, same adapted FEC6, staged blocks | 55.686 | 1.068 | 20,579,497 |

All character runs retained every sector and had complete exact certificates.
The two parity policies share precisely the same input chain; the other rows
are different recursive basis choices, not isolated elimination comparisons.
These are single pilots under shared-host load. The assembled parity pilot needed one
reconstruction prime per sector; staged parity needed two for the even sector
and one for the odd sector. Merely maximizing the number of sectors was
not beneficial here. The smaller assembled-parity FEC7 is selected as the next
input, with staged elimination for the larger weights.

A raw equation-basis control at weight six also passed, but raised condition
nonzeros from 9,128 to 13,140 and did not improve the measured time or memory.
The prepared adapted equations remain selected. This is an experimental
benchmark option, not a default change.

The driver now releases each certified kernel after checkpointing it.
`partitioned_tensor_view.hpp` joins those row partitions through the existing
WXF writer while loading only one saved partition at a time; the factored
operator is released before joining. Twelve encoding/readback controls passed
against unpartitioned matrices, including empty axes, empty rows/partitions,
rationals and large integers. This prevents already completed sectors from
adding to a later sector's elimination peak. It does not cap the largest
single sector, its retained fill, or its rational kernel.


## Third attempt — all parity sectors, sequential checkpoints

The frozen pipeline launched at 12:09:30 PDT on 4 October, controller PID
1091906, in `output_fec9_20261004/attempt3/`. It retains the original
5 October 09:22:13 PDT deadline and the 50 GiB address-space/RSS ceilings
with a 4 GiB available-host-memory floor. It first reconstructs full FEC8
from the smaller validated parity FEC7, then solves all FEC9 sectors. The
controller follows the active stage rather than assuming a fixed FEC9 log.
Input bundles, their full recursive checksums, source snapshots, the executable
and controller are frozen and hashed in the attempt manifest.

Before launch, the final executable passed the exact weight-six full-space
comparison. Its w1 through w6 recurrence files are byte-identical to the
preceding validated driver; `final-subgroup-controls.json` records those
hashes. The partitioned writer releases each completed sector before the next
elimination and releases the factored operator before final concatenation.
The saved final recurrence is the ordinary rank-three WXF format, with the
prepared alphabet change and all prior recurrences retained.

Successful FEC8 must also reproduce its known full dimension 17,471. FEC9 has
no supplied target nullity. Acceptance requires every sector certificate,
completed exact serialization, layout metadata and output hashes. A sector
checkpoint alone is useful retained work but is not the whole FEC9 basis.
The parity decomposition partitions computation; it is not claimed to supply
the irreducible decomposition of the full cyclic/flip/parity group.


## Third attempt — reconstruction memory stop

At 12:35:46 PDT attempt3 was stopped by the host-available-memory guard,
after 1,576.07 seconds at 22.05 GiB sampled peak RSS. Neither the nominal
50 GiB ceiling nor the deadline was reached. The unrelated Wolfram process
had grown to roughly 31 GiB by the follow-up inspection and remains untouched.
The controller and solver have exited. FEC9 was never started.

FEC8's first parity sector completed modular elimination at rank 75,839 of
89,754 reduced columns, with 13,915 modular kernel directions. Its triangular
front retained 116,953,719 entries, while the reduced rows exposed to kernel
reconstruction had 261,827,028 entries. Three prime eliminations reached the
same rank/pattern; each cost about eight minutes, and none yielded an accepted
rational kernel before the memory stop. This exposes coefficient growth and
CRT storage as a further bottleneck, beyond input assembly and triangular
fill. No sector checkpoint or complete exact FEC8 certificate was produced.
`attempt3-result.json` preserves the status, measurement and log tail.


## Reduced-row pivot reuse investigation

Inspection found that a staged complete pass already back-substitutes its
independent rows, but the common reconstruction layer subsequently runs
SparseRREF again. This can change free coordinates despite the available
identity pivot block. A generic fast path now scans for a private unit column
for every row; these columns certify independence and a reduced pivot block
directly. If any row lacks such a column, ordinary elimination is unchanged.
The fast path is a mathematical certificate, not an assumed rank or a
polygon-specific pivot rule. It also logs how many rational rows have
reconstructed at each CRT modulus. Its effect on the failed FEC8 sector has
not yet been established. The generic 303-case staged/affine/composition suite
and 120 coordinate-restriction cases passed with the change.


The private-unit control passed: a constructed exact large-coefficient matrix
needed three reconstruction primes with its existing unit chart, versus five
after an extra SparseRREF pivot pass. The heptagon FEC7 control, however,
produced byte-identical w1 through w7 files and unchanged prime counts. It took
61.17 seconds / 1.47 GiB including the new partitioned-file join, compared
with 55.69 seconds / 1.07 GiB for the earlier driver retaining its matrices.
These single wall/RSS observations include different serialization paths;
no heptagon coefficient-growth improvement is established by this shortcut.
The next full-weight trial enlarges the bounded initial front to 256 million
entries so the complete 170-million-entry FEC8 short band is pivoted globally,
rather than splitting it at 64 million entries. Existing generic front-limit
tests passed; a final exact lower-weight control precedes launch.


## Fourth attempt — larger global short-row front

Attempt4 launched at 12:58:31 PDT on 4 October, controller PID 1093977, with
the same original deadline and resource guards. Its frozen driver passed a
complete weight-six comparison against an independently assembled full
operator with the new explicit front option (5.53 seconds / 0.306 GiB).
`pivot-controls.json` records all associated generic and application controls.
The configured initial front is 268,435,456 nonzeros, versus 67,108,864 in
attempt3. The observed FEC8 initial short band has 170,612,163 entries, so
it should fit as one global pivot search unless a changed coordinate map
changes that count. Both parity sectors and all original equations remain
required. The same front cap is configured for FEC9; its actual cost is not
predicted from the lower-weight control.

The generic identity-block shortcut remains useful independently, but it did
not change the tested FEC7 basis. The larger front, rather than an assumed
shortcut improvement, is the substantive next trial for coefficient growth.
The diagnostics now distinguish rows reconstructed at each CRT modulus from
rows whose big-integer residues are still retained. The partitioned output
reader also has temporary parsing allocations; retaining only one partition
is not a guarantee that peak RSS equals its final tensor payload.


## Fourth-attempt reconstruction milestone

At the 13:32 PDT follow-up, attempt4 was active in the first FEC8 parity
sector, with its original controller, executable and unchanged configuration.
The larger front reduced the completed triangular basis from 116,953,719
entries in attempt3 to 78,337,287 (about one third less). Complete modular
passes took about 397–405 seconds, versus about 484–491 seconds previously.
The back-substituted row basis remained large: 264,569,288 entries.

After four reconstruction primes (241 modulus bits), 1,064 of the 13,915
sector rows had rational candidates; 264,368,361 CRT entries were still
retained for the remaining rows. These are partial reconstruction candidates,
not an accepted sector certificate. A fifth modular pass was in progress.
Sampled peak RSS was about 25.35 GiB. The unrelated user workload occupied
less memory at this inspection, leaving about 20 GiB available on the host,
so the guard was not near its floor. The productive computation was left
unchanged. No complete FEC8 sector or FEC9 basis was claimed.

## Attempt4 outcome and global-pivot retry

Attempt4 stopped at 14:17:40 PDT on 4 October, after 4,748.66 seconds
(79.14 minutes), with **38.47 GiB sampled peak RSS**. The host-available-memory
floor triggered; the 50 GiB process ceilings did not. Both owned processes
were confirmed absent before the next launch. Neither FEC8 parity sector was
certified and the FEC9 stage never started. Full measurements are preserved in
`attempt4-result.json` and the original attempt directory.

The enlarged front reduced retained triangular nonzeros from 116,953,719 to
78,337,287 and each modular pass by roughly 18%, but its back-substituted matrix
still contained 264,569,288 entries. After ten primes (601 modulus bits), only
2,585 of 13,915 kernel rows were tentative rational reconstructions, containing
3,097,225 rational nonzeros; 261,410,139 entries still needed CRT integer
storage. The eleventh modular elimination completed before the guard stopped
reconstruction. This demonstrates that reducing finite-field elimination fill
alone did not resolve this run's main memory growth. Tentative reconstructed
rows are not accepted sector certificates.

Attempt5 launched at 14:20:37 PDT with the same frozen, tested executable and
assembled-parity FEC7 chain. FEC8 now uses global assembled pivot selection
within each character; FEC9 retains the bounded staged policy with a
268,435,456-entry front. This tests whether the global choice of free
coordinates avoids the large rational coefficients selected by staged
elimination. It is motivated by the existing ordinary FEC8 two-prime result
and assembled-parity FEC7 one-prime results, but neither proves that the new
FEC8 run will fit. Its assembled input and Schur complement can also have
larger memory peaks. All exact full-operator, completeness, output-hash and
known FEC8 dimension checks remain active. No universal default was changed.

The attempt5 config SHA256 is
`f05945905ee5c644a633725f8376de3a9b62dd3881e3f8d7f1076416eed7d7a6`.
Controller PID at launch was 1101114 and solver PID 1101117. The active pointer
and heartbeat now follow attempt5. The original deadline remains 5 October
09:22:13 PDT; launching a retry did not replenish the budget.

## Full parity-adapted FEC8 completed; FEC9 started

At the 14:55 PDT inspection, attempt5 had successfully completed its FEC8
stage and advanced automatically to FEC9 at 14:48:46 PDT. The FEC8 process
exited zero with no memory or timeout flags. Its complete measured duration
was **1,686.978 seconds (28.116 minutes)** and sampled process-tree peak RSS
was **28.5539 GiB**, including checkpoint writing and full serialization.

| FEC8 parity sector | Dimension | Final nonzeros | Reconstruction primes | Modulus bits |
| --- | ---: | ---: | ---: | ---: |
| Even | 13,915 | 123,804,193 | 2 | 121 |
| Odd | 3,556 | 88,114,153 | 3 | 181 |
| Full basis | 17,471 | 211,918,346 | — | — |

Both sectors carry the complete factored-operator residual and height
certificate. Full serialization and sealing succeeded. The completed basis
file is `output_fec9_20261004/attempt5/weight8/w8.wxf`, 2,493,666,600 bytes,
SHA256 `f83dcd770ea3a25cbbc1cf87bb74c50558100fb5d9c8ed800626eba0f517ff31`.
All 22 files in the seal were independently hashed and checked, along with
the prepared-bundle binding and exact byte identity of the copied FEC1–FEC7
recurrences/layouts. Measurements, dimensions, certificates and hashes are
recorded in `attempt5-fec8-result.json`.

For the same even FEC8 input, global assembled pivot selection reduced the
core kernel from 264,507,364 entries in attempt4 (its recorded rational plus
CRT nonzeros after ten primes) to 118,730,228 entries in attempt5. More decisively, all even-sector rows reconstructed at
121 bits rather than remaining mostly unreconstructed at 601 bits. The even
sector's complete solve took 558.346 seconds, whereas attempt4 stopped after
79 minutes without certifying it. This supports optimizing the selection of
free coordinates as well as sparse fill. It does not establish that assembled
elimination fits the much larger FEC9 problem. It also is not a claim that
this adapted FEC8 construction beats the earlier ordinary-basis construction
in elapsed time; those runs use different input coordinates and conditions.

The FEC9 solver PID was verified as 1102936, under controller 1101114 and the
unchanged config hash. Its first sector begins with 439,404 columns. At the
14:54 inspection it was still performing exact short-relation presolve; peak
RSS so far was 12.43 GiB, including input loading. No FEC9 sector certificate
or full FEC9 output exists at this milestone. The productive run is unchanged.

## Attempt5 stopped; bounded investigation after five-loop publication

Attempt5's FEC9 stage stopped on the host-memory floor after 18,696.763 s,
with 39,809,060,864 bytes (37.08 GiB) sampled peak RSS and exit -15. No FEC9
sector was certified. The final logged first-sector rank was 190,108 of
234,379 reduced columns, with 2,233,371,682 retained triangular entries.
Their coefficient/index payload alone is about 24.96 GiB; the factored
input, temporary modular factors and other allocations are additional.
The original controller and FEC9 solver PIDs were confirmed absent at the
post-publication check. The completed parity FEC8 remains preserved, with
the provenance recorded above. It is not the ordinary FEC8 used for the
five-loop publication.

The clean five-loop publication is complete, so the publication hold no
longer blocks this investigation. The original FEC9 deadline remains
5 October 09:22:13 PDT. No six-loop calculation is being launched.
`output_fec9_20261004/policy-controls-v2/` contains sequential, bounded
full-space controls for larger streaming batches and a broader initial
equation band. Both characters are retained. The experimental driver only
exposes existing shared-library batch options; no default policy is changed.

A separate storage change is being prepared: the factored recurrence already
clears a common denominator for each older component, but retains rational
objects with unit denominators. Storing those factors as arbitrary-precision
integers removes the redundant denominator word. The completed validation
and resulting retry are recorded below.

## Validated integer factors and sixth attempt

The recurrence factors now use arbitrary-precision integers after clearing
one common denominator for all letters of each older component. The exact
operator, modular values and residual-height bounds are unchanged. The
106-case factored extension/sewing suite passed, including new large-rational
and prime-denominator cases in both directions and both solving strategies.
The 303 generic staged/affine/composition cases, 120 coordinate-restriction
cases, and 16 product-invariance cases also passed, including false-candidate
and extra-prime checks. All saved weight-one through weight-seven recurrence
files were byte-identical to the independently certified controls.

For the same weight-seven input, options and output, fresh serial runs gave:

| Factor representation | Wall s | Peak GiB |
| --- | ---: | ---: |
| Rational objects after denominator clearing | 63.044 | 1.362 |
| Arbitrary-precision integer objects | 63.551 | 1.298 |

These single runs show lower memory with essentially unchanged time, not a
speedup claim. For the 211,918,346-entry saved FEC8 input, omitting one
eight-byte denominator word per entry removes 1.579 GiB of factor payload.
This does not reduce retained modular elimination fill. The original input
tensor and parser can still set an earlier loading peak.

The policy controls retained both characters and passed independent complete
weight-six and weight-seven space comparisons. At weight seven, larger
batches reduced output nonzeros from 20,579,497 to 19,176,253; a broader
initial equation band reduced them to 18,170,443 but slowed the odd sector.
The row-only pivot policy was stopped as an unpromising diagnostic: after
five primes (301 bits), its even sector was still not fully reconstructed;
the default policy certified it with two. This stopped pilot is not reported
as a complete weight-seven result. All measurements and the stop reason
are retained in `policy-controls-v2/`.

Attempt6 launched on **5 October 2026 at 00:27:35 PDT**, controller PID
1176632 and solver PID 1176634 at the first verified observation. It uses
the validated integer representation, support bound 1024, a 536,870,912-entry
initial front, 16,777,216-entry / 16,384-row batches, and pivot method 0.
The larger initial front can hold the complete 304,506,965-entry short band
observed in attempt5. Both parity sectors remain required. No expected FEC9
nullity is supplied and no universal pivot-policy default is changed.

The source, executable, guards, prepared bundle and complete FEC8 chain are
frozen in `output_fec9_20261004/attempt6/`. All 22 files in the copied FEC8
seal were rehashed, as were the prepared bundle and the full basis checksum.
The new config SHA256 is
`1ed1c0aec791f9770b893e420c344710f4bac6c68b48016768dd824034f90bbd`.
The active pointer and existing heartbeat follow this attempt. Its deadline
is still **5 October 09:22:13 PDT**, with 50 GiB address-space/RSS ceilings
and a 4 GiB host-available-memory floor. No budget was replenished. At the
first process/config verification, the first sector had begun with 439,404
columns; peak RSS so far was 12.43 GiB from loading. No FEC9 output or sector
certificate is claimed at this milestone.

## Attempt6 first discovery milestone

At the 01:04 PDT follow-up on 5 October, attempt6's controller and solver
were active with matching frozen binary/config hashes and unchanged resource
limits/deadline. The full 304,506,965-entry initial band fitted in one global
front (`bounded=0`). The first unrestricted 1/16 discovery sample finished
at rank 131,250 of 234,379 columns, leaving 103,129 candidate directions.
It took 1,480.83 seconds from the start of modular discovery, compared with
1,842.12 seconds at the same band/rank in attempt5. This is a roughly 20%
lower intermediate elapsed time in separate shared-host runs, not a complete
FEC9 performance result or a controlled ablation of either change.

The solve had moved to the 1/4 discovery sample. Peak RSS at this check was
about 21.63 GiB, with about 30 GiB available on the host. The latest emitted
triangular count was 895,928,441 nonzeros at rank 129,894; emission batches
differ from attempt5, so that count is not a matched-rank memory comparison
or a peak fill measurement. No sector receipt or full FEC9 output existed.
The productive run was left unchanged. The process/resource observation is
retained in `attempt6-observations.jsonl`.

## Attempt6 complete-row discovery pass started

At 04:45 PDT on 5 October, the first sector's unrestricted 1/4 sample had
finished at rank 191,801 of 234,379, leaving 42,578 candidate directions,
at 14,433.6 seconds from modular discovery start. The solver had started
the unrestricted full-row pass (`sample_denominator=1`). This is still
first-prime discovery, not a certified rational sector or complete FEC9.
The latest emitted batch before that transition had rank 190,202 and
1,945,173,760 retained triangular entries; it is not a peak-fill measurement.

The run remained active after about 4 h 18 min, with sampled peak RSS about
33.7 GiB. Host available memory was about 8 GiB at the initial check, above
the unchanged 4 GiB floor. The frozen config and executable hashes, process
commands, 50 GiB limits and original 09:22:13 PDT deadline were verified.
No sector receipt or full FEC9 output existed, and the run was left unchanged.

## Attempt6 stopped; retained elimination fill is still dominant

The host-memory guard stopped attempt6 at **04:59:51 PDT on 5 October**.
Its measured wall time was **16,336.398944 s (4 h 32 min 16 s)** and peak
process-tree RSS was **36,151,386,112 bytes (33.669 GiB)**. Exit status was
-15, `memory_pressure=true`, `rss_limit_exceeded=false`, `timeout=false`.
Thus the host's available memory fell below 4 GiB before either the 50 GiB
process limit or original deadline was reached. The controller and solver
were confirmed absent, and their frozen config/executable hashes still match
the launch record. Full measurements, log text and hashes are preserved in
[`attempt6-result.json`](attempt6-result.json).

The completed 1/4 discovery sample had rank 191,801 of 234,379 reduced
columns. Its 42,578 remaining directions are only a candidate upper bound,
not a certified sector dimension. The full-row first-prime pass had started
but had not completed. There are no FEC9 sector checkpoints, rational
certificates, complete basis, or recurrence seal. The complete certified
parity FEC8 and all prior attempts remain preserved.

Inspection of `streaming_echelon.hpp` and `staged_kernel.hpp` identifies the
current bottleneck: bounded input batches do not bound the retained triangular
row basis. Its last emitted count was 1,945,173,760 entries at rank 190,202,
about **21.74 GiB** at 8-byte residues plus 4-byte indices per entry, before
factor storage, temporary modular factor copies, batches and allocator costs.
That log point is not the maximum rank or peak fill. `ConstraintRows::modular`
also retains a full modular copy of the integer recurrence factors during
discovery. Rational reconstruction and final output writing had not started;
they cannot explain this stop and their eventual cost is not yet measured.

The integer-factor change removes redundant denominator storage; it does not
address most of this retained fill. Further substantial changes should target
the retained triangular representation, pivot ordering, or elimination blocks
that can be discarded or stored externally, together with restartable modular
checkpoints. Avoiding a simultaneous full modular factor copy is a smaller
additional opportunity. These are prospective directions, not validated
performance claims or implemented changes in this attempt.

At the 05:16 PDT inspection, about four hours of the original budget remained,
less than this incomplete first-prime run consumed. The driver checkpoints
only completed rational sectors, so it cannot resume this elimination state.
An unchanged restart is not justified by the observed cost; none was launched.
A further retry would require a substantive change validated on exact smaller
full spaces within the same remaining budget. The monitor remains active until
the original deadline, with no budget extension or changes to unrelated jobs.
