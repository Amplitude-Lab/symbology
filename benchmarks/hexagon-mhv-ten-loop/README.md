# Ten-loop MHV hexagon: native result, physical validation and benchmarks

The **second fresh V27D ten-loop prime** also matches every native and
shifted reference coefficient and passes the recorded seam and physical
checks. It uses 23.329189 CPU-hours, 6,540.289 wall seconds and 47.601 GB
sampled peak RSS. Against the learning prime, CPU falls only 0.217%: final-seam
replay saves 49.27% in that small phase, while discovery remains essentially
unchanged. These are different primes without a same-prime replay-off control.
The service recorded soft-limit events, with no hard-limit/OOM events.
Prime 3 is active; the new chain's rational reconstruction and exact
certificates remain pending. See `evidence/v27-ten-prime2` and `VALIDATION.md`.

A fresh **V27D ten-loop learning prime** now reproduces all 275,246 native
and 518,048 shifted coefficients modulo the first prime. It uses 23.380
CPU-hours and 6,599.910 wall seconds, versus historical V16's 37.286 CPU-hours
and 12,711.891 seconds. Sampled peak RSS rises from 34.628 to **46.188 GB**,
so the combined CPU/RSS target fails despite the CPU saving. This compares
different storage/cache settings and one learning prime. The first replay
measurement is summarized above; rational reconstruction/exact certificates remain pending;
the existing certified references below remain intact. See `VALIDATION.md`
and `evidence/v27-ten-prime1` for the independent receipt review.

The current working nine-loop benchmark is **V27D with 16 workers**.
The qualification target is more than 15% CPU saving against **V16**, with
at most 5% RSS growth; V26 is an incremental comparator. The earlier
recommendation incorrectly reset that target against V26. Direct nine-loop
qualification against V16 remains **unverified**: the historical measurements
are V17 and retain that identity. The selection amendment is in
`evidence/v27-selection-amendment`; the original V27 receipts remain historical
measurement evidence.

All 27 V27-campaign native runs pass coefficient and physical checks. V27D/16
uses 0.74% fewer CPU-hours and 8.65% less wall time than fresh V26, with 1.65%
higher peak RSS: 0.867992 CPU-hours, 297.189 seconds and 6.010 GB per prime.
The small CPU difference comes from one learning/replay pair, not timing
repetitions. This working selection establishes neither the missing
nine-loop V16 qualification nor the performance of a complete ten-loop
multi-prime reconstruction.

Physical constraints must be available from the start as equations or lazy
block generators. Their assembly and elimination may be deferred according
to cost, expected rank gain, fill and memory, preserving every required
equation and its affine RHS/origin. This requirement does not mean solving
all physical blocks first. A generic pending-block scheduler is still a
future implementation comparison.

The historical optimization review, **V26**, is included in `VALIDATION.md` and
`evidence/v26-generic`. Its nine-loop learning-plus-replay pair averages
0.872271 CPU-hours and 324.935 wall seconds per prime, at 5.943 GB peak RSS
(the larger of sampled RSS and kernel high-water mark). This is 41.13% less
CPU time than historical V17 and 4.32% less than a fresh matched V22 pair.
Eight-loop controls average 0.097837 CPU-hours, 47.654 seconds and 1.883 GB.
These are two different primes, not timing repetitions. V26 ten-loop timing
has not been measured; the ten-loop results and costs below retain their
original producer identities.

The primary reference is the rational **17+3** symbol in `reference/17p3`.
The equivalent **16+4** symbol is in `reference/16p4`. Each consists of six
D3-sector matrices, together with the precise coordinate order and recurrence
basis hashes. These are factorized ordinary symbols, with the E-component
metric `(1,3)` and the package's physical-copy-primitive normalization.
They do not determine function-level constants.

| Representation | Coordinates | Nonzero entries | Maximum numerator / denominator bits |
| --- | ---: | ---: | ---: |
| Native 17+3 | 275,246 | 246,941 | 45 / 315 |
| Shifted 16+4 | 518,048 | 491,197 | 44 / 231 |

The native export was reconstructed directly from seven saved native prime
solutions and checked against the eighth. These are retrospective checks.
An exact rational 17+3 to 16+4 conversion then reproduced every coefficient
of the separately certified shifted result, verified using Wolfram's WXF
reader. The shifted candidate was frozen after six CRT primes and validated
against two subsequently computed primes. Prime 1 was computed by V16 and
imported with verified provenance; primes 2–8 were computed by V17.

Uniqueness remains conditional on extrapolated double-final-entry bonus
relations (Dixon–Liu, equations 76–79 and their D3 orbit). The exact checks do
not prove unrestricted uniqueness or the bonus extrapolation itself.

## Reproducing a complete native prime

This benchmark uses precomputed FEC1–17, MHV LEC1–3 and physical input data.
It starts a new solve directory with no saved solve checkpoints. All basis
hashes are verified before calculation. The reference coefficients are read
only after the solve, for an exhaustive comparison. A reference vector is
never a fitting input.

On the current workstation, run from the repository root:

```bash
python3 bench/hexagon_ten_loop_prime.py \
  --output output/hexagon-ten-prime2 --prime-index 2 --hours 8
```

The wrapper uses the frozen V17 executable and pinned recovered runtime.
It enforces 16 threads, a 48 GB cgroup pressure threshold, a **50,000,000,000
byte hard memory limit**, no swap, and a 3 GiB host reserve. It refuses to
replace an existing output directory and retains logs on failure. A plan can
be inspected without launching the solve by adding `--plan-only`.
The new wrapper's plan and input validation were checked locally; the archived
complete calculations were produced by the earlier hash-pinned supervisor.
A new four-hour complete prime is not included in the additional-validation
session merely to test this wrapper.

For relocated inputs, use `--fec`, `--lec`, `--data`, `--solver`, `--library`
and `--reference`. The supplied sources are under `solver-source`, including
the package's normalization instructions and original V17 driver. They need
a C++20 compiler with `<format>`, FLINT 3, MPFR, GMP and the allocator
libraries specified by the package. Use the recorded production executable
for a byte-identical executable comparison; a rebuild is a separate binary
and must retain its own build receipt.

This is a **warm-basis sewing benchmark**. Unlike the five-loop heptagon
benchmark, it does not yet offer a validated fresh-clone, cold regeneration
of every hexagon basis and external physical input. Those large inputs are
identified by hashes rather than silently treated as free computation.
A full cold ten-loop timing has not been measured.

## Controlled final-seam replay

The last sewing phase can be measured without rerunning discovery and the
physical fit. Both variants start from the same post-OPE checkpoint at prime
1152921504606847067, use the identical frozen finishing function, and check
all 275,246 final coefficients against the committed native result.

```bash
python3 bench/hexagon_ten_loop_benchmark.py replay \
  --output output/hexagon-replay-mapped --storage mapped
python3 bench/hexagon_ten_loop_benchmark.py replay \
  --output output/hexagon-replay-heap --storage heap
```

Run these serially on an otherwise idle workstation. Each run preserves its
command, source and checkpoint hashes, CPU model, systemd limits, log, wall
and CPU time, sampled peak process-tree RSS, and the timed finishing phase.
Input loading is included in process wall time but excluded from the reported
finishing-phase time. Compiler time is excluded. A successful replay is not
an additional CRT prime and is not a measurement of a complete new solve.

## Recorded production cost

The machine is an Intel Core Ultra 7 265; numerical jobs used 16 threads.
GB below means decimal bytes divided by 10^9. RSS was sampled every 100 ms.
See `production-measurements.json` for all eight individual measurements,
CPU times and source report paths.

| Measurement | Wall time | Sampled peak RSS |
| --- | ---: | ---: |
| V16 native pilot, prime 1 | 12,711.891 s (3.53 h) | 34.628 GB |
| V17 native primes 2–8, observed range | 12,627.650–13,382.666 s | 38.513–39.141 GB |
| V17 native mean | 12,991.948 s (3.61 h) | — |
| All eight successful native primes | 103,655.526 s (28.79 h) | 39.141 GB |
| Eight modular 17+3 → 16+4 conversions | 537.176 s (8.95 min) | approximately 12.02 GB |

These totals exclude basis/input preparation, failed or deliberately stopped
experiments, reboot downtime, reconstruction and later certification. CPU
time across the eight successful native solves was 1,116,756.594 s, about
10.77 CPU cores on average. Low utilization in an individual phase should
therefore be assessed separately from full-run utilization.

The additional physical checks, exact native certificates and controlled
backend measurements are recorded in `VALIDATION.md` and their JSON receipts.

## Completed validation and allocator comparison

Both rational representations pass their complete original seam exactly,
the 44 exact origin predictions, and exhaustive exact conversion between
splits. The primary native result additionally passed **9,072** finite-field
collinear/MRK/OPE checks at two new primes, including **912 distinct unused
OPE coefficients**. These physical checks do not certify every limit or
function-level constant. Detailed coverage and controls are in `VALIDATION.md`.

| Finishing replay | Phase seconds | Process seconds | Peak RSS, GB |
| --- | ---: | ---: | ---: |
| Mapped storage | 1276.920 | 1349.263 | 13.959 |
| Heap storage | 222.071 | 289.422 | 14.410 |

Heap storage is 5.75× faster for this phase and reproduces every native
coefficient. This is one serial pair on the same checkpoint with warm caches;
a whole-prime speedup and a globally optimal solver are not established.
Prefer `--storage heap` for reproducing this finishing benchmark. Large
full-prime discovery retains the frozen production memory policy.

`manifest.json` binds the reference matrices, source snapshots, scripts and
evidence. `additional-measurements.json` records the independent artifact
review and detailed times, CPU usage, peak RSS and validation coverage. The
technique report remains in `audits/hexagon-ten-loop-2026-10-06/TECHNIQUE.tex`
at the repository root.

A later **V18 lower-loop comparison** is now included in `VALIDATION.md` and
`evidence/v18-lower-loop`: matched complete-prime savings are about 3% at eight
and nine loops; nine-loop OPE-cache tuning raises the combined saving to 8.6%
at 5.984 GB peak RSS. These are separate O3/heap controls, not ten-loop V18
measurements or factors to multiply by the ten-loop allocator speedup.

The subsequent **V22 nine-loop study** adds repeated CPU measurements:
1.481645 → 0.918766 mean CPU-hours (37.99% less) and 508.864 → 347.391 mean
wall seconds, at 5.840–5.890 GB sampled peak RSS. See `VALIDATION.md` and
`evidence/v22-nine-loop`. This remains a nine-loop result; ten-loop V22
performance is unmeasured. The validation record also corrects earlier
reference-selection wording and verifies cache-free rational-reference
comparisons at both nine-loop primes.
