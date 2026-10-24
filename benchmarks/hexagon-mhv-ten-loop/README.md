# Ten-loop MHV hexagon: committed result and reproduction benchmark

## Result

The primary reference is the rational **17+3** symbol in
[reference/17p3](reference/17p3); the equivalent **16+4** symbol is in
[reference/16p4](reference/16p4). Each consists of six D3-sector matrices
together with the precise coordinate order and recurrence basis hashes.
These are factorized ordinary symbols with the E-component metric `(1,3)`
and the package's physical-copy-primitive normalization. They do not
determine function-level constants.

| Representation | Coordinates | Nonzero entries | Maximum numerator / denominator bits |
| --- | ---: | ---: | ---: |
| Native 17+3 | 275,246 | 246,941 | 45 / 315 |
| Shifted 16+4 | 518,048 | 491,197 | 44 / 231 |

The native export was reconstructed from seven saved native prime solutions
and checked against the eighth. An exact rational 17+3 → 16+4 conversion then
reproduced every coefficient of the separately certified shifted result,
verified using Wolfram's WXF reader. The shifted candidate was frozen after
six CRT primes and validated against two subsequently computed primes.
Prime 1 was computed by V16; primes 2–8 by V17.

**Current selection versus historical replay.** The later
[current-benchmark.json](current-benchmark.json) selects V27D/16, native 17+3,
and explicitly marks the older manifest state historical. Its coefficients are
identical to these committed references, but its first-prime RSS growth fails
the 5% qualification tolerance; selection is not a combined performance-gate
pass. The packaged complete-prime and final-seam commands below reproduce the
**historical V17 route**, not V27D. They are valid accepted replays, not
abandoned pipelines and not evidence of a fresh V27D run.

Completed exact validation: both representations pass their complete original
seam, the 44 exact origin predictions per split, and exhaustive exact
conversion between the splits. The native result additionally passed **9,072**
finite-field collinear/MRK/OPE checks at two new primes, including **912
distinct unused OPE coefficients**. These checks do not certify every limit or
function-level constant. Uniqueness remains conditional on extrapolated
double-final-entry bonus relations (Dixon–Liu, equations 76–79 and their D3
orbit); the exact checks do not prove unrestricted uniqueness or the bonus
extrapolation itself.

## Environment

All recorded production and replay measurements ran on one workstation:
Intel Core Ultra 7 265 (20 cores, no SMT), 62 GiB RAM, Ubuntu 22.04.5 LTS,
kernel 6.8.0-138-generic, NVMe storage. Numerical jobs used 16 threads under
a 50,000,000,000-byte hard memory cgroup limit with swap disabled and a 3 GiB
host reserve; RSS was sampled every 100 ms. The frozen solver binary and its
runtime libraries are kept locally under `inputs/` and are deliberately not
pushed (see below); rebuilding them from `solver-source/` needs the repository
[environment-linux.yml](../../environment-linux.yml) toolchain
(gxx_linux-64 14.4.0, libflint 3.2.2, gmp, tbb). No Wolfram installation
or Git LFS is required; the ten-loop warm-basis archives noted above are
the only non-repository dependency.

### Regenerating the bases from committed sources (cold path)

Everything that *generates* the warm-basis inputs is committed in this
repository:

- `solver-source/src/reduced_d3/` — the base and sewing programs
  (`irr_fec`, `irr_lec`, `irr_sew`, `irr_convert`, `irr_conditions`,
  `irr_pipeline`, `irr_cache`);
- `solver-source/Makefile` and `solver-source/run_workflow.sh` — the build
  and stage entries (`STAGE: FEC | LEC | SEW | CONDITIONS | SHIFT | RELEASE`);
- `solver-source/data/` — the dlogmat seed conditions and cyclic matrices
  (60 MB), with per-file checksums;
- [bench/hexagon_joint_build_v17.py](../../bench/hexagon_joint_build_v17.py)
  — builds the frozen V17 solver from the committed package sources with
  the environment-linux.yml toolchain;
- [bench/hexagon_seam_replay.cpp](../../bench/hexagon_seam_replay.cpp) —
  the seam-replay binary.

A recorded audit already re-ran `irr_fec 16` from the seed data (422 s wall,
1,169 CPU-s, 13.7 GB peak — `audits/hexagon-ten-loop-2026-10-06/fec16.json`),
which proves the regeneration mechanism; a complete cold rebuild of the
FEC1–17 chain has not been timed end-to-end and was excluded from the
production cost table. `solver-source/SOURCE_SHA256.tsv` pins the package's
build-time source state; five source files in the committed snapshot carry
later post-freeze edits (both states are documented by the TSV and the
manifest). The OPE small-s checkpoint and the per-prime joint caches remain
production outputs: a new prime computes its own cache via the solver, and
the committed `reference/` symbols are the certified endpoints compared
against.

## Current benchmark: V27D/16, native 17+3

V27D is the current selected benchmark, not V17. The selection receipt is
[current-benchmark.json](current-benchmark.json); its full audit records eight
fresh V27D native prime solves, seven CRT primes plus one retrospective check,
2,053,259 exact native seam rows, 44 exact origin checks, independent WXF
decoding and exact split equality. The selected V27D native WXF files are
byte-identical to `reference/17p3`, so the committed reference remains the
correct coefficient target. Selection does not erase the recorded first-prime
RSS-qualification failure.

Validate the current V27D installation without launching a new solve:

```bash
python3 bench/hexagon_ten_loop_runner_v27d.py --preflight-only
```

This checks the V27D/16 profile, solver SHA256
`493bbc1a5712a27060e3f3641b75d1820ff031fb4658d022e7c4df733f2d1091`,
frozen package, controls and runtime. It passed in the 2026-10-24 audit. The
full V27D producer is `bench/hexagon_ten_loop_runner_v27d.py` under its
installed 50 GB service; it is not the V17 wrapper below. No fresh V27D solve
was launched in this audit.

## Historical V17 warm-basis replays

Two historical reproduction levels are packaged. Both are **warm-basis** computations:
they consume the FEC1–17 and MHV LEC1–3 bases, physical input data, the
frozen V17 solver executable, its runtime libraries and the hash-bound
checkpoints — about 10 GB that is **not pushed to the remote** (too large
for the repository). Two local locations work interchangeably: the original
`output_hexagon_*` archives, or the split xz copies at
`inputs/` (4.1 GB in 48 parts plus the recorded seam-replay binary, pinned
in `inputs/SHA256SUMS`; the runners materialize them automatically into
`inputs/.materialized/`). On a machine without these archives the two
levels below cannot run at warm-basis speed; the cold regeneration path
from committed sources is described in its own subsection above. The bases
are never silently treated as free computation, and that dependency is
stated here explicitly rather than hidden. Everything else in this
benchmark — both committed rational symbols, all receipts, certificates and
evidence — is in the repository.

On the workstation the runners find the archives automatically; add
`--source bench` to force the `inputs/` split copies.

On the workstation with the original archives under `output_hexagon_*` the
runners use those directly; add `--source bench` to force the local split
inputs. **A fresh clone does not contain `inputs/` or `SHA256SUMS`.** Install
the split parts separately before using `--source bench`, or provide the
original archives; warm-basis commands cannot work from Git alone.

### 1. Reproduce one complete native prime (about 4 h, 16 threads, 50 GB cap)

From the repository root:

```bash
python3 bench/hexagon_ten_loop_prime.py \
  --output output/hexagon-ten-prime2 --prime-index 2 --hours 8
```

Every basis hash is verified before the solve; no solve checkpoints or
reference coefficients enter the calculation. After the solver exits, all
275,246 native coefficients are compared exhaustively against the committed
reference. **Success criterion:** exit code 0, `FULL_NATIVE_PRIME_BENCHMARK_PASS`
and the comparison receipt in the output directory. The wrapper refuses an
existing output directory, enforces 16 threads, the 50 GB hard limit, no
swap and the host reserve, and retains logs on failure. `--plan-only`
validates inputs and prints the plan without launching the solve.

### 2. Controlled final-seam replay (minutes)

```bash
python3 bench/hexagon_ten_loop_benchmark.py replay \
  --output output/hexagon-replay-heap --storage heap
```

Replays only the last sewing phase from the shared post-OPE checkpoint at
prime 1152921504606847067 with the identical frozen finishing function, then
checks all 275,246 final coefficients against the committed native result.
Prefer `--storage heap`: 222 s versus 1,277 s for `mapped` in the recorded
serial pair (a re-run under CPU-frequency throttling and a concurrent job
measured 324 s), reproducing every coefficient. Run variants serially on an
otherwise idle workstation. A successful replay is **not** an additional CRT
prime and not a measurement of a complete new solve.

Solver sources, the original V17 driver and the package's normalization
instructions are under [solver-source](solver-source). They need a C++20
compiler with `<format>`, FLINT 3, MPFR, GMP and the allocator libraries
specified by the package. A rebuild is a separate binary and must retain its
own build receipt; the local `inputs/` solver binary is the recorded
production executable for byte-identical comparison. The build helper accepts
`CXX`, `CONDA_PREFIX` and `LDLIBS`; it no longer depends on a removed `/tmp`
toolchain.

### Current-source replay audit (2026-10-24)

The current V27D installation passed `--preflight-only`; its previously
completed eight-prime selection and exact certificates were reviewed without
launching another solve. Separately, the historical V17 native-prime
`--plan-only` check passed with the local split inputs. The heap final-seam replay then passed with
`all_coefficients_equal=true` on all **275,246** columns and unchanged input
hashes: phase 290.212 s, whole measured process 364.334 s, sampled peak RSS
14.280 GB. See `output/hexagon-replay-heap-20261024/manifest.json`. The V17
build helper also compiled successfully with GCC 14.4.0 and the recovered
runtime. These checks validate the historical V17 replay only; they do not
claim a fresh complete prime, cold FEC1–17 regeneration, or a V27D run.

## Measured production cost

Measured on an Intel Core Ultra 7 265, 16 threads; GB means decimal bytes
divided by 10^9; RSS sampled every 100 ms. Full per-prime records are in
[production-measurements.json](production-measurements.json).

| Measurement | Wall time | Sampled peak RSS |
| --- | ---: | ---: |
| V16 native pilot, prime 1 | 12,711.891 s (3.53 h) | 34.628 GB |
| V17 native primes 2–8, observed range | 12,627.650–13,382.666 s | 38.513–39.141 GB |
| V17 native mean | 12,991.948 s (3.61 h) | — |
| All eight successful native primes | 103,655.526 s (28.79 h) | 39.141 GB |
| Eight modular 17+3 → 16+4 conversions | 537.176 s (8.95 min) | 12.035 GB |

These totals exclude basis/input preparation, failed or deliberately stopped
experiments, reboot downtime, reconstruction and later certification. Total
CPU time across the eight successful native solves was 1,116,756.594 s, about
10.77 CPU cores on average; low utilization in an individual phase should be
assessed separately from full-run utilization. A full cold ten-loop timing has
not been measured.

## Records

`manifest.json` binds the reference matrices, source snapshots, scripts and
evidence. [VALIDATION.md](VALIDATION.md) is the validation record **and the
changelog**: exact checks, physical-checker controls, adaptive rational
reconstruction, and the measured solver-optimization campaigns (V18, V22,
V26, and the V27/V27D ten-prime studies) with their receipts under
[evidence](evidence). The technique report is
[audits/hexagon-ten-loop-2026-10-06/TECHNIQUE.tex](../../audits/hexagon-ten-loop-2026-10-06/TECHNIQUE.tex);
`additional-measurements.json` records the independent artifact review.
