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

Completed exact validation: both representations pass their complete original
seam, the 44 exact origin predictions per split, and exhaustive exact
conversion between the splits. The native result additionally passed **9,072**
finite-field collinear/MRK/OPE checks at two new primes, including **912
distinct unused OPE coefficients**. These checks do not certify every limit or
function-level constant. Uniqueness remains conditional on extrapolated
double-final-entry bonus relations (Dixon–Liu, equations 76–79 and their D3
orbit); the exact checks do not prove unrestricted uniqueness or the bonus
extrapolation itself.

## How to repeat the calculation

Two reproduction levels are packaged. Both are **warm-basis** benchmarks:
they consume the archived hexagon packages on the workstation — the FEC1–17
and MHV LEC1–3 bases, physical input data, the frozen V17 solver executable
and its recovered runtime — identified by hashes in
[manifest.json](manifest.json). A validated cold regeneration of every
ten-loop basis is not part of this package; those inputs are never silently
treated as free computation.

Prerequisites: the repository, Python 3, and the local archives at their
default paths (`output_hexagon_16p4_20261006/`, `output_hexagon_17p3_20261007/`,
`output_hexagon_runtime_recovery_20261008/`; relocatable with `--fec`, `--lec`,
`--data`, `--solver`, `--library`, `--reference`). A complete prime needs
about 4 h on 16 threads under a 50,000,000,000-byte hard memory limit.

### 1. Reproduce one complete native prime (about 4 h, 16 threads, 50 GB cap)

From the repository root:

```bash
python3 bench/hexagon_ten_loop_prime.py \
  --output output/hexagon-ten-prime2 --prime-index 2 --hours 8
```

Every basis hash is verified before the solve; no solve checkpoints or
reference coefficients enter the calculation. After the solver exits, all
275,246 native coefficients are compared exhaustively against the committed
reference. **Success criterion:** exit code 0 with the comparison receipt in
the output directory. The wrapper refuses an existing output directory,
enforces 16 threads, the 50 GB hard limit, no swap and a 3 GiB host reserve,
and retains logs on failure. `--plan-only` validates inputs and prints the
plan without launching the solve.

### 2. Controlled final-seam replay (minutes)

```bash
python3 bench/hexagon_ten_loop_benchmark.py replay \
  --output output/hexagon-replay-heap --storage heap
```

Replays only the last sewing phase from the shared post-OPE checkpoint at
prime 1152921504606847067 with the identical frozen finishing function, then
checks all 275,246 final coefficients against the committed native result.
Prefer `--storage heap`: 222 s versus 1,277 s for `mapped` in the recorded
serial pair, reproducing every coefficient. Run variants serially on an
otherwise idle workstation. A successful replay is **not** an additional CRT
prime and not a measurement of a complete new solve.

Solver sources, the original V17 driver and the package's normalization
instructions are under [solver-source](solver-source). They need a C++20
compiler with `<format>`, FLINT 3, MPFR, GMP and the allocator libraries
specified by the package. A rebuild is a separate binary and must retain its
own build receipt; use the recorded production executable for a byte-identical
executable comparison.

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
| Eight modular 17+3 → 16+4 conversions | 537.176 s (8.95 min) | approximately 12.02 GB |

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
