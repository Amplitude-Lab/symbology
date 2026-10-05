# Further reductions in rational-kernel time and peak RAM

The larger heptagon chains are about 12% faster than the previous symmetry
solver in this comparison. Heptagon LEC peak memory falls from 334.2 to
258.0 MiB, a 22.8% reduction; FEC falls from 77.0 to 68.0 MiB. Hexagon MHV LEC
is 11.7% faster and uses 16.7% less RAM. Every saved rational basis remains
byte-identical to the previous exactly verified version.

The current heptagon LEC chain is comparable to the ordinary solver in time
(1.449 versus 1.495 s) and uses less RAM (258.0 versus 344.1 MiB). The ordinary
solver was rebuilt with the same shared reconstruction and output-verification
improvements, so this comparison does not withhold those gains from it.

The [research review](RESEARCH.md) covers representation condensation,
SpaSM/FFLAS-FFPACK, LinBox, block Wiedemann, fixed-prime arithmetic, recent
exact CPU/GPU GEMM, fraction-free methods and allocation strategies. It clearly
separates local experiments from unimplemented candidates.

## Complete-chain wall time (seconds)

| Chain | Previous symmetry | Current symmetry | Ordinary, updated library | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 0.099 | 0.112 | 0.073 | — |
| Heptagon LEC ≤3 | 0.097 | 0.093 | 0.065 | — |
| Heptagon FEC ≤5 | 0.412 | 0.363 | 0.431 | — |
| Heptagon LEC ≤4 | 1.646 | 1.449 | 1.495 | — |
| Hexagon FEC ≤8 | 0.086 | 0.093 | 0.154 | 0.050 |
| Hexagon MHV LEC ≤6 | 0.222 | 0.196 | 2.873 | 0.215 |

## Peak resident memory (MiB)

| Chain | Previous symmetry | Current symmetry | Ordinary, updated library | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 20.9 | 21.0 | 21.1 | — |
| Heptagon LEC ≤3 | 22.9 | 20.7 | 18.3 | — |
| Heptagon FEC ≤5 | 77.0 | 68.0 | 111.4 | — |
| Heptagon LEC ≤4 | 334.2 | 258.0 | 344.1 | — |
| Hexagon FEC ≤8 | 16.0 | 15.9 | 28.3 | 18.6 |
| Hexagon MHV LEC ≤6 | 39.4 | 32.8 | 354.7 | 37.3 |

Each cell is the median of five fresh sequential processes with two requested
workers and rotating method order. The previous binary is exactly the current
binary from the [streaming audit](../symrep-stream-2026-10-02/REPORT.md). Both
versions were rerun; the table does not compare timings taken in different
sessions. No build or test job ran concurrently with these benchmark processes.
OS caches were not flushed and unrelated user processes were not stopped.

Time includes startup, loading, recursive solving, writing, exact semantic WXF
readback and completion checksums. Prepared inputs are reused. Optional full
expansion and reference verification run separately and are excluded. Peak RAM
is GNU time's maximum resident set for the entire process, including temporary
allocations and allocator retention. It is not the sum of tensor payload sizes.

The current and previous project executables and the rebuilt ordinary reference
use GCC 14.4.0, FLINT 3.2.2, `-O3 -march=native -mtune=native -flto`, and no
mimalloc. The archive executables retain their own SparseRREF, pivot method 1
and forced mimalloc, with the same compiler optimization flags. Their hashes
match the preceding audit. This remains an end-to-end implementation comparison.

Small-case elapsed times were noisy. A separate nine-pair alternating check,
using the same GNU time wrapper and exact basis comparison after stopping the
timer, measured:

| Chain | Previous symmetry | Current symmetry |
| --- | ---: | ---: |
| Heptagon FEC ≤4 | 0.0892 s | 0.0884 s |
| Hexagon FEC ≤8 | 0.0635 s | 0.0556 s |

This did not confirm the small-case regressions in the main five-run sample.
The diagnostic is retained in `small-cases-trials.json`; it does not replace
the main table or establish a guaranteed improvement. Hexagon FEC remains
slower than the unpacked package, with lower peak RAM. The two smaller heptagon
chains also retain startup/representation overhead relative to ordinary solving.

## What changed and why

**Incremental exact output verification.** Previously the WXF readback retained
one parser token for every entry and decoded another rational tensor. Heptagon
LEC weight 4 has 1,774,777 carrier nonzeros. Instrumentation put the old peak
near 334 MiB after writing/readback, before orbit selection began; modular
elimination itself had a smaller peak. The parser now exposes `parse_each`,
and the writer checks the header, dimensions, row offsets, indices and each
exact rational value as tokens arrive. It retains the input byte buffer, but
neither a per-entry token list nor another complete rational tensor. Checksums
do not replace semantic validation.

**Direct final-carrier output.** A read-only tensor-shaped view serializes the
final matrix without copying its rational values into an intermediate tensor.
Its indices and row offsets reproduce the existing file format byte-for-byte.
Intermediate weights still retain the tensor required by the next contraction.
The final carrier is freed before orbit selection as in the preceding version;
bad-prime retries can reload it. No LEC component expansion is introduced.

**Cheaper exact reconstruction.** Input row heights now avoid redundant rational
construction and denominator-one LCM work. For sufficiently large matrices,
row-height calculation and first rational reconstruction run in parallel using
per-worker accumulators. Small signed residues within FLINT's reconstruction
uniqueness bound become integers directly. Other residues still use FLINT.
The original full-input height certificate and multi-prime CRT fallback remain.

**Earlier release of unused rational rows.** After the first modular forward
elimination, original non-pivot rows are freed before back-substitution and
reconstruction. Their contribution to the input height bound has already been
recorded. Subsequent CRT solves use the same original pivot-row subset as
before; the old implementation delayed this release until reconstruction failed.

These changes do not alter the equations, pivot choices, irrep conventions or
basis normalization. They apply generally; no D3/D7 dispatch was added. They
also do not turn heptagon elimination into a multiplicity-only solve: its default
still computes a full sparse rational kernel and an exact irreducible orbit
recipe. Hexagon retains scalar multiplicities throughout its normal recurrence.

## Where the time went

The final heptagon LEC weight has 140,538 equation rows before empty-row removal,
30,828 columns and 3,964,499 nonzeros. After removal, 67,463 rows remain.
Its 3,745-dimensional kernel is unchanged.

| Final weight-4 phase | Previous median | Current median |
| --- | ---: | ---: |
| Constraint assembly | 0.159 s | 0.153 s |
| QQ kernel and normalization | 0.730 s | 0.656 s |
| Inherited modular actions | 0.068 s | 0.075 s |
| Orbit selection/frame recipe | 0.349 s | 0.367 s |

The intended arithmetic gain appears in the kernel phase; orbit selection was
not accelerated by this pass. Phase medians need not add to the median total,
and the table excludes loading, output and earlier weights. For hexagon LEC6,
the streamed assembly/solve pipeline changes from 0.1254 to 0.1191 s. Summed
worker durations are not wall time because sector jobs overlap.

The one-off instrumented profile showed first-prime success in the large
heptagon solve, with about 0.09 s spent calculating original row heights and
0.04 s in initial reconstruction before optimization. This is why improving
those steps was more relevant than parallelizing several CRT primes that this
example never needs. Profile builds are diagnostic and are not used for the
main table.

## Experiments retained and rejected

- The reconstruction-only and incremental-verification candidates were tested
  separately on both polygons; `optimization-trials.json` retains the results.
  Its exploratory wall times include a post-run basis comparison, so they must
  not be mixed with the main benchmark table.
- The first direct-output candidate lowered heptagon LEC RAM from about 266 to
  252 MiB but was initially slower. Hot-loop error-string construction was
  removed before the production measurement. `direct-output-trials.json`
  records that earlier candidate, not the final executable.
- Earlier row release had mixed exploratory timing results: it helped LEC,
  while FEC5 was slower in that three-run sample and used less RAM. The final
  full comparison measures the complete adopted change, not a claim that each
  individual edit speeds up every case. See `pruning-trials.json`.
- Smaller primes, sparse versus dense FLINT RREF, glibc arena limits and
  mimalloc were tested. None justified a new production default.
  [RESEARCH.md](RESEARCH.md) explains the inputs, results and limitations;
  raw results are in `modular-trials.json` and `allocator-trials.json`.

## Exactness, continuation and reproducibility

- The public ten-group algebra suite passed, including C3/C5/C7, D3/D5/D7,
  D7×C2, S4, C2×C2 and Q8, non-split/quaternionic cases and direct comparison
  with SparseRREF kernel extraction.
- A new 22,000 × 22,002 system with known nullity two forces the parallel
  reconstruction path and multiple CRT primes using large rational
  coefficients. Its full rational residual is zero. Signed reconstruction
  boundary values are checked on both sides of the fast-path cutoff.
- New WXF tests compare direct-view output with the original tensor encoder
  and independent full reader. They cover large rationals, empty tensors,
  changed values, truncation, trailing bytes, invalid indices and zero
  denominators. Existing WXF regressions pass 27 valid cases, 354 invalid
  cases and 1,000 bounded mutations.
- The ordinary numerical suite passes 90 tests. Public CLI regressions pass
  recursive FEC/LEC, one/three-worker consistency, continued runs, vector
  terminals, empty spaces, factorized export and bundle integrity.
- At every saved weight of all six benchmark cases, current and ordinary
  basis files match the preceding exactly verified audit byte-for-byte; all
  five repetitions also agree. Prepared-bundle checksums must match before
  the runner reuses those proofs. This pass reuses the earlier exact change
  of basis and residual certificates; it does not claim to have recomputed
  every expensive reference transformation.
- All WXF files from all five repetitions of each hexagon archive chain also
  match the preceding archive run. Separate real-data heptagon continuations
  from FEC4 to FEC5 and LEC3 to LEC4 match fresh current runs in every saved
  basis file. See `additional-validation.json`.
- Applying the full five-patch stack to clean pinned SparseRREF reproduces
  every production vendor header exactly. Setup `--check`, idempotent setup
  and `git diff --check` pass. The new arithmetic and visitor patches are
  included in `scripts/setup-sparserref.sh`.

`source-manifest.json` records build settings, source/binary hashes, validation
logs and the exact patch-stack check. `solver-changes.patch` records the solver
changes relative to the start of this pass. Complete source snapshots, output
bundles and process logs are in `output_symrep_research_20261002/release/`.
The reference proof chain ultimately leads to
[the full reference validation](../symrep-compact-2026-10-02/reference-validation.json).

```bash
python3 bench/symrep_orbits.py \
  --native output_symrep_research_20261002/native-current \
  --previous output_symrep_research_20261002/before/symrep \
  --prepared-root output_symrep_sparse_20261002/validated \
  --archive /home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap \
  --archive-bin output_symrep_orbits_20261002/archive-release \
  --references output_symrep_validated_20261001 \
  --validated-bases output_symrep_stream_20261002/release \
  --output output_symrep_research_reproduction --repeats 5
```

The output directory must be new. Omit `--validated-bases` to recompute the
full exact reference transformations. Changes remain uncommitted on
`symrep-solver-dev`.
