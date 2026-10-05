# Symmetry kernel solver: implementation and benchmark

**Historical baseline.** The performance problems below were subsequently
investigated and the recurrence optimized. See the
[2026-10-02 analysis and measurements](../symrep-2026-10-02/REPORT.md) for current
results. This report and its raw measurements are retained unchanged below.

The new `symrep` solver produces complete rational irreducible copies during
kernel solving and recursive FEC/LEC extension. Exact comparisons pass for the
unpacked hexagon package and the repository's heptagon data. **This first general
implementation is much slower than ordinary solving on the heptagon benchmark.**
It establishes the representation-preserving workflow; it is not yet a
performance replacement for the existing solver at large weights.

Implementation and usage: [docs/symrep.md](../../docs/symrep.md).
The work is on `symrep-solver-dev`, started at `front-end-dev`
`ac664bba85f887c6255b940556c9952f6ed02b67`. SparseRREF was not modified.

## Verified spaces

| Case | Dimensions at successive weights starting at one |
| --- | --- |
| Hexagon FEC | 3, 6, 13, 26 |
| Hexagon MHV LEC | 6, 25, 83 |
| Heptagon FEC | 7, 28, 97, 308 |
| Heptagon LEC | 14, 118, 734 |

Every timed run is checked by exact invertible changes of basis at every weight,
not just by dimension. Both hexagon chains are also compared with the actual
archive executables after expanding their multiplicity recurrences. The heptagon
chains additionally match independently generated chains from the unchanged
native `bootstrap --extend` command. Each returned recurrence satisfies its
canonical generator intertwining identities exactly.

Public checks cover C3/C5/C7, D3/D5/D7, S4, C2 × C2 and Q8, including a non-scalar
quaternionic commutant coefficient. They cover dense rational input conjugations,
tensor-product changes and inverses, sequential kernels, saved-kernel reuse,
forward/backward chains, vector terminals, zero/unconstrained kernels,
zero/unconstrained recursive spaces, continuation, damaged bundle rejection,
and non-invariant constraint rejection. `make check-symrep` passes and is included
in CI. The existing numerical regression reports 90 passes, with all ten recorded
two-/three-loop tensors matching exactly.

See [public checks](public-checks.txt), [native checks](native-checks.txt),
[heptagon FEC comparison](heptagon-fec-native-verify.txt), and
[heptagon LEC comparison](heptagon-lec-native-verify.txt).

## Timing and memory

Seconds below are **median process wall times of three runs**, using two
SparseRREF worker threads. A chain time includes weights 2 through the stated
maximum, group/model loading, mandatory checks and WXF I/O; weight one is a seed.
Preparation is a separate, reusable conversion. Exact reference verification is
outside the timed chain process. Peak RSS is the maximum among the three adapted
chain processes, in MiB.

| Case | Maximum weight | Ordinary chain (s) | Adapted chain (s) | Preparation (s) | Adapted peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hexagon FEC | 4 | 0.008679 | 0.018832 | 0.013747 | 12.6 |
| Hexagon MHV LEC | 3 | 0.015252 | 0.026228 | 0.014510 | 12.8 |
| Heptagon FEC | 4 | 0.075555 | 92.093124 | 0.344505 | 246.2 |
| Heptagon LEC | 3 | 0.064096 | 71.987206 | 0.331077 | 2123.6 |

The ordinary baseline uses the same executable, compiler, SparseRREF checkout,
rational field and worker count, retaining the original input basis at each
weight. It uses `symrep ordinary`, not the separately installed native command
whose default thread count is different. The native command is an independent
correctness reference.

Ordinary times do not include a post-solve symmetry decomposition. The table
compares the existing basis output with the new representation-organized output.

The specialized archive executables took median wall times of
0.011815 s for FEC through weight 4 and 0.025514 s for
MHV LEC through weight 3. They use their own bundled SparseRREF and fixed
parity × D3 recurrence implementation, so these are contextual measurements,
not an isolated comparison of decomposition algorithms.

The host is an Intel Core Ultra 7 265 with 20 logical CPUs, Linux x86-64,
GCC 14.4.0 and FLINT 3.2.2. The build uses `-O2 -std=c++20`, without
`-march=native` or a mimalloc override. Benchmark solver processes run
sequentially in fresh output directories. OS caches are not flushed, CPUs are
not pinned, and the shared workstation is not isolated; the millisecond-scale
hexagon timings should not be overinterpreted.

| Case | Individual adapted chain times (s) |
| --- | --- |
| Hexagon FEC | 0.023121, 0.018832, 0.016045 |
| Hexagon MHV LEC | 0.020622, 0.030769, 0.026228 |
| Heptagon FEC | 92.423943, 92.093124, 90.842872 |
| Heptagon LEC | 71.987206, 71.910275, 72.003809 |

Internal phase measurements (medians, seconds):

| Case | Weight | Dimension | Tensor adaptation | Reduced assembly | Kernel and reconstruction |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hexagon FEC | 2 | 6 | 0.000423 | 0.000174 | 0.000357 |
| Hexagon FEC | 3 | 13 | 0.000492 | 0.000276 | 0.000547 |
| Hexagon FEC | 4 | 26 | 0.000669 | 0.000446 | 0.000950 |
| Hexagon MHV LEC | 2 | 25 | 0.000829 | 0.000266 | 0.000668 |
| Hexagon MHV LEC | 3 | 83 | 0.003054 | 0.000754 | 0.001806 |
| Heptagon FEC | 2 | 28 | 0.011103 | 0.025485 | 0.006040 |
| Heptagon FEC | 3 | 97 | 0.027956 | 7.916442 | 0.063437 |
| Heptagon FEC | 4 | 308 | 0.130102 | 82.025457 | 1.651808 |
| Heptagon LEC | 2 | 118 | 0.018186 | 0.030158 | 0.016228 |
| Heptagon LEC | 3 | 734 | 0.056967 | 57.488733 | 13.580883 |

`Kernel and reconstruction` includes SparseRREF, integer basis shortening,
selection of whole irreducible copies and lifting the recurrence. It is not
just RREF time. `Reduced assembly` includes the source/target corner changes
needed to build each equation block. Raw records, memory measurements, source
hashes and executable hash are in [results.json](results.json).

## What the heptagon output contains

Cyclic + flip + parity generate D7 × C2, order 28, for the supplied matrices.
Its QQ irreps have dimensions `1,1,1,1,6,6`, and endomorphism-algebra dimensions
`1,1,1,1,3,3`. The six-dimensional models retain the three conjugate complex
2D irreps together. No algebraic extension of QQ is needed for this output.

The final copy counts in the automatically chosen model convention are:

| Irrep ID | QQ dimension | Endomorphism dimension | FEC weight-4 copies | LEC weight-3 copies |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1 | 1 | 2 | 22 |
| 1 | 1 | 1 | 2 | 30 |
| 2 | 1 | 1 | 26 | 33 |
| 3 | 1 | 1 | 14 | 19 |
| 4 | 6 | 3 | 4 | 53 |
| 5 | 6 | 3 | 40 | 52 |

The ordered generator matrices and rational characters identify each numeric
irrep ID; IDs are not universal group-theory names. The weighted copy counts
sum to 308 and 734. Layouts are retained in
[heptagon-fec-copies.tsv](heptagon-fec-copies.tsv) and
[heptagon-lec-copies.tsv](heptagon-lec-copies.tsv).

The prepared data format carries the original group, canonical irrep models,
letter/equation/terminal basis changes, the adapted dlog tensor and seed, and
copy layouts. A saved chain is bound to that input by checksums. Reusing the
original group matters when an intermediate kernel is not faithful; missing
sectors can reappear after tensoring with the alphabet.

## Performance interpretation and remaining limits

At heptagon FEC weight 4, the adapted recurrence has 197,644 nonzero entries
and coefficient height 11 bits in the first timed run. The ordinary recurrence
has 13,109 entries and height 3 bits. The adapted reduced equation blocks also
become substantially denser. Most of the FEC time is spent assembling them.
LEC has larger kernel and integer-basis reduction costs and reaches about
2.1 GiB RSS. Smaller blocks alone do not guarantee a faster sparse computation.

The most useful next optimization is to cache contraction templates between
irreducible factors and keep the recurrence in multiplicity coordinates during
assembly, as the specialized hexagon code does. The present implementation
already caches tensor-product decomposition maps and solves primitive corners,
but reconstructing and contracting component tensors at every weight remains
expensive. That optimization must preserve the same exact reference checks.

Automatic model construction has explicit limits: default group order at most
256, at most 1,024 subgroup projectors, and certified Schur-index-one or the
supported quaternionic Schur-index-two models. Unsupported constructions fail
explicitly. This is not a claim to implement every finite-group rational-model
construction. The archive benchmark covers MHV LEC, not NMHV. Higher-weight
production use, sewing and GUI integration were not benchmarked or added.

## Reproduction and retained artifacts

Build and run the commands in [the guide](../../docs/symrep.md). The benchmark
entry point is [bench/symrep_benchmark.py](../../bench/symrep_benchmark.py), with
the archive exporter [hexagon_symrep_reference.cpp](../../bench/hexagon_symrep_reference.cpp).
It snapshots the tested solver and sources before running. The unpacked input
archive was `hexagon_irr_bootstrap_v2_20260928.tar.gz`, SHA-256
`7f25b6f754f4b9bf4d297eba577bf14c6412700d427272289fcb4e4accb51a97`.

The complete local datasets, per-process logs, prepared bundles, chains and
exact basis-change certificates are retained under
[`output_symrep_validated_20261001/`](../../output_symrep_validated_20261001/).
Earlier experimental runs are separate and are not used in these tables.
