# Whole-strategy audit and four-loop heptagon acceptance test

Time and RAM were **not** at a lower bound. The new general kernel path and
a change in projection order construct the unique four-loop, weight-eight
heptagon MHV symbol on this machine in **68.5 seconds**, with **5.05 GiB maximum
single-process RSS**. Sampling the parent and descendants together recorded
**5.27 GiB summed RSS**. The calculation starts from the supplied seed data in
a fresh directory, including lower-loop recursion, basis construction,
symmetry, collinear solving, expanded collinear outputs and output files.

The result is an exact recursive symbol over the original 42-letter alphabet,
not only its collinear limit. Its definition and dependencies are recorded in
`output_solver_strategy_20261002/amplitude-all-streamed/result-manifest.json`.
The final recursive tensor is about 196 KiB; its forward basis dependencies
are additional storage, including a 124 MiB FEC7 file. No claim is made that
this is the globally fastest solver or that every matrix benefits.

## Measurements and what is being compared

Intel Core Ultra 7 265, 20 logical CPUs, 62.1 GiB RAM, no swap; eight requested
workers. GCC 14.4.0, FLINT 3.2.2, native optimization and LTO, no mimalloc.
Large runs used a 28 GiB per-process virtual-address-space limit and a 900-second timeout.
Builds and other benchmarks did not run concurrently with these measurements.
OS caches were not flushed and unrelated user processes were not stopped.

| Complete FEC1–FEC7 chain | Wall time | Peak process RSS |
| --- | ---: | ---: |
| Existing locally patched SparseRREF v0.3.6, ordinary recurrence | 27.81 s | 9.08 GiB |
| Isolated upstream v0.4.2 trial, ordinary recurrence | 36.97 s | 12.45 GiB |
| New streamed factored kernel | 18.60 s | 2.04 GiB |

The chain probes load inputs, construct every weight, write outputs and seal
the result with exact readback/checksums. These are individual large runs,
not five-run medians. The new path changes basis; equality means equal symbol
spaces, not equal coefficient files. The original and upstream trials happen
to produce byte-identical files at every weight. The newer upstream trial uses
the compatible local WXF layer but does not include all local reconstruction
optimizations; it is a deployment comparison, not an isolated comparison of
their finite-field elimination kernels. It does **not** establish that all
upstream improvements are worse.

| Complete four-loop calculation | Wall time | Peak process RSS | Sampled summed RSS |
| --- | ---: | ---: | ---: |
| First successful revision: streamed sewing, ordinary extensions and full projection | 215.48 s | 8.72 GiB | 8.60 GiB |
| Streamed extensions and sewing, restricted projection | 68.45 s | 5.05 GiB | 5.27 GiB |

The first row already includes fixes and a streamed sewing implementation from
this investigation. It is **not** the untouched original four-loop program.
Its collinear solve subprocess used the old default worker count; the large
extension, sewing and projection subprocesses used eight. The final command
uses eight throughout. The process monitor adds about 0.08 s of polling latency,
hence its reported totals of 215.51 and 68.54 s. Summed RSS is sampled every
100 ms and double-counts shared pages; it can miss a short-lived peak, explaining
why it is below the GNU-time high-water mark in the first row.

In the final run the FEC7 extension took 16.76 s, sewing 32.45 s, and the prefix
projection subprocess 4.82 s. The earlier full projection subprocess took
116.09 s. These subprocess phases exclude some parent work and do not sum to
the full run time. Sewing remains the largest cost. A subsequently interrupted
measurement overlapped reference verification and is explicitly excluded.

A second fresh, sequential optimized calculation completed in 68.02 s
(68.15 s including monitor polling), with 5.26 GiB sampled summed RSS. It
passed the same physical checks. The result is therefore reproducible at
approximately 68 seconds in two fresh runs, not a cache-resume timing.

## What the original solver actually does

I inspected the local and current upstream implementations, including pivot
search, forward elimination, Schur updates, back substitution, reconstruction,
kernel extraction and tensor allocation. The original implementation is
already a sophisticated sparse modular solver: it removes singleton rows,
selects batches of sparse pivots using row/column incidence, updates the Schur
remainder in parallel, uses Shoup modular multiplication, back-substitutes,
and reconstructs rational RREF coefficients.

The avoidable cost is around and before that algorithm:

1. Tensor contractions materialize a large rational equation matrix.
2. A finite-field copy exists alongside rational input/reconstruction storage.
3. Many forced zeros and scaled-identical unknowns are discovered only after
   the large product has been assembled.
4. The pipeline asks for rational RREF although its required output is a kernel.
5. Collinear projection reduces the entire final forward space although the
   final amplitude occupies only a few combinations of it.

The ordinary FEC7 equation matrix has 1,085,001 rows, 107,310 columns and
157,271,724 nonzeros. With the new recursive basis and structural presolve,
the materialized finite-field core has 289,085 nonzero rows, 37,093 columns
and 53,753,676 nonzeros. Its nullity is still 6,826. This is a general reduction
in the actual problem submitted to elimination, not a faster parser alone.

An intermediate trial contracted the last-entry boundary first but still
materialized all rational sewing equations: 1,200,850 rows, 95,564 unknowns,
566,934,084 nonzeros. It failed with SIGSEGV under the 28 GiB virtual-memory
cap, after reaching about 26.7 GiB RSS. Unchecked allocation failure is a
plausible cause, but was not independently proved. It is recorded as a failed
trial, not as a valid timing or a proof that the original algorithm cannot
finish on a larger machine. Streaming presolve made that allocation unnecessary.

## Implemented general changes

**Exact structural presolve.** Singleton equations fix a variable to zero;
two-term equations identify variables with a rational scale. Weighted
union-find handles transitive relations and cycles. The row producer applies
the coordinate map before allocating the large product. Passes are bounded;
all remaining equations go to the core, so stopping presolve early does not
change the solution. A materialized-matrix entry point skips this work when
too few short rows exist to justify it.

**Finite-field assembly and kernel-only reconstruction.** Sparse factors are
converted once per prime and their product is assembled directly in the
finite field. Only the nullspace is reconstructed. Prime denominators,
rank changes, changed free-column patterns and vanishing coefficient support
are handled explicitly. CRT uses the union of supports. This removes the
need to retain a complete rational product just in case another prime is needed.

**Exact certification without another enormous product.** For integral
factors and a primitive integral kernel, the maximum local-row 1-norm times
the maximum factor-slice 1-norm times the maximum kernel coefficient bounds
every residual. If this is smaller than the reconstruction prime, the modular
kernel congruence implies an exactly zero rational residual. Private free
coordinates prove independence and the modular rank proves completeness.
The final FEC7 and sewing bounds were only 45 bits, below the 61-bit prime.
Nonintegral or insufficient-bound cases use an exact rational residual check.
This is a deterministic certificate, not a second random-prime spot check.

**Project after restricting the ansatz.** The new projection contracts the
small sewn space into the final forward factor, applies letter maps, and
uses the preceding prefix projection. It does not build or row-reduce the
projection of all 6,826 final forward vectors. The operation is a reassociation
of exact contractions, valid for arbitrary compatible tensors.

**Small recursive symmetry actions.** Private kernel coordinates recursively
identify word evaluations that form a dual chart. They determine action
matrices on the small complete sewn space without constructing every large
intermediate action matrix. Closure of the local conditions and both endpoint
spaces is checked first. All generator-invariant combinations are then solved
and checked exactly. This does not replace the existing general irrep frame
layer; the new kernel also plugs into that layer through `symrep extend`.

**Two dependency memory fixes.** COO-to-CSR conversion previously shrank the
index allocation but retained its old capacity; a subsequent copy could read
out of bounds. The new patch synchronizes capacities and handles empty
allocations. Separately, integer/rational move assignment leaked the old
destination allocation for large values; the cleanup was backported from
upstream. Address, undefined-behavior and leak sanitizer tests pass. All seven
dependency patches were replayed onto clean pinned sources and compared with
the working dependency.

## General matrix tests

Five fresh sequential processes per method, rotating order; eight workers.
These are synthetic matrix families, not a representative sample of every
application. Solver time excludes generation and exact validation. RSS
includes those phases and a retained validation input, equally for both paths.

| Family | Existing solve | Structural dispatch | Existing / new peak RSS |
| --- | ---: | ---: | ---: |
| Chain, 12,000 variables | 94.78 ms | 1.19 ms | 14.53 / 12.50 MiB |
| Graph relations, 12,000 variables | 55.35 ms | 1.79 ms | 18.57 / 12.50 MiB |
| Propagating forced zeros | 109.56 ms | 0.77 ms | 14.53 / 12.50 MiB |
| Sparse rational coefficients | 6.01 ms | 6.41 ms | 12.66 / 12.66 MiB |
| Wide, large-nullity system | 4.04 ms | 3.78 ms | 12.66 / 12.66 MiB |
| Dense tail | 14.74 ms | 14.16 ms | 21.64 / 21.72 MiB |
| Local block relations | 5.12 ms | 5.35 ms | 12.66 / 12.66 MiB |

The measured large gains are structural, not universal arithmetic speedups.
The small mixed differences on other families are why the old backend remains
available and presolve uses a structural dispatch. The first prototype paid
unnecessary overhead when only a handful of variables could be removed; that
measurement led to the dispatch rather than a polygon-specific exception.

## Four-loop result and independent validation

The final sewn space has dimension four. Cyclic, flip and parity, generating
`D7 × C2`, reduce it to three. Collinear constraints then give a unique solution.
In this saved invariant basis the coefficients are `(36, 1/16, 1/48)`; these
numbers refer only to the accompanying `SEW_7p1.wxf`, not to another run's basis.
The remainder has no divergent letters. All nine expanded E/R/boundary tensors
at loops two through four are byte-identical between the two different kernel
and projection routes.

The four-loop word coefficients `[1^7,8]`, `[1^7,9]`, `[1^7,12]` and `[1,8^7]`
are respectively `960, 960, -960, 120`, matching the published formulas; the
corresponding lower-loop checks also pass. These were validation only, never
inputs to the solve. The published three-dimensional invariant ansatz also
agrees with the computed dimension.
[He, Jiang, Li and Liu, Eq. (12), following paragraph and Table 1](https://arxiv.org/html/2511.09669v2).

Additional checks: 159 general exact kernel cases against independent exact
ranks/residuals, including unlucky primes, prime denominators, changing support
and multiprime reconstruction; 80 random forward/backward factored operators;
90 numerical-system tests; all ten recorded lower-loop reference tensors;
recursive word charts versus explicit tensor actions; restricted projection
versus an explicit five-factor rational sum; and the symmetry CLI suite,
including streamed recursion, continuation and generator checks.

On the actual prepared data, independent recursive reference comparisons
passed through heptagon FEC5 and LEC4, and hexagon FEC8 and LEC5. The optional
streamed **factorized** hexagon LEC6 solve completed, but its full recursive
change-of-basis comparison exceeded a separate 180-second validation timeout
after passing weight five. That weight-six reference comparison is not claimed
as passed. It does not affect the completed four-loop heptagon validation,
and is another reason to retain the existing fast hexagon multiplicity default.

## SOTA alternatives and remaining opportunities

The local trial used upstream SparseRREF v0.4.2 at commit
`8ca2c9764992677c5940816f4bdeff1320d11280`. Its changes include pivot bookkeeping,
column weights, Schur work buffers and back substitution. Its published
comparisons are useful leads, but the complete local chain did not improve,
so it was not promoted wholesale.
[Official implementation](https://github.com/munuxi/SparseRREF).

SpaSM combines sparse triangular solves, greedy pivot selection, dense modular
linear algebra and structural decomposition. A sparse-to-dense Schur switch
is a credible next experiment once a remainder becomes dense. Smaller primes
could also reduce storage, at the cost of more reconstruction passes. Neither
SpaSM nor a FFLAS-FFPACK hybrid was integrated or benchmarked in this pass.
[SpaSM](https://github.com/cbouilla/spasm),
[FFLAS-FFPACK](https://linbox-team.github.io/fflas-ffpack/).

Structured Gaussian elimination in CADO-NFS supports the broader direction
of shrinking the sparse problem before its expensive linear algebra. The
implemented rational two-term presolve is much simpler than that solver's
parallel filtering algorithm; their published speedups are not our results.
[Bouillaguet and Zimmermann](https://perso.lip6.fr/Charles.Bouillaguet/static/publis/merge.pdf).

Block Wiedemann or block Lanczos could avoid elimination fill for very large
systems with a small desired kernel. Full kernel recovery, output sparsity
and exact rational certification still matter; they are not automatically
better for recurrences with thousands of output basis vectors.
[LinBox](https://github.com/linbox-team/linbox),
[CADO-NFS](https://cado-nfs.gitlabpages.inria.fr/).

Compile-time finite-field specialization and modern exact matrix-product
kernels are additional engineering options, especially for dense panels.
SparseRREF already uses Shoup arithmetic in sparse updates, so a speedup must
be measured over that implementation. No GPU or specialized-prime speedup is
claimed here.
[FiniteFieldSolve](https://arxiv.org/abs/2311.01671),
[exact GEMM work](https://arxiv.org/abs/2601.07508).

The strongest next experiments are certified equation-subset selection for
small-nullity sewing, releasing/reusing Schur storage more aggressively, and
a measured sparse/dense crossover. There is no evidence that further gains
are impossible. The present result establishes a useful general strategy
and a completed laptop-scale four-loop acceptance test, not a universal optimum.

## Reproduction and evidence

See [CLI and algorithm documentation](../../docs/solver-strategy.md). Raw logs,
fresh output directories and the compatible upstream trial are under
`output_solver_strategy_20261002/`. Compact timing records, validation results,
source/binary hashes, patch replay and general benchmark samples are copied
next to this report. The heptagon acceptance command is recorded verbatim in
`amplitude-all-streamed.json`. Use `bench/measure_process.py` to measure the
process tree and `/usr/bin/time -v` for high-water process RSS.

The defaults remain unchanged. The existing split-rational multiplicity
backend remains the preferred measured path for the hexagon; the streamed
carrier kernel is an additional general backend, not a reason to replace
that faster representation-specific recurrence.
