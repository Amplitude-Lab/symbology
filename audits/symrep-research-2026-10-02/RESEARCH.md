# Further exact-solver optimization: research and local evidence

Research checked on 2026-10-02. This is a targeted survey of relevant methods,
including recent work, not a claim to have benchmarked every available solver.
Only the experiments explicitly marked as local were run here. The production
changes and complete-chain measurements are in [REPORT.md](REPORT.md).

The useful question is which method computes our **complete rational kernel,
preserves recursive symmetry information, and reduces total time and peak RAM**.
A fast modular rank calculation alone does not answer that question. Our
heptagon LEC weight-4 system has 30,828 unknowns, 67,463 nonzero equation rows,
3,964,499 input nonzeros, and a 3,745-dimensional kernel with 1,774,777 stored
carrier coefficients. The existing modular solve usually needs only one prime
on this example. These local facts determine the priorities below.

## Candidates and decisions

| Approach | What it offers | Evidence and decision here |
| --- | --- | --- |
| Incremental validation, shared views and shorter allocation lifetimes | Less simultaneous storage and less allocation work | Implemented; saves substantial RAM while preserving exact output readback. |
| Parallel rational reconstruction and integer fast path | Lower cost around the modular elimination | Implemented with the existing height certificate and CRT fallback. |
| Representation condensation and sparse multiplicity coordinates | Remove repeated equations and unknowns before elimination | Best mathematical direction for a larger heptagon gain; needs sparse recursive coupling maps over its endomorphism field. Existing general multiplicity backend is the starting point. |
| Sparse elimination with a dense Schur complement | Keep sparse early pivots, accelerate only the dense remainder | SpaSM/FFLAS-FFPACK are serious candidates; not integrated or benchmarked here. |
| Matching, block triangularization, COLAMD/Markowitz-type ordering | Reduce fill and separate independent work | Worth testing with full-chain basis sparsity included; no claim of a gain yet. Earlier simple column-weight tests were inconclusive. |
| Smaller primes and dense FLINT RREF | Cheaper arithmetic or better dense kernels | Tested locally. Small-prime orbit gains were variable; dense RREF was slower on the tested sparse eigenspace. Neither adopted. |
| Compile-time prime specialization | Eliminate runtime modular divisions | Relevant published technique, but SparseRREF already uses Shoup arithmetic in sparse update loops. A new gain needs measurement, and compilation must be amortized. |
| Block Wiedemann / block Lanczos / matrix-free contraction | Avoid elimination fill for much larger systems | Attractive when memory is dominated by fill; full basis recovery and rational certification must be benchmarked. Not a demonstrated win on these chains. |
| GPU and exact finite-field GEMM | Accelerate sufficiently large dense panels or block products | Recent exact methods exist; no GPU implementation tested here. Whole-matrix densification is a poor starting point. |
| p-adic lifting and multi-prime parallelism | Amortize factorization or parallelize reconstruction for large coefficients | Lower priority because the measured heptagon solve succeeds after one prime. |
| Fraction-free sparse LU / exact QR | Avoid modular reconstruction entirely | Useful alternatives under different coefficient growth; no evidence of an advantage here. |
| Allocator substitution or arena limits | Reduce fragmentation and retained pages | Tested; tradeoffs vary by chain. No new allocator default adopted. |

## The symmetry direction that matches the user's objective

The current hexagon backend already solves scalar multiplicity equations and
generates the component rows through fixed intertwiners. The heptagon default
still solves a full sparse QQ kernel and then certifies an irreducible orbit
recipe. Its equation tensor can be symmetry-adapted, but the default fast path
does not yet use multiplicity blocks for that elimination.

Condensation restricts a module and its maps to the image of an idempotent,
where smaller homomorphism spaces can be processed without constructing the
full tensor space. Lux, Neunhöffer and Noeske explicitly identify basis choice
and avoiding large intermediate maps as practical issues in tensor/Hom-space
condensation. Their setting is modular representation theory; adapting that
idea to our QQ recurrence is an engineering and mathematical proposal, not a
performance result established by their paper for this application.
[Primary paper, 2012](https://www.math.rwth-aachen.de/homes/neunhoef/Publications/pdf/homcond_final.pdf).

Our prepared heptagon models contain four one-dimensional rational irreps and
two six-dimensional irreps whose endomorphism algebras have rational dimension
three and center degree three. Thus, in a six-dimensional sector, a primitive
corner retains **three rational coordinates per copy**. Equivalently, one can
work with one coordinate over a cubic field. One cannot treat that coordinate
as an ordinary rational scalar or promise a sixfold reduction of rational
unknowns. The corner reduces 6m rational coordinates to 3m for m copies;
field arithmetic and changed sparsity still cost work.

The concrete next mathematical experiment should therefore be:

1. Reuse the certified models and reduced condition coefficients already saved
   during preparation.
2. Store commutant coefficients as field elements and solve the corresponding
   multiplicity maps directly, with sparse coupling maps that remain compact
   at the next weight. Keep the general division-algebra fallback.
3. Avoid constructing the six component rows during ordinary recursion.
   Retain an exact recipe for expansion and certify the field operations and
   descent to QQ.
4. Compare the complete chain, including template construction, field
   arithmetic, output and next-weight fill. Reject a smaller block layout if
   its basis is much denser or its conversion cost exceeds the savings.

This extends the existing general multiplicity backend; it is not accomplished
by merely replacing the small-group decomposition routine. GAP's MeatAxe
interface offers homogeneous components, induced actions and homomorphism
bases. It can help obtain or cross-check models, but that does not by itself
make the growing recurrence sparse.
[GAP documentation](https://gap-system.github.io/gap/doc/ref/chap69.html).

## Sparse/dense hybrids and structural ordering

SpaSM combines sparse modular elimination, sparse triangular solves, greedy
pivot selection and dense finite-field linear algebra. It also provides
Dulmage–Mendelsohn decomposition and kernel-basis tools, using odd 32-bit
primes. It is a more relevant alternative than changing every matrix to dense
storage. Integration would still need rational reconstruction, permutation
handling, exact certification and a comparison of the resulting recursive
basis. Its full PLUQ mode has additional time and memory costs according to
the project documentation.
[SpaSM](https://github.com/cbouilla/spasm).

FFLAS-FFPACK provides finite-field BLAS-style products and triangular solves,
as well as echelon forms and nullspace bases. A promising use is the dense
remainder after sparse pivots, with switching governed by measured density
and memory. This library was researched, not installed or benchmarked here.
[FFLAS-FFPACK](https://linbox-team.github.io/fflas-ffpack/).

The local dense experiment used **FLINT**, not SpaSM or FFLAS-FFPACK. On the
heptagon LEC weight-3 cyclic eigenspace equation (734 × 734, 8,058 nonzeros,
rank 630), sparse RREF took 0.0083 s versus dense FLINT's 0.0599 s for a prime
just above 2^60. For a prime just above 2^29 the respective times were 0.0085
and 0.0377 s. These are individual diagnostics, not full solver benchmarks.
They discourage unconditional densification; they do not rule out a hybrid.
The whole nonzero-row weight-4 constraint matrix would occupy about 15.5 GiB
as an array of 64-bit residues before any workspaces.

For ordering, the relevant target is both elimination fill and the sparsity of
the output basis used at the next weight. A lower RREF time at one weight can
lose overall by making the following contraction dense. The preceding audit's
simple column weights changed hexagon kernel nonzeros by less than 0.5%; this
was not a test of a complete matching/DM/COLAMD implementation.

## Arithmetic specialization and recent accelerator work

FiniteFieldSolve is directly relevant to high-energy bootstrap systems. Its
published approach specializes arithmetic by recompiling for a fixed prime to
reduce modular-division overhead. Its speedups concern its own comparisons;
they are not measured gains over this SparseRREF implementation.
[Mangan, 2023 preprint / 2024 publication](https://arxiv.org/abs/2311.01671).

Here, `SparseRREF/sparse_vec.h` already uses `n_mulmod_precomp_shoup` and
`n_mulmod_shoup` for sparse scalar updates. That limits how much additional
benefit can be inferred from replacing `% p` with constant-prime arithmetic.
Precompiled specializations would be preferable to compiling during a
subsecond chain. No prime-specialized production backend was added.

Berthomieu, Graillat, Lesnoff and Mary describe multiword decompositions for
**exact** finite-field matrix multiplication using floating-point GEMM,
including CPU/GPU measurements. The 2026 preprint uses controlled decompositions
and bounds; it does not justify approximate floating-point nullspace solving.
It is a candidate for future dense panels or batched modular products, not a
drop-in replacement for today's sparse QQ kernel.
[2026 preprint](https://arxiv.org/abs/2601.07508).

The local orbit-selector test used primes just above 2^60 and 2^29, both
congruent to 1 modulo the group order. Weight-4 orbit times were
0.391/0.340 s for the larger prime and 0.363/0.333 s for the smaller one.
These two observations per choice show a possible modest gain, but do not
establish a stable full-chain benefit. Smaller primes also affect bad-prime
frequency and rational reconstruction requirements. The default stays unchanged.

## Black-box, lifting and fraction-free alternatives

LinBox exposes elimination and black-box methods including Wiedemann,
BlockWiedemann, Lanczos and BlockLanczos; it also has a p-adic rational solver
using block Wiedemann. Such methods avoid materializing elimination fill and
are a natural option for larger future systems.
[Method documentation](https://linalg.org/linbox-html/group__solutions.html),
[rational block-Wiedemann solver](https://linalg.org/linbox-html/class_lin_box_1_1_rational_solver_3_01_ring_00_01_field_00_01_random_prime_00_01_method_1_1_block_wiedemann_01_4.html).

For this application the operator could apply the recursive contraction and
its transpose without assembling every equation row. The comparison must
include recovering all 3,745 kernel generators, reconstructing QQ values,
organizing the representation, and supporting continuation. A test returning
one kernel vector or a rank would be an incomplete comparison. This is a
proposal, not a tested implementation.

A 2026 ETH master's thesis by Pascal Skipness presents a hybrid CPU/GPU
block-Wiedemann rank method and reports a very large rank computation. It
supports investigating this family at much larger scales, but is a thesis
about rank/lower-bound algorithms, not evidence that our full rational basis
can be produced faster.
[Institutional record](https://www.research-collection.ethz.ch/entities/publication/25c65e49-5f06-4121-966f-7b626ab9b060).

SPEX implements sparse integer-preserving LU and bundles AMD/COLAMD ordering.
This is an exact alternative when coefficient growth and reconstruction costs
justify it. Our measured one-prime reconstruction and rank-deficient
rectangular kernel task provide no current reason to expect it to win; it
would need a suitable nullspace wrapper and measurements.
[SPEX and its primary citation](https://github.com/clouren/SPEX).

## Engineering experiments and priority

Limiting glibc arenas to two lowered the old heptagon LEC median peak from
333.3 to 291.6 MiB in three diagnostic runs, while time rose from 1.470 to
1.486 s. Hexagon LEC time rose from 0.174 to 0.194 s. Loading mimalloc did not
give consistent improvements. These trials used the pre-change binary and
are retained in `allocator-trials.json`; no allocator setting was imposed.

The adopted changes address measured costs without changing the field,
equations, pivots or saved basis. For the next major algorithm experiment,
prioritize sparse condensation over the cubic endomorphism fields. In parallel
with that line of investigation, an optional SpaSM-style sparse/dense backend
is worth a full rational-basis comparison. Reserve block-Wiedemann/GPU work
for sizes where fill and matrix storage actually dominate. These priorities
are inferences from this workload, not claims of a universally fastest method.
