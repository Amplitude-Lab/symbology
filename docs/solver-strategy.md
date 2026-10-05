# Exact kernels, staged intersection and restricted projection

The shared staged solver and arbitrary `FEC_(n-m) × LEC_m` sewing policy are
documented in [staged-kernel.md](staged-kernel.md). It exposes ordinary rational
matrices, affine systems, composed factored operators and optional product
symmetry through one intersection engine. Select `staged` explicitly for
small-kernel problems; retain global sparse elimination for large recursive
kernels unless measurements favor staging.
The shared policy now bounds the initial front, switches from discovery on
measured rank stagnation, and permits skipping a costly rational interface
reduction. The [staged-library audit](../audits/staged-kernel-2026-10-04/REPORT.md)
compares these choices in time, peak RAM and exact solution space.

The [certified five-loop MHV result](../audits/fiveloop-day-2026-10-04/REPORT.md)
completed in 2 h 10 min at 19.40 GiB for the solve and verification from cached
prerequisites; FEC8 construction previously peaked at 32.67 GiB. Those are
pre-refactor measurements and are not a cold end-to-end benchmark.
The [five-loop lessons and six-loop direction](#five-loop-lessons-and-six-loop-direction)
below are the retained project guidance for future large calculations.

The new backend solves homogeneous rational equations while retaining their
sparse factors. It works independently of a symmetry group. The symmetry
solver can use it for its carrier kernel and retain the existing irreducible
frame and recursive generator recipes.

For a four-loop heptagon MHV **symbol**, starting from the supplied seeds:

```bash
./compute_rhs --target SEW_7p1 --data-dir data --output-dir output_fourloop \
  --letter-projection output_fourloop/collinear/colprojdiv_w1.wxf --threads 8
```

Use a fresh output directory for a timing measurement. The projection file
in this command is generated automatically from the supplied data. The
separate `divergent` support-filter interface imposes all conditions with any
divergent letter; its cost is not the same measurement as the command above.

The resulting `4loop/hepMHV_4L_recursive.wxf` is a tensor with dimensions
`1 × dim(FEC7) × dim(LEC1)`. Together with `FEC_1` through `FEC_7` and `LEC_1`,
it specifies the full original-alphabet symbol. It is not an expanded list of
words. `E4.wxf` and `R4.wxf` are expanded **collinear-limit** tensors, with the
11-letter collinear alphabet. `solMHV_4L_original.wxf` gives coefficients in
the saved symmetry-invariant `SEW_7p1.wxf` basis. Coefficients depend on the
chosen basis and must not be transferred between runs with different kernels.

These options are separate so each part can be compared:

| Option | Effect | Scope |
| --- | --- | --- |
| `--kernel-strategy streamed` | Exact presolve, finite-field assembly and kernel-only reconstruction | Forward and backward extension; `bootstrap`, `compute_rhs`, and factorized `symrep extend` |
| `--sew-strategy streamed` | Contract the right recurrence before equation assembly and use the streamed kernel | `bootstrap` and `compute_rhs`; supports multiple terminal components and deeper right factors |
| `--projection-strategy restricted` | Contract the small sewn space before projecting its final forward factor | `compute_rhs`; one-step right boundary |
| `--threads N` | Set the worker count, including child bootstrap commands | `bootstrap`, `compute_rhs` |
| `--require-unique` | Fail if collinear constraints leave free coefficients | Single-target `bootstrap --solve-collinear`; automatically required by `compute_rhs` |

Following the [NMHV acceptance audit](../audits/nmhv-kernel-2026-10-02/REPORT.md),
`bootstrap` and `compute_rhs` default to `--kernel-strategy streamed`.
`--sew-strategy auto` contracts the side with fewer basis coefficients first:
it uses the streamed right contraction when the right basis is no larger,
and the existing left contraction otherwise. Small heptagon NMHV sewing was
faster with the latter. Explicit `original` and `streamed` choices remain.
The MHV-only `compute_rhs` workflow defaults to `--projection-strategy restricted`.
The factorized `symrep extend` backend also defaults to the streamed kernel;
the fast scalar multiplicity backend for split rational symmetry models is
retained as its own default. `make` now builds `symrep` too.

The dependency remains pinned SparseRREF `5bbee55` (v0.3.6) with all seven
validated local patches, installed by `scripts/setup-sparserref.sh`. Upstream
v0.4.2 is not silently adopted: its previous measured deployment cost was worse.
For historical basis-dependent coefficients and CRCs, select all three
`original` strategies. `make regression` does so explicitly; the public native
regression additionally checks the new defaults against the same physical
E/R tensors. Never transfer coefficients between different basis files.
For a streamed symmetry chain, use:

```bash
./symrep extend --input PREPARED --max-weight 5 --output CHAIN \
  --backend factorized --kernel-strategy streamed --threads 8
```

Preparation still uses the supplied generators, condition tensor and seed;
there is no polygon or group-specific kernel implementation. `--compare`
checks the streamed kernel against the independently assembled equations and
validates the recursive actions. It materializes the ordinary matrix, so use
it for correctness tests rather than memory measurements on large systems.

## Algorithm and certificate

For each homogeneous equation, a single term fixes a variable to zero. Two
terms identify two variables up to an exact rational scale. A weighted
union-find structure stores these relations and detects cycles forcing zero.
The remaining operator is assembled directly modulo a good prime, applying
those relations while generating rows. No complete rational equation matrix
is retained alongside it.

SparseRREF still performs the finite-field elimination. Only the nullspace
is reconstructed over the rationals. Bad denominators and rank-changing
primes are rejected. CRT combines the union of coefficient supports, and
changing free-column patterns restart reconstruction. The current retry
limit is 16 primes; failure is explicit, never an uncertified result.

CRT and rational reconstruction are fused by row. Once a row reconstructs,
its rational coefficients encode the same residues, so the large integer
copy can be released. A later prime decodes that row modulo the accumulated
modulus before continuing. This works even when different rows require
different numbers of primes or a final residual certificate rejects a
candidate. The shared CRT inverse is computed once per prime.

Factored modular equations are generated in parallel tiles and emitted in
their original deterministic order. File-producing extensions write the
kernel directly, without constructing another rational tensor. WXF writing
and exact byte readback use bounded buffers; this avoids allocating the
whole output file while a large kernel remains in memory.

Completeness follows from a modular rank lower bound and an independent
set of rational null vectors of the complementary dimension. For an integral
factored operator and a primitive integral candidate, a rigorous bound on
every residual is

```text
max_local_row_l1 × max_factor_slice_l1 × max_abs_kernel_entry.
```

If this bound is smaller than the accumulated CRT modulus, the already-established
modular kernel congruences prove the rational residual is exactly zero.
Rational local equations and complete older-component slices are first
scaled to clear denominators, without changing their kernels. It avoids
computing a large matrix/kernel product again. Otherwise the backend checks
the residual with exact rational arithmetic. This is not verification at a
random additional prime. Random rational/integer tensor tests independently
check the final residual, rank and nullspace independence with exact arithmetic.

The factored source represents
`sum_(b,l,j) F[b,o,l] C[q,l,j] X[b,j] = 0`; the same code handles either
extension direction and sewing with an arbitrary rank-three right recurrence. Its terminal/older index
is retained as an additional equation index. The generic
`structured_kernel::solve` also accepts a materialized sparse matrix, using
a cheap structural dispatch to skip presolve when too few short rows exist.

## Symmetry and projection

The MHV workflow checks closure of the local conditions and both endpoint
spaces under every available cyclic, flip and parity generator. Recursive
private word coordinates then determine action matrices on the small complete
sewn space. The invariant map is checked exactly against every generator
before collinear solving. A nonunique solve now fails instead of being
reported as an amplitude.

Restricted projection changes contraction order. For a small ansatz
`S[n,f,l]`, it contracts the `f` index into the final forward tensor first,
then applies the letter maps and the preceding prefix projection. It avoids
constructing and reducing a projection for every vector of the large final
forward space. Exact tests compare this operation with the explicit
five-factor sum. The original projection route remains available for comparison.

## Five-loop lessons and six-loop direction

The project direction is to preserve these general improvements and assess a
six-loop MHV heptagon symbol (weight twelve) on a supercomputer. Five loops
demonstrated that reformulating the problem can matter more than provisioning
RAM for the original matrix. It supports investigating six loops; it does not
yet establish a six-loop runtime, RAM requirement, or successful solve.

The completed five-loop run used `FEC_8 × LEC_2`, cyclic + flip + parity
(`D7 × C2`), and a recursive output. From cached prerequisites, solve and
verification took 7,796.38 s with 19.40 GiB sampled peak RSS. The separately
measured FEC8 construction peaked at 32.67 GiB. Keep both numbers when planning
resources: reporting only the small final kernel omits the expensive basis
and elimination stages. The accepted result, source/input hashes and complete
validation are recorded in the
[five-loop audit](../audits/fiveloop-day-2026-10-04/REPORT.md).

Preserve the following lessons when extending the solver:

1. **Change the space before enlarging the machine.** Early two-sided sewing
   avoided unrestricted FEC9, whose counted equation payload alone was
   82.55 GiB. For higher weight, compare `FEC_(n-m) × LEC_m` cuts before building
   additional unrestricted basis levels. Keep both sides recursive.
2. **Impose product symmetry early.** Retain nontrivial left/right sectors that
   pair into invariants. Five-loop success used supplied generator equations;
   it is not evidence that a fully irreducible basis alone caused the gain.
3. **Discover cheaply, then intersect completely.** Five-loop discovery streamed
   8,262,891 equations and 6,813,101,616 input nonzeros, retaining 211,386,644
   triangular nonzeros. The small candidate shrank from 128 to 5 to 4 under
   complete factored residuals. Streamed input volume is not resident memory;
   retained fill and candidate density are the quantities to monitor.
4. **Separate reconstruction from certificate enlargement.** Reconstruction
   succeeded at one 61-bit prime. A second good prime checked the small
   candidate to exceed the 76-bit residual bound, without repeating the large
   elimination. Sampling and discovery projections never replace complete
   mathematical or physical verification.
5. **Budget the entire pipeline.** Fuse CRT/reconstruction by row, release
   temporary copies, serialize with bounded buffers, and use compact projection
   and packed shuffle accumulation. Earlier attempts failed after certification
   during output: solving the equations is not the last memory-intensive step.
   Fit physical coefficients with small projections, then verify the complete
   collinear conditions. Preserve hashes, immutable inputs, stage results and
   an external shared time/RAM guard.
6. **Measure downstream effects and keep alternatives.** A smaller transient
   matrix can produce a denser recursive basis. Staging loses on some NMHV
   spaces, and raw versus reduced interface equations favor different cuts.
   The subsequent adaptive library improves the four-loop control and now
   passes the complete [cold five-loop reproduction](../benchmarks/heptagon-mhv-five-loop/VALIDATION.md):
   2 h 18 min including all prerequisites, with 31.48 GiB peak RSS and a
   byte-identical final symbol. The joint kernel took 6,743 s / 19.33 GiB.
   These separate measurements are not a controlled speedup comparison;
   do not extrapolate them to six loops.

For a six-loop preflight, compare at least `10+2`, `9+3` and `8+4`. The latter
can reuse FEC8 but requires building and measuring a larger right basis and
its actions; it is not automatically the best cut. Include prerequisite
construction, product columns, local-interface reduction, modular fill,
residual contractions, kernel density, boundary/projection costs and disk
usage. Use bounded pilots of the same exact operator before requesting a
production allocation; do not extrapolate linearly from loop order or the
four-dimensional five-loop kernel.

The present kernel uses shared-memory workers. A high-memory compute node is
the direct deployment route; additional nodes do not automatically combine
their RAM for this implementation. Independent cut pilots, modular-prime work
where needed, and verification blocks are potential cluster tasks, subject to
consistent coordinate/rank handling and measured I/O costs. Distributed
elimination and resumable in-progress elimination checkpoints would require
additional implementation. At five loops one reconstruction prime sufficed,
so parallel primes alone would not have removed the dominant elimination cost.

The physical validation drivers currently support loops two through five and
an explicit two-letter right cut. Before accepting six loops, generalize and
test the physical projection/boundary path for the selected cut, prepare the
required lower-loop E/R inputs, establish uniqueness rather than assume it,
and verify every original constraint plus the complete physical boundary.
Do not extend the existing published-word formula beyond its validated loop
range without independent justification. The deliverable remains a certified
recursive symbol with its basis chain, not a claim about function-level
constants or a fully expanded original-alphabet word list.

## Validation

An experimental higher-weight MHV route imposes invariance on the complete
product FEC × LEC before sewing, then discovers a small kernel and checks it
against every original factored equation. It uses deeper right factors to
avoid building the next unrestricted FEC space. The original implementation and
measured controls are documented in the
[early-sewing audit](../audits/early-sewing-2026-10-03/REPORT.md).
`make check-early-sew` checks the product-invariance solver and the compact
exact shuffle accumulator used for its collinear validation. This route is
available through the shared staged library and optional native sewing flags.
It does not change the public CLI defaults described above.

`make check-nmhv` checks the promoted CLI defaults on public heptagon NMHV
fixtures, including exact comparisons with the original solver and the
published pre-collinear dimensions 5 and 11. The external hexagon benchmark
and reproduction instructions are in the NMHV audit.

`make check-kernel` runs general rational kernel cases, random factored
operators, recursive word-chart comparisons, a restricted-projection test,
and tensor capacity/copy tests. `make check-symrep` also tests the streamed
factorized recursion. `python3 tests/native_regression.py` checks the numerical
solvers and recorded lower-loop physical outputs.

The separate `bench/heptagon_amplitude_check.cpp` validates published special
word coefficients after solving. These values are never supplied to the solver.
See the [strategy audit](../audits/solver-strategy-2026-10-02/REPORT.md) for
hardware, costs, limitations, alternative algorithms and raw evidence.

The [five-loop feasibility audit](../audits/fiveloop-2026-10-03/REPORT.md)
records larger-weight memory measurements and a matrix-count preflight.
An optional experimental batch elimination backend processes all equations
while retaining a triangular row basis. It is not the default: the initial
FEC7 trial produced a denser recursive basis and ran slower, despite a small
memory saving. Reducing temporary memory is insufficient if it makes the
next recursive problem larger.
