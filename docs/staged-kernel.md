# Staged exact kernels and arbitrary sewing cuts

The successful five-loop MHV calculation combines early FEC/LEC sewing,
product invariance, sparse discovery, complete kernel intersection, and
rational certification. The elimination and intersection machinery is now
shared by ordinary matrices, tensor recurrences and supplied symmetry
generators. `early_sew.hpp` is a compatibility adapter rather than a second
solver. SparseRREF remains the elimination and arithmetic dependency.

The [combined-constraints benchmark](combined-constraints-benchmark.md) derives
the sewing and product-invariance equations and compares three solving orders
with independent exact-space checks and separate time/peak-memory measurements.

## The identity behind staged solving

Let the rows of `K` be a basis for the solutions already constructed, so
every candidate has the form `x = y K`. For the next equations `B x^T = 0`,
solve the smaller system

```text
R = B K^T
T = a row basis of ker(R)
K_new = T K
```

This is an exact parametrization of the intersection. Every earlier equation
remains satisfied because every new row is a combination of old rows. It
applies to any exact linear system, with or without symmetry. Splitting
equations into stages does not itself guarantee faster computation: both the
candidate dimension and the cost of applying `B` matter. A large dense `K`
can be more expensive than ordinary elimination.

For rational systems, discovery uses finite fields. Exact short-row presolve
removes zeros and identifies proportional variables. Cheap equations establish
sparse pivots before larger equations are selected. Once the candidate is
small, every original operator is applied in complete blocks, preferably by
contracting its factors into `K` rather than assembling `B`. No expected
physical nullity is supplied. `target_nullity` is only a discovery threshold.

Discovery also observes completed elimination batches. Four consecutive batches
adding at most two pivots can trigger intersection early, provided at most 256
candidate coordinates remain. This avoids scanning increasingly expensive rows
for very little rank gain. The three thresholds are configurable and affect
work ordering only: the complete operator still determines the final kernel.

If discovery has already completed an unsampled, unrestricted pass, its
triangular basis spans every input equation and the redundant residual pass
is skipped. Otherwise all omitted rows are included in residual intersection.
Masks and sampling never apply to certification. Rational reconstruction is
accepted only with a complete modular rank bound, independent rational kernel
coordinates, and exact residual certification. Extra primes can enlarge a
rigorous residual-height certificate without repeating discovery.

## Library interfaces

For a materialized rational matrix with equations as rows:

```cpp
#include "staged_matrix.hpp"

rref_option_t opt;
opt->pool.reset(8);
staged_kernel::Options options;
options.target_nullity = 128;
auto K = staged_kernel::solve(A, opt, options);  // rows span ker(A)
auto affine = staged_kernel::solve_affine(A, b, opt, options);
// If consistent: x = affine.particular + y * affine.directions.
// If inconsistent: particular has zero rows; directions still span ker(A).
```

`solve_affine` takes one rational right-hand-side vector. It homogenizes the
system with one extra coordinate, computes the complete kernel, and intersects
with that coordinate equal to one. It distinguishes an inconsistent system
from a consistent system with free variables. This API is separate from the
existing collinear workflow; it does not change that workflow's RHS semantics.

Compose factored operators without a combined equation matrix:

```cpp
#include "tensor_kernel.hpp"
#include "product_invariance.hpp"

tensor_kernel::ConstraintRows local(FEC, local_conditions, LEC.dim(0));
staged_kernel::ProductInvarianceRows symmetry(left_actions, right_actions);
auto operators = staged_kernel::stack(symmetry, local);
auto K = staged_kernel::solve(operators, opt, options);
```

The stack and matrix adapters are non-owning; their inputs must remain alive
through the solve. Each source supplies `ncol`, `nrow`, exact short-row
generation, modular discovery, complete residual blocks and a rigorous integer
residual bound. See `MatrixRows` for the smallest complete adapter and
`ConstraintRows` for contraction of an implicit tensor operator. The height
bound is part of the correctness contract, not a performance estimate.
Sources use canonical sparse rows with valid indices and nonzero coefficients.

`column_restriction.hpp` wraps a replayable source on an arbitrary selection
of coordinates, retaining its factors and its complete residual certificate:

```cpp
#include "column_restriction.hpp"
staged_kernel::ColumnRestriction selected(local, selected_columns);
auto K = staged_kernel::solve(staged_kernel::stack(selected), opt, options);
selected.lift_in_place(K);  // restore original columns; omitted coordinates are zero
```

Selection order is respected, and exact short relations are transferred to
the original coordinates before generating rows. The adapter does not assert
that independent selections exhaust the full kernel. That requires a proved
block decomposition. The experimental `bench/split_character_extension`
checks a prepared one-dimensional rational character decomposition, solves
every character sequentially, and writes an ordinary adapted recurrence chain.
Each sector gets its own exact certificate and a hashed checkpoint. This also
allows using a subgroup to partition computation while retaining the entire
solution space; it does not label the result by the larger group's irreps.
Measure equation density, coefficient heights and recursive output size before
choosing a subgroup. More sectors do not necessarily mean lower cost.
The experimental driver accepts `front_nonzeros=N` to bound the initial
global pivot search separately from its support threshold. A front that splits
a large set of short rows can change the final coordinate chart, coefficient
heights and output size. Thus a smaller front need not minimize total memory.
`partitioned_tensor_view.hpp` lets the ordinary WXF writer concatenate saved
row partitions while retaining only one partition in memory. The experimental
driver releases each solved sector and then the factored operator before
writing the full recurrence. This reduces overlap of large allocations; the
largest single sector still needs to fit in memory.

The reconstruction backend recognizes an already reduced identity block by
requiring a private unit column for every row. This proves independence and
allows extracting the modular kernel without another SparseRREF pivot pass.
It falls back to ordinary elimination when that certificate is absent. This
avoids needless coordinate changes; it does not guarantee small coefficients
when the preceding elimination has already chosen a costly chart.

The public numeric API is over Q. Discovery/intersection internally operates
over good finite fields. This is not a floating-point rank solver and makes
no claims about conditioning or approximate numerical nullspaces.

## FEC_(n-m) × LEC_m

Any cut with positive left and right weight can use the same sewing solver.
The backward basis at weight `m` is a rank-three recurrence indexed by
`(basis, first letter, previous right basis)`. Its third index is retained as
an equation index; it is **not** expanded into all words. Internal equations
have already been enforced on both sides. Sewing imposes the interface
conditions between their adjacent letters.

```bash
./bootstrap --sew -c CONDITION.wxf -f FEC_n_minus_m.wxf -l LEC_m.wxf \
  -o SEW_n_minus_m_pm.wxf --sew-strategy staged --threads 8
```

The contracted interface itself may contain dependent equations. By default,
`--local-reduction reduced` first computes its rational row basis. Use
`--local-reduction raw` to pass the normalized, nonzero contracted rows directly
to the global solver instead. This removes a potentially expensive preliminary
rational elimination; dependencies are then handled together with the other
constraints. Both choices impose the same complete interface. Neither expands
`LEC_m` into words.

For the heptagon weight-six `3+3` cut, a pilot fell from 154.8 s / 740 MiB to
4.5 s / 288 MiB by choosing `raw`. For the weight-eight `6+2` cut, `raw` was
slower and used more memory. Keep the choice explicit; do not infer it just
from the cut weight. In C++, pass `BoundaryReduction::raw` to `sew_one` or
`right_boundary_conditions`. The CLI option belongs to `bootstrap --sew`;
the existing `compute_rhs` physical driver does not expose arbitrary cuts.

Optional symmetry acts on the product, not separately on invariant subspaces
of its two factors. Nontrivial irreps on either side can pair to make an
invariant. Supply matrices in the **actual saved FEC and LEC bases**, paired
in the same generator order:

```bash
./bootstrap --sew -c CONDITION.wxf -f FEC_LEFT.wxf -l LEC_RIGHT.wxf \
  -o SEW_INVARIANT.wxf --sew-strategy staged --threads 8 \
  --left-action F_cyclic.wxf --right-action L_cyclic.wxf \
  --left-action F_flip.wxf   --right-action L_flip.wxf
```

There is no hard-coded group or alphabet. Each generator imposes
`x (F_generator tensor L_generator) = x`. Matrix dimensions are checked
against the corresponding factor, and mismatched generator counts fail.
This builds an invariant product kernel; it is distinct from constructing
every irrep sector. The existing `symrep` representation and recursive-action
layers retain that responsibility.

An optimal cut need not be balanced in weight. Compare:

* the cost of constructing or loading both basis chains;
* `dim(FEC_(n-m)) * dim(LEC_m)` unknown coefficients;
* nonzeros after right-boundary contraction and exact short-row presolve;
* sparse elimination fill, generator-action density, and residual contraction;
* the size of the final kernel and subsequent projection/recursive actions.

The existing `ConstraintRows::profile` reports structural equation costs
without materializing their product. Use small certified pilots when estimates
are close; basis dimensions alone cannot predict fill. Choosing the smallest
coefficient product alone is not an automatic optimal-cut algorithm. The
public sewing interface supports arbitrary cuts, while the experimental MHV
physical-verification driver still uses its explicit two-letter right cut.

## Policy and resource controls

`--kernel-strategy staged` is available for forward/backward `bootstrap`
extension, its `compute_rhs` child commands, and the `symrep` factorized
backend. `--sew-strategy staged` is explicit for sewing. Existing defaults
remain `streamed` for extensions and `auto` for sewing: high-nullity recursive
spaces can favor global sparse elimination, so a universal staged default is
not justified by a small physical kernel.

`Options` controls the cheap initial support bound, complete cheap bands,
discovery nullity, equation masks, sampling, and batch sizes. The initial
front now also has a nonzero budget, after which additional seed equations
are streamed. No equations are truncated. These bounds control pending input,
not elimination fill, operator storage, or the required output basis. Use an
external process memory/deadline guard for an actual workstation RAM budget.

The adaptive controls are `stagnation_batches`, `stagnation_rank_gain` and
`max_residual_nullity`; setting `stagnation_batches = 0` disables this early
switch. The ordinary `target_nullity` stop still applies. These are general
policy knobs, with no stored polygon-specific dimensions or expected answers.

Two general implementation improvements avoid unnecessary work: unchanged
triangular rows retain their sorted allocation during short-relation
propagation, and an unchanged rank skips repeated propagation entirely.
Finite support-limited discovery also rejects long rows using indices before
allocating or multiplying coefficients. Complete residual checks are unchanged.

The next backend alternatives should be judged against the same operator and
certificate interface. LinBox documents sparse elimination and black-box
Wiedemann methods for exact systems; the latter can reduce storage but requires
many operator applications. It is a candidate for fill-dominated problems,
not an implemented replacement here. See the
[LinBox method overview](https://linalg.org/linbox-html/group__solutions.html).
Small dense residual systems could instead use
[FLINT's modular nullspace](https://flintlib.org/doc/nmod_mat.html), but a new
backend should be adopted only after measuring the conversion and reconstruction
costs as well as elimination time.

## Validation and measurements

`make check-kernel` includes exact generic matrix, stacked-operator and affine
tests, empty and inconsistent systems, difficult primes, a deliberately tiny
front budget and all five cuts of a weight-six commuting-alphabet example.
The cuts are compared after independent expansion to original words.
An adversarial plateau fixture places an essential equation after many
dependent rows, so an early switch must recover it through the full residual.
`make check-early-sew` includes nontrivial S3 pairings and false-candidate
rejection. NMHV and `symrep` regressions exercise staged CLI paths too.

`bench/staged_kernel_benchmark assembled|staged small-kernel|large-kernel N`
provides a deterministic non-polygon benchmark with exact-space checks. Measure
fresh processes and both time and peak RAM. Separate claims about reduced work
or storage from timing differences under unrelated CPU load. The completed
[five-loop audit](../audits/fiveloop-day-2026-10-04/REPORT.md) remains evidence
for the pre-refactor implementation, not a timing measurement of this library.
The [staged-library audit](../audits/staged-kernel-2026-10-04/REPORT.md) records
the new timings, peak RSS, arbitrary-cut comparisons, four-loop physical
verification, and the NMHV cases where staging loses to the existing solver.
