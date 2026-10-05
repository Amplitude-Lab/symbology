# Shared staged kernel and arbitrary-cut audit — 2026-10-04

The early-sewing implementation is now a shared exact-kernel library. It accepts
ordinary rational matrices, affine systems, factored tensor equations, and
paired symmetry generators. The new adaptive policy improves four-loop heptagon
sewing, and skipping an expensive preliminary interface reduction makes a
deeper right cut practical. Staging does not win on every problem: existing
recursive-extension defaults remain unchanged.

## General strategy

If rows of `K` span the current solutions, write `x = y K`. For a further block
of equations `B`, compute `R = B K^T`, solve for rows `T` spanning its kernel,
and replace `K` by `T K`. All previously imposed conditions remain satisfied.
This identity is independent of the alphabet, group, number of loops, and
whether the original operator is stored explicitly.

The library uses the identity only when a sufficiently small candidate space
has been found. Its complete policy is:

1. Eliminate exact one- and two-term relations before modular elimination.
2. Establish a sparse initial front with a configurable nonzero budget; stream
   further rows once that budget is reached.
3. Generate cheap support bands before sampled expensive bands. Observe rank
   gain after completed batches. A configurable stagnation rule can switch to
   residual intersection before reaching the ordinary nullity target.
4. Apply every original operator to the candidate in complete blocks, without
   sampling. Factor contractions avoid assembling the full equation matrix.
5. If an unrestricted, unsampled discovery pass already processed everything,
   reuse its row basis directly. Do not build and re-annihilate a large kernel.
6. Reconstruct only the kernel over Q and accept it only with a complete modular
   rank bound and exact residual-height certification. Extra residual primes
   can certify a candidate without repeating discovery.

The implementation also avoids scanning unchanged elimination rows, sorting
unchanged supports, copying discovery masks, and multiplying coefficients for
long rows that support-limited discovery will discard. These improvements do
not depend on heptagon dimensions. Pending-input budgets do not bound retained
fill or the required output basis; a hard workstation RAM limit still needs a
process guard.

Symmetry is an optional operator in this design. Product invariance imposes
`x (F_g tensor L_g) = x` for arbitrary supplied generators and keeps invariant
pairings of nontrivial irreps. It does not replace the existing `symrep` layer
that organizes the entire kernel by irreps and computes recursive actions.

## Library and repository integration

* `staged_kernel.hpp`: non-owning operator stacks, bounded discovery, adaptive
  intersection, complete-pass fallback and certification.
* `staged_matrix.hpp`: ordinary Q-matrix adapter and `solve_affine(A,b)`, returning
  consistency, a particular solution, and the homogeneous directions.
* `product_invariance.hpp`: factored product-generator operator.
* `tensor_kernel.hpp`: extension and arbitrary rank-three right-recurrence
  sewing adapters; raw or reduced contracted interface equations.
* `early_sew.hpp`: compatibility adapter to the shared implementation.

`bootstrap` exposes `--kernel-strategy staged`, `--sew-strategy staged`, paired
`--left-action` / `--right-action` matrices, and `--local-reduction reduced|raw`.
`compute_rhs` propagates the staged strategy to child commands; its existing
physical workflow still uses a one-letter right cut. The experimental physical
MHV verifier still uses a two-letter right cut. The general library and native
sewing command accept deeper cuts. The `symrep` factorized backend also accepts
the staged strategy and retains its representation frames and action checks.

For `FEC_(n-m) × LEC_m`, the backward recurrence retains its previous-basis index
as an equation index. No full right-word expansion is needed. The best cut is
selected by both basis-chain costs, product dimension, interface sparsity,
elimination fill, action density and output size. Balanced weight alone is not
a useful guarantee.

The contracted interface can itself be expensive to reduce over Q. `raw`
normalizes and removes zero rows but defers its dependencies to the global
solver. `reduced` computes its row basis first. Both impose the complete same
conditions. Raw interfaces can save substantial work, but can also increase
global fill and peak RAM; the default remains `reduced`.

## Time and peak memory

Fresh processes, eight workers, GCC 14.4.0, `-O3 -march=native -mtune=native`,
the existing patched SparseRREF/FLINT/GMP/TBB stack, no mimalloc. The host has
about 62 GiB RAM and an unrelated CPU-heavy process was left running throughout.
Methods were alternated; small differences are not reliable speedups under
this contention. Peak memory below is GNU time's peak process RSS, in MiB.
The process guard also recorded sampled tree RSS, with a 24 GiB RSS limit,
28 GiB address-space limit, 300 s trial deadline, and 4 GiB free-memory floor.
No reported trial hit those limits.

| Problem / policy | Wall seconds | Peak MiB | Repetitions |
| --- | ---: | ---: | ---: |
| Four-loop heptagon `6+2`, saved pre-refactor solver | 18.982 | 511.54 | 3, median |
| Four-loop heptagon `6+2`, shared adaptive solver | 15.590 | 520.74 | 3, median |
| Generic small kernel, assembled | 0.404 | 143.50 | 3, median |
| Generic small kernel, staged | 0.404 | 125.31 | 3, median |
| Generic large kernel, assembled | 0.404 | 136.17 | 3, median |
| Generic large kernel, staged | 0.404 | 136.28 | 3, median |
| Hexagon NMHV backward chain through weight 5, streamed | 1.460 | 85.50 | 2, median |
| Hexagon NMHV backward chain through weight 5, staged | 4.740 | 186.69 | 2, median |
| Heptagon NMHV backward chain through weight 3, streamed | 0.706 | 89.05 | 2, median |
| Heptagon NMHV backward chain through weight 3, staged | 0.956 | 86.98 | 2, median |

The four-loop sewing measurement starts from cached FEC and generator actions;
it includes right-basis construction, interface setup, kernel solving and
certification, but excludes building the forward chain and physical amplitude
verification. Its runtime improves by 17.9%; median RAM is 1.8% higher, not
lower. In the first paired run, discovery processes 991,896 equations instead
of 1,288,533 and 20,639,544 input nonzeros instead of 24,968,680. The adaptive
switch occurs at candidate dimension 155; complete intersections reduce it
through 96, 93, 22 and finally 3. This explains the reduced work without assuming
the physical answer's dimension.

The generic benchmark has 512 columns and about 12,000 dependent long rows;
its independently known nullities are 4 and 432. Staging reduces small-kernel
peak RSS by 12.7%, with no measurable end-to-end speedup. Its internal solve is
also not faster in these trials; input generation occupies much of the short
process duration. The large-kernel fallback avoids a major regression, but
does not establish a general speed advantage.

The NMHV results are negative evidence against making staged solving the
universal default. They measure recursive basis construction, not a complete
NMHV amplitude. Full exact recursive-space comparisons pass for every trial.
The initial extraction alone also showed no four-loop speedup; those trials
are retained in `evidence.json` as `initial_measurements`.

### Arbitrary cuts and interface reduction

| Heptagon cut | Interface | Wall seconds | Peak MiB | Kernel dimension |
| --- | --- | ---: | ---: | ---: |
| Weight 6, `4+2` | reduced | 1.311 | 42.56 | 2 |
| Weight 6, `4+2` | raw | 0.905 | 57.34 | 2 |
| Weight 6, `3+3` | reduced | 154.786 | 740.48 | 2 |
| Weight 6, `3+3` | raw | 4.537 | 288.35 | 2 |
| Weight 8, `6+2` | raw | 18.098 | 723.45 | 3 |

These are single pilots, not repeated medians. The `3+3` raw policy is about
34 times faster and uses 61% less peak RSS. Reduced-interface setup consumed
135.7 s of its 154.8 s total; reducing the small interface first was the main
bottleneck, and also changed downstream sparsity. The `6+2` raw pilot is worse
than the repeated adaptive reduced-interface result above. Neither policy is
universally preferred. Both cut kernels agree exactly after expanding only
the changed recursive interface, and every raw kernel agrees exactly with
its reduced counterpart.

## Correctness and limitations

Passed:

* 303 exact generic matrix, stacked-operator and affine cases, including
  inconsistent systems, bad primes, tiny front budgets, large kernels, and all
  five cuts of a weight-six commuting-alphabet example. An adversarial plateau
  test puts an essential condition after many dependent rows and checks that
  complete residual intersection recovers it.
* 102 rational factored extension/sewing cases against independent exact
  equations, including 20 additional staged/raw multi-terminal comparisons;
  160 structured-kernel cases; 16 product-invariance intersections with
  nontrivial S3 pairings and extra-prime/false-candidate checks.
* Native C3 generator CLI tests, staged NMHV recursive-space tests, and staged
  `symrep` forward/backward action verification.
* All adaptive/baseline and raw/reduced physical-kernel space comparisons.
* Public native regression: 90 numerical checks, condition export/reimport,
  all ten historical two-/three-loop CRCs, and all six default E/R/boundary
  tensors agree with legacy physical outputs.
* Independent four-loop amplitude verification using the frozen prior physical
  driver and the new adaptive kernel: unique coefficients, all 533,360
  divergent-word conditions, and published coefficients 960, 960, -960, 120.
  The complete 985,265-term collinear `E4_candidate.wxf` is byte-identical to
  the earlier ordinary-workflow `4loop/E4.wxf` (SHA-256
  `2fd9e43469ea16f5d975883b58e83499df0676c094823cb1e265acb8c6f36e6f`).

`make regression` was attempted and stopped before calculation because its
historical `output/collinear/` seed is absent. Its private front-end project
fixtures are also absent. This gate is **not reported as passed**. The public
native suite reconstructs and passes the two-/three-loop reference outputs
without those old directories; it does not substitute for the unavailable
private-project baselines.

At the time of this audit, the library had not been rerun through the complete
five-loop pipeline. The later [publication validation](../../benchmarks/heptagon-mhv-five-loop/VALIDATION.md)
completed that cold reproduction in 8,280.83 s / 31.48 GiB, including all
prerequisites, and recovered the same final symbol byte for byte. Its joint
kernel took 6,743.00 s / 19.33 GiB and passed complete exact certification.
The prior 2 h 10 min / 19.40 GiB five-loop result remains documented in the
[separate audit](../fiveloop-day-2026-10-04/REPORT.md), and is not a measurement
of this refactor. No claim is made that a particular cut or staging policy
minimizes all inputs, or that the current library implements black-box
Wiedemann or a new dense finite-field backend.

## Reproduction and evidence

See [the API and strategy guide](../../docs/staged-kernel.md). Standard build
and tests use `make check-kernel check-early-sew check-staged check-nmhv
check-symrep check-public`; supply the normal compiler/dependency flags for
the machine. The benchmark drivers are `bench/staged_kernel_benchmark.cpp`,
`bench/nmhv_kernel_benchmark.cpp` and `bench/early_sew_probe.cpp`.

`evidence.json` retains all trial commands, exit status, process measurements,
test results and log excerpts, source/input/executable hashes, and hashes of
the independently compared kernels. Local full logs and generated tensors are
under `output_staged_library_20261004/`. Baseline executables and the changed
pre-refactor headers are frozen in its `baseline/` directory; the initial
nonadaptive library executables are in `initial-library/`. Only comments and
whitespace were adjusted in calculation sources after the final build.
