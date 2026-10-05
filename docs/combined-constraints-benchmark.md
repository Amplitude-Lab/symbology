# Joint integrability and symmetry: a reusable bootstrap benchmark

The solver strategy is general. The current implementation can combine exact
linear operators, discover a manageable candidate space, and impose the
remaining equations in that space. Neither the number of letters, the symmetry
group, the sewing cut, nor an expected answer dimension is built into this
algorithm. The five-loop heptagon supplied the demanding application that
motivated it. Its choice of input bases, collinear boundary, and physical
normalization remain specific to that calculation.

This benchmark makes the distinction testable. It compares three solving
orders on identical saved inputs and checks their complete rational solution
spaces independently. The [measured audit](../audits/combined-constraints-2026-10-04/REPORT.md)
contains the results; the reusable runner is
[`bench/combined_constraints_benchmark.py`](../bench/combined_constraints_benchmark.py).

## What system are we solving?

Choose already integrable recursive bases
\(F_i\in\mathrm{FEC}_{n-m}\) and \(L_j\in\mathrm{LEC}_m\). Write

\[
 S=\sum_{i,j} C_{ij}\,F_i\otimes L_j.
\]

Here the product concatenates the two symbol segments. It is not a shuffle
product. The unknowns are the entries of \(C\). Both factors already obey their
internal adjacent-letter conditions; sewing imposes the remaining conditions
at their interface. These give a homogeneous linear operator \(M\).

The saved recurrences have the form
\(F_i=\sum_{a,\alpha}f_{ia\alpha}F'_a\otimes\ell_\alpha\) and
\(L_j=\sum_{\beta,b}l_{j\beta b}\ell_\beta\otimes L'_b\).
If \(D_{q\alpha\beta}\) encodes the supplied adjacent-letter conditions, then

\[
 (MC)_{aqb}
 =\sum_{i,j,\alpha,\beta}
    f_{ia\alpha}D_{q\alpha\beta}l_{j\beta b}C_{ij}=0.
\]

The older indices \(a,b\) stay as recursive basis indices. No expansion of
\(\mathrm{LEC}_m\) into all length-\(m\) words is required. The input condition
tensor can include adjacency restrictions in addition to integrability; the
solver enforces the tensor actually supplied.

For each generator \(g\), use the row-action convention

\[
 gF_i=\sum_a(A_g)_{ia}F_a,\qquad
 gL_j=\sum_b(B_g)_{jb}L_b.
\]

Invariance of the product is exactly

\[
 A_g^T C B_g-C=0.
\]

If \(x\) is the row vector obtained by flattening \(C\), with the right index
varying fastest, this is \(x(A_g\otimes B_g-I)=0\). Consequently the desired
space is

\[
 \ker\begin{pmatrix}
 M\\
 (A_{g_1}\otimes B_{g_1}-I)^T\\
 (A_{g_2}\otimes B_{g_2}-I)^T\\
 \vdots
 \end{pmatrix}.
\]

This display specifies the equations; the implementation retains their factors
instead of assembling this entire matrix. Generator constraints imply
invariance under the group they generate. For the heptagon benchmark, the
supplied generators are cyclic, flip, and parity, giving \(D_7\times C_2\).

We do **not** first solve the complete integrability system and subsequently
filter its kernel in the joint methods. Nor do we first construct the entire
invariant product space. Both types of equations participate from presolve
onward. In particular, nontrivial representations on the two factors can pair
to produce an invariant: retaining only the individually invariant factors
would lose solutions.

## How the staged method solves it

1. **Propagate exact short relations.** An equation \(a x_i=0\) eliminates a
   coordinate; \(a x_i+b x_j=0\) identifies two coordinates up to a rational
   scale. Both integrability and symmetry contribute these relations. This
   reduces the variables before the expensive elimination.
2. **Discover rank over a finite field.** Generate equations from their
   factors, beginning with small supports. Use bounded input batches and sparse
   elimination; select further support bands as needed. Expensive discovery
   bands may be sampled. These selections change the work order only.
3. **Switch to candidate coordinates when useful.** Suppose the rows of
   \(K\) span the current candidate space, so \(x=yK\). For a remaining
   equation block \(E\), solve
   \[
      R=EK^T,\qquad T=\operatorname{rowker}(R),\qquad K\leftarrow TK.
   \]
   This retains exactly the candidates satisfying the new block. Previously
   imposed equations remain satisfied because each new row is a linear
   combination of old rows. Factored contractions apply \(E\) directly to
   the small candidate without materializing \(E\).
4. **Complete every operator.** Residual intersection includes all original
   constraints, without sampling. An expected physical dimension is never an
   acceptance condition. If discovery has already made a complete unrestricted
   pass, the solver reuses that row space rather than constructing a large
   intermediate kernel solely to annihilate it again.
5. **Reconstruct and certify over the rationals.** Reconstruct rational kernel
   rows from modular calculations. Acceptance requires independent rows, the
   complete modular rank bound, and exact annihilation of every operator.
   Residuals are checked modulo enough primes that their combined modulus
   exceeds a rigorous bound on the absolute integer residuals after clearing
   denominators. A divisible integer smaller than this modulus is zero. Extra
   certification primes need not repeat the expensive rank calculation.

The identity in step 3 works for any linear system. Its efficiency depends on
candidate dimension, factor sparsity, elimination fill, and rational coefficient
size. A large or dense intermediate kernel can be costly. Input-batch limits
do not bound retained elimination fill or the output basis; the benchmark also
uses external process memory and time guards.

## What is standardized in the repository?

| Layer | Reusable implementation | Responsibility |
| --- | --- | --- |
| Operator composition and staged solving | [`staged_kernel.hpp`](../staged_kernel.hpp) | Discovery, intersections, complete checks, certification |
| Explicit rational matrices | [`staged_matrix.hpp`](../staged_matrix.hpp) | Homogeneous and affine matrix adapters |
| Recursive extension and arbitrary sewing cuts | [`tensor_kernel.hpp`](../tensor_kernel.hpp) | Factored local conditions and residual contractions |
| Supplied product generators | [`product_invariance.hpp`](../product_invariance.hpp) | Factored \(A_g^TCB_g-C\) equations |
| Exact short relations | [`structured_kernel.hpp`](../structured_kernel.hpp) | Coordinate elimination, proportionality, and lifting |
| Rational reconstruction | [`certified_kernel.hpp`](../certified_kernel.hpp) | Kernel reconstruction and acceptance certificates |

The shared C++ entry point is:

```cpp
auto local = tensor_kernel::right_boundary_conditions(
    condition, LEC, tensor_kernel::BoundaryReduction::raw);
tensor_kernel::ConstraintRows integrability(FEC, std::move(local), LEC.dim(0));
staged_kernel::ProductInvarianceRows symmetry(left_actions, right_actions);
auto equations = staged_kernel::stack(symmetry, integrability);
auto K = staged_kernel::solve(equations, opt); // rows are solution coefficients
```

The stack borrows its operators, which must remain alive during the solve.
Each operator supplies discovery rows, complete residual blocks, and a valid
integer residual bound. The public command supports any positive sewing cut:

```bash
./bootstrap --sew -c CONDITION.wxf -f FEC_LEFT.wxf -l LEC_RIGHT.wxf \
  -o INVARIANT_KERNEL.wxf --sew-strategy staged --threads 2 \
  --local-reduction raw \
  --left-action F_cyclic.wxf --right-action L_cyclic.wxf \
  --left-action F_flip.wxf   --right-action L_flip.wxf \
  --left-action F_parity.wxf --right-action L_parity.wxf
```

Actions must refer to the **actual saved bases**, in matching generator order.
This layer accepts supplied actions; it does not infer the group, prove all
group relations, or discover the correct action matrices from their dimensions.
The benchmark's preparation helper checks closure of the seed and condition
spaces and derives factor actions recursively. Preparing an input is separate
from solving the invariant kernel.

This is also distinct from producing the **entire kernel classified by all
irreducible representations**. The joint invariant solve selects the trivial
representation of the product. It works with ordinary or adapted factor bases.
The existing `symrep` layer handles representation frames and recursive actions;
its compact multiplicity backend currently requires scalar rational
endomorphism algebras. The full rational heptagon decomposition does not satisfy
that restriction. These benchmarks do not establish completion of that larger
recursive irrep solver.

See [the general library document](staged-kernel.md) for adapter contracts,
affine systems, subgroup restrictions, and policy options. Existing defaults
for large recursive spaces are unchanged by this benchmark.

## Benchmark contract and reproduction

Every case is solved by three methods:

| Method | Computation |
| --- | --- |
| `integrability-first` | Compute the complete integrability kernel, then impose generator invariance in its coordinates. |
| `joint-global` | Combine both families from presolve onward, generate all reduced modular equations, and use global elimination. |
| `joint-staged` | Combine both families and use adaptive discovery followed by complete residual intersections. |

All three use the same current SparseRREF dependency and shared optimizations.
The first method is an optimized baseline, not a frozen historical release.
The second comparison isolates the effect of staging after imposing symmetry
early. Generator restriction in the first method uses factor contractions too;
it is not penalized by constructing an unnecessary full Kronecker matrix.

On Linux, with the repository's build dependencies available:

```bash
make bench/combined_constraints_benchmark
python3 bench/combined_constraints_benchmark.py \
  --output output_combined_constraints_trial --threads 2 --repeats 3
# Equivalent default suite:
make bench-combined-constraints BENCH_OUTPUT=output_combined_constraints_trial2
```

The committed [default case file](../bench/combined_constraints_cases.json)
contains two three-letter \(S_3\) examples at cuts `1+5` and `3+3`, plus heptagon
`2+2` and `4+2`. The second toy uses rational changes of coordinates so its
generator actions are not limited to permutations or sign flips. All examples
solve the full invariant integrable product space, before additional amplitude
boundary conditions.

Use `--cases PATH.json` to add a different alphabet, group, or cut. Paths in
the file are relative to its directory. A scalar-terminal example is:

```json
[
  {
    "name": "my-example-4p3", "kind": "bootstrap",
    "left_weight": 4, "right_weight": 3,
    "condition": "data/condition.wxf",
    "first": "data/FEC_1.wxf", "last": "data/LEC_1.wxf",
    "generators": ["data/generator1.wxf", "data/generator2.wxf"]
  }
]
```

The helper builds the factors and their actions once, before timing. For
already prepared factors, including non-scalar terminal representations, use
`{"name":"my-bundle","kind":"bundle","directory":"prepared"}`. The sealed
directory contains `condition.wxf`, `left.wxf`, `right.wxf`, paired square
`left_g0.wxf` / `right_g0.wxf` matrices and so on, `input.tsv` containing
`combined-constraints-v1 NUMBER_OF_GENERATORS`, and the seal generated by
`symrep::seal`. The caller is responsible for correct recursive provenance,
internally integrable factors and actions in those coordinates. The generic
solve supports these inputs; the scalar seed preparation helper does not build
an NMHV terminal representation for the caller.

For the additional unpacked-hexagon case in the recorded run, we used its
`dlogmat_full_ca_pd3.wxf`, cyclic and flip matrices, and weight-one seeds
exported by `bench/hexagon_symrep_reference` from the package's MHV chains.
In those coordinates parity is the cube of cyclic, so the two generators
already impose it. The exact paths and input hashes are recorded in the audit.
This compares solving orders on the package's data; it does **not** time the
package's specialized irrep solver as a fourth method.

Each run directory freezes its executable, solver source files, raw inputs,
prepared inputs, commands, and hashes. `replay_cases.json` reuses its prepared
inputs without requiring the original archive paths. To repeat the recorded
five-case suite with the current executable:

```bash
python3 bench/combined_constraints_benchmark.py \
  --cases output_combined_constraints_20261004/replay_cases.json \
  --output output_combined_constraints_replay --repeats 3
```

Use `--program output_combined_constraints_20261004/snapshot/bin/combined_constraints_benchmark`
to run the frozen binary instead. Its shared-library dependencies must remain
available. Source snapshots and binary hashes are evidence of the measured
version; this is not a self-contained compiler/runtime container.

Methods run in fresh sequential processes, with their order rotated between
repeats. The default guards are 600 seconds, 6 GiB sampled process-tree RSS,
8 GiB address space, and a 4 GiB host-available-memory floor **per process**.
The suite's total duration can be longer. Set `--local-reduction reduced` for
a separate experiment that reduces the contracted interface before solving;
the default benchmark uses `raw` consistently for all methods.

`REPORT.md` reports median internal solve time, process wall time, and GNU
time's peak process RSS. Solve time includes interface construction, presolve,
elimination and required certificates. Process wall additionally includes
loading, input-seal checks, normalization, output and exact serialization
readback. The monitor polls every 100 ms, so its short process times have that
granularity. Millisecond toy timings are correctness controls, not reliable
performance rankings. Sampled tree RSS is also retained; neither RSS statistic
is total host memory or a measure of page-cache use.

Preparation and independent verification have separate measurements. After
the timed runs, the checker uses the original rational sewing assembler and
explicit product actions to form a reference. It compares canonical rational
row spaces for every output, establishing completeness as well as annihilation.
The reference uses the same exact arithmetic dependency but separate assembly
and symmetry application. If it cannot complete within its guard, that case
is not reported as verified. A failed method is never counted as a speed win.

## What the measurements establish

The recorded suite passed all 45 solves: five cases, three methods, three
repeats. On heptagon `4+2`, median solve time fell from **2.132 s to 1.063 s**
and peak process RAM from **154.40 MiB to 54.33 MiB**, comparing
`integrability-first` with `joint-staged`. Against `joint-global`, staging saved
about 13% of solve time and 45% of peak RAM in that case. The `2+2` and small
hexagon cases favor joint global elimination over staging. This supports a
general strategy with measured policy choices, not a universal speed claim.

There is a useful subtlety in the heptagon `4+2` example: the complete
integrability kernel is already two-dimensional and invariant. Early symmetry
still reduces the presolved coordinate count from 8,967 to 2,008 and exposes
cheaper elimination. The speedup there comes from useful redundant equations
and staged work reduction, rather than discarding a large final non-invariant
space. For the hexagon example, symmetry reduces the complete kernel from
41 dimensions to 9. The audit records both mechanisms explicitly.

The completed five-loop heptagon `8+2` calculation remains a separate
[large-scale reference](../audits/fiveloop-day-2026-10-04/REPORT.md): 7,796.38 s
and 19.40 GiB sampled peak RSS for solving and physical verification from
cached prerequisites. It used ordinary FEC8 coordinates with supplied symmetry
actions, and obtained a four-dimensional joint kernel before physical boundary
conditions selected unique coefficients. It did not construct all unrestricted
weight-ten integrable symbols first. That historical result predates this
benchmark and the shared solver refactor; it is neither a fresh run nor a
three-method timing comparison. Building the cached FEC8 had separately peaked
at 32.67 GiB. The later parity-adapted FEC8 is a different coordinate choice
for the known space, not an input required by the successful five-loop run.
