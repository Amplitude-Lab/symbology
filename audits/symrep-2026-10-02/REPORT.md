# Why the symmetry solver was slow, and what changed

The large regression was an implementation problem, not a requirement of
symmetry-adapted solving. The first implementation solved smaller blocks but
repeated expensive component contractions, lattice reduction and rank checks
around them. The archive avoids those operations by keeping a compact
multiplicity recurrence in a fixed sparse parity × D3 convention.

The optimized solver now follows that recurrence strategy with automatically
constructed rational models. It still accepts general finite matrix groups;
there is no D3-specific dispatch or ordinary-solver fallback.

## Controlled comparison

Median wall seconds, three runs, two SparseRREF workers. Preparation is a
separate reusable operation; chain times include group/model loading, template
construction, solving, mandatory checks and WXF output. Exact reference
comparison is outside the timed solve.

| Hexagon chain | Final dimension | Original native | Previous general | Optimized general | Unpacked archive |
| --- | ---: | ---: | ---: | ---: | ---: |
| FEC through weight 8 | 340 | 0.120926 | 0.350927 | 0.103262 | 0.037130 |
| MHV LEC through weight 5 | 696 | 0.213885 | 1.069730 | 0.093120 | 0.046655 |
| MHV LEC through weight 6 | 1,858 | 2.953775 | stopped at ~10 GiB | 0.564094 | 0.189850 |

The optimized FEC chain is close to the original; LEC is about 2.3× faster at
weight 5 and 5.2× faster at weight 6. Compared with the previous general solver,
FEC weight 8 improves about 3.4× and LEC weight 5 about 11.5×. These are measured
cases, not a promise of a speedup for every sparse matrix.

| Final-weight phase (seconds) | Old assembly | New assembly | Old kernel + reconstruction | New kernel + reconstruction |
| --- | ---: | ---: | ---: | ---: |
| FEC 8 | 0.038980 | 0.005883 | 0.163242 | 0.021457 |
| MHV LEC 5 | 0.044145 | 0.004327 | 0.821036 | 0.020889 |

Peak RSS for optimized hexagon FEC 8 / LEC 5 / LEC 6 is about
24 / 26 / 123 MiB. The corresponding native peaks are 28 / 55 / 428 MiB.

Reusable preparation medians (seconds): hexagon FEC 0.010568, hexagon LEC 0.013223, heptagon FEC 0.252109, heptagon LEC 0.246444.
Every newly prepared bundle is byte-identical in its checksummed content to
the earlier prepared input used by the benchmark.

The workstation is an Intel Core Ultra 7 265 running Linux. Native and general
solvers use GCC 14.4.0, FLINT 3.2.2, `-O2 -std=c++20`, no mimalloc override,
the same patched SparseRREF and two workers. Processes run sequentially in
fresh directories; caches are not flushed and the shared workstation is not
CPU-isolated. Archive binaries retain their own build and bundled SparseRREF.
The old-general executable SHA-256 exactly matches the earlier audited binary.

## Root causes

### 1. Shorter integer vectors were confused with cheaper sparse recurrences

The old `kernel_blocks` called lattice saturation and LLL on every sufficiently
large-coefficient corner kernel. These operations preserve the solution space,
but minimize neither nonzero count nor subsequent sparse elimination cost.
They can mix otherwise sparse free-coordinate basis rows.

On hexagon MHV LEC at weight 5, the old reduced equations had 105,784 nonzeros.
Omitting lattice reduction gives 52,519. An ablation changing only the LLL
condition reduced the complete weight-5 chain from about 1.09 s to 0.21 s.
The remaining optimizations reduce it further. Thus the slowdown cannot be
attributed simply to the cost of finding group irreps.

The old weight-6 LEC diagnostic reached approximately 10 GiB RSS and was stopped;
it is not counted as a completed timing. Disabling only LLL completed in a
median 4.14 s, showing that other costs also needed attention.

**Change:** when `dim_QQ End_G(U)=1`, use one primitive normalization per whole
copy and preserve the sparse kernel basis. Do not apply saturation/LLL.
Higher-dimensional endomorphism algebras retain the existing shortening policy.

### 2. Expanding and recontracting irrep components defeated the reduction

The old `assemble_reduced` evaluated a four-factor scalar contraction for
every sector at every weight. It repeatedly traversed expanded recurrence
components, dense adapted condition entries and both change-of-basis charts.
Many intermediate contributions cancelled only at the end. The heptagon
baseline spent roughly 82 s assembling FEC weight 4 and 57 s assembling LEC
weight 3, far more than its actual reduced linear solves.

The archive's `make_template` contracts its fixed CG maps and dlog tensor once;
`assemble` then multiplies the resulting routes by sparse multiplicity
coefficients. Its production recurrence does not expand every representation
component at each step.

**Change:** `RecurrenceTemplates` constructs a general rational basis for
`Hom_G(U_b, U_a tensor alphabet)` (reversed for LEC). For each local basis
intertwiner `H`, it caches the corresponding reduced equation matrix

```text
T[a,b,H,s] = left_chart[a,s] * contract(H, condition) * right_corner[b,s]^T.
```

At weight `w`, assembly only accumulates these fixed routes multiplied by
coefficients between older and previous copies. A solved copy's cyclic vector
is retained directly in these coordinates for the next step; it is not recovered
by expanding and decomposing the full tensor again. `wN_multiplicity.wxf` saves
the compact coefficients alongside compatible full exports.

The local template contraction itself uses a bounded sparse matrix contraction
followed by sparse products. Simply caching the old deeply nested scalar loop
was still too expensive for heptagon setup and was discarded.

This is not restricted to scalar multiplicities: every rational
`End_G(U_b)` coordinate is retained. Tests include degree-three commutants for
D7 and a non-scalar quaternionic coefficient for Q8.

### 3. The solver repeated rank calculations already certified by its kernel

For `End_G(U)=QQ`, independent primitive-corner kernel vectors already give
independent whole irreducible copies. The old code nevertheless projected each
orbit back into the corner, constructed coordinate charts and eliminated those
vectors again.

**Change:** lift these scalar-multiplicity copies directly. For larger
endomorphism algebras, use the original SparseRREF kernel's private free
coordinates to test independence. Restriction to those coordinates is
injective on the kernel, even after an invertible lattice basis change. Thus
rank tests use nullity columns rather than the larger ambient candidate space.

### 4. Canonical sparse products were unnecessarily rebuilt

The wrapper normalized every multiplication output through a `std::map` of
big rationals, including rows already sorted, merged and zero-free by
SparseRREF. That cost also appeared in mandatory intertwiner checks and export
reconstruction.

**Change:** detect canonical rows with a linear scan and return immediately;
retain the original merging path for unsorted, duplicated or zero entries.

## Why the archive remains faster

The archive and the generic solver solve the same spaces, as certified by
invertible rational changes of basis at every weight. Their arithmetic work
and outputs are not identical:

| Final constraint nonzeros | Native full matrix | Old general corner blocks | Optimized general corner blocks | Archive corner blocks |
| --- | ---: | ---: | ---: | ---: |
| Hexagon FEC weight 8 | 244,892 | 83,999 | 62,857 | 38,336 |
| Hexagon MHV LEC weight 5 | 433,088 | 105,784 | 52,519 | 46,653 |

The archive's fixed CG convention and primitive/sorted copy choices produce
sparser equations, especially for FEC. The generic solver also discovers and
checks models, exports full component tensors and candidate transformations,
and performs exact intertwiner and serialization checks. The archive normally
writes compact multiplicity tensors; expansion is an explicit separate export.
Its bundled SparseRREF version and build configuration also differ. These are
end-to-end implementation comparisons, not a claim that identical RREF calls
differ in speed.

The next performance work should preserve the multiplicity recurrence: persist
templates with prepared inputs, choose rational models/CG and copy bases for
sparsity as well as coefficient height, and make full component export an
explicit consumer requirement. For non-split rational irreps, further optimize
endomorphism-space arithmetic and basis shortening. Increasing thread counts
or returning to full component contraction would not address these causes.

There is no universal guarantee that adding symmetry information makes every
sparse problem faster. Tiny chains can be dominated by setup and the cost of
returning a richer answer; an adapted basis can also lose useful physical
sparsity. The meaningful goal is to eliminate avoidable work and compare
complete, equal solution spaces at representative sizes.

## Heptagon: improved, but not at native parity

| Chain | Original native (s) | Previous general (s) | Optimized general (s) | Optimized peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| FEC through weight 4 | 0.068247 | 92.093124 | 8.246404 | 561.0 |
| LEC through weight 3 | 0.064196 | 71.987206 | 15.763335 | 952.9 |

Old heptagon figures are from the retained three-run baseline audit; the native
and optimized figures were measured in this investigation. The new chains are
exactly the same spaces as the original native chains at every weight.

At FEC weight 4, assembly is 4.35 s, including 3.04 s of previously unseen
template construction; SparseRREF itself is about 0.43 s. At LEC weight 3,
template construction is 8.98 s, SparseRREF 1.61 s, lattice shortening 2.68 s,
and copy selection/lifting 1.14 s. The last operation previously took roughly
9.16 s before restricting rank tests to kernel free coordinates.

The rational D7 × C2 models have dimensions `1,1,1,1,6,6` and endomorphism
dimensions `1,1,1,1,3,3`. The six-dimensional sectors require three rational
coordinates per primitive corner, rather than the scalar multiplicities of
hexagon D3. Automatic adaptation also produces much denser conditions and
recurrences than the original physical basis. These remaining costs are real;
the current general heptagon implementation is still far slower than native.

Caching trades memory for repeated arithmetic. FEC peak RSS rises from the
old 246 MiB to about 561 MiB; LEC falls from about 2,124 MiB to 953 MiB after
the rank-test optimization. Persisting reusable templates with a bounded
resident cache is therefore a priority, alongside sparser model/CG choices.

## Verification and reproduction

The public checks exercise C3/C5/C7, D3/D5/D7, S4, C2 × C2 and Q8, dense rational
changes of input coordinates, sequential kernels, both recurrence directions,
vector terminals, empty spaces, continuation and corrupt bundle rejection.
New tests compare cached contractions with independently assembled full tensor
equations, including non-scalar endomorphism coefficients, changed recurrence
coefficients and reuse of compact solved coordinates. `extend --compare` also
checks each cached block against the previous component contraction.

The native baseline calls the unchanged `extend_forward`/`extend_backward`
functions in `bootstrap.hpp` through
[`native_symrep_reference.cpp`](../../bench/native_symrep_reference.cpp), with
the same compiler, SparseRREF checkout and two workers as `symrep`.
[`symrep_optimization.py`](../../bench/symrep_optimization.py) runs fresh,
sequential native, old-general, optimized-general and archive chains, snapshots
the optimized binary and sources, and checks exact references outside the
timed processes. The original archive and baseline artifacts are preserved.

Implementation:
[`symrep.hpp`](../../symrep.hpp),
[`symrep_bootstrap.hpp`](../../symrep_bootstrap.hpp),
[`tests/symrep_probe.cpp`](../../tests/symrep_probe.cpp), and the
[usage guide](../../docs/symrep.md).


The completed hexagon FEC 8 and LEC 5 runs each have direct exact changes of
basis against both the native and archive chains. For hexagon LEC 6 and
heptagon, the first extended reference chains were checked directly; final
runs are certified by byte-for-byte equality of **every** recurrence tensor
and copy layout to those verified chains. Their digests are retained. One
redundant repeat of the expensive LEC-6 reference comparison was interrupted
after the identical output was established; that partial directory is excluded
from completed results.

Evidence: [hexagon timings](hexagon-results.json),
[final follow-up timings and exact digests](followup-results.json),
[extended reference checks](extended-reference-results.json),
[LLL ablation](ablation-results.json),
[exact ablation identity checks](ablation-verification.json),
[preparation](preparation-results.json),
[public algebra checks](public-probe.txt),
[public CLI checks](public-regression.txt),
[continuation from a pre-optimization chain](legacy-resume.txt), and
[binary provenance](provenance.json).

All local chains, source/binary snapshots, phase logs, templates' compatible
component exports and basis-change certificates remain under
[`output_symrep_optimization_20261002/`](../../output_symrep_optimization_20261002/).
The final main benchmark is `validated-hexagon/`; final extended runs are named
`final-hexagon-*` and `final-heptagon-*`. Intermediate diagnostic runs are kept
separately and are not substituted for final measurements.
