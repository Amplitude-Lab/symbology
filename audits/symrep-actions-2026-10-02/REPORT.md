# Basis construction and recursive symmetry transformations

The best measured choice is to retain scalar multiplicity recurrences for the
hexagon, and a sparse carrier plus a certified irreducible recipe for the
heptagon. Symmetry transformations should then be constructed or applied in
the logical irreducible coordinates. Expanding the heptagon basis during
every solve remains much more expensive and provides no further reduction
in the cost of those block transformations.

This pass adds an exact transformation-export command, a blockwise application
routine, and faster recursive actions in carrier coordinates. It also reduces
heptagon LEC basis-construction peak RAM by about 5.4%. The complete solve time
is essentially unchanged; the large new speedups concern subsequent symmetry
calculations. No claim of a globally optimal algorithm is made.

## Basis construction: full process

Five fresh sequential processes per method, two requested workers, rotating
method order. Previous means the exact binary measured in the
[preceding research audit](../symrep-research-2026-10-02/REPORT.md), rerun in this
session. Times include startup, loading, recursive solving, output, exact WXF
readback and completion checksums. Prepared inputs are reused; expansion and
optional reference verification are separate.

| Complete chain | Previous symmetry, s | Current symmetry, s | Ordinary, s | Hexagon package, s |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 0.0749 | 0.0675 | 0.0639 | — |
| Heptagon LEC ≤3 | 0.0683 | 0.0700 | 0.0475 | — |
| Heptagon FEC ≤5 | 0.3036 | 0.3045 | 0.3626 | — |
| Heptagon LEC ≤4 | 1.2384 | 1.2351 | 1.2605 | — |
| Hexagon FEC ≤8 | 0.0412 | 0.0425 | 0.0999 | 0.0380 |
| Hexagon MHV LEC ≤6 | 0.1496 | 0.1497 | 2.4371 | 0.1766 |

| Complete chain | Previous symmetry, MiB | Current symmetry, MiB | Ordinary, MiB | Hexagon package, MiB |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 20.5 | 20.7 | 20.5 | — |
| Heptagon LEC ≤3 | 20.7 | 20.3 | 17.9 | — |
| Heptagon FEC ≤5 | 68.8 | 68.1 | 111.0 | — |
| Heptagon LEC ≤4 | 260.2 | 246.1 | 328.3 | — |
| Hexagon FEC ≤8 | 16.0 | 15.9 | 28.6 | 18.2 |
| Hexagon MHV LEC ≤6 | 32.3 | 32.5 | 354.1 | 37.0 |

Peak RAM is GNU time's maximum resident set for the whole process. The larger
heptagon timing differences are within ordinary run variation. This session
was faster for both versions than the preceding audit; subtracting its old
1.449 s result from this session's 1.235 s would falsely attribute machine/load
variation to this patch. The valid before/after comparison is 1.238 versus
1.235 s. Heptagon LEC RAM samples span 250.1–268.8 MiB before and 242.5–252.2
MiB after.

The current and ordinary programs were built with GCC 14.4.0, FLINT 3.2.2,
`-O3 -march=native -mtune=native -flto`, without mimalloc. The ordinary baseline
uses the original `bootstrap.hpp` recurrence and the same shared arithmetic
and I/O patches. The archive retains its own SparseRREF, pivot method and
mimalloc. No build, test or other benchmark ran concurrently with these
measurements. OS caches were not flushed and unrelated user processes were
not stopped. Raw samples are in `basis-results.json`.

## Are symmetry transformations easier in the new basis?

Yes, once coefficients are expressed in the logical irreducible basis. These
measurements start from completed saved chains. They construct **all supplied
generator matrices** at the final weight, exactly over QQ.

| Final weight | Dimension | Ordinary-basis construction | Irreducible-block construction | Ordinary / adapted matrix nonzeros |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC5 | 911 | 6.495 ms | 0.116 ms | 83,651 / 4,033 |
| Heptagon LEC4 | 3,745 | 82.057 ms | 0.730 ms | 1,910,976 / 16,585 |
| Hexagon FEC8 | 340 | 3.217 ms | 0.041 ms | 32,601 / 1,246 |
| Hexagon MHV LEC6 | 1,858 | 2,225.060 ms | 0.249 ms | 2,094,942 / 6,810 |

The ordinary-basis comparison already uses the **optimized** recursive action
algorithm from this pass, not an intentionally slow full Kronecker product
or a new unknown-coefficient solve. Ordinary and carrier bases are different
coordinate choices; the benchmark keeps them separate.

For heptagon LEC4, merely constructing the adapted matrices is about 112 times
faster than deriving the ordinary-basis matrices. This is not a 112-fold
speedup of basis construction. The irreducible layout was obtained as part
of the preceding symmetry solve.

The C++ helper `apply_irrep_action` avoids even the repeated block-matrix
allocation. Application measurements use a deterministic, mostly nonzero
rational coefficient row and apply every supplied generator:

| Final weight | Ordinary matrices | Explicit adapted matrices | Implicit irrep blocks |
| --- | ---: | ---: | ---: |
| Heptagon FEC5 | 1.171 ms | 0.103 ms | 0.060 ms |
| Heptagon LEC4 | 17.146 ms | 0.376 ms | 0.236 ms |
| Hexagon FEC8 | 0.591 ms | 0.044 ms | 0.024 ms |
| Hexagon MHV LEC6 | 80.192 ms | 0.286 ms | 0.116 ms |

These are five-process medians; each process records an eleven-sample median
after one warmup. Equal arrays of coefficients in different bases generally
represent different physical vectors. This comparison measures operator cost
for comparable dense coordinate inputs, not equality of those physical
vectors. Explicit and implicit actions **in the same adapted basis** are
checked for exact equality. Converting an existing ordinary/carrier vector
into adapted coordinates is extra work, excluded here. To obtain the benefit
repeatedly, retain coefficients in adapted coordinates between symmetry
operations rather than converting back and forth.

Across the whole transformation benchmark, peak memory was 116.7 MiB for
ordinary-basis heptagon LEC actions versus 13.3 MiB for adapted actions; for
hexagon LEC, 182.6 versus 12.3 MiB. These process peaks include input loading,
writing, exact readback and application trials. The probe's total wall time
also includes twelve application batches per weight and is therefore **not**
a construction-only timing. All observations and phase boundaries are in
`action-results.json`.

## Why the recursive action becomes simple

Let `rho_s(g)` be the fixed rational matrix of an irrep, and let `m(w,s)` be its
copy count at weight `w`. In the adapted basis,

```text
D_w(g) = direct sum_s [I_m(w,s) tensor rho_s(g)].
```

The recursion discovers the copy counts and the basis embedding, but the
small `rho_s(g)` never changes. Constructing `D_w` only repeats known blocks;
applying it directly needs a small accumulator for one copy at a time. A
product of group generators can likewise be evaluated in the small models
before applying its block action.

For a factorized heptagon basis, the exact frame identity is

```text
H_w R_carrier,w(g) = D_w(g) H_w.
```

The saved orbit recipe fixes `H_w` without requiring its dense entries during
normal solving or action export. Thus the logical irreducible basis and its
simple transformations remain available while the recursive constraint solve
continues to use sparse carrier coordinates. The default heptagon elimination
still solves a full sparse QQ kernel; this work does not claim to have turned
it into a multiplicity-only elimination.

For an explicit adapted recurrence tensor `T_w`, the corresponding check is
`T_w (D_(w-1)(g) tensor R_letter(g)) = D_w(g) T_w`, reversing tensor factors
for LEC and using the tensor's alphabet coordinates. Tests exercise these
identities directly in both recursive directions.

## Implemented changes

1. **Construct only selected tensor-action columns.** Coordinate recovery only
   needs columns selected by the carrier chart. The new implementation visits
   incoming nonzeros of those columns instead of scanning every nonzero in
   the full tensor-product action. Rows are sorted by construction. This
   applies to both QQ and finite-field action inheritance.
2. **Use private columns directly over QQ.** Kernel rows expose private
   free-variable columns. Their reciprocals recover coordinates without
   copying/transposing the whole rational basis through `Chart` or multiplying
   by an explicit diagonal inverse. A general chart remains available when
   private columns do not exist.
3. **Release modular echelon storage after extracting an inclusion.** Orbit
   selection no longer retains completed large RREF matrices while allocating
   the later pairing and sector workspaces. This is the principal additional
   full-solve memory reduction.
4. **Expose exact action export and implicit application.** `symrep actions`
   exports matrices at all requested saved weights, with generator names,
   coordinate metadata and integrity checksums. `apply_irrep_action` applies
   arbitrary group elements directly through the small models.

On the *same* heptagon LEC carrier, the old recursive rational action routine
took 197.9 ms at weight 4; the new one takes 75.7 ms, about 2.6 times faster.
Every saved action matrix is byte-identical. Peak process memory in this
comparison falls from 176.5 to 104.2 MiB. The new direct adapted action takes
0.730 ms, but acts in a different coordinate system.

Normal basis construction already inherited actions modulo a prime using a
private chart, so it did not have the same expensive QQ chart overhead. Its
weight-4 action phase changes only from 0.0659 to 0.0630 s. The QQ kernel and
orbit-selection phases still dominate. This explains why downstream exact
transformations improve much more than the complete solve.

## Explicit multiplicity-basis alternative

The general heptagon multiplicity backend was rerun with the current code,
followed by full exact reference verification at every weight:

| Complete chain | Explicit multiplicity backend | Peak RAM | Current factorized chain |
| --- | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 4.830 s | 747.4 MiB | 0.0675 s / 20.7 MiB |
| Heptagon LEC ≤3 | 10.653 s | 1,076.4 MiB | 0.0700 s / 20.3 MiB |

The explicit-backend entries are single diagnostic runs, not five-run medians.
They differ by orders of magnitude, so they suffice to reject changing the
default on these inputs. At LEC3, the corner systems retain 2,832 columns
instead of 4,956, but contain 1,364,502 nonzeros rather than the raw system's
102,551. Template construction, coefficient growth and shortening dominate.
The exact reference checks passed; the problem is efficiency, not validity.

Both representations already supply the same small irrep action blocks.
Paying for the explicit dense recurrence provides no corresponding advantage
for this transformation calculation. A compact recurrence over the heptagon's
cubic endomorphism fields, with sparse coupling maps, remains the important
unimplemented mathematical direction. Its benefit is not established by this
pass. The hexagon already benefits from its split scalar-multiplicity format.

## Validation and use

- The public ten-group exact algebra suite and CLI regressions pass.
- Private and general charts, zero-dimensional actions, and both tensor-factor
  orders agree with independent full tensor actions over QQ and modulo a prime.
- Implicit action application agrees with explicit blocks for all group
  elements in the public examples, including nonsplit and quaternionic models.
- New CLI tests verify exported matrices against actual recursive tensors or
  materialized frames, including vector terminals, empty kernels and C3.
- The ordinary heptagon and hexagon FEC/LEC transformations pass independent
  exact recursive intertwiner checks through weight 3.
- At every production weight, the new carrier action matrices match the old
  rational algorithm byte-for-byte, across all five repeats. All real-data
  action exports match the diagnostic matrices exactly.
- Every saved basis from the six main solver cases matches the preceding
  exactly verified audit, and all repetitions agree. Prepared checksums are
  checked before reusing the prior full reference proofs. This is not a claim
  that all expensive high-weight frame expansions were repeated in this pass.
- The explicit multiplicity alternatives received fresh full exact reference
  verification at every weight. Logs and `multiplicity-trials.json` retain it.

Example using the new heptagon LEC chain:

```bash
./symrep actions \
  --input output_symrep_sparse_20261002/validated/heptagon-lec-prepared \
  --chain output_symrep_actions_20261002/release/heptagon-lec-w4-symmetry-0 \
  --output heptagon_lec_irrep_actions --threads 2
```

The default output matrices act on logical irreducible row coefficients. Add
`--coordinates carrier` for the sparse working basis instead. Do not mix these
coordinate conventions. This command expands neither FEC nor LEC basis tensors.

The complete basis benchmark is reproduced with `bench/symrep_orbits.py`,
using `output_symrep_actions_20261002/native-current`, the preserved
`output_symrep_actions_20261002/before/symrep`, and
`--validated-bases output_symrep_research_20261002/release`. Other data/archive
arguments are unchanged from the preceding audit. To reproduce the action
benchmark after building `bench/symrep_actions_probe.cpp` with the same flags:

```bash
python3 bench/symrep_actions.py \
  --probe output_symrep_actions_20261002/actions-probe \
  --prepared-root output_symrep_sparse_20261002/validated \
  --chains output_symrep_actions_20261002/release \
  --output output_symrep_actions_reproduction --repeats 5
```

Use new output directories. `source-manifest.json` records the exact build and
hashes; `solver-changes.patch` records this pass's changes. Source snapshots,
raw logs and outputs are under `output_symrep_actions_20261002/`. Changes remain
uncommitted on `symrep-solver-dev`.
