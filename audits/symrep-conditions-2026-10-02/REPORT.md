# Reducing the condition map before recursive contraction

This audit follows the cached-template implementation in
[`../symrep-2026-10-02/REPORT.md`](../symrep-2026-10-02/REPORT.md). It evaluates
the next change on `symrep-solver-dev`: store and contract the independent
symmetry coefficients of `dlogmat`, including its equation representation.

## What was already reduced, and what changed

The preceding solver already adapted the equation space and restricted both
equation rows and candidate columns before SparseRREF. It did not solve a full
unorganized system and subsequently sort its answers. However, constructing
its cached contraction templates still expanded the adapted condition tensor
over the full alphabet and equation axes. Many terms were introduced only to
cancel after projection into the smaller blocks.

The new implementation reduces the condition map itself. Preparation saves
`dlogmat_multiplicity.wxf` and its coupling convention. Template construction
contracts these independent coefficients using cached maps between individual
irreps. These small component calculations do not depend on the number of
alphabet or equation copies. The growing recurrence combines sparse routes
between copies and passes compact solved coefficients to the next weight.

For the supplied heptagon, the adapted component tensor has 123,695 nonzero
entries. Its reduced description has **9,338 independent rational coefficients**,
a 13.25-fold reduction in stored nonzero coefficients. This is relative to the
adapted tensor: the original physical-coordinate input has only 5,450 nonzeros.
Symmetry adaptation still introduces fill-in compared with that original input.

This change preserves the matrices sent to SparseRREF, their dimensions,
their entries and the solved bases. Its improvement is in condition handling
and assembly, not a newly discovered reduction in the number of kernel variables.

## Why the condition compression is exact

Let `L` denote the letter representation and `E` the independent equation
representation. In the row-vector convention used by this repository,

```text
D (R_L(g) tensor R_L(g))^T = R_E(g)^T D.
```

Couple the two letter factors into irreducible copies using `C`. The matrix
`M = D C^T` has nonzero blocks only between equivalent irreps. For a rational
irrep `U_s`, a block from letter-pair copy `b` to equation copy `a` has the form

```text
M_s[a,b] = sum_e d_s[a,b,e] J_s[e]^T,
J_s[e] rho_s(g) = rho_s(g) J_s[e].
```

The matrices `J_s[e]` form a basis of `End_G(U_s)`. Only the rational numbers
`d_s[a,b,e]` depend on the supplied conditions; the rest follows from the
representation. The implementation reconstructs all of `M` from these
coefficients and requires exact equality during preparation. It also checks
that every `J_s[e]` commutes with the supplied generators. No numerical
tolerance or assumption of a scalar commutant is used.

The reduced file has one row per equation copy and columns indexed by
letter-pair copy and endomorphism coordinate, ordered within irrep sectors.
For heptagon its shape is 341 by 1,008, replacing a 1,191 by 1,764 component
map. The saved pair basis is checked on loading so a changed coupling
convention cannot silently reinterpret the coefficients. Legacy prepared
bundles are still accepted and compressed at process startup.

## Each recursive solve

1. Keep the preceding basis as irreducible copies. Represent its expansion
   using the image of one canonical generating vector per copy, retaining all
   endomorphism coordinates.
2. Couple the new candidate space and the new equation space with the fixed
   irreducible models. The factor order is reversed for LEC.
3. Combine preceding expansion coefficients with the reduced condition
   coefficients. Cached recoupling maps implement the contraction of the
   representation components. For the tested heptagon the raw matrices inside
   these local contractions are at most 36 by 36.
4. Assemble one rational primitive-corner matrix for each irrep sector, and
   compute its kernel with the unchanged SparseRREF backend.
5. Select generating kernel vectors and apply the fixed group actions to
   generate whole irreducible copies. Retain their compact coefficients for
   the next recursive step. Expanded tensors are exported for compatibility
   and checked against the requested generator actions.

The unknowns in step 4 are independent coefficients connecting equivalent
copies, rather than every component coefficient in the expanded basis.
Symmetry does not identify arbitrary different copies of the same irrep.
Consequently it cannot remove those multiplicity unknowns without solving
the constraints that connect them.

For example, coordinate exchange relates the equations `x1=0` and `x2=0`.
Dropping the latter while retaining two unrestricted unknowns changes the
kernel. Equations and unknowns must be reduced together; the intertwiner
description above is what makes generation of the remaining components valid.

## The remaining heptagon bottleneck

The selected group is D7 × C2, generated by cyclic, flip and parity. Its
rational irreps have dimensions `1,1,1,1,6,6`, and their endomorphism algebras
have rational dimensions `1,1,1,1,3,3`.

For a split irrep the familiar block is a scalar multiplicity matrix tensored
with an identity matrix. One solves only the scalar multiplicity matrix.
For either six-dimensional heptagon irrep, one independent coefficient instead
has three rational coordinates. The current QQ backend solves their rational
expansion. It reduces six component variables per candidate copy to three,
not to one rational scalar.

At LEC weight 3, the native system is 16,674 by 4,956 with 102,551 nonzeros.
The symmetry solver uses 2,832 rational corner variables across six blocks,
but those blocks collectively contain 1,364,502 nonzeros. The two largest
blocks have sizes 3,573 by 1,044 and 3,573 by 1,080. The reduced systems are
smaller but considerably denser; exact rational arithmetic and subsequent
lattice shortening dominate much of the remaining time.

An elimination backend operating directly in the cubic endomorphism field
could represent those two blocks as 1,191 by 348 and 1,191 by 360 matrices.
Its coefficients can still be stored as rational triples with rational
multiplication rules; this need not split the rational irreps into algebraic
two-dimensional models. It would compute kernel generators over the field
directly, avoiding extraction of generators from a threefold rational kernel.
That backend is **not implemented by this change**, and its speed must be
measured rather than inferred from the smaller dimensions.

The next priorities are preserving sparsity when choosing copy and coupling
bases, field-aware elimination for the heptagon sectors, and shortening only
the module generators that will actually be retained. Persisting reusable
recoupling tables can also reduce repeat-run setup; they currently live only
for one process. Arbitrarily dropping additional corner rows or coordinates
would not implement these optimizations correctly.

## Measured results

Median complete-chain wall seconds, three runs per method, two workers. Preparation is separate. “Before” means the immediately preceding cached-template implementation, not the initial slow prototype.

| Chain | Native | Before | Reduced conditions | Archive |
| --- | ---: | ---: | ---: | ---: |
| hexagon-fec through 8 | 0.126 | 0.115 | 0.113 | 0.042 |
| hexagon-lec through 6 | 2.940 | 0.608 | 0.599 | 0.189 |
| heptagon-fec through 4 | 0.076 | 8.574 | 5.242 | — |
| heptagon-lec through 3 | 0.058 | 16.585 | 11.116 | — |

Heptagon FEC improves by 38.9% and LEC by 33.0%. Hexagon changes are small enough that no substantial additional speedup is claimed. Heptagon still does not meet native performance; this optimization removes a measured assembly cost while retaining the same dense rational solve.

Final-weight phase medians:

| Phase | FEC4 before | FEC4 after | LEC3 before | LEC3 after |
| --- | ---: | ---: | ---: | ---: |
| Assembly, including templates | 4.522 | 2.508 | 10.376 | 4.794 |
| Template construction, part of assembly | 3.174 | 1.201 | 9.542 | 3.955 |
| SparseRREF | 0.432 | 0.440 | 1.721 | 1.790 |
| Lattice shortening | 0.764 | 0.770 | 2.739 | 2.744 |
| Copy lifting | 0.107 | 0.106 | 1.178 | 1.166 |

Single fresh preparation runs took hexagon-fec 0.014 s, hexagon-lec 0.016 s, heptagon-fec 0.304 s, heptagon-lec 0.315 s. These are separate setup measurements, not three-run medians.

Median peak RSS for heptagon FEC increases from 558 to 599 MiB because the new method retains additional small recoupling tables. LEC changes from 943 to 931 MiB. The improvement is primarily time, not a general memory reduction.

All 12 before/after chain pairs have byte-identical saved recurrence arrays. Exact recursive reference checks passed at every tested weight: native for all four cases, plus the archive for both hexagon cases. The probe and CLI regression also passed.

Detailed records: [results.json](results.json), [probe.log](probe.log), [regression.log](regression.log).

## Validation and reproduction

The algebra probe compares reduced contraction matrices entry for entry with
complete component contraction for C3, D3, D7 and Q8, in both directions.
It exercises non-scalar endomorphism coefficients, changed recursive
coefficients, template reuse and loading the compact condition representation.
The CLI regression checks kernels, FEC/LEC, continuation, vector terminals,
empty/full constraint spaces and saved-bundle integrity.

The benchmark runs fresh processes sequentially with two workers and three
repetitions. It compares every saved recurrence array byte for byte between
the preceding solver and the new implementation. It then independently
checks an exact invertible recursive basis change to the native bootstrap
at every weight, and to the unpacked package for both hexagon chains.

```bash
python3 bench/symrep_conditions.py \
  --before output_symrep_conditions_20261002/before/symrep \
  --native output_symrep_optimization_20261002/native_reference \
  --archive /home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap \
  --references output_symrep_validated_20261001 \
  --output NEW_OUTPUT_DIRECTORY
```

Compiler: GCC 14.4.0, C++20, `-O2`, FLINT 3.2.2, mimalloc disabled for both
general solvers and the native reference. The archive uses its own bundled
SparseRREF and executable configuration. Process wall times include startup,
mandatory validation, reconstruction and WXF exports; preparation and exact
reference comparisons are measured separately. OS file caches are not flushed.

Retained source snapshots, immutable bundles, arrays, full process logs and
timing files are under `output_symrep_conditions_20261002/validated/`.
