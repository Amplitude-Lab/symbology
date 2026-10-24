# Experimental: isotypic (symmetry-adapted) sewing optimization

This folder archives an **unused experimental approach to optimizing the
sewing step**. Nothing in the accepted four-loop chain depends on it: the
accepted sewing imposes the seven symmetry generators (dihedral D4
crossing + five root-sign Galois involutions) as explicit product-
invariance rows, and **no stage of the chain decomposes any basis into
irreducible representations**.

## The idea

The sewing system is a linear problem that commutes with the symmetry
group acting on the word bases. Decomposing the FEC/LEC bases into
isotypic components (via the `symrep` library's central idempotents and
`adapt`) would block-diagonalize the ~3.3 M-column constraint system into
independent blocks per irreducible representation, potentially cutting
both the presolve and the certified-kernel cost. This is a pure solver
optimization: solving within the trivial isotype — which is what the
invariance rows already enforce implicitly — is mathematically equivalent;
only the execution strategy changes. The abelian Z2^5 Galois part is
already exploited cheaply in the accepted chain as a parity grading
(diagonal action on letters; the `galois` helpers verify homogeneity of
every basis row and the masked sewing skips equations that would mix
grades). The non-abelian D4 part is where an adapted basis would add
value and was never built.

## What was attempted and recorded

| Receipt | Result |
| --- | --- |
| `symrep_prepare` (archived in `stage_receipts.json` here) | killed, exit 143 at 371.6 s; output `symrep_prepared/` read by no downstream stage and no longer exists |
| `symrep_cached_prepare` (archived here) | failed, exit 1 at 242.7 s |
| `symrep_separating_prepare.json` (copied from the research tree) | one-core, 5 GB-capped setup experiment: PASS in 26.7 s, full representation output written to `symrep_separating_prepared/` — also never consumed |

Timing caveat recorded on the original day: the 371.6 s `symrep_prepare`
receipt started four seconds before `interior.wxf` was rewritten by the
adjacency condition merge, so it computed against the older, smaller
pre-merge tensor; with the committed final tensor the same computation
takes hours. Its wall time is not a cost estimate for this approach.

## Where the sources live

- `solver-source/build_symrep_separating.py` — center-separation
  improvement (rejects provably colliding central elements before
  constructing projectors).
- `solver-source/build_symrep_variant.py` — builds the `symrep_cached`
  variant binary.
- `solver-source/run_separating_prepare.py` — the capped re-run driver
  that produced the PASS receipt above.
- The `symrep` binary itself builds from the repository root
  (`make symrep`); `symrep.cpp`/`symrep_bootstrap.hpp` implement
  `prepare` (group closure, rational irreps, `adapt`, central
  idempotents, the kron pair space).

Any future attempt should first measure whether the isotypic split
survives the seam/affine structure of the accepted sewing (the affine
RHS and the +1 homogenizing column are not group-equivariant objects),
and should be benchmarked against the accepted invariance-row baseline
under the CPU-hours-within-RSS-tolerance criterion.
