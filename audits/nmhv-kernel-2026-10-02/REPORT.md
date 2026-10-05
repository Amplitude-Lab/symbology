# NMHV kernel acceptance and default selection

This audit tests **basis construction and sewing kernels**, not newly solved
collinear-normalized NMHV amplitudes. The two are different acceptance tests.
The heptagon fixture is built from equations (7) and (8) of
[He et al., arXiv:2511.09669v2](https://arxiv.org/html/2511.09669v2).
It has 147 independent last-entry vectors with 15 physical terminal components.
The hexagon fixture is the archive's joint E/Et NMHV seed, with 18 vectors and
five terminal components. No MHV last-entry seed is reused as an NMHV test.

The machine is the same Intel Core Ultra 7 265 with 62.1 GiB RAM and no swap.
Each timed comparison uses eight workers, a fresh process and output directory,
and the **same executable** containing both kernels and the same patched
SparseRREF dependency. Three sequential repetitions rotate execution order.
Wall time includes input, output and readback. Peak RSS comes from GNU time;
there are no concurrent builds or correctness comparisons in these intervals.
Caches are not flushed. Sub-0.1-second sewing differences are too small for a
strong speed claim at GNU time's two-decimal reporting resolution.

## Results and defaults

Three-run medians for complete recursive LEC chains, including output:

| NMHV workload | Original time | New time | Original peak RSS | New peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Hexagon, weights 1–5 | 1.30 s | 1.06 s | 186 MiB | 121 MiB |
| Heptagon, weights 1–4 | 23.41 s | 20.21 s | 7.15 GiB | 4.65 GiB |

The hexagon dimensions are 18, 58, 166, 450, 1171. The heptagon dimensions are
147, 937, 4829, 21826. **Exact recursive changes of basis pass at every one of
these weights**, not just dimension comparisons. The full heptagon comparison
uses approximately 5.7 GiB and 124 s; validation is deliberately outside the
performance measurements. The hexagon reference comes from the archive's
independent irrep recurrence, the heptagon reference from the original kernel.

The unmodified archive's compressed NMHV recurrence through weight 5 took
about 0.16 s of its internally reported runtime in the exploratory reference
run. Its output is compact multiplicity data, unlike the expanded carrier
files above. Our retained default multiplicity backend constructs the same
1171-dimensional space in about 0.10 s of reported command time, with exact
recursive agreement. These single-run observations are not controlled
speedup claims between the two packages.

Both original and streamed sewing agree exactly. A separate weight-6 heptagon
NMHV sewing test (FEC5 × LEC1, 177 unrestricted solutions) takes 5.30 s and
1.36 GiB with the original method, versus 4.73 s and 0.55 GiB streamed. This is
one run per method. Its invariant space has the independently published
**24** dimensions. Together with the lower tests, weights 2, 4, 6 reproduce
**5, 11, 24**, respectively; no collinear normalization is claimed here.

Small weight-4 heptagon sewing is an exception: roughly 0.13 s originally
versus 0.34–0.36 s streamed, though the latter uses less RAM. Accordingly the
default is **adaptive**, not unconditional streamed sewing: reduce the side
with fewer basis coefficients first. This selects the original left
contraction for FEC3 × NMHV LEC1, and the streamed right contraction for
FEC5 × NMHV LEC1 and the large MHV test. This is a general dimension heuristic,
not a polygon/helicity test or a guarantee of universal optimality.

Promoted defaults:

- `bootstrap` / `compute_rhs`: streamed extension kernels.
- Sewing: `auto`, choosing the contraction order as above; either method can
  still be selected explicitly.
- `compute_rhs`: restricted projection for its MHV, one-step-tail workflow.
- `symrep extend`: streamed carrier kernels; split scalar multiplicity solving
  remains the default where applicable. NMHV terminal actions were prepared
  and checked recursively for hexagon through weight 5 and heptagon through
  weight 2, in addition to the larger ordinary-coordinate kernel tests.
- `make`: builds `symrep` alongside the existing executables.
- Dependency: pinned SparseRREF `5bbee55` (v0.3.6) plus the seven tested patches.
  The slower v0.4.2 deployment trial from the preceding audit is not promoted.

A fresh four-loop MHV calculation with **no strategy flags** completed in
**68.49 s**, with **5.21 GiB** maximum single-process RSS and **5.43 GiB**
sampled summed RSS. All **67 output WXF files** are byte-identical to the
previous explicitly optimized calculation, and all ten archived input tensors
match the repository seeds. Twelve independent published word coefficients
also pass. This validates the new defaults end to end; it is comparable to the
previous approximately 68-second MHV result, not a further MHV speedup claim.

Historical coefficient CRCs are explicitly tested with `original` strategies.
The new-default regression separately reproduces all six physical two-/three-
loop E/R/boundary tensors exactly. Public CI also runs the new NMHV CLI checks.

Final checks pass: 159 general exact kernel cases; 101 independently assembled
factored extension/sewing cases; the symmetry CLI suite (including vector
terminals and continuation); the new NMHV CLI suite; dependency patch checks;
and native numerical/physical regressions. The Python execution/export suite
reports 58 passed and 6 skipped. `git diff --check` is clean.

## What the NMHV tests exposed

The first streamed version was slower in both complete NMHV LEC chains despite
using less memory. These results were retained rather than reported as wins.
Three general changes address the measured work:

1. Clear denominators by scaling whole local equations, and all letter slices
   belonging to the same older component together. Independent letter scaling
   would change the equations and is never used.
2. Parallelize independent CRT, rational reconstruction and coordinate-lifting
   rows. Reuse the already validated lifted kernel instead of constructing it
   a second time. Disjoint relation components guarantee unique lifted indices.
3. Use the **complete CRT modulus** in the deterministic residual-height
   certificate. Hexagon NMHV weight 5 has a 101-bit residual bound: it exceeds
   one 61-bit prime but is below the accumulated 121-bit modulus. Falling back
   to a full rational residual product was unnecessary. Generic tests include
   a 161-bit bound certified only after accumulating a 181-bit modulus.

The exact rank lower bound, private free coordinates and zero-residual proof
are all retained. A failed reconstruction is never accepted on timing grounds.

Sewing was also generalized from a scalar terminal to arbitrary terminal
components. The local equations now carry both the condition and terminal
indices. This also supports a deeper backward recurrence tensor. Independent
random rational tests compare it with an explicit four-index contraction.

## References and scope

The archive reference is generated by its own `bin/irr_lec NMHV`, using its
own headers and dependency. The exporter now reads the non-scalar weight-zero
representation and converts its terminal coordinates to the physical basis.
The existing scalar multiplicity route remains the appropriate hexagon
symmetry backend; a carrier-kernel improvement is not a reason to discard it.

Heptagon cyclic and flip actions include their action on R-invariant terminal
components. They reproduce the published pre-collinear invariant dimensions
5, 11 and 24 at weights 2, 4 and 6. These counts are an independent physical check;
they are not inputs to the solve. Parity of a scalar MHV seed is not silently
imposed on NMHV's multi-component terminal space.

Large generated files and exploratory runs are under
`output_nmhv_kernel_20261002/`; compact measurements and validation logs are
retained alongside this report. Benchmark executables are generated from
`bench/nmhv_kernel_benchmark.cpp`, `bench/nmhv_kernel_check.cpp` and the archive
exporter. The full SparseRREF patch list is checked by
`scripts/setup-sparserref.sh --check`.

## Reproduction

Build the ordinary repository tools and the two NMHV probes, using the same
compiler/dependency options for both methods:

```bash
make all bench/nmhv_kernel_benchmark bench/nmhv_kernel_check
make check-kernel check-symrep check-nmhv check-public
```

Build the archive and `bench/hexagon_symrep_reference` as described in
[the symmetry benchmark instructions](../../docs/symrep.md). Set `HEXAGON` to
the unpacked archive directory, then use a fresh output directory:

```bash
NMHV_AUDIT="$PWD/output_nmhv_reproduce"
mkdir "$NMHV_AUDIT"
(cd "$HEXAGON" && THREADS=8 bin/irr_lec NMHV 5 "$NMHV_AUDIT/archive")
bench/hexagon_symrep_reference LEC "$NMHV_AUDIT/archive" 5 "$HEXAGON/data" "$NMHV_AUDIT/hex-reference"
bench/nmhv_kernel_check fixture "$NMHV_AUDIT/hept-fixture"
python3 bench/nmhv_kernel_benchmark.py --hexagon "$HEXAGON" \
  --hexagon-seed "$NMHV_AUDIT/hex-reference/w1.wxf" \
  --heptagon-seed "$NMHV_AUDIT/hept-fixture/LEC_1.wxf" \
  --output "$NMHV_AUDIT/timings" --threads 8 --repeats 3
bench/nmhv_kernel_check compare backward "$NMHV_AUDIT/timings/hexagon-streamed-0" "$NMHV_AUDIT/hex-reference" 5
bench/nmhv_kernel_check compare backward "$NMHV_AUDIT/timings/heptagon-streamed-0" "$NMHV_AUDIT/timings/heptagon-original-0" 4
```

The driver enforces a 28 GiB virtual-memory limit and a 300-second timeout per
run. Exact comparison time is separate from construction time. These fixtures
represent NMHV last-entry and bulk spaces, not complete physical amplitudes.

The old private multi-pair baseline manifests cannot all be replayed in this
checkout: their project data and previously generated `output/` seed files
are absent. This limitation is not counted as a passed test. The public native
regression reconstructs its inputs in a temporary directory and passes all ten
historical lower-loop MHV CRCs with the explicit legacy strategies, then checks
the new defaults against the same basis-independent physical tensors.
