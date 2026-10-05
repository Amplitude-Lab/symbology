# Streamed symmetry solving: further time and memory optimization

The new implementation reduces peak RAM without changing a saved basis, pivot
choice or mathematical algorithm. Hexagon MHV LEC drops from 56.0 to 39.7 MiB
(29% less), and heptagon LEC from 389.7 to 334.7 MiB (14% less). Their elapsed
times remain approximately unchanged. Hexagon FEC is 9% faster in the main
comparison and uses 12% less RAM.

Heptagon LEC now measures 1.628 s versus ordinary solving's 1.621 s, with lower
peak RAM: 334.7 versus 384.2 MiB. Hexagon MHV LEC measures 0.219 s versus the
unpacked package's 0.217 s, using 39.7 versus 37.1 MiB. Hexagon FEC still takes
0.062 s versus the package's 0.048 s, although it uses less RAM. These are
measured comparisons for these chains, not a guarantee for all representations.

## Full-process elapsed time (seconds)

| Complete chain | Previous symmetry | Current symmetry | Ordinary | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 0.092 | 0.106 | 0.081 | — |
| Heptagon LEC ≤3 | 0.095 | 0.091 | 0.066 | — |
| Heptagon FEC ≤5 | 0.411 | 0.406 | 0.470 | — |
| Heptagon LEC ≤4 | 1.643 | 1.628 | 1.621 | — |
| Hexagon FEC ≤8 | 0.067 | 0.062 | 0.146 | 0.048 |
| Hexagon MHV LEC ≤6 | 0.219 | 0.219 | 3.008 | 0.217 |

## Peak resident memory (MiB)

| Complete chain | Previous symmetry | Current symmetry | Ordinary | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 21.3 | 20.8 | 21.0 | — |
| Heptagon LEC ≤3 | 23.2 | 22.5 | 19.1 | — |
| Heptagon FEC ≤5 | 79.1 | 75.4 | 111.0 | — |
| Heptagon LEC ≤4 | 389.7 | 334.7 | 384.2 | — |
| Hexagon FEC ≤8 | 18.3 | 16.2 | 28.8 | 18.6 |
| Hexagon MHV LEC ≤6 | 56.0 | 39.7 | 361.2 | 37.1 |

The previous binary is exactly the current binary from the preceding
[compact-multiplicity audit](../symrep-compact-2026-10-02/REPORT.md). All six cases
were rerun with both versions, the ordinary reference, and the unpacked hexagon
package where applicable. Each cell is the median of five fresh sequential
processes, rotating method order, with two requested workers. No build or test
job ran concurrently with these measurements. OS caches were not flushed and
unrelated user processes were not stopped.

Elapsed time includes loading, recursive solving, writing, exact WXF readback
and final integrity checksums. Reusable prepared inputs, optional expanded
exports and optional full verification are excluded. Peak RAM is GNU time's
maximum resident set for the full process, including temporary allocations.
Both project solvers and the ordinary reference use GCC 14.4.0, FLINT 3.2.2,
`-O3 -march=native -mtune=native -flto`, and no mimalloc. The ordinary reference
was freshly rebuilt and calls the unchanged `bootstrap.hpp` recurrence.
The archive has the same optimization flags but retains its own SparseRREF,
pivot method 1 and forced mimalloc; this remains an end-to-end implementation
comparison.

The small heptagon FEC ≤4 result was noisy: the main five-run median increased
from 0.092 to 0.106 s. A separate nine-pair alternating diagnostic measured
0.093 s previous versus 0.089 s current, so it did not confirm a repeatable
regression. That diagnostic omits the GNU time wrapper and is kept separate
in `small-fec4-trials.json`; it does not replace the main table. The main
heptagon LEC ordinary RAM samples span 360.9–391.9 MiB, while current symmetry
samples span 330.2–339.4 MiB. All individual observations are in `results.json`.

## What changed

1. **Stream each reduced block into its solve.** Previously, assembly retained
   every constraint sector before elimination began. Now a worker assembles
   one sector, solves it immediately, and releases its constraint storage
   before taking another. Other workers can assemble or solve concurrently.
   Sectors are scheduled by estimated size, largest first. Kernels keep their
   original sector and row order. Each worker has a private, reusable
   single-worker SparseRREF workspace; a waiting assembly worker does not
   introduce extra simultaneous elimination work.
2. **Avoid retaining unused recurrence state.** The final split-model kernel
   is saved directly without constructing the coefficient-state copy used
   only by a subsequent weight. Continued runs reconstruct that state from
   the exact saved compact kernel, as they already did.
3. **Release the final rational carrier before orbit selection.** The factorized
   path derives the modular actions, writes and validates the exact carrier,
   then frees that matrix before the orbit-selection workspaces are allocated.
   The mathematical work and saved output are unchanged. A bad-prime retry
   reloads the saved exact carrier and rebuilds modular actions; it does not
   repeat the rational constraint solve. An incomplete selection still leaves
   an unsealed, non-resumable output. Intermediate weights retain the carrier
   tensor needed for the next recursive contraction.

These are general scheduling and ownership changes. They do not introduce a
hexagon/heptagon group dispatch, expand LEC, change the basis normalization,
or switch a rational calculation to an approximate result. The heptagon fast
path still solves a full sparse QQ kernel and represents its irreducible basis
with certified orbit recipes; this work does not make its elimination a
multiplicity-only solve.

## Where the remaining time goes

Assembly and elimination now overlap for split-model solving. `pipeline_s`
is the actual elapsed time of their joint pipeline, including serial planning.
`planning_s` is its serial template/mapping setup. `assembly_work_s` and
`rref_work_s` sum per-sector durations, which can overlap; adding them does not
give wall time. Primitive normalization remains separately timed in
`shorten_s`. Legacy separate wall-time columns are `nan` in the normal streamed
path; `--compare` retains separate assembly and elimination for inspection.
Historical timing rows are padded on resume.

| Final weight | Previous assembly + RREF wall time | Current pipeline wall time | Current summed assembly work | Current summed RREF work |
| --- | ---: | ---: | ---: | ---: |
| Hexagon FEC ≤8 | 0.021051 s | 0.018007 s | 0.003572 s | 0.029284 s |
| Hexagon MHV LEC ≤6 | 0.116425 s | 0.118696 s | 0.027269 s | 0.200835 s |

At heptagon LEC weight 4, median current phase times are 0.160 s assembly,
0.736 s QQ kernel, 0.070 s inherited modular actions and 0.330 s orbit selection.
These exclude earlier weights, loading and output costs. Eliminating the
remaining representation overhead requires a better sparse orbit-selection
algorithm or a better multiplicity-space basis; faster constraint assembly
alone cannot remove it.

The hexagon reduced matrices retain their previous nonzero counts (62,857 at
FEC weight 8 and 416,268 at LEC weight 6), still above the package's 38,336 and
267,104. Changing scheduling reduces simultaneous storage but cannot remove
this difference in recursive basis sparsity.

## Column-ordering experiments

`bench/symrep_reduced_probe.cpp` now tests column weights on identical dumped
reduced matrices: original index, column degree, reverse degree and sum of
incident row lengths. Each policy was tried twice with one and two workers on
the six hexagon LEC weight-6 blocks, using pivot method 0. Every result passed
an exact zero-residual check. Kernel nonzeros were:

| Column weight | Kernel nonzeros |
| --- | ---: |
| Original index | 142,635 |
| Degree | 144,114 |
| Reverse degree | 141,991 |
| Incident row lengths | 143,848 |

The runtime ranges overlap, and the best sparsity change is less than 0.5%.
No ordering policy was adopted on this evidence. These weights mainly affect
pivot tie-breaking in method 0; they are not a full sparse reordering or a
proof that better bases are impossible. The measurements are retained in
`ordering-trials.json`.

## Exactness and continuation

- The ten-group public algebra suite passed, including C3/C5/C7, D3/D5/D7,
  D7×C2, S4, C2×C2 and Q8, direct SparseRREF extraction comparisons and
  non-split/quaternionic coefficients.
- The template tests now compare streamed and materialized contractions in
  both recursive directions. They deliberately throw from a consumer and
  require all submitted jobs to finish before the error returns.
- Public CLI tests compare saved bases at one and three workers, continue
  streamed compact chains, check final-carrier release against `--compare`,
  and resume from the early-written carrier. They cover both directions,
  vector terminals, zero/free kernels and a public nonsplit C3 example.
- For every weight of all six production cases, every saved native and current
  symmetry `w*` basis file matches the preceding exactly verified audit
  byte-for-byte; all five repetitions agree as well. The reference proof is
  reused only after checking those bytes and the prepared-bundle checksum.
  The previous audit contains exact invertible changes of basis against
  ordinary chains and the unpacked hexagon package, plus full heptagon checks.
  The certificate paths/hashes are recorded in `results.json`, and the original
  full proof is retained in
  [the preceding reference validation](../symrep-compact-2026-10-02/reference-validation.json).
- Separate real-data continuation checks resumed previous-format hexagon
  chains and newly saved heptagon carriers to FEC8/LEC6 and FEC5/LEC4,
  respectively. All saved bases match fresh solves; old/new timing columns
  parse correctly. See `continuation-checks.json`.

The ordinary numerical core and SparseRREF patches are unchanged in this pass.
No new file format or export conversion is required. Only timing metadata gains
columns; basis output payload is unchanged. Optional full expansion costs
remain additional to compact solving.

## Reproduction

Build settings and source/binary hashes are recorded in `source-manifest.json`.
The production changes relative to the preceding audit are retained in
`solver-changes.patch`. Source snapshots, raw process logs and output bundles
are under `output_symrep_stream_20261002/release/`.

```bash
python3 bench/symrep_orbits.py \
  --native output_symrep_stream_20261002/native-current \
  --previous output_symrep_stream_20261002/before/symrep \
  --prepared-root output_symrep_sparse_20261002/validated \
  --archive /home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap \
  --archive-bin output_symrep_orbits_20261002/archive-release \
  --references output_symrep_validated_20261001 \
  --validated-bases output_symrep_compact_20261002/release \
  --output output_symrep_stream_reproduction --repeats 5
```

Omit `--validated-bases` to rerun all exact reference transformations instead
of requiring identical previously verified basis bytes. The command's output
directory must be new. Changes remain uncommitted on `symrep-solver-dev`.
