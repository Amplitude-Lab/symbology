# Compact multiplicity recurrence: time, memory and exactness audit

The current hexagon MHV LEC chain is 3.05 times faster than the previous symmetry
solver and uses 53% less peak resident memory. It is within 6% of the unpacked
package's elapsed time, but still uses 59% more RAM. Hexagon FEC uses approximately
the package's RAM, while remaining 54% slower (an additional 26 ms). The larger
heptagon cases retain approximately ordinary-solver runtime; LEC peak RAM falls
by 9.6% from the previous symmetry implementation. No claim of universal speed
or memory parity is made.

## Full-process elapsed time (seconds)

| Full chain | Previous symmetry | Current symmetry | Ordinary | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 0.100 | 0.090 | 0.077 | — |
| Heptagon LEC ≤3 | 0.098 | 0.099 | 0.063 | — |
| Heptagon FEC ≤5 | 0.410 | 0.403 | 0.471 | — |
| Heptagon LEC ≤4 | 1.726 | 1.684 | 1.625 | — |
| Hexagon FEC ≤8 | 0.127 | 0.073 | 0.138 | 0.048 |
| Hexagon MHV LEC ≤6 | 0.655 | 0.214 | 2.996 | 0.202 |

## Peak resident memory (MiB)

| Full chain | Previous symmetry | Current symmetry | Ordinary | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 22.0 | 21.2 | 21.1 | — |
| Heptagon LEC ≤3 | 24.7 | 23.3 | 19.0 | — |
| Heptagon FEC ≤5 | 83.6 | 81.0 | 112.5 | — |
| Heptagon LEC ≤4 | 431.4 | 389.8 | 360.7 | — |
| Hexagon FEC ≤8 | 23.4 | 18.1 | 28.8 | 18.6 |
| Hexagon MHV LEC ≤6 | 120.0 | 56.3 | 371.8 | 35.5 |

Each entry is a median of five fresh processes. Methods run sequentially with
rotated ordering and two requested workers, without concurrent test/build jobs.
Timings include loading, solving, writing, exact WXF readback and checksum
sealing. Prepared inputs are reused; preparation and optional expanded exports
are excluded. Peak RSS is GNU time's maximum resident set, not just matrix
payload or final live memory. OS file caches were not flushed and unrelated user
processes were not stopped. Small cases have visible scheduling noise: the
heptagon LEC ≤3 difference from the previous solver is only about 1 ms.

Both current project implementations use GCC 14.4.0, FLINT 3.2.2 and
`-O3 -march=native -mtune=native -flto`, with mimalloc disabled. The ordinary
helper still calls the unchanged `bootstrap.hpp` extension functions. It was
rebuilt with the same I/O improvements and SparseRREF ownership patch as the
current symmetry solver. The preserved previous symmetry binary is from the
preceding audit, with the same optimization flags and allocator setting.

The archive binaries have the same optimization flags but retain their own
SparseRREF, pivot method 1 and forced mimalloc. This remains an end-to-end
implementation comparison, including each implementation's output format.

## What removed the cost

1. **Keep scalar multiplicities through storage as well as solving.** For split
   rational irreps, each map between equivalent copies is `h I_dimension`.
   The driver now stores and recurses on `h`. It enumerates global coupling
   channels without constructing the large component coupling matrices, full
   kernel components, ambient generator matrices or expanded recurrence tensor.
   Small certified models prove the component symmetry once. The complete
   basis remains defined at every step; there is no deferred solve or basis
   selection. The previous implementation already reduced the equations, but
   expanded and checked these redundant component structures after each solve.
2. **Solve independent large irrep blocks concurrently.** Workers share the
   requested budget, with larger blocks scheduled first. The comparison with
   serial compact solving showed a runtime gain at the cost of a few MiB of
   additional peak RAM; the substantial savings from compact storage remain.
3. **Remove temporary matrix copies.** WXF readback compares CSR entries directly
   and releases the encoded buffer first. SparseRREF now moves the completed
   rational reconstruction instead of deep-copying it. The symmetry kernel
   writes its row basis directly, in column order, without constructing and
   transposing a second kernel matrix. The reconstruction patch is included in
   `patches/` and applied by the setup script, so a clean checkout reproduces it.
4. **Reduce assembly allocations.** Noncanonical sparse rows are sorted and
   merged in their existing storage instead of building a tree node for each
   coordinate. Recurrence assembly reads its stored coefficient state directly
   instead of copying the entire state for every weight.

These changes do not change pivot choices, rational solutions or basis order.
The final saved bases are byte-identical to the preceding implementation's
corresponding files. For hexagon, the redundant expanded and candidate/kernel
files are replaced by the compact multiplicity files that already specified
those same bases.

Final-weight phase medians help locate the removed work. Kernel times for the
current split path include concurrent RREF wall time and primitive normalization;
previous kernel times also include component lifting and contraction.

| Final weight | Implementation | Assembly | Kernel | Component checks |
| --- | --- | ---: | ---: | ---: |
| Hexagon FEC ≤8 | previous | 0.005023 s | 0.032380 s | 0.007701 s |
| Hexagon FEC ≤8 | symmetry | 0.003964 s | 0.021438 s | 0.000000 s |
| Hexagon MHV LEC ≤6 | previous | 0.031567 s | 0.237128 s | 0.098454 s |
| Hexagon MHV LEC ≤6 | symmetry | 0.021438 s | 0.101730 s | 0.000000 s |

The remaining package advantage is not explained by elapsed time alone. At the
last hexagon weight, the current reduced equations contain 62,857 nonzeros for
FEC and 416,268 for LEC. The package's corresponding totals are 38,336 and
267,104. Their sector dimensions agree, but their coupling conventions and
recursive basis choices yield different sparsity. A profitable next algorithmic
target is therefore a better choice of sparse intertwiners and basis within each
multiplicity space, without expensive global basis reduction.

## Controlled experiments that were not adopted

`bench/symrep_reduced_probe.cpp` dumped the six LEC weight-6 reduced blocks and
solved the identical matrices with pivot methods 0, 1 and 2. Method 1 produced
204,880 kernel nonzeros versus 142,635 for method 0, with no speed advantage in
that trial. Method 2 produced 143,529. Copying the archive's pivot setting alone
would not reproduce its sparser recurrence.

Sorting solved multiplicity rows by their nonzero count, as the archive does,
also failed to improve this general solver: FEC's final matrix increased from
62,857 to 64,061 nonzeros and LEC's from 416,268 to 416,589. Three interleaved trials
showed no consistent speed gain and slightly higher RAM. The existing basis
order is retained. Trial records are included beside this report. Improvements
in sparsity must account for the actual representation/coupling conventions;
a row ordering rule alone was insufficient.

## Output and compatibility

| Full chain | Previous symmetry | Current symmetry | Ordinary | Hexagon package |
| --- | ---: | ---: | ---: | ---: |
| Heptagon FEC ≤4 | 0.063 MiB | 0.063 MiB | 0.057 MiB | — |
| Heptagon LEC ≤3 | 0.184 MiB | 0.184 MiB | 0.165 MiB | — |
| Heptagon FEC ≤5 | 0.959 MiB | 0.959 MiB | 0.985 MiB | — |
| Heptagon LEC ≤4 | 10.383 MiB | 10.383 MiB | 11.758 MiB | — |
| Hexagon FEC ≤8 | 0.638 MiB | 0.087 MiB | 0.344 MiB | 0.359 MiB |
| Hexagon MHV LEC ≤6 | 6.359 MiB | 0.847 MiB | 16.213 MiB | 4.685 MiB |

Output sizes include each implementation's entire output directory. Archive
sizes also include its reconstruction-cache files and reports, so these are
storage totals, not a comparison of basis payload alone.

Split multiplicity chains now use `symbology-symrep-chain-v4 multiplicity`.
`w1.wxf` is explicit; later weights save `wN_multiplicity.wxf` and copy layouts.
The complete logical component basis is available through `symrep expand`.
Non-split groups retain the existing generic multiplicity and factorized paths;
heptagon's default remains the certified-orbit factorization of a full sparse
QQ kernel. This work does not turn heptagon elimination into a multiplicity-only
solve.

A separate one-run optional export took 0.032 s / 15.8 MiB for hexagon FEC ≤8,
and 0.150 s / 67.8 MiB for MHV LEC ≤6. These costs are additional to compact
solving, not hidden in the table. Every expanded weight is byte-identical to the
previous solver's expanded tensor. Genuine previous v1 chains were also resumed
in both directions; old files were preserved, new compact weights matched a
fresh solve, and `--compare` passed. Public tests cover compact seed-only resume,
zero/free kernels, vector terminals, export and continued expanded chains.

## Exact validation

The complete reference audit in `reference-validation.json` ran exact recursive
changes of basis at every weight against ordinary chains, and against the actual
unpacked package for both hexagon chains. It also ran full-matrix heptagon
comparisons. The final timing run reuses those proofs only after requiring the
same prepared-bundle checksum and byte equality of **every saved native and
symmetry basis file**, then checks all repeated runs against that basis. Its
certificate path and SHA-256 are recorded in `results.json`; these are exact
basis certificates, not numerical spot checks. Original basis data also matches
the preserved previous implementation at every weight.

The release passes the ten-group public algebra suite, direct-row versus
SparseRREF kernel extraction checks, the expanded CLI regression suite and all
90 ordinary-solver numerical cases. Patch application was checked against the
pinned upstream source and the setup script is idempotent. `git diff --check`,
Python syntax checks and shell syntax checks pass.

## Reproduction

Build `symrep`, `tests/symrep_probe` and `bench/native_symrep_reference.cpp` with
the compiler settings in `source-manifest.json`, after running
`scripts/setup-sparserref.sh`. The exact benchmark command is:

```bash
python3 bench/symrep_orbits.py \
  --native output_symrep_compact_20261002/native-current \
  --previous output_symrep_compact_20261002/before/symrep \
  --prepared-root output_symrep_sparse_20261002/validated \
  --archive /home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap \
  --archive-bin output_symrep_orbits_20261002/archive-release \
  --references output_symrep_validated_20261001 \
  --output output_symrep_compact_reproduction --repeats 5
```

This command runs all exact comparisons again. To reuse the completed basis
certificates, add `--validated-bases output_symrep_compact_20261002/final`;
any changed prepared checksum or saved basis byte causes failure. Full logs,
source/binary/header snapshots, intermediate trials and outputs are retained
under `output_symrep_compact_20261002/`. The report and compact JSON evidence
are retained in this audit directory. Changes remain uncommitted on
`symrep-solver-dev`.
