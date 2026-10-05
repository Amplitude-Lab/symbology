# Compact rational symmetry orbits: timing and exactness audit

The larger heptagon tests now run close to ordinary-solver speed while retaining
an exact irreducible basis description at every weight. FEC through weight 5
is 9.5% faster than ordinary solving; LEC through weight 4 is 2.1% slower.
Strictly equal or better timing on every tested chain has **not** been reached:
FEC through weight 4 costs an additional 13 ms, and LEC through weight 3 an
additional 23 ms. The smaller cases still expose symmetry setup and selection
costs. These are measured full-process results, not RREF-only timings.

| Full chain | Ordinary | Symmetry | Symmetry / ordinary | Archive |
| --- | ---: | ---: | ---: | ---: |
| heptagon-fec-w4 | 0.072 s | 0.085 s | 1.182× | — |
| heptagon-lec-w3 | 0.064 s | 0.087 s | 1.359× | — |
| heptagon-fec-w5 | 0.447 s | 0.405 s | 0.905× | — |
| heptagon-lec-w4 | 1.630 s | 1.665 s | 1.021× | — |
| hexagon-fec-w8 | 0.129 s | 0.104 s | 0.805× | 0.045 s |
| hexagon-lec-w6 | 2.967 s | 0.622 s | 0.210× | 0.211 s |

Each entry is the median of five fresh processes, run sequentially with two
workers and alternating native/symmetry order. All times include loading,
solving, the chosen output representation, WXF readback and checksum sealing.
Prepared inputs are reused. Optional exact comparison and expanded export are
measured separately. No existing result files are used to skip solving.
OS file caches were not flushed; unrelated user processes were not stopped.

Both project executables use GCC 14.4, FLINT 3.2.2, the same patched SparseRREF,
and `-O3 -march=native -mtune=native -flto`, with mimalloc disabled. The ordinary
helper calls the unchanged `bootstrap.hpp` extension functions. Its I/O helper
is the same as the symmetry solver's, including the improved tensor conversions.
The archive was rebuilt with the same optimization flags, but retains its own
SparseRREF, forced mimalloc allocator, pivot method and compact multiplicity
format. Its timings compare complete implementations, rather than isolated
elimination routines. The archive remains faster than the general hexagon path.

For context, the preceding audit measured heptagon LEC through weight 4 at
9.362 seconds, of which the last weight's explicit irreducible frame took
6.836 seconds. That earlier audit used `-O2`; it should not be treated as a
controlled compiler-identical speedup measurement. The current per-weight
phases identify where the work has actually been removed:

| Final weight | Assembly | QQ kernel | Modular actions | Irrep seed selection |
| --- | ---: | ---: | ---: | ---: |
| heptagon-fec-w4 | 0.0092 s | 0.0255 s | 0.0008 s | 0.0068 s |
| heptagon-lec-w3 | 0.0046 s | 0.0231 s | 0.0021 s | 0.0130 s |
| heptagon-fec-w5 | 0.0622 s | 0.1842 s | 0.0071 s | 0.0308 s |
| heptagon-lec-w4 | 0.1511 s | 0.7269 s | 0.0673 s | 0.3275 s |

The current LEC weight-4 carrier has 3,745 rows and 1,774,777 nonzero
coefficients. Its irreducible layout consists of 535 scalar copies and 535
six-dimensional rational copies: 1,070 certified generators specify all 3,745
rows. The normal solve stores the sparse carrier and an orbit recipe, without
materializing the rational frame or expanding the LEC recurrence into canonical
coordinates. Median peak RSS is 438 MiB versus 393 MiB for ordinary LEC through
weight 4; FEC through weight 5 uses 84 MiB versus 111 MiB.

## What changed

The logical basis is still `B_w = H_w C_w`, where `C_w` is sparse and `H_w`
orders rows by rational irrep, copy and component. The former implementation
computed and wrote large rational generator matrices and every entry of `H_w`.
The new path chooses all independent generators during the solve and stores
an exact formula for their group orbits.

A saved seed is

```text
v_(s,i) = e_i (R(h) - a I) P_s,
copy rows = v_(s,i) R(g_0), ..., v_(s,i) R(g_(d_s-1)).
```

`P_s` is the prepared rational primitive projector. The group element `h`,
integer `a`, carrier row `i`, and model words `g_j` are fully specified.
Every copy receives a common primitive integer normalization. This determines
an exact rational basis; evaluating it later requires no new kernel solve,
rank selection, basis choice or unknown coefficients.

The implementation inherits only modular generator actions during normal
recursion, directly from the exact kernel's free-coordinate chart. It checks
denominators and ranks, chooses primes congruent to 1 modulo the group order,
and retries bad primes using saved exact carriers. It never repeats the QQ
constraint solve merely because the modular certificate needs a new prime.
Rational generator actions remain exactly defined by the carrier recurrence
and are computed when an explicit consumer requests them.

The small certified models identify a generator whose eigenspace selects one
absolute component of each relevant rational irrep. For the heptagon, cyclic
symmetry provides this component. The two parity sectors share that large
eigenspace solve; parity is then resolved in a 535-dimensional space at LEC
weight 4. The scalar sectors likewise share the cyclic fixed space. These
two independent selections share the requested two-worker budget.

A nonzero minor in one residue-field embedding certifies independence over
the rational irrep's commutative endomorphism field. It therefore certifies
whole rational copies, rather than just the displayed finite-field component.
The small model checks establish the eigencomponent's dimension and a nonzero
pairing with the rational primitive projector. Character multiplicities and
the final dimension check establish completeness.

The initial free-coordinate rows do not always pair independently with the
primitive projector. When an involution `F` satisfies `F U F = U^-1`, changing
the seed to `e_i(U-a)P_s` gives the small matrix pencil

```text
Y_a = (lambda^-1 - a) Y + (lambda - lambda^-1) I.
```

Its rank is tested exactly modulo the selected prime. This reduced LEC weight
4's candidate set from 3,745 rows to 535. A full-row fallback remains when the
small candidate set cannot be certified. Groups without the spectral shortcut
use the fixed-corner selector; groups without a certified fixed-corner plan
retain the general rational frame construction. The logic checks the supplied
models and group relations, rather than dispatching on the name D7.

Additional improvements preallocate equation rows, assemble them in parallel,
remove empty equation rows before reconstruction, use direct CSR tensor
conversion, and resolve small rational coordinate charts with FLINT. Model
construction now checks generator intertwiners and derives all other actions
through group words; this preserves the same exact closure certificate.

Sharing a large flip-fixed space, combining joint eigenvalue equations, solving
a separate cyclic/parity spectral system for each sector, and changing the
SparseRREF pivot search all performed worse in the tested heptagon case. Their
exploratory snapshots are retained in the ignored local output directory; those
variants are not selected by the production implementation.

## Exact checks and compatibility

All six benchmark chains pass exact recursive changes of basis against the
ordinary solver at every weight. Both hexagon chains also pass exact checks
against the actual unpacked package. Every repeated basis is byte-identical
to its independently verified first run. Additional `--compare` runs for
heptagon FEC through 5 and LEC through 4 verify the complete rational residuals,
ambient carrier actions, frame intertwiners and frame invertibility.

Public algebra tests pass for ten groups: C3, C5, C7, D3, D5, D7, D7 × C2, S4,
C2 × C2 and Q8. They include dense rational conjugations, field and quaternionic
commutants, modular action recovery through private and general charts, empty
spaces, and a bad-prime fixture. Public CLI tests pass in both recurrence
directions, including vector terminals, continuation, zero/free kernels,
compact export, continued expanded chains and bundle-integrity rejection.
A separate check resumes old explicit-frame v2 chains into v3 for both FEC and
LEC, preserves the old basis files exactly, and verifies the resulting chains.

Separate real-data exports also passed exact verification against the native
chain. Materializing FEC through weight 4 took 0.558 seconds; LEC through
weight 3 took 1.158 seconds (one export each). These times are additional to
the compact solve and are recorded in `export-validation.json`. The fast
recursion itself does not require either export.

The new chain format is `symbology-symrep-chain-v3 factorized-orbits`.
`wN_orbits.tsv` uses orbit recipe v2 with carrier row, prefilter element and
shift. Orbit recipe v1 and full-frame chain v2 remain readable. The expanded
v1 chain format is unchanged. Higher-weight generator matrices are evaluated
recursively from `carrier_direction.tsv`, the exact seed actions and carriers.
`timings.tsv` distinguishes stored frame nonzeros from orbit seed counts.

The fast backend still computes the **full sparse QQ constraint kernel**.
It maintains the irreducible basis as a certified factorization. This result
does not establish equal performance for multiplicity-only elimination on the
heptagon; that backend remains available explicitly. Nor does a two-percent
large-case gap imply that all finite groups or future weights will be faster
than ordinary solving. Small-chain overhead and optional explicit expansion
remain material distinctions.

## Reproduction

Build `symrep`, `tests/symrep_probe`, `bench/native_symrep_reference.cpp`, and
the archive executables with the flags recorded in `source-manifest.json`,
then run:

```bash
python3 bench/symrep_orbits.py \
  --native output_symrep_orbits_20261002/native-release-current \
  --prepared-root output_symrep_sparse_20261002/validated \
  --archive /home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap \
  --archive-bin output_symrep_orbits_20261002/archive-release \
  --references output_symrep_validated_20261001 \
  --output output_symrep_orbits_reproduction --repeats 5
```

Full logs, bases and the timed source snapshot are in
`output_symrep_orbits_20261002/validated/`. `results.json` records every process
measurement and per-weight phase. `source-manifest.json` records source and
binary hashes, build flags and validation logs. Generated tensors and binaries
are excluded from the source changes. The user-facing format and commands are
in `docs/symrep.md`.
