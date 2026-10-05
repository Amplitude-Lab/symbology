# Combined integrability and symmetry benchmark — 2026-10-04

The joint staged strategy is implemented in the shared kernel library and is
independent of a polygon, loop count, or expected nullity. This audit adds a
reproducible comparison of solving orders. **All 45 timed outputs passed
independent exact full-space checks**, and no resource guard fired. A further
15 runs replayed the same inputs with reduced interface equations and also
passed. These are sewing-kernel benchmarks, before additional physical
boundary conditions.

Read the [mathematical and implementation explanation](../../docs/combined-constraints-benchmark.md)
first. The runner is [`bench/combined_constraints_benchmark.py`](../../bench/combined_constraints_benchmark.py),
with the three solving paths and independent checker in
[`bench/combined_constraints_benchmark.cpp`](../../bench/combined_constraints_benchmark.cpp).
The committed [default suite](../../bench/combined_constraints_cases.json)
requires only repository data. This recorded run also includes data from the
unpacked hexagon package.

## Controlled comparison

Each method loads identical sealed factors and paired generator matrices:

* **Integrability first:** compute its complete kernel, then restrict by
  generator invariance using factored rational contractions.
* **Joint global:** impose both constraint families from presolve onward,
  assemble all remaining modular equations, and eliminate globally.
* **Joint staged:** impose both families from presolve onward, discover rank
  adaptively, and intersect the candidate with complete residual operators.

All use the current patched SparseRREF stack, the same binary, two workers,
raw interface equations, and process guards. Thus the first method is an
optimized baseline; this is not a comparison with an old unmodified release.
The independent checker uses the original rational sewing assembler followed
by explicit Kronecker actions and compares canonical rational row spaces.
It shares exact arithmetic infrastructure, but not the factored assembly and
symmetry application used by the timed methods.

Three fresh sequential processes per method, rotated method order. Preparation
and independent validation are measured separately. Intel Core Ultra 7 265,
20 logical CPUs, approximately 62 GiB RAM; concurrent host activity was not
stopped. GCC 14.4.0, `-O3 -march=native -mtune=native -std=c++20`, FLINT/GMP/TBB,
no mimalloc. Compiler overrides used:

```bash
make bench/combined_constraints_benchmark \
  CXX=/tmp/symbology-audit-toolchain/env/bin/x86_64-conda-linux-gnu-c++ \
  CXXFLAGS='-O3 -march=native -mtune=native -std=c++20 -I. -I/tmp/symbology-audit-toolchain/env/include' \
  LDLIBS='-L/tmp/symbology-audit-toolchain/env/lib -Wl,-rpath,/tmp/symbology-audit-toolchain/env/lib -lflint -lgmp -ltbb' \
  MIMALLOC=0
```

Per-process guards: 600 seconds, 6 GiB sampled tree RSS, 8 GiB address space,
and a 4 GiB host-available-memory floor. GNU time provides peak single-process
RSS; the monitor separately records process-tree RSS at 100 ms intervals.
Reported values are medians, not maxima across repetitions.

## Results

| Case | Method | Solve seconds | Process wall seconds | Peak process MiB | Final dimension |
| --- | --- | ---: | ---: | ---: | ---: |
| S3, 1+5 | Integrability first | 0.002125 | 0.102 | 9.93 | 7 |
| S3, 1+5 | Joint global | 0.000581 | 0.102 | 9.61 | 7 |
| S3, 1+5 | Joint staged | 0.000606 | 0.102 | 9.93 | 7 |
| S3, 3+3, rational coordinates | Integrability first | 0.000933 | 0.102 | 9.77 | 7 |
| S3, 3+3, rational coordinates | Joint global | 0.002604 | 0.102 | 9.77 | 7 |
| S3, 3+3, rational coordinates | Joint staged | 0.004193 | 0.102 | 9.93 | 7 |
| Heptagon, 2+2 | Integrability first | 0.047857 | 0.102 | 16.48 | 1 |
| Heptagon, 2+2 | Joint global | 0.043755 | 0.102 | 15.39 | 1 |
| Heptagon, 2+2 | Joint staged | 0.058481 | 0.102 | 15.39 | 1 |
| Heptagon, 4+2 | Integrability first | **2.131670** | 2.215 | **154.40** | 2 |
| Heptagon, 4+2 | Joint global | 1.222600 | 1.309 | 98.82 | 2 |
| Heptagon, 4+2 | Joint staged | **1.062740** | 1.110 | **54.33** | 2 |
| Hexagon package, 4+2 | Integrability first | 0.016865 | 0.102 | 10.08 | 9 |
| Hexagon package, 4+2 | Joint global | 0.009080 | 0.102 | 9.93 | 9 |
| Hexagon package, 4+2 | Joint staged | 0.010884 | 0.102 | 9.77 | 9 |

Solve time includes interface construction, presolve, elimination, and required
certificates. Process wall additionally includes load/seal verification,
normalization, serialization, and exact byte readback. The monitor's wall time
has approximately 100 ms granularity, evident for the small cases. Millisecond
toy timings under contention should not be interpreted as stable speed ratios.

On heptagon 4+2, joint staging reduced median solve time by **50.1%** and peak
process RSS by **64.8%** relative to integrability first. Joint global already
captures much of the improvement: staging adds a **13.1%** time reduction and
**45.0%** RAM reduction relative to that method. Small cases can favor global
elimination. These measurements do not justify making staging the unconditional
default for every recursive basis problem.

## What changed inside the solve?

| Case | Product dimension | Presolved columns, integrability only | Presolved columns, joint | Complete integrability dimension | Invariant dimension |
| --- | ---: | ---: | ---: | ---: | ---: |
| S3, 1+5 | 3 × 21 = 63 | 28 | 7 | 28 | 7 |
| S3, 3+3, rational coordinates | 10 × 10 = 100 | 44 | 26 | 28 | 7 |
| Heptagon, 2+2 | 28 × 118 = 3,304 | 511 | 24 | 1 | 1 |
| Heptagon, 4+2 | 308 × 118 = 36,344 | 8,967 | 2,008 | 2 | 2 |
| Hexagon package, 4+2 | 26 × 25 = 650 | 435 | 238 | 41 | 9 |

For heptagon 4+2, the complete integrability kernel is already invariant.
Symmetry therefore does not reduce the final dimension in this fixture.
It supplies useful redundant short equations that expose a cheaper coordinate
system before elimination. This is a different mechanism from removing a large
non-invariant final space, which occurs in the S3 and hexagon examples.

The first heptagon 4+2 repetition records:

* Integrability first: 159,202 modular rows, 8,967 columns, 5,782,725 nonzeros.
* Joint global: 176,787 modular rows, 2,008 columns, 3,545,342 nonzeros.
* Joint staged: an initial 142,669-row, 1,314,990-nonzero front on 2,008
  columns, followed by complete residual intersection **7 → 2**. The final
  row basis passed to reconstruction has 2,006 rows and 3,402 nonzeros.

The complete-operator flag and exact residual-height certificate are present
in the staged log. The two-dimensional answer is a result of solving every
condition, not a requested stopping dimension. The logs do not identify which
individual operator caused the intersection, so no such attribution is made.

The S3 examples describe commuting three-letter symbols with the full symmetric
group generated by a cycle and a flip. The 3+3 fixture additionally changes
both factor bases using a rational shear. This tests non-monomial
actions and product-invariant pairings without polygon-specific data.

Heptagon inputs are the repository's `dlogmat_E6.wxf`, `FEC_1.wxf`, `LEC_1.wxf`,
and cyclic/flip/parity matrices. Hexagon uses the unpacked package's
`dlogmat_full_ca_pd3.wxf`, cyclic and flip matrices, and scalar weight-one
seeds exported from its MHV irrep chains into recurrence tensors. Its parity
is already generated by cyclic. The preparation helper uses the stated
condition tensor for both chains. We are timing three solving orders on those
inputs, not rerunning the package's specialized solver as a fourth competitor.
All paths and hashes are retained in [cases.json](cases.json) and
[results.json](results.json).

## Validation and replay

Every case's independent exact checker accepted all nine timed kernels. The
most expensive check, heptagon 4+2, took 2.013 s and 322.39 MiB peak process
RSS. Those costs are outside the solver table and are retained in the JSON.
Input hashes were unchanged after the benchmark.

The saved `replay_cases.json` was then used to run all three methods once
with `--local-reduction reduced`. All 15 additional kernels agreed with the
same independent exact reference. The results are in
[replay-reduced-results.json](replay-reduced-results.json); this single-repetition
run tests replay and policy equivalence, not a new performance ranking.

Additional checks, recorded in [validation.json](validation.json), verify that:

* unknown methods, altered sealed inputs, and action/basis dimension mismatches
  fail before writing an answer;
* a valid kernel from the same-shaped S3 product in different factor coordinates
  is rejected by the exact-space checker.

The new benchmark compiled successfully; Python bytecode compilation and
`git diff --check` passed. No production solver default or mathematical kernel
implementation was changed for this audit, so the existing full regression
suite was not rerun solely for the new benchmark and documentation.

The complete raw run, including logs, WXF outputs, frozen inputs, binary and
source snapshots, is
[`output_combined_constraints_20261004/`](../../output_combined_constraints_20261004/REPORT.md).
This audit keeps a copy of its machine-readable results. For a fresh benchmark
with the current binary and exactly those prepared inputs:

```bash
python3 bench/combined_constraints_benchmark.py \
  --cases output_combined_constraints_20261004/replay_cases.json \
  --output output_combined_constraints_review --threads 2 --repeats 3
```

Use a new output directory; the runner refuses to overwrite an existing one.
The full result JSON includes source and binary SHA256 hashes, git status,
host metadata, preparation/check costs, each command, each run, and all
input/output hashes. `integrability_dimension` in timed records is populated
only for `integrability-first`; the joint records use zero as an unused field.
The check records provide the independent complete integrability dimension.

The [five-loop reference](../fiveloop-day-2026-10-04/REPORT.md) remains the
previous 8+2 physical calculation from cached ordinary FEC8: 7,796.38 s and
19.40 GiB sampled peak RSS, including physical verification. It predates this
shared-library benchmark and was not rerun here. Its separate FEC8 preparation
peaked at 32.67 GiB. Neither these small-case comparisons nor that historical
success establishes a full recursive irrep solver or a successful FEC9 solve.
