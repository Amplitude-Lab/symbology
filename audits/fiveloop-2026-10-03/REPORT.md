# Five-loop feasibility and continued kernel optimization

Work on `symrep-solver-dev`, 3 October 2026. The task is the five-loop MHV
**symbol** bootstrap (weight ten), not a function-level amplitude including
constants invisible to symbols.

## Scope and machine

Intel Core Ultra 7 265, 20 logical CPUs, 62.1 GiB RAM, no swap; about 55 GiB
available initially. The existing workflow builds FEC9 and sews it with LEC1.
The immediate prerequisite tested here is FEC8, starting from the previously
certified FEC7. Benchmark outputs are in `output_fiveloop_20261003/`.

The published five-loop MHV computation reports about 24 hours on 160 cores
and 1.3 TB peak memory. Its expanded symbol contains 31,199,389,998 words.
Those are measurements of a different implementation, not lower bounds for
a recursive representation. At an illustrative 16 bytes per stored word,
an expanded list alone would use about 465 GiB. Recursive tensors avoid
storing that list.
[He et al., Appendix B and Section IV](https://arxiv.org/html/2511.09669v2).

## Changes

1. Kernel extraction reserves each modular output row exactly, avoiding
   geometric capacity growth.
2. CRT precomputes the shared inverse of the old modulus modulo the new
   prime. Equal supports update in place; changed supports merge their union.
3. CRT and rational reconstruction run together by row. Successful rows
   retain their rational representation and release their big-integer CRT
   storage. If another prime is needed, these residues are decoded exactly
   one row at a time. Different rows can finish at different primes. A failed
   final certificate also preserves the state needed to continue.
4. Factored modular equation assembly uses parallel, deterministic tiles.
   The tile budget bounds the number of terms before cancellation. Callbacks
   run serially after the worker barrier, so a callback can safely use the
   same pool for elimination.
5. File-producing extensions serialize the certified matrix directly,
   avoiding a duplicate rational CSR tensor. WXF encoding and exact byte
   readback use bounded chunks. The byte format is independently compared
   with the existing encoder and parser, including rational and big-integer
   coefficients, index widths, empty arrays, chunk boundaries and corruption.

All changes apply to general rational kernels or factored tensor equations.
There are no polygon dimensions, expected nullities or known amplitude
coefficients in these optimizations. The full exact residual/rank certificate
is retained.

## Rejected default: elimination in batches

`streaming_echelon.hpp` is an experimental finite-field backend. It reduces
every equation against an independent triangular row basis, compacts each
input batch with SparseRREF, and retains only independent rows. A heap visits
the active pivots in triangular order. This is an exact algorithm, not
uncertified equation sampling.

The first controlled FEC7 comparison used eight workers and the same compiled
probe, input FEC6, dependency, I/O and reconstruction revision:

| Backend | Wall time | Sampled peak RSS | Output nonzeros |
| --- | ---: | ---: | ---: |
| Complete finite-field equation matrix | 14.77 s | 2.36 GiB | 20,515,434 |
| Batches of at most 8,192 equations | 18.18 s | 2.18 GiB | 30,564,594 |

The batch ordering loses some global pivot choices. The resulting basis is
49% denser and would make subsequent recurrences more expensive. This backend
is therefore **not enabled by default**. Its optional argument is exercised
by exact general tests; the benchmark probe can select it explicitly.

Evidence: `fec7-current.{json,time,log}`, `fec7-batched.{json,time,log}`.

## Weight-eight measurements

The complete finite-field core has 957,631 nonzero rows, 102,389 columns and
659,216,640 nonzeros. Its rank is 84,918, giving nullity 17,471. Elimination
temporarily grows well beyond the input nonzero count. The certified
primitive integral output contains 230,509,403 nonzeros.

The original default, capped at 36 GiB virtual memory, terminated with signal
11 after 598.93 seconds, just after its two-prime exact certificate. Peak RSS
was approximately 32.91 GiB. It did not produce an output basis. This is a
bounded-run failure; it does not establish that the original program could
not finish if given the whole machine. The library uses unchecked raw tensor
allocations, and the process was near its virtual-memory ceiling. Allocation
exhaustion is consistent with this failure, but no captured stack trace proves
the precise fault site. The failure occurred after certification and before
any output basis was written.

An intermediate revision removed the CSR copy but still held large CRT
storage and read the entire WXF back into memory. It certified the kernel,
then failed with `std::bad_alloc` during output at the same 36 GiB cap.
This led to row-wise reconstruction and bounded output verification.

Evidence: `fec8-baseline.{json,time,log}` and
`fec8-engineering.{json,time,log}`. The later revision uses a 44 GiB virtual
cap to leave physical memory for the desktop and avoid repeatedly losing a
certified result solely to the earlier arbitrary cap.

The row-wise reconstruction and bounded-I/O revision **completed** in 629.34
seconds (10.49 minutes, contended), with sampled peak RSS **32.67 GiB**. The
1.6 GB `FEC_8.wxf` passed exact byte readback. Its two-prime modulus has 121
bits, exceeding the rigorous 69-bit residual bound. Together with modular
rank 84,918 and the private kernel coordinates, this certifies the full
17,471-dimensional rational kernel. This is a completed weight-eight FEC
basis, **not** a five-loop amplitude. Evidence: `fec8-final.{json,time,log}`.
The timing includes output and verification. The final scheduling refinement
balances tile tasks by estimated terms instead of row count, and an
additional first-prime early exit avoids reconstructing later rows once a
further prime is already required; the large run preceded those two refinements.

## Weight-nine preflight and the five-loop conclusion

Using the new FEC8 file, exact short-relation presolve reduces 733,782
coefficient variables to 266,452. A complete generation/count pass at prime
1,152,921,504,606,847,009 produces **2,769,885 nonzero equations and
7,386,376,131 nonzeros**, after merging and cancelling entries. It keeps
only a bounded tile, not the complete equation matrix.

At eight bytes per modular value and four bytes per column index, the matrix
payload alone is **88,636,513,572 bytes = 82.55 GiB**. The actual present
SparseRREF backend additionally needs row storage, input factors, transpose
indices and elimination fill. Consequently this route cannot fit in the
machine's 62.1 GiB RAM. This conclusion comes from an actual generated matrix
count, not extrapolation from four loops. The count pass took 471.24 seconds
under contention and peaked at 11.85 GiB. Evidence:
`fec9-count.{json,time,log}`.

This does **not** prove that five-loop bootstrap is mathematically impossible
on this machine. It rules out materializing this FEC9 equation matrix with
the current backend. The next substantial change should avoid constructing
the unrestricted FEC9 basis: impose both boundaries and symmetry on the
sewn coefficient space first, or retain an implicit/matrix-free operator
with a certified small-kernel method. Nontrivial FEC and LEC representations
must still be paired into invariants; discarding all nontrivial FEC sectors
would lose valid solutions. Neither alternative has been demonstrated here,
so no completed five-loop result or promised runtime is claimed.

## Timing caveat

At 08:34:23 a separate 16-thread CPU workload started on this machine. It
uses roughly fourteen to fifteen cores and was left untouched. The initial
baseline and FEC7 pair above preceded it. Later elapsed times must be marked
as contended; they are not clean speed comparisons with those earlier runs.
Per-process RSS and exact correctness remain measured independently.

## Accepted revision: NMHV comparisons

Three fresh runs per revision and case, eight workers, rotating before/after
order. Both revisions ran while the independent CPU workload remained active.
These are paired measurements under contention, not isolated machine timings.
The before executable was preserved before rebuilding. The driver records
both executable hashes and input hashes in `nmhv-comparison/results.json`.

| NMHV chain | Previous optimized kernel, median | New kernel, median | Previous peak RSS, median | New peak RSS, median |
| --- | ---: | ---: | ---: | ---: |
| Hexagon through weight 5 | 1.423 s | 1.227 s | 124.1 MiB | 84.2 MiB |
| Heptagon through weight 4 | 26.768 s | 24.255 s | 4.674 GiB | 3.645 GiB |

All 27 before/after WXF pairs (all weights, both cases, all repeats) are
byte-identical, including the bases previously compared exactly against the
unpacked hexagon package and the ordinary heptagon solver. Thus the new
revision preserves those benchmark spaces without another expensive basis
change calculation. This improves on the already-promoted streamed kernel,
not merely on the older ordinary solver.

The reproducible driver now accepts `--before-program` to compare two
streamed-kernel executables. Evidence: `nmhv-comparison.log`,
`nmhv-comparison/results.json`, `nmhv-comparison/exact-byte-comparison.json`,
and the per-run logs/RSS records in that directory.

## Validation of the final default

* 160 independent rational rank/residual/kernel-rank cases, including bad
  primes, vanishing modular supports, prime denominators, large multi-prime
  coefficients, mixed reconstruction heights and the optional batch backend.
  A further test deliberately rejects two correct candidates to exercise
  exact rational-to-CRT recovery across a multiword modulus.
* 102 factored tensor cases, with independent exact assembly for forward and
  backward extension and multiple-terminal sewing, plus equality of serial
  and parallel equation streams.
* Recursive private-word actions and restricted projections agree with
  explicit tensor contractions. WXF byte encoding agrees with SparseRREF's
  encoder and parser for all supported integer-index widths, rational and
  large-integer values, empty outputs and multiple chunks. Corruption is
  rejected.
* Public native, NMHV and symmetry CLI regressions all pass, including
  recursive actions, continuation, basis comparisons, lower-loop historical
  tensors and published NMHV invariant dimensions.
* A fresh full four-loop MHV workflow completed in **71.68 seconds** under
  contention, with **5.41 GiB sampled combined RSS**. All **67 WXF outputs**
  are byte-identical to the previously accepted default. All 12 independent
  published word-coefficient checks at two, three and four loops pass.
* The pinned SparseRREF revision and its seven existing patches pass the
  dependency check. No dependency-version change was made in this iteration.

Evidence is in `structured-accepted.log`, `tensor-accepted.log`,
`word-chart-accepted.log`, `native-regression.log`, `nmhv-regression.log`,
`symrep-regression.log`, `fourloop.{json,time,log}`,
`fourloop-compare.json`, `fourloop-words.log` and `provenance.json` under the
output directory. The private project-data regression remains unavailable,
as documented in the earlier audit; it is not claimed here.

The accepted memory/reconstruction, scheduling and output changes are now
used by the default streamed kernel. Experimental batch elimination remains
opt-in through the benchmark/API argument because its recursive basis cost
was worse in the first controlled comparison.

## Reproduction and preflight

`bench/streaming_kernel_benchmark.cpp` accepts:

```text
forward|backward|sew|profile|count|count-sew CONDITION SEED LAST_OR_DASH OUTPUT BATCH_ROWS METHOD THREADS
```

`BATCH_ROWS=0` uses the normal complete-core backend. Positive values enable
experimental elimination in batches. `profile` performs exact short-relation
presolve and counts terms before merging, without allocating the full matrix.
`count` additionally generates every normalized equation at the first modular
prime and counts its actual nonzeros, discarding each bounded tile. Its
reported payload is the coefficient/index storage alone, excluding row
metadata, factors, transpose structures and elimination fill. `count-sew`
does the same for a supplied right recurrence.

The preflight count is not a solve or a certificate for a new weight. Its
purpose is to distinguish the memory of a factored operator from the memory
needed by the current materialized finite-field elimination backend.
