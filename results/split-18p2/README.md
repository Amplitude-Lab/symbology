# Ten-loop 18+2 split experiment (2026-10-21/22)

Goal: bootstrap the ten-loop MHV hexagon with the split 18+2 (FEC head chain of
18 entries x LEC tail of 2 entries), with the Dixon-Liu double-final-entry
bonus relations imposed on the two-entry tail, and compare time and memory
against the preserved 17+3 first-prime receipts. No pre-existing file was
modified; all artifacts live in this folder. The working tree was committed
before the experiment started (commit 7675914 on 4pformfactor).

## Why this needed new engineering

The 18+2 split needs the FEC weight-18 chain, which did not exist. The
reference `irr_fec` materializes the assembled recurrence matrix over QQ and
extrapolates to ~90-130 GB peak RSS at weight 18 (measured: 13.7 GB at the
weight-16 step, 44.4 GB at weight 17), which does not fit this 62 GB machine.

`src/irr_fec_stream.cpp` is a streaming variant of the reference generator:
per CRT prime it scatters the previous chain through the template routes
directly into modular residues (no QQ assembly), reduces with the unchanged
SparseRREF kernels, saves byte-compatible k_crt caches, streams the caches
per bounded free-column slice to CRT-reconstruct kernel columns, and releases
the chain (loaded from the exported tensors per prime) during reductions and
the matrix rows during the WXF export. Natural row order everywhere keeps
cached pivots valid across process restarts, which the repeated memory-guard
restarts below relied on.

## Validation of the streaming variant

`FEC-validation/` reuses the existing w01..w16 (hardlinked tensors) and
regenerates w17 from scratch:

- all six `irrFEC_17_*.wxf` blocks byte-identical to
  `output_hexagon_16p4_20261006/FEC/w17` (dimension 38472, nnz 171,950,436,
  height 343);
- measurement receipt `audit/fec17-stream-validation3.json`:
  2,857 s wall, 39.8 GB peak RSS, 4,578 CPU-s user (16 threads), versus the
  reference `output_hexagon_17p3_20261007/audit/fec17.json`:
  1,596.6 s wall, 44.4 GB peak RSS (8 threads).

The variant trades ~1.8x wall time at weight 17 for a peak set by the modular
matrix instead of the QQ assembly, which is what makes weight 18 feasible.

## FEC weight-18 chain

`FEC/` reuses w01..w17 (hardlinked) and carries fresh w18 artifacts.
Final process receipt `audit/fec18-generation8.json` (caches from earlier
attempts were reused across restarts): exit 0, 14,434 s wall, 54.2 GB peak,
19.8 CPU-h. Restart receipts `audit/fec18-generation*.json` document the
memory-guard stops that were resolved by freeing slice buffers, per-prime
chain loading and releasing matrix rows during export.

Result (`FEC/w18`, summary row 18): dimension 61,165, nnz 363,613,184,
height 508 bits; multiplicities pe_A1 8778, pe_A2 7076, pe_E 15838,
po_A1 2299, po_A2 2256, po_E 4540.

## First native prime, heap mode: does not fit

`run_first_prime_18p2.sh` (V27D profile, 16 threads, same 50 GB decimal RSS
guard as the 17+3 receipts): stopped by the guard during chain loading at
854 s, peak 50.03 GB, before discovery (`audit/solve-18p2-first-prime.json`,
`rss_limit_exceeded: true`). The 363M-nnz FEC-18 chain alone pushes the
solver's heap baseline past the guard.

## First native prime, mapped mode

`run_first_prime_18p2_mapped.sh`: identical configuration except the
file-backed working storage of the validated replay evidence
(HEXAGON_MAPPED_FILE, GIB=192, SOFT_GIB=40), same 50 GB decimal guard.
Measured geometry (from the log):

- FACTORED_SEAM rows=477,136 columns=144,151 (17+3: 989,800 seam rows,
  275,246 columns);
- CONDITIONAL_MHV_DOUBLE_BONUS on the two-entry tail: 25,213 equations,
  rank gain 22,381, ansatz 144,151 -> 118,938 conditional columns
  (the audited pair-tail copy restriction 4,0,4,2,1,3 applied to the
  w18 multiplicities);
- discovery alone: 6,207 s wall / 44,731 CPU-s (12.4 CPU-h), already half
  again the 17+3 whole-prime CPU cost.

Receipts: `audit/solve-18p2-first-prime-mapped.json` (final numbers to be
appended when the prime completes), log `log/solve-18p2-first-prime-mapped.log`.

## Comparison against 17+3 (first native prime, V27D/16)

Receipt: `audit/solve-18p2-first-prime-mapped.json` — exit 0,
25,850.6 s wall (7.18 h), 115,964 s user + 6,782 s system (34.10 CPU-h),
48.67 GB peak sampled RSS under the shared 50 GB decimal guard. All checks
passed: V18_SEAM_VERIFIED full_operator=1 rows=477,136;
JOINT_PHYSICAL_CHECK collinear/mrk/ope_small_s each passed=1 remaining=0;
JOINT_PRIME_COMPLETE prime=1152921504606847009 directions=0 columns=144151.

| Quantity | 17+3 (V27D receipt) | 18+2 heap | 18+2 mapped (measured) |
| --- | --- | --- | --- |
| storage | heap | heap | file-backed |
| seam rows | 989,800 | - | 477,136 |
| solve columns | 275,246 (unrestricted) | - | 144,151 original / 118,938 after in-solve pair-tail bonus |
| wall seconds | 6,599.9 | killed at 854 | 25,850.6 (3.92x) |
| CPU hours | 23.38 | - | 34.10 (1.46x) |
| peak RSS GB | 46.19 (sampled) | > 50.03 at load | 48.67 |
| FEC chain prerequisite | existing w17 | w18 not generatable by reference irr_fec in 62 GB | w18 via streaming variant: 4.0 h wall / 19.8 CPU-h / 54.2 GB (final process; restarts reusing caches) |

Reading: the smaller 18+2 ansatz (2.3x fewer columns, 2.1x fewer seam rows)
does not pay for the 2.1x heavier FEC-18 chain. Every chain projection
(collinear/MRK/OPE fitting through 18 head weights) dominates, discovery
alone cost 12.4 CPU-h, and the heap mode does not even fit the shared 50 GB
envelope. On this machine and the current V27D solver, 17+3 remains the
cheaper split per native prime; 18+2 additionally required a one-time FEC-18
generation investment. CPU-hours (the scheduling criterion) rise 45.9% and
wall time 292%, the latter inflated by mapped-storage page traffic
(6,782 s of system time).

Caveats (corrected 2026-10-22 after re-reading the 17+3 log): the 17+3 V27D
first prime ALSO imposed the bonus in-solve (CONDITIONAL_MHV_DOUBLE_BONUS
equations=58,223, rank gain 47,440, 275,246 -> 217,024 conditional columns;
JOINT_PRIME_COMPLETE's "columns" field reports the original product on both
splits), so the two runs are directly comparable: both bonus-restricted, same
driver, threads, prime and guard. 18+2 heap exceeds the 50 GB decimal
envelope at chain load; 17+3 ran in heap mode, 18+2 in the validated
file-backed mode, so wall times are not storage-mode comparable (CPU-hours
are the fairer column).

Phase breakdown of the first native prime, from the PROFILE lines of both
logs (seconds; 17+3 = output_hexagon_ten_v27d_20261009 p01 receipt):

| phase | 17+3 | 18+2 |
| --- | ---: | ---: |
| bind/preflight | 23.0 | 90.9 |
| compute_discovery | 3,691.3 | 6,207.5 |
| compute_early_collinear | 1,266.4 | 13,592.5 |
| compute_early_mrk | 125.6 | 630.6 |
| compute_early_ope_small_s | 809.4 | 2,174.1 |
| compute_seam_complete | 42.5 | 15.5 |
| checks (col+mrk+ope) | 502.8 | 2,577.7 |

The 18+2 loss is localized in the collinear operator fitting (10.7x), i.e.
projection through the 18-entry head chain; the seam itself is cheaper
(15.5 s vs 42.5 s) because the ansatz is smaller. The adaptive scheduler
added an operator in the 18+2 run that 17+3 never needed
(ope_lowgrade4, 1,182 s; operator cost sums 17+3: mrk 338 s +
collinear_tail_exact 112 s).


## End-to-end pipeline comparison (added 2026-10-22, productive time only)

Path A = 17+3 with the existing w17 chain (V27D p01 receipt). Path B = build
FEC-18 from w17 (streaming generator; single-run equivalent = sum of
productive receipts 2,347.7 + 1,958.1 + 3,824.3 + 14,434.3 s; the 10,636 s
lost to memory-guard restarts and the export bug are excluded) plus the 18+2
first prime. "Free" = undetermined directions remaining.

| Step | A: 17+3 wall (CPU-h) -> free | B: 17->18 + 18+2 wall (CPU-h) -> free |
| --- | --- | --- |
| generate FEC head chain | / (w17 exists) | 22,564 s (31.2) -> head dim 61,165; ansatz 144,151 |
| preflight + bind | 23 s (0.02) -> 275,246 | 91 s (0.03) -> 144,151 |
| seam discovery | 3,691 s (13.4) -> 218,001 | 6,207 s (12.4) -> 111,525 |
| bonus seed (tail) | ~0 s (58,222 eq, +47,440) -> 170,555 | ~0 s (25,213 eq, +22,381) -> 89,120 |
| early physics COL | 1,266 s (4.9) | 13,593 s (5.0, near-serial) |
| early physics MRK | 126 s (0.5) | 631 s (1.3) |
| early physics OPE | 809 s (2.6) | 2,174 s (7.2) -> 9,806 |
| final seam replay | 43 s (0.1): 2,877 -> 0 | 16 s (0.03): 710 -> 0 |
| verification checks | 503 s (1.9) -> 0 | 2,578 s (8.0) -> 0 |
| total | 6,540 s / 23.33 CPU-h | 48,415 s / 65.31 CPU-h |

Marginal cost of every later prime: 23.3 CPU-h (17+3) vs 34.1 CPU-h (18+2),
so 18+2 loses both including and excluding the one-time chain build.
