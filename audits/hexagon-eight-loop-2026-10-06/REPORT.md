# Eight-loop hexagon with 12+4 sewing — 2026-10-06

The eight-loop **MHV symbol** (weight 16) was calculated at FEC12 × LEC4, uniquely fixed by collinear, MRK and OPE conditions, and exported as self-contained Mathematica-readable tensors. All 25 independent published coefficient checks pass. The final symbol has an exact certificate for every original interface equation. Every successful large run fits the user-authorized 30 GiB RAM ceiling.

**The best tested complete route is the package's direct irrep solver with left contraction, retaining finite-field objects through the physical stages and reconstructing only the unique final result.** It finished in **1172.73 s (19.55 min), peak 5.180 GiB**, starting from cached FEC/LEC prerequisites. Its six QQ origin blocks are byte-identical to the result obtained with the new factored staged sewing solver.

## Why symmetry still works

FEC12 has dimension 3224; MHV LEC4 has dimension 249. Their full product has 802776 coordinates. The complete parity-even D3 singlet embedding has 69988 multiplicity coordinates:

```
482*28 + 371*16 + 850*44 + 120*25 + 107*14 + 222*39 = 69988.
```

Every matching irrep pair is retained: A1, A2 and E in both parity blocks. Odd × odd pairs contribute to the even symbol. The E singlet uses the fixed metric diag(1,3). Taking invariants separately in the two factors would lose required pairings.

If P embeds these singlets into the full product, the system is M P c = 0, with solution image exactly ker(M) intersected with the invariant subspace. Symmetry therefore acts alongside interface integrability through an exact coordinate restriction. Both factors already satisfy their internal equations. Subsequent physical stages take affine combinations of this space and preserve its invariance.

FEC uses integrability, extended Steinmann and one-term cluster adjacency. LEC and the interface use the archive's full cluster-adjacency tensor. The fixed CG convention and physical-copy-primitive-v1 seed are preserved.

## What transferred from the heptagon strategy

Early product symmetry already exists in the hexagon archive's direct multiplicity ansatz. The new adapter adds bounded, factored row generation and the shared staged discovery/intersection engine. It keeps the six small right-side row bases and the left recurrence, generates rows in tiles of at most 4096 rows and approximately two million unmerged terms, and never stores the full contracted product matrix.

This route is correct, but the hexagon kernel has 1339 dimensions. The five-loop heptagon joint kernel had only four. At the default discovery thresholds, hexagon must finish a full pass, so the small-kernel intersection advantage does not appear. No expected physical nullity was supplied to any solver.

The most useful workflow optimization is **deferring QQ reconstruction until after the physical space becomes unique**. This avoids writing a 31.4-million-entry rational bulk basis and reduces the required bulk prime evaluations from eight to six. The output itself remains exact rational data.

## Matched first-prime contraction comparison

Same executable, input factors, compiler, dependency, four workers and first 60-bit prime; fresh sequential processes. One repetition per method. Input preparation, QQ reconstruction and export are excluded. Concurrent host activity was not stopped, so these are observed workload measurements rather than statistical timing claims.

| Method | Solve seconds | Peak process GiB | Materialized product rows | Materialized product nonzeros | Kernel dimension |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original left contraction (archive default for k=4) | 45.423 | 5.180 | 164213 | 247426838 | 1339 |
| Original right contraction | 98.473 | 7.954 | 246057 | 360791753 | 1339 |
| Factored right contraction with staged solving | 242.378 | 1.381 | 0 | 0 | 1339 |

The factored operator has 246057 implicit reduced rows. Its zero stored-row/nonzero counts mean the product is not materialized, not that equations are absent. All three return 31365569 kernel nonzeros at this prime.

The left contraction is the fastest tested sewing choice. Choosing the smaller factor first does not reliably predict the cheapest final elimination. Factoring lowers modular peak RAM substantially, but is slower in this case; it lowers the memory needed for modular elimination. A complete bulk QQ export still needs the larger memory shown below.

## Complete calculation costs

| Workload | Wall seconds | Peak process GiB | QQ output |
| --- | ---: | ---: | --- |
| Original right sewing, complete QQ basis | 997.04 | 8.093 | 1339 directions |
| New factored staged sewing, complete QQ basis | 2085.04 | 7.585 | Same 1339 directions, byte-identical |
| Physical solve from the QQ basis | 918.14 | 4.569 | Unique origin |
| Direct left sewing → physical stages → final QQ only | 1172.73 | 5.180 | Same unique origin, byte-identical |
| Self-contained release export | 3.01 | 0.484 | 77 MiB release |

Full QQ sewing requires six CRT plus two validation primes, with 361-bit reconstruction modulus. It produces 31365569 nonzeros and maximum coordinate coefficient height 193 bits. Rational writing/readback dominates its final peak RAM: the large modular-memory saving of factoring becomes only about 6.3% saving for a full QQ bulk export, while this run takes about 2.09 times longer than right contraction. This does not support changing the archive's default to staged sewing.

The direct final-only run genuinely recomputed all six bulk prime solves from factors, then applied physical constraints at each prime. It did not load the previously solved amplitude or a completed target-stage cache. It required four CRT plus two independent validation primes, modulus 241 bits. The final unique origin has 65344 nonzeros and coordinate height 91 bits.

The first physical chain has dimensions:

```
1339 sewing directions -> 109 after collinear -> 13 after MRK -> 0 after OPE.
```

The final-only route has the same dimensions. Its intermediate affine coordinate bases differ, but its final QQ origin is exactly equal. OPE uses 600 training and 900 held-out coefficient rows. Collinear/MRK checks use finite-field polarized slot samples; they are distinct from an exact symbolic identity proof of every physical limit. No function-level integration constants are claimed.

## Exact sewing and validation

Both complete QQ kernels pass a separate certificate against all 535680 original interface equations. The rigorous cleared-denominator residual bound is 294 bits; the product of eight fully solved and congruence-checked primes is 481 bits. Private rational coordinates prove independence, and complete finite-field rank bounds prove the full rational nullity. All six sewing WXF blocks have identical SHA256 hashes across the two solvers.

The unique physical origin also passes an exact original-interface certificate: bound 182 bits, verified modulus 361 bits across six complete physical caches, with the certified bulk ancestry and condition identity checked. Its singlet coordinates guarantee symmetry. The final-only origin is byte-identical and therefore has the same certificate.

Exact expanded assembly, QQ residuals, modular independence and complete-kernel checks passed for MHV 4+2, 2+4, 6+2, 4+4 and 8+4. Exact full-space comparisons passed for the first four. The first 8+4 run caught a zero-coordinate mapping bug in the experimental adapter; the failure is preserved and the repaired run passed before the full weight-16 calculation.

An initial staged trial that assembled the rational product hit the old 20 GiB RSS guard after 49.40 seconds during assembly. It produced no kernel. The user then authorized 30 GiB RAM. Subsequent large runs used 30 GiB sampled tree RSS, 40 GiB address space and a 6 GiB host-available-memory floor. Every failure, completed prime and old output is retained.

## Independent published checks

The [Dixon–Liu ancillary page](https://s3df.slac.stanford.edu/people/lance/Cosmic/index.html), associated with [arXiv:2308.08199](https://arxiv.org/abs/2308.08199), supplies `EZMHVg_u11_lin.txt`. It was used only for validation, never to fit physical parameters. The source URL, SHA256 and extracted cases are in the experiment's `published-validation/` directory.

Along v=w=1+epsilon, then epsilon -> 0, y_u ~ (1-u)^2/(u epsilon^2) and y_v=y_w ~ -epsilon/(1-u). Dropping dlog epsilon gives a pullback to {u,1-u}. The checker contracts the saved recurrences and singlet metric with these slot weights, without expanding the symbol. Published HPL words are reversed and receive (-1) for each index 1.

All three one-loop and nine two-loop control words pass exactly. All 25 weight-16 checks pass, including pure-u coefficient 2027025 and alternating-word coefficient 5017620. The weight-16 check took 7.93 seconds. This verifies a restriction of the symbol, not equality of the entire three-variable published amplitude.

## Deliverable and reproduction

The self-contained export is:

```
/home/ana/Documents/symbology/output_hexagon_12p4_20261006/release_MHV_w16_12p4/
```

Load `hexagon_symbols.wl` in Mathematica. It supplies `FCE[w]`, `MHVLCE[w]`, `MHVSEW[8]` and word-coefficient contractions. `hexagon_mhv_sewing_12p4.wxf` is the unique physical sewing matrix. All 17 WXF tensors pass exact write/readback checks; all 25 release inventory checksums pass. The loader needs only built-in SparseArray/WXF support and files within the release directory. No billion-word expansion is needed.

The full experiment, frozen package copies, source inventories, commands, caches and logs are in `output_hexagon_12p4_20261006/`. The original unpacked archive and earlier project modifications are untouched. Completed FEC1..8 were copied with their receipts; FEC9..12 were freshly calculated (3.32 seconds), and MHV LEC1..4 were separately generated. The FEC12 dimension is 3224 and its recurrence has 902542 nonzeros.

`bench/hexagon_staged_build.py` copies an explicit hexagon-v2 package to a fresh output directory, applies `patches/hexagon-v2-staged-sewing.patch`, and builds the sewing comparison, certificates and published-check helpers against the archive's coherent SparseRREF 23ceb2e dependency. It wraps the shared namespace lookup for the archive's physical drivers, without changing numerical algorithms. GCC 14.4.0, -O3 -std=c++20, FLINT/GMP/MPFR/TBB/mimalloc were used.

For the fastest tested complete route, use the archive's existing `irr_pipeline` from FEC12 and MHV LEC4, with `--prereduce left --conditions collinear+mrk+ope_small_s --primes auto`; omit `--input` so bulk is computed within each final prime. The exact measured command is in `mhv-12p4-final-only.json`. Use new output directories when changing strategy or input identity; preserve completed caches.
