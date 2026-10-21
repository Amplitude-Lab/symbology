# Nine-loop hexagon symbol with staged physical constraints

The calculation is complete. Early collinear, MRK and OPE restrictions, interleaved with sewing, determine a **unique rational nine-loop MHV hexagon symbol at weight 18 with split 14 + 4**, within the supplied parity × D3, FEC/LEC and physical ansatz. The selected production run peaked at **10.790 GiB (11.59 GB)**, below the 30 GiB RSS guard. The conventional left-precontracted assembly crossed that guard before completing its first sewing solve.

The portable [Wolfram release](/home/ana/Documents/symbology/output_hexagon_14p4_20261006/release_MHV_w18_14p4/README.md) is complete. It contains 19 checked WXF tensors, a loader and checksums. The compact rational object is in `/home/ana/Documents/symbology/output_hexagon_14p4_20261006/cached-nine/joint/`; its completion receipt, input identity and all seven native prime caches are retained.

## The affine schedule

Physical and integrability constraints can be imposed in either order when their complete intersection is retained. Let V be the symmetry-invariant product ansatz, M the complete sewing operator, and P the physical projection with right-hand side b. The target is {c in V : Mc = 0, Pc = b}. Maintain c = o + Nα at every stage. A sewing block gives MᵢNα = −Mᵢo; a physical block gives PⱼNα = bⱼ − Pⱼo. Solve for a particular solution and all null directions, then multiply both back into V. Every previously imposed equation survives.

The fixed CG convention builds V from every matching irrep pair, including nontrivial pairs and odd × odd contributions. Thus every updated origin and direction retains parity and D3 symmetry. Projecting the two factors separately onto their trivial representations would discard allowed invariant tensors; this calculation retains them.

Some physical maps require an integrable input for their usual interpretation. The implemented maps extend linearly to the partial candidate space, using the original ordered-integration prescription. Their values outside ker(M) can depend on that extension. The final intersection is valid because the prescription agrees on ker(M), every remaining sewing equation is imposed, and fresh physical checks are performed afterward. The six- and eight-loop controls also match the conventional ordering exactly modulo the control prime.

The economical schedule first discovers a sparse sewing core, applies physical restrictions in the resulting smaller affine space, and then contracts every omitted sewing equation. It reconstructs only the final rational physical object, with no presentation RREF or rational bulk kernel.

| Completed stage | Free directions |
| --- | ---: |
| Partial sewing discovery | 6,000 |
| Collinear | 1,555 |
| MRK | 541 |
| OPE | 501 |
| Remaining complete sewing operator | 0 |

The 6,000-direction target is a scheduling parameter. It does not prescribe the physical answer. The 501 directions after OPE violate integrability; the remaining sewing equations eliminate all of them. This is why the physical restriction of a partial candidate cannot be accepted as a completed symbol.

## Exact result and validation

The invariant coordinate space has dimension 195,155. The unique compact sewing representative has 184,573 nonzero rational coefficients and coefficient height 136 bits. Reconstruction succeeded with five CRT primes (301-bit modulus) and two independent validation primes. Smaller reconstruction budgets failed and were extended automatically, reusing completed native caches.

Each native solve passed a full, unsampled sweep of all 461,645 left-factored sewing rows. The independent certificate contracted the QQ representative directly through the original FEC, LEC and CG interface, without using that row reduction. All **1,541,565 original interface equations** vanish at all seven native primes, and the rational representative agrees with every native cache. After clearing denominators, the absolute integer residual bound has 286 bits, below the 421-bit certificate modulus. This proves that every original rational sewing residual is exactly zero. No residual-only extra primes were needed.

Every native solve also passed fresh collinear and MRK checks using 32 polarized samples for each map, and all 1,122 held-out OPE coefficients, including the held-out T¹ and higher-S coverage. These are the implemented finite-field physical checks; the report does not promote sampled limits into an exhaustive functional identity proof.

After reconstruction, an independent exact QQ checker passed all **37 independent pure-log origin count partitions**. This condition was not used in the fit. The expected coefficients follow from the ordinary-symbol origin limit in equation (42) of [Dixon and Liu, *An Eight Loop Amplitude via Antipodal Duality*](https://arxiv.org/html/2308.08199): E tends to exp(g²E¹), with one-loop symbol minus the sum of the three diagonal log-letter pairs. For even counts nᵢ, the ordered coefficient is (−1)ᴸ times the product of (nᵢ−1)!!; any odd count gives zero. At nine loops the all-u word gives −34,459,425 and counts (6,6,6) give −3,375, both matched exactly. Integrability and D3 symmetry make these 37 partitions sufficient for this full homogeneous pure-log origin polynomial. The checker first passed the known one-, two- and eight-loop controls (2, 4 and 30 cases).

The supplied MHV LEC chain uses its weight-one single-final-entry seed, integrability and cluster adjacency. Its dimensions 6, 25, 83 and 249 are the generic suffix spaces. The empirically observed higher-final-entry bonus relations were not added to select this result. The computation determines the ordinary symbol; function-level constants invisible to it are outside its scope.

## Resources and controls

FEC13 and FEC14 were generated from the retained validated chain. Their dimensions are 5,431 and 9,014; FEC14 has 8,064,127 compact nonzeros and 120-bit coefficient height. The extension took 26.42 seconds with about 0.98 GiB sampled peak RSS. LEC4 has dimension 249.

The conventional 14+4 trial reached 30.018 GiB sampled RSS and was stopped by the guard after 26.94 seconds during assembly. This is a memory failure, not a completed conventional solve or a total-runtime comparison.

The selected first native nine-loop prime used eight threads: 2,728.63 seconds (45.48 minutes), with 9.996 GiB sampled peak RSS. The six additional native primes used sixteen threads, sequentially. First-prime production plus reconstruction took 14984.98 seconds (4.162 hours), including a recorded 143.38-second suspension for sequential validator compilation and control checks. FEC extension, exploratory controls and export are recorded separately. Only one large native prime solve ran at a time.

At eight loops, the following single-prime controls all reproduced the known rational symbol modulo p = 1152921504606847009, with every sewing row and the independent physical checks passing:

| Native eight-loop control | Wall time | Sampled peak RSS |
| --- | ---: | ---: |
| Right contraction | 326.39 s | 1.322 GiB |
| Left contraction | 205.18 s | 1.360 GiB |
| Left contraction with bounded OPE cache | 176.84 s | 1.566 GiB |

The bounded OPE head cache stores projected input tensors, keyed by input SHA, endpoint/truncation, stride and prime. Three nine-loop entries total 3,045,794,544 bytes, below its 4 GiB limit. It is cleared before each new prime solve and reused in held-out checks; it contains no solution coefficients.

## Use the result

The physical sewing matrix is 9,014 × 249, with 307,781 nonzeros. All 19 WXF files were written, read back and checked for exact dimensions, indices and coefficients; all release checksums, including its manifest, were verified. The release is approximately 779 MiB and has no runtime dependency on SparseRREF or the build environment.

Load the helper in Wolfram Language:

```wolfram
Get["/home/ana/Documents/symbology/output_hexagon_14p4_20261006/release_MHV_w18_14p4/hexagon_symbols.wl"];
Dimensions[MHVSEW[9]]
```

Use `MHVcoeff[letters...]` with exactly 18 integer letters, numbered 1 through 9 in the alphabet (a,b,c,m_u,m_v,m_w,y_u,y_v,y_w). The definitions use a = u/(v w) and cyclic variants, not square roots of these ratios. Coefficients are obtained by contracting the saved chains; expansion of all symbol words is unnecessary.

The [reproduction instructions](/home/ana/Documents/symbology/audits/hexagon-nine-loop-2026-10-06/REPRODUCE.md), [machine-readable result](/home/ana/Documents/symbology/audits/hexagon-nine-loop-2026-10-06/result.json), [prime inventory](/home/ana/Documents/symbology/audits/hexagon-nine-loop-2026-10-06/completed-native-primes.json), source/executable hashes, completion receipts, resource reports and complete validation logs accompany this report. The driver is [hexagon_joint_physics.cpp](/home/ana/Documents/symbology/bench/hexagon_joint_physics.cpp); independent validators are [hexagon_joint_certificate.cpp](/home/ana/Documents/symbology/bench/hexagon_joint_certificate.cpp) and [hexagon_origin_check.cpp](/home/ana/Documents/symbology/bench/hexagon_origin_check.cpp). Completed stage caches and expensive input data are preserved in the experiment directory.
