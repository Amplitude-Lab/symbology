# Staged ten-loop hexagon calculation: working method

The target is the weight-20 MHV symbol at FEC16 × LEC4. This note records the method while the native computation runs; it does not assert a completed ten-loop result.

## Why the constraints can be interleaved

Write the current candidate as c = o + Nα in the invariant tensor space V. A block A c = b becomes (A N)α = b − A o. If α = a + Kβ solves that block, the updated candidate is c = (o + N a) + (N K)β. Thus both the particular solution and every surviving direction retain every earlier equation. Integrability uses b = 0; physical restrictions can have a nonzero b. The origin must therefore be present before the first physical restriction.

Let M be the complete sewing map. If a physical map P is defined on ker(M), an implemented linear extension P̃ to V is admissible when P̃ agrees with P on ker(M). The final intersection ker(M) ∩ P̃⁻¹(b) is then independent of the extension. Values outside the integrable locus do not define physical symbols, so all omitted sewing equations must still be imposed. The driver uses the original ordered projection prescription and tests the final integrable object with held-out physical data.

The parity-even D3 invariant tensor space includes every matching irrep pair. The multiplicities and coordinate counts at the requested split are:

| Irrep sector | FEC16 copies | LEC4 copies | Invariant coefficients |
| --- | ---: | ---: | ---: |
| pe_A1 | 3,445 | 28 | 96,460 |
| pe_A2 | 2,771 | 16 | 44,336 |
| pe_E | 6,213 | 44 | 273,372 |
| po_A1 | 902 | 25 | 22,550 |
| po_A2 | 873 | 14 | 12,222 |
| po_E | 1,772 | 39 | 69,108 |
| Total | | | 518,048 |

This preserves allowed A2 × A2, E × E and odd × odd contributions. Every affine update stays in this space and uses one shared solution ordering across all six blocks. The fixed CG and physical-copy-primitive-v1 conventions are retained; E components are not independently rescaled.

## Discovery before the physical restrictions

The FEC chain already enforces its internal integrability, extended Steinmann and one-term adjacency conditions. The generic MHV LEC chain enforces its internal integrability and full adjacency from the fixed single-final-entry seed. The joining interface uses the full-CA tensor. The requested invariant product has 518,048 coordinates; the left-factored joining operator has 1,232,602 rows. Its full product matrix is not stored.

The operator supplies rows in bounded tiles. Cheap sparse rows seed a finite-field triangular row basis; denser rows are admitted progressively. Exact singleton and doubleton substitutions prune redundant work. Deterministic nested equation samples guide discovery but never decide acceptance. The scheduling target is 12,000 directions, rather than an assumed physical dimension.

Once a partial candidate is available, the driver applies collinear, MRK and covered OPE equations. It then contracts every remaining sewing row into the smaller affine space, restricting and replaying until an entire unsampled sweep passes. Final physical checks use points distinct from the training stream and endpoint coefficients held out of the OPE fit.

## Keeping the computation within memory

The first native probe used sparse checkpoints and released temporary rational recurrences after constructing the factored operator. It was superseded during discovery to remove another avoidable candidate copy; no completed ten-loop native prime is claimed from that probe. Its logs and resource report remain available.

The left-contracted probe reached 37.9 GB at rank 245,424 while substantial discovery remained. The selected trial therefore compares contraction of the much smaller LEC4 factor; the mathematical operator and physical conditions are unchanged. This direction choice is a resource experiment, not a completed ten-loop result.

The V5 driver constructs the kernel while consuming the RREF pivot rows. Its free-coordinate convention agrees with the original implementation: each private free coordinate is −1 and pivot coordinates are the corresponding RREF coefficients. Storage for trivial field elements is created without eagerly clearing a second full array, and every entry is assigned before it is read. Each source pivot row is freed after its contribution is transferred.

Physical projection consumes the packed candidate row by row into the six irrep blocks. If C_b is an old block and T is the affine restriction matrix, the new block is T C_b. The same T acts on every block. Products are formed and packed one block at a time, with old blocks released after use. This avoids retaining a full global candidate and a second block copy during projection and lifting.

Private intermediate checkpoints store the complete sparse affine object, with origin last. Their file size is 80 + 8 r + 16 z bytes for r rows and z nonzeros; the conventional dense payload costs 72 + 8 r n bytes for n coordinates. Each write is streamed back and compared exactly with the live matrix, then protected by SHA256 and a stage-input fingerprint. Final native solve caches retain the package's original format. The probe tests zero and nonzero origins, chunk boundaries, wrong split identities and damaged data.

The run uses sixteen internal threads, one large native prime at a time, a 2 GiB modular recurrence cache and a 12 GiB OPE head cache. The process-tree RSS guard is 50,000,000,000 bytes, expressed as 46.5661287308 GiB, with a 3 GiB system-memory reserve. Completed input blocks, native caches and logs are retained.

## Coefficient ordering and observed resource limits

The V5 right-contracted trial reached rank 271,597 and a 34.9 GB process-tree RSS peak. It was stopped because system MemAvailable fell below the 3 GiB reserve, before the 50 GB RSS limit. It did not produce a completed native ten-loop prime. The complete log and resource history are retained.

Sparse elimination fill also depends on variable ordering. The V6 experiment orders the copy-matrix coefficients by last copy and then first copy within each irrep block. This is a permutation of the same 518,048 coordinates, not a restriction of the ansatz or a change of FEC/LEC normalization. Discovery occurs in that temporary order. Every candidate coefficient is then mapped back to the original first-copy/last-copy layout before physical projection, checkpointing, complete residual sewing and rational reconstruction. The permutation and inverse are checked explicitly; the established symbol provides a finite-field control. For a block with m first copies and n last copies, the temporary index map is f n + l → l m + f. Equivalently, if c′ = Πc, discovery uses MΠ⁻¹ and returns candidates to c = Π⁻¹c′ before every physical operation. This relabels variables while preserving the full constraint intersection.

## Physical data and independent checks

The ten-loop main and extended MHV OPE libraries both have complete, hash-bound generation receipts. The selected main library uses all its covered T1/T2 solve rows: 1,824 equations, with 456 higher-power coefficients held out. The maximal-helicity extraction cancels the unknown T2 F⁰ contribution; it never fills missing physical data with zero. MRK inputs are available through ten loops.

The known eight-loop symbol was reproduced exactly modulo p = 1152921504606847009 by both V4 and V5, using the ten-loop input library. These are ordering, storage and coefficient-convention controls, not ten-loop results.

Only the final affine product is reconstructed over the rationals, with two independent validation primes. The exact sewing certificate then contracts the reconstructed object directly through the original FEC/LEC/CG interface. Clear all denominators to obtain integer residuals Z with an explicit bound |Z| ≤ B. If every residual is zero modulo distinct good primes whose product exceeds B, the residuals vanish exactly over the rationals. This check is independent of the factored row reduction used in discovery.

An independent origin checker is prepared for all 44 symmetric pure-log count partitions at weight 20. It is not used in the initial fit. At ordinary symbol level the all-order origin limit gives exp(g² E1), where E1 has symbol minus the sum of diagonal log-letter pairs. Even counts n_i have coefficient (−1)^L times the product of (n_i−1)!!; an odd count gives zero. At ten loops the (20,0,0) and (8,6,6) cases predict 654,729,075 and 23,625. The all-order prediction is given in equation (42) of [Dixon and Liu, An Eight Loop Amplitude via Antipodal Duality](https://arxiv.org/html/2308.08199).

The final technical report will add the observed rank chain, native primes, rational reconstruction height, exact certificate bounds, independent checks, timings and memory measurements. A unique ordinary symbol does not determine function-level constants invisible to it.

## File-backed temporary arrays (V7)

The V6 trial stopped at a 43.1 GB peak because the system-memory reserve was exhausted. Its last reported rank was 296,591, with 2,125,869,115 retained nonzeros. It had not completed a native prime. Thus changing coefficient order reduced some fill and discovery time but did not remove the storage bottleneck.

V7 preserves the V6 numerical algorithm and maps SparseRREF arrays of at least 4,096 bytes into a private, fresh scratch file. Small arrays and GMP/FLINT limbs use their original allocator. The mapping reserves virtual addresses without reserving the corresponding RAM. Exact array bytes persist in shared file pages. At an RSS threshold, the allocator synchronizes those pages and releases their residency; subsequent accesses reload the same coefficients. This is temporary working storage, not an imported result or a completed-prime cache. The scratch filename is excluded from numerical stage identity, while the storage mode and executable hash are included.

The production mapping reserves 192 GiB of sparse-file address space under a 256 GiB virtual-address guard. The actual process-tree RSS guard remains exactly 50 GB, with a 3 GiB system reserve. The tentative residency threshold is 16 GiB. These settings do not make all working storage free: ordinary allocations, input caches, file-page faults, and output construction still consume memory and require monitoring.

The allocator probe tests expansion, shrinking, transition from heap to mapped storage, concurrent buffers, forced eviction and rational arrays with numerator/denominator values beyond the small-integer representation. It passes. The eight-loop control deliberately uses a 1 GiB residency threshold to exercise eviction during the full calculation; it reproduces the established eight-loop coefficient payload exactly modulo p = 1152921504606847009, with all sewing and physical checks passing. The control completed 28 evictions in 148.675 seconds, peaking at 2.203 GB process-tree RSS.

The first V7 ten-loop trial used a 16 GiB eviction trigger. It remained within the memory guard but page turnover increased: by 1,588 seconds its process counters recorded about 56.5 GB read and 238.5 GB written, while the scratch file occupied about 29.6 GB. The trial was manually superseded after 1,670.57 seconds with a 26.15 GB peak, before any completed native prime. Its last reported rank and counters are retained separately. The selected retry uses the same frozen executable and numerical algorithm, with a 32 GiB eviction trigger, the same 192 GiB sparse-file capacity, and the same exact 50 GB RSS cap. The eviction trigger is a storage policy, not a constraint or input change.
