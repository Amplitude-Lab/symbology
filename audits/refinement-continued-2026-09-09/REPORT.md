# Continued refinements — 9 September 2026

The second pass implements the next verified corrections from the original
audit: native WXF validation, dependency-aware intermediate caches, corrected
examples, safer standalone execution, installation fixes and narrow-screen
Materials layout. Changes remain uncommitted for review.

## Mathematical changes and evidence

No further changes were made to the shared `linear_solve.hpp` or
`incremental_solve.hpp` algorithms in this pass. The previously verified
contradiction-export correction remains in the working tree.

A separate error was reproduced in `tensor_ops`' symmetry wrapper. It treated
a SparseRREF pivot in an appended transformation column as proof that an
induced map did not exist. An independent rational calculation showed that
the transformed pentagon tensors satisfy the original integrability
conditions exactly. The wrapper now restricts pivots to the unknown columns,
following SparseRREF's own matrix-inverse convention, and checks reduced
residual rows for actual inconsistency. The temporary column preference is
restored even on exceptions. SparseRREF's elimination and reconstruction
implementations are unchanged.

A second reproduced symmetry error represented an empty invariant space as
`a × 0` instead of `0 × a`, causing a contraction failure. Empty results now
produce valid zero-row tensors and can be passed to another symmetry step.
A genuinely non-closed input space is still rejected.

The old pentagon example transformed the wrong axes and attempted to average
an integrability condition tensor with an incomplete group sum. The corrected
example starts with the unrestricted weight-two word basis, solves the
original conditions, then imposes cyclic and flip invariance on the resulting
symbols. The 4p example follows the same construction with its full conditions.
Its precomputed full tensor is checked against the exact concatenation of the
original integrability and extended-Steinmann tensors.

Independent verification establishes:

- **Pentagon:** 76 independent invariant symbols. Every original-condition,
  cyclic and flip residual is exactly zero. The supplied matrices satisfy
  the D5 relations. Constraint rank 885 plus nullity 76 equals 961.
- **4p:** 689 independent invariant symbols, with exact condition and symmetry
  residuals zero. Constraint rank 7960 plus nullity 689 equals 8649.
- The rank certificates use a separate Python modular elimination algorithm.
  Its lower bound on rational rank, together with exact rational residuals
  and independent solution rows, proves completeness over the rationals.
- **E6:** the two-loop template now wires E1 into the current RHS node and
  wires its boundary into the solver. All eight steps succeed; boundary and
  solution match the recorded reference exactly. The same comparison passes
  after exporting and relocating the project without Wolfram.

## Reliability and portability changes

- WXF decoding checks payload lengths, variable-length integers, expression
  completion, array dimensions, CSR row pointers/indices, value counts and
  exact rational denominators before tensor construction. Valid large
  rational values, empty tensors and fitting 64-bit encoded indices work.
  Unsupported compression and value types produce readable errors.
- A fresh dependency build without mimalloc exposed a small allocation leak
  when reading empty tensors. The reader now avoids that zero-capacity
  allocation. Both matrix and tensor readers pass AddressSanitizer,
  UndefinedBehaviorSanitizer and LeakSanitizer checks.
- The WXF changes are shipped as a reproducible dependency patch, integrated
  into `setup-sparserref.sh`. Application and reverse checks passed against a
  fresh copy of the pinned upstream source.
- Automatic FEC, SEW, collinear projection and RHS intermediates have SHA256
  dependency receipts, including executable and output contents. E1 is
  refreshed from its source. Legacy files are recomputed when their receipts
  cannot establish validity. A projection that becomes empty removes its
  obsolete nonempty output.
- Local and exported runs share content-based cache and isolated publication
  code, plus a per-project OS lock. Publication handles removed derived files
  and rolls them back if publication fails. Tests cover changed inputs with
  preserved timestamps, missing/partial output, timeout, relocation, executable
  lookup, Wolfram availability modes and concurrent admission.
- Exported scripts embed a Python standard-library runner. They handle spaces
  and apostrophes and no longer depend on GNU `timeout`. A missing Wolfram
  executable cannot make a missing or changed output silently acceptable.
- Mimalloc can actually be disabled, `make` builds all four required tools,
  `PORTABLE=1` omits host-specific CPU flags, and invocation through PATH finds
  the executable directory. Native subprocess quoting handles apostrophes.
  The Windows launcher now uses the correct virtual-environment path. The
  Unix launcher builds all required tools and uses the package lockfile.
- Materials navigation/cards fit 390-pixel screens; wide property tables stay
  within their card. Desktop and existing save/recovery behavior still pass.
- README, API contract, test instructions and CI configuration reflect these
  contracts. CI now specifies the tested GCC/FLINT dependency versions and
  includes allocator variants, parser diagnostics and exact public examples.
  Remote CI itself was not run in this session.

## Verification results

| Check | Result |
|---|---|
| Shared solver helpers | 90/90 exact systems |
| Contradiction export/reimport | 6/6 cases across both solvers |
| Public two-/three-loop references | All 10 tensors match |
| Native cache recovery | Changed seed/binary, corrupted output and stale empty projection pass; all 10 resumed outputs equal a fresh run |
| SHA256 implementation | 12 independent known-answer/padding-boundary checks |
| Fresh-build WXF memory checks | 27 valid files, 354 invalid files, 1,000 bounded mutations; no sanitizer report |
| Template mathematics | Both complete invariant spaces certified; E6 local/exported reference matches |
| API and standalone orchestration | 49/49 |
| Chromium | 14/14, including 390- and 768-pixel Materials views |
| Production web build, syntax and registry checks | Passed |

The native tests used GCC 14.4 and FLINT 3.2.2 on Linux, with mimalloc disabled.
The isolated dependency environment used for testing is under `/tmp`; the
built executables link against that environment. Rebuild with installed
matching dependencies before moving the repository or removing that runtime.
Python API tests used Python 3.12; browser tests used Chromium on Linux.

On this machine the public three-loop fixture took about **7.2 seconds from
scratch and 0.47 seconds on unchanged resume** in the final cache test. These
are workload-specific measurements, not a claim about high-weight performance.

Evidence: [public native gate](native-public-results.txt),
[cache and timing checks](cache-results.txt),
[template certificates](template-results.txt),
[API/export suite](api-export-results.txt),
[fresh-build memory checks](wxf-sanitizer-results.txt),
[browser checks](browser-results.json), [web build](web-build.txt).
Reproduction: [tests/README.md](../../tests/README.md).

## Applying the updated behavior

Run the dependency setup script before rebuilding. Re-export previously
exported scripts: existing scripts keep their embedded old implementation.
Native intermediates without valid new receipts will be regenerated once.
Existing saved user flows were preserved; create a fresh template project to
obtain the corrected example graphs.

## Remaining work

1. Measure high-weight workloads and replace conservative directory snapshots
   with precise, complete dependency lists. Current snapshots favor correctness
   and require additional temporary disk space and I/O.
2. Run Windows/WSL acceptance checks. The project owner subsequently confirmed
   macOS compatibility, closing that platform follow-up. See the
   [Windows audit](../windows-2026-09-09/REPORT.md) for prepared CI checks and
   remaining gaps. Full phone flow editing and Firefox/WebKit remain outside
   the recorded Chromium/Materials coverage.
3. Run private NMHV/integrability baselines and regenerate the Wolfram-dependent
   datasets. The corrected public examples use their shipped exact tensors;
   this does not certify every symbolic property-generation workflow.
4. Add durable external-process supervision if runs must survive forced server
   or machine crashes. Publication is atomic per file with rollback on ordinary
   errors, not across a power loss. Keep one server worker per project store.

Raw native CLI commands still write directly and do not acquire the frontend
project lock. Newly exported scripts share staging/cache behavior but do not
create server run records or execute server-only completion callbacks.
The WXF reader supports the documented exact CSR subset, not arbitrary WXF;
native encoding currently targets little-endian machines.
