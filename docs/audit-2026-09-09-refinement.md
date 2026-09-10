> **Archived audit report (2026-09-09, first refinements).** Kept as documentation of what was
> audited, fixed, and verified. The evidence files it references (result
> captures, screenshots, probe scripts) are not part of the user-facing
> repository; they live in the development history and on the author's
> machine. Findings described as open here were addressed by the following
> refinement reports and the regression suites in `tests/`.

# Verified first refinements — 9 September 2026

Follow-up: [continued refinements and latest verification](../refinement-continued-2026-09-09/REPORT.md). The results below record the first pass.

This implements the verified first corrections from the robustness audit,
with regression coverage. It does not close every item in the original six
work packages. The mathematical changes were kept small and independently
checked; elimination, pivot selection, reconstruction and null-space
algorithms are unchanged.

## What “exporting inconsistent equations” meant

Yes, this concerns the nonhomogeneous collinear constraint workflow. It was
specifically a loss of information when saving the combined `[M | r]` matrix.
For example, the equations `c0 = 2` and `0 = 1` have no solution. The CLI
correctly rejected them, but exported only `c0 = 2`; reimporting that file then
reported a solution. This was reproduced before editing the code.

Two row-collection sites now retain the contradiction: RHS-only support in a
seed/RHS pair, and RHS-only rows imported from a condition matrix. They still
trigger the existing early inconsistency verdict. Consistent systems retain
their original rows and numerical solution path.

A separate direct-call test found that the sampled helper in
`linear_solve.hpp` ignored zero-coefficient/nonzero-RHS rows during selection
and verification. All current collinear CLI call paths already guard that
case. The helper now performs the same elementary check itself. Before the
change, 15 of 90 exact tests failed for that case; afterward all 90 passed.
The incremental solver was not changed.

The FEC resume correction advances the predecessor when a cached weight is
reused. Previously, extending an existing two-loop run to three loops could
build the next weight from the wrong predecessor. This is a one-line state
update, not a change to the extension calculation.

## Implemented

- Lossless contradictory condition export and the sampled-helper guard.
- Correct FEC predecessor on partially cached runs.
- Strict integer parsing and tensor-add rational validation, including zero
  denominator rejection without restricting arbitrary-precision rationals.
- Alphabet-definition invalidation. Stale precomputed tensors require
  replacement; an old computation cannot mark a newly edited definition ready.
  Derived first/last-entry projection maps are scheduled for cache validation.
- Project mutations and completion callbacks are serialized per project.
  Flow saves require a matching revision and return 409 on conflicts.
- Serialized browser saves, pending-save flush on navigation, tab-local draft
  recovery, download/discard controls, and protection when draft storage fills.
  Run and export save the current flow first.
- Isolated local step attempts, verification of newly produced nonempty
  outputs, input-content stability checks, atomic replacement per output file,
  rollback on publication errors, and content-based input/output cache checks.
  Compiler steps expose file and conservative directory input manifests.
- Atomic job admission, one active run per project, cancellation coordinated
  with launch/cleanup, defined log decoding, and terminal-state cleanup after
  unexpected worker exceptions.
- A bounded event deque with monotonic IDs and reconnect cursors, streaming
  historical logs, reconciliation of unfinished run records, and opt-in server
  development reload.
- Readable unknown-run/network errors, bounded browser log batches and error
  reporting that consumes its own network failures.
- Structural graph validation, duplicate input/output ownership checks,
  numeric zero preservation and standard `FEC_2`/`LEC_2` weight recognition.
- Three verified exporter defects: malformed Wolfram preflight syntax,
  uninitialized dry-run description, and failure to notice rebuilt binaries.
- Updated setup/API/solver documentation, public regression tests and native
  plus API CI gates. CI configuration was added; remote CI was not run here.

## Verification

| Check | Result |
|---|---|
| Exact small systems, both solver helpers | 90/90; consistency, residuals, uniqueness and null-space checks |
| Export/reimport, both solvers | 6/6; imported contradictions, pair-support contradictions, consistent underdetermined systems |
| Integer/rational rejection and large exact fractions | Passed; rejected inputs publish no output |
| E6 two-loop then three-loop resume | All 10 recorded tensor CRC32 values match |
| API/concurrency/cache/output/SSE/export regression suite | 40/40 |
| Chromium save/recovery/network/768-pixel checks | 12/12 |
| Complete E6 example through the webpage | All 13 steps succeeded; results inspector worked; no page errors |
| Real Wolfram execution in an isolated attempt | Passed; input alphabet file unchanged; sparse tensor published |
| Production web build, Python compilation, RHS registry sync | Passed |

Evidence: [native tests](native-results.txt), [API tests](api-results.txt),
[browser checks](browser-results.json), [E6 webpage run](e6-browser-run.json),
[reference comparison](reference-results.json), [web build](web-build.txt).
Reproduction instructions: [tests/README.md](../../tests/README.md).

Native verification used GCC 14.4, FLINT 3.2.2 and the repository-pinned,
patched SparseRREF. API tests used Python 3.12. Browser tests used Chromium
on Linux. The checks used temporary projects, not existing user projects.

## Remaining work, in priority order

1. **Native file validation and internal cache provenance.** Malformed WXF
   inputs still need parser hardening and sanitizer/fuzz tests. The new local
   cache verifies step inputs/outputs; legacy directory-driven native commands
   also have internal existence-based caches. Their dependency/version rules
   need explicit treatment before claiming safe reuse after arbitrary source
   changes. Use fresh output directories when changing native dataset inputs.
2. **Mathematical validation of the remaining templates.** The pentagon
   symmetry example and outdated E6 RHS example need focused, independent
   mathematical checks before correction. The complete Wolfram-dependent 4p
   example and private NMHV/integrability baselines remain unverified.
3. **Finish standalone export and installation portability.** Share the local
   execution contract, validate skip/fail Wolfram modes, test relocation and
   paths with spaces, and check clean macOS/Windows installation. Exported
   scripts currently write directly and have a separate, weaker cache. The
   existing mimalloc switch and Windows launcher still need correction.
4. **Refine efficiency and broader usability.** Replace conservative directory
   snapshots with complete, precise dependency manifests after measuring real
   workloads. Profile high-volume logs/history; add browser/accessibility
   coverage. At 390 pixels the existing Materials layout still overflows;
   full phone editing has not been implemented.

Operational limits: one server worker per project store; one active run per
project. Directory/Wolfram staging adds disk I/O and temporary space. File
replacement is atomic individually, not across a power loss during multi-file
publication. A forced server/OS crash may leave an external process that needs
manual inspection; records are reconciled, but this is not a durable external
job scheduler. No performance claim is made for untested high-weight jobs.
