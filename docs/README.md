# Documentation

Audit and verification reports kept as documentation. They record what was
audited, what was fixed, and how each fix was verified; the current
regression suites that keep those guarantees live in [`tests/`](../tests/).

| Report | What it covers |
|---|---|
| [Core audit](audit-2026-09-09-core.md) | Calculation-core audit against the four design principles (robustness, efficiency, universality, transplantability); source of the regression-gate design. Includes the incremental-solver defect analysis behind the corrected collinear answers. |
| [Robustness audit](audit-2026-09-09-robustness.md) | Adversarial audit of the front-end, native I/O and solver contracts (findings F01–F13), with the refinement plan that the subsequent fixes implemented. |
| [First refinements](audit-2026-09-09-refinement.md) | Verified fixes from the robustness audit: contradiction-preserving condition export, alphabet invalidation, isolated step execution with content-checked caches. |
| [Continued refinements](audit-2026-09-09-refinement-continued.md) | WXF input validation, dependency-aware native caches, corrected pentagon/4p example construction with exact certificates, safer standalone export, narrow-screen layout. |
| [Windows compatibility](audit-2026-09-09-windows.md) | What was actually tested on native Windows (launcher, editor, API), and the WSL route for numerical work. |

The reports reference evidence files (result captures, screenshots, probe
scripts) that are not part of the user-facing repository.

For hands-on use of the web editor, see the separate
[illustrated front-end guide](front-end-guide.md); for a release-notes
style overview of the current build, see
[releases/2026-09.md](releases/2026-09.md).
