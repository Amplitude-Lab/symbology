> **Archived audit report (2026-09-09, robustness).** Kept as documentation of what was
> audited, fixed, and verified. The evidence files it references (result
> captures, screenshots, probe scripts) are not part of the user-facing
> repository; they live in the development history and on the author's
> machine. Findings described as open here were addressed by the following
> refinement reports and the regression suites in `tests/`.

# Symbology robustness audit and proposed refinement plan

Audited on 2026-09-09, at commit `7af08a7`. Application code, dependency manifests, and existing documentation were not changed. This directory contains the audit, executable reproductions, and evidence. Test projects and computational outputs were isolated in temporary directories.

The project can complete meaningful calculations, but I would address correctness and saved-work risks before investing in further numerical speedups. Several protections described in the previous audit are incomplete under concurrent requests, partially cached calculations, malformed inputs, or interrupted browser activity.

## What was actually tested

| Area | Result | Evidence |
|---|---|---|
| Native build | All six C++ executables built with GCC 14.4, FLINT 3.2.2, TBB, and mimalloc in a temporary environment. The original GCC 11 environment failed at `<format>`; a compatible environment without mimalloc failed at its missing header. | [Build log](core-build-results.txt) |
| Web production build | Passed: 209 modules; JavaScript 444.68 kB, 140.46 kB gzip. | [Build log](web-build-results.txt) |
| Backend/API/compiler/orchestration | 31 adversarial checks: **8 passed, 23 failed**. These tests deliberately target edge cases; this is not a statistical reliability score. | [Results](backend-results.txt), [reproduction suite](test_robustness.py) |
| Native CLI and file handling | 45 checks: **22 passed, 23 failed**. Seventeen failures are distinct malformed-file probes of the same parser failure class. | [Results](core-results.json), [runner](core_audit.py) |
| Exact linear systems | 90 deterministic cases, with exact residual, uniqueness, and null-space checks: **75 passed, 15 failed**. Incremental solver: 45/45; sampled helper: 30/45. Important scope qualification below. | [Results](numerical-results.txt), [C++ probe](numerical_probe.cpp) |
| E6 two-loop computation | Ten fresh runs all reproduced **all five** recorded two-loop CRC32 values. Observed wall times were roughly 0.25–0.31 seconds in the initial repeat set. | Repeated-run section in [native results](core-results.json) |
| Three-loop computation | Advancing from two loops failed. Explicitly building the FEC chain first allowed completion; **all ten** two-/three-loop reference CRC32 values then matched. | [Reference comparisons](prebuilt-three-loop.json) |
| Projection and symmetry | The standalone cyclic projection and symmetry solve for `SEW_2p2` completed; invariant dimension 1/1. Collinear projection is also exercised by the two-/three-loop tests. | Temporary core smoke logs; reference comparisons above |
| Built-in templates | Four of five example flows compile. E6 RHS example fails compilation. Pentagon example compiles but fails at its first tensor contraction. E6 smoke executes all 13 steps. | [Compile results](template-results.json), [execution results](template-run-results.json) |
| Browser failure scenarios | Nine checks: **4 passed, 5 failed**, including real lost edits and recursive error reporting. | [Results](browser-results.json), [runner](browser-audit.cjs) |
| Complete webpage workflow | Created an E6 project, compiled and ran the 13-step smoke flow, watched completion, opened Results, and inspected a tensor. All steps completed; 11 result files; no render errors. | [Results](browser-run-results.json), [screenshot](successful-run.png) |
| Large letter-support filter | Million-entry benchmark completed with `TOTAL: ALL PASS`; all compared filter implementations agreed. | [Full output](letter-filter-results.txt) |
| Simple graph compilation scale | 101 / 1,001 / 5,001 nodes compiled in about 1.8 / 13.9 / 73.7 ms. This measures simple matrix-power graphs, not complex custom-block expansion or browser routing. | [Measurements](compile-scale.json) |
| Existing repository gates | RHS registry sync passes. `make regression` cannot run from the clean checkout because recorded input directories are absent. SSR reports success with **zero renders**, because its private project fixtures are absent. | [Regression output](regression-results.txt), [SSR output](ssr-results.txt) |

Browser testing used Chromium 152 on Linux, at 1440-, 768-, and 390-pixel widths. Server tests used Python 3.12, FastAPI 0.141.1, and Pydantic 2.13.5; the web build used Node 24.20. Native audit builds used `-O2`, so timings should not be presented as optimized production benchmarks.

## Confirmed findings, in recommended priority order

Severity: **P0** means a demonstrated wrong-result or saved-work risk; **P1** means a major failure of an advertised workflow or robustness guarantee; **P2** means a more limited correctness, usability, portability, or efficiency gap. Severity describes impact, not frequency.

### F01 — P0: condition export changes the mathematical system

A two-row augmented matrix containing `c0 = 2` and `0 = 1` is correctly rejected as inconsistent. With condition export enabled, the generated file contains only `c0 = 2`. Reimporting that file succeeds and writes a solution with `c0 = 2`.

`ingest_cond_file` counts contradictory zero-coefficient rows and removes them before the export is constructed. The original solve retains the count, but the exported mathematical object does not. The same omission pattern exists in pair-row construction. This is a reproducible semantic corruption, not a numerical tolerance issue.

Source: [condition ingestion](/home/ana/Documents/symbology/solve_collinear.hpp:1266), [export](/home/ana/Documents/symbology/solve_collinear.hpp:1335). Reproduction: `condition export preserves inconsistency` in the native results.

**Proposed change:** preserve contradictory rows in the exported system, or refuse export with explicit failure semantics. Require equivalent consistency, solution space, and rank before and after every condition round trip.

### F02 — P0: changing an alphabet leaves obsolete tensors marked ready

Updating an alphabet's expressions from `{x,y}` to `{x,1-x}` leaves its existing integrability property `ready`, with its previous tensor and dimensions. Compilation trusts a ready property whose file exists, so it need not generate a new property step at all. Improving only the job fingerprint would not fix this path.

Source: [alphabet update](/home/ana/Documents/symbology/front-end/server/app/main.py:135), [property compilation](/home/ana/Documents/symbology/front-end/server/app/compile.py:803).

**Proposed change:** version alphabet definitions and dependent properties; invalidate affected property tensors, derived maps, and downstream computations whenever their definitions change. Give precomputed files an explicit association with the alphabet version they describe.

### F03 — P0: browser saves and concurrent project updates can lose work

Two browser reproductions succeeded:

- Rename a flow and navigate away before the 800 ms autosave delay: the edit disappears.
- Delay the first save, make another edit, let the newer save finish first: the older request overwrites the newer name. The browser still displays the newer name while the stored project contains the older one.

Separately, two concurrent flow-creation requests both return a created flow, but only one survives in `project.json`. Individual reads and writes are locked; the entire read–modify–write transaction is not.

Sources: [autosave](/home/ana/Documents/symbology/front-end/web/src/screens/FlowEditor.jsx:1791), [flow creation](/home/ana/Documents/symbology/front-end/server/app/main.py:571), [storage](/home/ana/Documents/symbology/front-end/server/app/storage.py:92).

**Proposed change:** add transactional project mutations and revision checks, serialize saves per flow, and flush or preserve pending edits when leaving the editor. The UI should report which revision is saved and offer recovery after a failed save. Test across two tabs as well as delayed responses.

### F04 — P0: cache checks can accept obsolete or damaged outputs

Confirmed in isolated engine tests:

- Truncating a cached output to zero bytes does not invalidate its signature.
- Changing an input inside an argv-referenced data directory leaves the fingerprint unchanged. Named projection workflows rely on such directories.
- A command that exits zero and produces nothing is accepted if a nonempty output from an earlier attempt already exists.
- Rebuilding a binary does not invalidate the standalone script's cache; the second exported run kept the first binary's output.

Sources: [engine fingerprints](/home/ana/Documents/symbology/front-end/server/app/jobs.py:185), [cache validation](/home/ana/Documents/symbology/front-end/server/app/jobs.py:226), [step output checks](/home/ana/Documents/symbology/front-end/server/app/jobs.py:314), [export fingerprints](/home/ana/Documents/symbology/front-end/server/app/compile.py:392).

**Proposed change:** make the compiler declare actual input files and executable identities, not infer all dependencies from argv. Stage outputs per attempt and publish them atomically only after validation. Record output integrity and reject empty/damaged cache entries. Share the dependency and cache contract between the local engine and exported execution.

### F05 — P1: cached FEC chains break progression to three loops

After a successful two-loop calculation, requesting `SEW_5p1` generates `FEC_4.wxf` with the **same CRC as FEC_2** (`49d4359a`). Projection then rejects dimensions 7 versus 97.

The loop skips existing FEC files without advancing `prev_fec`. The first missing file is therefore extended from an earlier tensor. Prebuilding `FEC_2` through `FEC_5` explicitly in order avoids the problem and produces all ten correct reference outputs.

Source: [FEC prerequisite generation](/home/ana/Documents/symbology/compute_rhs.hpp:295).

**Proposed change:** advance the predecessor whether a file is reused or generated; validate the chain's weight and dimensions. Cover fresh, partially cached, interrupted, and higher-target reruns in regression tests.

### F06 — P1: run registration, cancellation, and live reporting still race

- Two simultaneous runs targeting the same output both pass the overlap check. The lock is released between checking and registering the new run. The test forces this valid thread interleaving with a barrier.
- Cancelling while a subprocess is being launched can mark the run cancelled yet let that subprocess write an output afterward. A disposable real process reproduced this.
- After a streaming client consumes the 20,000-event buffer, later log/status/end events are lost: the buffer removes earlier items, but the client's index remains an absolute list position.
- Invalid UTF-8 from a subprocess raises out of `_execute` rather than settling the run reliably.
- A persisted `running` record remains `running` after the live engine has been recreated, with no corresponding controllable job.

Sources: [run admission](/home/ana/Documents/symbology/front-end/server/app/jobs.py:118), [subprocess handling](/home/ana/Documents/symbology/front-end/server/app/jobs.py:385), [cancellation](/home/ana/Documents/symbology/front-end/server/app/jobs.py:467), [event stream](/home/ana/Documents/symbology/front-end/server/app/main.py:787), [persisted history](/home/ana/Documents/symbology/front-end/server/app/main.py:705).

**Proposed change:** make admission and output reservations atomic; synchronize launch/cancel transitions; guarantee cleanup and final status in a top-level exception handler; decode logs with a defined error policy. Use monotonically numbered events with bounded storage and reconnect cursors. Reconcile interrupted records at startup.

### F07 — P1: malformed inputs are insufficiently validated

All 17 empty, truncated, or random WXF probes caused `tensor_ops dims` to terminate with **SIGSEGV** rather than a readable validation error. Other native input problems silently succeed: exponent `2junk` is treated as 2, `1.5` as 1, rational text `abc` becomes zero, and `1/0` is accepted by tensor addition.

Sources: [WXF reader](/home/ana/Documents/symbology/SparseRREF/wxf_support.h:387), [integer parsing](/home/ana/Documents/symbology/tensor_ops.cpp:199), [rational parsing](/home/ana/Documents/symbology/tensor_add.cpp:50).

**Proposed change:** validate the WXF envelope and token structure before indexing, then ranks, dimensions, sparse indices, and rational denominators. Require complete consumption of integer/rational strings. Add bounded fuzzing with AddressSanitizer/UndefinedBehaviorSanitizer builds. Rejected input must never publish output.

### F08 — P1: graph validation permits ambiguous or altered plans

The compiler accepts duplicate node IDs, dangling edges, multiple sources on one input port, and two distinct operations writing the same output path. Duplicate IDs are silently collapsed into a dictionary entry. A malformed graph containing a null node is accepted on save and produces HTTP 500 on compilation. Numeric project names and numeric alphabet letters also produce plain-text 500 responses.

Sources: [graph construction](/home/ana/Documents/symbology/front-end/server/app/compile.py:212), [edge selection](/home/ana/Documents/symbology/front-end/server/app/compile.py:746), [flow update](/home/ana/Documents/symbology/front-end/server/app/main.py:606).

**Proposed change:** introduce typed API/graph models and validate before persistence. Enforce unique IDs, valid handles, port cardinality, and unique or deliberately shared output ownership. Return structured errors that highlight the affected nodes and connections.

### F09 — P1: standalone export does not meet its advertised contract

- Every tested `--dry-run` reaches an uninitialized `desc` variable under `set -u` and exits 1.
- A generated plan containing a Wolfram step fails `bash -n`: its preflight contains `[[ -x "$WOLFRAMSCRIPT" (wolfram steps) ]]`.
- Executable changes do not invalidate exported caches, as described in F04.
- By source inspection, `WOLFRAM_MODE=fail` does not prevent skipping when existing outputs are available; `skip` is only considered when the executable is unavailable. These mode behaviors do not match the README's stated controls.

Sources: [Wolfram preflight](/home/ana/Documents/symbology/front-end/server/app/compile.py:429), [dry run](/home/ana/Documents/symbology/front-end/server/app/compile.py:484), [mode logic](/home/ana/Documents/symbology/front-end/server/app/compile.py:506).

**Proposed change:** test generated artifacts as executable products: shell syntax, dry run, fresh run, cached rerun, changed binary/input, missing Wolfram, explicit skip/fail modes, spaces in paths, timeout, and relocation. Make a reviewed execution manifest the common input to both runners.

### F10 — P1: browser error handling amplifies a network failure

A failed New Flow request triggers an unhandled rejection. The global reporter makes a fetch request whose rejection is also unhandled, triggering the reporter again. I observed 21 reporting attempts from one action and deliberately made attempt 21 succeed to stop the reproduction.

An unknown run URL also remains on “Loading run…” indefinitely. API errors are swallowed and the event stream has no user-visible recovery state.

Sources: [global error reporting](/home/ana/Documents/symbology/front-end/web/src/main.jsx:7), [flow actions](/home/ana/Documents/symbology/front-end/web/src/screens/Flows.jsx:14), [run loading](/home/ana/Documents/symbology/front-end/web/src/screens/RunDetail.jsx:14).

**Proposed change:** consume telemetry promise rejections, rate-limit/deduplicate reporting, and handle every user action's error locally. Provide explicit loading, disconnected, failed, not-found, and retry states.

### F11 — P1: built-in examples are not verified against current semantics

The E6 two-loop RHS example still sets an older `target`/`seed` contract; compilation now requires a target weight and the current RHS inputs. The Pentagon example applies 31-dimensional transformations to a `(31,31,361)` condition tensor using the general ternary operation's trailing axes, and fails its first contraction. It also supplies ten averaging weights for eight generated terms; the intended group average should be reviewed when repairing the example.

Sources: [E6 RHS template](/home/ana/Documents/symbology/front-end/server/app/templates.py:206), [Pentagon template](/home/ana/Documents/symbology/front-end/server/app/templates.py:219).

**Proposed change:** make every shipped example a checked fixture. Label tensor axes explicitly so an alphabet axis cannot be mistaken for a condition or basis axis. Verify the mathematical property an example advertises, rather than only successful compilation.

### F12 — P2: general-purpose contracts have small but real gaps

- `Matrix Power` rejects numeric JSON `0` but accepts string `"0"`, despite supporting non-negative powers. A truthiness default discards zero.
- Reusing a standard `FEC_2.wxf` as an extension input fails weight inference. The regular expression accepts `FEC2`, not the repository's `FEC_2` convention. File-only reuse also cannot always reconstruct the necessary provenance.
- All 45 incremental exact-system probes passed. The sampled **library helper** incorrectly accepted all 15 inconsistent systems containing a zero-coefficient/nonzero-RHS row. Current collinear CLI callers guard against this before calling it, so this test does **not** establish that those CLI paths currently return wrong solutions. It establishes that the advertised interchangeable/general solver contract is unsafe for direct reuse.

Sources: [zero exponent](/home/ana/Documents/symbology/front-end/server/app/compile.py:3042), [weight inference](/home/ana/Documents/symbology/front-end/server/app/compile.py:906), [sampled constraint selection](/home/ana/Documents/symbology/linear_solve.hpp:82).

**Proposed change:** preserve numeric zero, centralize target-name parsing, store versioned axis/weight provenance, and require both solver entry points to validate all constraints independently.

### F13 — P2: mobile layout and accessibility need a defined support target

The Materials page fits at 768 pixels. At a 390-pixel viewport it is 689 pixels wide, placing the project selector beyond the visible screen. The editor's fixed side panels also merit a separate narrow-screen design. Several controls use visual labels without explicit accessible associations; a keyboard/screen-reader audit remains needed.

Sources: [fixed project selector](/home/ana/Documents/symbology/front-end/web/src/styles.css:53), [phone screenshot](materials-390.png).

**Proposed change:** first support laptop windows and zoom reliably; then decide whether full phone editing is required. Collapsible panels, accessible names, focus management, and keyboard-operable graph controls are more useful than cosmetic redesign alone.

## Documentation and portability discrepancies

| Stated behavior | Actual behavior / required correction |
|---|---|
| README cluster prerequisite “GCC ≥ 12” | GCC 12 lacks the required standard formatting implementation. GCC's official table lists formatting at 13.1 and documents the evolving chrono support. Publish a tested compiler/version matrix; GCC 14 worked here. [GCC library status](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html). |
| mimalloc is optional and automatically detected | `bootstrap.hpp:16` defines `USE_MIMALLOC` unconditionally. The compatible temporary build failed without `mimalloc.h`. The Makefile can omit the link library while the headers still require it. |
| CI installs all native prerequisites | The workflow installs no mimalloc development package, then compiles code requiring its header. Its “Standalone-script generator smoke” step checks the setup script, not a generated standalone script. See [CI](/home/ana/Documents/symbology/.github/workflows/check.yml:23). |
| Front-end setup: build with `make` | `make`/`make all` build only `bootstrap`; working editor flows also need `tensor_ops` and sometimes `tensor_add`. This is correctly explained elsewhere in the README but not in the fresh front-end setup sequence. |
| Create a venv, install into it, then `python3 run.py` | The documented command does not activate or explicitly use that venv. On a clean system it runs a different interpreter. Ubuntu's `python3-venv` requirement is also absent from the package instructions. |
| Windows launcher starts the editor/server | After `cd server`, it runs `server\.venv\Scripts\python`, duplicating `server` in the path. This is a source-confirmed path defect; Windows execution was not available for testing. See [launcher](/home/ana/Documents/symbology/front-end/start.bat:45). |
| Export supports dry run, explicit Wolfram modes, and safe resume | F04/F09 demonstrate mismatches. These claims need executable examples attached to the documentation. |
| Wolfram is needed only for alphabet properties | Condition merging and the Results tensor-summary path also use Wolfram. The inspector worked on this host, but deployment requirements should distinguish these features. |
| API contract describes skips when outputs exist | The engine also requires fingerprints, and that contract needs to describe input identity, output validation, and failure semantics. Documented compiled-step kinds also omit actual utility kinds such as `tensor_ops` and `tensor_add`. |
| Per-module verification uses `wxf_roundtrip.wls` | `skills/01_bootstrap.md` and `data/DESCRIPTION.md` still reference an absent helper; the main README now acknowledges its absence. |
| Existing audit provides a current status list | Earlier limitations concerning inconsistent-solve exit codes and SSR aliases are marked fixed in a later appendix. The current inconsistent CLI test correctly exits 1, and SSR aliasing works. Consolidate current status so old entries are not mistaken for outstanding findings. |

References: [README](/home/ana/Documents/symbology/README.md:186), [API contract](/home/ana/Documents/symbology/front-end/API_CONTRACT.md), [module guide](/home/ana/Documents/symbology/skills/01_bootstrap.md:91), [data guide](/home/ana/Documents/symbology/data/DESCRIPTION.md:8).

Ubuntu 24.04's documented FLINT package is 3.0.1, so the earlier audit's uncertainty about whether that distribution supplies FLINT 3 can be resolved; the separate mimalloc problem remains. [Ubuntu package details](https://packages.ubuntu.com/noble/libflint-dev).

## Proposed work packages for your selection

| Order | Package | Changes | Completion criteria | Relative size |
|---|---|---|---|---|
| 1 | Preserve mathematical results | F01, F02, F04, F05: lossless condition export, alphabet invalidation, explicit dependencies, staged output publication, correct FEC resume. Include strict rational/integer parsing from F07. | Condition round trips preserve solution sets; changed inputs/binaries invalidate only dependent outputs; fresh and partially cached two-/three-loop runs match reference tensors. | Medium–large |
| 2 | Preserve saved work and control runs reliably | F03, F06, F10: transactional/versioned saves, pending-edit recovery, atomic job admission, cancellation synchronization, durable terminal states, bounded numbered event streams, error reporting that cannot recurse. | Slow/out-of-order saves never lose the newest edit; cancellation leaves no continuing computation; >20,000 events and reconnects deliver a correct final state; restart has no permanently “running” orphan records. | Medium–large |
| 3 | Validate graphs, files, and built-in examples | F07, F08, F11, F12: typed request models, handle/cardinality/output ownership checks, parser hardening, zero/weight fixes, axis metadata, correct templates. | Invalid input returns a readable error without a crash/output; all shipped examples compile and satisfy their advertised operation; both solver helpers pass exact adversarial tests. | Medium–large |
| 4 | Make installation and export reproducible | F09 and documentation table: complete build target, tested toolchain, real optional-mimalloc switch, correct venv/Windows launchers, explicit Wolfram/LibraryLink requirements, working exporter. | Follow the README on a clean Linux environment and a macOS host; relocate an exported project; pass syntax/dry-run/fresh/resume/no-Wolfram/path-with-spaces tests. Windows editor startup gets a separate test. | Medium |
| 5 | Make verification portable and continuous | Commit small, non-private fixtures; generate all necessary baseline prerequisites; run API/browser/native tests in CI; assert SSR rendered something; add bounded sanitizer/fuzz jobs and generated-script tests. | A fresh checkout runs meaningful tests without private `front-end/projects` data; CI reproduces the failures found here before fixes and passes after them. | Medium |
| 6 | Improve measured efficiency and usability | Replace shifted event lists with a deque/ring buffer; stream historical logs and CRC calculation; batch UI log updates; paginate run history; expose a bounded run queue/thread budget. Add keyboard/accessibility improvements and optional narrow-screen layouts. | Profile representative large flows/logs first; demonstrate lower memory/latency with unchanged exact results. Establish desktop/zoom requirements and optionally phone support. | Small–medium per item |

I recommend starting with packages **1 and 2**, with a minimal version of package **5** introduced alongside them so the fixes remain protected. Package **4** is the next priority if cluster transfer or onboarding another machine is imminent. Do not change the currently passing letter-filter algorithm solely for speed: its measured correctness and improvement are valuable. The simple graph compiler also showed no immediate scaling problem in the tested range.

## Limits and reproduction

This audit does not certify every high-weight calculation or every alphabet. The four original regression manifests were not all runnable because they depend on untracked data. The audit instead regenerated and checked the available two-/three-loop reference outputs. The 4p form-factor example compiled but its full Wolfram-dependent workflow was not executed. macOS, Windows, Firefox, Safari, high-concurrency throughput, disk-exhaustion recovery, and native sanitizer builds remain untested. A standalone Wolfram version probe printed Linux Wolfram 14.0 and then a shutdown segmentation-fault message; browser tensor inspection succeeded. That host/runtime behavior was not diagnosed as a repository defect.

The Python tests require the server requirements plus `pytest` and `httpx`; run from the repository root:

```bash
python -m pytest audits/robustness-2026-09-09/test_robustness.py -q
```

Compile `numerical_probe.cpp` with the same include/link flags as the application and run it. `core_audit.py` uses `/tmp/symbology-numerical-probe` by default, or the `NUMERICAL_PROBE` environment variable. The browser runners expect `serve_audit.py` at localhost port 18321 and use Playwright/Chromium; adjust their runtime paths for another machine. Generated project data is temporary, and the assertions intentionally remain failing where the application is defective.

Application source and existing docs remain unchanged so you can choose which work packages to apply. Built executables, downloaded dependencies, and the web bundle are ignored build artifacts; the only new repository deliverables are this audit directory.
