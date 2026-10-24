# Publication validation

The reference tensor was copied from the certified five-loop run, checked
against SHA256
`a06256013cbb03a7a1b67a18f930e02639cf2e777d48196ae36ba00355661691`,
and packaged with its coordinate manifest. Large local results are excluded.

The complete five-loop reproduction from a clean source export passed all
**11 stages** in **8,280.83 seconds (2 h 18 min 1 s)**, with **31.48 GiB sampled
peak process-tree RSS**. No old output directories, forward bases, generator
actions, lower-loop amplitudes or boundary caches were supplied. The pinned
SparseRREF dependency was fetched and patched from scratch. All 48 tracked
source files recorded by the run match published commit `b9ece6e`; the ten
dependency headers are also fingerprinted. The final recursive tensor is
byte-identical to the committed reference, with the SHA256 above.

| Cold-run stage | Wall seconds | Sampled peak RSS (GiB) |
| --- | ---: | ---: |
| Lower-loop amplitudes and FEC2…7 | 76.67 | 5.40 |
| FEC8 construction | 593.26 | 31.48 |
| Recursive generator actions | 70.90 | 12.52 |
| Independent five-loop boundary | 339.57 | 14.08 |
| Exact joint kernel | 6,743.00 | 19.33 |
| Unique physical coefficients | 67.46 | 12.04 |
| Complete collinear verification | 370.26 | 14.34 |
| Four published-word checks | 15.28 | 11.90 |
| Total, including bindings, controller and output hashes | 8,280.83 | 31.48 |

The adaptive discovery pass stopped at 151 parameters and the complete
operator intersection reduced this through 8 and 5 to **4**. Rational
reconstruction and an additional residual prime provided a 121-bit exact
certificate against a 76-bit height bound. Physical verification checked
all **55,728,602 divergent words** in the 88,620,356-term collinear expression.
The four independent published coefficients were −26,880, −26,880, +26,880
and −1,680. Every stage returned zero, with no time, RAM or host-memory guard
triggered. Compiler time is excluded. These are workstation measurements,
not a controlled performance comparison with the earlier cached run.

The compact evidence is retained in
[`cold-five-loop.json`](validation/cold-five-loop.json),
[`cold-five-loop-provenance.json`](validation/cold-five-loop-provenance.json),
the [kernel certificate](validation/cold-five-loop-kernel.txt),
[physical verification](validation/cold-five-loop-verify.txt), and
[published-word checks](validation/cold-five-loop-words.txt).
`$RUN` denotes the fresh local output directory; that directory is not uploaded.

Packaging the evidence exposed a receipt filename collision: a binding check
and its process monitor both wrote `<phase>-binding.json`. All reference
checks passed, but the monitor replaced their detailed receipts with timing
records. The driver now writes separate `<phase>-binding.reference.json`
receipts. All three five-loop binding checks were rerun successfully, and
the final detailed receipt is included in the provenance record. This changes
bookkeeping only; the numerical sources and checks are unchanged.
The corrected driver then passed another complete four-loop control in
**23.01 seconds**, with distinct successful reference receipts and measurement
files for all three binding stages; its
[record](validation/cold-four-loop-receipts.json) is retained.

The clean-start runner first passed the complete four-loop control using the
previously built executables: 24.66 seconds from seed binding through final
comparison. It reconstructed the lower-loop amplitudes and FEC6, certified the
joint kernel, checked the complete collinear divergent part and published
words, and reproduced the recorded E4 byte for byte. This initial control
validates the orchestration; clean-source build and publication checks are
recorded below. No receipt of this initial 24.66 s run is retained in the
package — the figure traces to a local control directory
(`output/five-loop-release-control4/`) that is not uploaded; the retained
four-loop control receipts are the 23.41 s and 23.01 s runs.

One packaging-time drift is disclosed: `cold-five-loop-provenance.json`
froze `reference/manifest.json` at sha256 `f263f0aa…f9ff` at launch, while
the published manifest hashes `b1464efd…b17c` (the manifest was edited
about 135 s after the run started — run start 21:04:25.1, manifest mtime
21:06:40.6, both Oct 4 — before the packaging commit). Every
hash field the binding checks consume — seeds, FEC_2..FEC_8, lower-loop
E/R files, LEC_2 and the solution — agrees between the run's recorded
artifacts and the published manifest; the edit touched no
binding-relevant field.

The selected Git files were exported into an empty directory, with no previous
output trees. The pinned SparseRREF dependency was fetched and patched from
scratch; its headers matched the working dependency byte for byte. All five
required executables compiled there. The complete four-loop pipeline then
passed in **23.41 seconds**, including every basis/lower-loop hash, complete
physical verification and the final E4 byte comparison. Details are in
[`validation/cold-four-loop.json`](validation/cold-four-loop.json); `$RUN`
denotes the fresh output directory.

The process-budget tests exposed a termination race: the monitor waited for
the group leader, which could exit before its descendants. It now waits for
all live members of the owned process group and escalates to SIGKILL after
the grace period even if the leader has exited. All eight checks passed,
including a descendant that ignores SIGTERM. The
[test transcript](validation/process-budget-tests.txt) is retained and these
checks are included in Linux CI. A changed seed was also rejected by the
new reference-binding checker.

The full local public C++ suite also passed: exact kernels, symmetry
representations and recursive actions, NMHV equivalence, staged solving,
collinear reconstruction, file decoding, cache invalidation and certificate
rejection. The [compact transcript](validation/public-checks.txt) records
the checks. CI exposed an old cache-test assumption that the default method
still generated the complete projection chain. The test now exercises both
the original and restricted strategies, checking each strategy's actual
outputs; no numerical algorithm was changed for this correction.

[GitHub Actions run 37265455121](https://github.com/Amplitude-Lab/symbology/actions/runs/37265455121)
passed on commit `043c346`: both C++ allocator configurations, process-budget
tests, WXF sanitizer checks, web build, API tests and Windows editor checks.

## Current-source replay audit (2026-10-24)

The accepted 11-stage runner was audited against the original and cold-run
receipts. It selects `refine-hybrid`, `coefficients`, `verify`, `words` and the
recorded prerequisite stages directly, with no exploratory receipt replay.
The current `tensor_kernel.hpp` (SHA256
`023042217792e325dd61d4ab09afbda948668424df7aff64b018ecc76b650058`)
was rebuilt and revalidated because it postdates the earlier cold receipts.

Targeted checks passed: 106 rational factored extension/sewing cases, 16
product-invariance intersections and certificate-rejection checks, 64 packed
shuffle cases, and all eight process-budget tests. A fresh four-loop control
passed all 11 stages and final E4 binding in 22.963 s at 0.533 GB sampled peak
RSS (`output/five-loop-control4-20261024/status.json`). The recovered runtime
lacked only the unversioned TBB linker alias; linking `libtbb.so.12` directly
resolved the build setup without changing numerical sources.

A fresh full five-loop cold replay then passed all 11 stages and final
reference binding. The recursive tensor is byte-identical to the accepted
reference and all four published-word checks passed. No timeout, memory guard
or host-memory guard triggered (`output/five-loop-heptagon-20261024/status.json`).

| Current five-loop stage | Wall seconds | Sampled peak RSS, GB |
| --- | ---: | ---: |
| Lower-loop amplitudes and FEC2…7 | 71.016 | 5.593 |
| FEC8 construction | 419.236 | 33.463 |
| Recursive generator actions | 49.768 | 13.445 |
| Independent five-loop boundary | 346.770 | 15.113 |
| Exact joint kernel | 4,551.720 | 19.064 |
| Unique physical coefficients | 48.419 | 12.934 |
| Complete collinear verification | 358.879 | 14.134 |
| Four published-word checks | 13.961 | 12.775 |
| Total including bindings and controller | **5,863.275** | **33.463** |

This validates the current source and fixed accepted pipeline, not an
alternative exploratory route or a controlled performance comparison with the
earlier cold run.
