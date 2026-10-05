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
recorded below.

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
