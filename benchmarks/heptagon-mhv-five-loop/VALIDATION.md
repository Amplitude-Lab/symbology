# Publication validation

The reference tensor was copied from the certified five-loop run, checked
against SHA256
`a06256013cbb03a7a1b67a18f930e02639cf2e777d48196ae36ba00355661691`,
and packaged with its coordinate manifest. Large local results are excluded.

The clean-start runner first passed the complete four-loop control using the
previously built executables: 24.66 seconds from seed binding through final
comparison. It reconstructed the lower-loop amplitudes and FEC6, certified the
joint kernel, checked the complete collinear divergent part and published
words, and reproduced the recorded E4 byte for byte. This initial control
validates the orchestration; clean-source build and publication checks are
recorded below once completed.

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
