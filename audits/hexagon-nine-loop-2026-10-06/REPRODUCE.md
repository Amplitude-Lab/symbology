# Reproduce the staged nine-loop hexagon calculation

The completed target is the nine-loop MHV weight-18 symbol at FEC14 × LEC4. The unique QQ result, all seven native caches and every completed intermediate stage remain in `output_hexagon_14p4_20261006/cached-nine/`. Run from `/home/ana/Documents/symbology`.

## Resume or inspect the exact production result

The preserved production executable is `output_hexagon_14p4_20261006/joint-cached-package/bin/joint_physics`. This executable includes the verified V3 algorithm and the bounded OPE head cache. The following invocation validates its existing completion receipt and returns without solving again:

```bash
env THREADS=16 PHYSICS_CACHE_MIB=512 OPE_HEAD_CACHE_MIB=4096 \
  output_hexagon_14p4_20261006/joint-cached-package/bin/joint_physics \
  14 4 output_hexagon_14p4_20261006/FEC \
  output_hexagon_14p4_20261006/LEC/MHV \
  output_hexagon_12p4_20261006/reference-package/data \
  output_hexagon_14p4_20261006/cached-nine 6000 reconstruct left
```

Output identity includes executable SHA256, algorithm, contraction direction, scheduling target, basis and physical-input hashes. Retain the exact executable when resuming a managed output directory. Changing code or inputs requires a fresh output root; partial candidate caches must not be relabelled as completed sewing solutions.

## Build and run a fresh calculation

Use the inspected hexagon v2 source package and its matching SparseRREF dependency (23ceb2e, including the WXF push-array fix):

```bash
python3 bench/hexagon_joint_build.py \
  --package /home/ana/Documents/symbology/output_hexagon_12p4_20261006/reference-package \
  --output /path/to/fresh-package \
  --cxx /tmp/symbology-audit-toolchain/env/bin/x86_64-conda-linux-gnu-c++ \
  --runtime-prefix /tmp/symbology-audit-toolchain/env
```

Keep the validated FEC chain through weight 14 and MHV LEC chain through weight 4. A guarded production invocation is:

```bash
python3 bench/measure_process.py \
  --report /path/to/fresh-audit/run.json \
  --log /path/to/fresh-log/run.log \
  --status /path/to/fresh-audit/run.status.json \
  --timeout 28800 --virtual-gib 42 --rss-gib 30 --min-available-gib 6 \
  -- env THREADS=16 PHYSICS_CACHE_MIB=512 OPE_HEAD_CACHE_MIB=4096 \
  /path/to/fresh-package/bin/joint_physics \
  14 4 /path/to/FEC /path/to/LEC/MHV /path/to/hexagon-data \
  /path/to/fresh-output 6000 reconstruct left
```

Use `cache` in place of `reconstruct` for one native prime. The target of 6,000 directions schedules partial sewing discovery; it is not an expected physical dimension. The actual completed production commands, timings and resources are in `cached-nine-first-prime.json` and `cached-nine-reconstruction.json`. The first prime used eight threads; additional primes used sixteen. The one recorded compilation/control suspension is in `reconstruction-pauses.json`.

Every checkpoint contains a full affine object, with origin last. Stages marked `partial-candidate` retain omitted sewing equations for later enforcement. The solver performs a full unsampled sewing sweep, fresh physical checks and two independent reconstruction validation primes. It writes the six-block rational format without an additional presentation RREF.

## Independent verification

The independent validators were compiled separately after the numerical driver was frozen. Their exact sources and executables are recorded in `final-source-manifest.json`; source snapshots are in the experiment's `validator-sources/` directory.

```bash
env THREADS=4 output_hexagon_14p4_20261006/bin/joint_certificate \
  14 4 output_hexagon_14p4_20261006/cached-nine/joint \
  output_hexagon_14p4_20261006/FEC \
  output_hexagon_14p4_20261006/LEC/MHV \
  output_hexagon_12p4_20261006/reference-package/data
```

This checks all original interface rows by direct factor contraction and proves their rational residuals vanish using the denominator-cleared integer height bound. It verifies QQ congruence with every native cache. Any additional residual-only primes are reported separately and are never written as native solve caches; this completed result needed none.

```bash
env THREADS=2 output_hexagon_14p4_20261006/bin/origin_check \
  14 4 output_hexagon_14p4_20261006/cached-nine/joint \
  output_hexagon_14p4_20261006/FEC \
  output_hexagon_14p4_20261006/LEC/MHV \
  output_hexagon_14p4_20261006/audit/origin-cases.tsv
```

All 37 independent count partitions pass exactly. `origin-provenance.json` records the all-order origin prediction and the successful one-, two- and eight-loop controls. This validator neither fits nor changes the symbol.

## Portable Wolfram export

The preserved converter uses the matching fixed SparseRREF writer and reads back every WXF after export:

```bash
env THREADS=4 output_hexagon_12p4_20261006/staged-package/bin/irr_convert \
  export MHV 14 4 output_hexagon_14p4_20261006/cached-nine/joint \
  /path/to/fresh-release \
  --fec-root output_hexagon_14p4_20261006/FEC \
  --lec-root output_hexagon_14p4_20261006/LEC/MHV
```

The completed self-contained release is `output_hexagon_14p4_20261006/release_MHV_w18_14p4/`. Its 19 exact WXF round trips and all file checksums are recorded in `release-validation.json`. Load its `hexagon_symbols.wl` and use `MHVSEW[9]` or `MHVcoeff[letters...]` with 18 letters numbered 1–9. FCE tensors use `{old,letter,new}`; LCE uses `{new,letter,old}`; the unique physical sewing matrix uses `{first,last}`. The physical-copy-primitive-v1 normalization and shared solution ordering are preserved.
