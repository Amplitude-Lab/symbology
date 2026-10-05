# Reproduce the five-loop MHV heptagon bootstrap

This benchmark rebuilds the weight-ten MHV heptagon symbol from the small
tracked inputs in `data/`. It uses FEC8 × LEC2 sewing, imposes integrability and
cyclic/flip/parity invariance together, and solves the combined equations in
stages. It then fixes the coefficients by the lower-loop collinear boundary,
checks the complete divergent part, and checks four independent published word
coefficients. No previously generated `output_*` directory is needed.

The accepted [recursive solution](reference/hepMHV_5L_recursive.wxf) is included
in Git: **6,825,548 bytes**, shape `1 × 17471 × 118`, SHA256
`a06256013cbb03a7a1b67a18f930e02639cf2e777d48196ae36ba00355661691`.
Its coefficients require the specific FEC1…8 and LEC1…2 bases recorded in
[manifest.json](manifest.json). The small [LEC2 tensor](reference/LEC_2.wxf) is
also included. The large forward bases are regenerated, not uploaded.
This is a recursive symbol, not an expanded list of original-alphabet words
or a determination of function-level constants.

## Build and run from a fresh clone

Use Linux or WSL, Python 3, a C++20 compiler, FLINT 3, GMP and TBB. The repository
[environment-linux.yml](../../environment-linux.yml) specifies the tested
compiler/library versions. For an activated conda/micromamba environment:

```bash
git clone --branch symrep-solver-dev https://github.com/Amplitude-Lab/symbology.git
cd symbology
./scripts/setup-sparserref.sh
make -j2 five-loop-tools MIMALLOC=0 \
  CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
  CXXFLAGS="-O3 -march=native -mtune=native -std=c++20 -I. -I$CONDA_PREFIX/include" \
  LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -lflint -lgmp -ltbb"

python3 bench/five_loop_heptagon.py run \
  --output output/five-loop-heptagon \
  --hours 24 --ram-gib 50 --min-available-gib 4
```

With system-installed dependencies, `make five-loop-tools MIMALLOC=0` is
sufficient. `make bench-five-loop BENCH_OUTPUT=output/five-loop-heptagon`
builds and runs with the default resource limits. The driver uses eight workers.
Use a new output path; it refuses to overwrite an existing run. No Git LFS,
external physics dataset, Wolfram installation, or local research cache is
required. Git is used to fetch the pinned SparseRREF source and apply the
committed patch stack.

A smaller complete control runs the same pipeline at four loops:

```bash
python3 bench/five_loop_heptagon.py run \
  --output output/five-loop-control4 --loops 4 --hours 1 --ram-gib 12
```

The expected solution is used only after construction, for comparison. Its
known dimension and published word coefficients are not fitting inputs or
discovery stopping criteria.

## What the runner regenerates

| Stage | Computation and acceptance |
| --- | --- |
| Seed binding | Verify the tracked alphabet, condition, projection and seed hashes. |
| Lower loops | Rebuild FEC2…7 and the two-, three- and four-loop E/R collinear tensors with `compute_rhs`. |
| Forward basis | Extend FEC7 to the complete ordinary FEC8 basis. |
| Basis binding | Check every recursive basis and lower-loop result against its recorded hash, before interpreting the archived coordinates. |
| Actions | Derive cyclic, flip and parity matrices recursively in these actual saved bases. |
| Boundary | Build the complete five-loop shuffle boundary from the reconstructed lower-loop amplitudes. |
| Joint kernel | Solve integrability and product invariance together with complete exact certification. |
| Coefficients | Determine a unique coefficient vector using small collinear projections. |
| Full verification | Check every divergent word against the independent full boundary before writing the accepted recursive symbol. |
| Published words | Check −26,880, −26,880, +26,880, −1,680 independently after construction. |
| Final binding | Check LEC2 coordinates and byte-identical agreement with the committed five-loop solution. |

The four-loop control uses FEC6 × LEC2, regenerates only the lower weights it
needs, and compares its full collinear E4 tensor with the recorded result.

The runner snapshots executables, seeds and controller scripts into its output
directory. `config.json` records every command. `manifest.json` records source,
binary and input hashes. `status.json` shows progress, acceptance and resource
usage; each stage has a log and measurement JSON. `final-binding.reference.json` is
written only after the final reference comparisons pass. A directory or
`candidate_*` file alone is not evidence of a successful run.

The single shared wall-time budget covers all numerical stages. The RAM setting
applies an address-space limit per process and a sampled aggregate process-tree
RSS limit; a separate host-available-memory floor protects other jobs. The
controller stops on a failed command, missing output, missing certificate,
timeout or memory guard. Keep the logs if a run stops; the runner deliberately
does not silently reuse partial or uncertified results. Compiler time is
outside the numerical benchmark.

## Cost and interpretation

The published driver passed a complete cold reproduction in **8,280.83 s
(2 h 18 min 1 s)** with **31.48 GiB sampled peak process-tree RSS**. This includes
rebuilding all forward bases, lower-loop amplitudes, generator actions and the
boundary, then solving and checking the final symbol. The final solution
matched the committed file byte for byte. The joint kernel alone took
6,743.00 s and peaked at 19.33 GiB; FEC8 construction set the overall memory
peak. See [VALIDATION.md](VALIDATION.md) for stage measurements and certificates.

The original successful calculation took **7,796.38 s (2 h 9 min 56 s)** and
**19.40 GiB sampled peak RSS** for the joint solve and physical verification,
starting from cached prerequisites. Its separately measured FEC8 construction
took 629.34 s and peaked at **32.67 GiB**. The large boundary was also prepared
separately. Those measurements are not a cold end-to-end timing. A machine
with about 64 GiB RAM and a 50 GiB calculation budget is the tested workstation
scale; timings depend on CPU and concurrent activity. Allow tens of GiB of
temporary disk space for bases, projections and verification tensors.

The new runner includes all prerequisites, and records their cost separately.
See [VALIDATION.md](VALIDATION.md) for the checks performed before publishing
this package and the [historical audit](../../audits/fiveloop-day-2026-10-04/REPORT.md)
for the original run. The
[general strategy document](../../docs/combined-constraints-benchmark.md)
explains the equations, staged intersection, certification, and the distinction
between product invariance and a full recursive irrep decomposition.

Only code, compact records, small seed inputs, the final solution and LEC2 are
tracked. FEC8, intermediate kernels, generator actions, expanded collinear
tensors, binaries and temporary run directories remain ignored. To use the
committed solution independently, regenerate and verify its recorded basis
chain first; matching tensor dimensions alone does not establish matching
coordinates.
