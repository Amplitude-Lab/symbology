# Symbology

Working prototype for symbol-space bootstrap experiments on the `E6` heptagon. The codebase covers the full bootstrap pipeline: recursive first-entry/last-entry growth and sewing (`bootstrap --extend` / `--sew`), universal projection matrices (`--project`), symmetry invariant-subspace solving (`--solve-symmetry`), collinear constraint solving (`--solve-collinear`), and recursive RHS (collinear boundary) computation (`compute_rhs`). All sparse rational linear algebra goes through [`SparseRREF`](https://github.com/munuxi/SparseRREF).

## What You Need

The executable is a C++20 program. A local build needs:

- a C++20 compiler with `<format>` and chrono time-zone support;

- `make`;

- FLINT and GMP;

- TBB;

- optional: mimalloc (`MIMALLOC=0` builds without it);

- Git, to fetch this repository and `SparseRREF`;

- optional: Wolfram/Mathematica, useful for inspecting or round-tripping WXF `SparseArray` data.

`SparseRREF` is not vendored as a submodule. Put a checkout at repository root, so the headers live under `SparseRREF/`. **Pin the tested commit and apply the bundled fixes** — `scripts/setup-sparserref.sh` does all of it (clone at the pinned commit, apply patches, verify):

```bash
git clone https://github.com/Amplitude-Lab/symbology.git
cd symbology
./scripts/setup-sparserref.sh
```

(Manually: `git clone https://github.com/munuxi/SparseRREF.git`, `git checkout 5bbee55`, `git apply ../patches/sparserref-5bbee55-race-and-init-fix.patch`. Use `--check` to verify an existing checkout.)

> **Required SparseRREF patches.** Upstream `sparse_mat_rref_forward/backward` (used by every
> `--project` / sewing-nullspace computation) publish rewritten rows through plain
> `int` flags with no memory ordering: consumer threads can observe the flag before the
> row data, index a phantom pivot, and dereference null — an intermittent `SIGSEGV`
> (\~10% of runs at weight ≥ 3 on multi-core machines; outputs are correct whenever a
> run survives). `patches/sparserref-5bbee55-race-and-init-fix.patch` fixes this with
> release/acquire atomics, exact-once flag consumption, and loud abort guards, and also
> adds a placement-new fix for `rat_t` values in `sparse_type.h` (present upstream as
> `9874061`). The patch is verified to reproduce the validated state byte-for-byte on
> `5bbee55`. If you track a newer upstream `master` instead, apply
> `patches/sparserref-race-fix.patch` (race fix only; the `sparse_type.h` part is
> already upstream since `9874061`) — it applies cleanly on `b5b3d3c` (v0.4.0).
> Without one of these patches the build works but **will crash intermittently**.
> A detailed bug report (root cause, reproduction data, and fix rationale) suitable for
> sharing upstream is in `patches/BUGREPORT-sparse-mat-rref-race.md`.

## Dependency Setup

### macOS

#### Why Homebrew GCC is required

`bootstrap` uses `std::chrono::zoned_time` (for timestamps) and `SparseRREF` uses `std::execution::par` (parallel algorithms) unconditionally. On macOS **every clang toolchain ships libc++**, and libc++ gates both of these behind build-time feature flags that no available distribution enables — Apple clang, Homebrew clang, and conda-forge clang all fail with `no member named 'par' in namespace 'std::execution'` and `no member named 'zoned_time' in namespace 'std::chrono'`. The `-D_LIBCPP_HAS_PARALLEL_ALGORITHMS` macro does not help (it is decided when libc++ itself is compiled, not at use-site). Conda-forge `gxx` on macOS is also a clang wrapper driving libc++, so it has the same problem.

The reliable route is **Homebrew GCC** (`g++-14`), which brings its own libstdc++ where `std::execution::par` and `zoned_time` are available unconditionally.

#### Setup

Install the libraries and Homebrew GCC:

```bash
xcode-select --install
brew install flint gmp tbb mimalloc gcc
```

Then build:

```bash
make
```

On macOS the Makefile handles two things automatically:

1. **Selects Homebrew GCC.** It picks a `g++-<version>` found under `/opt/homebrew/bin` or `/usr/local/bin` (override with `make CXX=...`, e.g. `make CXX=/opt/homebrew/bin/g++-14`).
2. **Adds the Homebrew include and library paths** (`-I$(brew --prefix)/include`, `-L$(brew --prefix)/lib`, plus an rpath) to the compile and link flags. Neither Apple Clang nor Homebrew GCC searches the Homebrew prefix by default, so without these flags the build fails with `flint/nmod.h: No such file or directory` (compile stage) or `ld: library 'flint' not found` (link stage).

Because of (2), plain `make` works in **any** shell — no `CPATH`/`LIBRARY_PATH`/`CPPFLAGS`/`LDFLAGS` exports are needed. If your login profile already exports them (some Homebrew mirror setup scripts add `export CPATH="$(brew --prefix)/include:$CPATH"` and the matching `LIBRARY_PATH` to `~/.zprofile`), the build used to work only because of those exports; the Makefile now provides the same paths itself, so CI runners, IDE build tasks, and non-login shells behave identically.

#### Non-standard Homebrew prefix

If `brew` is not on `PATH` in the build environment, the Makefile falls back to `/opt/homebrew/bin/brew` and then `/usr/local/bin/brew`. For a truly custom prefix, pass the paths explicitly:

```bash
make CXX=/opt/homebrew/bin/g++-14 \
  CXXFLAGS="-O3 -std=c++20 -I. -I<prefix>/include" \
  LDLIBS="-L<prefix>/lib -Wl,-rpath,<prefix>/lib -lflint -lgmp -lmimalloc -ltbb"
```

#### Linux note

On Linux (the archive's original target) libstdc++ is the default and the plain `make` rule works once a suitable GCC, FLINT 3, GMP and TBB are installed — none of the macOS libc++ issues arise.

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install build-essential make git libflint-dev libgmp-dev libtbb-dev libmimalloc-dev
make
```

Package names can vary slightly across distributions. The required libraries are FLINT 3, GMP and TBB; mimalloc is optional. Check distro versions before installing: an older FLINT package is insufficient.

## Build

There are six executables, all built with `make`:

- `bootstrap` — the main dispatcher (extension, sewing, projection, symmetry solving, collinear solving). Built from `bootstrap.cpp` plus the shared headers.

- `compute_rhs` — standalone recursive RHS (collinear boundary) computation. Built from `compute_rhs.cpp`.

- `inspect_tensors` — small diagnostic tool that prints tensor contents. Built from `inspect_tensors.cpp`.

- `tensor_add` — weighted sum of sparse tensors. Built from `tensor_add.cpp`.

- `tensor_ops` — the general tensor utility CLI the flow editor drives (ternary contraction, matrix power, join, squeeze, shuffle product, integrability conditioning, expansion, …). Built from `tensor_ops.cpp`.

- `letter_filter_bench` — semantics verification + benchmark for the divergent/finite letter support filters. Built from `letter_filter_bench.cpp`.

Each rule in the `Makefile` declares its own header dependencies, so `make` will only rebuild what has changed. The shared headers are: `bootstrap.hpp`, `projection.hpp`, `solve_symmetry.hpp`, `solve_collinear.hpp`, `linear_solve.hpp`, `incremental_solve.hpp`, `tensor_expand.hpp`, `tensor_shuffle.h`.

Build from the repository root:

```bash
make                  # builds bootstrap, compute_rhs, tensor_add and tensor_ops
make inspect_tensors  # optional diagnostic tool
make PORTABLE=1 MIMALLOC=0  # omit host-specific CPU flags and mimalloc
make regression       # re-run the recorded CRC32 baselines (gate for core changes)
```

Clean the binaries with:

```bash
make clean
```

## Design locally, run on the cluster (standalone script export)

A flow drawn in the web editor can be exported as a **single portable,
self-checking bash script** and carried to any machine where the C++ core is
built (a Linux x86_64 cluster, for example) — the front-end is not needed
there. The script is the exact plan the local run engine would execute,
wrapped in safety guards.

> **How to export (in the web UI):**
> 1. Open the flow in the **Flow Editor**.
> 2. Click **Compile** in the toolbar.
> 3. In the right-hand side panel, open the **Plan** tab — the compiled
>    command list appears.
> 4. Click **⇪ Export standalone script** (under "▶ Run this plan").
>
> The script is written to `front-end/projects/<project>/exported/<flow>.sh`,
> and the panel prints the invocation to run it elsewhere.

Then, on the target machine:

```bash
# one-time setup (clone the repo there, or copy it):
./scripts/setup-sparserref.sh && make bootstrap compute_rhs tensor_add tensor_ops     # SparseRREF (pinned) + build the core

# sync/copy the project dir, then run the exported script:
SYMBOLOGY_ROOT=/path/to/symbology \
PROJ_DIR=/path/to/synced/project-dir \
bash /path/to/project-dir/exported/<flow>.sh
```

The script embeds the same staging and cache implementation as the local
runner, using Bash and Python 3.10+ without the frontend dependencies. It
copies inputs into an isolated attempt, checks newly produced outputs and
unchanged inputs, and publishes only successful results. Cache hits require
matching SHA-256 input, executable and output contents. Local and exported
runs share a per-project operating-system lock. Logs go to `runs/export-*.log`.

- `--dry-run` prints the plan; `STEPS="3,5-9"` selects steps.
- `STEP_TIMEOUT=<seconds>` applies a timeout to each step, without requiring
  the GNU `timeout` utility. On POSIX, cancellation kills its process group.
- `SYMBOLOGY_ROOT` selects the built repository. `PROJ_DIR` defaults to the
  parent of the script's `exported/` directory. Paths may contain spaces or
  apostrophes. `PYTHON` and `WOLFRAMSCRIPT` override those executables.
- Ready property tensors require no Wolfram steps. If a plan includes Wolfram
  steps and Wolfram is unavailable, `WOLFRAM_MODE=auto` accepts only previously
  verified outputs with unchanged export-time inputs and output digests.
  `skip` explicitly permits precomputed export-time files with the same
  content checks; missing, empty or changed files still fail. `fail` requires
  Wolfram for those steps. Copy the complete computed project before exporting
  or install Wolfram when these checks cannot establish a usable result.

Cluster prerequisites: a C++20 GCC/libstdc++ toolchain (tested with GCC 14),
FLINT 3, GMP, TBB and Python 3.10+. Mimalloc is optional. Use `PORTABLE=1` for
binaries intended for other CPUs of the same architecture. Build separately
for different operating systems or architectures.

The tested Linux x86_64 dependency versions are in
[environment-linux.yml](environment-linux.yml). CI creates that environment
using the documented [setup-micromamba action](https://github.com/mamba-org/setup-micromamba).
An environment file pins the main package versions; it is not a complete
transitive dependency lockfile.

## The visual front-end (web editor)

The flow editor is a local web app: a FastAPI server that compiles graphs and
drives the C++ binaries, plus a React client it serves. Fresh-machine setup:

```bash
# 1) build the C++ core first — the server launches these binaries:
./scripts/setup-sparserref.sh && make bootstrap compute_rhs tensor_add tensor_ops

# 2) Python server (tested with Python 3.12):
cd front-end/server
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python run.py            # serves the app at http://127.0.0.1:8321

# 3) web client (second terminal, starting at the repository root):
cd front-end/web
npm ci
npm run build                      # builds dist/, served by the server at :8321
# or for development with hot reload:
npm run dev                        # http://localhost:5173 (API calls go to :8321)
```

Notes:

- Projects (graphs, tensors, run history, exported scripts) live under
  `front-end/projects/<name>/` — local data, gitignored; new projects can be
  seeded from the built-in templates.
- Wolfram-based alphabet properties need `wolframscript` on this machine;
  everything else (drawing, compiling, running flows, exporting scripts)
  does not.
- Platforms: **macOS and Linux are the execution targets** (the current public regression results were verified on Linux) (the run engine uses
  POSIX process groups for cancellation and timeouts). On Windows the editor
  and server run, but executing flows inside the web UI — and the exported
  bash scripts — need a POSIX environment such as WSL (the C++ core requires
  GCC/libstdc++ in any case).
- `make check` (regression baselines + registry sync) needs local project
  data and therefore runs meaningfully on the machine where the heptagon
  project lives; CI runs the machine-independent subset.

### Saving and local run guarantees

Flow updates carry a revision; an outdated tab receives HTTP 409. Saves from
one editor are serialized. Pending changes are kept in that tab's browser
session storage and flushed on navigation when auto save is enabled. After a
failed save, use **Download draft** to retain a copy or **Reload saved flow** to
discard it explicitly. Keep the page open if browser storage is unavailable.

Changing an alphabet definition marks dependent properties stale. Computable
properties must be regenerated; precomputed properties must be replaced for
the new definition. Renaming an alphabet does not invalidate its mathematics.

The local server permits one active run per project and computes each step in
an isolated attempt directory. Declared outputs must be newly produced and
nonempty; inputs must stay unchanged during the step. Validated files replace
saved outputs atomically per file, with rollback on publication errors. An
abrupt machine/power failure during a multi-file publication is not a database
transaction. Directory-driven and Wolfram steps conservatively copy their
input directories, requiring extra temporary disk space and I/O. Explicit
file-based tensor steps copy only their inputs.

Cache records check SHA-256 input/executable and output contents. Old cache
records are invalidated once. Cancellation is acknowledged after the child
process exits and the attempt is discarded. Run logs have numbered events and
reconnect cursors; older events remain in the disk log. After a restart,
unfinished records are marked failed. A forced server/OS crash can leave
external processes requiring manual inspection. Run one server worker per
project store. Development reload is opt-in with `SYMBOLOGY_RELOAD=1`.

Newly exported scripts use the same publication and cache protections.
Previously exported scripts retain the implementation embedded at export time;
export them again to obtain these corrections. The native directory-driven
RHS and collinear projection routines also verify `.native-cache` dependency
receipts. Legacy intermediates without receipts are recomputed once. Native
CLI commands themselves still write directly; use local/exported runs for
isolated publication, and do not run raw native commands concurrently against
the same output directory.

Run the portable correctness checks in [tests/README.md](tests/README.md).

## Multi-Project Layout

By default, every executable reads seed tensors from `data/` and writes outputs to `output/`. This is the `E6` problem that ships with the repository. For a different symmetry group or a different bootstrap project, use a per-project directory pair:

```
data_<PROJECT>/      # e.g. data_E7/, data_D5/
output_<PROJECT>/    # e.g. output_E7/, output_D5/
```

The default (no tag) is `data/` + `output/`, which is backward compatible.

### Conventions

- Inside `data_<PROJECT>/`, the seed files keep the same **roles** as in `data/` (see `data/DESCRIPTION.md`), but encode the symmetry group in the filename where it matters: `dlogmat_<group>.wxf` (e.g. `dlogmat_E7.wxf`), `FEC_1.wxf`, `LEC_1.wxf`, `colmat<N>.wxf` (where `<N>` is the FEC weight-1 dimension — `42` for `E6`), `colprojdiv.wxf`, `colprojfin.wxf`, `<group>repmat_*.wxf`, `E1.wxf`.

- All executables accept `--data-dir <dir>` and `--output-dir <dir>`. Relative paths are resolved against the executable directory. The flags are threaded through every pipeline mode: `bootstrap --project` / `--solve-symmetry` / `--solve-collinear`, `compute_rhs`, and `inspect_tensors`. When `compute_rhs` shells out to `./bootstrap --extend` / `--sew` / `--project` (to generate missing prerequisites) and `./bootstrap --solve-collinear` (to solve the collinear constraint at each loop order), it passes the absolute `--data-dir` / `--output-dir` to the subprocess so the same project directories are used end to end. The `--letter-projection <file|identity|divergent|finite>` flag is also threaded from `compute_rhs` to the `--solve-collinear` subprocess (file paths are made absolute first, since the subprocess resolves relative paths against its own executable directory; sentinels pass through verbatim).

- The driver scripts (`run_workflow.sh`, `run_projection.sh`, `run_solve.sh`) honor a `PROJECT=<name>` environment variable: setting `PROJECT=E7` makes them use `data_E7/` + `output_E7/`. With `PROJECT` unset they default to `data/` + `output/` (backward compatible).

- `run_workflow.sh` auto-detects the condition tensor as `dlogmat_*.wxf` inside the data directory, so it generalizes to other symmetry groups without editing the script.

## Minimal Smoke Test

After building, run one forward step, one backward step, and one sewing step:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_2.wxf -l output/LEC_2.wxf -o output/SEW_2p2.wxf
```

A successful run prints tensor ranks, dimensions, nonzero counts, RREF timing information, and CRC32 values for the files it reads and writes. `output/` is generated locally and is ignored by git.

For a longer but still local workflow, use:

```bash
./run_workflow.sh smoke
```

This runs:

- `FEC_1 -> FEC_6`;

- `LEC_1 -> LEC_4`;

- `SEW_2p2`, `SEW_3p1`, `SEW_4p2`, `SEW_5p1`.

The full workflow is available as:

```bash
./run_workflow.sh full
```

It lists the intended first-stage path through `FEC_8`, `LEC_5`, and `SEW_8p2`, but it is not expected to be practical on a normal laptop.

### Projection and solving smoke test

After running the bootstrap smoke test above (or `run_workflow.sh smoke`), the SEW tensors exist. Run the projection, symmetry solving, and RHS computation modules:

```bash
./bootstrap --project --symmetry collinear --target SEW_5p1
./bootstrap --solve-symmetry --symmetry cyclic --target SEW_5p1
./compute_rhs --target SEW_3p1 --letter-projection output/collinear/colprojdiv_w1.wxf    # 2-loop (requires only E1)
./compute_rhs --target SEW_5p1 --letter-projection output/collinear/colprojdiv_w1.wxf    # 3-loop (requires L=2 outputs)
```

Or via the driver scripts:

```bash
./run_projection.sh SEW_5p1   # all four symmetries
./run_solve.sh SEW_5p1        # cyclic, flip, parity
```

A successful `compute_rhs --target SEW_5p1 --letter-projection output/collinear/colprojdiv_w1.wxf` run prints the unique solution (`c[0] = -24, c[1] = 2`), verifies all 32 constraints, and confirms `R3` is divergent-free (no entries at letters `{0, 1}`). Note: `--letter-projection identity` enforces exact cancellation in the full letter space, which is inconsistent for the `E6` example at all `L ≥ 2` — use `colprojdiv_w1` for the standard workflow.

## CLI

### Core bootstrap (`./bootstrap`)

Forward extension:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf
```

Backward extension:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf
```

Sewing:

```bash
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_4.wxf -l output/LEC_2.wxf -o output/SEW_4p2.wxf
```

Projection matrix (collinear or symmetry):

```bash
./bootstrap --project --symmetry collinear --target SEW_5p1
./bootstrap --project --symmetry cyclic    --target SEW_5p1
```

Symmetry invariant subspace solver:

```bash
./bootstrap --solve-symmetry --symmetry cyclic --target SEW_5p1
```

Collinear constraint solver:

```bash
./bootstrap --solve-collinear --target SEW_5p1 --rhs output/3loop/boundary_3L.wxf --projection divergent \
    --letter-projection output/collinear/colprojdiv_w1.wxf
./bootstrap --solve-collinear --target SEW_3p1 --rhs 0 --projection divergent --letter-projection output/collinear/colprojdiv_w1.wxf   # empty RHS
```

Options:

- `--extend`: grow either forward (`-f/--first`) or backward (`-l/--last`) data by one weight;

- `--sew`: combine a forward tensor and a backward tensor into a sewing matrix;

- `--induce`: reserved for future induced-transformation workflows;

- `--project`: run the universal projection pipeline (requires `--symmetry`, `--target`);

- `--solve-symmetry`: compute the invariant subspace of a target's projection (requires `--symmetry`, `--target`);

- `--solve-collinear`: a general collinear-like constraint solver — finite/divergent split + expansion + linear solve (requires `--target`, `--rhs`, `--projection`, `--letter-projection`; `--basis` optional). The letter-space projection in the last step is **user-selectable** via `--letter-projection` (not hardcoded), so the same solver works for any collinear-like projection;

- `--symmetry <collinear|cyclic|flip|parity>`: symmetry name for `--project` / `--solve-symmetry`;

- `--target <SEW_FpL|FEC_W|LEC_W>`: target name (e.g. `SEW_5p1`, `FEC_3`, `LEC_2`);

- `--rhs <rhs.wxf>` or `--rhs 0`: RHS path for `--solve-collinear`; `"0"` means an all-zero RHS constructed in-memory. Missing → exit code 1;

- `--projection <finite|divergent|none>`: which projection to apply in `--solve-collinear` (required — no default). `none` is used with `--target-basis` (custom-seed mode, no auto-projection of the seed);

- `--letter-projection <file|identity|divergent|finite>`: letter-slot projection matrix for `--solve-collinear` (required — no default). A path (e.g. `output/collinear/colprojdiv_w1.wxf`) projects each 11-dim letter slot to a lower-dim subspace; the literal `identity` is the **do-nothing** value (no projection applied). The sentinels `divergent` / `finite` are **support filters**: `divergent` keeps entries whose letter key contains at least one divergent letter (the nonzero rows of `colprojdiv.wxf` in `--data-dir`), `finite` keeps entries where every letter is finite — they drop entries whole without changing dimensions, and apply to both sides of each pair. Under union matching, `c·A` must equal `boundary` exactly at every position where either is nonzero — `identity` enforces this in the full 11-dim space (typically inconsistent for `E6`), while a projection enforces it in the projected subspace (where supports coincide). Relative paths resolve against the executable directory;

- `--pair <seed.wxf> <rhs.wxf|0> <letter|identity|divergent|finite>`: multi-pair mode (repeatable). Each pair contributes its own seed tensor + RHS + per-pair letter projection; all rows are solved together in one system. Mutually exclusive with `--target` / `--target-basis` / `--rhs` / `--projection` / `--letter-projection` (each pair carries its own — combining them is an error);

- `--pair-cond <cond.wxf>`: extra constraint matrix re-ingested into a multi-pair solve (repeatable). Produced by `--export-conditions`; `[M | r]` rank-2 form, rhs = last column;

- `--export-conditions`: multi-pair mode. Write `output/collinear/cond_<stem>.wxf`: the combined non-homogeneous constraints as a rank-2 `[M | r]` matrix (n\_unknowns+1 columns, rhs = last column), one row per stacked constraint row;

- `--out-stem <name>`: override `sol_`/`cond_` output naming for multi-pair mode (default: the first pair's seed stem if there is exactly one pair/condition, else `<first-stem>_x<N>`);

- `--target-basis <seed.wxf>`: custom-seed mode for `--solve-collinear`. Use the given seed tensor directly (no `--target` auto-derivation, no automatic projection of the seed — pass `--projection none`); combined with `--rhs` and `--letter-projection` as usual;

- `--solver <incremental|sampled>`: linear solver used by `--solve-collinear` (default `incremental`);

- `--basis <basis.wxf>`: expansion basis file (repeatable; highest weight first). Auto-detected as `first_w{N}_basis.wxf` if omitted;

- `--data-dir <dir>`: data directory with seed files (default: `<exec_dir>/data`). Used by `--project`, `--solve-symmetry`, `--solve-collinear`; ignored by `--extend` / `--sew` (which use explicit `-c`/`-f`/`-l`/`-o` paths);

- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`). Same scope as `--data-dir`;

- `-c/--condition`: condition tensor, currently `data/dlogmat_E6.wxf`;

- `-f/--first`, `-l/--last`, `-o/--output`: input/output file paths;

- `-h/--help`: print usage.

Thread count is chosen automatically by `SparseRREF`.

### RHS computation (`./compute_rhs`)

```bash
./compute_rhs --target SEW_3p1 --letter-projection output/collinear/colprojdiv_w1.wxf    # 2-loop: computes E2, R2, boundary_2L
./compute_rhs --target SEW_5p1 --letter-projection output/collinear/colprojdiv_w1.wxf    # 3-loop: computes E3, R3, boundary_3L (requires L=2 outputs)
./compute_rhs --target SEW_5p1 --letter-projection output/collinear/colprojdiv_w1.wxf --data-dir data_E7 --output-dir output_E7   # multi-project
```

Options:

- `--target <SEW_FpL>`: target SEW name (required). Loop order `L = (F+L)/2`; supported `L = 2..5`;

- `--letter-projection <file|identity|divergent|finite>`: letter-slot projection matrix (required — no default). A path (e.g. `output/collinear/colprojdiv_w1.wxf`) projects each 11-dim letter slot to a lower-dim subspace; `identity` is the **do-nothing** value (no projection applied). The sentinels `divergent` / `finite` are **support filters** (see `bootstrap` options above). Under union matching, `c·A` must equal `boundary` exactly — `identity` enforces this in the full space (typically inconsistent for `E6`), while a projection enforces it in the projected subspace. Relative paths resolve against the executable directory. Threaded through to the `--solve-collinear` subprocess;

- `--data-dir <dir>`: data directory with seed files (default: `<exec_dir>/data`);

- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`);

- `-h/--help`: print usage.

### Inspection (`./inspect_tensors`)

```bash
./inspect_tensors                                  # default: reads from ./output/
./inspect_tensors --output-dir output_E7           # multi-project
```

Options:

- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`);

- `--data-dir <dir>`: accepted for symmetry with the other tools (not used by `inspect_tensors`);

- `-h/--help`: print usage.

Reads `<output-dir>/oneloop/E1.wxf` and `<output-dir>/2loop/boundary_2L.wxf`.

## Tensor Files

Tracked seed data (see `data/DESCRIPTION.md` for a complete listing with dimensions):

- `data/dlogmat_E6.wxf`: RREF-reduced `E6` adjacency/integrability condition tensor. It corresponds to `bootstrap_E6_archive/dlogmatE6RREF.wxf` and has dimensions `{42, 42, 1191}`.

- `data/FEC_1.wxf`: forward expansion coefficient seed, copied from `bootstrap_E6_archive/FCC_1.wxf`; layout `{basis_w, basis_{w-1}, letter}`.

- `data/LEC_1.wxf`: backward expansion coefficient seed, obtained from `bootstrap_E6_archive/LCC_1.wxf` by transposing to the new layout `{basis_w, letter, basis_{w-1}}`.

- `data/colmat42.wxf`: collinear seed on the FEC weight-1 space (`42 × 2`).

- `data/cycrepmat.wxf`, `data/fliprepmat.wxf`, `data/parityrepmat.wxf`: cyclic / flip / parity symmetry representation matrices on the FEC weight-1 space (`42 × 42`).

- `data/colprojdiv.wxf`: weight-1 colprojdiv seed (`11 × 2`); projects each letter slot to its 2-dim divergent subspace.

- `data/colprojfin.wxf`: weight-1 colprojfin seed (`11 × 9`); projects each letter slot to its 9-dim finite subspace.

- `data/E1.wxf`: one-loop collinear seed tensor (`11 × 11`, 5 nnz); renamed from the archive's `coloneloop.wxf`. Used by `compute_rhs`.

Generated files (under `output/`):

- `output/FEC_w.wxf`: forward expansion coefficients;

- `output/LEC_w.wxf`: backward expansion coefficients;

- `output/SEW_fpl.wxf`: sewing matrices with layout `{sew_basis, FEC_f_basis, LEC_l_basis}`;

- `output/collinear/`: collinear projection chain — `first_w{N}.wxf`, `last_w{N}.wxf`, `first_w{N}_basis.wxf`, `last_w{N}_basis.wxf`, `SEW_<name>_basis.wxf`, `colprojfin_w{N}.wxf`, `colprojdiv_w{N}.wxf`, `colprojfin_<sew_name>.wxf`, `colprojdiv_<sew_name>.wxf`, plus a `summary.txt`;

- `output/cyclic/`, `output/flip/`, `output/parity/`: symmetry projections — `first_w{N}.wxf`, `last_w{N}.wxf`, `SEW_<name>.wxf`, `<target>_invariant.wxf`, plus a `summary.txt`;

- `output/oneloop/E1.wxf`: copy of `data/E1.wxf` (written by `compute_rhs`);

- `output/{L}loop/` (digit prefix — `2loop`, `3loop`, `4loop`, `5loop`): per-loop results from `compute_rhs` — `solMHV_LL.wxf`, `hepMHV_LL.wxf`, `E_LL.wxf`, `R_LL.wxf`, `boundary_LL.wxf`;

- `logs/*.log`: stdout/stderr logs for each workflow step.

`output/`, `output_*/`, `logs/`, the compiled `bootstrap` / `compute_rhs` / `inspect_tensors` executables, `temp/`, and `tmp/` are ignored by git.

## SymbolBootstrap.wl

`SymbolBootstrap.wl` is the primary Wolfram Language package for generating symbol-level constraint tensors (dlogmat) from an alphabet of symbol letters. It depends on `SparseRREF/SparseRREF.wl` (clone SparseRREF into the repository root, or place it next to `SymbolBootstrap.wl`).

Load the package with:

```wolfram
Get[FileNameJoin[{<repo root>, "SymbolBootstrap.wl"}]];
```

### Workflow overview

1. **Provide an alphabet definition** (`alphabet.wl`): a Wolfram Language file defining the letter expressions as functions of kinematic variables, plus any square-root substitutions. See `data_pentagon/alphabet.wl` (variables `LetterRep`, `RootDef`) and `data_4pformfactor/alphabet.wl` (variables `alphabetf`, `sqrtrep`) for templates.
2. **Document the kinematics** (`Description.md`): record the variable definitions, letter classification, and any symmetry transformations (cyclic, flip, Galois) as kinematic substitution rules. See [data\_pentagon/Description.md](data_pentagon/Description.md) and [data\_4pformfactor/Description.md](data_4pformfactor/Description.md) for templates.
3. **Declare the alphabet** in Wolfram Language and **set its parametrized expressions** (letters as functions of kinematic variables, including any square roots).
4. **Set conditions** as needed: cluster adjacency, extended Steinmann, first/last entry, letter transformations (using the kinematic rules from `Description.md`).
5. **Request condition tensors** via the `Get*` functions. Results are cached per alphabet; re-setting a condition clears its cached tensor.
6. **Export** the returned `SparseArray` to WXF for use by the C++ pipeline, or combine multiple tensors with `GetAlphabetConditionTensor`.

### Alphabet management

| Function                                    | Description                                                                                                                                                                                                    |
| ------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `DeclareAlphabet[name, alphabet]`           | Register a new symbol alphabet (e.g. `{W[1], ..., W[n]}`). Fails if `name` is already declared.                                                                                                                |
| `ResetAlphabet[name, alphabet]`             | Overwrite an existing alphabet and clear its conditions/results.                                                                                                                                               |
| `ClearAlphabet[name]`                       | Remove an alphabet and all its conditions/results.                                                                                                                                                             |
| `SetAlphabetExpression[name, alphabetExpr]` | Set the parametrized letter expressions (a vector parallel to `alphabet`). Required for integrability and letter-transformation tensors. Re-setting clears the integrability and letter-transformation caches. |

### Condition setters

| Function                                            | Input                                               | Description                                                                                                                                                                |
| --------------------------------------------------- | --------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SetClusterAdjacency[name, adjpairs]`               | n×2 matrix of ordered adjacent pairs `{W[i], W[j]}` | Pairs of letters that **may** appear adjacent in symbol words.                                                                                                             |
| `SetExtendedSteinmann[name, nonadjpairs]`           | n×2 matrix of ordered non-adjacent pairs            | Pairs of letters that **may not** appear adjacent (Steinmann/extended-Steinmann relations). Internally converted to the complement adjacency set.                          |
| `SetFirstEntry[name, firstentry]`                   | vector of letters                                   | Letters allowed as the first entry of a symbol word.                                                                                                                       |
| `SetLastEntry[name, lastentry]`                     | vector of letters                                   | Letters allowed as the last entry of a symbol word.                                                                                                                        |
| `SetLetterTransformation[name, transName, kineMap]` | name + kinematic substitution rule                  | Register one named transformation (e.g. `"Cyclic"`, `"Flip"`, `"Galois1a"`) as a rule on the kinematic variables. Multiple transformations can coexist on one alphabet.    |
| `SetAlphabetCondition[name, key, content]`          | dispatch form of the above                          | Uniform setter; `key` is one of `"Expression"`, `"Cluster Adjacency"`, `"Extended Steinmann"`, `"First Entry"`, `"Last Entry"`, or `{"Letter Transformation", transName}`. |

### Tensor generators (Get\*)

All `Get*` functions require the corresponding condition to have been set and return a `SparseArray` (CSR format). Results are cached on the alphabet.

| Function                                               | Returns                             | Shape       | Notes                                                                                                                                                                                                             |
| ------------------------------------------------------ | ----------------------------------- | ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `GetIntegrabilityTensor[name, opts]`                   | Integrability tensor (dlog∧dlog=0)  | `{n, n, r}` | Uses `GenSqrtD` (sqrt separation + Z2 reduction + denominator rationalization) then `GenIntRelMat` (numeric sampling + SparseRREF). Supports square-root alphabets.                                               |
| `GetClusterAdjacencyTensor[name]`                      | Cluster adjacency condition tensor  | `{n, n, r}` | Built from the null space of the adjacency coefficient array.                                                                                                                                                     |
| `GetExtendedSteinmannTensor[name]`                     | Extended Steinmann condition tensor | `{n, n, r}` | Complement of cluster adjacency: takes the non-adjacent pairs and internally calls the cluster-adjacency generator on the complement.                                                                             |
| `GetFirstEntryTensor[name]`                            | First-entry seed tensor             | `{k, 1, n}` | One row per allowed first letter.                                                                                                                                                                                 |
| `GetLastEntryTensor[name]`                             | Last-entry seed tensor              | `{k, n, 1}` | One row per allowed last letter.                                                                                                                                                                                  |
| `GetLetterTransformationTensor[name, transName, opts]` | Transformation matrix               | `{n, n}`    | Square matrix mapping old dlog vector to new dlog vector under the kinematic substitution. Uses `GenSqrtD` on the joined old+new alphabet, then `GenLettRelMat`.                                                  |
| `GetAlphabetConditionTensor[name, key, opts]`          | Dispatch form                       | —           | Calls the matching `Get*` above. Also accepts a **list** of dlogmat-type conditions (`{"Integrability", "Extended Steinmann", ...}`) and returns their combined, row-reduced tensor via `CombineConditionTensor`. |

### Tensor ↔ expression conversion

| Function               | Input                    | Returns    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| ---------------------- | ------------------------ | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `SA2Exp[sarray]`       | `SparseArray` (any rank) | expression | Converts a sparse tensor to `Σ v · S[i, j, ...]`, one `S[...]` per nonzero position `{i, j, ...}` with its value `v` as coefficient. All-zero tensors give `0`. The `S` head resolves to ``Global`S``, matching the `Exp2SA` convention used by the front-end's symbol-tensor property scripts; `Exp2SA[SA2Exp[t]]` reproduces all nonzeros of `t`, and `Exp2SA[e, "dim" -> d]` restores the exact dimensions when each axis's max index is not tight. |
| `SA2Exp[sarray, head]` | + custom head symbol     | expression | Same with a different head (e.g. `W`), for displaying against a concrete alphabet.                                                                                                                                                                                                                                                                                                                                                                     |

`SA2Exp` is the inverse direction of the `Exp2SA` helper generated by the front-end server ([wolfram.py](front-end/server/app/wolfram.py)); typical use is turning a solved collinear combination `c . seed` back into a checkable expression in the symbol letters.

### Options

The integrability and letter-transformation generators accept these options (inherited from `GenIntRelMat` / `GenLettRelMat`):

- `"Samples" -> Automatic` (default): number of numeric sampling points. `Automatic` picks `10 + Ceiling[Binomial[n,2]/Binomial[v,2]]` for integrability and `10 + 2 Ceiling[n/v]` for letter transformations, where `n` is the letter count and `v` the variable count. Increase if reconstruction fails.

- `"Tries" -> 100`: max attempts to find a sampling point that avoids zero denominators.

- `"Threads" -> 0`: thread count for SparseRREF (0 = automatic).

- `"Verbose" -> True`: print progress and timing.

### Square-root alphabets

`GenSqrtD` handles alphabets containing square roots:

1. Each `Sqrt[x]` is rewritten to a custom `sqrt[x]` symbol.
2. The set of square roots is Z2-reduced (row reduction mod 2) to find the independent roots and express dependent ones as products.
3. Denominators containing square roots are rationalized by multiplying by all Galois conjugates of the sqrt-containing denominator part.
4. The resulting dlog expressions are rational functions of the kinematic variables only (no square roots in denominators).

The integrability sampling then expands each minor's coefficient in the basis of independent `sqrt[i]` symbols (2^k coefficients per minor, where k is the number of independent roots), so the reconstructed relations are exact rational identities. This has been verified for alphabets with up to 5 independent square roots (the 4-point form factor example in `data_4pformfactor/`).

### Example: 4-point form factor (`data_4pformfactor/`)

```wolfram
Get["SymbolBootstrap.wl"];
Get["data_4pformfactor/alphabet.wl"];  (* defines alphabetf, sqrtrep *)

DeclareAlphabet["4pFF", Table[W[i], {i, 93}]];
SetAlphabetExpression["4pFF", alphabetf[[All, 2]] /. sqrtrep];

(* Integrability tensor: {93, 93, 3774}, nnz=22092 *)
dlogmat = GetIntegrabilityTensor["4pFF"];
Export["data_4pformfactor/dlogmat_4pformfactor.wxf", dlogmat];

(* Extended Steinmann from non-adjacent pairs *)
SetExtendedSteinmann["4pFF", {
  {W[5], W[6]}, {W[5], W[7]}, {W[5], W[8]}, {W[6], W[5]},
  {W[6], W[7]}, {W[6], W[8]}, {W[7], W[5]}, {W[7], W[6]},
  {W[7], W[8]}, {W[8], W[5]}, {W[8], W[6]}, {W[8], W[7]}}];
esTensor = GetExtendedSteinmannTensor["4pFF"];  (* {93, 93, 12} *)
Export["data_4pformfactor/dlogmatES_4pformfactor.wxf", esTensor];

(* Joined integrability + ES tensor: {93, 93, 3786} *)
joined = GetAlphabetConditionTensor["4pFF", {"Integrability", "Extended Steinmann"}];
Export["data_4pformfactor/dlogmat_full_4pformfactor.wxf", joined];

(* Letter transformation matrices (cyclic, flip, 5 Galois) *)
SetLetterTransformation["4pFF", "Cyclic", {u1 -> u2, u2 -> u3, ...}];
cycmat = GetLetterTransformationTensor["4pFF", "Cyclic"];
Export["data_4pformfactor/cycmat.wxf", cycmat];
```

See [data\_4pformfactor/Description.md](data_4pformfactor/Description.md) for the full variable definitions, letter classification, transformation rules, and verified letter replacement rules.

### Example: pentagon (`data_pentagon/`)

```wolfram
Get["SymbolBootstrap.wl"];
Get["data_pentagon/alphabet.wl"];  (* defines LetterRep, RootDef *)

DeclareAlphabet["Pentagon", Table[W[i], {i, 31}]];
SetAlphabetExpression["Pentagon", LetterRep[[All, 2]] /. RootDef];

dlogmat = GetIntegrabilityTensor["Pentagon"];  (* {31, 31, 361}, nnz=1754 *)
```

## Skills and Changelog

- `skills/` holds per-module reference documents (concise, model-agnostic) for AI agents and new contributors. Start at `skills/README.md`.

- `CHANGELOG.md` records all notable changes (new files, modified files, new functionality) grouped by date.

- `data/DESCRIPTION.md` describes every seed file under `data/`.

## Format Notes

`bootstrap` writes SparseRREF-native WXF. Some archive files were exported through Mathematica, so byte-level CRC32 values can differ even when tensor contents are identical.

For forward files, the current workflow has been checked as follows:

1. Generate `output/FEC_2.wxf` through `output/FEC_6.wxf`.
2. Roundtrip each file through Mathematica (import the WXF `SparseArray`, re-export it) — the historical `wxf_roundtrip.wls` helper this step referred to is no longer in the tree.
3. Compare CRC32 with `bootstrap_E6_archive/FCC_2_rref.wxf` through `FCC_6_rref.wxf`.

After Mathematica roundtrip, all checked `FEC_2..FEC_6` CRC32 values match the archive. `LEC` files intentionally use a different axis order from the old `LCC` files.

## Troubleshooting

- If `make` cannot find `SparseRREF/sparse_mat.h`, clone `SparseRREF` into the repository root.

- If `git` fails on macOS with an `xcode-select` error, install the Apple command-line tools or use the conda-forge setup above.

- If compilation fails on macOS with `no member named 'par' in namespace 'std::execution'` or `no member named 'zoned_time' in namespace 'std::chrono'`, you are hitting the libc++ limitation described in the macOS section above. A newer clang will **not** fix it — build with Homebrew GCC instead: `make CXX=g++-14`.

- If the linker cannot find FLINT, GMP, TBB, or mimalloc, check that the matching include and library paths are visible to `make`.

- If byte-level WXF CRC32 values differ from archived Mathematica exports, compare after a Mathematica roundtrip rather than comparing raw SparseRREF-native WXF bytes.


## Verified example and WXF contracts

The pentagon template now builds **weight-two integrable symbols invariant
under cyclic and flip transformations**. The 4p template imposes the shipped
full integrability/extended-Steinmann conditions, then the same two symmetries.
They begin with the exact unrestricted word basis `B[n*i+j,i,j]=1` (zero-based
indices), solve the original conditions and transform the resulting symbols.
They do not average the condition tensor. The public tests independently
certify 76 pentagon and 689 4p invariant basis elements, including exact
residuals and completeness. The E6 two-loop template wires E1 to the current
RHS node and its output to the collinear solver. Existing user flows are not
rewritten when template definitions change; create a fresh template project
for the corrected examples.

Native tensor I/O accepts uncompressed WXF `8:` CSR `SparseArray` expressions
with rank at least two, implicit zero, exact integer/rational entries and
indices representable in the configured index type (the CLI uses signed
32-bit indices). Zero axes and 64-bit encoded indices whose values fit are
accepted. Compression, approximate/complex entries and other WXF expression
forms are rejected with readable errors. Lengths, shapes, row pointers,
indices and rational denominators are validated before tensor construction.
This is a checked subset of the [Wolfram WXF format](https://reference.wolfram.com/language/tutorial/WXFFormatDescription.html),
not a general Mathematica-expression importer. Current native encoding targets
little-endian machines. Run `make check-wxf-sanitized` for memory diagnostics.
