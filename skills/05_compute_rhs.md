# Skill 05 — Recursive RHS (boundary) computation

## Purpose

Recursively compute the collinear boundary `E[L]`, the remainder `R[L]`,
and the boundary tensor `boundary_LL` from loop 2 up to a target loop
order `L`. This is a standalone executable that **couples with** the main
`bootstrap` module: it invokes `./bootstrap --extend`, `--sew`,
`--project` (to generate missing SEW basis files) and
`--solve-collinear` (to solve the collinear constraint at each loop
order).

`--solve-collinear` is a general collinear-like constraint solver: the
letter-space projection used in its last step is **user-selectable via
`--letter-projection`** (not hardcoded). `compute_rhs` threads this
flag through to the subprocess. Pass `identity` to do nothing (solve
in the full letter space); pass a projection file (e.g.
`output/collinear/colprojdiv_w1.wxf`) to project each letter slot.
This choice is important — for the standard `E6` solve you must pass
a divergent projection (`colprojdiv_w1`) at every `L ≥ 2`, because
`E1` has divergent-letter entries and `identity` is inconsistent under
union matching.

The boundary at loop `L` is computed from the lower-loop results by the
master-equation recursion (**since 2026-09-10 the C++ implements this
directly**, replacing the historical per-L closed forms; at `L = 2..3` the
recursion reproduces the closed forms' outputs bit-identically — recorded as
the `compute_rhs` regression baseline):

```
boundary_L = (1/L) · Σ_{k=1}^{L-1} k · (R_k ⊗ E_{L-k}),   R_1 := E_1
```

Historical closed forms (kept for reference; the `L = 5` one was
**known-wrong** before the recursion replaced it — the Wolfram shuffle-word
check pinned its missing terms to `-10·R1^5 - 4·R1^3·R2 - 2·R1·R2^2` shuffle
combinations):

```
L=2: boundary = E1^2 / 2
L=3: boundary = E1^3 / 6 + E1 · R2
L=4: boundary = -E1^4 / 12 + E2^2 / 2 + E1 · R3
L=5: (retired — use the recursion)
```

The loop order is capped at `L ≤ 5` by the shuffle kernel (it refuses total
weight > 11, i.e. `2L` letter slots; see `tensor_shuffle.h`).

where `E_L` is the expanded collinear projection of `hepMHV_LL`, and
`R_L = E_L - boundary_L` is the remainder (which must be divergent-free).

## Flow-editor node (compile-time macro)

The visual flow editor exposes the same recursion as a single **Compute RHS**
node (`compute_rhs` in `flowdefs.js`; compiled by
`front-end/server/app/compile.py::_compile_compute_rhs`). At compile time it
expands into the equivalent sequence of `tensor_ops shuf` / `tensor_add` steps
(it does **not** call the `compute_rhs` binary).

The node's Inspector has an **Object type to generate** select (stored in
`data.mode`, default `mhv_boundary`). The choices come from the mirrored
`RHS_MODES` registries in `web/src/flowdefs.js` and
`server/app/compile.py` — each entry declares `label`, `requires` (seed ports
that must be wired) and `forbids` (ports the mode does not use). To add a new
object type: one `RHS_MODES` entry in each file + one generation branch in
`_compile_compute_rhs` (+ optional Inspector help in `FlowEditor.jsx`). The
node subtitle shows the chosen mode; the Inspector shows live warnings for
missing/stray seed wires; compile rejects unknown modes and wrong wiring with
mode-specific messages.

- **`mhv_boundary` — MHV boundary (boundary_L)** (needs `e1`; forbids `p1`;
  target weight `2L`, even, `L ≤ 5`):
  `boundary_L = (1/L) · Σ_{k=1}^{L-1} k · (R_k ⊗ E_{L-k})` with `R_1 := E1`.
  Master equations and any-order proof (Wolfram shuffle-word algebra,
  verified at every order `L ≤ 6`): `1 + E = exp_⊗(Σ R_k)` gives the
  log-derivative `E_L = (1/L)·Σ_m m·(R_m ⊗ E_{L-m})` (single sum, `E_0 = 1`);
  subtracting `R_L` yields the boundary sum. **The C++ `compute_boundary`
  implements this same recursion since 2026-09-10** (bit-identical to the
  historical closed forms at `L = 2, 3` — recorded in the compute_rhs
  regression baseline); the old closed-form branches, including the
  known-wrong `L = 5` one, are retired.
  **Divergent projection** (the hardcoded C++ logic, proven in the shuffle
  word algebra at every order `L ≤ 6`): the collinear solve matches only the
  *divergent* part of the boundary — the divergent letters are exactly those
  of `E1`; every `R_k` (`k ≥ 2`) is finite (finite remainder `R*`;
  `compute_rhs.hpp` Step 13 verifies `R_L` carries no divergent letter).
  Shuffle products preserve letter support, so `P_fin` is multiplicative
  and closes on the R's:
  `P_fin(boundary_L) = (1/L)·Σ_{k=2}^{L-1} k·(R_k ⊗ F_{L-k})` with
  `F_m := P_fin(E_m)`, `1 + F = exp_⊗(Σ_{k≥2} R_k)`. Finite parts of the
  C++ branches are `0, 0, R2²/2, R2⊗R3` (`L = 2..5`) — even the
  known-wrong L=5 branch has the correct finite part (its error is purely
  divergent). The node emits the **exact** formula; the pure-R finite
  content rides along inside the full `E_k` shuffles and the solve's letter
  filter (`--letter-projection divergent`) drops it downstream — `F_k` is
  never materialized. When only the divergent projection is needed, this
  is the suppressed part.
  The lower-loop `R_k` (`k ≥ 2`) and `E_k` tensors are **auto-loaded** from
  `<project>/output/` (`R2.wxf`, `E2.wxf`, …) — they must exist at compile
  time (same contract as `reuse_output`; the error names the missing file).
- **`e47me67_te` — NMHV E47mE67 (tE_L)** (needs `e1` + `p1`,
  `p1 = hep1LE47mE67`): master equation `tE = T ⊗ (1 + E)`, `T = Σ tP_k`.
  Grading gives `tE_L = tP_L + Σ_{k=1}^{L-1} tP_k ⊗ E_{L-k}`; the t-remainder
  `tP_L` is excluded exactly as `R_L` is from the MHV boundary, so the node
  emits the **sum** `Σ_{k=1}^{L-1} tP_k ⊗ E_{L-k}` (proven against the graded
  master product at every order `L ≤ 6`). `tP_1 = P1 := hep1LE47mE67`.
  For `k ≥ 2` the t-remainder is **not** zero — the weight-4 solve gives
  `tP_2` with 393 nonzeros (purely finite: its divergent part cancels
  against the boundary, exactly as the collinear condition demands) — so a
  missing `tP_k` must never be silently dropped. The remainder recursion
  `tP_k = tE_k − Σ_{j<k} tP_j ⊗ E_{k-j}` is **hardcoded into the node**:
  `tP_k` is taken from the cache `output/tP<k>.wxf` when present, otherwise
  derived from `output/tE<k>.wxf` (the solved E47mE67 tensor of the
  weight-`2k` flow — `sol` dotted with the E47−E67 seed through the
  collinear basis chain) and written back to `output/tP<k>.wxf` for later
  runs. When neither file exists the compile **fails loudly** — silently
  defaulting `tP_2 = 0` once truncated the weight-6 RHS and its solve
  came out inconsistent.

Fractions are folded into the shuffle weight (exact `Fraction(k, L)`), so all
`tensor_add` accumulations are plain `+1`. Shuffles are sequential-only (the
verified variant). The node's `out` port provides the concrete tensor
`output/<target>.wxf`, wireable into any downstream block. Ground truth:
weight-4 MHV output is byte-identical (`cmp`) to the hand-built
`boundary_2L.wxf` shuffle chain in the `heptagon` project.

## CLI entry point

```bash
./compute_rhs --target <SEW_FpL> --letter-projection <file|identity|divergent|finite> \
    [--data-dir <dir>] [--output-dir <dir>]
```

## Flags

| Flag | Description |
|------|-------------|
| `--target <SEW_FpL>` | Target SEW name (required). Loop order `L = (F+L)/2`. Supported `L = 2..5`. |
| `--letter-projection <file\|identity\|divergent\|finite>` | Letter-slot projection (required — user-selectable). A path (e.g. `output/collinear/colprojdiv_w1.wxf`) projects each 11-dim letter slot to a lower-dim subspace; `identity` is the do-nothing value (solve in full letter space); the sentinels `divergent` / `finite` are letter-space support filters (keep entries with any divergent letter / all letters finite, derived from the nonzero rows of `colprojdiv.wxf` in `--data-dir`; applied to both the seed and the RHS; dimensions unchanged). Relative paths resolve against the executable directory. |
| `--data-dir <dir>` | Data directory with seed files. Default: `<exec_dir>/data`. |
| `--output-dir <dir>` | Output directory. Default: `<exec_dir>/output`. |
| `-h` / `--help` | Print usage. |

## How it works (per loop order `L`)

1. **Load `E1`** from `data/E1.wxf` (rank 2, `11 × 11`, 5 nnz).
2. **Recursively compute `E[2..L-1]` and `R[2..L-1]`** by loading from
   `output/{l}loop/` if they already exist, or by invoking itself (this
   is the recursive part — it ensures the lower-loop prerequisites are
   present).
3. **Compute the boundary** (`compute_rhs.hpp::compute_boundary`):
   - `E1^n` via `shuffle_power` (sequential shuffle product).
   - `E1 · R_L` via `tensor_shuffle_product_parallel` (sequential).
   - Weighted sum via `tensor_add_weighted`.
4. **Invoke `--solve-collinear`** as a subprocess:
   `./bootstrap --solve-collinear --target SEW_<name> --rhs <boundary_file>
   --projection divergent --letter-projection <abs|identity|divergent|finite>
   --data-dir <abs> --output-dir <abs>`.
   This delegates the full collinear solve (projection chain + letter
   projection + matching + linear solve) to the collinear solver, which
   writes `solMHV_LL.wxf` to `output/<L>loop/`. See
   [04_collinear_solving.md](04_collinear_solving.md).
5. **Read `solMHV_LL.wxf`** and compute `hepMHV_LL` =
   `contract(solMHV_LL, SEW_basis, axis 1, 0)`. Per the design decision
   (Q3): use the SEW basis directly — do **not** apply `colprojdiv` to
   `hepMHV` or `E_L`.
6. **Expand `hepMHV_LL` to `E_L`** (rank `2L`, dims `11^2L`) via
   `expand_hepmhv`.
7. **Compute `R_L = E_L - boundary`** and save.
8. **Verify `R_L` is divergent-free** — three branches matching the
   `--letter-projection` value:
   - `identity`: skipped (no projection to test against);
   - `divergent` / `finite`: direct support check — collect the distinct
     letter indices of `R_L`'s nonzeros and test membership against the
     divergent-letter set from `colprojdiv.wxf` (`load_divergent_letters`,
     same derivation as the solver's support filter);
   - a file path: indicator-vector method — build an 11-dim indicator of
     `R_L`'s letter indices, contract with the `--letter-projection`
     matrix; zero means `R_L = R*` is divergent-free.

## Inputs

- `data/E1.wxf` (the one-loop seed).
- All other `data/` seeds (transitively, via `--project`).
- Lower-loop outputs from `output/{l}loop/` (if they exist; otherwise
  computed recursively).

## Outputs (to `output/<L>loop/`)

Note: directory names use a digit prefix (`2loop`, `3loop`, `4loop`,
`5loop`), **not** `twoloop`/`threeloop`.

| File | Dims | Description |
|------|------|-------------|
| `solMHV_LL.wxf` | `1 × n_unknowns` | Solution coefficient vector |
| `hepMHV_LL.wxf` | `FEC_F_basis × 11` | Contracted solution (rank 2) |
| `E_LL.wxf` | `11^2L` | Expanded collinear projection (rank `2L`) |
| `R_LL.wxf` | `11^2L` | Remainder `E_L - boundary` (rank `2L`) |
| `boundary_LL.wxf` | `11^2L` | The boundary tensor (rank `2L`) |

Also writes `output/oneloop/E1.wxf` (a copy of `data/E1.wxf`) on the
first run, and triggers writes to `output/collinear/` via `--project`.

## Key files

- `compute_rhs.cpp` — CLI parsing, loop-order derivation, path resolution.
- `compute_rhs.hpp` — `compute_rhs_for_loop`, `compute_boundary`,
  `shuffle_power`, `expand_hepmhv`, `ensure_fec_tensors`,
  `ensure_sew_basis`, `find_dlogmat`.
- `tensor_shuffle.h` — `tensor_shuffle_product_parallel` (sequential
  variant used for boundary computation; the parallel variant is
  incorrect).
- `solve_collinear.hpp` — the collinear solver invoked as a subprocess
  (see [04_collinear_solving.md](04_collinear_solving.md)).

## Conventions

- **`--data-dir` / `--output-dir`**: default to `<exec_dir>/data` and
  `<exec_dir>/output`. Relative paths resolve against the executable
  directory (same convention as `bootstrap`).
- **`--letter-projection` is required** — there is no default. Pass
  either a file path (e.g. `output/collinear/colprojdiv_w1.wxf`), the
  literal `identity` to skip projection (solve in full letter
  space), or the sentinels `divergent` / `finite` (letter-space
  support filters; divergent-letter set from the nonzero rows of
  `colprojdiv.wxf` in `--data-dir`). `identity` is the
  **do-nothing** value. This makes `--solve-collinear` reusable for
  any collinear-like projection: the user selects the letter subspace
  instead of the code hardcoding `colprojdiv_w1`. The value is
  threaded through to the `--solve-collinear` subprocess as an
  absolute path (file case) or verbatim (sentinel case).
- **Subprocess invocation**: `compute_rhs` resolves the `bootstrap` binary
  as a sibling of its own executable (any cwd) and shells out to
  `bootstrap --solve-collinear` to solve the collinear constraint at each
  loop order. It passes absolute, shell-quoted `--data-dir` /
  `--output-dir` and `--letter-projection` to the subprocess. It also shells
  out to `bootstrap --extend`, `--sew`, `--project` to generate missing SEW
  basis files. `find_dlogmat(data_dir)` scans for `dlogmat_*.wxf`
  instead of hardcoding `dlogmat_E6.wxf`.
- **Sequential shuffle product**: always pass `pool = nullptr` to
  `tensor_shuffle_product_parallel` for boundary computation. The
  parallel variant produces incorrect results.
- **Boundary = master-equation recursion** (any-L formula, `R_1 := E1`) in
  `compute_boundary` since 2026-09-10 — no per-L branches to extend. The
  loop order is capped at `L ≤ 5` solely by the shuffle kernel's weight-11
  limit (`tensor_shuffle.h`); raising it means widening that kernel.

## Smoke test

```bash
# L=2 (SEW_3p1): computes E2, R2, boundary_2L
./compute_rhs --target SEW_3p1 --letter-projection output/collinear/colprojdiv_w1.wxf

# L=3 (SEW_5p1): computes E3, R3, boundary_3L (requires L=2 outputs)
./compute_rhs --target SEW_5p1 --letter-projection output/collinear/colprojdiv_w1.wxf
```

Note: `--letter-projection identity` is inconsistent at all `L ≥ 2` for
the `E6` example (union matching requires exact cancellation; the
boundary has divergent-letter entries that A does not cover in the full
11-dim space). Use `colprojdiv_w1` for the standard `E6` workflow.

## Verified status

- **L=2 with `colprojdiv_w1`**: boundary `E1²/2` (68 nnz before scaling),
  `E2` (rank 4), `R2` (rank 4, divergent-free). Union matching: 8
  intersection, 0 homogeneous, 0 b-only. Solution `c[0] = 8`.
- **L=3 with `colprojdiv_w1`**: boundary `E1³/6 + E1·R2` (5894 nnz),
  `E3` (rank 6, 11606 nnz), `R3` (rank 6, 10461 nnz, divergent-free —
  contains only letters `{2,3,4,5,6,7,8,9,10}`). Union matching: 32
  intersection, 0 homogeneous, 0 b-only. Solution `c[0] = -24, c[1] = 2`
  (unique, all 32 constraints verified).
- **`identity` (both L=2 and L=3)**: union matching is inconsistent
  (24 b-only at L=2; 1857 b-only at L=3) — no solution written.

## Reusing this pattern (normalization-difference recursion)

The reason this module exists is a **normalization difference between two
BDS-like quantities**: the object `E_L` and the collinear boundary obey
`1 + E = exp_⊗(Σ R_k)`, so the boundary that the collinear constraint must
match is *not* `E_L` itself but a subtraction of lower-loop remainders. Any
future quantity pair with the same shape — a "normalized" object whose
logarithm mixes a recursion — is tackled the same way:

1. **Derive the recursion once, in the shuffle word algebra**, and prove it
   at low orders (`L ≤ 6` here). The master equations `1 + E = exp_⊗(Σ R_k)`
   and `E_L = (1/L)·Σ_m m·(R_m ⊗ E_{L-m})` are the template.
2. **Express the boundary as shuffle products of already-computed
   tensors only** (`R_k ⊗ E_{L-k}`) — never new primitives. This keeps the
   flow node a compile-time macro over `tensor_ops shuf` / `tensor_add`
   steps, so it inherits the engine's caching, fingerprints and logging.
3. **Cache the intermediate ladder** (`tP<k>` / `tE<k>`): each order's
   remainder is reused by every higher order. The compiler fails loudly
   when a ladder tensor it needs is missing — silently defaulting one to
   zero once produced a wrong weight-6 RHS.
4. **Feed the boundary to the general nonhomogeneous solver**
   (skills/04) — nothing about the solve is boundary-specific.

## Pitfalls

1. **`B.dims()` returns by value**: calling `B.dims().begin()` and
   `B.dims().end()` separately creates dangling iterators. Capture into
   a local first: `auto bdims = B.dims();`. (Fixed in `tensor_shuffle.h`.)
2. **`insert_add` corruption**: the COO `insert_add` (ordered insert)
   corrupts tensor dims for large results. Use `unordered_map` +
   `push_back` instead. (Fixed in `tensor_shuffle.h`.)
3. **Divergent-subspace projection is required at `L ≥ 2`** for `E6`:
   because `E1` has divergent-letter entries. Under union matching,
   `identity` is inconsistent at all `L ≥ 2`. See
   [04_collinear_solving.md](04_collinear_solving.md).
4. **Matching must preserve the sew axis**: do not collapse `A` into
   `std::map<key, T>`. See [04_collinear_solving.md](04_collinear_solving.md).
