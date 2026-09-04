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

The boundary at loop `L` is computed from the lower-loop results:

```
L=2: boundary = E1^2 / 2
L=3: boundary = E1^3 / 6 + E1 · R2
L=4: boundary = -E1^4 / 12 + E2^2 / 2 + E1 · R3
L=5: boundary = E1^5 / 20 - E1·E2^2 / 2 + E2·R3 + E1 · R4
```

where `E_L` is the expanded collinear projection of `hepMHV_LL`, and
`R_L = E_L - boundary_L` is the remainder (which must be divergent-free).

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
- **Subprocess invocation**: `compute_rhs` shells out to
  `./bootstrap --solve-collinear` to solve the collinear constraint at
  each loop order. It passes absolute `--data-dir` / `--output-dir` and
  `--letter-projection` to the subprocess. It also shells out to
  `./bootstrap --extend`, `--sew`, `--project` to generate missing SEW
  basis files. `find_dlogmat(data_dir)` scans for `dlogmat_*.wxf`
  instead of hardcoding `dlogmat_E6.wxf`.
- **Sequential shuffle product**: always pass `pool = nullptr` to
  `tensor_shuffle_product_parallel` for boundary computation. The
  parallel variant produces incorrect results.
- **Boundary formulas are hardcoded** for `L = 2..5` in
  `compute_boundary`. Adding `L = 6` requires extending this function.

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
