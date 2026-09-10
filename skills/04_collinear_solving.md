# Skill 04 — Collinear constraint solver

## Purpose

A general collinear-like constraint solver. It solves `c · A = boundary`
for any projection that is structurally similar to the collinear
projection (i.e. one that produces an expanded basis `A` and a boundary
of matching shape):

- `A` is the expanded SEW collinear basis (rank `2L+1`:
  `{sew_dim, 11, ..., 11}`).

- `boundary` is the expected collinear limit at loop order `L` (rank
  `2L`: `{11, ..., 11}`).

- `c` is the unknown coefficient vector (length `sew_dim`).

**The letter-space projection used in the last step is user-selectable
via** **`--letter-projection`, and this choice is important.** The flag is
required — it is not hardcoded. Pass either:

- a path to a projection matrix (e.g.
  `output/collinear/colprojdiv_w1.wxf`), which is applied to each
  11-dim letter slot of both `A` and the boundary via
  `apply_colprojdiv_slots` before matching; or

- the literal `identity`, which means "do nothing" — the solve runs
  in the full 11-dim letter space with no projection applied.

`identity` is the do-nothing value: no projection is applied, and the
constraint `c·A = boundary` is enforced exactly in the full 11-dim
letter space (union matching). For the standard `E6` collinear solve
this is **inconsistent at all** **`L ≥ 2`** because the boundary `E1^L/L!`
has entries at letter combinations involving the divergent letters
`{0, 1}` that the SEW collinear basis `A` does not cover. Use
`colprojdiv_w1` to project both sides to the divergent subspace, where
their supports coincide exactly and `c·A` cancels `boundary` exactly.

The finite part is `R* = c·A - boundary`, which must be divergent-free
when a divergent projection is used (verified by the indicator-vector
method, skipped for `identity`).

## CLI entry point

```bash
# Named target (SEW/FEC naming convention):
./bootstrap --solve-collinear --target <SEW_FpL|FEC_W> --rhs <rhs.wxf|0> \
    --projection <finite|divergent|none> --letter-projection <file|identity|divergent|finite> \
    [--basis <basis.wxf> ...] [--solver <incremental|sampled>] \
    [--data-dir <dir>] [--output-dir <dir>]

# Custom seed (any tensor, no naming convention):
./bootstrap --solve-collinear --target-basis <seed.wxf> --projection none \
    --rhs <rhs.wxf|0> --letter-projection <file|identity|divergent|finite> \
    [--basis <basis.wxf> ...] [--solver <incremental|sampled>] \
    [--data-dir <dir>] [--output-dir <dir>]

# Multi-pair mode: several {seed, rhs, letter projection} pairs, each with its
# own letter projection, plus optional pre-computed [M|r] conditions — all rows
# stacked into ONE linear solve:
./bootstrap --solve-collinear \
    (--pair <seed.wxf> <rhs.wxf|0> <letter>)... \
    [--pair-cond <cond.wxf>]... [--export-conditions] \
    [--out-stem <name>] [--basis <basis.wxf> ...] \
    [--solver <incremental|sampled>] [--data-dir <dir>] [--output-dir <dir>]
```

Multi-pair mode is mutually exclusive with the single-pair flags
(`--target`/`--target-basis`/`--rhs`/`--projection`/`--letter-projection`) —
each pair carries its own seed, rhs and letter projection. `--pair-cond`
only (no `--pair`) requires `--out-stem` for the output naming.

Two letter-projection modes:

- **Divergent projection** (`--letter-projection data/colprojdiv.wxf`): both
  sides are projected to the 2-dim divergent letter subspace first. This is
  the original MHV mode — a simplification that makes the calculation easier.

- **`identity`** **(finite part included)**: no projection at all; the constraint
  `c·A = boundary` is enforced in the full 11-dim letter space, so the finite
  part is also constrained. Use this when the divergent-only conditions are
  insufficient (e.g. rank 1/5 with null space 4). Both modes coexist;
  existing commands are unchanged.

Letter-projection values (the `--letter-projection` flag and the third
`--pair` argument) accept sentinels or a file path:

- `identity` — no projection; the solve runs in the full letter space.

- `divergent` — support filter: keep the tensor entries whose letter key
  contains **any** divergent letter (dimensions preserved). The backend
  derives the divergent-letter set itself from
  `<data-dir>/colprojdiv.wxf`; the front-end passes the string through.

- `finite` — keep the entries where **all** letters are finite (the exact
  complement of the `divergent` filter).

- any other string — a path to a `.wxf` projection matrix contracted per
  letter slot (LEGACY semantics: `data/colprojdiv.wxf` keeps only
  all-divergent keys, `data/colprojfin.wxf` only all-finite keys; the
  sentinels filter by support instead of contracting dimensions).

## Flags

| Flag                                                      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--solve-collinear`                                       | Run the collinear solver.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `--target <SEW_FpL\|FEC_W>`                               | Named target (e.g. `SEW_3p1` for 2-loop, `SEW_5p1` for 3-loop). The target basis is read from `output/collinear/<name>_basis.wxf` and the seed-space projection chain is computed/applied. Mutually exclusive with `--target-basis`.                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--target-basis <seed.wxf>`                               | **Custom seed mode**: use the given tensor file as-is (e.g. a summed NMHV expression like `E0+E23+E34`). No seed-space projection, no naming convention. Requires `--projection none`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `--rhs <rhs.wxf>` or `--rhs 0`                            | RHS path. Required — exits with code 1 if missing. `--rhs 0` means an all-zero RHS constructed in-memory.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `--projection <finite\|divergent\|none>`                  | Which projection to apply. `finite`/`divergent` select the seed-space projection for named targets; `none` is the custom-seed mode (requires `--target-basis`). Multi-pair mode does not use this flag — every `--pair` seed is used as-is (equivalent to `none`).                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `--letter-projection <file\|identity\|divergent\|finite>` | Letter-slot projection. Required (no default — user-selectable). `identity` = do nothing (solve in full letter space); `divergent` / `finite` = support filters (any-divergent / all-finite, dimensions preserved — see above); a path (e.g. `output/collinear/colprojdiv_w1.wxf`) is a projection matrix contracted into each letter slot. Relative paths resolve against the executable directory. Single-pair mode only — in multi-pair mode each `--pair` carries its own letter projection (third argument, same value grammar).                                                                                                                                                   |
| `--pair <seed> <rhs\|0> <letter>`                         | **Multi-pair mode** (repeatable). One `{seed, rhs, letter projection}` triple per pair: the seed is a concrete tensor used as-is (custom-seed semantics — no seed-space projection, no naming convention), optionally expanded with a shared `--basis` chain. Each pair may use a **different** letter projection (e.g. pair 1 `identity`, pair 2 `divergent` or `colprojdiv_w1.wxf`); the third argument accepts the same sentinels (`identity` / `divergent` / `finite`) or a file path as `--letter-projection`. Unwired/absent rhs may be given as `0` (homogeneous constraints). Mutually exclusive with `--target`/`--target-basis`/`--rhs`/`--projection`/`--letter-projection`. |
| `--pair-cond <cond.wxf>`                                  | **Multi-pair mode** (repeatable). Rank-2 `[M \| r]` condition matrix (rhs = last column) to stack as additional rows — re-ingest `cond_<stem>.wxf` files written by `--export-conditions` (from this or another flow) to combine constraints across flows. Contradictory `0 = nonzero` rows remain in exports, so reimporting cannot turn an inconsistent system into a consistent one. `--pair-cond` only (no `--pair`) requires `--out-stem`.                                                                                                                                                                                                                                                                                                                                                                      |
| `--export-conditions`                                     | **Multi-pair mode**. Write `output/collinear/cond_<stem>.wxf`: the combined non-homogeneous constraints as a rank-2 `[M \| r]` matrix (n\_unknowns+1 columns, rhs = last column), one row per stacked constraint row. Re-ingest via `--pair-cond` or the flow `cond` port.                                                                                                                                                                                                                                                                                                                                                                                                              |
| `--out-stem <name>`                                       | **Multi-pair mode**. Override the `sol_<stem>.wxf` / `cond_<stem>.wxf` naming. Auto: single pair/cond source → its stem; multiple → `stem1_xN`. **Required** when `--pair-cond` is given without any `--pair` (cond-only run: stack condition files and solve them directly).                                                                                                                                                                                                                                                                                                                                                                                                           |
| `--basis <basis.wxf>`                                     | Expansion basis file (repeatable; highest weight first). If omitted: auto-detected for named targets as `first_w{N}_basis.wxf` for weights `target_weight-1` down to 2; **no expansion** for custom seeds (use the tensor as-is — pass `--basis` explicitly to expand a compact custom seed).                                                                                                                                                                                                                                                                                                                                                                                           |
| `--solver <incremental\|sampled>`                         | Linear solver backend. `incremental` (default): samples constraints batch-by-batch with substitution sweeps and completes rank via full RREF when underdetermined. `sampled`: the legacy full-system modular RREF.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `--data-dir <dir>`                                        | Data directory with seed files (default: `<exec_dir>/data`). Resolved against the executable directory.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `--output-dir <dir>`                                      | Output directory (default: `<exec_dir>/output`). Same resolution as `--data-dir`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

## How it works

1. **Collinear projection chain** (`solve_collinear.hpp`):

   - Copy `data/colprojdiv.wxf` → `output/collinear/colprojdiv_w1.wxf`
     (and `colprojfin_w1.wxf`).

   - Recursively compute `colprojdiv_w{N}.wxf` and `colprojfin_w{N}.wxf`
     for `N = 2..F` by contracting the weight-`N-1` projection with the
     FEC weight-`N` basis.

   - Compute the SEW-level projection `colprojdiv_<sew_name>.wxf` (e.g.
     `colprojdiv_SEW_5p1.wxf`).

2. **Expand the SEW collinear basis** (`tensor_expand.hpp`):

   - Load `output/collinear/SEW_<name>_basis.wxf` (shape
     `{sew_dim, FEC_F_basis, 11}`).

   - Contract axis 1 with each `first_w{N}_basis.wxf` (highest weight
     first, down to weight 2).

   - Result: `A` with rank `2L+1`, dims `{sew_dim, 11, ..., 11}`.

3. **Project to the user-selected letter subspace** (`apply_colprojdiv_slots`):

   - The projection is **chosen by the user** via `--letter-projection`;
     it is not hardcoded. This makes the solver reusable for any
     collinear-like projection, not just the `E6` divergent one.

   - Apply the `--letter-projection` matrix to each 11-dim letter slot of
     both `A` and the boundary. The standard `E6` example is
     `colprojdiv_w1`, which projects to the 2-dim divergent subspace.

   - `A_proj` becomes `{sew_dim, 2, ..., 2}`; `boundary_proj` becomes
     `{2, ..., 2}`.

   - If `--letter-projection identity`, this step is **skipped** — the
     solve happens in the full 11-dim letter space ("do nothing").

   - A divergent projection is **required at** **`L ≥ 2`** for the `E6`
     example, because `E1` has divergent-letter entries
     (`E1[0,0] = -2`, `E1[1,1] = -2`) and the boundary `E1^L/L!`
     therefore has divergent components. Under union matching,
     `identity` is inconsistent at **all** `L ≥ 2` for `E6` because
     the boundary's entries at letter combinations involving `{0, 1}`
     are not covered by A in the full 11-dim space. Use
     `colprojdiv_w1` to project both sides to the divergent subspace,
     where the supports coincide exactly.

4. **Match positions — UNION of supports** (`A_match`, `b_match`):

   - Enforce `c·A = boundary` at **every** position where either `A` or
     `boundary` is nonzero (not just the intersection). This is the
     "exact match" semantics: after projecting both sides via
     `--letter-projection`, the projected supports should coincide and
     `c·A` cancels `boundary` exactly in the projected subspace.

   - Three cases per letter multi-index `key`:

     - **Both nonzero** (intersection): `c·A[key] = b[key]`.

     - **A nonzero, b zero** (homogeneous): `c·A[key] = 0` — enforces
       that A's extra support is annihilated by `c`.

     - **A zero, b nonzero** (b-only): `0 = b[key]` — trivially
       inconsistent. Detected before the linear solver; if any such
       position exists, the system is reported inconsistent and the
       solver is skipped.

   - Iterate `A_proj` directly (do **not** collapse into
     `std::map<key, T>` — that overwrites duplicate sew entries and
     breaks the multi-unknown system).

   - Preserve the sew axis in `A_match` so positions with both sew
     entries contribute two coefficients: `c0·A[0,key] + c1·A[1,key] = b[key]`.

   - When a divergent projection is used (`colprojdiv_w1`), the
     projected supports coincide exactly: intersection = union, with
     zero homogeneous and zero b-only positions. When `identity` is
     used, the boundary's finite components are not covered by A, so
     b-only positions appear and the system is inconsistent.

5. **Linear solve** (`linear_solve.hpp`):

   - Reshape `A_match` to `(n_unknowns, n_constraints)` and `b_match` to
     `(1, n_constraints)`.

   - Homogeneous constraints (A≠0, b=0) are included automatically:
     the solver collects "non-trivial" rows (where the A row has
     nonzero entries) and treats missing b entries as 0.

   - If too many constraints, sample `3 × n_unknowns` constraints, solve,
     then verify against all.

   - Uses modular RREF over `Z / 2^61` with reconstruction to `rat_t`.

6. **Indicator-vector verification** (skipped if `--letter-projection identity`):

   - Compute `R* = c·A - boundary` (in the **unprojected** full letter
     space — this is the finite remainder).

   - Collect the distinct letter indices appearing in `R*`'s non-zero
     entries.

   - Build an 11-dim indicator vector (1 at those indices, 0 elsewhere).

   - Project with the `--letter-projection` matrix. If the result is
     zero, `R*` is divergent-free (no entries at letters `{0, 1}`), so
     the collinear constraint is satisfied and `R*` is purely finite.

### Multi-pair mode (`--pair` / `--pair-cond` / `--export-conditions`)

Several `{seed, rhs, letter projection}` pairs, each with its **own**
letter projection (e.g. pair 1 `identity`, pair 2 the collinear divergent
projection `colprojdiv_w1.wxf`), plus optional pre-computed `[M|r]`
condition files. All rows are stacked into **one** linear system and
solved together (`solve_linear_system_incremental`; `incremental_solve.hpp`
is unchanged — the pairs only change row construction):

- **Per pair** (`build_pair_rows`): load the seed as-is (custom-seed
  semantics: no seed-space projection on axis 0 — named-target
  `--projection finite|divergent` does not apply), optionally expand with
  the shared `--basis` chain, load the rhs (or `0`), apply **that pair's**
  letter projection to every letter slot of both sides
  (`apply_colprojdiv_slots`, A first slot 1 / boundary first slot 0),
  then union-match and flatten to rank-2 rows `(n_unknowns+1)` —
  `[M_row | r]` with rhs in the last column.

- **Per cond file** (`--pair-cond`): read the rank-2 `[M|r]` matrix
  (rhs = last column) and stack its rows directly. Column count must
  match `n_unknowns + 1`.

- **Solve once**: all pair rows + all cond rows go into a single
  `(n_unknowns, R)` A + `(R,)` b. Consistency, rank and the particular
  solution are reported across the combined system.

- **Export** (`--export-conditions`): write
  `output/collinear/cond_<stem>.wxf` — exactly the stacked `[M|r]`
  (n\_unknowns+1 columns), re-ingestable later via `--pair-cond` (or the
  flow `cond` port) to combine constraints from this and other flows.

- **Naming** (`--out-stem`): `sol_<stem>.wxf` / `cond_<stem>.wxf`.
  Auto: single source → its stem; multiple sources → `stem1_xN`.
  Cond-only runs (no `--pair`) must pass `--out-stem` explicitly.

## Inputs

- `data/colprojdiv.wxf`, `data/colprojfin.wxf` (weight-1 seeds).

- `output/collinear/SEW_<name>_basis.wxf` (from `--project --symmetry collinear`).

- `output/collinear/first_w{N}_basis.wxf` (from the projection chain).

- The RHS file (passed via `--rhs`).

## Outputs

### To `output/collinear/`

| File                                                     | Description                                                                                                                                                                                                                                         |
| -------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `colprojfin_w1.wxf`, `colprojdiv_w1.wxf`                 | Copies of the `data/` seeds (skipped if present)                                                                                                                                                                                                    |
| `colprojfin_w{N}.wxf`, `colprojdiv_w{N}.wxf` (`N ≥ 2`)   | FEC-level projections                                                                                                                                                                                                                               |
| `colprojfin_<sew_name>.wxf`, `colprojdiv_<sew_name>.wxf` | SEW-level projections (e.g. `colprojdiv_SEW_5p1.wxf`)                                                                                                                                                                                               |
| `sol_<seed_name>.wxf`                                    | **Custom-seed mode only**: solution coefficient vector `(1 × n_unknowns)`, written only if the system is consistent. `<seed_name>` is the seed file's stem (e.g. `sol_cb198_79d057_E0E23E34tensor.wxf` for seed `cb198_79d057_E0E23E34tensor.wxf`). |
| `sol_<stem>.wxf`                                         | **Multi-pair mode**: same solution vector, `<stem>` = `--out-stem` or auto-derived (single source stem / `stem1_xN` / cond stems when cond-only).                                                                                                   |
| `cond_<stem>.wxf`                                        | **Multi-pair mode +** **`--export-conditions`**: the combined non-homogeneous constraints as a rank-2 `[M \| r]` matrix (rhs = last column). Re-ingest via `--pair-cond` or the flow `cond` port.                                                   |

Empty projections (0 rows) are not written.

### To `output/<L>loop/` (SEW targets only)

| File            | Description                                                                                                        |
| --------------- | ------------------------------------------------------------------------------------------------------------------ |
| `solMHV_LL.wxf` | Solution coefficient vector `(1 × n_unknowns)`, written only if the system is consistent. `L = target_weight / 2`. |

## Key files

- `solve_collinear.hpp` — `run_collinear_proj_chain`,
  `run_collinear_solver`, `apply_colprojdiv_slots`, `detect_chain_base_paths`.

- `linear_solve.hpp` — `solve_linear_system` (modular RREF + reconstruct).

- `tensor_expand.hpp` — `expand_tensor` (contracts with a basis chain).

## Conventions

- **`--rhs`** **is required.** Missing `--rhs` → exit code 1 with a helpful
  message. `--rhs 0` means an empty (all-zero) RHS tensor.

- **`--projection`** **is required** — there is no default. Pass `finite`,
  `divergent` (named targets), or `none` (custom seed via
  `--target-basis`) explicitly. `none` **requires** `--target-basis`;
  `finite`/`divergent` **require** `--target`.

- **Custom seed mode** (`--target-basis <seed.wxf> --projection none`):
  the seed tensor is used as-is — no seed-space projection chain, no
  `SEW_FpL`/`FEC_W` naming convention, no auto-detected expansion basis
  (pass `--basis` explicitly to expand a compact seed). The solver still
  applies `--letter-projection` to both the seed and the boundary, still
  uses union matching, and on success writes
  `output/collinear/sol_<seed stem>.wxf`. The seed may be any tensor of
  matching shape — e.g. a summed expression computed upstream in a flow
  (custom block output, tensor add, symmetry projection) and wired into
  the Solve Collinear `seed` port.

- **`--letter-projection`** **is required** — there is no default. Pass
  either a file path (e.g. `output/collinear/colprojdiv_w1.wxf`) or
  the literal `identity` to skip projection (solve in full letter
  space). `identity` is the **do-nothing** value: no projection is
  applied. This makes the solver reusable for any collinear-like
  projection — the user selects the letter subspace to project into,
  instead of the code hardcoding `colprojdiv_w1`. Note: under union
  matching, `identity` requires `c·A` to match `boundary` exactly in
  the full letter space, which is a **stronger** constraint than the
  projected solve — for the `E6` example it is inconsistent at all
  `L ≥ 2` because the boundary has entries at letter combinations
  that A does not cover.

- **`--data-dir`** **/** **`--output-dir`**: default to `<exec_dir>/data` and
  `<exec_dir>/output`. Relative paths resolve against the executable
  directory (same convention as `bootstrap`).

- **SEW-level naming**: `colprojdiv_<sew_name>.wxf` where `sew_name`
  already includes the `SEW_` prefix (e.g.
  `colprojdiv_SEW_5p1.wxf`). The old name `colprojdiv_w6.wxf` was
  renamed to this.

- **Slot order**: must remain in natural ascending weight order
  `[w1, w2, w3, w4]`. The `expansion_perm` function places new letters
  immediately after the FEC axis to maintain this; placing them at the
  end reverses to `[w1, w4, w3, w2]` and breaks verification.

- **`colprojdiv_w1`** **swaps** the two divergent letters: for `E6`,
  `colprojdiv[0,1] = colprojdiv[1,0] = 1`. Do not assume it is a
  projection onto the first coordinate.

- **Every run generates a log file** in `logs/` recording tensor
  dimensions and file CRC32.

## Smoke test

```bash
# L=2: solve for c such that c·A = E1²/2 (the boundary)
./bootstrap --solve-collinear --target SEW_3p1 --rhs output/2loop/boundary_2L.wxf \
    --projection divergent --letter-projection output/collinear/colprojdiv_w1.wxf

# L=3: boundary = E1³/6 + E1·R2 (requires L=2 outputs to exist)
./bootstrap --solve-collinear --target SEW_5p1 --rhs output/3loop/boundary_3L.wxf \
    --projection divergent --letter-projection output/collinear/colprojdiv_w1.wxf

# Custom seed (NMHV weight-2): solve c·(E0+E23+E34) = E1 for the divergent part.
# No naming convention — the summed seed comes from a flow (E0 + cyc·E23 + cyc·E34).
./bootstrap --solve-collinear --target-basis output/cb198_79d057_E0E23E34tensor.wxf \
    --projection none --rhs data/E1.wxf \
    --letter-projection data/colprojdiv.wxf --solver incremental

# Same seed, full letter space (identity = no projection): E0+E23+E34 equals
# the MHV symbol tensor in the collinear limit, so the finite part is also
# constrained — an NMHV-specific property fixing THREE coefficients.
./bootstrap --solve-collinear --target-basis output/cb198_79d057_E0E23E34tensor.wxf \
    --projection none --rhs data/E1.wxf \
    --letter-projection identity --solver incremental

# Multi-pair: the same two constraints solved together — pair 1 identity,
# pair 2 the divergent projection; export the combined [M|r] conditions.
./bootstrap --solve-collinear \
    --pair output/cb198_79d057_E0E23E34tensor.wxf data/E1.wxf identity \
    --pair output/cb198_79d057_E0E23E34tensor.wxf data/E1.wxf data/colprojdiv.wxf \
    --export-conditions --out-stem nmhv2pairs --solver incremental

# Cond-only: re-ingest exported conditions and solve them alone (requires --out-stem).
./bootstrap --solve-collinear --pair-cond output/collinear/cond_nmhv2pairs.wxf \
    --out-stem nmhv2cond --solver incremental
```

For the full recursive workflow (computing the boundary too), use
`compute_rhs` instead — see [05\_compute\_rhs.md](05_compute_rhs.md).
`compute_rhs` internally invokes `--solve-collinear` for each loop order.

## Verified status

- **L=2 (`SEW_3p1`)** with `colprojdiv_w1`: union matching reports
  8 intersection, 0 homogeneous, 0 b-only; unique solution `c[0] = 8`,
  all 8 constraints verified. `R2` is divergent-free.

- **L=3 (`SEW_5p1`)** with `colprojdiv_w1`: union matching reports
  32 intersection, 0 homogeneous, 0 b-only; unique solution
  `c[0] = -24, c[1] = 2`, all 32 constraints verified. `R3` contains
  only letters `{2,3,4,5,6,7,8,9,10}` (no divergent letters `{0,1}`),
  so `R3 = R*` is divergent-free.

- **`identity`** **at L=2**: union matching reports 44 intersection,
  87 homogeneous, 24 b-only → system is **inconsistent** (24 positions
  where `boundary ≠ 0` but `A = 0`). This is expected: the boundary
  `E1²/2` has divergent-letter entries that A does not cover in the
  full 11-dim space.

- **`identity`** **at L=3**: union matching reports 4037 intersection,
  7569 homogeneous, 1857 b-only → system is **inconsistent** (1857
  b-only positions).

- **Custom seed (NMHV weight-2,** **`E0+E23+E34`** **summed in a flow)** with
  `data/colprojdiv.wxf`: rank 1/5, unique leading coefficient
  `c[0] = -1`, null space 4, all constraints verified. Only seed
  component 0 has divergent support after letter projection — the
  solution picks exactly that component. Output
  `output/collinear/sol_cb198_79d057_E0E23E34tensor.wxf`.

- **Custom seed (NMHV weight-2) with** **`--letter-projection identity`**
  (full 11-dim letter space, finite part included): union matching
  reports 5 intersection, 4 homogeneous, **0 b-only — the system is
  consistent**, unlike the MHV `identity` case. 121 constraints
  (9 non-trivial), rank 3/5, particular solution
  `c[0] = 1, c[1] = 1, c[2] = 2`, null space 2. The full-letter-space
  conditions fix three of the five coefficients — strictly stronger
  than the divergent-only mode (rank 1/5, null space 4). This is the
  expected NMHV property: `E0+E23+E34` equals the MHV symbol tensor in
  the collinear limit, so the identity projection should (and does) fix
  three free parameters. Verified both via CLI and via the
  `NMHVw2collinear` flow in its earlier single-pair form (identity
  only) — identical results; the flow now runs the full two-pair
  system recorded below.

- **Multi-pair (NMHV weight-2, identity + divergent)**: pair 1
  (`E0+E23+E34`, `E1`, `identity`) and pair 2 (same seed/rhs,
  `data/colprojdiv.wxf`) stacked: 11 rows × 5 unknowns, rank 3/5,
  particular solution `c[0] = 1, c[1] = 1, c[2] = 2`, null space 2 —
  the union of both modes' conditions, consistent, matching the identity
  pair's fix of three coefficients. `--export-conditions` wrote
  `cond_nmhv2pairs.wxf` (11 × 6 `[M|r]`).

- **Cond round-trip**: re-ingesting `cond_nmhv2pairs.wxf` via
  `--pair-cond` alone (`--out-stem nmhv2cond`) reproduces the same
  solution `c = {1, 1, 2}`, rank 3/5 — the exported conditions carry
  the full combined information, and combining them with fresh pairs
  composes correctly.

- **`NMHVw2collinear`** **flow walkthrough (heptagonNMHV project, flow
  `b8f299fe`, run** **`2c7648e27513`, 2026-09-03)** — the full two-pair
  system, giving the **unique** NMHV weight-2 collinear solution:

  - Pair 1: seed `E0+E23+E34` (custom block: cyc-rep of `E12pre1`
    with S/S², colmat42 to each, sum with weights 1,1,1), RHS `E1`,
    letter projection `identity` (full 11-dim letter space) — 5
    intersection, 4 homogeneous, 0 b-only (consistent; the NMHV
    seed covers the full-space boundary, unlike the MHV case).

  - Pair 2: seed `E47−E67` (S³ of `E14pre1` − S⁵ of `E12pre1`, both
    through colmat42), RHS `hep1LE47mE67`, letter projection
    `divergent` (sentinel; divergent letters `{0, 1}` from
    `data/colprojdiv.wxf`) — the support filter applies to **both**
    sides of the pair: seed 21/24 entries kept, RHS 0/2 kept (both
    RHS entries are finite-letter keys, so after the divergent
    projection the pair-2 RHS is all-zero) → 10 homogeneous
    constraints `c·(E47−E67)₍div₎ = 0`, 0 b-only.

  - Stacked system 19 rows × 5 unknowns, incremental solver: rank
    **5/5**, null space 0 → **unique solution**
    `c = {1, 1, 2, −1, 1}` (`sol_nmhvsolw2.wxf` is 1×5, nnz=5).
    **Correction (2026-09-09):** the originally recorded solution
    `{1, 1, 2, 0, 1}` was wrong — the incremental solver silently
    mis-extracted it from a not-fully-reduced RREF basis (basis rows
    valid, read-off wrong; see the solver-bug note in
    `incremental_solve.hpp`). The corrected solver's answer matches
    Wolfram `LinearSolve` on the exported `cond` matrix exactly, and
    `tP_2 = tE_2 − (tP_1 ⊗ E_1)` now has a vanishing divergent part,
    the consistency the wrong `c[3] = 0` broke.

  - Physics check: neither pair is rank-5 alone — pair 1 (identity,
    inhomogeneous) fixes three coefficients `c[0]=1, c[1]=1,
    c[2]=2` (rank 3/5), pair 2 (divergent, homogeneous from a
    different seed combination) kills the remaining null space,
    fixing `c[3]=−1, c[4]=1`. Complementary constraint families
    intersecting trivially. The pair-2 RHS projecting to zero is
    correct, not a bug: the `hep1L` one-loop symbol tensor has only
    finite-letter support, so the equation `c·(E47−E67)₍div₎ = 0`
    carries no inhomogeneous information in the divergent subspace.

## What is E6-specific vs general (universality contract)

The solver is a **general nonhomogeneous-constraint solver**. Whenever a
new quantity produces constraints of the form `c · A = b` (coefficients
times known tensors against a fixed right-hand side — collinear limits,
soft limits, any normalization-difference boundary), it can be fed
through `--solve-collinear` unchanged:

| Piece | Status | Notes |
|---|---|---|
| `solve_linear_system_incremental` | fully general | any exact nonhomogeneous system; rows batched, arithmetic always exact |
| Union matching | fully general | enforces `c·A = b` on the union of supports; b-only ⇒ loud inconsistency |
| `--pair` / `--pair-cond` / `[M\|r]` round-trip | fully general | stack constraint families sharing the same unknown vector `c` |
| `--target-basis … --projection none` | fully general | custom seed = any tensor whose axis 0 indexes `c` |
| `--letter-projection` | general mechanism | the *mechanism* (file / `identity` / `divergent` / `finite`) is generic; which choice is consistent is problem-dependent |
| Projection chain (`colprojdiv/colprojfin`) | **E6-specific data** | lives in `data/colproj*.wxf`, derived from that project's alphabet; other projects supply their own seeds |
| Divergent-letter set | **project-specific data** | auto-derived from the nonzero rows of `data/colprojdiv.wxf`, never hardcoded |
| Alphabet size, letter count | general | derived from file dims (the historical hardcoded 11-dim indicator in `compute_rhs` was fixed 2026-09-09) |

Recipe for a new problem class: provide `A` (seed tensor, axis 0 = unknown
index), `b` (RHS, same letter slots), and a letter-projection choice; pick
`identity` first, and if union matching reports b-only positions, project
both sides to the subspace where the supports coincide (that is what made
`E6` solvable). Export the stacked `[M|r]` matrix (`--export-conditions`)
to cross-check in Wolfram `LinearSolve`.

## Pitfalls

1. **Map overwrite bug**: never collapse `A` into `std::map<key, T>`
   when `sew_dim > 1`. Iterate `A` directly and preserve the sew axis.
2. **Union matching requires exact cancellation**: the constraint
   `c·A = boundary` is enforced at the **union** of A's and boundary's
   supports. Positions where `boundary ≠ 0` but `A = 0` make the system
   trivially inconsistent (0 = nonzero). For the `E6` example, this
   means `--letter-projection identity` is inconsistent at all `L ≥ 2`
   — use a divergent projection (`colprojdiv_w1`) so both sides are
   projected to a subspace where their supports coincide.
3. **`expansion_perm`** **slot order**: new letters go immediately after
   the FEC axis, not at the end.

