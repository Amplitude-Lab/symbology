# Skill 02 — Universal projection matrix module

## Purpose

Compute the projection matrices that reduce a target tensor (FEC, LEC, or
SEW) onto a symmetry-invariant subspace. Supports two kinds of
projections:

- **Collinear** (dimension-shrinking): projects each 11-dim letter slot
  onto its 2-dim divergent subspace. Stores the reduced basis
  (`first_w{N}_basis.wxf`, `SEW_{name}_basis.wxf`) for use by the
  collinear solver.
- **Symmetry** (dimension-preserving): projects onto the invariant
  subspace of a discrete symmetry (cyclic, flip, parity). The basis is
  square, so no `_basis` file is stored.

## CLI entry point

```bash
./bootstrap --project --symmetry <collinear|cyclic|flip|parity> --target <name>
```

## Flags

| Flag | Description |
|------|-------------|
| `--project` | Run the universal projection pipeline. |
| `--symmetry <collinear\|cyclic\|flip\|parity>` | Symmetry name. |
| `--target <SEW_FpL \| FEC_W \| LEC_W>` | Target name (e.g. `SEW_5p1`, `FEC_3`, `LEC_2`). |

## How it works

1. **Read seeds** from `data/`:
   - `FEC_1.wxf`, `LEC_1.wxf` (the weight-1 expansion seeds).
   - `dlogmat_E6.wxf` (condition tensor — needed if the target's FEC/LEC
     chain must be generated).
   - The symmetry-specific seed:
     | Symmetry | Seed file | Seed shape |
     |----------|-----------|------------|
     | collinear | `colmat42.wxf` | `42 × 2` |
     | cyclic | `cycrepmat.wxf` | `42 × 42` |
     | flip | `fliprepmat.wxf` | `42 × 42` |
     | parity | `parityrepmat.wxf` | `42 × 42` |

2. **Compute weight-1 projections** `*first@w1` and `*last@w1` from the
   seeds:
   - Collinear: `*first@w1 = reshape(FEC_1, {7,42}) · colmat42`
   - Symmetry: `*first@w1 = reshape(FEC_1, {7,42}) · S · reshape(FEC_1, {7,42})^T`

3. **Recursively extend** the projections weight-by-weight up to the
   target's FEC/LEC weights (using the condition tensor to advance).

4. **For SEW targets**, additionally compute the SEW-level projection by
   contracting the FEC and LEC projections.

5. **For collinear only**, store the reduced basis at each weight
   (`first_w{N}_basis.wxf`, etc.) so the collinear solver can expand
   `hepMHV` and `A` back to full rank.

## Inputs (from `data/`)

`dlogmat_E6.wxf`, `FEC_1.wxf`, `LEC_1.wxf`, plus the symmetry-specific
seed listed above. Also reads previously-written `first_w{N-1}.wxf` /
`last_w{N-1}.wxf` from `output/<symmetry>/` (re-using prior work).

## Outputs (to `output/<symmetry>/`)

| File | When written |
|------|--------------|
| `first_w1.wxf`, `last_w1.wxf` | Always |
| `first_w{N}.wxf`, `last_w{N}.wxf` (`N = 2..F`) | Always |
| `first_w{N}_basis.wxf`, `last_w{N}_basis.wxf` | Collinear only |
| `SEW_{name}.wxf`, `SEW_{name}_basis.wxf` | SEW targets only |
| `summary.txt` | Always — tabular record of all files written this run, with dims/nnz/CRC32 |

## Key files

- `projection.hpp` — `run_projection_pipeline(symmetry, target, data_dir, output_dir, F, opt)`.
  Reads seed tensors from `data_dir` and writes to `output_dir/sym.name`.
  `bootstrap.cpp` resolves `--data-dir` / `--output-dir` (relative to the
  executable directory) and passes them in; defaults are `<exec_dir>/data` and
  `<exec_dir>/output`.
- `run_projection.sh` — driver script that loops over all four symmetries.
  Honors `PROJECT=<name>` (uses `data_<name>/` + `output_<name>/`).

## Conventions

- **Symmetry directory naming**: `output/collinear/`, `output/cyclic/`,
  `output/flip/`, `output/parity/`.
- **CSR-only**: all tensor functions return CSR. COO overloads were
  removed (see project memory).
- **Empty projections are not written**: if a projection has 0 rows, the
  WXF writer cannot serialize it (the WXF format requires `nnz > 0`).
- **Seed aliasing**: the collinear chain copies `data/colprojdiv.wxf`
  and `data/colprojfin.wxf` to `output/collinear/colprojdiv_w1.wxf` and
  `output/collinear/colprojfin_w1.wxf` (skipped if already present).
  The CLI also accepts the sentinels `divergent` / `finite` as
  letter-projection values (`--letter-projection` / the third `--pair`
  argument): they filter entries by support — any divergent letter /
  all letters finite — instead of contracting each letter slot with the
  matrix, while passing these file paths keeps the legacy per-slot
  contraction semantics (dims reduced).

## Smoke test

```bash
./bootstrap --project --symmetry collinear --target SEW_5p1
./bootstrap --project --symmetry cyclic --target SEW_5p1
```

Or via the driver script:

```bash
./run_projection.sh SEW_5p1   # runs all four symmetries
```

## Verified status

- `SEW_5p1` collinear: basis dims `2 × 546 × 11`, expansion to rank 7
  (`2 × 11^6`) verified.
- `SEW_5p1` cyclic/flip/parity: all three are identity `diag(1, 1)` for
  `E6` — the symmetry files are byte-identical. Only the collinear
  transformation is non-identity (`diag(-4, -2)`).
