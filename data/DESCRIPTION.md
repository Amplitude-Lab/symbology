# Description of `data/` seed files

This directory holds the **problem-specific seed tensors** for the `E6` heptagon
symbol-space bootstrap. Every run of `bootstrap`, `compute_rhs`, or
`inspect_tensors` reads from here; outputs are written to a sibling `output/`
(or `output_<PROJECT>/` — see `README.md`).

All files are SparseRREF-native WXF (CSR layout). Mathematica-exported WXF
may differ at the byte level even when tensor contents are identical; use the
Mathematica roundtrip script (`wxf_roundtrip.wls`) before comparing CRC32.

The current problem is the `E6` adjoint with 11 letters. For a different
symmetry group, replace the seed files in a fresh `data_<group>/` folder
following the same naming convention — the `42`/`11` dimensions below are
specific to `E6`; the file *roles* are universal.

## File listing

| File | Dims (rank) | Layout / shape | Role |
|------|-------------|----------------|------|
| `dlogmat_E6.wxf` | `42 × 42 × 1191` (3) | RREF-reduced `E6` adjacency/integrability condition tensor | The "dlog" condition applied at every extension step. Read by `--extend` and `--sew` via `-c/--condition`. Sourced from `bootstrap_E6_archive/dlogmatE6RREF.wxf`. |
| `FEC_1.wxf` | `7 × 1 × 42` (3) | `{basis_w, basis_{w-1}, letter}` — forward expansion coefficient seed | First-step forward basis. Read by `--extend -f` to grow `FEC_2, FEC_3, ...`. Sourced from `bootstrap_E6_archive/FCC_1.wxf`. |
| `LEC_1.wxf` | `14 × 42 × 1` (3) | `{basis_w, letter, basis_{w-1}}` — backward expansion coefficient seed | First-step backward basis. Read by `--extend -l` to grow `LEC_2, LEC_3, ...`. Sourced from `bootstrap_E6_archive/LCC_1.wxf` and transposed to the new layout. |
| `colmat42.wxf` | `42 × 2` (2) | collinear seed on the FEC weight-1 space | Collinear projection seed. Multiplied as `reshape(FEC_1, {7,42}) · colmat42` to obtain `*first@w1` (the weight-1 collinear projection). Used by `--project --symmetry collinear`. |
| `cycrepmat.wxf` | `42 × 42` (2) | cyclic symmetry representation on the FEC weight-1 space | Cyclic symmetry seed. Applied as `reshape(FEC_1, {7,42}) · S · reshape(FEC_1, {7,42})^T` to obtain the weight-1 cyclic projection. Used by `--project --symmetry cyclic`. |
| `fliprepmat.wxf` | `42 × 42` (2) | flip symmetry representation on the FEC weight-1 space | Same role as `cycrepmat.wxf` but for the flip (reflection) symmetry. Used by `--project --symmetry flip`. |
| `parityrepmat.wxf` | `42 × 42` (2) | parity symmetry representation on the FEC weight-1 space | Same role as `cycrepmat.wxf` but for the parity symmetry. Used by `--project --symmetry parity`. |
| `colprojdiv.wxf` | `11 × 2` (2) | `(input, output)` letter-space projection to the divergent subspace | Weight-1 colprojdiv seed. Projects each 11-dim letter slot down to its 2-dim divergent subspace (for `E6`: letters `{0, 1}`). Convention: `colprojdiv[0,1] = colprojdiv[1,0] = 1` (it swaps the two divergent letters). Used as the seed for the collinear projection chain. |
| `colprojfin.wxf` | `11 × 9` (2) | `(input, output)` letter-space projection to the finite subspace | Weight-1 colprojfin seed. Maps letter `i` (for `i = 2..10`) to column `i-2`, preserving order. Used as the seed for the collinear projection chain. |
| `E1.wxf` | `11 × 11` (2) | one-loop collinear seed tensor | One-loop boundary `E1` for the RHS computation module. Rank-2 tensor with 5 non-zero entries (`E1[0,0]=-2`, `E1[1,1]=-2`, `E1[2,6]=-2`, `E1[3,7]=-2`, `E1[4,5]=-2`). Renamed from the archive's `coloneloop.wxf`. Read by `compute_rhs --target SEW_3p1` (and higher loop orders, which depend on `E1` recursively). |

## Roles by module

| Module | Reads from `data/` |
|--------|--------------------|
| `bootstrap --extend` / `--sew` | `dlogmat_E6.wxf`, `FEC_1.wxf`, `LEC_1.wxf` |
| `bootstrap --project` (collinear) | `dlogmat_E6.wxf`, `FEC_1.wxf`, `LEC_1.wxf`, `colmat42.wxf`, `colprojdiv.wxf`, `colprojfin.wxf` |
| `bootstrap --project` (symmetry) | `dlogmat_E6.wxf`, `FEC_1.wxf`, `LEC_1.wxf`, and the matching `cycrepmat`/`fliprepmat`/`parityrepmat` |
| `bootstrap --solve-collinear` | `colprojdiv.wxf`, `colprojfin.wxf` (via the projection chain in `solve_collinear.hpp`); the last-step letter projection is **user-supplied via `--letter-projection`** (not read from `data/` — pass `identity` to do nothing, a path such as `output/collinear/colprojdiv_w1.wxf`, or the sentinels `divergent` / `finite`, whose divergent-letter set is the nonzero rows of `colprojdiv.wxf` from `--data-dir`) |
| `compute_rhs` | All of the above, plus `E1.wxf` (the one-loop seed) |

## Notes on dimensions

- The FEC weight-1 space for `E6` is **7**-dimensional; `FEC_1` has shape `(7, 1, 42)`.
- The condition tensor's first axis is **42**; this is the dimension of the
  FEC/LEC weight-1 "flattened" space that the seed maps (`colmat42`,
  `cycrepmat`, etc.) act on.
- The letter space for `E6` is **11**-dimensional, of which **2** letters are
  divergent and **9** are finite. These numbers appear in `colprojdiv.wxf`
  (shape `11 × 2`) and `colprojfin.wxf` (shape `11 × 9`).
- For a different symmetry group, these dimensions will differ; the seed
  file *names* should be parameterized accordingly (e.g. `colmat<N>.wxf`,
  `colprojdiv_<letter_dim>.wxf`) and documented in a per-project
  `data_<PROJECT>/DESCRIPTION.md`.
