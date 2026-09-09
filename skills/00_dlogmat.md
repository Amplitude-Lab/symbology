# Skill 00 — Alphabet → dlogmat generation

## Purpose

Produce the condition tensor `dlogmat_E6.wxf` (the dlog matrix) together with
the forward/backward seeds `FEC_1.wxf` and `LEC_1.wxf`. These are the **only
inputs the bootstrap pipeline needs from outside the C++ codebase** — every
downstream module (`01_bootstrap` onward) consumes them. This is therefore
skill 0: it runs once, by hand in Mathematica, before any C++ step.

## Entry point

### Primary: `SymbolBootstrap.wl` (project root)

The standard package [SymbolBootstrap.wl](../SymbolBootstrap.wl) provides a
clean, object-oriented API with caching. It auto-loads SparseRREF from the
project root. See [walkthrough.wl](../walkthrough.wl) for a full E6 example.

```mathematica
$ProjectRoot = "/path/to/symbology";
SetDirectory[$ProjectRoot];
Get[FileNameJoin[{$ProjectRoot, "SymbolBootstrap.wl"}]];

DeclareAlphabet["E6", alphabetE6];
SetAlphabetExpression["E6", alphabetExprE6];
SetExtendedSteinmann["E6", nonadjpairE6ES];
SetClusterAdjacency["E6", adjpairE6CA];
SetFirstEntry["E6", {a11,...,a17}];
SetLastEntry["E6", {a21,...,a37}];

dlogmatE6 = GetAlphabetConditionTensor["E6",
  {"Integrability", "Extended Steinmann", "Cluster Adjacency"}];
FEC[1] = GetFirstEntryTensor["E6"];
LEC[1] = GetLastEntryTensor["E6"];
```

Key advantages:
- **Algebraic alphabets handled natively**: `GenSqrtD` detects square roots,
  Z2-reduces them to an independent set, rationalizes denominators, and
  `GenIntRelMat` keeps `sq[i]` symbolic in the minors, extracting rational
  coefficients via `CoefficientArrays`. No float pipeline or kinematic
  sampler needed.
- **Automatic denominator non-vanishing**: `GenIntRelMat` computes
  `den$lcm = Times @@ Union @@ Flatten[...]` (the product of all distinct
  denominator factors) and retries sampling while `den$lcm /. numrule === 0`.
  No separate `"Constraint"` option needed.
- **No `NotebookDirectory[]`**: uses `$InputFileName` / `Directory[]`.

### Compiling SparseRREF (one-time, macOS)

`SparseRREF.wl` calls `sprreflink.dylib` via LibraryLink. The dylib must be
compiled once; it is **not** shipped in the repo. Build with Homebrew
FLINT/GMP/MPFR and the Mathematica C headers:

```bash
cd SparseRREF
g++-14 sprreflink.cpp -fPIC -shared -O3 -std=c++20 -o sprreflink.dylib \
  -I/opt/homebrew/opt/flint/include -I/opt/homebrew/opt/gmp/include \
  -I/opt/homebrew/opt/mpfr/include \
  -I/Applications/Mathematica.app/Contents/SystemFiles/IncludeFiles/C \
  -L/opt/homebrew/opt/flint/lib -L/opt/homebrew/opt/gmp/lib \
  -L/opt/homebrew/opt/mpfr/lib -lflint -lgmp -lmpfr -ltbb
```

`-ltbb` is required because `sparse_tensor.h` uses the C++20 parallel STL.
Without it, linking fails with `tbb::detail::r1::allocate` undefined
symbols.

Smoke check (returns a 2×3 RREF and a 3×1 kernel):

```mathematica
Needs["SparseRREF`", FileNameJoin[{$ProjectRoot, "SparseRREF", "SparseRREF.wl"}]];
SparseRREF[SparseArray @ {{1,0,2},{1/2,1/3,1/4}},
  "OutputMode" -> "RREF,Kernel", "Method" -> "Right",
  "BackwardSubstitution" -> True, "Threads" -> 0]
```

## Inputs

| Input | Symbol | Description |
|-------|--------|-------------|
| Alphabet | `alphabetE6` | List of `n` letters, each `Log[expr]` where `expr` is a rational function of kinematic variables (may contain square roots). For E6, `n = 42`. |
| Twistor matrix | `MatE6` | `k × 4` momentum-twistor matrix parameterized by face variables `f1..f6`. |
| Bracket rewrite | `rewrite` | Rules `aij -> (ratio of 4-brackets ab[...])`, with `cap1[...]` for square-root letters. |
| First-entry list | `{a11,...,a17}` | 7 letters allowed as the first symbol of a word. |
| Last-entry list | `{a21,...,a37}` | 14 letters allowed as the last symbol of a word. |
| ES non-adjacent pairs | (optional) | Pairs forbidden by extended Steinmann; fed to `GenDlogmatES`. |
| CA adjacent pairs | (optional) | Pairs allowed by cluster adjacency; fed to `GenDlogmatCA`. |

The alphabet is first converted to momentum-twistor rational functions:

```mathematica
alphabetExprE6 = Simplify[alphabetE6 /. rewrite /. abrules[MatE6]];
```

where `abrules[mat]` rewrites each `ab[a,b,c,d]` as `Det[mat[[{a,b,c,d}]]]`.

## Outputs (to `data/`)

| File | Dims (E6) | Contents |
|------|-----------|----------|
| `dlogmat_E6.wxf` | `{42, 42, 1191}` | RREF-reduced condition tensor: antisymmetric two-form coefficients in the free dlog basis. |
| `FEC_1.wxf` | `{7, 1, 42}` | Forward seed: one-hot rows over the 7 first-entry letters. |
| `LEC_1.wxf` | `{14, 42, 1}` | Backward seed: one-hot rows over the 14 last-entry letters. |

All three are exported with `Export[..., .wxf]` (WXF = CSR-native).

## Procedure (the core: `GenDlogmatInt`)

`GenDlogmatInt[alphabetExpr, opts]` in `SymbolBootstrap.wl` is the
integrability-based generator. It runs three stages:

### 1. Square-root-aware dlog vector — `GenSqrtD`

For each letter, compute `d log(letter) = D[letter]/letter` expanded in the
`D[var]` basis, with square roots handled algebraically:

- `vars = Variables[alphabetExpr]` — independent variables.
- Rewrite each `Power[a_, k_]` with `IntegerQ[k-1/2]` as
  `Power[Factor[a], Floor[k]] * sqrt[Factor[a]]`, then simplify products of
  `sqrt` via the `sqrt/:sqrt[x_]^k_Integer` UpValues.
- **Z2 row-reduce** the square-root set: find independent irrational
  polynomials, identify conjugate pairs (`p` and `-p`), and reduce the
  square-root basis via `NullSpace[..., Modulus->2]`.
- **Rationalize denominators**: multiply numerator and denominator by the
  conjugate of the `sqrt` factors in the denominator.
- Returns `{dlogexpr, vars, sqrtlist}` where `dlogexpr` is an `n × #vars`
  matrix rational in `vars` and `sqrtlist`.

### 2. Integrability relations — `GenIntRelMat`

The two-form `d log(W_a) ∧ d log(W_b)` must be closed; this gives linear
relations among the two-form coefficients.

- `subset = Subsets[Range[Length[vars]], {2}]` — all variable pairs.
- For each letter pair `{a,b}` and each variable pair, take the `2×2` minor
  `Det[dlogexpr[[{a,b}, {var_i, var_j}]]]`.
- **Algebraic extension**: `sq[i]` placeholders are kept symbolic in the
  minor. `CoefficientArrays[#, Array[sq, ...]]` extracts the rational
  coefficients — each algebraic minor becomes multiple rational rows (one
  per `sq[i]` power). This is the key mechanism for handling non-Galois-
  invariant alphabets like the pentagon.
- **Numeric sampling**: draws `RandomPrime[{2, 3*n}]` for the variables,
  retrying while `den$lcm /. numrule === 0` (where `den$lcm` is the product
  of all distinct denominator factors — prevents vanishing denominators
  automatically). Samples `n$samp ≈ 10 + Ceiling[Binomial[n,2]/Binomial[#var,2]]`
  points.
- `SparseRREF[mat, "Method"->"Right", "Threads"->...]` row-reduces the
  numeric system.

#### `GenIntRelMat` / `GenDlogmatInt` options

| Option | Default | Purpose |
|--------|---------|---------|
| `"Samples"` | `Automatic` (`10 + Ceiling[Binomial[n,2]/Binomial[#var,2]]`) | Number of numeric sampling points. Raise for large alphabets. |
| `"Tries"` | `100` | Max retries when `den$lcm` vanishes at a sampled point. |
| `"Threads"` | `0` | Thread count forwarded to `SparseRREF`. |
| `"Verbose"` | `True` | Print progress / show progress indicator. |

### 3. Assemble the dlog tensor

- `upper` = sparse `{n, n, Binomial[n,2]}` encoding the upper-triangular
  two-form basis `S[i,j]` (`i < j`).
- `result = (upper - Transpose[upper, {2,1,3}]) . Transpose[CanonSparseArray[rref]]`
  — contracts the RREF kernel with the antisymmetric two-form structure.
- Returns sparse `{n, n, n_basis}` where `n_basis` is the rank of the free
  two-form space.

### Combining condition tensors — `CombineConditionTensor`

When ES and/or CA constraints are added (via `GetAlphabetConditionTensor[name, {conditions}]`),
the tensors are combined:

- `Join[tensors, 3]` → `Transpose[{2,3,1}]` → flatten to `{1, n^2 * #tensors * dim3}`.
- `SortSparseArray` (by nonzero count) → `RowRescale` → `SparseRREF` →
  `CanonSparseArray` (drop zero rows) → `RowRescale` → `Transpose` → reshape
  back to `{n, n, _}`.

## Optional condition tensors

Two extra constraint families can be folded in before the final RREF. In
`SymbolBootstrap.wl`, set them via `SetExtendedSteinmann[name, pairs]` and
`SetClusterAdjacency[name, pairs]`, then retrieve via
`GetAlphabetConditionTensor[name, {"Integrability", "Extended Steinmann", "Cluster Adjacency"}]`.

- **Extended Steinmann** — `GenDlogmatES[alphabet, nonadjpairs]`: builds the
  complementary adjacent-pair list (`Complement[Tuples[alphabet,{2}],
  nonadjpairs]`) and delegates to `GenDlogmatCA`. Forbids the listed
  non-adjacent two-forms.
- **Cluster adjacency** — `GenDlogmatCA[alphabet, adjpairs]`: maps the
  adjacent pairs to letter ids, forms `S@@@adjpairIds`, and takes the null
  space of `CoefficientArrays[..., Flatten[Array[S,{n,n}],1]][[2]]`. Returns
  a `{n,n,_}` tensor.
- **Combine** — `CombineConditionTensor[t1, t2, ...]`: checks that all inputs
  share the same first two dims (`{n,n}`), stacks them along axis 3, flattens
  to `{1, n^2 * #tensors * dim3}`, sorts rows by nonzero count, rescales,
  RREFs, canonicalizes, rescales again, and reshapes back to `{n,n,_}`.

## FEC / LEC seeds

In `SymbolBootstrap.wl`, set via `SetFirstEntry[name, letters]` and
`SetLastEntry[name, letters]`, then retrieve via `GetFirstEntryTensor[name]`
and `GetLastEntryTensor[name]`.

```mathematica
SetFirstEntry["E6", {a11,a12,a13,a14,a15,a16,a17}];
SetLastEntry["E6", {a21,...,a37}];
FEC[1] = GetFirstEntryTensor["E6"];  (* {7, 1, n} *)
LEC[1] = GetLastEntryTensor["E6"];   (* {14, n, 1} *)
```

- `GenFEC[alphabet, firstentry]`: sparse `{len(firstentry), 1, n}` with
  `1` at `{i, 1, firstentryIds[[i]]}`. One-hot over the allowed first
  letters.
- `GenLEC[alphabet, lastentry]`: sparse `{len(lastentry), n, 1}` with
  `1` at `{i, lastentryIds[[i]], 1}`. One-hot over the allowed last letters.

Both use `Dispatch[Thread[alphabet -> Range[n]]]` for fast letter→id lookup.

## Conventions

- **Letter indexing**: `W[i]` runs `1..n` in the order letters appear in
  `alphabetE6`. This order is fixed for the whole pipeline — `FEC_1`,
  `LEC_1`, and `dlogmat_E6` must all use the same ordering.
- **Antisymmetry**: `dlog[{W[i],W[j]}] = -dlog[{W[j],W[i]}]`; the diagonal
  `dlog[{W[i],W[i]}] = 0`. Only the `i <= j` half is stored in `basis`.
- **Basis = free two-forms**: the third axis of `dlogmat` enumerates the
  `dlog[{W[a],W[b]}]` with `a <= b` that are **not** eliminated by the
  integrability/ES/CA relations. Its size is the rank of the free two-form
  space (1191 for E6, 361 for the pentagon).
- **Output format**: CSR WXF, RREF-reduced, primitive-integer rows. Matches
  the CSR-native convention every C++ module expects.
- **Path**: the notebook writes `dlogmat_E6.wxf` etc. into its own directory;
  copy them into `data/` (or `data_<PROJECT>/`) for the bootstrap to find.

## Algebraic alphabets

When the alphabet contains square-root letters (e.g. pentagon `W[26]-W[30]`
involving `eps5 = Sqrt[Delta5]`) whose **Galois conjugates are not also in
the alphabet**, the integrability minors are genuinely algebraic: they
contain `1/Sqrt[Delta5]` factors that do not cancel. `SymbolBootstrap.wl`
handles this natively via `QQ[sqrt]` coefficient extraction:

- `GenSqrtD` detects all square roots, Z2-reduces them to an independent
  set, and rationalizes denominators.
- `GenIntRelMat` keeps `sq[i]` symbolic in the minors and extracts rational
  coefficients via `CoefficientArrays[#, Array[sq, ...]]`. Each algebraic
  minor becomes multiple rational rows (one per `sq[i]` power).
- `den$lcm = Times @@ Union @@ Flatten[...]` (product of all distinct
  denominator factors) is used to reject sampled points where any
  denominator vanishes — no separate `"Constraint"` option needed.

No special handling is required from the caller. Verified on the pentagon:
`{31, 31, 361}`, nnz 1754, rank 361.

## Pitfalls

1. **Vanishing denominators at sampled points.** `GenIntRelMat` computes
   `den$lcm` (the product of all distinct denominator factors in the dlog
   vector) and retries sampling while `den$lcm /. numrule === 0`, up to
   `"Tries"` (default 100) attempts. If every draw hits a zero denominator,
   it prints an error and returns `$Failed`. Raise `"Tries"` or widen the
   prime range if this happens.
2. **Under-sampling.** Default sample count is
   `10 + Ceiling[Binomial[n,2]/Binomial[#var,2]]`. For large alphabets or
   high-degree letters, raise `"Samples"` — under-sampling produces a system
   whose RREF misses relations and inflates the basis dimension.
3. **`SmartArrayReshape` allows exactly one `_`.** More than one unknown
   dimension is an error; zero unknowns bypasses inference. The final
   reshape relies on `Times@@dims` being consistent — if `CanonSparseArray`
   dropped rows, the inferred third dim shrinks accordingly (this is
   intended).
4. **`CanonSparseArray` removes all-zero rows.** After RREF some letter-pair
   rows are identically zero (no constraint); they are dropped, so the
   returned tensor may have fewer than `n^2` letter-pair rows. The reshape
   back to `{n,n,_}` re-expands with implicit zero rows — consumers must
   treat absent rows as zero, not as missing.
5. **`RowRescale` assumes integer-valued rows.** It is applied after RREF;
   if the system has near-zero entries (from numeric noise), the GCD/LCM
   logic can break. This is rarely an issue with `SymbolBootstrap.wl` since
   all computations stay exact rational.
6. **Letter order is load-bearing.** The `alphabet` list order determines
   `W[i]` indexing. If you re-order the alphabet, you must regenerate
   `FEC_1`, `LEC_1`, **and** `dlogmat` together — mixing orderings
   between these files silently corrupts every downstream tensor.
7. **`CombineConditionTensor` requires matching `{n,n}` leading dims.** Mixing
   tensors built from different alphabets (different `n`) prints an error
   and returns `$Failed`. ES and CA tensors must be built from the same
   `alphabet` list as the integrability tensor.

## Smoke test (E6)

See [walkthrough.wl](../walkthrough.wl) for the full E6 example. Summary:

```mathematica
$ProjectRoot = "/path/to/symbology";
SetDirectory[$ProjectRoot];
Get[FileNameJoin[{$ProjectRoot, "SymbolBootstrap.wl"}]];

DeclareAlphabet["E6", alphabetE6];
SetAlphabetExpression["E6", alphabetExprE6];
SetExtendedSteinmann["E6", nonadjpairE6ES];
SetClusterAdjacency["E6", adjpairE6CA];
SetFirstEntry["E6", {a11,a12,a13,a14,a15,a16,a17}];
SetLastEntry["E6", {a21,...,a37}];

dlogmatE6 = GetAlphabetConditionTensor["E6",
  {"Integrability", "Extended Steinmann", "Cluster Adjacency"}];
FEC[1] = GetFirstEntryTensor["E6"];
LEC[1] = GetLastEntryTensor["E6"];

Export[FileNameJoin[{$ProjectRoot, "data", "dlogmat_E6.wxf"}], dlogmatE6];
Export[FileNameJoin[{$ProjectRoot, "data", "FEC_1.wxf"}], FEC[1]];
Export[FileNameJoin[{$ProjectRoot, "data", "LEC_1.wxf"}], LEC[1]];
```

A successful run prints: independent variables and square-root list from
`GenSqrtD`; sample-point count and timing from `GenIntRelMat`. The exported
`dlogmat_E6.wxf` should have dims `{42, 42, 1191}`.

## Smoke test (Pentagon)

The pentagon alphabet ([data_pentagon/alphabet.wl](../data_pentagon/alphabet.wl))
is stored as `LetterRep = {W[1] -> expr1, ..., W[31] -> expr31}` and
`RootDef = eps5 -> Sqrt[Delta5]`. The full pipeline from `alphabet.wl` to
`dlogmat_pentagon.wxf`:

```mathematica
$ProjectRoot = "/path/to/symbology";
SetDirectory[$ProjectRoot];
Get[FileNameJoin[{$ProjectRoot, "SymbolBootstrap.wl"}]];

(* 1. Load the raw alphabet definitions (LetterRep, RootDef) *)
Get[FileNameJoin[{$ProjectRoot, "data_pentagon", "alphabet.wl"}]];

(* 2. Declare the alphabet: 31 letters W[1]..W[31] *)
DeclareAlphabet["Pentagon", Table[W[i], {i, 1, 31}]];

(* 3. Set the parametrized expressions.
      Substitute eps5 -> Sqrt[Delta5] via RootDef so that the square root
      is inline (SymbolBootstrap.wl's GenSqrtD detects Power[_, 1/2]
      and handles it automatically). *)
SetAlphabetExpression["Pentagon", LetterRep[[All, 2]] /. RootDef];

(* 4. Generate the integrability condition tensor *)
dlogmat = GetIntegrabilityTensor["Pentagon"];

(* 5. Export *)
Export[FileNameJoin[{$ProjectRoot, "data_pentagon", "dlogmat_pentagon.wxf"}], dlogmat];
```

Expected output: dims `{31, 31, 361}`, nnz 1754, rank 361. The 5 algebraic
letters `W[26]-W[30]` (involving `eps5`) are handled automatically by
`GenSqrtD` + `GenIntRelMat` via `QQ[sqrt]` coefficient extraction — no
special options needed. See [data_pentagon/Description.md](../data_pentagon/Description.md)
for the full alphabet breakdown.

## Verification

1. Confirm `Dimensions[dlogmatE6]` is `{42, 42, 1191}` (E6).
2. Roundtrip `dlogmat_E6.wxf` through `wxf_roundtrip.wls` and compare CRC32
   against `data/dlogmat_E6.wxf`.
3. Sanity: the integrability tensor alone (`GetIntegrabilityTensor["E6"]`)
   must have rank consistent with `Binomial[42,2] - 1191 + (ES+CA
   contributions)`; an unexpectedly large basis dimension signals
   under-sampling in `GenIntRelMat`.
4. End-to-end: feed the regenerated `dlogmat_E6.wxf` + `FEC_1.wxf` into
   `./bootstrap --extend` (skill 01) and confirm `FEC_2.wxf` CRC32 matches
   the archive.
