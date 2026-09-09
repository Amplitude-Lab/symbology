# SparseRREF bug report: data race in `sparse_mat_rref_forward/backward` (phantom-pivot SIGSEGV)

**Upstream base tested:** `5bbee55` (merge of PR #20). The racy protocol is byte-identical on `master`@`b5b3d3c` — the fix patch applies cleanly there with no fuzz (all crash/repro data below is from the `5bbee55` build).
**Files affected:** `sparse_mat.h` — `sparse_mat_rref_forward` (\~L1344–1427) and the mirrored block in `sparse_mat_rref_backward` (\~L1682–1757)
**Severity:** intermittent crash (\~10% of runs on arm64) + formal data race (C++ UB); silent wrong-result risk in a rare interleaving
**Fix:** `sparserref-race-fix.patch` (race only, applies on `master`@`b5b3d3c`) / `sparserref-5bbee55-race-and-init-fix.patch` (race + the `sparse_type.h` placement-new fix, for `5bbee55`; the placement-new part is already upstream as `9874061`)

***

## 1. Summary

During overlapped Schur completion + transpose rebuild, worker threads publish "row rewritten"
through a plain `std::vector<int> flags` with no synchronization or memory ordering. The main
thread reads the flag and immediately re-reads the row's heap data to rebuild the transpose
index. Under the C++ memory model this is a data race (UB); on weakly-ordered arm64 hardware
it manifests as stale/torn row reads, which poison the transpose index, which yields a
**phantom pivot** `(r, c)` in the next pivot search, which null-derefs at the pivot rescale:

```cpp
T scalar = scalar_inv(*mat.find(r, c), F);   // sparse_mat.h:1259 — find() returns nullptr
```

## 2. Reproduction

Workload: symbolic-bootstrap projection (sparse RREF over QQ, `sparse_mat<unsigned long, int>`,
`rref_option.method = 2`, \~10-core machine, thread pool at hardware concurrency). Any run of
a large `--project`-style computation reproduces within dozens of runs.

Statistics on one machine (Apple Silicon, 10 cores; Homebrew GCC 14.2.0, `-O3 -flto`):

| build                                                                      | runs | crashes                              |
| -------------------------------------------------------------------------- | ---- | ------------------------------------ |
| unpatched `5bbee55`                                                        | 300  | **32 (10.7%)** — `SIGSEGV`, rc = −11 |
| patched                                                                    | 300  | **0**                                |
| patched, UBSan (`-fno-sanitize-recover`, vptr/alignment excluded — see §6) | 60   | **0**                                |

Patched runs are bit-identical to unpatched *surviving* runs (CRC32 of every output tensor
matches), i.e. the fix changes synchronization, not mathematics.

### Crash evidence

lldb (production binary, `thread #9`):

```
* thread #9, stop reason = EXC_BAD_ACCESS (code=1, address=0x0)
    frame #0: bootstrap`std::_Function_handler<void (), void BS::thread_pool<...>::detach_loop<
      ..., SparseRREF::sparse_mat_rref_forward<unsigned long, int>(...)::'lambda2'(unsigned long)>
      ...::_M_invoke(...) + 868
->  0x100029f64 <+868>: ldr    x0, [x0]     ; null load
```

The faulting frame is a `detach_loop` pool task inside `sparse_mat_rref_forward` — consistent
with the UBSan attribution below, which pinpoints the exact statement.

UBSan (independent confirmation, same line):

```
SparseRREF/sparse_mat.h:1259:27: runtime error: load of null pointer of type 'unsigned long'
SUMMARY: UndefinedBehaviorSanitizer: undefined-behavior SparseRREF/sparse_mat.h:1259:27
```

## 3. Root cause analysis

### 3.1 The racing protocol (upstream code)

In each RREF round, after new pivots are chosen, non-pivot rows are Schur-completed by pool
tasks while the main thread concurrently rebuilds the transpose index for finished rows:

```cpp
std::vector<int> flags(leftrows.size(), 0);          // plain ints
pool.detach_blocks<size_t>(0, leftrows.size(), [&](const size_t s, const size_t e) {
    for (size_t i = s; i < e; i++) {
        schur_complete_func(mat, leftrows[i], ...);  // (A) rewrites mat[row] heap data
        flags[i] = 1;                                // (B) plain store, no ordering
        ...
    }});
...
// main thread / serial loop:
if (flags[i]) {                                      // (C) plain load, no ordering
    auto row = leftrows[i];
    for (size_t j = 0; j < mat[row].nnz(); j++)      // (D) re-reads mat[row] heap data
        tranmat_vec[0][mat[row](j)].push_back(row);
    flags[i] = 0;
}
```

Problems:

1. **Data race / UB.** (B) and (C) are unsynchronized accesses to the same scalar from
   different threads — a race by definition. Moreover nothing orders (B) after (A), so the
   consumer can act on the flag while seeing stale or partially-visible row data from (A)/(D).
   On x86(TSO) this mostly survives in practice; on arm64 it fails readily (hence 10% here).

2. **Phantom pivot → null deref.** A stale row read makes the rebuilt `tranmat` disagree
   with `mat`. The next round's `pivots_search` then trusts `tranmat` and can return a pivot
   `(r, c)` for a coefficient that does not exist in `mat[r]`. The pivot-rescale loop does
   `*mat.find(r, c)` with no null check → UB → SIGSEGV at `sparse_mat.h:1259` (and the mirror
   site in `sparse_mat_rref_backward`).

### 3.2 A second protocol defect: double-indexed rows (silent-wrong-result risk)

The overlapped path also has a deterministic-logic hole independent of memory ordering. When
the main thread offloads a "second wave" of finished rows to the pool:

```cpp
std::vector<size_t> newleftrows;
for (size_t i = 0; i < leftrows.size(); i++)
    if (flags[i]) newleftrows.push_back(leftrows[i]);
pool.detach_loop(0, newleftrows.size(), [&](size_t i) {
    auto row = newleftrows[i];
    ... tranmat_vec[id][col].push_back(row, true);     // indexes row, but ...
});
pool.wait();
// serial loop then runs again:
if (flags[i]) { ...push_back(row); flags[i] = 0; localcount++; }
```

The second-wave tasks never clear the flags they consumed. After `pool.wait()`, the serial
loop re-examines those flags (still 1), pushes the same rows into `tranmat` **again**
(`sparse_vec::push_back` is a blind append — no dedup), and counts them **again** in
`localcount`. Consequences:

* duplicate `tranmat` entries are usually harmless (pivot search takes the first match per
  column), which is why surviving runs have always produced correct output;

* but the inflated `localcount` can satisfy the loop-exit condition
  (`localcount == leftrows.size()`) while some *other* rows were never indexed into `tranmat`
  at all — a latent missed-pivot / wrong-rank hazard in mixed interleavings.

The same two defects exist verbatim in `sparse_mat_rref_backward` (\~L1682–1757).

### 3.3 Related: `--verbose` wait-loop can hang

Inside the verbose branch of the overlapped path, the progress loop waits on
`localcount < leftrows.size()`, but at that point only the just-dispatched second wave
(`newleftrows.size()` of them) is in flight — if the wave completes fewer than all remaining
rows, the loop spins forever with the pool idle. It must wait on `newleftrows.size()` (the
serial loop afterwards finishes the rest).

## 4. The fix

`sparserref-race-fix.patch` (and its `5bbee55` companion) makes the flag protocol sound
while keeping the overlap. Four parts, mirrored in forward and backward passes:

1. **Atomic flags with release/acquire.** `flags` becomes
   `std::vector<std::atomic<int>>`; workers `store(1, memory_order_release)` **after**
   `schur_complete_func`, consumers `load(memory_order_acquire)` **before** reading
   `mat[row]`. This creates the happens-before edge that makes the rewritten row visible
   with the flag, and removes the race.

2. **Exact-once indexing.** Second-wave tasks now carry `(flag_index, row)` pairs and clear
   their own flag (`store(0, release)`) after indexing. Every row is therefore indexed
   exactly once and `localcount` counts each row exactly once — closing the duplicate /
   missed-row hole from §3.2.

3. **Loud pivot guards.** Both pivot-rescale sites now null-check `mat.find(r, c)` and
   `std::abort()` with a one-line diagnostic ("pivot missing — stale transpose index?").
   This converts any future regression from a bare `SIGSEGV` inside template soup into an
   actionable message. (Defensive only; never fired post-fix in 360 runs.)

4. **Verbose wait-loop** waits on `newleftrows.size()` (see §3.3).

The companion patch for `5bbee55` additionally carries the `sparse_type.h` placement-new fix
(construct `rat_t` values in `s_malloc`'ed storage before the parallel copy; destroy the old
array before `s_free`). This is independent of the race and is already upstream as `9874061`
("fix some UBs") — it is included only because our pinned base commit predates it.

## 5. How to apply

On upstream `master` (tested on `b5b3d3c`, v0.4.0):

```bash
cd SparseRREF
git apply sparserref-race-fix.patch
```

On the older `5bbee55` base (what our project pins; includes the `sparse_type.h` fix that
master already has):

```bash
cd SparseRREF
git apply sparserref-5bbee55-race-and-init-fix.patch
```

## 6. Reproducing the diagnosis yourself

* **Confirm the crash:** loop any large rref workload \~100–300 times on arm64; expect
  double-digit SIGSEGV rates. The crash log shows the fault in the pivot-rescale
  `detach_loop` of `sparse_mat_rref_forward`.

* **UBSan:** build with clang++ + `-fsanitize=undefined -fno-sanitize-recover=all`
  (exclude `vptr` and `alignment` — the vptr check misfires across the GCC-libstdc++/clang
  boundary, and one pre-existing benign misaligned read in the WXF reader would abort before
  reaching rref). The unpatched build reports
  `sparse_mat.h:1259: runtime error: load of null pointer of type 'unsigned long'`.

* **Verify the fix changes nothing mathematically:** CRC32 every output tensor before/after;
  they match (the fix only adds synchronization).

## 7. Test environment

* Apple Silicon (arm64), 10 cores (4P+6E), macOS

* Homebrew GCC 14.2.0 (`-O3 -flto -std=c++20`), linked `-lflint -lgmp -lmimalloc -ltbb`

* UBSan build: Homebrew llvm clang++ with GCC's libstdc++ (`-stdlib=libstdc++`)

* SparseRREF base `5bbee55` (race protocol identical on `master`@`b5b3d3c`)

