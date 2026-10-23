# Computation record index

Curated receipts and identities for the campaigns in this repository. Bulky
data (chains, solves, archives) stays out of git on disk; every entry below
records what exists, where, and how to restore or regenerate it. Workspace
cleanup of 2026-10-23 (manifest: `audits/workspace-cleanup-2026-10-23/`)
removed regenerable scratch (~428 GB) and moved regenerable-but-costly arrays
to `archives/` (tar.zst).

## Ten-loop 17+3 — certified result of record

- `ten-loop-17p3-certified/` — runner manifest/ledger/progress, final
  verification (exact 16+4 seam certificate, 44 origin holdouts), the
  native-17p3 chain stage identity, and the asymmetric QQ 16+4 export
  receipts. Data: `output_hexagon_ten_v27d_20261009/` (joint prime caches and
  QQ object kept; discovery checkpoints deleted — the certified object is
  final and does not resume).
- `ten-loop-17p3-campaign/` — V17 solver selection, FEC-17 generation cost
  receipt, first-native-prime completion, and the Dixon–Liu pair/triple-tail
  bonus audit (144,151→118,938 and 275,246→217,024 restrictions).
- Solver: frozen V27D/16 profile in `v27-profile/`; package
  `output_hexagon_v27_20261009/package-v27d` (kept on disk; binary hash in
  the profile).

## Shared data root

- `shared-data-root-16p4/` — FEC chain summary (w1..w17) and input identity,
  LEC/MHV summary. Tensors: `output_hexagon_16p4_20261006/FEC` (kept);
  chain-regeneration CRT caches deleted (regenerable by `irr_fec` /
  `bench/hexagon_fec_stream_v1.cpp`).

## Nine-loop 14+4 benchmark

- `nine-loop-14p4/` — cached-nine joint identity. Joint prime caches kept on
  disk; discovery/early-phase caches deleted.

## Split-height study

- `split-heights/` — rational heights of the exact symbols across splits
  ([17,3] vs [16,4]) at L=4..10 (decode jobs + comparison table).

## 18+2 split experiment (this branch, 2026-10-21/22)

- `split-18p2/` — complete record: measurement receipts (`audit/`), FEC-18
  chain identity + multiplicities + tensor checksums
  (`payload-checksums.sha256`), the first-prime stage identity, the streaming
  generator source, the comparison report (PDF + TeX), and the experiment
  README (full phase tables). Data on disk: `output_hexagon_18p2_20261021/`
  (FEC-18 tensors, first-prime joint payload, solver outputs; generation and
  solve scratch deleted). Status: exactly one 18+2 prime exists; no 18+2
  multi-prime reconstruction or certificates have been run.

## Archives (on disk only, not in git)

- `archives/hexagon-v18-v19-v23-packages.tar.zst` — superseded V18/V19/V23
  build-lineage packages. Restore: `tar -C . -I zstd -xf <archive>`.
- `archives/pentagon-last-entry-solve-20261013.tar.zst` — 4pformfactor
  pentagon last-entry modular solve arrays (105 GB raw; sources, code,
  figures, physical results and audit kept unarchived under `research/`).
