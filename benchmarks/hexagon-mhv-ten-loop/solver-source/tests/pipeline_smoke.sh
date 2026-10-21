#!/usr/bin/env bash
set -euo pipefail
# Existing compatible chains are inputs; all outputs go into a fresh directory.
[[ $# == 4 ]] || { echo 'pipeline_smoke.sh PIPELINE FEC_ROOT MHV_LEC_ROOT NEW_OUTPUT'; exit 2; }
exe=$1; fec=$2; lec=$3; out=$4
[[ ! -e $out ]] || { echo 'Use a fresh output directory'; exit 2; }
mkdir -p "$out"
args=(MHV 3 1 "$out/run" --fec-root "$fec" --lec-root "$lec" --data data)
"$exe" "${args[@]}" --cache-only --primes 1 > "$out/cold.log" 2>&1
"$exe" "${args[@]}" --cache-only --primes 1 > "$out/resume.log" 2>&1
grep -q PIPELINE_CACHE_HIT "$out/resume.log"
! grep -q BULK_MATRIX "$out/resume.log"
[[ ! -e "$out/run/SEW_3p1/MHV/bulk/qq_complete.tsv" ]]
if "$exe" "${args[@]}" --cache-only --primes 1 --prereduce left > "$out/identity_error.log" 2>&1; then echo 'Changed identity accepted'; exit 1; fi
"$exe" "${args[@]}" --export-all --conditions collinear+mrk --primes auto > "$out/export.log" 2>&1
for stage in bulk collinear collinear_mrk; do test -f "$out/run/SEW_3p1/MHV/$stage/qq_complete.tsv"; done
"$exe" "${args[@]}" --export-all --conditions collinear+mrk --primes auto > "$out/complete_resume.log" 2>&1
! grep -q BULK_MATRIX "$out/complete_resume.log"
echo 'PASS cache-only, unchanged resume, identity rejection, homogeneous/affine export and completed resume'
