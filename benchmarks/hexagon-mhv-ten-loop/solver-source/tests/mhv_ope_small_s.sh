#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${INPUT_SEW:?point to an MHV L6 collinear_mrk QQ object with basis identity}"
: "${FEC_ROOT:?}" "${LEC_ROOT:?point to the MHV chain}" "${REFERENCE_SEW:?point to the known unique L6 object}"
OUT=${TEST_OUTPUT:-output/test_mhv_small_s_$(date +%Y%m%d_%H%M%S)_$$}
test ! -e "$OUT"
THREADS=${THREADS:-4} bin/irr_conditions MHV 10 2 "$INPUT_SEW" "$OUT" --conditions ope_small_s --fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT" --data data --primes auto
for block in pe_A1 pe_A2 pe_E po_A1 po_A2 po_E; do
    cmp "$OUT/collinear_mrk_ope_small_s/irrSEW_10p2_${block}_origin.wxf" "$REFERENCE_SEW/irrSEW_10p2_${block}_origin.wxf"
done
echo 'MHV SMALL-S L6 EXACT COMPARISON PASS'
