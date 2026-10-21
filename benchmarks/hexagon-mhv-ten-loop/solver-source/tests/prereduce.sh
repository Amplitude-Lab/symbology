#!/usr/bin/env bash
# Existing low-weight chains only; never remove or overwrite results.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
: "${FEC_ROOT:?Set FEC_ROOT to an existing first chain}"
: "${LEC_ROOT:?Set LEC_ROOT to the parent of MHV and NMHV last chains}"
TEST_ROOT=${TEST_ROOT:-output/prereduce_test_$(date +%Y%m%d_%H%M%S)_$$}
[[ ! -e $TEST_ROOT && ! -L $TEST_ROOT ]] || { echo 'Fresh TEST_ROOT required' >&2; exit 2; }
mkdir -p "$TEST_ROOT"
export THREADS=${THREADS:-2}
for helicity in MHV NMHV; do
    for k in 1 2; do
        h=$((4-k)); reference=$TEST_ROOT/${helicity}_${h}p${k}_none
        for mode in none left right default; do
            target=$TEST_ROOT/${helicity}_${h}p${k}_${mode}
            extra=(); [[ $mode == default ]] || extra+=(--prereduce "$mode")
            ./bin/irr_sew "$helicity" "$h" "$k" "$target" --fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT/$helicity" --verify "${extra[@]}" > "$target.log" 2>&1
            for file in "$reference"/irrSEW_*.wxf "$reference"/irrSEW_*.wxf.shape; do
                [[ -f $file ]] || continue
                cmp "$file" "$target/$(basename "$file")"
            done
        done
    done
done
echo "PASS: $TEST_ROOT"
