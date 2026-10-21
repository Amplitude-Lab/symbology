#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OUTPUT_ROOT=${TEST_OUTPUT_ROOT:-$PWD/output/smoke_$(date +%Y%m%d_%H%M%S)_$$}
export THREADS=${THREADS:-2} BUILD=0 LOOP=3 TAIL_WEIGHT=2 VERIFY=1 EXPECT_DIM=4
STAGE=SEW bash run_workflow.sh
# Idempotent repeat must reuse blocks and the completed sewing.
STAGE=SEW bash run_workflow.sh
STAGE=SHIFT TARGET_HEAD=3 bash run_workflow.sh
STAGE=SHIFT HEAD_WEIGHT=3 TAIL_WEIGHT=3 TARGET_HEAD=4 INPUT_SEW="$OUTPUT_ROOT/SEW_3p3/NMHV/bulk" DESTINATION="$OUTPUT_ROOT/SEW_4p2/NMHV/roundtrip" bash run_workflow.sh
STAGE=RELEASE TARGET_HEAD=3 DESTINATION="$OUTPUT_ROOT/release_NMHV_w6" bash run_workflow.sh
HELICITY=MHV LOOP=2 TAIL_WEIGHT=2 EXPECT_DIM=3 STAGE=SEW bash run_workflow.sh
HELICITY=MHV LOOP=2 TAIL_WEIGHT=2 STAGE=RELEASE DESTINATION="$OUTPUT_ROOT/release_MHV_w4" bash run_workflow.sh
HELICITY=MHV LEC_WEIGHT=3 STAGE=LEC bash run_workflow.sh
for block in pe_A1 pe_A2 pe_E po_A1 po_A2 po_E; do
 for w in 1 2 3 4; do
  printf -v dir 'w%02d' "$w"
  bin/irr_tensor_compare "tests/reference/FEC/$dir/irrFEC_${w}_${block}.wxf" "$OUTPUT_ROOT/FEC/$dir/irrFEC_${w}_${block}.wxf"
 done
 for helicity in MHV NMHV; do
  for w in 1 2 3; do
   printf -v dir 'w%02d' "$w"
   bin/irr_tensor_compare "tests/reference/LEC/$helicity/$dir/irrLEC_${w}_${block}.wxf" "$OUTPUT_ROOT/LEC/$helicity/$dir/irrLEC_${w}_${block}.wxf"
  done
 done
 bin/irr_tensor_compare "$OUTPUT_ROOT/SEW_4p2/NMHV/bulk/irrSEW_4p2_${block}.wxf" "$OUTPUT_ROOT/SEW_4p2/NMHV/roundtrip/irrSEW_4p2_${block}.wxf"
done
while IFS= read -r sums; do (cd "$(dirname "$sums")"; tr -d '\r' < checksums.sha256 | sha256sum -c - >/dev/null); done < <(find "$OUTPUT_ROOT" -name checksums.sha256)
echo "SMOKE PASS: $OUTPUT_ROOT"
