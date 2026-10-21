#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
root=${1:?Pass a completed smoke output root}
export OUTPUT_ROOT=$root THREADS=${THREADS:-2} BUILD=0 LOOP=3 TAIL_WEIGHT=2 HELICITY=NMHV
STAGE=CONDITIONS CONDITIONS=collinear+mrk SAMPLES=64 EXPORT_INTERMEDIATE=1 bash run_workflow.sh
# Repeat a physical condition on the unique affine result: it must still be checked.
STAGE=CONDITIONS INPUT_SEW="$root/SEW_4p2/NMHV/collinear_mrk" CONDITIONS=mrk SAMPLES=64 bash run_workflow.sh
for b in pe_A1 pe_A2 pe_E po_A1 po_A2 po_E; do
 bin/irr_tensor_compare "$root/SEW_4p2/NMHV/collinear_mrk/irrSEW_4p2_${b}_origin.wxf" "$root/SEW_4p2/NMHV/collinear_mrk_mrk/irrSEW_4p2_${b}_origin.wxf"
done
STAGE=RELEASE INPUT_SEW="$root/SEW_4p2/NMHV/collinear_mrk" DESTINATION="$root/release_physics" bash run_workflow.sh
echo 'PHYSICAL STAGE RESUME, UNIQUE CHECK AND RELEASE PASS'
