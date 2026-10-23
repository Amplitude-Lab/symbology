#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
STAGE=${STAGE:-${1:-SEW}}; STAGE=${STAGE^^}
HELICITY=${HELICITY:-NMHV}; HELICITY=${HELICITY^^}
LOOP=${LOOP:-3}; TAIL_WEIGHT=${TAIL_WEIGHT:-2}
OUTPUT_ROOT=${OUTPUT_ROOT:-$ROOT/output}
LOG_ROOT=${LOG_ROOT:-$ROOT/log}
BUILD=${BUILD:-1}; BUILD_JOBS=${BUILD_JOBS:-2}
export THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN)}}
positive() { [[ $2 =~ ^[1-9][0-9]*$ ]] || { echo "$1 must be a positive integer" >&2; exit 2; }; }
case $STAGE in FEC|LEC|SEW|CONDITIONS|SHIFT|RELEASE) ;; *) echo 'STAGE: FEC, LEC, SEW, CONDITIONS, SHIFT, RELEASE' >&2; exit 2;; esac
[[ $HELICITY == MHV || $HELICITY == NMHV ]] || { echo 'HELICITY: MHV or NMHV' >&2; exit 2; }
positive LOOP "$LOOP"; positive TAIL_WEIGHT "$TAIL_WEIGHT"; positive THREADS "$THREADS"; positive BUILD_JOBS "$BUILD_JOBS"
HEAD_WEIGHT=${HEAD_WEIGHT:-$((2*LOOP-TAIL_WEIGHT))}; positive HEAD_WEIGHT "$HEAD_WEIGHT"
MAX_WEIGHT=${MAX_WEIGHT:-$HEAD_WEIGHT}; positive MAX_WEIGHT "$MAX_WEIGHT"
TARGET_HEAD=${TARGET_HEAD:-$HEAD_WEIGHT}; positive TARGET_HEAD "$TARGET_HEAD"
[[ $STAGE == FEC || $STAGE == LEC || $((HEAD_WEIGHT+TAIL_WEIGHT)) == $((2*LOOP)) ]] || { echo 'HEAD_WEIGHT+TAIL_WEIGHT must equal 2*LOOP' >&2; exit 2; }
TARGET_TAIL=$((HEAD_WEIGHT+TAIL_WEIGHT-TARGET_HEAD)); positive TARGET_TAIL "$TARGET_TAIL"
mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT"
OUTPUT_ROOT=$(cd "$OUTPUT_ROOT" && pwd)
exec > >(tee -a "$LOG_ROOT/${STAGE}_${HELICITY}_$(date +%Y%m%d_%H%M%S)_$$.log") 2>&1
echo "stage=$STAGE helicity=$HELICITY head=$HEAD_WEIGHT tail=$TAIL_WEIGHT threads=$THREADS"
if [[ $BUILD == 1 ]]; then "${MAKE:-make}" -j"$BUILD_JOBS"; elif [[ $BUILD != 0 ]]; then echo 'BUILD: 0 or 1'; exit 2; fi
FEC_ROOT=${FEC_ROOT:-$OUTPUT_ROOT/FEC}
LEC_ROOT=${LEC_ROOT:-$OUTPUT_ROOT/LEC/$HELICITY}
SEW_PARENT=$OUTPUT_ROOT/SEW_${HEAD_WEIGHT}p${TAIL_WEIGHT}/$HELICITY
SEW_ROOT=${SEW_ROOT:-$SEW_PARENT/bulk}
FEC_DLOGMAT=${FEC_DLOGMAT:-$ROOT/data/dlogmat_ca_1term_pd3.wxf}
FEC_BLOCKS=${FEC_BLOCKS:-$ROOT/data/dlogmat_ca_1term_pd3_blocks.tsv}
LEC_DLOGMAT=${LEC_DLOGMAT:-$ROOT/data/dlogmat_full_ca_pd3.wxf}
LEC_BLOCKS=${LEC_BLOCKS:-$ROOT/data/dlogmat_full_ca_pd3_blocks.tsv}
SEW_DLOGMAT=${SEW_DLOGMAT:-$LEC_DLOGMAT}; SEW_BLOCKS=${SEW_BLOCKS:-$LEC_BLOCKS}
# Locks live beside shared chains, so independent SEW jobs can share FEC/LEC.
# flock releases on exit/signals; no costly tensor/cache is removed.
locked() { local path=$1; shift; mkdir -p "$(dirname "$path")"; (flock -x 9; "$@") 9>"$path"; }
# Keep input metadata stable while another job resumes/extends a shared chain.
using_chains() { (flock -s 7; flock -s 8; "$@") 7>"$FEC_ROOT.lock" 8>"$LEC_ROOT.lock"; }
fec() {
 local args=()
 case ${FEC_CLEAN_BLOCK_CACHE:-0} in 0) ;; 1) args+=(--clean-block-cache) ;; *) echo 'FEC_CLEAN_BLOCK_CACHE: 0 or 1' >&2; return 2;; esac
 locked "$FEC_ROOT.lock" "$ROOT/bin/irr_fec" "$1" "$FEC_ROOT" --condition "$FEC_DLOGMAT" --condition-blocks "$FEC_BLOCKS" "${args[@]}"
}
lec() { locked "$LEC_ROOT.lock" "$ROOT/bin/irr_lec" "$HELICITY" "$1" "$LEC_ROOT" --data "$ROOT/data" --condition "$LEC_DLOGMAT" --condition-blocks "$LEC_BLOCKS"; }
conditions() {
 local args=(--fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT" --data "${PHYSICS_DATA:-$ROOT/data}" --primes "${MODULAR_PRIMES:-auto}")
 [[ -z ${CONDITION_PLAN:-} ]] || args+=(--plan "$CONDITION_PLAN")
 [[ -z ${CONDITIONS:-} ]] || args+=(--conditions "$CONDITIONS")
 [[ -z ${SAMPLES:-} ]] || args+=(--samples "$SAMPLES")
 [[ -z ${SAMPLE_BATCH:-} ]] || args+=(--batch "$SAMPLE_BATCH")
 [[ -z ${VALIDATION_SAMPLES:-} ]] || args+=(--validation-samples "$VALIDATION_SAMPLES")
 [[ ${EXPORT_INTERMEDIATE:-0} == 0 ]] || args+=(--export-all)
 using_chains "$ROOT/bin/irr_conditions" "$HELICITY" "$HEAD_WEIGHT" "$TAIL_WEIGHT" "${INPUT_SEW:-$SEW_ROOT}" "$SEW_PARENT" "${args[@]}"
 case ${KEEP_MODULAR_CACHE:-1} in
  1) ;; 0) "$ROOT/bin/irr_cache" "$SEW_PARENT" --apply ;;
  *) echo 'KEEP_MODULAR_CACHE: 0 or 1' >&2; return 2 ;;
 esac
}
case $STAGE in
 FEC) fec "$MAX_WEIGHT" ;;
 LEC) lec "${LEC_WEIGHT:-$TAIL_WEIGHT}" ;;
 SEW)
  if [[ ${FINITE_FIELD_PIPELINE:-0} == 1 ]]; then
   max_tail=$TAIL_WEIGHT
   if [[ -n ${CONDITION_PLAN:-} ]]; then
    while read -r label op value rest; do
     if [[ $op == shift ]]; then value=${value%$'\r'}; positive target_head "$value"; tail=$((HEAD_WEIGHT+TAIL_WEIGHT-value)); ((tail<=max_tail)) || max_tail=$tail; fi
    done < "$CONDITION_PLAN"
   fi
   fec "$HEAD_WEIGHT"; lec "$max_tail"
   args=(--fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT" --data "${PHYSICS_DATA:-$ROOT/data}" --condition "$SEW_DLOGMAT" --condition-blocks "$SEW_BLOCKS" --primes "${MODULAR_PRIMES:-auto}")
   [[ -z ${CONDITION_PLAN:-} ]] || args+=(--plan "$CONDITION_PLAN")
   [[ -z ${CONDITIONS:-} ]] || args+=(--conditions "$CONDITIONS")
   [[ -z ${INPUT_SEW:-} ]] || args+=(--input "$INPUT_SEW")
   [[ -z ${SEW_PREREDUCE:-} ]] || args+=(--prereduce "$SEW_PREREDUCE")
   [[ -z ${PROJECTION:-} ]] || args+=(--projection "$PROJECTION")
   args+=(--qq-cache-mib "${QQ_CACHE_MIB:-0}")
   [[ -z ${SAMPLES:-} ]] || args+=(--samples "$SAMPLES")
   [[ -z ${SAMPLE_BATCH:-} ]] || args+=(--batch "$SAMPLE_BATCH")
   [[ -z ${VALIDATION_SAMPLES:-} ]] || args+=(--validation-samples "$VALIDATION_SAMPLES")
   [[ ${CACHE_ONLY:-0} == 0 ]] || args+=(--cache-only)
   [[ ${EXPORT_INTERMEDIATE:-0} == 0 ]] || args+=(--export-all)
   "$ROOT/bin/irr_pipeline" "$HELICITY" "$HEAD_WEIGHT" "$TAIL_WEIGHT" "$OUTPUT_ROOT" "${args[@]}"
   case ${KEEP_MODULAR_CACHE:-1} in 1) ;; 0) "$ROOT/bin/irr_cache" "$OUTPUT_ROOT" --apply ;; *) echo 'KEEP_MODULAR_CACHE: 0 or 1'; exit 2;; esac
  else
  fec "$HEAD_WEIGHT"; lec "$TAIL_WEIGHT"
  extra=(); [[ ${VERIFY:-0} == 0 ]] || extra+=(--verify)
  [[ -z ${SEW_PREREDUCE:-} ]] || extra+=(--prereduce "$SEW_PREREDUCE")
  [[ -z ${EXPECT_DIM:-} ]] || extra+=(--expect "$EXPECT_DIM")
  using_chains locked "$SEW_ROOT.lock" "$ROOT/bin/irr_sew" "$HELICITY" "$HEAD_WEIGHT" "$TAIL_WEIGHT" "$SEW_ROOT" --fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT" --condition "$SEW_DLOGMAT" --condition-blocks "$SEW_BLOCKS" "${extra[@]}"
  [[ -z ${CONDITIONS:-}${CONDITION_PLAN:-} ]] || conditions
  fi
  ;;
 CONDITIONS) conditions ;;
 SHIFT|RELEASE)
  INPUT_SEW=${INPUT_SEW:-$SEW_ROOT}
  fec "$((HEAD_WEIGHT>TARGET_HEAD?HEAD_WEIGHT:TARGET_HEAD))"; lec "$((TAIL_WEIGHT>TARGET_TAIL?TAIL_WEIGHT:TARGET_TAIL))"
  args=(--head "$TARGET_HEAD" --fec-root "$FEC_ROOT" --lec-root "$LEC_ROOT")
  [[ -z ${DIRECTIONS:-} ]] || args+=(--directions "$DIRECTIONS")
  [[ -z ${ORIGIN:-} ]] || args+=(--origin "$ORIGIN")
  [[ ${UNIQUE:-0} == 0 ]] || args+=(--unique)
  if [[ $STAGE == SHIFT ]]; then
   history=bulk
   for meta in "$INPUT_SEW/condition_history.tsv" "$INPUT_SEW/stage.tsv"; do
    [[ ! -f $meta ]] || history=$(awk -F '\t' '$1=="history" {gsub(/\r/,"",$2); print $2}' "$meta")
   done
   DESTINATION=${DESTINATION:-$OUTPUT_ROOT/SEW_${TARGET_HEAD}p${TARGET_TAIL}/$HELICITY/$history}
   mode=shift
  else
   DESTINATION=${DESTINATION:-$ROOT/release_${HELICITY}_w$((2*LOOP))_${TARGET_HEAD}p${TARGET_TAIL}}
   mode=export
   [[ -z ${FIRST_NAME:-} ]] || args+=(--first-name "$FIRST_NAME")
   [[ -z ${LAST_NAME:-} ]] || args+=(--last-name "$LAST_NAME")
   [[ -z ${SEW_NAME:-} ]] || args+=(--sew-name "$SEW_NAME")
   [[ -z ${ORIGIN_NAME:-} ]] || args+=(--origin-name "$ORIGIN_NAME")
  fi
  using_chains locked "$DESTINATION.lock" "$ROOT/bin/irr_convert" "$mode" "$HELICITY" "$HEAD_WEIGHT" "$TAIL_WEIGHT" "$INPUT_SEW" "$DESTINATION" "${args[@]}"
  ;;
esac
echo 'WORKFLOW COMPLETE'
