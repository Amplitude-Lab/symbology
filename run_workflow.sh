#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Multi-project support: set PROJECT=<name> (e.g. PROJECT=E7) to use
# data_<name>/ and output_<name>/ instead of the default data/ and output/.
if [[ -n "${PROJECT:-}" ]]; then
  DATA_DIR="data_${PROJECT}"
  OUTPUT_DIR="output_${PROJECT}"
else
  DATA_DIR="data"
  OUTPUT_DIR="output"
fi

mode="${1:-smoke}"
case "$mode" in
  smoke|full) ;;
  *)
    echo "Usage: $0 [smoke|full]" >&2
    echo "  Set PROJECT=<name> to use data_<name>/ and output_<name>/." >&2
    exit 1
    ;;
esac

# Build. On macOS the Makefile auto-selects a Homebrew GCC (libc++ cannot
# handle std::execution::par and std::chrono::zoned_time used by the project)
# and adds the Homebrew include/lib/rpath flags itself, so a plain `make`
# works from any machine; override with `make CXX=...` if needed.
make

BOOTSTRAP="./bootstrap"
if [[ -x "./bootstrap.exe" ]]; then
  BOOTSTRAP="./bootstrap.exe"
fi

# Auto-detect dlogmat seed (dlogmat_<group>.wxf) in DATA_DIR.
# This avoids hardcoding the group name (e.g. E6) so the script generalizes
# to other symmetry groups via PROJECT=<name>.
shopt -s nullglob
dlogmat_files=("${DATA_DIR}"/dlogmat_*.wxf)
shopt -u nullglob
if [[ ${#dlogmat_files[@]} -eq 0 ]]; then
  echo "Error: no dlogmat_*.wxf found in ${DATA_DIR}/" >&2
  exit 1
fi
CONDITION="${dlogmat_files[0]}"
echo "Using condition: ${CONDITION}"

mkdir -p "${OUTPUT_DIR}" logs

run_step() {
  local name="$1"
  shift
  echo "==> $name"
  "$@" 2>&1 | tee "logs/${name}.log"
}

extend_fec_to_6() {
  run_step extend_FEC_2 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${DATA_DIR}/FEC_1.wxf"    -o "${OUTPUT_DIR}/FEC_2.wxf"
  run_step extend_FEC_3 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_2.wxf" -o "${OUTPUT_DIR}/FEC_3.wxf"
  run_step extend_FEC_4 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_3.wxf" -o "${OUTPUT_DIR}/FEC_4.wxf"
  run_step extend_FEC_5 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_4.wxf" -o "${OUTPUT_DIR}/FEC_5.wxf"
  run_step extend_FEC_6 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_5.wxf" -o "${OUTPUT_DIR}/FEC_6.wxf"
}

extend_lec_to_4() {
  run_step extend_LEC_2 "$BOOTSTRAP" --extend -c "$CONDITION" -l "${DATA_DIR}/LEC_1.wxf"    -o "${OUTPUT_DIR}/LEC_2.wxf"
  run_step extend_LEC_3 "$BOOTSTRAP" --extend -c "$CONDITION" -l "${OUTPUT_DIR}/LEC_2.wxf" -o "${OUTPUT_DIR}/LEC_3.wxf"
  run_step extend_LEC_4 "$BOOTSTRAP" --extend -c "$CONDITION" -l "${OUTPUT_DIR}/LEC_3.wxf" -o "${OUTPUT_DIR}/LEC_4.wxf"
}

if [[ "$mode" == "smoke" ]]; then
  extend_fec_to_6
  extend_lec_to_4
  run_step sew_2p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_2.wxf" -l "${OUTPUT_DIR}/LEC_2.wxf" -o "${OUTPUT_DIR}/SEW_2p2.wxf"
  run_step sew_3p1 "$BOOTSTRAP" --sew -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_3.wxf" -l "${DATA_DIR}/LEC_1.wxf"   -o "${OUTPUT_DIR}/SEW_3p1.wxf"
  run_step sew_4p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_4.wxf" -l "${OUTPUT_DIR}/LEC_2.wxf" -o "${OUTPUT_DIR}/SEW_4p2.wxf"
  run_step sew_5p1 "$BOOTSTRAP" --sew -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_5.wxf" -l "${DATA_DIR}/LEC_1.wxf"   -o "${OUTPUT_DIR}/SEW_5p1.wxf"
else
  extend_fec_to_6
  run_step extend_FEC_7 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_6.wxf" -o "${OUTPUT_DIR}/FEC_7.wxf"
  run_step extend_FEC_8 "$BOOTSTRAP" --extend -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_7.wxf" -o "${OUTPUT_DIR}/FEC_8.wxf"
  extend_lec_to_4
  run_step extend_LEC_5 "$BOOTSTRAP" --extend -c "$CONDITION" -l "${OUTPUT_DIR}/LEC_4.wxf" -o "${OUTPUT_DIR}/LEC_5.wxf"
  run_step sew_8p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f "${OUTPUT_DIR}/FEC_8.wxf" -l "${OUTPUT_DIR}/LEC_2.wxf" -o "${OUTPUT_DIR}/SEW_8p2.wxf"
fi

echo ""
echo "=== Workflow complete (${mode}) ==="
echo "Data dir:   ${DATA_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Logs:       logs/*.log"
