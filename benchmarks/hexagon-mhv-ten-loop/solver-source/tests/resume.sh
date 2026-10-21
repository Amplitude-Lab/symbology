#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
root=${1:?Pass a completed smoke output root}
destination="$root/resume_fec"
[[ ! -e $destination ]] || { echo 'Use a new resume test destination'; exit 2; }
mkdir -p "$destination/w02"
cp "$root/FEC/input_identity.txt" "$destination/"
cp -R "$root/FEC/w01" "$destination/"
# Simulate a job interrupted after saving primes, without deleting any result.
cp -R "$root/FEC/w02/cache" "$destination/w02/"
THREADS=2 bin/irr_fec 2 "$destination" > "$root/resume_fec.log" 2>&1
grep -q CRT_CACHE_REUSED "$root/resume_fec.log"
for block in pe_A1 pe_A2 pe_E po_A1 po_A2 po_E; do
 bin/irr_tensor_compare "$root/FEC/w02/irrFEC_2_${block}.wxf" "$destination/w02/irrFEC_2_${block}.wxf"
done
if THREADS=2 bin/irr_fec 2 "$destination" --condition data/dlogmat_full_ca_pd3.wxf --condition-blocks data/dlogmat_full_ca_pd3_blocks.tsv > "$root/reject_condition.log" 2>&1; then
 echo 'Changed condition was incorrectly accepted'; exit 1
fi
grep -q 'Input identity mismatch' "$root/reject_condition.log"
echo 'CACHE_RESUME_AND_CONDITION_REJECTION_PASS'
