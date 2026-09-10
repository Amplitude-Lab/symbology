#!/usr/bin/env bash
# Fetch SparseRREF at the pinned commit and apply the required local patches.
#
# The C++ core needs SparseRREF (sparse linear algebra + WXF I/O). Upstream
# has a data race in the parallel RREF (plain int flags publishing rewritten
# rows -> ~10% SIGSEGV at weight >= 3; see patches/BUGREPORT-sparse-mat-rref-race.md),
# and the symbology patch set additionally makes tensor_contract /
# file_to_ustr fail loudly and removes the dead >1GB mmap path in the WXF
# readers. This script reproduces the exact expected state of SparseRREF/
# on a fresh machine:
#
#   ./scripts/setup-sparserref.sh          # clone at PINNED_COMMIT + patch
#   ./scripts/setup-sparserref.sh --check  # only verify an existing checkout
#
# Idempotent: re-running on an already-patched tree verifies and exits 0.
set -euo pipefail

cd "$(dirname "$0")/.."

PINNED_COMMIT="5bbee55"
REPO_URL="https://github.com/munuxi/SparseRREF"
PATCH="patches/sparserref-5bbee55-race-and-init-fix.patch"
NEWER_PATCH="patches/sparserref-race-fix.patch"
WXF_PATCH="patches/sparserref-wxf-validation.patch"

apply_wxf_patch() {
  if git -C SparseRREF apply --reverse --check "../$WXF_PATCH" 2>/dev/null; then
    return
  fi
  git -C SparseRREF apply --check "../$WXF_PATCH" || die "WXF patch does not match this checkout; preserve local edits and use the pinned SparseRREF version."
  git -C SparseRREF apply "../$WXF_PATCH"
}

die() { echo "setup-sparserref: $*" >&2; exit 1; }

check_tree() {
  local n
  n=$(grep -c "std::atomic<int>" SparseRREF/sparse_mat.h || true)
  [[ "$n" == "2" ]] || die "SparseRREF/sparse_mat.h has $n (want 2) std::atomic<int> occurrences — the race fix is missing; intermittent SIGSEGV will follow. Re-run without --check to re-apply patches."
  grep -q "must not return an empty tensor" SparseRREF/sparse_tensor.h \
    || die "SparseRREF/sparse_tensor.h lacks the loud tensor_contract fix — patches out of date."
  grep -q 'Checked WXF decoding' SparseRREF/wxf_parser.h || die "WXF bounds patch is missing"
  grep -q 'wxf_checked::validate' SparseRREF/wxf_support.h || die "SparseArray validation patch is missing"
  echo "setup-sparserref: SparseRREF/ patches verified (race fix + loud-failure fixes present)."
}

[[ -f "$PATCH" ]] || die "patch file $PATCH not found — run from the repository root"

if [[ "${1:-}" == "--check" ]]; then
  [[ -f SparseRREF/sparse_mat.h ]] || die "SparseRREF/ not present"
  check_tree
  exit 0
fi

if [[ ! -d SparseRREF ]]; then
  echo "setup-sparserref: cloning SparseRREF at $PINNED_COMMIT"
  git clone "$REPO_URL" SparseRREF
  ( cd SparseRREF && git checkout -q "$PINNED_COMMIT" )
elif [[ -z "$(cd SparseRREF && git status --porcelain 2>/dev/null)" ]]; then
  # clean tree: make sure it is at the pinned commit
  ( cd SparseRREF \
    && [[ "$(git rev-parse --short HEAD)" == "$PINNED_COMMIT" ]] \
    || { echo "setup-sparserref: resetting SparseRREF to $PINNED_COMMIT"; git checkout -q "$PINNED_COMMIT"; } )
else
  echo "setup-sparserref: SparseRREF/ already patched (dirty tree) — verifying"
  apply_wxf_patch
  check_tree
  exit 0
fi

echo "setup-sparserref: applying $PATCH"
git -C SparseRREF apply "../$PATCH" || {
  echo "  (direct apply failed — trying newer-upstream race-only patch)" >&2
  git -C SparseRREF apply "../$NEWER_PATCH" || die "could not apply either patch; SparseRREF upstream may have moved — regenerate patches/ and update PINNED_COMMIT"
}

apply_wxf_patch
check_tree
echo "setup-sparserref: done. Build with: make"
