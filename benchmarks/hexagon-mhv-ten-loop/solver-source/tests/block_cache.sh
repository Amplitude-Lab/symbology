#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${TEST_BLOCK_ROOT:?Set TEST_BLOCK_ROOT to a fresh scratch directory}"
make -j2 bin/irr_fec${EXE:-} bin/irr_cache${EXE:-}
python3 tests/block_cache.py "$TEST_BLOCK_ROOT" --exe "bin/irr_fec${EXE:-}"
