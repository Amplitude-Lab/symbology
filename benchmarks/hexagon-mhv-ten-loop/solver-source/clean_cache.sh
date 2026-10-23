#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Default is a preview. Set OUTPUT_ROOT and pass --apply for one-click cleanup.
exec "$ROOT/bin/irr_cache" "${OUTPUT_ROOT:-$ROOT/output}" "$@"
