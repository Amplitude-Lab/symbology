#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bin/irr_conditions --self-test
bin/irr_cache --self-test
