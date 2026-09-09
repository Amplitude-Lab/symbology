#!/bin/zsh
# SSR render audit: bundles scripts/ssr-audit.mjs with the UI stubs applied
# (see ssr-audit-build.mjs — the stub routing is an esbuild plugin because
# relative alias names are not expressible on the CLI) and runs it.
set -e
cd "$(dirname "$0")/.."
node scripts/ssr-audit-build.mjs
node /tmp/ssr-audit.cjs
