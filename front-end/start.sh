#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
FRONTEND_DIR="$(pwd)"
REPO_ROOT="$(dirname "$FRONTEND_DIR")"
PORT=8321

echo "== Symbology Studio =="

if ! command -v python3 >/dev/null 2>&1; then
  echo "Error: python3 was not found. Please install Python 3.9+ first." >&2
  exit 1
fi

if [ ! -x "$REPO_ROOT/bootstrap" ]; then
  echo "The bootstrap binary is missing; building it (requires a C++23 compiler)..."
  (cd "$REPO_ROOT" && make bootstrap)
fi

if [ ! -d server/.venv ]; then
  echo "Setting up the Python environment (one-time)..."
  python3 -m venv server/.venv
  server/.venv/bin/pip install --quiet --upgrade pip
  server/.venv/bin/pip install --quiet -r server/requirements.txt
fi

if [ ! -d web/dist ]; then
  if command -v npm >/dev/null 2>&1; then
    echo "Building the web interface (one-time)..."
    (cd web && (npm install --no-audit --no-fund || npm install --no-audit --no-fund --strict-ssl=false || env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY npm install --no-audit --no-fund --strict-ssl=false) && npm run build)
  else
    echo "Note: npm not found; skipping the prebuilt UI."
    echo "You can run the UI in dev mode later with: cd web && npm install && npm run dev"
  fi
fi

URL="http://127.0.0.1:$PORT"
echo ""
echo "Starting the server at $URL"
echo "Press Ctrl+C to stop."
( sleep 2; (command -v open >/dev/null 2>&1 && open "$URL") || (command -v xdg-open >/dev/null 2>&1 && xdg-open "$URL") || true ) &

cd server
exec ./.venv/bin/python run.py
