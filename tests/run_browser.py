"""Run browser regressions with an isolated server on Linux or Windows."""
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def main():
    node = shutil.which('node')
    if not node:
        raise RuntimeError('Node.js is required for browser tests')
    output = Path(os.environ.get('TEST_OUTPUT_DIR', Path(tempfile.gettempdir()) / 'symbology-browser-results'))
    output.mkdir(parents=True, exist_ok=True)
    with socket.socket() as listener:
        listener.bind(('127.0.0.1', 0))
        port = listener.getsockname()[1]
    base = f'http://127.0.0.1:{port}'
    with tempfile.TemporaryDirectory(prefix='symbology-browser-') as projects, (output / 'server.log').open('wb') as log:
        env = {**os.environ, 'TEST_PORT': str(port), 'TEST_BASE_URL': base,
               'TEST_PROJECTS_DIR': projects, 'TEST_OUTPUT_DIR': str(output.resolve()),
               'PYTHONUTF8': '1'}
        server = subprocess.Popen([sys.executable, str(ROOT / 'tests/serve_browser.py')],
                                  cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            deadline = time.monotonic() + 30
            while True:
                if server.poll() is not None:
                    raise RuntimeError(f'Test server exited; see {output / "server.log"}')
                try:
                    with urllib.request.urlopen(base + '/api/projects', timeout=1) as response:
                        if response.status == 200:
                            break
                except (urllib.error.URLError, TimeoutError):
                    pass
                if time.monotonic() >= deadline:
                    raise RuntimeError(f'Test server did not start; see {output / "server.log"}')
                time.sleep(0.1)
            return subprocess.run([node, str(ROOT / 'tests/browser_regression.cjs')],
                                  cwd=ROOT, env=env, timeout=180).returncode
        finally:
            if server.poll() is None:
                server.terminate()
                try:
                    server.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    server.kill()
                    server.wait(timeout=10)


if __name__ == '__main__':
    raise SystemExit(main())
