"""Serve the real app against disposable project storage for browser tests."""
from pathlib import Path
import os
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'front-end/server'))
from app import config, storage
config.PROJECTS_DIR = storage.PROJECTS_DIR = Path(os.environ.get('TEST_PROJECTS_DIR') or tempfile.mkdtemp(prefix='symbology-browser-audit.'))
from app.main import app
import uvicorn
print('Isolated projects:', storage.PROJECTS_DIR, flush=True)
uvicorn.run(app, host='127.0.0.1', port=int(os.environ.get('TEST_PORT', '18321')), log_level='warning')
