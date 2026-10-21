#!/usr/bin/env python3
"""One-shot cleanup guard for this task's bounded validation operation."""
import json
from pathlib import Path
import subprocess
import sys
import time

validation, native, receipt = sys.argv[1:]
seen_running = False
start = time.time()
try:
    while time.time() - start < 3700:
        status = subprocess.run(['systemctl', '--user', 'show', validation,
            '--property=ActiveState', '--value'], text=True, capture_output=True).stdout.strip()
        if status in ('active', 'activating'):
            seen_running = True
        elif seen_running or time.time() - start > 60:
            break
        time.sleep(5)
    else:
        subprocess.run(['systemctl', '--user', 'stop', validation], timeout=90)
finally:
    result = subprocess.run(['systemctl', '--user', 'thaw', native], text=True, capture_output=True)
    Path(receipt).write_text(json.dumps({'started_unix_seconds': start,
        'finished_unix_seconds': time.time(), 'validation_seen_running': seen_running,
        'thaw_exit_code': result.returncode, 'stderr': result.stderr}, indent=2) + '\n')
