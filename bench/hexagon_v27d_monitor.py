#!/usr/bin/env python3
"""Read-only tmux dashboard; the solver belongs to a separate systemd cgroup."""
import json
from pathlib import Path
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]
CONTROL = REPO / 'output_hexagon_ten_v27d_20261009/runner'
UNIT = 'symbology-hexagon-ten-v27d-20261009.service'


def display():
    while True:
        print('\033[H\033[2J', end='')
        print('Ten-loop MHV hexagon — V27D / 16 workers / 50 GB maximum')
        print('Detach: Ctrl-b then d. The solver runs independently of this monitor.\n')
        subprocess.run(['systemctl', '--user', 'show', UNIT, '-p', 'ActiveState', '-p', 'SubState',
                        '-p', 'MemoryCurrent', '-p', 'MemoryPeak', '-p', 'CPUUsageNSec'])
        try:
            progress = json.loads((CONTROL / 'progress.json').read_text())
            print('\nState:', progress['state'], '| completed primes:', progress['native_primes_committed'])
            print('Previously accounted CPU-hours:', round(progress['accounted_cpu_hours'], 5))
            status_path = progress.get('resource_status')
            if status_path and Path(status_path).exists():
                status = json.loads(Path(status_path).read_text())
                print('Current phase: wall seconds', round(status.get('wall_seconds', 0)),
                      '| RSS GB', round(status.get('current_rss_bytes', 0) / 1e9, 3))
            log = progress.get('solver_log')
            if log and Path(log).exists():
                print('\nLatest solver output:', log, '\n', flush=True)
                subprocess.run(['tail', '-n', '18', log])
            if progress.get('error'):
                print('\nATTENTION:', progress['error'])
        except (OSError, ValueError, KeyError) as error:
            print('Waiting for progress:', error)
        sys.stdout.flush()
        time.sleep(10)


if __name__ == '__main__':
    if '--start' in sys.argv:
        if subprocess.run(['tmux', '-L', 'hexagon-v27d', 'has-session', '-t', 'hexagon-v27d'], capture_output=True).returncode:
            subprocess.run(['tmux', '-L', 'hexagon-v27d', 'new-session', '-d', '-s', 'hexagon-v27d',
                            '/usr/bin/python3', str(Path(__file__).resolve())], check=True)
    else:
        display()
