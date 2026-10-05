"""Run sequential, measured stages within one shared time and memory budget.

The JSON configuration contains a fresh run directory, an absolute deadline,
and stages with argv lists, required outputs and required success log markers.
This controller does not infer success from the existence of a candidate file.
"""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024**2), b''):
            h.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=h.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root = Path(config['run_directory']).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / 'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root / 'status.json').exists():
            raise RuntimeError('refusing to overwrite an existing pipeline run')
        started = time.time()
        monotonic_start = time.monotonic()
        allowed = max(0.0, config['deadline_unix_seconds'] - started)
        state = dict(state='running', controller_pid=os.getpid(),
                     started_unix_seconds=started,
                     deadline_unix_seconds=config['deadline_unix_seconds'],
                     rss_limit_gib=config['rss_gib'], stages=[])
        monitor = None
        cancelled = False

        def request_stop(signum, frame):
            nonlocal cancelled
            cancelled = True
            if monitor is not None and monitor.poll() is None:
                monitor.send_signal(signal.SIGTERM)

        signal.signal(signal.SIGTERM, request_stop)
        signal.signal(signal.SIGINT, request_stop)

        def publish():
            state['updated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            state['wall_seconds'] = time.monotonic() - monotonic_start
            atomic_json(root / 'status.json', state)

        publish()
        try:
            for stage in config['stages']:
                if cancelled:
                    state.update(state='cancelled', current_stage=None)
                    break
                remaining = min(config['deadline_unix_seconds'] - time.time(),
                                allowed - (time.monotonic() - monotonic_start))
                if remaining <= 0:
                    state.update(state='time_limit', current_stage=None)
                    break
                name = stage['name']
                if not name or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in name):
                    raise RuntimeError('invalid stage name')
                report, log = root / (name + '.json'), root / (name + '.log')
                live = root / (name + '-live.json')
                if report.exists() or log.exists():
                    raise RuntimeError('stage output already exists')
                state.update(current_stage=name, current_log=str(log), current_live=str(live))
                publish()
                command = [sys.executable, config['measure_program'],
                           '--report', str(report), '--log', str(log),
                           '--status', str(live), '--timeout', str(remaining),
                           '--virtual-gib', str(config['virtual_gib']),
                           '--rss-gib', str(config['rss_gib']),
                           '--min-available-gib', str(config.get('min_available_gib', 0)),
                           '--', *stage['command']]
                monitor = subprocess.Popen(command, cwd=config['working_directory'])
                if cancelled:
                    monitor.send_signal(signal.SIGTERM)
                code = monitor.wait()
                monitor = None
                measured = json.loads(report.read_text()) if report.exists() else {}
                entry = dict(name=name, measurement=measured, monitor_exit_code=code)
                state['stages'].append(entry)
                if code or cancelled:
                    state.update(state='cancelled' if cancelled else 'time_limit' if measured.get('timeout') else 'failed')
                    publish()
                    break
                missing = [p for p in stage.get('required_outputs', []) if not Path(p).is_file()]
                log_text = log.read_text(errors='replace')
                missing_markers = [m for m in stage.get('success_markers', []) if m not in log_text]
                if missing or missing_markers:
                    entry.update(missing_outputs=missing, missing_success_markers=missing_markers)
                    state.update(state='failed')
                    publish()
                    break
                entry['outputs'] = [digest(Path(p)) for p in stage.get('required_outputs', [])]
                publish()
            else:
                state.update(state='complete', current_stage=None)
            state['peak_stage_rss_bytes'] = max(
                (s['measurement'].get('peak_summed_rss_bytes', 0) for s in state['stages']), default=0)
            publish()
        except Exception as error:
            state.update(state='failed', error=str(error))
            publish()
            raise
        return 0 if state['state'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
