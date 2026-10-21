"""Measure a command and its descendants on Linux; no external dependencies.

RSS is sampled every 100 ms and summed across processes (shared pages can be
counted more than once). GNU time in the command can also record per-process
high-water RSS. Run benchmarks sequentially, without concurrent builds.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('--report', required=True, type=Path)
p.add_argument('--log', required=True, type=Path)
p.add_argument('--timeout', type=float, default=900)
p.add_argument('--virtual-gib', type=int, default=28)
p.add_argument('--rss-gib', type=float, default=0,
               help='stop the command when its sampled process-tree RSS exceeds this limit')
p.add_argument('--status', type=Path,
               help='atomically write live process/resource status every ten seconds')
p.add_argument('--min-available-gib', type=float, default=0,
               help='stop this command if Linux MemAvailable drops below this value')
p.add_argument('command', nargs=argparse.REMAINDER)
a = p.parse_args()
command = a.command[1:] if a.command[:1] == ['--'] else a.command

def limits():
    size = a.virtual_gib * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (size, size))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

def rss_tree(root):
    pending, total, count = [root], 0, 0
    while pending:
        pid = pending.pop()
        try:
            fields = Path(f'/proc/{pid}/statm').read_text().split()
            total += int(fields[1]) * os.sysconf('SC_PAGE_SIZE')
            count += 1
            pending.extend(map(int, Path(f'/proc/{pid}/task/{pid}/children').read_text().split()))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return total, count


def live_group(pgid):
    """Include orphaned descendants after their group leader has exited."""
    live = []
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            fields = path.read_text().rsplit(') ', 1)[1].split()
            if int(fields[2]) == pgid and fields[0] not in ('Z', 'X'):
                live.append(int(path.parent.name))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return live


def stop_group(child):
    # Waiting for the leader alone is insufficient: it may exit before its
    # children, including children that deliberately ignore SIGTERM.
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(child.pid, sig)
        except ProcessLookupError:
            break
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            child.poll()
            if not live_group(child.pid):
                child.wait()
                return
            time.sleep(.01)
    child.wait(timeout=5)
    if live_group(child.pid):
        raise RuntimeError('owned process group did not terminate')

start = time.monotonic()
peak = max_processes = 0
timed_out = False
memory_pressure = False
rss_limit_exceeded = False
cancel_signal = 0
last_status = -10.0
usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)

def write_status(value):
    if a.status:
        temporary = a.status.with_name(a.status.name + '.tmp')
        temporary.write_text(json.dumps(value, indent=2) + '\n')
        temporary.replace(a.status)

def request_stop(signum, frame):
    global cancel_signal
    cancel_signal = signum

signal.signal(signal.SIGTERM, request_stop)
signal.signal(signal.SIGINT, request_stop)

with a.log.open('w') as log:
    child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                             start_new_session=True, preexec_fn=limits)
    try:
        while child.poll() is None:
            rss, count = rss_tree(child.pid)
            peak = max(peak, rss)
            max_processes = max(max_processes, count)
            timed_out = time.monotonic() - start > a.timeout
            rss_limit_exceeded = bool(a.rss_gib and rss > a.rss_gib * 1024**3)
            if a.min_available_gib:
                available = next(int(line.split()[1]) * 1024 for line in
                                 Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:'))
                memory_pressure = available < a.min_available_gib * 1024**3
            elapsed = time.monotonic() - start
            if elapsed - last_status >= 10:
                write_status(dict(state='running', pid=child.pid, command=command,
                                  wall_seconds=elapsed, current_rss_bytes=rss,
                                  peak_summed_rss_bytes=peak,
                                  timeout_seconds=a.timeout, rss_limit_gib=a.rss_gib,
                                  virtual_limit_gib=a.virtual_gib,
                                  updated_unix_seconds=time.time()))
                last_status = elapsed
            if timed_out or memory_pressure or rss_limit_exceeded or cancel_signal:
                break
            time.sleep(.1)
    finally:
        if child.poll() is None or timed_out or memory_pressure or rss_limit_exceeded or cancel_signal:
            stop_group(child)
    code = child.wait()
usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
result = dict(command=command, exit_code=code, timeout=timed_out,
              memory_pressure=memory_pressure,
              rss_limit_exceeded=rss_limit_exceeded, rss_limit_gib=a.rss_gib,
              cancel_signal=cancel_signal,
              min_available_gib=a.min_available_gib,
              wall_seconds=time.monotonic()-start,
              peak_summed_rss_bytes=peak, sample_interval_seconds=.1,
              cpu_user_seconds=usage_after.ru_utime-usage_before.ru_utime,
              cpu_system_seconds=usage_after.ru_stime-usage_before.ru_stime,
              max_processes=max_processes, virtual_limit_gib=a.virtual_gib)
a.report.write_text(json.dumps(result, indent=2)+'\n')
write_status(dict(result, state='finished', pid=child.pid,
                  updated_unix_seconds=time.time()))
print(json.dumps(result), flush=True)
raise SystemExit(0 if code == 0 and not timed_out and not memory_pressure and not rss_limit_exceeded and not cancel_signal else 1)
