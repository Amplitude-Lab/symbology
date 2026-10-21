#!/usr/bin/env python3
"""Record a read-only live snapshot; reported rank may lag the current rank."""
import argparse
import json
import re
import time
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('root', type=Path)
p.add_argument('run')
p.add_argument('--scratch', type=Path)
a = p.parse_args()
status = json.loads((a.root / 'audit' / (a.run + '.status.json')).read_text())
log = (a.root / 'log' / (a.run + '.log')).read_text()
entry = {'snapshot_unix_seconds': time.time(), 'status': status}
reported = re.findall(r'streaming_batch=(\d+) input_rows=(\d+) rank=(\d+) retained_nnz=(\d+)', log)
if reported:
    entry['last_reported_batch'] = dict(zip(('batch', 'input_rows', 'rank', 'retained_nnz'), map(int, reported[-1])))
refine = re.findall(r'REFINE_DISCOVERY limit=(\d+) rank=(\d+) nullity=(\d+) elapsed=([\d.]+)', log)
if refine:
    q = refine[-1]
    entry['last_reported_discovery'] = {'limit': int(q[0]), 'rank': int(q[1]), 'nullity': int(q[2]), 'elapsed': float(q[3])}
entry['completed_evictions'] = log.count('MAPPED_EVICT end')
proc = Path('/proc') / str(status['pid']) / 'status'
if proc.exists():
    entry['process_memory'] = {k: int(v) * 1024 for k, v in re.findall(r'^(VmRSS|VmSize|RssAnon|RssFile|RssShmem):\s+(\d+) kB', proc.read_text(), re.M)}
    stat = proc.with_name('stat').read_text().rsplit(')', 1)[1].split()
    entry['process_counters'] = {'minor_faults': int(stat[7]), 'major_faults': int(stat[9]), 'cpu_ticks': int(stat[11]) + int(stat[12])}
    entry['process_io'] = {k: int(v) for k, v in re.findall(r'^(read_bytes|write_bytes):\s+(\d+)', proc.with_name('io').read_text(), re.M)}
if a.scratch and a.scratch.exists():
    s = a.scratch.stat()
    entry['scratch_logical_bytes'] = s.st_size
    entry['scratch_allocated_bytes'] = s.st_blocks * 512
with (a.root / 'audit' / (a.run + '.resource-history.jsonl')).open('a') as f:
    f.write(json.dumps(entry) + '\n')
print(json.dumps({k: v for k, v in entry.items() if k != 'status'}))
print('elapsed_seconds=', round(status['wall_seconds'], 1), 'rss_bytes=', status.get('current_rss_bytes'), 'state=', status['state'])
