#!/usr/bin/env python3
"""Record OOM/unclean service failures so boot cannot blindly repeat them."""
import os
from pathlib import Path
import time
from hexagon_boot_recovery import atomic_json

control = Path(__file__).resolve().parents[1] / 'output_hexagon_ten_v27d_20261009/runner'
result = {k: os.environ.get(k) for k in ('SERVICE_RESULT', 'EXIT_CODE', 'EXIT_STATUS')}
result['time'] = time.time()
relative = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
cgroup = Path('/sys/fs/cgroup') / relative.lstrip('/')
result['cgroup_counters'] = {name: (cgroup / name).read_text() for name in
                            ('memory.events', 'memory.peak', 'cpu.stat') if (cgroup / name).exists()}
atomic_json(control / 'last-service-exit.json', result)
atomic_json(control / 'audit' / ('service-exit-' + str(time.time_ns()) + '.json'), result)
if result['SERVICE_RESULT'] not in ('success', None) and not (control / 'attention-required.json').exists():
    atomic_json(control / 'attention-required.json', {
        'error': 'Service ended abnormally; inspect cgroup memory events and logs before retrying', **result})
