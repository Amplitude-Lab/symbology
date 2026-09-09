#!/usr/bin/env python3
"""One-off triage: categorize wire-audit findings by the router decision of
each implicated wire. Reads /tmp/wire-py.txt (WIRE_DEBUG dump) and the audit
stdout /tmp/wire-py-run.txt."""
import collections
import re

dec = {}
for ln in open('/tmp/wire-py.txt'):
    f = ln.rstrip('\n').split('|')
    if len(f) < 8:
        continue
    dec[(f[0].split('/')[-1], f[2], f[3], f[4], f[5])] = f[6]

cats = collections.Counter()
examples = collections.defaultdict(list)
flow = None
kinds = ('path passes through box', 'direct wires', 'detour runs share lane',
         'same-port outbound anchors', 'arrival not horizontal')
for ln in open('/tmp/wire-py-run.txt'):
    m = re.match(r'^([A-Za-z0-9_+\-]+): (\d+) problems', ln)
    if m:
        flow = m.group(1)
        continue
    if flow is None:
        continue
    kind = next((k for k in kinds if k in ln), None)
    if kind is None:
        continue
    tuples = re.findall(r'\(([^)]*)\)', ln)
    keys = []
    for t in tuples:
        for kk in re.findall(r"'([^']*)'", t):
            keys.append(kk)
    d1 = dec.get((flow, keys[0], keys[1], keys[2], keys[3]), '?')
    d2 = dec.get((flow, keys[4], keys[5], keys[6], keys[7]), '?') if len(keys) >= 8 else '-'
    cat = (kind, d1, d2)
    cats[cat] += 1
    if len(examples[cat]) < 2:
        examples[cat].append((flow, ln.strip()[:120]))

total = 0
for (kind, d1, d2), n in cats.most_common():
    total += n
    print('%4d  [%-8s|%-8s] %s' % (n, d1, d2, kind))
    for fl, e in examples[(kind, d1, d2)][:2]:
        print('        %s: %s' % (fl, e))
print('categorized:', total)
