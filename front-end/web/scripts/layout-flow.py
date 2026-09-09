#!/usr/bin/env python3
"""Auto-layout a flow in front-end/projects/heptagon/project.json.

House style: longest-path layering (sources left), barycenter ordering,
per-node true heights from the OpNode/AlphabetNode geometry, corridor
pitch DX wide enough for arrowheads + kind-colored wires to be readable.

Usage: python3 layout-flow.py <flow-name> [--dry-run]
Python 3.9 compatible. Writes back position-only changes.
"""
import json
import re
import sys

ROOT = '/Users/windfolgen/GitRepos/symbology'
PROJ = f'{ROOT}/front-end/projects/heptagon/project.json'
DEFS = open(f'{ROOT}/front-end/web/src/flowdefs.js').read()

DX = 260      # column pitch: node 190 + ~70 corridor for wires/arrowheads
GAP = 65      # vertical gap between nodes in a column
NODE_W = 190

# ---- node heights (measured from FlowEditor.jsx card geometry) ----
OP_BASE = 98          # title + padding + one row
SUB = 15              # subtitle line
ROW = 15
ALPHA_BASE = 62       # alphabet: body pad + name header + one row
IO_BASE = 60          # cb_out / reuse_output: title + one row


def defs_ports():
    """Parse static input/output counts per type from NODE_DEFS."""
    counts = {}
    for m in re.finditer(r'(\w+):\s*\{\s*title:[^{}]*inputs:\s*\[([^\]]*)\][^{}]*outputs:\s*\[([^\]]*)\]', DEFS, re.S):
        typ, ins, outs = m.group(1), m.group(2), m.group(3)
        counts[typ] = (len(re.findall(r"id:\s*'", ins)), len(re.findall(r"id:\s*'", outs)))
    return counts


STATIC = defs_ports()


def subtitle_types():
    return {'extend', 'project', 'solve_symmetry', 'symderive', 'sew', 'solve_collinear',
            'projection_chain', 'symmetry_invariant', 'compute_rhs', 'add_tensors',
            'ternary_contract', 'apply_symmetry', 'matrix_power', 'tensor_join',
            'tensor_dot', 'squeeze_tensor', 'shuffle_product', 'expand_tensor',
            'impose_integrability', 'integrability_condition', 'solve_conditions'}


def alphabet_rows(node, project):
    alpha = next((a for a in project['alphabets'] if a['id'] == node['data'].get('alphabet_id')), None)
    if not alpha:
        return 1
    kinds = {p['id']: p['type'] for p in alpha.get('properties', [])}
    n = 0
    for sp in node['data'].get('selected_properties', []):
        n += 1
        if kinds.get(sp) in ('first_entry', 'last_entry'):
            n += 1  # proj map slot
    return max(n, 1)


def node_rows(node, edges, project):
    t = node['type']
    if t == 'alphabet':
        return alphabet_rows(node, project), False
    if t in ('cb_out', 'cb_in', 'reuse_output'):
        return 1, False
    ins, outs = STATIC.get(t, (2, 1))
    wired_out = {e['sourceHandle'] for e in edges if e['source'] == node['id']}
    wired_in = {e.get('targetHandle') for e in edges if e['target'] == node['id']}
    if t in ('add_tensors', 'expand_tensor'):
        mx = 1
        for e in edges:
            if e['target'] == node['id']:
                m = re.match(r'in_(\d+)', e.get('targetHandle') or '')
                if m:
                    mx = max(mx, int(m.group(1)))
        ins = max(ins, mx + 2)
    if t == 'solve_collinear':
        d = node['data']
        pairs = d.get('pairs') or []
        extras = max(len(pairs) - 1, 0)
        mx = 0
        for e in edges:
            if e['target'] == node['id']:
                m = re.match(r'in_(?:seed|rhs)_(\d+)$', e.get('targetHandle') or '')
                if m:
                    mx = max(mx, int(m.group(1)))
        ins = 2 + (1 if (d.get('cond_enabled') or any(h == 'cond' for h in wired_in)) else 0) + max(extras, mx)
        outs = 2
    if t == 'projection_chain':
        basis = len([h for h in wired_out if re.match(r'basis_(last_)?w\d+$', h)])
        outs = 1 + basis
    if t == 'assemble':
        mx = -1
        for e in edges:
            if e['target'] == node['id']:
                m = re.match(r'in_(\d+)', e.get('targetHandle') or '')
                if m:
                    mx = max(mx, int(m.group(1)))
        ins = max(ins, mx + 2)
    return max(ins, outs, 1), t in subtitle_types()


def node_height(node, edges, project):
    rows, sub = node_rows(node, edges, project)
    t = node['type']
    if t == 'alphabet':
        return ALPHA_BASE + (rows - 1) * ROW
    if t in ('cb_out', 'cb_in', 'reuse_output'):
        return IO_BASE
    return OP_BASE + (rows - 1) * ROW + (SUB if sub else 0)


def layout(flow, project):
    nodes = flow['graph']['nodes']
    edges = flow['graph']['edges']
    ids = [n['id'] for n in nodes]
    succ = {i: [] for i in ids}
    pred = {i: [] for i in ids}
    for e in edges:
        succ[e['source']].append(e['target'])
        pred[e['target']].append(e['source'])
    # longest-path layering (iterative, cycle-safe fallback)
    layer = {i: 0 for i in ids}
    for _ in range(len(ids) + 1):
        changed = False
        for i in ids:
            l = max([layer[p] + 1 for p in pred[i]] + [0])
            if l > layer[i]:
                layer[i] = l
                changed = True
        if not changed:
            break
    # results rail: every cb_out sink goes to one terminal column BEYOND all
    # producers (stacked in flow order) so every wire strictly points
    # left→right into a single "outputs" rail on the right.
    producers = [i for i in ids if succ[i]]
    if producers:
        terminal = max(layer[i] for i in producers) + 1
        for n in nodes:
            if not succ[n['id']] and n['type'] == 'cb_out':
                layer[n['id']] = terminal
    cols = {}
    for i in ids:
        cols.setdefault(layer[i], []).append(i)
    # initial order: stable by original JSON order
    order = {c: sorted(members, key=ids.index) for c, members in cols.items()}
    # barycenter sweeps
    for direction in (1, 0, 1, 0):
        for c in sorted(cols, reverse=direction == 0):
            members = order[c]
            neigh = cols.get(c - 1 if direction else c + 1, [])
            if not neigh:
                continue
            pos = {n: k for k, n in enumerate(neigh)}
            bary = {}
            for n in members:
                if direction:
                    srcs = pred[n]
                else:
                    srcs = succ[n]
                ps = [pos.get(s) for s in srcs if s in pos]
                bary[n] = sum(ps) / len(ps) if ps else pos.get(members.index(n) and 0 or 0, 0)
            order[c] = sorted(members, key=lambda n: (bary[n], ids.index(n)))
    # positions: stack per column, then center columns on common midline
    H = {n['id']: node_height(n, edges, project) for n in nodes}
    pos = {}
    tops = {}
    for c, members in order.items():
        y = 0.0
        tops[c] = {}
        for n in members:
            tops[c][n] = y
            y += H[n] + GAP
    all_mid = []
    for c, members in order.items():
        first, last = members[0], members[-1]
        all_mid.append((tops[c][first] + tops[c][last] + H[last]) / 2)
    mid = sum(all_mid) / len(all_mid)
    node_by_id = {n['id']: n for n in nodes}
    for c, members in order.items():
        first, last = members[0], members[-1]
        cmid = (tops[c][first] + tops[c][last] + H[last]) / 2
        dy = mid - cmid
        for n in members:
            pos[n] = {'x': c * DX, 'y': round(tops[c][n] + dy)}
    return pos, H


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    dry = '--dry-run' in sys.argv
    if not args:
        print('usage: layout-flow.py <flow-name> [--dry-run]')
        sys.exit(2)
    project = json.load(open(PROJ))
    flow = next((f for f in project['flows'] if f['name'] == args[0]), None)
    if not flow:
        print('flow not found: ' + args[0])
        sys.exit(2)
    before = {n['id']: dict(n['position']) for n in flow['graph']['nodes']}
    pos, H = layout(flow, project)
    changed = 0
    for n in flow['graph']['nodes']:
        if before[n['id']] != pos[n['id']]:
            changed += 1
        n['position'] = pos[n['id']]
    xs = [p['x'] for p in pos.values()]
    ys = [p['y'] for p in pos.values()]
    ext = (max(xs) + NODE_W - min(xs), max(ys) + max(H.values()) - min(ys))
    print('%s: %d nodes moved, extent %dx%d' % (flow['name'], changed, ext[0], ext[1]))
    if not dry and changed:
        json.dump(project, open(PROJ, 'w'), indent=1)
        print('written (position-only changes)')
    # overlap self-check within columns
    bycol = {}
    for n in flow['graph']['nodes']:
        bycol.setdefault(n['position']['x'], []).append(n)
    bad = 0
    for x, ns in bycol.items():
        ns = sorted(ns, key=lambda n: n['position']['y'])
        for a, b in zip(ns, ns[1:]):
            if b['position']['y'] < a['position']['y'] + H[a['id']]:
                bad += 1
                print('OVERLAP x=%s: %s / %s' % (x, a['id'], b['id']))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
