import json, re, sys

ROOT = '/Users/windfolgen/GitRepos/symbology'
p = json.load(open(f'{ROOT}/front-end/projects/heptagon/project.json'))
compile_src = open(f'{ROOT}/front-end/server/app/compile.py').read()
ui_src = open(f'{ROOT}/front-end/web/src/screens/FlowEditor.jsx').read()
defs_src = open(f'{ROOT}/front-end/web/src/flowdefs.js').read()

problems = []
mine = {'MHVw4flow', 'MHVw6flow'}

def err(flow, node, msg):
    name = flow['name'] if isinstance(flow, dict) else str(flow)
    nid = node.get('id', '?')
    ntype = node.get('type', '?')
    problems.append(f'{name}/{nid} ({ntype}): {msg}')

# Layer 1: contract walk
for f in p['flows']:
    g = f['graph']
    ids = [n['id'] for n in g['nodes']]
    if len(ids) != len(set(ids)):
        err(f, {'id': '(flow)'}, 'duplicate node ids')
    idset = set(ids)
    for e in g['edges']:
        for end, handles in (('source', 'sourceHandle'), ('target', 'targetHandle')):
            if e[end] not in idset:
                err(f, {'id': e[end]}, f'dangling edge {end}')
    for n in g['nodes']:
        if 'data' not in n:
            err(f, n, 'missing data object')
        d = n.get('data', {})
        if n['type'] == 'alphabet':
            if not any(a['id'] == d.get('alphabet_id') for a in p['alphabets']):
                err(f, n, f"alphabet_id {d.get('alphabet_id')} not found")
            else:
                props = {pr['id'] for a in p['alphabets'] if a['id'] == d['alphabet_id'] for pr in a['properties']}
                for sp in d.get('selected_properties', []):
                    if sp not in props:
                        err(f, n, f'selected property {sp} not in alphabet')
        if n['type'] == 'reuse_output' and d.get('file'):
            import os
            if not os.path.exists(f"{ROOT}/front-end/projects/heptagon/{d['file']}"):
                err(f, n, f"reuse file {d['file']} missing on disk")

# Layer 2: data-key liveness for MY flows (typo detection)
for f in p['flows']:
    if f['name'] not in mine:
        continue
    for n in f['graph']['nodes']:
        d = n.get('data', {})
        for k in d:
            pat = f'"{k}"'
            dotted = re.search(rf'\.{k}\b|data\?\.{k}\b', ui_src + defs_src + compile_src)
            if pat not in compile_src and pat not in ui_src and pat not in defs_src and not dotted:
                err(f, n, f'data key "{k}" read by neither compiler nor UI (typo?)')

# Layer 2b: per-position sanity of my flows' semantic fields
for f in p['flows']:
    if f['name'] not in mine:
        continue
    for n in f['graph']['nodes']:
        d = n.get('data', {})
        if n['type'] == 'extend':
            tw = d.get('target_weight')
            if not isinstance(tw, int) or tw < 2:
                err(f, n, f'target_weight={tw!r} invalid')
        if n['type'] == 'tensor_dot':
            for a in ('axis_a', 'axis_b'):
                v = d.get(a)
                if not (isinstance(v, int) or (isinstance(v, str) and v.lstrip('-').isdigit())):
                    err(f, n, f'{a}={v!r} not an integer or numeric string')

print(f'contract+typo audit: {len(problems)} problems')
for x in problems:
    print(' ', x)
sys.exit(1 if problems else 0)
