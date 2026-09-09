import io, ast

P = 'scripts/wire-audit.py'
s = io.open(P, encoding='utf-8').read()

# 1. drop stale pre-ROW_STEP ordering + old port_y
old = """ROW_STEP_FROM_DUMP = 15


def port_y(n, handle, side):
    lists = PORT_LISTS.get(n['type'])
    if not lists:
        return n['position']['y'] + est_h(n) / 2
    ports = lists[1] if side == 'out' else lists[0]
    row = ports.index(handle) if handle in ports else 0
    return n['position']['y'] + ROW0 + (SUB_EXTRA if n['type'] in SUB else 0) + row * ROW_STEP


"""
assert s.count(old) == 1
s = s.replace(old, '')

# 2. constant before use (cosmetic, module-level)
old = "CLEAR = 6\n"
new = "CLEAR = 6\nROW_STEP_FROM_DUMP = 15\n"
assert s.count(old) == 1
s = s.replace(old, new)
old = "def port_y(flow_name, n, handle, side):\n    p = ports(flow_name, n['id'])\n    rows = p['out'] if side == 'out' else p['in']\n    row = rows.index(handle) if handle in rows else 0\n    return n['position']['y'] + p['baseY'] + row * ROW_STEP_FROM_DUMP\n"
assert s.count(old) == 1
s = s.replace(old, old)
s = s.replace("    return n['position']['y'] + p['baseY'] + row * ROW_STEP_FROM_DUMP\n\n\nROW_STEP_FROM_DUMP = 15", "    return n['position']['y'] + p['baseY'] + row * ROW_STEP_FROM_DUMP")

# 3. flow name threading in audit()
old = "def audit(flow):\n    nodes ="
new = "def audit(flow):\n    fn = flow['name']\n    nodes ="
assert s.count(old) == 1
s = s.replace(old, new)

pairs = [
    ("n['position']['y'] + est_h(n) + 26", "n['position']['y'] + est_h(fn, n) + 26"),
    ("sy = port_y(sn, e.get('sourceHandle'), 'out')", "sy = port_y(fn, sn, e.get('sourceHandle'), 'out')"),
    ("ty = port_y(tn, e.get('targetHandle'), 'in')", "ty = port_y(fn, tn, e.get('targetHandle'), 'in')"),
    ("y0 = port_y(sn, e.get('sourceHandle'), 'out') + spread_of[edge_key(e)]", "y0 = port_y(fn, sn, e.get('sourceHandle'), 'out') + spread_of[edge_key(e)]"),
    ("y1 = port_y(tn, e.get('targetHandle'), 'in')", "y1 = port_y(fn, tn, e.get('targetHandle'), 'in')"),
    ("y_run = port_y(tn, e.get('targetHandle'), 'in')", "y_run = port_y(fn, tn, e.get('targetHandle'), 'in')"),
]
for old, new in pairs:
    n = s.count(old)
    assert n >= 1, old
    s = s.replace(old, new)
    print('replaced %dx: %s' % (n, old[:50]))

# 4. drop now-unused SUB set + re import if unreferenced
if 're.' not in s and 'import re' in s:
    s = s.replace('import re\n', '')
    print('dropped unused import re')
sub_start = s.find('SUB = {')
if sub_start != -1 and s.count('SUB') == 1:
    sub_end = s.index('}\n', sub_start) + 2
    s = s[:sub_start] + s[sub_end:]
    print('dropped unused SUB set')

ast.parse(s)
io.open(P, 'w', encoding='utf-8').write(s)
print('wire-audit.py rewritten, %d lines' % (s.count('\n') + 1))
