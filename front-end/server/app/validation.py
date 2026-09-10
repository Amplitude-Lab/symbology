"""Portable file naming and validation of persisted graph structure."""

_WINDOWS_DEVICES = {'CON', 'PRN', 'AUX', 'NUL'} | {
    f'{prefix}{number}' for prefix in ('COM', 'LPT') for number in range(1, 10)
}


def portable_stem(name: str, prefix: str) -> str:
    """Protect an already sanitized ASCII name from Windows device aliases."""
    if name.split('.', 1)[0].upper() in _WINDOWS_DEVICES:
        return f'{prefix}-{name}'
    return name


def graph_errors(graph):
    if not isinstance(graph, dict) or not isinstance(graph.get('nodes'), list) or not isinstance(graph.get('edges'), list):
        return ['Graph must contain nodes and edges lists.']
    errors, ids, ports = [], set(), set()
    for n in graph['nodes']:
        if not isinstance(n, dict) or not isinstance(n.get('id'), str) or not n['id'] or not isinstance(n.get('type'), str) or not isinstance(n.get('data', {}), dict):
            errors.append('Every node needs a string id, type, and data object.')
            continue
        if n['id'] in ids: errors.append(f"Duplicate node id: {n['id']}")
        ids.add(n['id'])
    for e in graph['edges']:
        if not isinstance(e, dict) or any(not isinstance(e.get(k), str) for k in ('source', 'target', 'sourceHandle', 'targetHandle')):
            errors.append('Every connection needs source, target and port names.')
            continue
        if e['source'] not in ids or e['target'] not in ids:
            errors.append('A connection refers to a missing node.')
        port = (e['target'], e['targetHandle'])
        if port in ports: errors.append(f'Input port {port[0]}.{port[1]} has more than one connection.')
        ports.add(port)
    return errors
