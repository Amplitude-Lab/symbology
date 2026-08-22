from __future__ import annotations

import hashlib
import re
import shlex
import threading
from pathlib import Path

from . import storage
from .config import find_bootstrap, find_compute_rhs, find_tensor_add, find_tensor_ops, find_wolframscript
from .wolfram import PROP_KIND, merge_script, property_display_name, property_script, property_tensor_relpath

TENSOR_KINDS = {"dlogmat", "fec1", "fec", "lec1", "lec", "sew", "matrix", "basis", "solution", "boundary", "tensor"}
R3_KINDS = TENSOR_KINDS

EDGE_RULES = {
    ("merge_conditions", None): {"dlogmat"},
    ("extend", "condition"): {"dlogmat"},
    ("extend", "fec"): {"fec1", "fec"},
    ("extend", "lec"): {"lec1", "lec"},
    ("sew", "condition"): {"dlogmat"},
    ("sew", "fec"): {"fec1", "fec"},
    ("sew", "lec"): {"lec1", "lec"},
    ("project", "tensor"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("project", "rep"): {"matrix"},
    ("project", "map"): {"matrix"},
    ("solve_symmetry", "tensor"): R3_KINDS,
    ("solve_symmetry", "matrix"): {"matrix"},
    ("solve_symmetry", "matrix2"): {"matrix"},
    ("symderive", "tensor"): R3_KINDS,
    ("symderive", "matrix"): {"matrix"},
    ("symderive", "matrix2"): {"matrix"},
    ("solve_collinear", "seed"): {"fec1", "fec", "sew"},
    ("projection_chain", "seed"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("symmetry_invariant", "seed"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("compute_rhs", "seed"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("add_tensors", None): TENSOR_KINDS,
    ("ternary_contract", "tensor"): R3_KINDS,
    ("ternary_contract", "trans1"): {"matrix"},
    ("ternary_contract", "trans2"): {"matrix"},
    ("apply_symmetry", "tensor"): R3_KINDS,
    ("apply_symmetry", "trans1"): {"matrix"},
    ("apply_symmetry", "trans2"): {"matrix"},
    ("matrix_power", "matrix"): {"matrix"},
    ("tensor_join", "a"): TENSOR_KINDS,
    ("tensor_join", "b"): TENSOR_KINDS,
    ("tensor_dot", "a"): TENSOR_KINDS,
    ("tensor_dot", "b"): TENSOR_KINDS,
    ("impose_integrability", "tensor"): R3_KINDS,
    ("impose_integrability", "dlogmat"): {"dlogmat"},
    ("integrability_condition", "tensor"): R3_KINDS,
    ("integrability_condition", "dlogmat"): {"dlogmat"},
    ("solve_conditions", "cond"): {"matrix", "tensor"},
    ("assemble", None): TENSOR_KINDS,
}


RAT_RE = re.compile(r"^-?\d+(/\d+)?$")


def _abs(proj_dir: Path, rel: str) -> str:
    p = Path(rel)
    return str(p if p.is_absolute() else proj_dir / rel)


def _resolve_path_arg(proj_dir: Path, value: str) -> str:
    if value in ("0", "identity"):
        return value
    return _abs(proj_dir, value)


def _squeeze_dims(dims):
    return [d for d in dims if d != 1] if dims else None


def _join_dims(a, b, axis):
    if not a or not b or len(a) != len(b):
        return None
    k = axis - 1 if axis > 0 else len(a) + axis
    if k < 0 or k >= len(a) or any(x != y for i, (x, y) in enumerate(zip(a, b)) if i != k):
        return None
    out = list(a)
    out[k] += b[k]
    return out


def _write_gen_script(path: Path, content: str) -> None:
    try:
        if path.exists() and path.read_text(encoding="utf-8", errors="replace") == content:
            return
    except OSError:
        pass
    path.write_text(content, encoding="utf-8")


_out_tls = threading.local()


def _out() -> str:
    return getattr(_out_tls, "prefix", "output/")


def compile_flow(proj: dict, graph: dict, output_subdir: str | None = None) -> dict:
    errors: list[str] = []
    steps: list[dict] = []
    proj_dir = storage.project_dir(proj["id"])
    gen_dir = proj_dir / "wolfram_gen"
    gen_dir.mkdir(parents=True, exist_ok=True)

    subdir = (output_subdir or "").strip().strip("/")
    _out_tls.prefix = f"output/{subdir}/" if subdir else "output/"

    step_counter = [0]
    file_steps: set = set()
    provides: dict[tuple, dict] = {}

    def next_step_id() -> str:
        step_counter[0] += 1
        return f"step-{step_counter[0]}"

    def add_step(label, kind, argv, cwd, outputs, skip_if_exists, meta=None):
        steps.append({
            "id": next_step_id(),
            "label": label,
            "kind": kind,
            "command": shlex.join([str(a) for a in argv]),
            "argv": [str(a) for a in argv],
            "cwd": str(cwd) if cwd else None,
            "outputs": outputs,
            "skip_if_exists": skip_if_exists,
            "meta": meta or {},
        })

    ctx = {
        "proj": proj,
        "errors": errors,
        "provides": provides,
        "add_step": add_step,
        "file_steps": file_steps,
        "proj_dir": proj_dir,
        "gen_dir": gen_dir,
        "wolframscript": find_wolframscript(),
        "bootstrap": find_bootstrap(),
        "compute_rhs_bin": find_compute_rhs(),
        "tensor_add_bin": find_tensor_add(),
        "tensor_ops_bin": find_tensor_ops(),
        "stack": [],
        "sew_names": set(),
    }

    nodes = {n["id"]: n for n in graph.get("nodes", [])}
    if not nodes:
        return {"ok": False, "errors": ["The flow is empty. Add at least one node."], "steps": []}

    _compile_graph(ctx, list(nodes.values()), graph.get("edges", []), {})

    edges = graph.get("edges", [])
    for e in edges:
        src = (e.get("source"), e.get("sourceHandle"))
        if src not in provides:
            continue
        kind = provides[src]["kind"]
        tnode = nodes.get(e.get("target"))
        if not tnode:
            continue
        ttype = tnode.get("type")
        th = e.get("targetHandle") or ""
        if ttype in ("merge_conditions", "assemble", "add_tensors"):
            allowed = EDGE_RULES[(ttype, None)]
        else:
            allowed = EDGE_RULES.get((ttype, th))
            if allowed is None and th.startswith("seed"):
                allowed = EDGE_RULES.get((ttype, "seed"))
        if allowed is None:
            continue
        if kind not in allowed:
            errors.append(
                f"Type mismatch: a '{kind}' output cannot feed the '{th}' input of a {ttype} node."
            )

    if errors:
        return {"ok": False, "errors": errors, "steps": []}
    shapes = {f"{nid}|{handle}": info.get("dims") for (nid, handle), info in provides.items() if info.get("dims")}
    flow_outputs = []
    for n in nodes.values():
        if n.get("type") == "cb_out":
            info = provides.get((n["id"], "in"))
            if info and info.get("file"):
                flow_outputs.append({
                    "name": (n.get("data") or {}).get("name") or f"out{n['id'][:6]}",
                    "file": info["file"],
                    "kind": info.get("kind"),
                    "dims": info.get("dims"),
                })
    return {"ok": True, "errors": [], "shapes": shapes, "flow_outputs": flow_outputs,
            "steps": [{k: v for k, v in s.items() if k != "argv"} for s in steps], "_steps_full": steps}


def _compile_graph(ctx, node_list, edges, seed):
    """Compile one graph (main flow or an expanded custom block).

    `seed` maps (cb_in node id, 'out') -> provides entry resolved from the
    instance's incoming wires, so abstract inputs become concrete tensors.
    """
    errors = ctx["errors"]
    provides = ctx["provides"]
    add_step = ctx["add_step"]
    proj_dir = ctx["proj_dir"]
    gen_dir = ctx["gen_dir"]
    wolframscript = ctx["wolframscript"]
    bootstrap = ctx["bootstrap"]
    compute_rhs_bin = ctx["compute_rhs_bin"]
    tensor_add_bin = ctx["tensor_add_bin"]
    tensor_ops_bin = ctx["tensor_ops_bin"]
    file_steps = ctx["file_steps"]

    nodes = {n["id"]: n for n in node_list}
    incoming: dict[str, list] = {nid: [] for nid in nodes}
    outgoing: dict[str, list] = {nid: [] for nid in nodes}
    for e in edges:
        if e.get("source") in nodes and e.get("target") in nodes:
            incoming[e["target"]].append(e)
            outgoing[e["source"]].append(e)

    order = _topo_sort(nodes, incoming, outgoing)
    if order is None:
        errors.append("The flow contains a cycle. Break the loop and try again.")
        return

    for nid in order:
        node = nodes[nid]
        ntype = node.get("type")

        if ntype == "alphabet":
            _compile_alphabet(ctx["proj"], node, provides, add_step, errors, wolframscript, gen_dir, proj_dir, tensor_ops_bin, file_steps)
        elif ntype == "merge_conditions":
            _compile_merge(node, incoming, provides, add_step, errors, wolframscript, gen_dir, proj_dir)
        elif ntype == "extend":
            _compile_extend(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "sew":
            _compile_sew(node, incoming, provides, add_step, errors, bootstrap, proj_dir, ctx["sew_names"])
        elif ntype == "project":
            _compile_project(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "solve_symmetry":
            _compile_solve_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "symderive":
            _compile_symderive(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "solve_collinear":
            _compile_solve_collinear(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "projection_chain":
            _compile_projection_chain(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "symmetry_invariant":
            _compile_symmetry_invariant(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "compute_rhs":
            _compile_compute_rhs(node, incoming, provides, add_step, errors, compute_rhs_bin, proj_dir)
        elif ntype == "add_tensors":
            _compile_add_tensors(node, incoming, provides, add_step, errors, tensor_add_bin, proj_dir)
        elif ntype in ("ternary_contract", "apply_symmetry"):
            _compile_ternary_contract(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "matrix_power":
            _compile_matrix_power(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "tensor_join":
            _compile_tensor_join(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "tensor_dot":
            _compile_tensor_dot(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "impose_integrability":
            _compile_impose(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "integrability_condition":
            _compile_integrability_condition(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "solve_conditions":
            _compile_solve_conditions(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "assemble":
            _compile_assemble(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "customblock":
            _compile_custom_block(ctx, node, incoming)
        elif ntype == "cb_in":
            info = seed.get((nid, "out"))
            if info is None:
                if not ctx["stack"]:
                    errors.append(f"Input port '{node.get('data', {}).get('name') or nid}' belongs to a custom-block definition. Compile a flow that uses this block as an instance instead of compiling the definition itself.")
                else:
                    errors.append(f"Custom-block input port '{node.get('data', {}).get('name') or nid}' is not wired to the instance.")
            else:
                provides[(nid, "out")] = dict(info)
        elif ntype == "cb_out":
            src = None
            for e in incoming.get(nid, []):
                s = (e.get("source"), e.get("sourceHandle"))
                if s in provides:
                    src = provides[s]
            if src is None:
                errors.append(f"Custom-block output port '{node.get('data', {}).get('name') or nid}' has no wired input.")
            else:
                provides[(nid, "in")] = dict(src)
        else:
            errors.append(f"Unknown node type '{ntype}'.")


def _compile_custom_block(ctx, node, incoming):
    """Expand one custom-block instance: compile its definition flow inline."""
    errors = ctx["errors"]
    provides = ctx["provides"]
    block_id = (node.get("data") or {}).get("block")
    blk = next((f for f in ctx["proj"].get("flows", []) if f.get("id") == block_id and f.get("custom_block")), None)
    if blk is None:
        errors.append(f"Custom block '{node['id']}': definition flow not found (deleted?).")
        return
    if block_id in ctx["stack"]:
        errors.append(f"Custom block '{blk.get('name')}' contains itself (recursively) — that is not allowed.")
        return
    if len(ctx["stack"]) >= 8:
        errors.append("Custom blocks are nested more than 8 levels deep — refusing to expand.")
        return

    inner_nodes = (blk.get("graph") or {}).get("nodes", [])
    inner_edges = (blk.get("graph") or {}).get("edges", [])
    cb_ins = {n["id"]: n for n in inner_nodes if n.get("type") == "cb_in"}

    # instance wires -> abstract inputs
    seed = {}
    for e in incoming.get(node["id"], []):
        port = e.get("targetHandle")
        src = (e.get("source"), e.get("sourceHandle"))
        if port in cb_ins and src in provides:
            seed[(port, "out")] = provides[src]

    # namespace inner ids and target names so parallel instances never collide
    prefix = f"cb{node['id'].split('_')[-1]}_{hashlib.sha1(node['id'].encode()).hexdigest()[:6]}_"
    nodes_p = []
    for n in inner_nodes:
        m = {**n, "id": prefix + n["id"]}
        d = m.get("data")
        if isinstance(d, dict) and m.get("type") not in ("cb_in", "cb_out") and d.get("target"):
            m = {**m, "data": {**d, "target": f"{prefix}{d['target']}"}}
        nodes_p.append(m)
    edges_p = [{**e, "source": prefix + e["source"], "target": prefix + e["target"]} for e in inner_edges]
    seed_p = {(prefix + k[0], k[1]): v for k, v in seed.items()}

    ctx["stack"].append(block_id)
    _compile_graph(ctx, nodes_p, edges_p, seed_p)
    ctx["stack"].pop()

    # abstract outputs -> instance ports
    for n in inner_nodes:
        if n.get("type") == "cb_out":
            info = provides.get((prefix + n["id"], "in"))
            if info is not None:
                provides[(node["id"], n["id"])] = dict(info)


def _topo_sort(nodes, incoming, outgoing):
    indeg = {nid: len(incoming[nid]) for nid in nodes}
    queue = [nid for nid in nodes if indeg[nid] == 0]
    order = []
    while queue:
        nid = queue.pop(0)
        order.append(nid)
        for e in outgoing[nid]:
            t = e["target"]
            indeg[t] -= 1
            if indeg[t] == 0:
                queue.append(t)
    return order if len(order) == len(nodes) else None


def _edge_input(provides, incoming, nid, handle_prefix):
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        if th == handle_prefix or th.startswith(handle_prefix):
            src = (e.get("source"), e.get("sourceHandle"))
            if src in provides:
                return provides[src]
    return None


def _has_incoming(incoming, nid, handle_prefix) -> bool:
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        if th == handle_prefix or th.startswith(handle_prefix):
            return True
    return False


def _compile_alphabet(proj, node, provides, add_step, errors, wolframscript, gen_dir, proj_dir, tensor_ops_bin=None, file_steps=None):
    if file_steps is None:
        file_steps = set()
    nid = node["id"]
    data = node.get("data", {}) or {}
    alphabet = storage.find_alphabet(proj, data.get("alphabet_id", ""))
    if alphabet is None:
        errors.append("An alphabet node has no alphabet selected.")
        return
    selected = data.get("selected_properties") or []
    for prop_id in selected:
        prop = storage.find_property(alphabet, prop_id)
        if prop is None:
            errors.append(f"Alphabet '{alphabet['name']}': selected property not found.")
            continue
        kind = PROP_KIND.get(prop["type"])
        if kind is None:
            errors.append(f"Alphabet '{alphabet['name']}': property has unknown type '{prop['type']}'.")
            continue
        rel = prop.get("tensor_file") or property_tensor_relpath(alphabet, prop)
        if prop.get("status") != "ready" or not (proj_dir / rel).exists():
            if prop["type"] == "precomputed_tensor":
                errors.append(f"Alphabet '{alphabet['name']}': precomputed tensor file '{rel}' is missing from the project.")
                continue
            if wolframscript is None:
                errors.append("wolframscript was not found; cannot compute properties.")
                continue
            out_abs = _abs(proj_dir, rel)
            try:
                script = property_script(alphabet, prop, out_abs)
            except ValueError as exc:
                errors.append(f"Alphabet '{alphabet['name']}': {exc}")
                continue
            script_path = gen_dir / f"prop_{prop_id}.wl"
            _write_gen_script(script_path, script)

            add_step(
                f"Compute {property_display_name(prop)} ({prop['type'].replace('_', ' ')}) for alphabet '{alphabet['name']}'",
                "wolfram",
                [wolframscript, "-script", str(script_path)],
                proj_dir,
                [rel],
                True,
                {"type": "property", "alphabet_id": alphabet["id"], "property_id": prop_id, "tensor_file": rel},
            )
        base_info = {
            "kind": kind,
            "file": rel,
            "weight": 1 if kind in ("fec1", "lec1") else None,
            "name": Path(rel).stem,
            "dims": (prop.get("summary") or {}).get("dims"),
        }
        provides[(nid, f"prop_{prop_id}")] = dict(base_info)
        if prop["type"] in ("first_entry", "last_entry"):
            proj_rel = f"data/{Path(rel).stem}_proj.wxf"
            if not (proj_dir / proj_rel).exists() and proj_rel not in file_steps:
                file_steps.add(proj_rel)
                if tensor_ops_bin is None:
                    errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
                    continue
                add_step(
                    f"Derive projection map {Path(rel).stem} -> {Path(proj_rel).stem} (drop size-1 axis)",
                    "tensor_ops",
                    [tensor_ops_bin, "squeeze", _abs(proj_dir, rel), _abs(proj_dir, proj_rel)],
                    proj_dir,
                    [proj_rel],
                    True,
                    {"type": "squeeze", "tensor_file": proj_rel},
                )
            provides[(nid, f"proj_{prop_id}")] = {
                "kind": "matrix",
                "file": proj_rel,
                "weight": None,
                "name": Path(proj_rel).stem,
                "dims": _squeeze_dims((prop.get("summary") or {}).get("dims")),
            }


def _compile_merge(node, incoming, provides, add_step, errors, wolframscript, gen_dir, proj_dir):
    nid = node["id"]
    inputs = []
    for e in incoming.get(nid, []):
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            info = provides[src]
            if info["kind"] != "dlogmat":
                errors.append("Merge Conditions only accepts dlogmat inputs.")
                continue
            inputs.append(info["file"])
    if len(inputs) < 2:
        errors.append("Merge Conditions needs at least two dlogmat inputs.")
        return
    if wolframscript is None:
        errors.append("wolframscript was not found; cannot merge condition tensors.")
        return
    rel = f"data/merged_{nid[:8]}.wxf"
    script = merge_script([_abs(proj_dir, f) for f in inputs], _abs(proj_dir, rel))
    script_path = gen_dir / f"merge_{nid[:8]}.wl"
    _write_gen_script(script_path, script)
    add_step(
        f"Merge {len(inputs)} condition tensors",
        "wolfram",
        [wolframscript, "-script", str(script_path)],
        proj_dir,
        [rel],
        False,
        {"type": "merge", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "dlogmat", "file": rel, "weight": None, "name": Path(rel).stem}


def _compile_extend(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    cond = _edge_input(provides, incoming, nid, "condition")
    fec = _edge_input(provides, incoming, nid, "fec")
    lec = _edge_input(provides, incoming, nid, "lec")
    if cond is None:
        errors.append("An Extend node is missing its condition (dlogmat) input.")
        return
    if fec is not None and lec is not None:
        errors.append("An Extend node takes either a FEC input or a LEC input, not both.")
        return
    if fec is None and lec is None:
        errors.append("An Extend node is missing its FEC or LEC input.")
        return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found; build it with `make bootstrap`.")
        return
    seed = fec or lec
    direction = "FEC" if fec is not None else "LEC"
    flag = "-f" if fec is not None else "-l"
    w_in = seed.get("weight")
    if w_in is None:
        errors.append(f"Cannot determine the weight of the {direction} input to an Extend node.")
        return
    w_out = w_in + 1
    target = data.get("target_weight")
    if target not in (None, "", 0):
        try:
            target = int(target)
        except (TypeError, ValueError):
            errors.append(f"Extend node target weight '{target}' is not an integer.")
            return
        if target != w_out:
            errors.append(f"Extend node target weight {target} does not match input weight {w_in} + 1.")
            return
    rel = f"{_out()}{direction}_{w_out}.wxf"
    add_step(
        f"Extend {direction}_{w_in} -> {direction}_{w_out}",
        "bootstrap",
        [bootstrap, "--extend", "-c", _abs(proj_dir, cond["file"]), flag, _abs(proj_dir, seed["file"]), "-o", _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "extend", "tensor_file": rel},
    )
    provides[(nid, "fec" if fec is not None else "lec")] = {
        "kind": "fec" if fec is not None else "lec",
        "file": rel,
        "weight": w_out,
        "name": f"{direction}_{w_out}",
    }


def _compile_sew(node, incoming, provides, add_step, errors, bootstrap, proj_dir, sew_names):
    nid = node["id"]
    cond = _edge_input(provides, incoming, nid, "condition")
    fec = _edge_input(provides, incoming, nid, "fec")
    lec = _edge_input(provides, incoming, nid, "lec")
    for name, val in (("condition", cond), ("FEC", fec), ("LEC", lec)):
        if val is None:
            errors.append(f"A Sew node is missing its {name} input.")
            return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found; build it with `make bootstrap`.")
        return
    fw, lw = fec.get("weight"), lec.get("weight")
    if fw is None or lw is None:
        errors.append("Cannot determine FEC/LEC weights for a Sew node.")
        return
    base = (node.get("data") or {}).get("target") or f"SEW_{fw}p{lw}"
    name = base
    counter = 2
    while name in sew_names:
        name = f"{base}_{counter}"
        counter += 1
    sew_names.add(name)
    rel = f"{_out()}{name}.wxf"
    add_step(
        f"Sew FEC_{fw} + LEC_{lw} -> {name}",
        "bootstrap",
        [bootstrap, "--sew", "-c", _abs(proj_dir, cond["file"]), "-f", _abs(proj_dir, fec["file"]),
         "-l", _abs(proj_dir, lec["file"]), "-o", _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "sew", "tensor_file": rel},
    )
    provides[(nid, "sew")] = {"kind": "sew", "file": rel, "weight": fw + lw, "name": name}


def _derive_target(node, incoming, provides, data):
    if data.get("target"):
        return data["target"]
    seed = _edge_input(provides, incoming, node["id"], "seed")
    if seed and seed.get("name"):
        return seed["name"]
    return None


def _compile_project(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    rep = _edge_input(provides, incoming, nid, "rep")
    mapm = _edge_input(provides, incoming, nid, "map")
    if tensor is None or rep is None or mapm is None:
        missing = [h for h, v in (("tensor", tensor), ("rep", rep), ("map", mapm))
                   if v is None and not _has_incoming(incoming, nid, h)]
        if missing:
            errors.append(f"A Project node is missing its {' and '.join(missing)} input{'s' if len(missing) > 1 else ''}.")
        return
    if tensor.get("file") is None or rep.get("file") is None or mapm.get("file") is None:
        errors.append("Project: tensor/rep/map inputs must be concrete tensor files.")
        return
    target = _require_target(data, "Project", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Project {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "project", _abs(proj_dir, tensor["file"]), _abs(proj_dir, rep["file"]),
         _abs(proj_dir, mapm["file"]), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "project", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"), "name": target}


def _compile_solve_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    matrix = _edge_input(provides, incoming, nid, "matrix")
    matrix2 = _edge_input(provides, incoming, nid, "matrix2")
    if tensor is None or matrix is None:
        errors.append("A Symmetry Solve node is missing its tensor or sym M input.")
        return
    if tensor.get("file") is None or matrix.get("file") is None:
        errors.append("Symmetry Solve: tensor / sym M inputs must be concrete tensor files.")
        return
    sc = matrix2["file"] if (matrix2 is not None and matrix2.get("file") is not None) else "-"
    target = _require_target(data, "Symmetry Solve", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Symmetry solve {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "symsolve", _abs(proj_dir, tensor["file"]), _abs(proj_dir, matrix["file"]),
         sc if sc == "-" else _abs(proj_dir, sc), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "solve_symmetry", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"), "name": target}


def _compile_symderive(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    matrix = _edge_input(provides, incoming, nid, "matrix")
    matrix2 = _edge_input(provides, incoming, nid, "matrix2")
    if tensor is None or matrix is None:
        errors.append("A Symmetry Derive node is missing its tensor or sym M input.")
        return
    if tensor.get("file") is None or matrix.get("file") is None:
        errors.append("Symmetry Derive: tensor / sym M inputs must be concrete tensor files.")
        return
    sc = matrix2["file"] if (matrix2 is not None and matrix2.get("file") is not None) else "-"
    target = _require_target(data, "Symmetry Derive", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Symmetry derive R of {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "symderive", _abs(proj_dir, tensor["file"]), _abs(proj_dir, matrix["file"]),
         sc if sc == "-" else _abs(proj_dir, sc), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "symderive", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "name": target}


def _compile_solve_collinear(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    data = node.get("data", {}) or {}
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Solve Collinear node: provide a target (e.g. SEW_3p1).")
        return
    rhs = data.get("rhs", "")
    if not rhs:
        errors.append("Solve Collinear node: provide an RHS file (or '0').")
        return
    projection = data.get("projection", "finite")
    if projection not in ("finite", "divergent"):
        errors.append("Solve Collinear node: projection must be finite or divergent.")
        return
    letter_proj = data.get("letter_projection") or "identity"
    if bootstrap is None:
        errors.append("The bootstrap binary was not found.")
        return
    add_step(
        f"Solve collinear constraints on {target}",
        "bootstrap",
        [bootstrap, "--solve-collinear", "--target", target,
         "--rhs", _resolve_path_arg(proj_dir, rhs),
         "--projection", projection,
         "--letter-projection", _resolve_path_arg(proj_dir, letter_proj),
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [],
        False,
        {"type": "solve_collinear"},
    )
    provides[(node["id"], "solution")] = {"kind": "solution", "file": None, "weight": None, "name": target}


TARGET_RE = re.compile(r"^(SEW_\d+p\d+|FEC_(\d+)|LEC_(\d+))$", re.IGNORECASE)


def _parse_chain_target(target: str) -> tuple[str, int | None, int | None, int | None]:
    m = TARGET_RE.match(target or "")
    if not m:
        return (target, None, None, None)
    kind = target[:3].upper()
    if kind == "SEW":
        f, l = target[4:].split("p")
        return (target, int(f), int(l), int(f) + int(l))
    w = int(target[4:])
    return (target, w, None, w)


def _projection_out_rel(symmetry: str, target: str) -> str:
    """Tensor the projection chain reliably produces for this symmetry/target."""
    _, fw, _, _ = _parse_chain_target(target)
    kind = target[:3].upper()
    if symmetry == "collinear":
        if kind == "SEW":
            return f"output/collinear/{target}_basis.wxf"
        if kind == "FEC":
            return f"output/collinear/first_w{fw}_basis.wxf"
        return f"output/collinear/last_w{fw}_basis.wxf"
    return f"output/{symmetry}/{target}.wxf" if kind == "SEW" else f"output/{symmetry}/first_w{fw}.wxf" if kind == "FEC" else f"output/{symmetry}/last_w{fw}.wxf"


def _compile_projection_chain(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    data = node.get("data") or {}
    symmetry = data.get("symmetry") or ""
    if symmetry not in ("collinear", "cyclic", "flip", "parity"):
        errors.append("A Projection Chain node needs a symmetry (collinear, cyclic, flip or parity).")
        return
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Projection Chain node: provide a target (e.g. SEW_5p1).")
        return
    _, _, _, w = _parse_chain_target(target)
    if w is None:
        errors.append(f"Projection Chain node: invalid target '{target}' (expected SEW_FpL, FEC_W or LEC_W).")
        return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found; build it with `make bootstrap`.")
        return
    sym_dir = _abs(proj_dir, f"output/{symmetry}")
    summary = f"output/{symmetry}/summary.txt"
    add_step(
        f"Projection chain ({symmetry}) on {target}",
        "bootstrap",
        [bootstrap, "--project", "--symmetry", symmetry, "--target", target,
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [summary],
        False,
        {"type": "projection_chain", "symmetry": symmetry, "target": target, "sym_dir": sym_dir},
    )
    rel = _projection_out_rel(symmetry, target)
    provides[(node["id"], "out")] = {"kind": "basis" if symmetry == "collinear" else "tensor", "file": rel, "weight": None, "name": Path(rel).stem}


def _compile_symmetry_invariant(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    data = node.get("data") or {}
    symmetry = data.get("symmetry") or ""
    if symmetry not in ("cyclic", "flip", "parity"):
        errors.append("A Symmetry Invariant node needs a symmetry (cyclic, flip or parity; use Solve Collinear for collinear).")
        return
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Symmetry Invariant node: provide a target (e.g. SEW_5p1).")
        return
    _, _, _, w = _parse_chain_target(target)
    if w is None:
        errors.append(f"Symmetry Invariant node: invalid target '{target}' (expected SEW_FpL, FEC_W or LEC_W).")
        return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found; build it with `make bootstrap`.")
        return
    rel = f"output/{symmetry}/{target}_invariant.wxf"
    add_step(
        f"Symmetry invariant ({symmetry}) of {target}",
        "bootstrap",
        [bootstrap, "--solve-symmetry", "--symmetry", symmetry, "--target", target,
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [rel],
        False,
        {"type": "symmetry_invariant", "symmetry": symmetry, "target": target, "tensor_file": rel},
    )
    provides[(node["id"], "out")] = {"kind": "basis", "file": rel, "weight": None, "name": Path(rel).stem}


def _compile_compute_rhs(node, incoming, provides, add_step, errors, compute_rhs_bin, proj_dir):
    data = node.get("data", {}) or {}
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Compute RHS node: provide a target (e.g. SEW_3p1).")
        return
    letter_proj = data.get("letter_projection") or "identity"
    if compute_rhs_bin is None:
        errors.append("The compute_rhs binary was not found; build it with `make compute_rhs`.")
        return
    add_step(
        f"Compute collinear RHS for {target}",
        "bootstrap",
        [compute_rhs_bin, "--target", target,
         "--letter-projection", _resolve_path_arg(proj_dir, letter_proj),
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [],
        False,
        {"type": "compute_rhs"},
    )
    provides[(node["id"], "boundary")] = {"kind": "boundary", "file": None, "weight": None, "name": target}


def _compile_add_tensors(node, incoming, provides, add_step, errors, tensor_add_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    indexed = []
    seen_handles = set()
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        if th in ("a", "b"):
            th = "in_0" if th == "a" else "in_1"  # legacy two-input layout
        if not th.startswith("in_") or th in seen_handles:
            continue
        seen_handles.add(th)
        try:
            idx = int(th[3:].split("@")[0])
        except ValueError:
            idx = 0
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            indexed.append((idx, provides[src]))
    indexed.sort(key=lambda x: x[0])
    if len(indexed) < 2:
        errors.append("An Add Tensors node needs at least two wired inputs.")
        return
    tensors = [info for _, info in indexed]
    t0 = tensors[0]
    for t in tensors[1:]:
        if t["kind"] != t0["kind"]:
            errors.append(f"Add Tensors requires all inputs to have the same kind, got '{t0['kind']}' and '{t['kind']}'.")
            return
        if t.get("file") is None:
            errors.append("Add Tensors: all inputs must be concrete tensor files (the basis/solution/boundary outputs of Project/Solve/Compute RHS nodes are virtual and cannot be added).")
            return
        if t.get("weight") != t0.get("weight"):
            errors.append("Add Tensors inputs must have the same weight (tensors must have identical dimensions).")
            return

    weights = [w.strip() for w in str(data.get("weights") or "").split(",") if w.strip()]
    if not weights and (data.get("weight_a") or data.get("weight_b")):
        weights = [str(data.get("weight_a") or "1").strip(), str(data.get("weight_b") or "1").strip()]
    weights += ["1"] * (len(tensors) - len(weights))
    weights = weights[: len(tensors)]
    for w in weights:
        if not RAT_RE.match(w):
            errors.append(f"Add Tensors weight '{w}' is not a rational number (examples: 1, -2, 1/2).")
            return

    target = (data.get("target") or "").strip() or nid
    if tensor_add_bin is None:
        errors.append("The tensor_add binary was not found; build it with `make tensor_add`.")
        return

    def rel_path(name):
        return f"{_out()}{name}.wxf"

    def partial_name(k):
        return target if k == len(tensors) - 1 else f"{target}_p{k}"

    for k in range(1, len(tensors)):
        acc_file = tensors[0]["file"] if k == 1 else rel_path(partial_name(k - 1))
        acc_name = tensors[0]["name"] if k == 1 else partial_name(k - 1)
        rel = rel_path(partial_name(k))
        add_step(
            f"Add {acc_name} + {weights[k]}*{tensors[k]['name']} -> {partial_name(k)}",
            "tensor_add",
            [tensor_add_bin, _abs(proj_dir, acc_file), _abs(proj_dir, tensors[k]["file"]), "1" if k > 1 else weights[0], weights[k], _abs(proj_dir, rel)],
            proj_dir,
            [rel],
            True,
            {"type": "add_tensors", "tensor_file": rel},
        )
    provides[(nid, "out")] = {"kind": t0["kind"], "file": rel_path(target), "weight": t0.get("weight"), "name": target, "dims": t0.get("dims")}



SAFE_TARGET_RE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.\-]*$")


def _require_target(data, what, errors, nid=None):
    target = (data.get("target") or "").strip()
    if not target:
        target = nid or None
    if not target:
        errors.append(f"A {what} node needs a target name for the output tensor.")
        return None
    if target and (("/" in target) or (".." in target) or not SAFE_TARGET_RE.match(target)):
        errors.append(f"A {what} node has an invalid target name '{target}' (letters, digits, _, ., - only).")
        return None
    return target


def _compile_ternary_contract(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    m1 = _edge_input(provides, incoming, nid, "trans1")
    m2 = _edge_input(provides, incoming, nid, "trans2")
    if tensor is None or m1 is None or m2 is None:
        missing = [h for h, v in (("tensor", tensor), ("trans1", m1), ("trans2", m2))
                   if v is None and not _has_incoming(incoming, nid, h)]
        if missing:
            errors.append(f"A Ternary Contract node is missing its {' and '.join(missing)} input{'s' if len(missing) > 1 else ''}.")
        return
    if tensor.get("file") is None:
        errors.append("Ternary Contract: the tensor input must be a concrete rank-3 tensor file.")
        return
    if m1.get("file") is None or m2.get("file") is None:
        errors.append("Ternary Contract: trans1/trans2 must be concrete matrix tensors (a matrix property or a Matrix Power output).")
        return
    target = _require_target(data, "Ternary Contract", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Ternary contract {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]), _abs(proj_dir, m1["file"]), _abs(proj_dir, m2["file"]), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "ternary_contract", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"), "name": target}


def _compile_matrix_power(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    mat = _edge_input(provides, incoming, nid, "matrix")
    if mat is None:
        errors.append("A Matrix Power node needs a matrix input.")
        return
    if mat.get("file") is None:
        errors.append("Matrix Power: the matrix input must be a concrete matrix tensor (a symmetry matrix property or another Matrix Power output).")
        return
    try:
        n = int(str(data.get("n") or "").strip())
        if n < 0:
            raise ValueError
    except ValueError:
        errors.append("Matrix Power: n must be a non-negative integer.")
        return
    target = _require_target(data, "Matrix Power", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Matrix power {mat['name']}^{n} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "power", _abs(proj_dir, mat["file"]), str(n), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "matrix_power", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target, "dims": mat.get("dims")}


def _compile_assemble(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    indexed = []
    seen_handles = set()
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        if not th.startswith("in_") or th in seen_handles:
            continue
        seen_handles.add(th)
        try:
            idx = int(th[3:].split("@")[0])
        except ValueError:
            idx = 0
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            indexed.append((idx, provides[src]))
    indexed.sort(key=lambda x: x[0])
    elems = [info for _, info in indexed]
    if not elems:
        errors.append("An Assemble Solution Space node needs at least one element input.")
        return
    kind = elems[0]["kind"]
    for info in elems:
        if info["kind"] != kind:
            errors.append("Assemble Solution Space: all element inputs must have the same kind.")
            return
        if info.get("file") is None:
            errors.append("Assemble Solution Space: element inputs must be concrete tensor files.")
            return

    def _valid_coefs(tokens):
        for t in tokens:
            if not RAT_RE.match(t) or ("/" in t and int(t.split("/")[1]) == 0):
                return t
        return None

    outputs = data.get("outputs")
    if outputs is None:
        # legacy single-output migration: groups + coefs -> per-element coefs
        group_tokens = [t.strip() for t in str(data.get("groups") or "").split(",") if t.strip()]
        coef_tokens = [t.strip() for t in str(data.get("coefs") or "").split(",") if t.strip()]
        if len(group_tokens) == len(elems) and coef_tokens:
            try:
                expanded = [coef_tokens[int(g) - 1] for g in group_tokens]
            except (ValueError, IndexError):
                expanded = None
            if expanded is not None and _valid_coefs(expanded) is None:
                outputs = [{"name": data.get("target") or "out1", "coefs": ",".join(expanded)}]
        if outputs is None:
            outputs = [{"name": (data.get("target") or "").strip() or "out1", "coefs": ",".join(["1"] * len(elems))}]
    if not isinstance(outputs, list) or not outputs:
        errors.append("Assemble Solution Space: add at least one output in the inspector.")
        return
    target = _require_target(data, "Assemble Solution Space", errors, nid)
    if target is None:
        return
    for o in outputs:
        name = str(o.get("name") or "").strip() or target
        if not name:
            errors.append("Assemble Solution Space: every output needs a name.")
            return
        o["name"] = name
        tokens = [t.strip() for t in str(o.get("coefs") or "").split(",") if t.strip()]
        if len(tokens) != len(elems):
            errors.append(
                f"Assemble Solution Space output '{name}': needs one rational coefficient per element ({len(elems)} connected), got {len(tokens)}."
            )
            return
        bad = _valid_coefs(tokens)
        if bad is not None:
            errors.append(f"Assemble Solution Space output '{name}': coefficient '{bad}' is not a rational number (examples: 1, -2, 1/2).")
            return
    names = [str(o.get("name")).strip() for o in outputs]
    if len(set(names)) != len(names):
        errors.append("Assemble Solution Space: output names must be unique.")
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    elem_files = [_abs(proj_dir, info["file"]) for info in elems]
    for oi, o in enumerate(outputs):
        name = names[oi]
        rel = f"{_out()}{target}_{name}.wxf"
        add_step(
            f"Assemble output '{name}' of {target} from {len(elems)} elements",
            "tensor_ops",
            [tensor_ops_bin, "assemble", _abs(proj_dir, rel), "--elems", *elem_files,
             "--coefs", str(o.get("coefs"))],
            proj_dir,
            [rel],
            True,
            {"type": "assemble", "tensor_file": rel},
        )
        provides[(nid, f"out_{oi}")] = {"kind": kind, "file": rel, "weight": None, "name": f"{target}_{name}"}


def _compile_tensor_join(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    a = _edge_input(provides, incoming, nid, "a")
    b = _edge_input(provides, incoming, nid, "b")
    if a is None or b is None:
        errors.append("A Join Tensors node needs both A and B inputs.")
        return
    if a["kind"] != b["kind"]:
        errors.append(f"Join Tensors requires two tensors of the same kind, got '{a['kind']}' and '{b['kind']}'.")
        return
    if a.get("file") is None or b.get("file") is None:
        errors.append("Join Tensors: both inputs must be concrete tensor files (the basis/solution/boundary outputs of Project/Solve/Compute RHS nodes are virtual and cannot be joined).")
        return
    try:
        axis = int(str(data.get("axis") or "").strip())
        if axis == 0:
            raise ValueError
    except ValueError:
        errors.append("Join Tensors: axis must be a nonzero integer (1-based; negative counts from the end).")
        return
    target = _require_target(data, "Join Tensors", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Join {a['name']} + {b['name']} along axis {axis} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "join", _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), str(axis), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "tensor_join", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": a["kind"], "file": rel, "weight": None, "name": target, "dims": _join_dims(a.get("dims"), b.get("dims"), axis)}


def _compile_tensor_dot(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    a = _edge_input(provides, incoming, nid, "a")
    b = _edge_input(provides, incoming, nid, "b")
    if a is None or b is None:
        errors.append("A Tensor Dot node needs both A and B inputs.")
        return
    if a.get("file") is None or b.get("file") is None:
        errors.append("Tensor Dot: both inputs must be concrete tensor files (the basis/solution/boundary outputs of Project/Solve/Compute RHS nodes are virtual and cannot be dotted).")
        return
    def _axis(v, who):
        try:
            ax = int(str(v).strip())
            if ax == 0:
                raise ValueError
            return ax
        except (ValueError, TypeError):
            errors.append(f"Tensor Dot: axis of {who} must be a nonzero integer (1-based; negative counts from the end).")
            return None
    axis_a = _axis(data.get("axis_a", -1), "A")
    axis_b = _axis(data.get("axis_b", -1), "B")
    if axis_a is None or axis_b is None:
        return
    ad, bd = a.get("dims"), b.get("dims")
    if ad and bd:
        ra = axis_a - 1 if axis_a > 0 else len(ad) + axis_a
        rb = axis_b - 1 if axis_b > 0 else len(bd) + axis_b
        if not (0 <= ra < len(ad) and 0 <= rb < len(bd)):
            errors.append(f"Tensor Dot: axis out of range (A is {ad}, B is {bd}).")
            return
        if ad[ra] != bd[rb]:
            errors.append(f"Tensor Dot: A's axis has dimension {ad[ra]} but B's axis has dimension {bd[rb]} — they must match to contract.")
            return
        if len(ad) + len(bd) == 2:
            errors.append("Tensor Dot: this contraction reduces to a scalar (rank-0), which the WXF tensor format cannot store — keep at least one free axis.")
            return
    target = _require_target(data, "Tensor Dot", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Dot {a['name']}[{axis_a}] · {b['name']}[{axis_b}] -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "tdot", _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), str(axis_a), str(axis_b), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "tensor_dot", "tensor_file": rel},
    )
    dims_out = None
    if ad and bd:
        dims_out = [d for i, d in enumerate(ad) if i != ra] + [d for i, d in enumerate(bd) if i != rb]
    provides[(nid, "out")] = {"kind": "tensor", "file": rel, "weight": None, "name": target, "dims": dims_out}


def _check_icond_dims(tensor, dlog, errors, label):
    ts, ds = tensor.get("dims"), dlog.get("dims")
    if not ts or not ds or len(ts) < 2 or len(ds) < 2:
        return True  # dims unknown at compile time; the kernel checks at run time
    if len(ts) < 2 or len(ds) != 3:
        errors.append(f"{label}: tensor input has dims {ts} but needs rank >= 2 (dlogmat is {ds}).")
        return False
    if ts[-1] != ds[0] or ts[-2] != ds[1]:
        errors.append(
            f"{label}: tensor dims {ts} do not match dlogmat {ds} — "
            f"the tensor's last two entries must match the dlogmat's first two ({ds[0]} and {ds[1]})."
        )
        return False
    return True


def _compile_impose(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    dlog = _edge_input(provides, incoming, nid, "dlogmat")
    if tensor is None or dlog is None:
        errors.append("A Solve Integrability node needs tensor and integrability dlog inputs.")
        return
    if tensor.get("file") is None or dlog.get("file") is None:
        errors.append("Solve Integrability: both inputs must be concrete tensor files.")
        return
    if not _check_icond_dims(tensor, dlog, errors, "Solve Integrability"):
        return
    target = _require_target(data, "Solve Integrability", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    trans_flag = "1" if data.get("transpose", True) else "0"
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Solve integrability on {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "impose", _abs(proj_dir, tensor["file"]), _abs(proj_dir, dlog["file"]),
         trans_flag, _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "impose_integrability", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target}


def _compile_integrability_condition(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    dlog = _edge_input(provides, incoming, nid, "dlogmat")
    if tensor is None or dlog is None:
        errors.append("An Integrability Condition node needs tensor and integrability dlog inputs.")
        return
    if tensor.get("file") is None or dlog.get("file") is None:
        errors.append("Integrability Condition: both inputs must be concrete tensor files.")
        return
    if not _check_icond_dims(tensor, dlog, errors, "Integrability Condition"):
        return
    target = _require_target(data, "Integrability Condition", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Integrability conditions on {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "icond", _abs(proj_dir, tensor["file"]), _abs(proj_dir, dlog["file"]),
         _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "integrability_condition", "tensor_file": rel},
    )
    dims_out = None
    if tensor.get("dims") and dlog.get("dims") and len(tensor["dims"]) >= 2 and len(dlog["dims"]) == 3:
        dims_out = [tensor["dims"][0] if len(tensor["dims"]) > 2 else 1]
        inner = 1
        for d in tensor["dims"][1:-2]:
            inner *= d
        dims_out.append(inner * dlog["dims"][2])
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target, "dims": dims_out}


def _compile_solve_conditions(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    cond = _edge_input(provides, incoming, nid, "cond")
    if cond is None:
        errors.append("A Solve Conditions node needs a condition matrix input.")
        return
    if cond.get("file") is None:
        errors.append("Solve Conditions: the input must be a concrete matrix file.")
        return
    target = _require_target(data, "Solve Conditions", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    trans_flag = "1" if data.get("transpose", True) else "0"
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Solve conditions from {cond['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "isolve", _abs(proj_dir, cond["file"]),
         trans_flag, _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "solve_conditions", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target}
