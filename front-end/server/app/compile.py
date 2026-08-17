from __future__ import annotations

import re
import shlex
from pathlib import Path

from . import storage
from .config import find_bootstrap, find_compute_rhs, find_tensor_add, find_tensor_ops, find_wolframscript
from .wolfram import PROP_KIND, merge_script, property_display_name, property_script, property_tensor_relpath

TENSOR_KINDS = {"dlogmat", "fec1", "fec", "lec1", "lec", "sew", "matrix", "basis", "solution", "boundary"}

EDGE_RULES = {
    ("merge_conditions", None): {"dlogmat"},
    ("extend", "condition"): {"dlogmat"},
    ("extend", "fec"): {"fec1", "fec"},
    ("extend", "lec"): {"lec1", "lec"},
    ("sew", "condition"): {"dlogmat"},
    ("sew", "fec"): {"fec1", "fec"},
    ("sew", "lec"): {"lec1", "lec"},
    ("project", "seed"): {"fec1", "fec", "sew"},
    ("solve_symmetry", "seed"): {"fec1", "fec", "sew"},
    ("solve_collinear", "seed"): {"fec1", "fec", "sew"},
    ("add_tensors", "a"): TENSOR_KINDS,
    ("add_tensors", "b"): TENSOR_KINDS,
    ("apply_symmetry", "tensor"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("apply_symmetry", "trans1"): {"matrix"},
    ("apply_symmetry", "trans2"): {"matrix"},
    ("matrix_power", "matrix"): {"matrix"},
    ("tensor_join", "a"): TENSOR_KINDS,
    ("tensor_join", "b"): TENSOR_KINDS,
}

SYMMETRIES = {"collinear", "cyclic", "flip", "parity"}

RAT_RE = re.compile(r"^-?\d+(/\d+)?$")


def _abs(proj_dir: Path, rel: str) -> str:
    p = Path(rel)
    return str(p if p.is_absolute() else proj_dir / rel)


def _resolve_path_arg(proj_dir: Path, value: str) -> str:
    if value in ("0", "identity"):
        return value
    return _abs(proj_dir, value)


def compile_flow(proj: dict, graph: dict) -> dict:
    errors: list[str] = []
    steps: list[dict] = []
    nodes = {n["id"]: n for n in graph.get("nodes", [])}
    edges = graph.get("edges", [])
    proj_dir = storage.project_dir(proj["id"])
    gen_dir = proj_dir / "wolfram_gen"
    gen_dir.mkdir(exist_ok=True)

    wolframscript = find_wolframscript()
    bootstrap = find_bootstrap()
    compute_rhs_bin = find_compute_rhs()
    tensor_add_bin = find_tensor_add()
    tensor_ops_bin = find_tensor_ops()

    if not nodes:
        return {"ok": False, "errors": ["The flow is empty. Add at least one node."], "steps": []}

    incoming: dict[str, list] = {nid: [] for nid in nodes}
    outgoing: dict[str, list] = {nid: [] for nid in nodes}
    for e in edges:
        if e.get("source") in nodes and e.get("target") in nodes:
            incoming[e["target"]].append(e)
            outgoing[e["source"]].append(e)

    order = _topo_sort(nodes, incoming, outgoing)
    if order is None:
        return {"ok": False, "errors": ["The flow contains a cycle. Break the loop and try again."], "steps": []}

    provides: dict[tuple, dict] = {}
    step_counter = [0]

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

    for nid in order:
        node = nodes[nid]
        ntype = node.get("type")
        data = node.get("data", {}) or {}

        if ntype == "alphabet":
            _compile_alphabet(proj, node, provides, add_step, errors, wolframscript, gen_dir, proj_dir)
        elif ntype == "merge_conditions":
            _compile_merge(node, incoming, provides, add_step, errors, wolframscript, gen_dir, proj_dir)
        elif ntype == "extend":
            _compile_extend(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "sew":
            _compile_sew(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "project":
            _compile_project(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "solve_symmetry":
            _compile_solve_symmetry(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "solve_collinear":
            _compile_solve_collinear(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "compute_rhs":
            _compile_compute_rhs(node, incoming, provides, add_step, errors, compute_rhs_bin, proj_dir)
        elif ntype == "add_tensors":
            _compile_add_tensors(node, incoming, provides, add_step, errors, tensor_add_bin, proj_dir)
        elif ntype == "apply_symmetry":
            _compile_apply_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "matrix_power":
            _compile_matrix_power(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "tensor_join":
            _compile_tensor_join(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        else:
            errors.append(f"Unknown node type '{ntype}'.")

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
        key = (ttype, th.split("@")[0] if th.startswith("in_") else th)
        if ttype == "merge_conditions":
            allowed = EDGE_RULES[("merge_conditions", None)]
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
    return {"ok": True, "errors": [], "steps": [{k: v for k, v in s.items() if k != "argv"} for s in steps], "_steps_full": steps}


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


def _compile_alphabet(proj, node, provides, add_step, errors, wolframscript, gen_dir, proj_dir):
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
        rel = prop.get("tensor_file") or property_tensor_relpath(alphabet, prop)
        if prop.get("status") != "ready" or not (proj_dir / rel).exists():
            if wolframscript is None:
                errors.append("wolframscript was not found; cannot compute properties.")
                continue
            out_abs = _abs(proj_dir, rel)
            script = property_script(alphabet, prop, out_abs)
            script_path = gen_dir / f"prop_{prop_id}.wl"
            script_path.write_text(script)
            add_step(
                f"Compute {property_display_name(prop)} ({prop['type'].replace('_', ' ')}) for alphabet '{alphabet['name']}'",
                "wolfram",
                [wolframscript, "-script", str(script_path)],
                proj_dir,
                [rel],
                True,
                {"type": "property", "alphabet_id": alphabet["id"], "property_id": prop_id, "tensor_file": rel},
            )
        provides[(nid, f"prop_{prop_id}")] = {
            "kind": kind,
            "file": rel,
            "weight": 1 if kind in ("fec1", "lec1") else None,
            "name": Path(rel).stem,
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
    script_path.write_text(script)
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
    if target not in (None, "", 0) and int(target) != w_out:
        errors.append(f"Extend node target weight {target} does not match input weight {w_in} + 1.")
        return
    rel = f"output/{direction}_{w_out}.wxf"
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


def _compile_sew(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
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
    name = f"SEW_{fw}p{lw}"
    rel = f"output/{name}.wxf"
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


def _compile_project(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    data = node.get("data", {}) or {}
    symmetry = data.get("symmetry", "")
    if symmetry not in SYMMETRIES:
        errors.append(f"Project node: symmetry must be one of {sorted(SYMMETRIES)}.")
        return
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Project node: provide a target (e.g. SEW_5p1 or FEC_3) or connect a seed input.")
        return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found.")
        return
    add_step(
        f"Project {target} ({symmetry})",
        "bootstrap",
        [bootstrap, "--project", "--symmetry", symmetry, "--target", target,
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [],
        False,
        {"type": "project"},
    )
    provides[(node["id"], "basis")] = {"kind": "basis", "file": None, "weight": None, "name": target}


def _compile_solve_symmetry(node, incoming, provides, add_step, errors, bootstrap, proj_dir):
    data = node.get("data", {}) or {}
    symmetry = data.get("symmetry", "")
    if symmetry not in SYMMETRIES - {"collinear"}:
        errors.append("Solve Symmetry node: symmetry must be cyclic, flip or parity.")
        return
    target = _derive_target(node, incoming, provides, data)
    if not target:
        errors.append("Solve Symmetry node: provide a target or connect a seed input.")
        return
    if bootstrap is None:
        errors.append("The bootstrap binary was not found.")
        return
    add_step(
        f"Solve {symmetry} symmetry on {target}",
        "bootstrap",
        [bootstrap, "--solve-symmetry", "--symmetry", symmetry, "--target", target,
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [],
        False,
        {"type": "solve_symmetry"},
    )
    provides[(node["id"], "solution")] = {"kind": "solution", "file": None, "weight": None, "name": target}


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
    a = _edge_input(provides, incoming, nid, "a")
    b = _edge_input(provides, incoming, nid, "b")
    if a is None or b is None:
        errors.append("An Add Tensors node needs both A and B inputs.")
        return
    if a["kind"] != b["kind"]:
        errors.append(f"Add Tensors requires two tensors of the same kind, got '{a['kind']}' and '{b['kind']}'.")
        return
    if a.get("weight") != b.get("weight"):
        errors.append("Add Tensors inputs must have the same weight (tensors must have identical dimensions).")
        return
    wa = str(data.get("weight_a") or "1").strip()
    wb = str(data.get("weight_b") or "1").strip()
    for w in (wa, wb):
        if not RAT_RE.match(w):
            errors.append(f"Add Tensors weight '{w}' is not a rational number (examples: 1, -2, 1/2).")
            return
    target = (data.get("target") or "").strip()
    if not target:
        errors.append("An Add Tensors node needs a target name for the output tensor.")
        return
    if tensor_add_bin is None:
        errors.append("The tensor_add binary was not found; build it with `make tensor_add`.")
        return
    rel = f"output/{target}.wxf"
    add_step(
        f"Add {wa}*{a['name']} + {wb}*{b['name']} -> {target}",
        "tensor_add",
        [tensor_add_bin, _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), wa, wb, _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "add_tensors", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": a["kind"], "file": rel, "weight": a.get("weight"), "name": target}


def _require_target(data, what, errors):
    target = (data.get("target") or "").strip()
    if not target:
        errors.append(f"A {what} node needs a target name for the output tensor.")
        return None
    return target


def _compile_apply_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    m1 = _edge_input(provides, incoming, nid, "trans1")
    m2 = _edge_input(provides, incoming, nid, "trans2")
    if tensor is None or m1 is None or m2 is None:
        errors.append("An Apply Symmetry node needs tensor, trans1 and trans2 inputs.")
        return
    if m1.get("file") is None or m2.get("file") is None:
        errors.append("Apply Symmetry: trans1/trans2 must be concrete matrix tensors (a symmetry matrix property or a Matrix Power output).")
        return
    target = _require_target(data, "Apply Symmetry", errors)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"output/{target}.wxf"
    add_step(
        f"Apply symmetry to {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]), _abs(proj_dir, m1["file"]), _abs(proj_dir, m2["file"]), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "apply_symmetry", "tensor_file": rel},
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
    target = _require_target(data, "Matrix Power", errors)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"output/{target}.wxf"
    add_step(
        f"Matrix power {mat['name']}^{n} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "power", _abs(proj_dir, mat["file"]), str(n), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "matrix_power", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target}


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
    try:
        axis = int(str(data.get("axis") or "").strip())
        if axis == 0:
            raise ValueError
    except ValueError:
        errors.append("Join Tensors: axis must be a nonzero integer (1-based; negative counts from the end).")
        return
    target = _require_target(data, "Join Tensors", errors)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"output/{target}.wxf"
    add_step(
        f"Join {a['name']} + {b['name']} along axis {axis} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "join", _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), str(axis), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "tensor_join", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": a["kind"], "file": rel, "weight": None, "name": target}
