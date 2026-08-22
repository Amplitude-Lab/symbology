from __future__ import annotations

import shutil
import uuid

from . import storage
from .config import REPO_ROOT

TEMPLATES = {
    "4pformfactor": {
        "id": "4pformfactor",
        "name": "4-Point Form Factor (93 letters)",
        "description": "4-point form factor alphabet with 5 independent variables and 5 square roots (arXiv:2212.02410). Ships with precomputed dlogmat, extended-Steinmann dlogmat and symmetry matrices.",
        "n_letters": 93,
        "data_dir": "data_4pformfactor",
        "alphabet": {
            "name": "4pFF",
            "letters": [f"W[{i}]" for i in range(1, 94)],
            "variables": ["u1", "u2", "u3", "v1", "v2"],
            "expressions": [],
            "expr_loader": None,
            "roots": {},
        },
        "data_files": {
            "dlogmat_4pformfactor.wxf": "dlogmat_4pformfactor.wxf",
            "dlogmatES_4pformfactor.wxf": "dlogmatES_4pformfactor.wxf",
            "dlogmat_full_4pformfactor.wxf": "dlogmat_full_4pformfactor.wxf",
            "cycmat.wxf": "cycmat.wxf",
            "flipmat.wxf": "flipmat.wxf",
            "galois1a.wxf": "galois1a.wxf",
            "galois1b.wxf": "galois1b.wxf",
            "galois2a.wxf": "galois2a.wxf",
            "galois2b.wxf": "galois2b.wxf",
            "galois3.wxf": "galois3.wxf",
            "alphabet.wl": "alphabet.wl",
        },
        "properties": [
            {"type": "integrability", "params": {}, "tensor_file": "data/dlogmat_4pformfactor.wxf"},
            {"type": "extended_steinmann", "params": {}, "tensor_file": "data/dlogmatES_4pformfactor.wxf"},
            {"type": "transformation", "params": {"name": "Cyclic", "map": "{u1->u2, u2->u3, u3->1+u2-v1-v2, v1->v2, v2->1+u3-u1-v2}"}, "tensor_file": "data/cycmat.wxf"},
            {"type": "transformation", "params": {"name": "Flip", "map": "{u1->1+u2-v1-v2, u2->u3, u3->u2, v1->1+u3-u1-v2, v2->v2}"}, "tensor_file": "data/flipmat.wxf"},
            {"type": "transformation", "params": {"name": "Galois1a", "map": "{sqrt1a->-sqrt1a}"}, "tensor_file": "data/galois1a.wxf"},
            {"type": "transformation", "params": {"name": "Galois1b", "map": "{sqrt1b->-sqrt1b}"}, "tensor_file": "data/galois1b.wxf"},
            {"type": "transformation", "params": {"name": "Galois2a", "map": "{sqrt2a->-sqrt2a}"}, "tensor_file": "data/galois2a.wxf"},
            {"type": "transformation", "params": {"name": "Galois2b", "map": "{sqrt2b->-sqrt2b}"}, "tensor_file": "data/galois2b.wxf"},
            {"type": "transformation", "params": {"name": "Galois3", "map": "{sqrt3->-sqrt3}"}, "tensor_file": "data/galois3.wxf"},
        ],
    },
    "pentagon": {
        "id": "pentagon",
        "name": "Pentagon (31 letters)",
        "description": "Pentagon alphabet with one square root eps5 = Sqrt[Delta5]. Ships with precomputed dlogmat and cyclic/flip matrices.",
        "n_letters": 31,
        "data_dir": "data_pentagon",
        "alphabet": {
            "name": "pentagon",
            "letters": [f"W[{i}]" for i in range(1, 32)],
            "variables": ["s12", "s23", "s34", "s45", "s15"],
            "expressions": [],
            "expr_loader": None,
            "roots": {"eps5": "Sqrt[Delta5]"},
        },
        "data_files": {
            "dlogmat_pentagon.wxf": "dlogmat_pentagon.wxf",
            "cycmat.wxf": "cycmat.wxf",
            "flipmat.wxf": "flipmat.wxf",
            "alphabet.wl": "alphabet.wl",
        },
        "properties": [
            {"type": "integrability", "params": {}, "tensor_file": "data/dlogmat_pentagon.wxf"},
            {"type": "transformation", "params": {"name": "Cyclic", "map": "{s12->s23, s23->s34, s34->s45, s45->s15, s15->s12}"}, "tensor_file": "data/cycmat.wxf"},
            {"type": "transformation", "params": {"name": "Flip", "map": "{s12->s45, s23->s15, s34->s34, s45->s12, s15->s23}"}, "tensor_file": "data/flipmat.wxf"},
        ],
    },
    "e6": {
        "id": "e6",
        "name": "E6 Hexagon Wilson Loop (42 letters)",
        "description": "E6 adjoint hexagon Wilson loop seed tensors (11-letter collinear space). Ships with dlogmat, FEC_1/LEC_1 seeds, collinear and symmetry seeds, and the one-loop E1 tensor.",
        "n_letters": 42,
        "data_dir": "data",
        "alphabet": {
            "name": "E6",
            "letters": [f"W[{i}]" for i in range(1, 43)],
            "variables": [],
            "expressions": [],
            "expr_loader": None,
            "roots": {},
        },
        "data_files": {
            "dlogmat_E6.wxf": "dlogmat_E6.wxf",
            "FEC_1.wxf": "FEC_1.wxf",
            "LEC_1.wxf": "LEC_1.wxf",
            "colmat42.wxf": "colmat42.wxf",
            "colprojdiv.wxf": "colprojdiv.wxf",
            "colprojfin.wxf": "colprojfin.wxf",
            "cycrepmat.wxf": "cycrepmat.wxf",
            "fliprepmat.wxf": "fliprepmat.wxf",
            "parityrepmat.wxf": "parityrepmat.wxf",
            "E1.wxf": "E1.wxf",
        },
        "properties": [
            {"type": "integrability", "params": {}, "tensor_file": "data/dlogmat_E6.wxf"},
            {"type": "first_entry", "params": {}, "tensor_file": "data/FEC_1.wxf"},
            {"type": "last_entry", "params": {}, "tensor_file": "data/LEC_1.wxf"},
            {"type": "precomputed_tensor", "name": "colmat42 (collinear map)", "params": {"tensor_file": "data/colmat42.wxf"}, "tensor_file": "data/colmat42.wxf"},
            {"type": "precomputed_tensor", "name": "colprojdiv (divergent)", "params": {"tensor_file": "data/colprojdiv.wxf"}, "tensor_file": "data/colprojdiv.wxf"},
            {"type": "precomputed_tensor", "name": "colprojfin (finite)", "params": {"tensor_file": "data/colprojfin.wxf"}, "tensor_file": "data/colprojfin.wxf"},
            {"type": "precomputed_tensor", "name": "cycrepmat (cyclic rep)", "params": {"tensor_file": "data/cycrepmat.wxf"}, "tensor_file": "data/cycrepmat.wxf"},
            {"type": "precomputed_tensor", "name": "fliprepmat (flip rep)", "params": {"tensor_file": "data/fliprepmat.wxf"}, "tensor_file": "data/fliprepmat.wxf"},
            {"type": "precomputed_tensor", "name": "parityrepmat (parity rep)", "params": {"tensor_file": "data/parityrepmat.wxf"}, "tensor_file": "data/parityrepmat.wxf"},
            {"type": "precomputed_tensor", "name": "E1 (one-loop seed)", "params": {"tensor_file": "data/E1.wxf"}, "tensor_file": "data/E1.wxf"},
        ],
    },
}


def _expr_loader(template_id: str, data_dir_abs: Path) -> str | None:
    if template_id == "4pformfactor":
        return f'Get["{(data_dir_abs / "alphabet.wl").as_posix()}"];\nalphaExpr = alphabetf[[All, 2]] /. sqrtrep;'
    if template_id == "pentagon":
        return f'Get["{(data_dir_abs / "alphabet.wl").as_posix()}"];\nalphaExpr = LetterRep[[All, 2]] /. RootDef;'
    return None


def list_templates() -> list:
    return [
        {"id": t["id"], "name": t["name"], "description": t["description"], "n_letters": t["n_letters"]}
        for t in TEMPLATES.values()
    ]


def _flow(name, nodes, edges):
    return {"id": uuid.uuid4().hex[:8], "name": name, "graph": {"nodes": nodes, "edges": edges}}


def _n(nid, ntype, x, y, **data):
    return {"id": nid, "type": ntype, "position": {"x": x, "y": y}, "data": dict(data)}


def _e(s, sh, t, th):
    return {"id": f"{s}.{sh}->{t}.{th}", "source": s, "sourceHandle": sh, "target": t, "targetHandle": th}


def _sew_chain_nodes(prefix, pmap, fec_to=5, lec_to=4, sews=(("2p2", 2, 2, None),), x0=0, y0=0):
    """Alphabet + FEC/LEC extend chains + sew nodes, replicating run_workflow.sh (smoke).

    sews entries: (name, F, L, lec_from) where lec_from None means LEC_{L}
    from the extend chain, an int means the data/ seed LEC_1.
    """
    nodes = [
        _n("alpha", "alphabet", x0, y0 + 240,
           alphabet_id=pmap["__alphabet_id__"],
           selected_properties=[pmap["integrability"], pmap["first_entry"], pmap["last_entry"]]),
    ]
    edges = []
    prev_fec_handle = ("alpha", f"prop_{pmap['first_entry']}")
    x = x0 + 300
    for w in range(2, fec_to + 1):
        nid = f"ef{w}"
        nodes.append(_n(nid, "extend", x, y0 + (w - 2) * 130, target_weight=w))
        edges.append(_e(prev_fec_handle[0], prev_fec_handle[1], nid, "fec"))
        edges.append(_e("alpha", f"prop_{pmap['integrability']}", nid, "condition"))
        prev_fec_handle = (nid, "fec")
        x += 300
    prev_lec_handle = ("alpha", f"prop_{pmap['last_entry']}")
    for w in range(2, lec_to + 1):
        nid = f"el{w}"
        nodes.append(_n(nid, "extend", x0 + 300, y0 + 700 + (w - 2) * 130, target_weight=w))
        edges.append(_e(prev_lec_handle[0], prev_lec_handle[1], nid, "lec"))
        edges.append(_e("alpha", f"prop_{pmap['integrability']}", nid, "condition"))
        prev_lec_handle = (nid, "lec")
    lec_handles = {1: ("alpha", f"prop_{pmap['last_entry']}")}
    for w in range(2, lec_to + 1):
        lec_handles[w] = (f"el{w}", "lec")
    fec_handles = {w: (f"ef{w}", "fec") for w in range(2, fec_to + 1)}
    x = x0 + 300 + fec_to * 300
    for i, (name, f, l, lec_from) in enumerate(sews):
        nid = f"sew_{name}"
        nodes.append(_n(nid, "sew", x, y0 + i * 130, target=f"SEW_{name}"))
        edges.append(_e("alpha", f"prop_{pmap['integrability']}", nid, "condition"))
        edges.append(_e(*fec_handles[f], nid, "fec"))
        src = lec_handles[l if lec_from is None else lec_from]
        edges.append(_e(*src, nid, "lec"))
    return nodes, edges


def _e6_flows(pmap):
    flows = []
    nodes, edges = _sew_chain_nodes(
        "e6", pmap, fec_to=5, lec_to=4,
        sews=(("2p2", 2, 2, None), ("3p1", 3, 1, 1), ("4p2", 4, 2, None), ("5p1", 5, 1, 1)),
    )
    flows.append(_flow("E6 · extend & sew chain (smoke)", nodes, edges))

    nodes, edges = _sew_chain_nodes("e6", pmap, fec_to=5, lec_to=1, sews=(("5p1", 5, 1, 1),))
    base = list(nodes)
    base_e = list(edges)
    syms = [("pc_col", "projection_chain", "collinear"), ("pc_cyc", "projection_chain", "cyclic"),
            ("pc_flip", "projection_chain", "flip"), ("pc_par", "projection_chain", "parity"),
            ("si_cyc", "symmetry_invariant", "cyclic"), ("si_flip", "symmetry_invariant", "flip"),
            ("si_par", "symmetry_invariant", "parity")]
    for i, (nid, ntype, sym) in enumerate(syms):
        base.append(_n(nid, ntype, 2300, 200 + i * 150, symmetry=sym, target="SEW_5p1"))
        base_e.append(_e("sew_5p1", "sew", nid, "seed"))
    flows.append(_flow("E6 · symmetry projection & invariants", base, base_e))

    nodes, edges = _sew_chain_nodes("e6", pmap, fec_to=3, lec_to=1, sews=(("3p1", 3, 1, 1),))
    nodes.append(_n("pc_col3", "projection_chain", 1800, 200, symmetry="collinear", target="SEW_3p1"))
    edges.append(_e("sew_3p1", "sew", "pc_col3", "seed"))
    nodes.append(_n("rhs3", "compute_rhs", 2150, 380, target="SEW_3p1", letter_projection="data/colprojdiv.wxf"))
    edges.append(_e("sew_3p1", "sew", "rhs3", "seed"))
    nodes.append(_n("sc3", "solve_collinear", 2500, 560, target="SEW_3p1",
                    rhs="output/2loop/boundary_2L.wxf", projection="divergent",
                    letter_projection="data/colprojdiv.wxf"))
    edges.append(_e("sew_3p1", "sew", "sc3", "seed"))
    flows.append(_flow("E6 · collinear bootstrap & 2-loop RHS", nodes, edges))
    return flows


def _pentagon_flows(pmap):
    nodes = [
        _n("alpha", "alphabet", 0, 400, alphabet_id=pmap["__alphabet_id__"],
           selected_properties=[pmap["integrability"], pmap["Cyclic"], pmap["Flip"]]),
    ]
    edges = []
    for w, src in enumerate(("cyc2", "cyc3", "cyc4")):
        nodes.append(_n(src, "matrix_power", 300, 60 + w * 140, n=w + 2, target=f"cyc{w + 2}"))
        edges.append(_e("alpha", f"prop_{pmap['Cyclic']}", src, "matrix"))
    terms = []
    x = 620
    i = 0
    for cyc_src in ("alpha", "cyc2", "cyc3", "cyc4"):
        nid = f"t_c{i + 1}"
        nodes.append(_n(nid, "ternary_contract", x, 60 + i * 130, target=f"P5_cyc{i + 1}"))
        edges.append(_e("alpha", f"prop_{pmap['integrability']}", nid, "tensor"))
        cyc_handle = "prop_" + pmap["Cyclic"] if cyc_src == "alpha" else "out"
        edges.append(_e(cyc_src, cyc_handle, nid, "trans1"))
        edges.append(_e(cyc_src, cyc_handle, nid, "trans2"))
        terms.append(nid)
        i += 1
    for k, cyc_src in enumerate(("alpha", "cyc2", "cyc3", "cyc4")):
        nid = f"t_fc{k + 1}"
        nodes.append(_n(nid, "ternary_contract", x, 60 + i * 130, target=f"P5_fc{k + 1}"))
        edges.append(_e("alpha", f"prop_{pmap['integrability']}", nid, "tensor"))
        edges.append(_e("alpha", f"prop_{pmap['Flip']}", nid, "trans1"))
        cyc_handle = "prop_" + pmap["Cyclic"] if cyc_src == "alpha" else "out"
        edges.append(_e(cyc_src, cyc_handle, nid, "trans2"))
        terms.append(nid)
        i += 1
    nodes.append(_n("avg", "add_tensors", 1100, 500,
                    weights=",".join(["1/10"] * 10), target="P5_dlogmat_inv"))
    for k, nid in enumerate(terms):
        edges.append(_e(nid, "out", "avg", f"in_{k}"))
    return [_flow("Pentagon · D5-invariant integrability tensor", nodes, edges)]


def _4pff_flows(pmap):
    nodes = [
        _n("alpha", "alphabet", 0, 300, alphabet_id=pmap["__alphabet_id__"],
           selected_properties=[pmap["integrability"], pmap["extended_steinmann"],
                                 pmap["Cyclic"], pmap["Flip"]]),
        _n("merge", "merge_conditions", 320, 300),
        _n("ss_cyc", "solve_symmetry", 620, 200, target="S_cyc"),
        _n("ss_flip", "solve_symmetry", 920, 200, target="S_cycflip"),
    ]
    edges = [
        _e("alpha", f"prop_{pmap['integrability']}", "merge", "in_0"),
        _e("alpha", f"prop_{pmap['extended_steinmann']}", "merge", "in_1"),
        _e("merge", "out", "ss_cyc", "tensor"),
        _e("alpha", f"prop_{pmap['Cyclic']}", "ss_cyc", "matrix"),
        _e("ss_cyc", "out", "ss_flip", "tensor"),
        _e("alpha", f"prop_{pmap['Flip']}", "ss_flip", "matrix"),
    ]
    return [_flow("4pFF · full conditions + cyclic & flip invariants", nodes, edges)]


def _flows_for(template_id: str, pmap: dict) -> list:
    if template_id == "e6":
        return _e6_flows(pmap)
    if template_id == "pentagon":
        return _pentagon_flows(pmap)
    if template_id == "4pformfactor":
        return _4pff_flows(pmap)
    return []


def apply_template(proj: dict, template_id: str) -> None:
    t = TEMPLATES[template_id]
    src_dir = REPO_ROOT / t["data_dir"]
    data_dir = storage.project_dir(proj["id"]) / "data"
    for src_name, dst_name in t["data_files"].items():
        src = src_dir / src_name
        if src.exists():
            shutil.copy2(src, data_dir / dst_name)
    alpha = dict(t["alphabet"])
    alpha["id"] = uuid.uuid4().hex[:8]
    loader = _expr_loader(template_id, data_dir)
    if loader:
        alpha["expr_loader"] = loader
    props = []
    pmap = {"__alphabet_id__": alpha["id"]}
    for p in t["properties"]:
        pid = uuid.uuid4().hex[:8]
        pmap[p.get("name") or p["params"].get("name") or p["type"]] = pid
        props.append({
            "id": pid,
            "type": p["type"],
            "name": p.get("name", ""),
            "params": p["params"],
            "status": "ready",
            "error": None,
            "precomputed": True,
            "tensor_file": p["tensor_file"],
            "summary": None,
        })
    alpha["properties"] = props
    proj["alphabets"].append(alpha)
    proj["flows"].extend(_flows_for(template_id, pmap))
