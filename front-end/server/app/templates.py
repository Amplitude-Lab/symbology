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
    for p in t["properties"]:
        props.append({
            "id": uuid.uuid4().hex[:8],
            "type": p["type"],
            "params": p["params"],
            "status": "ready",
            "error": None,
            "precomputed": True,
            "tensor_file": p["tensor_file"],
            "summary": None,
        })
    alpha["properties"] = props
    proj["alphabets"].append(alpha)
