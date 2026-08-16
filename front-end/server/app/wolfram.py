from __future__ import annotations

import json
from pathlib import Path

from .config import REPO_ROOT

PACKAGE_PATH = (REPO_ROOT / "SymbolBootstrap.wl").as_posix()

DLOGMAT_TYPES = {"integrability", "extended_steinmann", "cluster_adjacency"}

PROP_KIND = {
    "integrability": "dlogmat",
    "extended_steinmann": "dlogmat",
    "cluster_adjacency": "dlogmat",
    "first_entry": "fec1",
    "last_entry": "lec1",
    "transformation": "matrix",
    "precomputed_tensor": "matrix",
}


def _q(s: str) -> str:
    return json.dumps(s)


def _wl_string_list(items: list) -> str:
    return "{" + ",".join(_q(str(i)) for i in items) + "}"


def _safe_name(name: str) -> str:
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in name)


def property_display_name(prop: dict) -> str:
    name = (prop.get("name") or "").strip()
    if name:
        return name
    if prop["type"] == "transformation":
        return prop.get("params", {}).get("name", "transformation")
    return prop["type"].replace("_", " ")


def property_tensor_relpath(alphabet: dict, prop: dict) -> str:
    aname = alphabet["name"]
    ptype = prop["type"]
    suffix = _safe_name(prop["name"].strip()) if (prop.get("name") or "").strip() else ""
    if ptype == "integrability":
        return f"data/dlogmat_{aname}{'_' + suffix if suffix else ''}.wxf"
    if ptype == "extended_steinmann":
        return f"data/dlogmatES_{aname}{'_' + suffix if suffix else ''}.wxf"
    if ptype == "cluster_adjacency":
        return f"data/dlogmatCA_{aname}{'_' + suffix if suffix else ''}.wxf"
    if ptype == "first_entry":
        return f"data/FEC_1{'_' + suffix if suffix else ''}.wxf"
    if ptype == "last_entry":
        return f"data/LEC_1{'_' + suffix if suffix else ''}.wxf"
    if ptype == "transformation":
        tname = prop.get("params", {}).get("name", "trans")
        return f"data/{_safe_name(tname)}.wxf"
    if ptype == "precomputed_tensor":
        rel = (prop.get("params", {}) or {}).get("tensor_file", "").strip()
        if rel:
            return rel
        return f"data/{_safe_name((prop.get('name') or 'tensor').strip())}.wxf"
    raise ValueError(f"unknown property type {ptype}")


def _preamble(alphabet: dict) -> str:
    lines = [
        f'Get[{_q(PACKAGE_PATH)}];',
        f'$alphaName = {_q(alphabet["name"])};',
        f'$letters = ToExpression[{_wl_string_list(alphabet["letters"])}];',
        'SymbolBootstrap`ResetAlphabet[$alphaName, $letters];',
    ]
    return "\n".join(lines)


def _expression_setup(alphabet: dict) -> str:
    loader = alphabet.get("expr_loader") or ""
    if loader.strip():
        return loader.rstrip().rstrip(";") + ";\nSymbolBootstrap`SetAlphabetExpression[$alphaName, alphaExpr];"
    exprs = alphabet.get("expressions") or []
    return f'SymbolBootstrap`SetAlphabetExpression[$alphaName, ToExpression[{_wl_string_list(exprs)}]];'


def _finish(out_abs: str) -> str:
    result_abs = out_abs + ".result.json"
    return f'''
If[FailureQ[$tensor] || !ArrayQ[$tensor],
  Print["@@RESULT@@FAIL"]; Exit[1]];
Export[{_q(out_abs)}, $tensor];
Export[{_q(result_abs)}, ExportString[<|"ok" -> True, "dims" -> Dimensions[$tensor], "nnz" -> Length[$tensor["NonzeroPositions"]]|>, "JSON"], "String"];
Print["@@RESULT@@OK"];
'''


def property_script(alphabet: dict, prop: dict, out_abs: str) -> str:
    ptype = prop["type"]
    params = prop.get("params", {}) or {}
    if ptype == "precomputed_tensor":
        raise ValueError("precomputed tensors are shipped files and cannot be computed")
    parts = [_preamble(alphabet)]

    if ptype in ("integrability", "transformation"):
        parts.append(_expression_setup(alphabet))

    if ptype == "integrability":
        parts.append('$tensor = SymbolBootstrap`GetIntegrabilityTensor[$alphaName];')
    elif ptype == "first_entry":
        parts.append(f'SymbolBootstrap`SetFirstEntry[$alphaName, ToExpression[{_wl_string_list(params.get("letters", []))}]];')
        parts.append('$tensor = SymbolBootstrap`GetFirstEntryTensor[$alphaName];')
    elif ptype == "last_entry":
        parts.append(f'SymbolBootstrap`SetLastEntry[$alphaName, ToExpression[{_wl_string_list(params.get("letters", []))}]];')
        parts.append('$tensor = SymbolBootstrap`GetLastEntryTensor[$alphaName];')
    elif ptype == "extended_steinmann":
        pairs = params.get("nonadjacent_pairs", [])
        pair_str = "{" + ",".join("{" + str(a) + "," + str(b) + "}" for a, b in pairs) + "}"
        parts.append(f'SymbolBootstrap`SetExtendedSteinmann[$alphaName, ToExpression[{_q(pair_str)}]];')
        parts.append('$tensor = SymbolBootstrap`GetExtendedSteinmannTensor[$alphaName];')
    elif ptype == "cluster_adjacency":
        pairs = params.get("adjacent_pairs", [])
        pair_str = "{" + ",".join("{" + str(a) + "," + str(b) + "}" for a, b in pairs) + "}"
        parts.append(f'SymbolBootstrap`SetClusterAdjacency[$alphaName, ToExpression[{_q(pair_str)}]];')
        parts.append('$tensor = SymbolBootstrap`GetClusterAdjacencyTensor[$alphaName];')
    elif ptype == "transformation":
        tname = params.get("name", "trans")
        kmap = params.get("map", "{}")
        parts.append(f'SymbolBootstrap`SetLetterTransformation[$alphaName, {_q(tname)}, ToExpression[{_q(kmap)}]];')
        parts.append(f'$tensor = SymbolBootstrap`GetLetterTransformationTensor[$alphaName, {_q(tname)}];')
    else:
        raise ValueError(f"unknown property type {ptype}")

    parts.append(_finish(out_abs))
    return "\n".join(parts) + "\n"


def merge_script(input_abs_paths: list, out_abs: str) -> str:
    files = "{" + ",".join(_q(p) for p in input_abs_paths) + "}"
    return f'''Get[{_q(PACKAGE_PATH)}];
$inputs = Import /@ {files};
If[MemberQ[$inputs, $Failed], Print["@@RESULT@@FAIL"]; Exit[1]];
$tensor = SymbolBootstrap`Private`CombineConditionTensor[Sequence @@ $inputs];
{_finish(out_abs)}'''


def summary_script(file_abs: str, result_abs: str) -> str:
    return f'''$t = Import[{_q(file_abs)}];
If[FailureQ[$t], Print["@@RESULT@@FAIL"]; Exit[1]];
$pos = $t["NonzeroPositions"];
$sample = Take[$pos, UpTo[20]];
$vals = Extract[$t, $sample];
Export[{_q(result_abs)}, ExportString[<|"ok" -> True, "dims" -> Dimensions[$t], "nnz" -> Length[$pos],
  "sample" -> MapThread[<|"index" -> #1, "value" -> ToString[#2, InputForm]|> &, {{$sample, $vals}}]|>, "JSON"], "String"];
Print["@@RESULT@@OK"];
'''


def parse_result_marker(line: str) -> dict | None:
    if "@@RESULT@@OK" in line:
        return {"ok": True}
    if "@@RESULT@@FAIL" in line:
        return {"ok": False}
    return None


def read_result_file(tensor_abs: str) -> dict | None:
    path = tensor_abs + ".result.json"
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None
