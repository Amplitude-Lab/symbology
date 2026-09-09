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
    "letter_symmetry": "matrix",
    "sparse_expression": "matrix",
}


def _q(s: str) -> str:
    return json.dumps(s)


EXP2SA_DEFS = """\
ClearAll[aCollect];
aCollect[exp_, name_] := Module[{nl},
   If[exp === 0, Return[{}]];
   nl = Cases[{exp // Expand}, _name, Infinity] // DeleteDuplicates;
   Return[{nl, Normal[CoefficientArrays[exp, nl]][[2]] // Expand}];
   ];
ClearAll[Exp2SA];
Options[Exp2SA] = {"dim" -> 0};
Exp2SA[exp_, OptionsPattern[]] := Module[{tem, result},
   tem = exp // aCollect[#, S] &;
   tem[[1]] = Identity @@@ (List @@ #) & /@ tem[[1]];
   result = SparseArray[Thread@Rule[tem[[1]], tem[[2]]]];
   If[OptionValue["dim"] === 0, result,
    SparseArray[result,
     ConstantArray[OptionValue["dim"], Length[Dimensions[result]]]]]
   ];
"""


def alphabet_file_expr_loader(letter_var: str, roots_var: str, file_abs: str) -> str:
    """A Get[...] + alphaExpr loader for an alphabet file in a known dialect."""
    return (
        f'Get[{_q(file_abs)}];\n'
        f'alphaExpr = {letter_var}[[All, 2]] /. {roots_var};'
    )


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
    if ptype == "letter_symmetry":
        sname = _safe_name((prop.get("name") or "").strip() or prop.get("params", {}).get("name", "sym"))
        return f"data/{sname}.wxf"
    if ptype == "sparse_expression":
        sname = _safe_name((prop.get("name") or "").strip() or "sparse")
        return f"data/{sname}.wxf"
    raise ValueError(f"unknown property type {ptype}")


def _symb_root_snippet() -> str:
    # The package path is relocatable: SYMBOLOGY_ROOT (if set) overrides the
    # baked-in absolute path, so generated scripts can travel to another
    # machine / checkout without rewriting. Lookup on the full environment
    # association avoids GetEnvironment's rule/None return-shape pitfalls.
    return (
        '$symbEnv = Lookup[GetEnvironment[], "SYMBOLOGY_ROOT", Null];\n'
        '$symbRoot = If[StringQ[$symbEnv] && $symbEnv =!= "", $symbEnv, '
        f'{_q(REPO_ROOT.as_posix())}];\n'
        f'Get[FileNameJoin[{{$symbRoot, {_q("SymbolBootstrap.wl")}}}]];'
    )


def _preamble(alphabet: dict) -> str:
    lines = [
        _symb_root_snippet(),
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
    elif ptype == "letter_symmetry":
        rule = (params.get("rule") or "").strip()
        defs = (params.get("defs_file") or "").strip()
        if defs:
            parts.append(f'If[Get[{_q(defs)}] === $Failed, Print["@@RESULT@@FAIL"]; Exit[1]];')
        parts.append(f'$rule = ToExpression[{_q(rule)}];')
        parts.append('$tensor = SparseArray[CoefficientArrays[$letters /. $rule, $letters][[2]]];')
        n = len(alphabet["letters"])
        parts.append(f'If[Dimensions[$tensor] =!= {{{n}, {n}}}, Print["@@RESULT@@FAIL"]; Exit[1]];')
    elif ptype == "sparse_expression":
        expr = (params.get("expression") or "").strip()
        dim = params.get("dim") or 0
        parts.append("ClearAll[S];")
        parts.append(EXP2SA_DEFS)
        parts.append(f'$exp = Quiet[Check[ToExpression[{_q(expr)}], $Failed]];')
        parts.append('If[$exp === $Failed, Print["could not parse the symbol tensor expression"]; Print["@@RESULT@@FAIL"]; Exit[1]];')
        parts.append('If[FreeQ[$exp, _S], Print["no S[...] terms found in the expression"]; Print["@@RESULT@@FAIL"]; Exit[1]];')
        parts.append('$sa = Exp2SA[$exp];')
        parts.append('If[!ArrayQ[$sa] || Length[$sa["NonzeroPositions"]] === 0, Print["the expression has no nonzero S[...] terms"]; Print["@@RESULT@@FAIL"]; Exit[1]];')
        if dim:
            parts.append(f'If[Max[Flatten[$sa["NonzeroPositions"]]] > {dim}, Print["an index exceeds the dimension {dim}"]; Print["@@RESULT@@FAIL"]; Exit[1]];')
            parts.append(f'$tensor = SparseArray[$sa, ConstantArray[{dim}, Length[Dimensions[$sa]]]];')
        else:
            parts.append('$tensor = $sa;')
    else:
        raise ValueError(f"unknown property type {ptype}")

    parts.append(_finish(out_abs))
    return "\n".join(parts) + "\n"


def merge_script(input_abs_paths: list, out_abs: str) -> str:
    files = "{" + ",".join(_q(p) for p in input_abs_paths) + "}"
    return f'''{_symb_root_snippet()}
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
