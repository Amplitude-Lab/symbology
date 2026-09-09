from __future__ import annotations

import hashlib
import json
import re
import shlex
import threading
from fractions import Fraction
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
    ("solve_collinear", "seed"): {"fec1", "fec", "sew", "basis", "solution", "boundary", "tensor", "matrix"},
    ("solve_collinear", "rhs"): {"boundary", "tensor", "basis", "solution", "matrix"},
    ("solve_collinear", "cond"): {"matrix", "tensor"},
    ("projection_chain", "seed"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("symmetry_invariant", "seed"): {"fec1", "fec", "lec1", "lec", "sew"},
    ("compute_rhs", "e1"): TENSOR_KINDS,
    ("compute_rhs", "p1"): TENSOR_KINDS,
    ("add_tensors", None): TENSOR_KINDS,
    ("ternary_contract", "tensor"): R3_KINDS,
    ("ternary_contract", "trans1"): TENSOR_KINDS,
    ("ternary_contract", "trans2"): TENSOR_KINDS,
    ("apply_symmetry", "tensor"): R3_KINDS,
    ("apply_symmetry", "trans1"): TENSOR_KINDS,
    ("apply_symmetry", "trans2"): TENSOR_KINDS,
    ("apply_symmetry", "sym"): {"matrix"},
    ("matrix_power", "matrix"): {"matrix"},
    ("tensor_join", "a"): TENSOR_KINDS,
    ("tensor_join", "b"): TENSOR_KINDS,
    ("tensor_dot", "a"): TENSOR_KINDS,
    ("tensor_dot", "b"): TENSOR_KINDS,
    ("squeeze_tensor", "in"): TENSOR_KINDS,
    ("shuffle_product", "a"): TENSOR_KINDS,
    ("shuffle_product", "b"): TENSOR_KINDS,
    ("impose_integrability", "tensor"): R3_KINDS,
    ("impose_integrability", "dlogmat"): {"dlogmat"},
    ("integrability_condition", "tensor"): R3_KINDS,
    ("integrability_condition", "dlogmat"): {"dlogmat"},
    ("solve_conditions", "cond"): {"matrix", "tensor"},
    ("assemble", None): TENSOR_KINDS,
    ("expand_tensor", None): TENSOR_KINDS,
}

# Compute RHS object types. Mirrored in web/src/flowdefs.js::RHS_MODES — keep
# the two registries in sync. requires/forbids constrain which seed ports the
# mode expects wired (e1 always, p1 for the pair modes); adding a new object
# type = one entry here + one entry in the JS registry + the generation branch
# in _compile_compute_rhs.
RHS_MODES = {
    "mhv_boundary": {
        "label": "MHV boundary (boundary_L)",
        "requires": ("e1",),
        "forbids": ("p1",),
    },
    "e47me67_te": {
        "label": "NMHV E47mE67 (tE_L)",
        "requires": ("e1", "p1"),
        "forbids": (),
    },
}


RAT_RE = re.compile(r"^-?\d+(/\d+)?$")


def _abs(proj_dir: Path, rel: str) -> str:
    p = Path(rel)
    return str(p if p.is_absolute() else proj_dir / rel)


def _resolve_path_arg(proj_dir: Path, value: str) -> str:
    if value in ("0", "identity", "divergent", "finite"):
        return value
    return _abs(proj_dir, value)


def _safe_file_arg(proj_dir: Path, value, what: str, errors) -> str | None:
    """Validate a free-text file reference from the graph before it reaches a
    CLI argv: reject anything flag-shaped (argument injection) and relative
    traversal ('..'). Existence is deliberately NOT required — the file may be
    produced by an earlier step of the same run; a genuinely missing file
    fails loudly inside the binary (file_to_ustr throws). Returns the absolute
    path string, or None after appending a compile error."""
    v = (value or "").strip()
    if not v:
        return None
    if v.startswith("-"):
        errors.append(f"{what}: '{value}' looks like a command-line flag, not a file path.")
        return None
    p = Path(v)
    if p.is_absolute():
        return str(p)
    if ".." in p.parts:
        errors.append(f"{what}: relative path may not contain '..' ({v}).")
        return None
    return str(proj_dir / p)


def _safe_flag_value(value, what: str, errors) -> str | None:
    """Validate a free-text value that becomes a bare argv element after a
    --flag (out stems, targets). It must not look like a flag itself."""
    v = (value or "").strip()
    if not v:
        return None
    if v.startswith("-"):
        errors.append(f"{what}: '{value}' may not start with '-' (it would be parsed as a flag).")
        return None
    return v


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
        if ttype in ("merge_conditions", "assemble", "add_tensors", "expand_tensor"):
            allowed = EDGE_RULES[(ttype, None)]
        else:
            allowed = EDGE_RULES.get((ttype, th))
            if allowed is None and th.startswith("seed"):
                allowed = EDGE_RULES.get((ttype, "seed"))
            if allowed is None and th.startswith("in_seed_"):
                allowed = EDGE_RULES.get((ttype, "seed"))
            if allowed is None and th.startswith("in_rhs_"):
                allowed = EDGE_RULES.get((ttype, "rhs"))
        if allowed is None:
            continue
        if kind not in allowed:
            errors.append(
                f"Type mismatch: a '{kind}' output cannot feed the '{th}' input of a {ttype} node."
            )

    if errors:
        return {"ok": False, "errors": errors, "steps": []}
    # Meaning persistence: annotate each producing step with the letters-axis
    # meaning / chain provenance of its outputs. The run engine writes
    # <output>.meaning.json sidecars when a step succeeds, and Reuse Output
    # nodes read them back — a leaf node has no upstream graph to walk, so
    # the provenance must travel with the file instead.
    meaning_by_file = {}
    for info in provides.values():
        f = (info or {}).get("file")
        if not f or any(part.startswith(".") for part in Path(f).parts):
            continue
        m = {k: info[k] for k in ("kind", "weight", "letters_axes", "axes_meaning", "chain")
             if info.get(k) is not None}
        if "letters_axes" in m or "chain" in m:
            meaning_by_file[f] = m
    for s in steps:
        tm = {o: meaning_by_file[o] for o in (s.get("outputs") or []) if o in meaning_by_file}
        if tm:
            s["tensor_meaning"] = tm
    shapes = {f"{nid}|{handle}": info.get("dims") for (nid, handle), info in provides.items() if info.get("dims")}
    flow_outputs = []
    seen_files = set()
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
                seen_files.add(info["file"])
    if not flow_outputs:
        out_handles = {}
        for n in nodes.values():
            ntype = n.get("type")
            if ntype in ("cb_in", "cb_out", "alphabet", "reuse_output"):
                continue
            out_handles[n["id"]] = ntype
        for (nid, handle), info in provides.items():
            if nid not in out_handles or not info or not info.get("file"):
                continue
            if info["file"] in seen_files:
                continue
            seen_files.add(info["file"])
            flow_outputs.append({
                "name": Path(info["file"]).stem,
                "file": info["file"],
                "kind": info.get("kind"),
                "dims": info.get("dims"),
            })
    return {"ok": True, "errors": [], "shapes": shapes, "flow_outputs": flow_outputs,
            "steps": [{k: v for k, v in s.items() if k not in ("argv", "tensor_meaning")} for s in steps],
            "_steps_full": steps}


# ============================================================================
# Standalone script export — "design on the laptop, run on the cluster".
#
# The run engine (jobs.py) executes compiled steps directly, which ties a run
# to the machine the front-end is on. export_flow_script renders the SAME
# compiled steps as a portable bash script with self-checks and safety guards:
#   - set -euo pipefail; loud, context-rich failures (never silent skips)
#   - preflight: binaries executable, seed inputs present with sizes, disk
#     space on the output dir, wolframscript presence (no-Mathematica mode
#     skips Wolfram steps with loud warnings when their outputs already exist)
#   - per-step banners, wall timings, full log tee'd to runs/exported-*.log
#   - post-step verification: every declared output exists and is nonempty,
#     with a CRC32 echo for traceability
#   - resume: a fingerprint cache (command + inputs) writes .sig.exported
#     sidecars; re-running skips completed steps (independent from the
#     server's .sig cache — a mismatch costs a recompute, never a stale skip)
#   - paths relativized to $SYMBOLOGY_ROOT / $PROJ_DIR / $WOLFRAMSCRIPT, so
#     the script + project dir + repo travel together
#   - --dry-run prints the plan; STEPS="3,5-9" runs a subset
# ============================================================================

def _portable_argv_element(arg: str, repo_root: str, proj_dir: str, kind: str, is_bin: bool) -> str:
    """Render one argv element as a bash word, relativized where possible."""
    if is_bin and kind == "wolfram":
        return '"$WOLFRAMSCRIPT"'
    for prefix, var in ((proj_dir + "/", '"$PROJ_DIR"/'),
                        (repo_root + "/", '"$SYMBOLOGY_ROOT"/')):
        if arg.startswith(prefix):
            rest = arg[len(prefix):]
            return var + shlex.quote(rest)
    return shlex.quote(arg)


def export_flow_script(proj: dict, graph: dict, flow_name: str, output_subdir: str | None = None) -> dict:
    from .config import REPO_ROOT

    result = compile_flow(proj, graph, output_subdir)
    if not result["ok"]:
        return {"ok": False, "errors": result["errors"]}
    steps = result["_steps_full"]
    if not steps:
        return {"ok": False, "errors": ["the flow compiled to zero steps — nothing to export"]}

    proj_dir = storage.project_dir(proj["id"])
    repo_root = REPO_ROOT.as_posix()
    wolfram = find_wolframscript()

    lines: list[str] = []
    w = lines.append
    w("#!/usr/bin/env bash")
    w(f"# Generated by the symbology flow editor: project '{proj['id']}', flow '{flow_name}'.")
    w("# Portable re-run of the compiled step list — works anywhere the C++ core is")
    w("# built (see README: build with make, set SYMBOLOGY_ROOT to the repo root).")
    w("# Self-checking: preflight + per-step output verification + CRC32 echoes.")
    w("# Resume: re-running skips steps whose outputs and .sig.exported fingerprints")
    w("# match. Env: SYMBOLOGY_ROOT, PROJ_DIR, WOLFRAMSCRIPT, STEP_TIMEOUT (s, 0=off),")
    w("#            STEPS='3,5-9' subset, WOLFRAM_MODE=skip|fail (default: skip only")
    w("#            when the outputs already exist, else fail). --dry-run prints the plan.")
    w("set -euo pipefail")
    w("")
    w(f'SYMBOLOGY_ROOT="${{SYMBOLOGY_ROOT:-{shlex.quote(repo_root)}}}"')
    w('PROJ_DIR="${PROJ_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"')
    w(f'WOLFRAMSCRIPT="${{WOLFRAMSCRIPT:-{shlex.quote(wolfram or "wolframscript")}}}"')
    w('STEP_TIMEOUT="${STEP_TIMEOUT:-0}"')
    w('WOLFRAM_MODE="${WOLFRAM_MODE:-auto}"')
    w('STEPS_FILTER="${STEPS:-}"')
    w('DRY_RUN=0; [[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1')
    w('RUN_TAG="$(date +%Y%m%d-%H%M%S)"')
    w('LOG_DIR="$PROJ_DIR/runs"; mkdir -p "$LOG_DIR"')
    w('LOG_FILE="$LOG_DIR/exported-${RUN_TAG}.log"')
    w('CUR_STEP="(preflight)"')
    w("")
    w('fail() { echo; echo "=== EXPORTED SCRIPT FAILURE ===" >&2; echo "step : $CUR_STEP" >&2;')
    w('      echo "log  : $LOG_FILE (tail below)" >&2; tail -40 "$LOG_FILE" >&2 || true;')
    w('      echo "disk : $(df -h "$PROJ_DIR" 2>/dev/null | tail -1)" >&2; exit 1; }')
    w('note() { echo "[export] $*"; }')
    w("")
    w('# ---- helpers ----')
    w('crc32() { python3 - "$1" <<\'PY\'')
    w("import sys, zlib")
    w("print(f'{zlib.crc32(open(sys.argv[1], \"rb\").read()):08x}')")
    w("PY")
    w("}")
    w("")
    w('# fingerprint: sha256 over the command and the (size, mtime) of every')
    w('# existing input file the step reads. Independent from the server-side')
    w('# .sig cache — a mismatch costs a recompute, never a stale skip.')
    w("fingerprint() {")
    w('  python3 - "$@" <<\'PY\'')
    w("import hashlib, sys, os")
    w("argv = sys.argv[1:]")
    w("outputs = set(os.environ.get('STEP_OUTPUTS', '').split())")
    w("h = hashlib.sha256()")
    w("h.update('\\x00'.join(argv).encode())")
    w("for a in argv[1:]:")
    w("    if a in outputs or not os.path.isfile(a):")
    w("        continue")
    w("    st = os.stat(a)")
    w("    h.update(f'{a}:{st.st_size}:{int(st.st_mtime)};'.encode())")
    w("print(h.hexdigest())")
    w("PY")
    w("}")
    w("")
    w('step_wanted() {')
    w('  [[ -z "$STEPS_FILTER" ]] && return 0')
    w('  local n="$1" part lo hi')
    w('  IFS="," read -ra parts <<< "$STEPS_FILTER"')
    w('  for part in "${parts[@]}"; do')
    w('    if [[ "$part" == *-* ]]; then lo=${part%-*}; hi=${part#*-}; else lo=$part; hi=$part; fi')
    w('    [[ $n -ge $lo && $n -le $hi ]] && return 0')
    w('  done')
    w('  return 1')
    w("}")
    w("")
    w("# ---- preflight ----")
    w('note "SYMBOLOGY_ROOT = $SYMBOLOGY_ROOT"')
    w('note "PROJ_DIR      = $PROJ_DIR"')
    w('note "log file      = $LOG_FILE"')
    w(": > \"$LOG_FILE\" || fail")
    binaries = sorted({s["argv"][0] for s in steps})
    bin_vars = []
    for b in binaries:
        kind = next(s["kind"] for s in steps if s["argv"][0] == b)
        if kind == "wolfram":
            bin_vars.append('"$WOLFRAMSCRIPT" (wolfram steps)')
            continue
        rel = b[len(repo_root):] if b.startswith(repo_root) else None
        if rel is None:
            bin_vars.append(shlex.quote(b))
        else:
            bin_vars.append(f'"$SYMBOLOGY_ROOT"/{shlex.quote(rel.lstrip("/"))}')
    for bv, raw in zip(bin_vars, binaries):
        w(f'[[ -x {bv} ]] || {{ note "MISSING BINARY: {raw}"; fail; }}')
        w(f'note "binary ok: {raw}"')
    # seed inputs that exist at export time (flow-produced files are checked
    # by the per-step output verification of their producing step)
    inputs_seen = []
    outputs_all = {o for s in steps for o in (s.get("outputs") or [])}
    for s in steps:
        for a in s["argv"][1:]:
            if a.endswith(".wxf") and not a.startswith("-"):
                p = Path(a)
                if p.is_file() and a not in inputs_seen:
                    inputs_seen.append(a)
    proj_prefix = proj_dir.as_posix() + "/"
    repo_prefix = repo_root + "/"

    def _shown_word(a: str) -> str:
        # keep the variable prefix OUTSIDE the shell-quoting so it expands
        if a.startswith(proj_prefix):
            return '"$PROJ_DIR"/' + shlex.quote(a[len(proj_prefix):])
        if a.startswith(repo_prefix):
            return '"$SYMBOLOGY_ROOT"/' + shlex.quote(a[len(repo_prefix):])
        return shlex.quote(a)

    for a in inputs_seen:
        rel_a = a[len(proj_prefix):] if a.startswith(proj_prefix) else a
        if rel_a in outputs_all:
            continue  # produced by an earlier step — verified per-step, not here
        shown = _shown_word(a)
        size = Path(a).stat().st_size
        w(f'[[ -s {shown} ]] || {{ note "MISSING INPUT: {shown}"; fail; }}')
        w(f'note "input ok : {shown} ({size} bytes)"')
    w('note "free space on PROJ_DIR: $(df -h "$PROJ_DIR" | tail -1)"')
    has_wolfram = any(s["kind"] == "wolfram" for s in steps)
    if has_wolfram:
        w('if ! command -v "$WOLFRAMSCRIPT" >/dev/null 2>&1 && ! [[ -x "$WOLFRAMSCRIPT" ]]; then')
        w('  note "WARNING: wolframscript not found — Wolfram steps will be handled per WOLFRAM_MODE=$WOLFRAM_MODE"')
        w("fi")
    w("")
    w("# ---- steps ----")
    w('# (bash-3.2-safe: ${arr[@]+...} guards for possibly-empty arrays under set -u)')
    w("run_step() {")
    w('  local n="$1" label="$2" kind="$3" skipflag="$4" outs_str="$5"; shift 5')
    w('  local -a argv=("$@")')
    w('  local -a outs=()')
    w('  [[ -n "$outs_str" ]] && read -ra outs <<< "$outs_str"')
    w('  CUR_STEP="[$n] $label"')
    w('  if ! step_wanted "$n"; then note "[$n] SKIPPED by STEPS filter"; return 0; fi')
    w('  if [[ $DRY_RUN == 1 ]]; then note "DRY-RUN $desc :: ${argv[*]}"; return 0; fi')
    w('  local desc="[$n/$N_STEPS] $label"')
    w('  local -a pre=()')
    w('  if [[ "$STEP_TIMEOUT" != "0" ]]; then pre=(timeout "$STEP_TIMEOUT"); fi')
    w('  # resume check: all outputs exist (parity with the engine: exists())')
    w('  # and the fingerprint matches')
    w('  local cached=1 fp=""')
    w('  for o in ${outs[@]+"${outs[@]}"}; do [[ -e "$PROJ_DIR/$o" ]] || cached=0; done')
    w('  if [[ "$skipflag" == "1" && $cached == 1 && ${#outs[@]} -gt 0 ]]; then')
    w('    STEP_OUTPUTS="$outs_str" fp="$(fingerprint "${argv[@]}")" || fail')
    w('    for o in ${outs[@]+"${outs[@]}"}; do')
    w('      sig="$PROJ_DIR/$o.sig.exported"')
    w('      [[ -f "$sig" && "$(cat "$sig")" == "$fp" ]] || cached=0')
    w("    done")
    w("  else")
    w('    cached=0')
    w("  fi")
    w('  if [[ $cached == 1 ]]; then note "$desc: CACHED (outputs + fingerprint match)"; return 0; fi')
    w('  # Wolfram steps without wolframscript')
    w('  if [[ "$kind" == "wolfram" ]] && ! command -v "$WOLFRAMSCRIPT" >/dev/null 2>&1 && ! [[ -x "$WOLFRAMSCRIPT" ]]; then')
    w('    local have_all=1; for o in ${outs[@]+"${outs[@]}"}; do [[ -s "$PROJ_DIR/$o" ]] || have_all=0; done')
    w('    if [[ $have_all == 1 || "$WOLFRAM_MODE" == "skip" ]]; then')
    w('      note "$desc: WARNING — wolframscript unavailable; SKIPPING the Wolfram step"')
    w('      note "   (outputs: ${outs_str:-none}; copy them from a machine that has Mathematica,"')
    w('      note "    or recompute there and sync this project dir)"')
    w("      return 0")
    w("    fi")
    w('    note "$desc: wolframscript unavailable and outputs missing (WOLFRAM_MODE=$WOLFRAM_MODE)"')
    w("    fail")
    w("  fi")
    w('  for o in ${outs[@]+"${outs[@]}"}; do mkdir -p "$PROJ_DIR/$(dirname "$o")"; done')
    w('  note "$desc: ${argv[0]} ..."')
    w('  local t0=$(date +%s) rc=0')
    w('  set +o pipefail')
    w('  { cd "$PROJ_DIR" && ${pre[@]+"${pre[@]}"} "${argv[@]}"; } 2>&1 | tee -a "$LOG_FILE"')
    w('  rc=${PIPESTATUS[0]}')
    w('  set -o pipefail')
    w('  note "$desc: exit=$rc elapsed=$(( $(date +%s) - t0 ))s"')
    w('  [[ $rc == 0 ]] || fail')
    w('  # post-step verification: every output must exist (files nonempty;')
    w('  # directory outputs are legal — scratch roots like output/.derived/...)')
    w('  for o in ${outs[@]+"${outs[@]}"}; do')
    w('    if [[ -d "$PROJ_DIR/$o" ]]; then')
    w('      note "   output $o  (directory)"')
    w('    elif [[ -s "$PROJ_DIR/$o" ]]; then')
    w('      note "   output $o  crc32=$(crc32 "$PROJ_DIR/$o")  $(stat -c %s "$PROJ_DIR/$o" 2>/dev/null || stat -f %z "$PROJ_DIR/$o") bytes"')
    w('    else')
    w('      note "step exited 0 but output $o is missing/empty"; fail')
    w('    fi')
    w("  done")
    w('  STEP_OUTPUTS="$outs_str" fp="$(fingerprint "${argv[@]}")" || true')
    w('  for o in ${outs[@]+"${outs[@]}"}; do echo "$fp" > "$PROJ_DIR/$o.sig.exported"; done')
    w("}")
    w("")
    w("N_STEPS=%d" % len(steps))
    for i, s in enumerate(steps, start=1):
        argv_words = [_portable_argv_element(a, repo_root, proj_dir.as_posix(), s["kind"], j == 0)
                      for j, a in enumerate(s["argv"])]
        outs_str = " ".join(shlex.quote(o) for o in (s.get("outputs") or []))
        w(f'run_step {i} {shlex.quote(s["label"])} {shlex.quote(s["kind"])} {1 if s.get("skip_if_exists") else 0} {shlex.quote(outs_str)} \\')
        w("  " + " ".join(argv_words))
        w("")
    w('note "ALL STEPS COMPLETE — final outputs:"')
    final_files = sorted({o for f in result.get("flow_outputs", []) for o in [f.get("file")] if o})
    for f in final_files:
        shown = _shown_word(str(proj_dir / f)) if not f.startswith("/") else shlex.quote(f)
        w(f'[[ -s {shown} ]] && note "  {f}  crc32=$(crc32 {shown})" || note "  {f}: MISSING"')
    w('note "full log: $LOG_FILE"')

    safe_flow = re.sub(r"[^A-Za-z0-9_.-]+", "_", flow_name).strip("_") or "flow"
    export_dir = proj_dir / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)
    script_path = export_dir / f"{safe_flow}.sh"
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)
    return {
        "ok": True,
        "path": str(script_path),
        "n_steps": len(steps),
        "n_wolfram_steps": sum(1 for s in steps if s["kind"] == "wolfram"),
        "flow_outputs": final_files,
    }


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
            _compile_solve_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir,
                                    nodes, edges, bootstrap, file_steps)
        elif ntype == "symderive":
            _compile_symderive(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "solve_collinear":
            _compile_solve_collinear(node, incoming, outgoing, provides, add_step, errors, bootstrap, proj_dir, tensor_ops_bin)
        elif ntype == "projection_chain":
            _compile_projection_chain(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "symmetry_invariant":
            _compile_symmetry_invariant(node, incoming, provides, add_step, errors, bootstrap, proj_dir)
        elif ntype == "compute_rhs":
            _compile_compute_rhs(node, incoming, provides, add_step, errors, tensor_ops_bin, tensor_add_bin, proj_dir)
        elif ntype == "add_tensors":
            _compile_add_tensors(node, incoming, provides, add_step, errors, tensor_add_bin, proj_dir)
        elif ntype == "ternary_contract":
            _compile_ternary_contract(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "apply_symmetry":
            _compile_apply_symmetry(node, incoming, provides, add_step, errors, bootstrap, tensor_ops_bin, proj_dir, nodes, edges, file_steps)
        elif ntype == "matrix_power":
            _compile_matrix_power(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "tensor_join":
            _compile_tensor_join(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "tensor_dot":
            _compile_tensor_dot(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "squeeze_tensor":
            _compile_squeeze_tensor(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "shuffle_product":
            _compile_shuffle_product(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "expand_tensor":
            _compile_expand_tensor(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "impose_integrability":
            _compile_impose(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "integrability_condition":
            _compile_integrability_condition(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "solve_conditions":
            _compile_solve_conditions(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "assemble":
            _compile_assemble(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        elif ntype == "reuse_output":
            _compile_reuse_output(node, provides, errors, proj_dir)
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


def _matrix_meaning(prop):
    """How a matrix property acts on the meaning of a letters axis.

    'preserve'  — the axis stays in the same alphabet basis (symmetry reps,
                  letter_symmetry maps): provenance survives the application.
    'collinear' — the axis moves to the collinear limit of the basis: the
                  meaning changes, further collinear projections are invalid.
    None        — unknown / custom matrix: no meaning claim.
    """
    action = prop.get("meaning_action")
    if action in ("preserve", "collinear"):
        return action
    name = (prop.get("name") or "").lower() + " " + (prop.get("id") or "")
    if "collinear" in name or "colmat" in name or "colproj" in name:
        return "collinear"
    if prop.get("type") in ("letter_symmetry",) or "repmat" in name or "invmap" in name:
        return "preserve"
    return None


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
            "meaning": _matrix_meaning(prop),
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
                "meaning": "preserve",
            }


def _compile_reuse_output(node, provides, errors, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    rel = (data.get("file") or "").strip().lstrip("/")
    kind = data.get("kind") or "tensor"
    if not rel:
        errors.append("A Reuse Output node has no output selected.")
        return
    if kind not in TENSOR_KINDS:
        errors.append(f"Reused output '{rel}': unknown tensor kind '{kind}'.")
        return
    if ".." in Path(rel).parts:
        errors.append(f"Reused output path '{rel}' must stay inside the project directory.")
        return
    if not (proj_dir / rel).exists():
        errors.append(f"Reused output '{rel}' does not exist yet — run the flow that produces it first.")
        return
    info = {"kind": kind, "file": rel, "weight": _infer_weight(rel), "name": Path(rel).stem, "dims": None}
    info.update(_load_meaning_sidecar(proj_dir, rel))
    provides[(nid, "out")] = info


def _load_meaning_sidecar(proj_dir, rel: str) -> dict:
    """Restore the letters-axis meanings / chain provenance of a tensor from
    its <file>.meaning.json sidecar (written by the run engine when the
    producing step succeeded). A Reuse Output node is a graph leaf, so the
    provenance cannot be recovered by walking upstream — it has to travel
    with the file. The size guard drops sidecars left stale by a tensor that
    was regenerated outside the flow engine."""
    sidecar = proj_dir / (rel + ".meaning.json")
    tensor = proj_dir / rel
    try:
        payload = json.loads(sidecar.read_text())
        meaning = payload.get("meaning") or {}
        if payload.get("size") != tensor.stat().st_size:
            return {}
        return {k: meaning[k] for k in ("letters_axes", "axes_meaning", "chain")
                if meaning.get(k) is not None}
    except (OSError, ValueError, AttributeError):
        return {}


def _infer_weight(rel: str):
    stem = Path(rel).stem
    m = re.match(r"^(?:FEC|LEC|first_w|last_w)(\d+)$", stem, re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.match(r"^SEW_\d+p(\d+)$", stem, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def flow_outputs_catalog(proj: dict) -> list[dict]:
    """Enumerate the declared outputs of every non-custom-block flow, for the
    'flow outputs' palette section. Uses compile-time info only; flows that do
    not compile are skipped silently (they have no usable outputs anyway)."""
    catalog = []
    for f in proj.get("flows", []):
        if f.get("custom_block"):
            continue
        try:
            result = compile_flow(proj, f.get("graph") or {}, output_subdir=f.get("output_subdir"))
        except Exception:
            continue
        if not result.get("ok"):
            continue
        outs = result.get("flow_outputs") or []
        if not outs:
            continue
        catalog.append({
            "flow_id": f["id"],
            "flow_name": f.get("name") or "Untitled flow",
            "outputs": [
                {"name": o["name"], "file": o["file"], "kind": o.get("kind") or "tensor", "dims": o.get("dims")}
                for o in outs
            ],
        })
    return catalog


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
        if target <= w_in:
            errors.append(f"Extend node target weight {target} must be greater than the input weight {w_in}.")
            return
        w_out = target
    cur_file, cur_w = seed["file"], w_in
    while cur_w < w_out:
        nxt = cur_w + 1
        rel = f"{_out()}{direction}_{nxt}.wxf"
        last = nxt == w_out
        add_step(
            f"Extend {direction}_{cur_w} -> {direction}_{nxt}",
            "bootstrap",
            [bootstrap, "--extend", "-c", _abs(proj_dir, cond["file"]), flag, _abs(proj_dir, cur_file), "-o", _abs(proj_dir, rel)],
            proj_dir,
            [rel],
            True,
            {"type": "extend", "tensor_file": rel, **({} if last else {"background": True})},
        )
        cur_file, cur_w = rel, nxt
    provides[(nid, "fec" if fec is not None else "lec")] = {
        "kind": "fec" if fec is not None else "lec",
        "file": cur_file,
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
    if base != f"SEW_{fw}p{lw}":
        base = _safe_flag_value(base, "A Sew node target", errors) or f"SEW_{fw}p{lw}"
        if "/" in base or not SAFE_TARGET_RE.match(base):
            errors.append(f"A Sew node has an invalid target name '{base}' (letters, digits, _, ., - only).")
            return
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


def _compile_solve_symmetry(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir,
                            nodes=None, edges=None, bootstrap=None, file_steps=None):
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

    sb_arg = _abs(proj_dir, matrix["file"])
    # A compressed middle axis (dim(1) in a chain basis, e.g. (229,97,42) with
    # 97 = weight-3 FEC basis): derive the INDUCED square map for that axis via
    # bootstrap --project and pass the letter matrix explicitly as Sc (instead
    # of "-", which would wrongly reuse the induced 97x97 map on the 42 axis).
    tensor_edge = None
    for e in incoming.get(nid, []):
        if (e.get("targetHandle") or "") == "tensor":
            tensor_edge = e
            break
    ch = (tensor or {}).get("chain")
    if ch and bootstrap is not None:
        # Chain-basis first axis (e.g. (229,97,42) = SEW_3p1 x FEC_3 x letters):
        # run_symsolve's R-derivation is ONLY valid on letter axes (it needs the
        # symmetry to act directly on axis 0's space). For a chain axis use the
        # projection pipeline's own solver: bootstrap --solve-symmetry computes the
        # induced action M (e.g. SEW_3p1.wxf 229x229) from the chain files and
        # returns K = ker(M^T - I); the invariant projection of T is then K . T
        # (tdot axis 2 / 1) — the exact analog of what symsolve does internally
        # for uniform letter axes.
        src = (tensor_edge or {}).get("source")
        head, sub = _expand_chain_provenance(provides, incoming, src, (tensor_edge or {}).get("sourceHandle") or "out",
                                             nodes or {}, edges or [])
        if head is not None:
            stem = Path(matrix["file"]).stem
            sym_name, usable = _symmetry_run_name(stem)
            if sym_name is None and not usable:
                errors.append(f"Symmetry Solve: cannot use '{stem}' as a symmetry name (unsafe characters).")
                return
            if sym_name is None:
                sym_name = stem
            need_bootstrap = head["kind"] == "sew" or any(
                e.get("weight", 1) > 1 for e in (sub.get("first") or []) + (sub.get("last") or []))
            maps = _emit_derived_projection(head, sub, _derived_matrix_sig(proj_dir, matrix["file"]),
                                            sym_name, matrix["file"], add_step, proj_dir, bootstrap,
                                            file_steps, errors, needed=need_bootstrap)
            if maps is not None:
                if head["kind"] == "sew":
                    fw = ((head.get("first") or [{}])[-1].get("weight") if head.get("first") else 1) or 1
                    lw = ((head.get("last") or [{}])[-1].get("weight") if head.get("last") else 1) or 1
                    target_name = f"SEW_{fw}p{lw}"
                else:
                    errors.append(
                        "Symmetry Solve: the chain head is not a sew tensor — "
                        "the invariant-space solver needs a SEW head (SEW_FpL)."
                    )
                    return
                inv_rel = maps["_root"] + f"/output/{sym_name}/{target_name}_invariant.wxf"
                if inv_rel not in file_steps:
                    file_steps.add(inv_rel)
                    add_step(
                        f"Solve invariant space of {target_name} under {sym_name} (shared cache)",
                        "bootstrap",
                        [bootstrap, "--solve-symmetry", "--symmetry", sym_name, "--target", target_name,
                         "--data-dir", _abs(proj_dir, f"{maps['_root']}/data"),
                         "--output-dir", _abs(proj_dir, f"{maps['_root']}/output")],
                        proj_dir,
                        [inv_rel],
                        True,
                        {"type": "symmetry_invariant", "symmetry": sym_name, "target": target_name},
                    )
                rel = f"{_out()}{target}.wxf"
                add_step(
                    f"Symmetry solve {tensor['name']} -> {target} (invariant space of the chain axis)",
                    "tensor_ops",
                    [tensor_ops_bin, "tdot", _abs(proj_dir, inv_rel), _abs(proj_dir, tensor["file"]),
                     "2", "1", _abs(proj_dir, rel)],
                    proj_dir,
                    [rel],
                    True,
                    {"type": "solve_symmetry", "tensor_file": rel},
                )
                provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"),
                                          "name": target, "letters_axes": tensor.get("letters_axes"),
                                          "axes_meaning": tensor.get("axes_meaning"),
                                          "chain": tensor.get("chain")}
                return
        errors.append(
            f"Symmetry Solve: the first axis of '{tensor['name']}' is a compressed chain basis; "
            "cannot derive the induced symmetry map for it (missing chain files or bootstrap)."
        )
        return

    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Symmetry solve {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "symsolve", _abs(proj_dir, tensor["file"]), sb_arg,
         sc if sc == "-" else _abs(proj_dir, sc), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "solve_symmetry", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"),
                              "name": target, "letters_axes": tensor.get("letters_axes"), "axes_meaning": tensor.get("axes_meaning"),
                              "chain": tensor.get("chain")}


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


def _collinear_pairs(data: dict) -> list[dict]:
    """Per-pair config rows [{projection, rhs}]; row 0 is pair 1 (the fixed
    seed/rhs ports). New graphs store data.pairs; legacy graphs stored
    data.rhs / data.letter_projection (pair 1) plus data.pairs_config (a
    flat list of letter projections for pairs 2+). The projection is a
    CLI sentinel ('identity' / 'divergent' / 'finite') passed through
    verbatim, or a .wxf file path (legacy per-slot contraction)."""
    def _proj(v) -> str:
        v = str(v or "identity").strip()
        return v or "identity"

    def _row(p) -> dict:
        p = p if isinstance(p, dict) else {}
        return {"projection": _proj(p.get("projection")), "rhs": str(p.get("rhs") or "").strip()}

    raw = data.get("pairs")
    if isinstance(raw, list) and raw:
        return [_row(p) for p in raw]
    rows = [{
        "projection": _proj(data.get("letter_projection")),
        "rhs": str(data.get("rhs") or "").strip(),
    }]
    for s in data.get("pairs_config") or []:
        if isinstance(s, dict):
            rows.append(_row(s))
        elif str(s or "").strip():
            rows.append({"projection": _proj(s), "rhs": ""})
    return rows


def _compile_solve_collinear(node, incoming, outgoing, provides, add_step, errors, bootstrap, proj_dir, tensor_ops_bin=None):
    data = node.get("data", {}) or {}
    nid = node["id"]
    seed = _edge_input(provides, incoming, nid, "seed")
    rhs_in = _edge_input(provides, incoming, nid, "rhs")
    cond_in = _edge_input(provides, incoming, nid, "cond")
    if bootstrap is None:
        errors.append("The bootstrap binary was not found.")
        return

    # Extra pairs arrive on dynamic in_seed_N / in_rhs_N ports (N >= 1; UI
    # pair N+1 — pair 1 is the fixed seed/rhs ports). Wired extra pairs force
    # multi-pair mode; so do a wired cond port (conditions re-ingest), pairs
    # configured in the Inspector's Pairs list, or — for custom seeds — a
    # request to export the combined conditions matrix.
    pair_edges: dict[int, dict] = {}
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        m = re.match(r"^in_seed_(\d+)$", th)
        if not m:
            continue
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            pair_edges.setdefault(int(m.group(1)), {})["seed"] = provides[src]
    for e in incoming.get(nid, []):
        th = e.get("targetHandle") or ""
        m = re.match(r"^in_rhs_(\d+)$", th)
        if not m:
            continue
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            pair_edges.setdefault(int(m.group(1)), {})["rhs"] = provides[src]
    wants_export = bool(data.get("export_conditions")) or any(
        (e.get("source") == nid) and (e.get("sourceHandle") == "conditions")
        for e in outgoing.get(nid, [])
    )

    # Seed selection: a wired fec/sew output keeps the named-target contract
    # (--target SEW_FpL / FEC_W); any other wired tensor switches to custom-seed
    # mode (--target-basis <file> --projection none), no naming convention needed.
    custom = False
    if seed is not None:
        if seed.get("kind") in ("fec1", "fec", "sew"):
            target = data.get("target") or seed.get("name")
        elif seed.get("file") is None:
            errors.append("Solve Collinear: the wired seed must be a concrete tensor file.")
            return
        else:
            custom = True
            target = None
    else:
        target = data.get("target")
    if target is not None:
        target = _safe_flag_value(target, "Solve Collinear target", errors)
        if target is None:
            return

    # Multi-pair is triggered by wiring (extra pairs or cond) or by adding
    # pairs in the Inspector's Pairs list (data.pairs rows beyond the first).
    # Cond-only (no seed at all) is allowed and re-solves ingested
    # conditions; it needs --out-stem, which we derive from the cond stems
    # when the user didn't set one.
    cfg_len = max(len(_collinear_pairs(data)) - 1, 0)
    multi = bool(pair_edges) or cond_in is not None or cfg_len > 0
    if multi:
        if cfg_len > 0 and seed is None and cond_in is None:
            errors.append(
                "Solve Collinear: pairs are configured in the Pairs list but the pair-1 seed port is "
                "not wired — wire a tensor into the seed port (all pairs are custom seeds)."
            )
            return
        if seed is not None and not custom:
            errors.append(
                "Solve Collinear: multi-pair mode needs concrete seed tensors — wire the pair-1 seed "
                "to a tensor (custom seed), not a named FEC/SEW target."
            )
            return
        if seed is None and pair_edges:
            errors.append(
                "Solve Collinear: multi-pair mode needs the pair-1 seed wired into the seed port "
                "(all pairs are custom seeds)."
            )
            return
        if seed is None and rhs_in is not None:
            errors.append("Solve Collinear: cond-only mode has no pairs — the rhs port is unused.")
            return
        if seed is None and (data.get("target") or "").strip():
            errors.append(
                "Solve Collinear: cond-only mode ignores the Target field — clear it (or wire a "
                "seed tensor for multi-pair mode)."
            )
            return
    elif wants_export:
        if not custom:
            errors.append(
                "Solve Collinear: exporting conditions needs a custom seed (wire a tensor into the "
                "seed port) or multi-pair/cond wiring."
            )
            return
        multi = True

    if not multi and not custom and not target:
        errors.append("Solve Collinear node: provide a target (e.g. SEW_3p1) or wire a tensor into the seed port.")
        return

    solver = data.get("solver") or "incremental"
    if solver not in ("incremental", "sampled"):
        errors.append("Solve Collinear node: solver must be incremental or sampled.")
        return

    if multi:
        _compile_solve_collinear_pairs(
            node, data, nid, seed, rhs_in, cond_in, pair_edges, wants_export,
            provides, add_step, errors, bootstrap, proj_dir, solver,
            incoming, tensor_ops_bin,
        )
        return

    # ---- single-pair mode (unchanged command line) ----
    # RHS: the wired rhs port wins over the text field; "0" means an all-zero RHS.
    pair_rows = _collinear_pairs(data)
    rhs = rhs_in["file"] if rhs_in is not None else pair_rows[0]["rhs"]
    if rhs_in is not None and rhs_in.get("file") is None:
        errors.append("Solve Collinear: the wired rhs must be a concrete tensor file.")
        return
    if not rhs:
        errors.append("Solve Collinear node: provide an RHS file (or '0'), or wire the rhs port.")
        return
    rhs_arg = "0"
    if rhs != "0":
        rhs_arg = _safe_file_arg(proj_dir, rhs, "Solve Collinear RHS", errors)
        if rhs_arg is None:
            return

    projection = data.get("projection") or "finite"
    if custom:
        projection = "none"  # custom seeds are used as-is; projection is meaningless
    elif projection == "none":
        errors.append("Solve Collinear: projection 'none' needs a custom seed — wire a tensor into the seed port.")
        return
    elif projection not in ("finite", "divergent"):
        errors.append("Solve Collinear node: projection must be finite or divergent.")
        return

    letter_proj = pair_rows[0]["projection"]
    if letter_proj not in ("identity", "divergent", "finite"):
        checked = _safe_file_arg(proj_dir, letter_proj, "Solve Collinear letter projection", errors)
        if checked is None:
            return
        letter_proj = checked

    if custom:
        seed_rel = seed["file"]
        label = f"Solve collinear constraints on custom seed {seed.get('name') or seed_rel}"
        target_flags = ["--target-basis", _abs(proj_dir, seed_rel), "--projection", "none"]
        # The solver writes output/collinear/sol_<seed stem>.wxf on success.
        sol_rel = f"output/collinear/sol_{Path(seed['file']).stem}.wxf"
        sol_weight = seed.get("weight")
    else:
        label = f"Solve collinear constraints on {target}"
        target_flags = ["--target", target, "--projection", projection]
        _, _, _, w = _parse_chain_target(target)
        sol_rel = None
        sol_weight = w
        if target[:3].upper() == "SEW" and w is not None and w % 2 == 0:
            # Step 6b of the named-target solver writes the MHV solution vector
            # as output/{L}loop/solMHV_{L}L.wxf with L = target_weight/2 (SEW
            # targets with an even total weight only; FEC/LEC targets keep a
            # virtual solution output).
            loops = w // 2
            sol_rel = f"output/{loops}loop/solMHV_{loops}L.wxf"

    # Declare the solution file when one is concretely written, so the run
    # engine (and exported scripts) verify it post-step.
    add_step(
        label,
        "bootstrap",
        [bootstrap, "--solve-collinear", *target_flags,
         "--rhs", rhs_arg,
         "--letter-projection", _resolve_path_arg(proj_dir, letter_proj),
         "--solver", solver,
         "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")],
        proj_dir,
        [sol_rel] if sol_rel else [],
        False,
        {"type": "solve_collinear", "custom_seed": custom},
    )
    if custom:
        provides[(nid, "solution")] = {"kind": "solution", "file": sol_rel, "weight": sol_weight, "name": Path(sol_rel).stem}
    else:
        provides[(nid, "solution")] = {"kind": "solution", "file": sol_rel, "weight": sol_weight, "name": Path(sol_rel).stem if sol_rel else target}


def _collinear_basis_flags(proj_dir, pairs, tensor_ops_bin):
    """When a seed's letter slots are compressed into a chain-basis axis
    (seed rank < rhs rank + 1, e.g. (11,76,11) vs rank-4 rhs), pass
    --basis <first_w3_basis> <first_w2_basis> so the solver expands the basis
    axis to its letter slots (11,25) and A.rank == rhs.rank + 1. When ranks
    are MEASURABLE (files from a previous run of the producing steps) the
    exact needed count is passed. When a seed/rhs file is not measurable at
    compile time — the first run schedules its producer in the same run —
    fall back to every contiguous first_w*_basis artifact: the solver's
    lenient expansion skips bases whose compressed side does not match, so
    over-passing is harmless while under-passing fails the run."""
    n_missing = None
    unmeasured = False
    for seed_rel, rhs, letter in pairs:
        sdim = _tensor_dims_or_none(tensor_ops_bin, proj_dir / seed_rel)
        if sdim is None:
            unmeasured = True
            continue
        if rhs == "0":
            continue
        rdim = _tensor_dims_or_none(tensor_ops_bin, proj_dir / rhs)
        if rdim is None:
            unmeasured = True
            continue
        miss = len(rdim) + 1 - len(sdim)
        if miss > 0 and (n_missing is None or miss > n_missing):
            n_missing = miss
    flags = []
    if n_missing is not None:
        w_range = range(n_missing + 1, 1, -1)
    elif unmeasured:
        # Optimistic: every contiguous first_w*_basis artifact, emitted
        # DESCENDING like the measured branch (the outermost compressed axis
        # carries the highest weight, so it must be contracted first). The
        # solver's lenient expansion skips the ones not needed.
        hi = 1
        while hi < 32 and (proj_dir / "output/collinear" / f"first_w{hi + 1}_basis.wxf").exists():
            hi += 1
        w_range = range(hi, 1, -1)
    else:
        return []
    for w in w_range:
        p = proj_dir / "output/collinear" / f"first_w{w}_basis.wxf"
        if not p.exists():
            break
        flags += ["--basis", str(p)]
    return flags


def _compile_solve_collinear_pairs(
    node, data, nid, seed, rhs_in, cond_in, pair_edges, wants_export,
    provides, add_step, errors, bootstrap, proj_dir, solver, incoming,
    tensor_ops_bin=None,
):
    """Multi-pair mode: every pair contributes one --pair triple; optional
    cond inputs append --pair-cond; the stacked system is solved once.
    Cond-only (seed is None) re-solves ingested [M|r] conditions."""
    # Per-pair config rows from the Inspector's Pairs list: pairs[i] = pair
    # i+1 (row 0 = the fixed seed/rhs ports). Each row carries its own letter
    # projection and optional RHS file; the wired rhs port wins over the row.
    pair_rows = _collinear_pairs(data)

    # Pair 0 comes from the fixed seed/rhs ports (seed is a concrete tensor
    # file — the dispatcher already guaranteed custom mode here).
    pairs = []
    if seed is not None:
        rhs0 = rhs_in["file"] if rhs_in is not None else pair_rows[0]["rhs"]
        if rhs_in is not None and rhs_in.get("file") is None:
            errors.append("Solve Collinear: the wired rhs must be a concrete tensor file.")
            return
        if not rhs0:
            rhs0 = "0"
        pairs.append((seed["file"], rhs0, pair_rows[0]["projection"]))

    for idx in sorted(pair_edges):
        entry = pair_edges[idx]
        ps = entry.get("seed")
        if ps is None or ps.get("file") is None:
            errors.append(f"Solve Collinear: pair {idx + 1} needs a concrete seed tensor wired to in_seed_{idx}.")
            return
        pr = entry.get("rhs")
        if pr is not None and pr.get("file") is None:
            errors.append(f"Solve Collinear: the rhs wired to in_rhs_{idx} must be a concrete tensor file.")
            return
        # Edge idx = in_seed_idx = pair idx+1 = config row idx (row 0 = pair 1).
        row = pair_rows[idx] if idx < len(pair_rows) else {"projection": "identity", "rhs": ""}
        rhs_field = row["rhs"] or "0"
        pairs.append((ps["file"], pr["file"] if pr is not None else rhs_field, row["projection"]))

    # A pair added in the Inspector's Pairs list must actually be wired —
    # silently dropping a configured constraint set would change the solve.
    for k in range(1, len(pair_rows)):
        if k not in pair_edges:
            errors.append(
                f"Solve Collinear: pair {k + 1} is configured in the Pairs list (letter projection "
                f"'{pair_rows[k]['projection']}') but its seed {k + 1} port is not wired — wire a tensor or remove the pair."
            )
            return

    cond_files = []
    for e in incoming.get(nid, []):
        if (e.get("targetHandle") or "") != "cond":
            continue
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            info = provides[src]
            if info.get("file") is None:
                errors.append("Solve Collinear: the wired cond must be a concrete matrix file.")
                return
            cond_files.append(info["file"])
    if not pairs and not cond_files:
        errors.append("Solve Collinear node: multi-pair mode needs at least one pair or cond input.")
        return

    # Output stem: --out-stem wins; else single source stem; else stem1_xN
    # (mirrors run_collinear_solver_pairs naming so provides matches the files
    # the backend writes). Cond-only always passes --out-stem (the backend
    # requires it when there are no --pair flags).
    out_stem = (data.get("out_stem") or "").strip()
    if out_stem:
        out_stem = _safe_flag_value(out_stem, "Solve Collinear out stem", errors)
        if out_stem is None:
            return
        if "/" in out_stem or not SAFE_TARGET_RE.match(out_stem):
            errors.append(f"Solve Collinear has an invalid out stem '{out_stem}' (letters, digits, _, ., - only).")
            return
    names = [Path(p[0]).stem for p in pairs] + [Path(c).stem for c in cond_files]
    stem = out_stem or (names[0] if len(names) == 1 else f"{names[0]}_x{len(names)}")
    need_stem_flag = bool(out_stem) or not pairs

    argv = [bootstrap, "--solve-collinear"]
    basis_flags = _collinear_basis_flags(proj_dir, pairs, tensor_ops_bin)
    for seed_rel, rhs, letter in pairs:
        rhs_arg = "0"
        if rhs != "0":
            rhs_arg = _safe_file_arg(proj_dir, rhs, "Solve Collinear pair RHS", errors)
            if rhs_arg is None:
                return
        if letter not in ("identity", "divergent", "finite"):
            letter = _safe_file_arg(proj_dir, letter, "Solve Collinear pair letter projection", errors)
            if letter is None:
                return
        argv += ["--pair", _abs(proj_dir, seed_rel), rhs_arg, letter]
    for c in cond_files:
        argv += ["--pair-cond", _abs(proj_dir, c)]
    argv += basis_flags
    if wants_export:
        argv += ["--export-conditions"]
    if need_stem_flag:
        argv += ["--out-stem", stem]
    argv += ["--solver", solver,
             "--data-dir", _abs(proj_dir, "data"), "--output-dir", _abs(proj_dir, "output")]

    label = f"Solve non-homogeneous constraints on {len(pairs)} pair{'s' if len(pairs) != 1 else ''}"
    if cond_files:
        label += f" + {len(cond_files)} cond"
    sol_rel = f"output/collinear/sol_{stem}.wxf"
    step_outputs = [sol_rel]
    if wants_export:
        step_outputs.append(f"output/collinear/cond_{stem}.wxf")
    # Declare the outputs the solver writes on success so the run engine (and
    # exported scripts) verify them — an exit-0 solve that wrote no solution
    # must fail loudly, not surface as a downstream file-not-found.
    add_step(
        label,
        "bootstrap",
        argv,
        proj_dir,
        step_outputs,
        False,
        {"type": "solve_collinear", "custom_seed": True, "multi_pair": True, "stem": stem},
    )
    provides[(nid, "solution")] = {"kind": "solution", "file": sol_rel, "name": Path(sol_rel).stem}
    if wants_export:
        cond_rel = f"output/collinear/cond_{stem}.wxf"
        provides[(nid, "conditions")] = {"kind": "matrix", "file": cond_rel, "name": Path(cond_rel).stem}


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
    if symmetry == "collinear":
        # The same --project run also materializes every per-weight expansion
        # basis of the chain next to the target basis (projection.hpp writes
        # first_w{w}_basis.wxf at each FEC-chain weight and last_w{w}_basis.wxf
        # at each LEC-chain weight). Expose them as extra outputs so Expand
        # Tensor nodes can wire them directly — one Projection Chain node then
        # supplies every basis the hepMHV expansion needs. Handle contract:
        # basis_w{w} = first-chain basis for FEC/SEW targets, last-chain basis
        # for LEC targets (the chain the main output rides on);
        # basis_last_w{w} = last-chain basis, only for SEW targets with L >= 2.
        _, fw, lw, _ = _parse_chain_target(target)
        kind = target[:3].upper()
        if kind == "LEC":
            firsts, lasts = 0, fw or 0
        elif kind == "SEW":
            firsts, lasts = fw or 0, lw or 0
        else:
            firsts, lasts = fw or 0, 0
        for w in range(2, firsts + 1):
            bw = f"output/collinear/first_w{w}_basis.wxf"
            provides[(node["id"], f"basis_w{w}")] = {"kind": "basis", "file": bw, "weight": None, "name": Path(bw).stem}
        for w in range(2, lasts + 1):
            bw = f"output/collinear/last_w{w}_basis.wxf"
            provides[(node["id"], f"basis_last_w{w}")] = {"kind": "basis", "file": bw, "weight": None, "name": Path(bw).stem}


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


def _compile_compute_rhs(node, incoming, provides, add_step, errors, tensor_ops_bin, tensor_add_bin, proj_dir):
    """One-block RHS generation, sealing the README-documented rules:

    The formulas below are the shuffle-algebra consequences of the master
    equations, valid for ANY loop order L (encoded up to L = 5 / weight 10;
    actually running the higher weights is a supercomputer-scale job - the
    compiled plan is still well-defined for every L).

    mhv_boundary:  master equation  1 + E = exp_⊗(Σ_k R_k)  (shuffle exp)
      Log-derivative of the exp:  E_L = (1/L)·Σ_{m=1}^{L} m·(R_m ⊗ E_{L-m})
      (single sum, E_0 = 1; proven symbolically in the shuffle word algebra
      against the graded exp at every order L <= 6). Hence
      boundary_L := E_L - R_L = (1/L)·Σ_{k=1}^{L-1} k·(R_k ⊗ E_{L-k}),
      R_1 := E1. Matches the compute_rhs binary branches at L = 2, 3, 4
      and ground-truth boundary tensors to L = 4; the C++ L = 5 branch is
      known-wrong, the recursion above is the master-equation value.

      Divergent projection (the hardcoded C++ logic, proven in the shuffle
      word algebra): the collinear solve matches only the DIVERGENT part of
      the boundary (bootstrap --solve-collinear --projection divergent; the
      divergent letters are exactly those of E1, and every R_k, k >= 2, is
      finite - compute_rhs.hpp Step 13 verifies R_L = R* carries no
      divergent letter). Shuffle products preserve letter support, so the
      finite projection P_fin is multiplicative and closes on the R's:
        P_fin(boundary_L) = (1/L)·Σ_{k=2}^{L-1} k·(R_k ⊗ F_{L-k}),
        F_m := P_fin(E_m),   1 + F = exp_⊗(Σ_{k≥2} R_k).
      Finite parts of the C++ branches: 0, 0, R2^2/2, R2⊗R3 (L = 2..5) -
      even the known-wrong L = 5 branch has the correct finite part; its
      error is purely divergent. When only the divergent projection is
      needed this pure-R content may be suppressed; the node keeps the
      exact formula (F_k is never materialized - the finite part rides
      along inside the full E_k shuffles) and the solve's letter filter
      drops it downstream.

    e47me67_te:  master equation  tE = T ⊗ (1 + E),  T = Σ_k tP_k
      Grading:  tE_L = tP_L + Σ_{k=1}^{L-1} tP_k ⊗ E_{L-k}  — the E_0 = 1
      term of (1+E) pairs with tP_L, the t-remainder, which is excluded
      here exactly as R_L is excluded from the MHV boundary. The node emits
      the sum Σ_{k=1}^{L-1} tP_k ⊗ E_{L-k} (proven symbolically against
      the graded master product at every order L <= 6).
      tP_1 = P1 := hep1LE47mE67. For k ≥ 2 the t-remainder is NOT zero —
      the weight-4 solve gives tP_2 with 393 nonzeros (its divergent part
      cancels against the boundary, as the collinear condition demands) —
      so a missing tP_k must never default to zero. Instead the recursion
        tP_k = tE_k − Σ_{j=1}^{k-1} tP_j ⊗ E_{k-j}
      is HARDCODED here: tP_k is taken from output/tP<k>.wxf when that
      cache exists, otherwise derived from output/tE<k>.wxf (the solved
      E47mE67 tensor of the weight-2k flow, expanded to letter slots;
      sol_<k> dotted with the E47−E67 seed through the collinear basis
      chain) and written back to output/tP<k>.wxf for later runs. When
      neither file exists the compile fails loudly — silently dropping the
      term produced a wrong weight-6 RHS once (the k=2 term was missing).

    The object type comes from data.mode (the Inspector select, backed by
    RHS_MODES — one registry entry per user-selectable object type). The
    mode's requires/forbids ports are validated against the actual wiring;
    lower-loop tensors E_2..E_{L-1} and R_2..R_{L-1} auto-load from output/
    (reuse_output semantics: they must exist on disk at compile time). Every
    product is a SEQUENTIAL shuffle (tensor_ops shuf) - the variant this
    codebase has verified correct; the coefficients k/L are exact fractions
    folded into the shuffle weight, so the adds are plain +1 accumulations.
    """
    nid = node["id"]
    data = node.get("data", {}) or {}
    mode_id = str(data.get("mode") or "mhv_boundary")
    spec = RHS_MODES.get(mode_id)
    if spec is None:
        errors.append(
            "Compute RHS: unknown object type '%s' - pick one of: %s."
            % (mode_id, ", ".join(RHS_MODES))
        )
        return
    weight_raw = str(data.get("weight") or "").strip()
    try:
        weight = int(weight_raw)
    except ValueError:
        errors.append("Compute RHS: choose a target weight (2, 4, 6, 8 or 10).")
        return
    if weight < 2 or weight % 2:
        errors.append("Compute RHS: weight must be an even integer >= 2 (loop order = weight/2).")
        return
    L = weight // 2

    e1 = _edge_input(provides, incoming, nid, "e1")
    p1 = _edge_input(provides, incoming, nid, "p1")
    missing = [p for p in spec["requires"] if (e1 if p == "e1" else p1) is None]
    stray = [p for p in spec["forbids"] if (e1 if p == "e1" else p1) is not None]
    if missing:
        errors.append(
            "Compute RHS (%s): wire the %s seed%s into the %s port%s."
            % (spec["label"], "/".join(missing), "s" if len(missing) > 1 else "",
               "/".join(missing), "s" if len(missing) > 1 else "")
        )
        return
    if stray:
        errors.append(
            "Compute RHS (%s) does not use the %s port%s - unwire %s or switch the object type."
            % (spec["label"], "/".join(stray), "s" if len(stray) > 1 else "",
               "/".join(stray))
        )
        return
    if e1.get("file") is None:
        errors.append("Compute RHS: the e1 seed must be a concrete tensor file.")
        return
    if p1 is not None and p1.get("file") is None:
        errors.append("Compute RHS: the p1 seed must be a concrete tensor file.")
        return

    target = _require_target(data, "Compute RHS", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    if tensor_add_bin is None:
        errors.append("The tensor_add binary was not found; build it with `make tensor_add`.")
        return

    def rel_path(name):
        return f"{_out()}{name}.wxf"

    def load_lower(basename, k):
        """Auto-load output/<basename><k>.wxf (must exist on disk at compile time)."""
        rel = f"{_out()}{basename}{k}.wxf"
        if not (proj_dir / rel).exists():
            errors.append(
                f"Compute RHS ({target}): lower-loop tensor '{rel}' not found - run the lower-loop "
                f"flow first (it produces {basename}{k})."
            )
            return None
        return {"kind": "tensor", "file": rel, "weight": None, "name": f"{basename}{k}"}

    def get_E(k):
        return e1 if k == 1 else load_lower("E", k)

    def shuf_step(a, b, w, name, desc):
        rel = rel_path(name)
        add_step(
            f"Shuffle {a['name']} (x) {b['name']} x {w} -> {name}  ({desc})",
            "tensor_ops",
            [tensor_ops_bin, "shuf", _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), w, _abs(proj_dir, rel)],
            proj_dir,
            [rel],
            True,
            {"type": "shuffle_product", "tensor_file": rel},
        )
        return rel

    def add_step_into(acc_file, acc_name, t, name):
        rel = rel_path(name)
        add_step(
            f"Add {acc_name} + {t['name']} -> {name}",
            "tensor_add",
            [tensor_add_bin, _abs(proj_dir, acc_file), _abs(proj_dir, t["file"]), "1", "1", _abs(proj_dir, rel)],
            proj_dir,
            [rel],
            True,
            {"type": "add_tensors", "tensor_file": rel},
        )
        return rel

    def load_lower_opt(basename, k):
        """Auto-load output/<basename><k>.wxf; missing file = the zero tensor."""
        rel = f"{_out()}{basename}{k}.wxf"
        if not (proj_dir / rel).exists():
            return None
        return {"kind": "tensor", "file": rel, "weight": None, "name": f"{basename}{k}"}

    if mode_id == "e47me67_te":
        if L == 1:
            provides[(nid, "out")] = {"kind": "tensor", "file": p1["file"], "weight": None,
                                      "name": "tE1", "dims": p1.get("dims")}
            return

        def sub_step(a_file, a_name, b_file, b_name, name):
            """Emit a − b (weighted tensor_add 1, −1)."""
            rel = rel_path(name)
            add_step(
                f"Subtract {a_name} - {b_name} -> {name}",
                "tensor_add",
                [tensor_add_bin, _abs(proj_dir, a_file), _abs(proj_dir, b_file), "1", "-1", _abs(proj_dir, rel)],
                proj_dir,
                [rel],
                True,
                {"type": "add_tensors", "tensor_file": rel},
            )
            return rel

        _tp_cache: dict[int, dict] = {}

        def ensure_tP(k):
            """tP_k for k >= 2: the cached output/tP<k>.wxf when present,
            otherwise DERIVED from the hardcoded remainder recursion
                tP_k = tE_k − Σ_{j=1}^{k-1} tP_j ⊗ E_{k-j},
            loading tE_k from output/tE<k>.wxf (the solved E47mE67 tensor of
            the weight-2k flow). The derived tensor is written back to
            output/tP<k>.wxf so later runs reuse it. Never defaults to zero:
            a missing term silently truncated the weight-6 RHS once."""
            if k in _tp_cache:
                return _tp_cache[k]
            rel = f"{_out()}tP{k}.wxf"
            if (proj_dir / rel).exists():
                _tp_cache[k] = {"kind": "tensor", "file": rel, "weight": None, "name": f"tP{k}"}
                return _tp_cache[k]
            te_rel = f"{_out()}tE{k}.wxf"
            if not (proj_dir / te_rel).exists():
                errors.append(
                    f"Compute RHS ({target}, {spec['label']}): neither output/tP{k}.wxf nor "
                    f"output/tE{k}.wxf exists. Produce tE{k} from the weight-{2 * k} collinear "
                    f"solve (sol dotted with the E47-E67 seed, expanded through the collinear "
                    "basis chain) — without it the tP recursion cannot close."
                )
                _tp_cache[k] = None
                return None
            # boundary_k = Σ_{j=1}^{k-1} tP_j ⊗ E_{k-j}, then tP_k = tE_k − boundary_k
            sub_terms = []
            for j in range(1, k):
                tpj = p1 if j == 1 else ensure_tP(j)
                if tpj is None:
                    _tp_cache[k] = None
                    return None
                e_right = get_E(k - j)
                if e_right is None:
                    _tp_cache[k] = None
                    return None
                sname = f"{target}_tp{k}_s{j}"
                sub_terms.append((shuf_step(tpj, e_right, "1", sname,
                                            f"tP_{k} recursion term j={j}: tP_{j} (x) E_{k-j}"), sname))
            if not sub_terms:
                errors.append(f"Compute RHS ({target}): tP_{k} recursion produced no terms.")
                _tp_cache[k] = None
                return None
            acc, acc_name = sub_terms[0]
            for i in range(1, len(sub_terms)):
                nname = f"{target}_tp{k}_p{i + 1}"
                acc = add_step_into(acc, acc_name, {"file": sub_terms[i][0], "name": sub_terms[i][1]}, nname)
                acc_name = nname
            sub_step(te_rel, f"tE{k}", acc, acc_name, f"tP{k}")
            _tp_cache[k] = {"kind": "tensor", "file": rel, "weight": None, "name": f"tP{k}"}
            return _tp_cache[k]

        terms = []
        for k in range(1, L):
            tp = p1 if k == 1 else ensure_tP(k)
            if tp is None:
                return
            e_right = get_E(L - k)
            if e_right is None:
                return
            name = target if L == 2 else f"{target}_t{k}"
            terms.append((shuf_step(tp, e_right, "1", name,
                                    f"tE_L term k={k}: tP_{k} (x) E_{L-k}"), name))
        if not terms:
            errors.append(
                "Compute RHS (%s): every tE_L term vanished - tP_1 seed is zero "
                "and no output/tP<k>.wxf files exist." % spec["label"]
            )
            return
        final_rel = terms[0][0]
        acc_name = terms[0][1]
        for i in range(1, len(terms)):
            final_rel = add_step_into(final_rel, acc_name, {"file": terms[i][0], "name": terms[i][1]},
                                      target if i == len(terms) - 1 else f"{target}_p{i + 1}")
            acc_name = target if i == len(terms) - 1 else f"{target}_p{i + 1}"
            terms[i] = (final_rel, terms[i][1])
    elif mode_id == "mhv_boundary":
        if L == 1:
            errors.append("Compute RHS: MHV weight 2 is the seed itself - the boundary starts at weight 4.")
            return
        r_map = {1: e1}
        for k in range(2, L):
            r_map[k] = load_lower("R", k)
            if r_map[k] is None:
                return
        terms = []
        for k in range(1, L):
            w = str(Fraction(k, L))
            e_right = get_E(L - k)
            if e_right is None:
                return
            name = target if L == 2 else f"{target}_t{k}"
            terms.append((shuf_step(r_map[k], e_right, w, name, f"boundary term k={k}: {w} R_{k} (x) E_{L-k}"), name))
        final_rel = terms[0][0]
        acc_name = terms[0][1]
        for i in range(1, len(terms)):
            final_rel = add_step_into(final_rel, acc_name, {"file": terms[i][0], "name": terms[i][1]},
                                      target if i == len(terms) - 1 else f"{target}_p{i + 1}")
            acc_name = target if i == len(terms) - 1 else f"{target}_p{i + 1}"
            terms[i] = (final_rel, terms[i][1])
    provides[(nid, "out")] = {"kind": "tensor", "file": final_rel, "weight": None, "name": target}


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
    letters_axes = all(t.get("letters_axes") for t in tensors) and any(t.get("letters_axes") for t in tensors)
    axis_meaning = None
    if letters_axes:
        for i, label in ((0, "second"), (1, "third")):
            ms = {t.get("axes_meaning", [None, None])[i] if t.get("axes_meaning") else t.get("axis_meaning") for t in tensors}
            if "collinear" in ms and ms - {"collinear"}:
                errors.append(
                    f"Add Tensors: the {label} entries of the operands mix collinear-basis and "
                    "non-collinear meanings — their coefficients live in different bases and cannot be summed."
                )
                return
    chain0 = t0.get("chain")
    if not all((t.get("chain") == chain0) for t in tensors[1:]):
        chain0 = None
    provides[(nid, "out")] = {"kind": t0["kind"], "file": rel_path(target), "weight": t0.get("weight"), "name": target, "dims": t0.get("dims"), "letters_axes": letters_axes, "axes_meaning": axis_meaning if letters_axes else None, "chain": chain0}


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
    if tensor is None:
        errors.append("A Ternary Contract node is missing its tensor input.")
        return
    if tensor.get("file") is None:
        errors.append("Ternary Contract: the tensor input must be a concrete rank-3 tensor file.")
        return
    # trans1/trans2 are OPTIONAL: an unwired trans port defaults to the identity
    # (that entry does not transform); a wired port may be a matrix or any
    # rank>=2 tensor (first axis contracted, remaining axes spliced in there).
    if m1 is not None and m1.get("file") is None:
        errors.append("Ternary Contract: trans1 must be a concrete matrix/tensor file.")
        return
    if m2 is not None and m2.get("file") is None:
        errors.append("Ternary Contract: trans2 must be a concrete matrix/tensor file.")
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
        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
         _abs(proj_dir, m1["file"]) if m1 is not None else "I",
         _abs(proj_dir, m2["file"]) if m2 is not None else "I",
         _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "ternary_contract", "tensor_file": rel},
    )
    chain_kind = tensor.get("kind") in ("fec1", "lec1", "fec", "lec", "sew")
    letters_axes = tensor.get("letters_axes") or chain_kind
    # The two trailing entries are INDEPENDENT bases: each one's meaning follows
    # the matrix that was applied to it (m1 on the second entry, m2 on the third).
    # Chain-kind inputs get per-axis meanings from their structure (FEC: third
    # entry is the alphabet, second is chain basis; LEC mirrored; SEW none).
    if tensor.get("letters_axes") and tensor.get("axes_meaning"):
        axes_meaning = [
            ((m1.get("meaning") if m1 is not None else None) or tensor["axes_meaning"][0]),
            ((m2.get("meaning") if m2 is not None else None) or tensor["axes_meaning"][1]),
        ]
    else:
        # After a ternary contraction both trailing axes are whatever m1/m2
        # project them to; only the matrices' own meanings are trustworthy.
        axes_meaning = [(m1.get("meaning") if m1 is not None else None),
                        (m2.get("meaning") if m2 is not None else None)]
    out_chain = None
    if chain_kind and m1 is None:
        # Keep the upstream chain provenance so downstream nodes (Symmetry Solve,
        # Apply Projection behind a custom-block cb_in port) can still derive
        # induced maps for the compressed chain-basis axis (axis 1 keeps the
        # basis when trans1 is unwired/identity; a wired trans1 consumed it).
        tensor_edge = next((e for e in incoming.get(nid, [])
                            if (e.get("targetHandle") or "") == "tensor"), None)
        if tensor_edge is not None:
            up_chain, _uerr = _chain_upstream(provides, incoming, tensor_edge["source"],
                                               tensor_edge.get("sourceHandle") or "out", {}, [])
            if up_chain:
                head = up_chain[-1]
                if head["kind"] == "sew" and not (head.get("first") or head.get("last")):
                    # The sew provides entry carries no sub-chains; walk the sew
                    # node's fec/lec inputs ourselves (same shape as _chain_upstream).
                    head = dict(head)
                    first = last = None
                    for e in incoming.get(tensor_edge["source"], []):
                        th = e.get("targetHandle") or ""
                        if th not in ("fec", "lec"):
                            continue
                        s = (e.get("source"), e.get("sourceHandle"))
                        if s not in provides:
                            continue
                        sub, _e2 = _chain_upstream(provides, incoming, e["source"], e["sourceHandle"], {}, [])
                        if th == "fec":
                            first = sub
                        else:
                            last = sub
                    if first is not None or last is not None:
                        head["first"], head["last"] = first or [], last or []
                if head["kind"] == "sew":
                    out_chain = {"kind": "sew", "file": head["file"], "weight": head.get("weight"), "name": head.get("name"),
                                 "first": [dict(x) for x in head.get("first") or []],
                                 "last": [dict(x) for x in head.get("last") or []]}
                else:
                    out_chain = {"kind": head["kind"], "file": head["file"], "weight": head.get("weight"), "name": head.get("name"),
                                 "first": up_chain if head["kind"].startswith("fec") else [],
                                 "last": up_chain if head["kind"].startswith("lec") else []}
        if out_chain is None and tensor.get("kind") in ("fec1", "lec1", "fec", "lec", "sew"):
            out_chain = tensor.get("chain")
    provides[(nid, "out")] = {"kind": tensor["kind"], "file": rel, "weight": tensor.get("weight"),
                              "name": target, "letters_axes": bool(letters_axes),
                              "axes_meaning": axes_meaning, "chain": out_chain}


def _chain_upstream(provides, incoming, nid, handle, nodes, edges, depth=0):
    """Walk a chain tensor input (fec/lec/sew handle) upstream through extend/sew
    nodes, collecting every chain tensor from the seed (weight 1) up to the given
    one. Returns (chain, error): chain is a list of provides-like dicts
    [{kind: 'fec1'|'lec1'|'fec'|'lec', file, weight}, ...] in weight order,
    or [sew_info] with .first/.last sub-chains attached for kind 'sew'."""
    if depth > 64:
        return None, "chain walk hit the recursion limit (cycle?)"
    info = provides.get((nid, handle))
    if info is None:
        for e in incoming.get(nid, []):
            if (e.get("targetHandle") or "").startswith(handle):
                s = (e.get("source"), e.get("sourceHandle"))
                if s in provides:
                    info = provides[s]
                    break
    if info is None:
        return None, None
    kind = info.get("kind")
    if kind in ("fec1", "lec1"):
        return [dict(info)], None
    if kind in ("fec", "lec"):
        node = nodes.get(nid)
        upstream = None
        for e in incoming.get(nid, []):
            th = e.get("targetHandle") or ""
            if th.startswith("fec") or th.startswith("lec") or th.startswith("seed"):
                s = (e.get("source"), e.get("sourceHandle"))
                if s in provides:
                    upstream = _chain_upstream(provides, incoming, e["source"], e["sourceHandle"] or "out", nodes, edges, depth + 1)
                    break
        if upstream is None:
            return None, None
        sub, err = upstream
        if err or sub is None:
            return None, err
        return sub + [dict(info)], None
    if kind == "sew":
        first = last = None
        for e in incoming.get(nid, []):
            th = e.get("targetHandle") or ""
            if th not in ("fec", "lec"):
                continue
            s = (e.get("source"), e.get("sourceHandle"))
            if s not in provides:
                continue
            sub, err = _chain_upstream(provides, incoming, e["source"], e["sourceHandle"], nodes, edges, depth + 1)
            if err:
                return None, err
            if th == "fec":
                first = sub
            else:
                last = sub
        if first is None or last is None:
            return None, None
        entry = dict(info)
        entry["first"], entry["last"] = first, last
        return [entry], None
    return None, None


def _tensor_dims_or_none(tensor_ops_bin, path):
    if tensor_ops_bin is None or path is None or not path.exists():
        return None
    import subprocess as _sp
    r = _sp.run([tensor_ops_bin, "dims", str(path)], capture_output=True, text=True)
    m = re.search(r"rank (\d+)((?: \d+)+) nnz", r.stdout or "")
    return [int(x) for x in m.group(2).split()] if m else None


def _symmetry_run_name(stem):
    """Map a rep-matrix stem to (bootstrap symmetry name or None, stem usable as a custom name)."""
    known = {"cycrepmat": "cyclic", "fliprepmat": "flip", "parityrepmat": "parity", "colmat42": "collinear"}
    return known.get(stem), SAFE_TARGET_RE.match(stem or "") is not None


def _expand_chain_provenance(provides, incoming, nid, handle, nodes, edges):
    """Inflate a stripped chain-basis provenance entry back into a chain-walk head.

    Tensor metadata carries {'chain': {...}, 'trailing_dims': [...]} through
    custom-block cb_in ports. Returns (head, sub) shaped exactly like the end of
    _chain_upstream, or (None, err-string) when the chain files are missing."""
    info = provides.get((nid, handle))
    if info is None:
        return None, None
    ch = info.get("chain") or {}
    kind = ch.get("kind")
    if kind not in ("fec", "lec", "sew", "fec1", "lec1"):
        return None, None

    def _mk(entry):
        return {"kind": entry["kind"], "file": entry["file"], "weight": entry["weight"],
                "name": entry.get("name")}

    first = [_mk(e) for e in (ch.get("first") or [])]
    last = [_mk(e) for e in (ch.get("last") or [])]
    if not first and not last:
        return None, None
    head = {"kind": kind, "file": ch.get("file") or (first + last)[-1]["file"],
            "weight": ch.get("weight"), "name": ch.get("name")}
    if kind == "sew":
        head["first"], head["last"] = first, last
    sub = {"first": first, "last": last}
    return head, sub


def _derive_basis_axis_maps(head, sub, sym, n, sig, add_step, proj_dir, bootstrap, tensor_ops_bin, file_steps, errors):
    """Derive the induced maps for a chain-BASIS axis (weight-w FEC/LEC basis)
    acting with symmetry matrix `sym`. Returns the SQUARE map path for the
    basis axis (stays compressed), or None on error. For a collinear matrix
    with n > 1 it errors (projections apply once)."""
    stem = Path(sym["file"]).stem
    sym_name, usable = _symmetry_run_name(stem)
    if sym_name is None and not usable:
        errors.append(f"cannot use '{stem}' as a symmetry name (unsafe characters).")
        return None
    if sym_name is None:
        sym_name = stem
    need_bootstrap = head["kind"] == "sew" or any(e.get("weight", 1) > 1 for e in (sub.get("first") or []) + (sub.get("last") or []))
    # For a sew head the compressed basis axis is the FEC side by the NMHV
    # convention (SEW, FEC_w, letters): derive the FEC chain only, so the
    # (possibly variant) LEC side is never staged (cyclic/flip/parity have no
    # action on a variant LEC basis; the full SEW derive would fail).
    maps = _emit_derived_projection(head, sub, sig, sym_name, sym["file"],
                                    add_step, proj_dir, bootstrap, file_steps, errors,
                                    needed=need_bootstrap,
                                    target_kind="first" if head["kind"] == "sew" else "auto")
    if maps is None:
        return None
    if head["kind"] == "sew":
        fw = ((head.get("first") or [{}])[-1].get("weight") if head.get("first") else 1) or 1
        cand = maps.get("first_w")
    else:
        fw = head.get("weight") or 1
        cand = maps.get("first_w") if head["kind"].startswith("fec") else maps.get("last_w")
    if cand is None:
        errors.append("could not derive the induced map for the chain head.")
        return None
    if n > 1:
        cand = _matrix_n_power(cand, n, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors)
    return cand


def _derived_matrix_sig(proj_dir: Path, rel: str) -> str:
    p = proj_dir / rel
    try:
        h = hashlib.sha1()
        h.update(rel.encode())
        h.update(str(p.stat().st_size).encode())
        h.update(str(int(p.stat().st_mtime)).encode())
        return h.hexdigest()[:10]
    except OSError:
        return "missing"


def _matrix_n_power(mat_rel, n, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors):
    """Return the rel path of M^n (cached under output/.derived/pow_<sig>/),
    or mat_rel unchanged when n == 1. Uses `tensor_ops power`."""
    if n == 1:
        return mat_rel
    stem = Path(mat_rel).stem
    pow_rel = f"output/.derived/pow_{sig}/{stem}_p{n}.wxf"
    if pow_rel not in file_steps:
        file_steps.add(pow_rel)
        (proj_dir / f"output/.derived/pow_{sig}").mkdir(parents=True, exist_ok=True)
        add_step(
            f"Raise {stem} to power {n} (shared cache)",
            "tensor_ops",
            [tensor_ops_bin, "power", _abs(proj_dir, mat_rel), str(n), _abs(proj_dir, pow_rel)],
            proj_dir,
            [pow_rel],
            True,
            {"type": "matrix_power_derived", "matrix": stem, "n": n},
        )
    return pow_rel


def _seed_composed_map(seed_file, sym_file, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors, n=1):
    """Compose E·S where E is the seed's embedding matrix ({stem}_proj.wxf, e.g. 7x42)
    and S the full letter symmetry matrix (42x42). The result (7x42) applies directly
    to an axis living in the short seed basis and lifts it to the full alphabet.
    Cached under output/.derived/composed_<sig>/ (shared across blocks/flows)."""
    stem = Path(seed_file).stem
    proj_rel = f"data/{stem}_proj.wxf"
    if not (proj_dir / proj_rel).exists() and proj_rel not in file_steps:
        if tensor_ops_bin is None:
            errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
            return None
        file_steps.add(proj_rel)
        add_step(
            f"Derive projection map {stem} -> {Path(proj_rel).stem} (drop size-1 axis)",
            "tensor_ops",
            [tensor_ops_bin, "squeeze", _abs(proj_dir, seed_file), _abs(proj_dir, proj_rel)],
            proj_dir,
            [proj_rel],
            True,
            {"type": "squeeze", "tensor_file": proj_rel},
        )
    sym_stem = Path(sym_file).stem
    sym_use = _matrix_n_power(sym_file, n, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors)
    suffix = "" if n == 1 else f"_p{n}"
    comp_rel = f"output/.derived/composed_{sig}/{stem}_x_{sym_stem}{suffix}.wxf"
    if comp_rel not in file_steps:
        if tensor_ops_bin is None:
            errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
            return None
        file_steps.add(comp_rel)
        (proj_dir / f"output/.derived/composed_{sig}").mkdir(parents=True, exist_ok=True)
        add_step(
            f"Compose seed map {stem}_proj · {sym_stem}{suffix} (shared cache)",
            "tensor_ops",
            [tensor_ops_bin, "matmul", _abs(proj_dir, proj_rel), _abs(proj_dir, sym_use),
             _abs(proj_dir, comp_rel)],
            proj_dir,
            [comp_rel],
            True,
            {"type": "seed_composed_map", "seed": stem, "symmetry": sym_stem},
        )
    return comp_rel


def _emit_derived_projection(head, sub, sig, symmetry_name, sym_file, add_step, proj_dir, bootstrap, file_steps, errors, needed=True, target_kind="auto"):
    """Stage a hidden scratch workspace (output/.derived/proj_<sig>/) holding the
    chain tensors and the symmetry rep matrix under their canonical names, then
    emit one `bootstrap --project` step that derives the induced per-weight
    transformation matrices (first_w*/last_w*, or SEW_FpL) exactly like the
    original projection pipeline — including the weight-1 special case and
    degeneracy handling. Returns the rel path of the map for the chain head,
    or None on error. Files under .derived/ are shared across blocks: the
    step is skipped when the .sig fingerprint already matches."""
    first = sub.get("first") or []
    last = sub.get("last") or []
    if not first and not last:
        errors.append("Apply Projection: empty chain — auto mode needs Extend/Sew-built tensors.")
        return None
    fw = (first[-1].get("weight") if first else 1) or 1
    lw = (last[-1].get("weight") if last else 1) or 1
    head_kind = head["kind"]
    if head_kind == "sew" and target_kind == "auto":
        target_kind = "sew"
    if head_kind == "sew" and target_kind == "first":
        # The axis being acted on is the FEC-side basis (e.g. the 97-dim w-3
        # basis of (116,97,42)): only the FEC chain's induced map is needed.
        # Deriving the full SEW target would also stage the (variant) LEC
        # side, where e.g. the cyclic symmetry has no action (last_w1 = 0)
        # and bootstrap --project fails before writing first_w*.
        head_kind = "fec"
        head = dict(head)
        head["kind"] = "fec"
    if head_kind == "sew":
        target = f"SEW_{fw}p{lw}"
    elif head_kind.startswith("fec"):
        target = f"FEC_{fw}"
    else:
        target = f"LEC_{lw}"

    scratch = f"output/.derived/proj_{sig}_{target}"
    sdata = proj_dir / scratch / "data"
    sout = proj_dir / scratch / "output"
    def _link(src_rel: str, dst: Path):
        src = proj_dir / src_rel
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        try:
            dst.symlink_to(src.resolve())
        except OSError:
            import shutil
            shutil.copy2(src, dst)
    try:
        sdata.mkdir(parents=True, exist_ok=True)
        sout.mkdir(parents=True, exist_ok=True)
        # Symmetry rep matrix under its expected data/<stem>.wxf name.
        _link(sym_file, sdata / Path(sym_file).name)
        # Chain tensors under their canonical names (FEC_w/LEC_w/SEW_FpL).
        # bootstrap --project walks the WHOLE chain (reads FEC_2..FEC_w when
        # projecting FEC_w), so stage the intermediate weights too — the walk
        # only reports the seed and the head. Dangling symlinks are fine: the
        # producing extend steps run before the derive step.
        for link in first + last:
            w = link.get("weight")
            if w is None:
                continue
            base = "FEC" if link["kind"].startswith("fec") else "LEC"
            if w == 1:
                _link(link["file"], sdata / f"{base}_1.wxf")
            else:
                tailname = f"{base}_{w}.wxf"
                prefix = link["file"][: len(link["file"]) - len(tailname)] if link["file"].endswith(tailname) else "output/"
                for k in range(2, w + 1):
                    rel_k = link["file"] if k == w else f"{prefix}{base}_{k}.wxf"
                    _link(rel_k, sout / f"{base}_{k}.wxf")
        # bootstrap --project reads BOTH weight-1 seeds at the seed level, even
        # for a single-sided chain — stage the opposite side from data/ if absent.
        for base in ("FEC", "LEC"):
            if not (sdata / f"{base}_1.wxf").exists():
                seed_rel = f"data/{base}_1.wxf"
                if (proj_dir / seed_rel).exists():
                    _link(seed_rel, sdata / f"{base}_1.wxf")
        if head_kind == "sew":
            _link(head["file"], sout / f"{target}.wxf")
    except OSError as e:
        errors.append(f"Apply Projection: cannot stage the derived-matrix scratch workspace: {e}")
        return None

    summary = f"{scratch}/output/{symmetry_name}/summary.txt"
    maps = {
        "_root": scratch,
        "first_w": f"{scratch}/output/{symmetry_name}/first_w{fw}.wxf" if first else None,
        "last_w": f"{scratch}/output/{symmetry_name}/last_w{lw}.wxf" if last else None,
        "first_w_prev": f"{scratch}/output/{symmetry_name}/first_w{max(fw - 1, 1)}.wxf" if first else None,
        "last_w_prev": f"{scratch}/output/{symmetry_name}/last_w{max(lw - 1, 1)}.wxf" if last else None,
    }
    if needed and summary not in file_steps:
        outs = [summary] + [m for m in maps.values() if m]
        add_step(
            f"Derive induced {symmetry_name} maps for {target} (shared cache)",
            "bootstrap",
            [bootstrap, "--project", "--symmetry", symmetry_name, "--target", target,
             "--data-dir", _abs(proj_dir, f"{scratch}/data"),
             "--output-dir", _abs(proj_dir, f"{scratch}/output")],
            proj_dir,
            outs,
            True,
            {"type": "derived_projection", "symmetry": symmetry_name, "target": target},
        )
        file_steps.add(summary)
    return maps


def _compile_apply_symmetry(node, incoming, provides, add_step, errors, bootstrap, tensor_ops_bin, proj_dir, nodes, edges, file_steps):
    nid = node["id"]
    data = node.get("data", {}) or {}
    tensor = _edge_input(provides, incoming, nid, "tensor")
    if tensor is None:
        errors.append("An Apply Projection node is missing its tensor input.")
        return
    sym = _edge_input(provides, incoming, nid, "sym")
    m1 = _edge_input(provides, incoming, nid, "trans1")
    m2 = _edge_input(provides, incoming, nid, "trans2")

    if sym is not None and (m1 is not None or m2 is not None):
        errors.append("Apply Projection: wire either the 'sym' input (auto-derived chain matrices) or trans1/trans2 (manual), not both.")
        return

    if sym is None:
        _compile_ternary_contract(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir)
        return

    if tensor.get("file") is None or sym.get("file") is None:
        errors.append("Apply Projection: the tensor and sym inputs must be concrete tensor/matrix files.")
        return
    if bootstrap is None or tensor_ops_bin is None:
        errors.append("The bootstrap/tensor_ops binaries were not found; build them with `make bootstrap tensor_ops`.")
        return

    sig = _derived_matrix_sig(proj_dir, sym["file"])

    raw_n = data.get("n")
    if raw_n in (None, ""):
        n = 1
    else:
        try:
            n = int(raw_n)
        except (TypeError, ValueError):
            errors.append("Apply Projection: n must be a positive integer.")
            return
    if n < 1:
        errors.append("Apply Projection: n must be a positive integer.")
        return
    if n > 1 and sym.get("meaning") == "collinear":
        errors.append(
            "Apply Projection: n > 1 is only valid for symmetry matrices; a collinear "
            "projection (basis change) can only be applied once (n = 1)."
        )
        return

    tensor_edge = None
    for e in incoming.get(nid, []):
        if (e.get("targetHandle") or "") == "tensor":
            tensor_edge = e
            break
    if tensor_edge is None:
        errors.append("Apply Projection: the tensor input must be wired (auto mode needs the upstream chain in the graph).")
        return

    if tensor.get("letters_axes"):
        # The tensor's trailing axes were already projected to the uniform full
        # alphabet basis (e.g. by a manual Ternary Contract with proj matrices,
        # possibly via solve_symmetry / custom-block ports): the symmetry acts
        # as plain S^n on both axes. The meaning of the axes carries through
        # only if the matrix preserves it (symmetry reps do; a collinear
        # projection changes the meaning to the collinear limit).
        tm = tensor.get("axes_meaning") or (
            [tensor.get("axis_meaning"), tensor.get("axis_meaning")] if tensor.get("axis_meaning") else [None, None])
        t_meaning = "collinear" if "collinear" in tm else None
        s_meaning = sym.get("meaning")
        if t_meaning == "collinear":
            # Axes are in the collinear limit basis. Only a matrix acting WITHIN
            # that basis is admissible; a full-alphabet matrix (preserve rep,
            # collinear projection of the alphabet, anything else) cannot act.
            errors.append(
                f"Apply Projection: the letters axes of '{tensor['name']}' are in the collinear "
                "limit basis; a full-alphabet matrix cannot act on them (no proper matrix for these axes)."
            )
            return
        sym_file_n = _matrix_n_power(sym["file"], n, sig, add_step, proj_dir,
                                     tensor_ops_bin, file_steps, errors)
        if sym_file_n is None:
            return
        # M0 fast-path guard: a compressed chain-basis axis (e.g. (5,97,42),
        # 97 = weight-3 FEC basis) cannot take the 42x42 letter matrix. Route
        # through the induced-map derivation instead (basis axis stays
        # compressed; S^n acts on the letters axis only). The file may not
        # exist yet on the first compile (pre-run), so provenance alone
        # (letters_axes + chain) decides; measured dims refine the check when
        # available. A collinear sym on the letters axis is handled by the
        # not-chain M2 branch below, not here.
        tpath = proj_dir / tensor["file"]
        tdims0 = _tensor_dims_or_none(tensor_ops_bin, tpath) or []
        sdims0 = sym.get("dims") or []
        compressed = False
        if tensor.get("chain"):
            if len(tdims0) >= 3 and len(sdims0) == 2:
                # The letters axis (last) must match the matrix's letters side;
                # for a projection like colmat42 (42->11) that side is sdims[0],
                # for a symmetry rep (42x42) both sides are 42. The basis axis
                # (second-to-last, e.g. 97) must match NEITHER side.
                letters_ok = tdims0[-1] in (sdims0[0], sdims0[1])
                compressed = letters_ok and tdims0[-2] != sdims0[1] and tdims0[-2] != sdims0[0]
            else:
                compressed = True
        if compressed and sym.get("meaning") != "collinear":
            prov_head, prov_sub = _expand_chain_provenance(provides, incoming, tensor_edge["source"],
                                                           tensor_edge.get("sourceHandle") or "out", nodes, edges)
            if prov_head is not None:
                basis_map = _derive_basis_axis_maps(prov_head, prov_sub, sym, n, sig, add_step, proj_dir,
                                                    bootstrap, tensor_ops_bin, file_steps, errors)
                if basis_map is not None:
                    target = _require_target(data, "Apply Projection", errors, nid)
                    if target is None:
                        return
                    rel = f"{_out()}{target}.wxf"
                    add_step(
                        f"Apply symmetry {sym.get('name') or 'S'}{'' if n == 1 else f' (S^{n})'} to {tensor['name']} -> {target} (induced map on the chain-basis axis)",
                        "tensor_ops",
                        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
                         _abs(proj_dir, basis_map), _abs(proj_dir, sym_file_n), _abs(proj_dir, rel)],
                        proj_dir,
                        [rel],
                        True,
                        {"type": "apply_symmetry"},
                    )
                    d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
                    provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                              "name": target, "dims": None, "letters_axes": True,
                                              "axes_meaning": [None, s_meaning], "chain": tensor.get("chain")}
                    return
        elif compressed:
            # Collinear projection on a compressed chain-basis tensor: derive
            # the induced COLLINER map for the basis axis (first_w of the
            # chain with the collinear symmetry) and pair it with the letter
            # matrix — mirrors the M2 branch below for the letters_axes path.
            prov_head, prov_sub = _expand_chain_provenance(provides, incoming, tensor_edge["source"],
                                                           tensor_edge.get("sourceHandle") or "out", nodes, edges)
            if prov_head is not None:
                cmap = _derive_basis_axis_maps(prov_head, prov_sub, sym, 1, sig, add_step, proj_dir,
                                               bootstrap, tensor_ops_bin, file_steps, errors)
                if cmap is not None:
                    target = _require_target(data, "Apply Projection", errors, nid)
                    if target is None:
                        return
                    rel = f"{_out()}{target}.wxf"
                    add_step(
                        f"Apply collinear projection {sym.get('name') or 'S'} to {tensor['name']} -> {target} (induced collinear map on the chain-basis axis)",
                        "tensor_ops",
                        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
                         _abs(proj_dir, cmap), _abs(proj_dir, sym["file"]), _abs(proj_dir, rel)],
                        proj_dir,
                        [rel],
                        True,
                        {"type": "apply_symmetry"},
                    )
                    d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
                    provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                              "name": target, "dims": None, "letters_axes": True,
                                              "axes_meaning": ["collinear", "collinear"], "chain": None}
                    return
        target = _require_target(data, "Apply Projection", errors, nid)
        if target is None:
            return
        rel = f"{_out()}{target}.wxf"
        add_step(
            f"Apply symmetry {sym.get('name') or 'S'}{'' if n == 1 else f' (S^{n})'} to {tensor['name']} -> {target} (uniform alphabet axes)",
            "tensor_ops",
            [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
             _abs(proj_dir, sym_file_n), _abs(proj_dir, sym_file_n), _abs(proj_dir, rel)],
            proj_dir,
            [rel],
            True,
            {"type": "apply_symmetry"},
        )
        d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
        provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                  "name": target, "dims": None, "letters_axes": True,
                                  "axes_meaning": [s_meaning, s_meaning]}
        return

    chain, err = _chain_upstream(provides, incoming, tensor_edge["source"], tensor_edge.get("sourceHandle") or "out", nodes, edges)
    if err:
        errors.append(f"Apply Projection: {err}")
        return
    if not chain:
        # Old-way compatibility fallback: for an arbitrary tensor, check the
        # trailing-axis dimensions directly and apply S^n if they match the
        # symmetry matrix dimension (uniform full-alphabet axes).
        tpath = proj_dir / tensor["file"]
        if not tpath.exists():
            errors.append(
                f"Apply Projection: tensor file '{tensor['file']}' not found; auto mode needs the "
                "tensor from a chain built in this graph, or an existing file whose trailing axes "
                "match the symmetry matrix dimension."
            )
            return
        import subprocess as _sp
        r = _sp.run([tensor_ops_bin, "dims", str(tpath)], capture_output=True, text=True)
        m = re.search(r"rank (\d+)((?: \d+)+) nnz", r.stdout or "")
        tdims = [int(x) for x in m.group(2).split()] if m else []
        sym_dims = sym.get("dims") or []
        if not sym_dims:
            spath = proj_dir / sym["file"]
            if spath.exists():
                r2 = _sp.run([tensor_ops_bin, "dims", str(spath)], capture_output=True, text=True)
                m2 = re.search(r"rank (\d+)((?: \d+)+) nnz", r2.stdout or "")
                if m2:
                    sym_dims = [int(x) for x in m2.group(2).split()][-2:]
        if len(tdims) >= 3 and len(sym_dims) == 2 and tdims[-1] == sym_dims[0] and tdims[-2] == sym_dims[0]:
            if n > 1 and sym_dims[0] != sym_dims[1]:
                errors.append(
                    f"Apply Projection: n > 1 needs a square matrix, but '{sym['file']}' is "
                    f"{sym_dims[0]}x{sym_dims[1]} (a projection can only be applied once, n = 1)."
                )
                return
            sym_file_n = _matrix_n_power(sym["file"], n, sig, add_step, proj_dir,
                                         tensor_ops_bin, file_steps, errors)
            if sym_file_n is None:
                return
            target = _require_target(data, "Apply Projection", errors, nid)
            if target is None:
                return
            rel = f"{_out()}{target}.wxf"
            add_step(
                f"Apply symmetry {sym.get('name') or 'S'}{'' if n == 1 else f' (S^{n})'} to {tensor['name']} -> {target} (dims-matched axes)",
                "tensor_ops",
                [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
                 _abs(proj_dir, sym_file_n), _abs(proj_dir, sym_file_n), _abs(proj_dir, rel)],
                proj_dir,
                [rel],
                True,
                {"type": "apply_symmetry"},
            )
            d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
            provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                      "name": target, "dims": None, "letters_axes": True,
                                      "axes_meaning": [sym.get("meaning"), sym.get("meaning")]}
            return
        # M0: the tensor's trailing axes are (chain-basis dim B, letter dim 42)
        # from a compressed NMHV-style ternary — act with the INDUCED map on
        # the basis axis (stays compressed) and S^n on the letters axis.
        prov_head, prov_sub = _expand_chain_provenance(provides, incoming, tensor_edge["source"],
                                                       tensor_edge.get("sourceHandle") or "out", nodes, edges)
        if prov_head is not None and (sym.get("meaning") != "collinear"):
            basis_map = _derive_basis_axis_maps(prov_head, prov_sub, sym, n, sig, add_step, proj_dir,
                                                bootstrap, tensor_ops_bin, file_steps, errors)
            if basis_map is None:
                return
            sym_file_n = _matrix_n_power(sym["file"], n, sig, add_step, proj_dir,
                                         tensor_ops_bin, file_steps, errors)
            if sym_file_n is None:
                return
            # (B, 42): basis map on axis-2, S^n on axis-3; (42, B): mirrored.
            m1_arg, m2_arg = basis_map, sym_file_n
            if len(tdims) >= 3 and tdims[-2] == sym_dims[0] and tdims[-1] != sym_dims[0]:
                m1_arg, m2_arg = sym_file_n, basis_map
            target = _require_target(data, "Apply Projection", errors, nid)
            if target is None:
                return
            rel = f"{_out()}{target}.wxf"
            add_step(
                f"Apply symmetry {sym.get('name') or 'S'}{'' if n == 1 else f' (S^{n})'} to {tensor['name']} -> {target} (induced map on the chain-basis axis)",
                "tensor_ops",
                [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
                 _abs(proj_dir, m1_arg), _abs(proj_dir, m2_arg), _abs(proj_dir, rel)],
                proj_dir,
                [rel],
                True,
                {"type": "apply_symmetry"},
            )
            d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
            prov_chain = dict(tensor.get("chain") or {})
            provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                      "name": target, "dims": None, "letters_axes": True,
                                      "axes_meaning": [None, sym.get("meaning")], "chain": prov_chain}
            return
        # M2: collinear block — first trailing axis is a chain basis that must
        # move to the collinear limit (76 = collinear rank of the w-3 basis):
        # bootstrap --project --symmetry collinear gives first_w3 (97x76).
        if prov_head is not None and sym.get("meaning") == "collinear":
            cmap = _derive_basis_axis_maps(prov_head, prov_sub, sym, 1, sig, add_step, proj_dir,
                                           bootstrap, tensor_ops_bin, file_steps, errors)
            if cmap is None:
                return
            target = _require_target(data, "Apply Projection", errors, nid)
            if target is None:
                return
            rel = f"{_out()}{target}.wxf"
            add_step(
                f"Apply collinear projection to {tensor['name']} -> {target} (chain-basis axis via induced collinear map)",
                "tensor_ops",
                [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
                 _abs(proj_dir, cmap), _abs(proj_dir, sym["file"]), _abs(proj_dir, rel)],
                proj_dir,
                [rel],
                True,
                {"type": "apply_symmetry"},
            )
            d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
            provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"),
                                      "name": target, "dims": None, "letters_axes": True,
                                      "axes_meaning": ["collinear", "collinear"], "chain": None}
            return
        got = f"{tdims[-2:]}" if len(tdims) >= 2 else "unknown"
        errors.append(
            f"Apply Projection: no proper symmetry matrix found for the trailing axes of '{tensor['file']}' "
            f"(axes {got} vs symmetry matrix {sym_dims or 'unknown'}). Wire trans1/trans2 manually, or feed "
            "a tensor from a chain built in this graph."
        )
        return

    head = chain[-1]
    fw = head.get("weight") or 1
    if head["kind"] == "sew":
        fw = ((head.get("first") or [{}])[-1].get("weight") if head.get("first") else 1) or 1
        lw = ((head.get("last") or [{}])[-1].get("weight") if head.get("last") else 1) or 1
    else:
        lw = 1
    if head["kind"] in ("fec1", "fec", "lec1", "lec"):
        sub = {"first": chain if head["kind"].startswith("fec") else [],
               "last": chain if head["kind"].startswith("lec") else []}
        if head["kind"].startswith("lec"):
            head = dict(head)
            head["kind"] = "lec"
    else:  # sew: sub-chains attached by the walker
        sub = {"first": head.get("first") or [], "last": head.get("last") or []}

    # Name the symmetry from the rep-matrix file (cycrepmat -> cyclic etc.).
    stem = Path(sym["file"]).stem
    known = {"cycrepmat": "cyclic", "fliprepmat": "flip", "parityrepmat": "parity",
             "colmat42": "collinear"}
    symmetry_name = known.get(stem)
    if symmetry_name is None:
        # Custom symmetry: bootstrap accepts data/<name>.wxf, name must be the stem.
        if not SAFE_TARGET_RE.match(stem or ""):
            errors.append(f"Apply Projection: cannot use '{stem}' as a symmetry name (unsafe characters).")
            return
        symmetry_name = stem

    # Does any trailing axis need a weight >= 2 induced map (bootstrap), or is
    # every trailing axis either the letters axis (full 42) or a weight-1 seed
    # basis (handled by the cheaper composed E·S map)?
    need_bootstrap = head["kind"] == "sew" or fw > 1 or lw > 1
    maps = _emit_derived_projection(head, sub, sig, symmetry_name, sym["file"],
                                    add_step, proj_dir, bootstrap, file_steps, errors,
                                    needed=need_bootstrap)
    if maps is None:
        return

    def _seed_axis_map(seed_entry):
        """Map for an axis living in a weight-1 seed basis: compose the seed
        embedding E ({stem}_proj) with the full symmetry matrix S, so the short
        axis is projected to the uniform full-alphabet dimension."""
        if seed_entry is None:
            return None
        return _seed_composed_map(seed_entry["file"], sym["file"], sig, add_step, proj_dir,
                                  tensor_ops_bin, file_steps, errors, n=n)

    def _pow(mat_arg):
        """Raise a derived square map to the n-th power (shared cache)."""
        if mat_arg is None or mat_arg == "I" or n == 1:
            return mat_arg
        return _matrix_n_power(mat_arg, n, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors)

    sym_file_n = _matrix_n_power(sym["file"], n, sig, add_step, proj_dir, tensor_ops_bin, file_steps, errors)

    if head["kind"] == "sew":
        # SEW dims (sew_basis, FEC_F_basis, LEC_L_basis): a weight-1 side has the
        # raw seed basis (e.g. 7 / 14) on that axis — lift it to the uniform full
        # alphabet dimension via the composed E·S map instead of a bootstrap map.
        m1_arg = _pow(maps.get("first_w")) if fw > 1 else _seed_axis_map((sub.get("first") or [None])[0])
        m2_arg = _pow(maps.get("last_w")) if lw > 1 else _seed_axis_map((sub.get("last") or [None])[0])
    elif head["kind"].startswith("fec"):
        # FEC_w dims (basis_w, basis_{w-1}, 42): axis-2 uses the letter map S,
        # axis-1 uses the induced map on the weight-(w-1) basis ('I' at w=1,
        # where the axis is literally dimension 1).
        m1_arg = "I" if fw == 1 else _pow(maps.get("first_w_prev"))
        m2_arg = sym_file_n
    else:
        # LEC_w dims (basis_w, 42, basis_{w-1}): axis-1 is the letter map S,
        # axis-2 uses the induced map ('I' at w=1, dimension-1 axis).
        m1_arg = sym_file_n
        m2_arg = "I" if lw == 1 else _pow(maps.get("last_w_prev"))
    if not m1_arg or not m2_arg:
        errors.append("Apply Projection: could not derive the induced maps for the chain head.")
        return

    target = _require_target(data, "Apply Projection", errors, nid)
    if target is None:
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Apply symmetry {sym.get('name') or 'S'}{'' if n == 1 else f' (S^{n})'} to {tensor['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "ternary", _abs(proj_dir, tensor["file"]),
         ("I" if m1_arg == "I" else _abs(proj_dir, m1_arg)),
         ("I" if m2_arg == "I" else _abs(proj_dir, m2_arg)),
         _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "apply_symmetry"},
    )
    d_kind = tensor.get("kind") if tensor.get("kind") in TENSOR_KINDS else "tensor"
    # No letters_axes tag: chain outputs have at most ONE alphabet axis (FEC/LEC
    # heads; SEW heads have none) — the other trailing axis is a chain-basis
    # axis, and the fast path applies S to BOTH trailing axes, so tagging here
    # would be wrong. A subsequent Apply Projection correctly re-walks the
    # chain or falls back to the dims check (which rejects the basis axis).
    provides[(nid, "out")] = {"kind": d_kind, "file": rel, "weight": tensor.get("weight"), "name": target, "dims": None}


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
    provides[(nid, "out")] = {"kind": "matrix", "file": rel, "weight": None, "name": target, "dims": mat.get("dims"), "meaning": mat.get("meaning")}


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
        provides[(nid, f"out_{oi}")] = {"kind": kind, "file": rel, "weight": None, "name": f"{target}_{name}", "letters_axes": elems[0].get("letters_axes"), "axes_meaning": elems[0].get("axes_meaning"), "chain": elems[0].get("chain")}


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
    rank = len(a.get("dims") or b.get("dims") or [0, 0, 0])
    trailing = axis in (-1, -2, rank, rank - 1)
    # Non-join axes must have matching dims, so both operands carry letters
    # axes or neither; require agreement to keep the claim trustworthy.
    agree = bool(a.get("letters_axes")) == bool(b.get("letters_axes"))
    letters_axes = agree and bool(a.get("letters_axes")) and not trailing
    out = {"kind": a["kind"], "file": rel, "weight": None, "name": target,
           "dims": _join_dims(a.get("dims"), b.get("dims"), axis), "letters_axes": letters_axes}
    if letters_axes:
        # Meaning inheritance: joining along the parameter axis keeps the
        # shared trailing axes, so their basis meanings carry over — but only
        # when both operands make the same claim (same guard as Add Tensors).
        for i, label in ((0, "second"), (1, "third")):
            ms = {(t.get("axes_meaning") or [None, None])[i] for t in (a, b)}
            if "collinear" in ms and ms - {"collinear"}:
                errors.append(
                    f"Join Tensors: the {label} entries of the operands mix collinear-basis and "
                    "non-collinear meanings — their dims match but their coefficients live in different bases."
                )
                return
        if a.get("axes_meaning") == b.get("axes_meaning"):
            out["axes_meaning"] = a.get("axes_meaning")
        chain_a, chain_b = a.get("chain"), b.get("chain")
        if chain_a is not None and chain_b is not None and chain_a != chain_b:
            errors.append(
                "Join Tensors: the operands carry different chain provenance — their basis axes "
                "come from different chains despite matching dims."
            )
            return
        out["chain"] = chain_a if (chain_a is not None and chain_a == chain_b) else None
    provides[(nid, "out")] = out


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
    out_meta = {"kind": "tensor", "file": rel, "weight": None, "name": target, "dims": dims_out}
    # The contraction removes one axis from A and one from B; the REMAINING
    # trailing axes of the result come from B's tail when B is the long side
    # (e.g. sol[-1]·E_pre1[1] with sol (11,345) and E_pre1 (11,97,42) ->
    # (11,97,42): the trailing (97,42) = B's chain-basis + letters axes).
    # Propagate letters_axes/chain provenance so downstream Apply Projection
    # nodes route to the induced-map branch instead of the dims-matched one
    # (which mis-fires on stale compile-time dims).
    if b.get("letters_axes"):
        out_meta["letters_axes"] = True
        out_meta["axes_meaning"] = b.get("axes_meaning")
        out_meta["chain"] = b.get("chain")
        out_meta["dims"] = out_meta.get("dims") or b.get("dims")
    elif a.get("letters_axes"):
        out_meta["letters_axes"] = True
        out_meta["axes_meaning"] = a.get("axes_meaning")
        out_meta["chain"] = a.get("chain")
    provides[(nid, "out")] = out_meta


def _compile_squeeze_tensor(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    src = _edge_input(provides, incoming, nid, "in")
    if src is None:
        errors.append("A Squeeze Tensor node needs an input.")
        return
    if src.get("file") is None:
        errors.append("Squeeze Tensor: the input must be a concrete tensor file (virtual outputs cannot be squeezed).")
        return
    dims_in = src.get("dims")
    if dims_in:
        dims_out = _squeeze_dims(dims_in)
        if dims_out is None or len(dims_out) < 2:
            errors.append(
                f"Squeeze Tensor: input dims {dims_in} would squeeze to "
                f"{'×'.join(str(d) for d in dims_out) if dims_out else 'a scalar/vector'}; "
                "the result must keep at least 2 axes (WXF cannot store rank-1 tensors)."
            )
            return
    else:
        dims_out = None
    target = _require_target(data, "Squeeze Tensor", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Squeeze {src['name']} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "squeeze", _abs(proj_dir, src["file"]), _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "squeeze_tensor", "tensor_file": rel},
    )
    provides[(nid, "out")] = {"kind": "tensor", "file": rel, "weight": None, "name": target, "dims": dims_out}


def _compile_shuffle_product(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    nid = node["id"]
    data = node.get("data", {}) or {}
    a = _edge_input(provides, incoming, nid, "a")
    b = _edge_input(provides, incoming, nid, "b")
    if a is None or b is None:
        errors.append("A Shuffle Product node needs both A and B inputs.")
        return
    if a.get("file") is None or b.get("file") is None:
        errors.append("Shuffle Product: both inputs must be concrete tensor files (the basis/solution/boundary outputs of Project/Solve/Compute RHS nodes are virtual and cannot be shuffled).")
        return
    weight = str(data.get("weight") or "1").strip() or "1"
    if not RAT_RE.match(weight):
        errors.append(f"Shuffle Product: weight '{weight}' is not a rational number (examples: 1, -1, 1/2, 1/6).")
        return
    target = _require_target(data, "Shuffle Product", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Shuffle {a['name']} ⊗ {b['name']} × {weight} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "shuf", _abs(proj_dir, a["file"]), _abs(proj_dir, b["file"]), weight, _abs(proj_dir, rel)],
        proj_dir,
        [rel],
        True,
        {"type": "shuffle_product", "tensor_file": rel},
    )
    ad, bd = a.get("dims"), b.get("dims")
    dims_out = None
    if ad and bd:
        dims_out = ad + bd
    provides[(nid, "out")] = {"kind": "tensor", "file": rel, "weight": None, "name": target, "dims": dims_out}


def _compile_expand_tensor(node, incoming, provides, add_step, errors, tensor_ops_bin, proj_dir):
    # Expand a (FEC, letter) or (1, FEC, letter) hepMHV-style tensor back to the
    # full 42-letter alphabet by contracting the FEC axis with a chain of
    # rank-3 (FEC, FEC', letter) bases, highest weight first.
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
            continue
        src = (e.get("source"), e.get("sourceHandle"))
        if src in provides:
            indexed.append((idx, provides[src]))
    indexed.sort(key=lambda x: x[0])
    if not indexed:
        errors.append("An Expand Tensor node needs a tensor input and at least one basis input.")
        return
    tensor = indexed[0][1]
    bases = [info for _, info in indexed[1:]]
    if not bases:
        errors.append("An Expand Tensor node needs at least one basis input (wire in_0 to the tensor and in_1.. to the bases, highest weight first).")
        return
    if tensor.get("file") is None:
        errors.append("Expand Tensor: the tensor input must be a concrete tensor file.")
        return
    for i, b in enumerate(bases, start=1):
        if b.get("file") is None:
            errors.append(f"Expand Tensor: basis input {i} must be a concrete tensor file (virtual basis/solution outputs cannot be expanded with).")
            return
    td = tensor.get("dims")
    if td and len(td) not in (2, 3):
        errors.append(f"Expand Tensor: the tensor input has dims {td} but must be rank 2 (FEC, letter) or rank 3 (1, FEC, letter).")
        return
    if td and len(td) == 3 and td[0] != 1:
        errors.append(f"Expand Tensor: the tensor input has dims {td} but its leading axis must be 1 (the (1, FEC, letter) convention).")
        return
    dims_out = None
    if td and all(b.get("dims") for b in bases):
        # Mirror the kernel contraction chain: the tensor's FEC axis is
        # consumed by each basis (dim basis[0]) and replaced by the basis's
        # FEC' axis (basis[1]); each basis appends one letter axis (basis[2]).
        dims_out = [td[-1]]
        fec = td[-2]
        for i, b in enumerate(bases, start=1):
            bd = b["dims"]
            if len(bd) != 3:
                errors.append(f"Expand Tensor: basis input {i} ({b.get('name')}) has dims {bd} but must be a rank-3 (FEC, FEC', letter) tensor.")
                return
            if bd[0] != fec:
                errors.append(f"Expand Tensor: basis input {i} ({b.get('name')}) expects FEC dim {fec} (axis 1) but provides {bd[0]}.")
                return
            dims_out.append(bd[2])
            fec = bd[1]
        dims_out = [fec] + dims_out
    target = _require_target(data, "Expand Tensor", errors, nid)
    if target is None:
        return
    if tensor_ops_bin is None:
        errors.append("The tensor_ops binary was not found; build it with `make tensor_ops`.")
        return
    rel = f"{_out()}{target}.wxf"
    add_step(
        f"Expand {tensor['name']} with {len(bases)} basis{'' if len(bases) == 1 else 'es'} -> {target}",
        "tensor_ops",
        [tensor_ops_bin, "expand", _abs(proj_dir, tensor["file"]), _abs(proj_dir, rel),
         *(_abs(proj_dir, b["file"]) for b in bases)],
        proj_dir,
        [rel],
        True,
        {"type": "expand_tensor", "tensor_file": rel},
    )
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
