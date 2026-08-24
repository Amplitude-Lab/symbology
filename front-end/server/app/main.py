from __future__ import annotations

import json
import subprocess
import uuid
from pathlib import Path

from fastapi import Body, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse

from . import storage, templates
from .compile import compile_flow, flow_outputs_catalog
from .config import WEB_DIST, env_status, find_tensor_ops, find_wolframscript
from .jobs import engine
from .wolfram import property_display_name, property_script, property_tensor_relpath, read_result_file, summary_script

app = FastAPI(title="Symbology Front-End")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["*"],
    allow_headers=["*"])


@app.on_event("startup")
def _startup() -> None:
    storage.init_dirs()


@app.post("/api/client-log")
def api_client_log(body: dict = Body(...)) -> dict:
    from datetime import datetime, timezone
    entry = {"ts": datetime.now(timezone.utc).isoformat(), **{k: body.get(k) for k in ("kind", "message", "stack", "url", "ua")}}
    log_file = storage.PROJECTS_DIR / "_client_errors.log"
    try:
        with open(log_file, "a") as f:
            f.write(json.dumps(entry) + "\n")
    except OSError:
        pass
    return {"ok": True}


def _get_project(pid: str) -> dict:
    proj = storage.load_project(pid)
    if proj is None:
        raise HTTPException(404, f"project '{pid}' not found")
    return proj


def _get_alphabet(proj: dict, aid: str) -> dict:
    alpha = storage.find_alphabet(proj, aid)
    if alpha is None:
        raise HTTPException(404, "alphabet not found")
    return alpha


@app.get("/api/env")
def api_env() -> dict:
    return env_status()


@app.get("/api/templates")
def api_templates() -> list:
    return templates.list_templates()


@app.get("/api/projects")
def api_list_projects() -> list:
    return storage.list_projects()


@app.post("/api/projects")
def api_create_project(body: dict = Body(...)) -> dict:
    name = (body.get("name") or "").strip()
    if not name:
        raise HTTPException(400, "project name is required")
    template_id = body.get("template_id")
    if template_id is not None and template_id not in templates.TEMPLATES:
        raise HTTPException(400, f"unknown template '{template_id}'")
    proj = storage.create_project(name)
    if template_id:
        templates.apply_template(proj, template_id)
        storage.save_project(proj)
    return proj


@app.get("/api/projects/{pid}")
def api_get_project(pid: str) -> dict:
    return _get_project(pid)


@app.delete("/api/projects/{pid}")
def api_delete_project(pid: str) -> dict:
    if not storage.delete_project(pid):
        raise HTTPException(404, "project not found")
    return {"ok": True}


@app.post("/api/projects/{pid}/alphabets")
def api_create_alphabet(pid: str, body: dict = Body(...)) -> dict:
    proj = _get_project(pid)
    letters = [s.strip() for s in (body.get("letters") or []) if str(s).strip()]
    if not body.get("name") or not letters:
        raise HTTPException(400, "alphabet name and letters are required")
    if len(set(letters)) != len(letters):
        raise HTTPException(400, "letters must be unique")
    if any(a["name"] == body["name"] for a in proj.get("alphabets", [])):
        raise HTTPException(400, f"an alphabet named '{body['name']}' already exists in this project")
    alpha = {
        "id": uuid.uuid4().hex[:8],
        "name": body["name"],
        "letters": letters,
        "variables": body.get("variables") or [],
        "expressions": body.get("expressions") or [],
        "expr_loader": body.get("expr_loader"),
        "roots": body.get("roots") or {},
        "properties": [],
    }
    proj["alphabets"].append(alpha)
    storage.save_project(proj)
    return alpha


@app.put("/api/projects/{pid}/alphabets/{aid}")
def api_update_alphabet(pid: str, aid: str, body: dict = Body(...)) -> dict:
    proj = _get_project(pid)
    alpha = _get_alphabet(proj, aid)
    for key in ("name", "letters", "variables", "expressions", "expr_loader", "roots"):
        if key in body:
            alpha[key] = body[key]
    storage.save_project(proj)
    return alpha


@app.delete("/api/projects/{pid}/alphabets/{aid}")
def api_delete_alphabet(pid: str, aid: str) -> dict:
    proj = _get_project(pid)
    before = len(proj["alphabets"])
    proj["alphabets"] = [a for a in proj["alphabets"] if a["id"] != aid]
    if len(proj["alphabets"]) == before:
        raise HTTPException(404, "alphabet not found")
    storage.save_project(proj)
    return {"ok": True}


VALID_PROP_TYPES = {"integrability", "first_entry", "last_entry", "extended_steinmann", "cluster_adjacency", "transformation", "precomputed_tensor", "letter_symmetry"}


@app.post("/api/projects/{pid}/alphabets/{aid}/properties")
def api_add_property(pid: str, aid: str, body: dict = Body(...)) -> dict:
    proj = _get_project(pid)
    alpha = _get_alphabet(proj, aid)
    ptype = body.get("type")
    if ptype not in VALID_PROP_TYPES:
        raise HTTPException(400, f"property type must be one of {sorted(VALID_PROP_TYPES)}")
    params = body.get("params") or {}
    letters = set(alpha["letters"])
    if ptype in ("first_entry", "last_entry"):
        unknown = [x for x in params.get("letters", []) if x not in letters]
        if unknown:
            raise HTTPException(400, f"unknown letters: {unknown}")
    if ptype in ("extended_steinmann", "cluster_adjacency"):
        key = "nonadjacent_pairs" if ptype == "extended_steinmann" else "adjacent_pairs"
        for pair in params.get(key, []):
            if len(pair) != 2 or pair[0] not in letters or pair[1] not in letters:
                raise HTTPException(400, f"invalid letter pair: {pair}")
    if ptype == "transformation" and not params.get("name"):
        raise HTTPException(400, "transformation property needs a name")
    pname = (body.get("name") or "").strip()
    if ptype == "transformation" and not pname:
        pname = params["name"]
    if pname:
        for other in alpha.get("properties", []):
            existing = (other.get("name") or "").strip() or (other.get("params", {}).get("name") if other["type"] == "transformation" else "") or other["type"]
            if existing == pname:
                raise HTTPException(400, f"a property named '{pname}' already exists on this alphabet")
    if ptype == "precomputed_tensor":
        if not pname:
            raise HTTPException(400, "a precomputed tensor property needs a name")
        rel = (params.get("tensor_file") or "").strip().lstrip("/")
        if not rel or ".." in rel.split("/"):
            raise HTTPException(400, "provide a valid relative tensor file path (e.g. data/cycrepmat.wxf)")
        target = (storage.project_dir(pid) / rel).resolve()
        if not str(target).startswith(str(storage.project_dir(pid).resolve())):
            raise HTTPException(400, "tensor file must live inside the project")
        if not target.exists():
            raise HTTPException(400, f"tensor file '{rel}' does not exist in this project")
        params["tensor_file"] = rel
    if ptype == "letter_symmetry":
        if not pname:
            raise HTTPException(400, "a letter-rule symmetry property needs a name")
        if not (params.get("rule") or "").strip():
            raise HTTPException(400, "provide a replacement rule (e.g. {W[1]->W[2], W[2]->W[1]})")
        defs = (params.get("defs_file") or "").strip().lstrip("/")
        if defs:
            if ".." in defs.split("/"):
                raise HTTPException(400, "invalid definitions file path")
            target = (storage.project_dir(pid) / defs).resolve()
            if not str(target).startswith(str(storage.project_dir(pid).resolve())) or not target.exists():
                raise HTTPException(400, f"definitions file '{defs}' does not exist in this project")
            params["defs_file"] = defs
        else:
            params.pop("defs_file", None)
    prop = {
        "id": uuid.uuid4().hex[:8],
        "type": ptype,
        "name": pname,
        "params": params,
        "status": "ready" if ptype == "precomputed_tensor" else "pending",
        "error": None,
        "precomputed": ptype == "precomputed_tensor",
        "tensor_file": params.get("tensor_file") if ptype == "precomputed_tensor" else None,
        "summary": None,
    }
    alpha["properties"].append(prop)
    storage.save_project(proj)
    return prop


@app.delete("/api/projects/{pid}/alphabets/{aid}/properties/{prop_id}")
def api_delete_property(pid: str, aid: str, prop_id: str) -> dict:
    proj = _get_project(pid)
    alpha = _get_alphabet(proj, aid)
    before = len(alpha["properties"])
    alpha["properties"] = [p for p in alpha["properties"] if p["id"] != prop_id]
    if len(alpha["properties"]) == before:
        raise HTTPException(404, "property not found")
    storage.save_project(proj)
    return {"ok": True}


def _property_step(proj: dict, alpha: dict, prop: dict) -> dict:
    proj_dir = storage.project_dir(proj["id"])
    rel = property_tensor_relpath(alpha, prop)
    out_abs = str(proj_dir / rel)
    script = property_script(alpha, prop, out_abs)
    gen_dir = proj_dir / "wolfram_gen"
    gen_dir.mkdir(exist_ok=True)
    script_path = gen_dir / f"prop_{prop['id']}.wl"
    script_path.write_text(script)
    ws = find_wolframscript()
    if ws is None:
        raise HTTPException(500, "wolframscript not found on this machine")
    return {
        "id": "step-1",
        "label": f"Compute {property_display_name(prop)} ({prop['type'].replace('_', ' ')}) for alphabet '{alpha['name']}'",
        "kind": "wolfram",
        "command": f"{ws} -script {script_path}",
        "argv": [ws, "-script", str(script_path)],
        "cwd": str(proj_dir),
        "outputs": [rel],
        "skip_if_exists": False,
        "meta": {"type": "property", "alphabet_id": alpha["id"], "property_id": prop["id"], "tensor_file": rel},
    }


@app.post("/api/projects/{pid}/alphabets/{aid}/properties/{prop_id}/compute")
def api_compute_property(pid: str, aid: str, prop_id: str) -> dict:
    proj = _get_project(pid)
    alpha = _get_alphabet(proj, aid)
    prop = storage.find_property(alpha, prop_id)
    if prop is None:
        raise HTTPException(404, "property not found")
    if prop["type"] == "precomputed_tensor":
        raise HTTPException(400, "a precomputed tensor is a shipped file and cannot be recomputed")
    step = _property_step(proj, alpha, prop)
    steps = [step]
    if prop["type"] in ("first_entry", "last_entry"):
        tensor_ops = find_tensor_ops()
        if tensor_ops is not None:
            proj_dir = storage.project_dir(proj["id"])
            rel = step["meta"]["tensor_file"]
            proj_rel = f"data/{Path(rel).stem}_proj.wxf"
            steps.append({
                "id": "step-2",
                "label": f"Derive projection map {Path(rel).stem} -> {Path(proj_rel).stem} (drop size-1 axis)",
                "kind": "tensor_ops",
                "command": f"{tensor_ops} squeeze {proj_dir / rel} {proj_dir / proj_rel}",
                "argv": [tensor_ops, "squeeze", str(proj_dir / rel), str(proj_dir / proj_rel)],
                "cwd": str(proj_dir),
                "outputs": [proj_rel],
                "skip_if_exists": False,
                "meta": {"type": "squeeze", "tensor_file": proj_rel},
            })
    prop["status"] = "computing"
    prop["error"] = None
    storage.save_project(proj)

    def on_step_done(run, st, result):
        if st.get("meta", {}).get("type") != "property":
            return
        p2 = storage.load_project(pid)
        if p2 is None:
            return
        a2 = storage.find_alphabet(p2, aid)
        pr2 = storage.find_property(a2, prop_id) if a2 else None
        if pr2 is None:
            return
        if result and result.get("ok"):
            data = read_result_file(str(storage.project_dir(pid) / st["meta"]["tensor_file"])) or {}
            pr2["status"] = "ready"
            pr2["tensor_file"] = st["meta"]["tensor_file"]
            pr2["summary"] = {"dims": data.get("dims"), "nnz": data.get("nnz")}
        else:
            pr2["status"] = "error"
            pr2["error"] = "computation did not report success"
        storage.save_project(p2)

    def on_done(run):
        if run.status == "failed":
            p2 = storage.load_project(pid)
            if p2 is None:
                return
            a2 = storage.find_alphabet(p2, aid)
            pr2 = storage.find_property(a2, prop_id) if a2 else None
            if pr2 is not None and pr2.get("status") == "computing":
                pr2["status"] = "error"
                pr2["error"] = "wolframscript exited with an error; see run log"
                storage.save_project(p2)

    run = engine.create_run(pid, None, step["label"], [step], on_step_done=on_step_done, on_done=on_done)
    return {"run_id": run.run_id}


@app.get("/api/projects/{pid}/tensors")
def api_list_tensors(pid: str, dir: str = "data") -> list:
    _get_project(pid)
    if dir not in ("data", "output"):
        raise HTTPException(400, "dir must be data or output")
    base = storage.project_dir(pid) / dir
    out = []
    if base.exists():
        for f in sorted(base.rglob("*.wxf")):
            rel = f.relative_to(storage.project_dir(pid))
            if any(part.startswith(".") for part in rel.parts):
                continue
            out.append({"file": rel.as_posix(), "size_bytes": f.stat().st_size})
    return out


def summarize_tensor_file(pid: str, rel_file: str) -> dict:
    proj_dir = storage.project_dir(pid)
    target = (proj_dir / rel_file).resolve()
    if not str(target).startswith(str(proj_dir.resolve())) or not target.exists():
        raise HTTPException(404, "tensor file not found")
    cache_dir = proj_dir / ".summary_cache"
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / (rel_file.replace("/", "_") + ".json")
    mtime = target.stat().st_mtime
    if cache_file.exists():
        try:
            cached = json.loads(cache_file.read_text())
            if cached.get("mtime") == mtime:
                return cached["summary"]
        except Exception:
            pass
    ws = find_wolframscript()
    if ws is None:
        raise HTTPException(500, "wolframscript not found")
    gen_dir = proj_dir / "wolfram_gen"
    gen_dir.mkdir(exist_ok=True)
    script_path = gen_dir / "summary.wl"
    result_path = cache_dir / (rel_file.replace("/", "_") + ".result.json")
    if result_path.exists():
        result_path.unlink()
    script_path.write_text(summary_script(str(target), str(result_path)))
    proc = subprocess.run([ws, "-script", str(script_path)], capture_output=True, text=True, timeout=600)
    result = None
    if result_path.exists():
        try:
            result = json.loads(result_path.read_text())
        except Exception:
            result = None
    if result is None or not result.get("ok"):
        raise HTTPException(500, f"failed to summarize tensor: {proc.stdout[-500:]} {proc.stderr[-500:]}")
    summary = {"dims": result.get("dims"), "nnz": result.get("nnz"), "sample": result.get("sample", [])}
    cache_file.write_text(json.dumps({"mtime": mtime, "summary": summary}))
    return summary


@app.get("/api/projects/{pid}/tensor_summary")
def api_tensor_summary(pid: str, file: str) -> dict:
    _get_project(pid)
    return summarize_tensor_file(pid, file)


@app.post("/api/projects/{pid}/alphabets/{aid}/properties/{prop_id}/summarize")
def api_summarize_property(pid: str, aid: str, prop_id: str) -> dict:
    proj = _get_project(pid)
    alpha = _get_alphabet(proj, aid)
    prop = storage.find_property(alpha, prop_id)
    if prop is None:
        raise HTTPException(404, "property not found")
    rel = prop.get("tensor_file") or property_tensor_relpath(alpha, prop)
    summary = summarize_tensor_file(pid, rel)
    prop["tensor_file"] = rel
    prop["summary"] = {"dims": summary.get("dims"), "nnz": summary.get("nnz")}
    storage.save_project(proj)
    return prop


@app.post("/api/projects/{pid}/flows")
def api_create_flow(pid: str, body: dict = Body(...)) -> dict:
    proj = _get_project(pid)
    flow = {
        "id": uuid.uuid4().hex[:8],
        "name": body.get("name") or "Untitled flow",
        "graph": {"nodes": [], "edges": []},
    }
    if body.get("custom_block"):
        flow["custom_block"] = True
    proj["flows"].append(flow)
    storage.save_project(proj)
    return flow


@app.post("/api/projects/{pid}/flows/{fid}/duplicate")
def api_duplicate_flow(pid: str, fid: str, body: dict = Body(default={})) -> dict:
    import copy as _copy
    proj = _get_project(pid)
    flow = storage.find_flow(proj, fid)
    if flow is None:
        raise HTTPException(404, "flow not found")
    dup = {
        "id": uuid.uuid4().hex[:8],
        "name": (body.get("name") or f"Copy of {flow['name']}")[:120],
        "graph": _copy.deepcopy(flow.get("graph") or {"nodes": [], "edges": []}),
        "output_subdir": f"flow_{uuid.uuid4().hex[:6]}",
    }
    if flow.get("custom_block"):
        dup["custom_block"] = True
    proj["flows"].append(dup)
    storage.save_project(proj)
    return dup


@app.put("/api/projects/{pid}/flows/{fid}")
def api_update_flow(pid: str, fid: str, body: dict = Body(...)) -> dict:
    proj = _get_project(pid)
    flow = storage.find_flow(proj, fid)
    if flow is None:
        raise HTTPException(404, "flow not found")
    if "name" in body:
        flow["name"] = body["name"]
    if "graph" in body:
        flow["graph"] = body["graph"]
    if "custom_block" in body:
        if body["custom_block"]:
            flow["custom_block"] = True
        else:
            flow.pop("custom_block", None)
    storage.save_project(proj)
    return flow


@app.delete("/api/projects/{pid}/flows/{fid}")
def api_delete_flow(pid: str, fid: str) -> dict:
    proj = _get_project(pid)
    before = len(proj["flows"])
    proj["flows"] = [f for f in proj["flows"] if f["id"] != fid]
    if len(proj["flows"]) == before:
        raise HTTPException(404, "flow not found")
    storage.save_project(proj)
    return {"ok": True}


@app.get("/api/projects/{pid}/flow_outputs")
def api_flow_outputs_catalog(pid: str) -> list:
    proj = _get_project(pid)
    return flow_outputs_catalog(proj)


@app.post("/api/projects/{pid}/flows/{fid}/compile")
def api_compile_flow(pid: str, fid: str) -> dict:
    proj = _get_project(pid)
    flow = storage.find_flow(proj, fid)
    if flow is None:
        raise HTTPException(404, "flow not found")
    result = compile_flow(proj, flow.get("graph") or {}, output_subdir=flow.get("output_subdir"))
    result.pop("_steps_full", None)
    return result


@app.post("/api/projects/{pid}/flows/{fid}/runs")
def api_run_flow(pid: str, fid: str) -> dict:
    proj = _get_project(pid)
    flow = storage.find_flow(proj, fid)
    if flow is None:
        raise HTTPException(404, "flow not found")
    result = compile_flow(proj, flow.get("graph") or {}, output_subdir=flow.get("output_subdir"))
    if not result["ok"]:
        raise HTTPException(400, {"message": "flow does not compile", "errors": result["errors"]})
    steps = result["_steps_full"]
    if not steps:
        raise HTTPException(400, "flow produced no steps")
    run = engine.create_run(pid, fid, f"Flow: {flow['name']}", steps, on_step_done=_flow_step_done(pid))
    return {"run_id": run.run_id}


def _flow_step_done(pid: str):
    def hook(run, st, result):
        meta = st.get("meta") or {}
        if meta.get("type") != "property" or not (result and result.get("ok")):
            return
        p2 = storage.load_project(pid)
        if p2 is None:
            return
        a2 = storage.find_alphabet(p2, meta.get("alphabet_id", ""))
        pr2 = storage.find_property(a2, meta.get("property_id", "")) if a2 else None
        if pr2 is not None:
            data = read_result_file(str(storage.project_dir(pid) / meta["tensor_file"])) or {}
            pr2["status"] = "ready"
            pr2["tensor_file"] = meta.get("tensor_file")
            pr2["summary"] = {"dims": data.get("dims"), "nnz": data.get("nnz")}
            storage.save_project(p2)
    return hook


@app.get("/api/projects/{pid}/runs")
def api_list_runs(pid: str) -> list:
    _get_project(pid)
    return [r.snapshot() for r in engine.list_for_project(pid)]


@app.get("/api/runs/{run_id}")
def api_get_run(run_id: str) -> dict:
    run = engine.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")
    return run.snapshot()


@app.post("/api/runs/{run_id}/cancel")
def api_cancel_run(run_id: str) -> dict:
    if not engine.cancel(run_id):
        raise HTTPException(400, "run not found or not running")
    return {"ok": True}


@app.get("/api/runs/{run_id}/events")
def api_run_events(run_id: str):
    run = engine.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")

    def gen():
        idx = 0
        while True:
            with run.cond:
                if idx >= len(run.events) and run.status in ("queued", "running"):
                    run.cond.wait(timeout=15)
                batch = run.events[idx:]
                idx += len(batch)
                done = run.status not in ("queued", "running") and idx >= len(run.events)
            if not batch and not done:
                yield ": keepalive\n\n"
                continue
            for ev in batch:
                yield f"event: {ev['event']}\ndata: {json.dumps(ev['data'])}\n\n"
            if done:
                return

    return StreamingResponse(gen(), media_type="text/event-stream")


if WEB_DIST.exists():
    @app.get("/{full_path:path}")
    def spa(full_path: str):
        candidate = (WEB_DIST / full_path).resolve()
        if full_path and candidate.exists() and candidate.is_file() and str(candidate).startswith(str(WEB_DIST.resolve())):
            headers = {"Cache-Control": "public, max-age=31536000, immutable"} if full_path.startswith("assets/") else None
            return FileResponse(candidate, headers=headers)
        return FileResponse(WEB_DIST / "index.html", headers={"Cache-Control": "no-cache"})
