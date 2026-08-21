# Front-End API Contract (backend <-> frontend)

Single source of truth for the Milestone-1 vertical slice.

- Backend: Python 3.9-compatible FastAPI + uvicorn, in `front-end/server/`.
- Frontend: React + Vite + React Flow, in `front-end/web/`, dev proxy to backend.
- Backend listens on `http://127.0.0.1:8321`. All API routes prefixed `/api`.
- Backend serves the built SPA from `front-end/web/dist/` at `/` (with SPA fallback to index.html) when the dist folder exists.
- All request/response bodies are JSON. Errors: `{"detail": "message"}` with proper HTTP status.

## Data model

Project (stored at `front-end/projects/<project_id>/project.json`; each project has `data/` and `output/` subdirs, and `runs/` for logs):

```json
{
  "id": "string (slug)",
  "name": "string",
  "created_at": "iso8601",
  "alphabets": [Alphabet],
  "flows": [Flow]
}
```

Alphabet (the "building block"):

```json
{
  "id": "string",
  "name": "string",
  "letters": ["x", "y", "z", ...],
  "variables": ["x", "y", ...],
  "expressions": ["x", "y", "1/x", ...],
  "roots": {"eps5": "Sqrt[Delta5]"},
  "properties": [Property]
}
```

Property types and params:

| type | params | tensor produced |
|---|---|---|
| `integrability` | `{}` | dlogmat (n,n,s) |
| `first_entry` | `{"letters": ["x","z",...]}` | FEC_1 seed (k,1,n) |
| `last_entry` | `{"letters": [...]}` | LEC_1 seed (k,n,1) |
| `extended_steinmann` | `{"nonadjacent_pairs": [["x","y"],...]}` (letter names) | ES dlogmat |
| `cluster_adjacency` | `{"adjacent_pairs": [["x","y"],...]}` | CA dlogmat |
| `transformation` | `{"name": "cyclic", "map": "{x->y, y->z, z->x}"}` (Wolfram rules string) | n x n rational matrix |
| `precomputed_tensor` | `{"tensor_file": "data/cycrepmat.wxf"}` | registers an existing project file (ready immediately, no compute path; matrix kind) |
| `letter_symmetry` | `{"rule": "{W[1]->W[2], ...}", "defs_file": "data/defs.wl"?}` | computes `CoefficientArrays[letters /. rule, letters][[2]]`, the n x n letter-space symmetry matrix (matrix kind) |

Property object:

```json
{
  "id": "string",
  "type": "integrability",
  "params": {},
  "status": "pending|computing|ready|error",
  "error": null,
  "tensor_file": "data/dlogmat_<name>.wxf (relative to project dir, null until ready)",
  "summary": {"dims": [7,7,21], "nnz": 126}
}
```

Flow:

```json
{
  "id": "string",
  "name": "string",
  "graph": {"nodes": [...], "edges": [...]}
}
```

Graph nodes (React Flow compatible): `{"id", "type", "position": {"x","y"}, "data": {...}}`.
Node `type` values and `data`:

| type | data | typed outputs |
|---|---|---|
| `alphabet` | `{"alphabet_id": "...", "selected_properties": ["prop_id", ...]}` | one output handle per selected property: `prop_<prop_id>` with tensor kind (`dlogmat`, `fec1`, `lec1`, `matrix`). Selected `first_entry`/`last_entry` properties additionally expose `proj_<prop_id>` (kind `matrix`): the projection map derived by dropping the size-1 axis of the rank-3 tensor (e.g. (a,1,b) → (a,b)), auto-computed as `data/<stem>_proj.wxf` via `tensor_ops squeeze` |
| `merge_conditions` | `{}` | input handles `in_0..in_n` (kind `dlogmat`), output `out` (kind `dlogmat`) |
| `extend` | `{"target_weight": 2}` | inputs `condition` (dlogmat) plus exactly one of `fec` (fec1/fec) or `lec` (lec1/lec); output `fec`/`lec` matching the input direction (weight +1) |
| `sew` | `{}` | inputs `condition`, `fec`, `lec`; output `sew` |
| `project` | `{"symmetry": "collinear", "target": "SEW_5p1"}` | input `seed`; output `basis` |
| `solve_symmetry` | `{"projection_file": "cycrepmat.wxf"}` | inputs; output `solution` |
| `solve_collinear` | `{"target": "SEW_3p1", "rhs": "boundary_2L.wxf", "projection": "finite"}` | output `solution` |
| `compute_rhs` | `{"loops": 2}` | output `boundary` |
| `add_tensors` | `{"weight_a": "1", "weight_b": "-1", "target": "SEW_3p1_total"}` | inputs `a`, `b` (any tensor kind, must match in kind and weight); output `out` carrying the inputs' kind/weight; runs `tensor_add A B wA wB out` |
| `ternary_contract` | `{"target": "FEC_2_sym"}` | inputs `tensor` (rank-3: fec1/fec/lec1/lec/sew), `trans1`, `trans2` (any matrices — symmetry transformations included); TernaryContract T'[a,b',c'] = Σ T[a,b,c]·M1[b,b']·M2[c,c'] on the last two entries; runs `tensor_ops ternary T M1 M2 out`. Legacy graphs saved with node type `apply_symmetry` are compiled identically |
| `matrix_power` | `{"n": 2, "target": "cyc2"}` | input `matrix`; output `matrix`; runs `tensor_ops power M n out` |
| `tensor_join` | `{"axis": -1, "target": "FEC_1_dup"}` | inputs `a`, `b` (same kind); output keeps the kind (weight unset); runs `tensor_ops join A B axis out`; axis 1-based, negative counts from the end |
| `impose_integrability` | `{"target": "NMHV_E14_w2f_integ", "transpose": true}` | inputs `tensor` (rank-3: fec1/fec/lec1/lec/sew), `dlogmat` (integrability); output `matrix` solution basis; runs `tensor_ops impose T D trans out`; contracts S[s,i,a]·D[a,i,c] and SparseRREF-solves — with `transpose` (default) the kernel of Mᵀ (tensor coefficients satisfying all conditions), without it the kernel of M (relations among conditions) |
| `assemble` | `{"groups": "1,1,2,3", "coefs": "1/2,2,0", "target": "sol_A"}` | inputs `in_0`, `in_1`, … (element tensors, same kind; any count — the canvas adds an empty slot as you connect; handle order = alignment order); output keeps the kind in the common aligned frame; runs `tensor_ops assemble out --elems … --groups … --coefs …`; assembles A = Σ cⱼ·Bⱼ as P_A·Join[c-scaled elements, 1] — groups with coefficient 0 keep their rows in the frame but zeroed |

Graphs may also carry a `groups` array (`{"id","name","node_ids","collapsed","position"}`) for collapsible canvas groups; groups are view-only and are expanded to their member nodes at compile time.

Edges: `{"id", "source", "sourceHandle", "target", "targetHandle"}`.
Tensor-kind compatibility (validated at compile time):
- `dlogmat` -> merge_conditions.in_*, extend.condition, sew.condition
- `fec1`/`fec` -> extend.fec, sew.fec
- `lec1` -> sew.lec
- `matrix` -> (stored as data file; used by solve_symmetry via projection_file)

Compiled step:

```json
{
  "id": "step-1",
  "label": "Compute integrability (dlogmat) for alphabet X",
  "kind": "wolfram|bootstrap",
  "command": "full display command string",
  "outputs": ["data/dlogmat_X.wxf"],
  "skip_if_exists": true
}
```

## REST endpoints

- `GET /api/env` -> `{"repo_root": "...", "wolframscript": {"found": true, "path": "..."}, "bootstrap": {"found": true, "path": "..."}, "projects_dir": "..."}`
- `GET /api/templates` -> `[{"id": "4pformfactor", "name": "...", "description": "...", "n_letters": 7}]` (also `pentagon`, `e6`)
- `GET /api/projects` -> `[{"id","name","created_at","n_alphabets","n_flows"}]`
- `POST /api/projects` body `{"name": "...", "template_id": "4pformfactor"|null}` -> Project (template copies the example's data files into project data/ and pre-creates its alphabet with ready properties)
- `GET /api/projects/{pid}` -> Project (full)
- `DELETE /api/projects/{pid}`
- `POST /api/projects/{pid}/alphabets` body Alphabet (without id/properties) -> Alphabet
- `PUT /api/projects/{pid}/alphabets/{aid}` body Alphabet -> Alphabet
- `DELETE /api/projects/{pid}/alphabets/{aid}`
- `POST /api/projects/{pid}/alphabets/{aid}/properties` body `{"type": "...", "params": {...}}` -> Property (status pending)
- `DELETE /api/projects/{pid}/alphabets/{aid}/properties/{prop_id}`
- `POST /api/projects/{pid}/alphabets/{aid}/properties/{prop_id}/compute` -> `{"run_id": "..."}` (async run; poll/SSE via run endpoints)
- `GET /api/projects/{pid}/tensors?dir=data|output` -> `[{"file": "data/x.wxf", "size_bytes": 123}]`
- `GET /api/projects/{pid}/tensor_summary?file=data/x.wxf` -> `{"dims": [...], "nnz": N, "sample": [{"index": [i,j,k], "value": "1/2"}, ...]}`
- `POST /api/projects/{pid}/flows` body `{"name": "..."}` -> Flow
- `PUT /api/projects/{pid}/flows/{fid}` body `{"name": "...", "graph": {...}}` -> Flow
- `DELETE /api/projects/{pid}/flows/{fid}`
- `POST /api/projects/{pid}/flows/{fid}/compile` -> `{"ok": true, "errors": [], "steps": [Step]}`
- `POST /api/projects/{pid}/flows/{fid}/runs` -> `{"run_id": "..."}` (compiles then executes)
- `GET /api/projects/{pid}/runs` -> `[RunSummary]` newest first
- `GET /api/runs/{run_id}` -> Run detail
- `POST /api/runs/{run_id}/cancel`
- `GET /api/runs/{run_id}/events` -> SSE stream

Run summary/detail:

```json
{
  "run_id": "string",
  "project_id": "string",
  "flow_id": "string|null",
  "label": "string",
  "status": "queued|running|done|failed|cancelled",
  "created_at": "iso8601",
  "steps": [{"id","label","kind","command","status","returncode"}],
  "current_step": "step id|null"
}
```

SSE events (`GET /api/runs/{run_id}/events`), `text/event-stream`; replay backlog on connect, then live; ends with `event: end`:
- `event: log`  data: `{"step_id": "...", "stream": "stdout|stderr", "line": "..."}`
- `event: step` data: `{"step_id": "...", "status": "running|done|failed|skipped"}`
- `event: status` data: `{"status": "running|done|failed|cancelled"}`
- `event: end`  data: `{}`

## Execution semantics

- Steps execute sequentially in topological order; a step with all `outputs` present and `skip_if_exists` is marked `skipped`.
- Subprocesses run in their own process group; cancel kills the whole group. Browser disconnect never affects a run.
- Logs persist to `projects/<pid>/runs/<run_id>.log`; run metadata to `<run_id>.json`.
- `wolfram` steps: backend generates a `.wl` script under `projects/<pid>/wolfram_gen/`, runs `wolframscript -script <file>`; script prints a marker line `@@RESULT@@<json>` with tensor dims/nnz.
- `bootstrap` steps run `<repo_root>/bootstrap <mode> ... --data-dir <project data> --output-dir <project output>`.
