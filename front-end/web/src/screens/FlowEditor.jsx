import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { useNavigate, useParams } from 'react-router-dom'
import ReactFlow, {
  Background, ConnectionMode, Controls, Handle, Position, ReactFlowProvider,
  addEdge, useEdges, useEdgesState, useNodesState, useUpdateNodeInternals,
} from 'reactflow'
import 'reactflow/dist/style.css'
import { api } from '../api'
import { useProject, useToast } from '../App'
import {
  NODE_DEFS, PALETTE_SECTIONS, PROP_KIND, propLabel,
  kindsCompatible, sourceKindFor, targetKindFor,
} from '../flowdefs'

const ROW0 = 29
const ROW_STEP = 15
const SUB_EXTRA = 15

function StatusDot({ status }) {
  return <span className={`dot ${status || 'pending'}`} />
}

function AlphabetNode({ id, data, selected }) {
  const { project } = useProject()
  const alphabet = project?.alphabets.find((a) => a.id === data.alphabet_id)
  const selectedProps = data.selected_properties || []
  const outRows = []
  if (alphabet) {
    for (const p of alphabet.properties.filter((x) => selectedProps.includes(x.id))) {
      outRows.push({ key: p.id, handle: `prop_${p.id}`, label: propLabel(p), kind: PROP_KIND[p.type], status: p.status })
      if (p.type === 'first_entry' || p.type === 'last_entry') {
        outRows.push({ key: `${p.id}_proj`, handle: `proj_${p.id}`, label: `${propLabel(p)} · proj map`, kind: 'matrix', status: p.status })
      }
    }
  }
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: NODE_DEFS.alphabet.color }}>
        Alphabet: {alphabet ? alphabet.name : '—'}
      </div>
      <div className="node-body">
        {!alphabet && <span>select in inspector →</span>}
        {alphabet && selectedProps.length === 0 && <span>no outputs selected</span>}
        {outRows.map((r, i) => (
          <div key={r.key} className="handle-row" style={{ textAlign: 'right' }}>
            <span className={`kind-${r.kind}`}>
              {r.label}
            </span>{' '}
            <StatusDot status={r.status} />
            <Handle
              type="source" position={Position.Right} id={r.handle}
              style={{ top: ROW0 + i * ROW_STEP, background: 'var(--accent)' }}
            />
          </div>
        ))}
      </div>
    </div>
  )
}

function OpNode({ type, data, selected }) {
  const def = NODE_DEFS[type]
  const inputs = def.inputs || []
  const outputs = def.outputs || []
  const rows = Math.max(inputs.length, outputs.length)
  const subtitle =
    type === 'extend' ? (data.target_weight ? `→ weight ${data.target_weight}` : 'weight +1')
    : type === 'project' ? `${data.symmetry || '…'} ${data.target || ''}`
    : type === 'solve_symmetry' ? `${data.symmetry || '…'} ${data.target || ''}`
    : type === 'solve_collinear' ? `${data.target || '…'}`
    : type === 'compute_rhs' ? `${data.target || '…'}`
    : type === 'add_tensors' ? `${data.weight_a || '1'}·A + ${data.weight_b || '1'}·B → ${data.target || '…'}`
    : (type === 'ternary_contract' || type === 'apply_symmetry') ? `→ ${data.target || '…'}`
    : type === 'matrix_power' ? `M^${data.n || '?'} → ${data.target || '…'}`
    : type === 'tensor_join' ? `axis ${data.axis || '?'} → ${data.target || '…'}`
    : type === 'impose_integrability' ? `${data.transpose === false ? 'relations among conditions' : 'solve coefficients'} → ${data.target || '…'}`
    : ''
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        {subtitle && <div style={{ marginBottom: 2 }}>{subtitle}</div>}
        {Array.from({ length: rows }).map((_, i) => (
          <div key={i} className="handle-row" style={{ display: 'flex', justifyContent: 'space-between' }}>
            <span className={inputs[i] ? `kind-${inputs[i].kind === 'seed' ? 'fec' : inputs[i].kind}` : ''}>
              {inputs[i]?.label || ''}
            </span>
            <span className={outputs[i] ? `kind-${outputs[i].kind}` : ''}>{outputs[i]?.label || ''}</span>
          </div>
        ))}
        {inputs.map((inp, i) => (
          <Handle
            key={inp.id} type="target" position={Position.Left} id={inp.id}
            style={{ top: (subtitle ? ROW0 + SUB_EXTRA : ROW0) + i * ROW_STEP }}
          />
        ))}
        {outputs.map((out, i) => (
          <Handle
            key={out.id} type="source" position={Position.Right} id={out.id}
            style={{ top: (subtitle ? ROW0 + SUB_EXTRA : ROW0) + i * ROW_STEP, background: 'var(--accent)' }}
          />
        ))}
      </div>
    </div>
  )
}

function AssembleNode({ id, data, selected }) {
  const edges = useEdges()
  const updateNodeInternals = useUpdateNodeInternals()
  const def = NODE_DEFS.assemble
  let maxConnected = -1
  for (const e of edges) {
    if (e.target === id && (e.targetHandle || '').startsWith('in_')) {
      const idx = parseInt(e.targetHandle.slice(3).split('@')[0], 10)
      if (!Number.isNaN(idx)) maxConnected = Math.max(maxConnected, idx)
    }
  }
  const slots = Math.max(def.inputs.length, maxConnected + 2)
  useEffect(() => { updateNodeInternals(id) }, [id, slots, updateNodeInternals])
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        <div style={{ marginBottom: 2 }}>{`A = Σc·B → ${data.target || '…'}`}</div>
        {Array.from({ length: slots }).map((_, i) => (
          <div key={i} className="handle-row" style={{ display: 'flex', justifyContent: 'space-between' }}>
            <span className="kind-tensor">elem {i + 1}</span>
            <span className={i === 0 ? 'kind-tensor' : ''}>{i === 0 ? def.outputs[0].label : ''}</span>
          </div>
        ))}
        {Array.from({ length: slots }).map((_, i) => (
          <Handle
            key={`in_${i}`} type="target" position={Position.Left} id={`in_${i}`}
            style={{ top: ROW0 + SUB_EXTRA + i * ROW_STEP }}
          />
        ))}
        <Handle
          type="source" position={Position.Right} id="out"
          style={{ top: ROW0 + SUB_EXTRA, background: 'var(--accent)' }}
        />
      </div>
    </div>
  )
}

const nodeTypes = {
  alphabet: AlphabetNode,
  merge_conditions: (p) => <OpNode {...p} type="merge_conditions" />,
  extend: (p) => <OpNode {...p} type="extend" />,
  sew: (p) => <OpNode {...p} type="sew" />,
  project: (p) => <OpNode {...p} type="project" />,
  solve_symmetry: (p) => <OpNode {...p} type="solve_symmetry" />,
  solve_collinear: (p) => <OpNode {...p} type="solve_collinear" />,
  compute_rhs: (p) => <OpNode {...p} type="compute_rhs" />,
  add_tensors: (p) => <OpNode {...p} type="add_tensors" />,
  ternary_contract: (p) => <OpNode {...p} type="ternary_contract" />,
  apply_symmetry: (p) => <OpNode {...p} type="apply_symmetry" />,
  matrix_power: (p) => <OpNode {...p} type="matrix_power" />,
  tensor_join: (p) => <OpNode {...p} type="tensor_join" />,
  impose_integrability: (p) => <OpNode {...p} type="impose_integrability" />,
  assemble: AssembleNode,
  groupBox: GroupNode,
}

function Inspector({ node, onChange, onDelete, groupOps }) {
  const { project } = useProject()
  if (!node) return (
    <div>
      <p className="muted">Select a node to edit its parameters.</p>
      <p className="muted" style={{ fontSize: 11 }}>Tips: Shift+drag or Cmd/Ctrl+click selects multiple blocks · Cmd/Ctrl+C / V / D copy, paste, duplicate · Delete removes the selection.</p>
    </div>
  )
  if (node.type === 'groupBox') {
    const g = groupOps.groups.find((x) => x.id === node.id)
    return (
      <div>
        <label>Group name</label>
        <input value={g?.name || ''} onChange={(e) => groupOps.onRename(node.id, e.target.value)} />
        <p className="muted">Collapsed group of {(g?.node_ids || []).length} blocks. Wires crossing the group boundary are shown dashed.</p>
        <div className="row" style={{ marginTop: 8 }}>
          <button className="small" onClick={() => groupOps.onExpand(node.id)}>Expand</button>
          <button className="small danger" onClick={() => groupOps.onUngroup(node.id)}>Ungroup</button>
        </div>
      </div>
    )
  }
  const d = node.data || {}
  const set = (patch) => onChange(node.id, patch)
  const myGroup = groupOps.groups.find((g) => !g.collapsed && g.node_ids.includes(node.id))

  if (node.type === 'alphabet') {
    const alphabet = project?.alphabets.find((a) => a.id === d.alphabet_id)
    return (
      <div>
        <h3>Alphabet node</h3>
        <label>Alphabet</label>
        <select value={d.alphabet_id || ''} onChange={(e) => set({ alphabet_id: e.target.value, selected_properties: [] })}>
          <option value="">— choose —</option>
          {project?.alphabets.map((a) => <option key={a.id} value={a.id}>{a.name}</option>)}
        </select>
        {alphabet && (
          <>
            <label>Exposed outputs (properties)</label>
            {alphabet.properties.map((p) => (
              <div key={p.id} style={{ padding: '2px 0' }}>
                <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, margin: 0, color: 'var(--text)' }}>
                  <input
                    type="checkbox" style={{ width: 'auto' }}
                    checked={(d.selected_properties || []).includes(p.id)}
                    onChange={(e) => {
                      const cur = d.selected_properties || []
                      set({ selected_properties: e.target.checked ? [...cur, p.id] : cur.filter((x) => x !== p.id) })
                    }}
                  />
                  <StatusDot status={p.status} />
                  <span className={`kind-${PROP_KIND[p.type]}`}>
                    {propLabel(p)}
                  </span>
                </label>
              </div>
            ))}
            {!alphabet.properties.length && <p className="muted">This alphabet has no properties yet — add them on the Materials screen.</p>}
            <p className="muted" style={{ fontSize: 11 }}>Properties that are not ready will be computed automatically when the flow runs. First/Last Entry outputs also expose a “proj map” matrix output: the rank-3 tensor with its size-1 axis dropped (e.g. (a,1,b) → (a,b)), derived automatically as a named matrix.</p>
          </>
        )}
        <div style={{ marginTop: 16, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
          {myGroup && (
            <button className="small" style={{ marginBottom: 8 }} onClick={() => groupOps.onCollapse(myGroup.id)}>
              Collapse group “{myGroup.name}”
            </button>
          )}
          <button className="danger" onClick={() => onDelete(node.id)}>Delete this node</button>
        </div>
      </div>
    )
  }

  return (
    <div>
      <h3>{NODE_DEFS[node.type]?.title || node.type}</h3>
      {node.type === 'extend' && (
        <>
          <label>Target weight</label>
          <input type="number" min="2" value={d.target_weight || ''} onChange={(e) => set({ target_weight: e.target.value ? Number(e.target.value) : null })} />
          <p className="muted" style={{ fontSize: 11 }}>Must equal input FEC/LEC weight + 1 (checked at compile time). Connect exactly one of FEC in / LEC in — LEC extends backward.</p>
        </>
      )}
      {node.type === 'project' && (
        <>
          <label>Symmetry</label>
          <select value={d.symmetry || ''} onChange={(e) => set({ symmetry: e.target.value })}>
            <option value="">— choose —</option>
            {['collinear', 'cyclic', 'flip', 'parity'].map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
          <label>Target (e.g. SEW_3p1 or FEC_3; derived from seed input if empty)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} />
        </>
      )}
      {node.type === 'solve_symmetry' && (
        <>
          <label>Symmetry</label>
          <select value={d.symmetry || ''} onChange={(e) => set({ symmetry: e.target.value })}>
            <option value="">— choose —</option>
            {['cyclic', 'flip', 'parity'].map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} />
        </>
      )}
      {node.type === 'solve_collinear' && (
        <>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="SEW_3p1" />
          <label>RHS file (or 0)</label>
          <input value={d.rhs || ''} onChange={(e) => set({ rhs: e.target.value })} placeholder="output/2loop/boundary_2L.wxf" />
          <label>Projection</label>
          <select value={d.projection || 'finite'} onChange={(e) => set({ projection: e.target.value })}>
            <option value="finite">finite</option>
            <option value="divergent">divergent</option>
          </select>
          <label>Letter projection (file or identity)</label>
          <input value={d.letter_projection || ''} onChange={(e) => set({ letter_projection: e.target.value })} placeholder="identity" />
        </>
      )}
      {node.type === 'compute_rhs' && (
        <>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="SEW_3p1" />
          <label>Letter projection (file or identity)</label>
          <input value={d.letter_projection || ''} onChange={(e) => set({ letter_projection: e.target.value })} placeholder="identity" />
        </>
      )}
      {node.type === 'add_tensors' && (
        <>
          <label>Weight of A (rational)</label>
          <input value={d.weight_a || ''} onChange={(e) => set({ weight_a: e.target.value })} placeholder="1" />
          <label>Weight of B (rational)</label>
          <input value={d.weight_b || ''} onChange={(e) => set({ weight_b: e.target.value })} placeholder="1" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. SEW_3p1_total" />
          <p className="muted" style={{ fontSize: 11 }}>Computes wA·A + wB·B with exact rational arithmetic (tensor_add). A and B must have identical dimensions (same kind and weight).</p>
        </>
      )}
      {(node.type === 'ternary_contract' || node.type === 'apply_symmetry') && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. FEC_2_sym" />
          <p className="muted" style={{ fontSize: 11 }}>TernaryContract: contracts trans1 with the 2nd-to-last axis and trans2 with the last axis of the rank-3 input tensor (T&apos;[a,b&apos;,c&apos;] = Σ T·M1·M2), exact rational arithmetic via tensor_ops. The transformations can be any matrices — symmetry transformations included.</p>
        </>
      )}
      {node.type === 'matrix_power' && (
        <>
          <label>Power n (non-negative integer)</label>
          <input value={d.n || ''} onChange={(e) => set({ n: e.target.value })} placeholder="2" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. cyc2" />
          <p className="muted" style={{ fontSize: 11 }}>M^n with exact rational arithmetic (binary exponentiation). n=0 gives the identity. Use it to generate group elements M, M², … feeding Ternary Contract.</p>
        </>
      )}
      {node.type === 'tensor_join' && (
        <>
          <label>Axis (1-based; negative counts from the end)</label>
          <input value={d.axis || ''} onChange={(e) => set({ axis: e.target.value })} placeholder="1 = first entry, -1 = last entry" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. FEC_1_joined" />
          <p className="muted" style={{ fontSize: 11 }}>Like Mathematica Join: all dimensions except the join axis must match; the join axis dimension grows.</p>
        </>
      )}
      {node.type === 'impose_integrability' && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. NMHV_E14_w2f_integ" />
          <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, marginTop: 8, color: 'var(--text)' }}>
            <input
              type="checkbox" style={{ width: 'auto' }}
              checked={d.transpose !== false}
              onChange={(e) => set({ transpose: e.target.checked })}
            />
            Transpose before solving (solve for tensor coefficients)
          </label>
          <p className="muted" style={{ fontSize: 11 }}>
            Contracts the last two axes of the rank-3 tensor S[s,i,a] with the integrability dlog D[a,i,c]
            (M[s,c] = Σ S[s,i,a]·D[a,i,c], like TensorContract[S . dlogmat, {'{2, 3}'}]), then row-reduces
            with SparseRREF. With transpose on, the output basis is the kernel of Mᵀ — combinations of the
            s tensors satisfying all c conditions. With transpose off, it is the kernel of M — linear
            relations among the conditions.
          </p>
        </>
      )}
      {node.type === 'assemble' && (
        <>
          <label>Groups (one group index per connected element, in elem order top to bottom)</label>
          <input value={d.groups || ''} onChange={(e) => set({ groups: e.target.value })} placeholder="e.g. 1,1,2,3 → B1={elem1,elem2}, B2={elem3}, B3={elem4}" />
          <label>Coefficients (one rational per group; 0 excludes the group)</label>
          <input value={d.coefs || ''} onChange={(e) => set({ coefs: e.target.value })} placeholder="e.g. 1/2,2,0 → A = (1/2)·B1 + 2·B2" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. sol_A" />
          <p className="muted" style={{ fontSize: 11 }}>
            Assembles A = Σ cⱼ·Bⱼ in the common aligned frame: the element solution spaces are scaled by
            their group&apos;s coefficient and joined along the first axis in the given order, then projected
            by P_A — rows of groups with coefficient 0 stay in the frame but are zeroed. Connect any number
            of elements (a new empty slot appears as you wire); all must have the same kind and matching
            trailing dimensions. Runs tensor_ops assemble.
          </p>
        </>
      )}
      {node.type === 'merge_conditions' && <p className="muted">Connect two or more dlogmat outputs (integrability, extended Steinmann, cluster adjacency) to merge them into a single condition tensor.</p>}
      {node.type === 'sew' && <p className="muted">Combines a condition tensor, an FEC tensor of weight F and an LEC tensor of weight L into SEW_FpL.</p>}
      <div style={{ marginTop: 16, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
        {myGroup && (
          <button className="small" style={{ marginBottom: 8 }} onClick={() => groupOps.onCollapse(myGroup.id)}>
            Collapse group “{myGroup.name}”
          </button>
        )}
        <button className="danger" onClick={() => onDelete(node.id)}>Delete this node</button>
      </div>
    </div>
  )
}

function CompilePanel({ result, onRun, running }) {
  if (!result) return null
  return (
    <div>
      <h3>Compiled plan</h3>
      {!result.ok && (
        <div>
          <p className="error-text">The flow does not compile:</p>
          <ul className="error-text">
            {result.errors.map((e, i) => <li key={i}>{e}</li>)}
          </ul>
        </div>
      )}
      {result.ok && (
        <>
          <ul className="step-list">
            {result.steps.map((s) => (
              <li key={s.id}>
                <span className={`kind-badge ${s.kind}`}>{s.kind}</span>
                {s.label}
                {s.skip_if_exists && <span className="muted"> (skips if output exists)</span>}
                <div className="step-cmd">$ {s.command}</div>
              </li>
            ))}
          </ul>
          <button className="primary" disabled={running} onClick={onRun} style={{ marginTop: 8, width: '100%' }}>
            {running ? 'Starting…' : '▶ Run this plan'}
          </button>
        </>
      )}
    </div>
  )
}

let nodeSeq = 1
let groupSeq = 1
let pasteCount = 0
let canvasClipboard = null

function bumpNodeSeqFromIds(ids) {
  for (const id of ids) {
    const m = /_(\d+)$/.exec(String(id))
    if (m) nodeSeq = Math.max(nodeSeq, parseInt(m[1], 10) + 1)
  }
}

function isOutputHandle(node, handleId) {
  if (!node || !handleId) return false
  if (node.type === 'groupBox') return handleId === 'out'
  if (node.type === 'alphabet') return /^(prop|proj)_/.test(handleId)
  const def = NODE_DEFS[node.type]
  return !!(def?.outputs || []).some((o) => o.id === handleId)
}

function GroupNode({ data }) {
  return (
    <div className="node-card group-node">
      <div className="node-title">▣ {data.name || 'Group'}</div>
      <div className="node-sub muted">{data.count || 0} block{(data.count || 0) === 1 ? '' : 's'} collapsed</div>
      <Handle type="target" position={Position.Left} id="in" style={{ top: '50%' }} />
      <Handle type="source" position={Position.Right} id="out" style={{ top: '50%' }} />
    </div>
  )
}

function groupBoxNode(g, extra = {}) {
  return { id: g.id, type: 'groupBox', position: g.position, data: { group_id: g.id, name: g.name, count: g.node_ids.length }, ...extra }
}

function FlowEditorInner() {
  const { fid } = useParams()
  const { project, refreshProject } = useProject()
  const toast = useToast()
  const navigate = useNavigate()
  const [nodes, setNodes, onNodesChange] = useNodesState([])
  const [edges, setEdges, onEdgesChange] = useEdgesState([])
  const [selectedId, setSelectedId] = useState(null)
  const [compileResult, setCompileResult] = useState(null)
  const [sideTab, setSideTab] = useState('inspector')
  const [running, setRunning] = useState(false)
  const [flowName, setFlowName] = useState('')
  const [saveState, setSaveState] = useState('saved')
  const [groups, setGroups] = useState([])
  const [selCount, setSelCount] = useState(0)
  const rf = useRef(null)
  const wrapper = useRef(null)
  const skipAutosave = useRef(true)

  const flow = project?.flows.find((f) => f.id === fid)

  useEffect(() => {
    if (flow) {
      skipAutosave.current = true
      setFlowName(flow.name)
      const gr0 = (flow.graph?.groups || []).map((g) => ({ collapsed: true, position: { x: 0, y: 0 }, ...g }))
      let ns = flow.graph?.nodes || []
      const seenIds = new Set()
      const renamedIds = new Map()
      ns = ns.map((n) => {
        if (!seenIds.has(n.id)) { seenIds.add(n.id); return n }
        bumpNodeSeqFromIds(seenIds)
        let nid = `${n.type}_${nodeSeq++}`
        while (seenIds.has(nid)) nid = `${n.type}_${nodeSeq++}`
        seenIds.add(nid)
        renamedIds.set(n.id, nid)
        return { ...n, id: nid }
      })
      const gr = renamedIds.size
        ? gr0.map((g) => ({ ...g, node_ids: g.node_ids.map((x) => renamedIds.get(x) || x) }))
        : gr0
      if (renamedIds.size) skipAutosave.current = false
      for (const g of gr) {
        const num = parseInt(String(g.id).replace(/\D/g, ''), 10)
        if (!Number.isNaN(num)) groupSeq = Math.max(groupSeq, num + 1)
        if (!g.collapsed) continue
        const ids = new Set(g.node_ids)
        ns = ns.map((n) => (ids.has(n.id) ? { ...n, hidden: true } : n))
        ns = ns.concat(groupBoxNode(g))
      }
      setGroups(gr)
      setNodes(ns)
      setEdges(flow.graph?.edges || [])
      bumpNodeSeqFromIds([...ns.map((n) => n.id), ...(flow.graph?.edges || []).map((e) => e.id)])
      setCompileResult(null)
      setSaveState('saved')
    }
  }, [flow?.id]) // eslint-disable-line

  const serializeGraph = useCallback(() => ({
    nodes: nodes.filter((n) => n.type !== 'groupBox').map((n) => ({ id: n.id, type: n.type, position: n.position, data: n.data })),
    edges: edges.map((e) => ({ id: e.id, source: e.source, sourceHandle: e.sourceHandle, target: e.target, targetHandle: e.targetHandle })),
    groups: groups.filter((g) => g.node_ids.some((id) => nodes.some((n) => n.id === id && n.type !== 'groupBox'))),
  }), [nodes, edges, groups])

  useEffect(() => {
    if (!flow) return undefined
    if (skipAutosave.current) { skipAutosave.current = false; return undefined }
    setSaveState('unsaved')
    const t = setTimeout(async () => {
      setSaveState('saving')
      try {
        await api.updateFlow(project.id, fid, { name: flowName, graph: serializeGraph() })
        setSaveState('saved')
      } catch {
        setSaveState('error')
      }
    }, 800)
    return () => clearTimeout(t)
  }, [nodes, edges, flowName, groups]) // eslint-disable-line

  const resolveSourceKind = useCallback((nodeId, handleId, depth = 0) => {
    const node = nodes.find((n) => n.id === nodeId)
    if (!node) return null
    if (handleId === 'out' && (node.type === 'add_tensors' || node.type === 'tensor_join' || node.type === 'ternary_contract' || node.type === 'apply_symmetry' || node.type === 'assemble')) {
      if (depth > 8) return 'tensor'
      const handle = (node.type === 'ternary_contract' || node.type === 'apply_symmetry') ? 'tensor' : 'a'
      const e = node.type === 'assemble'
        ? edges.find((ed) => ed.target === nodeId && (ed.targetHandle || '').startsWith('in_'))
        : edges.find((ed) => ed.target === nodeId && (ed.targetHandle === handle || (handle === 'a' && ed.targetHandle === 'b')))
      return e ? resolveSourceKind(e.source, e.sourceHandle, depth + 1) : 'tensor'
    }
    return sourceKindFor(node, handleId, project)
  }, [nodes, edges, project])

  const normalizeConn = useCallback((conn) => {
    const srcNode = nodes.find((n) => n.id === conn.source)
    if (isOutputHandle(srcNode, conn.sourceHandle)) return conn
    return { source: conn.target, sourceHandle: conn.targetHandle, target: conn.source, targetHandle: conn.sourceHandle }
  }, [nodes])

  const isValidConnection = useCallback((conn) => {
    const c = normalizeConn(conn)
    const tn = nodes.find((n) => n.id === c.target)
    const sk = resolveSourceKind(c.source, c.sourceHandle)
    const tk = targetKindFor(tn, c.targetHandle)
    return !!(sk && tk && kindsCompatible(sk, tk))
  }, [nodes, resolveSourceKind, normalizeConn])

  const onConnect = useCallback((conn) => {
    if (!isValidConnection(conn)) {
      toast('Those ports are not compatible.')
      return
    }
    const c = normalizeConn(conn)
    setEdges((eds) => addEdge(
      { ...c, animated: false },
      eds.filter((ed) => !(ed.target === c.target && (ed.targetHandle || null) === (c.targetHandle || null))),
    ))
    setCompileResult(null)
  }, [isValidConnection, normalizeConn, setEdges, toast])

  const onDrop = useCallback((e) => {
    e.preventDefault()
    const type = e.dataTransfer.getData('application/symbology-node')
    if (!type || !rf.current) return
    const pos = rf.current.screenToFlowPosition({ x: e.clientX, y: e.clientY })
    const defaults = type === 'alphabet' ? { alphabet_id: '', selected_properties: [] } : {}
    setNodes((ns) => {
      const used = new Set(ns.map((n) => n.id))
      let id = `${type}_${nodeSeq++}`
      while (used.has(id)) id = `${type}_${nodeSeq++}`
      return [...ns, { id, type, position: pos, data: defaults }]
    })
    setCompileResult(null)
  }, [setNodes])

  const save = async () => {
    setSaveState('saving')
    await api.updateFlow(project.id, fid, { name: flowName, graph: serializeGraph() })
    await refreshProject()
    setSaveState('saved')
    toast('Flow saved.')
  }

  const compile = async () => {
    try {
      await save()
    } catch (e) {
      toast('Save failed: ' + e.message)
      return
    }
    try {
      const res = await api.compileFlow(project.id, fid)
      setCompileResult(res)
      setSideTab('plan')
      if (!res.ok) toast('The flow does not compile — see the errors panel.')
    } catch (e) { toast(e.message) }
  }

  const run = async () => {
    setRunning(true)
    try {
      const { run_id } = await api.runFlow(project.id, fid)
      navigate(`/runs/${run_id}`)
    } catch (e) {
      setRunning(false)
      const errs = e.payload?.detail?.errors
      if (errs) { setCompileResult({ ok: false, errors: errs }); setSideTab('plan') }
      toast('Run failed to start: ' + (errs ? 'see errors panel' : e.message))
    }
  }

  const updateNodeData = useCallback((id, patch) => {
    setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, ...patch } } : n)))
    if ('selected_properties' in patch || 'alphabet_id' in patch) {
      const keep = new Set()
      for (const pid of patch.selected_properties || []) { keep.add(`prop_${pid}`); keep.add(`proj_${pid}`) }
      setEdges((es) => es.filter((e) => e.source !== id || !/^(prop|proj)_/.test(e.sourceHandle || '') || keep.has(e.sourceHandle)))
    }
    setCompileResult(null)
  }, [setNodes, setEdges])

  const groupOf = useMemo(() => {
    const m = new Map()
    for (const g of groups) if (g.collapsed) for (const id of g.node_ids) m.set(id, g.id)
    return m
  }, [groups])

  const visibleEdges = useMemo(() => {
    const out = []
    for (const e of edges) {
      const gs = groupOf.get(e.source)
      const gt = groupOf.get(e.target)
      if (gs && gt && gs === gt) continue
      if (!gs && !gt) { out.push(e); continue }
      out.push({
        ...e,
        id: `px_${e.id}`,
        source: gs || e.source,
        sourceHandle: gs ? 'out' : e.sourceHandle,
        target: gt || e.target,
        targetHandle: gt ? 'in' : e.targetHandle,
        selectable: false,
        focusable: false,
        style: { ...e.style, strokeDasharray: '6 3', opacity: 0.7 },
      })
    }
    return out
  }, [edges, groupOf])

  const groupSelection = useCallback(() => {
    const members = nodes.filter((n) => n.selected && n.type !== 'groupBox' && !n.hidden)
    if (members.length < 2) return
    const gid = `group_${groupSeq++}`
    const minX = Math.min(...members.map((n) => n.position.x))
    const minY = Math.min(...members.map((n) => n.position.y))
    const memberIds = members.map((n) => n.id)
    const g = { id: gid, name: `Group ${groupSeq - 1}`, node_ids: memberIds, collapsed: true, position: { x: minX, y: minY } }
    setNodes((ns) => ns
      .map((n) => (memberIds.includes(n.id) ? { ...n, hidden: true, selected: false } : n))
      .concat(groupBoxNode(g, { selected: true })))
    setGroups((gs) => [...gs, g])
    setSelectedId(gid)
    setCompileResult(null)
  }, [nodes, setNodes])

  const expandGroup = useCallback((gid, dissolve = false) => {
    const g = groups.find((x) => x.id === gid)
    if (!g) return
    const ids = new Set(g.node_ids)
    setNodes((ns) => ns.filter((n) => n.id !== gid).map((n) => (ids.has(n.id) ? { ...n, hidden: undefined } : n)))
    setGroups((gs) => (dissolve ? gs.filter((x) => x.id !== gid) : gs.map((x) => (x.id === gid ? { ...x, collapsed: false } : x))))
    setSelectedId(null)
    setCompileResult(null)
  }, [groups, setNodes])

  const collapseGroup = useCallback((gid) => {
    const g = groups.find((x) => x.id === gid)
    if (!g) return
    const ids = new Set(g.node_ids)
    const members = nodes.filter((n) => ids.has(n.id))
    if (!members.length) return
    const minX = Math.min(...members.map((n) => n.position.x))
    const minY = Math.min(...members.map((n) => n.position.y))
    setNodes((ns) => ns
      .map((n) => (ids.has(n.id) ? { ...n, hidden: true, selected: false } : n))
      .concat(groupBoxNode({ ...g, position: { x: minX, y: minY } })))
    setGroups((gs) => gs.map((x) => (x.id === gid ? { ...x, collapsed: true, position: { x: minX, y: minY } } : x)))
    setCompileResult(null)
  }, [groups, nodes, setNodes])

  const renameGroup = useCallback((gid, name) => {
    setGroups((gs) => gs.map((x) => (x.id === gid ? { ...x, name } : x)))
    setNodes((ns) => ns.map((n) => (n.id === gid && n.type === 'groupBox' ? { ...n, data: { ...n.data, name } } : n)))
  }, [setNodes])

  const onNodeDragStop = useCallback((e, node) => {
    if (node.type !== 'groupBox') return
    const g = groups.find((x) => x.id === node.id)
    if (!g || !g.position) return
    const dx = node.position.x - g.position.x
    const dy = node.position.y - g.position.y
    if (!dx && !dy) return
    const ids = new Set(g.node_ids)
    setNodes((ns) => ns.map((n) => (ids.has(n.id) ? { ...n, position: { x: n.position.x + dx, y: n.position.y + dy } } : n)))
    setGroups((gs) => gs.map((x) => (x.id === node.id ? { ...x, position: { x: node.position.x, y: node.position.y } } : x)))
  }, [groups, setNodes])

  const copySelection = useCallback(() => {
    const sel = nodes.filter((n) => n.selected && n.type !== 'groupBox' && !n.hidden)
    if (!sel.length) return false
    const ids = new Set(sel.map((n) => n.id))
    canvasClipboard = {
      nodes: sel.map((n) => ({ id: n.id, type: n.type, position: { ...n.position }, data: JSON.parse(JSON.stringify(n.data || {})) })),
      edges: edges.filter((ed) => ids.has(ed.source) && ids.has(ed.target))
        .map((ed) => ({ id: ed.id, source: ed.source, sourceHandle: ed.sourceHandle, target: ed.target, targetHandle: ed.targetHandle })),
    }
    return true
  }, [nodes, edges])

  const pasteClipboard = useCallback(() => {
    if (!canvasClipboard?.nodes?.length) return
    pasteCount += 1
    const off = 40 * pasteCount
    const used = new Set([...nodes.map((n) => n.id), ...edges.map((e) => e.id)])
    const fresh = (prefix) => {
      let id = `${prefix}_${nodeSeq++}`
      while (used.has(id)) id = `${prefix}_${nodeSeq++}`
      used.add(id)
      return id
    }
    const idMap = new Map()
    const newNodes = canvasClipboard.nodes.map((n) => {
      const nid = fresh(n.type)
      idMap.set(n.id, nid)
      return { id: nid, type: n.type, position: { x: n.position.x + off, y: n.position.y + off }, data: JSON.parse(JSON.stringify(n.data)), selected: true }
    })
    const newEdges = canvasClipboard.edges
      .filter((ed) => idMap.has(ed.source) && idMap.has(ed.target))
      .map((ed) => ({ id: fresh('edge'), source: idMap.get(ed.source), sourceHandle: ed.sourceHandle, target: idMap.get(ed.target), targetHandle: ed.targetHandle }))
    setNodes((ns) => ns.map((n) => ({ ...n, selected: false })).concat(newNodes))
    setEdges((es) => es.concat(newEdges))
    setCompileResult(null)
    toast(`Pasted ${newNodes.length} node${newNodes.length === 1 ? '' : 's'}.`)
  }, [nodes, edges, setNodes, setEdges, toast])

  const deleteNode = useCallback((id) => {
    const node = nodes.find((n) => n.id === id)
    if (node?.type === 'groupBox') { expandGroup(id); return }
    setNodes((ns) => ns.filter((n) => n.id !== id))
    setEdges((es) => es.filter((e) => e.source !== id && e.target !== id))
    setGroups((gs) => gs.map((g) => ({ ...g, node_ids: g.node_ids.filter((x) => x !== id) })))
    setSelectedId(null)
    setCompileResult(null)
    toast('Node deleted.')
  }, [nodes, setNodes, setEdges, expandGroup, toast])

  useEffect(() => {
    const handler = (e) => {
      const t = e.target
      if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable)) return
      const mod = e.metaKey || e.ctrlKey
      if (mod && (e.key === 'c' || e.key === 'C')) {
        if (copySelection()) toast('Selection copied.')
        e.preventDefault()
        return
      }
      if (mod && (e.key === 'v' || e.key === 'V')) { pasteClipboard(); e.preventDefault(); return }
      if (mod && (e.key === 'd' || e.key === 'D')) { if (copySelection()) pasteClipboard(); e.preventDefault(); return }
      if (e.key !== 'Delete' && e.key !== 'Backspace') return
      const selected = nodes.filter((n) => n.selected)
      const groupIds = selected.filter((n) => n.type === 'groupBox').map((n) => n.id)
      const nodeIds = new Set(selected.filter((n) => n.type !== 'groupBox').map((n) => n.id))
      const edgeIds = new Set(edges.filter((ed) => ed.selected).map((ed) => ed.id))
      if (!nodeIds.size && !edgeIds.size && !groupIds.length) return
      e.preventDefault()
      for (const gid of groupIds) expandGroup(gid)
      setNodes((ns) => ns.filter((n) => !nodeIds.has(n.id)))
      setEdges((es) => es.filter((ed) => !edgeIds.has(ed.id) && !nodeIds.has(ed.source) && !nodeIds.has(ed.target)))
      setGroups((gs) => gs.map((g) => ({ ...g, node_ids: g.node_ids.filter((x) => !nodeIds.has(x)) })))
      setSelectedId(null)
      setCompileResult(null)
    }
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [nodes, edges, setNodes, setEdges, copySelection, pasteClipboard, expandGroup, toast])

  const selectedNode = nodes.find((n) => n.id === selectedId)

  if (!flow) {
    if (project) return (
      <div className="page">
        <div className="empty-state">
          <div className="big">Flow not found</div>
          <p>This flow does not exist in project “{project.name}” — it may belong to another project.</p>
          <button className="primary" onClick={() => navigate('/flows')}>Back to flows</button>
        </div>
      </div>
    )
    return <div className="page"><p className="muted">Loading flow…</p></div>
  }

  return (
    <div className="flow-layout">
      <div className="flow-palette">
        {PALETTE_SECTIONS.map((sec) => (
          <div key={sec.title}>
            <div className="muted" style={{ fontSize: 10, margin: '10px 0 4px', textTransform: 'uppercase', letterSpacing: 0.5 }}>{sec.title}</div>
            {sec.items.map((p) => (
              <div
                key={p.type} className="palette-item" draggable
                onDragStart={(e) => e.dataTransfer.setData('application/symbology-node', p.type)}
              >
                {p.label}
                <div className="sub">{p.sub}</div>
              </div>
            ))}
          </div>
        ))}
        <div className="muted" style={{ fontSize: 11, marginTop: 12 }}>Drag a node onto the canvas. Wire colored ports of the same kind together.</div>
      </div>
      <div className="flow-canvas" ref={wrapper} onDrop={onDrop} onDragOver={(e) => e.preventDefault()}>
        <div className="flow-toolbar">
          <input style={{ width: 180 }} value={flowName} onChange={(e) => setFlowName(e.target.value)} />
          <span className="muted shrink" style={{ fontSize: 11, alignSelf: 'center', minWidth: 90 }}>
            {saveState === 'saved' && '✓ saved'}
            {saveState === 'saving' && 'saving…'}
            {saveState === 'unsaved' && 'unsaved changes'}
            {saveState === 'error' && <span className="error-text">save failed</span>}
          </span>
          <button onClick={save}>Save</button>
          <button onClick={compile}>Compile</button>
          <button className="primary" onClick={compile}>Run…</button>
        </div>
        <ReactFlow
          nodes={nodes} edges={visibleEdges}
          onNodesChange={(chg) => { onNodesChange(chg); setCompileResult(null) }}
          onEdgesChange={(chg) => { onEdgesChange(chg); setCompileResult(null) }}
          onConnect={onConnect}
          isValidConnection={isValidConnection}
          nodeTypes={nodeTypes}
          onInit={(inst) => { rf.current = inst }}
          onSelectionChange={(sel) => {
            setSelectedId(sel.nodes[0]?.id || null)
            setSelCount(sel.nodes.filter((n) => n.type !== 'groupBox').length)
          }}
          onNodeDragStop={onNodeDragStop}
          deleteKeyCode={null}
            selectionOnDrag
            panOnDrag={[1, 2]}
            selectionKeyCode="Shift"
            multiSelectionKeyCode={['Meta', 'Control']}
          fitView
        >
          <Background gap={18} />
          <Controls />
        </ReactFlow>
      </div>
      <div className="flow-side">
        <div className="row" style={{ marginBottom: 12 }}>
          <button className={sideTab === 'inspector' ? 'primary' : ''} onClick={() => setSideTab('inspector')}>Inspector</button>
          <button className={sideTab === 'plan' ? 'primary' : ''} onClick={() => setSideTab('plan')}>Plan</button>
        </div>
        {sideTab === 'inspector' && (
          <Inspector
            node={selectedNode}
            onChange={updateNodeData}
            onDelete={deleteNode}
            groupOps={{ groups, onExpand: (id) => expandGroup(id), onUngroup: (id) => expandGroup(id, true), onCollapse: collapseGroup, onRename: renameGroup }}
          />
        )}
        {sideTab === 'plan' && (
          compileResult
            ? <CompilePanel result={compileResult} onRun={run} running={running} />
            : <p className="muted">Press “Compile” to turn the diagram into the exact command sequence. You can inspect every command before running it.</p>
        )}
      </div>
    </div>
  )
}

export default function FlowEditor() {
  return (
    <ReactFlowProvider>
      <FlowEditorInner />
    </ReactFlowProvider>
  )
}
