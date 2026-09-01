import React, { useCallback, useContext, useEffect, useMemo, useRef, useState } from 'react'
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
  kindsCompatible, sourceKindFor, targetKindFor, alphabetOutSlots,
  registerCustomBlocks, customBlockPalette,
} from '../flowdefs'

const ROW0 = 29
const ROW_STEP = 15
const SUB_EXTRA = 15

const DimsContext = React.createContext(() => null)

function propDimsText(prop) {
  const dims = prop?.summary?.dims
  return Array.isArray(dims) && dims.length ? dims.join('×') : null
}

function projDimsText(prop) {
  const dims = (prop?.summary?.dims || []).filter((d) => d !== 1)
  return dims.length ? dims.join('×') : null
}

function alphabetHandleDims(project, node, handleId) {
  const alphabet = project?.alphabets.find((a) => a.id === node.data?.alphabet_id)
  const prop = alphabet?.properties.find((p) => handleId === `prop_${p.id}` || handleId === `proj_${p.id}`)
  if (!prop) return null
  return handleId.startsWith('proj_') ? projDimsText(prop) : propDimsText(prop)
}

function useHandleDims() {
  return useContext(DimsContext)
}

function StatusDot({ status }) {
  return <span className={`dot ${status || 'pending'}`} />
}

function AlphabetNode({ id, data, selected }) {
  const { project } = useProject()
  const alphabet = project?.alphabets.find((a) => a.id === data.alphabet_id)
  const selectedProps = data.selected_properties || []
  const outRows = alphabetOutSlots(alphabet, selectedProps).map((slot) => {
    const p = alphabet.properties.find((x) => x.id === slot.propId)
    return {
      ...slot,
      label: slot.proj ? `${propLabel(p)} · proj map` : propLabel(p),
      status: p.status,
      dims: slot.proj ? projDimsText(p) : propDimsText(p),
    }
  })
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-body">
        {!alphabet && <span>select property in inspector →</span>}
        {alphabet && selectedProps.length === 0 && <span>no output selected</span>}
        {alphabet && (
          <div style={{ fontWeight: 600, marginBottom: 4, paddingBottom: 4, borderBottom: '1px solid var(--border)' }}>
            {alphabet.name}
          </div>
        )}
        {outRows.map((r, i) => (
          <div key={r.key} className="handle-row" style={{ textAlign: 'right' }}>
            <span className={`kind-${r.kind}`}>
              {r.label}
            </span>{' '}
            {r.dims && <span className="dim-badge">{r.dims}</span>}{' '}
            <StatusDot status={r.status} />
            <Handle
              type="source" position={Position.Right} id={r.handle}
              style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -6, background: 'var(--accent)' }}
            />
          </div>
        ))}
      </div>
    </div>
  )
}

function OpNode({ id, type, data, selected }) {
  const def = NODE_DEFS[type]
  const handleDims = useHandleDims()
  const allEdges = useEdges()
  const updateNodeInternals = useUpdateNodeInternals()
  const dynamicInputs = type === 'add_tensors'
  const dynamicPairs = type === 'solve_collinear'
  let inputs = def.inputs || []
  if (dynamicInputs) {
    let maxConnected = 1
    for (const e of allEdges) {
      if (e.target === id && (e.targetHandle || '').startsWith('in_')) {
        const idx = parseInt(e.targetHandle.slice(3).split('@')[0], 10)
        if (!Number.isNaN(idx)) maxConnected = Math.max(maxConnected, idx)
      }
    }
    inputs = Array.from({ length: Math.max(inputs.length, maxConnected + 2) }, (_, i) => ({
      id: `in_${i}`, kind: 'tensor', label: String.fromCharCode(65 + i),
    }))
  }
  let pairCount = 0
  let hasCondEdge = false
  if (dynamicPairs) {
    // Extra {seed, rhs} pairs beyond the fixed seed/rhs ports (pair 0).
    // Each pair N >= 1 adds in_seed_N / in_rhs_N; one spare pair is always
    // shown so the next pair can be wired without Inspector round-trips.
    let maxPair = 0
    for (const e of allEdges) {
      if (e.target !== id) continue
      const m = (e.targetHandle || '').match(/^in_(?:seed|rhs)_(\d+)$/)
      if (m) maxPair = Math.max(maxPair, parseInt(m[1], 10))
    }
    hasCondEdge = allEdges.some((e) => e.target === id && (e.targetHandle || '') === 'cond')
    pairCount = maxPair + 1
    const nLetters = (data.pair_letters || '').split(',').map((s) => s.trim()).filter(Boolean).length
    const nPairs = Math.max(maxPair, nLetters - 1) + 1
    const extras = []
    for (let p = 1; p <= nPairs; p++) {
      extras.push({ id: `in_seed_${p}`, kind: 'seed_or_tensor', label: `seed ${p + 1}` })
      extras.push({ id: `in_rhs_${p}`, kind: 'boundary', label: `rhs ${p + 1}` })
    }
    inputs = [...inputs, ...extras]
  }
  useEffect(() => { if (dynamicInputs || dynamicPairs) updateNodeInternals(id) }, [id, dynamicInputs, dynamicPairs, inputs.length, updateNodeInternals])
  const outputs = def.outputs || []
  const rows = Math.max(inputs.length, outputs.length)
  const cbWeights = (data.weights || '').split(',').map((s) => s.trim())
  const subtitle =
    type === 'extend' ? (data.target_weight ? `→ weight ${data.target_weight}` : 'weight +1')
    : (type === 'project' || type === 'solve_symmetry' || type === 'symderive') ? `→ ${data.target || '…'}`
    : type === 'sew' ? `→ ${data.target || 'SEW_FpL'}`
    : type === 'solve_collinear' ? (pairCount > 0 || hasCondEdge
        ? `${pairCount || '—'} pair${pairCount === 1 ? '' : 's'}${hasCondEdge ? ' + cond' : ''} → ${data.out_stem || '…'}`
        : `${data.target || '…'}`)
    : type === 'projection_chain' ? `${data.symmetry || '?'} · ${data.target || '…'}`
    : type === 'symmetry_invariant' ? `${data.symmetry || '?'} · ${data.target || '…'}`
    : type === 'compute_rhs' ? `${data.target || '…'}`
    : type === 'add_tensors' ? `${inputs.map((_, i) => `${cbWeights[i] || '1'}·${String.fromCharCode(65 + i)}`).join(' + ')} → ${data.target || '…'}`
    : (type === 'ternary_contract' || type === 'apply_symmetry') ? `→ ${data.target || '…'}`
    : type === 'matrix_power' ? `M^${data.n || '?'} → ${data.target || '…'}`
    : type === 'tensor_join' ? `axis ${data.axis || '?'} → ${data.target || '…'}`
    : type === 'tensor_dot' ? `A[${data.axis_a ?? '?'}]·B[${data.axis_b ?? '?'}] → ${data.target || '…'}`
    : type === 'impose_integrability' ? `${data.transpose === false ? 'relations among conditions' : 'solve coefficients'} → ${data.target || '…'}`
    : type === 'integrability_condition' ? `M[(a), b·d] = Σ S·dlog → ${data.target || '…'}`
    : type === 'solve_conditions' ? `${data.transpose === false ? 'relations among conditions' : 'solve coefficients'} → ${data.target || '…'}`
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
              {inputs[i] && handleDims(id, inputs[i].id) && <span className="dim-badge"> {handleDims(id, inputs[i].id)}</span>}
            </span>
            <span className={outputs[i] ? `kind-${outputs[i].kind}` : ''}>
              {outputs[i]?.label || ''}
              {outputs[i] && handleDims(id, outputs[i].id) && <span className="dim-badge"> {handleDims(id, outputs[i].id)}</span>}
            </span>
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
  const handleDims = useHandleDims()
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
  const outputs = data.outputs || []
  const rows = Math.max(slots, outputs.length)
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        <div style={{ marginBottom: 2 }}>{`A = Σc·B → ${data.target || '…'}`}</div>
        {Array.from({ length: rows }).map((_, i) => (
          <div key={i} className="handle-row" style={{ display: 'flex', justifyContent: 'space-between', gap: 12 }}>
            {i < slots ? (
              <span className="kind-tensor">
                elem {i + 1}
                {handleDims(id, `in_${i}`) && <span className="dim-badge"> {handleDims(id, `in_${i}`)}</span>}
              </span>
            ) : <span />}
            {outputs[i] && (
              <span className="kind-tensor">
                {handleDims(id, `out_${i}`) && <span className="dim-badge">{handleDims(id, `out_${i}`)} </span>}
                {outputs[i].name || `out ${i + 1}`} →
              </span>
            )}
          </div>
        ))}
        {Array.from({ length: slots }).map((_, i) => (
          <Handle
            key={`in_${i}`} type="target" position={Position.Left} id={`in_${i}`}
            style={{ top: ROW0 + SUB_EXTRA + i * ROW_STEP }}
          />
        ))}
        {outputs.map((o, oi) => (
          <Handle
            key={`out_${oi}`} type="source" position={Position.Right} id={`out_${oi}`}
            style={{ top: ROW0 + SUB_EXTRA + oi * ROW_STEP, background: 'var(--accent)' }}
          />
        ))}
      </div>
    </div>
  )
}

function CBIONode({ id, data, selected, type }) {
  const def = NODE_DEFS[type]
  const isIn = type === 'cb_in'
  const label = (data?.name || '').trim() || (isIn ? 'in' : 'out')
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined, borderStyle: 'dashed' }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        <div className="handle-row" style={{ position: 'relative', display: 'flex', justifyContent: isIn ? 'flex-end' : 'flex-start' }}>
          <span className={`kind-${data?.kind || 'tensor'}`}>
            {label}
          </span>
          {isIn ? (
            <Handle type="source" position={Position.Right} id="out"
              style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -18, background: 'var(--accent)' }} />
          ) : (
            <Handle type="target" position={Position.Left} id="in"
              style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', left: -18 }} />
          )}
        </div>
      </div>
    </div>
  )
}

function CustomBlockNode({ id, data, selected }) {
  const { project } = useProject()
  const def = NODE_DEFS[`cb_${data?.block}`]
  if (!def) {
    return (
      <div className="node-card" style={{ borderColor: '#c33', opacity: 0.8 }}>
        <div className="node-title" style={{ background: '#fdd' }}>Custom block</div>
        <div className="node-body"><span className="error-text">missing block — deleted?</span></div>
      </div>
    )
  }
  const rows = Math.max(def.inputs.length, def.outputs.length, 1)
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        {Array.from({ length: rows }).map((_, i) => (
          <div key={i} className="handle-row" style={{ display: 'flex', justifyContent: 'space-between' }}>
            <span className={def.inputs[i] ? `kind-${def.inputs[i].kind}` : ''}>{def.inputs[i]?.label || ''}</span>
            <span className={def.outputs[i] ? `kind-${def.outputs[i].kind}` : ''}>
              {def.outputs[i]?.label || ''}
              {' →'}
            </span>
          </div>
        ))}
        {def.inputs.map((inp, i) => (
          <Handle key={inp.id} type="target" position={Position.Left} id={inp.id}
            style={{ top: ROW0 + i * ROW_STEP }} />
        ))}
        {def.outputs.map((out, i) => (
          <Handle key={out.id} type="source" position={Position.Right} id={out.id}
            style={{ top: ROW0 + i * ROW_STEP, background: 'var(--accent)' }} />
        ))}
      </div>
    </div>
  )
}

function ReuseOutputNode({ data, selected }) {
  const def = NODE_DEFS.reuse_output
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined, borderStyle: 'dashed' }}>
      <div className="node-title" style={{ background: def.color }}>{def.title}</div>
      <div className="node-body">
        <div className="handle-row" style={{ position: 'relative', display: 'flex', justifyContent: 'flex-end' }}>
          <span className={`kind-${data?.kind || 'tensor'}`}>
            {data?.name || data?.file || 'select output in inspector →'}
          </span>
          <Handle
            type="source" position={Position.Right} id="out"
            style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -6, background: 'var(--accent)' }}
          />
        </div>
      </div>
    </div>
  )
}

const nodeTypes = {
  alphabet: AlphabetNode,
  cb_in: CBIONode,
  cb_out: CBIONode,
  reuse_output: ReuseOutputNode,
  customblock: CustomBlockNode,
  merge_conditions: (p) => <OpNode {...p} type="merge_conditions" />,
  extend: (p) => <OpNode {...p} type="extend" />,
  sew: (p) => <OpNode {...p} type="sew" />,
  project: (p) => <OpNode {...p} type="project" />,
  solve_symmetry: (p) => <OpNode {...p} type="solve_symmetry" />,
  symderive: (p) => <OpNode {...p} type="symderive" />,
  solve_collinear: (p) => <OpNode {...p} type="solve_collinear" />,
  projection_chain: (p) => <OpNode {...p} type="projection_chain" />,
  symmetry_invariant: (p) => <OpNode {...p} type="symmetry_invariant" />,
  compute_rhs: (p) => <OpNode {...p} type="compute_rhs" />,
  add_tensors: (p) => <OpNode {...p} type="add_tensors" />,
  ternary_contract: (p) => <OpNode {...p} type="ternary_contract" />,
  apply_symmetry: (p) => <OpNode {...p} type="apply_symmetry" />,
  matrix_power: (p) => <OpNode {...p} type="matrix_power" />,
  tensor_join: (p) => <OpNode {...p} type="tensor_join" />,
  tensor_dot: (p) => <OpNode {...p} type="tensor_dot" />,
  impose_integrability: (p) => <OpNode {...p} type="impose_integrability" />,
  integrability_condition: (p) => <OpNode {...p} type="integrability_condition" />,
  solve_conditions: (p) => <OpNode {...p} type="solve_conditions" />,
  assemble: AssembleNode,
  groupBox: GroupNode,
}

function ReuseOutputInspector({ node, set, onDelete, project, currentFlowId }) {
  const [catalog, setCatalog] = useState(null)
  const [err, setErr] = useState(null)
  useEffect(() => {
    let alive = true
    setCatalog(null)
    setErr(null)
    api.flowOutputs(project.id)
      .then((c) => { if (alive) setCatalog(c) })
      .catch((e) => { if (alive) setErr(String(e.message || e)) })
    return () => { alive = false }
  }, [project.id, project.flows])
  const d = node.data || {}
  return (
    <div>
      <h3>Reuse Output</h3>
      <p className="muted" style={{ fontSize: 11 }}>
        Pick a calculated output of another flow in this project. The tensor file must already exist
        (run that flow first); it is then usable here exactly like any other tensor source.
      </p>
      {err && <p className="error-text">{err}</p>}
      {catalog === null && !err && <p className="muted">Loading outputs…</p>}
      {catalog?.length === 0 && (
        <p className="muted">No compiled flow with outputs yet. Compile or run another flow first.</p>
      )}
      {(catalog || []).filter((g) => g.flow_id !== currentFlowId).map((g) => (
        <div key={g.flow_id} style={{ marginTop: 10 }}>
          <div className="muted" style={{ fontSize: 11 }}>{g.flow_name}</div>
          {g.outputs.map((o) => (
            <label key={o.file} style={{ display: 'flex', alignItems: 'center', gap: 6, margin: '2px 0', cursor: 'pointer' }}>
              <input
                type="radio" name={`reuse_${node.id}`} style={{ width: 'auto' }}
                checked={d.file === o.file}
                onChange={() => set({ file: o.file, kind: o.kind, name: o.name })}
              />
              <span className={`kind-${o.kind}`}>{o.name}</span>
              <span className="muted" style={{ fontSize: 10 }}>{o.kind}</span>
            </label>
          ))}
        </div>
      ))}
      {d.file && (
        <>
          <label style={{ marginTop: 12 }}>Tensor kind override</label>
          <select value={d.kind || 'tensor'} onChange={(e) => set({ kind: e.target.value })}>
            {['tensor', 'matrix', 'dlogmat', 'fec', 'lec', 'fec1', 'lec1', 'sew', 'basis', 'solution', 'boundary'].map((k) => (
              <option key={k} value={k}>{k}</option>
            ))}
          </select>
          <p className="muted" style={{ fontSize: 11 }}>File: {d.file}</p>
        </>
      )}
      <div style={{ marginTop: 16, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
        <button className="danger" onClick={() => onDelete(node.id)}>Delete this node</button>
      </div>
    </div>
  )
}

function Inspector({ node, onChange, onDelete, groupOps, onOpenBlock, flowId }) {
  const { project } = useProject()
  const edges = useEdges()
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

  if (node.type === 'cb_in' || node.type === 'cb_out') {
    return (
      <div>
        <h3>{node.type === 'cb_in' ? 'Block Input port' : 'Block Output port'}</h3>
        <p className="muted" style={{ fontSize: 11 }}>
          Abstract port of the custom block being defined. Name and kind become the port of the sealed block.
        </p>
        <label>Port name</label>
        <input value={d.name || ''} onChange={(e) => set({ name: e.target.value })} placeholder="e.g. M / tensor / sym" />
        <label>Port kind</label>
        <select value={d.kind || 'tensor'} onChange={(e) => set({ kind: e.target.value })}>
          {['tensor', 'matrix', 'dlogmat', 'fec', 'lec', 'seed'].map((k) => <option key={k} value={k}>{k}</option>)}
        </select>
        <p className="muted" style={{ fontSize: 11 }}>
          Ports are ordered by placement on the card; the block exposes them in that order.
        </p>
        <button className="danger" onClick={() => onDelete(node.id)}>Delete node</button>
      </div>
    )
  }
  if (node.type === 'customblock') {
    const cb = project?.flows.find((f) => f.id === d.block)
    return (
      <div>
        <h3>Custom block: {cb?.name || '(missing)'}</h3>
        <p className="muted" style={{ fontSize: 11 }}>
          An instance of the sealed block “{cb?.name}”. Its internals are defined in that block&apos;s own
          diagram — edit them there; every instance updates automatically. Definition-level parameters
          (shared by all instances) are set inside the block.
        </p>
        <button onClick={() => onOpenBlock(d.block)}>Open block definition</button>
        <p className="muted" style={{ fontSize: 11, marginTop: 12 }}>
          Inputs: {(NODE_DEFS[`cb_${d.block}`]?.inputs || []).map((i) => `${i.label} (${i.kind})`).join(', ') || 'none'}
          <br />
          Outputs: {(NODE_DEFS[`cb_${d.block}`]?.outputs || []).map((o) => `${o.label} (${o.kind})`).join(', ') || 'none'}
        </p>
        <button className="danger" onClick={() => onDelete(node.id)}>Delete node</button>
      </div>
    )
  }
  if (node.type === 'reuse_output') {
    return <ReuseOutputInspector node={node} set={set} onDelete={onDelete} project={project} currentFlowId={flowId} />
  }
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
            <label>Exposed output (property)</label>
            {alphabet.properties.map((p) => (
              <div key={p.id} style={{ padding: '2px 0' }}>
                <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, margin: 0, color: 'var(--text)', cursor: 'pointer' }}>
                  <input
                    type="radio" name={`prop_${node.id}`} style={{ width: 'auto' }}
                    checked={(d.selected_properties || []).includes(p.id)}
                    onChange={() => set({ selected_properties: [p.id] })}
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
      {node.type === 'sew' && (
        <>
          <label>Target name (optional)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="auto: SEW_FpL" />
          <p className="muted" style={{ fontSize: 11 }}>Defaults to SEW_FpL from the FEC/LEC weights; a counter suffix (_2, _3, …) is appended automatically when multiple Sew nodes share the same name.</p>
        </>
      )}
      {(node.type === 'project' || node.type === 'solve_symmetry' || node.type === 'symderive') && (
        <>
          <label>Output name</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. FEC_proj" />
          <p className="muted" style={{ fontSize: 11 }}>
            {node.type === 'project'
              ? 'out[d,b\',c\'] = Σ P[d,a]·T[a,b,c]·S[b\',b]·S[c\',c] — wire the symmetry rep (square, axes 2,3) and the projection map (d×a, axis 1).'
              : node.type === 'symderive'
              ? 'Two steps: ① impose sym M on the last two entries of the extended tensor (s,a,b), ② solve R·T = T\' for the induced transformation matrix R (s×s). Feed R back as the sym M input of the next Symmetry Derive level. Wire sym M (c) only if the two axes use different matrices.'
              : 'Three steps: ① impose sym M on the last two entries (ternary contract), ② derive the induced R (a×a) from R·T = T\' and take K = ker(Rᵀ−I), ③ contract K with the original tensor → (e,b,c). Wire sym M (c) only if the two axes use different matrices.'}
          </p>
        </>
      )}
      {node.type === 'solve_collinear' && (
        <>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="SEW_3p1 (blank when seed is wired)" />
          <p className="muted" style={{ fontSize: 11 }}>
            Wire a tensor into the <b>seed</b> port to solve a custom seed (projection = none): the tensor is used
            as-is (no seed-space projection, no expansion unless --basis). Leave blank for named SEW/FEC targets.
          </p>
          <label>RHS file (or 0)</label>
          <input value={d.rhs || ''} onChange={(e) => set({ rhs: e.target.value })} placeholder="data/E1.wxf (or wire the rhs port)" />
          <label>Projection</label>
          <select value={d.projection || 'finite'} onChange={(e) => set({ projection: e.target.value })}>
            <option value="finite">finite</option>
            <option value="divergent">divergent</option>
            <option value="none">none (custom seed)</option>
          </select>
          <label>Letter projection (file or identity)</label>
          <input value={d.letter_projection || ''} onChange={(e) => set({ letter_projection: e.target.value })} placeholder="identity" />
          <p className="muted" style={{ fontSize: 11 }}>
            Applies to pair 0 (the seed/rhs ports). In multi-pair mode, entry 0 of “Pair letters” below
            overrides this field.
          </p>
          <label>Pair letters (per pair, comma-separated)</label>
          <input value={d.pair_letters || ''} onChange={(e) => set({ pair_letters: e.target.value })} placeholder="e.g. identity, output/collinear/colprojdiv_w1.wxf" />
          <p className="muted" style={{ fontSize: 11 }}>
            Each pair may use a different letter projection: entry 0 = pair 0 (seed/rhs ports above), entry N =
            seed N / rhs N (extra ports appear as pairs are wired). Each entry is a letter-projection file or
            “identity”. An unwired rhs N defaults to “0” (homogeneous constraints). Wiring any seed N port or the
            cond port switches the node to multi-pair mode, where every seed is used as-is (custom seed, no
            seed-space projection).
          </p>
          <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, marginTop: 8, color: 'var(--text)' }}>
            <input
              type="checkbox" style={{ width: 'auto' }}
              checked={!!d.export_conditions}
              onChange={(e) => set({ export_conditions: e.target.checked })}
            />
            Export conditions [M|r]
          </label>
          <p className="muted" style={{ fontSize: 11 }}>
            Writes the combined non-homogeneous constraints as a rank-2 [M | r] matrix (rhs = last column) next
            to the solution, and enables the conditions output port. Re-ingest it later via the cond port to
            combine constraints from different flows.
          </p>
          <label>Out stem (output naming)</label>
          <input value={d.out_stem || ''} onChange={(e) => set({ out_stem: e.target.value })} placeholder="auto: seed stem / stem1_xN" />
          <label>Solver</label>
          <select value={d.solver || 'incremental'} onChange={(e) => set({ solver: e.target.value })}>
            <option value="incremental">incremental</option>
            <option value="sampled">sampled</option>
          </select>
        </>
      )}
      {node.type === 'projection_chain' && (
        <>
          <label>Symmetry</label>
          <select value={d.symmetry || 'collinear'} onChange={(e) => set({ symmetry: e.target.value })}>
            {['collinear', 'cyclic', 'flip', 'parity'].map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="SEW_5p1" />
          <p className="muted" style={{ fontSize: 11 }}>
            Runs the full <code>bootstrap --project</code> pipeline for the chosen symmetry: builds the projection
            matrices from the seed chain tensor (or a previously produced SEW/FEC/LEC chain tensor already in output/)
            and applies them. Output goes to output/&lt;symmetry&gt;/ — collinear yields a basis
            (&lt;T&gt;_basis.wxf), cyclic/flip/parity yield the projected tensor.
          </p>
        </>
      )}
      {node.type === 'symmetry_invariant' && (
        <>
          <label>Symmetry</label>
          <select value={d.symmetry || 'cyclic'} onChange={(e) => set({ symmetry: e.target.value })}>
            {['cyclic', 'flip', 'parity'].map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
          <label>Target</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="SEW_5p1" />
          <p className="muted" style={{ fontSize: 11 }}>
            Runs <code>bootstrap --solve-symmetry</code>: solves for the symmetry-invariant solution basis of the
            seed chain tensor and writes output/&lt;symmetry&gt;/&lt;T&gt;_invariant.wxf. Not valid for collinear
            (use the Projection Chain node or Solve Collinear instead).
          </p>
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
          <label>Weights (comma-separated, one per input)</label>
          <input value={d.weights || ''} onChange={(e) => set({ weights: e.target.value })} placeholder="1, -1, 1/2" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. SEW_3p1_total" />
          <p className="muted" style={{ fontSize: 11 }}>Computes Σ wᵢ·Aᵢ with exact rational arithmetic (tensor_add). All inputs must have identical dimensions (same kind and weight). Weights default to 1 when omitted; extra entries are ignored.</p>
        </>
      )}
      {(node.type === 'ternary_contract' || node.type === 'apply_symmetry') && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. FEC_2_sym" />
          {node.type === 'apply_symmetry' && (
            <>
              <label>Apply n times (Sⁿ, n ≥ 1)</label>
              <input value={d.n ?? ''} onChange={(e) => set({ n: e.target.value.replace(/[^\d]/g, '') })} placeholder="1" />
              <p className="muted" style={{ fontSize: 11 }}>
                Auto mode: wire only the <b>sym</b> input with the letter-representation matrix S (e.g. cycrepmat).
                The flow walks the chain upstream (Extend/Sew → weight-1 seed), runs the recursive
                projection pipeline (weight-1 special case + per-weight induced maps, degeneracy-safe
                kernel extraction) to derive the induced basis transformations for every chain weight,
                caches them as shared matrices under <code>output/.derived/</code> (reused across blocks/flows), and
                applies the correct pair to the tensor&apos;s two trailing axes. When a trailing axis lives in a short
                seed basis (dim &lt; alphabet), the flow first composes the seed&apos;s proj matrix with S
                (E·S) to project that entry to the uniform full-alphabet dimension. With n &gt; 1 every
                matrix (S and each induced map) is first raised to the n-th power (Sⁿ) — the induced maps
                form a representation, so R(S)ⁿ = R(Sⁿ) — and E·Sⁿ is used for short seed axes.
                For FEC the seed axis comes second and
                letters last; for LEC they are swapped; a SEW applies R&#7432; on axis 2 and R&#7460; on axis 3.
                Wiring trans1/trans2 instead gives the manual mode (plain ternary contraction).
                Provenance is inherited: sums of tensors sharing the same trailing-axis meaning
                (e.g. e12 + e13 built the same way) keep it, so the same matrices apply to the sum.
              </p>
            </>
          )}
          <p className="muted" style={{ fontSize: 11 }}>TernaryContract: contracts trans1 with the 2nd-to-last axis and trans2 with the last axis of the rank-3 input tensor (T&apos;[a,b&apos;,c&apos;] = Σ T·M1·M2), exact rational arithmetic via tensor_ops. The transformations can be any matrices — symmetry transformations included.</p>
        </>
      )}
      {node.type === 'matrix_power' && (
        <>
          <label>Power n (non-negative integer)</label>
          <input value={d.n || ''} onChange={(e) => set({ n: e.target.value.replace(/[^\d-]/g, '') })} placeholder="2" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. cyc2" />
          <p className="muted" style={{ fontSize: 11 }}>M^n with exact rational arithmetic (binary exponentiation). n=0 gives the identity. Use it to generate group elements M, M², … feeding Ternary Contract.</p>
        </>
      )}
      {node.type === 'tensor_join' && (
        <>
          <label>Axis (1-based; negative counts from the end)</label>
          <input value={d.axis || ''} onChange={(e) => set({ axis: e.target.value.replace(/[^\d-]/g, '') })} placeholder="1 = first entry, -1 = last entry" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. FEC_1_joined" />
          <p className="muted" style={{ fontSize: 11 }}>Like Mathematica Join: all dimensions except the join axis must match; the join axis dimension grows.</p>
        </>
      )}
      {node.type === 'tensor_dot' && (
        <>
          <label>Axis of A (1-based; negative counts from the end)</label>
          <input value={d.axis_a ?? ''} onChange={(e) => set({ axis_a: e.target.value.replace(/[^\d-]/g, '') })} placeholder="-1 = last entry" />
          <label>Axis of B (1-based; negative counts from the end)</label>
          <input value={d.axis_b ?? ''} onChange={(e) => set({ axis_b: e.target.value.replace(/[^\d-]/g, '') })} placeholder="-1 = last entry" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. contracted" />
          <p className="muted" style={{ fontSize: 11 }}>Contraction Σ A[...,i,...]·B[...,i,...]: the two chosen axes must have equal dimension. Result = A's remaining axes followed by B's remaining axes. Dotting the last axis of A with the first of a matrix is ordinary matrix–tensor multiplication.</p>
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
      {node.type === 'integrability_condition' && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. w2f_icond" />
          <p className="muted" style={{ fontSize: 11 }}>
            Step 1 of Solve Integrability: contracts the last two axes of the rank-3 tensor S[s,i,a] with
            the integrability dlog D[a,i,c] and writes the condition matrix M[s,c] = Σ S[s,i,a]·D[a,i,c].
            Feed it to Solve Conditions (or manipulate it with matrix operations first).
          </p>
        </>
      )}
      {node.type === 'solve_conditions' && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. w2f_sol" />
          <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, marginTop: 8, color: 'var(--text)' }}>
            <input
              type="checkbox" style={{ width: 'auto' }}
              checked={d.transpose !== false}
              onChange={(e) => set({ transpose: e.target.checked })}
            />
            Transpose before solving (solve for tensor coefficients)
          </label>
          <p className="muted" style={{ fontSize: 11 }}>
            Step 2 of Solve Integrability: computes the kernel of the condition matrix M. With transpose on,
            the output basis is the kernel of Mᵀ — combinations of the s tensors satisfying all c conditions.
            With transpose off, it is the kernel of M — linear relations among the conditions.
          </p>
        </>
      )}
      {node.type === 'assemble' && (() => {
        const outs = d.outputs || []
        const nElems = edges.filter((e) => e.target === node.id && (e.targetHandle || '').startsWith('in_')).length
        const setOuts = (next) => set({ outputs: next })
        return (
          <>
            <label>Target name (output files &lt;target&gt;_&lt;output&gt;.wxf)</label>
            <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. sol_A" />
            <label>
              Outputs — each is a projection of the combined frame ({nElems || '…'} element{nElems === 1 ? '' : 's'} connected;
              {' '}one rational coefficient per element, 0 keeps zero rows)
            </label>
            {outs.map((o, oi) => (
              <div key={oi} style={{ display: 'flex', gap: 4, alignItems: 'center', marginBottom: 4 }}>
                <input
                  style={{ width: 90 }} value={o.name || ''}
                  onChange={(e) => setOuts(outs.map((x, i) => (i === oi ? { ...x, name: e.target.value } : x)))}
                  placeholder={`out ${oi + 1}`}
                />
                <input
                  value={o.coefs || ''}
                  onChange={(e) => setOuts(outs.map((x, i) => (i === oi ? { ...x, coefs: e.target.value } : x)))}
                  placeholder={nElems ? Array(nElems).fill('1').join(',') : 'e.g. 1,0,-1/2'}
                />
                <button
                  className="btn ghost" style={{ padding: '2px 7px' }}
                  onClick={() => setOuts(outs.filter((_, i) => i !== oi))}
                >✕</button>
              </div>
            ))}
            <button className="btn ghost" style={{ width: '100%' }}
              onClick={() => setOuts([...outs, { name: '', coefs: nElems ? Array(nElems).fill('1').join(',') : '' }])}
            >+ add output</button>
            <p className="muted" style={{ fontSize: 11 }}>
              Elements first combine into one aligned frame (stacked along the first axis, in elem order).
              Each output then takes its own linear combination: coefficient 0 leaves that element&apos;s block
              as zero rows, so every output has the same shape as the frame. Runs one tensor_ops assemble per output.
            </p>
          </>
        )
      })()}
      {node.type === 'merge_conditions' && <p className="muted">Connect two or more dlogmat outputs (integrability, extended Steinmann, cluster adjacency) to merge them into a single condition tensor.</p>}
      {node.type === 'sew' && <p className="muted">Combines a condition tensor, an FEC tensor of weight F and an LEC tensor of weight L into SEW_FpL (or a custom target name). Multiple Sew nodes with the same auto-name get _2, _3… suffixes.</p>}
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
          {(result.flow_outputs || []).length > 0 && (
            <div style={{ marginBottom: 8 }}>
              <div className="muted" style={{ fontSize: 10, textTransform: 'uppercase', letterSpacing: 0.5 }}>Flow outputs</div>
              {result.flow_outputs.map((o, i) => (
                <div key={i} className="step-cmd" style={{ display: 'flex', gap: 6, alignItems: 'baseline' }}>
                  <span className={`kind-badge ${o.kind || 'tensor'}`}>{o.kind || '?'}</span>
                  <strong>{o.name}</strong>
                  <span className="muted">{o.dims ? `(${o.dims.join('×')})` : ''} {o.file}</span>
                </div>
              ))}
            </div>
          )}
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
  if (node.type === 'assemble') return /^out_\d+$/.test(handleId)
  if (node.type === 'customblock') {
    const def = NODE_DEFS[`cb_${node.data?.block}`]
    if (!def) return true
    return !!(def?.outputs || []).some((o) => o.id === handleId)
  }
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
  const [autoSave, setAutoSave] = useState(() => localStorage.getItem('symbology.autosave') !== '0')
  const [groups, setGroups] = useState([])
  const [selCount, setSelCount] = useState(0)
  const [outCatalog, setOutCatalog] = useState(null)
  const [calcResultKey, setCalcResultKey] = useState('')
  const calcItems = (outCatalog || []).filter((g) => g.flow_id !== fid).flatMap((g) =>
    g.outputs.map((o) => ({ o, key: g.flow_id + o.file, label: `${g.flow_name} · ${o.name}` }))
  )
  const calcSel = calcItems.find((it) => it.key === calcResultKey) || (calcItems.length ? calcItems[0] : null)
  const rf = useRef(null)
  const wrapper = useRef(null)
  const skipAutosave = useRef(true)
  const loadedSig = useRef('')
  const stashTimer = useRef(null)

  const flow = project?.flows.find((f) => f.id === fid)
  registerCustomBlocks(project)

  useEffect(() => {
    let alive = true
    if (project?.id) {
      api.flowOutputs(project.id)
        .then((c) => { if (alive) setOutCatalog(Array.isArray(c) ? c : []) })
        .catch(() => { if (alive) setOutCatalog([]) })
    }
    return () => { alive = false }
  }, [project?.id, project?.flows])

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
      // migrate legacy assemble nodes: groups/coefs -> outputs list, out -> out_0
      const legacyEdges = flow.graph?.edges || []
      const migNodes = new Set()
      ns = ns.map((n) => {
        if (n.type !== 'assemble' || (n.data?.outputs?.length)) return n
        migNodes.add(n.id)
        const g = (n.data?.groups || '').split(',').map((s) => s.trim()).filter(Boolean)
        const c = (n.data?.coefs || '').split(',').map((s) => s.trim()).filter(Boolean)
        let coefs = null
        if (g.length && c.length) {
          try {
            const ex = g.map((gi) => c[parseInt(gi, 10) - 1])
            if (ex.every(Boolean)) coefs = ex.join(',')
          } catch { coefs = null }
        }
        const outputs = coefs
          ? [{ name: 'out1', coefs }]
          : [{ name: 'out1', coefs: '' }]
        const { groups, coefs: _c, ...rest } = n.data || {}
        return { ...n, data: { ...rest, outputs } }
      })
      const migEdges = migNodes.size
        ? legacyEdges.map((e) => (e.sourceHandle === 'out' && migNodes.has(e.source) ? { ...e, sourceHandle: 'out_0' } : e))
        : legacyEdges
      // migrate legacy add_tensors edges: a -> in_0, b -> in_1 (same shape, so ports stay aligned)
      const addMig = migEdges.map((e) => {
        if (e.targetHandle === 'a' || e.targetHandle === 'b') {
          const t = ns.find((n) => n.id === e.target)
          if (t?.type === 'add_tensors') return { ...e, targetHandle: e.targetHandle === 'a' ? 'in_0' : 'in_1' }
        }
        return e
      })
      setGroups(gr)
      setNodes(ns)
      setEdges(addMig)
      bumpNodeSeqFromIds([...ns.map((n) => n.id), ...(flow.graph?.edges || []).map((e) => e.id)])
      setCompileResult(null)
      loadedSig.current = JSON.stringify({
        name: flow.name || '',
        nodes: ns.filter((n) => n.type !== 'groupBox').map((n) => ({ id: n.id, type: n.type, position: n.position, data: n.data })),
        edges: addMig.map((e) => ({ id: e.id, source: e.source, sourceHandle: e.sourceHandle, target: e.target, targetHandle: e.targetHandle })),
        groups: gr.filter((g) => g.node_ids.some((id) => ns.some((n) => n.id === id && n.type !== 'groupBox'))),
      })
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
    if (JSON.stringify({ name: flowName, ...serializeGraph() }) === loadedSig.current) return undefined
    if (!autoSave) { setSaveState('unsaved'); return undefined }
    setSaveState('unsaved')
    const t = setTimeout(async () => {
      setSaveState('saving')
      try {
        const graph = serializeGraph()
        await api.updateFlow(project.id, fid, { name: flowName, graph })
        loadedSig.current = JSON.stringify({ name: flowName, ...graph })
        setSaveState('saved')
      } catch {
        setSaveState('error')
      }
    }, 800)
    return () => clearTimeout(t)
  }, [nodes, edges, flowName, groups, autoSave, flow]) // eslint-disable-line

  const resolveSourceKind = useCallback((nodeId, handleId, depth = 0) => {
    const node = nodes.find((n) => n.id === nodeId)
    if (!node) return null
    if ((handleId === 'out' || /^out_\d+$/.test(handleId)) && (node.type === 'add_tensors' || node.type === 'tensor_join' || node.type === 'ternary_contract' || node.type === 'apply_symmetry' || node.type === 'assemble')) {
      if (depth > 8) return 'tensor'
      const e = (node.type === 'assemble' || node.type === 'add_tensors')
        ? edges.find((ed) => ed.target === nodeId && (ed.targetHandle || '').startsWith('in_'))
        : edges.find((ed) => {
            const h = node.type === 'ternary_contract' || node.type === 'apply_symmetry' ? 'tensor' : 'a'
            return ed.target === nodeId && (ed.targetHandle === h || (h === 'a' && ed.targetHandle === 'b'))
          })
      return e ? resolveSourceKind(e.source, e.sourceHandle, depth + 1) : 'tensor'
    }
    return sourceKindFor(node, handleId, project)
  }, [nodes, edges, project])

  const compileShapes = compileResult?.ok ? compileResult.shapes : null

  const inputDimsFor = useCallback((nodeId, handleId, depth = 0) => {
    const node = nodes.find((n) => n.id === nodeId)
    if (!node || depth > 8) return null
    if (node.type === 'alphabet') return alphabetHandleDims(project, node, handleId)
    const inbound = edges.find((ed) => ed.target === nodeId && ed.targetHandle === handleId)
    if (inbound) {
      const key = `${inbound.source}|${inbound.sourceHandle}`
      const shape = compileShapes?.[key]
      if (shape) return shape.join('×')
      return inputDimsFor(inbound.source, inbound.sourceHandle, depth + 1)
    }
    return null
  }, [nodes, edges, project, compileShapes])

  const outputDimsFor = useCallback((nodeId, handleId, depth = 0) => {
    const node = nodes.find((n) => n.id === nodeId)
    if (!node || depth > 8) return null
    if (node.type === 'alphabet') return alphabetHandleDims(project, node, handleId)
    const shape = compileShapes?.[`${nodeId}|${handleId}`]
    if (shape) return shape.join('×')
    const inDims = (hid) => {
      const e = edges.find((ed) => ed.target === nodeId && ed.targetHandle === hid)
      return e ? inputDimsFor(e.source, e.sourceHandle, depth + 1) : null
    }
    const inDimsList = (hid) => {
      const d = inDims(hid)
      if (!d) return null
      return d.split('×').map((x) => parseInt(x, 10))
    }
    switch (node.type) {
      case 'add_tensors': {
        const d = inDims('in_0') || inDims('in_1') || inDims('a') || inDims('b')
        if (d) return d
        break
      }
      case 'matrix_power': {
        const d = inDims('matrix')
        if (d) return d
        break
      }
      case 'tensor_join': {
        const a = inDimsList('a')
        const b = inDimsList('b')
        const axis = Number(node.data?.axis)
        if (a && b && Number.isInteger(axis) && axis !== 0) {
          const k = axis > 0 ? axis - 1 : a.length + axis
          if (k >= 0 && k < a.length && a.length === b.length && a.every((x, i) => i === k || x === b[i])) {
            const out = [...a]
            out[k] += b[k]
            return out.join('×')
          }
        }
        break
      }
      case 'tensor_dot': {
        const a = inDimsList('a')
        const b = inDimsList('b')
        const aa = Number(node.data?.axis_a ?? -1)
        const ab = Number(node.data?.axis_b ?? -1)
        if (a && b && Number.isInteger(aa) && Number.isInteger(ab)) {
          const ka = aa > 0 ? aa - 1 : a.length + aa
          const kb = ab > 0 ? ab - 1 : b.length + ab
          if (ka >= 0 && ka < a.length && kb >= 0 && kb < b.length && a[ka] === b[kb])
            return [...a.filter((_, i) => i !== ka), ...b.filter((_, i) => i !== kb)].join('×')
        }
        break
      }
      case 'integrability_condition': {
        const t = inDimsList('tensor')
        const d = inDimsList('dlog')
        if (t && t.length >= 2 && d && d.length === 3) {
          let inner = 1
          for (let i = 1; i < t.length - 2; i++) inner *= t[i]
          return [t.length > 2 ? t[0] : 1, inner * d[2]].join('×')
        }
        break
      }
      case 'ternary_contract': {
        const t = inDimsList('tensor')
        if (t && t.length === 3) return [t[0], t[2]].join('×')
        break
      }
      case 'assemble': {
        const first = edges.find((ed) => ed.target === nodeId && (ed.targetHandle || '').startsWith('in_'))
        if (first) return inputDimsFor(first.source, first.sourceHandle, depth + 1)
        break
      }
      default:
        break
    }
    return null
  }, [nodes, edges, project, compileShapes])

  const dimsFor = useCallback((nodeId, handleId) => {
    const node = nodes.find((n) => n.id === nodeId)
    if (!node) return null
    const inbound = edges.find((ed) => ed.target === nodeId && ed.targetHandle === handleId)
    if (inbound) return inputDimsFor(inbound.source, inbound.sourceHandle)
    if (isOutputHandle(node, handleId)) return outputDimsFor(nodeId, handleId)
    return null
  }, [nodes, edges, project, inputDimsFor, outputDimsFor])

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
    const raw = e.dataTransfer.getData('application/symbology-node')
    if (!raw || !rf.current) return
    const pos = rf.current.screenToFlowPosition({ x: e.clientX, y: e.clientY })
    let type = raw
    let defaults = type === 'alphabet' ? { alphabet_id: '', selected_properties: [] } : {}
    if (raw.startsWith('customblock:')) {
      type = 'customblock'
      defaults = { block: raw.slice('customblock:'.length) }
    }
    if (raw.startsWith('reuse_output:')) {
      const [file, kind, name] = raw.slice('reuse_output:'.length).split('|')
      type = 'reuse_output'
      defaults = { file, kind, name }
    }
    if (type === 'cb_in' || type === 'cb_out') defaults = { name: '', kind: 'tensor', order: Date.now() % 100000 }
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
    try {
      const graph = serializeGraph()
      await api.updateFlow(project.id, fid, { name: flowName, graph })
      loadedSig.current = JSON.stringify({ name: flowName, ...graph })
      await refreshProject()
      setSaveState('saved')
      toast('Flow saved.')
    } catch (e) {
      setSaveState('error')
      throw e
    }
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
    if ('selected_properties' in patch && !('alphabet_id' in patch)) {
      // Defer pruning so the uncheck-A / check-B sequence can migrate the
      // dangling edges to the replacement property within a short window.
      if (stashTimer.current) clearTimeout(stashTimer.current)
      const node = nodes.find((n) => n.id === id)
      if (node?.type === 'alphabet') {
        const prevSel = node.data?.selected_properties || []
        const sel = patch.selected_properties || []
        setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, ...patch } } : n)))
        const keep = new Set()
        for (const pid of sel) { keep.add(`prop_${pid}`); keep.add(`proj_${pid}`) }
        // migrate dangling edges to newly added same-type properties
        const alphabet = project?.alphabets?.find((a) => a.id === node.data?.alphabet_id)
        const added = sel.filter((p) => !prevSel.includes(p))
        const used = new Set()
        const typeOf = (pid) => alphabet?.properties.find((p) => p.id === pid)?.type
        const migrate = new Map()
        for (const e of edges) {
          if (e.source !== id || keep.has(e.sourceHandle)) continue
          const m = /^(prop|proj)_(.+)$/.exec(e.sourceHandle || '')
          if (!m) continue
          const cand = added.find((a) => typeOf(a) === typeOf(m[2]) && !used.has(a))
          if (cand) { used.add(cand); migrate.set(e.id, `${m[1]}_${cand}`) }
        }
        if (migrate.size) setEdges((es) => es.map((e) => (migrate.has(e.id) ? { ...e, sourceHandle: migrate.get(e.id) } : e)))
        // prune anything still dangling after a grace period
        if (stashTimer.current) clearTimeout(stashTimer.current)
        stashTimer.current = setTimeout(() => {
          setEdges((es) => es.filter((e) => e.source !== id || !/^(prop|proj)_/.test(e.sourceHandle || '') || keep.has(e.sourceHandle)))
        }, 2500)
        setCompileResult(null)
        return
      }
    }
    setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, ...patch } } : n)))
    if ('selected_properties' in patch || 'alphabet_id' in patch) {
      const keep = new Set()
      for (const pid of patch.selected_properties || []) { keep.add(`prop_${pid}`); keep.add(`proj_${pid}`) }
      setEdges((es) => es.filter((e) => e.source !== id || !/^(prop|proj)_/.test(e.sourceHandle || '') || keep.has(e.sourceHandle)))
    }
    setCompileResult(null)
  }, [nodes, edges, project, setNodes, setEdges])

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
        {flow.custom_block && (
          <div>
            <div className="muted" style={{ fontSize: 10, margin: '10px 0 4px', textTransform: 'uppercase', letterSpacing: 0.5 }}>Block ports</div>
            {[
              { type: 'cb_in', label: 'Block Input', sub: 'abstract input port' },
            ].map((p) => (
              <div
                key={p.type} className="palette-item" draggable
                onDragStart={(e) => e.dataTransfer.setData('application/symbology-node', p.type)}
              >
                {p.label}
                <div className="sub">{p.sub}</div>
              </div>
            ))}
          </div>
        )}
        <div>
          <div className="muted" style={{ fontSize: 10, margin: '10px 0 4px', textTransform: 'uppercase', letterSpacing: 0.5 }}>Flow outputs</div>
          <div
            className="palette-item" draggable
            onDragStart={(e) => e.dataTransfer.setData('application/symbology-node', 'cb_out')}
          >
            Flow Output
            <div className="sub">named result of this flow</div>
          </div>
        </div>
        {Array.isArray(outCatalog) && outCatalog.filter((g) => g.flow_id !== fid && Array.isArray(g.outputs)).length > 0 && (
          <div>
            <div className="muted" style={{ fontSize: 10, margin: '10px 0 4px', textTransform: 'uppercase', letterSpacing: 0.5 }}>Calculated results</div>
            <div
              className="palette-item" draggable={Boolean(calcSel)}
              title={calcSel ? `Drag to reuse ${calcSel.label}` : 'Select an output first'}
              onDragStart={(e) => calcSel && e.dataTransfer.setData('application/symbology-node', `reuse_output:${calcSel.o.file}|${calcSel.o.kind}|${calcSel.o.name}`)}
            >
              Calculated Result
              <div className="sub">{calcSel ? calcSel.label : 'pick an output below'}</div>
            </div>
            <select
              value={calcResultKey}
              onChange={(e) => setCalcResultKey(e.target.value)}
              style={{ width: '100%', marginTop: 4, fontSize: 11 }}
            >
              {(outCatalog || []).filter((g) => g.flow_id !== fid).flatMap((g) =>
                g.outputs.map((o) => ({ o, key: g.flow_id + o.file, label: `${g.flow_name} · ${o.name}` }))
              ).map((it) => (
                <option key={it.key} value={it.key}>{it.label}</option>
              ))}
            </select>
          </div>
        )}
        {customBlockPalette(project).length > 0 && (
          <div>
            <div className="muted" style={{ fontSize: 10, margin: '10px 0 4px', textTransform: 'uppercase', letterSpacing: 0.5 }}>Custom blocks</div>
            {customBlockPalette(project).map((p) => (
              <div
                key={p.block} className="palette-item" draggable
                onDragStart={(e) => e.dataTransfer.setData('application/symbology-node', `customblock:${p.block}`)}
              >
                {p.label}
                <div className="sub">{p.sub}</div>
              </div>
            ))}
          </div>
        )}
        <div className="muted" style={{ fontSize: 11, marginTop: 12 }}>Drag a node onto the canvas. Wire colored ports of the same kind together.</div>
      </div>
      <div className="flow-canvas" ref={wrapper} onDrop={onDrop} onDragOver={(e) => e.preventDefault()}>
        <div className="flow-toolbar">
          <input style={{ width: 180 }} value={flowName} onChange={(e) => setFlowName(e.target.value)} />
          {flow.custom_block && (
            <span className="muted shrink" style={{ fontSize: 11, alignSelf: 'center' }} title="This flow is sealed as a reusable custom block (see Flows to unseal)">
              ⬢ custom block
            </span>
          )}
          <span className="muted shrink" style={{ fontSize: 11, alignSelf: 'center', minWidth: 90 }}>
            {saveState === 'saved' && '✓ saved'}
            {saveState === 'saving' && 'saving…'}
            {saveState === 'unsaved' && 'unsaved changes'}
            {saveState === 'error' && <span className="error-text">save failed</span>}
          </span>
          <label className="muted shrink" style={{ fontSize: 11, display: 'inline-flex', alignItems: 'center', gap: 4, alignSelf: 'center', cursor: 'pointer' }}>
            <input
              type="checkbox" style={{ width: 'auto' }}
              checked={autoSave}
              onChange={(e) => { setAutoSave(e.target.checked); localStorage.setItem('symbology.autosave', e.target.checked ? '1' : '0') }}
            />
            auto save
          </label>
          <button onClick={save}>Save</button>
          <button onClick={compile}>Compile</button>
          <button className="primary" onClick={compile}>Run…</button>
        </div>
        <DimsContext.Provider value={dimsFor}>
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
        </DimsContext.Provider>
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
            onOpenBlock={(bid) => navigate(`/flows/${bid}`)}
            flowId={fid}
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
