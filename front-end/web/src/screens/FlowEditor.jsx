import { openFlowSave } from '../flowSave'
import React, { useCallback, useContext, useEffect, useMemo, useRef, useState } from 'react'
import { useNavigate, useParams } from 'react-router-dom'
import ReactFlow, {
  Background, BaseEdge, ConnectionMode, Controls, Handle, MarkerType, Position, ReactFlowProvider,
  addEdge, getBezierPath, useEdges, useEdgesState, useNodesState, useUpdateNodeInternals,
} from 'reactflow'
import 'reactflow/dist/style.css'
import { api } from '../api'
import { useProject, useToast } from '../App'
import {
  NODE_DEFS, PALETTE_SECTIONS, PROP_KIND, propLabel,
  kindsCompatible, sourceKindFor, targetKindFor, alphabetOutSlots,
  registerCustomBlocks, customBlockPalette,
  normalizeCollinearPairs, PAIR_PRESETS, projectionChainBasisOutputs,
  portsOf, portBaseY, estHeight, SUBTITLE_TYPES, NODE_W,
  ROW0, ROW_STEP, SUB_EXTRA, RHS_MODES,
} from '../flowdefs'

// Wire colors follow the port-kind palette (styles.css --dlogmat/--fec/…)
// so a wire is visually the continuation of the port label it comes from.
// Output ports render their dot in the same color (see OpNode/…): a wire
// and the dot it leaves/enters never disagree.
const KIND_COLORS = {
  dlogmat: '#7c3aed', fec: '#0d9488', lec: '#ea580c', fec1: '#0d9488', lec1: '#ea580c',
  matrix: '#92400e', tensor: '#475569', basis: '#4f46e5', seed: '#0d9488',
  xtrans: '#92400e',
  matrix_or_tensor: '#92400e',
  seed_or_tensor: '#475569', solution: '#059669', boundary: '#d97706', any: '#64748b',
  sew: '#7c3aed',
}

function kindColor(kind) {
  return KIND_COLORS[kind] || '#64748b'
}

// Output-handle dot tint (input dots stay default gray: they accept several
// kinds, so no single color would be truthful).
function outDotStyle(kind) {
  return { background: kindColor(kind) }
}

// Direct port-to-port wiring ("follow the line"): a wire is ONE straight
// segment from the output port to the input port wherever possible — the
// simplest, most traceable shape. Wires that would run on top of a neighbor
// get a gentle quadratic bow (kappa, assigned in routeEdges) that keeps the
// ports pinned and the take-off/arrival directions unchanged; wires whose
// direct path would cut through node boxes route around the diagram through
// an external lane, re-entering their port horizontally. All knobs are
// derived at render time; project.json stores positions/semantics only.
function KindEdge({ id, sourceX, sourceY, targetX, targetY,
  sourcePosition, targetPosition, style = {}, markerEnd, data: edata }) {
  if (!edata || edata.fallback) {
    const fy = sourceY + ((edata && edata.spread) || 0)
    const [d] = getBezierPath({ sourceX, sourceY: fy, sourcePosition, targetX, targetY, targetPosition })
    return <BaseEdge id={id} path={d} style={style} markerEnd={markerEnd} />
  }
  const y0 = sourceY + (edata.spread || 0)
  if (typeof edata.detourY === 'number') {
    const ox = typeof edata.outboundX === 'number' ? edata.outboundX : sourceX + 30
    const inX = typeof edata.inboundX === 'number' ? edata.inboundX : targetX - 30
    const ly = edata.detourY
    // Port tangents stay horizontal (arrowheads land like on a straight
    // wire); the vertical transition happens along the lane leg. Every join
    // is G1 by construction: the rise cubic ENDS with control point
    // (ox-k0, ly) — tangent (k0,0) = the lane's run direction — and the
    // descent cubic STARTS with (inX+k1, ly). Generous corner radii
    // (k = 10..24 by leg length) keep the elbows round, and the descent
    // finishes its turn `lead` px BEFORE the port: one straight horizontal
    // glide carries the arrowhead in.
    const k0 = Math.max(10, Math.min(24, 0.6 * Math.abs(ox - sourceX)))
    const k1 = Math.max(10, Math.min(24, 0.6 * Math.abs(targetX - inX)))
    const lead = Math.max(8, Math.min(16, 0.4 * (targetX - inX)))
    const d = [
      `M ${sourceX} ${y0}`,
      `C ${sourceX + k0} ${y0} ${ox - k0} ${ly} ${ox} ${ly}`,
      `L ${inX} ${ly}`,
      `C ${inX + k1} ${ly} ${targetX - k1 - lead} ${targetY} ${targetX - lead} ${targetY}`,
      `L ${targetX} ${targetY}`,
    ].join(' ')
    return <BaseEdge id={id} path={d} style={style} markerEnd={markerEnd} />
  }
  if (edata.kappa) {
    // Bow perpendicular to the chord so it separates neighbors whatever
    // the wire's direction (vertical hops included); the apex rises
    // 0.75*b off the chord (b = 4 + 1.75|kappa|), capped by a quarter of
    // the chord LENGTH so short hops stay subtle.
    const rx = targetX - sourceX
    const ry = targetY - y0
    const len = Math.hypot(rx, ry) || 1
    const b = (4 + 1.75 * Math.abs(edata.kappa)) * Math.sign(edata.kappa)
    const bb = Math.max(-len / 4, Math.min(len / 4, b))
    const cx = (sourceX + targetX) / 2 + 1.5 * bb * (ry / len)
    const cy = (y0 + targetY) / 2 + 1.5 * bb * (-rx / len)
    const d = `M ${sourceX} ${y0} Q ${cx} ${cy} ${targetX} ${targetY}`
    return <BaseEdge id={id} path={d} style={style} markerEnd={markerEnd} />
  }
  return <BaseEdge id={id} path={`M ${sourceX} ${y0} L ${targetX} ${targetY}`} style={style} markerEnd={markerEnd} />
}

const edgeTypes = { kind: KindEdge }

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
              style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -6, ...outDotStyle(r.kind) }}
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
  const dynamicInputs = type === 'add_tensors' || type === 'expand_tensor'
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
    const labels = type === 'expand_tensor'
      ? (i) => (i === 0 ? 'T (FEC,letter)' : `basis ${i}`)
      : (i) => String.fromCharCode(65 + i)
    inputs = Array.from({ length: Math.max(inputs.length, maxConnected + 2) }, (_, i) => ({
      id: `in_${i}`, kind: 'tensor', label: labels(i),
    }))
  }
  let pairCount = 0
  let multiPair = false
  let hasCondEdge = false
  if (dynamicPairs) {
    // Pair 1 = the fixed seed/rhs ports; every entry of data.pairs beyond
    // the first adds one in_seed_N / in_rhs_N port pair (N = pair index).
    // maxPair from wired edges keeps ports alive if a saved graph wires a
    // pair that the (possibly legacy) config no longer lists.
    let maxPair = 0
    for (const e of allEdges) {
      if (e.target !== id) continue
      const m = (e.targetHandle || '').match(/^in_(?:seed|rhs)_(\d+)$/)
      if (m) maxPair = Math.max(maxPair, parseInt(m[1], 10))
    }
    hasCondEdge = allEdges.some((e) => e.target === id && (e.targetHandle || '') === 'cond')
    const pairs = normalizeCollinearPairs(data)
    const cfgLen = pairs.length - 1
    pairCount = Math.max(maxPair + 1, cfgLen + 1)
    multiPair = maxPair > 0 || cfgLen > 0 || hasCondEdge
    const nPairs = Math.max(maxPair, cfgLen)
    const extras = []
    for (let p = 1; p <= nPairs; p++) {
      extras.push({ id: `in_seed_${p}`, kind: 'seed_or_tensor', label: `seed ${p + 1}` })
      extras.push({ id: `in_rhs_${p}`, kind: 'boundary', label: `rhs ${p + 1}` })
    }
    // cond is opt-in: shown only when enabled in the Inspector or wired.
    inputs = inputs.filter((i) => i.id !== 'cond')
    if (data.cond_enabled || hasCondEdge) inputs = [...inputs, { id: 'cond', kind: 'matrix', label: 'cond [M|r] (opt)' }]
    inputs = [...inputs, ...extras]
  }
  const dynamicChainOutputs = type === 'projection_chain'
  let outputs = def.outputs || []
  if (dynamicChainOutputs) {
    // Collinear targets: one extra output port per expansion basis the same
    // --project run materializes (contract shared with compile.py).
    const basis = projectionChainBasisOutputs(data)
    // Keep ports alive for saved graphs that wire a basis handle even when
    // data.target is empty (the server derives the target from the seed).
    for (const e of allEdges) {
      if (e.source === id && /^(basis_w\d+|basis_last_w\d+)$/.test(e.sourceHandle || '')) {
        if (!basis.some((b) => b.id === e.sourceHandle)) basis.push({ id: e.sourceHandle, kind: 'basis', label: e.sourceHandle })
      }
    }
    outputs = [...outputs, ...basis]
  }
  useEffect(() => { if (dynamicInputs || dynamicPairs || dynamicChainOutputs) updateNodeInternals(id) }, [id, dynamicInputs, dynamicPairs, dynamicChainOutputs, inputs.length, outputs.length, updateNodeInternals])
  const rows = Math.max(inputs.length, outputs.length)
  const cbWeights = (data.weights || '').split(',').map((s) => s.trim())
  const rhsMode = (RHS_MODES[data.mode] || RHS_MODES.mhv_boundary).short
  const subtitle =
    type === 'extend' ? (data.target_weight ? `→ weight ${data.target_weight}` : 'weight +1')
    : (type === 'project' || type === 'solve_symmetry' || type === 'symderive') ? `→ ${data.target || '…'}`
    : type === 'sew' ? `→ ${data.target || 'SEW_FpL'}`
    : type === 'solve_collinear' ? (multiPair
        ? `${pairCount} pair${pairCount === 1 ? '' : 's'}${hasCondEdge ? ' + cond' : ''} → ${data.out_stem || '…'}`
        : `${data.target || '…'}`)
    : type === 'projection_chain' ? `${data.symmetry || '?'} · ${data.target || '…'}`
    : type === 'symmetry_invariant' ? `${data.symmetry || '?'} · ${data.target || '…'}`
    : type === 'compute_rhs' ? `${rhsMode} w${data.weight || '?'} → ${data.target || '…'}`
    : type === 'add_tensors' ? `${inputs.map((_, i) => `${cbWeights[i] || '1'}·${String.fromCharCode(65 + i)}`).join(' + ')} → ${data.target || '…'}`
    : (type === 'ternary_contract' || type === 'apply_symmetry') ? `→ ${data.target || '…'}`
    : type === 'matrix_power' ? `M^${data.n || '?'} → ${data.target || '…'}`
    : type === 'tensor_join' ? `axis ${data.axis || '?'} → ${data.target || '…'}`
    : type === 'tensor_dot' ? `A[${data.axis_a ?? '?'}]·B[${data.axis_b ?? '?'}] → ${data.target || '…'}`
    : type === 'squeeze_tensor' ? `${handleDims(id, 'in') || '?'} → ${data.target || '…'}`
    : type === 'shuffle_product' ? `${data.weight || '1'}·A⊗B → ${data.target || '…'}`
    : type === 'expand_tensor' ? `${inputs.length - 1 || '?'} basis${(inputs.length - 1) === 1 ? '' : 'es'} → ${data.target || '…'}`
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
            style={{ top: (subtitle ? ROW0 + SUB_EXTRA : ROW0) + i * ROW_STEP, ...outDotStyle(out.kind) }}
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
          <Handle key={`out_${oi}`} type="source" position={Position.Right} id={`out_${oi}`}
            style={{ top: ROW0 + SUB_EXTRA + oi * ROW_STEP, ...outDotStyle(o.kind || 'tensor') }} />
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
              style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -18, ...outDotStyle(data?.kind || 'any') }} />
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
            style={{ top: ROW0 + i * ROW_STEP, ...outDotStyle(out.kind) }} />
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
            style={{ position: 'absolute', top: '50%', transform: 'translateY(-50%)', right: -6, ...outDotStyle(data?.kind || 'tensor') }}
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
  squeeze_tensor: (p) => <OpNode {...p} type="squeeze_tensor" />,
  shuffle_product: (p) => <OpNode {...p} type="shuffle_product" />,
  expand_tensor: (p) => <OpNode {...p} type="expand_tensor" />,
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
          <p className="muted" style={{ fontSize: 11 }}>Any weight above the input FEC/LEC weight — intermediate extends run automatically in the background (each +1 step is a cached step). Connect exactly one of FEC in / LEC in — LEC extends backward.</p>
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
            Wire a tensor into the <b>seed</b> port to solve a custom seed (used as-is): the tensor is used
            as-is (no seed-space projection, no expansion unless --basis). Leave blank for named SEW/FEC targets.
            In multi-pair mode every seed — pair 1's port included — must be a wired tensor.
          </p>
          <label>Seed projection (single-pair, named targets)</label>
          <select value={d.projection || 'finite'} onChange={(e) => set({ projection: e.target.value })}>
            <option value="finite">finite</option>
            <option value="divergent">divergent</option>
            <option value="none">none (custom seed)</option>
          </select>
          <p className="muted" style={{ fontSize: 11 }}>
            Seed-space projection for the single-pair named-target contract. Custom / multi-pair seeds are always
            used as-is; the per-pair option below is the <b>letter</b> projection.
          </p>
          <label>Pairs — pair 1 = the seed / rhs ports; each added pair grows its own seed N / rhs N ports</label>
          {(() => {
            const pairs = normalizeCollinearPairs(d)
            const modeOf = (v) => (PAIR_PRESETS.some(([p]) => p === v) ? v : 'custom')
            const setRow = (i, patch) => {
              set({ pairs: pairs.map((p, j) => (j === i ? { ...p, ...patch } : p)) })
            }
            return (
              <>
                {pairs.map((p, i) => {
                  const mode = modeOf(p.projection)
                  return (
                    <div key={i} style={{ border: '1px solid var(--border)', borderRadius: 6, padding: '4px 6px', marginBottom: 6 }}>
                      <div style={{ display: 'flex', gap: 4, alignItems: 'center', marginBottom: 4 }}>
                        <span className="muted" style={{ fontSize: 11, width: 20 }}>P{i + 1}</span>
                        <span className="muted" style={{ fontSize: 10 }}>proj</span>
                        <select value={mode} onChange={(e) => setRow(i, { projection: e.target.value === 'custom' ? '' : e.target.value })}>
                          {PAIR_PRESETS.map(([pv, label]) => <option key={pv} value={pv}>{label}</option>)}
                          <option value="custom">custom file…</option>
                        </select>
                        {mode === 'custom' && (
                          <input style={{ flex: 1 }} value={p.projection} onChange={(e) => setRow(i, { projection: e.target.value })} placeholder="letterproj.wxf" />
                        )}
                        {i > 0 && (
                          <button
                            className="btn ghost" style={{ padding: '2px 7px' }}
                            onClick={() => set({ pairs: pairs.filter((_, j) => j !== i) })}
                          >✕</button>
                        )}
                      </div>
                      <div style={{ display: 'flex', gap: 4, alignItems: 'center' }}>
                        <span className="muted" style={{ fontSize: 11, width: 20 }} />
                        <span className="muted" style={{ fontSize: 10 }}>rhs</span>
                        <input style={{ flex: 1 }} value={p.rhs} onChange={(e) => setRow(i, { rhs: e.target.value })}
                          placeholder={i === 0 ? 'data/E1.wxf or 0 (or wire rhs)' : 'rhs.wxf or 0 (or wire rhs N)'} />
                      </div>
                    </div>
                  )
                })}
                <button
                  className="btn ghost" style={{ width: '100%' }}
                  onClick={() => set({ pairs: [...pairs, { projection: 'identity', rhs: '' }] })}
                >+ add pair</button>
              </>
            )
          })()}
          <p className="muted" style={{ fontSize: 11 }}>
            Each pair is one { '{seed, rhs}' } constraint set: its own letter projection — a sentinel
            (identity / divergent = keep entries with any divergent letter / finite = keep all-finite
            entries) or a custom .wxf file path (legacy per-slot contraction, e.g. data/colprojdiv.wxf)
            — and its own RHS — a file path, “0” for homogeneous constraints,
            or a wired rhs port (the wire wins over the field). Pair 1 uses the fixed seed / rhs ports; adding a pair
            grows the seed N / rhs N ports and switches the node to multi-pair mode: all rows are stacked and solved
            together, every seed is used as-is.
          </p>
          <label style={{ display: 'inline-flex', alignItems: 'center', gap: 6, marginTop: 8, color: 'var(--text)' }}>
            <input
              type="checkbox" style={{ width: 'auto' }}
              checked={!!d.cond_enabled}
              onChange={(e) => set({ cond_enabled: e.target.checked })}
            />
            Conditions input (cond port)
          </label>
          <p className="muted" style={{ fontSize: 11 }}>
            Shows a third input port for rank-2 [M|r] conditions matrices (e.g. the conditions output of another Solve
            Collinear node); every wired cond file stacks extra rows into the same solve.
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
            (&lt;T&gt;_basis.wxf), cyclic/flip/parity yield the projected tensor. Collinear targets also expose one
            extra output per chain-weight expansion basis (first_w&#123;w&#125;_basis / last_w&#123;w&#125;_basis) written by the same
            run — wire them into Expand Tensor inputs.
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
        (() => {
          const mode = RHS_MODES[d.mode] || RHS_MODES.mhv_boundary
          const wired = (p) => edges.some((e) => e.target === node.id && e.targetHandle === p)
          const missingWires = mode.requires.filter((p) => !wired(p))
          const strayWires = mode.forbids.filter(wired)
          return (
            <>
              <label>Object type to generate</label>
              <select value={d.mode || 'mhv_boundary'} onChange={(e) => set({ mode: e.target.value })}>
                {Object.entries(RHS_MODES).map(([k, m]) => (
                  <option key={k} value={k}>{m.label}</option>
                ))}
              </select>
              <label>Target weight (even integer; loop order = weight/2)</label>
              <input value={d.weight || ''} onChange={(e) => set({ weight: e.target.value })} placeholder="4" />
              <label>Target name (output tensor)</label>
              <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder={mode.targetPlaceholder} />
              {missingWires.length > 0 && (
                <p style={{ fontSize: 11, color: 'var(--warn)' }}>
                  This mode needs {missingWires.map((p) => <code key={p}>{p}</code>)} wired — compile will reject it until then.
                </p>
              )}
              {strayWires.length > 0 && (
                <p style={{ fontSize: 11, color: 'var(--warn)' }}>
                  <code>{strayWires.join(', ')}</code> is not used by this mode — unwire it or switch mode.
                </p>
              )}
              {d.mode === 'e47me67_te' ? (
                <p className="muted" style={{ fontSize: 11 }}>
                  NMHV E47mE67: tE_L boundary = Σ tP_k ⊗ E_(L−k) with tP_1 = P1 = hep1LE47mE67
                  wired into <b>p1</b>; tP_k (k ≥ 2) loads from output/tP&lt;k&gt;.wxf or is derived
                  from output/tE&lt;k&gt;.wxf via the hardcoded recursion tP_k = tE_k − Σ tP_j ⊗ E_(k−j)
                  (missing both = compile error — never silently zero). Lower-loop E_k auto-load
                  from output/. Output is a concrete tensor (output/&lt;target&gt;.wxf). Shuffles are
                  sequential (the verified variant).
                </p>
              ) : (
                <p className="muted" style={{ fontSize: 11 }}>
                  MHV boundary: boundary_L = (1/L)·Σ k·(R_k ⊗ E_(L−k)), R_1 = E1; lower-loop
                  R_k/E_k auto-load from output/. Output is a concrete tensor
                  (output/&lt;target&gt;.wxf) wireable anywhere. Shuffles are sequential
                  (the verified variant).
                </p>
              )}
            </>
          )
        })()
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
          <p className="muted" style={{ fontSize: 11 }}>TernaryContract: contracts trans1 with the second axis and trans2 with the last axis. A trans port left unlinked defaults to the identity (that entry does not transform). Besides matrices, a trans port accepts any tensor of rank 2 or more: its first axis is contracted with the tensor axis and its remaining axes are inserted at the contracted position (e.g. expansion with a chain basis).</p>
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
      {node.type === 'squeeze_tensor' && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. hepMHV_2L" />
          <p className="muted" style={{ fontSize: 11 }}>
            Drops every size-1 axis of the input tensor (e.g. the (1, 76, 11) Tensor Dot of solMHV with
            the SEW basis becomes (76, 11), matching the ground-truth hepMHV files). The result must
            keep at least 2 axes — vectors are stored as (1, n) matrices in WXF, so squeeze cannot
            produce them.
          </p>
        </>
      )}
      {node.type === 'shuffle_product' && (
        <>
          <label>Weight (rational, applied to the product)</label>
          <input value={d.weight || ''} onChange={(e) => set({ weight: e.target.value })} placeholder="1/2 (for E1 ⊗ E1 → boundary_2L)" />
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. boundary_2L" />
          <p className="muted" style={{ fontSize: 11 }}>
            Shuffle product of two word tensors (letters of A shuffle with letters of B): the result&apos;s
            letter axes are A&apos;s axes followed by B&apos;s, every entry scaled by the rational weight.
            This is the boundary term building block: E1 ⊗ E1 with weight 1/2 gives E1²/2,
            E1 ⊗ R2 with weight 1 gives E1·R2, etc.
          </p>
        </>
      )}
      {node.type === 'expand_tensor' && (
        <>
          <label>Target name (output file)</label>
          <input value={d.target || ''} onChange={(e) => set({ target: e.target.value })} placeholder="e.g. E2" />
          <p className="muted" style={{ fontSize: 11 }}>
            Expands a compressed (FEC, letter) or (1, FEC, letter) tensor back to the full 42-letter
            alphabet: in_0 is the tensor (e.g. the Tensor Dot of solMHV with the SEW basis), in_1.. are
            the rank-3 (FEC, FEC&apos;, letter) bases from Projection Chain nodes — highest weight first
            (e.g. first_w3_basis then first_w2_basis for the 2-loop E2). Each basis input consumes the
            current FEC axis and appends one letter axis; extra inputs can be added by wiring more
            basis ports.
          </p>
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

function CompilePanel({ result, onRun, onExport, running, exporting, exportInfo }) {
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
          <button disabled={exporting} onClick={onExport} style={{ marginTop: 6, width: '100%' }}
                  title="Compile this flow into a portable, self-checking bash script you can copy to a cluster (see README: design locally, run anywhere the C++ core is built)">
            {exporting ? 'Exporting…' : '⇪ Export standalone script'}
          </button>
          <div className="muted" style={{ marginTop: 4, fontSize: 11 }}>
            Export writes this plan as a portable, self-checking bash script (preflight checks, per-step
            output verification, CRC32 echoes, resume) to <code>exported/&lt;flow&gt;.sh</code> inside the
            project. Copy the repo + project dir to another machine, build with <code>make</code>, and run
            it there — no front-end needed. See README → “Design locally, run on the cluster”.
          </div>
          {exportInfo && (
            <div className="step-cmd" style={{ marginTop: 6 }}>
              {exportInfo.error
                ? <span className="error-text">{exportInfo.error}</span>
                : <div>Written to <code>{exportInfo.path}</code> ({exportInfo.n_steps} steps,
                  {exportInfo.n_wolfram_steps > 0
                    ? ` ${exportInfo.n_wolfram_steps} Wolfram step${exportInfo.n_wolfram_steps > 1 ? 's' : ''} — needs Mathematica or pre-copied outputs`
                    : ' no Wolfram steps'}).<br />
                  Run anywhere the C++ core is built:
                  <code> SYMBOLOGY_ROOT=&lt;repo&gt; PROJ_DIR=&lt;project&gt; bash {exportInfo.path}</code>
                </div>}
            </div>
          )}
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


// Pure wire-routing pipeline over explicit inputs, extracted from the old
// component-local visibleEdges so the SSR audit can execute the REAL
// routing against every flow (component-only code was exactly how the
// alphabet def.outputs crash passed every audit and shipped).
export function routeEdges({ edges, nodes, project, groupOf, resolveSourceKind }) {
    // Direct port-to-port routing ("follow the line"). Most wires are ONE
    // straight segment between their two ports; the interesting cases get
    // small explicit knobs, all derived here at render time:
    //  - fan-out spread: wires leaving the same port take off a few px apart;
    //  - obstacle detour: a wire whose direct path would cut through node
    //    boxes routes around via a LOCAL external lane and re-enters its
    //    port horizontally (port tangents stay clean);
    //  - kappa bow: wires that would run on top of a neighbor (locally
    //    parallel, <7px apart) bow apart symmetrically by bundle rank.
    const boxes = nodes.map((n) => ({
      x0: n.position.x, x1: n.position.x + NODE_W,
      y0: n.position.y - 26, y1: n.position.y + estHeight(n, project, edges) + 26,
    }))
    // True handle y from the TRUE rendered port lists (portsOf) — dynamic
    // types (alphabet, chain basis, collinear pairs, in_N slots) included;
    // unknown handles fall back to row 0 instead of crashing.
    const portY = (n, handle, side) => {
      const { inputs, outputs } = portsOf(n, project, edges)
      const list = side === 'out' ? outputs : inputs
      const row = Math.max(0, list.findIndex((p) => p.id === handle))
      return n.position.y + portBaseY(n) + row * ROW_STEP
    }
    const nodeById = new Map(nodes.map((n) => [n.id, n]))
    const wires = []
    for (const e of edges) {
      const gs = groupOf.get(e.source)
      const gt = groupOf.get(e.target)
      if (gs || gt) continue
      const sn = nodeById.get(e.source)
      const tn = nodeById.get(e.target)
      if (!sn || !tn) continue
      wires.push({
        e, sn, tn,
        x0: sn.position.x + NODE_W, x1: tn.position.x,
        y0: portY(sn, e.sourceHandle, 'out'), y1: portY(tn, e.targetHandle, 'in'),
        spread: 0, kappa: 0, detourY: null, outboundX: 0, inboundX: 0, fallback: false,
      })
    }
    // 1) Fan-out spread per source port (rank by target y, then target x, so
    // sibling wires part ways immediately at the take-off and never re-cross).
    const portGrp = new Map()
    for (const w of wires) {
      const k = `${w.e.source}|${w.e.sourceHandle || ''}`
      if (!portGrp.has(k)) portGrp.set(k, [])
      portGrp.get(k).push(w)
    }
    for (const group of portGrp.values()) {
      if (group.length < 2) continue
      group.sort((a, b) => (a.y1 - b.y1) || (a.x1 - b.x1))
      group.forEach((w, i) => { w.spread = (i - (group.length - 1) / 2) * 7 })
    }
    for (const w of wires) w.y0 += w.spread
    // 2) Obstacle classification. The direct path is the STRAIGHT chord
    // (what renders unless bowed): a wire is blocked when the chord passes
    // through a box lying strictly between the two ports. The source and
    // target boxes are excluded BY CONSTRUCTION: the chord starts at the
    // source's right edge (x0 = box.x1) and ends at the target's left edge
    // (x1 = box.x0), and the strict test interval (x0+1, x1-1) never
    // includes either. (Never replace this with a generic segment-rect
    // intersection — that would flag every wire on its own two boxes.)
    // Exact segment-vs-rect interior test (Liang-Barsky). The old 9-point
    // sampler could step over thin corner slivers and miss real hits.
    const segRect = (ax, ay, bx, by, b) => {
      const dx = bx - ax
      const dy = by - ay
      let t0 = 0
      let t1 = 1
      const clip = (p, q) => {
        if (p === 0) return q >= 0
        const r = q / p
        if (p < 0) {
          if (r > t1) return false
          if (r > t0) t0 = r
        } else {
          if (r < t0) return false
          if (r < t1) t1 = r
        }
        return true
      }
      if (!clip(-dx, ax - b.x0)) return false
      if (!clip(dx, b.x1 - ax)) return false
      if (!clip(-dy, ay - b.y0)) return false
      if (!clip(dy, b.y1 - ay)) return false
      return t1 > t0
    }
    const detourSet = new Set()
    for (const w of wires) {
      // A wire doubling back more than 24px would under-run its own source
      // box and re-emerge from its LEFT edge — always detour. Everything
      // else gets the full chord-vs-box test: a dx>=90 shortcut would be
      // unsound because jittered/stacked columns overlap horizontally and
      // even short chords can cut a sibling box.
      if (w.x1 - w.x0 < -24) { detourSet.add(w); continue }
      for (const b of boxes) {
        if (b.x1 <= w.x0 + 1 || b.x0 >= w.x1 - 1) continue
        if (segRect(w.x0, w.y0, w.x1, w.y1, b)) { detourSet.add(w); break }
      }
    }
    let directs = wires.filter((w) => !detourSet.has(w))
    // Chord polylines + true segment-segment distance (parallel case via
    // point-segment mins) — point-pair distance underestimates closeness of
    // near-parallel neighbors between sample points.
    const samples = (x0, y0, x1, y1) => {
      const out = []
      for (let i = 0; i <= 16; i++) {
        const t = i / 16
        out.push([x0 + (x1 - x0) * t, y0 + (y1 - y0) * t])
      }
      return out
    }
    const segSeg = (ax0, ay0, ax1, ay1, bx0, by0, bx1, by1) => {
      const dx = ax1 - ax0, dy = ay1 - ay0
      const ex = bx1 - bx0, ey = by1 - by0
      const denom = dx * ey - dy * ex
      const ptSeg = (px, py, qx0, qy0, qx1, qy1) => {
        const l2 = (qx1 - qx0) ** 2 + (qy1 - qy0) ** 2
        if (!l2) return Math.hypot(px - qx0, py - qy0)
        let t = ((px - qx0) * (qx1 - qx0) + (py - qy0) * (qy1 - qy0)) / l2
        t = Math.max(0, Math.min(1, t))
        return Math.hypot(px - (qx0 + t * (qx1 - qx0)), py - (qy0 + t * (qy1 - qy0)))
      }
      if (Math.abs(denom) < 1e-9) {
        return Math.min(ptSeg(ax0, ay0, bx0, by0, bx1, by1), ptSeg(ax1, ay1, bx0, by0, bx1, by1),
          ptSeg(bx0, by0, ax0, ay0, ax1, ay1), ptSeg(bx1, by1, ax0, ay0, ax1, ay1))
      }
      const t0 = ((bx0 - ax0) * ey - (by0 - ay0) * ex) / denom
      const t1 = ((bx0 - ax0) * dy - (by0 - ay0) * dx) / denom
      if (t0 >= 0 && t0 <= 1 && t1 >= 0 && t1 <= 1) return 0
      return Math.min(ptSeg(ax0, ay0, bx0, by0, bx1, by1), ptSeg(ax1, ay1, bx0, by0, bx1, by1),
        ptSeg(bx0, by0, ax0, ay0, ax1, ay1), ptSeg(bx1, by1, ax0, ay0, ax1, ay1))
    }
    const chordAngle = (a, b) => {
      const dax = a.x1 - a.x0, day = a.y1 - a.y0
      const dbx = b.x1 - b.x0, dby = b.y1 - b.y0
      const dot = dax * dbx + day * dby
      const cos = dot / (Math.hypot(dax, day) * Math.hypot(dbx, dby) || 1)
      return { dot, ang: Math.acos(Math.max(-1, Math.min(1, cos))) }
    }
    // 2b) Narrow-X promotion. Antiparallel wires crossing at a SHALLOW
    // angle (<20°) closer than 2.5px tangle near-vertically at the port —
    // bows cannot undo a crossing. Promote the longer-chord member to a
    // detour (lane + horizontal descent) so it arrives cleanly from a lane.
    // Steep crossings are readable X shapes and stay.
    const narrowX = []
    for (let i = 0; i < directs.length; i++) {
      for (let j = i + 1; j < directs.length; j++) {
        const a = directs[i]
        const b = directs[j]
        const pa = samples(a.x0, a.y0, a.x1, a.y1)
        const pb = samples(b.x0, b.y0, b.x1, b.y1)
        let best = Infinity
        for (let m = 0; m <= 15; m++) {
          for (let n2 = 0; n2 <= 15; n2++) {
            const dd = segSeg(pa[m][0], pa[m][1], pa[m + 1][0], pa[m + 1][1],
              pb[n2][0], pb[n2][1], pb[n2 + 1][0], pb[n2 + 1][1])
            if (dd < best) best = dd
          }
        }
        if (best >= 2.5) continue
        const { dot, ang } = chordAngle(a, b)
        const acute = dot < 0 ? Math.PI - ang : ang
        if (acute >= 0.35) continue
        // Bows cannot undo a crossing, and a crossing this shallow is an
        // unreadable tangle: promote when the chords actually CROSS
        // (proper intersection) or run antiparallel overlapping.
        const cross = segSeg(a.x0, a.y0, a.x1, a.y1, b.x0, b.y0, b.x1, b.y1) === 0
        if (cross || dot < 0) narrowX.push([a, b])
      }
    }
    for (const [a, b] of narrowX) {
      const promote = Math.hypot(a.x1 - a.x0, a.y1 - a.y0) >= Math.hypot(b.x1 - b.x0, b.y1 - b.y0) ? a : b
      detourSet.add(promote)
    }
    if (narrowX.length) directs = directs.filter((w) => !detourSet.has(w))
    // 3) Local external lanes for blocked wires. A lane sits only outside
    // the boxes the horizontal run actually passes over (keeps lanes close
    // to their wires and the fitView zoom unchanged); stacked lanes are
    // 22px apart, and wires whose lane x-spans overlap never share a lane y.
    const topLanes = []
    const botLanes = []
    for (const w of detourSet) {
      // Scan margin mirrors the anchor geometry: the lane run extends up to
      // maxAnchorOffset (79) beyond each end — left of inX, and for backward
      // wires also past the target box to the outbound anchor. Forward
      // detours (blocked chord, target right of source) run BETWEEN the
      // boxes: ox (right of source) to inX (left of target) — span [x0, x1].
      const lo = Math.min(w.x0, w.x1) - 80
      const hi = w.x1 >= w.x0 ? w.x1 + 80 : Math.max(w.x0, w.x1 + NODE_W) + 80
      let topY = 0
      let botY = 0
      for (const b of boxes) {
        if (b.x1 <= lo || b.x0 >= hi) continue
        if (!topY || b.y0 < topY) topY = b.y0
        if (!botY || b.y1 > botY) botY = b.y1
      }
      const firstTop = (topY || 0) - 40
      const firstBot = (botY || 0) + 40
      const takeTop = w.y0 <= w.y1
      const overlaps = (a) => Math.min(a.hi, hi) - Math.max(a.lo, lo) > 0
      const lanes = takeTop ? topLanes : botLanes
      let y = takeTop ? firstTop : firstBot
      let stack = 1
      while (lanes.some((a) => overlaps(a) && a.y === y)) {
        y = takeTop ? firstTop - 22 * stack : firstBot + 22 * stack
        stack += 1
      }
      lanes.push({ lo, hi, y })
      w.detourY = y
    }
    // 3b) Anchor search: rise/drop legs must not cut through boxes. Try
    // offsets 30, 37, … (7px steps fan the rises apart) and take the first
    // pair whose legs clear the boxes. No pair clears → plain fallback
    // bezier (rare, e.g. two nodes drawn almost on top of each other).
    // A leg's DRAWN shape is a cubic whose control points can pull the
    // curve inside a box the straight chord clears — sample the real curve
    // (control points exactly as KindEdge builds them for each leg).
    const cubicHits = (p0, c1, c2, p1) => {
      const hx0 = Math.min(p0[0], c1[0], c2[0], p1[0])
      const hx1 = Math.max(p0[0], c1[0], c2[0], p1[0])
      for (const b of boxes) {
        if (b.x1 <= hx0 || b.x0 >= hx1) continue
        for (let i = 0; i <= 20; i++) {
          const t = i / 20
          const u = 1 - t
          const px = u * u * u * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t * t * t * p1[0]
          const py = u * u * u * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t * t * t * p1[1]
          if (px > b.x0 && px < b.x1 && py > b.y0 - 2 && py < b.y1 + 2) return true
        }
      }
      return false
    }
    const riseClear = (fromX, portYv, toX, laneY) => {
      const k = Math.max(10, Math.min(24, 0.6 * Math.abs(toX - fromX)))
      return !cubicHits([fromX, portYv], [fromX + k, portYv],
        [toX - k, laneY], [toX, laneY])
    }
    const descClear = (toX, portYv, fromX, laneY) => {
      // EXACT drawn shape: cubic from the lane down to the glide start
      // (already horizontal), then a straight glide into the port.
      const leg = toX - fromX
      const k = Math.max(10, Math.min(24, 0.6 * Math.abs(leg)))
      const lead = Math.max(8, Math.min(16, 0.4 * Math.abs(leg)))
      if (cubicHits([fromX, laneY], [fromX + k, laneY],
        [toX - k - lead, portYv], [toX - lead, portYv])) return false
      for (const b of boxes) {
        if (segRect(toX - lead, portYv, toX, portYv, b)) return false
      }
      return true
    }
    // Two detours leaving the SAME source port must not share a riser:
    // siblings rotate the anchor ladder so their first choice differs.
    const detourRank = new Map()
    for (const w of detourSet) {
      const k = `${w.e.source}|${w.e.sourceHandle || ''}`
      const r = detourRank.get(k) || 0
      detourRank.set(k, r + 1)
      w.anchorShift = r
    }
    // Anchor ladder, FAR-FIRST: prefer o≈44 so the descent elbow and its
    // horizontal lead-in clear the target block early; walk IN to 30 for
    // tight corridors, then OUT to 79 for wide obstacles. Siblings start
    // r steps into the rotation, keeping riser xs ≥7px apart.
    const ANCHOR_OFFS = [44, 37, 30, 51, 58, 65, 72, 79]
    for (const w of detourSet) {
      let ok = false
      for (let i = 0; i < 8 && !ok; i++) {
        const o = ANCHOR_OFFS[(i + w.anchorShift) % 8]
        // Outbound anchor bases: right of the source box, plus — for
        // backward wires whose target box overlaps the source column —
        // right of the target box, so the rise leg can clear it.
        const bases = [w.x0]
        if (w.x1 < w.x0 && w.x1 + NODE_W > w.x0) bases.push(w.x1 + NODE_W)
        // Inbound anchor sits LEFT of the target port in BOTH directions —
        // the arrowhead approaches from the left like every other input,
        // and the descent leg stays outside the target box's x-span. For
        // backward wires the lane run then passes over the boxes at lane-y
        // (outside all boxes), leftward from ox to inX.
        const inX = w.x1 - o
        if (Math.abs(inX - w.x0) < 24) continue
        for (const base of bases) {
          const ox = base + o
          if (riseClear(w.x0, w.y0, ox, w.detourY) && descClear(w.x1, w.y1, inX, w.detourY)) {
            w.outboundX = ox
            w.inboundX = inX
            ok = true
            break
          }
        }
      }
      if (!ok) { w.detourY = null; w.fallback = true }
    }
    // 4) Kappa bows for near-parallel overlapping direct wires (helpers
    // samples/segSeg/chordAngle hoisted above the lane pass). Shared-port
    // exclusion zones shrink the comparison window near the shared
    // take-off/arrival (fan-out spread separates them there).
    const tang = []
    for (let i = 0; i < directs.length; i++) {
      for (let j = i + 1; j < directs.length; j++) {
        const a = directs[i]
        const b = directs[j]
        const sharedSrc = a.e.source === b.e.source && (a.e.sourceHandle || '') === (b.e.sourceHandle || '')
        const sharedTgt = a.e.target === b.e.target && (a.e.targetHandle || '') === (b.e.targetHandle || '')
        if (sharedSrc && sharedTgt) continue
        const pa = samples(a.x0, a.y0, a.x1, a.y1)
        const pb = samples(b.x0, b.y0, b.x1, b.y1)
        let best = Infinity
        for (let m = 0; m <= 15; m++) {
          const jLo = sharedSrc ? 5 : 0
          const jHi = sharedTgt ? 10 : 15
          for (let n2 = jLo; n2 <= jHi; n2++) {
            const dd = segSeg(pa[m][0], pa[m][1], pa[m + 1][0], pa[m + 1][1],
              pb[n2][0], pb[n2][1], pb[n2 + 1][0], pb[n2 + 1][1])
            if (dd < best) best = dd
          }
        }
        if (best >= 7) continue
        // SIGNED angle test: antiparallel wires (dot < 0) never bundle —
        // they cross (shallow X pairs were detoured in 2b); near-parallel
        // same-direction wires are the only bow candidates.
        const { dot, ang } = chordAngle(a, b)
        if (dot < 0 || ang >= 0.35) continue
        tang.push([a, b])
      }
    }
    if (tang.length) {
      const adj = new Map()
      for (const w of directs) adj.set(w, [])
      for (const [a, b] of tang) { adj.get(a).push(b); adj.get(b).push(a) }
      const seen = new Set()
      for (const w of directs) {
        if (seen.has(w) || !adj.get(w).length) continue
        const bundle = []
        const queue = [w]
        seen.add(w)
        while (queue.length) {
          const cur = queue.shift()
          bundle.push(cur)
          for (const nb of adj.get(cur)) if (!seen.has(nb)) { seen.add(nb); queue.push(nb) }
        }
        bundle.sort((p, q) => ((p.y0 + p.y1) - (q.y0 + q.y1)) || (p.x0 - q.x0))
        // A bow must never push a wire into a node box: shrink each member's
        // effective bow until its DRAWN quad path clears all boxes (dropping
        // to kappa 0 keeps the wire straight and safe).
        const quadHitsBox = (w, kappa) => {
          const rx = w.x1 - w.x0
          const ry = w.y1 - w.y0
          const len = Math.hypot(rx, ry) || 1
          const b = (4 + 1.75 * Math.abs(kappa)) * Math.sign(kappa)
          const bb = Math.max(-len / 4, Math.min(len / 4, b))
          const cx = (w.x0 + w.x1) / 2 + 1.5 * bb * (ry / len)
          const cy = (w.y0 + w.y1) / 2 + 1.5 * bb * (-rx / len)
          for (const b2 of boxes) {
            if (b2.x1 <= Math.min(w.x0, w.x1, cx) || b2.x0 >= Math.max(w.x0, w.x1, cx)) continue
            for (let i = 0; i <= 24; i++) {
              const t = i / 24
              const u = 1 - t
              const px = u * u * w.x0 + 2 * u * t * cx + t * t * w.x1
              const py = u * u * w.y0 + 2 * u * t * cy + t * t * w.y1
              if (px > b2.x0 && px < b2.x1 && py > b2.y0 && py < b2.y1) return true
            }
          }
          return false
        }
        // A bow must also not swing into a NON-BUNDLE wire's corridor: a
        // bowed chord can cross a neighbor the straight chord cleared
        // (guards re-check the drawn shape against the neighbor's chord;
        // bundle-mates are excluded — the fan itself separates them).
        const mates = new Set(bundle)
        const bowPts = (w, kappa) => {
          const rx = w.x1 - w.x0
          const ry = w.y1 - w.y0
          const len = Math.hypot(rx, ry) || 1
          const b = (4 + 1.75 * Math.abs(kappa)) * Math.sign(kappa)
          const bb = Math.max(-len / 4, Math.min(len / 4, b))
          const cx = (w.x0 + w.x1) / 2 + 1.5 * bb * (ry / len)
          const cy = (w.y0 + w.y1) / 2 + 1.5 * bb * (-rx / len)
          const pts = []
          for (let i = 0; i <= 24; i++) {
            const t = i / 24
            const u = 1 - t
            pts.push([u * u * w.x0 + 2 * u * t * cx + t * t * w.x1,
              u * u * w.y0 + 2 * u * t * cy + t * t * w.y1])
          }
          return pts
        }
        const bowHitsWire = (w, kappa) => {
          const pw = bowPts(w, kappa)
          for (const d of directs) {
            if (mates.has(d)) continue
            const pd = samples(d.x0, d.y0, d.x1, d.y1)
            const sharedSrc = w.e.source === d.e.source && (w.e.sourceHandle || '') === (d.e.sourceHandle || '')
            const sharedTgt = w.e.target === d.e.target && (w.e.targetHandle || '') === (d.e.targetHandle || '')
            for (let m = 0; m <= 23; m++) {
              const jLo = sharedSrc ? 5 : 0
              const jHi = sharedTgt ? 10 : 15
              for (let n2 = jLo; n2 <= jHi; n2++) {
                if (segSeg(pw[m][0], pw[m][1], pw[m + 1][0], pw[m + 1][1],
                  pd[n2][0], pd[n2][1], pd[n2 + 1][0], pd[n2 + 1][1]) < 7) return true
              }
            }
          }
          return false
        }
        bundle.forEach((c, i) => {
          // Rank 0 stays flat; every later member bows AWAY from its
          // previous-rank neighbor, stacking on the opposite side of the
          // corridor. Fixed perp sign (no normal flip) keeps the fan
          // monotone — it never pushes a bow back onto rank i-1.
          let kappa
          if (i === 0) {
            kappa = 0
          } else {
            // side = cross(u, prev->c start): which side of THIS chord the
            // previous rank's start sits on. Bow to the opposite side.
            const prev = bundle[i - 1]
            const side = (c.x1 - c.x0) * (prev.y0 - c.y0) - (c.y1 - c.y0) * (prev.x0 - c.x0)
            kappa = (side > 0 ? 4 : -4) * i
          }
          const shrink = (k) => (k >= 4 ? k - 4 : (k <= -4 ? k + 4 : 0))
          while (kappa !== 0 && (quadHitsBox(c, kappa) || bowHitsWire(c, kappa))) {
            kappa = shrink(kappa)
          }
          c.kappa = kappa
        })
      }
    }
    // 5) Emit. Styles/markers by kind as before; the geometry knobs ride
    // in edge.data and are read only by KindEdge. Keyed by edge OBJECT —
    // some hand-authored flows (heptagon MHVw4/MHVw6) have edges with no
    // id, which would collapse to one shared undefined key and stamp every
    // such wire with the LAST wire's routing knobs.
    const wireByEdge = new Map(wires.map((w) => [w.e, w]))
    const eidOf = new Map(edges.map((e, i) => [e, e.id || `e${i}`]))
    const out = []
    for (const e of edges) {
      const gs = groupOf.get(e.source)
      const gt = groupOf.get(e.target)
      if (gs && gt && gs === gt) continue
      const kind = (gt && !gs ? null : resolveSourceKind(e.source, e.sourceHandle)) || 'any'
      const color = kindColor(kind)
      const w = wireByEdge.get(e)
      if (!gs && !gt && w) {
        out.push({
          ...e,
          id: eidOf.get(e),
          type: 'kind',
          data: {
            spread: w.spread,
            kappa: w.kappa,
            ...(w.detourY !== null ? {
              detourY: w.detourY, outboundX: w.outboundX, inboundX: w.inboundX,
            } : {}),
            ...(w.fallback ? { fallback: true } : {}),
          },
          style: { stroke: color, strokeWidth: 1.6, ...e.style },
          markerEnd: { type: MarkerType.ArrowClosed, color, width: 16, height: 16 },
        })
        continue
      }
      out.push({
        ...e,
        type: 'kind',
        id: `px_${eidOf.get(e)}`,
        source: gs || e.source,
        sourceHandle: gs ? 'out' : e.sourceHandle,
        target: gt || e.target,
        targetHandle: gt ? 'in' : e.targetHandle,
        selectable: false,
        focusable: false,
        style: { stroke: color, strokeWidth: 1.6, ...e.style, strokeDasharray: '6 3', opacity: 0.7 },
        markerEnd: { type: MarkerType.ArrowClosed, color, width: 16, height: 16 },
      })
    }
    return out
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
  const [exporting, setExporting] = useState(false)
  const [exportInfo, setExportInfo] = useState(null)
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
  const saveSession = useRef(null)
  const autosaveRef = useRef(autoSave)
  autosaveRef.current = autoSave
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
      const session = openFlowSave(project.id, flow)
      saveSession.current = session
      const recovered = session.draft
      const loadedFlow = recovered ? { ...flow, ...recovered } : flow
      if (recovered && !session.isSaved()) toast('Recovered an unsaved draft from this tab. Save to retry; a conflicting server edit will be protected.')
      skipAutosave.current = true
      setFlowName(loadedFlow.name)
      const gr0 = (loadedFlow.graph?.groups || []).map((g) => ({ collapsed: true, position: { x: 0, y: 0 }, ...g }))
      // default data: hand-authored or legacy graphs may omit it; OpNode/Inspector
      // dereference node.data unconditionally.
      let ns = (loadedFlow.graph?.nodes || []).map((n) => ({ ...n, data: n.data || {} }))
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
      const legacyEdges = loadedFlow.graph?.edges || []
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
      bumpNodeSeqFromIds([...ns.map((n) => n.id), ...(loadedFlow.graph?.edges || []).map((e) => e.id)])
      setCompileResult(null)
      loadedSig.current = JSON.stringify({
        name: loadedFlow.name || '',
        nodes: ns.filter((n) => n.type !== 'groupBox').map((n) => ({ id: n.id, type: n.type, position: n.position, data: n.data })),
        edges: addMig.map((e) => ({ id: e.id, source: e.source, sourceHandle: e.sourceHandle, target: e.target, targetHandle: e.targetHandle })),
        groups: gr.filter((g) => g.node_ids.some((id) => ns.some((n) => n.id === id && n.type !== 'groupBox'))),
      })
      setSaveState(session.isSaved() ? 'saved' : 'unsaved')
    }
  }, [project?.id, flow?.id]) // eslint-disable-line

  const serializeGraph = useCallback(() => ({
    nodes: nodes.filter((n) => n.type !== 'groupBox').map((n) => ({ id: n.id, type: n.type, position: n.position, data: n.data })),
    edges: edges.map((e) => ({ id: e.id, source: e.source, sourceHandle: e.sourceHandle, target: e.target, targetHandle: e.targetHandle })),
    groups: groups.filter((g) => g.node_ids.some((id) => nodes.some((n) => n.id === id && n.type !== 'groupBox'))),
  }), [nodes, edges, groups])

  useEffect(() => {
    const session = saveSession.current
    return () => {
      if (autosaveRef.current) session?.flush().catch((e) => toast('Draft preserved; save failed: ' + e.message))
    }
  }, [project?.id, fid]) // eslint-disable-line

  useEffect(() => {
    if (!flow || !saveSession.current) return undefined
    if (skipAutosave.current) { skipAutosave.current = false; return undefined }
    const graph = serializeGraph()
    const session = saveSession.current
    try { session.remember({ name: flowName, graph }) }
    catch { toast('Browser draft storage is full or unavailable. Keep this page open until saved.') }
    if (session.isSaved()) { setSaveState('saved'); return undefined }
    setSaveState('unsaved')
    if (!autoSave) return undefined
    let alive = true
    const t = setTimeout(async () => {
      setSaveState('saving')
      try {
        await session.flush()
        if (alive) setSaveState(session.isSaved() ? 'saved' : 'unsaved')
      } catch (e) {
        if (alive) { setSaveState('error'); toast('Draft preserved; save failed: ' + e.message) }
      }
    }, 800)
    return () => { alive = false; clearTimeout(t) }
  }, [nodes, edges, flowName, groups, autoSave, flow]) // eslint-disable-line

  useEffect(() => {
    const guard = (e) => {
      if (saveSession.current && !saveSession.current.isSaved()) {
        e.preventDefault()
        e.returnValue = ''
      }
    }
    window.addEventListener('beforeunload', guard)
    return () => window.removeEventListener('beforeunload', guard)
  }, [])

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
      case 'squeeze_tensor': {
        const t = inDimsList('in')
        if (t && t.length >= 2) {
          const out = t.filter((x) => x !== 1)
          if (out.length >= 2) return out.join('×')
        }
        break
      }
      case 'shuffle_product': {
        const a = inDimsList('a')
        const b = inDimsList('b')
        if (a && b) return [...a, ...b].join('×')
        break
      }
      case 'expand_tensor': {
        const bases = edges
          .filter((ed) => ed.target === nodeId && /^in_\d+$/.test(ed.targetHandle || '') && ed.targetHandle !== 'in_0')
          .map((ed) => ({
            idx: parseInt((ed.targetHandle || '').slice(3), 10),
            dims: (inputDimsFor(ed.source, ed.sourceHandle, depth + 1) || '').split('×').map((x) => parseInt(x, 10)),
          }))
          .filter((b) => b.dims.length > 0 && b.dims.every((x) => Number.isFinite(x)))
          .sort((x, y) => x.idx - y.idx)
        const t = inDimsList('in_0')
        if (t && t.length >= 2 && bases.length) {
          let fec = t[t.length - 2]
          const letters = [t[t.length - 1]]
          for (const b of bases) {
            if (b.dims.length !== 3 || b.dims[0] !== fec) return null
            letters.push(b.dims[2])
            fec = b.dims[1]
          }
          return [fec, ...letters].join('×')
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
      try { saveSession.current.remember({ name: flowName, graph }) }
      catch { toast('Browser draft recovery is unavailable; saving to the server.') }
      await saveSession.current.flush()
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
      await save()
      const { run_id } = await api.runFlow(project.id, fid)
      navigate(`/runs/${run_id}`)
    } catch (e) {
      setRunning(false)
      const errs = e.payload?.detail?.errors
      if (errs) { setCompileResult({ ok: false, errors: errs }); setSideTab('plan') }
      toast('Run failed to start: ' + (errs ? 'see errors panel' : e.message))
    }
  }

  const exportScript = async () => {
    setExporting(true)
    try {
      // export recompiles server-side; make sure the latest graph is saved first
      await save()
      const info = await api.exportFlowScript(project.id, fid)
      setExportInfo(info)
      toast(`Standalone script written (${info.n_steps} steps)`)
    } catch (e) {
      const errs = e.payload?.detail?.errors
      if (errs) { setCompileResult({ ok: false, errors: errs }); setSideTab('plan') }
      setExportInfo({ error: errs ? 'flow does not compile — see errors panel' : (e.message || 'export failed') })
    } finally {
      setExporting(false)
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
    const node = nodes.find((n) => n.id === id)
    if (node?.type === 'solve_collinear') {
      // Pair rows shifted (a middle row removed): renumber in_seed_N /
      // in_rhs_N edges to follow their config row; wires of removed rows
      // are pruned so the ports really disappear. Only a shrinking array
      // means a removal — edits/adds keep every wire.
      if (Array.isArray(patch.pairs) && patch.pairs.length < normalizeCollinearPairs(node.data).length) {
        const sig = (r) => `${String(r?.projection ?? '').trim()}||${String(r?.rhs ?? '').trim()}`
        const oldRows = normalizeCollinearPairs(node.data).slice(1)
        const newRows = patch.pairs.slice(1)
        const mapping = new Map()
        const removed = []
        let cursor = 0
        for (let o = 0; o < oldRows.length; o++) {
          let matched = -1
          for (let n = cursor; n < newRows.length; n++) {
            if (sig(newRows[n]) === sig(oldRows[o])) { matched = n; break }
          }
          if (matched >= 0) { mapping.set(o + 1, matched + 1); cursor = matched + 1 }
          else removed.push(o + 1)
        }
        if (removed.length || [...mapping.entries()].some(([o, n]) => o !== n)) {
          setEdges((es) => es.map((e) => {
            if (e.target !== id) return e
            const m = (e.targetHandle || '').match(/^in_(seed|rhs)_(\d+)$/)
            if (!m) return e
            const oldN = parseInt(m[2], 10)
            const newN = mapping.get(oldN)
            if (newN === undefined) return null
            return newN === oldN ? e : { ...e, targetHandle: `in_${m[1]}_${newN}` }
          }).filter(Boolean))
        }
      }
      // cond toggle-off drops the port: remove its wires too, otherwise the
      // wired edge would keep the port alive and the toggle would look broken.
      if (patch.cond_enabled === false) {
        setEdges((es) => es.filter((e) => !(e.target === id && (e.targetHandle || '') === 'cond')))
      }
    }
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

  const visibleEdges = useMemo(() => routeEdges({
    edges, nodes, project, groupOf, resolveSourceKind,
  }), [edges, nodes, project, groupOf, resolveSourceKind])

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
          <input aria-label="Flow name" style={{ width: 180 }} value={flowName} onChange={(e) => setFlowName(e.target.value)} />
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
          <button onClick={() => save().catch((e) => toast('Draft preserved; save failed: ' + e.message))}>Save</button>
          {saveState === 'error' && <>
            <button onClick={() => saveSession.current.download()}>Download draft</button>
            <button onClick={() => {
              if (window.confirm('Discard this tab’s draft and load the currently saved flow? Download the draft first to keep a copy.')) {
                saveSession.current.discard()
                window.location.reload()
              }
            }}>Reload saved flow</button>
          </>}
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
          edgeTypes={edgeTypes}
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
            ? <CompilePanel result={compileResult} onRun={run} onExport={exportScript} running={running} exporting={exporting} exportInfo={exportInfo} />
            : <div>
              <p className="muted">Press “Compile” to turn the diagram into the exact command sequence. You can inspect every command before running it.</p>
              <p className="muted">After compiling you can <strong>▶ Run</strong> the plan here, or <strong>⇪ Export standalone script</strong> to get a portable bash script you can carry to another machine (see README → “Design locally, run on the cluster”).</p>
            </div>
        )}
      </div>
    </div>
  )
}

export { OpNode, Inspector, nodeTypes, KindEdge }
export default function FlowEditor() {
  return (
    <ReactFlowProvider>
      <FlowEditorInner />
    </ReactFlowProvider>
  )
}
