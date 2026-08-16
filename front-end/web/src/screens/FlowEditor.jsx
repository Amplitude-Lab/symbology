import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { useNavigate, useParams } from 'react-router-dom'
import ReactFlow, {
  Background, Controls, Handle, Position, ReactFlowProvider,
  addEdge, useEdgesState, useNodesState,
} from 'reactflow'
import 'reactflow/dist/style.css'
import { api } from '../api'
import { useProject, useToast } from '../App'
import {
  NODE_DEFS, PALETTE, PROP_KIND, propLabel,
  kindsCompatible, sourceKindFor, targetKindFor,
} from '../flowdefs'

function StatusDot({ status }) {
  return <span className={`dot ${status || 'pending'}`} />
}

function AlphabetNode({ id, data, selected }) {
  const { project } = useProject()
  const alphabet = project?.alphabets.find((a) => a.id === data.alphabet_id)
  const selectedProps = data.selected_properties || []
  return (
    <div className="node-card" style={{ borderColor: selected ? 'var(--accent)' : undefined }}>
      <div className="node-title" style={{ background: NODE_DEFS.alphabet.color }}>
        Alphabet: {alphabet ? alphabet.name : '—'}
      </div>
      <div className="node-body">
        {!alphabet && <span>select in inspector →</span>}
        {alphabet && selectedProps.length === 0 && <span>no outputs selected</span>}
        {alphabet && alphabet.properties.filter((p) => selectedProps.includes(p.id)).map((p, i) => (
          <div key={p.id} className="handle-row" style={{ textAlign: 'right' }}>
            <span className={`kind-${PROP_KIND[p.type]}`}>
              {propLabel(p)}
            </span>{' '}
            <StatusDot status={p.status} />
            <Handle
              type="source" position={Position.Right} id={`prop_${p.id}`}
              style={{ top: 52 + i * 22, background: 'var(--accent)' }}
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
            style={{ top: 58 + i * 20 }}
          />
        ))}
        {outputs.map((out, i) => (
          <Handle
            key={out.id} type="source" position={Position.Right} id={out.id}
            style={{ top: 58 + i * 20, background: 'var(--accent)' }}
          />
        ))}
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
}

function Inspector({ node, onChange, onDelete }) {
  const { project } = useProject()
  if (!node) return (
    <div>
      <p className="muted">Select a node to edit its parameters.</p>
      <p className="muted" style={{ fontSize: 11 }}>Tip: click a node or edge, then press Delete/Backspace — or use the Delete button here — to remove it.</p>
    </div>
  )
  const d = node.data || {}
  const set = (patch) => onChange(node.id, patch)

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
            <p className="muted" style={{ fontSize: 11 }}>Properties that are not ready will be computed automatically when the flow runs.</p>
          </>
        )}
        <div style={{ marginTop: 16, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
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
      {node.type === 'merge_conditions' && <p className="muted">Connect two or more dlogmat outputs (integrability, extended Steinmann, cluster adjacency) to merge them into a single condition tensor.</p>}
      {node.type === 'sew' && <p className="muted">Combines a condition tensor, an FEC tensor of weight F and an LEC tensor of weight L into SEW_FpL.</p>}
      <div style={{ marginTop: 16, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
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
  const rf = useRef(null)
  const wrapper = useRef(null)
  const skipAutosave = useRef(true)

  const flow = project?.flows.find((f) => f.id === fid)

  useEffect(() => {
    if (flow) {
      skipAutosave.current = true
      setFlowName(flow.name)
      setNodes(flow.graph?.nodes || [])
      setEdges(flow.graph?.edges || [])
      setCompileResult(null)
      setSaveState('saved')
    }
  }, [flow?.id]) // eslint-disable-line

  const serializeGraph = useCallback(() => ({
    nodes: nodes.map((n) => ({ id: n.id, type: n.type, position: n.position, data: n.data })),
    edges: edges.map((e) => ({ id: e.id, source: e.source, sourceHandle: e.sourceHandle, target: e.target, targetHandle: e.targetHandle })),
  }), [nodes, edges])

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
  }, [nodes, edges, flowName]) // eslint-disable-line

  const isValidConnection = useCallback((conn) => {
    const sn = nodes.find((n) => n.id === conn.source)
    const tn = nodes.find((n) => n.id === conn.target)
    const sk = sourceKindFor(sn, conn.sourceHandle, project)
    const tk = targetKindFor(tn, conn.targetHandle)
    return !!(sk && tk && kindsCompatible(sk, tk))
  }, [nodes, project])

  const onConnect = useCallback((conn) => {
    if (!isValidConnection(conn)) {
      toast('Those ports are not compatible.')
      return
    }
    setEdges((eds) => addEdge({ ...conn, animated: false }, eds))
    setCompileResult(null)
  }, [isValidConnection, setEdges, toast])

  const onDrop = useCallback((e) => {
    e.preventDefault()
    const type = e.dataTransfer.getData('application/symbology-node')
    if (!type || !rf.current) return
    const pos = rf.current.screenToFlowPosition({ x: e.clientX, y: e.clientY })
    const id = `${type}_${nodeSeq++}`
    const defaults = type === 'alphabet' ? { alphabet_id: '', selected_properties: [] } : {}
    setNodes((ns) => [...ns, { id, type, position: pos, data: defaults }])
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
    await save()
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
    setCompileResult(null)
  }, [setNodes])

  const deleteNode = useCallback((id) => {
    setNodes((ns) => ns.filter((n) => n.id !== id))
    setEdges((es) => es.filter((e) => e.source !== id && e.target !== id))
    setSelectedId(null)
    setCompileResult(null)
    toast('Node deleted.')
  }, [setNodes, setEdges, toast])

  useEffect(() => {
    const handler = (e) => {
      if (e.key !== 'Delete' && e.key !== 'Backspace') return
      const t = e.target
      if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable)) return
      const nodeIds = new Set(nodes.filter((n) => n.selected).map((n) => n.id))
      const edgeIds = new Set(edges.filter((ed) => ed.selected).map((ed) => ed.id))
      if (!nodeIds.size && !edgeIds.size) return
      e.preventDefault()
      setNodes((ns) => ns.filter((n) => !nodeIds.has(n.id)))
      setEdges((es) => es.filter((ed) => !edgeIds.has(ed.id) && !nodeIds.has(ed.source) && !nodeIds.has(ed.target)))
      setSelectedId(null)
      setCompileResult(null)
    }
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [nodes, edges, setNodes, setEdges])

  const selectedNode = nodes.find((n) => n.id === selectedId)

  if (!flow) return <div className="page"><p className="muted">Loading flow…</p></div>

  return (
    <div className="flow-layout">
      <div className="flow-palette">
        <div className="muted" style={{ fontSize: 11, marginBottom: 8, textTransform: 'uppercase', letterSpacing: 0.5 }}>Nodes</div>
        {PALETTE.map((p) => (
          <div
            key={p.type} className="palette-item" draggable
            onDragStart={(e) => e.dataTransfer.setData('application/symbology-node', p.type)}
          >
            {p.label}
            <div className="sub">{p.sub}</div>
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
          nodes={nodes} edges={edges}
          onNodesChange={(chg) => { onNodesChange(chg); setCompileResult(null) }}
          onEdgesChange={(chg) => { onEdgesChange(chg); setCompileResult(null) }}
          onConnect={onConnect}
          isValidConnection={isValidConnection}
          nodeTypes={nodeTypes}
          onInit={(inst) => { rf.current = inst }}
          onSelectionChange={(sel) => setSelectedId(sel.nodes[0]?.id || null)}
          deleteKeyCode={['Backspace', 'Delete']}
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
        {sideTab === 'inspector' && <Inspector node={selectedNode} onChange={updateNodeData} onDelete={deleteNode} />}
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
