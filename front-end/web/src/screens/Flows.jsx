import React, { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { api } from '../api'
import { useProject, useToast } from '../App'

export default function Flows() {
  const { project, refreshProject } = useProject()
  const toast = useToast()
  const navigate = useNavigate()
  const [name, setName] = useState('')

  if (!project) return null

  const create = async () => {
    const flow = await api.createFlow(project.id, name.trim() || 'Untitled flow')
    await refreshProject()
    navigate(`/flows/${flow.id}`)
  }

  const createBlock = async () => {
    const flow = await api.createFlow(project.id, name.trim() || 'New custom block', { custom_block: true })
    await refreshProject()
    navigate(`/flows/${flow.id}`)
  }

  const toggleBlock = async (f) => {
    if (!f.custom_block) {
      const nodes = f.graph?.nodes || []
      const edges = f.graph?.edges || []
      const ins = nodes.filter((n) => n.type === 'cb_in')
      const outs = nodes.filter((n) => n.type === 'cb_out')
      const problems = []
      if (!ins.length) problems.push('at least one Input port')
      if (!outs.length) problems.push('at least one Output port')
      const srcs = new Set(edges.map((e) => e.source))
      const tgts = new Set(edges.map((e) => e.target))
      if (ins.some((n) => !srcs.has(n.id))) problems.push('every Input port wired onwards')
      if (outs.some((n) => !tgts.has(n.id))) problems.push('every Output port wired from something')
      const unnamed = [...ins, ...outs].filter((n) => !(n.data?.name || '').trim())
      if (unnamed.length) problems.push('all ports named')
      if (problems.length) {
        toast(`Cannot seal yet — the block needs: ${problems.join('; ')}.`)
        return
      }
    }
    await api.updateFlow(project.id, f.id, { custom_block: !f.custom_block })
    await refreshProject()
    toast(f.custom_block ? 'Converted back to a normal flow.' : 'Sealed as a custom block — it now appears in the flow editor palette.')
  }

  const duplicate = async (f) => {
    const dup = await api.duplicateFlow(project.id, f.id)
    await refreshProject()
    navigate(`/flows/${dup.id}`)
  }

  const remove = async (fid) => {
    await api.deleteFlow(project.id, fid)
    await refreshProject()
    toast('Flow deleted.')
  }

  return (
    <div className="page">
      <h1>Flows</h1>
      <p className="muted">
        A flow is a visual program: wire alphabet building blocks into operations, compile it into the
        exact command sequence, then run it. Open a flow to draw it on the canvas — the graph autosaves
        as you work.
      </p>
      <p className="muted" style={{ fontSize: 11 }}>
        How to run a flow: in the editor click <strong>Compile</strong>, then open the <strong>Plan</strong> tab
        on the right — it shows the exact command list the flow compiles to. From there you can
        <strong> ▶ Run this plan</strong> on this machine, or <strong>⇪ Export standalone script</strong> to
        write the same plan as a portable, self-checking bash script (<code>exported/&lt;flow&gt;.sh</code> inside
        this project). Copy that script together with the project directory to any machine where the C++ core
        is built and run it there — no front-end needed. See README → “Design locally, run on the cluster”.
      </p>
      <div className="card" style={{ marginBottom: 16 }}>
        <div className="row">
          <input placeholder="New flow name…" value={name} onChange={(e) => setName(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && create()} />
          <button className="primary shrink" onClick={create}>+ New flow</button>
          <button className="shrink" onClick={createBlock}>+ New custom block</button>
        </div>
        <p className="muted" style={{ marginBottom: 0, fontSize: 11 }}>
          A custom block is a sub-diagram with abstract Input/Output ports. Once it contains wired Input and
          Output nodes, seal it — it becomes a single reusable block in every flow&apos;s palette.
        </p>
      </div>
      {!project.flows.length && (
        <div className="empty-state"><div className="big">No flows yet</div><p>Create your first flow to wire operations together.</p></div>
      )}
      <div className="grid">
        {project.flows.map((f) => (
          <div key={f.id} className="card clickable" onClick={() => navigate(`/flows/${f.id}`)}>
            <h2 style={{ marginTop: 0 }}>{f.name}{f.custom_block && <span className="muted" style={{ fontSize: 11, fontWeight: 400 }}> · custom block</span>}</h2>
            <p className="muted">{(f.graph?.nodes || []).length} nodes, {(f.graph?.edges || []).length} connections</p>
            <div className="row">
              <button className="small" onClick={(e) => { e.stopPropagation(); duplicate(f) }}>Duplicate</button>
              <button className="small" onClick={(e) => { e.stopPropagation(); toggleBlock(f) }}>
                {f.custom_block ? 'Unseal (normal flow)' : 'Seal as custom block'}
              </button>
              <button className="small danger" onClick={(e) => { e.stopPropagation(); remove(f.id) }}>Delete</button>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
