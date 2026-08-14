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
        exact command sequence, then run it.
      </p>
      <div className="card" style={{ marginBottom: 16 }}>
        <div className="row">
          <input placeholder="New flow name…" value={name} onChange={(e) => setName(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && create()} />
          <button className="primary shrink" onClick={create}>+ New flow</button>
        </div>
      </div>
      {!project.flows.length && (
        <div className="empty-state"><div className="big">No flows yet</div><p>Create your first flow to wire operations together.</p></div>
      )}
      <div className="grid">
        {project.flows.map((f) => (
          <div key={f.id} className="card clickable" onClick={() => navigate(`/flows/${f.id}`)}>
            <h2 style={{ marginTop: 0 }}>{f.name}</h2>
            <p className="muted">{(f.graph?.nodes || []).length} nodes, {(f.graph?.edges || []).length} connections</p>
            <button className="small danger" onClick={(e) => { e.stopPropagation(); remove(f.id) }}>Delete</button>
          </div>
        ))}
      </div>
    </div>
  )
}
