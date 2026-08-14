import React, { useEffect, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { api } from '../api'
import { useProject } from '../App'

export default function Runs() {
  const { project } = useProject()
  const navigate = useNavigate()
  const [runs, setRuns] = useState([])

  useEffect(() => {
    if (!project) return
    let alive = true
    const load = () => api.listRuns(project.id).then((r) => alive && setRuns(r)).catch(() => {})
    load()
    const t = setInterval(load, 3000)
    return () => { alive = false; clearInterval(t) }
  }, [project])

  if (!project) return null

  return (
    <div className="page">
      <h1>Runs</h1>
      {!runs.length && <div className="empty-state"><div className="big">No runs yet</div><p>Compile a flow and run it — progress and logs will appear here.</p></div>}
      {!!runs.length && (
        <div className="card">
          <table>
            <thead>
              <tr><th>Run</th><th>Status</th><th>Started</th><th>Current step</th></tr>
            </thead>
            <tbody>
              {runs.map((r) => (
                <tr key={r.run_id} style={{ cursor: 'pointer' }} onClick={() => navigate(`/runs/${r.run_id}`)}>
                  <td>{r.label}</td>
                  <td className={`status-${r.status}`}>{r.status}</td>
                  <td className="muted">{new Date(r.created_at).toLocaleString()}</td>
                  <td className="muted">{r.current_step || '—'}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}
