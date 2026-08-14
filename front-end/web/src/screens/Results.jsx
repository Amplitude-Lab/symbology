import React, { useEffect, useState } from 'react'
import { api } from '../api'
import { useProject, useToast } from '../App'

function TensorInspector({ file, onClose }) {
  const { project } = useProject()
  const [summary, setSummary] = useState(null)
  const [error, setError] = useState(null)

  useEffect(() => {
    setSummary(null); setError(null)
    api.tensorSummary(project.id, file).then(setSummary).catch((e) => setError(e.message))
  }, [file, project.id])

  return (
    <div className="card" style={{ marginTop: 16 }}>
      <div className="row">
        <h2 style={{ margin: 0, flex: 1 }} className="mono">{file}</h2>
        <button className="small shrink" onClick={onClose}>Close</button>
      </div>
      {!summary && !error && <p className="muted">Loading tensor summary (this calls Wolfram, may take a few seconds)…</p>}
      {error && <p className="error-text">{error}</p>}
      {summary && (
        <>
          <p>
            <span className="badge">dims: {summary.dims?.join(' × ')}</span>{' '}
            <span className="badge">nnz: {summary.nnz}</span>
          </p>
          {!!summary.sample?.length && (
            <table>
              <thead><tr><th>Index</th><th>Value</th></tr></thead>
              <tbody>
                {summary.sample.map((s, i) => (
                  <tr key={i}><td className="mono">[{s.index.join(', ')}]</td><td className="mono">{s.value}</td></tr>
                ))}
              </tbody>
            </table>
          )}
          {!summary.sample?.length && <p className="muted">The tensor is empty.</p>}
        </>
      )}
    </div>
  )
}

export default function Results() {
  const { project } = useProject()
  const toast = useToast()
  const [tab, setTab] = useState('output')
  const [files, setFiles] = useState([])
  const [inspect, setInspect] = useState(null)

  useEffect(() => {
    if (!project) return
    api.listTensors(project.id, tab).then(setFiles).catch(() => setFiles([]))
  }, [project, tab])

  if (!project) return null

  return (
    <div className="page">
      <h1>Results</h1>
      <div className="row" style={{ marginBottom: 12, maxWidth: 420 }}>
        <button className={tab === 'output' ? 'primary' : ''} onClick={() => { setTab('output'); setInspect(null) }}>output/</button>
        <button className={tab === 'data' ? 'primary' : ''} onClick={() => { setTab('data'); setInspect(null) }}>data/</button>
      </div>
      <div className="card">
        {!files.length && <p className="muted">No tensor files in {tab}/ yet.</p>}
        {!!files.length && (
          <table>
            <thead><tr><th>File</th><th>Size</th><th></th></tr></thead>
            <tbody>
              {files.map((f) => (
                <tr key={f.file}>
                  <td className="mono">{f.file}</td>
                  <td className="muted">{(f.size_bytes / 1024).toFixed(1)} KB</td>
                  <td style={{ textAlign: 'right' }}>
                    <button className="small" onClick={() => setInspect(f.file)}>Inspect</button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </div>
      {inspect && <TensorInspector file={inspect} onClose={() => setInspect(null)} />}
    </div>
  )
}
