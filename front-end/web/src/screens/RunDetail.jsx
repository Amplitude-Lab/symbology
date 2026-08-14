import React, { useEffect, useRef, useState } from 'react'
import { useParams } from 'react-router-dom'
import { api } from '../api'
import { useToast } from '../App'

export default function RunDetail() {
  const { rid } = useParams()
  const toast = useToast()
  const [run, setRun] = useState(null)
  const [logs, setLogs] = useState([])
  const logRef = useRef(null)
  const esRef = useRef(null)

  useEffect(() => {
    let alive = true
    api.getRun(rid).then((r) => alive && setRun(r)).catch(() => {})
    const es = new EventSource(api.runEventsUrl(rid))
    esRef.current = es
    es.addEventListener('log', (e) => {
      const d = JSON.parse(e.data)
      setLogs((ls) => [...ls.slice(-4000), d.line])
    })
    es.addEventListener('step', () => { api.getRun(rid).then(setRun).catch(() => {}) })
    es.addEventListener('status', () => { api.getRun(rid).then(setRun).catch(() => {}) })
    es.addEventListener('end', () => { es.close() })
    return () => { alive = false; es.close() }
  }, [rid])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [logs])

  const cancel = async () => {
    try { await api.cancelRun(rid); toast('Run cancelled.') } catch (e) { toast(e.message) }
  }

  if (!run) return <div className="page"><p className="muted">Loading run…</p></div>
  const active = run.status === 'queued' || run.status === 'running'

  return (
    <div className="page">
      <div className="row" style={{ marginBottom: 16 }}>
        <h1 style={{ margin: 0, flex: 1 }}>{run.label} <span className={`status-${run.status}`} style={{ fontSize: 14 }}>● {run.status}</span></h1>
        {active && <button className="danger shrink" onClick={cancel}>Cancel run</button>}
      </div>
      <div className="card" style={{ marginBottom: 16 }}>
        <ul className="step-list">
          {run.steps.map((s) => (
            <li key={s.id}>
              <span className={`kind-badge ${s.kind}`}>{s.kind}</span>
              <span className={`status-${s.status}`}>{s.status === 'running' ? '▶' : s.status === 'done' ? '✓' : s.status === 'failed' ? '✗' : s.status === 'skipped' ? '↷' : '○'}</span>{' '}
              {s.label}
              <div className="step-cmd">$ {s.command}</div>
            </li>
          ))}
        </ul>
      </div>
      <h2>Log</h2>
      <div className="log-console" ref={logRef}>
        {logs.length ? logs.join('\n') : <span className="muted">waiting for output…</span>}
      </div>
    </div>
  )
}
