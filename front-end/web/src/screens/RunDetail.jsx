import React, { useEffect, useRef, useState } from 'react'
import { useParams } from 'react-router-dom'
import { api } from '../api'
import { useToast } from '../App'

export default function RunDetail() {
  const { rid } = useParams()
  const toast = useToast()
  const [run, setRun] = useState(null)
  const [logs, setLogs] = useState([])
  const [error, setError] = useState('')
  const [connection, setConnection] = useState('connecting')
  const [retry, setRetry] = useState(0)
  const logRef = useRef(null)
  const esRef = useRef(null)

  useEffect(() => {
    let alive = true
    let refreshing = false
    let refreshAgain = false
    let pending = []
    setRun(null)
    setLogs([])
    setError('')
    setConnection('connecting')
    const es = new EventSource(api.runEventsUrl(rid))
    esRef.current = es
    const refresh = async () => {
      if (refreshing) { refreshAgain = true; return }
      refreshing = true
      try {
        const r = await api.getRun(rid)
        if (alive) { setRun(r); setError('') }
      } catch (e) {
        if (alive) { setError(e.message); if (e.status === 404) es.close() }
      } finally {
        refreshing = false
        if (alive && refreshAgain) { refreshAgain = false; refresh() }
      }
    }
    refresh()
    const flush = () => {
      if (alive && pending.length) {
        const batch = pending
        pending = []
        setLogs((ls) => [...ls, ...batch].slice(-4000))
      }
    }
    const timer = setInterval(flush, 100)
    es.onopen = () => { if (alive) setConnection('connected') }
    es.onerror = () => { if (alive) { setConnection('disconnected; reconnecting…'); refresh() } }
    es.addEventListener('log', (e) => { pending.push(JSON.parse(e.data).line); if (pending.length > 4000) pending.shift() })
    es.addEventListener('gap', () => { pending.push('[Older output is available in the saved run log.]') })
    es.addEventListener('step', refresh)
    es.addEventListener('status', refresh)
    es.addEventListener('end', () => { es.close(); flush(); refresh(); if (alive) setConnection('finished') })
    return () => { alive = false; clearInterval(timer); es.close() }
  }, [rid, retry])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [logs])

  const cancel = async () => {
    try { await api.cancelRun(rid); toast('Run cancelled.') } catch (e) { toast(e.message) }
  }

  if (!run) return <div className="page"><p className={error ? 'error-text' : 'muted'}>{error || 'Loading run…'}</p>{error && <button onClick={() => setRetry((v) => v + 1)}>Retry</button>}</div>
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
      {error && <p role="alert" className="error-text">{error}</p>}
      <p className="muted">Log connection: {connection}</p>
      <h2>Log</h2>
      <div className="log-console" ref={logRef}>
        {logs.length ? logs.join('\n') : <span className="muted">waiting for output…</span>}
      </div>
    </div>
  )
}
