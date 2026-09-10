import React, { createContext, useCallback, useContext, useEffect, useState } from 'react'
import { NavLink, Navigate, Route, Routes, useNavigate } from 'react-router-dom'
import { api } from './api'
import Materials from './screens/Materials'
import Flows from './screens/Flows'
import FlowEditor from './screens/FlowEditor'
import Runs from './screens/Runs'
import RunDetail from './screens/RunDetail'
import Results from './screens/Results'

const ProjectCtx = createContext(null)
export const useProject = () => useContext(ProjectCtx)

const ToastCtx = createContext(() => {})
export const useToast = () => useContext(ToastCtx)

class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props)
    this.state = { error: null }
  }
  static getDerivedStateFromError(error) {
    return { error }
  }
  componentDidCatch(error, info) {
    try {
      fetch('/api/client-log', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          kind: 'render',
          message: String(error?.message || error).slice(0, 4000),
          stack: String(error?.stack || '').slice(0, 4000) + '\n' + String(info?.componentStack || '').slice(0, 2000),
          url: window.location.href,
          ua: navigator.userAgent,
        }),
      }).catch(() => {})
    } catch { /* best effort */ }
  }
  render() {
    if (this.state.error) {
      return (
        <div className="page">
          <div className="empty-state">
            <div className="big">Something went wrong rendering this page</div>
            <p className="error-text mono" style={{ whiteSpace: 'pre-wrap', maxWidth: 640 }}>{String(this.state.error?.message || this.state.error)}</p>
            <button className="primary" onClick={() => this.setState({ error: null })}>Try again</button>
          </div>
        </div>
      )
    }
    return this.props.children
  }
}

function NewProjectModal({ templates, onClose, onCreated }) {
  const [name, setName] = useState('')
  const [templateId, setTemplateId] = useState('')
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState(null)
  const create = async () => {
    if (!name.trim()) return setError('Please give the project a name.')
    setBusy(true); setError(null)
    try {
      const p = await api.createProject(name.trim(), templateId || null)
      onCreated(p)
    } catch (e) { setError(e.message); setBusy(false) }
  }
  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>New project</h3>
        <label>Project name</label>
        <input autoFocus value={name} onChange={(e) => setName(e.target.value)} placeholder="e.g. 4p form factor study" />
        <label>Start from template (optional)</label>
        <select value={templateId} onChange={(e) => setTemplateId(e.target.value)}>
          <option value="">Empty project</option>
          {templates.map((t) => (
            <option key={t.id} value={t.id}>{t.name}</option>
          ))}
        </select>
        {templateId && <p className="muted" style={{ fontSize: 12 }}>{templates.find((t) => t.id === templateId)?.description}</p>}
        {error && <p className="error-text">{error}</p>}
        <div className="actions">
          <button onClick={onClose}>Cancel</button>
          <button className="primary" disabled={busy} onClick={create}>{busy ? 'Creating…' : 'Create'}</button>
        </div>
      </div>
    </div>
  )
}

export default function App() {
  const [env, setEnv] = useState(null)
  const [projects, setProjects] = useState([])
  const [templates, setTemplates] = useState([])
  const [currentId, setCurrentId] = useState(() => localStorage.getItem('currentProject') || '')
  const [project, setProject] = useState(null)
  const [showNew, setShowNew] = useState(false)
  const [toast, setToast] = useState(null)
  const navigate = useNavigate()

  const showToast = useCallback((msg) => {
    setToast(msg)
    setTimeout(() => setToast(null), 3500)
  }, [])

  const refreshProjects = useCallback(async () => {
    const list = await api.listProjects()
    setProjects(list)
    return list
  }, [])

  const refreshProject = useCallback(async (pid) => {
    const id = pid || currentId
    if (!id) { setProject(null); return null }
    const p = await api.getProject(id)
    setProject(p)
    return p
  }, [currentId])

  useEffect(() => {
    api.env().then(setEnv).catch(() => setEnv(null))
    api.templates().then(setTemplates).catch(() => {})
    refreshProjects().catch(() => {})
  }, [refreshProjects])

  useEffect(() => {
    if (currentId) {
      localStorage.setItem('currentProject', currentId)
      refreshProject(currentId).catch(() => { setProject(null); setCurrentId('') })
    } else {
      localStorage.removeItem('currentProject')
      setProject(null)
    }
  }, [currentId]) // eslint-disable-line

  const envProblem = env && (!env.wolframscript.found || !env.bootstrap.found)

  return (
    <ProjectCtx.Provider value={{ project, currentId, setCurrentId, refreshProject, refreshProjects }}>
      <ToastCtx.Provider value={showToast}>
        <div className="topbar">
          <span className="brand">Symbology Studio</span>
          <nav>
            <NavLink to="/materials">Materials</NavLink>
            <NavLink to="/flows">Flows</NavLink>
            <NavLink to="/runs">Runs</NavLink>
            <NavLink to="/results">Results</NavLink>
          </nav>
          <div className="project-picker">
            <select value={currentId} onChange={(e) => setCurrentId(e.target.value)}>
              <option value="">— select project —</option>
              {projects.map((p) => <option key={p.id} value={p.id}>{p.name}</option>)}
            </select>
            <button className="primary small" onClick={() => setShowNew(true)}>+ New</button>
          </div>
        </div>
        {envProblem && (
          <div className="env-banner">
            Environment incomplete:
            {!env.wolframscript.found && ' wolframscript was not found.'}
            {!env.bootstrap.found && ' the bootstrap binary was not found (run `make bootstrap` in the repository root).'}
            {' '}Property computation and C++ steps will fail until this is fixed.
          </div>
        )}
        {!currentId ? (
          <div className="page">
            <div className="empty-state">
              <div className="big">Welcome to Symbology Studio</div>
              <p>Create a project to start. The 4-point form factor template is a good first choice — it ships with a ready-made alphabet and precomputed tensors.</p>
              <button className="primary" onClick={() => setShowNew(true)}>Create your first project</button>
            </div>
          </div>
        ) : (
          <ErrorBoundary>
            <Routes>
              <Route path="/materials" element={<Materials />} />
              <Route path="/flows" element={<Flows />} />
              <Route path="/flows/:fid" element={<FlowEditor />} />
              <Route path="/runs" element={<Runs />} />
              <Route path="/runs/:rid" element={<RunDetail />} />
              <Route path="/results" element={<Results />} />
              <Route path="*" element={<Navigate to="/materials" replace />} />
            </Routes>
          </ErrorBoundary>
        )}
        {showNew && (
          <NewProjectModal
            templates={templates}
            onClose={() => setShowNew(false)}
            onCreated={(p) => { setShowNew(false); refreshProjects(); setCurrentId(p.id); navigate('/materials') }}
          />
        )}
        {toast && <div className="toast">{toast}</div>}
      </ToastCtx.Provider>
    </ProjectCtx.Provider>
  )
}
