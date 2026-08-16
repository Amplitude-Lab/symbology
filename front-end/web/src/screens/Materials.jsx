import React, { useState } from 'react'
import { api } from '../api'
import { useProject, useToast } from '../App'

const PROP_LABEL = {
  integrability: 'Integrability (dlogmat)',
  first_entry: 'First Entry',
  last_entry: 'Last Entry',
  extended_steinmann: 'Extended Steinmann',
  cluster_adjacency: 'Cluster Adjacency',
  transformation: 'Transformation',
}

function StatusDot({ status }) {
  return <span className={`dot ${status || 'pending'}`} />
}

function LetterMultiSelect({ letters, value, onChange }) {
  const toggle = (l) => onChange(value.includes(l) ? value.filter((x) => x !== l) : [...value, l])
  return (
    <div style={{ maxHeight: 180, overflow: 'auto', border: '1px solid var(--border)', borderRadius: 6, padding: 6 }}>
      {letters.map((l) => (
        <span key={l} className={`chip selectable ${value.includes(l) ? 'selected' : ''}`} onClick={() => toggle(l)}>{l}</span>
      ))}
    </div>
  )
}

function parsePairList(text, letters) {
  const letterSet = new Set(letters)
  const pairs = []
  const errors = []
  const body = (text || '').trim()
  if (!body) return { pairs, errors: ['Nothing to import.'] }
  const groups = body.match(/\{[^{}]+\}/g)
  if (!groups) return { pairs, errors: ['No {letter1,letter2} groups found. Expected e.g. {{W[1],W[2]},{W[3],W[5]}}'] }
  for (const g of groups) {
    const items = g.slice(1, -1).split(',').map((s) => s.trim()).filter(Boolean)
    if (items.length !== 2) {
      errors.push(`${g} is not a pair`)
      continue
    }
    const bad = items.filter((x) => !letterSet.has(x))
    if (bad.length) {
      errors.push(`unknown letter(s): ${bad.join(', ')}`)
      continue
    }
    pairs.push([items[0], items[1]])
  }
  return { pairs, errors }
}

function PairEditor({ letters, value, onChange }) {
  const [a, setA] = useState('')
  const [b, setB] = useState('')
  const [showPaste, setShowPaste] = useState(false)
  const [pasteText, setPasteText] = useState('')
  const [pasteErrors, setPasteErrors] = useState([])
  const add = () => {
    if (!a || !b) return
    onChange([...value, [a, b]])
    setA(''); setB('')
  }
  const importPairs = () => {
    const { pairs, errors } = parsePairList(pasteText, letters)
    setPasteErrors(errors)
    if (pairs.length) {
      const seen = new Set(value.map((p) => `${p[0]}|${p[1]}`))
      const fresh = pairs.filter((p) => !seen.has(`${p[0]}|${p[1]}`))
      onChange([...value, ...fresh])
      if (!errors.length) {
        setShowPaste(false)
        setPasteText('')
      }
    }
  }
  return (
    <div>
      <div className="row">
        <select value={a} onChange={(e) => setA(e.target.value)}>
          <option value="">letter A…</option>
          {letters.map((l) => <option key={l} value={l}>{l}</option>)}
        </select>
        <select value={b} onChange={(e) => setB(e.target.value)}>
          <option value="">letter B…</option>
          {letters.map((l) => <option key={l} value={l}>{l}</option>)}
        </select>
        <button className="shrink" onClick={add}>Add pair</button>
      </div>
      <div style={{ marginTop: 6 }}>
        {value.map((p, i) => (
          <span key={i} className="chip">
            ({p[0]}, {p[1]})
            <button className="small danger" style={{ marginLeft: 6, border: 'none', background: 'none' }} onClick={() => onChange(value.filter((_, j) => j !== i))}>×</button>
          </span>
        ))}
      </div>
      <div style={{ marginTop: 8 }}>
        <button className="small" onClick={() => { setShowPaste(!showPaste); setPasteErrors([]) }}>
          {showPaste ? 'Hide list input' : 'Paste a pair list…'}
        </button>
      </div>
      {showPaste && (
        <div style={{ marginTop: 6 }}>
          <textarea
            rows={4} value={pasteText}
            onChange={(e) => setPasteText(e.target.value)}
            placeholder="{{W[1],W[2]},{W[3],W[5]},...}"
          />
          <div className="row" style={{ marginTop: 6 }}>
            <span className="muted" style={{ fontSize: 11 }}>Wolfram-style list; imported pairs are added to the ones above.</span>
            <button className="small primary shrink" onClick={importPairs}>Import</button>
          </div>
          {pasteErrors.map((e, i) => <p key={i} className="error-text" style={{ margin: '4px 0' }}>{e}</p>)}
        </div>
      )}
    </div>
  )
}

function AddPropertyModal({ alphabet, onClose, onAdded }) {
  const { project } = useProject()
  const [type, setType] = useState('')
  const [propName, setPropName] = useState('')
  const [letters, setLetters] = useState([])
  const [pairs, setPairs] = useState([])
  const [transName, setTransName] = useState('')
  const [transMap, setTransMap] = useState('')
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState(null)

  const submit = async () => {
    let params = {}
    if (type === 'first_entry' || type === 'last_entry') {
      if (!letters.length) return setError('Select at least one letter.')
      params = { letters }
    } else if (type === 'extended_steinmann') {
      if (!pairs.length) return setError('Add at least one non-adjacent pair.')
      params = { nonadjacent_pairs: pairs }
    } else if (type === 'cluster_adjacency') {
      if (!pairs.length) return setError('Add at least one adjacent pair.')
      params = { adjacent_pairs: pairs }
    } else if (type === 'transformation') {
      if (!transName.trim() || !transMap.trim()) return setError('Provide a name and a kinematic map.')
      params = { name: transName.trim(), map: transMap.trim() }
    }
    const name = type === 'transformation' ? transName.trim() : propName.trim()
    setBusy(true); setError(null)
    try {
      await api.addProperty(project.id, alphabet.id, { type, name, params })
      onAdded()
    } catch (e) { setError(e.message); setBusy(false) }
  }

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>Add property to “{alphabet.name}”</h3>
        <label>Property type</label>
        <select value={type} onChange={(e) => { setType(e.target.value); setError(null) }}>
          <option value="">— choose —</option>
          {Object.entries(PROP_LABEL).map(([k, v]) => <option key={k} value={k}>{v}</option>)}
        </select>
        {(type === 'first_entry' || type === 'last_entry') && (
          <>
            <label>{type === 'first_entry' ? 'Allowed first-entry letters' : 'Allowed last-entry letters'}</label>
            <LetterMultiSelect letters={alphabet.letters} value={letters} onChange={setLetters} />
          </>
        )}
        {type === 'extended_steinmann' && (
          <>
            <label>Non-adjacent letter pairs (Steinmann)</label>
            <PairEditor letters={alphabet.letters} value={pairs} onChange={setPairs} />
          </>
        )}
        {type === 'cluster_adjacency' && (
          <>
            <label>Adjacent letter pairs (cluster)</label>
            <PairEditor letters={alphabet.letters} value={pairs} onChange={setPairs} />
          </>
        )}
        {type === 'transformation' && (
          <>
            <label>Transformation name</label>
            <input value={transName} onChange={(e) => setTransName(e.target.value)} placeholder="e.g. Cyclic" />
            <label>Kinematic map (Wolfram rules)</label>
            <textarea rows={4} value={transMap} onChange={(e) => setTransMap(e.target.value)} placeholder="{u1->u2, u2->u3, u3->1+u2-v1-v2, v1->v2, v2->1+u3-u1-v2}" />
          </>
        )}
        {type && type !== 'transformation' && (
          <>
            <label>Property name (optional — needed if you add several {PROP_LABEL[type] || type} properties, e.g. for different physical objects)</label>
            <input value={propName} onChange={(e) => setPropName(e.target.value)} placeholder="e.g. e12FirstEntry, mhvFirstEntry" />
          </>
        )}
        {type === 'integrability' && (
          <p className="muted" style={{ fontSize: 12 }}>
            The integrability condition (dlogmat) follows from the alphabet and its parameterization — no extra input needed.
          </p>
        )}
        {error && <p className="error-text">{error}</p>}
        <div className="actions">
          <button onClick={onClose}>Cancel</button>
          <button className="primary" disabled={!type || busy} onClick={submit}>{busy ? 'Adding…' : 'Add property'}</button>
        </div>
      </div>
    </div>
  )
}

function NewAlphabetModal({ onClose, onAdded }) {
  const { project } = useProject()
  const [name, setName] = useState('')
  const [letters, setLetters] = useState('')
  const [variables, setVariables] = useState('')
  const [expressions, setExpressions] = useState('')
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState(null)

  const submit = async () => {
    const ls = letters.split(',').map((s) => s.trim()).filter(Boolean)
    if (!name.trim() || !ls.length) return setError('Name and a comma-separated letter list are required.')
    const exprs = expressions.split('\n').map((s) => s.trim()).filter(Boolean)
    if (exprs.length && exprs.length !== ls.length) return setError(`Got ${exprs.length} expressions for ${ls.length} letters. Provide one expression per line, or leave empty.`)
    setBusy(true); setError(null)
    try {
      await api.createAlphabet(project.id, {
        name: name.trim(),
        letters: ls,
        variables: variables.split(',').map((s) => s.trim()).filter(Boolean),
        expressions: exprs,
        roots: {},
      })
      onAdded()
    } catch (e) { setError(e.message); setBusy(false) }
  }

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>New alphabet</h3>
        <label>Name</label>
        <input value={name} onChange={(e) => setName(e.target.value)} placeholder="e.g. myAlphabet" />
        <label>Letters (comma-separated)</label>
        <input value={letters} onChange={(e) => setLetters(e.target.value)} placeholder="W[1], W[2], W[3]" />
        <label>Independent variables (comma-separated)</label>
        <input value={variables} onChange={(e) => setVariables(e.target.value)} placeholder="u1, u2, v1" />
        <label>Letter expressions (one per line, in letter order)</label>
        <textarea rows={6} value={expressions} onChange={(e) => setExpressions(e.target.value)} placeholder={'u1\nu2\n1-u1-u2'} />
        {error && <p className="error-text">{error}</p>}
        <div className="actions">
          <button onClick={onClose}>Cancel</button>
          <button className="primary" disabled={busy} onClick={submit}>{busy ? 'Creating…' : 'Create'}</button>
        </div>
      </div>
    </div>
  )
}

function PropertyRow({ alphabet, prop, onChanged }) {
  const { project, refreshProject } = useProject()
  const toast = useToast()
  const [busy, setBusy] = useState(false)

  const compute = async () => {
    setBusy(true)
    try {
      const { run_id } = await api.computeProperty(project.id, alphabet.id, prop.id)
      toast('Computation started…')
      const poll = setInterval(async () => {
        try {
          const run = await api.getRun(run_id)
          if (run.status !== 'queued' && run.status !== 'running') {
            clearInterval(poll)
            setBusy(false)
            await refreshProject()
            toast(run.status === 'done' ? 'Property computed.' : `Computation ${run.status}.`)
          }
        } catch { clearInterval(poll); setBusy(false) }
      }, 2000)
    } catch (e) { setBusy(false); toast(e.message) }
  }

  const remove = async () => {
    await api.deleteProperty(project.id, alphabet.id, prop.id)
    onChanged()
  }

  const display = (prop.name || '').trim() || (prop.type === 'transformation' ? prop.params?.name : '') || PROP_LABEL[prop.type] || prop.type
  const showType = display !== (PROP_LABEL[prop.type] || prop.type)
  return (
    <tr>
      <td><StatusDot status={busy ? 'computing' : prop.status} /> <strong>{display}</strong>{showType && <span className="muted"> ({PROP_LABEL[prop.type] || prop.type})</span>}{prop.precomputed ? ' (precomputed)' : ''}</td>
      <td className="mono">
        {prop.type === 'transformation' && (prop.params?.map || '').length > 60 ? `${prop.params.map.slice(0, 60)}…` : prop.params?.map}
        {(prop.type === 'first_entry' || prop.type === 'last_entry') && (prop.params?.letters || []).join(', ')}
        {prop.type === 'extended_steinmann' && `${(prop.params?.nonadjacent_pairs || []).length} pairs`}
        {prop.type === 'cluster_adjacency' && `${(prop.params?.adjacent_pairs || []).length} pairs`}
        {prop.summary?.dims && (
          <div className="muted" style={{ fontSize: 11 }}>dims {prop.summary.dims.join('×')}, nnz {prop.summary.nnz}</div>
        )}
      </td>
      <td className="mono">{prop.tensor_file || '—'}</td>
      <td style={{ textAlign: 'right' }}>
        {!prop.precomputed && prop.status !== 'ready' && (
          <button className="small primary" disabled={busy} onClick={compute}>{busy ? 'Computing…' : 'Compute'}</button>
        )}
        {' '}
        {!prop.precomputed && <button className="small danger" onClick={remove}>Delete</button>}
      </td>
    </tr>
  )
}

function AlphabetDetail({ alphabet, onBack }) {
  const { refreshProject } = useProject()
  const [showAdd, setShowAdd] = useState(false)
  const refresh = () => refreshProject()

  return (
    <div className="page">
      <button className="small" onClick={onBack}>← All alphabets</button>
      <h1 style={{ marginTop: 10 }}>{alphabet.name}</h1>
      <div className="card" style={{ marginBottom: 16 }}>
        <h2 style={{ marginTop: 0 }}>Letters ({alphabet.letters.length})</h2>
        <div style={{ maxHeight: 120, overflow: 'auto' }}>
          {alphabet.letters.map((l) => <span key={l} className="chip">{l}</span>)}
        </div>
        {alphabet.variables?.length > 0 && (
          <p className="muted" style={{ fontSize: 12 }}>Variables: {alphabet.variables.join(', ')}</p>
        )}
        {alphabet.expressions?.length > 0 && (
          <details>
            <summary className="muted" style={{ fontSize: 12, cursor: 'pointer' }}>Parameterization ({alphabet.expressions.length} expressions)</summary>
            <table style={{ marginTop: 8 }}>
              <tbody>
                {alphabet.letters.map((l, i) => (
                  <tr key={l}><td className="mono">{l}</td><td className="mono">{alphabet.expressions[i]}</td></tr>
                ))}
              </tbody>
            </table>
          </details>
        )}
        {alphabet.expr_loader && <p className="muted" style={{ fontSize: 12 }}>Parameterization loaded from the template alphabet file.</p>}
      </div>
      <div className="card">
        <div className="row" style={{ marginBottom: 8 }}>
          <h2 style={{ margin: 0, flex: 1 }}>Properties</h2>
          <button className="primary small shrink" onClick={() => setShowAdd(true)}>+ Add property</button>
        </div>
        <table>
          <thead>
            <tr><th>Property</th><th>Details</th><th>Tensor</th><th></th></tr>
          </thead>
          <tbody>
            {alphabet.properties.map((p) => (
              <PropertyRow key={p.id} alphabet={alphabet} prop={p} onChanged={refresh} />
            ))}
            {!alphabet.properties.length && (
              <tr><td colSpan={4} className="muted">No properties yet. Add one, or compute the integrability dlogmat.</td></tr>
            )}
          </tbody>
        </table>
      </div>
      {showAdd && <AddPropertyModal alphabet={alphabet} onClose={() => setShowAdd(false)} onAdded={() => { setShowAdd(false); refresh() }} />}
    </div>
  )
}

export default function Materials() {
  const { project, refreshProject, setCurrentId } = useProject()
  const [selected, setSelected] = useState(null)
  const [showNew, setShowNew] = useState(false)

  if (!project) return null
  const alphabet = selected ? project.alphabets.find((a) => a.id === selected) : null
  if (alphabet) return <AlphabetDetail alphabet={alphabet} onBack={() => setSelected(null)} />

  return (
    <div className="page">
      <div className="row" style={{ marginBottom: 16 }}>
        <h1 style={{ margin: 0, flex: 1 }}>Materials — alphabet building blocks</h1>
        <button className="primary shrink" onClick={() => setShowNew(true)}>+ New alphabet</button>
      </div>
      {!project.alphabets.length && (
        <div className="empty-state">
          <div className="big">No alphabets yet</div>
          <p>An alphabet is the basic building block of a symbol bootstrap problem. Create one, or start a new project from a template.</p>
        </div>
      )}
      <div className="grid">
        {project.alphabets.map((a) => (
          <div key={a.id} className="card clickable" onClick={() => setSelected(a.id)}>
            <h2 style={{ marginTop: 0 }}>{a.name}</h2>
            <div style={{ maxHeight: 66, overflow: 'hidden', marginBottom: 8 }}>
              {a.letters.slice(0, 24).map((l) => <span key={l} className="chip">{l}</span>)}
              {a.letters.length > 24 && <span className="chip">+{a.letters.length - 24} more</span>}
            </div>
            <div>
              {a.properties.map((p) => (
                <span key={p.id} className="badge" style={{ marginRight: 4, marginBottom: 4 }}>
                  <StatusDot status={p.status} /> {(p.name || '').trim() || (p.type === 'transformation' ? p.params?.name : '') || PROP_LABEL[p.type] || p.type}
                </span>
              ))}
            </div>
          </div>
        ))}
      </div>
      {showNew && <NewAlphabetModal onClose={() => setShowNew(false)} onAdded={() => { setShowNew(false); refreshProject() }} />}
    </div>
  )
}
