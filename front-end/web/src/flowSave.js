import { api } from './api'

const sessions = new Map()
const keyFor = (pid, fid) => `symbology.draft.${pid}.${fid}`

export function openFlowSave(pid, flow) {
  const key = keyFor(pid, flow.id)
  if (sessions.has(key)) {
    const current = sessions.get(key)
    if (current.isSaved() && (flow.revision ?? 0) >= current.revision) {
      current.revision = flow.revision ?? 0
      current.saved = JSON.stringify({ name: flow.name, graph: flow.graph })
      current.draft = null
    }
    return current
  }
  let recovered = null
  try { recovered = JSON.parse(sessionStorage.getItem(key) || 'null') } catch { /* unavailable */ }
  const session = {
    revision: recovered?.revision ?? flow.revision ?? 0,
    draft: recovered?.draft || null,
    saved: JSON.stringify({ name: flow.name, graph: flow.graph }),
    chain: Promise.resolve(),
    pending: false,
    remember(draft) {
      this.draft = draft
      // A tab owns its recovery draft; another tab cannot overwrite it.
      sessionStorage.setItem(key, JSON.stringify({ revision: this.revision, draft }))
    },
    flush() {
      const work = this.chain.catch(() => {}).then(async () => {
        this.pending = true
        try {
          while (this.draft && JSON.stringify(this.draft) !== this.saved) {
            const draft = this.draft
            const result = await api.updateFlow(pid, flow.id, { ...draft, revision: this.revision })
            this.revision = result.revision
            this.saved = JSON.stringify(draft)
            try {
              if (JSON.stringify(this.draft) === this.saved) sessionStorage.removeItem(key)
              else this.remember(this.draft)
            } catch { /* the acknowledged server save is still valid */ }
          }
        } finally { this.pending = false }
      })
      this.chain = work
      return work
    },
    download() {
      const url = URL.createObjectURL(new Blob([JSON.stringify(this.draft, null, 2)], { type: 'application/json' }))
      const link = document.createElement('a')
      link.href = url; link.download = `${flow.id}-draft.json`; link.click()
      setTimeout(() => URL.revokeObjectURL(url), 1000)
    },
    discard() { this.draft = null; sessionStorage.removeItem(key); sessions.delete(key) },
    isSaved() { return !this.pending && (!this.draft || JSON.stringify(this.draft) === this.saved) },
  }
  sessions.set(key, session)
  return session
}
