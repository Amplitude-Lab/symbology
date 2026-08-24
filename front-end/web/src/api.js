const BASE = '/api'

async function request(method, path, body) {
  const opts = { method, headers: {} }
  if (body !== undefined) {
    opts.headers['Content-Type'] = 'application/json'
    opts.body = JSON.stringify(body)
  }
  const res = await fetch(BASE + path, opts)
  const text = await res.text()
  let data = null
  try { data = text ? JSON.parse(text) : null } catch { data = text }
  if (!res.ok) {
    const detail = data && data.detail ? data.detail : res.statusText
    const err = new Error(typeof detail === 'string' ? detail : JSON.stringify(detail))
    err.status = res.status
    err.payload = data
    throw err
  }
  return data
}

export const api = {
  env: () => request('GET', '/env'),
  templates: () => request('GET', '/templates'),
  listProjects: () => request('GET', '/projects'),
  createProject: (name, template_id) => request('POST', '/projects', { name, template_id }),
  getProject: (pid) => request('GET', `/projects/${pid}`),
  deleteProject: (pid) => request('DELETE', `/projects/${pid}`),
  createAlphabet: (pid, a) => request('POST', `/projects/${pid}/alphabets`, a),
  updateAlphabet: (pid, aid, a) => request('PUT', `/projects/${pid}/alphabets/${aid}`, a),
  deleteAlphabet: (pid, aid) => request('DELETE', `/projects/${pid}/alphabets/${aid}`),
  addProperty: (pid, aid, p) => request('POST', `/projects/${pid}/alphabets/${aid}/properties`, p),
  deleteProperty: (pid, aid, propId) => request('DELETE', `/projects/${pid}/alphabets/${aid}/properties/${propId}`),
  computeProperty: (pid, aid, propId) => request('POST', `/projects/${pid}/alphabets/${aid}/properties/${propId}/compute`),
  summarizeProperty: (pid, aid, propId) => request('POST', `/projects/${pid}/alphabets/${aid}/properties/${propId}/summarize`),
  listTensors: (pid, dir) => request('GET', `/projects/${pid}/tensors?dir=${dir}`),
  tensorSummary: (pid, file) => request('GET', `/projects/${pid}/tensor_summary?file=${encodeURIComponent(file)}`),
  createFlow: (pid, name, extra) => request('POST', `/projects/${pid}/flows`, { name, ...extra }),
  duplicateFlow: (pid, fid) => request('POST', `/projects/${pid}/flows/${fid}/duplicate`, {}),
  updateFlow: (pid, fid, flow) => request('PUT', `/projects/${pid}/flows/${fid}`, flow),
  deleteFlow: (pid, fid) => request('DELETE', `/projects/${pid}/flows/${fid}`),
  compileFlow: (pid, fid) => request('POST', `/projects/${pid}/flows/${fid}/compile`, {}),
  flowOutputs: (pid) => request('GET', `/projects/${pid}/flow_outputs`),
  runFlow: (pid, fid) => request('POST', `/projects/${pid}/flows/${fid}/runs`, {}),
  listRuns: (pid) => request('GET', `/projects/${pid}/runs`),
  getRun: (rid) => request('GET', `/runs/${rid}`),
  cancelRun: (rid) => request('POST', `/runs/${rid}/cancel`),
  runEventsUrl: (rid) => `${BASE}/runs/${rid}/events`,
}
