export const PROP_KIND = {
  integrability: 'dlogmat',
  extended_steinmann: 'dlogmat',
  cluster_adjacency: 'dlogmat',
  first_entry: 'fec1',
  last_entry: 'lec1',
  transformation: 'matrix',
  precomputed_tensor: 'matrix',
  letter_symmetry: 'matrix',
}

export const PROP_LABEL = {
  integrability: 'Integrability',
  first_entry: 'First Entry',
  last_entry: 'Last Entry',
  extended_steinmann: 'Ext. Steinmann',
  cluster_adjacency: 'Cluster Adj.',
  transformation: 'Transformation',
  precomputed_tensor: 'Tensor file',
  letter_symmetry: 'Symmetry (rule)',
}

export function propLabel(prop) {
  const base = PROP_LABEL[prop.type] || prop.type
  const custom = (prop.name || '').trim() || (prop.type === 'transformation' ? prop.params?.name : '')
  return custom ? `${custom} · ${base}` : base
}

export const NODE_DEFS = {
  cb_in: { title: 'Block Input', color: '#e8f4ea', inputs: [], outputs: [{ id: 'out', kind: 'any', label: '' }] },
  cb_out: { title: 'Block Output', color: '#fdeeea', inputs: [{ id: 'in', kind: 'any', label: '' }], outputs: [] },
  reuse_output: {
    title: 'Reuse Output', color: '#e8f0fe',
    inputs: [],
    outputs: [{ id: 'out', kind: 'tensor', label: 'reused tensor' }],
  },
  alphabet: { title: 'Alphabet', color: '#eef0fe', inputs: [] },
  merge_conditions: {
    title: 'Merge Conditions', color: '#f3e8ff',
    inputs: [
      { id: 'in_0', kind: 'dlogmat', label: 'dlog 1' },
      { id: 'in_1', kind: 'dlogmat', label: 'dlog 2' },
      { id: 'in_2', kind: 'dlogmat', label: 'dlog 3' },
      { id: 'in_3', kind: 'dlogmat', label: 'dlog 4' },
    ],
    outputs: [{ id: 'out', kind: 'dlogmat', label: 'merged' }],
  },
  extend: {
    title: 'Extend', color: '#ccfbf1',
    inputs: [
      { id: 'condition', kind: 'dlogmat', label: 'condition' },
      { id: 'fec', kind: 'fec', label: 'FEC in' },
      { id: 'lec', kind: 'lec', label: 'LEC in' },
    ],
    outputs: [
      { id: 'fec', kind: 'fec', label: 'FEC out' },
      { id: 'lec', kind: 'lec', label: 'LEC out' },
    ],
  },
  sew: {
    title: 'Sew', color: '#ffedd5',
    inputs: [
      { id: 'condition', kind: 'dlogmat', label: 'condition' },
      { id: 'fec', kind: 'fec', label: 'FEC' },
      { id: 'lec', kind: 'lec', label: 'LEC' },
    ],
    outputs: [{ id: 'sew', kind: 'sew', label: 'SEW' }],
  },
  project: {
    title: 'Project', color: '#e2e8f0',
    inputs: [
      { id: 'tensor', kind: 'seed', label: 'tensor' },
      { id: 'rep', kind: 'matrix', label: 'rep (axes 2,3)' },
      { id: 'map', kind: 'matrix', label: 'map (axis 1)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'projected' }],
  },
  solve_symmetry: {
    title: 'Symmetry Solve', color: '#e2e8f0',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'matrix', kind: 'matrix', label: 'sym M (b)' },
      { id: 'matrix2', kind: 'matrix', label: 'sym M (c, opt)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'invariant part' }],
  },
  symderive: {
    title: 'Symmetry Derive', color: '#e2e8f0',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'matrix', kind: 'matrix', label: 'sym M (b)' },
      { id: 'matrix2', kind: 'matrix', label: 'sym M (c, opt)' },
    ],
    outputs: [{ id: 'out', kind: 'matrix', label: 'R (a,a)' }],
  },
  solve_collinear: {
    title: 'Solve Collinear', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'solution', kind: 'solution', label: 'solution' }],
  },
  projection_chain: {
    title: 'Projection Chain', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'chain seed (opt)' }],
    outputs: [{ id: 'out', kind: 'basis', label: 'basis / projection' }],
  },
  symmetry_invariant: {
    title: 'Symmetry Invariant', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'chain seed (opt)' }],
    outputs: [{ id: 'out', kind: 'basis', label: 'invariant basis' }],
  },
  compute_rhs: {
    title: 'Compute RHS', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'boundary', kind: 'boundary', label: 'boundary' }],
  },
  add_tensors: {
    title: 'Add Tensors', color: '#dcfce7',
    inputs: [
      { id: 'in_0', kind: 'tensor', label: 'A' },
      { id: 'in_1', kind: 'tensor', label: 'B' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'Σ wᵢ·Aᵢ' }],
  },
  ternary_contract: {
    title: 'Ternary Contract', color: '#fef3c7',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'trans1', kind: 'matrix', label: 'trans1 (2nd axis)' },
      { id: 'trans2', kind: 'matrix', label: 'trans2 (last axis)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'transformed' }],
  },
  apply_symmetry: {
    title: 'Apply Projection', color: '#fef3c7',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'sym', kind: 'matrix', label: 'sym (auto chain)' },
      { id: 'trans1', kind: 'matrix', label: 'trans1 (2nd axis)' },
      { id: 'trans2', kind: 'matrix', label: 'trans2 (last axis)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'transformed' }],
  },
  matrix_power: {
    title: 'Matrix Power', color: '#fef3c7',
    inputs: [{ id: 'matrix', kind: 'matrix', label: 'matrix M' }],
    outputs: [{ id: 'out', kind: 'matrix', label: 'M^n' }],
  },
  tensor_join: {
    title: 'Join Tensors', color: '#dcfce7',
    inputs: [
      { id: 'a', kind: 'tensor', label: 'A' },
      { id: 'b', kind: 'tensor', label: 'B' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'joined' }],
  },
  tensor_dot: {
    title: 'Tensor Dot', color: '#dcfce7',
    inputs: [
      { id: 'a', kind: 'tensor', label: 'A' },
      { id: 'b', kind: 'tensor', label: 'B' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'dot' }],
  },
  impose_integrability: {
    title: 'Solve Integrability', color: '#ede9fe',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'dlogmat', kind: 'dlogmat', label: 'integrability dlog' },
    ],
    outputs: [{ id: 'out', kind: 'matrix', label: 'solution basis' }],
  },
  integrability_condition: {
    title: 'Integrability Condition', color: '#ddd6fe',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'dlogmat', kind: 'dlogmat', label: 'integrability dlog' },
    ],
    outputs: [{ id: 'out', kind: 'matrix', label: 'condition matrix' }],
    hint: 'M[(a), b·d] = Σ S[(a),b,i,j]·D[j,i,d]; first entry kept as rows, middle entries × d flattened to columns',
  },
  solve_conditions: {
    title: 'Solve Conditions', color: '#ddd6fe',
    inputs: [
      { id: 'cond', kind: 'any', label: 'condition matrix' },
    ],
    outputs: [{ id: 'out', kind: 'matrix', label: 'solution basis' }],
  },
  assemble: {
    title: 'Assemble Solution Space', color: '#cffafe',
    inputs: [
      { id: 'in_0', kind: 'tensor', label: 'elem 1' },
      { id: 'in_1', kind: 'tensor', label: 'elem 2' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'aligned space' }],
  },
}

export const PALETTE_SECTIONS = [
  {
    title: 'Building blocks',
    items: [
      { type: 'alphabet', label: 'Alphabet', sub: 'letters, properties, tensors' },
    ],
  },
  {
    title: 'Basic tensor operations',
    items: [
      { type: 'merge_conditions', label: 'Merge Conditions', sub: 'combine dlogmats' },
      { type: 'add_tensors', label: 'Add Tensors', sub: 'Σ wᵢ·Aᵢ' },
      { type: 'ternary_contract', label: 'Ternary Contract', sub: 'T·M1·M2 contraction' },
      { type: 'apply_symmetry', label: 'Apply Projection', sub: 'auto-derive chain basis maps from sym' },
      { type: 'matrix_power', label: 'Matrix Power', sub: 'M^n' },
      { type: 'tensor_join', label: 'Join Tensors', sub: 'join along an axis' },
      { type: 'tensor_dot', label: 'Tensor Dot', sub: 'contract one axis of A with one of B' },
      { type: 'impose_integrability', label: 'Solve Integrability', sub: 'contract with dlog & solve' },
      { type: 'integrability_condition', label: 'Integrability Condition', sub: 'contract with dlog → conditions' },
      { type: 'solve_conditions', label: 'Solve Conditions', sub: 'kernel of condition matrix' },
      { type: 'assemble', label: 'Assemble Solution Space', sub: 'align & project element spaces' },
    ],
  },
  {
    title: 'Composite operations',
    items: [
      { type: 'extend', label: 'Extend', sub: 'FEC/LEC weight +1' },
      { type: 'sew', label: 'Sew', sub: 'FEC + LEC → SEW' },
      { type: 'project', label: 'Project', sub: 'apply rep + contract map' },
      { type: 'solve_symmetry', label: 'Symmetry Solve', sub: 'invariant part of tensor' },
      { type: 'symderive', label: 'Symmetry Derive', sub: 'derive R (a,a) transformation' },
      { type: 'solve_collinear', label: 'Solve Collinear', sub: 'collinear constraints' },
      { type: 'projection_chain', label: 'Projection Chain', sub: 'bootstrap --project (basis chain)' },
      { type: 'symmetry_invariant', label: 'Symmetry Invariant', sub: 'bootstrap --solve-symmetry' },
      { type: 'compute_rhs', label: 'Compute RHS', sub: 'collinear RHS / boundary' },
    ],
  },
]

const ACCEPTS = {
  dlogmat: ['dlogmat'],
  fec: ['fec1', 'fec'],
  lec: ['lec1', 'lec'],
  seed: ['seed', 'fec1', 'fec', 'lec1', 'lec', 'sew'],
  tensor: ['dlogmat', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'matrix', 'basis', 'solution', 'boundary', 'tensor'],
  matrix: ['matrix'],
  any: ['dlogmat', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'matrix', 'basis', 'solution', 'boundary', 'tensor'],
}

export function kindsCompatible(sourceKind, targetKind) {
  if (targetKind === 'any') return true
  if (sourceKind === 'any') return true // abstract IO inside a custom-block canvas
  const acc = ACCEPTS[targetKind]
  return acc ? acc.includes(sourceKind) : false
}

export function alphabetOutSlots(alphabet, selectedProps) {
  const rows = []
  if (!alphabet) return rows
  for (const p of (alphabet.properties || []).filter((x) => (selectedProps || []).includes(x.id))) {
    rows.push({ key: p.id, handle: `prop_${p.id}`, propId: p.id, proj: false, kind: PROP_KIND[p.type] })
    if (p.type === 'first_entry' || p.type === 'last_entry') {
      rows.push({ key: `${p.id}_proj`, handle: `proj_${p.id}`, propId: p.id, proj: true, kind: 'matrix' })
    }
  }
  return rows
}

export function sourceKindFor(node, handleId, project) {
  if (!node) return null
  if (node.type === 'cb_in') return 'any'
  if (node.type === 'reuse_output') {
    const k = node.data?.kind
    return k || 'tensor'
  }
  if (node.type === 'customblock') {
    const cb = customBlockDef(node, project)
    const out = cb?.outputs?.find((o) => o.id === handleId)
    return out ? out.kind : null
  }
  if (node.type === 'alphabet') {
    const alphabet = project?.alphabets.find((a) => a.id === node.data?.alphabet_id)
    if ((handleId || '').startsWith('proj_')) {
      const propId = handleId.replace(/^proj_/, '')
      const prop = alphabet?.properties.find((p) => p.id === propId)
      return prop && (prop.type === 'first_entry' || prop.type === 'last_entry') ? 'matrix' : null
    }
    const propId = (handleId || '').replace(/^prop_/, '')
    const prop = alphabet?.properties.find((p) => p.id === propId)
    return prop ? PROP_KIND[prop.type] : null
  }
  if (node.type === 'assemble' && /^out_\d+$/.test(handleId || '')) {
    return 'tensor'
  }
  const def = NODE_DEFS[node.type]
  const out = def?.outputs?.find((o) => o.id === handleId)
  return out ? out.kind : null
}

export function targetKindFor(node, handleId) {
  if (!node) return null
  if (node.type === 'cb_out') return 'any'
  if (node.type === 'customblock') {
    const cb = NODE_DEFS[`cb_${node.data?.block}`]
    const inp = cb?.inputs?.find((i) => i.id === handleId)
    return inp ? inp.kind : null
  }
  const def = NODE_DEFS[node.type]
  if ((node.type === 'assemble' || node.type === 'add_tensors') && (handleId || '').startsWith('in_')) return 'tensor'
  const inp = def.inputs.find((i) => i.id === handleId)
  return inp ? inp.kind : null
}

// ---- custom composite blocks ----
// A custom block is a flow flagged `custom_block: true` whose graph contains
// abstract cb_in / cb_out nodes. registerCustomBlocks() mirrors each such flow
// into NODE_DEFS under `cb_<flowId>` so instances render and type-check like
// any built-in op.

const CB_KINDS = ['tensor', 'matrix', 'dlogmat', 'fec', 'lec', 'seed']

export function cbPortKind(v) {
  return CB_KINDS.includes(v) ? v : 'tensor'
}

export function registerCustomBlocks(project) {
  const cbs = (project?.flows || []).filter((f) => f.custom_block)
  for (const f of cbs) delete NODE_DEFS[`cb_${f.id}`]
  for (const f of cbs) {
    const nodes = f.graph?.nodes || []
    const inputs = nodes
      .filter((n) => n.type === 'cb_in')
      .sort((a, b) => (a.data?.order ?? 0) - (b.data?.order ?? 0))
      .map((n) => ({ id: n.id, kind: cbPortKind(n.data?.kind), label: n.data?.name || n.id }))
    const outputs = nodes
      .filter((n) => n.type === 'cb_out')
      .sort((a, b) => (a.data?.order ?? 0) - (b.data?.order ?? 0))
      .map((n) => ({ id: n.id, kind: cbPortKind(n.data?.kind), label: n.data?.name || n.id }))
    NODE_DEFS[`cb_${f.id}`] = {
      title: f.name || 'Custom block',
      color: '#efe6f7',
      custom: true,
      blockId: f.id,
      inputs,
      outputs,
    }
  }
}

export function customBlockDef(node, project) {
  if (!node || node.type !== 'customblock') return null
  return NODE_DEFS[`cb_${node.data?.block}`] || null
}

export function customBlockPalette(project) {
  return (project?.flows || []).filter((f) => f.custom_block).map((f) => ({
    type: 'customblock', block: f.id, label: f.name || 'Custom block',
    sub: `${(f.graph?.nodes || []).filter((n) => n.type === 'cb_out').length} output(s)`,
  }))
}
