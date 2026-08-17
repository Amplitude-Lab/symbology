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
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'basis', kind: 'basis', label: 'basis' }],
  },
  solve_symmetry: {
    title: 'Solve Symmetry', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'solution', kind: 'solution', label: 'solution' }],
  },
  solve_collinear: {
    title: 'Solve Collinear', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'solution', kind: 'solution', label: 'solution' }],
  },
  compute_rhs: {
    title: 'Compute RHS', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'seed' }],
    outputs: [{ id: 'boundary', kind: 'boundary', label: 'boundary' }],
  },
  add_tensors: {
    title: 'Add Tensors', color: '#dcfce7',
    inputs: [
      { id: 'a', kind: 'tensor', label: 'A' },
      { id: 'b', kind: 'tensor', label: 'B' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'wA·A + wB·B' }],
  },
  apply_symmetry: {
    title: 'Apply Symmetry', color: '#fef3c7',
    inputs: [
      { id: 'tensor', kind: 'r3tensor', label: 'tensor' },
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
}

export const PALETTE = [
  { type: 'alphabet', label: 'Alphabet', sub: 'building block' },
  { type: 'merge_conditions', label: 'Merge Conditions', sub: 'join dlogmats' },
  { type: 'extend', label: 'Extend', sub: 'grow weight by 1' },
  { type: 'sew', label: 'Sew', sub: 'FEC + LEC → SEW' },
  { type: 'project', label: 'Project', sub: 'symmetry / collinear' },
  { type: 'solve_symmetry', label: 'Solve Symmetry', sub: 'constrain by symmetry' },
  { type: 'solve_collinear', label: 'Solve Collinear', sub: 'collinear constraints' },
  { type: 'compute_rhs', label: 'Compute RHS', sub: 'boundary terms' },
  { type: 'add_tensors', label: 'Add Tensors', sub: 'weighted sum wA·A + wB·B' },
  { type: 'apply_symmetry', label: 'Apply Symmetry', sub: 'contract last two entries' },
  { type: 'matrix_power', label: 'Matrix Power', sub: 'M^n group elements' },
  { type: 'tensor_join', label: 'Join Tensors', sub: 'concat along an axis' },
]

const ACCEPTS = {
  dlogmat: ['dlogmat'],
  fec: ['fec1', 'fec'],
  lec: ['lec1', 'lec'],
  seed: ['fec1', 'fec', 'lec1', 'sew'],
  tensor: ['dlogmat', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'matrix', 'basis', 'solution', 'boundary'],
  r3tensor: ['fec1', 'fec', 'lec1', 'lec', 'sew'],
  matrix: ['matrix'],
}

export function kindsCompatible(sourceKind, targetKind) {
  const acc = ACCEPTS[targetKind]
  return acc ? acc.includes(sourceKind) : false
}

export function sourceKindFor(node, handleId, project) {
  if (!node) return null
  if (node.type === 'alphabet') {
    const propId = (handleId || '').replace(/^prop_/, '')
    const alphabet = project?.alphabets.find((a) => a.id === node.data?.alphabet_id)
    const prop = alphabet?.properties.find((p) => p.id === propId)
    return prop ? PROP_KIND[prop.type] : null
  }
  const def = NODE_DEFS[node.type]
  const out = def?.outputs?.find((o) => o.id === handleId)
  return out ? out.kind : null
}

export function targetKindFor(node, handleId) {
  const def = node ? NODE_DEFS[node.type] : null
  if (!def) return null
  const inp = def.inputs.find((i) => i.id === handleId)
  return inp ? inp.kind : null
}
