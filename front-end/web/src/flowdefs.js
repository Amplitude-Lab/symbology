export const PROP_KIND = {
  integrability: 'dlogmat',
  extended_steinmann: 'dlogmat',
  cluster_adjacency: 'dlogmat',
  first_entry: 'fec1',
  last_entry: 'lec1',
  transformation: 'matrix',
  precomputed_tensor: 'matrix',
  letter_symmetry: 'matrix',
  sparse_expression: 'matrix',
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
  sparse_expression: 'Symbol Tensor',
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
    inputs: [
      { id: 'seed', kind: 'seed_or_tensor', label: 'seed (target or custom tensor)' },
      { id: 'rhs', kind: 'boundary', label: 'rhs (boundary, opt)' },
      { id: 'cond', kind: 'matrix', label: 'cond [M|r] (opt)' },
    ],
    outputs: [
      { id: 'solution', kind: 'solution', label: 'solution' },
      { id: 'conditions', kind: 'matrix', label: 'conditions [M|r]' },
    ],
  },
  // solve_collinear port model: the node always shows the two pair-1 ports
  // (seed / rhs). Each entry of data.pairs beyond the first adds one
  // in_seed_N / in_rhs_N pair of ports (N = pair index). The cond port is
  // opt-in via data.cond_enabled (a wired cond edge also keeps it visible).
  // Legacy graphs stored rhs / letter_projection / pairs_config —
  // normalizeCollinearData() upgrades them on load.
  projection_chain: {
    title: 'Projection Chain', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'chain seed (opt)' }],
    // Collinear runs additionally expose one dynamic output per expansion
    // basis the same --project run materializes (basis_w{w} / basis_last_w{w})
    // — see projectionChainBasisOutputs(); OpNode renders them from data.target.
    outputs: [{ id: 'out', kind: 'basis', label: 'basis / projection' }],
  },
  symmetry_invariant: {
    title: 'Symmetry Invariant', color: '#e2e8f0',
    inputs: [{ id: 'seed', kind: 'seed', label: 'chain seed (opt)' }],
    outputs: [{ id: 'out', kind: 'basis', label: 'invariant basis' }],
  },
  compute_rhs: {
    title: 'Compute RHS', color: '#e2e8f0',
    inputs: [
      { id: 'e1', kind: 'tensor', label: 'E1 seed' },
      { id: 'p1', kind: 'tensor', label: 'P1 seed (per-mode)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'rhs tensor' }],
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
      { id: 'trans1', kind: 'matrix_or_tensor', label: 'trans1 (2nd axis)' },
      { id: 'trans2', kind: 'matrix_or_tensor', label: 'trans2 (last axis)' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'transformed' }],
  },
  apply_symmetry: {
    title: 'Apply Projection', color: '#fef3c7',
    inputs: [
      { id: 'tensor', kind: 'tensor', label: 'tensor' },
      { id: 'sym', kind: 'matrix', label: 'sym (auto chain)' },
      { id: 'trans1', kind: 'matrix_or_tensor', label: 'trans1 (2nd axis)' },
      { id: 'trans2', kind: 'matrix_or_tensor', label: 'trans2 (last axis)' },
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
  shuffle_product: {
    title: 'Shuffle Product', color: '#dcfce7',
    inputs: [
      { id: 'a', kind: 'tensor', label: 'A' },
      { id: 'b', kind: 'tensor', label: 'B' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'w·(A ⊗ B)' }],
  },
  squeeze_tensor: {
    title: 'Squeeze Tensor', color: '#dcfce7',
    inputs: [{ id: 'in', kind: 'tensor', label: 'T' }],
    outputs: [{ id: 'out', kind: 'tensor', label: 'squeezed' }],
  },
  expand_tensor: {
    title: 'Expand Tensor', color: '#dcfce7',
    inputs: [
      { id: 'in_0', kind: 'tensor', label: 'T (FEC,letter)' },
      { id: 'in_1', kind: 'tensor', label: 'basis 1' },
    ],
    outputs: [{ id: 'out', kind: 'tensor', label: 'expanded' }],
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

// ---- geometry (shared single source of truth) ----
// port-row geometry constants and the TRUE rendered port lists for every
// node type (dynamic in_N slots, collinear pairs, chain basis outputs,
// alphabet proj slots, custom-block port lists). FlowEditor renders and
// routes wires from these; scripts/ports-dump.mjs + wire-audit.py and the
// SSR audit consume them directly, so the audits can never silently
// diverge from the renderer again (the alphabet def.outputs crash class).

export const NODE_W = 190
export const ROW0 = 29
export const ROW_STEP = 15
export const SUB_EXTRA = 15

export const SUBTITLE_TYPES = new Set([
  'extend', 'project', 'solve_symmetry', 'symderive', 'sew', 'solve_collinear',
  'projection_chain', 'symmetry_invariant', 'compute_rhs', 'add_tensors',
  'ternary_contract', 'apply_symmetry', 'matrix_power', 'tensor_join',
  'tensor_dot', 'squeeze_tensor', 'shuffle_product', 'expand_tensor',
  'impose_integrability', 'integrability_condition', 'solve_conditions',
])

// Compute RHS object types. Mirrored in server/app/compile.py::RHS_MODES —
// keep the two registries in sync. Adding a new object type = one entry here
// (palette/Inspector select + subtitle) + one entry there (the generation
// rule) + optional Inspector help text in FlowEditor.jsx.
export const RHS_MODES = {
  mhv_boundary: {
    label: 'MHV boundary (boundary_L)',
    short: 'MHV',
    requires: ['e1'], forbids: ['p1'],
    targetPlaceholder: 'boundary_2L',
  },
  e47me67_te: {
    label: 'NMHV E47mE67 (tE_L)',
    short: 'E47mE67',
    requires: ['e1', 'p1'], forbids: [],
    targetPlaceholder: 'tE2',
  },
}

// TRUE rendered port lists per node — the single source of truth mirroring
// what OpNode / AlphabetNode / CustomBlockNode / AssembleNode / CBIONode
// actually draw (dynamic in_N slots, collinear pairs, chain basis outputs,
// alphabet proj slots, custom-block port lists). estHeight and wire routing
// both derive from it; never read NODE_DEFS[t].outputs directly for geometry
// (alphabet has no outputs key; dynamic types exceed the static lists).
export function portsOf(node, project, allEdges) {
  const t = node.type
  const d = node.data || {}
  if (t === 'alphabet') {
    const alphabet = project?.alphabets?.find((a) => a.id === d.alphabet_id)
    const slots = alphabetOutSlots(alphabet, d.selected_properties || [])
    return { inputs: [], outputs: slots.map((s) => ({ id: s.handle })) }
  }
  if (t === 'customblock') {
    const def = NODE_DEFS[`cb_${d.block}`]
    return { inputs: def?.inputs || [], outputs: def?.outputs || [] }
  }
  if (t === 'assemble') {
    let maxConnected = -1
    for (const e of allEdges) {
      if (e.target === node.id && (e.targetHandle || '').startsWith('in_')) {
        const idx = parseInt((e.targetHandle || '').slice(3).split('@')[0], 10)
        if (!Number.isNaN(idx)) maxConnected = Math.max(maxConnected, idx)
      }
    }
    const def = NODE_DEFS.assemble
    const slots = Math.max(def.inputs.length, maxConnected + 2)
    return {
      inputs: Array.from({ length: slots }, (_, i) => ({ id: `in_${i}` })),
      outputs: (d.outputs || []).map((o, i) => ({ id: `out_${i}` })),
    }
  }
  if (t === 'cb_in') return { inputs: [], outputs: [{ id: 'out' }] }
  if (t === 'cb_out') return { inputs: [{ id: 'in' }], outputs: [] }
  if (t === 'reuse_output') return { inputs: [], outputs: [{ id: 'out' }] }
  if (t === 'groupBox') return { inputs: [], outputs: [] }
  const def = NODE_DEFS[t]
  if (!def) return { inputs: [], outputs: [] }
  let inputs = def.inputs || []
  let outputs = def.outputs || []
  if (t === 'add_tensors' || t === 'expand_tensor') {
    let maxConnected = 1
    for (const e of allEdges) {
      if (e.target === node.id && (e.targetHandle || '').startsWith('in_')) {
        const idx = parseInt((e.targetHandle || '').slice(3).split('@')[0], 10)
        if (!Number.isNaN(idx)) maxConnected = Math.max(maxConnected, idx)
      }
    }
    inputs = Array.from({ length: Math.max(inputs.length, maxConnected + 2) }, (_, i) => ({ id: `in_${i}` }))
  }
  if (t === 'solve_collinear') {
    let maxPair = 0
    let hasCondEdge = false
    for (const e of allEdges) {
      if (e.target !== node.id) continue
      const m = (e.targetHandle || '').match(/^in_(?:seed|rhs)_(\d+)$/)
      if (m) maxPair = Math.max(maxPair, parseInt(m[1], 10))
      if ((e.targetHandle || '') === 'cond') hasCondEdge = true
    }
    const nPairs = Math.max(maxPair, normalizeCollinearPairs(d).length - 1)
    inputs = inputs.filter((i) => i.id !== 'cond')
    if (d.cond_enabled || hasCondEdge) inputs = [...inputs, { id: 'cond' }]
    const extras = []
    for (let p = 1; p <= nPairs; p++) extras.push({ id: `in_seed_${p}` }, { id: `in_rhs_${p}` })
    inputs = [...inputs, ...extras]
  }
  if (t === 'projection_chain') {
    const basis = projectionChainBasisOutputs(d)
    for (const e of allEdges) {
      if (e.source === node.id && /^(basis_w\d+|basis_last_w\d+)$/.test(e.sourceHandle || '')) {
        if (!basis.some((b) => b.id === e.sourceHandle)) basis.push({ id: e.sourceHandle })
      }
    }
    outputs = [...outputs, ...basis]
  }
  return { inputs, outputs }
}

// Y of port row 0 inside a node card (matches the Handle top styles).
export function portBaseY(n) {
  if (n.type === 'alphabet') return 31
  if (n.type === 'customblock' || n.type === 'cb_in' || n.type === 'cb_out' || n.type === 'reuse_output') return ROW0
  if (n.type === 'assemble') return ROW0 + SUB_EXTRA
  return ROW0 + (SUBTITLE_TYPES.has(n.type) ? SUB_EXTRA : 0)
}

// Estimated rendered node height (mirror of scripts/layout-flow.py's model;
// used only for wire detour routing — small errors are absorbed by margins).
export function estHeight(n, project, allEdges) {
  if (n.type === 'groupBox' || n.type === 'cb_in' || n.type === 'cb_out' || n.type === 'reuse_output') return 60
  if (n.type === 'alphabet') {
    const { outputs } = portsOf(n, project, allEdges)
    return 62 + (Math.max(outputs.length, 1) - 1) * 15
  }
  const { inputs, outputs } = portsOf(n, project, allEdges)
  const sub = (SUBTITLE_TYPES.has(n.type) || n.type === 'assemble') ? SUB_EXTRA : 0
  return 98 + (Math.max(inputs.length, outputs.length, 1) - 1) * 15 + sub
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
      { type: 'squeeze_tensor', label: 'Squeeze Tensor', sub: 'drop all size-1 axes' },
      { type: 'shuffle_product', label: 'Shuffle Product', sub: 'w·(A ⊗ B) word shuffle' },
      { type: 'expand_tensor', label: 'Expand Tensor', sub: 'contract FEC axis with basis chain' },
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
      { type: 'solve_collinear', label: 'Solve Collinear', sub: 'non-homogeneous constraints' },
      { type: 'projection_chain', label: 'Projection Chain', sub: 'bootstrap --project (basis chain)' },
      { type: 'symmetry_invariant', label: 'Symmetry Invariant', sub: 'bootstrap --solve-symmetry' },
      { type: 'compute_rhs', label: 'Compute RHS', sub: 'one block: boundary_L / tE_L (shuffle recursion)' },
    ],
  },
]

const ACCEPTS = {
  dlogmat: ['dlogmat'],
  fec: ['fec1', 'fec'],
  lec: ['lec1', 'lec'],
  seed: ['seed', 'fec1', 'fec', 'lec1', 'lec', 'sew'],
  // Solve Collinear seed: a named target (SEW/FEC chain head) OR any tensor
  // produced by other flows (add_tensors / apply_symmetry / custom blocks).
  seed_or_tensor: ['seed', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'basis', 'solution', 'boundary', 'tensor', 'matrix'],
  boundary: ['boundary', 'tensor', 'basis', 'solution', 'matrix'],
  tensor: ['dlogmat', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'matrix', 'basis', 'solution', 'boundary', 'tensor'],
  matrix: ['matrix'],
  matrix_or_tensor: ['matrix', 'tensor', 'basis', 'solution', 'boundary', 'fec1', 'fec', 'lec1', 'lec', 'sew'],
  xtrans: ['matrix', 'fec1', 'fec', 'lec1', 'lec', 'sew', 'basis', 'solution', 'boundary', 'tensor'],
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

// ---- projection_chain dynamic basis outputs ----
// A collinear `bootstrap --project` run materializes one expansion basis per
// chain weight next to the target basis (first_w{w}_basis / last_w{w}_basis).
// The compile server exposes each as an extra output handle (basis_w{w} /
// basis_last_w{w}); this mirrors that contract for the UI so OpNode can
// render the ports and sourceKindFor can type-check wires. Must stay in sync
// with _compile_projection_chain in server/app/compile.py.

export function projectionChainBasisOutputs(data) {
  const rows = []
  if ((data?.symmetry || 'collinear') !== 'collinear') return rows
  const t = String(data?.target || '').trim()
  const m = t.match(/^(SEW|FEC|LEC)_(\d+)(?:p(\d+))?$/i)
  if (!m) return rows
  const kind = m[1].toUpperCase()
  const fw = kind === 'LEC' ? 1 : parseInt(m[2], 10)
  const lw = kind === 'FEC' ? 1 : kind === 'LEC' ? parseInt(m[2], 10) : parseInt(m[3], 10)
  if (kind === 'LEC') {
    for (let w = 2; w <= lw; w++) rows.push({ id: `basis_w${w}`, kind: 'basis', label: `w${w} basis (last chain)` })
  } else {
    for (let w = 2; w <= fw; w++) rows.push({ id: `basis_w${w}`, kind: 'basis', label: `w${w} basis` })
    for (let w = 2; w <= lw; w++) rows.push({ id: `basis_last_w${w}`, kind: 'basis', label: `w${w} basis (last chain)` })
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
  if (node.type === 'projection_chain' && /^(basis_w\d+|basis_last_w\d+)$/.test(handleId || '')) {
    return 'basis'
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
  if ((node.type === 'assemble' || node.type === 'add_tensors' || node.type === 'expand_tensor') && (handleId || '').startsWith('in_')) return 'tensor'
  if (node.type === 'solve_collinear') {
    // Dynamic extra pairs: in_seed_N / in_rhs_N (N >= 1). Pair 0 uses the
    // fixed seed/rhs ports (backwards compatible).
    if ((handleId || '').startsWith('in_seed_')) return 'seed_or_tensor'
    if ((handleId || '').startsWith('in_rhs_')) return 'boundary'
  }
  const inp = def.inputs.find((i) => i.id === handleId)
  return inp ? inp.kind : null
}

// ---- solve_collinear pair normalization ----
// Current shape: data.pairs = [{ projection, rhs }] with pair 1 at index 0
// (the fixed seed/rhs ports); entries beyond index 0 drive the dynamic
// in_seed_N / in_rhs_N ports. Legacy graphs instead stored data.rhs /
// data.letter_projection (pair 1) plus data.pairs_config (flat string list
// of letter projections for pairs 2+). This returns the canonical list,
// trimmed to the wired/configured pair count, without mutating the input.

export const PAIR_PRESETS = [
  ['identity', 'identity'],
  ['divergent', 'divergent (any letter)'],
  ['finite', 'finite (all letters)'],
]

export function normalizeCollinearPairs(data) {
  const proj = (v) => {
    let s = String(v ?? 'identity').trim()
    return s || 'identity'
  }
  if (Array.isArray(data?.pairs) && data.pairs.length) {
    return data.pairs.map((p) => ({
      projection: proj(p?.projection),
      rhs: typeof p?.rhs === 'string' ? p.rhs.trim() : '',
    }))
  }
  const rows = [{ projection: proj(data?.letter_projection), rhs: String(data?.rhs ?? '').trim() }]
  for (const s of Array.isArray(data?.pairs_config) ? data.pairs_config : []) {
    if (s && typeof s === 'object') rows.push({ projection: proj(s.projection), rhs: String(s.rhs ?? '').trim() })
    else if (String(s ?? '').trim()) rows.push({ projection: proj(s), rhs: '' })
  }
  return rows
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
