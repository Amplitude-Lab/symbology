import { renderToStaticMarkup } from 'react-dom/server'
import { existsSync, readFileSync, writeFileSync } from 'node:fs'
import { pathToFileURL } from 'node:url'
import path from 'node:path'
import { registerCustomBlocks, NODE_DEFS, sourceKindFor, portsOf, portBaseY, ROW_STEP, NODE_W } from '../src/flowdefs'
import { OpNode, Inspector, nodeTypes, routeEdges, KindEdge } from '../src/screens/FlowEditor'
import { Position } from './stubs/reactflow.js'

// Fixture projects live outside the bundle (gitignored); resolve them
// relative to this script (or the web root when bundled by ssr-audit.sh)
// and skip whatever is not present locally.
const baseDir = import.meta.url && String(import.meta.url).startsWith('file:')
  ? path.resolve(path.dirname(new URL(import.meta.url).pathname), '..', '..')
  : path.resolve(process.cwd(), '..')  // ssr-audit.sh runs from the web root; projects/ is one level up
const readProject = (name) => {
  const p = path.join(baseDir, 'projects', name, 'project.json')
  return existsSync(p) ? JSON.parse(readFileSync(p, 'utf8')) : null
}
const projects = [
  ['heptagon', readProject('heptagon')],
  ['test4p', readProject('test4p')],
].filter(([, p]) => p != null)

const projectGlobal = globalThis

let pass = 0
let fail = 0
const problems = []

function record(kind, flowName, nodeId, type, message) {
  fail++
  problems.push(`[${kind}] ${flowName}/${nodeId} (${type}): ${message}`)
}

function tryRender(label, flowName, nodeId, type, element) {
  try {
    const html = renderToStaticMarkup(element)
    if (typeof html !== 'string' || html.length === 0) {
      record('empty-render', flowName, nodeId, type, label + ' produced empty HTML')
    } else {
      pass++
    }
  } catch (e) {
    record('render-crash', flowName, nodeId, type, `${label}: ${e && e.message}`)
  }
  return null
}

// Path endpoints as KindEdge draws them: starts at (sourceX, sourceY+spread),
// ends exactly at (targetX, targetY). Returns {ok, why}.
function checkPathEndpoints(html, sx, sy, tx, ty, label, flowName, edgeId) {
  const dm = html.match(/d="([^"]+)"/)
  if (!dm) {
    record('kind-edge', flowName, edgeId, '-', `${label}: rendered markup has no path d attribute`)
    return
  }
  const nums = [...dm[1].matchAll(/-?\d+(?:\.\d+)?/g)].map((m) => parseFloat(m[0]))
  if (nums.length < 4) {
    record('kind-edge', flowName, edgeId, '-', `${label}: path too short: ${dm[1]}`)
    return
  }
  const EPS = 0.51
  const bad = []
  if (Math.abs(nums[0] - sx) > EPS) bad.push(`start x ${nums[0]} != ${sx}`)
  if (Math.abs(nums[1] - sy) > EPS) bad.push(`start y ${nums[1]} != ${sy}`)
  const ex = nums[nums.length - 2]
  const ey = nums[nums.length - 1]
  if (Math.abs(ex - tx) > EPS) bad.push(`end x ${ex} != ${tx}`)
  if (Math.abs(ey - ty) > EPS) bad.push(`end y ${ey} != ${ty}`)
  if (bad.length) {
    record('kind-edge', flowName, edgeId, '-', `${label}: ${bad.join('; ')}`)
  } else {
    pass++
  }
}

const groupOpsStub = {
  groups: [],
  onRename: () => {},
  onExpand: () => {},
  onUngroup: () => {},
}

const stats = { wires: 0, detours: 0, fallbacks: 0, kappa: 0 }
const dbgLines = []
const DBG = !!process.env.WIRE_DEBUG

for (const [projName, project] of projects) {
  projectGlobal.__PROJECT = project
  registerCustomBlocks(project)
  const flows = project.flows || []
  for (const flow of flows) {
    const flowLabel = `${projName}/${flow.name}`
    const nodes = flow.graph?.nodes || []
    const edges = flow.graph?.edges || []
    projectGlobal.__EDGES = edges
    for (const node of nodes) {
      const data = node.data || {}
      if (!nodeTypes[node.type]) {
        record('no-renderer', flowLabel, node.id, node.type, 'no nodeTypes entry')
        continue
      }
      const Renderer = nodeTypes[node.type]
      tryRender('RegistryNode', flowLabel, node.id, node.type, <Renderer id={node.id} type={node.type} data={data} selected={false} />)
      if (NODE_DEFS[node.type]) {
        tryRender('OpNode', flowLabel, node.id, node.type, <OpNode id={node.id} type={node.type} data={data} selected={false} />)
      }
      tryRender('Inspector', flowLabel, node.id, node.type, (
        <Inspector
          node={node}
          onChange={() => {}}
          onDelete={() => {}}
          groupOps={groupOpsStub}
          onOpenBlock={() => {}}
          flowId={flow.id}
        />
      ))
    }

    // Execute the REAL wire-routing pipeline (routeEdges) on every flow — the
    // exact component-computed visibleEdges the browser renders. The alphabet
    // def.outputs crash shipped precisely because no audit ran this code path.
    const groupOf = new Map()
    try {
      const routed = routeEdges({
        edges,
        nodes,
        project,
        groupOf,
        resolveSourceKind: (nodeId, handleId) => sourceKindFor(
          nodes.find((n) => n.id === nodeId) || {}, handleId, project),
      })
      if (!Array.isArray(routed) || routed.length !== edges.length) {
        record('route-edges', flowLabel, '-', '-', `expected ${edges.length} routed edges, got ${routed.length}`)
      } else {
        // Port coordinates exactly as routeEdges derives them (portsOf is
        // the shared truth incl. dynamic ports); unknown handles → row 0.
        const nodeById = new Map(nodes.map((n) => [n.id, n]))
        const portY = (n, handle, side) => {
          const { inputs, outputs } = portsOf(n, project, edges)
          const list = side === 'out' ? outputs : inputs
          const row = Math.max(0, list.findIndex((p) => p.id === handle))
          return n.position.y + portBaseY(n) + row * ROW_STEP
        }
        for (const re of routed) {
          const edgeId = re.id || '-'
          if (!re.style || !re.style.stroke) {
            record('route-edges', flowLabel, edgeId, '-', 'routed edge missing style.stroke')
            continue
          }
          if (String(edgeId).startsWith('px_')) {
            // With an empty groupOf map the proxy branch is only reachable
            // for dangling source/target node refs — contract-audit's turf,
            // but flag it here too rather than silently passing.
            record('route-edges', flowLabel, edgeId, '-', 'proxy px_ edge emitted with no groups (dangling node ref?)')
            continue
          }
          const d = re.data
          if (!d || typeof d.spread !== 'number' || typeof d.kappa !== 'number') {
            record('route-edges', flowLabel, edgeId, '-', 'wire data missing spread/kappa numbers')
            continue
          }
          if (d.detourY !== undefined) {
            if (typeof d.detourY !== 'number' || typeof d.outboundX !== 'number' || typeof d.inboundX !== 'number') {
              record('route-edges', flowLabel, edgeId, '-', 'detour wire has non-numeric detourY/outboundX/inboundX')
              continue
            }
            stats.detours++
          }
          if (d.fallback !== undefined && d.fallback !== true) {
            record('route-edges', flowLabel, edgeId, '-', 'fallback flag must be true when present')
            continue
          }
          if (d.fallback) stats.fallbacks++
          if (d.kappa) stats.kappa++
          stats.wires++
          // Render each wire through the REAL KindEdge with true port
          // coordinates and verify the drawn path pins both endpoints.
          const sn = nodeById.get(re.source)
          const tn = nodeById.get(re.target)
          if (!sn || !tn) {
            record('route-edges', flowLabel, edgeId, '-', 'routed wire endpoints not in nodes')
            continue
          }
          const sourceX = sn.position.x + NODE_W
          const sourceY = portY(sn, re.sourceHandle, 'out')
          const targetX = tn.position.x
          const targetY = portY(tn, re.targetHandle, 'in')
          if (DBG) {
            const dec = d.fallback ? 'fallback' : d.detourY !== undefined ? 'detour' : d.kappa ? 'kappa' : 'direct'
            dbgLines.push([
              flowLabel, edgeId, re.source, re.sourceHandle || '', re.target, re.targetHandle || '', dec,
              `${sourceX.toFixed(1)},${(sourceY + (d.spread || 0)).toFixed(1)},${targetX.toFixed(1)},${targetY.toFixed(1)}`,
              `s=${d.spread}`, `k=${d.kappa}`,
            ].join('|'))
          }
          const el = (
            <KindEdge
              id={edgeId}
              sourceX={sourceX}
              sourceY={sourceY}
              targetX={targetX}
              targetY={targetY}
              sourcePosition={Position.Right}
              targetPosition={Position.Left}
              style={re.style}
              data={d}
            />
          )
          try {
            const html = renderToStaticMarkup(el)
            if (typeof html !== 'string' || html.length === 0) {
              record('kind-edge', flowLabel, edgeId, '-', 'KindEdge produced empty HTML')
            } else {
              checkPathEndpoints(html, sourceX, sourceY + (d.spread || 0), targetX, targetY,
                d.fallback ? 'fallback' : d.detourY !== undefined ? 'detour' : d.kappa ? 'kappa' : 'straight',
                flowLabel, edgeId)
            }
          } catch (e) {
            record('kind-edge', flowLabel, edgeId, '-', `KindEdge render: ${e && e.message}`)
          }
        }
        pass++
      }
    } catch (e) {
      record('route-edges', flowLabel, '-', '-', `${e && e.message}`)
    }
  }
}
console.log(`SSR audit: ${pass} renders OK, ${fail} problems`)
console.log(`wires: ${stats.wires} (detours ${stats.detours}, kappa bows ${stats.kappa}, fallbacks ${stats.fallbacks})`)
if (dbgLines.length) {
  writeFileSync(process.env.WIRE_DEBUG, dbgLines.join('\n') + '\n')
}
if (problems.length) {
  console.log('--- problems ---')
  for (const p of problems) console.log(p)
  process.exit(1)
}
