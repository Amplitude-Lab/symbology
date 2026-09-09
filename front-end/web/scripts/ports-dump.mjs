// Dump the TRUE rendered port rows for every node of every flow, straight
// from the shared geometry model in src/flowdefs.js (portsOf/portBaseY/
// estHeight — the same functions FlowEditor renders and routes with).
// Consumed by wire-audit.py so the Python audit always mirrors the renderer
// exactly; a mismatch here is a code bug, never a silent fallback.
//
// usage: node scripts/ports-dump.mjs <project.json>
// stdout: {"<flow name>": {"<node id>": {in: [...], out: [...], baseY, h}}}
import { readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { portsOf, portBaseY, estHeight, registerCustomBlocks } from '../src/flowdefs.js'

const projPath = process.argv[2]
if (!projPath) {
  console.error('usage: node ports-dump.mjs <project.json>')
  process.exit(2)
}
const project = JSON.parse(readFileSync(projPath, 'utf8'))
registerCustomBlocks(project)
const dump = {}
for (const flow of project.flows || []) {
  const nodes = flow.graph?.nodes || []
  const edges = flow.graph?.edges || []
  const fd = {}
  for (const n of nodes) {
    const p = portsOf(n, project, edges)
    fd[n.id] = {
      in: p.inputs.map((x) => x.id),
      out: p.outputs.map((x) => x.id),
      baseY: portBaseY(n),
      h: estHeight(n, project, edges),
    }
  }
  dump[flow.name] = fd
}
process.stdout.write(JSON.stringify(dump))
