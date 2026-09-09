// Bundle scripts/ssr-audit.mjs for Node with the SSR stubs applied.
//
// esbuild's CLI --alias cannot express relative alias names ('../api',
// '../App' — rejected as invalid on current esbuild), so the stub routing
// lives here as an onResolve plugin: exact import specifiers are redirected
// to scripts/stubs/ regardless of the importing file. 'reactflow' and
// 'react-router-dom' keep working as plain aliases.
import { build } from 'esbuild'
import { fileURLToPath } from 'node:url'
import path from 'node:path'

const webRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..')
const stubs = path.join(webRoot, 'scripts', 'stubs')

const stubMap = {
  reactflow: path.join(stubs, 'reactflow.js'),
  '../api': path.join(stubs, 'api.js'),
  '../App': path.join(stubs, 'app.js'),
  'react-router-dom': path.join(stubs, 'app.js'),
}

await build({
  entryPoints: [path.join(webRoot, 'scripts', 'ssr-audit.mjs')],
  bundle: true,
  format: 'cjs',
  platform: 'node',
  jsx: 'automatic',
  // the audit entry itself is .mjs but embeds JSX elements
  loader: { '.mjs': 'jsx' },
  outfile: '/tmp/ssr-audit.cjs',
  logLevel: 'warning',
  plugins: [{
    name: 'ssr-stubs',
    setup(builder) {
      builder.onResolve(
        { filter: /^(reactflow|react-router-dom|\.\.\/api|\.\.\/App)$/ },
        (args) => {
          const target = stubMap[args.path]
          return target ? { path: target } : null
        },
      )
    },
  }],
})
