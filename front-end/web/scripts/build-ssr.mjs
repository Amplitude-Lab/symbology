import { build } from 'esbuild'
import { fileURLToPath } from 'node:url'

const stubs = fileURLToPath(new URL('./stubs/', import.meta.url))

const stubPlugin = {
  name: 'ssr-stubs',
  setup(build) {
    build.onResolve({ filter: /\.css$/ }, () => ({ path: stubs + 'empty.css', namespace: 'stub' }))
    build.onLoad({ filter: /.*/, namespace: 'stub' }, (args) => ({ contents: '', loader: 'css' }))
    build.onResolve({ filter: /^reactflow$/ }, () => ({ path: stubs + 'reactflow.js' }))
    build.onResolve({ filter: /(\.\.\/api|\.\.\/App|react-router-dom)$/ }, (args) => ({ path: stubs + 'app.js' }))
  },
}

await build({
  entryPoints: [new URL('./ssr-audit.mjs', import.meta.url).pathname],
  bundle: true,
  format: 'cjs',
  platform: 'node',
  jsx: 'automatic',
  loader: { '.mjs': 'jsx' },
  plugins: [stubPlugin],
  outfile: '/tmp/ssr-audit.cjs',
  logLevel: 'warning',
})
console.log('built /tmp/ssr-audit.cjs')
