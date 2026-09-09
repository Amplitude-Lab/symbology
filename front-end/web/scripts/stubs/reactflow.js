import React from 'react'

export const useEdges = () => globalThis.__EDGES || []
export const useUpdateNodeInternals = () => () => {}
export const useEdgesState = (v) => [v || [], () => {}]
export const useNodesState = (v) => [v || [], () => {}, () => {}]
export const useReactFlow = () => ({ fitView: () => {}, getNodes: () => [] })
export const addEdge = (p) => p
export const ReactFlow = () => null
export const ReactFlowProvider = ({ children }) => children
export const Background = () => null
export const Controls = () => null
export const MiniMap = () => null
export const Handle = () => null
export const Position = { Top: 'top', Bottom: 'bottom', Left: 'left', Right: 'right' }
export const ConnectionMode = { Loose: 'loose', Strict: 'strict' }
export const MarkerType = { Arrow: 'arrow', ArrowClosed: 'arrowclosed' }

// Faithful-enough stand-ins so SSR exercises the real KindEdge path math
// (the browser renders <path d={path}/>; a null-returning BaseEdge would
// silently skip that code path in audits).
export const BaseEdge = ({ path }) => React.createElement('path', { d: path })
export const getBezierPath = ({ sourceX, sourceY, targetX, targetY }) => [
  `M ${sourceX} ${sourceY} C ${(sourceX + targetX) / 2} ${sourceY} ${(sourceX + targetX) / 2} ${targetY} ${targetX} ${targetY}`,
  0,
  0,
]
export default ReactFlow
