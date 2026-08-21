import React from 'react'
import ReactDOM from 'react-dom/client'
import { BrowserRouter } from 'react-router-dom'
import App from './App'
import './styles.css'

function reportClientError(kind, message, stack) {
  try {
    fetch('/api/client-log', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        kind,
        message: String(message).slice(0, 4000),
        stack: String(stack || '').slice(0, 4000),
        url: window.location.href,
        ua: navigator.userAgent,
      }),
    })
  } catch { /* best effort */ }
}

window.addEventListener('error', (e) => {
  reportClientError('error', e.message, e.error?.stack || `${e.filename}:${e.lineno}:${e.colno}`)
})
window.addEventListener('unhandledrejection', (e) => {
  reportClientError('unhandledrejection', e.reason?.message || String(e.reason), e.reason?.stack)
})

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <BrowserRouter>
      <App />
    </BrowserRouter>
  </React.StrictMode>
)
