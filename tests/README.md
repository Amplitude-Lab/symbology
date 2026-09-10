# Portable regression checks

These tests create disposable projects and small exact tensors. They do not
read private `front-end/projects` data. The original dated audit under
`audits/robustness-2026-09-09` is historical evidence; its original probes
intentionally expose additional issues outside this refinement.

From the repository root, after installing the documented native dependencies:

```bash
./scripts/setup-sparserref.sh
make check-public
make check-wxf-sanitized # optional local memory diagnostics

python3 -m venv .venv-tests
.venv-tests/bin/pip install -r tests/requirements.txt
.venv-tests/bin/python -m pytest tests/test_robustness.py tests/test_export_runtime.py -q
```

The native test checks both solver helpers on 90 exact rational systems,
including residuals, uniqueness, and null spaces. It tests six condition
export/reimport cases, strict scalar parsing, and regenerates all ten public
E6 two-/three-loop reference outputs, including a partial-cache resume.
Compiler/toolchain overrides use the normal Makefile variables.

The API tests cover transactional edits, revision conflicts, cancellation
while launching, output isolation/rollback, content cache integrity, SSE buffer
rollover/reconnect, restart records, graph validation, and generated-script
tests for relocation, spaces/apostrophes, changed input contents, Wolfram
availability modes, shared project locks, missing outputs and timeout cleanup.

For browser tests, first build `front-end/web`, install Playwright and its
Chromium browser in a test environment, then run these in separate terminals:

```bash
.venv-tests/bin/python tests/serve_browser.py
PLAYWRIGHT_MODULE=/absolute/path/to/playwright node tests/browser_regression.cjs
```

`TEST_BASE_URL` defaults to localhost port 18321. `CHROME_PATH` can select an
installed Chrome instead of Playwright Chromium. `TEST_OUTPUT_DIR` defaults to
`/tmp/symbology-browser-results`. The server uses disposable project storage.
The browser test covers pending/delayed saves, two-tab conflicts, draft
recovery, connection errors, and 768-/390-pixel Materials views. Full phone flow editing and other
browser engines remain separate acceptance targets.

`check-public` also verifies native cache invalidation after seed/output
changes, SHA256 padding boundaries, malformed WXF inputs, and the corrected
examples. The template test checks exact rational residuals and certifies the
complete invariant spaces using an independent modular rank calculation.
Sanitizer tests cover both matrix and tensor readers; they require a runtime
that permits AddressSanitizer and LeakSanitizer to inspect the process.
