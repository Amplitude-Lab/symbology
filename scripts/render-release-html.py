"""Render the release note (and the guide) to standalone styled HTML.

Usage:  python3 scripts/render-release-html.py [docs-dir]   # default: ./docs

The HTML keeps the GIF paths relative (../gifs/...), so the exact same file
works opened over file:// by double-click AND when the server serves
docs/ under /docs/ — no build step, no CDN, no JS. Needs the `markdown`
package (pip install markdown). The Markdown files stay the source of
truth; re-run this after editing them.
"""
import re
import sys
from pathlib import Path

import markdown

DOCS = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "docs").resolve()

STYLE = """
:root { color-scheme: light; }
* { box-sizing: border-box; }
body {
  margin: 0; padding: 0 16px 96px;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", sans-serif;
  color: #1f2328; background: #ffffff; line-height: 1.6; font-size: 16px;
}
main { max-width: 780px; margin: 0 auto; }
h1 { font-size: 34px; line-height: 1.2; margin: 48px 0 8px; letter-spacing: -0.5px; }
h1 + p strong:first-child { color: #0969da; }
h2 {
  font-size: 22px; margin: 44px 0 12px; padding-bottom: 8px;
  border-bottom: 1px solid #d8dee4; letter-spacing: -0.2px;
}
h3 { font-size: 17px; margin: 28px 0 8px; }
a { color: #0969da; text-decoration: none; }
a:hover { text-decoration: underline; }
blockquote {
  margin: 16px 0; padding: 10px 16px; border-left: 4px solid #0969da;
  background: #f6f8fa; border-radius: 0 6px 6px 0; color: #57606a;
}
blockquote a { color: #0969da; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 87%;
       background: #f6f8fa; padding: 2px 5px; border-radius: 4px; }
pre { background: #f6f8fa; padding: 12px 16px; border-radius: 8px; overflow-x: auto; }
pre code { background: none; padding: 0; }
img {
  display: block; max-width: 100%; margin: 14px auto 22px;
  border: 1px solid #d8dee4; border-radius: 10px;
  box-shadow: 0 1px 6px rgba(31,35,40,.08);
}
table { border-collapse: collapse; margin: 14px 0; width: 100%; }
th, td { border: 1px solid #d8dee4; padding: 6px 12px; text-align: left; }
th { background: #f6f8fa; }
tr:nth-child(even) td { background: #fbfcfd; }
ol li, ul li { margin: 6px 0; }
hr { border: none; border-top: 1px solid #d8dee4; margin: 32px 0; }
footer { max-width: 780px; margin: 48px auto 0; color: #57606a; font-size: 13px;
         border-top: 1px solid #d8dee4; padding-top: 12px; }
"""

TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{style}</style>
</head>
<body>
<main>
{body}
</main>
<footer>Rendered from <code>{source}</code> in the symbology repository — the
Markdown file is the source of truth; regenerate the HTML with
<code>scripts/render-release-html.py</code>.</footer>
</body>
</html>
"""


def render(md_path: Path, out_path: Path, title: str, rewrites=None):
    text = md_path.read_text()
    for old, new in (rewrites or {}).items():
        text = text.replace(old, new)
    body = markdown.markdown(text, extensions=["tables", "toc", "sane_lists", "fenced_code"])
    out_path.write_text(TEMPLATE.format(title=title, style=STYLE, body=body, source=md_path.name))
    print("wrote", out_path)


render(
    DOCS / "releases" / "2026-09.md",
    DOCS / "releases" / "2026-09.html",
    "Symbology Studio — September 2026 update",
    # The in-app page and the double-clicked file both resolve relative
    # links; point the guide reference at its rendered twin.
    {"(../front-end-guide.md)": "(../front-end-guide.html)"},
)
render(
    DOCS / "front-end-guide.md",
    DOCS / "front-end-guide.html",
    "Operating the front-end — an illustrated guide",
    {"../README.md#the-visual-front-end-web-editor": "https://github.com/Amplitude-Lab/symbology/blob/front-end-dev/README.md",
     "../README.md#design-locally-run-on-the-cluster-standalone-script-export": "https://github.com/Amplitude-Lab/symbology/blob/front-end-dev/README.md"},
)
