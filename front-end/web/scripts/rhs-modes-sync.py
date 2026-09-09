#!/usr/bin/env python3
"""Check that the two RHS_MODES registries agree.

Compute RHS object types are registered twice by design (the web UI needs them
for the palette/Inspector, the compiler needs them for the generation rule):
  - front-end/web/src/flowdefs.js   (RHS_MODES: label/short/requires/forbids)
  - front-end/server/app/compile.py (RHS_MODES: label/requires/forbids)

A drifted registry means the UI offers an object type the compiler rejects (or
vice versa). This script parses both and exits nonzero on any mismatch in the
key set, labels, or requires/forbids fields. Run it after editing either side;
it is also suitable for CI.
"""
import re
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parents[1] / "src" / "flowdefs.js"
SERVER = Path(__file__).resolve().parents[2] / "server" / "app" / "compile.py"


def parse_js(text):
    m = re.search(r"export const RHS_MODES = \{(.*?)\n\}", text, re.S)
    if not m:
        raise SystemExit("flowdefs.js: RHS_MODES block not found")
    out = {}
    for key, body in re.findall(r"^\s{2}(\w+): \{(.*?)\},?$", m.group(1), re.S | re.M):
        entry = {}
        if lm := re.search(r"label:\s*'([^']*)'", body):
            entry["label"] = lm.group(1)
        if rm := re.search(r"requires:\s*\[(.*?)\]", body):
            entry["requires"] = tuple(sorted(x.strip("'\" ") for x in rm.group(1).split(",") if x.strip()))
        if fm := re.search(r"forbids:\s*\[(.*?)\]", body):
            entry["forbids"] = tuple(sorted(x.strip("'\" ") for x in fm.group(1).split(",") if x.strip()))
        out[key] = entry
    return out


def parse_py(text):
    m = re.search(r"^RHS_MODES = \{(.*?)^\}", text, re.S | re.M)
    if not m:
        raise SystemExit("compile.py: RHS_MODES block not found")
    out = {}
    for key, body in re.findall(r'^\s{4}"(\w+)": \{(.*?)\},?$', m.group(1), re.S | re.M):
        entry = {}
        if lm := re.search(r'"label":\s*"([^"]*)"', body):
            entry["label"] = lm.group(1)
        if rm := re.search(r'"requires":\s*\((.*?)\)', body):
            entry["requires"] = tuple(sorted(x.strip('"\' ') for x in rm.group(1).split(",") if x.strip()))
        if fm := re.search(r'"forbids":\s*\((.*?)\)', body):
            entry["forbids"] = tuple(sorted(x.strip('"\' ') for x in fm.group(1).split(",") if x.strip()))
        out[key] = entry
    return out


def main():
    js = parse_js(WEB.read_text())
    py = parse_py(SERVER.read_text())
    bad = []
    for key in sorted(set(js) | set(py)):
        if key not in js:
            bad.append(f"  '{key}' exists in compile.py but is missing from flowdefs.js")
            continue
        if key not in py:
            bad.append(f"  '{key}' exists in flowdefs.js but is missing from compile.py")
            continue
        for field in ("label", "requires", "forbids"):
            if js[key].get(field) != py[key].get(field):
                bad.append(f"  '{key}'.{field}: flowdefs.js has {js[key].get(field)!r} vs compile.py {py[key].get(field)!r}")
    if bad:
        print("RHS_MODES registries are OUT OF SYNC:")
        print("\n".join(bad))
        sys.exit(1)
    print(f"RHS_MODES registries in sync ({len(js)} object types: {', '.join(sorted(js))})")


if __name__ == "__main__":
    main()
