#!/usr/bin/env python3
"""Record or verify ignored large files without involving Git's content scanner.

From the repository root:
  python3 scripts/large_file_checksums.py status
  python3 scripts/large_file_checksums.py verify [repository/relative/path ...]
  sha256sum --check docs/large_file_checksums.sha256

Status checks metadata only; verify reads all selected contents. Neither changes
the saved baseline. Use record only when deliberately accepting new contents.
The initial inventory covers previously unignored files of at least 10 MiB;
already ignored output directories retain their existing provenance records.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "docs/large_file_checksums.json"


def signature(stat):
    return {name: getattr(stat, "st_" + name)
            for name in ("size", "mtime_ns", "ctime_ns", "dev", "ino")}


def local_path(relative):
    path = ROOT / relative
    if Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ValueError(f"Expected a repository-relative path: {relative}")
    if path.is_symlink() or not path.resolve().is_relative_to(ROOT):
        raise ValueError(f"Refusing a symlink or external path: {relative}")
    return path


def digest(relative):
    path = local_path(relative)
    before = signature(path.stat())
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            sha.update(block)
    after = signature(path.stat())
    if before != after:
        raise RuntimeError(f"File changed while hashing: {relative}")
    return {"path": relative, "sha256": sha.hexdigest(), **after}


def atomic_write(path, text):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("record", "verify", "status"))
    parser.add_argument("paths", nargs="*", help="Optional subset of recorded paths")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--paths-file", type=Path,
                        help="record only: one repository-relative path per line")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error("Use between one and four workers")
    if args.paths_file and args.command != "record":
        parser.error("--paths-file is only valid with record")
    manifest = args.manifest.resolve()
    old = json.loads(manifest.read_text()) if manifest.exists() else {"files": []}
    baseline = {entry["path"]: entry for entry in old["files"]}
    paths = args.paths or (args.paths_file.read_text().splitlines()
                          if args.paths_file else list(baseline))
    paths = sorted(set(paths))
    if not paths:
        parser.error("No files selected")
    if args.command != "record" and any(p not in baseline for p in paths):
        parser.error("All verification paths must already be in the manifest")
    failures = []
    if args.command == "status":
        for relative in paths:
            try:
                now = signature(local_path(relative).stat())
                if any(now[key] != baseline[relative][key] for key in now):
                    failures.append(relative + ": metadata changed; run verify")
            except (OSError, ValueError) as exc:
                failures.append(f"{relative}: {exc}")
        print("Metadata check only; this does not verify file contents.")
    else:
        started = last_report = time.monotonic()
        total_bytes = 0
        results = {}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            pending = {pool.submit(digest, path): path for path in paths}
            for done, future in enumerate(as_completed(pending), 1):
                relative = pending[future]
                try:
                    result = future.result()
                    results[relative] = result
                    total_bytes += result["size"]
                    if (args.command == "verify" and
                            result["sha256"] != baseline[relative]["sha256"]):
                        failures.append(relative + ": SHA-256 MISMATCH")
                except (OSError, RuntimeError, ValueError) as exc:
                    failures.append(f"{relative}: {exc}")
                if time.monotonic() - last_report >= 15 or done == len(paths):
                    print(f"{done}/{len(paths)} files; {total_bytes / 1e9:.3f} GB; "
                          f"{time.monotonic() - started:.1f} s", flush=True)
                    last_report = time.monotonic()
        if args.command == "record" and not failures:
            baseline.update(results)
            entries = [baseline[p] for p in sorted(baseline)]
            payload = {"schema": 1, "algorithm": "SHA-256",
                       "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
                       "scope": "Explicitly ignored large local research files; "
                                "source code and small provenance records stay visible to Git.",
                       "files": entries}
            atomic_write(manifest, json.dumps(payload, indent=2) + "\n")
            # Standard sha256sum --check format for the repository's plain paths.
            if all(not any(c in e["path"] for c in "\n\r\\") for e in entries):
                atomic_write(manifest.with_suffix(".sha256"), "".join(
                    f'{e["sha256"]}  {e["path"]}\n' for e in entries))
    for failure in failures:
        print(failure, file=sys.stderr)
    print(f"{args.command}: {len(paths)} files, {len(failures)} failures")
    return bool(failures)


if __name__ == "__main__":
    sys.exit(main())
