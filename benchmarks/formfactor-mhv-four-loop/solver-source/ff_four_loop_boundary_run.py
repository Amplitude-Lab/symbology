#!/usr/bin/env python3
"""Prepare full collinear image bases independently of the weight-eight solve."""
from pathlib import Path
import hashlib
import json
import os
import signal
import subprocess
import sys
import time

base = Path(__file__).resolve().parents[1]
root = base.parents[1]
experiment = base / "results/ff_four_loop"
kind, target = sys.argv[1:3]
memory_cap = int(sys.argv[3]) if len(sys.argv) > 3 else 4_000_000_000
assert kind in ("DC", "TC")
assert target in ("FEC_4", "FEC_6", "LEC_2")
out = experiment / "boundary" / kind
data = out / "data"
data.mkdir(parents=True, exist_ok=True)
source_map = base / "results" / ("ff_boot_boundary_map.wxf" if kind == "DC" else "ff_boot_triple_map.wxf")
links = {data / "colmat42.wxf": source_map,
         data / "FEC_1.wxf": experiment / "data/FEC_1.wxf",
         data / "LEC_1.wxf": experiment / "data/LEC_1.wxf"}
links.update({out / f"FEC_{w}.wxf": experiment / f"FEC_{w}.wxf" for w in range(2, 7)})
links[out / "LEC_2.wxf"] = experiment / "LEC_2.wxf"
for dest, source in links.items():
    assert source.is_file(), source
    if dest.is_symlink() or dest.exists():
        assert dest.resolve() == source.resolve(), dest
    else:
        dest.symlink_to(source)

name = f"project_{target}" + (f"_cap{memory_cap // 1_000_000_000}GB" if len(sys.argv) > 3 else "")
resources = out / f"{name}.resources.json"
command = [str(root / "bootstrap"), "--project", "--symmetry", "collinear",
           "--target", target, "--data-dir", str(data), "--output-dir", str(out),
           "--threads", "1"]
env = dict(os.environ, LD_LIBRARY_PATH=str(root / "output_hexagon_runtime_recovery_20261008/env/lib"),
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
start = time.monotonic()
peak = 0
stop_reason = None
samples = []
with (out / f"{name}.log").open("w") as log:
    process = subprocess.Popen(["/usr/bin/time", "-f", '{"user":%U,"system":%S,"maxrss_kib":%M,"exit":%x}',
                                "-o", str(resources), *command], env=env,
                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    while process.poll() is None:
        family = [process.pid]
        for pid in family:
            children = Path(f"/proc/{pid}/task/{pid}/children")
            try:
                family.extend(map(int, children.read_text().split()))
            except FileNotFoundError:
                pass
        rss = 0
        for pid in family:
            try:
                rss += int(Path(f"/proc/{pid}/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE")
            except FileNotFoundError:
                pass
        peak = max(peak, rss)
        elapsed = time.monotonic() - start
        if rss > memory_cap:
            stop_reason = f"{memory_cap} byte independent-preparation process cap"
        main_rss = 0
        try:
            main = json.loads((experiment / "progress.json").read_text())
            if main.get("state") == "running" and main.get("pid"):
                main_family = [main["pid"]]
                for pid in main_family:
                    children = Path(f"/proc/{pid}/task/{pid}/children")
                    if children.exists(): main_family.extend(map(int, children.read_text().split()))
                    main_rss += int(Path(f"/proc/{pid}/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE")
        except (FileNotFoundError, json.JSONDecodeError):
            pass
        if rss + main_rss > 40_000_000_000:
            stop_reason = "main plus boundary40GB guard, reserving10GB of shared50GB for audits/other processes"
        if target == "FEC_4" and elapsed > 180:
            stop_reason = "bounded low-weight profile exceeded180seconds"
        if stop_reason:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait()
            break
        if not samples or elapsed - samples[-1][0] >= 10:
            samples.append([round(elapsed, 2), rss])
        (out / "progress.json").write_text(json.dumps({"state":"running", "target":target,
                "seconds":elapsed, "rss_peak_bytes":peak, "pid":process.pid}, indent=2) + "\n")
        time.sleep(.2)
receipt = {"state":"pass" if process.returncode == 0 else "stopped" if stop_reason else "failed",
           "boundary":kind, "target":target, "command":command,
           "wall_seconds":time.monotonic()-start, "peak_sampled_rss_bytes":peak,
           "returncode":process.returncode, "stop_reason":stop_reason, "samples":samples,
           "boundary_map_sha256":hashlib.sha256(source_map.read_bytes()).hexdigest(),
           "bootstrap_sha256":hashlib.sha256((root / "bootstrap").read_bytes()).hexdigest(),
           "boundary_alphabet_size":22 if kind == "DC" else 29,
           "own_memory_cap_bytes":memory_cap, "combined_main_boundary_cap_bytes":40_000_000_000,
           "all_spurious_directions_retained":True, "finite_targets_read":False,
           "resources":resources.read_text() if resources.exists() else None}
(out / f"{name}.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
(out / "progress.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps({k:receipt[k] for k in ("state","boundary","target","wall_seconds","peak_sampled_rss_bytes","stop_reason")}))
raise SystemExit(process.returncode != 0)
