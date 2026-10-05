#!/usr/bin/env python3
"""Compare three exact solving orders using fresh, sequential processes on Linux.

Preparation and independent exact validation are measured separately. Per-run
budgets apply equally to each method; a failed method is never a timing win.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("integrability-first", "joint-global", "joint-staged")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def read_record(log, prefix):
    records = [json.loads(line[len(prefix):]) for line in log.read_text().splitlines()
               if line.startswith(prefix)]
    if len(records) != 1:
        raise RuntimeError(f"Expected one {prefix.strip()} record in {log}")
    return records[0]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases", type=Path, default=ROOT/"bench/combined_constraints_cases.json")
    p.add_argument("--program", type=Path, default=ROOT/"bench/combined_constraints_benchmark")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=float, default=600)
    p.add_argument("--virtual-gib", type=int, default=8)
    p.add_argument("--rss-gib", type=float, default=6)
    p.add_argument("--min-available-gib", type=float, default=4)
    p.add_argument("--local-reduction", choices=("raw", "reduced"), default="raw")
    a = p.parse_args()
    if min(a.threads, a.repeats, a.timeout, a.virtual_gib, a.rss_gib) <= 0:
        p.error("thread/repeat counts and process budgets must be positive")
    if a.min_available_gib < 0:
        p.error("host memory floor must be nonnegative")
    cases_file = a.cases.resolve()
    cases = json.loads(cases_file.read_text())
    if not isinstance(cases, list) or not cases:
        p.error("cases must be a nonempty JSON list")
    names = [c["name"] for c in cases]
    if len(set(names)) != len(names) or any(not n or any(x not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for x in n) for n in names):
        p.error("case names must be unique lowercase letters, digits, hyphens or underscores")
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    snapshot = out/"snapshot"; snapshot.mkdir(); (snapshot/"bin").mkdir()
    program = snapshot/"bin/combined_constraints_benchmark"
    shutil.copy2(a.program.resolve(), program)
    monitor = snapshot/"measure_process.py"
    shutil.copy2(ROOT/"bench/measure_process.py", monitor)
    sources = [ROOT/"Makefile", Path(__file__).resolve(), ROOT/"bench/combined_constraints_benchmark.cpp",
               *ROOT.glob("*.hpp"), *ROOT.glob("*.h"), *ROOT.glob("SparseRREF/*.hpp"), *ROOT.glob("SparseRREF/*.h")]
    for source in sources:
        target = snapshot/"sources"/source.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    def host_command(command):
        r = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        return r.stdout.strip() if not r.returncode else None
    result = {"schema":"combined-constraints-benchmark-v1", "state":"running",
              "started_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "settings":{k:v for k,v in vars(a).items() if not isinstance(v, Path)},
              "host":{"platform":platform.platform(),"logical_cpus":os.cpu_count(),
                      "cpu":host_command(["lscpu"]),"git_head":host_command(["git","rev-parse","HEAD"]),
                      "git_status":host_command(["git","status","--short"]),
                      "meminfo_at_start":Path("/proc/meminfo").read_text(),
                      "load_average_at_start":os.getloadavg()},
              "cases":cases,"preparations":[],"runs":[],"checks":[],"failures":[],
              "snapshot_hashes":{str(f.relative_to(snapshot)):digest(f) for f in snapshot.rglob("*") if f.is_file()}}
    save(out/"results.json",result); save(out/"cases.json",cases)

    def measured(label, command):
        log, report, rss = out/(label+".log"), out/(label+".json"), out/(label+".rss")
        command = [str(x) for x in command]
        wrapped = [sys.executable,str(monitor),"--report",str(report),"--log",str(log),
                   "--timeout",str(a.timeout),"--virtual-gib",str(a.virtual_gib),"--rss-gib",str(a.rss_gib),
                   "--min-available-gib",str(a.min_available_gib),"--","/usr/bin/time","-f","%M","-o",str(rss),*command]
        r = subprocess.run(wrapped, cwd=ROOT, capture_output=True, text=True)
        value = json.loads(report.read_text()) if report.exists() else {"monitor_error":r.stderr,"exit_code":r.returncode}
        value.update(label=label,command=command,monitor_exit_code=r.returncode)
        if rss.exists():
            last = rss.read_text().splitlines()
            if last and last[-1].isdigit():value["peak_process_rss_kib"] = int(last[-1])
        return value,log

    prepared = {}
    for case in cases:
        name = case["name"]; target = out/"inputs"/name; target.parent.mkdir(exist_ok=True)
        resolve = lambda key: (cases_file.parent/Path(case[key])).resolve()
        input_paths = []
        if case["kind"] == "toy":
            command = [program,"toy",case["left_weight"],case["right_weight"],case.get("coordinates","permutation"),target,a.threads]
        elif case["kind"] == "bootstrap":
            generators = [(cases_file.parent/Path(g)).resolve() for g in case["generators"]]
            input_paths = [resolve(k) for k in ("condition","first","last")]+generators
            command = [program,"prepare",*input_paths[:3],case["left_weight"],case["right_weight"],"scalar",target,a.threads,*generators]
        elif case["kind"] == "bundle":
            shutil.copytree(resolve("directory"),target)
            prepared[name]=target
            result["preparations"].append({"case":name,"copied_bundle":str(resolve("directory")),"excluded_from_timing":True})
            continue
        else:raise ValueError("unknown case kind: "+case["kind"])
        raw = out/"raw-inputs"/name; raw.mkdir(parents=True)
        recorded = []
        for i,path in enumerate(input_paths):
            copy = raw/(str(i)+"-"+path.name);shutil.copy2(path,copy)
            recorded.append({"original":str(path),"snapshot":str(copy.relative_to(out)),"sha256":digest(copy)})
        # Use frozen copies for preparation, not mutable external data.
        if input_paths:
            replacements = {str(path):str(out/r["snapshot"]) for path,r in zip(input_paths,recorded)}
            command = [replacements.get(str(x),x) for x in command]
        measurement,log = measured(name+"-prepare",command)
        measurement.update(case=name,raw_inputs=recorded)
        result["preparations"].append(measurement)
        if measurement["monitor_exit_code"] or "PREPARE_PASS" not in log.read_text():
            result["failures"].append(name+": preparation failed")
        else:prepared[name]=target
        save(out/"results.json",result)
    result["prepared_input_hashes"]={str(f.relative_to(out)):digest(f) for target in prepared.values() for f in target.iterdir() if f.is_file()}
    save(out/"replay_cases.json", [{"name":name,"kind":"bundle","directory":str(path.relative_to(out))} for name,path in prepared.items()])
    for repeat in range(a.repeats):
        order = METHODS[repeat%len(METHODS):]+METHODS[:repeat%len(METHODS)]
        for name,target in prepared.items():
            for method in order:
                label=f"{name}-{method}-{repeat}";output=out/(label+".wxf")
                measurement,log=measured(label,[program,"run",method,target,output,a.threads,a.local_reduction])
                measurement.update(case=name,method=method,repeat=repeat,output=str(output.relative_to(out)))
                if measurement["monitor_exit_code"]:
                    result["failures"].append(label+": solve failed")
                else:
                    try:measurement.update(read_record(log,"BENCH_RESULT "));measurement["output_sha256"]=digest(output)
                    except (RuntimeError,ValueError,OSError) as error:
                        result["failures"].append(label+": "+str(error));measurement["record_error"]=str(error)
                result["runs"].append(measurement);save(out/"results.json",result)
                print(label,measurement.get("solve_seconds"),measurement.get("peak_process_rss_kib"),flush=True)
    for name,target in prepared.items():
        outputs=[out/r["output"] for r in result["runs"] if r["case"]==name and "output_sha256" in r]
        if not outputs:continue
        measurement,log=measured(name+"-independent-check",[program,"check",target,a.threads,*outputs])
        measurement["case"]=name
        if measurement["monitor_exit_code"]:result["failures"].append(name+": independent check failed")
        else:
            try:measurement.update(read_record(log,"CHECK_RESULT "))
            except (RuntimeError,ValueError) as error:result["failures"].append(name+": "+str(error))
        result["checks"].append(measurement);save(out/"results.json",result)
    for path,expected in result["prepared_input_hashes"].items():
        if digest(out/path)!=expected:result["failures"].append("input changed: "+path)
    result["state"]="failed" if result["failures"] else "complete"
    result["finished_utc"]=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime());save(out/"results.json",result)
    lines=["# Combined integrability and symmetry benchmark", "", "State: **"+result["state"]+"**.", "",
           "Prepared inputs are identical across methods. Preparation and independent verification are excluded from the solve table. All methods use the same binary, workers, interface policy and per-process guards. Values are medians of fresh sequential processes; method order rotates. Host load is recorded and timings are not isolated-machine measurements.","",
           "| Case | Method | Repeats | Solve seconds | Process wall seconds | Peak process MiB | Kernel dimension | Exact reference |",
           "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    for name in prepared:
        checked=any(c.get("exact_space")=="pass" for c in result["checks"] if c["case"]==name)
        for method in METHODS:
            rows=[r for r in result["runs"] if r["case"]==name and r["method"]==method and "output_sha256" in r]
            if not rows:lines.append(f"| {name} | {method} | 0 | — | — | — | — | failed | ");continue
            median=lambda key:statistics.median(r[key] for r in rows)
            lines.append(f"| {name} | {method} | {len(rows)} | {median('solve_seconds'):.6f} | {median('wall_seconds'):.3f} | {median('peak_process_rss_kib')/1024:.2f} | {int(median('dimension'))} | {'pass' if checked else 'unverified'} |")
    lines += ["", "Solve time includes interface construction, presolve, elimination and required certificates; process wall time also includes loading, input-seal verification, normalization, serialization and exact byte readback. GNU time reports peak single-process RSS; sampled process-tree RSS is separately retained in JSON. Independent checks construct the original rational sewing kernel and use explicit product actions, then compare complete rational row spaces for every output.", "", "See [results.json](results.json) for every run, preparation/verification costs, resource failures, input/output hashes, host metadata and source/binary snapshots."]
    if result["failures"]:lines += ["", "Failures:", *["- "+f for f in result["failures"]]]
    (out/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(out/"REPORT.md");return int(bool(result["failures"]))


if __name__ == "__main__":
    raise SystemExit(main())
