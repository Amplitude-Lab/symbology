#!/usr/bin/env python3
"""Preserve a named calculation's log, exit status and measured process family."""
from pathlib import Path
import hashlib,json,os,signal,subprocess,sys,time

base=Path(__file__).resolve().parents[2]
out=base/'results/ff_four_loop'
name=sys.argv[1]; cap=int(sys.argv[2]); command=sys.argv[3:]
assert command and '/' not in name
logpath=out/(name+'.log');receiptpath=out/(name+'.receipt.json')
assert not logpath.exists() and not receiptpath.exists(), 'preserve prior run'
def family(pid):
    pids=[pid]
    for p in pids:
        try:pids.extend(map(int,Path(f'/proc/{p}/task/{p}/children').read_text().split()))
        except FileNotFoundError:pass
    return pids
def rss(pids):
    total=0
    for p in set(pids):
        try:total+=int(Path(f'/proc/{p}/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
        except FileNotFoundError:pass
    return total
start=time.monotonic();peak=0;stop=None;samples=[]
with logpath.open('w') as log:
    process=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/(name+'.resources.txt')),*command],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    while process.poll() is None:
        own=family(process.pid); use=rss(own);peak=max(peak,use);elapsed=time.monotonic()-start
        other=[]
        for path in [out/'progress.json',out/'boundary/DC/progress.json',out/'boundary/TC/progress.json']:
            try:
                state=json.loads(path.read_text())
                if state.get('state')=='running' and state.get('pid'):other.extend(family(state['pid']))
            except (FileNotFoundError,json.JSONDecodeError):pass
        if use>cap:stop=f'own{cap}bytecap'
        if rss(own+other)>45_000_000_000:stop='45GB combined recorded process guard, reserve5GB'
        if stop:os.killpg(process.pid,signal.SIGTERM);process.wait();break
        if not samples or elapsed-samples[-1][0]>10:samples.append([round(elapsed,2),use])
        time.sleep(.2)
receipt={'command':command,'returncode':process.returncode,'wall_seconds':time.monotonic()-start,'sampled_peak_rss_bytes':peak,'own_memory_cap_bytes':cap,'shared_accounted_cap_bytes':45_000_000_000,'stop_reason':stop,'samples':samples,'binary_sha256':hashlib.sha256(Path(command[0]).read_bytes()).hexdigest()}
receiptpath.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k]for k in ('returncode','wall_seconds','sampled_peak_rss_bytes','stop_reason')}))
raise SystemExit(process.returncode!=0)
