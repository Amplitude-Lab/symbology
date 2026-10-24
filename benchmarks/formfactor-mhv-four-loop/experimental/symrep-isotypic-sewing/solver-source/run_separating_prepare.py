"""One-core,5GB-capped setup experiment; does not extend symbol weights."""
from pathlib import Path
import json,subprocess,os,time,signal,hashlib
B=Path(__file__).resolve().parents[2];R=B.parents[1];O=B/'results/ff_four_loop'
old=next(x for x in reversed(json.loads((O/'stage_receipts.json').read_text())) if x.get('stage')=='symrep_cached_prepare')
cmd=list(old['command']);cmd[0]=str(B/'code/ff_four_loop/symrep_separating')
cmd[cmd.index('--output')+1]=str(O/'symrep_separating_prepared');cmd[cmd.index('--threads')+1]='1'
env=dict(os.environ,LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
t=time.monotonic();peak=0;samples=[]
with (O/'symrep_separating_prepare.log').open('w') as log:
 p=subprocess.Popen(['taskset','-c','3','/usr/bin/time','-f','{"user_seconds":%U,"system_seconds":%S,"maxrss_kib":%M}','-o',str(O/'symrep_separating_prepare.time'),*cmd],stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
 while p.poll() is None:
  pids=[p.pid]
  try:
   for pid in pids:
    cp=Path(f'/proc/{pid}/task/{pid}/children')
    if cp.exists():pids.extend(map(int,cp.read_text().split()))
   rss=sum(int(Path(f'/proc/{pid}/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE') for pid in pids)
   peak=max(peak,rss)
   if rss>5_000_000_000:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('5GB isolated setup cap exceeded')
   if not samples or time.monotonic()-t-samples[-1][0]>=5:samples.append([time.monotonic()-t,rss])
  except FileNotFoundError:pass
  (O/'symrep_separating_prepare.progress.json').write_text(json.dumps({'state':'running','pid':p.pid,'seconds':time.monotonic()-t,'peak_RSS_bytes':peak,'samples':samples})+'\n');time.sleep(.2)
receipt={'state':'PASS' if p.returncode==0 else 'FAILED','command':cmd,'exit_code':p.returncode,'wall_seconds':time.monotonic()-t,'peak_sampled_RSS_bytes':peak,'resources':(O/'symrep_separating_prepare.time').read_text(),'samples':samples,'threads':1,'affinity':[3],'memory_limit_bytes':5_000_000_000,'binary_sha256':hashlib.sha256(Path(cmd[0]).read_bytes()).hexdigest(),'scope':'Representation and seed/condition preparation only. Does not construct or certify a four-loop symbol.'}
(O/'symrep_separating_prepare.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));raise SystemExit(p.returncode)
