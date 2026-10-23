#!/usr/bin/env python3
"""One-core, 6 GB isolated holdout test with a shared 45 GB guard."""
import hashlib,json,os,pathlib,signal,subprocess,sys,time
base=pathlib.Path(__file__).resolve().parents[1]
exp=base/'results/ff_four_loop'
boundary=sys.argv[1]
assert boundary in ('DC','TC')
stem=exp/'holdouts'/(boundary+'_prefix_probe')
binary=base/'code/ff_four_loop/boundary_holdout_probe'
command=['/usr/bin/time','-v','taskset','-c','3','nice','-n','10',str(binary),str(exp),boundary]
log=open(str(stem)+'.log','w');resources=open(str(stem)+'.resources.txt','w')
start=time.time();p=subprocess.Popen(command,stdout=log,stderr=resources,start_new_session=True)
peak=sharedpeak=0;stop=None;samples=[]
while p.poll() is None:
 own=shared=0
 for f in pathlib.Path('/proc').glob('[0-9]*/status'):
  try:
   pid=int(f.parent.name);s=f.read_text();rss=int(next(x for x in s.splitlines() if x.startswith('VmRSS:')).split()[1])*1024
   if os.getpgid(pid)==p.pid:own+=rss
   cmd=(f.parent/'cmdline').read_bytes().decode(errors='replace')
   if ('research/4pformfactor' in cmd or '\0bootstrap\0' in cmd) and pid!=os.getpid():shared+=rss
  except (OSError,ProcessLookupError,StopIteration):pass
 peak=max(peak,own);sharedpeak=max(sharedpeak,shared)
 if not samples or time.time()-start-samples[-1][0]>5:samples.append([round(time.time()-start,3),own,shared])
 if own>6_000_000_000:stop='own process group exceeds 6 GB decimal'
 if shared>45_000_000_000:stop='shared calculation exceeds 45 GB decimal'
 if stop:os.killpg(p.pid,signal.SIGTERM);break
 time.sleep(.2)
p.wait();log.close();resources.close()
inputs=[binary,base/'code/ff_four_loop/boundary_holdout_probe.cpp',base/'code/ff_four_loop/boundary_tensor_train.hpp']
receipt={'command':command,'returncode':p.returncode,'stop_reason':stop,'wall_seconds':time.time()-start,'peak_RSS_bytes':peak,'peak_shared_RSS_bytes':sharedpeak,'own_cap_bytes':6_000_000_000,'shared_cap_bytes':45_000_000_000,'workers':1,'affinity':[3],'finite_targets_used_only_as_holdouts':True,'sha256':{str(x.relative_to(base)):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},'samples':samples}
pathlib.Path(str(stem)+'.guard.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='samples'},indent=2));sys.exit(p.returncode)
