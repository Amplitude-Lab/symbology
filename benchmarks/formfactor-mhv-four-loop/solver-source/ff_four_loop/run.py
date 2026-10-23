#!/usr/bin/env python3
"""Measured, checkpoint-preserving native stages for the four-loop experiment."""
from pathlib import Path
import subprocess,os,json,time,sys,signal
BASE=Path(__file__).resolve().parents[2];ROOT=BASE.parents[1];OUT=BASE/'results/ff_four_loop';ENV=dict(os.environ,LD_LIBRARY_PATH=str(ROOT/'output_hexagon_runtime_recovery_20261008/env/lib'),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1')
def run(name,args,executable="bootstrap"):
 dest=OUT/(name+'.resources.json');cmd=[str(ROOT/executable),*map(str,args),'--threads','4'];t=time.monotonic();peak=0;samples=[]
 with(OUT/(name+'.log')).open('w')as f:
  p=subprocess.Popen(['/usr/bin/time','-f','{"user_seconds":%U,"system_seconds":%S,"maxrss_kib":%M}','-o',str(dest),*cmd],env=ENV,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while p.poll()is None:
   family=[p.pid]
   try:
    for pid in family:
     c=Path(f'/proc/{pid}/task/{pid}/children')
     if c.exists():family.extend(map(int,c.read_text().split()))
    rss=sum(int(Path(f'/proc/{pid}/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')for pid in family);peak=max(peak,rss)
    if rss>44_000_000_000:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('44GB own-process cap')
    if len(samples)==0 or time.monotonic()-t-samples[-1][0]>=10:samples.append([round(time.monotonic()-t,2),rss])
   except FileNotFoundError:pass
   (OUT/'progress.json').write_text(json.dumps({'state':'running','stage':name,'seconds':time.monotonic()-t,'rss_peak_bytes':peak,'samples':samples,'pid':p.pid},indent=2)+'\n');time.sleep(.2)
 receipt={'stage':name,'command':cmd,'exit_code':p.returncode,'wall_seconds':time.monotonic()-t,'peak_sampled_rss_bytes':peak,'samples':samples,'resources':json.loads(dest.read_text().splitlines()[-1]) if dest.exists() and dest.read_text().strip() else None};path=OUT/'stage_receipts.json';rs=json.loads(path.read_text())if path.exists()else[];rs.append(receipt);path.write_text(json.dumps(rs,indent=2)+'\n');(OUT/'progress.json').write_text(json.dumps({'state':'completed'if p.returncode==0 else'failed',**receipt},indent=2)+'\n');print(receipt,flush=True)
 if p.returncode:raise RuntimeError(name)
if __name__=='__main__':
 if sys.argv[1]=='prefix':
  for w in range(3,int(sys.argv[2])+1):
   if (OUT/f'FEC_{w}.wxf').exists():continue
   run(f'prefix_w{w}',['--extend','-c',OUT/'data/interior.wxf','-f',OUT/f'FEC_{w-1}.wxf','-o',OUT/f'FEC_{w}.wxf'])
 elif sys.argv[1]=='actions':
  w=int(sys.argv[2])
  for s in ['cyclic','flip']:run(f'{s}_w{w}',['--project','--symmetry',s,'--target',f'FEC_{w}','--data-dir',OUT/'data','--output-dir',OUT])
 elif sys.argv[1]=='sew':
  w=int(sys.argv[2]);cmd=['--sew','-c',OUT/'data/interior.wxf','-f',OUT/f'FEC_{w-1}.wxf','-l',OUT/'data/LEC_1.wxf','-o',OUT/f'SEW_{w-1}p1_invariant.wxf']
  for k,s in enumerate(['galois1a','galois1b','galois2a','galois2b','galois3'],1):cmd+=['--left-action',OUT/f'galois{k}_w{w-1}.wxf','--right-action',OUT/'data'/f'{s}.wxf']
  for s,t in [('cyclic','cycmat'),('flip','flipmat')]:cmd+=['--left-action',OUT/s/f'first_w{w-1}.wxf','--right-action',OUT/'data'/f'{t}.wxf']
  run(f'sew_w{w}',cmd)
 else:raise ValueError(sys.argv)
