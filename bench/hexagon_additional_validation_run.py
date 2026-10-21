#!/usr/bin/env python3
"""Run extra read-only validation in series, with receipts and bounded resources."""
import json,os,subprocess,time,hashlib,re
from pathlib import Path
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009';B=R/'output_hexagon_17p3_20261007';P=R/'output_hexagon_16p4_20261006'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
plan=json.loads((V/'physical-validation-plan.json').read_text());assert all(sha(p)==h for p,h in plan['inputs'].items())
env=os.environ.copy();env.update(THREADS='16',PHYSICS_CACHE_MIB='8192',OPE_HEAD_CACHE_MIB='16384',LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib'))
state={'schema':'hexagon-additional-validation-v1','started':time.time(),'plan_sha256':sha(V/'physical-validation-plan.json'),'attempts':[],'status':'running'}
def run(name,cmd,expect=0,timeout=2400):
 report=V/(name+'.resources.json');log=V/(name+'.log');assert not report.exists(),name
 state['stage']=name;save(V/'progress.json',state)
 measured=['python3','bench/measure_process.py','--report',str(report),'--log',str(log),'--status',str(V/'resource-status.json'),'--timeout',str(timeout),'--virtual-gib','128','--rss-gib','46.566128730773926','--min-available-gib','3','--',*map(str,cmd)]
 subprocess.run(measured,env=env,cwd=R,check=False)
 r=json.loads(report.read_text());a={'name':name,'report':str(report),'report_sha256':sha(report),'log':str(log),'log_sha256':sha(log),'measurement':r,'expected_exit_code':expect};state['attempts'].append(a);save(V/'progress.json',state)
 assert r['exit_code']==expect and not any(r[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal']),name
 return log.read_text()
try:
 # Wait only on the owned exact conversion; never overlap two large jobs.
 while not (V/'native-exact-shift.resources.json').exists():time.sleep(2)
 assert json.loads((V/'native-exact-shift.resources.json').read_text())['exit_code']==0
 run('native-independent-wxf',['wolframscript','-file','bench/hexagon_candidate_export_controls.wls','verify',V/'native-qq-17p3',V/'native-reconstruction-7p96.tsv'])
 run('exact-split-equality',['wolframscript','-file','bench/hexagon_candidate_export_controls.wls','verify',V/'native-exact-shift-16p4',B/'audit/crt-six-prime-speculative-candidate-20261008T1532Z.candidate.tsv'])
 run('native-exact-seam',[V/'hexagon_joint_certificate','17','3',V/'native-qq-17p3',P/'FEC',P/'LEC/MHV',P/'data'])
 run('native-exact-origin',[P/'joint-v8-package/bin/origin_check','17','3',V/'native-qq-17p3',P/'FEC',P/'LEC/MHV',R/'audits/hexagon-ten-loop-2026-10-06/origin-cases.tsv'])
 roots=[(12,4,R/'output_hexagon_12p4_20261006/modular_direct/SEW_12p4/MHV/collinear_mrk_ope_small_s'),(14,4,R/'output_hexagon_14p4_20261006/cached-nine/joint')]
 ext=P/'joint-v8-package/data/ope_L10_extended'
 def phys(name,h,k,root,data,prime,mode,samples=64,multiplier=1):
  return run(name,[V/'hexagon_physical_holdouts',h,k,root,P/'FEC',P/'LEC/MHV',data,V/name,prime,plan['seed'],samples,mode,multiplier],expect=2 if multiplier!=1 else 0)
 for h,k,root in roots:
  label='control-L'+str((h+k)//2)
  for mode in ['collinear','mrk']:phys(label+'-'+mode,h,k,root,P/'data',plan['primes'][0],mode,8)
  phys(label+'-extended-ope',h,k,root,ext,plan['primes'][0],'ope_small_s',8)
 phys('negative-control-scaled-symbol',12,4,roots[0][2],P/'data',plan['primes'][0],'mrk',4,2)
 for i,prime in enumerate(plan['primes'],1):
  for mode in ['collinear','mrk','ope_small_s']:
   phys('native-fresh'+str(i)+'-'+mode,17,3,V/'native-qq-17p3',ext if mode=='ope_small_s' else P/'data',prime,mode)
 state['status']='passed';state['finished']=time.time();state['physics_scope']=plan['holdout_scope'];state['conditional_uniqueness']=True
 save(V/'progress.json',state);save(V/'additional-verification.json',state)
except BaseException as e:
 state['status']='failed';state['error']=repr(e);save(V/'progress.json',state);raise
