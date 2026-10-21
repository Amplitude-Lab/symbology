#!/usr/bin/env python3
"""Nine-loop switch ablations; replay-cache ablation includes both primes."""
import json,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v27_20261009';B=O/'ablations';B.mkdir(exist_ok=False)
P=O/'package-v27b';state={'state':'waiting_for_combined_controls','steps':[],'started':time.time()}
def save():
    q=B/'progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(B/'progress.json')
save()
try:
    while True:
        previous=json.loads((O/'suite-b/progress.json').read_text());assert previous['state']!='stopped','Combined candidate failed'
        if previous['state']=='passed':break
        time.sleep(3)
    for flag in ('compression','schedule','routes','kernels'):
        label='L9-without-'+flag
        jobs='V27-p1,V27-p2' if flag=='routes' else 'V27-p1'
        cmd=['python3',str(R/'bench/hexagon_v27_controls.py'),'--package',str(P),'--output',str(B/label),'--loop','9','--jobs',jobs,'--'+flag,'0']
        step={'label':label,'command':cmd,'started':time.time(),'status':'running'};state['steps'].append(step);state['state']=label;save()
        with (B/(label+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT)
        step.update(returncode=r.returncode,finished=time.time(),status='passed' if not r.returncode else 'failed');save();assert not r.returncode,label
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
