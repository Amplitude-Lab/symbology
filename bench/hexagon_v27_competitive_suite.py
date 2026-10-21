#!/usr/bin/env python3
"""Finish the running references, then build and test the reviewed revision."""
import json,os,signal,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v27_20261009';B=O/'suite-c';B.mkdir(exist_ok=False)
P=O/'package-v27c';state={'state':'waiting_for_scheduler_comparison','steps':[],'started':time.time()}
def save():
    q=B/'progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(B/'progress.json')
def run(label,cmd):
    step={'label':label,'command':list(map(str,cmd)),'started':time.time(),'status':'running'};state['steps'].append(step);state['state']=label;save()
    with (B/(label+'.log')).open('x') as f:r=subprocess.run(list(map(str,cmd)),cwd=R,stdout=f,stderr=subprocess.STDOUT)
    step.update(returncode=r.returncode,finished=time.time(),status='passed' if not r.returncode else 'failed');save();assert not r.returncode,label
save()
try:
    while True:
        previous=json.loads((O/'suite-b/progress.json').read_text());assert previous['state']!='stopped'
        if previous['state']=='passed':break
        time.sleep(3)
    run('build',['python3',R/'bench/hexagon_v27/build.py','--output',P])
    run('generic-new',[P/'bin/v27_generic_control'])
    run('generic-legacy',['env','SOLVER_WORK_KERNELS=0','SOLVER_WORK_SCHEDULE=0','SOLVER_EXACT_COMPRESSION=0','SOLVER_PACKED_ROUTES=0',P/'bin/v27_generic_control'])
    run('additional-controls',[P/'bin/v27_additional_controls'])
    for loop in (8,9):
        run(f'L{loop}-V27',['python3',R/'bench/hexagon_v27_controls.py','--package',P,'--output',B/f'L{loop}-V27','--loop',loop,'--jobs','V27-p1,V27-p2'])
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
