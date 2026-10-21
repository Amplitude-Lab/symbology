#!/usr/bin/env python3
import json,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v23_20261009'
state={'state':'waiting_for_v24_timings','started':time.time(),'steps':[]}
def save():
    q=O/'v25-suite-progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(O/'v25-suite-progress.json')
def run(label,cmd):
    state['state']=label;step={'label':label,'command':cmd,'started':time.time()};state['steps'].append(step);save()
    with (O/(label+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT)
    step.update(returncode=r.returncode,finished=time.time());save();assert r.returncode==0,label
save()
try:
    while True:
        previous=json.loads((O/'v24-suite-progress.json').read_text());assert previous['state']!='stopped','Earlier suite failed; inspect first'
        if previous['state']=='passed':break
        time.sleep(5)
    run('v24-affine-families',[str(O/'package-v24/bin/v24_affine_benchmark')])
    P=O/'package-v25'
    run('v25-build',['python3','bench/hexagon_joint_build_v25.py','--output',str(P)])
    run('v25-L8',['python3','bench/hexagon_v23_controls.py','--package',str(P),'--output',str(O/'L8-v25'),'--loop','8','--jobs','V25-p1,V25-p2'])
    run('v25-L9',['python3','bench/hexagon_v23_controls.py','--package',str(P),'--output',str(O/'L9-v25'),'--jobs','V25-p1,V25-p2'])
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
