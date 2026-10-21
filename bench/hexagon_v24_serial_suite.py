#!/usr/bin/env python3
"""Do not compile or start new numerical work during the existing timing suite."""
import json,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v23_20261009'
state={'state':'waiting_for_v23_and_v22_timings','started':time.time(),'steps':[]}
def save():
    q=O/'v24-suite-progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(O/'v24-suite-progress.json')
def run(label,cmd):
    state['state']=label;step={'label':label,'command':cmd,'started':time.time()};state['steps'].append(step);save()
    with (O/(label+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT)
    step.update(returncode=r.returncode,finished=time.time());save();assert r.returncode==0,label
save()
try:
    while True:
        prior=json.loads((O/'L9-combined/progress.json').read_text())
        assert prior['status']!='stopped','Earlier suite failed; inspect first'
        if prior['status']=='passed':break
        time.sleep(5)
    P=O/'package-v24'
    run('v24-build',['python3','bench/hexagon_joint_build_v24.py','--output',str(P)])
    run('v24-generic',[str(P/'bin/v24_generic_control')])
    run('v24-affine-families',[str(P/'bin/v24_affine_benchmark')])
    run('v24-L8',['python3','bench/hexagon_v23_controls.py','--package',str(P),'--output',str(O/'L8-v24'),'--loop','8','--jobs','V24-p1,V24-p2,V22-p1,V22-p2'])
    run('v24-L9',['python3','bench/hexagon_v23_controls.py','--package',str(P),'--output',str(O/'L9-v24'),'--jobs','V24-p1,V24-p2'])
    run('v23-L9-factors',['python3','bench/hexagon_v23_controls.py','--package',str(O/'package-v23b'),'--output',str(O/'L9-factors-only'),'--replay','0','--factors','1','--jobs','V23-p1'])
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
