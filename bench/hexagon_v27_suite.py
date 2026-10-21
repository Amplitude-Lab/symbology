#!/usr/bin/env python3
"""Serialized fresh V26/candidate comparisons, preserving every result."""
import argparse,json,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();O=a.output.resolve();P=a.package.resolve()
assert (P/'v27-build-receipt.json').exists();O.mkdir(exist_ok=False)
state={'state':'starting','started':time.time(),'steps':[],'baseline':'V26'}
def save():
    q=O/'progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(O/'progress.json')
def run(label,cmd):
    step={'label':label,'command':list(map(str,cmd)),'started':time.time(),'status':'running'};state['steps'].append(step);state['state']=label;save()
    with (O/(label+'.log')).open('x') as f:r=subprocess.run(list(map(str,cmd)),cwd=R,stdout=f,stderr=subprocess.STDOUT)
    step.update(returncode=r.returncode,finished=time.time(),status='passed' if not r.returncode else 'failed');save();assert not r.returncode,label
try:
    run('generic-new',[P/'bin/v27_generic_control'])
    run('generic-legacy',['env','SOLVER_WORK_KERNELS=0','SOLVER_WORK_SCHEDULE=0','SOLVER_EXACT_COMPRESSION=0','SOLVER_PACKED_ROUTES=0',P/'bin/v27_generic_control'])
    run('additional-controls',[P/'bin/v27_additional_controls'])
    for loop in (8,9):
        for version in ('V26','V27'):
            label=f'L{loop}-{version}'
            run(label,['python3',R/'bench/hexagon_v27_controls.py','--package',P,'--output',O/label,'--loop',loop,'--jobs',version+'-p1,'+version+'-p2'])
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
