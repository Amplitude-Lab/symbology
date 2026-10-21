#!/usr/bin/env python3
"""Bounded, serial V17/V18 warm-basis controls at eight and nine loops."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--package',type=Path,required=True)
p.add_argument('--resume',action='store_true')
p.add_argument('--math-only',action='store_true')
p.add_argument('--tune-cache',action='store_true',help='After standard controls, test the measured three-head cache capacity')
a=p.parse_args();O=a.output.resolve();P=a.package.resolve()
O.mkdir(parents=True,exist_ok=a.resume)
F=R/'output_hexagon_16p4_20261006/FEC';L=R/'output_hexagon_16p4_20261006/LEC/MHV';D=R/'output_hexagon_16p4_20261006/data'
B=R/'output_hexagon_17p3_20261007';V=R/'output_hexagon_validation_20261009'
runtime=R/'output_hexagon_runtime_recovery_20261008/env/lib'
cases={8:(13,2048,B/'v17-schedule-controls/eight-t4'),9:(15,6000,B/'v17-nine-control/nine')}
known={8:R/'output_hexagon_12p4_20261006/modular_direct/SEW_12p4/MHV/collinear_mrk_ope_small_s',
       9:R/'output_hexagon_14p4_20261006/cached-nine/joint'}
prime=1152921504606847009;prime2=1152921504606847067
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
state={'schema':'hexagon-v18-validation-v1','status':'starting','started':time.time(),
       'threads':16,'hard_memory_bytes':20000000000,'memory_high_bytes':18000000000,
       'production_memory_limit_bytes':50000000000,'storage':'native heap, identical for both variants',
       'scope':'Warm-basis complete native prime and same-checkpoint final-seam comparison. Excludes cold basis/input generation. No production replacement.',
       'steps':[],'binary_sha256':{str(P/'bin'/n):sha(P/'bin'/n) for n in
           ['joint_physics_v17','joint_physics_v18','hexagon_v18_seam_control','pipeline_compare']}}
if a.resume:
    previous=json.loads((O/'progress.json').read_text());assert previous['binary_sha256']==state['binary_sha256']
    state=previous;state['resumed']=time.time()
def save():
    t=O/'progress.tmp';t.write_text(json.dumps(state,indent=2)+'\n');t.replace(O/'progress.json')
def run(label,command,timeout=1800,ope_cache_mib=2048):
    previous=[x for x in state['steps'] if x['label']==label]
    if previous and previous[-1]['status']=='passed':return
    assert not previous,'Failed work retained; inspect before a new attempt'
    unit='symbology-v18-'+uuid.uuid4().hex[:10]
    env=['env','-u','HEXAGON_MAPPED_FILE','-u','HEXAGON_MAPPED_GIB','-u','HEXAGON_MAPPED_SOFT_GIB',
         'THREADS=16','PHYSICS_CACHE_MIB=4096','OPE_HEAD_CACHE_MIB='+str(ope_cache_mib),'HEXAGON_MHV_DOUBLE_BONUS=1',
         'LD_LIBRARY_PATH='+str(runtime),'MIMALLOC_PURGE_DELAY=0']
    cmd=['systemd-run','--user','--quiet','--wait','--collect','--unit='+unit,
         '--property=MemoryMax=20000000000','--property=MemoryHigh=18000000000',
         '--property=MemorySwapMax=0','--property=TasksMax=128',
         '--property=WorkingDirectory='+str(R),sys.executable,str(R/'bench/measure_process.py'),
         '--report',str(O/(label+'.json')),'--log',str(O/(label+'.log')),
         '--status',str(O/(label+'.status.json')),'--timeout',str(timeout),
         '--virtual-gib','96','--rss-gib',str(19_000_000_000/1024**3),'--min-available-gib','3',
         '--',*env,*map(str,command)]
    step={'label':label,'status':'running','command':cmd};state['steps'].append(step)
    state['status']='running_'+label;save();print('RUN',label,flush=True)
    cg_peak=0;cg_events={};cg_path=None
    with (O/(label+'.manager.log')).open('x') as stream:
        proc=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            if cg_path is None:
                q=subprocess.run(['systemctl','--user','show',unit,'--property=ControlGroup','--value'],
                    capture_output=True,text=True)
                value=q.stdout.strip()
                if value:cg_path=Path('/sys/fs/cgroup'+value)
            if cg_path:
                try:
                    cg_peak=max(cg_peak,int((cg_path/'memory.peak').read_text()))
                    cg_events={k:int(v) for k,v in (line.split() for line in (cg_path/'memory.events').read_text().splitlines())}
                except FileNotFoundError:pass
            time.sleep(.5)
    report=json.loads((O/(label+'.json')).read_text())
    text=(O/(label+'.log')).read_text()
    step.update(status='passed' if proc.returncode==0 and report['exit_code']==0 else 'failed',
                measurement=report,log_sha256=sha(O/(label+'.log')),
                cgroup_peak_bytes=cg_peak,cgroup_memory_events=cg_events,
                phase_results=[json.loads(x.split(' ',1)[1]) for x in text.splitlines() if x.startswith('V18_CONTROL_PASS ')],
                seam_results=[x for x in text.splitlines() if x.startswith('V18_SEAM_COMPLETE ')])
    save()
    if step['status']!='passed':raise RuntimeError(label+' failed; all outputs retained')
    print('PASS',label,round(report['wall_seconds'],3),report['peak_summed_rss_bytes'],flush=True)

def gate():
    deadline=time.time()+7200;last=None
    while True:
        statuses={}
        for mode in ['mapped','heap']:
            f=V/('replay-'+mode)/'manifest.json'
            statuses[mode]=json.loads(f.read_text()).get('status') if f.exists() else 'not_started'
        if all(s=='passed' for s in statuses.values()):
            state['isolation_gate']={'time':time.time(),'backend_receipts':statuses,
                'receipt_hashes':{mode:sha(V/('replay-'+mode)/'manifest.json') for mode in statuses}}
            save();return
        if any(s=='failed' for s in statuses.values()):
            raise RuntimeError('Other backend benchmark failed; inspect before isolated timing')
        state['status']='waiting_for_backend_benchmarks';state['backend_status']=statuses;save()
        if statuses!=last:print('WAIT',statuses,flush=True);last=statuses
        if time.time()>deadline:raise RuntimeError('Backend isolation gate timed out; no timed calculation started')
        time.sleep(10)

try:
    if not a.math_only:gate()
    for loop,(h,target,ref) in cases.items():
        cp=ref/'early_ope_small_s';native=ref/'joint';plan=O/f'L{loop}-phase-plan'
        if a.math_only:
            run(f'L{loop}-contraction-audit',[P/'bin/hexagon_v18_seam_control',h,3,F,L,D,cp,native,prime,'audit',plan])
            continue
        common=[P/'bin/hexagon_v18_seam_control',h,3,F,L,D,cp,native,prime]
        # Same prime/checkpoint/binary, with independent original-row replay
        # outside the measured region in every case.
        run(f'L{loop}-phase-v17',common+['v17',plan])
        run(f'L{loop}-phase-v18-learn',common+['v18',plan])
        run(f'L{loop}-phase-v18-replay',common+['v18',plan])
        for version in ([17,18] if loop==8 else [18,17]):
            root=O/f'L{loop}-V{version}'
            run(f'L{loop}-V{version}-prime1',[P/'bin'/f'joint_physics_v{version}',h,3,F,L,D,root,target,'cache','right','--prime-index',1])
            path=root/'joint/cache'/f'p{prime}.bin';reference=native/'cache'/path.name
            actual,expected=path.read_bytes(),reference.read_bytes()
            assert actual[:64]==expected[:64] and actual[72:]==expected[72:],'Every native reference coefficient differs'
            equality={'loop':loop,'version':version,'prime':prime,'columns':(38989 if loop==8 else 105898),'payload_sha256':hashlib.sha256(actual[72:]).hexdigest()}
            if equality not in state.setdefault('native_equalities',[]):state['native_equalities'].append(equality)
            save()
        root=O/f'L{loop}-V18'
        run(f'L{loop}-V18-prime2',[P/'bin/joint_physics_v18',h,3,F,L,D,root,target,'cache','right','--prime-index',2])
        # Compare both primes at the common split against saved exact QQ data.
        shifted=O/f'L{loop}-V18-shifted'
        converter=R/'output_hexagon_16p4_20261006/joint-v8-package/bin/irr_convert'
        run(f'L{loop}-shift-both',[converter,'shift-cache','MHV',h,3,root/'joint',shifted,
            '--head',h-1,'--fec-root',F,'--lec-root',L,'--kind','unique'])
        for pp in [prime,prime2]:
            run(f'L{loop}-exact-reference-{pp}',[P/'bin/pipeline_compare',shifted,known[loop],h-1,4,F,L,pp])
    if a.tune_cache and not a.math_only:
        root=O/'L9-V18-cache3g'
        run('L9-V18-cache3g-prime1',[P/'bin/joint_physics_v18',15,3,F,L,D,root,6000,'cache','right','--prime-index',1],ope_cache_mib=3072)
        actual=(root/'joint/cache'/f'p{prime}.bin').read_bytes()
        expected=(cases[9][2]/'joint/cache'/f'p{prime}.bin').read_bytes()
        assert actual[:64]==expected[:64] and actual[72:]==expected[72:],'Tuned complete native reference differs'
        state['cache_tuning']={'budget_mib':3072,'every_native_coefficient_equal':True,
            'payload_sha256':hashlib.sha256(actual[72:]).hexdigest(),'hard_memory_bytes':20000000000}
    state['status']='mathematical_audits_passed' if a.math_only else 'passed'
    state['finished']=time.time();save()
except BaseException as exc:
    state['status']='stopped';state['error']=repr(exc);state['finished']=time.time();save();raise
