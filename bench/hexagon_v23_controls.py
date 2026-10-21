#!/usr/bin/env python3
"""Serialized nine-loop CPU-hour benchmark, retaining all fresh outputs."""
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
p.add_argument('--package',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--jobs',default='V17-p1,V19-p1,V19-p2,V17-p2,V19-p1-r2,V17-p1-r2')
p.add_argument('--resume',action='store_true')
p.add_argument('--threads',type=int,default=16)
p.add_argument('--loop',type=int,default=9,choices=[8,9])
p.add_argument('--replay',type=int,default=1,choices=[0,1])
p.add_argument('--factors',type=int,default=1,choices=[0,1])
a=p.parse_args();O=a.output.resolve();P=a.package.resolve()
O.mkdir(parents=True,exist_ok=a.resume)
F=R/'output_hexagon_16p4_20261006/FEC';L=R/'output_hexagon_16p4_20261006/LEC/MHV';D=R/'output_hexagon_16p4_20261006/data'
native=R/('output_hexagon_17p3_20261007/v17-nine-control/nine/joint' if a.loop==9 else 'output_hexagon_17p3_20261007/v17-schedule-controls/eight-t4/joint')
runtime=R/'output_hexagon_runtime_recovery_20261008/env/lib'
def sha(q):
    h=hashlib.sha256()
    with Path(q).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def producer_for(version):
    for base in (P,R/'output_hexagon_v19_20261009/package-v22',R/'output_hexagon_v19_20261009/package-v21',R/'output_hexagon_v19_20261009/package-v20',R/'output_hexagon_v19_20261009/package'):
        q=base/'bin'/('joint_physics_'+version.lower())
        if q.exists():return q
    raise FileNotFoundError(version)
state=dict(schema='hexagon-v23-cpu-controls-v1',loop=a.loop,replay=a.replay,factors=a.factors,status='starting',started=time.time(),threads=a.threads,
    acceptance=dict(baseline='original V17 with 2 GiB OPE cache',cpu_saving_strictly_greater_than=.15,
                    rss_ratio_max=1.05,metric='process user+system CPU seconds / 3600',
                    scope='complete fresh prime including all current independent physical and complete seam checks; warm precomputed basis inputs'),
    hard_cgroup_bytes=20000000000,high_cgroup_bytes=18000000000,swap_bytes=0,
    production_limit_bytes=50000000000,production_modified=False,steps=[],
    binaries={v:sha(producer_for(v)) for v in ['V17','V18','V19']+
              (['V20'] if (P/'bin/joint_physics_v20').exists() or (P/'bin/joint_physics_v21').exists() else [])+
              (['V21'] if (P/'bin/joint_physics_v21').exists() else [])+
              (['V22'] if (P/'bin/joint_physics_v22').exists() else [])})
for requested in a.jobs.split(','):
    version=requested.split('-')[0]
    state['binaries'][version]=sha(producer_for(version))
if a.resume:
    prior=json.loads((O/'progress.json').read_text());assert prior['binaries']==state['binaries'] and prior['threads']==a.threads
    state=prior;state['resumed']=time.time()
def save():
    q=O/'progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(O/'progress.json')
try:
    # Do not overlap numerical timings with the older independent backend tests.
    gate=json.loads((R/'output_hexagon_validation_20261009/replay-comparison.json').read_text())
    assert gate['runs']['mapped']['status']=='passed' and gate['runs']['heap']['status']=='passed'
    for label in a.jobs.split(','):
        old=[s for s in state['steps'] if s['label']==label]
        if old:
            assert old[-1]['status']=='passed','Failed outputs retained; inspect before retrying'
            continue
        m=re.fullmatch(r'(V17|V18|V19|V20|V21|V22|V23|V24|V25|V26)-p([12])(?:-r(\d+))?',label);assert m,label
        version,index,repetition=m.groups();index=int(index)
        prime=[1152921504606847009,1152921504606847067][index-1]
        root=O/(version+('-r'+repetition if repetition else ''))
        producer=producer_for(version)
        command=[producer,2*a.loop-3,3,F,L,D,root,6000 if a.loop==9 else 2048,'cache','right','--prime-index',index]
        cache=2048 if version=='V17' else 3072
        env=['env','-u','HEXAGON_MAPPED_FILE','-u','HEXAGON_MAPPED_GIB','-u','HEXAGON_MAPPED_SOFT_GIB',
             'THREADS='+str(a.threads),'PHYSICS_CACHE_MIB=4096','OPE_HEAD_CACHE_MIB='+str(cache),
             'HEXAGON_MHV_DOUBLE_BONUS=1','HEXAGON_DISCOVERY_REPLAY='+str(a.replay),'HEXAGON_AFFINE_FACTORS='+str(a.factors),'LD_LIBRARY_PATH='+str(runtime),'MIMALLOC_PURGE_DELAY=0']
        unit='symbology-v23-'+uuid.uuid4().hex[:10]
        cmd=['systemd-run','--user','--quiet','--wait','--collect','--unit='+unit,
             '--property=MemoryMax=20000000000','--property=MemoryHigh=18000000000',
             '--property=MemorySwapMax=0','--property=TasksMax=128','--property=WorkingDirectory='+str(R),
             sys.executable,str(R/'bench/measure_process.py'),'--report',str(O/(label+'.json')),
             '--log',str(O/(label+'.log')),'--status',str(O/(label+'.status.json')),
             '--timeout','2400','--virtual-gib','96','--rss-gib',str(19000000000/1024**3),
             '--min-available-gib','3','--',*env,*map(str,command)]
        step=dict(label=label,status='running',started=time.time(),command=cmd,ope_cache_mib=cache,prime=prime)
        state['steps'].append(step);state['status']='running_'+label;save();print('RUN',label,flush=True)
        cg_path=None;peak=0;events={}
        with (O/(label+'.manager.log')).open('x') as stream:
            proc=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                if cg_path is None:
                    q=subprocess.run(['systemctl','--user','show',unit,'--property=ControlGroup','--value'],capture_output=True,text=True)
                    if q.stdout.strip():cg_path=Path('/sys/fs/cgroup'+q.stdout.strip())
                if cg_path:
                    try:
                        peak=max(peak,int((cg_path/'memory.peak').read_text()))
                        events={k:int(v) for k,v in (r.split() for r in (cg_path/'memory.events').read_text().splitlines())}
                    except FileNotFoundError:pass
                time.sleep(.5)
        measurement=json.loads((O/(label+'.json')).read_text());assert proc.returncode==measurement['exit_code']==0
        log=(O/(label+'.log')).read_text()
        phases=[]
        for line in log.splitlines():
            if line.startswith('PROFILE '):
                phases.append({k:(v if k=='phase' else float(v)) for k,v in re.findall(r'(\w+)=([^ ]+)',line)})
        actual=root/'joint/cache'/f'p{prime}.bin';reference=native/'cache'/actual.name
        reference_scope='Original V17 native cache'
        if not reference.exists():
            reference=R/('output_hexagon_v18_20261009/validation/L'+str(a.loop)+'-V18/joint/cache')/actual.name
            reference_scope='Previously validated V18 second native prime; separate forced-QQ comparison required'
        x=actual.read_bytes();y=reference.read_bytes()
        assert x[:64]==y[:64] and x[72:]==y[72:],'Native coefficient difference'
        assert all('JOINT_PHYSICAL_CHECK operator='+name+' passed=1' in log for name in ('collinear','mrk','ope_small_s'))
        assert 'JOINT_PRIME_COMPLETE' in log
        cpu=measurement['cpu_user_seconds']+measurement['cpu_system_seconds']
        step.update(status='passed',finished=time.time(),measurement=measurement,cpu_seconds=cpu,cpu_hours=cpu/3600,
                    phase_profiles=phases,cgroup_peak_bytes=peak,cgroup_events=events,
                    all_native_coefficients_equal=True,cache_sha256=sha(actual),reference=str(reference),reference_scope=reference_scope,
                    reference_sha256=sha(reference),log_sha256=sha(O/(label+'.log')))
        save();print('PASS',label,'CPU_HOURS',cpu/3600,'RSS',measurement['peak_summed_rss_bytes'],flush=True)
    state['status']='passed';state['finished']=time.time();save()
except BaseException as exc:
    state['status']='stopped';state['error']=repr(exc);state['finished']=time.time();save();raise
