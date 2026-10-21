#!/usr/bin/env python3
"""Resolve worker efficiency, then ablate each implemented direction."""
import json,subprocess,time,statistics
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v27_20261009';B=O/'final-comparisons';B.mkdir(exist_ok=False)
P=O/'package-v27d';state={'state':'waiting_for_wide_controls','steps':[],'started':time.time()}
def save():
    q=B/'progress.tmp';q.write_text(json.dumps(state,indent=2)+'\n');q.replace(B/'progress.json')
def run(label,loop=9,threads=16,jobs='V27-p1,V27-p2',extra=()):
    cmd=['python3',str(R/'bench/hexagon_v27_controls.py'),'--package',str(P),'--output',str(B/label),'--loop',str(loop),'--threads',str(threads),'--jobs',jobs,*extra]
    step={'label':label,'command':cmd,'started':time.time(),'status':'running'};state['steps'].append(step);state['state']=label;save()
    with (B/(label+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT)
    step.update(returncode=r.returncode,finished=time.time(),status='passed' if not r.returncode else 'failed');save();assert not r.returncode,label
    return json.loads((B/label/'progress.json').read_text())
def metrics(d):
    return statistics.mean(s['cpu_hours'] for s in d['steps']),max(max(s['measurement']['peak_summed_rss_bytes'],max(x.get('process_peak_bytes',0) for x in s['phase_profiles'])) for s in d['steps'])
save()
try:
    while True:
        previous=json.loads((O/'suite-d/progress.json').read_text());assert previous['state']!='stopped','Combined candidate failed'
        if previous['state']=='passed':break
        time.sleep(3)
    fewer=run('L9-eight-workers',threads=8)
    cpu8,rss8=metrics(fewer);cpu16,rss16=metrics(json.loads((O/'suite-d/L9-V27/progress.json').read_text()))
    _,baseline_rss=metrics(json.loads((O/'suite-a/L9-V26/progress.json').read_text()))
    threads=8 if cpu8<cpu16 and rss8<=1.05*baseline_rss else 16
    state['selected_workers']=threads;state['worker_comparison']={'eight_cpu_hours':cpu8,'sixteen_cpu_hours':cpu16,'eight_rss_bytes':rss8,'sixteen_rss_bytes':rss16};save()
    if threads==8:
        run('L8-eight-workers',loop=8,threads=8)
        run('L9-V26-eight-workers',threads=8,jobs='V26-p1,V26-p2')
    for flag in ('kernels','compression','schedule','routes'):
        jobs='V27-p1,V27-p2' if flag=='routes' else 'V27-p1'
        run('L9-without-'+flag,threads=threads,jobs=jobs,extra=('--'+flag,'0'))
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
