#!/usr/bin/env python3
"""Resolve worker efficiency, then ablate each implemented direction."""
import json,subprocess,time,statistics,os,signal
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v27_20261009';B=O/'final-comparisons';B.mkdir(exist_ok=True)
P=O/'package-v27d';state=json.loads((B/'progress.json').read_text())
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
try:
    while True:
        fewer=json.loads((B/'L9-eight-workers/progress.json').read_text());assert fewer['status']!='stopped'
        if fewer['status']=='passed':break
        time.sleep(3)
    receipt=json.loads((O/'final-comparisons-dispatch.json').read_text());pid=receipt['dispatcher_pid']
    assert Path(f'/proc/{pid}/cmdline').read_bytes().split(b'\0')[1]==b'bench/hexagon_v27_final_comparisons.py'
    os.kill(pid,signal.SIGTERM);os.kill(pid,signal.SIGCONT)
    state['steps'][-1].update(status='passed',returncode=0,finished=time.time())
    cpu8,rss8=metrics(fewer);cpu16,rss16=metrics(json.loads((O/'suite-d/L9-V27/progress.json').read_text()))
    baseline_cpu,baseline_rss=metrics(json.loads((O/'suite-a/L9-V26/progress.json').read_text()))
    accepted8=cpu8<.85*baseline_cpu and rss8<=1.05*baseline_rss
    threads=8 if accepted8 else 16
    state['selected_workers_for_ablations']=threads;state['worker_comparison']={'eight_cpu_hours':cpu8,'sixteen_cpu_hours':cpu16,'eight_rss_bytes':rss8,'sixteen_rss_bytes':rss16,'eight_passes_adoption_gate':accepted8};save()
    if accepted8:
        run('L8-eight-workers',loop=8,threads=8)
        run('L9-V26-eight-workers',threads=8,jobs='V26-p1,V26-p2')
    else:
        # A rejected resource setting needs a matched first-prime comparison,
        # not a separate adoption/cross-loop validation campaign.
        run('L9-V26-eight-workers',threads=8,jobs='V26-p1')
    for flag in ('kernels','compression','schedule','routes'):
        jobs='V27-p1,V27-p2' if flag=='routes' else 'V27-p1'
        run('L9-without-'+flag,threads=threads,jobs=jobs,extra=('--'+flag,'0'))
    state['state']='passed';state['finished']=time.time();save()
except BaseException as e:
    state['state']='stopped';state['error']=repr(e);state['finished']=time.time();save();raise
