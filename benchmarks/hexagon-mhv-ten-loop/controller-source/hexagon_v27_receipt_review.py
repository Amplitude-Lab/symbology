#!/usr/bin/env python3
"""Review archived V27 controls without launching a numerical calculation."""
from pathlib import Path
import hashlib,json,math,re,shutil,statistics,time

R=Path(__file__).resolve().parents[1]
N=R/'output_hexagon_v27_20261009'
V=R/'output_hexagon_validation_20261009'
D=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v27-experimental'
D.mkdir(exist_ok=True)
bindings={}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def bind(p,expected=None):
    p=Path(p).resolve();h=sha(p)
    assert expected is None or h==expected,('hash',str(p))
    bindings[str(p)]=h
    return h
def keep(p,name=None):
    p=Path(p);bind(p);target=D/(name or p.name)
    if target.exists():assert sha(target)==sha(p),('archive collision',str(target))
    else:shutil.copy2(p,target)
def same(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)

x=json.loads((N/'RESULTS.json').read_text())
c=json.loads((N/'completion-receipt.json').read_text())
profile=json.loads((N/'experimental-profile.json').read_text())
gate=json.loads((N/'execution-gate.json').read_text())
assert x['status']==c['state']=='complete'
assert x['selected']==c['selected']==gate['selected']=='V26'
assert not x['candidate_accepted'] and not c['candidate_accepted'] and not gate['acceptance']['passed']
assert not c['production_modified'] and c['production_memory_limit_bytes']==50000000000
for p,h in c['input_hashes'].items():bind(p,h)
source_count=0
for item in c['builds']:
    receipt=Path(item['receipt']);bind(receipt,item['receipt_sha256'])
    build=json.loads(receipt.read_text());assert len(build['hashes'])==item['files_verified']
    for p,h in build['hashes'].items():bind(receipt.parent/p,h)
    source_count+=len(build['hashes']);keep(receipt,receipt.parent.name+'-'+receipt.name)
bind(profile['binary'],profile['binary_sha256'])
for plan in c['plans']:
    p=Path(plan['path']);bind(p,plan['sha256']);assert p.stat().st_size==plan['bytes']

intervals=[]
for z in x['runs']:
    log=Path(z['log']);receipt=log.parent/'progress.json'
    state=json.loads(receipt.read_text());a=next(a for a in state['steps'] if a['label']==z['label'])
    m=json.loads(Path(z['measurement']).read_text())
    assert a['status']=='passed' and m==a['measurement']
    assert m['exit_code']==0 and not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
    same((m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,z['cpu_hours'])
    assert m['wall_seconds']==z['wall_seconds'] and m['peak_summed_rss_bytes']==z['sampled_peak_bytes']
    assert a['phase_profiles']==z['phases']
    kernel=max(q.get('process_peak_bytes',0) for q in z['phases'])
    assert z['peak_rss_bytes']==max(kernel,z['sampled_peak_bytes'])
    assert all(n==0 for n in z['cgroup_events'].values()) and z['cgroup_events']==a['cgroup_events']
    assert state['hard_cgroup_bytes']==20000000000 and state['high_cgroup_bytes']==18000000000 and state['swap_bytes']==0
    bind(log,a['log_sha256']);binary=Path(a['command'][-12])
    assert binary.name=='joint_physics_'+z['version'].lower()
    assert bind(binary,state['binaries'][z['version']])==z['binary_sha256']
    assert state['threads']==z['threads'] and f"THREADS={z['threads']}" in a['command']
    assert state['flags']==z['flags']
    for flag,env in [('kernels','SOLVER_WORK_KERNELS'),('schedule','SOLVER_WORK_SCHEDULE'),('compression','SOLVER_EXACT_COMPRESSION'),('routes','SOLVER_PACKED_ROUTES')]:
        assert f"{env}={z['flags'][flag]}" in a['command']
    root=Path(a['command'][-6]);cache=root/'joint/cache'/f"p{z['prime']}.bin"
    bind(cache,z['cache_sha256']);bind(z['reference'],z['reference_sha256'])
    actual=cache.read_bytes();ref=Path(z['reference']).read_bytes()
    assert actual[:64]==ref[:64] and actual[72:]==ref[72:]
    t=log.read_text()
    for op in ['collinear','mrk','ope_small_s']:assert 'JOINT_PHYSICAL_CHECK operator='+op+' passed=1' in t
    assert 'V18_SEAM_VERIFIED full_operator=1' in t
    assert z['native_coefficients_equal'] and z['all_physical_checks'] and z['complete_seam_check']
    parsed=[]
    for line in t.splitlines():
        if line.startswith('PROFILE phase='):
            q=dict(re.findall(r'(\w+)=([^ ]+)',line))
            parsed.append({k:(value if k=='phase' else float(value)) for k,value in q.items()})
    assert parsed==z['phases'],('phase parser',str(log))
    intervals.append((a['started'],a['finished'],z['suite'],z['label']))
    prefix=z['suite'].replace('/','-')+'-'
    for p in [receipt,log,Path(z['measurement'])]:keep(p,prefix+p.name)
for previous,current in zip(sorted(intervals),sorted(intervals)[1:]):
    assert previous[1]<=current[0],('overlapping numerical jobs',previous,current)
for g in x['groups']:
    runs=[z for z in x['runs'] if z['suite']==g['suite'] and z['version']==g['version']]
    assert len(runs)==g['runs']
    same(sum(z['cpu_hours'] for z in runs),g['cpu_hours_total'])
    same(statistics.mean(z['cpu_hours'] for z in runs),g['cpu_hours_mean'])
    same(statistics.mean(z['wall_seconds'] for z in runs),g['wall_seconds_mean'])
    assert max(z['peak_rss_bytes'] for z in runs)==g['peak_rss_bytes']
assert len(x['runs'])==c['total_fresh_native_runs']==27
same(sum(z['cpu_hours'] for z in x['runs']),c['total_benchmark_cpu_hours'])
groups={g['suite']:g for g in x['groups']}
for key,suite in [('baseline','suite-a/L9-V26'),('nine_loop','suite-d/L9-V27'),('eight_loop','suite-d/L8-V27'),('nine_loop_eight_workers','final-comparisons/L9-eight-workers')]:
    assert groups[suite]==profile[key]
for label,new,old in [('comparison_nine','suite-d/L9-V27','suite-a/L9-V26'),('comparison_nine_eight_workers','final-comparisons/L9-eight-workers','suite-a/L9-V26'),('comparison_eight','suite-d/L8-V27','suite-a/L8-V26')]:
    saving=1-groups[new]['cpu_hours_mean']/groups[old]['cpu_hours_mean']
    ratio=groups[new]['peak_rss_bytes']/groups[old]['peak_rss_bytes']
    same(saving,profile[label]['cpu_saving']);same(ratio,profile[label]['rss_ratio'])
    assert (saving>0.15 and ratio<=1.05)==profile[label]['passed']
p1=1152921504606847009
full=next(z for z in x['runs'] if z['suite']=='suite-d/L9-V27' and z['prime']==p1)
for ablation in profile['ablations']:
    off=next(z for z in x['runs'] if z['suite']==ablation['off_suite'] and z['prime']==p1)
    same(full['cpu_hours'],ablation['full_cpu_hours']);same(off['cpu_hours'],ablation['off_cpu_hours'])
    same(1-full['cpu_hours']/off['cpu_hours'],ablation['enabling_feature_cpu_saving'])
matched=profile['matched_eight_worker_first_prime']
base8=next(z for z in x['runs'] if z['suite']=='final-comparisons/L9-V26-eight-workers')
new8=next(z for z in x['runs'] if z['suite']=='final-comparisons/L9-eight-workers' and z['prime']==p1)
same(base8['cpu_hours'],matched['baseline_cpu_hours']);same(new8['cpu_hours'],matched['candidate_cpu_hours'])
same(1-new8['cpu_hours']/base8['cpu_hours'],matched['cpu_saving'])
same(1-groups['suite-d/L9-V27']['cpu_hours_mean']/groups['final-comparisons/L9-without-routes']['cpu_hours_mean'],profile['route_pair']['enabling_packed_routes_cpu_saving'])

for name in ['generic-new.log','generic-legacy.log']:
    p=N/'suite-d'/name;t=p.read_text()
    assert 'GENERIC_CONTROL passed=1 field_cases=96 successful_replays=48 guarded_fallbacks=2 factor_cases=72' in t
    keep(p)
p=N/'suite-d/additional-controls.log';t=p.read_text()
assert 'WIDE_CONTROL overflow_cases=5077 canonical_reference_equal=1' in t
assert 'ADDITIONAL_CONTROL passed=1 queue_operations=756718 compression_cases=144 deterministic_schedules=64 packed_budget_and_damage_controls=1 exact_density_and_persistent_cycle_controls=1' in t
keep(p);keep(N/'suite-d/progress.json','suite-d-controls-progress.json')
# Recheck the prior reference-chain review, including genuine QQ-only comparisons.
prior=V/'v26-generic-review.json';q=json.loads(prior.read_text());assert q['passed']
for p,h in q['input_hashes'].items():bind(p,h)
bind(prior)
for name in ['RESULTS.json','REPORT.txt','experimental-profile.json','completion-receipt.json','execution-gate.json','reference-binding-audit.json','v26-to-v27d-review.patch']:keep(N/name)
bind(__file__)
out={'schema':'hexagon-v27-experimental-receipt-review-v1','created':time.time(),'passed':True,'native_runs_reviewed':27,'source_build_files_reviewed':source_count,'native_payload_equalities_rechecked':True,'CPU_means_conservative_RSS_ablations_and_rejection_recomputed':True,'serialized_native_intervals_verified':True,'prior_QQ_reference_chain_verified':True,'selected_nine_loop_profile':'V26','V27_accepted':False,'experimental_profile':profile,'scope':'Archived eight/nine-loop controls only. Learning/replay pairs are not timing repetitions. No new numerical solve or ten-loop performance inference.','input_hashes':bindings}
p=V/'v27-experimental-review.json';p.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(p,D/p.name)
print(json.dumps({k:v for k,v in out.items() if k not in ['input_hashes','experimental_profile']},indent=2));print('bindings',len(bindings))
