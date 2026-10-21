#!/usr/bin/env python3
"""Summarize complete, fresh-prime measurements without hiding learning cost."""
import hashlib,json,re,statistics
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v23_20261009'
baseline=json.loads((R/'output_hexagon_v19_20261009/RESULTS.json').read_text())
# Original V17 two matched first-prime repeats, already frozen and reviewed.
v17_cpu=1.4816453901388888;v17_rss=5903253504
records=[]
for directory in sorted(O.iterdir()):
    p=directory/'progress.json'
    if not directory.is_dir() or not p.exists():continue
    state=json.loads(p.read_text())
    for step in state.get('steps',[]):
        if step.get('status')!='passed':continue
        log=directory/(step['label']+'.log');text=log.read_text();m=step['measurement']
        def lines(prefix):
            return [{k:(float(v) if re.fullmatch(r'[-+0-9.eE]+',v) else v) for k,v in re.findall(r'(\w+)=([^ ]+)',line)} for line in text.splitlines() if line.startswith(prefix+' ')]
        records.append({'suite':directory.name,'loop':state['loop'],'version':step['label'].split('-')[0],'label':step['label'],
            'replay_enabled':state['replay'],'factors_enabled':state['factors'],'prime':step['prime'],
            'cpu_hours':step['cpu_hours'],'wall_seconds':m['wall_seconds'],'sampled_peak_rss_bytes':m['peak_summed_rss_bytes'],'kernel_peak_rss_bytes':max([x.get('process_peak_bytes',0) for x in step['phase_profiles']],default=0),'peak_rss_bytes':max(m['peak_summed_rss_bytes'],max([x.get('process_peak_bytes',0) for x in step['phase_profiles']],default=0)),
            'cgroup_peak_bytes':step['cgroup_peak_bytes'],'cgroup_events':step['cgroup_events'],
            'phases':step['phase_profiles'],'structural_profile':lines('STRUCTURAL_PROFILE'),'dependency_profile':lines('DEPENDENCY_PROFILE'),
            'discovery_structure':lines('DISCOVERY_STRUCTURE'),'fallbacks':[s for s in text.splitlines() if 'STRUCTURE_FALLBACK' in s],
            'reference_coefficients_equal':step['all_native_coefficients_equal'],'all_physical_checks':all('JOINT_PHYSICAL_CHECK operator='+name+' passed=1' in text for name in ['collinear','mrk','ope_small_s']),
            'full_seam_check':'V18_SEAM_VERIFIED full_operator=1' in text,'log':str(log.resolve()),'measurement':str((directory/(step['label']+'.json')).resolve()),
            'reference':step['reference'],'reference_sha256':step['reference_sha256'],'cache_sha256':step['cache_sha256']})
groups=[]
for key in sorted({(x['suite'],x['version']) for x in records}):
    xs=[x for x in records if (x['suite'],x['version'])==key];peak=max(x['peak_rss_bytes'] for x in xs);mean=statistics.mean(x['cpu_hours'] for x in xs)
    g={'suite':key[0],'version':key[1],'loop':xs[0]['loop'],'runs':len(xs),'primes':[x['prime'] for x in xs],
       'cpu_hours_total':sum(x['cpu_hours'] for x in xs),'cpu_hours_mean':mean,'wall_seconds_mean':statistics.mean(x['wall_seconds'] for x in xs),
       'peak_rss_bytes':peak,'sampled_peak_rss_bytes':max(x['sampled_peak_rss_bytes'] for x in xs),'kernel_peak_rss_bytes':max(x['kernel_peak_rss_bytes'] for x in xs),'correctness_passed':all(x['reference_coefficients_equal'] and x['all_physical_checks'] and x['full_seam_check'] for x in xs)}
    if g['loop']==9:g.update(cpu_saving_vs_original_v17=1-mean/v17_cpu,rss_ratio_vs_conservative_original_v17=peak/v17_rss,
        original_v17_numeric_gate=(mean<.85*v17_cpu and peak<=1.05*v17_rss))
    groups.append(g)
result={'status':'interim','scope':'Fresh native primes on warm basis/RHS inputs, including all physical and complete original seam checks; CPU=user+system. Excludes compilation/cold input generation/new QQ reconstruction. First prime includes structural learning and plan serialization; subsequent prime includes plan read and validation. Two-prime groups are learning+replay amortization, not same-prime repetitions.',
        'baseline':{'version':'original V17','cpu_hours_mean':v17_cpu,'conservative_peak_rss_bytes':v17_rss,'cpu_saving_required_strictly_greater_than':.15,'rss_ratio_max':1.05},'groups':groups,'runs':records,'rss_metric':'Conservative maximum of 100ms sampled process RSS and kernel process high-water mark in phase logs; one process per native solve. Historical V17 reference uses its sampled minimum, which is conservative for the new gate.','production_modified':False,'production_memory_limit_bytes':50000000000}
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
for g in groups:print(g['suite'],g['version'],g['runs'],f"{g['cpu_hours_mean']:.6f} CPUh",f"{g['peak_rss_bytes']/1e9:.3f} GB")
