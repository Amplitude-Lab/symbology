#!/usr/bin/env python3
import json,re,statistics,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v27_20261009'
records=[]
for q in sorted(O.rglob('progress.json')):
    state=json.loads(q.read_text())
    if state.get('schema')!='hexagon-v27-cpu-controls-v1':continue
    for step in state['steps']:
        if step['status']!='passed':continue
        log=q.parent/(step['label']+'.log');txt=log.read_text();m=step['measurement']
        def lines(prefix):
            return [{k:(float(v) if re.fullmatch(r'[-+0-9.eE]+',v) else v) for k,v in re.findall(r'(\w+)=([^ ]+)',line)} for line in txt.splitlines() if line.startswith(prefix+' ')]
        records.append(dict(suite=str(q.parent.relative_to(O)),loop=state['loop'],label=step['label'],version=step['label'].split('-')[0],threads=state['threads'],binary_sha256=state['binaries'][step['label'].split('-')[0]],flags=state['flags'],prime=step['prime'],cpu_hours=step['cpu_hours'],wall_seconds=m['wall_seconds'],
            peak_rss_bytes=max(m['peak_summed_rss_bytes'],max(x.get('process_peak_bytes',0) for x in step['phase_profiles'])),
            sampled_peak_bytes=m['peak_summed_rss_bytes'],cgroup_events=step['cgroup_events'],
            phases=step['phase_profiles'],structural_profile=lines('STRUCTURAL_PROFILE'),dependency_profile=lines('DEPENDENCY_PROFILE'),kernel_work=lines('KERNEL_WORK'),wide_kernel=lines('WIDE_KERNEL'),compression=lines('OPERATOR_COMPRESSION'),
            work_schedule=lines('WORK_PHYSICS_SELECT'),fallbacks=[x for x in txt.splitlines() if 'STRUCTURE_FALLBACK' in x],
            native_coefficients_equal=step['all_native_coefficients_equal'],all_physical_checks=all('JOINT_PHYSICAL_CHECK operator='+name+' passed=1' in txt for name in ['collinear','mrk','ope_small_s']),
            complete_seam_check='V18_SEAM_VERIFIED full_operator=1' in txt,reference=step['reference'],reference_sha256=step['reference_sha256'],cache_sha256=step['cache_sha256'],log=str(log),measurement=str(q.parent/(step['label']+'.json'))))
groups=[]
for key in sorted({(r['suite'],r['version']) for r in records}):
    xs=[r for r in records if (r['suite'],r['version'])==key]
    groups.append(dict(suite=key[0],version=key[1],loop=xs[0]['loop'],runs=len(xs),threads=xs[0]['threads'],flags=xs[0]['flags'],cpu_hours_mean=statistics.mean(r['cpu_hours'] for r in xs),cpu_hours_total=sum(r['cpu_hours'] for r in xs),peak_rss_bytes=max(r['peak_rss_bytes'] for r in xs),correctness_passed=all(r['native_coefficients_equal'] and r['all_physical_checks'] and r['complete_seam_check'] for r in xs)))
result=dict(status='interim',baseline='V26',acceptance=dict(cpu_saving_required_strictly_greater_than=.15,rss_ratio_max=1.05),scope='Fresh native finite-field solves on warm precomputed inputs, including input loading, learning/replay plan I/O, checkpoint verification, final full seam and all existing physical checks. CPU=user+system. Excludes compilation, cold inputs and new QQ reconstruction.',groups=groups,runs=records)
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
for g in groups:print(g['suite'],g['runs'],f"{g['cpu_hours_mean']:.6f} CPUh",f"{g['peak_rss_bytes']/1e9:.3f} GB",g['correctness_passed'])
