#!/usr/bin/env python3
"""Read-only review of two completed primes; never touches active prime 3."""
from pathlib import Path
import hashlib,json,math,re,shutil,time
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009';N=R/'output_hexagon_ten_v27d_20261009/runner'
D=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v27-ten-prime2';D.mkdir(exist_ok=True)
bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def bind(p,expected=None):
 p=Path(p).resolve();h=sha(p);assert expected is None or h==expected,('hash',str(p));bindings[str(p)]=h;return h
def keep(p):
 p=Path(p);bind(p);dst=D/p.name
 if dst.exists():assert sha(dst)==sha(p),('archive collision',str(dst))
 else:shutil.copy2(p,dst)
def same(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
def records(text,prefix):
 return [dict(re.findall(r'(\w+)=([^ ]+)',line)) for line in text.splitlines() if line.startswith(prefix+' ')]
p=N/'audit/prime2-complete-comparison.json';c=json.loads(p.read_text());keep(p)
assert c['prime']==1152921504606847067 and c['correctness']['directions']==0
for path,h in c['sources'].items():bind(path,h)
logs={};derived={}
for index in [1,2]:
 key=str(index);metric=c['metrics'][key];assert metric['version']=='V27D/16'
 for phase in ['native','convert']:
  report=next(Path(s) for s in c['sources'] if f'/audit/p{index:02}-{phase}-' in s)
  m=json.loads(report.read_text());assert m['exit_code']==0
  assert not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
  label='native' if phase=='native' else 'conversion'
  same((m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,metric[label+'_cpu_hours'])
  same(m['wall_seconds'],metric[label+'_wall_seconds'])
  if phase=='native':assert m['peak_summed_rss_bytes']==metric['native_peak_sampled_rss_bytes']
  keep(report);log=next(Path(s) for s in c['sources'] if f'/log/p{index:02}-{phase}-' in s);keep(log)
  if phase=='native':logs[key]=log.read_text()
 t=logs[key];assert 'V18_SEAM_VERIFIED full_operator=1 rows=989800 directions=0' in t
 for op in ['collinear','mrk','ope_small_s']:assert 'JOINT_PHYSICAL_CHECK operator='+op+' passed=1' in t
 phase_rows=records(t,'PROFILE');parsed={q['phase']:{k:float(v) for k,v in q.items() if k!='phase'} for q in phase_rows}
 for name,q in c['phases'][key].items():assert parsed[name]==q
 assert max(float(q.get('process_peak_bytes',0)) for q in phase_rows)==metric['native_kernel_peak_rss_bytes']
 assert records(t,'V18_SEAM_COMPLETE')[-1]==c['seam_profiles'][key]
 structural=records(t,'STRUCTURAL_PROFILE')[0]
 assert int(structural['replay'])==int(metric['discovery_replay'])
 derived[key]={'dependency':{k:int(v) for k,v in records(t,'DEPENDENCY_PROFILE')[0].items()},'kernel':{k:int(v) for k,v in records(t,'KERNEL_WORK')[0].items()},'structural':{k:float(v) for k,v in structural.items()}}
for kind,q in c['correctness']['coefficient_comparisons'].items():
 bind(q['path'],q['sha256']);bind(q['reference'],q['reference_sha256'])
 a=Path(q['path']).read_bytes();b=Path(q['reference']).read_bytes();assert len(a)==72+8*q['columns'] and a[:64]==b[:64] and a[72:]==b[72:]
 assert q['every_coefficient_equal']
a,b=c['metrics']['1'],c['metrics']['2'];o=c['observed_comparison']
same(o['native_cpu_saving_replay_vs_learning'],1-b['native_cpu_hours']/a['native_cpu_hours'])
same(o['native_wall_saving_replay_vs_learning'],1-b['native_wall_seconds']/a['native_wall_seconds'])
same(o['sampled_peak_rss_ratio_replay_vs_learning'],b['native_peak_sampled_rss_bytes']/a['native_peak_sampled_rss_bytes'])
s1=c['phases']['1']['compute_seam_complete']['cpu_seconds'];s2=c['phases']['2']['compute_seam_complete']['cpu_seconds']
same(o['seam_phase_cpu_saving'],1-s2/s1)
same(o['seam_phase_cpu_hours_learning'],s1/3600);same(o['seam_phase_cpu_hours_replay'],s2/3600)
same(o['seam_cpu_saving_fraction_of_whole_learning_native'],(s1-s2)/3600/a['native_cpu_hours'])
same(c['completed_native_and_conversion_cpu_hours'],sum(q['native_cpu_hours']+q['conversion_cpu_hours'] for q in c['metrics'].values()))
prod=R/'benchmarks/hexagon-mhv-ten-loop/production-measurements.json';history=json.loads(prod.read_text());old=next(q for q in history['rows'] if q['prime_index']==2)
assert old['producer']=='V17';h=c['historical_same_prime_reference'];same(h['native_cpu_hours'],(old['native']['cpu_user_seconds']+old['native']['cpu_system_seconds'])/3600)
same(h['native_wall_seconds'],old['native']['wall_seconds']);assert h['native_peak_sampled_rss_bytes']==old['native']['peak_summed_rss_bytes']
same(h['CPU_saving_V27D_vs_V17_same_prime'],1-b['native_cpu_hours']/h['native_cpu_hours'])
same(h['RSS_ratio_V27D_vs_V17_same_prime'],b['native_peak_sampled_rss_bytes']/h['native_peak_sampled_rss_bytes'])
assert c['memory_events']['high']>0 and all(c['memory_events'][k]==0 for k in ['max','oom','oom_kill','oom_group_kill'])
assert c['configured_hard_limit_bytes']==50000000000 and c['swap_limit_bytes']==0 and not c['active_producer_modified']
dependency=derived['2']['dependency'];coverage=dependency['route_steps']/(dependency['route_steps']+dependency['heap_steps'])
assert dependency['route_hits']==4402 and dependency['local_fallbacks']==3739 and 0.019<coverage<0.021
discovery_saving=1-c['phases']['2']['compute_discovery']['cpu_seconds']/c['phases']['1']['compute_discovery']['cpu_seconds']
source=R/'output_hexagon_v27_20261009/package-v27d/src/shared/streaming_echelon.hpp'
prior=V/'v27-ten-prime1-review.json';q=json.loads(prior.read_text());assert q['passed']
for path,h in q['input_hashes'].items():bind(path,h)
bind(prior);bind(source);bind(D/'active-run-observation.json');bind(__file__)
out={'schema':'hexagon-v27d-ten-replay-review-v1','created':time.time(),'passed':True,'native_primes_measured':2,'prime2_native_coefficients_recompared':275246,'prime2_shifted_coefficients_recompared':518048,'complete_factored_seam_rows_per_prime':989800,'CPU_RSS_phase_and_route_evidence_recomputed':True,'discovery_CPU_saving_fraction':discovery_saving,'discovery_route_step_fraction':coverage,'discovery_profiles':derived,'new_QQ_certificates_complete':False,'same_prime_replay_off_control':False,'active_producer_unchanged':True,'comparison':c,'scope':'Two different primes, learning and replay, not statistical repeats or a same-prime replay-off experiment. Exact original QQ certificates for this fresh chain remain pending. Historical prime2 producer is V17, not V16.','input_hashes':bindings}
dest=V/'v27-ten-prime2-review.json';dest.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(dest,D/dest.name)
print(json.dumps({k:v for k,v in out.items() if k not in ['input_hashes','comparison','discovery_profiles']},indent=2));print('bindings',len(bindings))
