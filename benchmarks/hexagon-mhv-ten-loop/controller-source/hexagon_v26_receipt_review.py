#!/usr/bin/env python3
"""Read-only numerical evidence review; copies receipts into the benchmark."""
from pathlib import Path
import json,hashlib,statistics,math,shutil,time
R=Path(__file__).resolve().parents[1];r=R/'output_hexagon_v23_20261009';v=R/'output_hexagon_validation_20261009';d=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v26-generic';d.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
bindings={}
def bind(p,expected=None):
 p=Path(p).resolve();h=sha(p);assert expected is None or h==expected,('hash',str(p));bindings[str(p)]=h;return h
def keep(p,name=None):
 p=Path(p);bind(p);shutil.copy2(p,d/(name or p.name))
x=json.loads((r/'RESULTS.json').read_text());c=json.loads((r/'completion-receipt.json').read_text());profile=json.loads((r/'accepted-nine-loop-profile.json').read_text());assert x['status']=='complete' and x['selected']=='V26' and c['state']=='complete'
build=Path(profile['build_receipt']);b=json.loads(build.read_text());assert len(b['hashes'])==c['source_binary_files_verified']
for name,h in b['hashes'].items():bind(build.parent/name,h)
bind(profile['binary'],profile['binary_sha256']);keep(build)
for plan in c['plans']:
 p=Path(plan['path']);bind(p,plan['sha256']);assert p.stat().st_size==plan['bytes']
 assert plan['learning_prime']==1152921504606847009
for z in x['runs']:
 log=Path(z['log']);receipt=log.parent/'progress.json';state=json.loads(receipt.read_text());a=next(a for a in state['steps'] if a['label']==z['label']);m=json.loads(Path(z['measurement']).read_text());assert a['status']=='passed' and m==a['measurement']
 assert m['exit_code']==0 and not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
 assert math.isclose((m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,z['cpu_hours'],rel_tol=1e-12)
 assert m['wall_seconds']==z['wall_seconds'] and m['peak_summed_rss_bytes']==z['sampled_peak_rss_bytes']
 assert a['phase_profiles']==z['phases'];kernel=max(q.get('process_peak_bytes',0) for q in z['phases']);assert kernel==z['kernel_peak_rss_bytes'];assert z['peak_rss_bytes']==max(kernel,z['sampled_peak_rss_bytes'])
 assert all(n==0 for n in z['cgroup_events'].values());assert z['cgroup_events']==a['cgroup_events'];bind(log,a['log_sha256'])
 binary=Path(a['command'][-12]);assert binary.name=='joint_physics_'+z['version'].lower();bind(binary,state['binaries'][z['version']])
 root=Path(a['command'][-6]);cache=root/'joint/cache'/f"p{z['prime']}.bin";bind(cache,z['cache_sha256']);bind(z['reference'],z['reference_sha256']);actual=cache.read_bytes();ref=Path(z['reference']).read_bytes();assert actual[:64]==ref[:64] and actual[72:]==ref[72:]
 t=log.read_text()
 for op in ['collinear','mrk','ope_small_s']:assert 'JOINT_PHYSICAL_CHECK operator='+op+' passed=1' in t
 assert z['reference_coefficients_equal'] and z['all_physical_checks'] and z['full_seam_check']
 for p in [receipt,log,Path(z['measurement'])]:keep(p,log.parent.name+'-'+p.name)
for g in x['groups']:
 runs=[z for z in x['runs'] if z['suite']==g['suite'] and z['version']==g['version']];assert len(runs)==g['runs'];assert math.isclose(sum(z['cpu_hours'] for z in runs),g['cpu_hours_total'],rel_tol=1e-12);assert math.isclose(statistics.mean(z['cpu_hours'] for z in runs),g['cpu_hours_mean'],rel_tol=1e-12);assert math.isclose(statistics.mean(z['wall_seconds'] for z in runs),g['wall_seconds_mean'],rel_tol=1e-12);assert max(z['peak_rss_bytes'] for z in runs)==g['peak_rss_bytes']
assert len(x['runs'])==c['total_fresh_native_runs']==21
selected9=next(g for g in x['groups'] if g['suite']=='L9-v26');selected8=next(g for g in x['groups'] if g['suite']=='L8-v26');fresh22=next(g for g in x['groups'] if g['suite']=='L9-combined' and g['version']=='V22')
for key,g in [('nine_loop',selected9),('eight_loop',selected8)]:
 assert g['cpu_hours_mean']==profile[key]['two_prime_mean_cpu_hours'];assert g['peak_rss_bytes']==profile[key]['peak_rss_bytes']
assert math.isclose(1-selected9['cpu_hours_mean']/fresh22['cpu_hours_mean'],profile['nine_loop']['saving_vs_matched_v22'],rel_tol=1e-12)
assert 1-selected9['cpu_hours_mean']/x['baseline']['cpu_hours_mean']>0.15
assert selected9['peak_rss_bytes']<=1.05*x['baseline']['conservative_peak_rss_bytes']
old8=json.loads(Path(profile['eight_loop']['original_v17_measurement']).read_text());assert math.isclose((old8['cpu_user_seconds']+old8['cpu_system_seconds'])/3600,profile['eight_loop']['original_v17_cpu_hours'],rel_tol=1e-12);bind(profile['eight_loop']['original_v17_measurement'])
assert selected8['peak_rss_bytes']<=1.05*old8['peak_summed_rss_bytes'] and 1-selected8['cpu_hours_mean']/profile['eight_loop']['original_v17_cpu_hours']>0.15
link=json.loads((r/'reference-binding-audit.json').read_text());bind(link['prior_forced_QQ_binding_audit'],link['prior_forced_QQ_binding_audit_sha256'])
for z in link['native_comparisons']:
 bind(z['candidate'],z['candidate_sha256']);bind(z['reference'],z['reference_sha256']);a=Path(z['candidate']).read_bytes();b=Path(z['reference']).read_bytes();assert a[:64]==b[:64] and a[72:]==b[72:]
# Preserve the stronger nine-loop QQ-only chain established by the prior review.
old=v/'v22-nine-loop-review.json';q=json.loads(old.read_text());assert q['passed']
for p,h in q['input_hashes'].items():bind(p,h)
bind(old)
for field in ['generic_controls','affine_family_controls']:bind(c[field],c[field+'_sha256']);keep(c[field])
assert 'GENERIC_CONTROL passed=1 field_cases=96 successful_replays=48 guarded_fallbacks=2 factor_cases=72' in Path(c['generic_controls']).read_text()
lines=Path(c['affine_family_controls']).read_text().splitlines();assert len(lines)==5 and all('coefficient_equality=1 original_equations=1' in line for line in lines)
for name in ['RESULTS.json','REPORT.txt','accepted-nine-loop-profile.json','completion-receipt.json','execution-gate.json','reference-binding-audit.json']:keep(r/name)
bind(__file__)
out={'schema':'hexagon-v26-generic-receipt-review-v1','passed':True,'created':time.time(),'native_runs_reviewed':21,'selected_build_files_reviewed':c['source_binary_files_verified'],'native_payload_equalities_rechecked':True,'CPU_means_and_conservative_RSS_recomputed':True,'nine_loop_QQ_reference_chain_verified':True,'selected_profile':profile,'scope':'Review of archived controls. Two-prime CPU means include learning and replay, not statistical repetitions. RSS=max(100ms sample,kernel process high-water mark). No new numerical solve or ten-loop performance test.','input_hashes':bindings}
p=v/'v26-generic-review.json';p.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(p,d/p.name)
print(json.dumps({k:v for k,v in out.items() if k not in ['input_hashes','selected_profile']},indent=2));print('bindings',len(bindings))
