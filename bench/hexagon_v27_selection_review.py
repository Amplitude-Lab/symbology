#!/usr/bin/env python3
"""Bind the selection amendment without rewriting historical measurement receipts."""
from pathlib import Path
import hashlib,json,math,shutil,time
R=Path(__file__).resolve().parents[1];N=R/'output_hexagon_v27_20261009';V=R/'output_hexagon_validation_20261009'
D=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v27-selection-amendment';D.mkdir(exist_ok=True)
bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def bind(p,expected=None):
 p=Path(p).resolve();h=sha(p);assert expected is None or h==expected,('hash',str(p));bindings[str(p)]=h;return h
def same(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
a=json.loads((N/'selection-amendment.json').read_text());p=json.loads((N/'current-nine-loop-profile.json').read_text())
prior=V/'v27-experimental-review.json';old=json.loads(prior.read_text());assert old['passed']
relocations={}
for path,h in old['input_hashes'].items():
 if sha(path)==h:bind(path,h)
 else:
  # The active gate was explicitly superseded. Require its exact old bytes
  # in both the archived benchmark and the amendment's historical snapshot.
  assert Path(path)==N/'execution-gate.json',('unexpected evidence change',path)
  archived=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v27-experimental/execution-gate.json'
  superseded=N/'execution-gate-superseded-v26-threshold.json'
  bind(archived,h);bind(superseded,h);current_hash=bind(path)
  g=json.loads(Path(path).read_text())
  assert g['selected_working_benchmark']=='V27D' and g['qualification_baseline']=='V16'
  assert g['qualification_status']=='not_verified_against_V16_at_nine_loops'
  assert Path(g['amendment'])==N/'selection-amendment.json' and Path(g['superseded_gate'])==superseded
  relocations[path]={'archived_path':str(archived),'historical_sha256':h,'amended_sha256':current_hash,'reason':'Selection gate superseded by explicit V16-baseline amendment; original bytes preserved.'}
bind(prior)
for path,h in a['evidence'].items():bind(path,h)
assert a['selected_working_benchmark']==p['selected']==p['current_nine_loop_default']=='V27D'
assert a['threads']==p['threads']==16 and p['state']=='selected_working_benchmark'
assert a['qualification']==p['qualification'] and a['qualification']['nine_loop_status']=='not_verified_against_V16'
assert a['qualification']['cpu_saving_required_strictly_greater_than']==0.15 and a['qualification']['rss_ratio_max']==1.05
assert a['qualification']['incremental_comparator']==p['comparison_reference']=='V26'
assert a['qualification']['baseline'].startswith('V16')
assert not a['production_modified'] and a['production_memory_limit_bytes']==50000000000
q=old['experimental_profile']
for key in ['package','binary','binary_sha256','build_receipt','baseline','nine_loop','eight_loop','nine_loop_eight_workers','matched_eight_worker_first_prime','comparison_nine','comparison_nine_eight_workers','comparison_eight','environment','hardware','unset_environment','ablations','route_pair']:
 assert p[key]==q[key],key
assert p['superseded_acceptance_decision']==q['acceptance']
bind(p['binary'],p['binary_sha256']);bind(p['build_receipt'])
c=a['comparison_with_v26'];b=p['baseline'];n=p['nine_loop']
for kind,field in [('cpu_hours','cpu_hours_mean'),('wall_seconds','wall_seconds_mean'),('rss_bytes','peak_rss_bytes')]:
 same(c[kind+'_v26'],b[field]);same(c[kind+'_v27d'],n[field])
same(c['cpu_saving'],1-n['cpu_hours_mean']/b['cpu_hours_mean'])
same(c['wall_saving'],1-n['wall_seconds_mean']/b['wall_seconds_mean'])
same(c['rss_ratio'],n['peak_rss_bytes']/b['peak_rss_bytes'])
e=a['qualification']['eight_loop_historical_reference'];bind(e['path'],e['sha256']);m=json.loads(Path(e['path']).read_text())
same(e['cpu_hours'],(m['cpu_user_seconds']+m['cpu_system_seconds'])/3600)
assert e['peak_sampled_rss_bytes']==m['peak_summed_rss_bytes']
same(e['v27d_cpu_saving'],1-p['eight_loop']['cpu_hours_mean']/e['cpu_hours'])
same(e['v27d_rss_ratio_to_sampled_reference'],p['eight_loop']['peak_rss_bytes']/e['peak_sampled_rss_bytes'])
for name in ['selection-amendment.json','current-nine-loop-profile.json','SELECTION_CORRECTION.txt','execution-gate.json','execution-gate-superseded-v26-threshold.json']:
 src=N/name;bind(src);dst=D/name
 if dst.exists():assert sha(src)==sha(dst)
 else:shutil.copy2(src,dst)
bind(__file__)
out={'schema':'hexagon-v27-selection-amendment-review-v1','created':time.time(),'passed':True,'current_working_benchmark':'V27D/16','qualification_baseline':'V16','qualification_status':'not_verified_against_V16','incremental_comparator':'V26','historical_measurements_and_review_unchanged':True,'historical_input_relocations':relocations,'incremental_comparison_recomputed':True,'eight_loop_historical_V16_comparison_recomputed':True,'early_physics_policy':'All required physical equations or lazy generators available from the start; assembly and elimination deferred by cost, rank gain, fill and memory while preserving affine RHS/origin and completion checks. A generic pending-block scheduler is not established by these measurements.','scope':'Selection and documentation amendment only. No new solve; no nine-loop V16 qualification or ten-loop V27 timing claimed.','amendment':a,'input_hashes':bindings}
dest=V/'v27-selection-amendment-review.json';dest.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(dest,D/dest.name)
print(json.dumps({k:v for k,v in out.items() if k not in ['amendment','input_hashes']},indent=2));print('bindings',len(bindings))
