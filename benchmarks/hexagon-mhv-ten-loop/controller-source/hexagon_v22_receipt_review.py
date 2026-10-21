#!/usr/bin/env python3
"""Review saved V19--V22 measurements; no numerical jobs or producer changes."""
from pathlib import Path
import json,hashlib,statistics,math,time,shutil
R=Path(__file__).resolve().parents[1];r=R/'output_hexagon_v19_20261009';v=R/'output_hexagon_validation_20261009';d=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v22-nine-loop';d.mkdir(exist_ok=True)
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
x=json.loads((r/'RESULTS.json').read_text());c=json.loads((r/'completion-receipt.json').read_text());profile=json.loads((r/'accepted-nine-loop-profile.json').read_text());assert c['status']=='passed' and x['selected']=='V22/16'
for group in ['artifacts','mathematical_receipts','dynamic_libraries_sha256']:
 for p,h in c[group].items():bind(p,h)
source_files=0
for z in x['source_audits']:
 p=Path(z['receipt']);bind(p,z['sha256']);q=json.loads(p.read_text());assert len(q['hashes'])==z['files_verified']
 for name,h in q['hashes'].items():bind(p.parent/name,h);source_files+=1
 keep(p,p.parent.name+'-'+p.name)
for z in x['runs']:
 receipt=Path(z['receipt']);s=json.loads(receipt.read_text());a=next(a for a in s['steps'] if a['label']==z['label']);assert a['status']=='passed';m=a['measurement'];assert json.loads((receipt.parent/(z['label']+'.json')).read_text())==m
 assert m['exit_code']==0 and not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
 assert math.isclose(z['cpu_hours'],(m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,rel_tol=1e-12)
 assert z['wall_seconds']==m['wall_seconds'] and z['peak_rss_bytes']==m['peak_summed_rss_bytes'];assert all(n==0 for n in a['cgroup_events'].values())
 # The oldest first-prime receipt records only the reference hash. Its
 # historical path is known and must match that hash before use.
 reference_path=Path(a.get('reference',str(R/'output_hexagon_17p3_20261007/v17-nine-control/nine/joint/cache'/f"p{z['prime']}.bin")))
 bind(z['cache'],a['cache_sha256']);bind(reference_path,a['reference_sha256']);bind(z['binary'],z['binary_sha256']);actual=Path(z['cache']).read_bytes();reference=reference_path.read_bytes();assert actual[:64]==reference[:64] and actual[72:]==reference[72:]
 log=receipt.parent/(z['label']+'.log');bind(log,a['log_sha256']);t=log.read_text()
 for operator in ['collinear','mrk','ope_small_s']:assert 'JOINT_PHYSICAL_CHECK operator='+operator+' passed=1' in t
 assert z['phases']==a['phase_profiles']
 for p in [receipt,log,log.with_suffix('.json')]:keep(p,receipt.parent.name+'-'+p.name)
for group,reported in x['groups'].items():
 selected=[z for z in x['runs'] if z['version']+'/'+str(z['threads'])==group];first=[z for z in selected if z['prime']==1152921504606847009];assert len(first)==reported['first_prime_repeats'];assert math.isclose(statistics.mean(z['cpu_hours'] for z in first),reported['cpu_hours_mean'],rel_tol=1e-12);assert math.isclose(statistics.mean(z['wall_seconds'] for z in first),reported['wall_seconds_mean'],rel_tol=1e-12);assert max(z['peak_rss_bytes'] for z in selected)==reported['peak_rss_bytes_max'] and min(z['peak_rss_bytes'] for z in selected)==reported['peak_rss_bytes_min']
g=x['groups'];assert 1-g['V22/16']['cpu_hours_mean']/g['V17/16']['cpu_hours_mean']>0.15;assert g['V22/16']['peak_rss_bytes_max']<=1.05*g['V17/16']['peak_rss_bytes_min'];assert profile['acceptance']==g['V22/16'];bind(profile['binary'],profile['binary_sha256'])
qroot=r/'qq-only-references/L9';q=json.loads((qroot/'snapshot-provenance.json').read_text());assert not (qroot/'cache').exists()
for name,h in q['files'].items():bind(qroot/name,h);bind(Path(q['source'])/name,h)
keep(qroot/'snapshot-provenance.json','L9-QQ-only-snapshot-provenance.json')
shift=r/'V20-shifted-forced-qq';link=json.loads((r/'reference-binding-audit.json').read_text());assert link['passed']
for z in link['shift_bindings']:
 prime=z['prime'];native=r/'controls-v20/V20/joint/cache'/f'p{prime}.bin';shifted=shift/'cache'/native.name;bind(native,z['native_sha256']);bind(shifted,z['shifted_sha256']);assert (shift/f'source_p{prime}.sha256').read_text().strip()==sha(native)
 target=native.read_bytes()
 for a in x['runs']:
  if a['version']=='V22' and a['prime']==prime:
   actual=Path(a['cache']).read_bytes();assert actual[:64]==target[:64] and actual[72:]==target[72:]
for tag in ['forced-qq-shift','forced-qq-compare-p1','forced-qq-compare-p2']:
 m=json.loads((r/(tag+'.json')).read_text());assert m['exit_code']==0 and not m['timeout'];log=(r/(tag+'.log')).read_text()
 if 'compare' in tag:assert 'PASS affine spaces directions=0 columns=195155' in log and any(str(qroot.relative_to(R))==a for a in m['command'])
 for ext in ['json','log']:keep(r/(tag+'.'+ext))
mathroot=R/'bonusrelation/solver_audit_20261009';triple=json.loads((mathroot/'completed-ten-triple-receipt.json').read_text());assert triple['passed']
for z in triple['checks']:
 assert z['Eq92_nonzeros']==0;bind(R/'output_hexagon_17p3_20261007/ten-v17-main/joint/cache'/f"p{z['prime']}.bin",z['cache_sha256'])
for name in ['REPORT.txt','independent-verification.json','native-triple-independent-verification.json','completed-ten-triple-receipt.json','triple-control-v3.log']:
 keep(mathroot/name,'math-'+name)
for name in ['RESULTS.json','REPORT.txt','completion-receipt.json','accepted-nine-loop-profile.json','execution-gate.json','reference-binding-audit.json','ope-weight-control.json','ope-weight-control.log']:keep(r/name)
t=(r/'ope-weight-control.log').read_text();assert 'cases=2160' in t and t.count('coefficients=1178100 equal=1')==2
bind(__file__)
out={'schema':'v22-nine-loop-receipt-review-v1','passed':True,'created':time.time(),'numeric_runs_reviewed':len(x['runs']),'source_files_reviewed':source_files,'all_native_payload_comparisons_repeated':True,'CPU_means_and_RSS_gate_recomputed':True,'QQ_only_snapshot_verified_against_original':True,'V22_to_V20_payload_binding_at_both_primes':True,'ten_loop_triple_receipt_and_cache_bindings_reviewed':True,'new_calculations_launched':False,'prior_validation_correction':'The V18 pipeline_compare command could use saved Fp reference caches. The new nine-loop comparison forces actual QQ reduction with a cache-free hash-matched snapshot. This correction changes the stated validation method, not the accepted ten-loop result.','mathematical_scope':'Triple checks are a posteriori at two primes; no QQ/all-loop triple proof or new constraints in V22.','review_repair':'Initial review stopped on a missing reference-path field in the earliest first-prime receipts. Explicit historical reference paths are now verified against the recorded hashes; no assertion was bypassed.','input_hashes':bindings}
f=v/'v22-nine-loop-review.json';f.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(f,d/f.name);print(json.dumps({k:v for k,v in out.items() if k!='input_hashes'},indent=2));print('hash bindings',len(bindings))
