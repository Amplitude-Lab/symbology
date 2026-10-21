#!/usr/bin/env python3
"""Verify and package completed additional checks. Does not claim cold timing."""
import json,re,shutil,sys,time,hashlib
from pathlib import Path
from hexagon_ten_loop_benchmark import sha
from hexagon_v17_cache_migration import verify_receipt
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009';D=R/'benchmarks/hexagon-mhv-ten-loop';B=R/'output_hexagon_17p3_20261007'
def main():
 x=json.loads((V/'additional-verification.json').read_text());assert x['status']=='passed'
 comparison=json.loads((V/'replay-comparison.json').read_text());assert comparison['passed']
 plan=json.loads((V/'physical-validation-plan.json').read_text());assert sha(V/'physical-validation-plan.json')==x['plan_sha256']
 bindings={};physical=[]
 for a in x['attempts']:
  for key in ['report','log']:
   assert sha(a[key])==a[key+'_sha256'];bindings[a[key]]=a[key+'_sha256']
  m=a['measurement'];assert m==json.loads(Path(a['report']).read_text());assert m['exit_code']==a['expected_exit_code'];assert not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
  lines=Path(a['log']).read_text().splitlines();rows=[json.loads(l.split(' ',1)[1]) for l in lines if l.startswith('PHYSICAL_SUMMARY ')];parts=[json.loads(l.split(' ',1)[1]) for l in lines if l.startswith('PHYSICAL_HOLDOUT ')]
  if a['name'].startswith('native-fresh'):
   assert len(rows)==1 and rows[0]['passed'] and rows[0]['mismatches']==0
   physical.append({'name':a['name'],'rows':rows[0]['rows'],'mismatches':0,'new_ope_coefficients':sum(r['S_power_above_10'] for r in parts),'measurement':m})
  if a['name']=='negative-control-scaled-symbol':assert rows==[{'rows':64,'mismatches':64,'multiplier':2,'passed':False}]
 assert len(physical)==6 and sorted(r['rows'] for r in physical)==[64,64,1280,1280,3192,3192]
 assert sum(r['new_ope_coefficients'] for r in physical)==1824
 for p,h in plan['inputs'].items():assert sha(p)==h;bindings[p]=h
 code=json.loads((V/'validation-build-binding.json').read_text())
 for p,h in code['artifacts'].items():assert sha(p)==h['sha256'];bindings[p]=h['sha256']
 old=json.loads((B/'ten-asymmetric-runner/final-verification.json').read_text())
 for root in [V/'native-qq-17p3',B/'ten-asymmetric-qq-16p4']:
  seal=dict(line.split('\t',1) for line in (root/'qq_complete.tsv').read_text().splitlines())
  assert seal['schema']=='irr-physics-qq-v1' and sha(root/'input_identity.txt')==seal['identity']
  for name,h in seal.items():
   if name not in ['schema','identity']:
    assert Path(name).name==name and sha(root/name)==h;bindings[str(root/name)]=h
  manifest=json.loads((root/'manifest.json').read_text())
  for p,h in manifest['input_hashes'].items():assert sha(p)==h;bindings[p]=h
 for p,h in old['artifact_hashes'].items():assert sha(p)==h;bindings[p]=h
 verify_receipt(B/'ten-v17-runner/v16-prime-import.json')
 native=json.loads((V/'native-export-install.json').read_text());assert native['passed']
 for p,h in native['inputs'].items():assert sha(p)==h;bindings[p]=h
 assert sha(V/'native-qq-17p3/qq_complete.tsv')==native['qq_complete_sha256']
 log=(V/'native-exact-seam.log').read_text();seam=json.loads(next(l.split(' ',1)[1] for l in log.splitlines() if l.startswith('JOINT_CERTIFICATE ')))
 assert seam['exact_QQ_original_seam'] and seam['original_rows']==2053259 and seam['bound_bits']<seam['modulus_bits'] and seam['directions']==0
 assert 'ORIGIN_CHECK PASS words=44' in (V/'native-exact-origin.log').read_text()
 assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=518048' in (V/'exact-split-equality.log').read_text()
 assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=275246' in (V/'native-independent-wxf.log').read_text()
 replay={};live_resources={}
 for mode in ['mapped','heap']:
  d=comparison['runs'][mode];assert d['status']=='passed' and d['phase_result'][0]['all_coefficients_equal']
  for p,h in {**d['input_hashes'],**d['artifacts']}.items():assert sha(p)==h;bindings[p]=h
  replay[mode]={'measurement':d['measurement'],'phase_result':d['phase_result'][0]}
  snapshots=[json.loads(l) for l in (V/('replay-'+mode)/'process-snapshots.jsonl').read_text().splitlines()]
  resource_rows=[r for r in snapshots if 'resources' in r];assert resource_rows
  for r in resource_rows:
   q=r['resources'];assert 'hexagon_seam_replay 17 3 ' in r['cmdline']
   assert int(q['memory.high'])==48000000000 and 50000000000-4096<int(q['memory.max'])<=50000000000 and int(q['memory.swap.max'])==0
   events=dict(l.split() for l in q['memory.events'].splitlines());assert all(int(v)==0 for v in events.values())
  live_resources[mode]={'verified_live_samples':len(resource_rows),'hard_memory_bytes':int(resource_rows[-1]['resources']['memory.max']),'high_bytes':48000000000,'swap_max_bytes':0,'sampled_events_all_zero':True,'post_exit_status_note':'The transient unit was collected before the post-exit systemctl query, which returned default infinity limits. Live cgroup snapshots and the launch command verify the actual bounded job.'}
 assert comparison['runs']['mapped']['input_hashes']==comparison['runs']['heap']['input_hashes']
 (D/'evidence').mkdir(exist_ok=True)
 for f in [V/'additional-verification.json',V/'replay-comparison.json',V/'physical-validation-plan.json',V/'validation-build-binding.json',V/'minheight-experiment.json',V/'minheight.resources.json']:
  shutil.copy2(f,D/'evidence'/f.name);bindings[str(f)]=sha(f)
 for a in x['attempts']:
  for k in ['report','log']:
   f=Path(a[k]);shutil.copy2(f,D/'evidence'/f.name)
 for mode in ['mapped','heap']:
  for name in ['manifest.json','resources.json','run.log','systemd-status.txt']:
   f=V/('replay-'+mode)/name;shutil.copy2(f,D/'evidence'/(mode+'-'+name));bindings[str(f)]=sha(f)
  for name in ['process-snapshots.jsonl','resource-snapshots.jsonl']:
   f=V/('replay-'+mode)/name
   if f.exists():shutil.copy2(f,D/'evidence'/(mode+'-'+name));bindings[str(f)]=sha(f)
 for split,root in [('17p3',V/'native-qq-17p3'),('16p4',B/'ten-asymmetric-qq-16p4')]:
  for f in (D/'reference'/split).iterdir():assert sha(f)==sha(root/f.name)
 result={'schema':'hexagon-native-additional-verification-and-benchmark-v1','passed':True,'created':time.time(),'primary_native_root':str(V/'native-qq-17p3'),'equivalent_shifted_root':str(B/'ten-asymmetric-qq-16p4'),'native_original_seam':seam,'origin_predictions_per_representation':44,'exact_all_coordinate_split_equality':True,'physical_checks':physical,'physical_rows_total':sum(r['rows'] for r in physical),'new_OPE_distinct_coefficients':912,'physical_scope':plan['holdout_scope'],'ope_scope':plan['ope_scope'],'replays':replay,'phase_speed_ratio_mapped_over_heap':replay['mapped']['phase_result']['phase_seconds']/replay['heap']['phase_result']['phase_seconds'],'conditional_bonus_assumption':True,'unrestricted_uniqueness_claimed':False,'cold_end_to_end_measured':False,'input_hashes':bindings}
 result['replay_live_resource_verification']=live_resources
 (V/'final-additional-review.json').write_text(json.dumps(result,indent=2)+'\n');(D/'additional-measurements.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['input_hashes','physical_checks','replays']},indent=2))
if __name__=='__main__':main()
