#!/usr/bin/env python3
"""Bind final benchmark results, sources, evidence and physical scope."""
import json,time
from pathlib import Path
from hexagon_ten_loop_benchmark import sha
R=Path(__file__).resolve().parents[1];D=R/'benchmarks/hexagon-mhv-ten-loop';V=R/'output_hexagon_validation_20261009'
def main():
 review=json.loads((V/'final-additional-review.json').read_text());assert review['passed']
 for marker in ['PHYSICAL_RESULTS_PENDING','REPLAY_RESULTS_PENDING']:
  assert marker not in (D/'VALIDATION.md').read_text()
 source=R/'audits/hexagon-ten-loop-2026-10-06/TECHNIQUE.tex';assert 'REPLAY_REPORT_PENDING' not in source.read_text()
 scripts=['hexagon_result_reconstruct.py','hexagon_native_export_install.py','hexagon_physical_holdouts.cpp','hexagon_validation_build.py','hexagon_seam_replay.cpp','hexagon_ten_loop_benchmark.py','hexagon_ten_loop_prime.py','hexagon_replay_suite.py','hexagon_minheight_reconstruction.py','hexagon_additional_validation_run.py','hexagon_validation_finalize.py','hexagon_benchmark_manifest.py','hexagon_joint_certificate.cpp','hexagon_candidate_qq_export.cpp','hexagon_candidate_export_controls.wls','measure_process.py']
 scripts+=['hexagon_crt_asymmetric_probe.py','hexagon_crt_progress_audit.py','hexagon_candidate_export_install.py','hexagon_v17_cache_migration.py','hexagon_runtime_recovery_validate.py','hexagon_ten_loop_runner_v17.py']
 if (V/'v22-nine-loop-review.json').exists():scripts.append('hexagon_v22_receipt_review.py')
 if (V/'v26-generic-review.json').exists():scripts.append('hexagon_v26_receipt_review.py')
 if (V/'v27-experimental-review.json').exists():scripts.append('hexagon_v27_receipt_review.py')
 if (V/'v27-ten-prime1-review.json').exists():scripts.append('hexagon_v27_ten_prime_review.py')
 if (V/'v27-ten-prime2-review.json').exists():scripts.append('hexagon_v27_ten_replay_review.py')
 selection_review=V/'v27-selection-amendment-review.json';selection=None
 if selection_review.exists():
  scripts.append('hexagon_v27_selection_review.py')
  selection=json.loads(selection_review.read_text());assert selection['passed']
  assert all(sha(p)==h for p,h in selection['input_hashes'].items())
 source_dir=D/'controller-source';source_dir.mkdir(exist_ok=True)
 import shutil
 for name in scripts:shutil.copy2(R/'bench'/name,source_dir/name)
 files={str(p.relative_to(D)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(D.rglob('*')) if p.is_file() and p!=D/'manifest.json'}
 r={'schema':'hexagon-ten-native-and-shifted-benchmark-v1','created':time.time(),'status':'verified prepared-basis benchmark','primary':'reference/17p3','secondary':'reference/16p4','native_coordinates':275246,'native_nonzero':246941,'shifted_coordinates':518048,'shifted_nonzero':491197,'native_reconstruction':'Seven CRT primes and one retrospective check; supplementary minimum-product experiment also succeeds with six plus two retrospective checks.','shifted_reconstruction':'Six CRT primes with candidate frozen before two fresh validation primes finished.','prime_provenance':'Prime1 V16; validated identity import. Primes2-8 V17. No committed prime recomputed.','physical_checks':9072,'distinct_unused_OPE_coefficients':912,'assumptions':'Extrapolated MHV double-final-entry bonus relations, fixed FEC/LEC/CG normalization.','scope':review['physical_scope'],'ope_scope':review['ope_scope'],'exact_native_seam':True,'exact_shifted_seam':True,'exact_split_equality':True,'exact_origin_predictions_per_split':44,'unrestricted_uniqueness_claimed':False,'fully_optimized_claimed':False,'cold_end_to_end_measured':False,'new_prime_wrapper_full_run_tested':False,'phase_benchmark':'One sequential mapped and heap replay of the identical frozen V17 finishing function, same prime2 post-OPE checkpoint. Warm input/OS caches; shared workstation; fixed order.','memory_units':'GB = 10^9 bytes; sampled process-tree RSS every0.1s; per-job hard cgroup cap50,000,000,000 bytes.','evidence_review':{'path':str(V/'final-additional-review.json'),'sha256':sha(V/'final-additional-review.json')},'report':{'path':str(source),'sha256':sha(source)},'files':files}
 r['report']['compilation']=json.loads((V/'report-progress.json').read_text())['status']
 lower_review=V/'v18-lower-loop-review.json'
 if lower_review.exists():
  lower=json.loads(lower_review.read_text());assert lower['passed']
  assert all(sha(p)==h for p,h in lower['input_hashes'].items())
  r['later_optimization_controls']={'version':'V18','validated_loops':[8,9],'review':str(lower_review),'review_sha256':sha(lower_review),'ten_loop_V18_measured':False,'ten_loop_producer_changed':False}
 v22_review=V/'v22-nine-loop-review.json'
 if v22_review.exists():
  v22=json.loads(v22_review.read_text());assert v22['passed']
  assert all(sha(p)==h for p,h in v22['input_hashes'].items())
  r['nine_loop_CPU_optimization']={'selected_measured_profile':'V22/16','review':str(v22_review),'review_sha256':sha(v22_review),'ten_loop_V22_measured':False,'ten_loop_producer_changed':False}
  r['reference_comparison_correction']=v22['prior_validation_correction']
  r['additional_triple_checks']='Complete published Eq92 and cyclic images pass at two saved ten-loop primes a posteriori; no new rational triple certificate or all-loop proof.'
 v26_review=V/'v26-generic-review.json'
 if v26_review.exists():
  v26=json.loads(v26_review.read_text());assert v26['passed']
  assert all(sha(p)==h for p,h in v26['input_hashes'].items())
  r['generic_discovery_affine_optimization']={'selected_measured_profile':'V26/16','validated_loops':[8,9],'review':str(v26_review),'review_sha256':sha(v26_review),'native_runs_reviewed':v26['native_runs_reviewed'],'CPU_scope':'One learning and one replay prime per selected workload; not same-prime timing repetitions.','RSS_metric':'Maximum of 100ms sampled RSS and kernel process high-water marks; historical V17 gate reference remains the conservative sampled minimum.','ten_loop_V26_measured':False,'ten_loop_producer_changed':False,'fully_implicit_cross_stage_maps_measured':False,'profile':v26['selected_profile']}
 v27_review=V/'v27-experimental-review.json'
 if v27_review.exists():
  v27=json.loads(v27_review.read_text());assert v27['passed'] and not v27['V27_accepted']
  for p,h in v27['input_hashes'].items():
   relocation=(selection or {}).get('historical_input_relocations',{}).get(p)
   if relocation:
    assert relocation['historical_sha256']==h and sha(relocation['archived_path'])==h
    assert sha(p)==relocation['amended_sha256']
   else:assert sha(p)==h
  r['further_generic_experiments']={'candidate':'V27D','candidate_accepted':False,'retained_nine_loop_profile':'V26','validated_loops':[8,9],'review':str(v27_review),'review_sha256':sha(v27_review),'native_runs_reviewed':v27['native_runs_reviewed'],'baseline':'Fresh V26; previous V17 savings not counted again.','ten_loop_V27_measured':False,'ten_loop_producer_changed':False,'profile':v27['experimental_profile']}
 if selection:
  historical=r['further_generic_experiments']
  r['further_generic_experiments']={'historical_campaign_decision':historical,'recommendation_superseded':True,'current_working_benchmark':'V27D/16','amendment_review':str(selection_review),'amendment_review_sha256':sha(selection_review)}
  r['current_nine_loop_benchmark']={'version':'V27D','threads':16,'state':'selected_working_benchmark','qualification_baseline':'V16','qualification_status':'not_verified_against_V16','incremental_comparator':'V26','amendment':selection['amendment'],'historical_input_relocations':selection['historical_input_relocations'],'physical_constraint_policy':selection['early_physics_policy'],'ten_loop_measured':False,'ten_loop_producer_changed':False}
 ten_review=V/'v27-ten-prime1-review.json'
 if ten_review.exists():
  ten=json.loads(ten_review.read_text());assert ten['passed'] and not ten['combined_CPU_RSS_gate_passed']
  assert all(sha(p)==h for p,h in ten['input_hashes'].items())
  r['fresh_ten_loop_V27D_first_prime']={'review':str(ten_review),'review_sha256':sha(ten_review),'measured_native_primes':1,'replay_measured_at_review':False,'new_rational_reconstruction_complete':False,'new_exact_QQ_certificates_complete':False,'combined_CPU_RSS_gate_passed':False,'scope':ten['scope'],'comparison':ten['comparison'],'optimization_hypothesis':ten['optimization_hypothesis']}
  if 'current_nine_loop_benchmark' in r:
   r['current_nine_loop_benchmark']['ten_loop_measured']=True
   r['current_nine_loop_benchmark']['ten_loop_measurement_scope']='One fresh learning prime only; full multi-prime reconstruction and exact certificates pending.'
 replay_review=V/'v27-ten-prime2-review.json'
 if replay_review.exists():
  replay=json.loads(replay_review.read_text());assert replay['passed'] and replay['native_primes_measured']==2
  assert all(sha(p)==h for p,h in replay['input_hashes'].items())
  r['fresh_ten_loop_V27D_replay']={'review':str(replay_review),'review_sha256':sha(replay_review),'measured_native_primes':2,'new_rational_reconstruction_complete':False,'new_exact_QQ_certificates_complete':False,'same_prime_replay_off_control':False,'scope':replay['scope'],'comparison':replay['comparison'],'discovery_profiles':replay['discovery_profiles'],'discovery_CPU_saving_fraction':replay['discovery_CPU_saving_fraction'],'discovery_route_step_fraction':replay['discovery_route_step_fraction']}
  if 'current_nine_loop_benchmark' in r:r['current_nine_loop_benchmark']['ten_loop_measurement_scope']='Two distinct fresh primes, learning and replay; full multi-prime reconstruction and exact certificates pending.'
 (D/'manifest.json').write_text(json.dumps(r,indent=2)+'\n')
 assert all(sha(D/p)==f['sha256'] and (D/p).stat().st_size==f['bytes'] for p,f in files.items())
 print(json.dumps({'passed':True,'files':len(files),'bundle_bytes':sum(f['bytes'] for f in files.values()),'manifest_sha256':sha(D/'manifest.json')},indent=2))
if __name__=='__main__':main()
