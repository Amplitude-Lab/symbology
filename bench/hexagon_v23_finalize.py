#!/usr/bin/env python3
"""Freeze the selected measured profile; do not modify production outputs."""
import argparse,fcntl,hashlib,json,statistics,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v23_20261009'
p=argparse.ArgumentParser();p.add_argument('--selected',required=True,choices=['V23','V26']);a=p.parse_args()
subprocess.run(['python3',str(R/'bench/hexagon_v23_report.py')],check=True)
d=json.loads((O/'RESULTS.json').read_text())
folder9='L9-v26' if a.selected=='V26' else 'L9-combined'
folder8='L8-v26' if a.selected=='V26' else 'L8-v23b-confirmation'
P=O/('package-v26' if a.selected=='V26' else 'package-v23b')
xs=[x for x in d['runs'] if x['suite']==folder9 and x['version']==a.selected]
lower=[x for x in d['runs'] if x['suite']==folder8 and x['version']==a.selected]
assert len(xs)==len(lower)==2 and len({x['prime'] for x in xs})==2
assert all(x['reference_coefficients_equal'] and x['all_physical_checks'] and x['full_seam_check'] and not x['fallbacks'] for x in xs+lower)
assert all(x['cgroup_events'].get('oom',0)==x['cgroup_events'].get('oom_kill',0)==0 for x in xs+lower)
old8_path=R/'output_hexagon_v18_20261009/validation/L8-V17-prime1.json';old8=json.loads(old8_path.read_text());old8cpu=(old8['cpu_user_seconds']+old8['cpu_system_seconds'])/3600
old9=d['baseline']['cpu_hours_mean'];old9rss=d['baseline']['conservative_peak_rss_bytes']
base22=next(g for g in d['groups'] if g['suite']=='L9-combined' and g['version']=='V22')
cpu=statistics.mean(x['cpu_hours'] for x in xs);rss=max(x['peak_rss_bytes'] for x in xs)
cpu8=statistics.mean(x['cpu_hours'] for x in lower);rss8=max(x['peak_rss_bytes'] for x in lower)
assert cpu<.85*old9 and rss<=1.05*old9rss
assert cpu8<.85*old8cpu and rss8<=1.05*old8['peak_summed_rss_bytes']
# The learning prime alone also clears the original gate; amortization does
# not hide a first-prime regression.
assert xs[0]['cpu_hours']<.85*old9 and lower[0]['cpu_hours']<.85*old8cpu
receipt=P/('v26-build-receipt.json' if a.selected=='V26' else 'v23-build-receipt.json')
build=json.loads(receipt.read_text())
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
for name,digest in build['hashes'].items():assert sha(P/name)==digest,name
refs=[]
for x in xs+lower:
 assert sha(x['reference'])==x['reference_sha256']
 root=O/x['suite']/x['version'];actual=root/'joint/cache'/('p'+str(x['prime'])+'.bin')
 got=actual.read_bytes();expected=Path(x['reference']).read_bytes()
 assert got[:64]==expected[:64] and got[72:]==expected[72:]
 assert sha(actual)==x['cache_sha256']
 refs.append({'candidate':str(actual),'candidate_sha256':sha(actual),'reference':x['reference'],'reference_sha256':x['reference_sha256'],'comparison':'all header/payload bytes equal except the algorithm/input identity word at bytes64:72'})
plans=[]
for folder in [folder8,folder9]:
 root=O/folder/a.selected
 for path in sorted((root/'discovery_plans').glob('*.tsv')):
  digest=sha(path);assert path.with_suffix(path.suffix+'.sha256').read_text()==digest+'\n'
  plans.append({'path':str(path),'sha256':digest,'bytes':path.stat().st_size,'learning_prime':1152921504606847009,'input_identity_sha256':sha(root/'joint/input_identity.txt')})
actual_qq=R/'output_hexagon_v19_20261009/reference-binding-audit.json'
reference_audit={'native_comparisons':refs,'prior_forced_QQ_binding_audit':str(actual_qq),'prior_forced_QQ_binding_audit_sha256':sha(actual_qq),'scope':'New native outputs match historical V17/V18 finite-field references. Nine-loop references are also bound transitively to the prior forced-QQ-only comparison recorded in the linked audit. No new QQ reconstruction was performed.'}
(O/'reference-binding-audit.json').write_text(json.dumps(reference_audit,indent=2)+'\n')
profile={'schema':'generic-replay-affine-tested-profile-v1','selected':a.selected,'package':str(P),'binary':str(P/'bin'/('joint_physics_'+a.selected.lower())),
 'binary_sha256':sha(P/'bin'/('joint_physics_'+a.selected.lower())),'build_receipt':str(receipt),
 'environment':{'THREADS':16,'PHYSICS_CACHE_MIB':4096,'OPE_HEAD_CACHE_MIB':3072,'HEXAGON_MHV_DOUBLE_BONUS':1,'HEXAGON_DISCOVERY_REPLAY':1,'HEXAGON_AFFINE_FACTORS':1,'MIMALLOC_PURGE_DELAY':0,'LD_LIBRARY_PATH':str(R/'output_hexagon_runtime_recovery_20261008/env/lib')},
 'unset_environment':['HEXAGON_MAPPED_FILE','HEXAGON_MAPPED_GIB','HEXAGON_MAPPED_SOFT_GIB'],
 'dependency_route_budget_bytes':67108864 if a.selected=='V26' else 0,
 'execution':'One prime process at a time, 16 workers. Use fresh output roots. Plans contain structure, never old field values.',
 'nine_loop':{'learning_cpu_hours':xs[0]['cpu_hours'],'replay_cpu_hours':xs[1]['cpu_hours'],'two_prime_total_cpu_hours':sum(x['cpu_hours'] for x in xs),'two_prime_mean_cpu_hours':cpu,'peak_rss_bytes':rss,'peak_sampled_rss_bytes':max(x['sampled_peak_rss_bytes'] for x in xs),'saving_vs_original_v17':1-cpu/old9,'saving_vs_matched_v22':1-cpu/base22['cpu_hours_mean'],'conservative_rss_ratio_vs_original_v17':rss/old9rss},
 'eight_loop':{'learning_cpu_hours':lower[0]['cpu_hours'],'replay_cpu_hours':lower[1]['cpu_hours'],'two_prime_mean_cpu_hours':cpu8,'peak_rss_bytes':rss8,'original_v17_cpu_hours':old8cpu,'original_v17_peak_rss_bytes':old8['peak_summed_rss_bytes'],'saving_vs_original_v17':1-cpu8/old8cpu,'rss_ratio_vs_original_v17':rss8/old8['peak_summed_rss_bytes'],'original_v17_measurement':str(old8_path)},
 'acceptance':{'baseline':'Original V17; matched V22 is a separate incremental comparator','cpu_saving_strictly_greater_than':.15,'rss_ratio_max':1.05,'passed':True,'rss_metric':d['rss_metric']},
 'scope':'Generic finite-field sparse/affine kernels, tested on arbitrary matrix families and two primes at both eight and nine loops. This does not establish universal speedups or global optimality. No ten-loop performance claim.',
 'production_modified':False,'production_memory_limit_bytes':50000000000,'remaining':['Ten-loop pilot under the existing 50 GB decimal limit','Dependency-route compression and cost-aware selection','Deterministic work-based early-physics scheduling','Fully implicit cross-stage affine projection; current checkpoint boundaries still materialize complete affine maps']}
(O/'accepted-nine-loop-profile.json').write_text(json.dumps(profile,indent=2)+'\n')
genericlog=O/('v26-generic.log' if a.selected=='V26' else 'generic-control-v23b.log');assert 'GENERIC_CONTROL passed=1' in genericlog.read_text()
micro=O/'v24-affine-families.log';assert micro.read_text().count('coefficient_equality=1 original_equations=1')==5
lines=[f"Selected {a.selected}: generic structural discovery replay and factored triangular affine maps.",
 f"Nine-loop mean: {cpu:.6f} CPU-hours per prime, including one learning prime and one replay prime. This is {100*(1-cpu/old9):.2f}% below original V17 and {100*(1-cpu/base22['cpu_hours_mean']):.2f}% below the fresh matched V22 pair. Conservative peak RSS is {rss/1e9:.3f} GB; the original-V17 reference for the 5% gate is {old9rss/1e9:.3f} GB. The learning prime alone also passes the CPU gate.",
 "The earlier V18-V22 structural replay covered the final seam only. Discovery now records source positions and pivot groups, regenerates all coefficients at each prime, and verifies every discarded equation against the new basis. Dependency routes accelerate reduction against retained pivots; absent or incomplete routes fall back locally to ordinary heap elimination. A changed source stream, pivot/rank or plan identity triggers a fresh discovery solve. All original final seam and physical checks remain mandatory.",
 "The affine factor represents triangular equations x_j = b_j - sum_{k>j} U_jk*x_k. Free-coordinate maps and one shared constant coordinate are composed directly, avoiding an intermediate full reduced matrix. The generic API supports lift, operator projection and materialization. The current driver uses materialization at its existing complete-affine checkpoint boundaries; it does not yet carry implicit factors through every physical projection stage.",
 "The core algorithms contain no amplitude, loop, alphabet, helicity, expected rank or known solution. The adapter supplies the existing equations and retains all normalization conventions and conditional bonus assumptions.",
 "Nine-loop measurements (CPU-hours; RSS uses the larger of sampled RSS and the kernel high-water mark):",
 "Variant/suite                  learning       replay       mean       peak GB"]
for g in d['groups']:
 if g['loop']!=9:continue
 ys=[x for x in d['runs'] if x['suite']==g['suite'] and x['version']==g['version']]
 vals=[f"{x['cpu_hours']:.6f}" for x in ys]
 lines.append(f"{g['version']+' / '+g['suite']:<30} {vals[0]:>9} {vals[1] if len(vals)>1 else 'n/a':>12} {g['cpu_hours_mean']:>10.6f} {g['peak_rss_bytes']/1e9:>12.3f}")
lines += [f"Original V17: {old9:.6f} CPU-hours, mean of two historical first-prime repetitions, 16 workers and 2 GiB OPE cache. The two-prime means above include learning plus replay; they are not same-prime statistical repetitions. All optimized runs use 3 GiB OPE cache. Small differences between neighboring candidates should not be treated as established universal speedups.",
 f"Eight-loop cross-instance check: {cpu8:.6f} CPU-hours per prime versus {old8cpu:.6f} for original V17 ({100*(1-cpu8/old8cpu):.2f}% saving). Peak RSS is {rss8/1e9:.3f} GB versus {old8['peak_summed_rss_bytes']/1e9:.3f} GB. Both fields match all 38,989 native reference coefficients; nine-loop outputs match all 105,898. Complete original seam sweeps and fresh collinear, MRK and OPE checks pass in every selected run.",
 "Generic controls passed 96 field cases, 48 successful structural replays and 72 consistent affine-map comparisons against independent dense Gaussian elimination. They include sparse/dense and deficient/inconsistent systems, two genuine exceptional-prime rank cases, identity and triangular-shape guards. The dependency variant also passed 48 deliberately emptied-route controls using exact local fallback. Five additional affine matrix families matched both the previous implementation and every original equation; their short timing measurements are supporting examples, not broad performance guarantees.",
 "The dominant measured discovery cost was incoming-row reduction (825 CPU seconds), not new pivot search (27 CPU seconds). This motivated dependency-path replay. A 256 MiB route cache was correct but raised eight-loop kernel RSS to 2.186 GB; reclaiming free allocator pages reduced that to 2.060 GB, still above the desired tolerance. The final tested route budget is 64 MiB for every input. Unused glibc pages are reclaimed only after temporary owners are destroyed at joined phase boundaries; no live coefficient/input object is removed. The 100 ms RSS sampler missed short allocation peaks, so the new acceptance audit also uses kernel high-water marks.",
 "Validation against the available rational result is transitive through the prior forced-QQ-only reference audit. New timings include input loading, discovery, physical fitting, checkpoint writes, structural-plan I/O and all current final checks. Compilation, cold basis/RHS generation and a new rational reconstruction are excluded. No ten-loop job was changed or benchmarked here; the 50 GB decimal production limit remains in force.",
 "Replay remains sensitive to the existing time-based physics schedule: a changed emitted row stream safely falls back, but can erase the speed advantage. Future generic directions are packed/delta-encoded dependency routes, selecting cached routes by measured work per byte, and work-based scheduling. Full cross-stage implicit maps are still untested. These results do not establish a globally optimal solver.",
 "Reproducible profile: "+str(O/'accepted-nine-loop-profile.json'),"All measurements: "+str(O/'RESULTS.json'),"Reference bindings: "+str(O/'reference-binding-audit.json'),"Frozen source/build receipt: "+str(receipt)]
(O/'REPORT.txt').write_text('\n\n'.join(lines)+'\n')
d['status']='complete';d['selected']=a.selected;d['scope_notes']=profile['scope'];(O/'RESULTS.json').write_text(json.dumps(d,indent=2)+'\n')
completion={'state':'complete','selected':a.selected,'finished':time.time(),'source_binary_files_verified':len(build['hashes']),'plans':plans,'generic_controls':str(genericlog),'generic_controls_sha256':sha(genericlog),'affine_family_controls':str(micro),'affine_family_controls_sha256':sha(micro),'native_reference_binding_audit':str(O/'reference-binding-audit.json'),'total_fresh_native_runs':len(d['runs']),'total_benchmark_cpu_hours':sum(x['cpu_hours'] for x in d['runs']),'production_modified':False}
(O/'completion-receipt.json').write_text(json.dumps(completion,indent=2)+'\n')
gate={'state':'complete_accepted_for_eight_and_nine_loop_measured_scope','selected':a.selected,'acceptance':profile['acceptance'],'profile':str(O/'accepted-nine-loop-profile.json'),'results':str(O/'RESULTS.json'),'completion':str(O/'completion-receipt.json'),'production_modified':False,'production_memory_limit_bytes':50000000000};(O/'execution-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
with (R/'docs/hexagon_optimization_lessons.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);p=R/'docs/hexagon_optimization_lessons.json';lessons=json.loads(p.read_text())
 entry={'id':'v23-v24-universal-replay-affine-factors','state':'complete','selected':a.selected,'finding':'Implemented generic guarded discovery pivot/group and bounded dependency-route replay plus triangular affine factors. Accepted the measured CPU/RSS profile on original V17 baseline at both8 and9 loops; no ten-loop performance inference. A 256 MiB dependency cache failed the lower-loop RSS requirement;64MiB plus joined-boundary reclamation is the tested policy. Kernel RSS peaks are required because100ms samples missed brief spikes.','profile':profile,'evidence':[str(O/'REPORT.txt'),str(O/'RESULTS.json'),str(O/'completion-receipt.json')],'failed_alternatives':['256MiB route cache without reclamation: eight-loop kernel RSS2.186GB','256MiB with glibc reclamation: eight-loop kernel RSS2.060GB'],'remaining_comparisons':profile['remaining']}
 lessons['lessons']=[x for x in lessons['lessons'] if x.get('id')!=entry['id']]+[entry];lessons['current_nine_loop_profile']=str(O/'accepted-nine-loop-profile.json');lessons['updated']='2026-10-09'
 lessons['scope']='Conditional ten-loop result remains complete and unchanged. '+a.selected+' generic replay/affine profile accepted in the measured eight/nine-loop scope; ten-loop speed and broad universal optimality remain unmeasured.'
 q=p.with_suffix('.tmp');q.write_text(json.dumps(lessons,indent=2)+'\n');q.replace(p)
print(json.dumps(profile['nine_loop'],indent=2));print(json.dumps(profile['eight_loop'],indent=2))
