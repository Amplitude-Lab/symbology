#!/usr/bin/env python3
"""Read-only first-prime measurement review; never changes the active V27D run."""
from pathlib import Path
import hashlib,json,math,re,shutil,time
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009'
N=R/'output_hexagon_ten_v27d_20261009/runner'
D=R/'benchmarks/hexagon-mhv-ten-loop/evidence/v27-ten-prime1';D.mkdir(exist_ok=True)
bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def bind(p,expected=None):
 p=Path(p).resolve();h=sha(p);assert expected is None or h==expected,('hash',str(p));bindings[str(p)]=h;return h
def keep(p,name=None):
 p=Path(p);bind(p);dst=D/(name or p.name)
 if dst.exists():assert sha(dst)==sha(p),('archive collision',str(dst))
 else:shutil.copy2(p,dst)
def same(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
p=N/'audit/prime1-complete-comparison.json';c=json.loads(p.read_text());bind(p)
assert c['prime']==1152921504606847009 and c['candidate_version']=='V27D/16'
for src,h in c['sources'].items():bind(src,h)
paths=list(c['sources'])
new_report=next(Path(s) for s in paths if '/audit/p01-native-' in s)
new_log=next(Path(s) for s in paths if '/log/p01-native-' in s)
convert_report=next(Path(s) for s in paths if '/audit/p01-convert-' in s)
old_report=next(Path(s) for s in paths if '/audit/pindex01-native-' in s)
old_log=next(Path(s) for s in paths if '/log/pindex01-native-' in s)
for report,label in [(new_report,'candidate'),(old_report,'baseline')]:
 m=json.loads(report.read_text());assert m['exit_code']==0
 assert not any(m[k] for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
 same((m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,c[label]['native_cpu_hours'])
 same(m['wall_seconds'],c[label]['native_wall_seconds'])
 assert m['peak_summed_rss_bytes']==c[label]['native_peak_sampled_rss_bytes']
 keep(report)
m=json.loads(convert_report.read_text());assert m['exit_code']==0
same((m['cpu_user_seconds']+m['cpu_system_seconds'])/3600,c['candidate']['conversion_cpu_hours'])
same(m['wall_seconds'],c['candidate']['conversion_wall_seconds']);keep(convert_report)
for key,field in [('cpu_saving','native_cpu_hours'),('wall_saving','native_wall_seconds')]:
 same(c['comparison'][key],1-c['candidate'][field]/c['baseline'][field])
same(c['comparison']['rss_ratio'],c['candidate']['native_peak_sampled_rss_bytes']/c['baseline']['native_peak_sampled_rss_bytes'])
assert c['comparison']['cpu_saving']>0.15 and c['comparison']['rss_ratio']>1.05 and not c['comparison']['combined_gate_passed']
for kind,z in c['correctness']['coefficient_comparisons'].items():
 bind(z['path'],z['sha256']);bind(z['reference'],z['reference_sha256'])
 a=Path(z['path']).read_bytes();b=Path(z['reference']).read_bytes()
 assert len(a)==72+8*z['columns'] and a[:64]==b[:64] and a[72:]==b[72:]
 assert z['every_coefficient_equal']
t=new_log.read_text();assert 'V18_SEAM_VERIFIED full_operator=1 rows=989800 directions=0' in t
for op in ['collinear','mrk','ope_small_s']:assert 'JOINT_PHYSICAL_CHECK operator='+op+' passed=1' in t
parsed={};peaks=[]
for line in t.splitlines():
 if line.startswith('PROFILE phase='):
  z=dict(re.findall(r'(\w+)=([^ ]+)',line));name=z.pop('phase');parsed[name]={k:float(v) for k,v in z.items()};peaks.append(float(z.get('process_peak_bytes',0)))
for name,z in c['phases'].items():assert z==parsed[name]
assert max(peaks)==c['candidate']['native_kernel_peak_rss_bytes']
assert 'retained_nnz=1649199276 retained_alloc=1649199276 retained_array_bytes=19790391312' in t
manifest=json.loads((N/'manifest.json').read_text());assert manifest['memory_max_bytes']==50000000000 and not manifest['imports_old_coefficients']
assert manifest['profile_environment']['THREADS']==16 and manifest['profile_environment']['OPE_HEAD_CACHE_MIB']==3072
bind(manifest['driver'],manifest['driver_sha256'])
assert Path(json.loads(new_report.read_text())['command'][0])==Path(manifest['driver'])
old_manifest=R/'output_hexagon_17p3_20261007/ten-v16-runner/manifest.json';old=json.loads(old_manifest.read_text())
assert old['ope_head_cache_mib']==8192 and old['threads']==16 and old['mapped_capacity_gib']>0
assert 'OPE_HEAD_CACHE budget_bytes=8589934592' in old_log.read_text()
keep(old_manifest,'v16-manifest.json');keep(N/'manifest.json','v27d-frozen-run-manifest.json')
# Read only the small frozen source receipts, not live checkpoints or growing plans.
package=Path(manifest['driver']).parents[1];build=package/'v27-build-receipt.json';build_data=json.loads(build.read_text())
for name in ['src/shared/streaming_echelon.hpp','src/SparseRREF/sparse_mat.h']:
 bind(package/name,build_data['hashes'][name]);keep(package/name)
keep(build)
sample_path=N/'audit/heartbeat-20261009T172507Z.json';sample=json.loads(sample_path.read_text())['sample']
same(sample['sample_wall_seconds'],sample['end']['time']-sample['start']['time'])
same(sample['active_core_equivalents'],(sample['end']['cpu_seconds']-sample['start']['cpu_seconds'])/sample['sample_wall_seconds'])
assert 15.05<sample['active_core_equivalents']<15.07
keep(sample_path)
for source in [new_log,old_log,p]:keep(source)
convert_log=next(Path(s) for s in paths if '/log/p01-convert-' in s);keep(convert_log)
# Active status is a dated snapshot, not an immutable binding to the live file.
bind(D/'active-run-observation.json');bind(__file__)
out={'schema':'hexagon-v27d-ten-loop-prime1-review-v1','passed':True,'created':time.time(),'native_coordinates_checked':275246,'shifted_coordinates_checked':518048,'coefficient_equality_rechecked':True,'timing_memory_and_phase_receipts_recomputed':True,'complete_factored_seam_rows_mod_prime':989800,'new_QQ_reconstruction_complete':False,'new_exact_QQ_certificates_complete':False,'combined_CPU_RSS_gate_passed':False,'comparison':c,'scope':'One fresh ten-loop learning prime. Historical V16 mapped/8GiB OPE versus V27D heap/3GiB OPE. Finite-field checks and reference equality do not replace this new run\'s pending rational reconstruction and exact certificates. Replay timing remains pending at review. Existing certified ten-loop references remain unchanged.','optimization_hypothesis':'Count and allocate exact-capacity pivot-column incidence in the backsolve transpose; prove consumer coverage and validate generic finite-field, rational, affine, rank and exceptional cases before any replacement. No implementation change or competing benchmark made.','input_hashes':bindings}
dest=V/'v27-ten-prime1-review.json';dest.write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(dest,D/dest.name)
print(json.dumps({k:v for k,v in out.items() if k not in ['comparison','input_hashes']},indent=2));print('bindings',len(bindings))
