#!/usr/bin/env python3
"""Record a completed, independently audited optimization without changing production."""
import fcntl, hashlib, json, subprocess, time
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'output_hexagon_v19_20261009'
d=json.loads((O/'RESULTS.json').read_text());selected=d['groups'][d['selected']]
assert selected['numeric_gate_pass'] and selected['first_prime_repeats']>=2
assert d['validation']['all_coefficients_equal'] and d['validation']['forced_QQ_reduction']
control=json.loads((O/'ope-weight-control.json').read_text());assert control['exit_code']==0
assert 'OPE_WEIGHT_CONTROL PASS cases=2160 guards=12 complete_heads=2 square_reference=1' in (O/'ope-weight-control.log').read_text()
checks={}
for name in ['independent-verification.json','native-triple-independent-verification.json','completed-ten-triple-receipt.json']:
    q=R/'bonusrelation/solver_audit_20261009'/name;x=json.loads(q.read_text());assert x['passed']
    checks[str(q)]=hashlib.sha256(q.read_bytes()).hexdigest()
gate=dict(state='complete_accepted_for_nine_loop_measured_scope',selected=d['selected'],
 baseline=d['acceptance']['baseline'],metric=d['acceptance']['metric'],acceptance=d['acceptance'],
 measurements=selected,results=str(O/'RESULTS.json'),profile=str(O/'accepted-nine-loop-profile.json'),
 source_freeze=str(O/'package-v22/v22-build-receipt.json'),
 proof_scope='Conditional four double-final-entry ordinary-symbol bonuses; no additional triple bonus imposed.',
 current_progress=[str(p) for p in sorted(O.glob('controls*/progress.json'))],
 production_modified=False,production_memory_limit_bytes=50000000000,
 remaining=['bounded ten-loop pilot','persistent discovery replay','precomposed bonus copy basis','deterministic early-physics schedule'],
 completed=time.time())
(O/'execution-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
entry=dict(id='nine-loop-affine-suffix-and-weight-dependent-ope-cpu-gate',
 finding=f"{d['selected']} uses {selected['cpu_hours_mean']:.6f} mean CPU-hours: {selected['cpu_saving_vs_original_v17']:.2%} below original V17, {selected['cpu_saving_vs_v21_16']:.2%} below V21/16. All selected RSS peaks satisfy the conservative +5% gate. Complete native outputs and physical/seam checks pass.",
 practice='Profile full process user+system CPU, not only wall time or ranks. Compress the affine suffix operator with every RHS column and explicit row-span certificates; recompute exact zero-head support per field/batch. During OPE head recurrence store only a+b<=current weight, then restore the external/cache layout. Preserve complete residual and physical checks. Fewer workers can increase CPU-hours and can alter wall-time-driven discovery scheduling.',
 measurements=d['groups'],
 unsuccessful_alternatives=['Raw V19 exact suffix projection misses the 15% CPU gate.','V21 at eight threads passes correctness but is slower in CPU and wall time than V21/16; its early-physics admission schedule and discovery fill differ.'],
 validation=dict(**d['validation'],ope_square_reference_cases=2160,ope_complete_head_fields=2,
   ope_complete_head_coefficients_per_field=1178100,source_audits=d['source_audits'],
   prior_QQ_comparison_caveat='Earlier pipeline_compare used available Fp cache when present; current snapshots omit cache directories and force actual rational reduction.'),
 resource_settings=json.loads((O/'accepted-nine-loop-profile.json').read_text()),
 higher_bonus_scope='Three additional full Eq92 triple rows have exact rank3, reducing65 to62, and pass all known L8/L9 QQ contractions. Certified local source does not imply them; exact local witnesses are not physical amplitude counterexamples. Ten-loop completed caches pass at two primes a posteriori. No new all-loop triple theorem or new constraints in solver. Quadruple hierarchy calibration169 vs published166 remains unresolved.',
 evidence=[str(O/'REPORT.txt'),str(O/'RESULTS.json'),str(R/'bonusrelation/solver_audit_20261009/REPORT.txt')],
 scope=d['scope'],production_modified=False,production_memory_limit_bytes=50000000000,
 remaining_comparisons=gate['remaining'])
lock=R/'docs/hexagon_optimization_lessons.lock'
with lock.open('a') as f:
    fcntl.flock(f,fcntl.LOCK_EX)
    p=R/'docs/hexagon_optimization_lessons.json';lessons=json.loads(p.read_text())
    assert not any(x.get('id')==entry['id'] for x in lessons['lessons']),'Already finalized; inspect before updating'
    lessons['lessons'].append(entry);lessons['updated']='2026-10-09'
    temp=p.with_suffix('.json.tmp');temp.write_text(json.dumps(lessons,indent=2,ensure_ascii=False)+'\n');temp.replace(p)
files=['RESULTS.json','REPORT.txt','accepted-nine-loop-profile.json','execution-gate.json','ope-weight-control.json','ope-weight-control.log']
receipt=dict(status='passed',selected=d['selected'],completed=time.time(),
  artifacts={str(O/n):hashlib.sha256((O/n).read_bytes()).hexdigest() for n in files},
  mathematical_receipts=checks,production_modified=False,ten_loop_candidate_performance_measured=False)
profile=json.loads((O/'accepted-nine-loop-profile.json').read_text())
linked=subprocess.run(['ldd',profile['binary']],check=True,capture_output=True,text=True).stdout
libraries={}
for line in linked.splitlines():
    for token in line.split():
        if token.startswith('/') and Path(token).is_file():
            q=Path(token).resolve();libraries[str(q)]=hashlib.sha256(q.read_bytes()).hexdigest()
receipt['dynamic_libraries_sha256']=libraries
hardware=subprocess.run(['lscpu','--json'],check=True,capture_output=True,text=True)
receipt['hardware']=json.loads(hardware.stdout)
(O/'completion-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
handoff=f'''The completed nine-loop comparison selects {d['selected']}. Mean CPU-hours {selected['cpu_hours_mean']:.6f}, saving {selected['cpu_saving_vs_original_v17']:.2%} against original V17 and {selected['cpu_saving_vs_v21_16']:.2%} against V21/16. All selected RSS peaks meet the user's +5% cap. Exact affine suffix compression and weight-dependent OPE head storage preserve all native coefficients and fresh physical/full seam checks. Two-prime actual QQ reductions are bound through matching V20 native payloads. The eight-thread alternative is slower than sixteen threads.

Please treat this as a completed result for review, not a request to interrupt physical checks or launch a ten-loop run. No production modification or new ten-loop performance claim. The production cap remains 50,000,000,000 bytes.

Report: {O/'REPORT.txt'}
Machine-readable results: {O/'RESULTS.json'}
Pinned profile: {O/'accepted-nine-loop-profile.json'}
Source freeze: {O/'package-v22/v22-build-receipt.json'}
Reproducers: {R/'bench/hexagon_joint_build_v22.py'}, {R/'bench/hexagon_v19_controls.py'}
Math audit: {R/'bonusrelation/solver_audit_20261009/REPORT.txt'}

Full published Eq92 and its cyclic images add rank3 beyond the four double bonuses and pass exact L8/L9 and completed L10 two-prime tests. No new all-loop triple proof: independent rational local witnesses show the currently certified local source does not imply these targets. They were not imposed in this optimization. Earlier incomplete Eq92 probes are marked invalid; full 21-term v3 probes pass. An optional quadruple hierarchy169 vs published166 remains an unresolved diagnostic.
'''
(O/'HANDOFF.txt').write_text(handoff)
print(json.dumps(receipt,indent=2))
