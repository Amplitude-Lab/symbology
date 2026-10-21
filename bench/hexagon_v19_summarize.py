#!/usr/bin/env python3
"""Audit completed nine-loop controls and consolidate their measured scope."""
import argparse, hashlib, json, statistics, time
from pathlib import Path
R=Path(__file__).resolve().parents[1]
O=R/'output_hexagon_v19_20261009'
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--selected',default='V21/16')
a=p.parse_args()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
rows=[];groups={};audits=[]
for path in sorted(O.glob('controls*/progress.json')):
    d=json.loads(path.read_text())
    assert d['status']=='passed',(str(path),d['status'])
    for s in d['steps']:
        assert s['status']=='passed' and s['all_native_coefficients_equal']
        m=s['measurement'];assert m['exit_code']==0
        assert not any(m.get(k) for k in ['timeout','memory_pressure','rss_limit_exceeded','cancel_signal'])
        assert not any(s['cgroup_events'].values())
        cmd=m['command'];i=next(i for i,x in enumerate(cmd) if '/bin/joint_physics_' in x)
        binary=Path(cmd[i]);root=Path(cmd[i+6]);cache=root/'joint/cache'/f"p{s['prime']}.bin"
        if 'reference' in s:
            reference=Path(s['reference'])
        else:
            # The first runner revision recorded the fixed V17 reference hash
            # but not its path. Resolve that documented path and verify its hash.
            assert s['prime']==1152921504606847009
            reference=R/'output_hexagon_17p3_20261007/v17-nine-control/nine/joint/cache'/cache.name
        x=cache.read_bytes();y=reference.read_bytes()
        assert x[:64]==y[:64] and x[72:]==y[72:]
        assert sha(cache)==s['cache_sha256'] and sha(reference)==s['reference_sha256']
        logpath=path.parent/(s['label']+'.log')
        assert sha(logpath)==s['log_sha256']
        log=logpath.read_text()
        assert 'JOINT_PRIME_COMPLETE' in log
        for operator in ('collinear','mrk','ope_small_s'):
            assert 'JOINT_PHYSICAL_CHECK operator='+operator+' passed=1' in log
        assert ('V18_SEAM_VERIFIED full_operator=1 rows=372388 directions=0' in log or
                'JOINT_SEAM_COMPLETE full_operator=1 rows=372388' in log)
        version=s['label'].split('-')[0]
        assert sha(binary)==d['binaries'][version]
        row=dict(label=s['label'],version=version,threads=d['threads'],prime=s['prime'],
          cpu_hours=s['cpu_hours'],wall_seconds=m['wall_seconds'],peak_rss_bytes=m['peak_summed_rss_bytes'],
          cgroup_peak_bytes=s['cgroup_peak_bytes'],receipt=str(path),cache=str(cache),binary=str(binary),
          binary_sha256=d['binaries'][version],native_coefficients_equal=True,
          cgroup_limit_events_zero=True,phases=s['phase_profiles'])
        rows.append(row)
        if s['prime']==1152921504606847009:groups.setdefault(version+'/'+str(d['threads']),[]).append(row)
for package,receipt in [('package','v19-build-receipt.json'),('package-v20','v20-build-receipt.json'),('package-v21','v21-build-receipt.json'),('package-v22','v22-build-receipt.json')]:
    q=O/package/receipt
    if not q.exists():continue
    d=json.loads(q.read_text())
    for name,digest in d['hashes'].items():assert sha(q.parent/name)==digest,(package,name)
    audits.append(dict(receipt=str(q),sha256=sha(q),files_verified=len(d['hashes'])))
baseline=groups['V17/16'];best_current=groups['V21/16']
base_cpu=statistics.mean(x['cpu_hours'] for x in baseline)
current_cpu=statistics.mean(x['cpu_hours'] for x in best_current)
# A conservative common RSS gate: compare all candidate peaks to the LOWER
# peak of the two original V17 trials, rather than selecting a favorable pair.
rss_reference=min(x['peak_rss_bytes'] for x in baseline)
summary={}
for name,items in groups.items():
    cpu=statistics.mean(x['cpu_hours'] for x in items)
    all_fields=[x for x in rows if x['version']+'/'+str(x['threads'])==name]
    peak=max(x['peak_rss_bytes'] for x in all_fields)
    summary[name]=dict(first_prime_repeats=len(items),cpu_hours_mean=cpu,
      wall_seconds_mean=statistics.mean(x['wall_seconds'] for x in items),
      peak_rss_bytes_min=min(x['peak_rss_bytes'] for x in all_fields),peak_rss_bytes_max=peak,
      rss_scope='all tested fields and repetitions',
      cpu_saving_vs_original_v17=1-cpu/base_cpu,cpu_saving_vs_v21_16=1-cpu/current_cpu,
      conservative_rss_ratio=peak/rss_reference,
      numeric_gate_pass=cpu<.85*base_cpu and peak<=1.05*rss_reference)
    summary[name]['phase_cpu_seconds_mean']={phase:statistics.mean(
        next(p['cpu_seconds'] for p in x['phases'] if p['phase']==phase) for x in items)
        for phase in ['compute_discovery','compute_early_collinear','compute_early_mrk',
                      'compute_early_ope_small_s','compute_seam_complete']}
selected=summary[a.selected]
assert selected['numeric_gate_pass'] and selected['first_prime_repeats']>=2
sv,st=a.selected.split('/');chosen=[x for x in rows if x['version']==sv and x['threads']==int(st)]
assert {x['prime'] for x in chosen}=={1152921504606847009,1152921504606847067}
for index in [1,2]:
    d=json.loads((O/f'forced-qq-compare-p{index}.json').read_text());assert d['exit_code']==0
    log=(O/f'forced-qq-compare-p{index}.log').read_text();assert 'PASS affine spaces directions=0 columns=195155' in log
    # Bind selected outputs to the candidate that was converted to 14+4 and
    # compared to actual QQ-only references (whose directories contain no cache).
    prime=[1152921504606847009,1152921504606847067][index-1]
    v20=(O/'controls-v20/V20/joint/cache'/f'p{prime}.bin').read_bytes()
    for x in chosen:
        if x['prime']==prime:
            z=Path(x['cache']).read_bytes();assert z[:64]==v20[:64] and z[72:]==v20[72:]
assert not list((O/'qq-only-references/L9').rglob('cache'))
result=dict(schema='nine-loop-cpu-optimization-summary-v1',updated=time.time(),
 acceptance=dict(baseline='original V17,16 threads,OPE2GiB',metric='process user+system CPU-hours',
   cpu_saving_strictly_greater_than=.15,rss_ratio_max=1.05,conservative_rss_reference_bytes=rss_reference),
 selected=a.selected,groups=summary,runs=rows,source_audits=audits,
 validation=dict(native_columns=105898,shifted_QQ_reference_columns=195155,
   primes=[1152921504606847009,1152921504606847067],all_coefficients_equal=True,
   all_current_fresh_physical_and_full_seam_checks=True,forced_QQ_reduction=True,
   new_QQ_reconstruction=False),
 scope='Serialized complete fresh primes on warm precomputed basis inputs; compilation, cold basis generation and new rational reconstruction excluded. RSS is sampled at 100ms. Shared workstation; no confidence interval or global-optimality claim.',
 production_modified=False,production_memory_limit_bytes=50000000000)
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
profile=dict(schema='accepted-nine-loop-profile-v1',selection=a.selected,
 binary=chosen[0]['binary'],binary_sha256=chosen[0]['binary_sha256'],
 loop=9,head=15,tail=3,target=6000,mode='cache',prereduce='right',
 inputs=dict(fec_root=str(R/'output_hexagon_16p4_20261006/FEC'),
   lec_root=str(R/'output_hexagon_16p4_20261006/LEC/MHV'),
   data_root=str(R/'output_hexagon_16p4_20261006/data')),
 environment=dict(THREADS=int(st),PHYSICS_CACHE_MIB=4096,OPE_HEAD_CACHE_MIB=3072,
    HEXAGON_MHV_DOUBLE_BONUS=1,MIMALLOC_PURGE_DELAY=0,
    LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib')),
 unset_environment=['HEXAGON_MAPPED_FILE','HEXAGON_MAPPED_GIB','HEXAGON_MAPPED_SOFT_GIB'],
 limits=dict(memory_max_bytes=20000000000,memory_high_bytes=18000000000,swap_bytes=0,
    virtual_GiB=96,sampled_RSS_limit_bytes=19000000000,tasks=128),
 results=str(O/'RESULTS.json'),production_memory_limit_bytes=50000000000,
 acceptance=selected,scope='Nine-loop measured profile only; use fresh output roots. Ten-loop pilot remains required.')
(O/'accepted-nine-loop-profile.json').write_text(json.dumps(profile,indent=2)+'\n')
lines=['Nine-loop solver optimization — measured CPU-hours, 2026-10-09','',
 f"Selected measured profile: {a.selected}. Mean CPU saving against original V17: {selected['cpu_saving_vs_original_v17']:.2%}; against V21/16: {selected['cpu_saving_vs_v21_16']:.2%}.",
 'Acceptance: strictly more than 15% fewer user+system CPU-hours than original V17, with at most 5% higher peak RSS. Candidate peaks are conservatively compared to the LOWER of the two original V17 peaks.',
 '', 'Version/threads | p1 repetitions | mean CPU-hours | mean wall seconds | peak RSS range (decimal GB) | CPU saving vs original V17']
for name,s in summary.items():
    lines.append(f"{name:14} | {s['first_prime_repeats']:2} | {s['cpu_hours_mean']:.6f} | {s['wall_seconds_mean']:.3f} | {s['peak_rss_bytes_min']/1e9:.3f}–{s['peak_rss_bytes_max']/1e9:.3f} | {s['cpu_saving_vs_original_v17']:.2%}")
lines+=['',
 'The p2 runs are second-field correctness/replay controls, not repetitions of the p1 timing; the RSS ranges include every tested field. The selected profile has two complete p1 timings plus a p2 control. V17 uses its original algorithm/settings in the matched build with phase instrumentation. All numerical timings were serialized, without concurrent builds or other numerical probes. Build time, cold basis generation and new rational reconstruction are excluded. Shared-workstation timing has no confidence interval; these tests do not establish global optimality.',
 '', 'What changed and why:',
 'V18 uses direct final-seam contraction, forward continuation after restriction and guarded source-row replay across primes. A 3 GiB OPE-head cache holds all three nine-loop terminal heads; V17 used 2 GiB and rebuilt one during validation.',
 'V19 extends exact collinear suffix projection to the post-discovery fit. Its 12.29% CPU saving in the first matched trial fails the requested 15% gate.',
 'V20 compresses the affine suffix operator, including every sampled RHS column. It selects independent ORIGINAL suffix rows and explicitly verifies source row-span certificates. The 36 suffix rows reduce to 14, removing repeated large projected sparse products.',
 'V21 additionally checks every relevant head copy at every current sample for exact zero, separately in each field and batch. Its exact support mask reduces 36 suffix rows to nine in this case. RHS-only constraints remain. Final scalar physical checks are unchanged.',
 'The V21 eight-thread control passes correctness but costs more CPU-hours than sixteen threads. Wall-time-based early-physics scheduling changes its discovery space (5,962 directions/48,535,187 nonzeros versus 5,977/47,127,363), so this is a complete configuration comparison rather than isolated hardware scaling.',
 'V22 stores only a+b<=w logarithm slots at intermediate head weight w. It keeps the full T/S jet and all kernels, and restores the existing final/cache layout before projection. At head weight 15, 136 logarithm slots replace 190 final-weight slots. Two-prime square-reference controls cover all letters and weights 1–20 (2,160 cases), plus complete FEC8 heads (1,178,100 coefficients per field). No new physical relation is imposed.',
 'The V22 build also changes the shared recurrence step to accept separate input/output strides. Collinear timing improved in both V22 trials even though its mathematical operator is unchanged; no separate compiler/layout attribution experiment was performed. Phase CPU means are reported independently in RESULTS.json, and the full improvement must not all be assigned to OPE storage.',
 '', 'Correctness and provenance:',
 'Every completed native output has all 105,898 coefficients equal to the baseline reference at its prime (only the producer fingerprint differs). Fresh COL, MRK and OPE checks and the complete 372,388-row factored seam are retained. All cgroup high/max/OOM counters are zero.',
 'V20 was converted from 15+3 to 14+4 with exact shift-back verification. At both primes all 195,155 shifted coordinates match reductions of the existing rational nine-loop solution. The reference snapshots deliberately contain no cache directory, forcing actual QQ reduction. Selected candidate payloads match V20 at both primes. This repairs a validation-selection caveat: earlier pipeline_compare calls could choose saved Fp caches when present.',
 'The final source/binary hash audit and per-run cache/reference/log hash checks are recorded in RESULTS.json. Synthetic suffix tests preserve RHS-only constraints, changing affine rank, zero operators and exact head-zero support. Existing production inputs, outputs and checkpoints are preserved.',
 '', 'Mathematical strategy:',
 'The current four double-final-entry bonus proof is conditional on its heptagon/collinear/adjacency/Qbar hypotheses at ordinary-symbol level. This optimization does not strengthen those assumptions. The complete published triple Eq. (92) and cyclic images add three independent local constraints (65 to 62 triple directions) and pass exact known eight/nine-loop checks. The certified local source alone does not derive them; rational local witnesses certify this non-implication, without being physical-amplitude counterexamples. New triple constraints were not used in the numerical candidates.',
 'The original truncated Eq. (92) diagnostic was invalid and has been marked accordingly; corrected v3 checks pass. A separate quadruple-dimension calibration returned 169 instead of the published 166 and remains unresolved. See bonusrelation/solver_audit_20261009/REPORT.txt for the precise scope and certificates.',
 '', 'Remaining useful comparisons:',
 'Persistent discovery row/pivot replay with numerical coefficient/RHS recomputation, exceptional-rank guards and full verification; precomposition of the proven double-bonus tail copy basis; deterministic early-physics schedules; alternative splits with certified compressed suffix operators; a bounded ten-loop pilot. These remain proposals, not measured improvements.',
 '', 'Use accepted-nine-loop-profile.json for the exact binary, hash, resource settings and measured scope. Sixteen-thread controls use heap allocation, 4 GiB physical cache, 3 GiB OPE cache (except original V17: 2 GiB), 20 GB hard cgroup limit, 18 GB high threshold, no swap, 96 GiB virtual limit. The production ten-loop limit remains 50,000,000,000 decimal bytes. No new ten-loop candidate has been selected or deployed.',
 '', 'Reproduce from the workspace root: python3 bench/hexagon_v19_controls.py --package output_hexagon_v19_20261009/package-v22 --output /path/to/new-output-directory --threads 16 --jobs V22-p1,V22-p2,V22-p1-r2. Replace the output placeholder with a NEW directory. Do not reuse a completed output directory as a fresh timing.',
 '', 'Related reports: output_hexagon_algorithm_comparison_20261008/COMPARISON.txt; output_hexagon_v18_20261009/REPORT.txt; bonusrelation/solver_audit_20261009/REPORT.txt.',
 'Primary higher-entry source: https://arxiv.org/html/2308.08199 (complete Eq. 92).']
(O/'REPORT.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(selected=a.selected,groups=summary,verified_runs=len(rows),source_audits=audits),indent=2))
