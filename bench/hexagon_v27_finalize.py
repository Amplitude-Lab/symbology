#!/usr/bin/env python3
"""Audit the completed V27 campaign and retain V26 if the new gate fails."""
import fcntl
import hashlib
import json
import math
import statistics
import subprocess
import time
from pathlib import Path

R = Path(__file__).resolve().parents[1]
O = R / 'output_hexagon_v27_20261009'
P = O / 'package-v27d'


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


assert read(O / 'final-comparisons/progress.json')['state'] == 'passed'
subprocess.run(['python3', str(R / 'bench/hexagon_v27_report.py')], check=True)
d = read(O / 'RESULTS.json')
bindings = {}


def bind(path, expected=None):
    path = Path(path).resolve()
    digest = sha(path)
    assert expected is None or digest == expected, ('hash mismatch', str(path))
    bindings[str(path)] = digest
    return digest


builds = []
for revision in 'abcd':
    package = O / ('package-v27' + revision)
    receipt = package / 'v27-build-receipt.json'
    build = read(receipt)
    for name, digest in build['hashes'].items():
        bind(package / name, digest)
    bind(receipt)
    builds.append(dict(package=str(package), files_verified=len(build['hashes']),
                       receipt=str(receipt), receipt_sha256=sha(receipt)))

native = []
roots = set()
for run in d['runs']:
    folder = O / run['suite']
    progress = read(folder / 'progress.json')
    assert progress['status'] == 'passed'
    step = next(x for x in progress['steps'] if x['label'] == run['label'])
    measurement = read(run['measurement'])
    assert step['status'] == 'passed' and step['measurement'] == measurement
    assert measurement['exit_code'] == 0
    assert not any(measurement[k] for k in
                   ['timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal'])
    assert math.isclose((measurement['cpu_user_seconds'] +
                         measurement['cpu_system_seconds']) / 3600,
                        run['cpu_hours'], rel_tol=1e-12)
    assert run['wall_seconds'] == measurement['wall_seconds']
    assert run['peak_rss_bytes'] == max(measurement['peak_summed_rss_bytes'],
                                       max(x.get('process_peak_bytes', 0)
                                           for x in step['phase_profiles']))
    assert run['phases'] == step['phase_profiles']
    assert run['cgroup_events'] == step['cgroup_events']
    assert run['cgroup_events'] and all(v == 0 for v in run['cgroup_events'].values())
    bind(run['log'], step['log_sha256'])
    bind(run['measurement'])
    bind(folder / 'progress.json')
    command = measurement['command']
    i = next(i for i, x in enumerate(command)
             if Path(x).name == 'joint_physics_' + run['version'].lower())
    bind(command[i], run['binary_sha256'])
    root = Path(command[i + 6])
    roots.add(root)
    cache = root / 'joint/cache' / ('p' + str(run['prime']) + '.bin')
    bind(cache, run['cache_sha256'])
    bind(run['reference'], run['reference_sha256'])
    actual = cache.read_bytes()
    reference = Path(run['reference']).read_bytes()
    assert actual[:64] == reference[:64] and actual[72:] == reference[72:]
    assert run['native_coefficients_equal'] and run['all_physical_checks']
    assert run['complete_seam_check']
    log = Path(run['log']).read_text()
    assert 'JOINT_PRIME_COMPLETE' in log and 'V18_SEAM_VERIFIED full_operator=1' in log
    for operator in ['collinear', 'mrk', 'ope_small_s']:
        assert 'JOINT_PHYSICAL_CHECK operator=' + operator + ' passed=1' in log
    native.append(dict(candidate=str(cache), candidate_sha256=sha(cache),
                       reference=run['reference'], reference_sha256=sha(run['reference']),
                       comparison='All header and coefficient bytes equal except identity at bytes 64:72'))

plans = []
for root in sorted(roots):
    bind(root / 'joint/input_identity.txt')
    for plan in sorted((root / 'discovery_plans').glob('*.tsv')):
        digest = bind(plan)
        assert plan.with_suffix(plan.suffix + '.sha256').read_text() == digest + '\n'
        bind(plan.with_suffix(plan.suffix + '.sha256'))
        plans.append(dict(path=str(plan), sha256=digest, bytes=plan.stat().st_size))

for name in ['generic-new', 'generic-legacy', 'additional-controls']:
    log = O / 'suite-d' / (name + '.log')
    bind(log)
    expected = ('ADDITIONAL_CONTROL passed=1' if name == 'additional-controls'
                else 'GENERIC_CONTROL passed=1 field_cases=96 successful_replays=48 guarded_fallbacks=2 factor_cases=72')
    assert expected in log.read_text()
assert 'WIDE_CONTROL overflow_cases=5077 canonical_reference_equal=1' in (O / 'suite-d/additional-controls.log').read_text()
assert read(O / 'suite-d/progress.json')['state'] == 'passed'
bind(O / 'suite-d/progress.json')

prior = R / 'output_hexagon_v23_20261009/reference-binding-audit.json'
prior_data = read(prior)
bind(prior)
bind(prior_data['prior_forced_QQ_binding_audit'], prior_data['prior_forced_QQ_binding_audit_sha256'])
reference_audit = dict(native_comparisons=native, prior_reference_audit=str(prior),
                       prior_reference_audit_sha256=sha(prior),
                       scope='Every new finite-field output was rechecked against the historical V17/V18 native references. The prior audit binds the nine-loop references to a forced-QQ-only comparison. No new QQ reconstruction or independent new physics proof was performed.')
write(O / 'reference-binding-audit.json', reference_audit)


def group(suite):
    return next(x for x in d['groups'] if x['suite'] == suite)


def runs(suite):
    return [x for x in d['runs'] if x['suite'] == suite]


for g in d['groups']:
    rows = runs(g['suite'])
    assert len(rows) == g['runs'] and g['correctness_passed']
    assert statistics.mean(x['cpu_hours'] for x in rows) == g['cpu_hours_mean']
    assert max(x['peak_rss_bytes'] for x in rows) == g['peak_rss_bytes']
    g['wall_seconds_mean'] = statistics.mean(x['wall_seconds'] for x in rows)

baseline = group('suite-a/L9-V26')
baseline8 = group('suite-a/L8-V26')
accepted_v26_path = R / 'output_hexagon_v23_20261009/accepted-nine-loop-profile.json'
accepted_v26 = read(accepted_v26_path)
assert accepted_v26['selected'] == 'V26'
bind(accepted_v26_path)
bind(accepted_v26['binary'], accepted_v26['binary_sha256'])
assert all(x['binary_sha256'] == accepted_v26['binary_sha256']
           for x in d['runs'] if x['version'] == 'V26')
candidate = group('suite-d/L9-V27')
candidate8 = group('suite-d/L8-V27')
worker8 = group('final-comparisons/L9-eight-workers')
worker8_baseline_first = runs('final-comparisons/L9-V26-eight-workers')[0]
worker8_candidate_first = runs('final-comparisons/L9-eight-workers')[0]
matched_worker8 = dict(
    scope='Matched first prime only; eight workers for both executables',
    baseline_cpu_hours=worker8_baseline_first['cpu_hours'],
    candidate_cpu_hours=worker8_candidate_first['cpu_hours'],
    cpu_saving=1-worker8_candidate_first['cpu_hours']/worker8_baseline_first['cpu_hours'],
    baseline_peak_rss_bytes=worker8_baseline_first['peak_rss_bytes'],
    candidate_peak_rss_bytes=worker8_candidate_first['peak_rss_bytes'])


def compare(g, b):
    return dict(cpu_saving=1 - g['cpu_hours_mean'] / b['cpu_hours_mean'],
                rss_ratio=g['peak_rss_bytes'] / b['peak_rss_bytes'],
                passed=(g['cpu_hours_mean'] < .85 * b['cpu_hours_mean'] and
                        g['peak_rss_bytes'] <= 1.05 * b['peak_rss_bytes']))


paired9 = [g for g in d['groups'] if g['loop'] == 9 and g['runs'] == 2 and g['version'] == 'V27']
for g in paired9:
    assert not compare(g, baseline)['passed'], 'An eligible new profile needs an explicit promotion review'
# Do not silently dismiss a potentially successful first-prime ablation.
baseline_first = runs('suite-a/L9-V26')[0]
for g in d['groups']:
    if g['loop'] == 9 and g['version'] == 'V27' and g['runs'] == 1:
        assert not (g['cpu_hours_mean'] < .85 * baseline_first['cpu_hours'] and
                    g['peak_rss_bytes'] <= 1.05 * baseline['peak_rss_bytes']), 'Complete second-prime and eight-loop validation of the promising ablation'

ablations = []
full_first = runs('suite-d/L9-V27')[0]
for flag in ['kernels', 'compression', 'schedule', 'routes']:
    suite = 'final-comparisons/L9-without-' + flag
    off = runs(suite)[0]
    ablations.append(dict(feature=flag, off_suite=suite, matching_prime=off['prime'],
                          full_cpu_hours=full_first['cpu_hours'], off_cpu_hours=off['cpu_hours'],
                          enabling_feature_cpu_saving=1-full_first['cpu_hours']/off['cpu_hours'],
                          off_peak_rss_bytes=off['peak_rss_bytes'],
                          scope='One-at-a-time ablation; interactions remain. Disabled features use their V26 implementation.'))
legacy_routes = group('final-comparisons/L9-without-routes')
route_pair = dict(full_cpu_hours_mean=candidate['cpu_hours_mean'],
                  legacy_routes_cpu_hours_mean=legacy_routes['cpu_hours_mean'],
                  enabling_packed_routes_cpu_saving=1-candidate['cpu_hours_mean']/legacy_routes['cpu_hours_mean'],
                  scope='Both first-prime learning and second-prime replay included, 16 workers')

profile = dict(schema='hexagon-v27-experimental-profile-v1', state='tested_not_accepted',
               package=str(P), binary=str(P / 'bin/joint_physics_v27'),
               binary_sha256=sha(P / 'bin/joint_physics_v27'),
               build_receipt=str(P / 'v27-build-receipt.json'),
               baseline=baseline, nine_loop=candidate, eight_loop=candidate8,
               nine_loop_eight_workers=worker8,
               matched_eight_worker_first_prime=matched_worker8,
               comparison_nine=compare(candidate, baseline),
               comparison_nine_eight_workers=compare(worker8, baseline),
               comparison_eight=compare(candidate8, baseline8),
               acceptance=dict(baseline='Fresh V26, 16 workers, 3 GiB OPE cache',
                               cpu_saving_required_strictly_greater_than=.15, rss_ratio_max=1.05,
                               passed=False),
               environment=dict(THREADS=16, PHYSICS_CACHE_MIB=4096, OPE_HEAD_CACHE_MIB=3072,
                                HEXAGON_MHV_DOUBLE_BONUS=1, HEXAGON_DISCOVERY_REPLAY=1,
                                HEXAGON_AFFINE_FACTORS=1, MIMALLOC_PURGE_DELAY=0,
                                SOLVER_WORK_KERNELS=1, SOLVER_WORK_SCHEDULE=1,
                                SOLVER_EXACT_COMPRESSION=1, SOLVER_PACKED_ROUTES=1,
                                LD_LIBRARY_PATH=str(R / 'output_hexagon_runtime_recovery_20261008/env/lib')),
               hardware=dict(cpu_model='Intel Core Ultra 7 265', logical_cpus=20,
                             perf_event_paranoid=4, hardware_counters_available=False),
               unset_environment=['HEXAGON_MAPPED_FILE','HEXAGON_MAPPED_GIB','HEXAGON_MAPPED_SOFT_GIB'],
               ablations=ablations, route_pair=route_pair,
               current_nine_loop_default='V26', production_modified=False,
               production_memory_limit_bytes=50000000000)
write(O / 'experimental-profile.json', profile)

lines = [
    'All four directions implemented and tested in experimental V27; V26 remains the accepted nine-loop solver.',
    f"The latest complete V27D profile averages {candidate['cpu_hours_mean']:.6f} CPU-hours per prime versus {baseline['cpu_hours_mean']:.6f} for fresh V26: {100*compare(candidate,baseline)['cpu_saving']:.2f}% less CPU. Peak RSS is {candidate['peak_rss_bytes']/1e9:.3f} versus {baseline['peak_rss_bytes']/1e9:.3f} decimal GB. With eight workers, V27D averages {worker8['cpu_hours_mean']:.6f} CPU-hours and peaks at {worker8['peak_rss_bytes']/1e9:.3f} GB. Neither passes the required strictly-greater-than-15% CPU saving, although both satisfy the 5% RSS tolerance. The new benchmark is V26, as explicitly requested; earlier savings versus V17 are not credited again.",
    'Measurements below use total process user+system CPU divided by 3600. Two-prime means combine a learning prime and a replay prime, not independent statistical repetitions. RSS is the larger of the 100 ms process sample and kernel process high-water mark. One native solve ran at a time; compilation and other numerical tests did not overlap these timings.',
    f"Matched eight-worker first-prime control: V26 {matched_worker8['baseline_cpu_hours']:.6f}, V27D {matched_worker8['candidate_cpu_hours']:.6f} CPU-hours ({100*matched_worker8['cpu_saving']:.2f}% saving). Peak RSS is {matched_worker8['baseline_peak_rss_bytes']/1e9:.3f} and {matched_worker8['candidate_peak_rss_bytes']/1e9:.3f} GB, respectively. This is a single-prime comparison, separate from the main two-prime/16-worker baseline.",
    'Nine-loop measurements: suite, workers, primes, mean CPU-hours, peak GB, mean wall seconds',
]
for g in d['groups']:
    if g['loop'] == 9:
        lines.append(f"{g['suite']:<46} {g['threads']:>2} {g['runs']:>2} {g['cpu_hours_mean']:>10.6f} {g['peak_rss_bytes']/1e9:>8.3f} {g['wall_seconds_mean']:>9.1f}")
lines += [
    'Ablations, same first prime and 16 workers. Positive values mean enabling the named feature reduced CPU; negative values mean it increased CPU. These are local effects with other features enabled, and must not be added together.',
]
for a in ablations:
    lines.append(f"{a['feature']:<14} enabled {a['full_cpu_hours']:.6f}; disabled {a['off_cpu_hours']:.6f}; enabling saving {100*a['enabling_feature_cpu_saving']:+.2f}%")
lines.append(f"Packed routes, learning plus replay: {route_pair['full_cpu_hours_mean']:.6f} CPU-hours/prime versus legacy routes {route_pair['legacy_routes_cpu_hours_mean']:.6f}; enabling packed routes saves {100*route_pair['enabling_packed_routes_cpu_saving']:.2f}% in this pair.")
phase_names = ['compute_discovery', 'compute_early_collinear', 'compute_early_mrk',
               'compute_early_ope_small_s', 'compute_seam_complete']
lines.append('First-prime phase CPU seconds (phase timers are nested; do not sum all PROFILE records): V26/16, V27D/16, V27D/8')
for phase in phase_names:
    values = [next(p['cpu_seconds'] for p in x['phases'] if p['phase'] == phase)
              for x in [baseline_first, full_first, worker8_candidate_first]]
    lines.append(f"{phase:<30} " + ' '.join(f'{v:>10.3f}' for v in values))
lines.append('Discovery illustrates the CPU/wall tradeoff: V26/16 used 930.056 CPU seconds in 126.644 wall seconds; V27D/16 used 1112.65 CPU seconds in 104.756 wall seconds. V27D/8 used 868.333 CPU seconds in 140.247 wall seconds. Greater worker utilization and shorter elapsed time do not by themselves meet the requested CPU-hour objective. A phase-dependent worker policy remains an untested follow-up; these observations do not establish the hardware-level cause of the scaling difference.')
lines += [
    'Numerical kernels: exact hierarchical bit queue, work-balanced row ranges, reusable per-worker projection scratch, exact density-based sparse/dense dispatch, and delayed 128-bit modular accumulation with exact overflow correction. Fields larger than 61 bits use the existing canonical reduction path. Dependency hints always finish by checking all residual pivots. The arithmetic core does not depend on amplitude, loop order, expected rank or a known answer.',
    'Physical scheduling: input-nonzero triggers and measured rank gain per arithmetic work replace clock deadlines when enabled. Physical work competes with seam reduction, including retained fill. Persistent cycles and recorded task/stop decisions make structural replay reproducible. The initial work-budget-only V27B policy admitted expensive MRK/OPE work too early and regressed badly at nine loops. V27C/D added the competition with seam work. It remains a heuristic, not a proof of an optimal schedule.',
    'Exact compression: a reusable augmented feature-space row-span compressor validates each discarded equation with an explicit exact combination of source rows. Every RHS/origin column participates. It is connected to direct collinear discovery, the existing collinear suffix projection, and a conservative MRK polynomial-convolution adapter. MRK targets stayed full rank in these workloads, so that extension saved no projection rows. Independent physical checks still generate their complete original target lists. This removes redundant equations; it introduces no new physical constraints or bonus-relation assumptions.',
    'Dependency routes: delta-varint encoding and estimated coefficient updates saved per stored byte select retained batches. A rejected candidate does not evict useful routes. The retained accounting budget is 64 MiB, with an 8 MiB logical learning budget; allocator and vector-capacity overhead are included in measured RSS. Incomplete routes are valid hints and must finish with exact residual elimination. Turning this feature off restores the legacy prefix-admission cache, not an absence of structural replay.',
    f"Eight-loop cross-instance result: {candidate8['cpu_hours_mean']:.6f} CPU-hours versus V26 {baseline8['cpu_hours_mean']:.6f}, a {100*compare(candidate8,baseline8)['cpu_saving']:.2f}% reduction, with {candidate8['peak_rss_bytes']/1e9:.3f} versus {baseline8['peak_rss_bytes']/1e9:.3f} GB peak RSS. This useful lower-loop result does not satisfy the separate nine-loop adoption criterion.",
    'Correctness: all new complete native runs reproduce all 38,989 eight-loop or 105,898 nine-loop reference coefficients. Complete original seam sweeps and fresh collinear, MRK and OPE checks pass. Source/binary hashes, coefficients, plan checksums, timing receipts and reference bindings were rechecked when sealing this report. No cgroup pressure or OOM event occurred. The prior forced-QQ reference binding is preserved; no new rational reconstruction was timed.',
    'Generic controls: 96 field cases, 48 successful structural replays, two exceptional-prime guarded fallbacks and 72 affine-map cases passed with both new and legacy switches. Additional controls cover 144 augmented compression cases, 64 deterministic schedules, 756,718 ordered-queue operations, worker ownership/exception propagation, packed-budget eviction and malformed-plan checks, exact density aliasing, persistent schedule cycles and 5,077 deliberately exercised wide-accumulator overflow corrections against canonical arithmetic.',
    'Failed alternatives are retained: A used uniform density sampling that aliased polynomial layouts; B corrected density but scheduled costly physical work too early; C added competition with seam work; D added exact delayed modular accumulation. All completed variants were numerically correct. Corrections were built into new package directories; timed binaries and outputs were not overwritten.',
    'Scope and resources: warm precomputed basis/RHS inputs; 16 workers unless explicitly labeled eight; 4 GiB physics and 3 GiB OPE-head caches; heap backend; one prime process at a time; 20 GB decimal benchmark cgroup maximum, 18 GB high threshold and no swap. Full input loading, checks, checkpoints and learning/replay I/O are included. Cold input generation, compilation and a new QQ reconstruction are excluded. Hardware counters were unavailable, so phase scaling does not establish measured cache-miss or memory-bandwidth rates. No ten-loop job was modified or benchmarked; its 50 GB decimal limit remains in force.',
    'These measurements do not establish universal speedups or global optimality. The reusable implementations are retained for further work, but the current nine-loop profile stays V26. Further experiments should target the profiled incoming-row reduction and physical projection costs, with complete CPU/RSS comparisons before adoption.',
    'Experimental profile: ' + str(O / 'experimental-profile.json'),
    'All measurements: ' + str(O / 'RESULTS.json'),
    'Reference bindings: ' + str(O / 'reference-binding-audit.json'),
    'Frozen source/build receipt: ' + str(P / 'v27-build-receipt.json'),
]
(O / 'REPORT.txt').write_text('\n\n'.join(lines) + '\n')
d.update(status='complete', selected='V26', experimental_candidate='V27D',
         candidate_accepted=False, comparisons=profile, ablations=ablations,
         sampling_note='Learning/replay pairs are not statistical repetitions.')
write(O / 'RESULTS.json', d)

completion = dict(state='complete', selected='V26', candidate_accepted=False,
                  finished=time.time(), builds=builds, plans=plans,
                  total_fresh_native_runs=len(d['runs']),
                  total_benchmark_cpu_hours=sum(x['cpu_hours'] for x in d['runs']),
                  native_coefficients_and_receipts_rechecked=True,
                  source_and_binary_hashes_rechecked=True, input_hashes=bindings,
                  production_modified=False, production_memory_limit_bytes=50000000000)
for name in ['REPORT.txt', 'RESULTS.json', 'reference-binding-audit.json', 'experimental-profile.json']:
    bind(O / name)
bind(O / 'v26-to-v27d-review.patch')
for source in sorted((R / 'bench/hexagon_v27').glob('*')):
    if source.is_file():
        bind(source)
for name in ['hexagon_v27_controls.py', 'hexagon_v27_report.py', 'hexagon_v27_gated_comparisons.py']:
    bind(R / 'bench' / name)
bind(__file__)
write(O / 'completion-receipt.json', completion)
write(O / 'execution-gate.json', dict(state='complete_candidate_rejected_retained_V26',
      selected='V26', candidate='V27D', acceptance=profile['acceptance'],
      results=str(O / 'RESULTS.json'), report=str(O / 'REPORT.txt'),
      production_modified=False, production_memory_limit_bytes=50000000000))

with (R / 'docs/hexagon_optimization_lessons.lock').open('a') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    path = R / 'docs/hexagon_optimization_lessons.json'
    lessons = read(path)
    entry = dict(id='v27-generic-workload-compression-campaign', state='complete',
                 baseline='Fresh V26, per current explicit user request',
                 candidate_accepted=False, selected='V26',
                 finding='All four generic directions implemented and validated. Eight-loop gains did not transfer to a >15% nine-loop CPU saving. Preserve V26 as the accepted nine-loop default. Compare physical gain/work against seam work including fill; exact density checks avoid polynomial-layout sampling aliases; reduced thread count has opposite CPU effects in different phases; packed coverage alone is not a speedup guarantee.',
                 profile=profile, failed_alternatives=[
                     'Uniform sampled density aliases polynomial recurrence layouts.',
                     'V27B work-only physical scheduling over-admits expensive early MRK/OPE constraints.',
                     'MRK augmented feature-space rows were full rank in these examples.',
                     'V27D/8 workers lowered elimination CPU and RSS but increased collinear/OPE CPU; the full nine-loop gate still failed.'],
                 evidence=[str(O / x) for x in ['REPORT.txt', 'RESULTS.json', 'completion-receipt.json']],
                 production_modified=False, production_limit_bytes=50000000000)
    lessons['lessons'] = [x for x in lessons['lessons'] if x.get('id') != entry['id']] + [entry]
    lessons['universal_optimization_policy']['existing_acceptance'] = (
        'V26 retains its historical acceptance against original V17. The current V27 campaign explicitly uses fresh V26 as the baseline and requires >15% fewer CPU-hours with <=5% higher peak RSS. V27 did not pass; do not credit earlier V17 savings again.')
    lessons['updated'] = '2026-10-09'
    lessons['scope'] = ('Conditional ten-loop result remains complete and unchanged. V26 remains the accepted nine-loop profile. V27 implements four generic directions and passes correctness controls, but fails the new nine-loop CPU gate against V26. Ten-loop speed and broad universal optimality remain unmeasured.')
    tmp = path.with_suffix('.tmp')
    write(tmp, lessons)
    tmp.replace(path)
print(json.dumps(dict(state='complete', runs=len(d['runs']),
                      candidate_accepted=False, comparisons={k: profile[k] for k in
                          ['comparison_nine', 'comparison_nine_eight_workers', 'comparison_eight']}), indent=2))
