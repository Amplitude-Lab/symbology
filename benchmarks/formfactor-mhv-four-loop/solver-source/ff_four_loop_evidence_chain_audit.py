#!/usr/bin/env python3
"""Audit the completed four-loop certificate chain without resolving it."""
import hashlib
import json
import math
from pathlib import Path
import re
import time

base = Path(__file__).resolve().parents[1]
out = base / 'results/ff_four_loop'
started = time.monotonic()

def load(p):
    return json.loads((out / p).read_text())

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

sewing = load('affine_TR_sewing_certificate.json')
assert sewing['augmented_kernel_dimension'] == 4
assert sewing['modular_rank'] + 4 == sewing['presolved_columns']
assert sewing['certified_residue_product_bits'] > sewing['integer_residual_bound_bits']
assert sewing['finite_boundary_targets_read'] is False
for name, digest in sewing['hashes'].items():
    assert sha(out / name) == digest, name
extract = load('affine_extraction_independent_audit.json')
coordinate = load('final_coordinate_chain_independent_audit.json')
assert extract['all_checks_pass'] and coordinate['all_checks_pass']

# Independently replay exact Lucas certificates.  Factorization is not trusted:
# every product, recursive prime hypothesis, modular power and gcd is checked.
old = load('projection_cfirst/TC_w6_modcert.prime_audit.json')
certs = {c['prime']: c for c in old['certificates']}
proved = set()

def prove(n):
    if n in proved:
        return
    if n == 2:
        assert certs[n]['base_case']
    else:
        c = certs[n]
        factors = c['complete_n_minus_one_factors']
        assert math.prod(factors) == n - 1
        for q in set(factors):
            prove(q)
        a = c['witness']
        assert pow(a, n - 1, n) == 1
        assert all(math.gcd(pow(a, (n - 1) // q, n) - 1, n) == 1
                   for q in set(factors))
    proved.add(n)

boundaries = {}
for name, q, initial in [('DC', 22, 3), ('TC', 29, 0)]:
    folder = out / 'boundary' / name / 'affine_w8_mixed'
    receipt = json.loads((folder / 'receipt.json').read_text())
    blocks = [json.loads(x) for x in (folder / 'blocks.jsonl').read_text().splitlines()]
    assert receipt['consistent'] and receipt['initial_affine_dimension'] == initial
    assert receipt['final_augmented_dimension'] == 1
    assert receipt['all_suffix_pairs_checked'] and receipt['all_spurious_words_retained']
    assert receipt['finite_targets_read'] is False
    assert len(blocks) == receipt['checked_pair_blocks'] == receipt['nonzero_input_pair_blocks']
    assert len(blocks) + receipt['exact_zero_input_pair_blocks'] == q*q
    assert len({tuple(b['pair']) for b in blocks}) == len(blocks)
    assert all(b['exact_kernel_residual_zero'] for b in blocks)
    assert sum(b['components'] for b in blocks) == receipt['full_coordinate_components_checked']
    dimension = initial + 1
    for i, b in enumerate(blocks, 1):
        assert b['block'] == i and b['before'] == dimension
        assert 1 <= b['after'] <= dimension
        dimension = b['after']
    assert dimension == 1
    compare = load(f'holdouts/{name}_final_candidate_check/receipt.json')
    guard = load(f'holdouts/{name}_final_candidate_check.guard.json')
    assert guard['returncode'] == 0 and guard['stop_reason'] is None
    assert compare['exact_boundary_equality'] and compare['all_output_components_certified']
    assert compare['output_component_count'] == compare['boundary_prefix_dimension'] * q*q
    primes = compare['accepted_primes']
    assert len(primes) == len(set(primes)) == compare['certificate_primes']
    for p in primes:
        prove(p)
    assert math.prod(primes) == int(compare['certificate_modulus'])
    assert int(compare['certificate_modulus']) > int(compare['certificate_integer_bound'])
    log = (out / f'holdouts/{name}_final_candidate_check.log').read_text()
    passes = re.findall(r'CERT_PRIME prime=(\d+) all_components=(\d+)', log)
    assert [int(p) for p, _ in passes] == primes
    assert all(int(n) == compare['output_component_count'] for _, n in passes)
    for path, digest in guard['sha256'].items():
        assert sha(path) == digest, path
    boundaries[name] = {
        'initial_affine_dimension': initial,
        'final_affine_dimension': 0,
        'checked_nonzero_pair_blocks': len(blocks),
        'verified_zero_input_blocks': receipt['exact_zero_input_pair_blocks'],
        'regularity_components': receipt['full_coordinate_components_checked'],
        'holdout_components': compare['output_component_count'],
        'accepted_proven_primes': primes,
        'certificate_modulus': compare['certificate_modulus'],
        'certificate_integer_bound': compare['certificate_integer_bound'],
        'all_logged_prime_passes_complete': True,
        'all_recorded_hashes_match_current_inputs': True,
        'finite_target_not_used_to_solve': True,
    }

dc = out / 'boundary/DC/affine_w8_mixed/origin.wxf'
tc = out / 'boundary/TC/affine_w8_mixed/origin.wxf'
assert sha(dc) == sha(tc)
result = {
    'all_checks_pass': True,
    'sewing': {'complete_augmented_dimension': 4,
               'sampled_intermediate_nullity_not_used_as_final': 256,
               'rank_plus_nullity_equals_presolved_columns': True,
               'full_operator_residue_product_bits': 121,
               'integer_residual_bound_bits': 91,
               'full_operator_not_rerun_by_this_audit': True},
    'independent_coordinate_extraction': extract['all_checks_pass'],
    'independent_final_coordinate_transport': coordinate['all_checks_pass'],
    'boundaries': boundaries,
    'final_DC_TC_origin_sha256': sha(tc),
    'prime_proof': 'Recursive exact Lucas order certificates independently replayed with Python integers.',
    'replayed_prime_certificates': [certs[n] for n in sorted(proved)],
    'scope': 'Unique four-loop candidate in the proposed adjacency and fixed TR12 class, using regularity only. Full normalized DC and TC boundaries equal independently published targets. This does not certify larger-ES uniqueness, local adjacency/ES equivalence, or equality to an independently known full four-point symbol.',
    'review_limit': 'Boundary block receipts, prior independent algorithm audits, exact saved coordinate-map replay, complete comparator prime logs, heights and hashes were checked. The full native sewing and boundary equation operators were not independently recomputed in this light audit.',
    'source_sha256': sha(__file__),
    'seconds': time.monotonic() - started,
}
(out / 'four_loop_evidence_chain_independent_audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('replayed_prime_certificates',)}, indent=2))
