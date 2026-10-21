#!/usr/bin/env python3
"""Verify an explicit scheduling-only import against native and shifted controls."""
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import hexagon_v17_cache_migration as m

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output_hexagon_17p3_20261007'
SHARED = ROOT / 'output_hexagon_16p4_20261006'
OUT = BASE / 'v17-migration-controls'
PRIME = 1152921504606847009
DRIVER = BASE / 'joint-v17-schedule-package/bin/joint_physics_v17'
CONVERTER = SHARED / 'joint-v8-package/bin/irr_convert'


def main():
    m.require(not OUT.exists(), 'Preserve the existing validation; do not overwrite it')
    OUT.mkdir()
    result = {'schema': 'hexagon-v17-migration-controls-v1', 'started': time.time(), 'status': 'running',
              'migration_module_sha256': m.sha(Path(m.__file__)), 'driver_sha256': m.sha(DRIVER), 'steps': []}
    def save():
        (OUT / 'progress.json').write_text(json.dumps(result, indent=2) + '\n')
    def run(label, arguments):
        command = [sys.executable, str(ROOT / 'bench/measure_process.py'), '--report', str(OUT / (label + '.json')),
            '--status', str(OUT / (label + '.status.json')), '--log', str(OUT / (label + '.log')),
            '--timeout', '600', '--rss-gib', '5', '--min-available-gib', '3', '--virtual-gib', '96', '--',
            'env', 'THREADS=4', *map(str, arguments)]
        subprocess.run(command, cwd=ROOT, stdout=subprocess.DEVNULL, check=True)
        report = json.loads((OUT / (label + '.json')).read_text())
        m.require(report['exit_code'] == 0 and not any(report.get(k) for k in ['timeout', 'memory_pressure', 'rss_limit_exceeded', 'cancel_signal']), label + ' failed')
        result['steps'].append({'label': label, 'report': report});save()
    save()
    try:
        source = BASE / 'v16-controls/v16-eight-t16/joint'
        native = OUT / 'imported-eight/joint'
        expected = BASE / 'v17-schedule-controls/eight-t4/joint'
        imported = m.import_native(source, native, DRIVER, PRIME)
        for name in ['input_identity.txt', 'stage.tsv', 'basis_identity.tsv', f'cache/p{PRIME}.bin']:
            m.require((native / name).read_bytes() == (expected / name).read_bytes(), 'Independent native control mismatch: ' + name)
        result['native_matches_independently_computed_v17_byte_for_byte'] = True
        run('native-bind-resume', ['HEXAGON_MHV_DOUBLE_BONUS=1', 'HEXAGON_MAPPED_GIB=32', 'HEXAGON_MAPPED_SOFT_GIB=1', 'HEXAGON_MAPPED_FILE=' + str(OUT / 'unused.mmap'),
            DRIVER, 13, 3, SHARED / 'FEC', SHARED / 'LEC/MHV', SHARED / 'data', native.parent,
            2048, 'cache', 'right', '--prime-index', 1])
        log = (OUT / 'native-bind-resume.log').read_text()
        m.require('JOINT_PRIME_COMPLETE' not in log and m.sha(native / 'cache' / f'p{PRIME}.bin') == imported['destination_sha256'],
                  'Imported prime unexpectedly recomputed')
        for label, parent in [('old-shift', source), ('new-shift', native)]:
            run(label, [CONVERTER, 'shift-cache', 'MHV', 13, 3, parent, OUT / label, '--head', 12,
                '--fec-root', SHARED / 'FEC', '--lec-root', SHARED / 'LEC/MHV', '--kind', 'unique'])
        migration = m.import_shifted(OUT / 'old-shift', OUT / 'imported-shift', source, native, PRIME)
        for name in ['input_identity.txt', 'stage.tsv', 'basis_identity.tsv', f'source_p{PRIME}.sha256', f'cache/p{PRIME}.bin']:
            m.require((OUT / 'imported-shift' / name).read_bytes() == (OUT / 'new-shift' / name).read_bytes(), 'Independent shift mismatch: ' + name)
        result['shifted_import_matches_independent_conversion_byte_for_byte'] = True
        # Fault injection uses disposable copies, never the original controls.
        rejected = []
        with tempfile.TemporaryDirectory(prefix='hexagon-migration-', dir=OUT) as temporary:
            tmp = Path(temporary)
            payload = (source / 'cache' / f'p{PRIME}.bin').read_bytes()
            identity = tmp / 'input_identity.txt';identity.write_bytes((source / 'input_identity.txt').read_bytes())
            for label, data, keep_hash in [
                ('checksum', payload[:-1] + bytes([payload[-1] ^ 1]), True),
                ('identity', payload[:64] + struct.pack('<Q', 1) + payload[72:], False),
                ('nonunique', payload[:40] + struct.pack('<Q', 2) + payload[48:], False),
                ('dimension', payload[:48] + struct.pack('<Q', 7) + payload[56:], False),
                ('truncated', payload[:-8], False),
                ('noncanonical', payload[:72] + struct.pack('<Q', PRIME) + payload[80:], False),
            ]:
                path = tmp / (label + '.bin');path.write_bytes(data)
                Path(str(path) + '.sha256').write_text((m.sha(source / 'cache' / f'p{PRIME}.bin') if keep_hash else m.sha(path)) + '\n')
                try:
                    m.cache(path, identity, 13, 3, PRIME)
                except RuntimeError:
                    rejected.append(label)
                else:
                    raise RuntimeError('Fault was not rejected: ' + label)
        result.update(status='migration_controls_passed', finished=time.time(), rejected_corruptions=rejected,
                      native_import=imported, shifted_import=migration)
    except Exception as error:
        result.update(status='failed', error=str(error), finished=time.time());save();raise
    save()


if __name__ == '__main__':
    main()
