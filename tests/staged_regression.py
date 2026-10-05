"""Public native product-invariance contract, compared with explicit equations."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(command, success=True):
    result = subprocess.run(list(map(str, command)), cwd=ROOT, capture_output=True,
                            text=True, timeout=120)
    assert (result.returncode == 0) == success, (command, result.stdout[-2000:], result.stderr[-2000:])
    return result.stdout + result.stderr


with tempfile.TemporaryDirectory(prefix='symbology-staged-') as directory:
    work = Path(directory)
    run([ROOT/'tests/staged_kernel_probe', work])
    base = [ROOT/'bootstrap', '--sew', '-c', work/'condition.wxf',
            '-f', work/'FEC_2.wxf', '-l', work/'LEC_2.wxf', '--threads', 2]
    generators = ['--left-action', work/'left_action.wxf',
                  '--right-action', work/'right_action.wxf']
    for strategy in ('auto', 'staged'):
        output = work/f'{strategy}.wxf'
        run([*base, '-o', output, '--sew-strategy', strategy, '--local-reduction', 'raw' if strategy == 'staged' else 'reduced', *generators])
        run([ROOT/'bench/nmhv_kernel_check', 'sew-compare', output, work/'reference.wxf'])
    rejected = work/'rejected.wxf'
    assert 'counts differ' in run([*base, '-o', rejected, '--left-action', work/'left_action.wxf'], False)
    assert 'require --sew-strategy' in run([*base, '-o', rejected, '--sew-strategy', 'original', *generators], False)
    run([*base, '-o', rejected, '--left-action', work/'condition.wxf',
         '--right-action', work/'right_action.wxf'], False)
    assert not rejected.exists()
print('PASS native staged/auto product invariance against independent C3 equations; invalid action combinations rejected')
