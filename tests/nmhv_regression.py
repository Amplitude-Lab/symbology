"""Public heptagon NMHV acceptance: exact spaces and published symmetry counts."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
check = ROOT/'bench/nmhv_kernel_check'
bootstrap = ROOT/'bootstrap'

def run(*args):
    result = subprocess.run(list(map(str, args)), cwd=ROOT, capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, (args, result.stdout[-2000:], result.stderr[-2000:])
    return result.stdout

with tempfile.TemporaryDirectory(prefix='symbology-nmhv-') as tmp:
    work = Path(tmp); fixture = work/'fixture'
    run(check, 'fixture', fixture)
    chains = {}
    for mode in ('original', 'default', 'staged'):
        chain = work/mode; chain.mkdir(); chains[mode] = chain
        shutil.copy2(fixture/'LEC_1.wxf', chain/'w1.wxf')
        for w in (2, 3):
            options = [] if mode == 'default' else ['--kernel-strategy', mode]
            run(bootstrap, '--extend', '-c', ROOT/'data/dlogmat_E6.wxf', '-l',
                chain/f'w{w-1}.wxf', '-o', chain/f'w{w}.wxf', '--threads', 2, *options)
    for mode in ('default','staged'):
        print(run(check, 'compare', 'backward', chains[mode], chains['original'], 3).strip())
    fec = work/'fec'; fec.mkdir(); shutil.copy2(ROOT/'data/FEC_1.wxf', fec/'w1.wxf')
    for w in (2, 3):
        run(bootstrap, '--extend', '-c', ROOT/'data/dlogmat_E6.wxf', '-f', fec/f'w{w-1}.wxf',
            '-o', fec/f'w{w}.wxf', '--threads', 2)
    for weight in (1, 3):
        for strategy in ('auto', 'original', 'streamed', 'staged'):
            options = [] if strategy == 'auto' else ['--sew-strategy', strategy]
            output = work/f'sew-{weight}-{strategy}.wxf'
            log = run(bootstrap, '--sew', '-c', ROOT/'data/dlogmat_E6.wxf', '-f', fec/f'w{weight}.wxf',
                      '-l', fixture/'LEC_1.wxf', '-o', output, '--threads', 2, *options)
            if strategy == 'auto':
                assert 'sew_strategy=original' in log, log  # smaller left coefficient space
        for strategy in ('auto', 'streamed', 'staged'):
            print(run(check, 'sew-compare', work/f'sew-{weight}-{strategy}.wxf',
                      work/f'sew-{weight}-original.wxf').strip())
        print(run(check, 'hept-symmetry', fixture, fec, weight, work/f'sew-{weight}-auto.wxf').strip())
print('PASS NMHV default CLI, exact original/streamed/staged spaces, and independent published dimensions 5/11')
