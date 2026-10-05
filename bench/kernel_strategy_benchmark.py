"""Sequential, rotated fresh-process benchmark of general exact kernels."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--probe',type=Path,default=Path('bench/kernel_strategy_benchmark'))
p.add_argument('--output',type=Path,required=True)
p.add_argument('--repeats',type=int,default=3)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
result=dict(repeats=a.repeats,threads=8,cases=[],notes=[
    'Synthetic general matrix families, not polygon fixtures; not a representative sample of all applications.',
    'Solve time excludes matrix generation and independent exact validation.',
    'Peak process RSS includes input generation, retained validation input and exact validation.',
    'Fresh sequential processes, rotated methods; no OS cache flushing.'])
for name in ['chain','graph','zero-propagation','rational','wide','dense-tail','blocks']:
    case=dict(name=name,runs=[]);result['cases'].append(case)
    for repeat in range(a.repeats):
        for method in (['baseline','structured'] if repeat%2==0 else ['structured','baseline']):
            command=[str(a.probe.resolve()),name,method,'8']
            run=subprocess.run(command,capture_output=True,text=True,timeout=180)
            (a.output/f'{name}-{method}-{repeat}.log').write_text(run.stdout+run.stderr)
            assert run.returncode==0,(command,run.returncode,run.stdout[-1000:],run.stderr[-1000:])
            line=next(x for x in run.stdout.splitlines() if x.startswith('RESULT '))
            fields=dict(x.split('=',1) for x in line.split()[1:]);fields['repeat']=repeat
            case['runs'].append(fields)
            print(line,flush=True)
    case['medians']={m:{k:statistics.median(float(x[k]) for x in case['runs'] if x['method']==m)
                       for k in ['seconds','rss_kib']} for m in ['baseline','structured']}
    (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
