import subprocess,time,json,hashlib
from pathlib import Path
R=Path('/home/ana/Documents/symbology');O=R/'output_symrep_optimization_20261002';B=R/'output_symrep_validated_20261001'
report={'solver_sha256':hashlib.sha256((R/'symrep').read_bytes()).hexdigest(),'runs':[]}
for name,kind,w in [('hexagon','lec',6),('heptagon','fec',4),('heptagon','lec',3)]:
 for i in range(3):
  label=f'final-{name}-{kind}-{i}';dest=O/label;start=time.perf_counter()
  with (O/(label+'.log')).open('w') as log:
   subprocess.run(['/usr/bin/time','-f','{"peak_rss_kib":%M}','-o',str(O/(label+'.usage.json')),str(R/'symrep'),'extend','--input',str(B/f'{name}-{kind}-0-prepared'),'--max-weight',str(w),'--output',str(dest),'--threads','2'],stdout=log,stderr=subprocess.STDOUT,check=True)
  record=dict(label=label,wall_s=time.perf_counter()-start,**json.loads((O/(label+'.usage.json')).read_text()))
  reference=O/f'extra-{name}-{kind}-0-after';digests={}
  for weight in range(1,w+1):
   for suffix in ('.wxf','_copies.tsv'):
    file=f'w{weight}{suffix}';actual=(dest/file).read_bytes();assert actual==(reference/file).read_bytes(),f'different certified recurrence: {file}';digests[file]=hashlib.sha256(actual).hexdigest()
  record['exactly_identical_to_verified_chain']=str(reference);record['sha256']=digests
  report['runs'].append(record);(O/'final-followup-results.json').write_text(json.dumps(report,indent=2)+'\n');print(record['label'],record['wall_s'],record['peak_rss_kib'],'EXACT PASS',flush=True)
