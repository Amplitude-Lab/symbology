import subprocess,time,json,statistics,os,hashlib
from pathlib import Path
R=Path('/home/ana/Documents/symbology');O=R/'output_symrep_optimization_20261002';B=R/'output_symrep_validated_20261001';A=Path('/home/ana/Downloads/hexagon_irr_bootstrap_v2_20260928_inspected/hexagon_irr_bootstrap');records=json.loads((O/"extra-results.json").read_text()) if (O/"extra-results.json").exists() else []
def run(label,cmd,cwd=R):
 if any(r["label"]==label for r in records):return
 t=time.perf_counter()
 with (O/(label+'.log')).open('w') as log:
  subprocess.run(['/usr/bin/time','-f','{"peak_rss_kib":%M}','-o',str(O/(label+'.usage.json')),*map(str,cmd)],cwd=cwd,env=dict(os.environ,THREADS='2'),stdout=log,stderr=subprocess.STDOUT,check=True)
 r=dict(label=label,wall_s=time.perf_counter()-t,**json.loads((O/(label+'.usage.json')).read_text()));records.append(r);print(r,flush=True);(O/'extra-results.json').write_text(json.dumps(records,indent=2)+'\n')
for name,kind,w in [('hexagon','lec',6),('heptagon','fec',4),('heptagon','lec',3)]:
 prepared=B/f'{name}-{kind}-0-prepared';condition=A/'data/dlogmat_full_ca_pd3.wxf' if name=='hexagon' else R/'data/dlogmat_E6.wxf';seed=B/'hexagon-reference-lec/w1.wxf' if name=='hexagon' else R/f'data/{kind.upper()}_1.wxf'
 for i in range(3):
  label=f'extra-{name}-{kind}-{i}';native=O/(label+'-native');after=O/(label+'-after')
  run(label+'-native',[O/'native_reference','forward' if kind=='fec' else 'backward',condition,seed,w,native,2])
  run(label+'-after',[R/'symrep','extend','--input',prepared,'--max-weight',w,'--output',after,'--threads',2])
  if name=='hexagon':
   archive=O/(label+'-archive');ref=O/(label+'-export');run(label+'-archive',[A/'bin/irr_lec','MHV',w,archive],A)
   run(label+'-export',[R/'bench/hexagon_symrep_reference','LEC',archive,w,A/'data',ref])
   if i==0:run(label+'-archive-verify',[R/'symrep','verify','--input',prepared,'--chain',after,'--reference',ref,'--max-weight',w,'--output',O/(label+'-archive-verify'),'--threads',2])
  if i==0:run(label+'-native-verify',[R/'symrep','verify','--input',prepared,'--chain',after,'--reference',native,'--max-weight',w,'--output',O/(label+'-native-verify'),'--threads',2])
  if i:
   first=O/f'extra-{name}-{kind}-0-after'
   matches=[]
   for weight in range(1,w+1):
    for suffix in ('.wxf','_copies.tsv'):
     file=f'w{weight}{suffix}'
     assert (after/file).read_bytes()==(first/file).read_bytes(),f'nonidentical repeated result: {after/file}'
     matches.append(file)
   (O/(label+'-identical.json')).write_text(json.dumps({'exactly_identical_to':str(first),'files':matches},indent=2)+'\n')
   print(label+' exactly identical to verified repeat 0',flush=True)
print('ALL EXTRA REFERENCES PASS',flush=True)
