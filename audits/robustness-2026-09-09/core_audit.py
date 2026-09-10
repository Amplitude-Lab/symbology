"""Real C++ executable probes, bounded to temporary files and short timeouts."""
import json, os, random, resource, subprocess, tempfile, time, zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
work=Path(tempfile.mkdtemp(prefix='symbology-cpp-audit.'))
probe=Path(os.environ.get('NUMERICAL_PROBE','/tmp/symbology-numerical-probe'))
records=[]
def limits(): resource.setrlimit(resource.RLIMIT_CORE,(0,0))
def call(argv,timeout=15):
    start=time.monotonic()
    try:
        p=subprocess.run(list(map(str,argv)),cwd=work,capture_output=True,text=True,errors='replace',timeout=timeout,preexec_fn=limits)
        return dict(rc=p.returncode,seconds=round(time.monotonic()-start,4),stdout=p.stdout,stderr=p.stderr)
    except subprocess.TimeoutExpired as e:
        return dict(rc=124,seconds=timeout,stdout=str(e.stdout or '')[-2000:],stderr='timeout')
def record(name,passed,r):
    item=dict(name=name,pass_=bool(passed),**r);records.append(item)
    print(name,'PASS' if passed else 'FAIL',r.get('rc'),flush=True)
def dims(file):
    r=call([probe,'dump',file]);return r
def positive(name,args):
    r=call(args);record(name,r['rc']==0,r);return r
def negative(name,args):
    r=call(args);record(name,r['rc']==1,r);return r
call([probe,'gen',work/'fixtures'])
f=work/'fixtures';o=work/'out';o.mkdir()
ops=ROOT/'tensor_ops';bs=ROOT/'bootstrap'
for p in [ROOT/'data/FEC_1.wxf', ROOT/'data_pentagon/dlogmat_pentagon.wxf',ROOT/'data_4pformfactor/dlogmat_full_4pformfactor.wxf']:
    positive('read bundled '+p.name,[ops,'dims',p])
for n in ['0','2','17']:
    r=positive('matrix power '+n,[ops,'power',f/'swap.wxf',n,o/f'power{n}.wxf'])
    if r['rc']==0:
        got=dims(o/f'power{n}.wxf')
        expected=['ENTRY 0 0 = 1','ENTRY 1 1 = 1'] if int(n)%2==0 else ['ENTRY 0 1 = 1','ENTRY 1 0 = 1']
        record('matrix power exact values '+n,all(e in got['stdout'] for e in expected),got)
for n in ['-1','2junk','1.5','999999999999999999999999999']:
    negative('reject invalid power '+n,[ops,'power',f/'swap.wxf',n,o/'bad-power.wxf'])
negative('reject nonsquare power',[ops,'power',f/'rect.wxf','2',o/'bad-square.wxf'])
negative('reject missing input',[ops,'dims',f/'missing.wxf'])
negative('reject missing output directory',[ops,'power',f/'swap.wxf','2',o/'missing'/'out.wxf'])
negative('reject full output device',[ops,'power',f/'swap.wxf','2','/dev/full'])
positive('ternary identity',[ops,'ternary',f/'tensor.wxf','I','I',o/'ternary.wxf'])
negative('reject contraction dimension mismatch',[ops,'tdot',f/'swap.wxf',f/'rect.wxf','1','2',o/'bad-dot.wxf'])
negative('reject invalid dot axis',[ops,'tdot',f/'swap.wxf',f/'swap.wxf','0','1',o/'bad-axis.wxf'])
negative('reject unsupported rank-one squeeze',[ops,'squeeze',f/'vector.wxf',o/'rank-one.wxf'])
positive('zero from exact cancellation',[ROOT/'tensor_add',f/'swap.wxf',f/'swap.wxf','1','-1',o/'zero.wxf'])
negative('reject zero denominator',[ROOT/'tensor_add',f/'swap.wxf',f/'swap.wxf','1/0','1',o/'bad-fraction.wxf'])
negative('reject malformed rational',[ROOT/'tensor_add',f/'swap.wxf',f/'swap.wxf','abc','1',o/'bad-fraction.wxf'])
raw=(f/'swap.wxf').read_bytes()
rng=random.Random(11)
corrupt=[b'',b'8:',raw[:12],raw[:len(raw)//2],b'not wxf']+[rng.randbytes(i*7+1) for i in range(12)]
for i,b in enumerate(corrupt):
    p=f/f'corrupt-{i}.wxf';p.write_bytes(b)
    negative('reject malformed WXF '+str(i),[ops,'dims',p])
r=call([bs,'--solve-collinear','--pair-cond',f/'inconsistent.wxf','--export-conditions','--out-stem','bad','--output-dir',o])
record('inconsistent solve returns nonzero',r['rc']!=0,r)
cond=o/'collinear/cond_bad.wxf'
if cond.exists():
    r=call([bs,'--solve-collinear','--pair-cond',cond,'--out-stem','roundtrip','--output-dir',o])
    record('condition export preserves inconsistency','INCONSISTENT' in r['stdout'] and not (o/'collinear/sol_roundtrip.wxf').exists(),r)
baseline=json.loads((ROOT/'audits/baseline-2026-09-09/manifest-computerhs.json').read_text())['crc32']
repeated=[]
for i in range(10):
    output=work/f'fresh {i}'/'output'
    r=call([ROOT/'compute_rhs','--target','SEW_3p1','--letter-projection',output/'collinear/colprojdiv_w1.wxf','--data-dir',ROOT/'data','--output-dir',output],timeout=30)
    actual={rel:f'{zlib.crc32((output.parent/rel).read_bytes()):08x}' for rel in baseline if '/2loop/' in rel and (output.parent/rel).exists()}
    ok=r['rc']==0 and len(actual)==5 and all(baseline[k]==v for k,v in actual.items())
    repeated.append({'repeat':i,'pass':ok,'seconds':r['seconds'],'crc32':actual})
    if i==0:
        r3=call([ROOT/'compute_rhs','--target','SEW_5p1','--letter-projection',output/'collinear/colprojdiv_w1.wxf','--data-dir',ROOT/'data','--output-dir',output],timeout=30)
        record('resume from two loops to three loops',r3['rc']==0,r3)
record('10 fresh two-loop runs match all 5 baseline files',all(r['pass'] for r in repeated),dict(repeats=repeated))
OUT.joinpath('core-results.json').write_text(json.dumps({'work':str(work),'records':records},indent=2))
print('WORK',work)
print('SUMMARY',sum(r['pass_'] for r in records),'passed',sum(not r['pass_'] for r in records),'failed')
