#!/usr/bin/env python3
"""Reproduce one complete native 17+3 prime from precomputed, hash-bound bases.

No checkpoints or reference coefficients enter the solve. Reference comparison
occurs only after the solver exits successfully. This excludes basis generation.
"""
import argparse,hashlib,json,subprocess,time,uuid,struct
from pathlib import Path
from hexagon_ten_loop_benchmark import sha
from hexagon_bench_inputs import roots as input_roots
R=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--source',choices=['auto','local','bench'],default='auto',help='local: workstation archives; bench: local benchmarks/hexagon-mhv-ten-loop/inputs parts (not pushed)')
 a0=p.parse_known_args()[0];roots=input_roots(a0.source)
 for name,path in [('fec',roots['fec']),('lec',roots['lec']),('data',roots['data']),('solver',roots['solver']),('reference',roots['reference']),('library',roots['library'])]:p.add_argument('--'+name,type=Path,default=path)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--prime-index',type=int,choices=range(1,9),default=2);p.add_argument('--hours',type=float,default=8);p.add_argument('--plan-only',action='store_true');a=p.parse_args();out=a.output.resolve();assert not out.exists()
 n=a.prime_index;primes=[1152921504606847009,1152921504606847067,1152921504606847081,1152921504606847123,1152921504606847127,1152921504606847189,1152921504606847201,1152921504606847229];prime=primes[n-1]
 # Verify the coordinate system against the archived result before solving.
 basis=a.reference/'basis_identity.tsv';inputs={str(basis):sha(basis),str(a.solver):sha(a.solver)}
 for line in basis.read_text().splitlines():
  path,digest=line.split('\t');parts=Path(path).parts;root=a.fec if parts[0]=='FEC' else a.lec;f=root.joinpath(*parts[1:]);f=Path(str(f)+'.shape') if not f.exists() and f.suffix=='.wxf' else f;assert sha(f)==digest,str(f);inputs[str(f)]=digest
 for f in a.data.rglob('*'):
  if f.is_file():inputs[str(f)]=sha(f)
 env={'THREADS':'16','PHYSICS_CACHE_MIB':'4096','OPE_HEAD_CACHE_MIB':'8192','MIMALLOC_PURGE_DELAY':'0','HEXAGON_MHV_DOUBLE_BONUS':'1','HEXAGON_MAPPED_FILE':str(out/'scratch.mmap'),'HEXAGON_MAPPED_GIB':'192','HEXAGON_MAPPED_SOFT_GIB':'40','LD_LIBRARY_PATH':str(a.library.resolve())}
 solver=[str(a.solver.resolve()),'17','3',str(a.fec.resolve()),str(a.lec.resolve()),str(a.data.resolve()),str(out/'solve'),'9000','cache','right','--prime-index',str(n)]
 cmd=['systemd-run','--user','--wait','--pipe','--unit=symbology-hexagon-prime-benchmark-'+uuid.uuid4().hex[:10],'--property=SuccessExitStatus=2','--property=MemoryHigh=48000000000','--property=MemoryMax=50000000000','--property=MemorySwapMax=0','--property=TasksMax=128','--working-directory='+str(R)]
 cmd+=['--setenv='+k+'='+v for k,v in env.items()]
 cmd+=['/usr/bin/python3',str(R/'bench/measure_process.py'),'--report',str(out/'resources.json'),'--log',str(out/'run.log'),'--status',str(out/'resource-status.json'),'--timeout',str(a.hours*3600),'--virtual-gib','256','--rss-gib','46.566128730773926','--min-available-gib','3','--',*solver]
 manifest={'schema':'hexagon-complete-native-prime-benchmark-v1','prime_index':n,'prime':prime,'status':'planned','scope':'Full native solve from precomputed bases; no precomputed solve checkpoints consumed. Fitting assumptions identical to V17.','input_hashes':inputs,'command':cmd,'environment':env,'threads':16,'hard_memory_bytes':50000000000}
 if a.plan_only:print(json.dumps(manifest,indent=2));return
 out.mkdir(parents=True);manifest['status']='running';(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 code=subprocess.run(cmd,cwd=R).returncode;r=json.loads((out/'resources.json').read_text());assert code==0 and r['exit_code']==0
 actual=out/f'solve/joint/cache/p{prime}.bin';ref=a.reference/f'cache/p{prime}.bin'
 for f in [actual,ref]:assert sha(f)==Path(str(f)+'.sha256').read_text().strip()
 x=actual.read_bytes();y=ref.read_bytes();assert struct.unpack('<8s7Q',x[:64])==struct.unpack('<8s7Q',y[:64]);assert x[72:]==y[72:]
 assert all(sha(f)==h for f,h in inputs.items())
 manifest.update(status='passed',measurement=r,all_275246_coefficients_equal=True,actual_sha256=sha(actual),reference_sha256=sha(ref),completed=time.time());(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('FULL_NATIVE_PRIME_BENCHMARK_PASS')
if __name__=='__main__':main()
