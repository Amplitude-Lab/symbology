#!/usr/bin/env python3
"""Measured replay benchmark using hash-bound saved ten-loop inputs.

This is a warm-basis phase benchmark, not cold basis generation or a new prime.
The same frozen numerical function is used for both memory backends.
"""
import argparse,hashlib,json,os,subprocess,time,uuid
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['replay']);p.add_argument('--output',type=Path,required=True);p.add_argument('--storage',choices=['heap','mapped'],required=True);p.add_argument('--threads',type=int,default=16);p.add_argument('--timeout',type=int,default=3600);a=p.parse_args()
 assert 1<=a.threads<=16
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 P=R/'output_hexagon_16p4_20261006';B=R/'output_hexagon_17p3_20261007';V=R/'output_hexagon_validation_20261009';pkg=B/'joint-v17-schedule-package';prime=1152921504606847067
 cp=B/'ten-v17-main/early_ope_small_s';ref=B/'ten-v17-main/joint'
 files=[V/'hexagon_seam_replay',R/'bench/hexagon_seam_replay.cpp',R/'bench/hexagon_ten_loop_benchmark.py',R/'bench/measure_process.py',V/'hexagon_seam_replay.build.json',pkg/'tests/joint_physics_v17.cpp',cp/'basis_identity.tsv',cp/'input_identity.txt',cp/f'cache/p{prime}.csr',cp/f'cache/p{prime}.csr.sha256',ref/'basis_identity.tsv',ref/'input_identity.txt',ref/f'cache/p{prime}.bin',ref/f'cache/p{prime}.bin.sha256',P/'data/dlogmat_full_ca_pd3.wxf',P/'data/dlogmat_full_ca_pd3_blocks.tsv']
 files += [f for f in (pkg/'src').rglob('*') if f.is_file()]
 files += [R/'bench/hexagon_runtime_recovery_validate.py',R/'output_hexagon_runtime_recovery_20261008/validation-passed.json']
 subprocess.run(['python3',str(R/'bench/hexagon_runtime_recovery_validate.py'),'--verify-ready'],check=True,cwd=R,stdout=subprocess.DEVNULL)
 binding={str(f):sha(f) for f in files};assert sha(cp/f'cache/p{prime}.csr')==(cp/f'cache/p{prime}.csr.sha256').read_text().strip();assert sha(ref/f'cache/p{prime}.bin')==(ref/f'cache/p{prime}.bin.sha256').read_text().strip()
 command=[str(V/'hexagon_seam_replay'),'17','3',str(P/'FEC'),str(P/'LEC/MHV'),str(P/'data'),str(cp),str(ref),str(prime)]
 env={'HEXAGON_MAPPED_FILE':'','THREADS':str(a.threads),'LD_LIBRARY_PATH':str(R/'output_hexagon_runtime_recovery_20261008/env/lib'),'MIMALLOC_PURGE_DELAY':'0'}
 if a.storage=='mapped':env.update(HEXAGON_MAPPED_FILE=str(out/'scratch.mmap'),HEXAGON_MAPPED_GIB='192',HEXAGON_MAPPED_SOFT_GIB='40')
 unit='symbology-hexagon-replay-'+a.storage+'-'+uuid.uuid4().hex[:10]
 manager=['systemd-run','--user','--wait','--pipe','--unit='+unit,'--property=MemoryHigh=48000000000','--property=MemoryMax=50000000000','--property=MemorySwapMax=0','--property=TasksMax=128','--working-directory='+str(R)]
 for k,v in env.items():manager.append('--setenv='+k+'='+v)
 manager+=['/usr/bin/python3',str(R/'bench/measure_process.py'),'--report',str(out/'resources.json'),'--log',str(out/'run.log'),'--status',str(out/'resource-status.json'),'--timeout',str(a.timeout),'--virtual-gib','256','--rss-gib','46.566128730773926','--min-available-gib','3','--',*command]
 manifest={'schema':'hexagon-ten-seam-replay-benchmark-v1','created':time.time(),'scope':'Post-OPE final-seam replay with precomputed bases and checkpoints. Not cold generation, not an additional CRT prime.','storage':a.storage,'threads':a.threads,'prime':prime,'hard_memory_bytes':50000000000,'input_hashes':binding,'command':manager,'environment':env,'host':subprocess.check_output(['lscpu'],text=True),'status':'running'}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 code=subprocess.run(manager,cwd=R).returncode
 stat=subprocess.run(['systemctl','--user','show',unit,'-p','MemoryPeak','-p','CPUUsageNSec','-p','Result','-p','MemoryHigh','-p','MemoryMax','-p','MemorySwapMax'],capture_output=True,text=True)
 (out/'systemd-status.txt').write_text(stat.stdout+stat.stderr)
 report=json.loads((out/'resources.json').read_text());lines=(out/'run.log').read_text().splitlines();matches=[json.loads(l.split(' ',1)[1]) for l in lines if l.startswith('REPLAY_PASS ')]
 mapped_enabled=any(l.startswith('MAPPED_BUFFERS capacity_bytes=') for l in lines)
 passed=code==0 and report['exit_code']==0 and len(matches)==1 and matches[0]['all_coefficients_equal'] and mapped_enabled==(a.storage=='mapped')
 manifest.update(status='passed' if passed else 'failed',finished=time.time(),measurement=report,phase_result=matches,artifacts={str(f):sha(f) for f in [out/'resources.json',out/'run.log']});assert all(sha(f)==h for f,h in binding.items())
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');assert passed
 print(json.dumps({'passed':True,'storage':a.storage,'measurement':report,'phase_result':matches},indent=2))
if __name__=='__main__':main()
