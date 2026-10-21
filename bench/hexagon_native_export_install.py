#!/usr/bin/env python3
"""Explicitly export reconstructed native QQ data without changing producers."""
import json,os,shutil,subprocess,time
from pathlib import Path
from hexagon_candidate_export_install import sha,write
R=Path(__file__).resolve().parents[1];B=R/'output_hexagon_17p3_20261007';V=R/'output_hexagon_validation_20261009'
source=B/'ten-v17-main/joint';target=V/'native-qq-17p3';receipt=V/'native-reconstruction-7p96.json'
r=json.loads(receipt.read_text());assert r['passed'] and r['holdout_mismatches']==[0]
vector=Path(r['vector']);assert sha(vector)==r['vector_sha256']
export=R/'output_hexagon_candidate_export_20261009/candidate_qq_export'
vpath=export.parent/'validation-passed.json';v=json.loads(vpath.read_text());assert v['passed'] and sha(export)==v['binary_sha256']
env=os.environ.copy();env.update(THREADS='16',LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib'))
inputs={str(p):sha(p) for p in [receipt,vector,vpath,export,Path(__file__)]}
subprocess.run([str(export),'17','3',str(R/'output_hexagon_16p4_20261006/FEC'),str(R/'output_hexagon_16p4_20261006/LEC/MHV'),str(vector),str(target)],env=env,check=True)
for n in ['input_identity.txt','basis_identity.tsv','stage.tsv']:
 inputs[str(source/n)]=sha(source/n);shutil.copy2(source/n,target/n)
(target/'cache').mkdir()
roles=[]
for c in r['cache_receipts']:
 p=Path(c['path']);assert sha(p)==c['sha256'];roles.append(f"p{c['prime']}\t{c['role']}\n")
 for path in [p,Path(str(p)+'.sha256')]:inputs[str(path)]=sha(path);shutil.copy2(path,target/'cache'/path.name)
write(target/'reconstruction_primes.tsv',''.join(roles))
write(target/'summary.tsv','crt_primes\t7\nretrospective_holdout_primes\t1\nmodulus_bits\t421\ndirections\t0\n')
manifest={'schema':'hexagon-native-QQ-explicit-export-v1','created':time.time(),'source_native_root':str(source),'reconstruction_receipt':str(receipt),'input_hashes':inputs,'producer_unchanged':True,'primes_1_to_7_role':'CRT','prime_8_role':'retrospective check; was not fresh after candidate selection','physics_certificate':False,'required_acceptance':'Exact rational conversion equality with independently certified 16+4; new physical checks recorded separately.','assumption':'Extrapolated MHV double-final-entry bonus relations; conditional uniqueness.'}
write(target/'manifest.json',json.dumps(manifest,indent=2)+'\n')
subprocess.run([str(export),'--seal',str(target)],env=env,check=True)
assert all(sha(p)==h for p,h in inputs.items())
write(V/'native-export-install.json',json.dumps({'passed':True,'root':str(target),'inputs':inputs,'qq_complete_sha256':sha(target/'qq_complete.tsv'),'manifest_sha256':sha(target/'manifest.json')},indent=2)+'\n')
