#!/usr/bin/env python3
"""Install an explicitly provenanced, separate QQ candidate; never alter production."""
import hashlib, importlib.util, json, os, shutil, struct, subprocess, time
from pathlib import Path

R=Path(__file__).resolve().parents[1]
B=R/'output_hexagon_17p3_20261007'
E=R/'output_hexagon_candidate_export_20261009'
TARGET=B/'ten-asymmetric-qq-16p4'
SHIFT=B/'ten-v17-shifted-16p4'
NATIVE=B/'ten-v17-main/joint'
ENV=os.environ.copy()
ENV.update(THREADS='2',LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib'))

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def write(p,data):
    with p.open('x') as f:
        f.write(data);f.flush();os.fsync(f.fileno())

def main():
    assert not TARGET.exists(), 'Refuse to replace any candidate'
    validation=E/'validation-passed.json';v=json.loads(validation.read_text())
    assert v['passed'] and sha(E/'candidate_qq_export')==v['binary_sha256']
    assert sha(R/'bench/hexagon_candidate_qq_export.cpp')==v['source_sha256']
    for c in v['known_result_controls']:
        for f,h in c['files'].items():
            assert sha(Path(c['reference'])/f)==h==sha(Path(c['export'])/f)
    frozen=B/'audit/crt-six-prime-speculative-candidate-20261008T1532Z.json'
    vector=frozen.with_suffix('.candidate.tsv')
    assert sha(frozen)=='5f8a78259503f24d08f614d1aea65a1ca2083a09ecea0869b62ea772b6c4b2ca'
    assert sha(vector)=='97355630801f12642cc284e46862af5239fe2dd5d98f3ca330752e511fca07d8'
    fresh=B/'audit/crt-six-prime-candidate-fresh-primes7-8-20261009T0534Z.json'
    f=json.loads(fresh.read_text());assert f['both_fresh_checks_pass'] and f['candidate_sha256']==sha(vector)
    assert [c['prime_index'] for c in f['checks']]==[7,8]
    spec=importlib.util.spec_from_file_location('migration',R/'bench/hexagon_v17_cache_migration.py')
    migration=importlib.util.module_from_spec(spec);spec.loader.exec_module(migration)
    migration_receipt=B/'ten-v17-runner/v16-prime-import.json'
    migration.verify_receipt(migration_receipt)
    ledger=json.loads((B/'ten-v17-runner/ledger.json').read_text())
    inputs={str(p):sha(p) for p in [frozen,vector,fresh,validation,migration_receipt,R/'bench/hexagon_v17_cache_migration.py']}
    for path in (R/'output_hexagon_16p4_20261006/joint-v8-package/src').rglob('*'):
        if path.is_file():inputs[str(path)]=sha(path)
    caches=[];roles={};modulus=1
    for i in range(1,9):
        entry=ledger['primes'][str(i)]
        for key,h,k,n in [('native',17,3,275246),('conversion',16,4,518048)]:
            a=entry[key];path=Path(a['path']);assert sha(path)==a['sha256']
            identity=path.parent.parent/'input_identity.txt'
            data=migration.cache(path,identity,h,k,a['prime']);assert len(data)==72+8*n
            inputs[str(path)]=sha(path);inputs[str(Path(str(path)+'.sha256'))]=sha(str(path)+'.sha256')
        prime=entry['conversion']['prime'];roles['p'+str(prime)]='CRT' if i<=6 else 'validation'
        if i<=6:modulus*=prime
        source_receipt=SHIFT/('source_p'+str(prime)+'.sha256')
        assert source_receipt.read_text()==entry['native']['sha256']+'\n'
        inputs[str(source_receipt)]=sha(source_receipt)
        if i>=7:
            check=f['checks'][i-7];assert check['prime']==prime and check['sha256']==entry['conversion']['sha256'] and check['mismatches']==check['noninvertible']==0
        caches.append({'prime_index':i,**entry})
    for name in ['input_identity.txt','basis_identity.tsv','stage.tsv']:
        inputs[str(SHIFT/name)]=sha(SHIFT/name)
    command=[str(E/'candidate_qq_export'),'16','4',str(R/'output_hexagon_16p4_20261006/FEC'),str(R/'output_hexagon_16p4_20261006/LEC/MHV'),str(vector),str(TARGET)]
    with (E/'ten-export.log').open('x') as log:subprocess.run(command,env=ENV,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    # Keep the exact shift-map identity and original native parent. No rebinding
    # or claim that the native/shifted producer generated this rational vector.
    for name in ['input_identity.txt','basis_identity.tsv','stage.tsv']:shutil.copy2(SHIFT/name,TARGET/name)
    (TARGET/'cache').mkdir()
    for c in caches:
        p=Path(c['conversion']['path'])
        for src in [p,Path(str(p)+'.sha256')]:shutil.copy2(src,TARGET/'cache'/src.name)
        name='source_p'+str(c['conversion']['prime'])+'.sha256';shutil.copy2(SHIFT/name,TARGET/name)
    write(TARGET/'reconstruction_primes.tsv',''.join(k+'\t'+v+'\n' for k,v in sorted(roles.items())))
    write(TARGET/'summary.tsv',f'crt_primes\t6\nvalidation_primes\t2\nmodulus_bits\t{modulus.bit_length()}\ndirections\t0\n')
    manifest={'schema':'hexagon-external-asymmetric-QQ-export-v1','created_unix_seconds':time.time(),'source_shifted_root':str(SHIFT),'source_native_root':str(NATIVE),'candidate_receipt':str(frozen),'candidate_sha256':sha(vector),'fresh_check_receipt':str(fresh),'export_validation':str(validation),'exporter_binary_sha256':sha(E/'candidate_qq_export'),'exporter_source_sha256':v['source_sha256'],'installer_sha256':sha(__file__),'prime_provenance':caches,'input_hashes':inputs,'production_modified':False,'qq_export_not_physics_certificate':True,'unrestricted_uniqueness_claimed':False,'assumption':'Extrapolated double-final-entry bonus relations remain explicit and conditional.'}
    write(TARGET/'manifest.json',json.dumps(manifest,indent=2)+'\n')
    subprocess.run([str(E/'candidate_qq_export'),'--seal',str(TARGET)],env=ENV,check=True,timeout=120)
    for p,h in inputs.items():assert sha(p)==h, 'Input changed during export: '+p
    for p in TARGET.rglob('*'):
        if p.is_file():
            with p.open('rb') as fd:os.fsync(fd.fileno())
    for p in [TARGET/'cache',TARGET,TARGET.parent]:
        fd=os.open(p,os.O_RDONLY|os.O_DIRECTORY);os.fsync(fd);os.close(fd)
    receipt={'schema':'hexagon-external-QQ-install-receipt-v1','created_unix_seconds':time.time(),'qq_root':str(TARGET),'manifest_sha256':sha(TARGET/'manifest.json'),'qq_receipt_sha256':sha(TARGET/'qq_complete.tsv'),'installer_sha256':sha(__file__),'candidate_unchanged':True,'six_CRT_plus_two_fresh_primes':True,'production_modified':False,'physics_certified':False}
    write(E/'ten-install-receipt.json',json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
