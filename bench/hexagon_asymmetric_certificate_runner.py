#!/usr/bin/env python3
"""Bounded independent certification of the frozen six-plus-two QQ export.

Production is frozen only while this process runs and is thawed in finally and
by the service's ExecStopPost. Its sources, ledger and coefficient files are
never written here. A success is conditional on the recorded bonus assumptions.
"""
import hashlib, importlib.util, json, os, re, subprocess, time
from pathlib import Path

R=Path(__file__).resolve().parents[1];B=R/'output_hexagon_17p3_20261007'
C=B/'ten-asymmetric-runner';Q=B/'ten-asymmetric-qq-16p4'
E=R/'output_hexagon_candidate_export_20261009';S=R/'output_hexagon_16p4_20261006'
UNIT='symbology-hexagon-ten-v17-20261007.service'
CERT=R/'output_hexagon_certificate_provenance_20261007/bin/hexagon_shifted_joint_certificate'
ORIGIN=S/'joint-v8-package/bin/origin_check'
CASES=R/'audits/hexagon-ten-loop-2026-10-06/origin-cases.tsv'
MIG=B/'ten-v17-runner/v16-prime-import.json'
LIMIT=8_000_000_000

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for a in iter(lambda:f.read(1<<20),b''):h.update(a)
    return h.hexdigest()

def put(name,obj):
    p=C/name;t=p.with_suffix(p.suffix+'.tmp')
    with t.open('w') as f:json.dump(obj,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(t,p);fd=os.open(C,os.O_RDONLY|os.O_DIRECTORY);os.fsync(fd);os.close(fd)

def check():
    assert sha(CERT)=='d9ad435ac40f527ad65bf31c15f2058f498215645d232f6b17d79727e204b34c'
    assert sha(ORIGIN)=='8a32180f0d60d83fe3fa59ca02c77edc5a5ab596e9e9b73c6514d43f54c6b6ac'
    assert sha(CASES)=='7f193ed3d554778af10671db02c4d0b42cf2fe60a2e2e7d310cf2c860f158791'
    install=json.loads((E/'ten-install-receipt.json').read_text())
    assert sha(Q/'manifest.json')==install['manifest_sha256']
    assert sha(Q/'qq_complete.tsv')==install['qq_receipt_sha256']
    for name,h in (line.split('\t') for line in (Q/'qq_complete.tsv').read_text().splitlines()):
        if name not in ['schema','identity']:assert sha(Q/name)==h
    for name,h in (line.split('\t') for line in (Q/'cache_files.tsv').read_text().splitlines()):assert sha(Q/'cache'/name)==h
    manifest=json.loads((Q/'manifest.json').read_text())
    for p,h in manifest['input_hashes'].items():assert sha(p)==h,p
    assert 'INDEPENDENT_WOLFRAM_VERIFY_PASS columns=518048 nnz=491197' in (E/'ten-independent-verify.log').read_text()
    spec=importlib.util.spec_from_file_location('migration',R/'bench/hexagon_v17_cache_migration.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.verify_receipt(MIG)
    subprocess.run(['/usr/bin/python3',str(R/'bench/hexagon_runtime_recovery_validate.py'),'--verify-ready'],check=True,timeout=120)
    return manifest

def measured(name,args):
    report=C/(name+'.json');log=C/(name+'.log');status=C/(name+'.status.json')
    command=['/usr/bin/python3',str(R/'bench/measure_process.py'),'--report',str(report),'--log',str(log),'--status',str(status),'--timeout','21600','--virtual-gib','64','--rss-gib',str(7_600_000_000/1024**3),'--min-available-gib','2','--',*map(str,args)]
    put('progress.json',{'state':'running_'+name,'updated':time.time(),'command':command,'resource_status':str(status),'solver_log':str(log),'production_frozen':True,'final_symbol_certified':False})
    env=os.environ.copy();env.update(THREADS='16',LD_LIBRARY_PATH=str(R/'output_hexagon_runtime_recovery_20261008/env/lib'),MIMALLOC_PURGE_DELAY='0')
    result=subprocess.run(command,env=env,timeout=21700)
    assert result.returncode==0,f'{name} manager exit {result.returncode}'
    data=json.loads(report.read_text());assert data['exit_code']==0 and not data['timeout'] and not data['memory_pressure'] and not data['rss_limit_exceeded']
    return {'report':str(report),'report_sha256':sha(report),'log':str(log),'log_sha256':sha(log),'measurement':data}

def main():
    C.mkdir(exist_ok=True);assert not (C/'progress.json').exists(),'No unchanged retries/overwrites'
    manifest=check();frozen=False
    try:
        pid=subprocess.check_output(['systemctl','--user','show',UNIT,'--property=MainPID','--value'],text=True).strip()
        assert pid=='15259','Production identity changed; reassess before freezing'
        assert 'hexagon_ten_loop_runner_v17.py' in Path('/proc',pid,'cmdline').read_bytes().decode()
        subprocess.run(['systemctl','--user','freeze',UNIT],check=True);frozen=True
        cg=Path('/sys/fs/cgroup')/subprocess.check_output(['systemctl','--user','show',UNIT,'--property=ControlGroup','--value'],text=True).strip().lstrip('/')
        assert 'frozen 1' in (cg/'cgroup.events').read_text()
        memory=int((cg/'memory.current').read_text());assert memory+LIMIT+300_000_000<=50_000_000_000,'Insufficient safe certificate budget'
        put('execution-gate.json',{'state':'frozen_production_for_exact_certificate','production_pid':int(pid),'production_unit':UNIT,'production_memory_at_freeze':memory,'certificate_memory_max':LIMIT,'combined_bound_with_reserve':memory+LIMIT+300_000_000,'threads':16,'controller_sha256':sha(__file__),'selected_export_manifest_sha256':sha(Q/'manifest.json'),'unmodified_candidate_six_plus_two':True,'automatic_thaw_on_all_exits':True,'created':time.time()})
        seam=measured('exact-target-seam',[CERT,'16','4',Q,S/'FEC',S/'LEC/MHV',S/'data'])
        records=re.findall(r'^SHIFTED_JOINT_CERTIFICATE (.+)$',Path(seam['log']).read_text(),re.M);assert len(records)==1
        certificate=json.loads(records[0]);assert certificate['exact_QQ_original_target_seam'] and certificate['shift_source_hash_receipts_verified']
        assert certificate['directions']==0 and certificate['columns']==518048 and certificate['derived_cache_congruence_primes']==8
        assert [certificate[x] for x in ['source_h','source_k','target_h','target_k']]==[17,3,16,4]
        origin=measured('exact-origin-holdouts',[ORIGIN,'16','4',Q,S/'FEC',S/'LEC/MHV',CASES])
        text=Path(origin['log']).read_text();assert re.findall(r'^ORIGIN_CHECK PASS words=(\d+)$',text,re.M)==['44']
        assert len(re.findall(r'^ORIGIN_WORD loop=10 ',text,re.M))==44
        check()
        pinned=[Q/'manifest.json',Q/'qq_complete.tsv',Q/'reconstruction_primes.tsv',E/'ten-install-receipt.json',E/'validation-passed.json',E/'ten-independent-verify.log',MIG,CERT,ORIGIN,CASES,Path(__file__)]
        final={'schema':'hexagon-ten-conditional-asymmetric-certificate-v1','qq_root':str(Q),'created':time.time(),'rational_reconstruction':{'CRT_primes':6,'fresh_validation_primes':2,'candidate_sha256':manifest['candidate_sha256'],'candidate_frozen_before_validation_primes_completed':True},'exact_original_target_seam':certificate,'exact_origin_holdouts':44,'seam_attempt':seam,'origin_attempt':origin,'artifact_hashes':{str(p):sha(p) for p in pinned},'imported_prime_provenance':'Prime1 computed by V16; explicit validated V17 import. Primes2-8 computed by V17. Frozen asymmetric QQ export is a separate producer.','migration_receipt':str(MIG),'bonus_assumptions':'MHV double-final bonus relations, Dixon-Liu 2308.08199 eq76-79, full D3 orbit; extrapolated to L10','unrestricted_uniqueness_claimed':False,'physical_scope':'Exact original 16+4 seam and44 exact origin holdouts. Native primes passed supplied finite sampled collinear/MRK/OPE checks; not a proof of every physical identity over QQ.','report_update_pending':True}
        put('final-verification.json',final)
        put('progress.json',{'state':'conditional_candidate_certified_report_pending','updated':time.time(),'final_symbol_certified':True,'final_verification':str(C/'final-verification.json'),'report_update_pending':True})
    except BaseException as error:
        put('failure.json',{'time':time.time(),'error':repr(error),'final_symbol_certified':False,'action':'Preserve evidence; diagnose before retry. Production auto-thaws.'});raise
    finally:
        if frozen:subprocess.run(['systemctl','--user','thaw',UNIT],check=True)

if __name__=='__main__':main()
