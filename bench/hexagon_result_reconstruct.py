#!/usr/bin/env python3
"""Reconstruct a saved native or shifted result without running new solves.

The requested final caches are excluded from CRT. They are retrospective holdouts,
not newly generated primes. Independent exact split comparison is required.
"""
from array import array
import argparse, json, math, os, resource, struct, sys, time
from pathlib import Path
from hexagon_crt_asymmetric_probe import reconstruct, controls
from hexagon_crt_progress_audit import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True);p.add_argument('--kind',choices=['native','conversion'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--crt-primes',type=int,default=6)
    p.add_argument('--numerator-bits',type=int,default=48);p.add_argument('--check-primes',type=int,default=2);a=p.parse_args()
    assert a.crt_primes>=1 and a.check_primes>=0
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30));resource.setrlimit(resource.RLIMIT_CPU,(180,180));os.nice(5)
    start=time.monotonic();checked=controls();ledger=json.loads((a.control/'ledger.json').read_text())
    h,k=(17,3) if a.kind=='native' else (16,4);cols=275246 if a.kind=='native' else 518048
    caches=[];primes=[];receipts=[]
    for i in range(1,a.crt_primes+a.check_primes+1):
        e=ledger['primes'][str(i)][a.kind];path=Path(e['path']);digest=sha(path)
        assert digest==e['sha256']==Path(str(path)+'.sha256').read_text().strip()
        with path.open('rb') as f:
            header=struct.unpack('<8s8Q',f.read(72));assert header==(b'NMHVSEW\0',1,10,h,k,1,cols,e['prime'],int(sha(path.parent.parent/'input_identity.txt')[:16],16))
            v=array('Q');v.fromfile(f,cols);assert not f.read(1)
        if sys.byteorder!='little':v.byteswap()
        assert all(x<header[7] for x in v);caches.append(v);primes.append(header[7]);receipts.append({'path':str(path),'sha256':digest,'prime':header[7],'role':'CRT' if i<=a.crt_primes else 'retrospective_holdout'})
    residues=list(caches[0]);modulus=primes[0]
    for prime,v in zip(primes[1:a.crt_primes],caches[1:a.crt_primes]):
        inv=pow(modulus,-1,prime);residues=[x+modulus*((y-x)*inv%prime) for x,y in zip(residues,v)];modulus*=prime
    N=(1<<a.numerator_bits)-1;D=(modulus-1)//(2*N);vector=a.output.with_suffix('.tsv')
    failed=[];bad=[0]*a.check_primes;nnz=hn=hd=0
    with vector.open('x') as out:
        out.write(f'# native {h}+{k} reconstruction; zero-based column, numerator, denominator\n')
        for j,r in enumerate(residues):
            q=reconstruct(r,modulus,N,D)
            if q is None:failed.append(j);continue
            n,d=q;assert d>0 and math.gcd(n,d)==1
            assert all(d%prime and (n-c[j]*d)%prime==0 for c,prime in zip(caches[:a.crt_primes],primes[:a.crt_primes]))
            for t,(c,prime) in enumerate(zip(caches[a.crt_primes:],primes[a.crt_primes:])):bad[t]+=int(d%prime==0 or (n-c[j]*d)%prime!=0)
            nnz+=n!=0;hn=max(hn,abs(n).bit_length());hd=max(hd,d.bit_length());out.write(f'{j}\t{n}\t{d}\n')
        out.flush();os.fsync(out.fileno())
    receipt={'schema':'hexagon-saved-cache-rational-reconstruction-v1','created':time.time(),'source_sha256':sha(Path(__file__)),'arithmetic_controls_passed':checked,'split':[h,k],'columns':cols,'nonzero':nnz,'numerator_bits':hn,'denominator_bits':hd,'numerator_bound_bits':a.numerator_bits,'denominator_bound_bits':D.bit_length(),'modulus_bits':modulus.bit_length(),'cache_receipts':receipts,'failed_columns':failed,'holdout_mismatches':bad,'vector':str(vector.resolve()),'vector_sha256':sha(vector),'passed':not failed and not any(bad),'scope':'Reconstruction and retrospective modular holdouts only; exact split equality, original sewing and physical checks separate. No claim of newly generated fresh primes.','wall_seconds':time.monotonic()-start}
    with a.output.open('x') as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps({k:receipt[k] for k in ['passed','columns','nonzero','numerator_bits','denominator_bits','failed_columns','holdout_mismatches','wall_seconds']},indent=2));assert receipt['passed']

if __name__=='__main__':main()
