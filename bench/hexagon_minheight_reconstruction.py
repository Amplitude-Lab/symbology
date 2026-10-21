#!/usr/bin/env python3
"""Retrospective minimum-product reconstruction experiment, never a fresh check.

Select the least |numerator|*denominator Euclidean convergent with 2ND<M.
The rule is heuristic absent coefficient bounds; exhaustive exact-reference
and withheld-cache comparisons below assess this data set, not a theorem
of physical uniqueness. Zero is handled canonically.
"""
import json,math,time
from pathlib import Path
from array import array
from hexagon_ten_loop_benchmark import sha
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009'
def reconstruct(r,m):
 if r==0:return (0,1)
 a,b=m,r;u,v=0,1;best=None;score=None;tie=False
 while b:
  if v:
   n,d=(b,v) if v>0 else (-b,-v);height=abs(n)*d
   if 2*height<m and (score is None or height<=score):
    if height==score and best!=(n,d):tie=True
    elif score is None or height<score:best=(n,d);score=height;tie=False
  q=a//b;a,b=b,a-q*b;u,v=v,u-q*v
 return None if tie else best
def vector(path):
 out=[]
 for line in path.read_text().splitlines():
  if not line.startswith('#'):
   j,n,d=map(int,line.split());assert j==len(out);out.append((n,d))
 return out
def main():
 start=time.monotonic();receipt=json.loads((V/'native-reconstruction-7p96.json').read_text());refs=vector(V/'native-reconstruction-7p96.tsv');ps=[x['prime'] for x in receipt['cache_receipts']];caches=[]
 for c in receipt['cache_receipts']:
  p=Path(c['path']);assert sha(p)==c['sha256'];a=array('Q');a.frombytes(p.read_bytes()[72:]);caches.append(a)
 m=ps[0];res=list(caches[0])
 for p,c in zip(ps[1:6],caches[1:6]):
  inv=pow(m,-1,p);res=[a+m*((b-a)*inv%p) for a,b in zip(res,c)];m*=p
 result=[reconstruct(r,m) for r in res];mismatches=sum(a!=b for a,b in zip(result,refs));checks=[sum(q is None or q[1]%p==0 or (q[0]-c[j]*q[1])%p for j,q in enumerate(result)) for p,c in zip(ps[6:],caches[6:])]
 assert mismatches==0 and not any(checks)
 controls=[]
 for name in ['eight','nine']:
  path=R/f'output_hexagon_candidate_export_20261009/{name}-reference.canonical.tsv';v=vector(path);bad=0
  for n,d in v:
   r=(n*pow(d,-1,m))%m;bad+=reconstruct(r,m)!=(n,d)
  assert bad==0;controls.append({'reference':str(path),'sha256':sha(path),'coefficients':len(v),'mismatches':bad})
 out=V/'native-minheight-six-prime.tsv'
 with out.open('x') as f:
  f.write('# native17+3 minimum-product six-prime retrospective experiment\n')
  for j,(n,d) in enumerate(result):f.write(f'{j}\t{n}\t{d}\n')
 record={'schema':'hexagon-minheight-reconstruction-experiment-v1','passed':True,'source_sha256':sha(Path(__file__)),'modulus_bits':m.bit_length(),'CRT_primes':6,'retrospective_selection_or_check_primes':2,'native_coefficients':len(result),'native_reference_mismatches':mismatches,'holdout_mismatches':checks,'controls':controls,'maximum_individual_2ND_bits':max((2*abs(n)*d).bit_length() for n,d in result),'maximum_global_box_2ND_bits':(2*max(abs(n) for n,d in result)*max(d for n,d in result)).bit_length(),'wall_seconds':time.monotonic()-start,'vector_sha256':sha(out),'required_future_protocol':'Freeze this algorithm before future independent validation primes. Retrospective success does not relabel previously computed primes as fresh. This experiment does not replace the accepted seven-prime native export provenance.'}
 (V/'minheight-experiment.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
