#!/usr/bin/env python3
"""Build exact factorized four-loop regularity RHS from reconstructed lower loops."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from itertools import combinations
from functools import lru_cache
import hashlib,json,random,time,resource
from ff_four_loop_H8_oracle import H8Oracle,splits,take
base=Path(__file__).resolve().parents[1];out=base/'results/ff_four_loop'

def clean(a):return {w:c for w,c in a.items() if c}
def add(a,b,s=Q(1)):
    d=defaultdict(Q,a)
    for w,c in b.items():d[w]+=s*c
    return clean(d)
def shuffle(a,b):
    d=defaultdict(Q);n=len(next(iter(a)));m=len(next(iter(b)))
    for u,cu in a.items():
        for v,cv in b.items():
            c=cu*cv
            for s,t in splits(n+m,n):
                w=[0]*(n+m)
                for i,x in zip(s,u):w[i]=x
                for i,x in zip(t,v):w[i]=x
                d[tuple(w)]+=c
    return clean(d)
def serial(a,n):
    return {'dimensions':[n]*len(next(iter(a))),'rules':[[list(w),str(c)] for w,c in sorted(a.items())]}
def shuffled_word(a,b,rng):
    s=rng.choice(splits(len(a)+len(b),len(a)));w=[0]*(len(a)+len(b))
    for i,x in zip(s[0],a):w[i]=x
    for i,x in zip(s[1],b):w[i]=x
    return tuple(w)

def run():
    start=time.monotonic();reports=[]
    for tag in ['DC','TC']:
        inp=out/f'H8_{tag}_lower_boundary.json';j=json.loads(inp.read_text());n=j['boundary_size'];finite=set(j['finite_letters'])
        factors={k:{tuple(w):Q(c) for w,c in v['rules']} for k,v in j['factors'].items()}
        e1,f2,f3=(factors[k] for k in ['E1','F2','F3'])
        square=shuffle(e1,e1);cube=shuffle(e1,square);e1f2=shuffle(e1,f2)
        r2=add(f2,square,-Q(1,2));r3=add(add(f3,e1f2,-1),cube,Q(1,3))
        assert r2 and r3
        assert all(all(i in finite for i in w) for w in r2)
        assert all(all(i in finite for i in w) for w in r3)
        name=out/f'H8_{tag}_factors.json'
        record={'boundary':tag,'boundary_size':n,'finite_letters':sorted(finite),
          'factor_origin':'Analytic F1 and independently reconstructed F2/F3, freshly projected; no known finite targets read.',
          'factors':{k:serial(a,n) for k,a in [('E1',e1),('R2',r2),('R3',r3)]},
          'spurious_formula':[{'coefficient':'1','shuffle_factors':['E1','R3']},{'coefficient':'1/2','shuffle_factors':['E1','E1','R2']},{'coefficient':'1/24','shuffle_factors':['E1','E1','E1','E1']}],
          'wholly_finite_correction':{'coefficient':'1/2','shuffle_factors':['R2','R2']},
          'quotient':'Only words with every letter in finite_letters are quotiented out; every mixed spurious word remains.',
          'lower_boundary_input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest()}
        name.write_text(json.dumps(record,separators=(',',':'))+'\n')
        oracle=H8Oracle(name)
        # Independent labeled-factor shuffle recursion for the original H8 formula.
        originals={'E1':e1,'F2':f2,'F3':f3}
        @lru_cache(maxsize=100000)
        def coeff(names,w):
            if len(names)==1:return originals[names[0]].get(w,Q(0))
            a=originals[names[0]];wt=len(next(iter(a)))
            return sum((a.get(take(w,s),Q(0))*coeff(names[1:],take(w,t)) for s,t in splits(len(w),wt)),Q(0))
        def direct(w):
            return (coeff(('E1','F3'),w)+coeff(('F2','F2'),w)/2
                    -coeff(('E1','E1','F2'),w)+coeff(('E1','E1','E1','E1'),w)/4)
        rng=random.Random(20261010);words=set();e=list(e1);r=list(r3);ff=list(f2)
        for _ in range(80):
            words.add(shuffled_word(rng.choice(e),rng.choice(r),rng))
            words.add(shuffled_word(rng.choice(ff),rng.choice(ff),rng))
            w=()
            for _ in range(4):w=shuffled_word(w,rng.choice(e),rng)
            words.add(w)
        sampled=[];nonzero=0;spurious=0
        for w in sorted(words):
            a=direct(w);b=oracle.coefficient(w,quotient=False);assert a==b,(tag,w,a,b)
            q=oracle.coefficient(w,quotient=True);is_spur=not all(i in finite for i in w)
            assert q==(a if is_spur else 0)
            nonzero+=bool(a);spurious+=is_spur
            if a and len(sampled)<12:sampled.append({'word':w,'full_H8':str(a),'spurious_H8':str(q)})
        # Exhaustive lower-weight control: matching formula versus the fully expanded square/cube.
        assert all(oracle.pairings(w)==c/2 for w,c in square.items())
        assert all(oracle.pairings(w)==c/6 for w,c in cube.items())
        reports.append({'boundary':tag,'factor_word_counts':{k:len(a) for k,a in [('E1',e1),('R2',r2),('R3',r3)]},
          'all_R2_R3_mixed_spurious_coefficients_cancel_exactly':True,
          'independent_original_formula_words_checked':len(words),'nonzero_checked':nonzero,'mixed_spurious_checked':spurious,
          'matching_formula_square_cube_nonzero_words_checked':len(square)+len(cube),'sample_coefficients':sampled,
          'factor_json_sha256':hashlib.sha256(name.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start})
        print(json.dumps(reports[-1]),flush=True)
        coeff.cache_clear();oracle.pairings.cache_clear()
    report={'status':'PASS','reports':reports,'known_finite_targets_read':False,'published_four_point_symbols_read':False,
      'weight8_tensor_materialized':False,'seconds':time.monotonic()-start,
      'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
      'scope':'Exact factor construction and exhaustive R2/R3 boundary regularity; sampled independent weight8 coefficient audit. Algebraic shuffle identity proves the factorization for every word.'}
    (out/'H8_factor_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
