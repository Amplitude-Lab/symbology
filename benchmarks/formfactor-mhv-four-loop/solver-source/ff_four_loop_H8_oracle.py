#!/usr/bin/env python3
"""Exact lazy four-loop affine boundary coefficients, with no weight-eight tensor.

Input tensors have one-based boundary-letter indices. The spurious quotient is
by the span of words all of whose letters are finite; a mixed word is retained.
Only reconstructed lower loops and analytic F1 occur in the factor file.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse,json

@lru_cache(None)
def splits(n,k):
    return tuple((s,tuple(i for i in range(n) if i not in s)) for s in combinations(range(n),k))

def take(w,s): return tuple(w[i] for i in s)

class H8Oracle:
    def __init__(self,path):
        self.path=Path(path)
        self.data=json.loads(self.path.read_text())
        self.n=self.data['boundary_size']
        self.finite=frozenset(self.data['finite_letters'])
        self.factors={k:{tuple(w):Q(c) for w,c in v['rules']} for k,v in self.data['factors'].items()}
        self.e1=self.factors['E1'];self.r2=self.factors['R2'];self.r3=self.factors['R3']

    @lru_cache(maxsize=200000)
    def pairings(self,w):
        """Coefficient of E1^shuffle(n/2)/(n/2)! on the even word w.

        A perfect matching gives an unordered collection of identical factors;
        every factor reads its two slots in their original order. No symmetry
        of the E1 tensor is assumed.
        """
        if not w:return Q(1)
        ans=Q(0)
        for j in range(1,len(w)):
            c=self.e1.get((w[0],w[j]),Q(0))
            if c:ans+=c*self.pairings(w[1:j]+w[j+1:])
        return ans

    def remainder_product(self,w):
        """Coefficient of R2 shuffle R2 / 2, including finite words."""
        return sum((self.r2.get(take(w,s),Q(0))*self.r2.get(take(w,t),Q(0))
                    for s,t in splits(8,4)),Q(0))/2

    def coefficient(self,word,quotient=True):
        w=tuple(word)
        if len(w)!=8 or any(i<1 or i>self.n for i in w):raise ValueError('expected eight valid one-based boundary letters')
        wholly_finite=all(i in self.finite for i in w)
        if quotient and wholly_finite:return Q(0)
        ans=sum((self.e1.get(take(w,s),Q(0))*self.r3.get(take(w,t),Q(0))
                 for s,t in splits(8,2)),Q(0))
        ans+=sum((self.r2.get(take(w,s),Q(0))*self.pairings(take(w,t))
                  for s,t in splits(8,4)),Q(0))
        ans+=self.pairings(w)
        if not quotient:ans+=self.remainder_product(w)
        return ans

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('boundary',choices=['DC','TC']);p.add_argument('letters',type=int,nargs=8)
    p.add_argument('--full',action='store_true',help='Include wholly finite words and R2 shuffle R2 / 2')
    a=p.parse_args();base=Path(__file__).resolve().parents[1]/'results/ff_four_loop'
    o=H8Oracle(base/f'H8_{a.boundary}_factors.json')
    print(o.coefficient(a.letters,quotient=not a.full))

if __name__=='__main__':main()
