from pathlib import Path
from fractions import Fraction as Q
import json
O=Path(__file__).resolve().parents[2]/'results/ff_four_loop';j=json.load(open(O/'H8_tail_factors.json'))
A={(i,k):Q(v)for(i,k),v in j['annTR']['rules']};C={(i,k):Q(v)for(i,k),v in j['prefix_log_factors']['rules']};B={(i,k):Q(v)for(i,k),v in j['annihilator_component_factors']['rules']}
D=[[sum(C.get((i,z),0)*A.get((r,i),0)for i in range(1,89))for r in range(1,77)]for z in range(1,4)]
def rank(rows):
 piv={}
 for raw in rows:
  a=raw.copy()
  for p,b in sorted(piv.items()):
   q=a[p];a=[x-q*y for x,y in zip(a,b)]
  if any(a):
   p=next(i for i,x in enumerate(a)if x);q=a[p];piv[p]=[x/q for x in a]
 return len(piv)
br=[[B.get((i,k),Q(0))for k in range(1,77)]for i in range(1,4)]
r={'universal_log_space_dimension':3,'quotient_map_on_log_space_rank':rank(D),'projected_one_loop_tail_rank':rank(br),'projected_log_rows':[[[k+1,str(v)]for k,v in enumerate(row)if v]for row in D],'no_absorbing_N_in_Sym2L':rank(D)<rank(br),'proof':'If N has symmetric matrix C M C^T, projectedtail is C M D with D=C^T annTR^T. Its rank is at mostrankD=1. ProjectedS1tail has rank3; equality is impossible.','scope':'Rules out only a quadraticnormalization polynomial in the three universallyallowedlogs. No broadernormalizationno-go asserted.','four_loop_target_read':False};(O/'tail_normalization_reduction.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
u=[[A.get((r,i),Q(0))for r in range(1,77)]for i in range(1,5)]
r['larger_geometric_first_entry_kernel']={'basis':['W1','W2','W3','W4'],'dimension':4,'quotient_rank':rank(u),'no_absorbing_N_in_Sym2U':rank(u)<rank(br),'scope':'The full first-entry intersection ker(q) is the span of these four logs, certified in tail_prefix_gate.json. Any N in its symmetric square has projected tail rank at most2, whereas the analytic one-loop source has rank3.'}
(O/'tail_normalization_reduction.json').write_text(json.dumps(r,indent=2)+'\n');print(r['larger_geometric_first_entry_kernel'])
