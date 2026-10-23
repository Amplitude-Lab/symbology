#!/usr/bin/env python3
"""Exact local test of the affine last-entry RHS, with no weight-eight target."""
from pathlib import Path
from fractions import Fraction as Q
import json,sys,time,hashlib
B=Path(__file__).resolve().parents[2];O=B/'results/ff_four_loop';sys.path.insert(0,str(B/'code'));import ff_es_cluster_compare as m
start=time.monotonic();a=json.loads((O/'data/adjacency_symmetric.json').read_text());A={}
for (i,j,k),v in a['rules']:A.setdefault((i,j),{})[k]=Q(v)
f=json.loads((O/'data/first_pair.json').read_text());fr=[{}for _ in range(50)]
for (r,i,j),v in f['rules']:fr[r-1][88*(i-1)+j-1]=Q(v)
fb=m.basis(fr);logs=[{1:Q(1),2:Q(-1)},{1:Q(-1),3:Q(1)},{1:Q(1),4:Q(-1)}]
def tensor(l,i):
 d={}
 for k,v in l.items():m.add(d,{88*(k-1)+i-1:v});m.add(d,{88*(i-1)+k-1:v})
 return d
R=[];first=[]
for l in logs:
 rows=[]
 for i in range(1,89):
  r={}
  for k,v in l.items():
   for pair in [(k,i),(i,k)]:m.add(r,{q:v*x for q,x in A.get(pair,{}).items()})
  rows.append(r)
 R.append(rows);first.append([m.fullreduce(tensor(l,i),fb)for i in range(1,9)])
result={'identity':'Pi_TR d_last H8=sum_i(F3 shuffle li) tensor Bi, li=W1-W2,-W1+W3,W1-W4','local_reduction':'Insertion outside a tested pair is annihilated by the corresponding F3 pair condition. Remaining contribution is sym(li tensor a), where a is a contextual single-letter slice of F3.','interior_universal_log_compatibility':[all(not r for r in rows)for rows in R],'first_pair_universal_on_first8':[all(not r for r in rows)for rows in first],'interior_nonzero_singleletter_maps':[[i+1 for i,r in enumerate(rows)if r]for rows in R],'first_nonzero_maps':[[i+1 for i,r in enumerate(rows)if r]for rows in first],'candidate_first_pair_dimension':50,'candidate_interior_symmetric_annihilator':2047,'F3_word_tensor_read':False,'F4_target_read':False,'seconds':time.monotonic()-start,'input_sha256':{str(p.relative_to(B)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [O/'H8_tail_factors.json',O/'data/adjacency_symmetric.json',O/'data/first_pair.json']}}
result['necessary_tail_prefix_gate_pass']=all(result['interior_universal_log_compatibility'])and all(result['first_pair_universal_on_first8']);(O/'tail_prefix_gate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

nd=json.loads((B/'results/ff_endpoint_node_dictionary.json').read_text());geo=json.loads((B/'results/ff_ca_extras_classify.json').read_text());short={i-1 for i in geo['distinguished_node_coordinates']};keep=[i for i in range(64) if i not in short]
wr=[{i:Q(nd['W_in_cluster_coordinates'][k][i]) for i in keep if Q(nd['W_in_cluster_coordinates'][k][i])} for k in range(8)]
u=m.kernel(wr);zero=[not m.combine({k-1:v for k,v in l.items()},wr) for l in logs]
result['geometric_explanation']={'q_first_entry_rank':len(m.basis(wr)),'first_entry_kernel_dimension':len(u),'kernel_basis':[[[i+1,str(v)] for i,v in sorted(row.items())] for row in u],'three_logs_map_zero_under_q':zero,'conclusion':'The three tail logs span a 3-dimensional subspace of the 4-dimensional first-entry kernel of q. Both ordered tensors with one of these logs as a factor lie in each pullback adjacency space.'}
(O/'tail_prefix_gate.json').write_text(json.dumps(result,indent=2)+'\n');print(result['geometric_explanation'])
