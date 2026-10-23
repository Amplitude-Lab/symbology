#!/usr/bin/env python3
"""Freeze an adjacency-only generation problem before weight-eight discovery."""
from pathlib import Path
import sys,json,hashlib,time,resource
from fractions import Fraction as Q
BASE=Path(__file__).resolve().parents[2];sys.path.insert(0,str(BASE/'code'));import ff_es_cluster_compare as m
out=BASE/'results/ff_four_loop';t=time.monotonic();N=100
cs=json.loads((BASE/'results/ff_ca_extras_candidate_spaces.json').read_text());sp=json.loads((BASE/'results/ff_ca_extras_colored_limits_spaces.json').read_text());nd=json.loads((BASE/'results/ff_endpoint_node_dictionary.json').read_text());maps=json.loads((BASE/'results/ff_exchange_limit_pair_maps.json').read_text());allowed=set(sp['allowed']);W={idx:{i:Q(v)for i,v in enumerate(row)if Q(v)}for idx,row in zip(m.inp['rational_indices'],nd['W_in_cluster_coordinates'])}
for k,i in enumerate([45,46,55,56]):W[i]={76+k:Q(1)}
for k,i in enumerate(maps['q_indices']):W[i]={80+k:Q(1)}
def fold(row):
 z={}
 for k,v in row.items():
  if k not in allowed:m.add(z,{N*min(k//N,k%N)+max(k//N,k%N):v})
 return z
sym=m.basis(fold({int(k):Q(v)for k,v in row})for row in cs['interior_generators']);assert len(sym)==281
ids=sorted(W);assert len(ids)==88
rows=[]
for i in ids:
 for j in ids:
  rows.append(m.fullreduce(fold({N*a+b:x*y for a,x in W[i].items()for b,y in W[j].items()}),sym))
cols=sorted({k for r in rows for k in r});assert len(cols)==2047
ci={k:i+1 for i,k in enumerate(cols)}
rules=[[[i//88+1,i%88+1,ci[k]],str(v)]for i,row in enumerate(rows)for k,v in sorted(row.items())]
(out/'data/adjacency_symmetric.json').write_text(json.dumps({'dimensions':[88,88,len(cols)],'rules':rules},separators=(',',':'))+'\n')
pure=json.loads((BASE/'results/ff_ca_extras_pure60.json').read_text());pr=[{}for _ in range(60)]
for(a,b,c),v in pure['rules']:
 if Q(v):pr[a-1][b,c]=Q(v)
fec=[m.combine({int(k):Q(v)for k,v in rr},pr)for rr in cs['firstpair_prefix_kernel']];assert len(fec)==50
(out/'data/first_pair.json').write_text(json.dumps({'dimensions':[50,8,88],'rules':[[[r+1,b,c],str(v)]for r,row in enumerate(fec)for(b,c),v in row.items()]},separators=(',',':'))+'\n')
manifest={'objective':'Weight8 finite four-point form-factor bootstrap and prospective test of proposed adjacency; fixed retained88alphabet.','registered_constraints':{'first_entries':'W1..W8','first_pair':'exact proposed50-dimensional integrable prefix','all_later_pairs':'proposed symmetric adjacency annihilator2047 plus exact pure wedge tensor','symmetries':'physical D4,parity,fourquadratic Galois','parent_last_entry':'unrestricted88; no literature tail hypothesis','normalization':'analytic F1; logarithm subtraction H8=F1 shuffle F3+1/2 F2 shuffle F2-F1 shuffle F1 shuffle F2+1/4 F1 shuffle F1 shuffle F1 shuffle F1','boundaries':'all slotwise spurious double/triple-collinear components of F4-H8 vanish; finite targets not imposed','external_ES':'not added; implied by proposed pairspaces'},'limits':{'shared_memory_bytes':50000000000,'initial_workers':4,'own_process_stop_bytes':44000000000},'schedule':'Generate exact adjacency prefixes; profile cost/fill. Keep complete DC/TC lazy generators registered while scheduling algebraic symmetry/seam/regularity blocks by measured cost. Preserve affine origins and all checkpoints.','inputs':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [BASE/'results/ff_ca_extras_candidate_spaces.json',BASE/'results/ff_exchange_complete_generators.json',BASE/'results/ff_freeze_w6/solution_w6.wxf',BASE/'results/ff_boot_boundary_map.wxf',BASE/'results/ff_boot_triple_map.wxf'] if p.exists()},'preparation':{'adjacency_symmetric_rank':2047,'adjacency_nonzeros':len(rules),'firstpair_dimension':50,'seconds':time.monotonic()-t,'maximum_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}}
(out/'execution-gate.json').write_text(json.dumps(manifest,indent=2)+'\n');print(manifest['preparation'])
