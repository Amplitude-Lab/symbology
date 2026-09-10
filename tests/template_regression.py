"""Run bundled examples and independently certify their exact output spaces."""
from collections import defaultdict
from fractions import Fraction
import json
import shutil
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import zlib

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'front-end/server'))
from app import storage, templates, compile as compiler, jobs


def read(path):
    text=subprocess.check_output([str(ROOT/'tests/numerical_probe'),'dump',str(path)],text=True)
    lines=text.splitlines();dims=tuple(map(int,lines[0].split()[1:]));entries={}
    for line in lines[1:]:
        left,right=line.split(' = ')
        entries[tuple(map(int,left.split()[1:]))]=Fraction(right)
    return dims,entries


def rank_mod(rows, prime=1000000007):
    # Independent sparse Gaussian elimination over a prime field. A modular
    # rank is a lower bound on rational rank. Combined with exact residuals
    # and independent solution rows below it certifies completeness over QQ.
    pivots={}
    for values in rows:
        row={i:int(v.numerator)*pow(int(v.denominator),-1,prime)%prime for i,v in values.items() if v}
        row={i:v for i,v in row.items() if v}
        while row:
            lead=min(row)
            if lead not in pivots:
                inv=pow(row[lead],-1,prime)
                pivots[lead]={i:v*inv%prime for i,v in row.items()}
                break
            scale=row[lead]
            for i,v in pivots[lead].items():
                value=(row.get(i,0)-scale*v)%prime
                if value:row[i]=value
                else:row.pop(i,None)
    return len(pivots)


def certify(root,tid):
    n=templates.TEMPLATES[tid]['n_letters']
    dfile='dlogmat_pentagon.wxf' if tid=='pentagon' else 'dlogmat_full_4pformfactor.wxf'
    _,D=read(root/'data'/dfile)
    if tid == '4pformfactor':
        dd,base=read(root/'data/dlogmat_4pformfactor.wxf')
        _,steinmann=read(root/'data/dlogmatES_4pformfactor.wxf')
        combined={**base, **{(i,j,c+dd[2]):v for (i,j,c),v in steinmann.items()}}
        assert D == combined, 'Precomputed full conditions differ from the original two tensors'
    dims,T=read(root/'output/cyclic_flip_symbols.wxf')
    assert dims[1:]==(n,n)
    by_pair=defaultdict(list);conditions=defaultdict(dict)
    for (i,j,c),v in D.items():
        by_pair[i,j].append((c,v));conditions[c][i*n+j]=v
    residual=defaultdict(Fraction);basis=[{} for _ in range(dims[0])]
    for (k,i,j),v in T.items():
        basis[k][i*n+j]=v
        for c,d in by_pair[i,j]:residual[k,c]+=v*d
    assert not any(residual.values()),'Original integrability/Steinmann residual is nonzero'
    all_conditions=list(conditions.values())
    matrices=[]
    for filename in ('cycmat.wxf','flipmat.wxf'):
        md,M=read(root/'data'/filename);assert md==(n,n);matrices.append(M)
        m=defaultdict(list)
        for (i,j),v in M.items():m[i].append((j,v))
        transformed=defaultdict(Fraction)
        for (k,i,j),v in T.items():
            for a,x in m[i]:
                for b,y in m[j]:transformed[k,a,b]+=v*x*y
        assert {i:v for i,v in transformed.items() if v}==T,filename+' invariance failed'
        symrows=[defaultdict(Fraction) for _ in range(n*n)]
        for i in range(n):
            for j in range(n):
                symrows[i*n+j][i*n+j]-=1
                for a,x in m[i]:
                    for b,y in m[j]:symrows[a*n+b][i*n+j]+=x*y
        all_conditions.extend(symrows)
    if tid == 'pentagon':
        identity={(i,i):Fraction(1) for i in range(n)}
        def multiply(A,B):
            by_row=defaultdict(list);out=defaultdict(Fraction)
            for (i,j),v in B.items():by_row[i].append((j,v))
            for (i,j),v in A.items():
                for k,w in by_row[j]:out[i,k]+=v*w
            return {p:v for p,v in out.items() if v}
        C,F=matrices;power=identity
        for _ in range(5):power=multiply(power,C)
        assert power==identity and multiply(F,F)==identity
        inverse=multiply(multiply(multiply(C,C),C),C)
        assert multiply(multiply(F,C),F)==inverse
    r=rank_mod(all_conditions);k=rank_mod(basis)
    assert k==dims[0] and r+k==n*n,(r,k,dims)
    print(f'PASS {tid}: {k} independent invariant symbols; exact original-condition and cyclic/flip residuals zero; full space certified ({r}+{k}={n*n})',flush=True)


def symmetry_edge_cases(work):
    from app.template_seeds import write_word_basis
    def call(*argv, code=0):
        r=subprocess.run(list(map(str,argv)),capture_output=True,text=True,timeout=20)
        assert r.returncode==code,r.stdout+r.stderr
        return r.stdout
    call(ROOT/'tests/numerical_probe','gen',work)
    write_word_basis(work/'words.wxf',2)
    call(ROOT/'tensor_ops','power',work/'swap.wxf',0,work/'identity.wxf')
    call(ROOT/'tensor_add',work/'identity.wxf',work/'identity.wxf',-1,0,work/'minus.wxf')
    call(ROOT/'tensor_ops','symsolve',work/'words.wxf',work/'identity.wxf',work/'minus.wxf',work/'empty.wxf')
    assert 'rank 3 0 2 2 nnz 0' in call(ROOT/'tensor_ops','dims',work/'empty.wxf')
    call(ROOT/'tensor_ops','symsolve',work/'empty.wxf',work/'identity.wxf',work/'minus.wxf',work/'empty-again.wxf')
    assert (work/'empty.wxf').read_bytes()==(work/'empty-again.wxf').read_bytes()
    call(ROOT/'tensor_ops','symsolve',work/'tensor.wxf',work/'identity.wxf',work/'swap.wxf',work/'bad.wxf',code=1)
    assert not (work/'bad.wxf').exists()
    print('PASS symmetry edge cases: empty invariant space, empty input, and genuine failure of closure',flush=True)


with tempfile.TemporaryDirectory(prefix='symbology-templates ') as tmp:
    symmetry_edge_cases(Path(tmp))
    storage.PROJECTS_DIR=Path(tmp);storage.init_dirs()
    for tid in ('pentagon','4pformfactor','e6'):
        project=storage.create_project(tid);templates.apply_template(project,tid);storage.save_project(project)
        flow=project['flows'][-1];plan=compiler.compile_flow(project,flow['graph']);assert plan['ok'],plan['errors']
        engine=jobs.Engine();run=engine.create_run(project['id'],flow['id'],flow['name'],plan['_steps_full'])
        until=time.monotonic()+120
        while run.worker_active and time.monotonic()<until:time.sleep(.05)
        root=storage.project_dir(project['id'])
        if run.worker_active:engine.cancel(run.run_id)
        assert run.status=='done',(root/'runs'/f'{run.run_id}.log').read_text()[-2000:]
        if tid!='e6':certify(root,tid)
        else:
            ref=json.loads((ROOT/'tests/baselines/2026-09-09/manifest-computerhs.json').read_text())['crc32']
            for actual,expected in [('output/boundary_2L.wxf','output/2loop/boundary_2L.wxf'),('output/2loop/solMHV_2L.wxf','output/2loop/solMHV_2L.wxf')]:
                assert f'{zlib.crc32((root/actual).read_bytes()):08x}'==ref[expected]
            print('PASS E6 two-loop example: all 8 steps succeed; boundary and solution match recorded exact reference',flush=True)
            exported=compiler.export_flow_script(project,flow['graph'],"E6 portable")
            assert exported['ok'],exported
            relocated=Path(tmp)/"relocated E6's project"
            shutil.copytree(root,relocated)
            result=subprocess.run(['bash',str(relocated/'exported'/Path(exported['path']).name)],
                                  env={**os.environ,'PROJ_DIR':str(relocated),'WOLFRAMSCRIPT':'/no/such/wolfram'},
                                  capture_output=True,text=True,timeout=120)
            assert result.returncode==0,result.stdout[-2000:]+result.stderr
            for actual,expected in [('output/boundary_2L.wxf','output/2loop/boundary_2L.wxf'),('output/2loop/solMHV_2L.wxf','output/2loop/solMHV_2L.wxf')]:
                assert f'{zlib.crc32((relocated/actual).read_bytes()):08x}'==ref[expected]
            print('PASS relocated standalone E6 flow: spaces/apostrophes, no Wolfram, exact boundary and solution',flush=True)
