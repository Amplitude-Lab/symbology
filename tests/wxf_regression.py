"""Bounded malformed-input coverage; compatible with the ASan/UBSan reader."""
import argparse
from pathlib import Path
import random
import struct
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(); p.add_argument('--probe',type=Path,default=ROOT/'tests/wxf_probe'); args=p.parse_args()

def vi(n):
    out=bytearray()
    while n>127: out.append((n&127)|128); n >>=7
    out.append(n); return bytes(out)
def sym(s): return b's'+vi(len(s))+s.encode()
def fun(s,n): return b'f'+vi(n)+sym(s)
def arr(shape,values,width=8):
    return bytes([193,{1:0,2:1,4:2,8:3}[width]])+vi(len(shape))+b''.join(vi(n) for n in shape)+b''.join(int(v).to_bytes(width,'little',signed=True) for v in values)
def integer(n):return b'C'+struct.pack('b',n)
def sparse(dims=(2,2),rows=(0,1,2),cols=(1,2),vals=None,nnz=2,width=4):
    vals=vals or fun('List',2)+integer(3)+fun('Rational',2)+integer(-2)+integer(3)
    return b'8:'+fun('SparseArray',4)+sym('Automatic')+arr([len(dims)],dims)+integer(0)+fun('List',3)+integer(1)+fun('List',2)+arr([len(rows)],rows)+arr([nnz,len(dims)-1],cols,width)+vals
big = b'I'+vi(101)+b'1'+b'0'*100
valid=[sparse(vals=fun('List',2)+big+fun('Rational',2)+big+integer(3)),sparse(),sparse(width=8),sparse(dims=(2,0),rows=(0,0,0),cols=(),nnz=0,vals=fun('List',0))]
invalid=[b'',b'8:',b'not wxf',b'8C:invalid compressed',b'8:L'+b'\0'*2,sparse()+integer(4),
 sparse(rows=(1,1,2)),sparse(rows=(0,2,1)),sparse(rows=(0,1)),sparse(cols=(0,2)),sparse(cols=(3,2)),sparse(dims=(-1,2)),
 sparse(vals=fun('List',2)+integer(3)+fun('Rational',2)+integer(1)+integer(0)),
 sparse(vals=fun('List',2)+integer(3)+b'I'+vi(3)+b'1x2'),
 sparse(vals=fun('List',1)+integer(3)), b'8:'+bytes([193,3])+vi(2)+vi(2**63)+vi(8), b'8:s'+b'\xff'*10]
raw=valid[0]; invalid += [raw[:i] for i in range(len(raw))]
rng=random.Random(73)
fuzz=[rng.randbytes(rng.randrange(1,256)) for _ in range(500)]
for _ in range(500):
    b=bytearray(raw)
    for _ in range(rng.randrange(1,5)): b[rng.randrange(len(b))]=rng.randrange(256)
    fuzz.append(bytes(b))
with tempfile.TemporaryDirectory(prefix='symbology-wxf ') as tmp:
    paths=[]
    for i,b in enumerate(valid+invalid+fuzz):
        f=Path(tmp)/f'{i}.wxf';f.write_bytes(b);paths.append(f)
    known=list((ROOT/'data').glob('*.wxf'))+list((ROOT/'data_pentagon').glob('*.wxf'))+list((ROOT/'data_4pformfactor').glob('*.wxf'))
    result=subprocess.run([str(args.probe.resolve()),*map(str,paths+known)],capture_output=True,text=True,timeout=60)
    assert result.returncode==0,(result.returncode,result.stdout[-2000:],result.stderr[-3000:])
    lines=result.stdout.splitlines()
    assert len(lines)==len(paths)+len(known),result.stdout[-2000:]
    assert all(s.startswith('OK ') for s in lines[:len(valid)]),lines[:len(valid)]
    assert all(s.startswith('REJECT ') for s in lines[len(valid):len(valid)+len(invalid)]),lines[len(valid):len(valid)+len(invalid)]
    assert all(s.startswith('OK ') for s in lines[-len(known):]),lines[-len(known):]
    assert not result.stderr,result.stderr
    print(f'PASS {len(valid)+len(known)} valid files, {len(invalid)} invalid files, {len(fuzz)} bounded mutations; no crash or sanitizer report')
