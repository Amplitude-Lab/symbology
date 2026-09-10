"""Native cache receipts must track content, including hidden seed dependencies."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--probe',type=Path,default=ROOT/'tests/cache_probe');args=p.parse_args()
with tempfile.TemporaryDirectory(prefix="symbology cache's ") as tmp:
    work=Path(tmp);files=[]
    for n in (0,1,3,55,56,63,64,65,127,128,129,1000000):
        f=work/str(n);f.write_bytes(b'a'*n);files.append(f)
    got=subprocess.check_output([str(args.probe.resolve()),*map(str,files)],text=True).splitlines()
    assert got==[hashlib.sha256(f.read_bytes()).hexdigest() for f in files]
    print('PASS 12 independent SHA256 checks')
    data=work/'data';out=work/'output';shutil.copytree(ROOT/'data',data)
    def compute(destination, executable=ROOT/'compute_rhs'):
        start=time.monotonic()
        r=subprocess.run([str(executable),'--target','SEW_5p1','--data-dir',str(data),
                          '--output-dir',str(destination),'--letter-projection',str(destination/'collinear/colprojdiv_w1.wxf')],
                         capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout[-2000:]+r.stderr
        return time.monotonic()-start,r.stdout
    cold,_=compute(out)
    protected=out/'FEC_4.wxf';before=protected.stat().st_mtime_ns
    warm,log=compute(out)
    assert protected.stat().st_mtime_ns==before
    assert 'Verified cached loop 3' in log
    # A rebuilt sibling executable invalidates native receipts even when the
    # current output happens to have the correct mathematical contents.
    binaries=work/'binaries';binaries.mkdir()
    for name in ('bootstrap','compute_rhs'):shutil.copy2(ROOT/name,binaries/name)
    with (binaries/'bootstrap').open('ab') as f:f.write(b'\0')
    compute(out,binaries/'compute_rhs')
    assert protected.stat().st_mtime_ns!=before
    # A valid WXF file with wrong content must not be accepted just because it exists.
    expected=protected.read_bytes();protected.write_bytes((out/'FEC_3.wxf').read_bytes())
    compute(out);assert protected.read_bytes()==expected
    # A formerly nonempty projection must disappear when recomputation is empty.
    fin=out/'collinear/colprojfin_SEW_3p1.wxf';assert not fin.exists()
    fin.write_bytes((data/'colprojfin.wxf').read_bytes())
    compute(out);assert not fin.exists()
    # Change a hidden data-directory input without trusting its timestamp.
    original=(data/'E1.wxf').stat()
    changed=work/'changed.wxf'
    r=subprocess.run([str(ROOT/'tensor_add'),str(data/'E1.wxf'),str(data/'E1.wxf'),'2','0',str(changed)],capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    (data/'E1.wxf').write_bytes(changed.read_bytes())
    os.utime(data/'E1.wxf',ns=(original.st_atime_ns,original.st_mtime_ns))
    compute(out)
    fresh=work/'fresh';compute(fresh)
    refs=json.loads((ROOT/'audits/baseline-2026-09-09/manifest-computerhs.json').read_text())['crc32']
    for rel in refs:
        rel=Path(rel).relative_to('output')
        assert (out/rel).read_bytes()==(fresh/rel).read_bytes(),rel
    assert (out/'oneloop/E1.wxf').read_bytes()==changed.read_bytes()
    print(f'PASS changed seeds/binary, corrupted cached tensor, stale empty projection, apostrophes/spaces; all {len(refs)} resumed outputs equal a fresh run')
    print(f'Measured cold {cold:.3f}s, unchanged resume {warm:.3f}s (this machine, public 3-loop fixture)')
