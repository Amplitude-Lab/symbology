"""Real FEC integration: cleanup must preserve results and reject damaged data."""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import stat
import subprocess

p=argparse.ArgumentParser()
p.add_argument('root',type=Path)
p.add_argument('--exe',default='bin/irr_fec.exe')
a=p.parse_args()
assert not a.root.exists(), 'Use a fresh test directory'
a.root.mkdir(parents=True)
env=dict(os.environ,THREADS='2')

def run(root,clean=False,ok=True):
    cmd=[a.exe,'4',str(root)] + (['--clean-block-cache'] if clean else [])
    result=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (a.root/(root.name+('_clean' if clean else '')+'.log')).write_bytes(result.stdout)
    assert (result.returncode==0)==ok, result.stdout[-2500:].decode(errors='replace')
    return result.stdout

def copy(source,target):
    assert not any(x.is_symlink() for x in source.rglob('*'))
    shutil.copytree(source,target)

def outputs(root):
    return {str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest()
            for x in root.rglob('*') if x.is_file() and (x.suffix=='.wxf' or x.name.endswith('.wxf.shape'))}

def unlink(file):
    for part in (file,*file.parents):
        assert not part.is_symlink()
        assert not (getattr(part.lstat(),'st_file_attributes',0)&getattr(stat,'FILE_ATTRIBUTE_REPARSE_POINT',1024))
    assert file.is_file() and file.stat().st_nlink==1
    file.unlink()

baseline=a.root/'keep'
run(baseline)
assert list(baseline.rglob('prime_*.bin'))
fresh=a.root/'fresh'
out=run(fresh,True)
assert outputs(fresh)==outputs(baseline)
assert not list(fresh.rglob('prime_*.bin'))
assert list(fresh.rglob('cleanup_done.tsv'))
assert b'BLOCK_CACHE_CLEAN' in out
run(fresh,True)
run(fresh)

resume=a.root/'resume'
copy(baseline,resume)
run(resume,True)
assert outputs(resume)==outputs(baseline)
assert not list(resume.rglob('prime_*.bin'))

# The cloud w01..w18 were cleaned by the old whole-stage tool, not this flag.
legacy=a.root/'legacy_cleaned'
copy(baseline,legacy)
cache_exe=Path(a.exe).with_name('irr_cache'+Path(a.exe).suffix)
result=subprocess.run([str(cache_exe),str(legacy),'--apply'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
assert result.returncode==0 and b'SKIP' not in result.stdout,result.stdout
assert not list(legacy.rglob('prime_*.bin'))
run(legacy,True)
result=subprocess.run([str(cache_exe),str(fresh),'--apply'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
assert result.returncode==0 and b'SKIP' not in result.stdout,result.stdout

# Recreate interruption after authorization and one unlink, but before the done receipt.
partial=a.root/'partial'
copy(fresh,partial)
first=next(baseline.rglob('prime_*.bin'))
relative=first.parent.relative_to(baseline)
for file in first.parent.glob('prime_*'):
    shutil.copy2(file,partial/relative/file.name)
unlink(partial/relative/'cleanup_done.tsv')
unlink(partial/relative/first.name)
run(partial,True)
assert outputs(partial)==outputs(baseline)
assert not list(partial.rglob('prime_*.bin'))
assert (partial/relative/'cleanup_done.tsv').exists()

for kind in ('damaged_output','damaged_cache','foreign','hardlink','missing_cache','linked','missing_marker','active'):
    fixture=a.root/kind
    copy(baseline,fixture)
    cache=next(fixture.rglob('prime_*.bin'))
    if kind=='damaged_output':
        file=next((fixture/'w01').glob('*.wxf'))
        with file.open('ab') as f:f.write(b'bad')
    elif kind=='damaged_cache':
        with cache.open('ab') as f:f.write(b'bad')
    elif kind=='foreign':
        (cache.parent/'foreign.bin').write_bytes(b'preserve')
    elif kind=='hardlink':
        os.link(cache,a.root/'hardlink_sentinel.bin')
    elif kind=='missing_cache':
        unlink(cache)
    elif kind=='missing_marker':
        unlink(next((fixture/'w01').glob('*.ready')))
    elif kind=='linked':
        outside=a.root/'outside';outside.mkdir()
        (outside/'sentinel').write_bytes(b'must survive')
        if os.name=='nt':
            subprocess.run(['cmd','/c','mklink','/J',str((cache.parent/'outside').resolve()),str(outside.resolve())],check=True,stdout=subprocess.DEVNULL)
        else:
            (cache.parent/'outside').symlink_to(outside.resolve(),target_is_directory=True)
    handle=None
    if kind=='active':
        if os.name=='nt':
            import ctypes
            create=ctypes.windll.kernel32.CreateFileW;create.restype=ctypes.c_void_p
            handle=create(str((fixture/'.active.lock').resolve()),0xC0000000,0,None,3,0x80,None)
            assert handle!=ctypes.c_void_p(-1).value
        else:
            import fcntl
            handle=(fixture/'.active.lock').open('r+');fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    before={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in cache.parent.glob('*.bin')}
    try:run(fixture,True,False)
    finally:
        if handle is not None:
            if os.name=='nt':ctypes.windll.kernel32.CloseHandle(ctypes.c_void_p(handle))
            else:handle.close()
    after={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in cache.parent.glob('*.bin')}
    assert before==after,kind
    if kind=='linked':assert (outside/'sentinel').read_bytes()==b'must survive'
print('FEC BLOCK CACHE INTEGRATION PASS: fresh/reuse/default/partial/corruption/foreign/hardlink/missing/link/lock')
