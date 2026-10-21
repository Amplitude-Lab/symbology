"""Check SHA256 padding/read boundaries against Python's independent hashlib."""
import hashlib
from pathlib import Path
import subprocess
import sys

exe,root=sys.argv[1],Path(sys.argv[2])
root.mkdir(parents=True,exist_ok=False)
for size in (0,1,55,56,63,64,65,65535,65536,65537,4194304):
    data=bytes(range(256))*(size//256)+bytes(range(size%256))
    file=root/str(size)
    file.write_bytes(data)
    actual=subprocess.check_output([exe,str(file)]).decode().strip()
    assert actual==hashlib.sha256(data).hexdigest(),size
print('STREAM SHA256 PASS: empty/padding/block/read-buffer/4MiB boundaries')
