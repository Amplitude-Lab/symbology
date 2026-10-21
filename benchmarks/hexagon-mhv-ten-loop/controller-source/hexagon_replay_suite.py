#!/usr/bin/env python3
"""Serialize validation and allocator comparisons; no overlapping large jobs."""
import json,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009'
while True:
 x=json.loads((V/'progress.json').read_text())
 if x['status']=='failed':raise RuntimeError('Physical validation failed; no performance selection')
 if x['status']=='passed':break
 time.sleep(5)
for mode in ['mapped','heap']:
 subprocess.run(['python3',str(R/'bench/hexagon_ten_loop_benchmark.py'),'replay','--output',str(V/('replay-'+mode)),'--storage',mode],check=True,cwd=R)
rows={mode:json.loads((V/('replay-'+mode)/'manifest.json').read_text()) for mode in ['mapped','heap']}
(V/'replay-comparison.json').write_text(json.dumps({'schema':'hexagon-post-OPE-backend-comparison-v1','passed':True,'runs':rows,'scope':'One serial warm-basis replay per backend, fixed O2 binary and checkpoint. No end-to-end speedup or global optimum claimed.'},indent=2)+'\n')
