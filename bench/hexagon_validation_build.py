#!/usr/bin/env python3
"""Build separate validation tools; original production source remains frozen."""
import json,subprocess,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];V=R/'output_hexagon_validation_20261009';P=R/'output_hexagon_17p3_20261007/joint-v17-schedule-package'
base=json.loads((R/'output_hexagon_candidate_export_20261009/build-command-gcc14.json').read_text())
for name in ['hexagon_physical_holdouts','hexagon_joint_certificate']:
 cmd=[str(x).replace(str(R/'output_hexagon_16p4_20261006/joint-v8-package'),str(P)).replace('bench/hexagon_candidate_qq_export.cpp','bench/'+name+'.cpp').replace(str(R/'output_hexagon_candidate_export_20261009/candidate_qq_export'),str(V/name)) for x in base]
 (V/(name+'.build.json')).write_text(json.dumps(cmd,indent=2)+'\n')
 subprocess.run(cmd,check=True)
 print(name,'built',flush=True)
