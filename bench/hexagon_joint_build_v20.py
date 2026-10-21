#!/usr/bin/env python3
"""Build an isolated affine exact-suffix basis compression variant."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();assert not out.exists()
parent=R/'output_hexagon_v19_20261009/package'
shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
projection=out/'src/reduced_d3/irr_projection.hpp';code=projection.read_text()
start=code.index('inline equations collinear_tail_exact(context& c){')
end=code.index('// OPE propagation goes',start)
body=code[start:end]
marker='    c.record("collinear_rhs",rows.size(),1,0,rows.size()*8,start);start=c.now();'
assert body.count(marker)==1
extension='''
    auto selected=independent_collinear_suffix(c,la,rows,n);
    std::vector<size_t> selected_words;for(auto beta:selected)selected_words.push_back(c.tail_words[beta]);
    equations compressed;compressed.reserve(n*selected.size());
    for(size_t s=0;s<n;++s)for(auto beta:selected)compressed.push_back(std::move(rows[s*t+beta]));
    rows=std::move(compressed);t=selected.size();stride=n*t;
'''
body=body.replace(marker,marker+extension)
body=body.replace('c.exact_tail[l*q+c.tail_words[beta]]','c.exact_tail[l*q+selected_words[beta]]')
code=code[:start]+'} // namespace irr_physics\n#include "shared/hexagon_collinear_tail_compression.hpp"\nnamespace irr_physics {\n'+body+code[end:]
projection.write_text(code)
shutil.copy2(R/'bench/hexagon_collinear_tail_compression.hpp',out/'src/shared/hexagon_collinear_tail_compression.hpp')
driver=(parent/'tests/joint_physics_v19.cpp').read_text()
assert driver.count('joint-discovery-physical-affine-v19-exact-collinear-tail')==1
driver=driver.replace('joint-discovery-physical-affine-v19-exact-collinear-tail','joint-discovery-physical-affine-v20-affine-suffix-compression')
driver=driver.replace('exact-all-tail-slots; scalar independent validation','exact affine suffix row basis; all coefficient/RHS row-span certificates; scalar independent validation')
(out/'tests/joint_physics_v20.cpp').write_text(driver)
runtime=R/'output_hexagon_runtime_recovery_20261008/env'
compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
     '-I'+str(runtime/'include'),'tests/joint_physics_v20.cpp','-o','bin/joint_physics_v20',
     '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*')
        if q.is_file() and ('src' in q.parts or 'tests' in q.parts or 'bin' in q.parts)}
(out/'v20-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),command=cmd,hashes=hashes,
 change='Exact local compression of the affine collinear suffix operator, including RHS-only constraints. Scalar independent checks and physical/bonus inputs preserved.'),indent=2)+'\n')
print('BUILT V20',flush=True)
