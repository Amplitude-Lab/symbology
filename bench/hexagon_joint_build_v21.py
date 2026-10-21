#!/usr/bin/env python3
"""Build exact head-support-aware affine suffix compression."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();assert not out.exists()
parent=R/'output_hexagon_v19_20261009/package-v20'
shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'));(out/'bin').mkdir()
header=out/'src/shared/hexagon_collinear_tail_compression.hpp';code=header.read_text()
code=code.replace('const equations& raw,size_t samples){','const equations& raw,size_t samples,const std::vector<unsigned char>* active=nullptr){')
code=code.replace('const size_t q=size_t(1)<<(2*c.k);const auto& F=c.F;',
                  'const size_t q=size_t(1)<<(2*c.k);const auto& F=c.F;require(!active||active->size()==last.dim,"Head support shape");')
old='source[beta][i]=c.exact_tail[i*q+c.tail_words[beta]];';assert code.count(old)==1
code=code.replace(old,'source[beta][i]=(!active||(*active)[i])?c.exact_tail[i*q+c.tail_words[beta]]:0;')
header.write_text(code)
projection=out/'src/reduced_d3/irr_projection.hpp';code=projection.read_text()
old='    auto selected=independent_collinear_suffix(c,la,rows,n);';assert code.count(old)==1
new='''    // A tail component can be omitted from the local operator only if
    // EVERY corresponding head copy at EVERY current sample is exactly zero.
    // Recompute this support at each field/batch; no parity assumption is used.
    std::vector<unsigned char> active(la.dim);size_t active_count=0;
    for(size_t b=0;b<6;++b)for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){
        bool used=false;
        for(size_t f=0;f<fa.m[b]&&!used;++f)for(size_t s=0;s<n;++s)
            if(head[(fa.off[b]+mu*fa.m[b]+f)*n+s]){used=true;break;}
        if(used)for(size_t l=0;l<la.m[b];++l){active[la.off[b]+mu*la.m[b]+l]=1;++active_count;}
    }
    std::cout<<"COL_HEAD_SUPPORT active_tail_components="<<active_count<<" total_tail_components="<<la.dim
        <<" every_skipped_head_copy_checked_zero=1"<<std::endl;
    auto selected=independent_collinear_suffix(c,la,rows,n,&active);'''
projection.write_text(code.replace(old,new))
driver=(parent/'tests/joint_physics_v20.cpp').read_text()
driver=driver.replace('joint-discovery-physical-affine-v20-affine-suffix-compression',
                      'joint-discovery-physical-affine-v21-head-supported-affine-suffix')
driver=driver.replace('exact affine suffix row basis; all coefficient/RHS row-span certificates; scalar independent validation',
                      'exact head-supported affine suffix row basis; all coefficient/RHS row-span certificates; scalar independent validation')
(out/'tests/joint_physics_v21.cpp').write_text(driver)
runtime=R/'output_hexagon_runtime_recovery_20261008/env';compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
 '-I'+str(runtime/'include'),'tests/joint_physics_v21.cpp','-o','bin/joint_physics_v21',
 '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*')
 if q.is_file() and ('src' in q.parts or 'tests' in q.parts or 'bin' in q.parts)}
(out/'v21-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),command=cmd,hashes=hashes,
 change='Exact batch/field-specific head zero-support masks before affine suffix row compression. All RHS columns retained. No extra physics relation or omitted final verification.'),indent=2)+'\n')
print('BUILT V21',flush=True)
