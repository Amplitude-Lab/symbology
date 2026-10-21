#!/usr/bin/env python3
"""Build a separate weight-dependent OPE-head storage candidate."""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args(); out=a.output.resolve(); assert not out.exists()
parent=R/'output_hexagon_v19_20261009/package-v21'
shutil.copytree(parent,out,ignore=shutil.ignore_patterns('bin','.git','output','log'))
(out/'bin').mkdir()
g=out/'src/nmhv_ope_germ.hpp'; code=g.read_text()
start=code.index('    void integrate('); end=code.rindex('\n};')
method=code[start:end]
method=method.replace('void integrate(const ulong* previous,ulong* next,size_t letter,size_t w)const {',
'''void integrate_from(const nmhv_ope_germ_t& previous_layout,const ulong* previous,
        ulong* next,size_t letter,size_t w)const {
        if(w<1||w>weight||previous_layout.weight<w-1||previous_layout.max_t!=max_t||
           previous_layout.grade!=grade||previous_layout.F.mod.n!=F.mod.n)
            throw std::runtime_error("Incompatible compact OPE integration layouts");''')
method=method.replace('previous[index(', 'previous[previous_layout.index(')
code=code[:end]+'''
    // At weight w, a+b<=w. Preserve the complete T/S jet and every kernel;
    // only logarithmic slots that are identically zero are absent.
    nmhv_ope_germ_t with_log_weight(size_t w)const {
        if(w>weight)throw std::runtime_error("Cannot enlarge compact OPE source layout");
        auto copy=*this;copy.weight=w;copy.logs=w+1;
        copy.stride=size_t(max_t+1)*size_t(grade+1)*copy.logs*(copy.logs+1)/2;
        return copy;
    }
'''+method+code[end:]
g.write_text(code)
q=out/'src/reduced_d3/irr_projection.hpp'; code=q.read_text()
start=code.index('template<class Transfer>'); end=code.index('inline std::vector<ulong> terminal(',start)
step=code[start:end]
old='const std::array<mat,6>* modular=nullptr){'
assert step.count(old)==1
step=step.replace(old,'const std::array<mat,6>* modular=nullptr,size_t input_stride=0){\n    if(!input_stride)input_stride=stride;\n    require(previous.size()==old.dim*input_stride,"Recurrence input stride mismatch");')
assert step.count('previous.data()+r*stride')==2
step=step.replace('previous.data()+r*stride','previous.data()+r*input_stride')
helper='''
inline std::vector<ulong> compact_ope_head(context& c,const nmhv_ope_germ_t& germ){
    auto previous_layout=germ.with_log_weight(0);
    std::vector<ulong> head(previous_layout.stride);
    head[previous_layout.index(0,0,0,0)]=1;
    for(size_t w=1;w<=c.h;++w){
        auto current_layout=germ.with_log_weight(w);auto K=c.modular_recurrence(w,true);
        axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
        head=step(c,irr::matrices{},old,next,head,current_layout.stride,
            [&](const ulong* a,ulong* b,size_t l){current_layout.integrate_from(previous_layout,a,b,l,w);},
            "ope_head_"+std::to_string(w),K.get(),previous_layout.stride);
        previous_layout=std::move(current_layout);
    }
    // Restore the existing external/cache layout. Descending copies permit
    // in-place expansion; omitted logarithmic coordinates are explicitly zero.
    size_t rows=head.size()/previous_layout.stride;
    size_t compact_size=head.size();head.resize(rows*germ.stride);
    for(size_t r=rows;r-->0;)for(int t=germ.max_t;t>=0;--t)
        for(int z=germ.grade-t;z>=-t;--z)for(size_t a=c.h+1;a-->0;)
            for(size_t b=c.h-a+1;b-->0;)
                head[r*germ.stride+germ.index(t,z,a,b)]=
                    head[r*previous_layout.stride+previous_layout.index(t,z,a,b)];
    for(size_t r=0;r<rows;++r)for(int t=0;t<=germ.max_t;++t)
        for(int z=-t;z<=germ.grade-t;++z)for(size_t a=0;a<=germ.weight;++a)
            for(size_t b=0;b<=germ.weight-a;++b)if(a+b>c.h)
                head[r*germ.stride+germ.index(t,z,a,b)]=0;
    std::cout<<"OPE_WEIGHT_STORAGE terminal_weight="<<c.h<<" compact_bytes="<<compact_size*sizeof(ulong)
        <<" restored_bytes="<<head.size()*sizeof(ulong)<<" omitted_slots_exact_zero=1"<<std::endl;
    return head;
}

'''
code=code[:start]+step+helper+code[end:]
old='''        std::vector<ulong> calculated(stride);calculated[germ.index(0,0,0,0)]=1;
        for(size_t w=1;w<=c.h;++w){auto K=c.modular_recurrence(w,true);axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
            irr::matrices unused;calculated=step(c,unused,old,next,calculated,stride,[&](const ulong* a,ulong* b,size_t l){germ.integrate(a,b,l,w);},"ope_head_"+std::to_string(w),K.get());}'''
assert code.count(old)==1
new='''        std::vector<ulong> calculated;
        if constexpr(std::is_same_v<Germ,nmhv_ope_germ_t>)calculated=compact_ope_head(c,germ);
        else {
'''+old+'''
        head_storage=std::make_shared<const std::vector<ulong>>(std::move(calculated));
        }
'''
# Keep the generic legacy MHV branch intact without shadowing calculated.
new=new.replace('        std::vector<ulong> calculated(stride);','        calculated.resize(stride);')
new=new.replace('        head_storage=std::make_shared<const std::vector<ulong>>(std::move(calculated));\n        }','        }')
q.write_text(code.replace(old,new))
q=out/'src/shared/direct_ope_rows.hpp';code=q.read_text()
start=code.index('    std::vector<ulong> head(germ.stride);');end=code.index('    auto result=',start)
code=code[:start]+'    auto head=compact_ope_head(c,germ);\n'+code[end:];q.write_text(code)
driver=(parent/'tests/joint_physics_v21.cpp').read_text().replace(
 'joint-discovery-physical-affine-v21-head-supported-affine-suffix',
 'joint-discovery-physical-affine-v22-weight-dependent-ope-head')
(out/'tests/joint_physics_v22.cpp').write_text(driver)
shutil.copy2(R/'bench/hexagon_ope_weight_control.cpp',out/'tests/ope_weight_control.cpp')
runtime=R/'output_hexagon_runtime_recovery_20261008/env'
compiler=R/'output_hexagon_candidate_export_20261009/build-env/bin/x86_64-conda-linux-gnu-c++'
commands=[]
for target in ['ope_weight_control','joint_physics_v22']:
    cmd=[str(compiler),'-O3','-std=c++20','-Isrc','-Isrc/shared','-Isrc/reduced_d3','-Isrc/SparseRREF',
      '-I'+str(runtime/'include'),'tests/'+target+'.cpp','-o','bin/'+target,
      '-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lflint','-lmpfr','-lgmp','-lmimalloc','-l:libtbb.so.12','-pthread']
    commands.append(cmd);subprocess.run(cmd,cwd=out,check=True)
hashes={str(q.relative_to(out)):hashlib.sha256(q.read_bytes()).hexdigest() for q in out.rglob('*')
 if q.is_file() and ('src' in q.parts or 'tests' in q.parts or 'bin' in q.parts)}
(out/'v22-build-receipt.json').write_text(json.dumps(dict(parent=str(parent),commands=commands,hashes=hashes,
 change='Weight-dependent logarithmic triangle during OPE head recurrence; complete existing cache/projection layout restored; no physics/jet-order change.'),indent=2)+'\n')
print('BUILT V22',flush=True)
