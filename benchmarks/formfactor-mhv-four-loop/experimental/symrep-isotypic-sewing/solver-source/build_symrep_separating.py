"""Isolated center-separation improvement; no repository headers are edited.

For a semisimple rational group algebra the center has dimension equal to
the number of conjugacy classes. A single central element can separate it
only when its squarefree regular characteristic polynomial has that degree.
Reject provably colliding elements before constructing their projectors.
"""
from pathlib import Path
import shutil,hashlib,json
B=Path(__file__).resolve().parents[2];R=B.parents[1]
S=B/'results/ff_four_loop/symrep_cached_build';O=B/'results/ff_four_loop/symrep_separating_build';O.mkdir(exist_ok=True)
for p in S.glob('*.hpp'):shutil.copy2(p,O/p.name)
s=(S/'symrep.hpp').read_text()
needle='long weight=attempt==0?long(c+1):long(rng%101)-50;'
assert s.count(needle)==1
s=s.replace(needle,'long weight=attempt==0?long(c+1):long(rng%1000003)-500001;')
needle='        auto power_start=std::chrono::steady_clock::now();'
assert s.count(needle)==1
s=s.replace(needle,'''        // A central element generates the whole rational center only if its
        // squarefree regular characteristic polynomial has full center degree.
        // Reject collisions before any projector polynomial is constructed.
        if(size_t(fmpq_poly_degree(total.p))!=g.classes.size()) {
            std::cout<<"PROFILE central_early_reject attempt="<<attempt<<std::endl;
            continue;
        }
        auto power_start=std::chrono::steady_clock::now();''')
(O/'symrep.hpp').write_text(s);shutil.copy2(S/'symrep.cpp',O/'symrep.cpp')
c=(S/'cache_control.cpp').read_text()
c=c.replace('require(a==b&&x==y,"central decomposition changed");','require(x.size()==y.size(),"central decomposition count changed"); for(const auto& v:x)require(std::find(y.begin(),y.end(),v)!=y.end(),"central projector set changed");')
(O/'separating_control.cpp').write_text(c)
manifest={'source_variant':'symrep_cached_build','changes':['reject central elements whose squarefree polynomial degree is smaller than the rational center dimension','deterministic random class weights widened from101 to1000003 values'],'unchanged_checks':['exact central-factor coprimality','each projected center dimension equals field degree','every central projector is idempotent','sum of projectors is1','rational representation and group-action certificates'],'source_hash':hashlib.sha256((S/'symrep.hpp').read_bytes()).hexdigest(),'new_hash':hashlib.sha256(s.encode()).hexdigest(),'scope':'Setup strategy only. May reorder central projectors after a changed separating element; controls compare exact projector sets. No physics constraints altered and no performance qualification claimed.'}
(O/'variant_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(O)
