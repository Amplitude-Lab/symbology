#!/usr/bin/env python3
"""Isolated universal central-polynomial power cache; repository headers untouched."""
from pathlib import Path
import shutil,hashlib,json
B=Path(__file__).resolve().parents[2];R=B.parents[1];O=B/'results/ff_four_loop/symrep_cached_build';O.mkdir(exist_ok=True)
for p in R.glob('*.hpp'):shutil.copy2(p,O/p.name)
s=(R/'symrep.hpp').read_text();needle='        std::vector<std::vector<Q>> result;degrees.clear();bool separated=true;'
repl='''        std::cout<<"PROFILE central_attempt="<<attempt<<" group_order="<<g.order()<<" factors="<<fs.size()<<" distinct_degree="<<fmpq_poly_degree(total.p)<<std::endl;
        auto power_start=std::chrono::steady_clock::now();
        std::vector<std::vector<Q>> central_powers(1,std::vector<Q>(g.order(),Q(0)));central_powers[0][0]=Q(1);
        for(slong k=1;k<fmpq_poly_degree(total.p);++k)central_powers.push_back(g.product(central_powers.back(),z));
        auto cached_evaluate=[&](const Poly& polynomial){std::vector<Q> value(g.order(),Q(0));
            for(slong k=0;k<=fmpq_poly_degree(polynomial.p);++k){Q c;fmpq_poly_get_coeff_fmpq(c.data(),polynomial.p,k);if(c==0)continue;
                for(size_t h=0;h<g.order();++h)if(central_powers.at(k)[h]!=0)value[h]+=c*central_powers[k][h];}
            return value;};
        std::cout<<"PROFILE central_power_cache_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-power_start).count()<<std::endl;
        std::vector<std::vector<Q>> result;degrees.clear();bool separated=true;'''
assert s.count(needle)==1;s=s.replace(needle,repl);assert s.count('auto coeff=evaluate(g,e,z);')==1;s=s.replace('auto coeff=evaluate(g,e,z);','auto coeff=cached_evaluate(e);')
s=s.replace('auto centers=central_idempotents(g,degrees);','auto centers=central_idempotents(g,degrees);std::cout<<"PROFILE central_idempotents_complete count="<<centers.size()<<std::endl;')
(O/'symrep.hpp').write_text(s)
p=(R/'symrep.cpp').read_text().replace(';auto ms=saved?load_models', ';std::cout<<"PROFILE group_closed order="<<group.order()<<" classes="<<group.classes.size()<<" elapsed="<<seconds(start)<<std::endl;auto ms=saved?load_models')
(O/'symrep.cpp').write_text(p)
(O/'variant_manifest.json').write_text(json.dumps({'change':'Cache shared powers of a separating central element once per attempt; evaluate all projector polynomials linearly in these cached powers. Retain exact central separation,idempotent,completeness and subsequent representation checks.','original_symrep_hpp_sha256':hashlib.sha256((R/'symrep.hpp').read_bytes()).hexdigest(),'modified_symrep_hpp_sha256':hashlib.sha256(s.encode()).hexdigest(),'profile_scope':'Representation setup only; scoped comparison and low-group independent equality controls required. No V16 qualification claim.'},indent=2)+'\n')
print(O)
