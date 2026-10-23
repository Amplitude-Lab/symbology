from pathlib import Path
import shutil,json,hashlib
B=Path(__file__).resolve().parents[2];R=B.parents[1];D=B/'results/ff_four_loop/sew_masked_build';D.mkdir(exist_ok=True)
for p in R.glob('*.hpp'):shutil.copy2(p,D/p.name)
p=D/'product_invariance.hpp';s=p.read_text().replace('namespace staged_kernel {','namespace staged_kernel {\ninline std::vector<unsigned char> ff_active_output;\ninline size_t ff_diagonal_generators=5;',1)
s=s.replace('size_t diagonal=l*right_dimension+r;','size_t diagonal=l*right_dimension+r;\n            if(g>=ff_diagonal_generators&&!ff_active_output.empty()&&!ff_active_output[diagonal])continue;')
s=s.replace('auto build=[&](size_t c){size_t id=g*ncol()+c;', 'auto build=[&](size_t c){if(g>=ff_diagonal_generators&&!ff_active_output.empty()&&!ff_active_output[c])return MV{};size_t id=g*ncol()+c;')
p.write_text(s)
p=D/'tensor_kernel.hpp';s=p.read_text();mark='class ConstraintRows {';s=s.replace(mark,'inline std::vector<unsigned> ff_old_grades,ff_local_column_grades;\n'+mark,1).replace('    bool ordered;','    bool ordered;\n    std::vector<unsigned> ff_row_grades;',1)
needle='        require(local.ncol==letters*next,"local condition dimensions mismatch");'
s=s.replace(needle,needle+'''\n        if(!ff_old_grades.empty()){require(ff_old_grades.size()==older&&ff_local_column_grades.size()==local.ncol,"grade filter dimensions");
            for(const auto& row:local.rows){require(row.nnz()>0,"zero local grade row");unsigned g=ff_local_column_grades[row(0)];for(auto c:row.index_span())require(ff_local_column_grades[c]==g,"local equation mixes Galois grades");ff_row_grades.push_back(g);}}''',1)
s=s.replace('for(size_t old=0;old<older;++old)for(const auto& equation:local.rows){\n            if(limit!=SIZE_MAX)', 'for(size_t old=0;old<older;++old)for(const auto& equation:local.rows){\n            if(!ff_row_grades.empty()&&ff_old_grades[old]!=ff_row_grades[&equation-local.rows.data()])continue;\n            if(limit!=SIZE_MAX)',1)
s=s.replace('old=id/cm.size();if((seen&&', 'old=id/cm.size();if(!ff_row_grades.empty()&&ff_old_grades[old]!=ff_row_grades[id%cm.size()])return ModVec{};if((seen&&',1)
s=s.replace('for(auto [lj,a]:cm[id%cm.size()])cost+=fm[old*letters+lj/next].nnz();','if(ff_row_grades.empty()||ff_old_grades[old]==ff_row_grades[id%cm.size()])for(auto [lj,a]:cm[id%cm.size()])cost+=fm[old*letters+lj/next].nnz();',1)
p.write_text(s)
s=(B/'code/ff_four_loop/sew_selective_presolve.cpp').read_text();s=s.replace('require(argc>=5&&argc<=6,','require(argc>=5,').replace('require(mode=="full"||mode=="galois","bad mode");rref_option_t opt;opt->pool.reset(1);','require(mode=="full"||mode=="galois"||mode=="masked","bad mode");size_t threads=1;std::string referencePath;for(int i=5;i<argc;++i){if(std::string(argv[i])=="--threads"){require(i+1<argc,"missing threads");threads=std::stoul(argv[++i]);}else referencePath=argv[i];}rref_option_t opt;opt->pool.reset(threads);')
needle=' staged_kernel::ProductInvarianceRows invariants(std::move(left),std::move(right)),galois(std::move(short_left),std::move(short_right));'
repl=''' if(mode=="masked"){
 auto grades=[](const std::vector<Mat>&a){std::vector<unsigned>g(a[0].nrow);for(size_t k=0;k<5;++k)for(size_t i=0;i<g.size();++i){require(a[k][i].nnz()==1&&a[k][i](0)==index(i)&&(a[k][i][0]==1||a[k][i][0]==-1),"non-diagonal Galois action");if(a[k][i][0]==-1)g[i]|=1u<<k;}return g;};
 auto lg=grades(left),rg=grades(right);std::vector<Mat>old,letters;std::array<std::string,5>gn{"galois1a","galois1b","galois2a","galois2b","galois3"};for(size_t k=1;k<=5;++k){old.push_back(read_matrix(out/("galois"+std::to_string(k)+"_w"+std::to_string(w-1)+".wxf")));letters.push_back(read_matrix(out/"data"/(gn[k-1]+".wxf")));}
 tensor_kernel::ff_old_grades=grades(old);auto ag=grades(letters);tensor_kernel::ff_local_column_grades.resize(88*nl);for(size_t a=0;a<88;++a)for(size_t b=0;b<nl;++b)tensor_kernel::ff_local_column_grades[a*nl+b]=ag[a]^rg[b];
 for(size_t i=0;i<nf;++i)for(size_t p=f.rowptr()[i];p<f.rowptr()[i+1];++p){auto c=f.index(p);require(lg[i]==(tensor_kernel::ff_old_grades[c[0]]^ag[c[1]]),"prefix recurrence grade mismatch");}
 for(size_t h=5;h<7;++h){std::array<int,32>trans;trans.fill(-1);size_t terms=0;auto check=[&](const Mat&a,const std::vector<unsigned>&g){for(size_t i=0;i<a.nrow;++i)for(auto[j,x]:a[i]){auto&v=trans[g[i]];if(v<0)v=g[j];require(v==int(g[j]),"D4 action not grade homogeneous");++terms;}};check(left[h],lg);check(right[h],rg);std::set<int>images;for(int v:trans)if(v>=0)require(images.insert(v).second,"D4 grade map not injective");std::cout<<"GRADE_NORMALITY generator="<<h<<" exact_support_terms="<<terms<<" injective=1"<<std::endl;}
 staged_kernel::ff_active_output.resize(nf*nl);size_t count=0;for(size_t a=0;a<nf;++a)for(size_t b=0;b<nl;++b)count+=(staged_kernel::ff_active_output[a*nl+b]=lg[a]==rg[b]);std::cout<<"GRADE_MASK active_product="<<count<<" of="<<nf*nl<<std::endl;
 }
'''+needle
assert s.count(needle)==1;s=s.replace(needle,repl).replace('if(mode=="full")k=','if(mode=="full"||mode=="masked")k=').replace('if(argc==6){auto reference=read_matrix(argv[5]);','if(!referencePath.empty()){auto reference=read_matrix(referencePath);')
(D/'sew_masked.cpp').write_text(s)
(D/'manifest.json').write_text(json.dumps({'scope':'Exact static totalGaloiseven output masks only. Five diagonal equations remain intact. D4 grade permutation is checked on every action entry and injective on all occupied grades; every prefix coefficient and local RREF row is grade-checked. Full final residual operators are unchanged. Global repository headers unchanged.','files_sha256':{n:hashlib.sha256((D/n).read_bytes()).hexdigest()for n in ['sew_masked.cpp','product_invariance.hpp','tensor_kernel.hpp']}},indent=2)+'\n');print(D)
