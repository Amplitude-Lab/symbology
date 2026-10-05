#include "symrep_bootstrap.hpp"
#include <iostream>
using namespace symrep;
void test_incremental_wxf() {
    auto root=fs::temp_directory_path()/("symrep-wxf-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));fresh(root);
    Mat m(3,8);m[0].push_back(0,Q(1,3));
    m[0].push_back(7,Q(int_t("123456789012345678901234567890123456789"),int_t("100000000000000000000000000003")));
    m[2].push_back(1,Q(-7));auto t=tensor(m,{3,2,4});
    write(root/"tensor.wxf",t);MatrixTensorView view(m,{3,2,4});write_tensor_view(root/"view.wxf",view);
    auto bytes=file_to_ustr(root/"tensor.wxf");
    require(bytes==file_to_ustr(root/"view.wxf"),"direct output differs from tensor serialization");
    auto ordinary=read(root/"view.wxf");require(equal(flatten(ordinary),m),"full reader disagrees with incremental readback");
    // The visitor must reject valid-but-different values as well as malformed
    // syntax. It must not replace exact semantic checking with a checksum.
    auto rejects=[&](const std::vector<uint8_t>& bad){std::ofstream out(root/"bad.wxf",std::ios::binary);out.write(reinterpret_cast<const char*>(bad.data()),bad.size());out.close();
        bool caught=false;try{verify_written_tensor(root/"bad.wxf",t);}catch(const std::exception&){caught=true;}require(caught,"incremental verifier accepted corrupted output");};
    auto bad=bytes;bad.pop_back();rejects(bad);bad=bytes;bad.back()^=1;rejects(bad);bad=bytes;bad.push_back(0);rejects(bad);
    WXF_PARSER::Parser parser(bytes);parser.parse();
    bad=bytes;bad[parser.tokens[11].data-bytes.data()]=0;rejects(bad);
    for(size_t i=14;i+3<parser.tokens.size();++i)if(parser.tokens[i].type==WXF_PARSER::WXF_HEAD::func){
        bad=bytes;const auto& denominator=parser.tokens[i+3];
        if(denominator.type==WXF_PARSER::WXF_HEAD::i8){bad[denominator.data-bytes.data()]=0;rejects(bad);break;}}
    for(size_t rows:{0,3}){Mat empty(rows,0);MatrixTensorView v(empty,{rows,0,4});auto path=root/("empty"+std::to_string(rows)+".wxf");
        write_tensor_view(path,v);auto check=read(path);require(check.dims()==v.dims()&&!check.nnz(),"empty view roundtrip");}
    fs::remove_all(root);std::cout<<"PASS incremental WXF exact values, large rationals, empty tensors and corruption rejection"<<std::endl;
}
void test_reconstruction_paths(rref_option_t opt) {
    // I | a | b has known rank and nullity. Its size exercises parallel first
    // reconstruction, and the large rational needs multiple CRT primes.
    const size_t n=22000;Mat a(n,n+2);Q big(int_t("1208925819614629174706177"),int_t("1180591620717411303425"));
    for(size_t i=0;i<n;++i){a[i].push_back(index(i),Q(1));a[i].push_back(index(n),i%2?big:-big);a[i].push_back(index(n+1),Q(5,7));}
    auto k=kernel(a,opt);require(k.nrow==2&&!mul(a,transpose(k),&opt->pool).nnz(),"parallel multi-prime reconstruction failed");
    ulong prime=n_nextprime(1ULL<<60,0),bound=n_sqrt((prime-1)/2);
    for(long delta:{-1L,0L,1L})for(long sign:{-1L,1L}){Mat b(1,2);b[0].push_back(0,Q(1));b[0].push_back(1,Q(sign*(long(bound)+delta)));
        auto basis=kernel(b,opt);require(basis.nrow==1&&!mul(b,transpose(basis)).nnz(),"integer reconstruction bound failed");}
    std::cout<<"PASS parallel CRT and signed integer reconstruction boundaries"<<std::endl;
}
void test_direct_kernel(rref_option_t opt) {
    for(size_t trial=0;trial<24;++trial){size_t nr=trial%9,nc=trial%13;Mat a(nr,nc);
        for(size_t i=0;i<nr;++i)for(size_t j=0;j<nc;++j){long numerator=long((i*17+j*11+trial*7)%19)-9;
            if((i+j+trial)%3&&numerator)a[i].push_back(index(j),Q(numerator,long((j+trial)%5)+1));}
        if(nr>2)a[nr-1]=a[0]; // dependent rows, zero columns and rectangular cases
        auto actual=kernel(a,opt);require(actual.ncol==nc,"direct kernel shape mismatch");
        require(!mul(a,transpose(actual),&opt->pool).nnz(),"direct kernel residual");
        if(!a.nnz()){require(equal(actual,identity(nc)),"free direct kernel");continue;}
        auto reduced=a;std::erase_if(reduced.rows,[](const Vec& row){return !row.nnz();});reduced.nrow=reduced.rows.size();reduced.sort_rows_by_nnz();
        auto pivots=sparse_mat_rref_reconstruct(reduced,opt);size_t rank=0;for(const auto& batch:pivots)rank+=batch.size();
        auto expected=rank==nc?Mat(0,nc):sparse_mat_rref_kernel(reduced,pivots,QQ,opt).transpose();
        for(auto& row:expected.rows)normalize(row);require(equal(actual,expected),"direct kernel differs from SparseRREF column kernel");
    }
    std::cout<<"PASS direct row kernel versus SparseRREF extraction"<<std::endl;
}
Mat permutation(const std::vector<size_t>& p){Mat m(p.size(),p.size());for(size_t i=0;i<p.size();++i)m[i].push_back(index(p[i]),Q(1));return m;}
void test_carrier_coordinates(rref_option_t opt){
    std::vector<Mat> older{permutation({1,0})},letters{permutation({1,2,0})};
    field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));
    for(bool backward:{false,true})for(size_t mode=0;mode<3;++mode){
        Mat basis(mode==2?0:6,6);
        for(size_t i=0;i<basis.nrow;++i)for(size_t j=0;j<6;++j)
            if(mode==1||i==j)basis[i].push_back(index(j),Q((i==j?2:1)*(i+1),i+2));
        auto ambient=backward?kron(letters[0],older[0]):kron(older[0],letters[0]);
        auto expected=induced(basis,{ambient},opt),actual=carrier_actions(basis,older,letters,backward,opt);
        require(equal(actual[0],expected[0]),"selected-column exact carrier action differs from full tensor action");
        require(equal(mul(actual[0],basis),mul(basis,ambient)),"recursive exact carrier intertwiner failed");
        auto modular=carrier_actions_mod(basis,{modular_matrix(older[0],fp)},{modular_matrix(letters[0],fp)},backward,fp,opt);
        auto reduction=modular_matrix(expected[0],fp);
        require(modular[0].nrow==reduction.nrow&&modular[0].ncol==reduction.ncol,"modular carrier shape mismatch");
        for(size_t i=0;i<reduction.nrow;++i){auto row=modular[0][i];sparse_vec_sub_mul(row,reduction[i],ulong(1),fp);require(!row.nnz(),"modular carrier differs from exact action");}
    }
    std::cout<<"PASS selected tensor columns, private/general charts, both directions and zero-dimensional actions"<<std::endl;
}
std::vector<Mat> dihedral(size_t n,bool reflection=true){std::vector<size_t> c(n),f(n);for(size_t i=0;i<n;++i){c[i]=(i+1)%n;f[i]=(n-i)%n;}auto out=std::vector<Mat>{permutation(c)};if(reflection)out.push_back(permutation(f));return out;}
void test_orbits(const Group& g,const std::vector<Model>& ms,const std::vector<Mat>& gs,rref_option_t opt){
    Layout copies;copies.copies.assign(ms.size(),3);size_t n=copies.dimension(ms);Vec vector;
    for(size_t j=0;j<n;++j)if(j%3)vector.push_back(index(j),Q(long(j%7)-3,long(j%5)+1));normalize(vector);
    for(size_t h=0;h<g.order();++h){auto actual=apply_irrep_action(vector,copies,ms,h),expected=mul(vector,action(copies,ms,h));
        require(equal(rows({actual},n),rows({expected},n)),"implicit irrep action differs from explicit block matrix");}
    if(!fixed_corner_plans(g,ms))return;
    auto recipe=select_fixed_orbits(g,ms,gs,opt);
    auto frame=materialize_orbits(recipe,g,ms,gs,opt);
    verify_adaptation(frame,g,ms,gs);
    require(row_basis(frame.basis).nrow==gs[0].nrow,"compact orbit certificate lost QQ rank");
}
void test_recurrence(const Group& g,const std::vector<Model>& ms,const Layout& letters,rref_option_t opt) {
    TensorAdapter adapter(g,ms,opt);ProductCopies seed_channels;auto cg=adapter.product(letters,letters,&seed_channels);size_t n=letters.dimension(ms);
    Mat wedge(n*(n-1)/2,n*n);size_t r=0;for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j){wedge[r].push_back(index(i*n+j),Q(1));wedge[r++].push_back(index(j*n+i),Q(-1));}
    auto ag=generators(g,ms,letters);std::vector<Mat> pair;for(auto& a:ag)pair.push_back(kron(a,a));
    auto eq=adapt(g,ms,equation_actions(wedge,pair,opt),opt);
    auto condition=condition_tensor(mul(transpose(sparse_mat_inverse(eq.basis,QQ,opt)),wedge),n);
    // Use a nontrivial End(U) coefficient, including the quaternionic case.
    Mat commute(cg.basis.nrow,cg.basis.nrow);auto offsets=cg.layout.offsets(ms);
    for(size_t s=0;s<ms.size();++s){auto corner=row_basis(g.evaluate(ms[s].primitive,ms[s].actions));
        for(size_t copy=0;copy<cg.layout.copies[s];++copy)for(size_t i=0;i<ms[s].dimension;++i)
            for(auto [j,x]:mul(corner[corner.nrow-1],ms[s].actions[ms[s].orbit[i]]))commute[offsets[s]+copy*ms[s].dimension+i].push_back(index(offsets[s]+copy*ms[s].dimension+j),Q(copy+1)*x);}
    auto expansion=mul(commute,cg.basis);
    for(bool backward:{false,true}){
        Prepared p{letters,eq.layout,cg.layout,letters,condition,tensor(expansion,{expansion.nrow,n,n}),backward};
        EquivariantConditions reduced(p,g,ms,adapter,opt);
        // Exercise the saved reduced-data path too, including number-field and
        // quaternionic coefficients; the forward path still compresses on demand.
        if(backward){p.reduced_condition=reduced.coefficients;p.condition_pair_basis=reduced.pair.basis;}
        RecurrenceTemplates cache(g,ms,p,adapter,opt);ProductCopies source_channels,target_channels;
        auto domain=backward?adapter.product(letters,cg.layout,&source_channels):adapter.product(cg.layout,letters,&source_channels);
        auto target=backward?adapter.product(eq.layout,letters,&target_channels):adapter.product(letters,eq.layout,&target_channels);
        auto chart=transpose(block_inverse(target.basis,opt));auto so=domain.layout.offsets(ms),to=target.layout.offsets(ms);
        size_t count=0;
        for(int pass=0;pass<2;++pass){auto seed=p.expansion;if(pass){auto scaled=expansion;for(auto& row:scaled.rows)sparse_vec_rescale(row,Q(2,3),QQ);seed=tensor(scaled,{scaled.nrow,n,n});
                auto compact=commute;for(auto& row:compact.rows)sparse_vec_rescale(row,Q(2,3),QQ);cache.remember({compact,cg.layout},cg.layout,seed_channels);}
            auto blocks=cache.assemble(seed,letters,cg.layout,source_channels,target_channels);auto full=assemble(seed,condition,backward);
            std::atomic<size_t> consumed=0;
            cache.assemble(seed,letters,cg.layout,source_channels,target_channels,[&](size_t s,Mat&& block,double){
                require(equal(blocks[s],block),"streamed recurrence differs from materialized blocks");++consumed;});
            require(consumed==ms.size(),"streamed recurrence missed a sector");
            // A consumer failure must propagate only after all workers finish,
            // so they cannot keep using the assembler's destroyed local maps.
            consumed=0;bool failed=false;
            try{cache.assemble(seed,letters,cg.layout,source_channels,target_channels,[&](size_t s,Mat&&,double){
                ++consumed;if(s==0)throw std::runtime_error("deliberate sector failure");});}
            catch(const std::runtime_error& error){failed=std::string(error.what())=="deliberate sector failure";}
            require(failed&&consumed==ms.size(),"streamed recurrence did not drain tasks on failure");
            for(size_t s=0;s<ms.size();++s){auto local=row_basis(g.evaluate(ms[s].primitive,ms[s].actions));Mat q(domain.layout.copies[s]*local.nrow,domain.basis.nrow),left(target.layout.copies[s]*local.nrow,target.basis.nrow);
                for(size_t c=0;c<domain.layout.copies[s];++c)for(size_t j=0;j<local.nrow;++j)for(auto [k,x]:local[j])q[c*local.nrow+j].push_back(index(so[s]+c*ms[s].dimension+k),x);
                for(size_t c=0;c<target.layout.copies[s];++c)for(size_t j=0;j<local.nrow;++j)left[c*local.nrow+j]=chart[to[s]+c*ms[s].dimension+local[j](0)];
                auto expected=mul(left,mul(full,transpose(mul(q,domain.basis))));require(equal(blocks[s],expected),"template differs from complete tensor contraction");}
            if(pass)require(cache.templates_built==count,"template cache was rebuilt for new coefficients");count=cache.templates_built;
        }
    }
}
void test(const std::string& name,std::vector<Mat> gens,size_t order,rref_option_t opt){
    std::vector<std::string> labels;for(size_t i=0;i<gens.size();++i)labels.push_back("g"+std::to_string(i));
    Group g(labels,gens);require(g.order()==order,"wrong group order in "+name);auto ms=models(g,opt);
    auto a=adapt(g,ms,gens,opt);verify_adaptation(a,g,ms,gens);require(row_basis(a.basis).nrow==a.basis.nrow,"adaptation lost rank");
    auto frame=carrier_frame(g,ms,gens,opt);verify_adaptation(frame,g,ms,gens);
    require(row_basis(frame.basis).nrow==frame.basis.nrow,"carrier frame lost rank");
    test_orbits(g,ms,gens,opt);
    test_orbits(g,ms,std::vector<Mat>(gens.size(),Mat(0,0)),opt);
    // Recover the induced action from a sparse kernel chart and compare it
    // with the complete ambient product, independently of frame selection.
    size_t n=gens[0].nrow;Mat wedge(n*(n-1)/2,n*n);size_t r=0;
    for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j){wedge[r].push_back(index(i*n+j),Q(1));wedge[r++].push_back(index(j*n+i),Q(-1));}
    auto carrier=kernel(wedge,opt);std::vector<Mat> ambient;for(const auto& x:gens)ambient.push_back(kron(x,x));
    auto recovered=carrier_actions(carrier,gens,gens,false,opt),full=induced(carrier,ambient,opt);
    for(size_t j=0;j<full.size();++j)require(equal(recovered[j],full[j]),"carrier action differs from full tensor action");
    field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));std::vector<ModMat> modular;
    for(const auto& x:gens)modular.push_back(modular_matrix(x,fp));
    for(bool dense:{false,true}){
        auto c=carrier;
        if(dense){Mat change(c.nrow,c.nrow);for(size_t i=0;i<c.nrow;++i)for(size_t j=0;j<c.nrow;++j)change[i].push_back(index(j),Q(i==j?2:1));c=mul(change,c);}
        auto exact=carrier_actions(c,gens,gens,false,opt);auto actual=carrier_actions_mod(c,modular,modular,false,fp,opt);
        for(size_t h=0;h<actual.size();++h){auto expected=modular_matrix(exact[h],fp);require(actual[h].nrow==expected.nrow&&actual[h].ncol==expected.ncol,"modular action shape");
            for(size_t i=0;i<expected.nrow;++i){require(actual[h][i].nnz()==expected[i].nnz(),"modular action nnz");
                for(size_t k=0;k<expected[i].nnz();++k)require(actual[h][i](k)==expected[i](k)&&actual[h][i][k]==expected[i][k],"modular action coefficient");}}
        test_orbits(g,ms,exact,opt);
    }
    auto kernel_frame=carrier_frame(g,ms,recovered,opt);
    require(row_basis(kernel_frame.basis).nrow==carrier.nrow,"kernel carrier frame lost rank");
    verify_adaptation({mul(kernel_frame.basis,carrier),kernel_frame.layout},g,ms,ambient);
    TensorAdapter cg(g,ms,opt);auto t=cg.product(a.layout,a.layout);std::vector<Mat> tg;
    for(auto h:g.generator_ids)tg.push_back(kron(action(a.layout,ms,h),action(a.layout,ms,h)));
    verify_adaptation(t,g,ms,tg);require(equal(mul(t.basis,block_inverse(t.basis,opt)),identity(t.basis.nrow)),"tensor inverse");
    for(size_t s=0;s<ms.size();++s){
        size_t d=ms[s].dimension;Layout source;source.copies.assign(ms.size(),0);source.copies[s]=3;
        Layout target;target.copies.assign(ms.size(),0);target.copies[s]=1;
        Mat center(d,d);
        for(size_t i=0;i<d;++i){std::map<I,Q> row;for(auto h:g.classes.back())for(auto [j,v]:ms[s].actions[h][i])row[j]+=v;center[i]=packed(row);}
        if(name=="Q8"&&ms[s].division_dimension==4) {
            // Reynolds-project E_01 into the commutant. This exercises an
            // actual quaternionic coefficient, beyond scalar copy coupling.
            for(size_t i=0;i<d;++i){std::map<I,Q> row;
                for(size_t h=0;h<g.order();++h)if(auto v=ms[s].actions[g.inverses[h]][i].find(0))
                    for(auto [j,x]:ms[s].actions[h][1])row[j]+=*v*x;
                center[i]=packed(row);}
            require(center.nnz()>0,"zero quaternionic commutant fixture");Q trace=0;
            for(size_t i=0;i<d;++i)if(auto v=center[i].find(index(i)))trace+=*v;
            require(trace==0,"quaternionic fixture should be non-scalar");
        }
        center=transpose(center);
        Mat condition(d,3*d);for(size_t i=0;i<d;++i){condition[i].push_back(index(i),Q(1));for(auto [j,v]:center[i])condition[i].push_back(index(d+j),v);condition[i].push_back(index(2*d+i),Q(-2));}
        auto sg=generators(g,ms,source),qg=generators(g,ms,target);verify_equivariance(condition,sg,qg);
        auto k=kernel_adapted(condition,g,ms,source,opt);require(k.solution.basis.nrow==2*d&&k.solution.layout.copies[s]==2,"wrong irrep kernel multiplicity");
        verify_adaptation(k.solution,g,ms,sg);Chart baseline(kernel(condition,opt),opt);baseline.coordinates(k.solution.basis);
        // A second solve consumes the previous canonical representation directly.
        Mat second(d,2*d);for(size_t i=0;i<d;++i){second[i].push_back(index(i),Q(1));second[i].push_back(index(d+i),Q(-3));}
        auto kg=generators(g,ms,k.solution.layout);verify_equivariance(second,kg,qg);
        auto k2=kernel_adapted(second,g,ms,k.solution.layout,opt);require(k2.solution.basis.nrow==d,"recursive kernel");
        auto full=mul(k2.solution.basis,k.solution.basis);require(!mul(condition,transpose(full)).nnz(),"recursive original residual");
        auto zero=kernel_adapted(identity(3*d),g,ms,source,opt);require(zero.solution.basis.nrow==0&&zero.solution.basis.ncol==3*d,"zero kernel shape");
        auto free=kernel_adapted(Mat(0,3*d),g,ms,source,opt);require(free.solution.basis.nrow==3*d,"unconstrained space");
        // Dense rational conjugation prevents reliance on permutation/signed-permutation inputs.
        Mat change=identity(3*d);for(size_t i=0;i+1<3*d;++i)change[i].push_back(index(i+1),Q(long(i%3)+1,2));auto inv=sparse_mat_inverse(change,QQ,opt);
        std::vector<Mat> scrambled;for(auto& r:sg)scrambled.push_back(mul(mul(inv,r),change));
        auto dec=adapt(g,ms,scrambled,opt);verify_adaptation(dec,g,ms,scrambled);
        auto sparse_frame=carrier_frame(g,ms,scrambled,opt);verify_adaptation(sparse_frame,g,ms,scrambled);
        require(row_basis(sparse_frame.basis).nrow==3*d,"scrambled carrier frame lost rank");
        test_orbits(g,ms,scrambled,opt);
        auto original=mul(condition,transpose(inv));auto independent=row_basis(original);auto eq=equation_actions(independent,scrambled,opt);verify_equivariance(independent,scrambled,eq);
        auto kk=kernel_adapted(mul(original,transpose(dec.basis)),g,ms,dec.layout,opt);
        auto physical=mul(kk.solution.basis,dec.basis);require(physical.nrow==2*d&&!mul(original,transpose(physical)).nnz(),"scrambled kernel");
    }
    if(name=="C3"||name=="D3"||name=="D7"||name=="Q8")test_recurrence(g,ms,a.layout,opt);
    std::cout<<"PASS "<<name<<" order="<<order<<" rational_irreps="<<ms.size()<<std::endl;
}
int main(int argc,char** argv){try{
    rref_option_t opt;opt->pool.reset(2);opt->verbose=false;
    if(argc==6&&std::string(argv[1])=="verify-actions"){
        fs::path prepared=argv[2],chain=argv[3],exported=argv[4];bool carrier=std::string(argv[5])=="carrier";
        auto g=load_group(prepared);auto ms=load_models(prepared,g,opt);auto p=load_prepared(prepared,ms);
        bool factorized=false,compact=false;size_t maximum=chain_weight(chain,prepared,&factorized,&compact);verify_seal(exported);
        TensorAdapter adapter(g,ms,opt);std::vector<Mat> previous=generators(g,ms,p.terminal),letters=generators(g,ms,p.alphabet);
        for(size_t w=1;w<=maximum;++w){auto layout=load_layout(chain/(weight_name(w)+"_copies.tsv"),ms);
            std::vector<Mat> actual;for(size_t j=0;j<g.names.size();++j)actual.push_back(read_matrix(exported/(weight_name(w)+"_g"+std::to_string(j)+".wxf")));
            if(factorized){auto gs=load_carrier_actions(chain,w,g,layout.dimension(ms),opt);
                if(carrier){for(size_t j=0;j<gs.size();++j)require(equal(actual[j],gs[j]),"exported carrier action mismatch");}
                else{auto frame=load_carrier_frame(chain,w,g,ms,gs,opt);
                    for(size_t j=0;j<gs.size();++j)require(equal(mul(actual[j],frame.basis),mul(frame.basis,gs[j])),"exported adapted action fails frame intertwiner");}
            }else{auto tensor=compact?materialize_multiplicity(chain,w,p,ms,adapter,opt):read(chain/(weight_name(w)+".wxf"));auto basis=flatten(tensor);
                for(size_t j=0;j<letters.size();++j){auto ambient=p.backward?kron(letters[j],previous[j]):kron(previous[j],letters[j]);
                    require(equal(mul(actual[j],basis),mul(basis,ambient)),"exported action fails recursive tensor intertwiner");}previous=actual;
            }
        }std::cout<<"PASS exported exact recursive symmetry actions"<<std::endl;return 0;
    }
    if(argc==3&&std::string(argv[1])=="fixtures"){
        fs::path root=argv[2];fresh(root);auto gens=dihedral(3);write(root/"cycle.wxf",gens[0]);write(root/"flip.wxf",gens[1]);
        Mat a(1,3);for(size_t i=0;i<3;++i)a[0].push_back(index(i),Q(1));write(root/"sum.wxf",a);
        Mat bad(1,3);bad[0].push_back(0,Q(1));write(root/"bad.wxf",bad);write(root/"identity.wxf",identity(3));write(root/"empty.wxf",Mat(0,3));
        write(root/"identity2.wxf",identity(2));write(root/"empty2.wxf",Mat(0,2));write(root/"empty0.wxf",Mat(0,0));
        write(root/"terminal.wxf",identity(1));
        write(root/"fec.wxf",tensor(identity(3),{3,1,3}));write(root/"lec.wxf",tensor(identity(3),{3,3,1}));
        Mat wedge(3,9);size_t q=0;for(size_t i=0;i<3;++i)for(size_t j=i+1;j<3;++j){wedge[q].push_back(index(i*3+j),Q(1));wedge[q++].push_back(index(j*3+i),Q(-1));}
        write(root/"condition.wxf",condition_tensor(wedge,3));write(root/"full_condition.wxf",condition_tensor(identity(9),3));
        write(root/"no_condition.wxf",condition_tensor(Mat(0,9),3));
        write(root/"vector_seed.wxf",tensor(identity(9),{9,3,3}));return 0;
    }
    test_incremental_wxf();test_reconstruction_paths(opt);test_direct_kernel(opt);test_carrier_coordinates(opt);
    for(size_t n:{3,5,7}){test("C"+std::to_string(n),dihedral(n,false),n,opt);test("D"+std::to_string(n),dihedral(n),2*n,opt);}
    {auto gs=dihedral(7);for(auto& x:gs)x=kron(x,identity(2));gs.push_back(kron(identity(7),permutation({1,0})));
        test("D7xC2",std::move(gs),28,opt);}
    test("S4",{permutation({1,2,3,0}),permutation({1,0,2,3})},24,opt);
    test("C2xC2",{permutation({1,0,3,2}),permutation({2,3,0,1})},4,opt);
    Mat qi(4,4),qj(4,4);qi[0].push_back(1,Q(1));qi[1].push_back(0,Q(-1));qi[2].push_back(3,Q(1));qi[3].push_back(2,Q(-1));
    qj[0].push_back(2,Q(1));qj[1].push_back(3,Q(-1));qj[2].push_back(0,Q(-1));qj[3].push_back(1,Q(1));test("Q8",{qi,qj},8,opt);
    bool caught=false;try{Mat two(1,1);two[0].push_back(0,Q(2));Group infinite({"scale"},{two},16);}catch(const std::exception&){caught=true;}require(caught,"infinite group not rejected");
    caught=false;try{Group singular({"zero"},{Mat(2,2)});}catch(const std::exception&){caught=true;}require(caught,"singular generator not rejected");
    caught=false;try{Mat bad(1,3);bad[0].push_back(0,Q(1));equation_actions(bad,dihedral(3),opt);}catch(const std::exception&){caught=true;}require(caught,"symmetry-breaking constraints not rejected");
    {auto gs=dihedral(3);Group group({"cycle","flip"},gs);auto model=models(group,opt);auto a=adapt(group,model,gs,opt);
        Mat bad(1,3);bad[0].push_back(0,Q(1));caught=false;
        try{kernel_adapted(mul(bad,transpose(a.basis)),group,model,a.layout,opt);}catch(const std::exception&){caught=true;}
        require(caught,"adapted API accepted a non-invariant kernel");
        Mat change=identity(3);change[0].push_back(1,Q(1,long(n_nextprime(1ULL<<60,0))));
        auto inv=identity(3);inv[0].push_back(1,-change[0][1]);std::vector<Mat> bad_prime;
        for(const auto& x:gs)bad_prime.push_back(mul(mul(inv,x),change));
        auto frame=carrier_frame(group,model,bad_prime,opt);
        require(row_basis(frame.basis).nrow==3,"bad-prime retry lost rank");
        verify_adaptation(frame,group,model,bad_prime);
        test_orbits(group,model,bad_prime,opt);
        field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));auto c=identity(3);sparse_vec_rescale(c[0],Q(long(fp.mod.n)),QQ);
        std::vector<ModMat> ambient,terminal;for(const auto& x:gs){ambient.push_back(modular_matrix(x,fp));terminal.push_back(modular_matrix(identity(1),fp));}
        caught=false;try{carrier_actions_mod(c,terminal,ambient,false,fp,opt);}catch(const BadCarrierPrime&){caught=true;}
        require(caught,"singular modular private chart was not rejected");
    }
    Mat long_basis(2,3);long_basis[0].push_back(0,Q(1));long_basis[0].push_back(2,Q(1,1009));long_basis[1].push_back(1,Q(1));long_basis[1].push_back(2,Q(101,1009));
    auto shortened=short_kernel_basis(long_basis,opt);Chart(long_basis,opt).coordinates(shortened);
    require(row_basis(shortened).nrow==2,"basis reduction lost rank");for(auto& r:shortened.rows)for(auto [j,x]:r)require(x.den()==1,"nonintegral reduced lattice");
    auto large=kron(identity(65),long_basis);auto batched=short_kernel_basis(large,opt);
    Chart(large,opt).coordinates(batched);require(row_basis(batched).nrow==large.nrow,"batched basis reduction lost rank");
    std::cout<<"ALL SYMREP CHECKS PASSED"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
