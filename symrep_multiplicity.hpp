#pragma once
// Included inside namespace symrep, after recurrence templates and chain I/O.
// For split QQ models, a map between copies is scalar times the identity.
// Keep precisely these scalars throughout the recurrence; expand only on request.
inline bool split_models(const std::vector<Model>& ms) {
    return std::all_of(ms.begin(),ms.end(),[](const Model& m){return m.division_dimension==1;});
}
inline void certify_scalar_corners(const Group& g,const std::vector<Model>& ms) {
    require(split_models(ms),"compact scalar recurrence requires split rational irreps");
    for(const auto& m:ms){auto local=row_basis(g.evaluate(m.primitive,m.actions));Mat orbit(m.dimension,m.dimension);
        require(local.nrow==1,"scalar corner has wrong dimension");
        for(size_t i=0;i<m.dimension;++i)orbit[i]=mul(local[0],m.actions[m.orbit[i]]);
        require(equal(orbit,identity(m.dimension)),"scalar corner does not reproduce the saved model convention");}
}
inline void primitive_integer(Vec& row) {
    int_t den=1,divisor=0;
    for(auto [j,x]:row)den=LCM(den,x.den());
    for(auto [j,x]:row){int_t v=den;v/=x.den();v*=x.num();divisor=GCD(divisor,v);}
    require(divisor!=0,"zero kernel basis vector");sparse_vec_rescale(row,Q(den,divisor),QQ);
}
inline Mat lift_scalar_kernel(const Mat& compact,const Layout& solution,const Layout& domain,const std::vector<Model>& ms) {
    Mat out(solution.dimension(ms),domain.dimension(ms));auto so=solution.offsets(ms),d=domain.offsets(ms);size_t rb=0,cb=0;
    for(size_t s=0;s<ms.size();++s){size_t n=ms[s].dimension;
        require(rb+solution.copies[s]<=compact.nrow,"compact kernel row mismatch");
        for(size_t c=0;c<solution.copies[s];++c)for(size_t i=0;i<n;++i){auto& row=out[so[s]+c*n+i];row.reserve(compact[rb+c].nnz());
            for(auto [j,x]:compact[rb+c]){require(size_t(j)>=cb&&size_t(j)<cb+domain.copies[s],"compact kernel mixes irreps");row.push_back(index(d[s]+(j-cb)*n+i),x);}}
        rb+=solution.copies[s];cb+=domain.copies[s];}
    require(rb==compact.nrow&&cb==compact.ncol,"compact kernel dimensions disagree with layouts");return out;
}
inline Tensor materialize_multiplicity(const fs::path& chain,size_t w,const Prepared& p,const std::vector<Model>& ms,TensorAdapter& adapter,rref_option_t opt) {
    auto stem=weight_name(w);if(fs::exists(chain/(stem+".wxf")))return read(chain/(stem+".wxf"));
    require(w>1,"missing compact chain seed");auto previous=load_layout(chain/(weight_name(w-1)+"_copies.tsv"),ms);
    auto solution=load_layout(chain/(stem+"_copies.tsv"),ms);
    auto domain=p.backward?adapter.product(p.alphabet,previous):adapter.product(previous,p.alphabet);
    auto kernel=lift_scalar_kernel(read_matrix(chain/(stem+"_multiplicity.wxf")),solution,domain.layout,ms);
    auto coefficients=mul(kernel,domain.basis,&opt->pool);size_t letters=p.alphabet.dimension(ms);
    return tensor(coefficients,p.backward?std::vector<size_t>{coefficients.nrow,letters,previous.dimension(ms)}:
        std::vector<size_t>{coefficients.nrow,previous.dimension(ms),letters});
}
inline void extend_multiplicity(const Prepared& p,const Group& g,const std::vector<Model>& ms,size_t maximum,
    const fs::path& root,rref_option_t opt,bool compare,const fs::path& prepared,const fs::path& resume={}) {
    certify_scalar_corners(g,ms);TensorAdapter adapter(g,ms,opt);RecurrenceTemplates recurrence(g,ms,p,adapter,opt);
    auto seed=p.expansion;auto previous=p.seed,older=p.terminal;size_t first=2,letters=p.alphabet.dimension(ms);
    std::ofstream stats(root/"timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tadapt_s\tkernel_s\tordinary_kernel_s\tverification_s\tcorner_columns\ttemplate_s\tcompression_s\trref_s\tshorten_s\tlift_s\tpipeline_s\tplanning_s\tassembly_work_s\trref_work_s\n";
    if(resume.empty()){write(root/"w1.wxf",seed);save_layout(root/"w1_copies.tsv",previous,ms);}
    else {bool factored=false,compact=false;auto last=chain_weight(resume,prepared,&factored,&compact);
        require(!factored,"expand the factorized chain before continuing with the multiplicity backend");
        require(maximum>last,"--max-weight must exceed the saved chain weight");first=last+1;
        require(equal(flatten(read(resume/"w1.wxf")),flatten(p.expansion)),"saved seed differs from prepared seed");
        for(const auto& e:fs::directory_iterator(resume))if(e.is_regular_file()&&e.path().filename().string().starts_with("w"))fs::copy_file(e.path(),root/e.path().filename());
        std::ifstream oldstats(resume/"timings.tsv");require(bool(oldstats),"missing saved timings");std::string line;std::getline(oldstats,line);
        while(std::getline(oldstats,line)){if(line.empty())continue;auto fields=1+std::count(line.begin(),line.end(),'\t');stats<<line;for(size_t j=fields;j<20;++j)stats<<"\tnan";stats<<'\n';}
        previous=load_layout(resume/(weight_name(last)+"_copies.tsv"),ms);
        if(last>1)older=load_layout(resume/(weight_name(last-1)+"_copies.tsv"),ms);
        if(compact&&last>1){ProductCopies channels;auto domain=p.backward?adapter.product(p.alphabet,older,&channels,false):adapter.product(older,p.alphabet,&channels,false);
            recurrence.remember_compact(read_matrix(resume/(weight_name(last)+"_multiplicity.wxf")),previous,domain.layout,channels);}
        if(compare||!compact)seed=materialize_multiplicity(resume,last,p,ms,adapter,opt);
        if(compare||!compact)require(seed.rank()==3&&seed.dim(0)==previous.dimension(ms)&&seed.dim(p.backward?1:2)==letters&&seed.dim(p.backward?2:1)==older.dimension(ms),"saved recurrence dimensions disagree with layouts");
    }
    std::vector<std::unique_ptr<rref_option>> sector_options;
    if(!compare&&first<=maximum)for(size_t i=0;i<opt->pool.get_thread_count();++i){auto local=std::make_unique<rref_option>();
        local->method=opt->method;local->col_weight=opt->col_weight;local->eliminate_one_nnz=opt->eliminate_one_nnz;
        local->shrink_memory=opt->shrink_memory;local->is_back_sub=opt->is_back_sub;sector_options.push_back(std::move(local));}
    for(size_t w=first;w<=maximum;++w){auto tick=std::chrono::steady_clock::now();ProductCopies source_channels,target_channels;
        auto domain=p.backward?adapter.product(p.alphabet,previous,&source_channels,compare):adapter.product(previous,p.alphabet,&source_channels,compare);
        auto target=p.backward?adapter.product(p.equations,older,&target_channels,compare):adapter.product(older,p.equations,&target_channels,compare);
        double adapt_s=seconds(tick);tick=std::chrono::steady_clock::now();
        Layout solution;solution.copies.resize(ms.size());std::vector<Mat> kernels(ms.size());
        std::vector<size_t> block_nnz(ms.size());std::vector<double> assembly_work(ms.size()),rref_work(ms.size());
        size_t nnz=0,nr=0,nc=std::accumulate(domain.layout.copies.begin(),domain.layout.copies.end(),size_t(0));
        double assembly_s=0,rref_s=0,normalization_s=0,pipeline_s=0;
        auto solve=[&](size_t s,Mat&& block,double assembly_time){
            require(block.ncol==domain.layout.copies[s]&&block.nrow==target.layout.copies[s],"compact block/layout mismatch");
            block_nnz[s]=block.nnz();assembly_work[s]=assembly_time;auto start=std::chrono::steady_clock::now();
            kernels[s]=kernel(std::move(block),sector_options.at(SparseRREF::thread_id()).get());rref_work[s]=seconds(start);
        };
        if(!compare){recurrence.assemble(seed,older,previous,source_channels,target_channels,solve);pipeline_s=seconds(tick);
            // Assembly and elimination overlap: report the joint wall time and
            // summed per-sector work separately, rather than adding overlapping times.
            assembly_s=rref_s=std::numeric_limits<double>::quiet_NaN();
        }else{
            auto blocks=recurrence.assemble(seed,older,previous,source_channels,target_channels);assembly_s=seconds(tick);
            auto target_chart=transpose(block_inverse(target.basis,opt));auto to=target.layout.offsets(ms),d=domain.layout.offsets(ms);
            for(size_t s=0;s<ms.size();++s){auto& block=blocks[s];block_nnz[s]=block.nnz();
                require(block.ncol==domain.layout.copies[s]&&block.nrow==target.layout.copies[s],"compact block/layout mismatch");
                Mat left(block.nrow,target_chart.ncol),right(block.ncol,domain.basis.ncol);
                for(size_t i=0;i<left.nrow;++i)left[i]=target_chart[to[s]+i*ms[s].dimension];
                for(size_t i=0;i<right.nrow;++i)right[i]=domain.basis[d[s]+i*ms[s].dimension];
                require(equal(block,assemble_reduced(seed,p.condition,p.backward,left,right,opt)),"compact contraction differs from component contraction");
                auto start=std::chrono::steady_clock::now();kernels[s]=kernel(std::move(block),opt);rref_work[s]=seconds(start);
            }
            rref_s=std::accumulate(rref_work.begin(),rref_work.end(),0.0);pipeline_s=assembly_s+rref_s;
        }
        for(size_t s=0;s<ms.size();++s){nnz+=block_nnz[s];std::cout<<"  irrep="<<s<<" rows="<<target.layout.copies[s]<<" cols="<<domain.layout.copies[s]<<" nnz="<<block_nnz[s]<<std::endl;}
        tick=std::chrono::steady_clock::now();
        for(size_t s=0;s<ms.size();++s){for(auto& row:kernels[s].rows)primitive_integer(row);solution.copies[s]=kernels[s].nrow;nr+=kernels[s].nrow;}
        normalization_s=seconds(tick);
        Mat compact(nr,nc);size_t rb=0,cb=0;
        for(size_t s=0;s<ms.size();++s){for(auto& row:kernels[s].rows){for(size_t j=0;j<row.nnz();++j)row(j)+=index(cb);compact[rb++]=std::move(row);}cb+=kernels[s].ncol;}
        if(w<maximum)recurrence.remember_compact(compact,solution,domain.layout,source_channels);
        double verification_s=0,ordinary_s=0;
        if(compare){tick=std::chrono::steady_clock::now();auto lifted=lift_scalar_kernel(compact,solution,domain.layout,ms);auto coefficients=mul(lifted,domain.basis,&opt->pool);
            auto lg=generators(g,ms,p.alphabet),pg=generators(g,ms,previous);std::vector<Mat> source;
            for(size_t j=0;j<g.names.size();++j)source.push_back(p.backward?kron(lg[j],pg[j]):kron(pg[j],lg[j]));
            verify_adaptation({coefficients,solution},g,ms,source);auto a=assemble(seed,p.condition,p.backward);
            require(!mul(a,transpose(coefficients),&opt->pool).nnz(),"compact recurrence residual is nonzero");verification_s=seconds(tick);
            tick=std::chrono::steady_clock::now();auto baseline=kernel(std::move(a),opt);ordinary_s=seconds(tick);
            tick=std::chrono::steady_clock::now();require(baseline.nrow==coefficients.nrow,"ordinary and compact kernel dimensions differ");Chart chart(baseline,opt);chart.coordinates(coefficients);verification_s+=seconds(tick);
            seed=tensor(coefficients,p.backward?std::vector<size_t>{coefficients.nrow,letters,previous.dimension(ms)}:std::vector<size_t>{coefficients.nrow,previous.dimension(ms),letters});
        }
        auto stem=weight_name(w);write(root/(stem+"_multiplicity.wxf"),compact);save_layout(root/(stem+"_copies.tsv"),solution,ms);
        stats<<std::setprecision(9)<<w<<'\t'<<target.layout.dimension(ms)<<'\t'<<domain.layout.dimension(ms)<<'\t'<<nnz<<'\t'<<solution.dimension(ms)<<'\t'<<assembly_s<<'\t'<<adapt_s<<'\t'<<rref_s+normalization_s<<'\t'<<ordinary_s<<'\t'<<verification_s<<'\t';
        for(size_t s=0;s<ms.size();++s){if(s)stats<<',';stats<<domain.layout.copies[s];}
        stats<<'\t'<<recurrence.template_s<<'\t'<<recurrence.compression_s<<'\t'<<rref_s<<'\t'<<normalization_s<<"\t0\t"<<pipeline_s<<'\t'<<recurrence.compression_s+recurrence.template_s
             <<'\t'<<std::accumulate(assembly_work.begin(),assembly_work.end(),0.0)<<'\t'<<std::accumulate(rref_work.begin(),rref_work.end(),0.0)<<std::endl;
        std::cout<<"weight="<<w<<" dimension="<<solution.dimension(ms)<<" multiplicity_nnz="<<compact.nnz()<<std::endl;
        older=previous;previous=std::move(solution);
    }
    stats.close();require(bool(stats),"cannot write compact timing report");std::ofstream meta(root/"chain.tsv");
    meta<<"symbology-symrep-chain-v4 prepared_sha256 "<<native_cache::digest(prepared/"checksums.tsv")<<' '<<maximum<<" multiplicity\n";
    meta.close();require(bool(meta),"cannot write compact chain metadata");seal(root);
}
inline void expand_multiplicity_chain(const fs::path& prepared,const fs::path& chain,size_t maximum,const Group& g,
    const std::vector<Model>& ms,rref_option_t opt,const fs::path& output) {
    bool compact=false;require(maximum<=chain_weight(chain,prepared,nullptr,&compact)&&compact,"expand requires a compact multiplicity chain covering the requested weights");
    certify_scalar_corners(g,ms);auto p=load_prepared(prepared,ms);TensorAdapter adapter(g,ms,opt);
    std::ofstream stats(output/"timings.tsv"),exports(output/"expansion_timings.tsv");
    stats<<"weight\trows\tcolumns\tnnz\tdimension\tassembly_s\tadapt_s\tkernel_s\tordinary_kernel_s\tverification_s\tcorner_columns\ttemplate_s\tcompression_s\trref_s\tshorten_s\tlift_s\n";
    exports<<"weight\tdimension\texpansion_s\texpanded_nnz\n";
    for(size_t w=1;w<=maximum;++w){auto tick=std::chrono::steady_clock::now();auto t=materialize_multiplicity(chain,w,p,ms,adapter,opt);auto stem=weight_name(w);
        write(output/(stem+".wxf"),t);save_layout(output/(stem+"_copies.tsv"),load_layout(chain/(stem+"_copies.tsv"),ms),ms);
        if(w>1){stats<<w<<"\tnan\t"<<t.dim(1)*t.dim(2)<<"\tnan\t"<<t.dim(0);for(size_t j=5;j<16;++j)stats<<"\tnan";stats<<'\n';}
        exports<<w<<'\t'<<t.dim(0)<<'\t'<<seconds(tick)<<'\t'<<t.nnz()<<'\n';}
    stats.close();exports.close();require(bool(stats)&&bool(exports),"cannot write expansion timings");save_chain_meta(output,prepared,maximum,false);
}
