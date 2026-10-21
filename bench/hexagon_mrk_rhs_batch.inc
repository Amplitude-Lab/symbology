// Included inside namespace irr_physics, after context and symbols().
inline std::vector<std::vector<ulong>> mrk_rhs_reuse_batch(context& c,size_t degree){
    require(!BS::this_thread::get_pool().has_value(),"MRK RHS batch must run outside a worker");
    const auto& data=symbols(c);const size_t loop=(c.h+c.k)/2;
    size_t first=c.sample_offset,count=c.samples;
    if(c.samples==1){first-=first%4;count=4;}
    using key=std::tuple<ulong,uint64_t,size_t,size_t,size_t>;
    auto make_key=[&](size_t sample){return key{c.F.mod.n,c.seed,loop,degree,sample};};
    const bool cacheable=count<=128;
    if(c.mrk_rhs_cache.size()+count>128)c.mrk_rhs_cache.clear();
    std::vector<size_t> missing;
    for(size_t i=0;i<count;++i)if(!c.mrk_rhs_cache.count(make_key(first+i)))missing.push_back(first+i);
    std::vector<std::vector<ulong>> generated(missing.size());
    auto started=c.now();
    c.opt->pool.detach_loop(size_t(0),missing.size(),[&](size_t i){
        auto slots=nmhv_physical_sample_mrk_slots(c.seed,missing[i],2*loop-1,c.F);
        generated[i]=hexagon_mrk_reuse::rhs(loop,slots,data,degree,c.F);
    });
    c.opt->pool.wait();
    if(!cacheable){
        require(first==c.sample_offset&&generated.size()==c.samples,"MRK uncached batch layout");
        std::cout<<"MRK_RHS_REUSE requested="<<c.samples<<" offset="<<c.sample_offset
                 <<" generated="<<missing.size()<<" cache_entries=0 bypass=1 seconds="
                 <<std::chrono::duration<double>(c.now()-started).count()<<std::endl;
        return generated;
    }
    for(size_t i=0;i<missing.size();++i)c.mrk_rhs_cache.emplace(make_key(missing[i]),std::move(generated[i]));
    std::vector<std::vector<ulong>> result(c.samples);
    for(size_t i=0;i<c.samples;++i)result[i]=c.mrk_rhs_cache.at(make_key(c.sample_offset+i));
    std::cout<<"MRK_RHS_REUSE requested="<<c.samples<<" offset="<<c.sample_offset
             <<" generated="<<missing.size()<<" prefetch_width="<<(c.samples==1?4:c.samples)
             <<" cache_entries="<<c.mrk_rhs_cache.size()<<" seconds="
             <<std::chrono::duration<double>(c.now()-started).count()<<std::endl;
    return result;
}
