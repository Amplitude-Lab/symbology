#pragma once
#include "irr_conversion.hpp"
#include "irr_profile.hpp"
#include <memory>

namespace irr {
// One store is shared by all pipeline splits. Cached tensors are immutable;
// eviction never invalidates a caller's shared ownership.
struct recurrence_store {
    using residues=std::array<sparse_mat<ulong,int32_t>,6>;
    template<class T> struct entry {std::shared_ptr<const T> value;size_t bytes,tick;};
    std::map<std::string,entry<matrices>> qq;
    std::map<std::string,entry<residues>> fp;
    struct identity {uintmax_t size;std::filesystem::file_time_type time;std::string hash;};
    std::map<std::string,identity> identities;
    size_t qq_limit,fp_limit,qq_bytes=0,fp_bytes=0,tick=0;
    size_t qq_hits=0,qq_misses=0,qq_bypasses=0,qq_blocks_read=0,fp_hits=0,fp_misses=0,fp_bypasses=0,evictions=0;
    ulong prime=0;
    recurrence_store(size_t rational_budget,size_t modular_budget):qq_limit(rational_budget),fp_limit(modular_budget){}
    std::string file_hash(const std::filesystem::path& file){
        auto path=std::filesystem::weakly_canonical(file).generic_string();auto size=std::filesystem::file_size(file);auto time=std::filesystem::last_write_time(file);
        auto it=identities.find(path);if(it==identities.end()||it->second.size!=size||it->second.time!=time)identities[path]={size,time,hex_sha256_file(file)};
        return identities.at(path).hash;
    }
    std::string key(const std::filesystem::path& root,size_t w,bool first){
        auto dir=root/joint_eet_weight_name(w);std::string result=std::filesystem::weakly_canonical(root).generic_string()+":"+std::to_string(first)+":"+std::to_string(w);
        std::vector<std::filesystem::path> files{dir/(first?"multiplicities.tsv":"blocks.tsv")};
        if(w>1)files.push_back(root/joint_eet_weight_name(w-1)/(first?"multiplicities.tsv":"blocks.tsv"));
        else if(!first)files.push_back(root/"LEC_0_blocks.tsv");
        for(size_t b=0;b<6;++b)files.push_back(tensor_path(dir/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+specs()[b].name+".wxf")));
        for(const auto& file:files)result+=':'+file_hash(file);
        return result;
    }
    template<class T> void retain(std::map<std::string,entry<T>>& cache,const std::string& key,std::shared_ptr<const T> value,size_t bytes,size_t limit,size_t& used){
        if(bytes>limit)return;
        while(bytes>limit-used){auto victim=cache.end();for(auto it=cache.begin();it!=cache.end();++it)
            if(it->second.value.use_count()==1&&(victim==cache.end()||it->second.tick<victim->second.tick))victim=it;
            if(victim==cache.end())return;used-=victim->second.bytes;cache.erase(victim);++evictions;
        }
        cache.emplace(key,entry<T>{std::move(value),bytes,++tick});used+=bytes;
    }
    std::shared_ptr<const matrices> rational(const std::filesystem::path& root,size_t w,bool first,const channels_t& cg,thread_pool* pool){
        auto id=key(root,w,first);auto it=qq.find(id);if(it!=qq.end()){++qq_hits;it->second.tick=++tick;return it->second.value;}
        ++qq_misses;irr_run::phase timer("recurrence_read_decode_CSR");
        auto value=std::make_shared<matrices>(read_recurrence(root,w,first,cg,pool,true));qq_blocks_read+=6;size_t bytes=sizeof(matrices);
        for(const auto& block:*value){bytes+=block.rows.capacity()*sizeof(block.rows[0])+block.alloc()*(sizeof(rat_t)+sizeof(int32_t));
            // Conservative GMP limb payload estimate; allocator/RSS overhead is
            // reported separately, not promised to fit this cache budget.
            for(const auto& row:block.rows)for(size_t p=0;p<row.nnz();++p){const auto& x=row[p];bytes+=((x.num().bits()+63)/64+(x.den().bits()+63)/64)*sizeof(ulong);}
        }
        retain<matrices>(qq,id,value,bytes,qq_limit,qq_bytes);if(!qq.count(id))++qq_bypasses;return value;
    }
    void select_prime(ulong p){if(prime!=p){fp.clear();fp_bytes=0;prime=p;}}
    std::shared_ptr<const residues> modular(const std::filesystem::path& root,size_t w,bool first,const channels_t& cg,const field_t& F,thread_pool* pool){
        select_prime(F.mod.n);auto id=key(root,w,first);auto it=fp.find(id);if(it!=fp.end()){++fp_hits;it->second.tick=++tick;return it->second.value;}
        ++fp_misses;std::shared_ptr<const matrices> input;
        if(qq_limit)input=rational(root,w,first,cg,pool);
        irr_run::phase timer("recurrence_modular_total");
        auto value=std::make_shared<residues>();size_t bytes=sizeof(residues);double read_seconds=0,csr_seconds=0,mod_seconds=0;
        auto old=chain_copies(root,w-1,first),m=chain_copies(root,w,first);
        if(!input){++qq_misses;++qq_bypasses;}
        for(size_t b=0;b<6;++b){
            sparse_mat<rat_t,int32_t> temporary;
            if(!input){irr_run::phase read("recurrence_block_read_decode_CSR");
                auto path=root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+specs()[b].name+".wxf");
                auto start=std::chrono::steady_clock::now();auto tensor=read_tensor(path,pool);require(tensor.rank()==2,"Recurrence must be rank two");
                read_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
                temporary=joint_eet_flatten_tensor(std::move(tensor),pool);csr_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();++qq_blocks_read;}
            const auto& matrix=input?(*input)[b]:temporary;
            size_t columns=0;for(size_t a=0;a<6;++a)columns+=old[a]*cg.channels[b][a].size();
            require(matrix.nrow==m[b]&&matrix.ncol==columns,"Recurrence dimensions differ");
            for(const auto& row:matrix.rows)for(size_t j=0;j<row.nnz();++j)require(row[j].den()%F.mod!=0,"Bad prime: recurrence denominator vanishes");
            auto start=std::chrono::steady_clock::now();(*value)[b]=nmhv_sew_modular_reduce(matrix,F,pool);
            mod_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            bytes+=(*value)[b].rows.capacity()*sizeof((*value)[b].rows[0])+(*value)[b].alloc()*(sizeof(ulong)+sizeof(int32_t));}
        if(!input)std::cout<<"RECURRENCE_IO first="<<first<<" weight="<<w<<" read_decode_seconds="<<read_seconds<<" csr_seconds="<<csr_seconds<<std::endl;
        std::cout<<"RECURRENCE_MOD first="<<first<<" weight="<<w<<" prime="<<F.mod.n<<" seconds="<<mod_seconds<<std::endl;
        retain<residues>(fp,id,value,bytes,fp_limit,fp_bytes);if(!fp.count(id))++fp_bypasses;return value;
    }
    void report()const{
        auto mem=hex_platform::current_process_memory();
        std::cout<<"SHARED_RECURRENCE prime="<<prime<<" qq_entries="<<qq.size()<<" qq_bytes="<<qq_bytes<<" qq_limit="<<qq_limit
            <<" qq_hits="<<qq_hits<<" qq_reads="<<qq_misses<<" qq_blocks_read="<<qq_blocks_read<<" qq_bypasses="<<qq_bypasses<<" fp_entries="<<fp.size()<<" fp_bytes="<<fp_bytes
            <<" qq_remaining="<<qq_limit-qq_bytes<<" fp_limit="<<fp_limit<<" fp_remaining="<<fp_limit-fp_bytes
            <<" fp_hits="<<fp_hits<<" fp_misses="<<fp_misses<<" fp_bypasses="<<fp_bypasses<<" evictions="<<evictions
            <<" rss_bytes="<<mem.working_set_bytes<<" process_peak_rss_bytes="<<mem.peak_working_set_bytes<<std::endl;
    }
};
}
