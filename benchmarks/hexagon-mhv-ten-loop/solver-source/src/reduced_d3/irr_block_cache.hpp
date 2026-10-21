#pragma once
#include "irr_stage.hpp"

namespace irr_block_cache {
namespace fs=std::filesystem;
using irr::require;

inline void sync(const fs::path& file,bool directory=false) {
#ifdef _WIN32
    if(directory)return;
    auto h=CreateFileW(file.c_str(),GENERIC_WRITE,FILE_SHARE_READ,nullptr,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    require(h!=INVALID_HANDLE_VALUE,"Cannot open for flush: "+file.string());bool ok=FlushFileBuffers(h);CloseHandle(h);require(ok,"File flush failed");
#else
    int fd=::open(file.c_str(),O_RDONLY|(directory?O_DIRECTORY:0));require(fd>=0,"Cannot open for fsync: "+file.string());
    bool ok=::fsync(fd)==0;::close(fd);require(ok,"fsync failed: "+file.string());
#endif
}
inline void regular(const fs::path& file) {
    irr_stage::no_links(file);require(fs::is_regular_file(file)&&fs::hard_link_count(file)==1,"Not an exclusive regular file: "+file.string());
}
inline void publish(const fs::path& file,const irr_stage::fields& values) {
    // An interrupted temporary write may be replaced; the committed receipt may not.
    auto tmp=fs::path(file.string()+".tmp");irr_stage::no_links(tmp);if(fs::exists(tmp))regular(tmp);
    require(!fs::exists(file),"Receipt already exists: "+file.string());
    irr_stage::write(tmp,values);sync(tmp);fs::rename(tmp,file);sync(file.parent_path(),true);
}
inline bool prime_name(const std::string& name) {
    if(!name.starts_with("prime_")||!name.ends_with(".bin"))return false;
    auto middle=name.substr(6,name.size()-10);auto sep=middle.find('_');
    return sep!=std::string::npos&&sep>0&&sep+1<middle.size()&&middle.find('_',sep+1)==std::string::npos&&
        std::all_of(middle.begin(),middle.end(),[](char c){return (c>='0'&&c<='9')||c=='_';});
}
inline void clean(const fs::path& root,const std::string& block,const fs::path& logical,const irr_stage::lock&,bool& legacy_verified) {
    auto dir=root/"cache"/block;irr_stage::no_links(dir);if(!fs::exists(dir))return;
    std::cout<<"BLOCK_CACHE_VERIFY "<<dir<<std::endl;
    auto auth=dir/"qq_block.tsv";bool authorized=fs::exists(auth);
    if(authorized)regular(auth);
    auto allowed=authorized?irr_stage::verify_block_cache(root,block):irr_stage::read(dir/"used_cache.tsv");
    size_t primes=0;for(auto [name,hash]:allowed){
        require(irr_stage::plain(name),"Unsafe block cache name");
        if(prime_name(name)){++primes;require(allowed.count(name+".ready"),"Missing prime marker in inventory");}
        else require(name.ends_with(".ready")&&prime_name(name.substr(0,name.size()-6))&&allowed.count(name.substr(0,name.size()-6)),"Unrecognized producer cache file");
    }
    std::ifstream reconstruction(dir/"reconstruction.txt");std::string key;size_t count=0;
    reconstruction>>key>>count;require(bool(reconstruction)&&key=="primes"&&count==primes&&primes>0,"Incomplete reconstruction record");
    // Check the entire local block before deleting any file. Never recurse or
    // follow directory links; a foreign file or a damaged prime preserves it all.
    for(const auto& entry:fs::directory_iterator(dir)){
        regular(entry.path());auto name=entry.path().filename().string();
        require(allowed.count(name)||name=="input.txt"||name=="reconstruction.txt"||name=="used_cache.tsv"||
            name=="qq_block.tsv"||name=="qq_block.tsv.tmp"||name=="cleanup_done.tsv"||name=="cleanup_done.tsv.tmp","Foreign file in block cache: "+name);
    }
    bool empty=true;for(auto [name,hash]:allowed){irr_stage::no_links(dir/name);if(fs::exists(dir/name))empty=false;}
    if(!authorized&&empty&&fs::exists(root/"qq_cache_safe.tsv")){
        // Legacy whole-stage cleanup leaves producer records but no binaries.
        // Verify all six exported results; this branch deletes nothing.
        if(!legacy_verified){irr_stage::verify_chain(root);legacy_verified=true;}
        std::cout<<"BLOCK_CACHE_ALREADY_CLEAN "<<dir<<std::endl;return;
    }
    std::vector<fs::path> disposable;uintmax_t bytes=0;
    for(auto [name,hash]:allowed){auto file=dir/name;irr_stage::no_links(file);
        if(!fs::exists(file)){require(authorized,"Missing uncommitted prime cache; retained");continue;}
        regular(file);require(hex_sha256_file(file)==hash,"Changed prime cache; retained: "+file.string());
        bytes+=fs::file_size(file);disposable.push_back(file);
    }
    if(!authorized)require(irr_run::ready(logical),"Block not ready; cache retained");
    auto actual=irr::tensor_path(logical),marker=fs::path(logical.string()+".ready"),report=root/(block+".tsv");
    require(actual.parent_path()==root,"Block output outside stage");regular(actual);regular(marker);regular(report);
    if(!authorized){
        std::ifstream marker_in(marker);std::string output_hash;marker_in>>output_hash;
        irr_stage::fields receipt{{"schema","irr-fec-block-cache-v1"},{"chain_identity",hex_sha256_file(root.parent_path()/"input_identity.txt")},
            {"output",actual.filename().string()},{"output_sha256",output_hash},
            {"marker",marker.filename().string()},{"marker_sha256",hex_sha256_file(marker)},
            {"report",report.filename().string()},{"report_sha256",hex_sha256_file(report)},
            {"files",std::to_string(disposable.size())},{"bytes",std::to_string(bytes)}};
        for(auto file:{actual,marker,report})sync(file);
        sync(root,true);
        for(auto name:{"input.txt","reconstruction.txt","used_cache.tsv"}){regular(dir/name);receipt[name]=hex_sha256_file(dir/name);sync(dir/name);}
        publish(auth,receipt);
        require(irr_stage::read(auth)==receipt,"Authorization readback mismatch");
    }
    auto receipt=irr_stage::read(auth);
    require(receipt.at("output")==actual.filename().string()&&receipt.at("marker")==marker.filename().string()&&receipt.at("report")==report.filename().string(),"Wrong block receipt");
    auto done=dir/"cleanup_done.tsv";
    if(fs::exists(done))require(irr_stage::read(done).at("authorization_sha256")==hex_sha256_file(auth),"Changed cleanup receipt");
    for(auto file:disposable){regular(file);require(fs::remove(file),"Cannot unlink prime cache");}
    sync(dir,true);
    if(!fs::exists(done))publish(done,{{"schema","irr-fec-block-cleaned-v1"},{"authorization_sha256",hex_sha256_file(auth)},
        {"files",receipt.at("files")},{"bytes",receipt.at("bytes")}});
    std::cout<<"BLOCK_CACHE_CLEAN "<<dir<<" files="<<disposable.size()<<" bytes="<<bytes<<std::endl;
}
}
