#pragma once
#include "irr_object.hpp"
#include "irr_runtime.hpp"
#include "platform_resources.hpp"
#include <map>
#ifndef _WIN32
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#endif

namespace irr_stage {
namespace fs=std::filesystem;
using irr::require;
using fields=std::map<std::string,std::string>;
inline fields read(const fs::path& path) {
    std::ifstream in(path);require(bool(in),"Missing metadata: "+path.string());fields out;std::string line;
    while(std::getline(in,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();auto tab=line.find('\t');require(tab!=std::string::npos,"Malformed metadata: "+path.string());require(out.emplace(line.substr(0,tab),line.substr(tab+1)).second,"Duplicate metadata key");}return out;
}
inline void write(const fs::path& path,const fields& values) {
    std::ofstream out(path,std::ios::binary);for(const auto& [k,v]:values){require(v.find_first_of("\r\n\t")==std::string::npos,"Invalid metadata value");out<<k<<'\t'<<v<<'\n';}out.close();require(bool(out),"Cannot write "+path.string());
}
inline std::string text(const fs::path& path) {std::ifstream in(path,std::ios::binary);require(bool(in),"Cannot read "+path.string());return {std::istreambuf_iterator<char>(in),{}};}
inline void fixed_text(const fs::path& path,const std::string& value) {
    if(fs::exists(path)){require(text(path)==value,"Identity/profile mismatch: "+path.string());return;}
    std::ofstream out(path,std::ios::binary);out<<value;out.close();require(bool(out),"Cannot write "+path.string());
}
inline bool linked(const fs::path& path) {
#ifdef _WIN32
    auto a=GetFileAttributesW(path.c_str());if(a!=INVALID_FILE_ATTRIBUTES&&(a&FILE_ATTRIBUTE_REPARSE_POINT))return true;
#endif
    return fs::is_symlink(fs::symlink_status(path));
}
inline void no_links(fs::path path) {
    path=fs::absolute(path).lexically_normal();for(auto p=path;!p.empty();p=p.parent_path()){require(!linked(p),"Refusing link/reparse point: "+p.string());if(p==p.parent_path())break;}
}
// An OS-held lock also excludes cleanup. A stale pathname is harmless after exit.
struct lock {
#ifdef _WIN32
    HANDLE handle=INVALID_HANDLE_VALUE;
    explicit lock(const fs::path& p,bool shared=false){handle=CreateFileW(p.c_str(),shared?GENERIC_READ:GENERIC_READ|GENERIC_WRITE,shared?FILE_SHARE_READ:0,nullptr,OPEN_ALWAYS,FILE_ATTRIBUTE_NORMAL,nullptr);require(handle!=INVALID_HANDLE_VALUE,"Stage busy: "+p.string());}
    ~lock(){if(handle!=INVALID_HANDLE_VALUE)CloseHandle(handle);}
#else
    int fd=-1;
    explicit lock(const fs::path& p,bool shared=false){fd=::open(p.c_str(),O_CREAT|O_RDWR,0600);if(fd<0||flock(fd,(shared?LOCK_SH:LOCK_EX)|LOCK_NB)){if(fd>=0)::close(fd);throw std::runtime_error("Stage busy: "+p.string());}}
    ~lock(){if(fd>=0){flock(fd,LOCK_UN);::close(fd);}}
#endif
    lock(const lock&)=delete;lock& operator=(const lock&)=delete;
};
inline std::string history(const fs::path& root) {
    if(fs::exists(root/"stage.tsv"))return read(root/"stage.tsv").at("history");
    if(fs::exists(root/"condition_history.tsv"))return read(root/"condition_history.tsv").at("history");return "";
}
inline std::string append(const std::string& old,const std::string& op){return old.empty()?op:old+"_"+op;}
inline fields basis_identity(const fs::path& fec,const fs::path& lec,size_t h,size_t k) {
    fields out;
    for(bool first:{true,false}){auto root=first?fec:lec;std::string prefix=first?"FEC/":"LEC/";
        for(size_t w=1;w<=(first?h:k);++w){auto dir=joint_eet_weight_name(w);auto meta=dir+"/"+(first?"multiplicities.tsv":"blocks.tsv");out[prefix+meta]=hex_sha256_file(root/meta);
            for(size_t b=0;b<6;++b){auto name=dir+"/"+(first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf";out[prefix+name]=hex_sha256_file(irr::tensor_path(root/name));}
        }
    }
    for(auto name:{"LEC_0_components.wxf","LEC_0_blocks.tsv","LEC_1_physical.wxf","common/CG.tsv"})out[std::string("LEC/")+name]=hex_sha256_file(lec/name);return out;
}
inline void verify_basis(const fs::path& root,const fields& basis,size_t h,size_t k) {
    if(fs::exists(root/"basis_identity.tsv")){require(read(root/"basis_identity.tsv")==basis,"Source recurrence basis changed");return;}
    require(!fs::exists(root/"stage.tsv"),"Legacy physics QQ lacks recurrence identity");
    std::map<std::string,std::string> expected;
    for(size_t b=0;b<6;++b){auto f="irrFEC_"+std::to_string(h)+"_"+irr::specs()[b].name+".wxf";auto l="irrLEC_"+std::to_string(k)+"_"+irr::specs()[b].name+".wxf";
        expected[f]=basis.at("FEC/"+joint_eet_weight_name(h)+"/"+f);expected[l]=basis.at("LEC/"+joint_eet_weight_name(k)+"/"+l);}
    std::ifstream in(root/"inputs.sha256");require(bool(in),"Legacy bulk lacks boundary provenance");std::string line;std::set<std::string> matched;
    while(std::getline(in,line))if(line.size()>66){if(line.back()=='\r')line.pop_back();auto name=fs::path(line.substr(66)).filename().string();if(name.ends_with(".shape"))name.resize(name.size()-6);
        if(expected.count(name)){require(expected.at(name)==line.substr(0,64),"Legacy boundary recurrence basis changed");matched.insert(name);}}
    require(matched.size()==12,"Legacy bulk must identify all twelve distinct boundary recurrences");
}
inline bool plain(const std::string& name){return !name.empty()&&name!="."&&name!=".."&&name.find_first_of("/\\:\r\n\t")==std::string::npos;}
inline void verify_complete(const fs::path& root) {
    auto receipt=read(root/"qq_complete.tsv");require(receipt.at("schema")=="irr-physics-qq-v1","Unknown QQ receipt");
    require(receipt.at("identity")==hex_sha256_file(root/"input_identity.txt"),"QQ identity mismatch");
    size_t files=0;for(const auto& [name,hash]:receipt)if(name!="schema"&&name!="identity"){
        require(plain(name),"Unsafe receipt path");no_links(root/name);require(hex_sha256_file(root/name)==hash,"QQ checksum mismatch: "+name);++files;
    }
    require(files>=5&&receipt.count("object.tsv")&&receipt.count("columns.tsv")&&receipt.count("manifest.json")&&receipt.count("stage.tsv"),"Incomplete QQ receipt");
    auto obj=read(root/"object.tsv");size_t n=irr_io::number(obj.at("directions"));auto kind=obj.at("kind");
    require(kind=="affine"||kind=="unique"||kind=="homogeneous","Unknown sewing object kind");
    if(kind!="homogeneous")require((n==0)==(kind=="unique"),"Invalid affine object kind");
    for(size_t b=0;b<6;++b){auto stem=irr_io::stem(irr_io::number(obj.at("h")),irr_io::number(obj.at("k")),b);
        auto check=[&](const std::string& name){auto actual=irr::tensor_path(root/name).filename().string();require(receipt.count(actual),"QQ receipt omits "+actual);};
        if(n||kind=="homogeneous")check(stem+".wxf");if(kind!="homogeneous")check(stem+"_origin.wxf");
    }
}
// A per-block deletion authorization pins both the exported tensor and the
// original producer inventory. Missing cache files are allowed only afterwards.
inline fields verify_block_cache(const fs::path& root,const std::string& block) {
    bool known=false;for(const auto& spec:irr::specs())if(spec.name==block)known=true;
    require(known,"Unknown cache block");auto dir=root/"cache"/block;no_links(dir);
    auto receipt=read(dir/"qq_block.tsv");require(receipt.at("schema")=="irr-fec-block-cache-v1","Unknown block receipt");
    require(receipt.at("chain_identity")==hex_sha256_file(root.parent_path()/"input_identity.txt"),"Block chain identity mismatch");
    for(auto key:{"output","marker","report"}){
        const auto& name=receipt.at(key);require(plain(name),"Unsafe block output path");no_links(root/name);
        require(hex_sha256_file(root/name)==receipt.at(std::string(key)+"_sha256"),"Changed block output: "+name);
    }
    for(auto name:{"input.txt","reconstruction.txt","used_cache.tsv"}){no_links(dir/name);require(hex_sha256_file(dir/name)==receipt.at(name),"Changed block provenance");}
    return read(dir/"used_cache.tsv");
}
inline void seal_cache(const fs::path& root,bool chain=false) {
    fields files;auto cache=root/"cache";if(!fs::exists(cache)){write(root/"cache_files.tsv",files);return;}
    no_links(cache);
    bool binaries=false;for(auto it=fs::recursive_directory_iterator(cache);it!=fs::recursive_directory_iterator();++it){require(!linked(it->path()),"Cannot seal a linked cache");if(it->is_regular_file()&&it->path().extension()==".bin")binaries=true;}
    if(!binaries){write(root/"cache_files.tsv",files);return;}
    for(auto it=fs::recursive_directory_iterator(cache);it!=fs::recursive_directory_iterator();++it){
        require(!linked(it->path()),"Cannot seal a linked cache");if(!it->is_regular_file())continue;
        if(it->path().parent_path()==cache&&it->path().filename()=="initial_prime.tsv"){
            auto sha=fs::path(it->path().string()+".sha256");auto hash=hex_sha256_file(it->path());
            require(text(sha)==hash+"\n","Damaged initial-prime receipt; retained");
            files["initial_prime.tsv"]=hash;files["initial_prime.tsv.sha256"]=hex_sha256_file(sha);
        }
        bool block_directory=false;for(auto& spec:irr::specs())if(it->path().parent_path()==cache/spec.name)block_directory=true;
        if(chain&&block_directory&&it->path().filename()=="used_cache.tsv"){
            bool cleaned=fs::exists(it->path().parent_path()/"qq_block.tsv");
            if(cleaned)verify_block_cache(root,it->path().parent_path().filename().string());
            for(auto [name,hash]:read(it->path())){require(plain(name),"Invalid producer cache list");auto file=it->path().parent_path()/name;
                if(cleaned&&!fs::exists(file))continue;
                require(hex_sha256_file(file)==hash,"Changed producer cache");files[fs::relative(file,cache).generic_string()]=hash;}
        }else if(!chain&&it->path().extension()==".bin"){
            std::ifstream in(it->path(),std::ios::binary);nmhv_sew_cache_header_t h{};in.read(reinterpret_cast<char*>(&h),sizeof(h));
            auto id=std::stoull(hex_sha256_file(root/"input_identity.txt").substr(0,16),nullptr,16);
            require(bool(in)&&std::memcmp(h.magic,"NMHVSEW",7)==0&&h.version==1&&h.fingerprint==id,"Foreign modular cache; retained");
            require(fs::file_size(it->path())==sizeof(h)+h.nrow*h.ncol*sizeof(ulong),"Truncated modular cache; retained");
            auto hash=hex_sha256_file(it->path());auto sha=fs::path(it->path().string()+".sha256");
            require(fs::exists(sha)&&text(sha)==hash+"\n","Missing/damaged producer cache checksum; retained");
            files[fs::relative(it->path(),cache).generic_string()]=hash;files[fs::relative(sha,cache).generic_string()]=hex_sha256_file(sha);
        }
    }
    write(root/"cache_files.tsv",files);
}
inline void complete(const fs::path& root) {
    seal_cache(root);
    fields receipt{{"schema","irr-physics-qq-v1"},{"identity",hex_sha256_file(root/"input_identity.txt")}};
    for(const auto& e:fs::directory_iterator(root))if(e.is_regular_file()){
        auto name=e.path().filename().string();if(name=="manifest.json"||name=="stage.tsv"||name=="object.tsv"||name=="columns.tsv"||name=="summary.tsv"||name=="cache_files.tsv"||name=="basis_identity.tsv"||name.ends_with(".wxf")||name.ends_with(".wxf.shape"))receipt[name]=hex_sha256_file(e.path());
    }
    write(root/"qq_complete.tsv",receipt);verify_complete(root);
}
// Published only after all six FEC/LEC/bulk blocks have been exported. Keep the
// receipt separate from legacy manifests, so cleanup never guesses completion.
inline void complete_chain(const fs::path& root,bool parent_lock=true) {
    seal_cache(root,parent_lock);
    fields receipt{{"schema","irr-exact-cache-safe-v1"},{"lock_scope",parent_lock?"parent":"self"},{"cache_files.tsv",hex_sha256_file(root/"cache_files.tsv")}};size_t blocks=0;
    for(const auto& e:fs::directory_iterator(root))if(e.is_regular_file()){
        auto name=e.path().filename().string();if(name=="manifest.json"||name=="checksums.sha256"||name=="basis_identity.tsv"||name.ends_with(".wxf")||name.ends_with(".wxf.shape")||name=="columns.tsv"||name=="blocks.tsv"||name=="multiplicities.tsv"){
            receipt[name]=hex_sha256_file(e.path());if(name.ends_with(".wxf")||name.ends_with(".wxf.shape"))++blocks;
        }
    }
    require(blocks==6&&receipt.count("manifest.json")&&receipt.count("checksums.sha256"),"Incomplete chain export");
    write(root/"qq_cache_safe.tsv",receipt);
}
inline void verify_chain(const fs::path& root) {
    auto r=read(root/"qq_cache_safe.tsv");require(r.at("schema")=="irr-exact-cache-safe-v1","Unknown chain receipt");size_t blocks=0;
    require(r.at("lock_scope")=="parent"||r.at("lock_scope")=="self","Unknown lock scope");
    for(auto [name,hash]:r)if(name!="schema"&&name!="lock_scope"){
        require(plain(name),"Unsafe receipt path");no_links(root/name);require(hex_sha256_file(root/name)==hash,"Chain result checksum mismatch");
        if(name.ends_with(".wxf")||name.ends_with(".wxf.shape"))++blocks;
    }
    require(blocks==6&&r.count("manifest.json")&&r.count("checksums.sha256"),"Incomplete chain receipt");
}
}
