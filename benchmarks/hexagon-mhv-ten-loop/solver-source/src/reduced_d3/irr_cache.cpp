#include "irr_stage.hpp"
#include <chrono>

namespace fs=std::filesystem;
using irr::require;
// Return only regular files in a verified, link-free cache tree. Do not use
// recursive remove_all: enumerate, check, then unlink files individually.
std::vector<fs::path> inventory(const fs::path& root) {
    irr_stage::no_links(root);std::vector<fs::path> out;
    for(auto it=fs::recursive_directory_iterator(root);it!=fs::recursive_directory_iterator();++it){
        require(!irr_stage::linked(it->path()),"Cache contains a link: "+it->path().string());
        if(it->is_regular_file()){require(fs::hard_link_count(it->path())==1,"Cache contains shared hard links");out.push_back(it->path());}
        else require(it->is_directory(),"Cache contains a special file");
    }return out;
}
bool clean(const fs::path& stage,bool apply) {
    if(!fs::exists(stage/"cache"))return false;
    irr_stage::no_links(stage);bool chain=fs::exists(stage/"qq_cache_safe.tsv");
    auto lockroot=chain&&irr_stage::read(stage/"qq_cache_safe.tsv").at("lock_scope")=="parent"?stage.parent_path():stage;
    irr_stage::lock lock(lockroot/".active.lock");
    auto verify=[&](){if(chain)irr_stage::verify_chain(stage);else irr_stage::verify_complete(stage);};verify();
    if(chain)for(const auto& spec:irr::specs())if(fs::exists(stage/"cache"/spec.name/"qq_block.tsv"))irr_stage::verify_block_cache(stage,spec.name);
    auto files=inventory(stage/"cache");uintmax_t bytes=0;
    auto receipt=irr_stage::read(stage/(chain?"qq_cache_safe.tsv":"qq_complete.tsv"));
    require(receipt.count("cache_files.tsv"),"No sealed producer cache inventory; retained");
    auto allowed=irr_stage::read(stage/"cache_files.tsv");
    std::vector<fs::path> disposable;
    for(auto file:files){auto name=file.filename().string();
        if(name=="input.txt"||name=="reconstruction.txt"||name=="used_cache.tsv")continue;
        if(chain&&(name=="qq_block.tsv"||name=="cleanup_done.tsv"))continue;
        auto relative=fs::relative(file,stage/"cache").generic_string();
        require(allowed.count(relative),"File not owned by this completed computation; retained: "+file.string());
        require(allowed.at(relative)==hex_sha256_file(file),"Changed/truncated cache; retained: "+file.string());bytes+=fs::file_size(file);disposable.push_back(file);
    }
    // Partial deletion or missing completion markers also require inspection.
    if(!disposable.empty())for(auto [relative,hash]:allowed){fs::path p=relative;require(p.is_relative()&&p.lexically_normal()==p&&relative.find("..") == std::string::npos,"Unsafe cache inventory path");require(fs::is_regular_file(stage/"cache"/p),"Incomplete cache inventory; retained");}
    std::cout<<(apply?"CLEAN ":"DRY-RUN ")<<stage<<" files="<<disposable.size()<<" bytes="<<bytes<<std::endl;
    if(apply){verify();for(auto file:disposable){irr_stage::no_links(file);require(fs::remove(file),"Cannot remove cache file");}
        // Keep the empty cache directory and a record. No scientific output is removed.
        irr_stage::write(stage/"cache_cleanup.tsv",{{"files",std::to_string(disposable.size())},{"bytes",std::to_string(bytes)},{"verified_receipt",hex_sha256_file(stage/(chain?"qq_cache_safe.tsv":"qq_complete.tsv"))}});
    }return true;
}
void self_test(){
    require(!irr_stage::plain("../bad")&&irr_stage::plain("irrSEW_2p2_pe_A1.wxf"),"Receipt path validation");
    require(irr_stage::append("mrk","collinear")=="mrk_collinear","Ordered paths");
    std::cout<<"CACHE GUARD UNIT TESTS PASS\n";
}
int main(int argc,char** argv)try {
    if(argc==2&&std::string(argv[1])=="--self-test"){self_test();return 0;}
    require(argc>=2,"irr_cache OUTPUT_ROOT [--apply]");bool apply=false;
    for(int i=2;i<argc;++i){require(std::string(argv[i])=="--apply","Unknown option");apply=true;}
    fs::path root=fs::absolute(argv[1]).lexically_normal();irr_stage::no_links(root);require(fs::is_directory(root),"Output root is not a directory");
    std::vector<fs::path> stages;
    for(auto it=fs::recursive_directory_iterator(root);it!=fs::recursive_directory_iterator();++it){
        if(irr_stage::linked(it->path())){it.disable_recursion_pending();std::cout<<"SKIP link "<<it->path()<<std::endl;continue;}
        if(it->is_directory()&&(it->path().filename()=="cache"||it->path().filename()==".git"))it.disable_recursion_pending();
        if(it->is_regular_file()&&(it->path().filename()=="qq_complete.tsv"||it->path().filename()=="qq_cache_safe.tsv"))stages.push_back(it->path().parent_path());
    }
    size_t count=0;for(auto stage:stages)try{count+=clean(stage,apply);}catch(const std::exception& e){std::cout<<"SKIP "<<stage<<": "<<e.what()<<std::endl;}
    std::cout<<"Verified completed stages: "<<count<<"; partial/cache-only/unrecognized stages retained.\n";return 0;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}
