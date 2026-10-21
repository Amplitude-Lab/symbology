#pragma once
#include "hex_hash.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace irr_run {
inline void threads(rref_option_t& opt) {
	const char* value=std::getenv("THREADS");if(!value||!*value)value=std::getenv("SLURM_CPUS_PER_TASK");
	if(value&&*value) {
		std::string s=value;size_t used=0,n=std::stoull(s,&used);
		if(used!=s.size()||!n||s[0]=='-')throw std::runtime_error("Invalid THREADS/SLURM_CPUS_PER_TASK");
		opt->pool.reset(n);
	}else opt->pool.reset();
	std::cout<<"threads: "<<opt->pool.get_thread_count()<<std::endl;
}
// Content identity is independent of mount point. Configuration files are editable.
inline void bind(const std::filesystem::path& root,const std::string& label,
	const std::vector<std::filesystem::path>& inputs) {
	std::string identity="irr-bulk-v1 "+label+"\n";
	for(const auto& path:inputs)identity+=hex_sha256_file(path)+"\n";
	auto file=root/"input_identity.txt";
	if(std::filesystem::exists(root)) {
		std::ifstream in(file);std::string old,line;
		while(std::getline(in,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();old+=line+"\n";}
		if(old!=identity)throw std::runtime_error("Input identity mismatch or unmanaged output: "+root.string());
	}else {
		std::filesystem::create_directories(root);std::ofstream out(file,std::ios::binary);out<<identity;out.close();
		if(!out)throw std::runtime_error("Cannot write input identity");
	}
}
inline bool ready(const std::filesystem::path& logical) {
	auto actual=irr::tensor_path(logical),marker=std::filesystem::path(logical.string()+".ready");
	if(!std::filesystem::exists(actual)&&!std::filesystem::exists(marker))return false;
	std::ifstream in(marker);std::string hash;in>>hash;
	if(!in||!std::filesystem::exists(actual)||hash!=hex_sha256_file(actual))
		throw std::runtime_error("Incomplete or changed block (retained for inspection): "+logical.string());
	return true;
}
inline void mark(const std::filesystem::path& logical) {
	std::ofstream out(logical.string()+".ready");out<<hex_sha256_file(irr::tensor_path(logical))<<'\n';out.close();
	if(!out)throw std::runtime_error("Cannot publish block completion");
}
inline void verify_checksums(const std::filesystem::path& root) {
	std::ifstream in(root/"checksums.sha256");if(!in)throw std::runtime_error("Missing output checksums");
	std::string line;while(std::getline(in,line)) {
		if(!line.empty()&&line.back()=='\r')line.pop_back();
		if(line.size()<67||line.substr(64,2)!="  "||hex_sha256_file(root/line.substr(66))!=line.substr(0,64))
			throw std::runtime_error("Output checksum mismatch: "+root.string());
	}
}
}
