#pragma once

// Experimental disk-backed CRT. The installed SparseRREF remains untouched.
#include <fstream>
#include <map>
#include "hex_hash.hpp"
#include "irr_parallel.hpp"

namespace k_crt {
using Matrix=sparse_mat<ulong,int32_t>;
using Pivots=std::vector<std::vector<pivot_t<int32_t>>>;
constexpr uint64_t magic=0x4b43525400000001ULL;
inline bool same(const Pivots& a,const Pivots& b){
	if(a.size()!=b.size())return false;
	for(size_t g=0;g<a.size();++g){if(a[g].size()!=b[g].size())return false;for(size_t i=0;i<a[g].size();++i)if(a[g][i].r!=b[g][i].r||a[g][i].c!=b[g][i].c)return false;}
	return true;
}

struct Stream {
	std::fstream f;uint64_t hash=14695981039346656037ULL;
	Stream(const std::filesystem::path& p,bool write):f(p,std::ios::binary|(write?std::ios::out:std::ios::in)){
		if(!f)throw std::runtime_error("Cannot open CRT cache: "+p.string());
	}
	void bytes(char* data,size_t n,bool write){
		if(write)f.write(data,n);else f.read(data,n);
		if(!f)throw std::runtime_error("Incomplete CRT cache or disk I/O error");
		for(size_t i=0;i<n;++i){hash^=static_cast<unsigned char>(data[i]);hash*=1099511628211ULL;}
	}
	uint64_t number(uint64_t x=0,bool write=false){bytes(reinterpret_cast<char*>(&x),sizeof(x),write);return x;}
	void finish(bool write){
		uint64_t expected=hash;
		if(write){f.write(reinterpret_cast<char*>(&expected),8);f.flush();if(!f)throw std::runtime_error("CRT cache flush failed");}
		else{uint64_t actual=0;f.read(reinterpret_cast<char*>(&actual),8);if(!f||actual!=expected||f.peek()!=EOF)throw std::runtime_error("CRT cache checksum mismatch");}
	}
};

inline Pivots header(Stream& s,ulong prime,size_t nr,size_t nc,const Pivots& pivots={},bool write=false){
	auto field=[&](uint64_t x){auto y=s.number(x,write);if(!write&&x!=y)throw std::runtime_error("CRT header mismatch");};
	field(magic);field(prime);field(nr);field(nc);
	uint64_t groups=s.number(pivots.size(),write);if(groups>nr+1)throw std::runtime_error("Invalid pivot groups");
	Pivots result(groups);size_t rank=0;
	for(size_t g=0;g<groups;++g){auto n=s.number(write?pivots[g].size():0,write);if((rank+=n)>std::min(nr,nc))throw std::runtime_error("Invalid cached rank");
		for(size_t i=0;i<n;++i){auto r=s.number(write?pivots[g][i].r:0,write),c=s.number(write?pivots[g][i].c:0,write);
			if(r>=nr||c>=nc)throw std::runtime_error("Invalid cached pivot");result[g].push_back({static_cast<int32_t>(r),static_cast<int32_t>(c)});}}
	return result;
}

inline void save(const std::filesystem::path& path,const Matrix& mat,ulong prime,const Pivots& pivots){
	Stream s(path,true);header(s,prime,mat.nrow,mat.ncol,pivots,true);
	for(size_t r=0;r<mat.nrow;++r){s.number(mat[r].nnz(),true);for(size_t j=0;j<mat[r].nnz();++j){s.number(mat[r](j),true);s.number(mat[r][j],true);}}
	s.finish(true);s.f.close();
	// The marker is published only after the complete file has been closed.
	std::ofstream marker(path.string()+".ready");marker<<"complete\n";marker.flush();
	if(!marker)throw std::runtime_error("Cannot publish CRT cache marker");
}

inline bool reconstruct(const std::vector<std::pair<ulong,std::filesystem::path>>& files,
	sparse_mat<rat_t,int32_t>& result,size_t nr,size_t nc,const Pivots& pivots,const int_t& input_height,thread_pool* pool=nullptr){
	irr_parallel::phase timing("recurrence_crt");
	std::vector<std::unique_ptr<Stream>> streams;
	for(auto& [prime,path]:files){streams.push_back(std::make_unique<Stream>(path,false));if(!same(header(*streams.back(),prime,nr,nc),pivots))throw std::runtime_error("Cached pivots changed");}
	result=sparse_mat<rat_t,int32_t>(nr,nc);int_t modulus=1;
	for(auto& f:files)modulus*=f.first;
	bool ok=true;
	// Sequential checked I/O, bounded batches, parallel independent rows. The
	// 16 MiB work estimate may be exceeded by ONE indivisible row, not per thread.
	using row_values=std::vector<std::vector<std::pair<int32_t,ulong>>>;
	for(size_t first=0;first<nr;) {
		std::vector<row_values> batch;size_t bytes=0;
		do {
			row_values values(files.size());
			for(size_t p=0;p<files.size();++p){auto prime=files[p].first;auto& s=*streams[p];auto n=s.number();
				if(n>nc)throw std::runtime_error("Invalid CRT row size");
				auto& row=values[p];row.reserve(n);
				for(size_t j=0;j<n;++j){auto c=s.number(),v=s.number();if(c>=nc||v>=prime)throw std::runtime_error("Invalid CRT row entry");row.emplace_back(int32_t(c),v);}
				if(!std::is_sorted(row.begin(),row.end()))std::sort(row.begin(),row.end());
				for(size_t j=1;j<row.size();++j)if(row[j-1].first==row[j].first)throw std::runtime_error("Duplicate CRT column");
				bytes+=n*(128+2*((modulus.bits()+7)/8));
			}
			batch.push_back(std::move(values));
		}while(first+batch.size()<nr&&batch.size()<64&&bytes<(16ULL<<20));
		std::vector<unsigned char> success(batch.size(),1);
		if(ok)irr_parallel::rows(batch.size(),pool,[&](size_t r){
			std::map<int32_t,int_t> residues;int_t mod=1;
			for(size_t p=0;p<files.size();++p){const auto& values=batch[r][p];
				for(const auto& [c,v]:values)residues.try_emplace(c,0);
				auto it=values.begin();
				for(auto& [c,v]:residues){while(it!=values.end()&&it->first<c)++it;v=CRT(v,mod,it!=values.end()&&it->first==c?it->second:0,files[p].first);}
				mod*=files[p].first;
			}
			for(auto& [c,v]:residues){rat_t q;if(!rational_reconstruct(q,v,modulus)){success[r]=0;return;}if(q!=0)result[first+r].push_back(c,std::move(q));}
		});
		if(std::find(success.begin(),success.end(),0)!=success.end()){ok=false;result.clear();}
		first+=batch.size();
	}
	for(auto& s:streams)s->finish(false);
	if(!ok)return false;
	size_t parts=std::min(nr,size_t(pool?pool->get_thread_count():1));std::vector<int_t> denominators(parts,1),heights(parts,1);
	irr_parallel::rows(parts,pool,[&](size_t part){for(size_t r=part;r<nr;r+=parts)for(auto [j,x]:result[r])denominators[part]=LCM(denominators[part],x.den());});
	int_t d=1,h=1;for(const auto& x:denominators)d=LCM(d,x);
	irr_parallel::rows(parts,pool,[&](size_t part){for(size_t r=part;r<nr;r+=parts)for(auto [j,x]:result[r]){int_t v=(x*d).num().abs();if(v>heights[part])heights[part]=v;}});
	for(const auto& x:heights)if(x>h)h=x;
	if(input_height*nc*h>=modulus){result.clear();return false;}
	return true;
}

inline Pivots rref(sparse_mat<rat_t,int32_t>& mat,rref_option_t opt,const std::filesystem::path& dir){
	std::filesystem::create_directories(dir);
	// Protect reuse even when invoked outside the Python runner. Hash canonical input values.
	uint64_t hash=14695981039346656037ULL;
	auto feed=[&](const std::string& x){for(unsigned char c:x){hash^=c;hash*=1099511628211ULL;}};
	feed(std::to_string(mat.nrow)+","+std::to_string(mat.ncol)+";");
	for(size_t r=0;r<mat.nrow;++r){mat[r].compress();feed(std::to_string(r)+":");for(size_t j=0;j<mat[r].nnz();++j)feed(std::to_string(mat[r](j))+"="+mat[r][j].get_str()+";");}
	std::string identity="KCRT1 method=1 backward=true input="+std::to_string(hash);
	auto identity_path=dir/"input.txt";
	if(std::filesystem::exists(identity_path)){std::ifstream f(identity_path);std::string old;std::getline(f,old);if(old!=identity)throw std::runtime_error("CRT input identity mismatch");}
	else{std::ofstream f(identity_path);f<<identity<<'\n';f.flush();if(!f)throw std::runtime_error("CRT identity write failed");}
	auto& pool=opt->pool;int_t input_height=sparse_mat_height(mat),modulus=1;
	Pivots pivots;std::vector<std::pair<ulong,std::filesystem::path>> files;
	ulong prime=n_nextprime(1ULL<<60,0);
	for(size_t iteration=1;;++iteration,prime=n_nextprime(prime,0)){
		std::filesystem::path path;
		const std::string prefix="prime_"+std::to_string(prime)+"_";
		for(auto& entry:std::filesystem::directory_iterator(dir))if(entry.path().filename().string().starts_with(prefix)&&entry.path().extension()==".bin"&&std::filesystem::exists(entry.path().string()+".ready")){path=entry.path();break;}
		if(path.empty()){
			field_t F(FIELD_Fp,prime);Matrix modular(mat.nrow,mat.ncol);
			pool.detach_loop(0,mat.nrow,[&](size_t r){modular[r]=mat[r]%F.mod;});pool.wait();
			if(files.empty())pivots=sparse_mat_rref_forward(modular,F,opt);else sparse_mat_direct_rref(modular,pivots,F,opt);
			triangular_solver_2(modular,pivots,F,opt);
			size_t attempt=1;do{path=dir/(prefix+std::to_string(attempt++)+".bin");}while(std::filesystem::exists(path));
			save(path,modular,prime,pivots);modular.clear();
			std::cout<<"CRT_CACHE_SAVED "<<path<<std::endl;
		}else{
			Stream s(path,false);auto cached=header(s,prime,mat.nrow,mat.ncol);if(files.empty())pivots=std::move(cached);else if(!same(cached,pivots))throw std::runtime_error("CRT pivot mismatch");
			std::cout<<"CRT_CACHE_REUSED "<<path<<std::endl;
		}
		files.push_back({prime,path});modulus*=prime;
		// As upstream, retain only the independent original rows for later primes.
		if(files.size()==1){std::vector<bool> selected(mat.nrow,false);for(auto& group:pivots)for(auto [r,c]:group)selected[r]=true;for(size_t r=0;r<mat.nrow;++r)if(!selected[r])mat[r].clear();}
		sparse_mat<rat_t,int32_t> rational;
		if(reconstruct(files,rational,mat.nrow,mat.ncol,pivots,input_height,&pool)){
			mat=std::move(rational);std::ofstream done(dir/"reconstruction.txt");done<<"primes "<<files.size()<<"\nmodulus_bits "<<modulus.bits()<<'\n';
			std::ofstream used(dir/"used_cache.tsv");for(auto& [p,file]:files){used<<file.filename().string()<<'\t'<<hex_sha256_file(file)<<'\n';used<<file.filename().string()<<".ready\t"<<hex_sha256_file(file.string()+".ready")<<'\n';}
			std::cout<<"** Reconstruct success! Using mod ~ 2^"<<modulus.bits()<<std::endl;return pivots;
		}
		std::cout<<">> Reconstruct failed, now mod ~ 2^"<<modulus.bits()<<std::endl;
	}
}
}
