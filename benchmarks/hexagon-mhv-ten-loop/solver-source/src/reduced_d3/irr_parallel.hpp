#pragma once
#include <atomic>
#include <exception>
#include <chrono>
#include "platform_resources.hpp"

namespace irr_parallel {
// Independent rows only. No nested submissions to this pool; exceptions must
// reach the caller rather than disappear in a detached worker.
template<class Function>
void rows(size_t count,thread_pool* pool,Function function,size_t maximum_workers=SIZE_MAX) {
	if(!count)return;
	size_t workers=pool?std::min({count,size_t(pool->get_thread_count()),maximum_workers}):1;
	if(workers<=1){for(size_t r=0;r<count;++r)function(r);return;}
	std::atomic<size_t> next=0;std::atomic<bool> failed=false;
	std::vector<std::exception_ptr> errors(workers);
	pool->detach_loop(size_t(0),workers,[&](size_t worker){
		try {while(!failed){auto r=next.fetch_add(1);if(r>=count)break;function(r);}}
		catch(...){errors[worker]=std::current_exception();failed=true;}
	});pool->wait();
	for(const auto& error:errors)if(error)std::rethrow_exception(error);
}

// Stable, non-owning references grouped by destination copy. Eight bytes per
// source nonzero, for ONE input irrep at a time; never duplicate rational data.
struct grouped_entries {
	struct entry { uint32_t row,position; };
	std::vector<size_t> offsets;
	std::vector<entry> entries;

	template<class T>
	grouped_entries(const sparse_mat<T,int32_t>& matrix,const std::vector<size_t>& groups,size_t count):offsets(count+1) {
		if(matrix.nrow>UINT32_MAX)throw std::runtime_error("Source rows exceed grouped index range");
		for(const auto& row:matrix.rows) {
			if(row.nnz()>UINT32_MAX)throw std::runtime_error("Source row exceeds grouped index range");
			for(size_t p=0;p<row.nnz();++p)++offsets[groups.at(row(p))+1];
		}
		for(size_t i=1;i<offsets.size();++i)offsets[i]+=offsets[i-1];
		entries.resize(offsets.back());auto cursor=offsets;
		for(size_t r=0;r<matrix.nrow;++r)for(size_t p=0;p<matrix[r].nnz();++p)
			entries[cursor[groups[matrix[r](p)]]++]={uint32_t(r),uint32_t(p)};
	}
};

struct phase {
	const char* name;
	std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
	explicit phase(const char* n):name(n){}
	~phase(){std::cout<<"PERF phase="<<name<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
		<<" peak_rss_bytes="<<hex_platform::current_process_memory().peak_working_set_bytes<<std::endl;}
};
}
