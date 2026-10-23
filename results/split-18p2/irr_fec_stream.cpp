// Streaming FEC chain generator for weights whose rational assembly exceeds RAM.
//
// The reference irr_fec materializes the assembled recurrence matrix over QQ
// (weight-17 step peak: 44.4 GB, dominated by that assembly), then reduces it
// prime by prime through the disk-backed k_crt cache. This variant keeps the
// k_crt cache format, the SparseRREF reduction, the deterministic prime
// sequence and every downstream artifact of the reference pipeline, but never
// materializes the QQ matrix:
//
//   * per prime, the previous chain K is scattered through the template routes
//     directly into modular residues (nmod arithmetic); rows are merged and
//     sorted exactly as the QQ assembly + merge + "%mod" chain would produce
//     (a QQ entry vanishing mod p has probability ~2^-60 per entry and is
//     handled identically by both paths when it does occur);
//   * the row permutation of the first prime's sort is reused for every later
//     prime, matching the reference where the QQ matrix is sorted once;
//   * after the first prime, rows outside the pivot set are exactly zero over
//     QQ and are no longer assembled (the reference clears the same rows of
//     its QQ matrix before reducing further primes);
//   * reconstruction streams the saved per-prime caches once per attempt and
//     extracts kernel columns in bounded free-column slices; each kernel row
//     (one per free column) is canonicalized as soon as its slice completes.
//
// Acceptance keeps the reference's global sufficiency test shape
// (input_height * nc * h < modulus) with h bounded above by
// numerator_max * LCM(denominators); the bound is monotone during an attempt,
// so an insufficient modulus aborts early instead of after a full pass. The
// input height itself is the exact per-entry product bound
// height(K entry) * height(route) * height(row scale), which is never smaller
// than the reference's assembled-matrix height.
//
// The kernel construction reproduces sparse_mat_rref_kernel + transpose:
// kernel row i (free column f = nonpivs[i]) carries -1 at f and the RREF
// entry M[r][f] at every pivot column c(r). Sign, layout and the final
// sort_rows_by_nnz follow the reference, so the written chain files match the
// reference byte for byte for the same inputs.
#include "irr_recurrence.hpp"
#include "joint_eet_io.hpp"
#include "irr_runtime.hpp"
#include "k_crt_cache.hpp"
#include "irr_stage.hpp"
#include "irr_block_cache.hpp"
#include <atomic>
#include <memory>
#include <mutex>

namespace fs=std::filesystem;
const field_t QQ(FIELD_QQ);
using Matrix=sparse_mat<ulong,int32_t>;
using Pivots=k_crt::Pivots;

namespace fec_stream {

struct layout_t {
	size_t nr=0,nc=0;
	std::vector<size_t> rowbase;
	irr::copies colbase{};
};

inline layout_t layout(const irr::channels_t& cg,const irr::template_t& t,size_t target,
	const std::vector<fec_irrep_axis_block_t>& condition,const irr::copies& old,const irr::copies& previous){
	layout_t L;L.rowbase.resize(t.rows.size());
	for(size_t r=0;r<t.rows.size();){
		size_t end=r+1;
		while(end<t.rows.size()&&t.rows[end].left==t.rows[r].left&&t.rows[end].condition==t.rows[r].condition)++end;
		for(size_t s=r;s<end;++s)L.rowbase[s]=L.nr+t.rows[s].copy;
		L.nr+=old[t.rows[r].left]*condition[t.rows[r].condition].multiplicity;
		r=end;
	}
	for(size_t b=0;b<6;++b){L.colbase[b]=L.nc;L.nc+=previous[b]*cg.channels[target][b].size();}
	if(L.nc>size_t(INT32_MAX))throw std::runtime_error("Candidate columns exceed int32");
	return L;
}

inline void merge_mod_row(sparse_vec<ulong,int32_t>& row,nmod_t mod){
	row.sort_indices();
	size_t k=0;
	for(size_t i=0;i<row._nnz;){
		size_t j=i;ulong s=0;
		while(j<row._nnz&&row.indices[j]==row.indices[i]){s=nmod_add(s,row.entries[j],mod);++j;}
		if(s){row.indices[k]=row.indices[i];row.entries[k]=s;++k;}
		i=j;
	}
	row._nnz=k;
}

// Modular image of the reference QQ assembly for one prime, in natural row
// order. skip_mask (natural row indices) drops rows outside the pivot set
// after the first prime; they are exactly zero over QQ by then. Row order is
// never permuted: the reduced echelon form is canonical, so the kernel and
// every written artifact are unchanged, while every prime's assembly stays
// consistent with the cached pivots even across process restarts.
inline Matrix assemble_modular(const irr::channels_t& cg,const irr::template_t& t,size_t target,
	const std::vector<fec_irrep_axis_block_t>& condition,
	const std::vector<std::vector<std::pair<size_t,size_t>>>& decode,
	const std::vector<std::unique_ptr<irr_parallel::grouped_entries>>& index,
	const irr::matrices& K,const layout_t& L,const std::vector<char>& skip_mask,
	rref_option_t opt,const field_t& F){
	irr_parallel::phase timing("modular_assembly");
	const auto& mod=F.mod;
	Matrix modular(L.nr,L.nc);
	std::vector<ulong> scale_mod(t.row_scale.size());
	for(size_t s=0;s<t.row_scale.size();++s)scale_mod[s]=t.row_scale[s]%mod;
	for(size_t b=0;b<6;++b){
		const auto& routes=t.routes[b];
		std::vector<ulong> route_mod;std::vector<size_t> route_off(routes.size()+1,0);
		for(size_t g=0;g<routes.size();++g){
			for(const auto& route:routes[g])route_mod.push_back(route.coefficient%mod);
			route_off[g+1]=route_mod.size();
		}
		auto scatter=[&](size_t owner){
			for(size_t p=index[b]->offsets[owner];p<index[b]->offsets[owner+1];++p){
				auto [a,position]=index[b]->entries[p];
				auto col=K[b][a](position);
				ulong v=K[b][a][position]%mod;
				if(!v)continue;
				auto [g,c]=decode[b].at(col);
				for(size_t j=0;j<routes[g].size();++j){
					const auto& route=routes[g][j];
					const auto& row=t.rows[route.row];
					size_t at=L.rowbase[route.row]+c*condition[row.condition].multiplicity;
					if(!skip_mask.empty()&&!skip_mask[at])continue;
					ulong x=nmod_mul(v,route_mod[route_off[g]+j],mod);
					x=nmod_mul(x,scale_mod[route.row],mod);
					if(x)modular[at].push_back(
						int32_t(L.colbase[b]+a*cg.channels[target][b].size()+route.channel),x);
				}
			}
		};
		irr_parallel::rows(index[b]->offsets.size()-1,&opt->pool,scatter);
		// Merge after each source block: the raw scatter carries duplicates
		// from all six blocks at once, which would spike memory ~6x otherwise.
		opt->pool.detach_loop(0,modular.nrow,[&](size_t r){
			if(modular[r].nnz())merge_mod_row(modular[r],mod);
		});
		opt->pool.wait();
	}
	return modular;
}

inline size_t chain_height_bits(const irr::matrices& K,thread_pool* pool){
	std::atomic<size_t> best{1};
	for(size_t b=0;b<6;++b)
		pool->detach_loop(0,K[b].nrow,[&,b](size_t r){
			size_t local=1;
			for(size_t p=0;p<K[b][r].nnz();++p)local=std::max(local,K[b][r][p].height_bits());
			size_t cur=best.load(std::memory_order_relaxed);
			while(local>cur&&!best.compare_exchange_weak(cur,local,std::memory_order_relaxed)){}
		});
	pool->wait();
	return best.load();
}

// The previous chain is only needed while assembling a prime's modular image.
// Loading it from the exported tensors per prime (instead of retaining ~12 GB
// of rationals through every reduction) keeps the peak at the RREF working
// set itself; the deterministic files make each reload bit-identical.
inline irr::matrices load_chain(const std::vector<fs::path>& chain_files,thread_pool* pool){
	irr::matrices K;
	for(size_t b=0;b<6;++b)K[b]=joint_eet_flatten_tensor(irr::read_tensor(chain_files[b],pool),pool);
	return K;
}
inline void unload_chain(irr::matrices& K){
	for(size_t b=0;b<6;++b)K[b]=sparse_mat<rat_t,int32_t>();
}

struct tagged_entry { int32_t ipos=0,row=0; uint32_t prime=0; ulong value=0; };

// One reconstruction attempt: stream every saved prime cache once, CRT and
// rationally reconstruct the free-column entries slice by slice, and build the
// kernel rows directly. Returns false when the modulus is still insufficient;
// next is then meaningless.
inline bool reconstruct_kernel(const std::vector<std::pair<ulong,fs::path>>& files,
	const Pivots& pivots,size_t nr,size_t nc,const int_t& input_height,const int_t& modulus,
	rref_option_t opt,sparse_mat<rat_t,int32_t>& next,size_t& max_2nd_bits){
	std::vector<pivot_t<int32_t>> flat;
	for(const auto& group:pivots)for(auto p:group)flat.push_back(p);
	std::vector<char> is_pivot_column(nc,0),is_pivot_row(nr,0);
	std::vector<int32_t> pivot_column_of_row(nr,-1),nonpivs_pos(nc,-1);
	for(auto [r,c]:flat){is_pivot_column[c]=1;is_pivot_row[r]=1;pivot_column_of_row[r]=int32_t(c);}
	std::vector<int32_t> nonpivs;
	for(size_t c=0;c<nc;++c)if(!is_pivot_column[c]){nonpivs_pos[c]=int32_t(nonpivs.size());nonpivs.push_back(int32_t(c));}
	const size_t kdim=nonpivs.size();
	next=sparse_mat<rat_t,int32_t>(kdim,nc);
	if(!kdim)return true;

	// Bounded slices over kernel rows: each slice re-streams the saved prime
	// caches and keeps only its own free-column entries, so reconstruction
	// memory stays at one slice regardless of matrix size. The extra NVMe
	// passes replace a full-entry-set buffer that would not fit in RAM.
	const size_t per_slice=std::max<size_t>(1,std::min(kdim,size_t(256)));
	const size_t nslices=(kdim+per_slice-1)/per_slice;

	int_t denominator_lcm=1,numerator_max=1;
	bool ok=true;
	for(size_t s=0;s<nslices;++s){
		const size_t i0=s*per_slice,i1=std::min(kdim,i0+per_slice);
		// One sequential read per prime file for this slice only.
		std::vector<std::vector<tagged_entry>> by_prime(files.size());
		std::atomic<bool> io_ok{true};
		std::mutex io_mutex;
		std::string io_error;
		opt->pool.detach_loop(0,files.size(),[&](size_t p){
			try{
				k_crt::Stream stream(files[p].second,false);
				if(!k_crt::same(k_crt::header(stream,files[p].first,nr,nc),pivots))
					throw std::runtime_error("CRT pivot mismatch during streaming");
				for(size_t r=0;r<nr;++r){
					auto n=stream.number();
					if(n>nc)throw std::runtime_error("Invalid CRT row size");
					if(!is_pivot_row[r]||!n){
						for(size_t j=0;j<2*n;++j)stream.number();
						continue;
					}
					for(size_t j=0;j<n;++j){
						auto c=stream.number(),v=stream.number();
						if(v>=files[p].first)throw std::runtime_error("Invalid CRT row entry");
						if(is_pivot_column[c])continue; // identity part; the kernel needs only free columns
						auto ipos=nonpivs_pos[c];
						if(size_t(ipos)<i0||size_t(ipos)>=i1)continue;
						by_prime[p].push_back({ipos,int32_t(r),uint32_t(p),v});
					}
				}
				stream.finish(false);
			}catch(const std::exception& e){
				io_ok=false;
				std::lock_guard<std::mutex> guard(io_mutex);
				if(io_error.empty())io_error=e.what();
			}
		});
		opt->pool.wait();
		if(!io_ok)throw std::runtime_error(io_error);
		size_t total=0;
		for(size_t p=0;p<files.size();++p)total+=by_prime[p].size();
		std::vector<tagged_entry> entries;
		entries.reserve(total);
		for(size_t p=0;p<files.size();++p){
			entries.insert(entries.end(),by_prime[p].begin(),by_prime[p].end());
			std::vector<tagged_entry>().swap(by_prime[p]);   // free before sorting
		}
		std::sort(entries.begin(),entries.end(),[](const tagged_entry& a,const tagged_entry& b){
			return a.ipos!=b.ipos?a.ipos<b.ipos:(a.row!=b.row?a.row<b.row:a.prime<b.prime); });

		// Group spans of (kernel row, pivot row) in sorted order.
		std::vector<std::pair<size_t,size_t>> span;
		for(size_t e=0;e<entries.size();){
			size_t k=e+1;
			while(k<entries.size()&&entries[k].ipos==entries[e].ipos&&entries[k].row==entries[e].row)++k;
			span.push_back({e,k});
			e=k;
		}
		// Groups sharing a kernel row must stay on one thread: parallelize
		// over kernel rows, never over groups.
		const size_t row_count=i1-i0;
		std::vector<std::pair<size_t,size_t>> row_groups(row_count,{span.size(),span.size()});
		for(size_t i=0,g=0;i<row_count;++i){
			size_t row_ipos=i0+i;
			while(g<span.size()&&size_t(entries[span[g].first].ipos)<row_ipos)++g;
			size_t start=g;
			while(g<span.size()&&size_t(entries[span[g].first].ipos)==row_ipos)++g;
			row_groups[i]={start,g};
		}
		std::vector<rat_t> value(span.size());
		std::atomic<size_t> slice_max2nd{0};
		std::atomic<bool> slice_ok{true};
		opt->pool.detach_loop(0,i1-i0,[&](size_t i){
			auto free_column=nonpivs[i0+i];   // kernel row i0+i carries -1 at this column
			for(size_t g=row_groups[i].first;g<row_groups[i].second;++g){
				int_t residue=0,running=1;
				for(size_t e=span[g].first;e<span[g].second;++e){
					residue=CRT(residue,running,entries[e].value,files[entries[e].prime].first);
					running*=files[entries[e].prime].first;
				}
				// A QQ-nonzero entry vanishes mod some prime with probability
				// ~2^-60 per entry; treat any such gap as a failed attempt.
				if(running!=modulus){slice_ok=false;return;}
				rat_t q;
				// Per-entry uniqueness window 2|num*den| < modulus; a failure
				// here or in the reconstruction rejects the whole attempt.
				if(!rational_reconstruct(q,residue,modulus)||
					2*(q.num().abs()*q.den()).bits()>=modulus.bits()){slice_ok=false;return;}
				size_t bits=(2*(q.num().abs()*q.den())).bits();
				size_t cur=slice_max2nd.load(std::memory_order_relaxed);
				while(bits>cur&&!slice_max2nd.compare_exchange_weak(cur,bits,std::memory_order_relaxed)){}
				value[g]=q;
				next[i0+i].push_back(pivot_column_of_row[size_t(entries[span[g].first].row)],value[g]);
			}
			auto& row=next[i0+i];
			row.push_back(free_column,rat_t(-1));
			row.sort_indices();
			row.canonicalize();
			fec_irrep_cancel_divisor(row);
		});
		opt->pool.wait();
		if(!slice_ok.load()){ok=false;break;}
		// Sequential slot aggregation for the global height bound.
		const size_t parts=std::max<size_t>(1,std::min(span.size(),size_t(opt->pool.get_thread_count())));
		std::vector<int_t> dens(parts,1),nums(parts,1);
		for(size_t g=0;g<span.size();++g){
			size_t part=g%parts;
			dens[part]=LCM(dens[part],value[g].den());
			int_t nabs=value[g].num().abs();
			if(nabs>nums[part])nums[part]=nabs;
		}
		for(size_t part=0;part<parts;++part){
			denominator_lcm=LCM(denominator_lcm,dens[part]);
			if(nums[part]>numerator_max)numerator_max=nums[part];
		}
		const size_t slice_bits=slice_max2nd.load();
		if(slice_bits>max_2nd_bits)max_2nd_bits=slice_bits;
		// Monotone global sufficiency test; aborting early is always safe.
		if(input_height*nc*(numerator_max*denominator_lcm)>=modulus){ok=false;break;}
		std::cout<<"KERNEL_SLICE slice="<<(s+1)<<"/"<<nslices<<" rows=["<<i0<<","<<i1
			<<") entries="<<total<<" max_2nd_bits="<<slice_bits<<std::endl;
	}
	return ok;
}

struct element_result { sparse_mat<rat_t,int32_t> next; size_t nr=0,nc=0,assembled_nnz=0; };

inline element_result chain_element(const irr::channels_t& cg,const irr::template_t& t,size_t target,
	const std::vector<fec_irrep_axis_block_t>& condition,const irr::copies& old,const irr::copies& previous,
	const std::vector<fs::path>& chain_files,rref_option_t opt,const fs::path& cache_dir){
	const auto L=layout(cg,t,target,condition,old,previous);
	irr::matrices K;                       // loaded per prime from chain_files
	bool K_loaded=false;
	irr::copies oldbase{};size_t groups=0;
	for(size_t b=0;b<6;++b){oldbase[b]=groups;groups+=old[b];}
	std::vector<std::vector<std::pair<size_t,size_t>>> decode(6);
	std::vector<std::unique_ptr<irr_parallel::grouped_entries>> index(6);
	{
		K=load_chain(chain_files,&opt->pool);K_loaded=true;
		for(size_t b=0;b<6;++b){
			std::vector<size_t> group;
			for(size_t left=0;left<6;++left)for(size_t a=0;a<old[left];++a)
				for(size_t g:cg.channels[b][left]){decode[b].push_back({g,a});group.push_back(oldbase[left]+a);}
			if(K[b].ncol!=decode[b].size()||K[b].nrow!=previous[b])
				throw std::runtime_error("Copy/channel layout mismatch");
			index[b]=std::make_unique<irr_parallel::grouped_entries>(K[b],group,groups);
		}
	}
	if(!L.nr||!L.nc){ // reference edge path over the (empty-dimension) QQ assembly
		auto matrix=irr::assemble(cg,t,target,condition,old,previous,K,opt);
		const size_t nnz=matrix.nnz();
		auto kernel=irr::solve(matrix,opt);
		return {std::move(kernel),L.nr,L.nc,nnz};
	}

	uint64_t hash=14695981039346656037ULL;
	auto feed=[&](const std::string& x){for(unsigned char ch:x){hash^=ch;hash*=1099511628217ULL;}};
	feed("nr="+std::to_string(L.nr)+",nc="+std::to_string(L.nc)+";");
	for(const auto& file:chain_files)feed(hex_sha256_file(irr::tensor_path(file))+";");
	for(size_t b=0;b<6;++b)for(size_t g=0;g<t.routes[b].size();++g)
		for(const auto& route:t.routes[b][g])
			feed(std::to_string(b)+":"+std::to_string(g)+":"+std::to_string(route.row)+":"
				+std::to_string(route.channel)+"="+route.coefficient.get_str()+";");
	for(const auto& scale:t.row_scale)feed("s="+scale.get_str()+";");

	fs::create_directories(cache_dir);
	const std::string identity="FECSTREAM2 method=1 backward=true input="+std::to_string(hash);
	const auto identity_path=cache_dir/"input.txt";
	if(fs::exists(identity_path)){
		std::ifstream f(identity_path);std::string old_line;std::getline(f,old_line);
		if(old_line!=identity)throw std::runtime_error("CRT input identity mismatch");
	}else{
		std::ofstream f(identity_path);f<<identity<<'\n';f.flush();
		if(!f)throw std::runtime_error("CRT identity write failed");
	}

	const size_t kheight=chain_height_bits(K,&opt->pool);
	size_t cheight=1,sheight=1;
	for(size_t b=0;b<6;++b)for(size_t g=0;g<t.routes[b].size();++g)
		for(const auto& route:t.routes[b][g])cheight=std::max(cheight,size_t(route.coefficient.height_bits()));
	for(const auto& scale:t.row_scale)sheight=std::max(sheight,size_t(scale.height_bits()));
	const int_t input_height=int_t(1)<<(kheight+cheight+sheight+2);

	Pivots pivots;
	std::vector<std::pair<ulong,fs::path>> files;
	int_t modulus=1;
	std::vector<char> pivot_rows;          // natural row indices; empty until first prime
	size_t assembled_nnz=0;
	ulong prime=n_nextprime(1ULL<<60,0);
	for(size_t iteration=1;;++iteration,prime=n_nextprime(prime,0)){
		(void)iteration;
		fs::path path;
		const std::string prefix="prime_"+std::to_string(prime)+"_";
		for(auto& entry:fs::directory_iterator(cache_dir))
			if(entry.path().filename().string().starts_with(prefix)&&entry.path().extension()==".bin"
				&&fs::exists(entry.path().string()+".ready")){path=entry.path();break;}
		if(path.empty()){
			const field_t F(FIELD_Fp,prime);
			if(!K_loaded){K=load_chain(chain_files,&opt->pool);K_loaded=true;}
			Matrix modular=assemble_modular(cg,t,target,condition,decode,index,K,L,pivot_rows,opt,F);
			if(!assembled_nnz)assembled_nnz=modular.nnz();
			unload_chain(K);K_loaded=false;   // the reduction never touches K
			if(files.empty())pivots=sparse_mat_rref_forward(modular,F,opt);
			else sparse_mat_direct_rref(modular,pivots,F,opt);
			triangular_solver_2(modular,pivots,F,opt);
			size_t attempt=1;fs::path candidate;
			do{candidate=cache_dir/(prefix+std::to_string(attempt++)+".bin");}while(fs::exists(candidate));
			k_crt::save(candidate,modular,prime,pivots);
			modular.clear();
			std::cout<<"CRT_CACHE_SAVED "<<candidate<<std::endl;
			path=candidate;
		}else{
			k_crt::Stream stream(path,false);
			auto cached=k_crt::header(stream,prime,L.nr,L.nc);
			if(files.empty())pivots=std::move(cached);
			else if(!k_crt::same(cached,pivots))throw std::runtime_error("CRT pivot mismatch");
			std::cout<<"CRT_CACHE_REUSED "<<path<<std::endl;
		}
		files.push_back({prime,path});
		modulus*=prime;
		if(files.size()==1){
			pivot_rows.assign(L.nr,0);
			for(const auto& group:pivots)for(auto [r,c]:group)pivot_rows[r]=1;
			size_t rank=0;for(const auto& group:pivots)rank+=group.size();
			if(rank==L.nc){ // empty kernel, as in the reference
				auto next=sparse_mat<rat_t,int32_t>(0,L.nc);
				std::ofstream done(cache_dir/"reconstruction.txt");
				done<<"primes "<<files.size()<<"\nmodulus_bits "<<modulus.bits()<<"\nkernel_rows 0\n";
				std::ofstream used(cache_dir/"used_cache.tsv");
				for(auto& [p,f]:files){
					used<<f.filename().string()<<'\t'<<hex_sha256_file(f)<<'\n';
					used<<f.filename().string()<<".ready\t"<<hex_sha256_file(f.string()+".ready")<<'\n';
				}
				std::cout<<"** Reconstruct success! Using mod ~ 2^"<<modulus.bits()<<" kernel_rows=0"<<std::endl;
				return {std::move(next),L.nr,L.nc,assembled_nnz};
			}
		}
		sparse_mat<rat_t,int32_t> next;
		size_t max_2nd_bits=0;
		if(reconstruct_kernel(files,pivots,L.nr,L.nc,input_height,modulus,opt,next,max_2nd_bits)){
			next.sort_rows_by_nnz();
			std::ofstream done(cache_dir/"reconstruction.txt");
			done<<"primes "<<files.size()<<"\nmodulus_bits "<<modulus.bits()
				<<"\nmax_individual_2nd_bits "<<max_2nd_bits<<"\n";
			std::ofstream used(cache_dir/"used_cache.tsv");
			for(auto& [p,f]:files){
				used<<f.filename().string()<<'\t'<<hex_sha256_file(f)<<'\n';
				used<<f.filename().string()<<".ready\t"<<hex_sha256_file(f.string()+".ready")<<'\n';
			}
			std::cout<<"** Reconstruct success! Using mod ~ 2^"<<modulus.bits()
				<<" max_2nd_bits="<<max_2nd_bits<<std::endl;
			return {std::move(next),L.nr,L.nc,assembled_nnz};
		}
		std::cout<<">> Reconstruct failed, now mod ~ 2^"<<modulus.bits()<<std::endl;
	}
}

} // namespace fec_stream

int main(int argc,char** argv){try {
	if(argc<3)throw std::runtime_error("irr_fec_stream MAX_WEIGHT OUTPUT [--condition WXF] [--condition-blocks TSV] [--clean-block-cache]");
	size_t maximum=std::stoul(argv[1]);if(!maximum)throw std::runtime_error("MAX_WEIGHT must be positive");
	fs::path root=argv[2],condition="data/dlogmat_ca_1term_pd3.wxf",blocks="data/dlogmat_ca_1term_pd3_blocks.tsv";
	bool clean_blocks=false;
	for(int i=3;i<argc;++i){std::string key=argv[i];if(key=="--clean-block-cache"){clean_blocks=true;continue;}if(i+1==argc)throw std::runtime_error("Missing option value");
		if(key=="--condition")condition=argv[++i];else if(key=="--condition-blocks")blocks=argv[++i];else throw std::runtime_error("Unknown option "+key);}
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=true;irr_run::threads(opt);
	if(clean_blocks)irr_stage::no_links(root);
	irr_run::bind(root,"FEC primitive method1",{condition,blocks});
	if(clean_blocks)irr_stage::no_links(root/".active.lock");
	irr_stage::lock active(root/".active.lock");
	auto* pool=&opt->pool;irr::channels_t cg(pool);
	std::ofstream terms(root/"CG.tsv");terms<<"target\tinput\tchannel\ttarget_component\tinput_component\tletter\tcoefficient\n";
	for(size_t b=0;b<6;++b)for(size_t a=0;a<6;++a)for(size_t g:cg.channels[b][a])
		for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)for(const auto& term:cg.cg[b].vectors[g].component[s])
			terms<<irr::specs()[b].name<<'\t'<<irr::specs()[a].name<<'\t'<<g<<'\t'<<s<<'\t'<<term.coordinate/9-cg.axis[a].offset<<'\t'<<term.coordinate%9<<'\t'<<term.coefficient<<'\n';
	terms.close();
	auto W=irr::read_tensor(condition,pool);auto cb=nmhv_sew_read_irrep_axis(blocks);
	std::array<irr::template_t,6> templates;for(size_t b=0;b<6;++b)templates[b]=irr::make_template(cg,b,W,cb,opt);
	irr::copies older{},previous{1,0,0,0,0,0};
	std::ofstream summary(root/"summary.tsv");summary<<"weight\tdimension\tirr_nnz\theight_bits\n";
	for(size_t w=1;w<=maximum;++w) {
		auto dir=root/joint_eet_weight_name(w);fs::create_directories(dir);irr::matrices next;bool legacy_verified=false;
		std::array<size_t,6> block_nnz{},block_height{};
		std::vector<fs::path> chain_files;   // logical .wxf paths; read_tensor falls back to .shape
		if(w>=2){auto previous_dir=root/joint_eet_weight_name(w-1);
			for(size_t b=0;b<6;++b)
				chain_files.push_back(previous_dir/("irrFEC_"+std::to_string(w-1)+"_"+irr::specs()[b].name+".wxf"));}
		for(size_t b=0;b<6;++b) {
			auto file=dir/("irrFEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf");
			if(clean_blocks){irr_stage::no_links(file);irr_stage::no_links(fs::path(file.string()+".ready"));irr_stage::no_links(dir/"cache"/irr::specs()[b].name);}
			if(irr_run::ready(file)) {
				next[b]=joint_eet_flatten_tensor(irr::read_tensor(file,pool),pool);
				if(clean_blocks)irr_block_cache::clean(dir,irr::specs()[b].name,file,active,legacy_verified);
				std::cout<<"REUSE weight="<<w<<" block="<<irr::specs()[b].name<<std::endl;continue;
			}
			size_t nr=0,nc=0,nnz=0;
			if(w==1) {
				auto seed=irr::skeleton(previous,pool);auto c=make_hhex_candidates(fec_symmetry_t::parity_dihedral,seed,irr::specs()[b].name);
				next[b]=sparse_mat<rat_t,int32_t>(b==0||b==2?1:0,c.size());if(next[b].nrow)next[b][0].push_back(0,rat_t(1));
			}else {
				auto result=fec_stream::chain_element(cg,templates[b],b,cb,older,previous,chain_files,opt,dir/"cache"/irr::specs()[b].name);
				nr=result.nr;nc=result.nc;nnz=result.assembled_nnz;
				next[b]=std::move(result.next);
				std::cout<<"IRR_FEC_STREAM weight="<<w<<" block="<<irr::specs()[b].name
					<<" rows="<<nr<<" cols="<<nc<<" assembly_nnz="<<nnz
					<<" kernel_rows="<<next[b].nrow<<" kernel_nnz="<<next[b].nnz()<<std::endl;
			}
			// The WXF serialization buffers the whole tensor, so the QQ rows
			// are released right after the tensor copy; only the captured
			// statistics and dimensions are needed afterwards.
			const size_t kept_rows=next[b].nrow,kept_cols=next[b].ncol;
			block_nnz[b]=next[b].nnz();
			auto stats=joint_eet_rational_stats(next[b]);
			block_height[b]=stats.height_bits;
			const auto peak_bytes=hex_platform::current_process_memory().peak_working_set_bytes;
			{
				sparse_tensor<rat_t,int32_t,SPARSE_CSR> tensor(next[b],pool);
				next[b]=sparse_mat<rat_t,int32_t>(kept_rows,kept_cols);
				irr::write_tensor(file,tensor,false,clean_blocks);
			}
			std::ofstream report(dir/(irr::specs()[b].name+".tsv"));
			report<<"rows\tcolumns\tmatrix_nnz\tcopies\tirr_nnz\theight_bits\tpeak_bytes\n"<<nr<<'\t'<<nc<<'\t'<<nnz<<'\t'<<kept_rows<<'\t'<<block_nnz[b]<<'\t'<<block_height[b]<<'\t'<<peak_bytes<<'\n';
			report.close();if(!report)throw std::runtime_error("Cannot write block report");
			irr_run::mark(file);
			if(clean_blocks)irr_block_cache::clean(dir,irr::specs()[b].name,file,active,legacy_verified);
		}
		irr::copies current{};size_t dim=0,nnz=0,height=0;std::ofstream mult(dir/"multiplicities.tsv"),hashes(dir/"checksums.sha256"),cols(dir/"columns.tsv");
		cols<<"target\tinput\told_copies\tcolumn_offset\tchannels_per_copy\n";
		for(size_t b=0;b<6;++b) {
			current[b]=next[b].nrow;dim+=current[b]*irr::specs()[b].irrep_dimension;nnz+=block_nnz[b];height=std::max(height,block_height[b]);
			mult<<irr::specs()[b].name<<'\t'<<current[b]<<'\n';size_t offset=0;
			for(size_t a=0;a<6;++a){cols<<irr::specs()[b].name<<'\t'<<irr::specs()[a].name<<'\t'<<previous[a]<<'\t'<<offset<<'\t'<<cg.channels[b][a].size()<<'\n';offset+=previous[a]*cg.channels[b][a].size();}
			if(next[b].ncol!=offset)throw std::runtime_error("FEC channel dimensions mismatch");
			auto file=irr::tensor_path(dir/("irrFEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"));hashes<<hex_sha256_file(file)<<"  "<<file.filename().string()<<'\n';
		}
		mult.close();cols.close();hashes<<hex_sha256_file(dir/"multiplicities.tsv")<<"  multiplicities.tsv\n"<<hex_sha256_file(dir/"columns.tsv")<<"  columns.tsv\n";hashes.close();
		std::ofstream manifest(dir/"manifest.json");manifest<<"{\"schema\":\"hexagon-irr-fec-v1\",\"complete\":true,\"weight\":"<<w<<",\"dimension\":"<<dim<<",\"irr_nnz\":"<<nnz<<",\"height_bits\":"<<height<<",\"method\":1,\"checksums_sha256\":\""<<hex_sha256_file(dir/"checksums.sha256")<<"\"}\n";
		summary<<w<<'\t'<<dim<<'\t'<<nnz<<'\t'<<height<<std::endl;
		manifest.close();irr_stage::complete_chain(dir);
		older=previous;previous=current;std::cout<<"COMPLETE FEC weight="<<w<<" dimension="<<dim<<" nnz="<<nnz<<" height_bits="<<height<<std::endl;
	}
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}}
