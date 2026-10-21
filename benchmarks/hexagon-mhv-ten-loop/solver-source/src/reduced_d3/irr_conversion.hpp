#pragma once
#include "irr_recurrence.hpp"
#include "joint_eet_io.hpp"
#include "nmhv_irrep_sew.hpp"
#include "hex_hash.hpp"
#include "irr_profile.hpp"
#include <set>
#include <map>

namespace irr {
inline void conversion_progress(const std::string& stage,bool detail=false) {
	static const auto start=std::chrono::steady_clock::now();
	static auto last=start;
	static const bool detailed=[](){auto v=std::getenv("SHIFT_VERBOSE_DETAIL");return v&&std::string(v)=="1";}();
	auto now=std::chrono::steady_clock::now();
	if(detail&&!detailed&&now-last<std::chrono::seconds(30))return;
	last=now;
	auto memory=hex_platform::current_process_memory();
	std::cout<<"CONVERT "<<stage<<" elapsed_seconds="
		<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
		<<" memory_valid="<<memory.valid<<" rss_bytes="<<memory.working_set_bytes
		<<" process_peak_rss_bytes="<<memory.peak_working_set_bytes<<std::endl;
}
inline void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
inline copies chain_copies(const std::filesystem::path& root,size_t w,bool first) {
	copies m{};if(!w&&first){m[0]=1;return m;}
	if(first) {
		std::ifstream in(root/joint_eet_weight_name(w)/"multiplicities.tsv");
		for(size_t b=0;b<6;++b){std::string name;in>>name>>m[b];require(bool(in)&&name==specs()[b].name,"Invalid irrFEC metadata");}
	}else {
		auto a=nmhv_sew_read_irrep_axis(w?root/joint_eet_weight_name(w)/"blocks.tsv":root/"LEC_0_blocks.tsv");
		require(a.size()==6,"Expected six irrLEC blocks");
		for(size_t b=0;b<6;++b){require(a[b].spec.name==specs()[b].name,"Invalid irrLEC block order");m[b]=a[b].multiplicity;}
	}return m;
}
inline matrices read_recurrence(const std::filesystem::path& root,size_t w,bool first,const channels_t& cg,thread_pool* pool,bool profile=false) {
	const field_t F(FIELD_QQ);matrices K;auto old=chain_copies(root,w-1,first),m=chain_copies(root,w,first);
	double read_seconds=0,csr_seconds=0;
	for(size_t b=0;b<6;++b) {
		auto file=root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+specs()[b].name+".wxf");
		auto start=std::chrono::steady_clock::now();
		auto t=read_tensor(file,pool);require(t.rank()==2,"Recurrence must be a matrix");
		read_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
		K[b]=joint_eet_flatten_tensor(std::move(t),pool);csr_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
		size_t columns=0;for(size_t a=0;a<6;++a)columns+=old[a]*cg.channels[b][a].size();
		require(K[b].nrow==m[b]&&K[b].ncol==columns,"Recurrence copy/channel dimensions differ");
	}if(profile)std::cout<<"RECURRENCE_IO first="<<first<<" weight="<<w<<" read_decode_seconds="<<read_seconds<<" csr_seconds="<<csr_seconds<<std::endl;
	return K;
}

// A sparse exact right inverse, preferring private columns of kernel rows.
// These are normally diagonal charts even after primitive row normalization.
inline sparse_mat<rat_t,int32_t> right_inverse(const sparse_mat<rat_t,int32_t>& K,rref_option_t opt) {
	const field_t F(FIELD_QQ);sparse_mat<rat_t,int32_t> inverse(K.ncol,K.nrow);
	if(!K.nrow)return inverse;
	std::vector<size_t> count(K.ncol,0),pivot(K.nrow,SIZE_MAX);
	for(const auto& row:K.rows)for(auto [j,x]:row)if(x!=0)++count[j];
	bool monomial=true;
	for(size_t i=0;i<K.nrow;++i) {
		for(auto [j,x]:K[i])if(x!=0&&count[j]==1){pivot[i]=j;inverse[j].push_back(int32_t(i),rat_t(1)/x);break;}
		if(pivot[i]==SIZE_MAX){monomial=false;break;}
	}
	if(!monomial) {
		pivot=dlogmat_irrep_pivot_columns(K,opt);require(pivot.size()==K.nrow,"Target recurrence rows are dependent");
		auto chart=sparse_mat_inverse(dlogmat_irrep_select_columns(K,pivot),F,opt);
		inverse=sparse_mat<rat_t,int32_t>(K.ncol,K.nrow);
		for(size_t i=0;i<pivot.size();++i)inverse[pivot[i]]=std::move(chart[i]);
	}
	require(nmhv_sew_matrices_equal(sparse_mat_mul(K,inverse,F,&opt->pool),joint_eet_identity<rat_t,int32_t>(K.nrow)),"Invalid exact recurrence chart");
	return inverse;
}

struct recouple_t { size_t block,channel; rat_t value; };
struct seam_recoupling_t {
	std::array<std::vector<std::vector<recouple_t>>,6> routes;
	seam_recoupling_t(const channels_t& cg,rref_option_t opt) {
		const field_t F(FIELD_QQ);
		for(size_t b=0;b<6;++b) {
			routes[b].resize(cg.cg[b].size());
			for(size_t a=0;a<6;++a) {
				const auto& left=cg.channels[b][a];const auto& right=cg.channels[a][b];
				size_t dim=specs()[a].irrep_dimension*9*specs()[b].irrep_dimension;
				sparse_mat<rat_t,int32_t> A(left.size(),dim),B(right.size(),dim);
				for(size_t c=0;c<left.size();++c)for(size_t s=0;s<specs()[b].irrep_dimension;++s)
					for(const auto& term:cg.cg[b].vectors[left[c]].component[s]) {
						size_t t=term.coordinate/9-cg.axis[a].offset,j=term.coordinate%9;
						A[c].push_back(int32_t((t*9+j)*specs()[b].irrep_dimension+s),term.coefficient*(s?rat_t(3):rat_t(1)));
					}
				for(size_t c=0;c<right.size();++c)for(size_t t=0;t<specs()[a].irrep_dimension;++t)
					for(const auto& term:cg.cg[a].vectors[right[c]].component[t]) {
						size_t s=term.coordinate/9-cg.axis[b].offset,j=term.coordinate%9;
						B[c].push_back(int32_t((t*9+j)*specs()[b].irrep_dimension+s),term.coefficient*(t?rat_t(3):rat_t(1)));
					}
				for(auto& row:A.rows)fec_irrep_merge_entries(row);for(auto& row:B.rows)fec_irrep_merge_entries(row);
				auto R=sparse_mat_mul(A,right_inverse(B,opt),F,&opt->pool);
				require(nmhv_sew_matrices_equal(sparse_mat_mul(R,B,F,&opt->pool),A),"Fixed CG recoupling failed");
				for(size_t c=0;c<left.size();++c)for(auto [d,x]:R[c])routes[b][left[c]].push_back({a,size_t(d),x});
			}
		}
	}
};

// Internally an affine origin is a final computational row, never a basis
// vector. Only the declared direction rows carry a public solution index.
struct sewing_object_t {
	matrices S;
	size_t h=0,k=0,directions=0;
	std::string kind="homogeneous";
	copies first{},last{};
	size_t rows()const{return directions+(kind=="homogeneous"?0:1);}
};

inline sewing_object_t shift_one(const sewing_object_t& input,bool left,const std::filesystem::path& fec,
	const std::filesystem::path& lec,const channels_t& cg,const seam_recoupling_t& recouple,rref_option_t opt) {
	irr_parallel::phase timing("shift");
	const field_t F(FIELD_QQ);auto* pool=&opt->pool;sewing_object_t out;
	out.h=left?input.h-1:input.h+1;out.k=left?input.k+1:input.k-1;
	require(out.h&&out.k,"Moving through weight zero is not supported");out.kind=input.kind;out.directions=input.directions;
	out.first=chain_copies(fec,out.h,true);out.last=chain_copies(lec,out.k,false);
	conversion_progress("phase=read_moving begin h="+std::to_string(input.h)+" k="+std::to_string(input.k));
	auto moving=read_recurrence(left?fec:lec,left?input.h:input.k,left,cg,pool);
	conversion_progress("phase=read_target begin");
	auto target=read_recurrence(left?lec:fec,left?out.k:out.h,!left,cg,pool);
	conversion_progress("phase=build_index begin");
	const auto& old=left?out.first:out.last;const auto& stationary=left?input.last:input.first;
	std::array<std::vector<std::tuple<size_t,size_t,size_t>>,6> decode;
	copies oldbase{};size_t groups=0;for(size_t a=0;a<6;++a){oldbase[a]=groups;groups+=old[a];}
	for(size_t b=0;b<6;++b)for(size_t a=0;a<6;++a)for(size_t copy=0;copy<old[a];++copy)
		for(size_t g:cg.channels[b][a])decode[b].push_back({a,copy,g});
	std::array<std::unique_ptr<irr_parallel::grouped_entries>,6> index;
	for(size_t b=0;b<6;++b){std::vector<size_t> group;
		for(auto [a,copy,g]:decode[b])group.push_back(oldbase[a]+copy);
		index[b]=std::make_unique<irr_parallel::grouped_entries>(moving[b],group,groups);
	}
	matrices inverse;std::array<copies,6> base{};copies columns{};
	for(size_t a=0;a<6;++a) {
		conversion_progress("phase=right_inverse begin block="+specs()[a].name);
		require(out.first[a]*out.last[a]<=size_t(INT32_MAX),"Sewing block exceeds int32 range");
		out.S[a]=sparse_mat<rat_t,int32_t>(input.rows(),out.first[a]*out.last[a]);
		inverse[a]=right_inverse(target[a],opt);
		for(size_t b=0;b<6;++b){base[a][b]=columns[a];columns[a]+=stationary[b]*cg.channels[a][b].size();}
		require(columns[a]==target[a].ncol,"Target channel ordering mismatch");
	}
	// One solution and one ambient irrep at a time. Workers own destination
	// rows and immediately merge them; neither rationals nor matrices are cloned
	// per thread. The source index stores only positions in the input solution.
	for(size_t sol=0;sol<input.rows();++sol) {
		const auto progress=" solution="+std::to_string(sol+1)+"/"+std::to_string(input.rows());
		conversion_progress("phase=source_index begin"+progress,true);
		std::array<std::vector<std::vector<std::pair<size_t,size_t>>>,6> source;
		for(size_t b=0;b<6;++b){source[b].resize(moving[b].nrow);
			for(size_t p=0;p<input.S[b][sol].nnz();++p){size_t ij=input.S[b][sol](p),i=ij/input.last[b],j=ij%input.last[b];
				source[b][left?i:j].push_back({left?j:i,p});
			}
		}
		for(size_t a=0;a<6;++a) {
			const auto context=progress+" block="+specs()[a].name;
			conversion_progress("phase=accumulate begin"+context,true);
			sparse_mat<rat_t,int32_t> ambient(old[a],columns[a]);
			irr_parallel::rows(old[a],pool,[&](size_t copy){
				// Accumulate by destination, not by multiplication path: many
				// recurrence terms contribute to the same small tail coordinate.
				std::map<int32_t,rat_t> sums;
				for(size_t b=0;b<6;++b){const auto& idx=*index[b];size_t group=oldbase[a]+copy;
					for(size_t p=idx.offsets[group];p<idx.offsets[group+1];++p){auto [r,pos]=idx.entries[p];
						auto [block,c,g]=decode[b][moving[b][r](pos)];const auto& y=moving[b][r][pos];
						for(auto [stationary_copy,s]:source[b][r])for(const auto& route:recouple.routes[b][g]) {
							require(route.block==a,"Recoupling changes old irrep");
							sums[int32_t(base[a][b]+stationary_copy*cg.channels[a][b].size()+route.channel)]+=input.S[b][sol][s]*y*route.value;
						}
					}
				}
				for(auto& [column,value]:sums)if(value!=0)
					ambient[copy].push_back(column,std::move(value));
			});
			conversion_progress("phase=target_coordinates begin"+context+" ambient_nnz="+std::to_string(ambient.nnz()),true);
			auto coefficients=sparse_mat_mul(ambient,inverse[a],F,pool);
			conversion_progress("phase=verify begin"+context+" coefficients_nnz="+std::to_string(coefficients.nnz()),true);
			require(nmhv_sew_matrices_equal(sparse_mat_mul(coefficients,target[a],F,pool),ambient),"Target recurrence does not span the shifted tensor");
			conversion_progress("phase=store_row begin"+context,true);
			std::vector<size_t> offsets(coefficients.nrow+1);
			for(size_t r=0;r<coefficients.nrow;++r)offsets[r+1]=offsets[r]+coefficients[r].nnz();
			auto& row=out.S[a][sol];row.reserve(offsets.back());row.resize(offsets.back());
			irr_parallel::rows(coefficients.nrow,pool,[&](size_t copy){for(size_t p=0;p<coefficients[copy].nnz();++p){
				size_t j=coefficients[copy](p),dst=offsets[copy]+p;
				row(dst)=int32_t(left?copy*out.last[a]+j:j*out.last[a]+copy);row[dst]=std::move(coefficients[copy][p]);
			}});
			if(!left)row.sort_indices();
			conversion_progress("phase=store_row end"+context,true);
		}
		if(sol==0||(sol+1)%std::max(size_t(1),(input.rows()+99)/100)==0||sol+1==input.rows())
			conversion_progress("phase=solution_complete"+progress);
	}
	return out;
}
}
