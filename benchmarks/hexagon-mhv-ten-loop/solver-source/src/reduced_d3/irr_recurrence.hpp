#pragma once
#include "fec_irrep_reduced.hpp"
#include "k_crt_cache.hpp"
#include "irr_parallel.hpp"
#include <array>
#include <tuple>

// Multiplicity recurrence in the existing rational parity x D3 convention.
// LEC uses the same coordinate CG, but its new letter is physically on the left.
namespace irr {
// Zero-dimensional blocks are metadata, not Mathematica SparseArrays. Keep
// their shape separately; nonzero shapes (even nnz=0) remain ordinary WXF.
inline std::filesystem::path tensor_path(const std::filesystem::path& path) {
	return std::filesystem::exists(path)?path:std::filesystem::path(path.string()+".shape");
}
inline sparse_tensor<rat_t,int32_t,SPARSE_CSR> read_tensor(const std::filesystem::path& path,thread_pool* pool) {
	if(std::filesystem::exists(path))return sparse_tensor_read_wxf<rat_t,int32_t>(path,field_t(FIELD_QQ),pool);
	std::ifstream in(path.string()+".shape");size_t rank=0;in>>rank;
	if(!in||rank<2||rank>8)throw std::runtime_error("Missing tensor or invalid empty shape: "+path.string());
	std::vector<size_t> dims(rank);for(auto& d:dims)in>>d;
	if(!in||std::find(dims.begin(),dims.end(),0)==dims.end())throw std::runtime_error("Invalid empty tensor shape: "+path.string());
	return sparse_tensor<rat_t,int32_t,SPARSE_CSR>(dims);
}
inline std::filesystem::path write_tensor(const std::filesystem::path& path,
	const sparse_tensor<rat_t,int32_t,SPARSE_CSR>& tensor,bool mma_layout=false,bool verify_bytes=false) {
	if(std::filesystem::exists(path)||std::filesystem::exists(path.string()+".shape"))throw std::runtime_error("Refusing to overwrite "+path.string());
	if(std::find(tensor.dims().begin(),tensor.dims().end(),0)!=tensor.dims().end()) {
		if(tensor.nnz())throw std::runtime_error("Nonzero entries in zero-dimensional tensor");
		std::ofstream out(path.string()+".shape");out<<tensor.rank();for(auto d:tensor.dims())out<<'\t'<<d;out<<'\n';out.close();
		if(!out)throw std::runtime_error("Cannot write tensor shape");return path.string()+".shape";
	}
	auto bytes=sparse_tensor_write_wxf(tensor,true,mma_layout);
	if(bytes.empty())throw std::runtime_error("Empty WXF serialization: "+path.string());
	std::ofstream out(path,std::ios::binary);out.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());out.close();
	if(!out)throw std::runtime_error("Cannot write WXF: "+path.string());
	if(verify_bytes){
		// Compare the serialized tensor with a bounded readback, not a second tensor.
		std::ifstream in(path,std::ios::binary);std::array<char,65536> buffer;
		for(size_t offset=0;offset<bytes.size();){auto n=std::min(buffer.size(),bytes.size()-offset);in.read(buffer.data(),n);
			if(!in||std::memcmp(buffer.data(),bytes.data()+offset,n))throw std::runtime_error("WXF readback mismatch: "+path.string());offset+=n;}
		if(in.peek()!=EOF)throw std::runtime_error("Trailing WXF data: "+path.string());
	}
	return path;
}
inline const auto& specs() { return hhex_block_specs(fec_symmetry_t::parity_dihedral); }
using matrices = std::array<sparse_mat<rat_t, int32_t>, 6>;
using copies = std::array<size_t, 6>;

inline std::vector<fec_block_t> skeleton(const copies& m, thread_pool* pool) {
	std::vector<fec_block_t> blocks;
	for (size_t b = 0; b < 6; ++b) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> t({m[b]*specs()[b].irrep_dimension,1,9});
		blocks.push_back({specs()[b],m[b],sparse_tensor<rat_t,int32_t,SPARSE_CSR>(std::move(t),pool)});
	}
	return blocks;
}

inline size_t axis_block(const std::vector<fec_irrep_axis_block_t>& axis, size_t i) {
	for (size_t b = 0; b < axis.size(); ++b)
		if (i >= axis[b].offset && i < axis[b].offset+axis[b].dimension) return b;
	throw std::runtime_error("Coordinate outside irrep axis");
}

struct channels_t {
	std::vector<fec_block_t> prototype;
	std::vector<fec_irrep_axis_block_t> axis;
	std::array<fec_candidate_basis_t,6> cg;
	std::array<std::array<std::vector<size_t>,6>,6> channels;
	explicit channels_t(thread_pool* pool) {
		prototype=skeleton({1,1,1,1,1,1},pool); axis=fec_irrep_axis_from_blocks(prototype);
		for (size_t b=0;b<6;++b) {
			cg[b]=make_hhex_candidates(fec_symmetry_t::parity_dihedral,prototype,specs()[b].name);
			for (size_t c=0;c<cg[b].size();++c) {
				auto input=axis_block(axis,cg[b].vectors[c].component[0][0].coordinate/9);
				channels[b][input].push_back(c);
			}
		}
	}
};
struct route_t { size_t row, channel; rat_t coefficient; };
struct row_t { size_t left, condition, copy; };
struct template_t {
	std::vector<row_t> rows;
	std::array<std::vector<std::vector<route_t>>,6> routes;
	std::vector<rat_t> row_scale;
};

inline template_t make_template(const channels_t& cg, size_t target,
	const sparse_tensor<rat_t,int32_t,SPARSE_CSR>& W,
	const std::vector<fec_irrep_axis_block_t>& condition, rref_option_t opt) {
	const field_t F(FIELD_QQ); template_t t;
	if (W.rank()!=3 || W.dim(0)!=9 || W.dim(1)!=9 || fec_irrep_axis_dimension(condition)!=W.dim(2))
		throw std::runtime_error("Condition shape/metadata mismatch");
	auto codomain=make_hhex_codomain_basis(cg.axis,condition,specs()[target],F,opt);
	if (!codomain.monomial_chart) throw std::runtime_error("Direct recurrence requires a monomial codomain chart");
	std::vector<int32_t> selected(codomain.raw_dimension,-1);
	for (size_t r=0;r<codomain.multiplicity;++r) {
		size_t raw=codomain.pivot_coordinates[r],q=raw%W.dim(2),b=axis_block(condition,q);
		selected[raw]=int32_t(r);
		t.rows.push_back({axis_block(cg.axis,raw/W.dim(2)),b,(q-condition[b].offset)%condition[b].multiplicity});
	}
	std::array<std::array<std::vector<std::pair<size_t,rat_t>>,9>,9> lookup;
	for(size_t i=0;i<9;++i)for(size_t p=W.rowptr()[i];p<W.rowptr()[i+1];++p)
		lookup[i][W.index(p)[0]].push_back({size_t(W.index(p)[1]),W.val(p)});
	for(size_t b=0;b<6;++b) {
		t.routes[b].resize(cg.cg[b].size());
		for(size_t g=0;g<cg.cg[b].size();++g) {
			std::map<std::pair<size_t,size_t>,rat_t> sum;
			for(size_t h=0;h<cg.channels[target][b].size();++h) {
				auto c=cg.channels[target][b][h];
				for(const auto& term:cg.cg[target].vectors[c].component[0]) {
					size_t component=term.coordinate/9-cg.axis[b].offset,j=term.coordinate%9;
					for(const auto& prev:cg.cg[b].vectors[g].component[component])
						for(const auto& [q,value]:lookup[prev.coordinate%9][j]) {
							int32_t r=selected[(prev.coordinate/9)*W.dim(2)+q];
							if(r>=0)sum[{size_t(r),h}]+=prev.coefficient*term.coefficient*value*codomain.transpose_inverse_chart[r][0];
						}
				}
			}
			for(const auto& [key,value]:sum)if(value!=0)t.routes[b][g].push_back({key.first,key.second,value});
		}
	}
	// One scale per COMPLETE equation, shared by every b/g/h storage fragment.
	t.row_scale.assign(t.rows.size(),rat_t(1));
	std::vector<int_t> den(t.rows.size(),int_t(1)),gcd(t.rows.size(),int_t(0));
	std::vector<int> sign(t.rows.size(),0);
	for(const auto& block:t.routes)for(const auto& channel:block)for(const auto& r:channel)
		fmpz_lcm(den[r.row].data(),den[r.row].data(),fmpq_denref(r.coefficient.data()));
	for(const auto& block:t.routes)for(const auto& channel:block)for(const auto& r:channel) {
		int_t value; fmpz_divexact(value.data(),den[r.row].data(),fmpq_denref(r.coefficient.data()));
		fmpz_mul(value.data(),value.data(),fmpq_numref(r.coefficient.data()));
		if(!sign[r.row])sign[r.row]=fmpz_sgn(value.data());
		fmpz_gcd(gcd[r.row].data(),gcd[r.row].data(),value.data());
	}
	for(size_t r=0;r<t.rows.size();++r)if(sign[r]) {
		fmpq_set_fmpz_frac(t.row_scale[r].data(),den[r].data(),gcd[r].data());
		if(sign[r]<0)t.row_scale[r]=-t.row_scale[r];
	}
	return t;
}

inline sparse_mat<rat_t,int32_t> assemble(const channels_t& cg,const template_t& t,size_t target,
	const std::vector<fec_irrep_axis_block_t>& condition,const copies& old,const copies& previous,
	const matrices& K,rref_option_t opt,bool primitive=true) {
	irr_parallel::phase timing("recurrence_assembly");
	std::vector<size_t> rowbase(t.rows.size());size_t nr=0,nc=0;copies colbase{};
	for(size_t r=0;r<t.rows.size();) {
		size_t end=r+1;while(end<t.rows.size()&&t.rows[end].left==t.rows[r].left&&t.rows[end].condition==t.rows[r].condition)++end;
		for(size_t s=r;s<end;++s)rowbase[s]=nr+t.rows[s].copy;
		nr+=old[t.rows[r].left]*condition[t.rows[r].condition].multiplicity;r=end;
	}
	for(size_t b=0;b<6;++b){colbase[b]=nc;nc+=previous[b]*cg.channels[target][b].size();}
	if(nc>size_t(INT32_MAX))throw std::runtime_error("Candidate columns exceed int32");
	sparse_mat<rat_t,int32_t> matrix(nr,nc);
	copies oldbase{};size_t groups=0;for(size_t b=0;b<6;++b){oldbase[b]=groups;groups+=old[b];}
	for(size_t b=0;b<6;++b) {
		std::vector<std::pair<size_t,size_t>> decode;
		std::vector<size_t> group;
		for(size_t left=0;left<6;++left)for(size_t a=0;a<old[left];++a)
			for(size_t g:cg.channels[b][left]){decode.push_back({g,a});group.push_back(oldbase[left]+a);}
		if(K[b].ncol!=decode.size()||K[b].nrow!=previous[b])throw std::runtime_error("Copy/channel layout mismatch");
		irr_parallel::grouped_entries index(K[b],group,groups);
		irr_parallel::rows(groups,&opt->pool,[&](size_t owner){
		for(size_t p=index.offsets[owner];p<index.offsets[owner+1];++p) {
			auto [a,position]=index.entries[p];auto col=K[b][a](position);const auto& value=K[b][a][position];
			auto [g,c]=decode.at(col);
			for(const auto& route:t.routes[b][g]) {
				const auto& row=t.rows[route.row];
				matrix[rowbase[route.row]+c*condition[row.condition].multiplicity].push_back(
					int32_t(colbase[b]+a*cg.channels[target][b].size()+route.channel),
					value*route.coefficient*(primitive?t.row_scale[route.row]:rat_t(1)));
			}
		}
		});
	}
	opt->pool.detach_loop(0,matrix.nrow,[&](size_t r){fec_irrep_merge_entries(matrix[r]);});opt->pool.wait();
	return matrix;
}

inline sparse_mat<rat_t,int32_t> solve(sparse_mat<rat_t,int32_t>& matrix,rref_option_t opt,const std::filesystem::path& cache={}) {
	const field_t F(FIELD_QQ);matrix.sort_rows_by_nnz();sparse_mat<rat_t,int32_t> kernel;
	if(matrix.ncol==0)kernel=sparse_mat<rat_t,int32_t>(0,0);
	else if(matrix.nrow==0) {
		kernel=sparse_mat<rat_t,int32_t>(matrix.ncol,matrix.ncol);
		for(size_t i=0;i<matrix.ncol;++i)kernel[i].push_back(int32_t(i),rat_t(1));
	}else{
		auto pivots=cache.empty()?sparse_mat_rref_reconstruct(matrix,opt):k_crt::rref(matrix,opt,cache);size_t rank=0;for(const auto& g:pivots)rank+=g.size();
		if(rank==matrix.ncol)kernel=sparse_mat<rat_t,int32_t>(0,matrix.ncol);
		else {kernel=sparse_mat_rref_kernel(matrix,pivots,F,opt);kernel.transpose_replace();}
	}
	matrix.clear();opt->pool.detach_loop(0,kernel.nrow,[&](size_t r){kernel[r].canonicalize();fec_irrep_cancel_divisor(kernel[r]);});
	opt->pool.wait();kernel.sort_rows_by_nnz();return kernel;
}

// Only seed conversion, explicit export and regression tests use this adapter.
inline sparse_mat<rat_t,int32_t> expand(const channels_t& cg,size_t target,const copies& old,
	const sparse_mat<rat_t,int32_t>& K) {
	size_t dim=0;copies offsets{};for(size_t b=0;b<6;++b){offsets[b]=dim;dim+=old[b]*specs()[b].irrep_dimension;}
	std::vector<std::tuple<size_t,size_t,size_t>> decode;
	for(size_t b=0;b<6;++b)for(size_t a=0;a<old[b];++a)for(size_t g:cg.channels[target][b])decode.push_back({b,a,g});
	if(K.ncol!=decode.size())throw std::runtime_error("Expansion column mismatch");
	sparse_mat<rat_t,int32_t> out(K.nrow*specs()[target].irrep_dimension,dim*9);
	for(size_t s=0;s<specs()[target].irrep_dimension;++s)for(size_t a=0;a<K.nrow;++a) {
		for(auto [col,value]:K[a]) {
			auto [b,copy,g]=decode.at(col);
			for(const auto& term:cg.cg[target].vectors[g].component[s]) {
				size_t component=term.coordinate/9-cg.axis[b].offset;
				out[s*K.nrow+a].push_back(int32_t((offsets[b]+component*old[b]+copy)*9+term.coordinate%9),value*term.coefficient);
			}
		}
		fec_irrep_merge_entries(out[s*K.nrow+a]);
	}
	return out;
}

// Row-basis action: B_i -> sum_j D_ij B_j; D is d x d, not (d*m) x (d*m).
inline sparse_mat<rat_t,int32_t> action(size_t b,const std::string& name) {
	const auto& s=specs().at(b);sparse_mat<rat_t,int32_t> D(s.irrep_dimension,s.irrep_dimension);
	if(name=="C"||name=="Cyclic") {
		rat_t p=name=="Cyclic"?rat_t(s.parity):rat_t(1);
		if(s.irrep==fec_irrep_t::E) {
			D[0].push_back(0,p*rat_t(-1,2));D[0].push_back(1,p*rat_t(3,2));
			D[1].push_back(0,p*rat_t(-1,2));D[1].push_back(1,p*rat_t(-1,2));
		}else D[0].push_back(0,p);
	}else if(name=="Flip"||name=="Parity") {
		for(size_t i=0;i<s.irrep_dimension;++i)D[i].push_back(int32_t(i),
			rat_t(name=="Parity"?s.parity:(s.irrep==fec_irrep_t::A2||i==1?-1:1)));
	}else throw std::runtime_error("Unknown compact action");
	return D;
}

inline sparse_mat<rat_t,int32_t> apply_action(size_t b,const std::string& name,
	const sparse_mat<rat_t,int32_t>& rows) {
	auto D=action(b,name);if(rows.nrow%D.nrow)throw std::runtime_error("Action component dimensions");
	size_t m=rows.nrow/D.nrow;sparse_mat<rat_t,int32_t> out(rows.nrow,rows.ncol);
	for(size_t s=0;s<D.nrow;++s)for(size_t a=0;a<m;++a) {
		for(auto [t,x]:D[s])for(auto [c,y]:rows[t*m+a])out[s*m+a].push_back(c,x*y);
		fec_irrep_merge_entries(out[s*m+a]);
	}
	return out;
}
} // namespace irr
