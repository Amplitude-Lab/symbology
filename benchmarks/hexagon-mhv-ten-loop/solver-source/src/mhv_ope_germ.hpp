#pragma once
#include "mhv_sewing.hpp"

// Ordered T->0 then z=1/S->0 expansion. Logarithms remain formal variables.
// A last-letter integral is a linear operator, not multiplication by a tail function.
struct mhv_ope_germ_t {
	size_t weight, logs, stride;
	static constexpr int max_z = 4, z_size = 9;
	field_t F;
	struct term_t { int t, z; ulong value; };
	std::array<std::vector<term_t>, 9> kernels;
	std::array<ulong, 9> at{}, az{};
	std::vector<ulong> inverse;

	mhv_ope_germ_t(size_t w, ulong prime, const std::filesystem::path& file)
		: weight(w), logs(w + 1), stride(3 * z_size * logs * logs), F(FIELD_Fp, prime), inverse(w + 8) {
		for (size_t i = 1; i < inverse.size(); ++i) inverse[i] = n_invmod(i, prime);
		std::ifstream stream(file); int letter, t, z; std::string coefficient;
		while (stream >> letter >> t >> z >> coefficient) {
			if (coefficient.front() == '"') coefficient = coefficient.substr(1, coefficient.size()-2);
			if (letter < 0 || letter >= 9 || t < -1 || t > 2) throw std::runtime_error("Invalid OPE kernel index");
			const ulong value = rat_t(coefficient) % F.mod;
			kernels[letter].push_back({t,z,value});
			if (t == 0 && z == 0) at[letter] = value;
			if (t == -1 && z == 0) az[letter] = value;
		}
		if (!stream.eof() || kernels[0].empty()) throw std::runtime_error("Missing or invalid OPE kernel: " + file.string());
	}
	ulong add(ulong a, ulong b) const { return nmod_add(a,b,F.mod); }
	ulong mul(ulong a, ulong b) const { return nmod_mul(a,b,F.mod); }
	size_t index(int t, int z, size_t a, size_t b) const { return ((t*z_size+z+2)*logs+a)*logs+b; }
	bool valid(int t, int z) const { return t>=0 && t<=2 && z>=-t && z<=6-t; }

	void integrate(const ulong* previous, ulong* next, size_t letter, size_t w) const {
		std::fill(next,next+stride,0);
		for (auto term : kernels[letter])
			for (int t=0;t<=2;++t) for (int z=-t;z<=6-t;++z) {
				if ((t==0 && z==0) || (t==0)!=(term.t==-1)) continue;
				int old_t=t-(term.t==-1?0:term.t), old_z=z-term.z;
				if (!valid(old_t,old_z)) continue;
				for (size_t a=0;a<w;++a) for (size_t b=0;b<w-a;++b) {
					auto value=previous[index(old_t,old_z,a,b)];
					if (value) { auto& v=next[index(t,z,a,b)]; v=add(v,mul(value,term.value)); }
				}
			}
		for (int t=0;t<=2;++t) for (int z=-t;z<=6-t;++z) {
			if (t==0 && z==0) continue;
			for (int a=static_cast<int>(w);a>=0;--a) for (int b=static_cast<int>(w)-a;b>=0;--b) {
				auto& value=next[index(t,z,a,b)];
				if (t && a+1<logs) value=nmod_sub(value,mul(a+1,next[index(t,z,a+1,b)]),F.mod);
				if (!t && b+1<logs) value=nmod_sub(value,mul(b+1,next[index(t,z,a,b+1)]),F.mod);
				value=mul(value,inverse[t?t:z]);
			}
		}
		for (size_t a=0;a<w;++a) for (size_t b=0;b<w-a;++b) {
			const auto value=mul(previous[index(0,0,a,b)],inverse[w]);
			next[index(0,0,a+1,b)]=add(next[index(0,0,a+1,b)],mul(value,at[letter]));
			next[index(0,0,a,b+1)]=add(next[index(0,0,a,b+1)],mul(value,az[letter]));
		}
	}

	std::vector<ulong> head(const std::filesystem::path& root, size_t maximum, thread_pool* pool) const {
		std::vector<ulong> previous(stride); previous[index(0,0,0,0)]=1;
		for (size_t w=1;w<=maximum;++w) {
			auto dir=nmhv_limit_weight_directory(root,w);
			auto blocks=nmhv_limit_read_irrep_axis(dir/("FEC_"+std::to_string(w)+"_blocks.tsv"));
			std::vector<ulong> current(fec_irrep_axis_dimension(blocks)*stride);
			// Only one letter image is resident at a time; no ninefold permanent buffer.
			std::vector<ulong> image(previous.size());
			for (size_t letter=0;letter<9;++letter) {
				pool->detach_loop(0,previous.size()/stride,[&](size_t r){integrate(previous.data()+r*stride,image.data()+r*stride,letter,w);});
				pool->wait();
				for (auto& block:blocks) if (block.dimension) {
					auto tensor=sparse_tensor_read_wxf<ulong,int32_t>(dir/("FEC_"+std::to_string(w)+"_"+block.spec.name+".wxf"),F,pool);
					if (tensor.rank()!=3 || tensor.dim(1)!=previous.size()/stride || tensor.dim(2)!=9) throw std::runtime_error("OPE FEC shape mismatch");
					pool->detach_loop(0,block.dimension,[&](size_t r){
						auto target=current.data()+(block.offset+r)*stride;
						for(size_t p=tensor.rowptr()[r];p<tensor.rowptr()[r+1];++p) {
							auto i=tensor.index(p); if(i[1]!=letter) continue;
							auto source=image.data()+i[0]*stride;
							for(size_t j=0;j<stride;++j) if(source[j]) target[j]=add(target[j],mul(source[j],tensor.val(p)));
						}
					}); pool->wait();
				}
			}
			previous=std::move(current);
		}
		return previous;
	}

	template<typename index_t>
	std::vector<ulong> project(const std::vector<ulong>& head_values,
		const sparse_mat<ulong,index_t>& sewing, const std::vector<joint_eet_lec_t<index_t>>& chain,
		size_t h, thread_pool* pool) const {
		const size_t tail_dimension=chain.empty()?1:chain.back().dimension();
		if(sewing.ncol!=(head_values.size()/stride)*tail_dimension) throw std::runtime_error("OPE sewing shape mismatch");
		std::vector<ulong> result(sewing.nrow*stride);
		// Process one solution at a time. Do not materialize head x tail x all solutions.
		for(size_t solution=0;solution<sewing.nrow;++solution) {
			std::vector<ulong> values(tail_dimension*stride);
			for(auto [column,coefficient]:sewing[solution]) {
				auto source=head_values.data()+(column/tail_dimension)*stride;
				auto target=values.data()+(column%tail_dimension)*stride;
				for(size_t j=0;j<stride;++j) if(source[j]) target[j]=add(target[j],mul(coefficient,source[j]));
			}
			for(size_t k=chain.size();k>0;--k) {
				const auto& lec=chain[k-1].LEC;
				std::vector<ulong> next(lec.dim(2)*stride), image(stride);
				for(size_t r=0;r<lec.dim(0);++r) for(size_t letter=0;letter<9;++letter) {
					bool used=false;
					for(size_t p=lec.rowptr()[r];p<lec.rowptr()[r+1];++p) if(lec.index(p)[0]==letter) {used=true;break;}
					if(!used) continue;
					integrate(values.data()+r*stride,image.data(),letter,h+chain.size()-k+1);
					for(size_t p=lec.rowptr()[r];p<lec.rowptr()[r+1];++p) {
						auto i=lec.index(p);if(i[0]!=letter) continue;
						auto target=next.data()+i[1]*stride;const ulong coefficient=lec.val(p)%F.mod;
						for(size_t j=0;j<stride;++j) if(image[j]) target[j]=add(target[j],mul(coefficient,image[j]));
					}
				}
				values=std::move(next);
			}
			if(values.size()!=stride) throw std::runtime_error("MHV tail must end in one scalar component");
			std::copy(values.begin(),values.end(),result.begin()+solution*stride);
		}
		return result;
	}
};
