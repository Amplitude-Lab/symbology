#pragma once
#include "irr_recurrence.hpp"

namespace irr {
struct seam_row_t {
	size_t left, right, condition, copy, pivot;
	rat_t pivot_value;
	std::vector<std::pair<size_t,rat_t>> embedding;
};
struct seam_route_t { size_t row; rat_t value; };
struct seam_template_t {
	size_t condition_dimension;
	std::vector<seam_row_t> rows;
	std::array<std::vector<std::vector<std::vector<seam_route_t>>>,6> routes;
};

// Couple (old first x old last) to the condition irrep. The condition axis
// consists of equation covectors, so its pairing is evaluation, not diag(1,3).
// The two SYMBOL factors below do use the E singlet metric diag(1,3).
inline seam_template_t make_seam_template(const channels_t& cg,
	const sparse_tensor<rat_t,int32_t,SPARSE_CSR>& W,
	const std::vector<fec_irrep_axis_block_t>& condition,rref_option_t opt) {
	if(W.rank()!=3||W.dim(0)!=9||W.dim(1)!=9||W.dim(2)!=fec_irrep_axis_dimension(condition))
		throw std::runtime_error("Sewing condition shape/metadata mismatch");
	const field_t F(FIELD_QQ);seam_template_t out;out.condition_dimension=W.dim(2);
	size_t dim=fec_irrep_axis_dimension(cg.axis);
	for(size_t q=0;q<condition.size();++q) {
		auto product=make_hhex_codomain_basis(cg.axis,cg.axis,condition[q].spec,F,opt);
		for(const auto& v:product.vectors)for(size_t a=0;a<condition[q].multiplicity;++a) {
			size_t coord=v.component[0][0].coordinate;
			seam_row_t row{axis_block(cg.axis,coord/dim),axis_block(cg.axis,coord%dim),q,a,0,rat_t(0),{}};
			for(size_t s=0;s<condition[q].spec.irrep_dimension;++s)
				for(const auto& term:v.component[s]) {
					size_t raw=((term.coordinate/dim)*W.dim(2)+condition[q].offset+s*condition[q].multiplicity+a)*dim+term.coordinate%dim;
					row.embedding.push_back({raw,term.coefficient});
				}
			std::sort(row.embedding.begin(),row.embedding.end(),[](const auto& a,const auto& b){return a.first<b.first;});
			row.pivot=row.embedding.front().first;row.pivot_value=row.embedding.front().second;
			out.rows.push_back(std::move(row));
		}
	}
	std::array<std::array<std::vector<std::pair<size_t,rat_t>>,9>,9> lookup;
	for(size_t i=0;i<9;++i)for(size_t p=W.rowptr()[i];p<W.rowptr()[i+1];++p)
		lookup[i][W.index(p)[0]].push_back({size_t(W.index(p)[1]),W.val(p)});
	for(size_t b=0;b<6;++b) {
		out.routes[b].resize(cg.cg[b].size());
		for(size_t l=0;l<cg.cg[b].size();++l) {
			out.routes[b][l].resize(cg.cg[b].size());
			for(size_t r=0;r<cg.cg[b].size();++r) {
				std::map<size_t,rat_t> raw,restored;
				for(size_t s=0;s<specs()[b].irrep_dimension;++s)
					for(const auto& first:cg.cg[b].vectors[l].component[s])
						for(const auto& last:cg.cg[b].vectors[r].component[s])
							for(const auto& [q,value]:lookup[first.coordinate%9][last.coordinate%9])
								raw[((first.coordinate/9)*W.dim(2)+q)*dim+last.coordinate/9]
									+=first.coefficient*last.coefficient*value*(s==1?rat_t(3):rat_t(1));
				for(size_t i=0;i<out.rows.size();++i) {
					const auto& row=out.rows[i];auto pos=raw.find(row.pivot);
					if(pos==raw.end()||pos->second==0)continue;
					rat_t value=pos->second/row.pivot_value;
					out.routes[b][l][r].push_back({i,value});
					for(const auto& [coord,c]:row.embedding)restored[coord]+=value*c;
				}
				// A small exact check on EVERY channel, not on sampled large matrices.
				for(const auto& [coord,c]:raw)if(restored[coord]!=c) {
					std::ostringstream message;message<<"Sewing condition chart mismatch b="<<b<<" l="<<l<<" r="<<r<<" coordinate="<<coord<<" actual="<<c<<" restored="<<restored[coord];
					throw std::runtime_error(message.str());
				}
				for(const auto& [coord,c]:restored)if(raw[coord]!=c)
					throw std::runtime_error("Sewing row chart loses a representation component");
			}
		}
	}
	return out;
}

struct seam_layout_t {
	std::vector<size_t> rowbase;
	copies colbase{};
	size_t rows=0,columns=0;
};
inline seam_layout_t seam_layout(const seam_template_t& t,const copies& oldfirst,
	const copies& oldlast,const copies& first,const copies& last) {
	seam_layout_t layout;
	for(const auto& r:t.rows){layout.rowbase.push_back(layout.rows);layout.rows+=oldfirst[r.left]*oldlast[r.right];}
	for(size_t b=0;b<6;++b){layout.colbase[b]=layout.columns;layout.columns+=first[b]*last[b];}
	if(layout.rows>size_t(INT32_MAX)||layout.columns>size_t(INT32_MAX))throw std::runtime_error("Sewing matrix exceeds int32 index range");
	return layout;
}
struct decoded_channel_t { size_t block,copy,channel; };
inline std::vector<decoded_channel_t> decode_channels(const channels_t& cg,size_t target,const copies& old) {
	std::vector<decoded_channel_t> result;
	for(size_t b=0;b<6;++b)for(size_t a=0;a<old[b];++a)for(size_t g:cg.channels[target][b])result.push_back({b,a,g});
	return result;
}

// Only multiplicity coefficients enter this assembly. Each worker owns one
// candidate row, avoiding concurrent sparse-vector insertions; transpose once.
inline sparse_mat<rat_t,int32_t> assemble_seam(const channels_t& cg,const seam_template_t& t,
	const seam_layout_t& layout,const copies& oldfirst,const copies& oldlast,
	const matrices& first,const matrices& last,rref_option_t opt) {
	irr_parallel::phase timing("sewing_assembly");
	sparse_mat<rat_t,int32_t> image(layout.columns,layout.rows);
	for(size_t b=0;b<6;++b) {
		auto left=decode_channels(cg,b,oldfirst),right=decode_channels(cg,b,oldlast);
		if(first[b].ncol!=left.size()||last[b].ncol!=right.size())throw std::runtime_error("Sewing input copy/channel layout mismatch");
		opt->pool.detach_loop(0,first[b].nrow*last[b].nrow,[&](size_t i) {
			auto& row=image[layout.colbase[b]+i];
			for(auto [j,x]:first[b][i/last[b].nrow])for(auto [k,y]:last[b][i%last[b].nrow]) {
				const auto& l=left[j];const auto& r=right[k];
				for(const auto& route:t.routes[b][l.channel][r.channel])
					row.push_back(int32_t(layout.rowbase[route.row]+l.copy*oldlast[r.block]+r.copy),x*y*route.value);
			}
			fec_irrep_merge_entries(row);
		});
		opt->pool.wait();
	}
	image.transpose_replace();
	return image;
}
}
