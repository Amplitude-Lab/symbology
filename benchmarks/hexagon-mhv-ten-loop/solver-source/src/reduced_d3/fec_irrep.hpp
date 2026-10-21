#pragma once

#ifndef USE_MIMALLOC
#define USE_MIMALLOC 1
#endif

#include "sparse_mat.h"
#include "sparse_tensor.h"
#include "wxf_support.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace SparseRREF;

enum class fec_symmetry_t {
	parity_flip,
	parity_dihedral
};

enum class fec_irrep_t {
	reflection_even,
	reflection_odd,
	A1,
	A2,
	E
};

struct fec_block_spec_t {
	std::string name;
	int parity;
	fec_irrep_t irrep;
	size_t irrep_dimension;
};

struct fec_block_t {
	fec_block_spec_t spec;
	size_t multiplicity = 0;
	sparse_tensor<rat_t, int32_t, SPARSE_CSR> tensor;
};

struct fec_raw_term_t {
	int32_t coordinate;
	rat_t coefficient;
};

struct fec_candidate_t {
	std::array<std::vector<fec_raw_term_t>, 2> component;
};

struct fec_candidate_lookup_entry_t {
	int32_t column;
	rat_t coefficient;
};

struct fec_candidate_basis_t {
	std::string block;
	int parity = 0;
	fec_irrep_t irrep = fec_irrep_t::A1;
	size_t irrep_dimension = 1;
	size_t previous_dimension = 0;
	std::vector<fec_candidate_t> vectors;
	std::vector<std::vector<fec_candidate_lookup_entry_t>> first_component_lookup;
	std::vector<bool> active_previous_blocks;

	size_t size() const { return vectors.size(); }
};

struct fec_extension_metrics_t {
	size_t candidate_columns = 0;
	size_t matrix_rows = 0;
	size_t matrix_nnz = 0;
	size_t rank = 0;
	size_t output_nnz = 0;
};

inline bool validate_hhex_candidate_symmetry(
	const fec_candidate_basis_t& candidates,
	const std::vector<fec_block_t>& blocks);

inline const std::vector<fec_block_spec_t>& hhex_block_specs(fec_symmetry_t symmetry) {
	static const std::vector<fec_block_spec_t> parity_flip{
		{"pe_re", 1, fec_irrep_t::reflection_even, 1},
		{"pe_ro", 1, fec_irrep_t::reflection_odd, 1},
		{"po_re", -1, fec_irrep_t::reflection_even, 1},
		{"po_ro", -1, fec_irrep_t::reflection_odd, 1}
	};
	static const std::vector<fec_block_spec_t> parity_dihedral{
		{"pe_A1", 1, fec_irrep_t::A1, 1},
		{"pe_A2", 1, fec_irrep_t::A2, 1},
		{"pe_E", 1, fec_irrep_t::E, 2},
		{"po_A1", -1, fec_irrep_t::A1, 1},
		{"po_A2", -1, fec_irrep_t::A2, 1},
		{"po_E", -1, fec_irrep_t::E, 2}
	};
	return symmetry == fec_symmetry_t::parity_flip
		? parity_flip : parity_dihedral;
}

inline size_t hhex_global_dimension(const std::vector<fec_block_t>& blocks) {
	size_t result = 0;
	for (const auto& block : blocks)
		result += block.tensor.dim(0);
	return result;
}

inline void hhex_validate_blocks(
	fec_symmetry_t symmetry, const std::vector<fec_block_t>& blocks) {
	const auto& specs = hhex_block_specs(symmetry);
	if (blocks.size() != specs.size())
		throw std::runtime_error("Wrong number of FEC symmetry blocks.");
	size_t predecessor_dimension = blocks.empty() ? 0 : blocks[0].tensor.dim(1);
	for (size_t i = 0; i < blocks.size(); ++i) {
		const auto& block = blocks[i];
		if (block.spec.name != specs[i].name
			|| block.spec.parity != specs[i].parity
			|| block.spec.irrep != specs[i].irrep
			|| block.spec.irrep_dimension != specs[i].irrep_dimension)
			throw std::runtime_error("FEC symmetry blocks are not in canonical order.");
		if (block.tensor.rank() != 3 || block.tensor.dim(2) != 9
			|| block.tensor.dim(1) != predecessor_dimension)
			throw std::runtime_error("FEC block must have layout {basis,previous_basis,9}.");
		if (block.tensor.dim(0) != block.multiplicity * block.spec.irrep_dimension)
			throw std::runtime_error("FEC block multiplicity does not match its first dimension.");
	}
}

inline fec_block_t make_hhex_block(
	fec_symmetry_t symmetry, const std::string& name,
	sparse_tensor<rat_t, int32_t, SPARSE_CSR>&& tensor) {
	for (const auto& spec : hhex_block_specs(symmetry)) {
		if (spec.name != name)
			continue;
		if (tensor.rank() != 3 || tensor.dim(2) != 9
			|| tensor.dim(0) % spec.irrep_dimension != 0)
			throw std::runtime_error("Invalid FEC block dimensions for " + name + ".");
		return {spec, tensor.dim(0) / spec.irrep_dimension, std::move(tensor)};
	}
	throw std::runtime_error("Unknown FEC symmetry block: " + name);
}

inline void fec_irrep_merge_entries(sparse_vec<rat_t, int32_t>& row) {
	row.sort_indices();
	size_t result = 0;
	for (size_t i = 0; i < row.nnz(); ++i) {
		if (result != 0 && row(result - 1) == row(i))
			row[result - 1] += row[i];
		else {
			row(result) = row(i);
			row[result] = row[i];
			++result;
		}
	}
	row.resize(result);
	row.compress();
}

inline void fec_irrep_cancel_divisor(sparse_vec<rat_t, int32_t>& row) {
	if (row.nnz() == 0)
		return;
	int_t denominator = 1;
	for (size_t i = 0; i < row.nnz(); ++i)
		denominator = LCM(denominator, row[i].den());
	sparse_vec_rescale(row, rat_t(denominator), field_t(FIELD_QQ));
	int_t divisor = 0;
	for (size_t i = 0; i < row.nnz(); ++i)
		divisor = GCD(divisor, row[i].num());
	if (divisor != 0)
		sparse_vec_rescale(row, rat_t(1, divisor), field_t(FIELD_QQ));
}

inline void fec_candidate_add_term(
	fec_candidate_t& candidate, size_t component,
	size_t previous, size_t letter, const rat_t& coefficient) {
	if (coefficient == 0)
		return;
	const int32_t coordinate = static_cast<int32_t>(previous * 9 + letter);
	for (auto& term : candidate.component[component]) {
		if (term.coordinate == coordinate) {
			term.coefficient += coefficient;
			return;
		}
	}
	candidate.component[component].push_back({coordinate, coefficient});
}

inline void fec_candidate_add_letter_vector(
	fec_candidate_t& candidate, size_t component, size_t previous,
	size_t triplet, const std::array<int, 3>& coefficients,
	const rat_t& scale = rat_t(1)) {
	for (size_t i = 0; i < 3; ++i)
		fec_candidate_add_term(
			candidate, component, previous, 3 * triplet + i,
			scale * rat_t(coefficients[i]));
}

inline fec_candidate_basis_t make_hhex_flip_candidates(
	const std::vector<fec_block_t>& blocks, const fec_block_spec_t& target) {
	fec_candidate_basis_t result;
	result.block = target.name;
	result.parity = target.parity;
	result.irrep = target.irrep;
	result.irrep_dimension = 1;
	result.previous_dimension = hhex_global_dimension(blocks);
	result.active_previous_blocks.assign(blocks.size(), false);
	const int target_reflection = target.irrep == fec_irrep_t::reflection_even ? 1 : -1;
	const std::array<int32_t, 9> reflection{2, 1, 0, 5, 4, 3, 8, 7, 6};

	size_t offset = 0;
	for (size_t block_index = 0; block_index < blocks.size(); ++block_index) {
		const auto& block = blocks[block_index];
		const int old_reflection = block.spec.irrep == fec_irrep_t::reflection_even ? 1 : -1;
		for (size_t previous = 0; previous < block.tensor.dim(0); ++previous) {
			for (int32_t letter = 0; letter < 9; ++letter) {
				const int32_t image = reflection[letter];
				if (letter > image)
					continue;
				const int letter_parity = letter < 6 ? 1 : -1;
				if (block.spec.parity * letter_parity != target.parity)
					continue;
				if (letter == image && old_reflection != target_reflection)
					continue;
				fec_candidate_t candidate;
				fec_candidate_add_term(candidate, 0, offset + previous, letter, rat_t(1));
				if (letter != image)
					fec_candidate_add_term(
						candidate, 0, offset + previous, image,
						rat_t(target_reflection * old_reflection));
				result.vectors.push_back(std::move(candidate));
				result.active_previous_blocks[block_index] = true;
			}
		}
		offset += block.tensor.dim(0);
	}
	return result;
}

inline fec_candidate_basis_t make_hhex_dihedral_candidates(
	const std::vector<fec_block_t>& blocks, const fec_block_spec_t& target) {
	fec_candidate_basis_t result;
	result.block = target.name;
	result.parity = target.parity;
	result.irrep = target.irrep;
	result.irrep_dimension = target.irrep_dimension;
	result.previous_dimension = hhex_global_dimension(blocks);
	result.active_previous_blocks.assign(blocks.size(), false);
	const std::array<int, 3> letter_A1{1, 1, 1};
	const std::array<int, 3> letter_Ep{1, -2, 1};
	const std::array<int, 3> letter_Em{1, 0, -1};

	size_t offset = 0;
	for (size_t block_index = 0; block_index < blocks.size(); ++block_index) {
		const auto& block = blocks[block_index];
		const size_t candidates_before_block = result.vectors.size();
		const int letter_parity = target.parity * block.spec.parity;
		const size_t triplet_begin = letter_parity == 1 ? 0 : 2;
		const size_t triplet_end = letter_parity == 1 ? 2 : 3;
		for (size_t multiplicity = 0; multiplicity < block.multiplicity; ++multiplicity) {
			for (size_t triplet = triplet_begin; triplet < triplet_end; ++triplet) {
				const size_t old_single = offset + multiplicity;
				const size_t old_Ep = offset + multiplicity;
				const size_t old_Em = offset + block.multiplicity + multiplicity;

				if (target.irrep == fec_irrep_t::A1) {
					if (block.spec.irrep == fec_irrep_t::A1) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_single, triplet, letter_A1);
						result.vectors.push_back(std::move(candidate));
					} else if (block.spec.irrep == fec_irrep_t::E) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_Ep, triplet, letter_Ep);
						fec_candidate_add_letter_vector(candidate, 0, old_Em, triplet, letter_Em, rat_t(3));
						result.vectors.push_back(std::move(candidate));
					}
				} else if (target.irrep == fec_irrep_t::A2) {
					if (block.spec.irrep == fec_irrep_t::A2) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_single, triplet, letter_A1);
						result.vectors.push_back(std::move(candidate));
					} else if (block.spec.irrep == fec_irrep_t::E) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_Ep, triplet, letter_Em, rat_t(-1));
						fec_candidate_add_letter_vector(candidate, 0, old_Em, triplet, letter_Ep);
						result.vectors.push_back(std::move(candidate));
					}
				} else {
					if (block.spec.irrep == fec_irrep_t::A1) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_single, triplet, letter_Ep);
						fec_candidate_add_letter_vector(candidate, 1, old_single, triplet, letter_Em);
						result.vectors.push_back(std::move(candidate));
					} else if (block.spec.irrep == fec_irrep_t::A2) {
						fec_candidate_t candidate;
						fec_candidate_add_letter_vector(candidate, 0, old_single, triplet, letter_Em, rat_t(3));
						fec_candidate_add_letter_vector(candidate, 1, old_single, triplet, letter_Ep, rat_t(-1));
						result.vectors.push_back(std::move(candidate));
					} else if (block.spec.irrep == fec_irrep_t::E) {
						fec_candidate_t candidate_A1;
						fec_candidate_add_letter_vector(candidate_A1, 0, old_Ep, triplet, letter_A1);
						fec_candidate_add_letter_vector(candidate_A1, 1, old_Em, triplet, letter_A1);
						result.vectors.push_back(std::move(candidate_A1));

						fec_candidate_t candidate_E;
						fec_candidate_add_letter_vector(candidate_E, 0, old_Ep, triplet, letter_Ep, rat_t(-1));
						fec_candidate_add_letter_vector(candidate_E, 0, old_Em, triplet, letter_Em, rat_t(3));
						fec_candidate_add_letter_vector(candidate_E, 1, old_Ep, triplet, letter_Em);
						fec_candidate_add_letter_vector(candidate_E, 1, old_Em, triplet, letter_Ep);
						result.vectors.push_back(std::move(candidate_E));
					}
				}
			}
		}
		if (result.vectors.size() != candidates_before_block)
			result.active_previous_blocks[block_index] = true;
		offset += block.tensor.dim(0);
	}
	return result;
}

inline fec_candidate_basis_t make_hhex_candidates(
	fec_symmetry_t symmetry, const std::vector<fec_block_t>& blocks,
	const std::string& block_name) {
	hhex_validate_blocks(symmetry, blocks);
	for (const auto& target : hhex_block_specs(symmetry)) {
		if (target.name != block_name)
			continue;
		auto result = symmetry == fec_symmetry_t::parity_flip
			? make_hhex_flip_candidates(blocks, target)
			: make_hhex_dihedral_candidates(blocks, target);
		result.first_component_lookup.resize(result.previous_dimension * 9);
		for (size_t column = 0; column < result.vectors.size(); ++column) {
			for (const auto& term : result.vectors[column].component[0])
				result.first_component_lookup[term.coordinate].push_back(
					{static_cast<int32_t>(column), term.coefficient});
		}
		return result;
	}
	throw std::runtime_error("Unknown FEC symmetry block: " + block_name);
}

inline std::vector<fec_block_t> make_hhex_seed_blocks(
	fec_symmetry_t symmetry, thread_pool* pool) {
	std::vector<fec_block_t> result;
	if (symmetry == fec_symmetry_t::parity_flip) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> even({2, 1, 9});
		even.push_back({0, 0, 0}, rat_t(1));
		even.push_back({0, 0, 2}, rat_t(1));
		even.push_back({1, 0, 1}, rat_t(1));
		sparse_tensor<rat_t, int32_t, SPARSE_COO> odd({1, 1, 9});
		odd.push_back({0, 0, 0}, rat_t(1));
		odd.push_back({0, 0, 2}, rat_t(-1));
		result.push_back(make_hhex_block(symmetry, "pe_re",
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(even), pool)));
		result.push_back(make_hhex_block(symmetry, "pe_ro",
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(odd), pool)));
		for (const auto& name : {"po_re", "po_ro"}) {
			sparse_tensor<rat_t, int32_t, SPARSE_COO> empty({0, 1, 9});
			result.push_back(make_hhex_block(symmetry, name,
				sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(empty), pool)));
		}
	} else {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> A1({1, 1, 9});
		for (int32_t letter = 0; letter < 3; ++letter)
			A1.push_back({0, 0, letter}, rat_t(1));
		sparse_tensor<rat_t, int32_t, SPARSE_COO> E({2, 1, 9});
		for (int32_t letter = 0; letter < 3; ++letter)
			E.push_back({0, 0, letter}, rat_t(std::array<int, 3>{1, -2, 1}[letter]));
		E.push_back({1, 0, 0}, rat_t(1));
		E.push_back({1, 0, 2}, rat_t(-1));
		result.push_back(make_hhex_block(symmetry, "pe_A1",
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(A1), pool)));
		for (const auto& name : {"pe_A2"}) {
			sparse_tensor<rat_t, int32_t, SPARSE_COO> empty({0, 1, 9});
			result.push_back(make_hhex_block(symmetry, name,
				sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(empty), pool)));
		}
		result.push_back(make_hhex_block(symmetry, "pe_E",
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(E), pool)));
		for (const auto& name : {"po_A1", "po_A2", "po_E"}) {
			sparse_tensor<rat_t, int32_t, SPARSE_COO> empty({0, 1, 9});
			result.push_back(make_hhex_block(symmetry, name,
				sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(empty), pool)));
		}
	}
	hhex_validate_blocks(symmetry, result);
	return result;
}

inline fec_block_t extend_hhex_block(
	fec_symmetry_t symmetry,
	const sparse_tensor<rat_t, int32_t, SPARSE_CSR>& dlogmat,
	const std::vector<fec_block_t>& blocks, const std::string& block_name,
	const field_t& F, rref_option_t& opt,
	fec_extension_metrics_t* metrics = nullptr) {
	hhex_validate_blocks(symmetry, blocks);
	if (dlogmat.rank() != 3 || dlogmat.dim(0) != 9 || dlogmat.dim(1) != 9)
		throw std::runtime_error("Condition tensor must have layout {9,9,condition}.");
	auto candidates = make_hhex_candidates(symmetry, blocks, block_name);
	if (!validate_hhex_candidate_symmetry(candidates, blocks))
		throw std::runtime_error("FEC candidate basis violates its symmetry label.");
	fec_extension_metrics_t local_metrics;
	local_metrics.candidate_columns = candidates.size();
	local_metrics.matrix_rows = blocks[0].tensor.dim(1) * dlogmat.dim(2);
	if (candidates.size() == 0) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> empty(
			{0, candidates.previous_dimension, 9});
		if (metrics != nullptr)
			*metrics = local_metrics;
		return make_hhex_block(symmetry, block_name,
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(empty), &(opt->pool)));
	}

	struct dlog_entry_t {
		int32_t letter;
		int32_t relation;
		rat_t coefficient;
	};
	std::array<std::vector<dlog_entry_t>, 9> dlog_lookup;
	for (size_t letter = 0; letter < 9; ++letter) {
		for (size_t pos = dlogmat.rowptr()[letter]; pos < dlogmat.rowptr()[letter + 1]; ++pos) {
			const auto index = dlogmat.index(pos);
			dlog_lookup[letter].push_back({index[0], index[1], dlogmat.val(pos)});
		}
	}

	sparse_mat<rat_t, int32_t> matrix(local_metrics.matrix_rows, candidates.size());
	size_t offset = 0;
	for (size_t block_index = 0; block_index < blocks.size(); ++block_index) {
		const auto& block = blocks[block_index];
		if (!candidates.active_previous_blocks[block_index]) {
			offset += block.tensor.dim(0);
			continue;
		}
		for (size_t previous = 0; previous < block.tensor.dim(0); ++previous) {
			for (size_t pos = block.tensor.rowptr()[previous];
				pos < block.tensor.rowptr()[previous + 1]; ++pos) {
				const auto index = block.tensor.index(pos);
				for (const auto& dlog_entry : dlog_lookup[index[1]]) {
					const auto& lookup = candidates.first_component_lookup[
						(offset + previous) * 9 + dlog_entry.letter];
					for (const auto& candidate_entry : lookup) {
						matrix[index[0] * dlogmat.dim(2) + dlog_entry.relation].push_back(
							candidate_entry.column,
							block.tensor.val(pos) * dlog_entry.coefficient
								* candidate_entry.coefficient);
					}
				}
			}
		}
		offset += block.tensor.dim(0);
	}
	opt->pool.detach_loop(0, matrix.nrow, [&](size_t row) {
		fec_irrep_merge_entries(matrix[row]);
	});
	opt->pool.wait();
	matrix.sort_rows_by_nnz();
	local_metrics.matrix_nnz = matrix.nnz();
	auto pivots = sparse_mat_rref_reconstruct(matrix, opt);
	for (const auto& group : pivots)
		local_metrics.rank += group.size();
	sparse_mat<rat_t, int32_t> kernel;
	if (local_metrics.rank == candidates.size())
		kernel = sparse_mat<rat_t, int32_t>(0, candidates.size());
	else
		kernel = sparse_mat_rref_kernel(matrix, pivots, F, opt).transpose();
	matrix.clear();
	opt->pool.detach_loop(0, kernel.nrow, [&](size_t row) {
		kernel[row].canonicalize();
		fec_irrep_cancel_divisor(kernel[row]);
	});
	opt->pool.wait();
	kernel.sort_rows_by_nnz();
	if (kernel.nrow == 0) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> empty(
			{0, candidates.previous_dimension, 9});
		if (metrics != nullptr)
			*metrics = local_metrics;
		return make_hhex_block(symmetry, block_name,
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(empty), &(opt->pool)));
	}

	// The two E components are expanded with the same multiplicity kernel.
	// Their component-major row order is fixed here and is never RREF'ed separately.
	sparse_mat<rat_t, int32_t> expanded(
		kernel.nrow * candidates.irrep_dimension, candidates.previous_dimension * 9);
	for (size_t component = 0; component < candidates.irrep_dimension; ++component) {
		for (size_t row = 0; row < kernel.nrow; ++row) {
			auto& output_row = expanded[component * kernel.nrow + row];
			for (size_t pos = 0; pos < kernel[row].nnz(); ++pos) {
				const auto& candidate = candidates.vectors[kernel[row](pos)];
				for (const auto& term : candidate.component[component])
					output_row.push_back(term.coordinate, kernel[row][pos] * term.coefficient);
			}
			fec_irrep_merge_entries(output_row);
		}
	}
	local_metrics.output_nnz = expanded.nnz();
	sparse_tensor<rat_t, int32_t, SPARSE_COO> output(expanded, &(opt->pool));
	output.reshape({expanded.nrow, candidates.previous_dimension, 9});
	if (metrics != nullptr)
		*metrics = local_metrics;
	return make_hhex_block(symmetry, block_name,
		sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(output), &(opt->pool)));
}

inline std::vector<fec_block_t> extend_hhex_blocks(
	fec_symmetry_t symmetry,
	const sparse_tensor<rat_t, int32_t, SPARSE_CSR>& dlogmat,
	const std::vector<fec_block_t>& blocks,
	const field_t& F, rref_option_t& opt) {
	hhex_validate_blocks(symmetry, blocks);
	size_t candidate_dimension = 0;
	for (const auto& spec : hhex_block_specs(symmetry)) {
		auto candidates = make_hhex_candidates(symmetry, blocks, spec.name);
		candidate_dimension += candidates.size() * candidates.irrep_dimension;
	}
	if (candidate_dimension != 9 * hhex_global_dimension(blocks))
		throw std::runtime_error("Symmetry-adapted candidate decomposition is incomplete.");
	std::vector<fec_block_t> result;
	for (const auto& spec : hhex_block_specs(symmetry))
		result.push_back(extend_hhex_block(
			symmetry, dlogmat, blocks, spec.name, F, opt));
	return result;
}

inline std::map<size_t, rat_t> fec_irrep_component_map(
	const std::vector<fec_raw_term_t>& component) {
	std::map<size_t, rat_t> result;
	for (const auto& term : component)
		result[term.coordinate] += term.coefficient;
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline std::vector<std::pair<size_t, rat_t>> fec_irrep_old_cyclic_image(
	const std::vector<fec_block_t>& blocks, size_t previous) {
	size_t offset = 0;
	for (const auto& block : blocks) {
		if (previous >= offset + block.tensor.dim(0)) {
			offset += block.tensor.dim(0);
			continue;
		}
		const size_t local = previous - offset;
		if (block.spec.irrep != fec_irrep_t::E)
			return {{previous, rat_t(1)}};
		const size_t m = block.multiplicity;
		if (local < m)
			return {{offset + local, rat_t(-1, 2)}, {offset + m + local, rat_t(3, 2)}};
		return {{offset + local - m, rat_t(-1, 2)}, {offset + local, rat_t(-1, 2)}};
	}
	throw std::runtime_error("Raw candidate predecessor is out of range.");
}

inline int fec_irrep_old_reflection_sign(
	const std::vector<fec_block_t>& blocks, size_t previous) {
	size_t offset = 0;
	for (const auto& block : blocks) {
		if (previous >= offset + block.tensor.dim(0)) {
			offset += block.tensor.dim(0);
			continue;
		}
		if (block.spec.irrep == fec_irrep_t::reflection_odd
			|| block.spec.irrep == fec_irrep_t::A2)
			return -1;
		if (block.spec.irrep == fec_irrep_t::E
			&& previous - offset >= block.multiplicity)
			return -1;
		return 1;
	}
	throw std::runtime_error("Raw candidate predecessor is out of range.");
}

inline int fec_irrep_old_parity(
	const std::vector<fec_block_t>& blocks, size_t previous) {
	size_t offset = 0;
	for (const auto& block : blocks) {
		if (previous < offset + block.tensor.dim(0))
			return block.spec.parity;
		offset += block.tensor.dim(0);
	}
	throw std::runtime_error("Raw candidate predecessor is out of range.");
}

inline std::map<size_t, rat_t> fec_irrep_cyclic_image(
	const std::vector<fec_raw_term_t>& component,
	const std::vector<fec_block_t>& blocks) {
	std::map<size_t, rat_t> result;
	for (const auto& term : component) {
		const size_t previous = static_cast<size_t>(term.coordinate) / 9;
		const size_t letter = static_cast<size_t>(term.coordinate) % 9;
		const size_t image = 3 * (letter / 3) + (letter + 1) % 3;
		for (const auto& [old_image, coefficient] :
			fec_irrep_old_cyclic_image(blocks, previous))
			result[old_image * 9 + image] += term.coefficient * coefficient;
	}
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline std::map<size_t, rat_t> fec_irrep_reflection_image(
	const std::vector<fec_raw_term_t>& component,
	const std::vector<fec_block_t>& blocks) {
	const std::array<size_t, 9> reflection{2, 1, 0, 5, 4, 3, 8, 7, 6};
	std::map<size_t, rat_t> result;
	for (const auto& term : component) {
		const size_t previous = static_cast<size_t>(term.coordinate) / 9;
		const size_t letter = static_cast<size_t>(term.coordinate) % 9;
		result[previous * 9 + reflection[letter]] +=
			term.coefficient * rat_t(fec_irrep_old_reflection_sign(blocks, previous));
	}
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline std::map<size_t, rat_t> fec_irrep_linear_combination(
	const std::vector<fec_raw_term_t>& left, const rat_t& left_coefficient,
	const std::vector<fec_raw_term_t>& right, const rat_t& right_coefficient) {
	std::map<size_t, rat_t> result;
	for (const auto& term : left)
		result[term.coordinate] += left_coefficient * term.coefficient;
	for (const auto& term : right)
		result[term.coordinate] += right_coefficient * term.coefficient;
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline bool validate_hhex_candidate_symmetry(
	const fec_candidate_basis_t& candidates,
	const std::vector<fec_block_t>& blocks) {
	for (const auto& candidate : candidates.vectors) {
		for (size_t component = 0; component < candidates.irrep_dimension; ++component) {
			for (const auto& term : candidate.component[component]) {
				const size_t previous = static_cast<size_t>(term.coordinate) / 9;
				const size_t letter = static_cast<size_t>(term.coordinate) % 9;
				if (fec_irrep_old_parity(blocks, previous)
					* (letter < 6 ? 1 : -1) != candidates.parity)
					return false;
			}
		}
		const auto reflection_first = fec_irrep_reflection_image(candidate.component[0], blocks);
		const rat_t first_sign = candidates.irrep == fec_irrep_t::A2
			|| candidates.irrep == fec_irrep_t::reflection_odd ? rat_t(-1) : rat_t(1);
		if (reflection_first != fec_irrep_linear_combination(
			candidate.component[0], first_sign, {}, rat_t(0)))
			return false;
		if (candidates.irrep_dimension == 1)
			continue;
		const auto reflection_second = fec_irrep_reflection_image(candidate.component[1], blocks);
		if (reflection_second != fec_irrep_linear_combination(
			candidate.component[1], rat_t(-1), {}, rat_t(0)))
			return false;
		if (fec_irrep_cyclic_image(candidate.component[0], blocks)
			!= fec_irrep_linear_combination(
				candidate.component[0], rat_t(-1, 2),
				candidate.component[1], rat_t(3, 2)))
			return false;
		if (fec_irrep_cyclic_image(candidate.component[1], blocks)
			!= fec_irrep_linear_combination(
				candidate.component[0], rat_t(-1, 2),
				candidate.component[1], rat_t(-1, 2)))
			return false;
	}
	return true;
}

inline bool validate_hhex_d3_clebsch_gordan() {
	rref_option_t opt;
	opt->verbose = false;
	opt->pool.reset(1);
	std::vector<fec_block_t> blocks;
	for (const auto& spec : hhex_block_specs(fec_symmetry_t::parity_dihedral)) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> tensor(
			{spec.irrep_dimension, 1, 9});
		blocks.push_back(make_hhex_block(
			fec_symmetry_t::parity_dihedral, spec.name,
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(tensor), &(opt->pool))));
	}
	for (const auto& spec : hhex_block_specs(fec_symmetry_t::parity_dihedral)) {
		auto candidates = make_hhex_candidates(
			fec_symmetry_t::parity_dihedral, blocks, spec.name);
		if (!validate_hhex_candidate_symmetry(candidates, blocks)) {
			std::cerr << "Invalid Hhex Clebsch-Gordan block: " << spec.name << std::endl;
			return false;
		}
	}
	return true;
}

inline sparse_mat<rat_t, int32_t> expanded_hhex_block_matrix(
	const std::vector<fec_block_t>& blocks) {
	const size_t dimension = hhex_global_dimension(blocks);
	const size_t previous_dimension = blocks[0].tensor.dim(1);
	sparse_mat<rat_t, int32_t> result(dimension, previous_dimension * 9);
	size_t offset = 0;
	for (const auto& block : blocks) {
		for (size_t row = 0; row < block.tensor.dim(0); ++row) {
			for (size_t pos = block.tensor.rowptr()[row];
				pos < block.tensor.rowptr()[row + 1]; ++pos) {
				const auto index = block.tensor.index(pos);
				result[offset + row].push_back(
					static_cast<int32_t>(index[0] * 9 + index[1]), block.tensor.val(pos));
			}
			fec_irrep_merge_entries(result[offset + row]);
		}
		offset += block.tensor.dim(0);
	}
	return result;
}

inline sparse_mat<rat_t, int32_t> expand_hhex_blocks(
	const std::vector<fec_block_t>& blocks,
	const sparse_mat<rat_t, int32_t>& previous) {
	if (blocks.empty() || blocks[0].tensor.dim(1) != previous.nrow)
		throw std::runtime_error("FEC predecessor and physical basis dimensions differ.");
	sparse_mat<rat_t, int32_t> result(
		hhex_global_dimension(blocks), previous.ncol * 9);
	size_t offset = 0;
	for (const auto& block : blocks) {
		for (size_t row = 0; row < block.tensor.dim(0); ++row) {
			for (size_t pos = block.tensor.rowptr()[row];
				pos < block.tensor.rowptr()[row + 1]; ++pos) {
				const auto index = block.tensor.index(pos);
				for (size_t previous_pos = 0;
					previous_pos < previous[index[0]].nnz(); ++previous_pos) {
					result[offset + row].push_back(
						static_cast<int32_t>(previous[index[0]](previous_pos) * 9 + index[1]),
						previous[index[0]][previous_pos] * block.tensor.val(pos));
				}
			}
			fec_irrep_merge_entries(result[offset + row]);
		}
		offset += block.tensor.dim(0);
	}
	return result;
}

inline sparse_mat<rat_t, int32_t> expanded_hhex_weight_two_matrix(
	const std::vector<fec_block_t>& blocks, thread_pool* pool) {
	const auto symmetry = blocks.size() == 4
		? fec_symmetry_t::parity_flip : fec_symmetry_t::parity_dihedral;
	auto seed = make_hhex_seed_blocks(symmetry, pool);
	return expand_hhex_blocks(blocks, expanded_hhex_block_matrix(seed));
}

inline size_t fec_irrep_matrix_rank(
	sparse_mat<rat_t, int32_t>& matrix, rref_option_t opt) {
	matrix.sort_rows_by_nnz();
	auto pivots = sparse_mat_rref_reconstruct(matrix, opt);
	size_t rank = 0;
	for (const auto& group : pivots)
		rank += group.size();
	return rank;
}

inline bool equivalent_expanded_hhex_spaces(
	const std::vector<fec_block_t>& left,
	const std::vector<fec_block_t>& right,
	const field_t&, rref_option_t opt) {
	auto left_matrix = expanded_hhex_weight_two_matrix(left, &(opt->pool));
	auto right_matrix = expanded_hhex_weight_two_matrix(right, &(opt->pool));
	if (left_matrix.ncol != right_matrix.ncol
		|| left_matrix.nrow != right_matrix.nrow)
		return false;
	auto left_copy = left_matrix;
	auto right_copy = right_matrix;
	const size_t left_rank = fec_irrep_matrix_rank(left_copy, opt);
	const size_t right_rank = fec_irrep_matrix_rank(right_copy, opt);
	if (left_rank != left_matrix.nrow || right_rank != right_matrix.nrow
		|| left_rank != right_rank) {
		std::cerr << "Expanded FEC rank mismatch: left=" << left_rank
			<< "/" << left_matrix.nrow << ", right=" << right_rank
			<< "/" << right_matrix.nrow << std::endl;
		return false;
	}
	sparse_mat<rat_t, int32_t> combined(
		left_matrix.nrow + right_matrix.nrow, left_matrix.ncol);
	for (size_t row = 0; row < left_matrix.nrow; ++row)
		combined[row] = std::move(left_matrix[row]);
	for (size_t row = 0; row < right_matrix.nrow; ++row)
		combined[left_matrix.nrow + row] = std::move(right_matrix[row]);
	const size_t combined_rank = fec_irrep_matrix_rank(combined, opt);
	if (combined_rank != left_rank)
		std::cerr << "Expanded FEC row-space mismatch: combined="
			<< combined_rank << ", expected=" << left_rank << std::endl;
	return combined_rank == left_rank;
}
