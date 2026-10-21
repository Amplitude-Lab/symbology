#pragma once

#include "dlogmat_irrep.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace SparseRREF;

struct fec_irrep_axis_block_t {
	fec_block_spec_t spec;
	size_t multiplicity = 0;
	size_t offset = 0;
	size_t dimension = 0;
};

struct fec_codomain_basis_t {
	fec_block_spec_t target;
	size_t left_dimension = 0;
	size_t right_dimension = 0;
	size_t raw_dimension = 0;
	size_t multiplicity = 0;
	std::vector<fec_candidate_t> vectors;
	sparse_mat<rat_t, int32_t> embedding;
	std::vector<size_t> pivot_coordinates;
	sparse_mat<rat_t, int32_t> transpose_inverse_chart;
	bool monomial_chart = false;
};

struct fec_reduced_extension_metrics_t {
	size_t candidate_columns = 0;
	size_t raw_matrix_rows = 0;
	size_t selected_matrix_rows = 0;
	size_t reduced_matrix_rows = 0;
	size_t selected_matrix_nnz = 0;
	size_t matrix_nnz = 0;
	size_t rank = 0;
	size_t output_nnz = 0;
	bool monomial_chart = false;
};

inline size_t fec_irrep_axis_dimension(
	const std::vector<fec_irrep_axis_block_t>& blocks) {
	size_t result = 0;
	for (const auto& block : blocks) {
		if (block.offset != result
			|| block.dimension != block.multiplicity * block.spec.irrep_dimension)
			throw std::runtime_error("Irrep axis blocks are not contiguous or consistent.");
		result += block.dimension;
	}
	return result;
}

inline fec_irrep_axis_block_t fec_irrep_axis_block(
	const dlogmat_irrep_block_t& block) {
	return {block.spec, block.multiplicity, block.offset, block.dimension};
}

inline std::vector<fec_irrep_axis_block_t> fec_irrep_condition_axis(
	const std::vector<dlogmat_irrep_block_t>& blocks) {
	std::vector<fec_irrep_axis_block_t> result;
	result.reserve(blocks.size());
	for (const auto& block : blocks)
		result.push_back(fec_irrep_axis_block(block));
	return result;
}

inline std::vector<fec_irrep_axis_block_t> fec_irrep_axis_from_blocks(
	const std::vector<fec_block_t>& blocks) {
	std::vector<fec_irrep_axis_block_t> result;
	result.reserve(blocks.size());
	size_t offset = 0;
	for (const auto& block : blocks) {
		result.push_back({
			block.spec, block.multiplicity, offset, block.tensor.dim(0)});
		offset += block.tensor.dim(0);
	}
	return result;
}

inline void fec_irrep_product_add_term(
	fec_candidate_t& candidate, size_t component,
	size_t left, size_t right, size_t right_dimension,
	const rat_t& coefficient) {
	if (coefficient == 0)
		return;
	const size_t coordinate = left * right_dimension + right;
	if (coordinate > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
		throw std::overflow_error("Irrep product coordinate exceeds int32_t.");
	for (auto& term : candidate.component[component]) {
		if (term.coordinate == static_cast<int32_t>(coordinate)) {
			term.coefficient += coefficient;
			return;
		}
	}
	candidate.component[component].push_back(
		{static_cast<int32_t>(coordinate), coefficient});
}

inline std::vector<std::pair<size_t, rat_t>> fec_irrep_axis_cyclic_image(
	const std::vector<fec_irrep_axis_block_t>& blocks, size_t coordinate) {
	for (const auto& block : blocks) {
		if (coordinate < block.offset || coordinate >= block.offset + block.dimension)
			continue;
		const size_t local = coordinate - block.offset;
		if (block.spec.irrep != fec_irrep_t::E)
			return {{coordinate, rat_t(1)}};
		if (local < block.multiplicity)
			return {{coordinate, rat_t(-1, 2)},
				{coordinate + block.multiplicity, rat_t(3, 2)}};
		return {{coordinate - block.multiplicity, rat_t(-1, 2)},
			{coordinate, rat_t(-1, 2)}};
	}
	throw std::runtime_error("Irrep axis coordinate is out of range.");
}

inline std::pair<size_t, rat_t> fec_irrep_axis_reflection_image(
	const std::vector<fec_irrep_axis_block_t>& blocks, size_t coordinate) {
	for (const auto& block : blocks) {
		if (coordinate < block.offset || coordinate >= block.offset + block.dimension)
			continue;
		const size_t local = coordinate - block.offset;
		if (block.spec.irrep == fec_irrep_t::A2
			|| (block.spec.irrep == fec_irrep_t::E
				&& local >= block.multiplicity))
			return {coordinate, rat_t(-1)};
		return {coordinate, rat_t(1)};
	}
	throw std::runtime_error("Irrep axis coordinate is out of range.");
}

inline int fec_irrep_axis_parity(
	const std::vector<fec_irrep_axis_block_t>& blocks, size_t coordinate) {
	for (const auto& block : blocks)
		if (coordinate >= block.offset && coordinate < block.offset + block.dimension)
			return block.spec.parity;
	throw std::runtime_error("Irrep axis coordinate is out of range.");
}

inline std::map<size_t, rat_t> fec_irrep_product_cyclic_image(
	const std::vector<fec_raw_term_t>& component,
	const std::vector<fec_irrep_axis_block_t>& left,
	const std::vector<fec_irrep_axis_block_t>& right,
	size_t right_dimension) {
	std::map<size_t, rat_t> result;
	for (const auto& term : component) {
		const size_t left_coordinate = static_cast<size_t>(term.coordinate) / right_dimension;
		const size_t right_coordinate = static_cast<size_t>(term.coordinate) % right_dimension;
		for (const auto& [left_image, left_coefficient] :
			fec_irrep_axis_cyclic_image(left, left_coordinate))
			for (const auto& [right_image, right_coefficient] :
				fec_irrep_axis_cyclic_image(right, right_coordinate))
				result[left_image * right_dimension + right_image] +=
					term.coefficient * left_coefficient * right_coefficient;
	}
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline std::map<size_t, rat_t> fec_irrep_product_reflection_image(
	const std::vector<fec_raw_term_t>& component,
	const std::vector<fec_irrep_axis_block_t>& left,
	const std::vector<fec_irrep_axis_block_t>& right,
	size_t right_dimension) {
	std::map<size_t, rat_t> result;
	for (const auto& term : component) {
		const size_t left_coordinate = static_cast<size_t>(term.coordinate) / right_dimension;
		const size_t right_coordinate = static_cast<size_t>(term.coordinate) % right_dimension;
		const auto [left_image, left_coefficient] =
			fec_irrep_axis_reflection_image(left, left_coordinate);
		const auto [right_image, right_coefficient] =
			fec_irrep_axis_reflection_image(right, right_coordinate);
		result[left_image * right_dimension + right_image] +=
			term.coefficient * left_coefficient * right_coefficient;
	}
	std::erase_if(result, [](const auto& item) { return item.second == 0; });
	return result;
}

inline bool validate_hhex_codomain_candidates(
	const fec_codomain_basis_t& basis,
	const std::vector<fec_irrep_axis_block_t>& left,
	const std::vector<fec_irrep_axis_block_t>& right) {
	for (const auto& candidate : basis.vectors) {
		for (size_t component = 0; component < basis.target.irrep_dimension; ++component) {
			for (const auto& term : candidate.component[component]) {
				const size_t left_coordinate =
					static_cast<size_t>(term.coordinate) / basis.right_dimension;
				const size_t right_coordinate =
					static_cast<size_t>(term.coordinate) % basis.right_dimension;
				if (fec_irrep_axis_parity(left, left_coordinate)
					* fec_irrep_axis_parity(right, right_coordinate)
					!= basis.target.parity)
					return false;
			}
		}
		const auto reflection_first = fec_irrep_product_reflection_image(
			candidate.component[0], left, right, basis.right_dimension);
		const rat_t first_sign = basis.target.irrep == fec_irrep_t::A2
			? rat_t(-1) : rat_t(1);
		if (reflection_first != fec_irrep_linear_combination(
			candidate.component[0], first_sign, {}, rat_t(0)))
			return false;
		if (basis.target.irrep_dimension == 1)
			continue;
		const auto reflection_second = fec_irrep_product_reflection_image(
			candidate.component[1], left, right, basis.right_dimension);
		if (reflection_second != fec_irrep_linear_combination(
			candidate.component[1], rat_t(-1), {}, rat_t(0)))
			return false;
		if (fec_irrep_product_cyclic_image(
			candidate.component[0], left, right, basis.right_dimension)
			!= fec_irrep_linear_combination(
				candidate.component[0], rat_t(-1, 2),
				candidate.component[1], rat_t(3, 2)))
			return false;
		if (fec_irrep_product_cyclic_image(
			candidate.component[1], left, right, basis.right_dimension)
			!= fec_irrep_linear_combination(
				candidate.component[0], rat_t(-1, 2),
				candidate.component[1], rat_t(-1, 2)))
			return false;
	}
	return true;
}

inline void fec_irrep_add_product_candidate(
	std::vector<fec_candidate_t>& result,
	const fec_irrep_axis_block_t& left,
	const fec_irrep_axis_block_t& right,
	const fec_block_spec_t& target,
	size_t left_multiplicity, size_t right_multiplicity,
	size_t right_dimension) {
	const size_t left_single = left.offset + left_multiplicity;
	const size_t left_plus = left.offset + left_multiplicity;
	const size_t left_minus = left.offset + left.multiplicity + left_multiplicity;
	const size_t right_single = right.offset + right_multiplicity;
	const size_t right_plus = right.offset + right_multiplicity;
	const size_t right_minus = right.offset + right.multiplicity + right_multiplicity;
	fec_candidate_t candidate;
	if (target.irrep == fec_irrep_t::A1) {
		if ((left.spec.irrep == fec_irrep_t::A1 && right.spec.irrep == fec_irrep_t::A1)
			|| (left.spec.irrep == fec_irrep_t::A2 && right.spec.irrep == fec_irrep_t::A2)) {
			fec_irrep_product_add_term(
				candidate, 0, left_single, right_single, right_dimension, rat_t(1));
		} else if (left.spec.irrep == fec_irrep_t::E
			&& right.spec.irrep == fec_irrep_t::E) {
			fec_irrep_product_add_term(
				candidate, 0, left_plus, right_plus, right_dimension, rat_t(1));
			fec_irrep_product_add_term(
				candidate, 0, left_minus, right_minus, right_dimension, rat_t(3));
		} else
			return;
	} else if (target.irrep == fec_irrep_t::A2) {
		if ((left.spec.irrep == fec_irrep_t::A1 && right.spec.irrep == fec_irrep_t::A2)
			|| (left.spec.irrep == fec_irrep_t::A2 && right.spec.irrep == fec_irrep_t::A1)) {
			fec_irrep_product_add_term(
				candidate, 0, left_single, right_single, right_dimension, rat_t(1));
		} else if (left.spec.irrep == fec_irrep_t::E
			&& right.spec.irrep == fec_irrep_t::E) {
			fec_irrep_product_add_term(
				candidate, 0, left_plus, right_minus, right_dimension, rat_t(-1));
			fec_irrep_product_add_term(
				candidate, 0, left_minus, right_plus, right_dimension, rat_t(1));
		} else
			return;
	} else if (left.spec.irrep == fec_irrep_t::A1
		&& right.spec.irrep == fec_irrep_t::E) {
		fec_irrep_product_add_term(
			candidate, 0, left_single, right_plus, right_dimension, rat_t(1));
		fec_irrep_product_add_term(
			candidate, 1, left_single, right_minus, right_dimension, rat_t(1));
	} else if (left.spec.irrep == fec_irrep_t::E
		&& right.spec.irrep == fec_irrep_t::A1) {
		fec_irrep_product_add_term(
			candidate, 0, left_plus, right_single, right_dimension, rat_t(1));
		fec_irrep_product_add_term(
			candidate, 1, left_minus, right_single, right_dimension, rat_t(1));
	} else if (left.spec.irrep == fec_irrep_t::A2
		&& right.spec.irrep == fec_irrep_t::E) {
		fec_irrep_product_add_term(
			candidate, 0, left_single, right_minus, right_dimension, rat_t(3));
		fec_irrep_product_add_term(
			candidate, 1, left_single, right_plus, right_dimension, rat_t(-1));
	} else if (left.spec.irrep == fec_irrep_t::E
		&& right.spec.irrep == fec_irrep_t::A2) {
		fec_irrep_product_add_term(
			candidate, 0, left_minus, right_single, right_dimension, rat_t(3));
		fec_irrep_product_add_term(
			candidate, 1, left_plus, right_single, right_dimension, rat_t(-1));
	} else if (left.spec.irrep == fec_irrep_t::E
		&& right.spec.irrep == fec_irrep_t::E) {
		fec_irrep_product_add_term(
			candidate, 0, left_plus, right_plus, right_dimension, rat_t(-1));
		fec_irrep_product_add_term(
			candidate, 0, left_minus, right_minus, right_dimension, rat_t(3));
		fec_irrep_product_add_term(
			candidate, 1, left_plus, right_minus, right_dimension, rat_t(1));
		fec_irrep_product_add_term(
			candidate, 1, left_minus, right_plus, right_dimension, rat_t(1));
	} else
		return;
	result.push_back(std::move(candidate));
}

inline fec_codomain_basis_t make_hhex_codomain_basis(
	const std::vector<fec_irrep_axis_block_t>& left,
	const std::vector<fec_irrep_axis_block_t>& right,
	const fec_block_spec_t& target,
	const field_t& F, rref_option_t opt) {
	fec_codomain_basis_t result;
	result.target = target;
	result.left_dimension = fec_irrep_axis_dimension(left);
	result.right_dimension = fec_irrep_axis_dimension(right);
	result.raw_dimension = result.left_dimension * result.right_dimension;
	if (result.raw_dimension > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
		throw std::overflow_error("Codomain dimension exceeds int32_t.");
	for (const auto& left_block : left) {
		for (const auto& right_block : right) {
			if (left_block.spec.parity * right_block.spec.parity != target.parity)
				continue;
			for (size_t i = 0; i < left_block.multiplicity; ++i)
				for (size_t j = 0; j < right_block.multiplicity; ++j)
					fec_irrep_add_product_candidate(
						result.vectors, left_block, right_block, target,
						i, j, result.right_dimension);
		}
	}
	result.multiplicity = result.vectors.size();
	result.embedding = sparse_mat<rat_t, int32_t>(
		result.multiplicity, result.raw_dimension);
	for (size_t row = 0; row < result.vectors.size(); ++row) {
		for (const auto& term : result.vectors[row].component[0])
			result.embedding[row].push_back(term.coordinate, term.coefficient);
		fec_irrep_merge_entries(result.embedding[row]);
	}

	std::vector<size_t> support_count(result.raw_dimension, 0);
	for (size_t row = 0; row < result.embedding.nrow; ++row)
		for (size_t pos = 0; pos < result.embedding[row].nnz(); ++pos)
			++support_count[result.embedding[row](pos)];
	result.pivot_coordinates.resize(result.multiplicity);
	result.monomial_chart = true;
	for (size_t row = 0; row < result.multiplicity; ++row) {
		auto position = std::find_if(
			result.embedding[row].indices,
			result.embedding[row].indices + result.embedding[row].nnz(),
			[&](int32_t coordinate) { return support_count[coordinate] == 1; });
		if (position == result.embedding[row].indices + result.embedding[row].nnz()) {
			result.monomial_chart = false;
			break;
		}
		result.pivot_coordinates[row] = static_cast<size_t>(*position);
	}
	if (!result.monomial_chart)
		result.pivot_coordinates = dlogmat_irrep_pivot_columns(result.embedding, opt);
	auto selected = dlogmat_irrep_select_columns(
		result.embedding, result.pivot_coordinates);
	if (result.monomial_chart) {
		result.transpose_inverse_chart = sparse_mat<rat_t, int32_t>(
			result.multiplicity, result.multiplicity);
		for (size_t row = 0; row < result.multiplicity; ++row) {
			const auto value = *selected[row].find(static_cast<int32_t>(row));
			result.transpose_inverse_chart[row].push_back(
				static_cast<int32_t>(row), rat_t(1) / value);
		}
	} else {
		auto inverse = sparse_mat_inverse(selected, F, opt);
		if (inverse.nrow != result.multiplicity || inverse.ncol != result.multiplicity)
			throw std::runtime_error("Codomain coordinate chart is not invertible.");
		result.transpose_inverse_chart = inverse.transpose();
	}
	return result;
}

inline bool validate_hhex_codomain_basis(
	const fec_codomain_basis_t& basis,
	const std::vector<fec_irrep_axis_block_t>& left,
	const std::vector<fec_irrep_axis_block_t>& right,
	const field_t& F, rref_option_t opt) {
	if (basis.left_dimension != fec_irrep_axis_dimension(left)
		|| basis.right_dimension != fec_irrep_axis_dimension(right)
		|| basis.raw_dimension != basis.left_dimension * basis.right_dimension
		|| basis.embedding.nrow != basis.multiplicity
		|| basis.embedding.ncol != basis.raw_dimension
		|| basis.pivot_coordinates.size() != basis.multiplicity
		|| basis.transpose_inverse_chart.nrow != basis.multiplicity
		|| basis.transpose_inverse_chart.ncol != basis.multiplicity
		|| !validate_hhex_codomain_candidates(basis, left, right))
		return false;
	auto selected = dlogmat_irrep_select_columns(
		basis.embedding, basis.pivot_coordinates);
	auto identity = sparse_mat_mul(
		basis.transpose_inverse_chart, selected.transpose(), F, &(opt->pool));
	return dlogmat_irrep_matrices_equal(
		identity, dlogmat_irrep_identity(basis.multiplicity));
}

inline fec_block_t extend_hhex_block_reduced(
	const sparse_tensor<rat_t, int32_t, SPARSE_CSR>& dlogmat,
	const std::vector<fec_irrep_axis_block_t>& condition_blocks,
	const std::vector<fec_irrep_axis_block_t>& row_blocks,
	const std::vector<fec_block_t>& blocks, const std::string& block_name,
	const field_t& F, rref_option_t& opt,
	fec_reduced_extension_metrics_t* metrics = nullptr) {
	hhex_validate_blocks(fec_symmetry_t::parity_dihedral, blocks);
	if (dlogmat.rank() != 3 || dlogmat.dim(0) != 9 || dlogmat.dim(1) != 9)
		throw std::runtime_error("Adapted condition tensor must have dimensions {9,9,q}.");
	if (fec_irrep_axis_dimension(condition_blocks) != dlogmat.dim(2))
		throw std::runtime_error("Condition block metadata does not match the dlogmat third axis.");
	if (fec_irrep_axis_dimension(row_blocks) != blocks[0].tensor.dim(1))
		throw std::runtime_error("FEC row-axis metadata does not match the predecessor second axis.");
	auto target = std::find_if(
		hhex_block_specs(fec_symmetry_t::parity_dihedral).begin(),
		hhex_block_specs(fec_symmetry_t::parity_dihedral).end(),
		[&](const auto& spec) { return spec.name == block_name; });
	if (target == hhex_block_specs(fec_symmetry_t::parity_dihedral).end())
		throw std::runtime_error("Unknown parity-dihedral FEC block: " + block_name + ".");
	auto candidates = make_hhex_candidates(
		fec_symmetry_t::parity_dihedral, blocks, block_name);
	if (!validate_hhex_candidate_symmetry(candidates, blocks))
		throw std::runtime_error("FEC candidate basis violates its symmetry label.");
	auto codomain = make_hhex_codomain_basis(
		row_blocks, condition_blocks, *target, F, opt);
	if (!validate_hhex_codomain_basis(codomain, row_blocks, condition_blocks, F, opt))
		throw std::runtime_error("FEC codomain CG basis failed exact validation.");

	fec_reduced_extension_metrics_t local_metrics;
	local_metrics.candidate_columns = candidates.size();
	local_metrics.raw_matrix_rows = blocks[0].tensor.dim(1) * dlogmat.dim(2);
	local_metrics.selected_matrix_rows = codomain.pivot_coordinates.size();
	local_metrics.reduced_matrix_rows = codomain.multiplicity;
	local_metrics.monomial_chart = codomain.monomial_chart;
	if (candidates.size() == 0) {
		sparse_tensor<rat_t, int32_t, SPARSE_COO> empty(
			{0, candidates.previous_dimension, 9});
		if (metrics != nullptr)
			*metrics = local_metrics;
		return make_hhex_block(
			fec_symmetry_t::parity_dihedral, block_name,
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(
				std::move(empty), &(opt->pool)));
	}

	std::vector<int32_t> selected_row(local_metrics.raw_matrix_rows, -1);
	for (size_t row = 0; row < codomain.pivot_coordinates.size(); ++row) {
		const size_t coordinate = codomain.pivot_coordinates[row];
		if (coordinate >= selected_row.size() || selected_row[coordinate] != -1)
			throw std::runtime_error("Invalid or repeated codomain pivot coordinate.");
		selected_row[coordinate] = static_cast<int32_t>(row);
	}
	struct dlog_entry_t {
		int32_t letter;
		int32_t relation;
		rat_t coefficient;
	};
	std::array<std::vector<dlog_entry_t>, 9> dlog_lookup;
	for (size_t letter = 0; letter < 9; ++letter) {
		for (size_t pos = dlogmat.rowptr()[letter];
			pos < dlogmat.rowptr()[letter + 1]; ++pos) {
			const auto index = dlogmat.index(pos);
			dlog_lookup[letter].push_back(
				{index[0], index[1], dlogmat.val(pos)});
		}
	}

	sparse_mat<rat_t, int32_t> selected_matrix(
		local_metrics.selected_matrix_rows, candidates.size());
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
					const size_t raw_row = index[0] * dlogmat.dim(2)
						+ dlog_entry.relation;
					const int32_t row = selected_row[raw_row];
					if (row < 0)
						continue;
					const auto& lookup = candidates.first_component_lookup[
						(offset + previous) * 9 + dlog_entry.letter];
					for (const auto& candidate_entry : lookup)
						selected_matrix[row].push_back(
							candidate_entry.column,
							block.tensor.val(pos) * dlog_entry.coefficient
								* candidate_entry.coefficient);
				}
			}
		}
		offset += block.tensor.dim(0);
	}
	opt->pool.detach_loop(0, selected_matrix.nrow, [&](size_t row) {
		fec_irrep_merge_entries(selected_matrix[row]);
	});
	opt->pool.wait();
	local_metrics.selected_matrix_nnz = selected_matrix.nnz();
	auto matrix = sparse_mat_mul(
		codomain.transpose_inverse_chart, selected_matrix, F, &(opt->pool));
	selected_matrix.clear();
	opt->pool.detach_loop(0, matrix.nrow, [&](size_t row) {
		fec_irrep_merge_entries(matrix[row]);
	});
	opt->pool.wait();
	local_metrics.matrix_nnz = matrix.nnz();
	matrix.sort_rows_by_nnz();

	sparse_mat<rat_t, int32_t> kernel;
	if (matrix.nrow == 0) {
		kernel = sparse_mat<rat_t, int32_t>(candidates.size(), candidates.size());
		for (size_t row = 0; row < candidates.size(); ++row)
			kernel[row].push_back(static_cast<int32_t>(row), rat_t(1));
	} else {
		auto pivots = sparse_mat_rref_reconstruct(matrix, opt);
		for (const auto& group : pivots)
			local_metrics.rank += group.size();
		if (local_metrics.rank == candidates.size())
			kernel = sparse_mat<rat_t, int32_t>(0, candidates.size());
		else
			kernel = sparse_mat_rref_kernel(matrix, pivots, F, opt).transpose();
	}
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
		return make_hhex_block(
			fec_symmetry_t::parity_dihedral, block_name,
			sparse_tensor<rat_t, int32_t, SPARSE_CSR>(
				std::move(empty), &(opt->pool)));
	}

	sparse_mat<rat_t, int32_t> expanded(
		kernel.nrow * candidates.irrep_dimension,
		candidates.previous_dimension * 9);
	for (size_t component = 0; component < candidates.irrep_dimension; ++component) {
		for (size_t row = 0; row < kernel.nrow; ++row) {
			auto& output_row = expanded[component * kernel.nrow + row];
			for (size_t pos = 0; pos < kernel[row].nnz(); ++pos) {
				const auto& candidate = candidates.vectors[kernel[row](pos)];
				for (const auto& term : candidate.component[component])
					output_row.push_back(
						term.coordinate, kernel[row][pos] * term.coefficient);
			}
			fec_irrep_merge_entries(output_row);
		}
	}
	local_metrics.output_nnz = expanded.nnz();
	sparse_tensor<rat_t, int32_t, SPARSE_COO> output(expanded, &(opt->pool));
	output.reshape({expanded.nrow, candidates.previous_dimension, 9});
	if (metrics != nullptr)
		*metrics = local_metrics;
	return make_hhex_block(
		fec_symmetry_t::parity_dihedral, block_name,
		sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(output), &(opt->pool)));
}
