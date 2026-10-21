#pragma once

#include "fec_irrep.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace SparseRREF;

enum class hhex_dlogmat_action_t {
	identity,
	cyclic,
	reflection,
	parity
};

struct dlogmat_irrep_block_t {
	fec_block_spec_t spec;
	size_t multiplicity = 0;
	size_t offset = 0;
	size_t dimension = 0;
};

struct dlogmat_irrep_decomposition_t {
	sparse_tensor<rat_t, int32_t, SPARSE_CSR> tensor;
	sparse_mat<rat_t, int32_t> basis;
	std::vector<dlogmat_irrep_block_t> blocks;
	sparse_mat<rat_t, int32_t> source_rows;
	sparse_mat<rat_t, int32_t> source_cyclic;
	sparse_mat<rat_t, int32_t> source_reflection;
	sparse_mat<rat_t, int32_t> source_parity;
	sparse_mat<rat_t, int32_t> adapted_cyclic;
	sparse_mat<rat_t, int32_t> adapted_reflection;
	sparse_mat<rat_t, int32_t> adapted_parity;
};

inline bool dlogmat_irrep_matrices_equal(
	sparse_mat<rat_t, int32_t> left,
	sparse_mat<rat_t, int32_t> right) {
	if (left.nrow != right.nrow || left.ncol != right.ncol)
		return false;
	for (size_t row = 0; row < left.nrow; ++row) {
		fec_irrep_merge_entries(left[row]);
		fec_irrep_merge_entries(right[row]);
		if (left[row].nnz() != right[row].nnz())
			return false;
		for (size_t pos = 0; pos < left[row].nnz(); ++pos) {
			if (left[row](pos) != right[row](pos)
				|| left[row][pos] != right[row][pos])
				return false;
		}
	}
	return true;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_identity(size_t dimension) {
	sparse_mat<rat_t, int32_t> result(dimension, dimension);
	for (size_t i = 0; i < dimension; ++i)
		result[i].push_back(static_cast<int32_t>(i), rat_t(1));
	return result;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_linear_combination(
	const std::vector<std::pair<const sparse_mat<rat_t, int32_t>*, rat_t>>& terms) {
	if (terms.empty())
		return {};
	const size_t nrow = terms[0].first->nrow;
	const size_t ncol = terms[0].first->ncol;
	sparse_mat<rat_t, int32_t> result(nrow, ncol);
	for (const auto& [matrix, coefficient] : terms) {
		if (matrix->nrow != nrow || matrix->ncol != ncol)
			throw std::runtime_error("Cannot combine matrices with different dimensions.");
		if (coefficient == 0)
			continue;
		for (size_t row = 0; row < nrow; ++row)
			for (size_t pos = 0; pos < (*matrix)[row].nnz(); ++pos)
				result[row].push_back(
					(*matrix)[row](pos), coefficient * (*matrix)[row][pos]);
	}
	for (size_t row = 0; row < nrow; ++row)
		fec_irrep_merge_entries(result[row]);
	return result;
}

inline sparse_mat<rat_t, int32_t> hhex_dlogmat_rows(
	const sparse_tensor<rat_t, int32_t, SPARSE_CSR>& dlogmat,
	hhex_dlogmat_action_t action) {
	if (dlogmat.rank() != 3 || dlogmat.dim(0) != 9 || dlogmat.dim(1) != 9)
		throw std::runtime_error("Hhex condition tensor must have dimensions {9,9,q}.");
	const std::array<size_t, 9> reflection{2, 1, 0, 5, 4, 3, 8, 7, 6};
	sparse_mat<rat_t, int32_t> result(dlogmat.dim(2), 81);
	for (size_t left = 0; left < 9; ++left) {
		for (size_t pos = dlogmat.rowptr()[left]; pos < dlogmat.rowptr()[left + 1]; ++pos) {
			const auto index = dlogmat.index(pos);
			const size_t right = index[0];
			size_t image_left = left;
			size_t image_right = right;
			rat_t coefficient = dlogmat.val(pos);
			if (action == hhex_dlogmat_action_t::cyclic) {
				image_left = 3 * (left / 3) + (left + 1) % 3;
				image_right = 3 * (right / 3) + (right + 1) % 3;
			} else if (action == hhex_dlogmat_action_t::reflection) {
				image_left = reflection[left];
				image_right = reflection[right];
			} else if (action == hhex_dlogmat_action_t::parity) {
				coefficient *= rat_t(
					(left < 6 ? 1 : -1) * (right < 6 ? 1 : -1));
			}
			result[index[1]].push_back(
				static_cast<int32_t>(9 * image_left + image_right), coefficient);
		}
	}
	for (size_t row = 0; row < result.nrow; ++row)
		fec_irrep_merge_entries(result[row]);
	return result;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_select_columns(
	const sparse_mat<rat_t, int32_t>& matrix,
	const std::vector<size_t>& columns) {
	std::vector<int32_t> positions(matrix.ncol, -1);
	for (size_t i = 0; i < columns.size(); ++i) {
		if (columns[i] >= matrix.ncol)
			throw std::runtime_error("Selected dlogmat column is out of range.");
		positions[columns[i]] = static_cast<int32_t>(i);
	}
	sparse_mat<rat_t, int32_t> result(matrix.nrow, columns.size());
	for (size_t row = 0; row < matrix.nrow; ++row) {
		for (size_t pos = 0; pos < matrix[row].nnz(); ++pos) {
			const auto target = positions[matrix[row](pos)];
			if (target >= 0)
				result[row].push_back(target, matrix[row][pos]);
		}
	}
	return result;
}

inline std::vector<size_t> dlogmat_irrep_pivot_columns(
	const sparse_mat<rat_t, int32_t>& matrix, rref_option_t opt) {
	auto reduced = matrix;
	reduced.sort_rows_by_nnz();
	auto pivot_groups = sparse_mat_rref_reconstruct(reduced, opt, true);
	std::vector<size_t> result;
	for (const auto& group : pivot_groups)
		for (const auto& pivot : group)
			result.push_back(static_cast<size_t>(pivot.c));
	std::ranges::sort(result);
	result.erase(std::unique(result.begin(), result.end()), result.end());
	if (result.size() != matrix.nrow)
		throw std::runtime_error("Condition slices are not linearly independent.");
	return result;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_recover_action(
	const sparse_mat<rat_t, int32_t>& basis,
	const sparse_mat<rat_t, int32_t>& transformed,
	const field_t& F, rref_option_t opt) {
	if (basis.nrow != transformed.nrow || basis.ncol != transformed.ncol)
		throw std::runtime_error("Condition basis and transformed basis dimensions differ.");
	const auto pivots = dlogmat_irrep_pivot_columns(basis, opt);
	auto pivot_block = dlogmat_irrep_select_columns(basis, pivots);
	auto transformed_block = dlogmat_irrep_select_columns(transformed, pivots);
	auto inverse = sparse_mat_inverse(pivot_block, F, opt);
	if (inverse.nrow != basis.nrow || inverse.ncol != basis.nrow)
		throw std::runtime_error("Condition pivot block is not invertible.");
	auto action = sparse_mat_mul(transformed_block, inverse, F, &(opt->pool));
	auto reconstructed = sparse_mat_mul(action, basis, F, &(opt->pool));
	if (!dlogmat_irrep_matrices_equal(reconstructed, transformed))
		throw std::runtime_error("Recovered action does not reconstruct transformed dlogmat rows.");
	return action;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_row_basis(
	sparse_mat<rat_t, int32_t> matrix, rref_option_t opt) {
	matrix.sort_rows_by_nnz();
	sparse_mat_rref_reconstruct(matrix, opt);
	std::vector<size_t> rows;
	for (size_t row = 0; row < matrix.nrow; ++row) {
		if (matrix[row].nnz() == 0)
			continue;
		fec_irrep_cancel_divisor(matrix[row]);
		rows.push_back(row);
	}
	std::ranges::stable_sort(rows, [&](size_t left, size_t right) {
		return matrix[left](0) < matrix[right](0);
	});
	sparse_mat<rat_t, int32_t> result(rows.size(), matrix.ncol);
	for (size_t i = 0; i < rows.size(); ++i)
		result[i] = std::move(matrix[rows[i]]);
	return result;
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_projector_basis(
	int parity, fec_irrep_t irrep,
	const sparse_mat<rat_t, int32_t>& cyclic,
	const sparse_mat<rat_t, int32_t>& reflection,
	const sparse_mat<rat_t, int32_t>& parity_action,
	const field_t& F, rref_option_t opt) {
	const size_t dimension = cyclic.nrow;
	auto identity = dlogmat_irrep_identity(dimension);
	auto cyclic2 = sparse_mat_mul(cyclic, cyclic, F, &(opt->pool));
	auto parity_projector = dlogmat_irrep_linear_combination({
		{&identity, rat_t(1, 2)}, {&parity_action, rat_t(parity, 2)}});
	auto cyclic_fixed = dlogmat_irrep_linear_combination({
		{&identity, rat_t(1, 3)}, {&cyclic, rat_t(1, 3)},
		{&cyclic2, rat_t(1, 3)}});
	auto cyclic_E = dlogmat_irrep_linear_combination({
		{&identity, rat_t(1)}, {&cyclic_fixed, rat_t(-1)}});
	const int reflection_sign = irrep == fec_irrep_t::A2 ? -1 : 1;
	auto reflection_projector = dlogmat_irrep_linear_combination({
		{&identity, rat_t(1, 2)},
		{&reflection, rat_t(reflection_sign, 2)}});
	auto d3_projector = irrep == fec_irrep_t::E ? cyclic_E : cyclic_fixed;
	auto combined = sparse_mat_mul(
		parity_projector, d3_projector, F, &(opt->pool));
	combined = sparse_mat_mul(combined, reflection_projector, F, &(opt->pool));
	return dlogmat_irrep_row_basis(std::move(combined), opt);
}

inline sparse_mat<rat_t, int32_t> dlogmat_irrep_expected_action(
	const std::vector<dlogmat_irrep_block_t>& blocks,
	hhex_dlogmat_action_t action) {
	size_t dimension = 0;
	for (const auto& block : blocks)
		dimension += block.dimension;
	sparse_mat<rat_t, int32_t> result(dimension, dimension);
	for (const auto& block : blocks) {
		if (block.spec.irrep != fec_irrep_t::E) {
			int sign = 1;
			if (action == hhex_dlogmat_action_t::reflection
				&& block.spec.irrep == fec_irrep_t::A2)
				sign = -1;
			if (action == hhex_dlogmat_action_t::parity)
				sign = block.spec.parity;
			for (size_t i = 0; i < block.dimension; ++i)
				result[block.offset + i].push_back(
					static_cast<int32_t>(block.offset + i), rat_t(sign));
			continue;
		}
		for (size_t i = 0; i < block.multiplicity; ++i) {
			const size_t plus = block.offset + i;
			const size_t minus = block.offset + block.multiplicity + i;
			if (action == hhex_dlogmat_action_t::cyclic) {
				result[plus].push_back(static_cast<int32_t>(plus), rat_t(-1, 2));
				result[plus].push_back(static_cast<int32_t>(minus), rat_t(3, 2));
				result[minus].push_back(static_cast<int32_t>(plus), rat_t(-1, 2));
				result[minus].push_back(static_cast<int32_t>(minus), rat_t(-1, 2));
			} else {
				const int plus_sign = action == hhex_dlogmat_action_t::parity
					? block.spec.parity : 1;
				const int minus_sign = action == hhex_dlogmat_action_t::reflection
					? -1 : plus_sign;
				result[plus].push_back(static_cast<int32_t>(plus), rat_t(plus_sign));
				result[minus].push_back(static_cast<int32_t>(minus), rat_t(minus_sign));
			}
		}
	}
	return result;
}

inline sparse_tensor<rat_t, int32_t, SPARSE_CSR> dlogmat_irrep_tensor_from_rows(
	const sparse_mat<rat_t, int32_t>& rows, thread_pool* pool) {
	if (rows.ncol != 81)
		throw std::runtime_error("Adapted dlogmat rows must have 81 columns.");
	sparse_tensor<rat_t, int32_t, SPARSE_COO> tensor({9, 9, rows.nrow});
	tensor.reserve(rows.nnz());
	for (size_t relation = 0; relation < rows.nrow; ++relation) {
		for (size_t pos = 0; pos < rows[relation].nnz(); ++pos) {
			const size_t pair = rows[relation](pos);
			tensor.push_back({static_cast<int32_t>(pair / 9),
				static_cast<int32_t>(pair % 9), static_cast<int32_t>(relation)},
				rows[relation][pos]);
		}
	}
	tensor.sort_indices(pool);
	return sparse_tensor<rat_t, int32_t, SPARSE_CSR>(std::move(tensor), pool);
}

inline bool validate_hhex_dlogmat_decomposition(
	const dlogmat_irrep_decomposition_t& decomposition,
	const field_t& F, rref_option_t opt) {
	const size_t dimension = decomposition.source_rows.nrow;
	if (decomposition.source_rows.ncol != 81
		|| decomposition.basis.nrow != dimension
		|| decomposition.basis.ncol != dimension
		|| decomposition.tensor.rank() != 3
		|| decomposition.tensor.dim(0) != 9
		|| decomposition.tensor.dim(1) != 9
		|| decomposition.tensor.dim(2) != dimension)
		return false;
	auto identity = dlogmat_irrep_identity(dimension);
	auto cyclic2 = sparse_mat_mul(
		decomposition.source_cyclic, decomposition.source_cyclic, F, &(opt->pool));
	auto cyclic3 = sparse_mat_mul(
		cyclic2, decomposition.source_cyclic, F, &(opt->pool));
	auto reflection2 = sparse_mat_mul(
		decomposition.source_reflection, decomposition.source_reflection,
		F, &(opt->pool));
	auto parity2 = sparse_mat_mul(
		decomposition.source_parity, decomposition.source_parity, F, &(opt->pool));
	auto reflection_cyclic = sparse_mat_mul(
		decomposition.source_reflection, decomposition.source_cyclic,
		F, &(opt->pool));
	auto reflection_cyclic_reflection = sparse_mat_mul(
		reflection_cyclic, decomposition.source_reflection, F, &(opt->pool));
	auto parity_cyclic = sparse_mat_mul(
		decomposition.source_parity, decomposition.source_cyclic, F, &(opt->pool));
	auto cyclic_parity = sparse_mat_mul(
		decomposition.source_cyclic, decomposition.source_parity, F, &(opt->pool));
	auto parity_reflection = sparse_mat_mul(
		decomposition.source_parity, decomposition.source_reflection, F, &(opt->pool));
	auto reflection_parity = sparse_mat_mul(
		decomposition.source_reflection, decomposition.source_parity, F, &(opt->pool));
	if (!dlogmat_irrep_matrices_equal(cyclic3, identity)
		|| !dlogmat_irrep_matrices_equal(reflection2, identity)
		|| !dlogmat_irrep_matrices_equal(parity2, identity)
		|| !dlogmat_irrep_matrices_equal(reflection_cyclic_reflection, cyclic2)
		|| !dlogmat_irrep_matrices_equal(parity_cyclic, cyclic_parity)
		|| !dlogmat_irrep_matrices_equal(parity_reflection, reflection_parity))
		return false;
	auto expected_cyclic = dlogmat_irrep_expected_action(
		decomposition.blocks, hhex_dlogmat_action_t::cyclic);
	auto expected_reflection = dlogmat_irrep_expected_action(
		decomposition.blocks, hhex_dlogmat_action_t::reflection);
	auto expected_parity = dlogmat_irrep_expected_action(
		decomposition.blocks, hhex_dlogmat_action_t::parity);
	if (!dlogmat_irrep_matrices_equal(decomposition.adapted_cyclic, expected_cyclic)
		|| !dlogmat_irrep_matrices_equal(
			decomposition.adapted_reflection, expected_reflection)
		|| !dlogmat_irrep_matrices_equal(decomposition.adapted_parity, expected_parity))
		return false;
	auto adapted_rows = hhex_dlogmat_rows(
		decomposition.tensor, hhex_dlogmat_action_t::identity);
	auto reconstructed = sparse_mat_mul(
		decomposition.basis, decomposition.source_rows, F, &(opt->pool));
	return dlogmat_irrep_matrices_equal(adapted_rows, reconstructed);
}

inline dlogmat_irrep_decomposition_t decompose_hhex_dlogmat(
	const sparse_tensor<rat_t, int32_t, SPARSE_CSR>& dlogmat,
	const field_t& F, rref_option_t opt) {
	dlogmat_irrep_decomposition_t result;
	result.source_rows = hhex_dlogmat_rows(
		dlogmat, hhex_dlogmat_action_t::identity);
	result.source_cyclic = dlogmat_irrep_recover_action(
		result.source_rows,
		hhex_dlogmat_rows(dlogmat, hhex_dlogmat_action_t::cyclic), F, opt);
	result.source_reflection = dlogmat_irrep_recover_action(
		result.source_rows,
		hhex_dlogmat_rows(dlogmat, hhex_dlogmat_action_t::reflection), F, opt);
	result.source_parity = dlogmat_irrep_recover_action(
		result.source_rows,
		hhex_dlogmat_rows(dlogmat, hhex_dlogmat_action_t::parity), F, opt);

	for (const auto& spec : hhex_block_specs(fec_symmetry_t::parity_dihedral)) {
		auto basis = dlogmat_irrep_projector_basis(
			spec.parity, spec.irrep, result.source_cyclic,
			result.source_reflection, result.source_parity, F, opt);
		if (spec.irrep == fec_irrep_t::E) {
			auto cyclic_image = sparse_mat_mul(
				basis, result.source_cyclic, F, &(opt->pool));
			auto paired = dlogmat_irrep_linear_combination({
				{&cyclic_image, rat_t(2, 3)}, {&basis, rat_t(1, 3)}});
			const size_t offset = result.basis.nrow;
			result.blocks.push_back({spec, basis.nrow, offset, 2 * basis.nrow});
			result.basis.append(std::move(basis));
			result.basis.append(std::move(paired));
		} else {
			const size_t offset = result.basis.nrow;
			result.blocks.push_back({spec, basis.nrow, offset, basis.nrow});
			result.basis.append(std::move(basis));
		}
	}
	if (result.basis.nrow != result.source_rows.nrow
		|| result.basis.ncol != result.source_rows.nrow)
		throw std::runtime_error("Parity-dihedral dlogmat decomposition is incomplete.");
	auto inverse = sparse_mat_inverse(result.basis, F, opt);
	if (inverse.nrow != result.basis.nrow || inverse.ncol != result.basis.ncol)
		throw std::runtime_error("Adapted dlogmat basis is not invertible.");
	result.adapted_cyclic = sparse_mat_mul(
		sparse_mat_mul(result.basis, result.source_cyclic, F, &(opt->pool)),
		inverse, F, &(opt->pool));
	result.adapted_reflection = sparse_mat_mul(
		sparse_mat_mul(result.basis, result.source_reflection, F, &(opt->pool)),
		inverse, F, &(opt->pool));
	result.adapted_parity = sparse_mat_mul(
		sparse_mat_mul(result.basis, result.source_parity, F, &(opt->pool)),
		inverse, F, &(opt->pool));
	auto adapted_rows = sparse_mat_mul(
		result.basis, result.source_rows, F, &(opt->pool));
	result.tensor = dlogmat_irrep_tensor_from_rows(adapted_rows, &(opt->pool));
	if (!validate_hhex_dlogmat_decomposition(result, F, opt))
		throw std::runtime_error("Adapted dlogmat decomposition failed exact validation.");
	return result;
}
