#pragma once

#include "../nmhv_limit.hpp"

#include <filesystem>
#include <string>
#include <vector>

template <typename scalar_t, typename index_t>
struct nmhv_l0c_generator_maps_t {
	sparse_mat<scalar_t, index_t> E;
	sparse_mat<scalar_t, index_t> Et_before;
	sparse_mat<scalar_t, index_t> Et;
};

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_l0c_coordinate_basis(
	size_t dimension, const std::vector<size_t>& coordinates,
	const field_t& F) {
	sparse_mat<scalar_t, index_t> result(coordinates.size(), dimension);
	for (size_t row = 0; row < coordinates.size(); ++row)
		result[row].push_back(static_cast<index_t>(coordinates[row]),
			nmhv_limit_scalar<scalar_t>(1, 1, F));
	return result;
}

template <typename scalar_t, typename index_t>
void nmhv_l0c_add_scaled(sparse_mat<scalar_t, index_t>& target,
	const sparse_mat<scalar_t, index_t>& source, const scalar_t& scale,
	const field_t& F, thread_pool* pool) {
	if (target.nrow != source.nrow || target.ncol != source.ncol)
		throw std::runtime_error("Cannot add matrices with different dimensions.");
	pool->detach_loop(0, target.nrow, [&](size_t row) {
		sparse_vec_add_mul(target[row], source[row], scale, F);
	});
	pool->wait();
}

template <typename scalar_t, typename index_t>
void nmhv_l0c_shift_columns(
	sparse_mat<scalar_t, index_t>& matrix, size_t offset) {
	for (size_t row = 0; row < matrix.nrow; ++row)
		for (size_t pos = 0; pos < matrix[row].nnz(); ++pos)
			matrix[row](pos) += static_cast<index_t>(offset);
	matrix.ncol += offset;
}

template <typename scalar_t, typename index_t>
void nmhv_l0c_accumulate_variable(
	sparse_mat<scalar_t, index_t>& constraints,
	const sparse_mat<scalar_t, index_t>& variable_map,
	const sparse_vec<scalar_t, index_t>& relation_coefficients,
	size_t prefix_count, const field_t& F, thread_pool* pool) {
	pool->detach_loop(0, prefix_count, [&](size_t prefix) {
		for (size_t pos = 0; pos < relation_coefficients.nnz(); ++pos) {
			const size_t relation = relation_coefficients(pos);
			sparse_vec_add_mul(
				constraints[relation * prefix_count + prefix],
				variable_map[prefix], relation_coefficients[pos], F);
		}
	});
	pool->wait();
}

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_l0c_final_entry_constraints(
	const std::vector<sparse_mat<scalar_t, index_t>>& letters,
	const sparse_mat<scalar_t, index_t>& cyclic,
	const sparse_mat<scalar_t, index_t>& cyclic2,
	const sparse_mat<scalar_t, index_t>& E,
	const sparse_mat<scalar_t, index_t>& Et,
	const sparse_mat<scalar_t, index_t>& relations,
	const field_t& F, thread_pool* pool) {
	if (letters.size() != 9 || relations.nrow != 27 || relations.ncol != 54)
		throw std::runtime_error(
			"NMHV last-entry inputs require nine letters and a 27x54 relation matrix.");
	const size_t prefix_count = letters.front().ncol;
	sparse_mat<scalar_t, index_t> constraints(
		relations.nrow * prefix_count, E.nrow + Et.nrow);
	auto relations_copy = relations;
	auto relation_by_variable = relations_copy.transpose();
	for (size_t component = 0; component < 3; ++component) {
		sparse_mat<scalar_t, index_t> E_component;
		sparse_mat<scalar_t, index_t> Et_component;
		if (component == 0) {
			E_component = E;
			Et_component = Et;
		}
		else if (component == 1) {
			E_component = sparse_mat_mul(E, cyclic, F, pool);
			Et_component = sparse_mat_mul(Et, cyclic, F, pool);
		}
		else {
			E_component = sparse_mat_mul(E, cyclic2, F, pool);
			Et_component = sparse_mat_mul(Et, cyclic2, F, pool);
		}
		for (size_t letter = 0; letter < letters.size(); ++letter) {
			auto map = sparse_mat_mul(
				E_component, letters[letter], F, pool).transpose();
			nmhv_l0c_accumulate_variable(
				constraints, map,
				relation_by_variable[component * 9 + letter],
				prefix_count, F, pool);
		}
		for (size_t letter = 0; letter < letters.size(); ++letter) {
			auto map = sparse_mat_mul(
				Et_component, letters[letter], F, pool).transpose();
			nmhv_l0c_shift_columns(map, E.nrow);
			nmhv_l0c_accumulate_variable(
				constraints, map,
				relation_by_variable[27 + component * 9 + letter],
				prefix_count, F, pool);
		}
	}
	constraints.compress(pool);
	return constraints;
}

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_l0c_collinear_generator_map(
	const sparse_mat<scalar_t, index_t>& E,
	const sparse_mat<scalar_t, index_t>& Et,
	const sparse_mat<scalar_t, index_t>& cyclic2,
	const field_t& F, thread_pool* pool) {
	auto cyclic_E = sparse_mat_mul(E, cyclic2, F, pool);
	auto result = E;
	nmhv_l0c_add_scaled(result, cyclic_E,
		nmhv_limit_scalar<scalar_t>(1, 1, F), F, pool);
	sparse_mat<scalar_t, index_t> zero_Et(Et.nrow, Et.ncol);
	result.append(std::move(zero_Et));
	result.compress(pool);
	return result;
}

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_l0c_project_constraints(
	const sparse_mat<scalar_t, index_t>& generator_map,
	const nmhv_limit_chain_values_t<scalar_t>& chain,
	const field_t& F, thread_pool* pool) {
	if (generator_map.ncol != chain.basis_count || chain.stride != 1)
		throw std::runtime_error("Collinear projection dimensions are inconsistent.");
	sparse_mat<scalar_t, index_t> result(
		chain.sample_count, generator_map.nrow);
	pool->detach_loop(0, chain.sample_count, [&](size_t sample) {
		for (size_t generator = 0; generator < generator_map.nrow; ++generator) {
			scalar_t value = 0;
			for (size_t pos = 0; pos < generator_map[generator].nnz(); ++pos) {
				value = scalar_add(value, scalar_mul(
					generator_map[generator][pos],
					chain(sample, generator_map[generator](pos), 0), F), F);
			}
			if (value != 0)
				result[sample].push_back(static_cast<index_t>(generator), value);
		}
	});
	pool->wait();
	return result;
}
