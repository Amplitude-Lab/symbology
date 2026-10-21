#pragma once

#include "joint_eet_lec.hpp"
#include "nmhv/exact_collinear.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

template <typename index_t>
struct joint_eet_collinear_prime_result_t {
	sparse_mat<ulong, index_t> basis;
	size_t source_nnz = 0;
	size_t image_ncol = 0;
	size_t image_nnz = 0;
	size_t collinear_rank = 0;
	size_t raw_kernel_nnz = 0;
	size_t kernel_nnz = 0;
	size_t restricted_nnz = 0;
};

template <typename T, typename index_t>
sparse_mat<T, index_t> joint_eet_tail_collinear_projection(
	const sparse_tensor<T, index_t, SPARSE_CSR>& tail,
	const field_t& F, thread_pool* pool) {
	if (tail.rank() < 3 || tail.dim(tail.rank() - 1) != 5)
		throw std::runtime_error(
			"Joint E/Et tail must have layout {tail,letters...,component}.");
	const size_t weight = tail.rank() - 2;
	for (size_t axis = 1; axis <= weight; ++axis)
		if (tail.dim(axis) != 9)
			throw std::runtime_error("Joint E/Et tail uses a non-Hhex alphabet.");
	const size_t target_count = nmhv_sew_checked_power(
		4, weight, "Collinear tail word count");
	// Strict collinear projects {E,CE,C^2E,Et,CEt} to E+C^2E;
	// the parity-odd Et components vanish on the line.
	static constexpr std::array<int, 5> component_projection{{1, 0, 1, 0, 0}};
	sparse_mat<T, index_t> result(tail.dim(0), target_count);
	pool->detach_loop(0, tail.dim(0), [&](size_t row) {
		std::vector<T> values(target_count);
		for (size_t pos = tail.rowptr()[row]; pos < tail.rowptr()[row + 1]; ++pos) {
			const auto index = tail.index(pos);
			const int component_coefficient = component_projection[index[weight]];
			if (component_coefficient == 0)
				continue;
			std::vector<std::pair<size_t, int>> terms{{0, component_coefficient}};
			for (size_t depth = 0; depth < weight; ++depth) {
				std::vector<std::pair<size_t, int>> next;
				for (const auto& [target_word, coefficient] : terms) {
					for (size_t target = 0; target < 4; ++target) {
						const int letter_coefficient =
							nmhv_limit_collinear_map()[index[depth]][target];
						if (letter_coefficient != 0)
							next.emplace_back(target_word * 4 + target,
								coefficient * letter_coefficient);
					}
				}
				terms = std::move(next);
			}
			for (const auto& [target_word, coefficient] : terms)
				values[target_word] = scalar_add(values[target_word],
					scalar_mul(tail.val(pos),
						nmhv_limit_scalar<T>(coefficient, 1, F), F), F);
		}
		for (size_t target = 0; target < values.size(); ++target)
			if (values[target] != 0)
				result[row].push_back(static_cast<index_t>(target), values[target]);
	});
	pool->wait();
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> joint_eet_project_natural(
	const sparse_mat<ulong, index_t>& basis,
	const sparse_mat<ulong, index_t>& head_projection,
	const sparse_mat<ulong, index_t>& tail_projection,
	const field_t& F, thread_pool* pool) {
	if (basis.ncol != head_projection.nrow * tail_projection.nrow)
		throw std::runtime_error(
			"Natural sewing and collinear projection dimensions are incompatible: "
			+ std::to_string(basis.nrow) + "x" + std::to_string(basis.ncol)
			+ ", head " + std::to_string(head_projection.nrow) + "x"
			+ std::to_string(head_projection.ncol) + ", tail "
			+ std::to_string(tail_projection.nrow) + "x"
			+ std::to_string(tail_projection.ncol) + ".");
	std::vector<sparse_mat<ulong, index_t>> head_slices;
	for (size_t target = 0; target < tail_projection.ncol; ++target)
		head_slices.emplace_back(basis.nrow, head_projection.nrow);

	pool->detach_loop(0, basis.nrow, [&](size_t row) {
		std::vector<std::vector<ulong>> values(
			tail_projection.ncol, std::vector<ulong>(head_projection.nrow));
		for (size_t pos = 0; pos < basis[row].nnz(); ++pos) {
			const size_t column = static_cast<size_t>(basis[row](pos));
			const size_t head = column / tail_projection.nrow;
			const size_t tail = column % tail_projection.nrow;
			for (size_t tail_pos = 0;
				tail_pos < tail_projection[tail].nnz(); ++tail_pos) {
				const size_t target = tail_projection[tail](tail_pos);
				values[target][head] = scalar_add(values[target][head],
					scalar_mul(basis[row][pos],
						tail_projection[tail][tail_pos], F), F);
			}
		}
		for (size_t target = 0; target < values.size(); ++target)
			for (size_t head = 0; head < values[target].size(); ++head)
				if (values[target][head] != 0)
					head_slices[target][row].push_back(
						static_cast<index_t>(head), values[target][head]);
	});
	pool->wait();

	sparse_mat<ulong, index_t> result(
		basis.nrow, head_projection.ncol * tail_projection.ncol);
	for (size_t target = 0; target < head_slices.size(); ++target) {
		auto projected = sparse_mat_mul(
			head_slices[target], head_projection, F, pool);
		pool->detach_loop(0, projected.nrow, [&](size_t row) {
			for (size_t pos = 0; pos < projected[row].nnz(); ++pos)
				result[row].push_back(static_cast<index_t>(
					static_cast<size_t>(projected[row](pos))
						* tail_projection.ncol + target), projected[row][pos]);
		});
		pool->wait();
	}
	result.compress(pool);
	return result;
}

template <typename index_t>
joint_eet_collinear_prime_result_t<index_t> joint_eet_restrict_basis(
	sparse_mat<ulong, index_t> basis,
	const sparse_mat<ulong, index_t>& head_projection,
	const sparse_mat<ulong, index_t>& tail_projection,
	nmhv_sew_pivot_stage_t<index_t>& collinear_stage,
	nmhv_sew_pivot_stage_t<index_t>& final_stage,
	nmhv_sew_pivot_stage_t<index_t>& second_kernel_stage,
	int presentation_mode,
	const field_t& F, rref_option_t& opt) {
	joint_eet_collinear_prime_result_t<index_t> result;
	result.source_nnz = basis.nnz();
	auto image = joint_eet_project_natural(
		basis, head_projection, tail_projection, F, &(opt->pool));
	result.image_ncol = image.ncol;
	result.image_nnz = image.nnz();
	auto kernel = nmhv_sew_modular_kernel_rows(
		image.transpose(), collinear_stage, F, opt);
	result.collinear_rank = nmhv_sew_modular_rank(collinear_stage.pivots);
	result.raw_kernel_nnz = kernel.nnz();
	if (presentation_mode == 2) {
		auto complement = nmhv_sew_modular_kernel_rows(
			std::move(kernel), final_stage, F, opt);
		kernel = nmhv_sew_modular_kernel_rows(
			std::move(complement), second_kernel_stage, F, opt);
	}
	else if (presentation_mode < 0 || presentation_mode > 2)
		throw std::runtime_error("Unknown collinear presentation mode.");
	result.kernel_nnz = kernel.nnz();
	basis = sparse_mat_mul(kernel, basis, F, &(opt->pool));
	nmhv_sew_modular_canonicalize(basis, &(opt->pool));
	result.restricted_nnz = basis.nnz();
	if (presentation_mode == 1) {
		nmhv_sew_modular_rref(basis, final_stage, F, opt);
		result.basis = nmhv_sew_modular_pivot_rows(
			std::move(basis), final_stage.pivots);
		nmhv_sew_modular_canonicalize(result.basis, &(opt->pool));
	}
	else
		result.basis = std::move(basis);
	return result;
}
