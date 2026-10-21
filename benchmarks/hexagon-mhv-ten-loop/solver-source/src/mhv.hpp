#pragma once

#include "joint_eet_mrk.hpp"
#include "joint_eet_solver.hpp"
#include "nmhv/exact_collinear.hpp"

// All MHV unknowns live in the pe_A1 block. The exported axis is nevertheless
// the complete FEC axis, so ordinary FEC contraction needs no generator map.
template <typename index_t>
sparse_mat<ulong, index_t> mhv_bulk(const std::filesystem::path& root,
	size_t weight, const field_t& F, rref_option_t& opt) {
	const auto directory = nmhv_limit_weight_directory(root, weight);
	auto FEC = sparse_tensor_read_wxf<ulong, index_t>(directory /
		("FEC_" + std::to_string(weight) + "_pe_A1.wxf"), F, &(opt->pool));
	if (FEC.rank() != 3 || FEC.dim(2) != 9)
		throw std::runtime_error("MHV requires rank-3 Hhex FEC tensors.");
	sparse_mat<ulong, index_t> mat(3 * FEC.dim(1), FEC.dim(0));
	for (size_t row = 0; row < FEC.dim(0); ++row)
		for (size_t pos = FEC.rowptr()[row]; pos < FEC.rowptr()[row + 1]; ++pos) {
			const auto index = FEC.index(pos);
			if (index[1] < 3)
				mat[3 * index[0] + index[1]].push_back(row, FEC.val(pos));
		}
	nmhv_sew_pivot_stage_t<index_t> stage;
	return nmhv_sew_modular_kernel_rows(std::move(mat), stage, F, opt);
}

template <typename index_t>
std::vector<nmhv_physical_equation_t> mhv_limit_equations(
	const std::filesystem::path& root, size_t loop, bool mrk,
	size_t samples, uint64_t seed, size_t offset, size_t dimension,
	const nmhv_physical_symbol_table_t& symbols, const field_t& F, thread_pool* pool) {
	const size_t degree = mrk ? 2 * loop - 1 : 0;
	auto values = nmhv_limit_evaluate_fec_chain<ulong, index_t>(root, 2 * loop,
		mrk ? nmhv_limit_kind_t::mrk : nmhv_limit_kind_t::collinear,
		samples, degree, seed, F, pool);
	std::vector<nmhv_physical_equation_t> equations(samples * (degree + 1));
	pool->detach_loop(0, samples, [&](size_t sample) {
		// The shared NMHV RHS is for the sum of two E components. MHV has
		// half that normalization; its MRK input has no helicity-flip factor.
		const auto rhs = mrk ? nmhv_physical_mrk_rhs(loop,
			nmhv_physical_sample_mrk_slots(seed, sample, 2 * loop - 1, F),
			symbols, degree, F) : std::vector<ulong>{nmhv_physical_collinear_rhs(loop,
			nmhv_physical_sample_collinear_slots(seed, sample, 2 * loop, F), F)};
		for (size_t k = 0; k <= degree; ++k) {
			auto& equation = equations[sample * (degree + 1) + k];
			equation.condition = mrk ? "MHV MRK" : "MHV collinear";
			equation.sample = sample;
			equation.rhs = scalar_div(rhs[k], ulong(2), F);
			equation.coefficients.resize(dimension);
			for (size_t column = 0; column < dimension; ++column)
				equation.coefficients[column] = values(sample, offset + column, k);
		}
	});
	pool->wait();
	return equations;
}

template <typename index_t>
void mhv_check_affine(const std::vector<nmhv_physical_equation_t>& equations,
	const nmhv_physical_affine_solution_t<index_t>& space,
	const field_t& F, thread_pool* pool) {
	const auto residual = nmhv_physical_restrict_equations(
		equations, space.origin, space.directions, F, pool);
	for (const auto& row : residual) {
		if (row.rhs != 0 || std::any_of(row.coefficients.begin(), row.coefficients.end(),
			[](ulong value) { return value != 0; }))
			throw std::runtime_error("Independent MHV limit validation failed.");
	}
}

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> mhv_collinear(
	const std::filesystem::path& root, size_t loop,
	const sparse_mat<ulong, index_t>& bulk, size_t offset,
	size_t samples, uint64_t seed, const field_t& F, rref_option_t& opt) {
	// Determine the homogeneous kernel by exact line projection, without
	// identifying sampled zeros with exact line-symbol relations.
	auto line = nmhv_l0c_v2_d3_exact_collinear<index_t>(root, 2 * loop, F, opt);
	sparse_mat<ulong, index_t> local(bulk.ncol, line.coordinates.ncol);
	for (size_t row = 0; row < bulk.ncol; ++row)
		local[row] = std::move(line.coordinates[offset + row]);
	auto image = sparse_mat_mul(bulk, local, F, &(opt->pool));
	nmhv_sew_pivot_stage_t<index_t> stage;
	auto kernel = nmhv_sew_modular_kernel_rows(image.transpose(), stage, F, opt);
	// SparseRREF represents a full-rank kernel as 0 x 0.
	if (kernel.nrow == 0) kernel = sparse_mat<ulong, index_t>(0, bulk.nrow);
	nmhv_physical_affine_solution_t<index_t> parent;
	parent.origin.resize(bulk.ncol);
	parent.directions = bulk;
	const auto equations = mhv_limit_equations<index_t>(root, loop, false,
		samples, seed, offset, bulk.ncol, {}, F, &(opt->pool));
	auto restricted = nmhv_physical_restrict_equations(
		equations, parent.origin, bulk, F, &(opt->pool));
	auto solution = nmhv_physical_solve_affine_space<index_t>(restricted, bulk.nrow, F);
	if (solution.directions.nrow != kernel.nrow)
		throw std::runtime_error("Collinear samples did not reach the exact rank; increase --samples.");
	solution.directions = std::move(kernel);
	auto result = nmhv_physical_compose_affine(parent, solution, F, &(opt->pool));
	if (result.directions.nrow == 0) result.directions = sparse_mat<ulong, index_t>(0, bulk.ncol);
	mhv_check_affine(equations, result, F, &(opt->pool));
	mhv_check_affine(mhv_limit_equations<index_t>(root, loop, false, 16,
		seed ^ 0x713ab19ULL, offset, bulk.ncol, {}, F, &(opt->pool)), result, F, &(opt->pool));
	return result;
}

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> mhv_mrk(
	const std::filesystem::path& root, size_t loop,
	const nmhv_physical_affine_solution_t<index_t>& collinear, size_t offset,
	size_t samples, uint64_t seed, const nmhv_physical_symbol_table_t& symbols,
	const field_t& F, rref_option_t& opt) {
	const auto equations = mhv_limit_equations<index_t>(root, loop, true,
		samples, seed, offset, collinear.origin.size(), symbols, F, &(opt->pool));
	const auto restricted = nmhv_physical_restrict_equations(equations,
		collinear.origin, collinear.directions, F, &(opt->pool));
	const auto solution = nmhv_physical_solve_affine_space<index_t>(
		restricted, collinear.directions.nrow, F);
	auto result = nmhv_physical_compose_affine(collinear, solution, F, &(opt->pool));
	if (result.directions.nrow == 0)
		result.directions = sparse_mat<ulong, index_t>(0, collinear.origin.size());
	mhv_check_affine(mhv_limit_equations<index_t>(root, loop, true, 16,
		seed ^ 0x592be17ULL, offset, collinear.origin.size(), symbols, F, &(opt->pool)),
		result, F, &(opt->pool));
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> mhv_pack(
	const nmhv_physical_affine_solution_t<index_t>& space, size_t offset, size_t dimension) {
	sparse_mat<ulong, index_t> result(space.directions.nrow + 1, dimension);
	for (size_t column = 0; column < space.origin.size(); ++column)
		if (space.origin[column] != 0)
			result[0].push_back(offset + column, space.origin[column]);
	for (size_t row = 0; row < space.directions.nrow; ++row)
		for (size_t pos = 0; pos < space.directions[row].nnz(); ++pos)
			result[row + 1].push_back(offset + space.directions[row](pos), space.directions[row][pos]);
	return result;
}

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> mhv_unpack(
	const sparse_mat<ulong, index_t>& packed, size_t offset, size_t dimension) {
	if (packed.nrow == 0 || packed.ncol < offset + dimension)
		throw std::runtime_error("Invalid MHV affine cache dimensions.");
	nmhv_physical_affine_solution_t<index_t> result;
	result.origin.resize(dimension);
	result.directions = sparse_mat<ulong, index_t>(packed.nrow - 1, dimension);
	for (size_t row = 0; row < packed.nrow; ++row)
		for (size_t pos = 0; pos < packed[row].nnz(); ++pos) {
			const size_t column = packed[row](pos);
			if (column < offset || column >= offset + dimension)
				throw std::runtime_error("MHV cache has coefficients outside pe_A1.");
			if (row == 0) result.origin[column - offset] = packed[row][pos];
			else result.directions[row - 1].push_back(column - offset, packed[row][pos]);
		}
	return result;
}
