#pragma once

#include "fec_provider_d3.hpp"

#include <array>
#include <filesystem>
#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

struct nmhv_l0c_v2_line_stage_stats_t {
	size_t weight = 0;
	size_t fec_dimension = 0;
	size_t ambient_dimension = 0;
	size_t image_nnz = 0;
	size_t line_rank = 0;
	size_t coordinate_nnz = 0;
	size_t basis_nnz = 0;
};

template <typename index_t>
struct nmhv_l0c_v2_line_reduction_t {
	sparse_mat<ulong, index_t> coordinates;
	sparse_mat<ulong, index_t> basis;
	nmhv_l0c_v2_line_stage_stats_t stats;
};

template <typename index_t>
struct nmhv_l0c_v2_line_projection_t {
	sparse_mat<ulong, index_t> coordinates;
	std::vector<nmhv_l0c_v2_line_stage_stats_t> stages;
};

// Select the RREF pivot rows as an abstract line-symbol basis. Since their
// pivot columns form an identity block, the same columns of the unreduced
// image are exactly its coordinates in that basis.
template <typename index_t>
nmhv_l0c_v2_line_reduction_t<index_t> nmhv_l0c_v2_reduce_line_image(
	sparse_mat<ulong, index_t>&& image, size_t weight,
	const field_t& F, rref_option_t& opt) {
	image.compress(&(opt->pool));
	const size_t image_nnz = image.nnz();
	auto reduced = image;
	auto groups = sparse_mat_rref(reduced, F, opt);
	std::vector<pivot_t<index_t>> pivots;
	for (const auto& group : groups)
		pivots.insert(pivots.end(), group.begin(), group.end());
	std::ranges::sort(pivots, {}, &pivot_t<index_t>::c);

	const index_t missing = index_sval<index_t>();
	std::vector<index_t> pivot_owner(image.ncol, missing);
	sparse_mat<ulong, index_t> basis(pivots.size(), image.ncol);
	for (size_t row = 0; row < pivots.size(); ++row) {
		if (pivots[row].r < 0 || pivots[row].c < 0
			|| static_cast<size_t>(pivots[row].r) >= reduced.nrow
			|| static_cast<size_t>(pivots[row].c) >= reduced.ncol)
			throw std::runtime_error("Line-image RREF returned an invalid pivot.");
		pivot_owner[pivots[row].c] = static_cast<index_t>(row);
		basis[row] = std::move(reduced[pivots[row].r]);
		basis[row].canonicalize();
		basis[row].sort_indices();
	}
	for (size_t row = 0; row < basis.nrow; ++row) {
		const auto diagonal = basis[row].find(pivots[row].c, false);
		if (diagonal == nullptr || *diagonal != 1)
			throw std::runtime_error(
				"Line-image RREF pivot block is not the identity matrix.");
		for (size_t pos = 0; pos < basis[row].nnz(); ++pos) {
			const index_t owner = pivot_owner[basis[row](pos)];
			if (owner != missing && owner != static_cast<index_t>(row))
				throw std::runtime_error(
					"Line-image RREF did not complete backward substitution.");
		}
	}

	sparse_mat<ulong, index_t> coordinates(image.nrow, pivots.size());
	opt->pool.detach_loop(0, image.nrow, [&](size_t row) {
		for (size_t pos = 0; pos < image[row].nnz(); ++pos) {
			const index_t column = pivot_owner[image[row](pos)];
			if (column != missing)
				coordinates[row].push_back(column, image[row][pos]);
		}
	});
	opt->pool.wait();

	const nmhv_l0c_v2_line_stage_stats_t stats{
		weight, image.nrow, image.ncol, image_nnz, pivots.size(),
		coordinates.nnz(), basis.nnz()};
	return {std::move(coordinates), std::move(basis), stats};
}

template <typename index_t>
void nmhv_l0c_v2_accumulate_collinear_slices(
	std::array<sparse_mat<ulong, index_t>, 4>& slices,
	const sparse_tensor<ulong, index_t, SPARSE_CSR>& tensor,
	size_t row_offset, const field_t& F) {
	for (size_t row = 0; row < tensor.dim(0); ++row) {
		for (size_t pos = tensor.rowptr()[row]; pos < tensor.rowptr()[row + 1]; ++pos) {
			const auto index = tensor.index(pos);
			for (size_t target = 0; target < slices.size(); ++target) {
				const int coefficient =
					nmhv_limit_collinear_map()[index[1]][target];
				if (coefficient != 0)
					slices[target][row_offset + row].push_back(index[0],
						scalar_mul(tensor.val(pos),
							nmhv_limit_scalar<ulong>(coefficient, 1, F), F));
			}
		}
	}
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_l0c_v2_build_line_image(
	std::array<sparse_mat<ulong, index_t>, 4>& slices,
	const sparse_mat<ulong, index_t>& previous,
	const field_t& F, thread_pool* pool) {
	const size_t dimension = slices.front().nrow;
	sparse_mat<ulong, index_t> image(dimension, previous.ncol * slices.size());
	for (size_t target = 0; target < slices.size(); ++target) {
		slices[target].compress(pool);
		auto mapped = sparse_mat_mul(slices[target], previous, F, pool);
		pool->detach_loop(0, mapped.nrow, [&](size_t row) {
			for (size_t pos = 0; pos < mapped[row].nnz(); ++pos)
				image[row].push_back(static_cast<index_t>(
					static_cast<size_t>(mapped[row](pos)) * slices.size() + target),
					mapped[row][pos]);
		});
		pool->wait();
	}
	return image;
}

template <typename index_t>
nmhv_l0c_v2_line_projection_t<index_t> nmhv_l0c_v2_d3_exact_collinear(
	const std::filesystem::path& basis_root, size_t total_weight,
	const field_t& F, rref_option_t& opt) {
	sparse_mat<ulong, index_t> previous(1, 1);
	previous[0].push_back(0, 1);
	std::vector<nmhv_l0c_v2_line_stage_stats_t> stages;
	stages.reserve(total_weight);
	for (size_t weight = 1; weight <= total_weight; ++weight) {
		const auto directory = nmhv_limit_weight_directory(basis_root, weight);
		const auto blocks = nmhv_limit_read_irrep_axis(directory /
			("FEC_" + std::to_string(weight) + "_blocks.tsv"));
		const size_t dimension = fec_irrep_axis_dimension(blocks);
		std::array<sparse_mat<ulong, index_t>, 4> slices{
			sparse_mat<ulong, index_t>(dimension, previous.nrow),
			sparse_mat<ulong, index_t>(dimension, previous.nrow),
			sparse_mat<ulong, index_t>(dimension, previous.nrow),
			sparse_mat<ulong, index_t>(dimension, previous.nrow)};
		for (const auto& block : blocks) {
			if (block.dimension == 0)
				continue;
			const auto path = directory / ("FEC_" + std::to_string(weight)
				+ "_" + block.spec.name + ".wxf");
			auto tensor = sparse_tensor_read_wxf<ulong, index_t>(
				path, F, &(opt->pool));
			if (tensor.rank() != 3 || tensor.dim(0) != block.dimension
				|| tensor.dim(1) != previous.nrow || tensor.dim(2) != 9)
				throw std::runtime_error(
					"D3 FEC block is incompatible with the exact collinear chain.");
			nmhv_l0c_v2_accumulate_collinear_slices(
				slices, tensor, block.offset, F);
		}
		auto image = nmhv_l0c_v2_build_line_image(
			slices, previous, F, &(opt->pool));
		auto current = nmhv_l0c_v2_reduce_line_image(
			std::move(image), weight, F, opt);
		if (opt->verbose)
			std::cout << "Exact collinear weight " << weight << "/"
				<< total_weight << ": FEC " << current.stats.fec_dimension
				<< ", line rank " << current.stats.line_rank
				<< ", coordinates nnz " << current.stats.coordinate_nnz
				<< std::endl;
		previous = std::move(current.coordinates);
		stages.push_back(current.stats);
	}
	return {std::move(previous), std::move(stages)};
}
