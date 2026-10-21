#pragma once

#include "fec_provider.hpp"

#include <string_view>

template <typename scalar_t, typename index_t>
struct nmhv_l0c_d3_final_entry_data_t {
	std::vector<fec_irrep_axis_block_t> blocks;
	nmhv_l0c_generator_maps_t<scalar_t, index_t> generators;
	sparse_mat<scalar_t, index_t> cyclic;
	sparse_mat<scalar_t, index_t> cyclic2;
	std::vector<sparse_mat<scalar_t, index_t>> letters;
	size_t prefix_count = 0;
};

inline const fec_irrep_axis_block_t& nmhv_l0c_d3_find_block(
	const std::vector<fec_irrep_axis_block_t>& blocks,
	std::string_view name) {
	for (const auto& block : blocks)
		if (block.spec.name == name)
			return block;
	throw std::runtime_error(
		"Missing parity-dihedral FEC block: " + std::string(name) + ".");
}

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_l0c_d3_raw_cyclic_action(
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, thread_pool* pool) {
	return nmhv_limit_raw_cyclic_action<scalar_t, index_t>(blocks, F, pool);
}

template <typename scalar_t, typename index_t>
nmhv_l0c_generator_maps_t<scalar_t, index_t>
nmhv_l0c_d3_generator_maps(
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F) {
	const size_t dimension = fec_irrep_axis_dimension(blocks);
	const auto& pe_A1 = nmhv_l0c_d3_find_block(blocks, "pe_A1");
	const auto& pe_E = nmhv_l0c_d3_find_block(blocks, "pe_E");
	const auto& po_A2 = nmhv_l0c_d3_find_block(blocks, "po_A2");
	const auto& po_E = nmhv_l0c_d3_find_block(blocks, "po_E");
	std::vector<size_t> E_coordinates;
	std::vector<size_t> Et_before_coordinates;
	std::vector<size_t> Et_coordinates;
	for (size_t i = 0; i < pe_A1.dimension; ++i)
		E_coordinates.push_back(pe_A1.offset + i);
	for (size_t i = 0; i < pe_E.multiplicity; ++i)
		E_coordinates.push_back(pe_E.offset + i);
	for (size_t i = 0; i < po_A2.dimension; ++i)
		Et_before_coordinates.push_back(po_A2.offset + i);
	for (size_t i = 0; i < po_E.multiplicity; ++i) {
		Et_before_coordinates.push_back(po_E.offset + po_E.multiplicity + i);
		Et_coordinates.push_back(po_E.offset + po_E.multiplicity + i);
	}
	return {
		nmhv_l0c_coordinate_basis<scalar_t, index_t>(
			dimension, E_coordinates, F),
		nmhv_l0c_coordinate_basis<scalar_t, index_t>(
			dimension, Et_before_coordinates, F),
		nmhv_l0c_coordinate_basis<scalar_t, index_t>(
			dimension, Et_coordinates, F)};
}

template <typename scalar_t, typename index_t>
nmhv_l0c_d3_final_entry_data_t<scalar_t, index_t>
nmhv_l0c_d3_final_entry_data(
	const std::filesystem::path& basis_root, size_t weight,
	const field_t& F, thread_pool* pool) {
	const auto directory = nmhv_limit_weight_directory(basis_root, weight);
	auto blocks = nmhv_limit_read_irrep_axis(
		directory / ("FEC_" + std::to_string(weight) + "_blocks.tsv"));
	const size_t dimension = fec_irrep_axis_dimension(blocks);
	auto generators = nmhv_l0c_d3_generator_maps<scalar_t, index_t>(blocks, F);
	auto cyclic = nmhv_l0c_d3_raw_cyclic_action<scalar_t, index_t>(
		blocks, F, pool);
	auto cyclic2 = sparse_mat_mul(cyclic, cyclic, F, pool);
	std::vector<sparse_mat<scalar_t, index_t>> letters;
	size_t prefix_count = 0;
	for (const auto& name : {"pe_A1", "pe_E", "po_E"}) {
		const auto& block = nmhv_l0c_d3_find_block(blocks, name);
		if (block.dimension == 0)
			continue;
		const auto path = directory / ("FEC_" + std::to_string(weight)
			+ "_" + std::string(name) + ".wxf");
		if (!std::filesystem::is_regular_file(path))
			throw std::runtime_error("Missing D3 FEC block: " + path.string());
		auto tensor = sparse_tensor_read_wxf<scalar_t, index_t>(path, F, pool);
		if (tensor.rank() != 3 || tensor.dim(0) != block.dimension
			|| tensor.dim(2) != 9)
			throw std::runtime_error("D3 FEC block has incompatible dimensions.");
		if (letters.empty()) {
			prefix_count = tensor.dim(1);
			for (size_t letter = 0; letter < 9; ++letter)
				letters.emplace_back(dimension, prefix_count);
		}
		else if (tensor.dim(1) != prefix_count) {
			throw std::runtime_error("D3 FEC blocks use different predecessor axes.");
		}
		for (size_t row = 0; row < tensor.dim(0); ++row) {
			for (size_t pos = tensor.rowptr()[row]; pos < tensor.rowptr()[row + 1]; ++pos) {
				const auto index = tensor.index(pos);
				letters[index[1]][block.offset + row].push_back(index[0], tensor.val(pos));
			}
		}
	}
	if (letters.empty())
		throw std::runtime_error("D3 final-entry provider found no generator blocks.");
	for (auto& letter : letters)
		letter.compress(pool);
	return {std::move(blocks), std::move(generators), std::move(cyclic),
		std::move(cyclic2), std::move(letters), prefix_count};
}

template <typename scalar_t, typename index_t>
nmhv_limit_chain_values_t<scalar_t> nmhv_l0c_d3_evaluate_collinear(
	const std::filesystem::path& basis_root, size_t weight,
	size_t samples, uint64_t seed, const field_t& F, thread_pool* pool) {
	return nmhv_limit_evaluate_fec_chain<scalar_t, index_t>(
		basis_root, weight, nmhv_limit_kind_t::collinear,
		samples, 0, seed, F, pool);
}
