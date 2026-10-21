#pragma once

#include "fec_irrep_reduced.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

inline const std::array<std::array<int, 4>, 9>&
nmhv_limit_collinear_map() {
	// Target letters are {T,S,1+S^2,F} in the standard w -> 0 limit.
	static constexpr std::array<std::array<int, 4>, 9> map{{
		{{-1, -1, 0, 0}}, {{-1, 1, 0, 0}}, {{1, -1, 1, 0}},
		{{0, 1, 0, 0}}, {{0, -1, 0, 0}}, {{-1, 0, 0, 0}},
		{{0, 0, 0, 0}}, {{0, 0, 0, 0}}, {{0, 0, 0, 1}}
	}};
	return map;
}

inline const std::array<std::array<int, 5>, 9>&
nmhv_limit_mrk_map() {
	// Target letters are {delta,w,1+w,wbar,1+wbar} in u -> 1 MRK.
	static constexpr std::array<std::array<int, 5>, 9> map{{
		{{-2, -1, 2, -1, 2}}, {{0, -1, 0, -1, 0}},
		{{0, 1, 0, 1, 0}}, {{1, 0, 0, 0, 0}},
		{{-1, 0, 1, 0, 1}}, {{-1, -1, 1, -1, 1}},
		{{0, 0, 0, 0, 0}}, {{0, 0, -1, 0, 1}},
		{{0, -1, 1, 1, -1}}
	}};
	return map;
}

template <typename scalar_t>
scalar_t nmhv_limit_scalar(int numerator, int denominator, const field_t& F) {
	if constexpr (std::is_same_v<scalar_t, rat_t>)
		return rat_t(numerator, denominator);
	else {
		const auto absolute = static_cast<scalar_t>(
			numerator < 0 ? -static_cast<int64_t>(numerator) : numerator) % F.mod.n;
		const auto signed_numerator = numerator < 0 && absolute != 0
			? F.mod.n - absolute : absolute;
		return scalar_div(signed_numerator,
			static_cast<scalar_t>(denominator) % F.mod.n, F);
	}
}

template <typename scalar_t, typename index_t>
sparse_mat<scalar_t, index_t> nmhv_limit_raw_cyclic_action(
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, thread_pool* pool) {
	const size_t dimension = fec_irrep_axis_dimension(blocks);
	sparse_mat<scalar_t, index_t> cyclic(dimension, dimension);
	pool->detach_loop(0, blocks.size(), [&](size_t block_index) {
		const auto& block = blocks[block_index];
		for (size_t local = 0; local < block.dimension; ++local) {
			auto& row = cyclic[block.offset + local];
			if (block.spec.irrep != fec_irrep_t::E) {
				row.push_back(static_cast<index_t>(block.offset + local),
					nmhv_limit_scalar<scalar_t>(block.spec.parity, 1, F));
			}
			else if (local < block.multiplicity) {
				row.push_back(static_cast<index_t>(block.offset + local),
					nmhv_limit_scalar<scalar_t>(-block.spec.parity, 2, F));
				row.push_back(static_cast<index_t>(
					block.offset + block.multiplicity + local),
					nmhv_limit_scalar<scalar_t>(3 * block.spec.parity, 2, F));
			}
			else {
				row.push_back(static_cast<index_t>(
					block.offset + local - block.multiplicity),
					nmhv_limit_scalar<scalar_t>(-block.spec.parity, 2, F));
				row.push_back(static_cast<index_t>(block.offset + local),
					nmhv_limit_scalar<scalar_t>(-block.spec.parity, 2, F));
			}
		}
	});
	pool->wait();
	return cyclic;
}

template <typename scalar_t>
size_t nmhv_limit_dense_rank(
	std::vector<std::vector<scalar_t>> matrix, const field_t& F) {
	if (matrix.empty())
		return 0;
	const size_t nrow = matrix.size();
	const size_t ncol = matrix.front().size();
	size_t rank = 0;
	for (size_t column = 0; column < ncol && rank < nrow; ++column) {
		size_t pivot = rank;
		while (pivot < nrow && matrix[pivot][column] == 0)
			++pivot;
		if (pivot == nrow)
			continue;
		std::swap(matrix[rank], matrix[pivot]);
		const auto inverse = scalar_inv(matrix[rank][column], F);
		for (size_t j = column; j < ncol; ++j)
			matrix[rank][j] = scalar_mul(matrix[rank][j], inverse, F);
		for (size_t row = rank + 1; row < nrow; ++row) {
			const auto factor = matrix[row][column];
			if (factor == 0)
				continue;
			for (size_t j = column; j < ncol; ++j)
				matrix[row][j] = scalar_sub(matrix[row][j],
					scalar_mul(factor, matrix[rank][j], F), F);
		}
		++rank;
	}
	return rank;
}

inline std::vector<std::string> nmhv_limit_split_tabs(const std::string& line) {
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, '\t'))
		fields.push_back(field);
	return fields;
}

inline std::vector<fec_irrep_axis_block_t> nmhv_limit_read_irrep_axis(
	const std::filesystem::path& path) {
	std::ifstream input(path);
	if (!input)
		throw std::runtime_error("Cannot read FEC block metadata: " + path.string());
	std::string line;
	std::getline(input, line);
	if (line != "block\tparity\tirrep\tirrep_dimension\tmultiplicity\toffset\tdimension")
		throw std::runtime_error("Unexpected FEC block metadata header.");

	std::vector<fec_irrep_axis_block_t> blocks;
	const auto& specs = hhex_block_specs(fec_symmetry_t::parity_dihedral);
	while (std::getline(input, line)) {
		if (line.empty())
			continue;
		auto fields = nmhv_limit_split_tabs(line);
		if (fields.size() != 7 || blocks.size() >= specs.size())
			throw std::runtime_error("Malformed FEC block metadata row.");
		const auto& spec = specs[blocks.size()];
		const std::string irrep = spec.irrep == fec_irrep_t::A1 ? "A1"
			: spec.irrep == fec_irrep_t::A2 ? "A2" : "E";
		if (fields[0] != spec.name || std::stoi(fields[1]) != spec.parity
			|| fields[2] != irrep
			|| std::stoull(fields[3]) != spec.irrep_dimension)
			throw std::runtime_error("FEC block metadata uses a noncanonical order.");
		blocks.push_back({spec, std::stoull(fields[4]),
			std::stoull(fields[5]), std::stoull(fields[6])});
	}
	if (blocks.size() != specs.size())
		throw std::runtime_error("FEC block metadata must contain six blocks.");
	(void)fec_irrep_axis_dimension(blocks);
	return blocks;
}

inline std::filesystem::path nmhv_limit_weight_directory(
	const std::filesystem::path& basis_root, size_t weight) {
	std::ostringstream name;
	name << 'w' << std::setw(2) << std::setfill('0') << weight;
	return basis_root / name.str();
}

inline uint64_t nmhv_limit_splitmix64(uint64_t& state) {
	uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
	value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
	value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
	return value ^ (value >> 31);
}

enum class nmhv_limit_kind_t {
	collinear,
	mrk
};

template <typename scalar_t>
struct nmhv_limit_slot_weights_t {
	std::array<scalar_t, 9> regular{};
	std::array<scalar_t, 9> logarithmic{};
};

template <typename scalar_t>
nmhv_limit_slot_weights_t<scalar_t> nmhv_limit_source_weights(
	nmhv_limit_kind_t kind, size_t weight, uint64_t& state,
	const field_t& F) {
	nmhv_limit_slot_weights_t<scalar_t> source;
	if (kind == nmhv_limit_kind_t::mrk && weight == 1) {
		for (size_t letter = 0; letter < 3; ++letter)
			source.regular[letter] = nmhv_limit_scalar<scalar_t>(
				letter == 0 ? 1 : -1, 1, F);
		return source;
	}

	if (kind == nmhv_limit_kind_t::collinear) {
		std::array<scalar_t, 4> target{};
		for (auto& value : target)
			value = static_cast<scalar_t>(
				nmhv_limit_splitmix64(state) % (F.mod.n - 1)) + 1;
		for (size_t letter = 0; letter < source.regular.size(); ++letter) {
			for (size_t target_letter = 0; target_letter < target.size(); ++target_letter) {
				const int coefficient = nmhv_limit_collinear_map()[letter][target_letter];
				if (coefficient != 0)
					source.regular[letter] = scalar_add(source.regular[letter],
						scalar_mul(nmhv_limit_scalar<scalar_t>(coefficient, 1, F),
							target[target_letter], F), F);
			}
		}
		return source;
	}

	std::array<scalar_t, 5> target{};
	for (auto& value : target)
		value = static_cast<scalar_t>(
			nmhv_limit_splitmix64(state) % (F.mod.n - 1)) + 1;
	// Single-valuedness permits matching only the holomorphic coproducts.
	target[3] = target[4] = 0;
	for (size_t letter = 0; letter < source.regular.size(); ++letter) {
		source.logarithmic[letter] = scalar_mul(
			nmhv_limit_scalar<scalar_t>(nmhv_limit_mrk_map()[letter][0], 1, F),
			target[0], F);
		for (size_t target_letter = 1; target_letter < target.size(); ++target_letter) {
			const int coefficient = nmhv_limit_mrk_map()[letter][target_letter];
			if (coefficient != 0)
				source.regular[letter] = scalar_add(source.regular[letter],
					scalar_mul(nmhv_limit_scalar<scalar_t>(coefficient, 1, F),
						target[target_letter], F), F);
		}
	}
	return source;
}

template <typename scalar_t>
struct nmhv_limit_chain_values_t {
	size_t sample_count = 0;
	size_t basis_count = 0;
	size_t stride = 0;
	std::vector<scalar_t> values;

	scalar_t& operator()(size_t sample, size_t basis, size_t degree) {
		return values[(sample * basis_count + basis) * stride + degree];
	}
	const scalar_t& operator()(size_t sample, size_t basis, size_t degree) const {
		return values[(sample * basis_count + basis) * stride + degree];
	}
};

template <typename scalar_t, typename index_t>
nmhv_limit_chain_values_t<scalar_t> nmhv_limit_evaluate_fec_chain(
	const std::filesystem::path& basis_root, size_t total_weight,
	nmhv_limit_kind_t kind, size_t sample_count, size_t max_log_degree,
	uint64_t seed, const field_t& F, thread_pool* pool) {
	const size_t stride = max_log_degree + 1;
	nmhv_limit_chain_values_t<scalar_t> previous{
		sample_count, 1, stride, std::vector<scalar_t>(sample_count * stride)};
	for (size_t sample = 0; sample < sample_count; ++sample)
		previous(sample, 0, 0) = nmhv_limit_scalar<scalar_t>(1, 1, F);
	std::vector<uint64_t> states(sample_count);
	for (size_t sample = 0; sample < sample_count; ++sample)
		states[sample] = seed ^ (0xd1b54a32d192ed03ULL * (sample + 1));

	for (size_t weight = 1; weight <= total_weight; ++weight) {
		const auto directory = nmhv_limit_weight_directory(basis_root, weight);
		const auto blocks = nmhv_limit_read_irrep_axis(directory
			/ ("FEC_" + std::to_string(weight) + "_blocks.tsv"));
		const size_t basis_count = fec_irrep_axis_dimension(blocks);
		nmhv_limit_chain_values_t<scalar_t> current{sample_count, basis_count,
			stride, std::vector<scalar_t>(sample_count * basis_count * stride)};
		std::vector<nmhv_limit_slot_weights_t<scalar_t>> weights(sample_count);
		for (size_t sample = 0; sample < sample_count; ++sample)
			weights[sample] = nmhv_limit_source_weights<scalar_t>(
				kind, weight, states[sample], F);

		for (const auto& block : blocks) {
			if (block.dimension == 0)
				continue;
			const auto path = directory / ("FEC_" + std::to_string(weight)
				+ "_" + block.spec.name + ".wxf");
			auto tensor = sparse_tensor_read_wxf<scalar_t, index_t>(path, F, pool);
			if (tensor.rank() != 3 || tensor.dim(0) != block.dimension
				|| tensor.dim(1) != previous.basis_count || tensor.dim(2) != 9)
				throw std::runtime_error("FEC block dimensions do not match the chain.");
			pool->detach_loop(0, sample_count * block.dimension, [&](size_t task) {
				const size_t sample = task / block.dimension;
				const size_t row = task % block.dimension;
				for (size_t pos = tensor.rowptr()[row];
					pos < tensor.rowptr()[row + 1]; ++pos) {
					auto index = tensor.index(pos);
					const size_t old = index[0];
					const size_t letter = index[1];
					for (size_t degree = 0; degree < stride; ++degree) {
						const auto old_value = previous(sample, old, degree);
						if (old_value == 0)
							continue;
						const auto coefficient = scalar_mul(tensor.val(pos), old_value, F);
						if (weights[sample].regular[letter] != 0)
							current(sample, block.offset + row, degree) = scalar_add(
								current(sample, block.offset + row, degree),
								scalar_mul(coefficient,
									weights[sample].regular[letter], F), F);
						if (degree + 1 < stride && weights[sample].logarithmic[letter] != 0)
							current(sample, block.offset + row, degree + 1) = scalar_add(
								current(sample, block.offset + row, degree + 1),
								scalar_mul(coefficient,
									weights[sample].logarithmic[letter], F), F);
					}
				}
			});
			pool->wait();
		}
		previous = std::move(current);
		std::cout << "Projected FEC weight " << weight << "/" << total_weight
			<< ", basis " << previous.basis_count << std::endl;
	}
	return previous;
}

template <typename scalar_t, typename index_t>
void nmhv_limit_add_scaled(sparse_mat<scalar_t, index_t>& target,
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
struct nmhv_limit_amplitude_maps_t {
	sparse_mat<scalar_t, index_t> collinear;
	sparse_mat<scalar_t, index_t> mrk;
};

template <typename scalar_t, typename index_t>
nmhv_limit_amplitude_maps_t<scalar_t, index_t> nmhv_limit_build_amplitude_maps(
	const std::filesystem::path& last0_dir,
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, thread_pool* pool) {
	auto read_matrix = [&](const std::string& name) {
		const auto path = last0_dir / name;
		if (!std::filesystem::exists(path))
			throw std::runtime_error("Input file does not exist: " + path.string());
		return sparse_mat_read_wxf<scalar_t, index_t>(path, F);
	};
	auto legal = read_matrix("Legal.wxf");
	auto E = read_matrix("E.wxf");
	auto Et = read_matrix("Et.wxf");
	if (legal.ncol != E.nrow + Et.nrow)
		throw std::runtime_error("Legal/E/Et dimensions are inconsistent.");

	auto cyclic = nmhv_limit_raw_cyclic_action<scalar_t, index_t>(blocks, F, pool);
	auto cyclic2 = sparse_mat_mul(cyclic, cyclic, F, pool);
	auto E_cyclic2 = sparse_mat_mul(E, cyclic2, F, pool);
	auto Et_cyclic2 = sparse_mat_mul(Et, cyclic2, F, pool);

	// Homogeneous variation of (5.6), cyclically relabeled to the w -> 0 limit.
	auto collinear_E = E;
	nmhv_limit_add_scaled(collinear_E, E_cyclic2,
		nmhv_limit_scalar<scalar_t>(1, 1, F), F, pool);
	sparse_mat<scalar_t, index_t> collinear_Et(Et.nrow, Et.ncol);
	collinear_E.append(std::move(collinear_Et));

	// Homogeneous variation of the E/Et combination in equation (5.16).
	auto mrk_E = E;
	nmhv_limit_add_scaled(mrk_E, E_cyclic2,
		nmhv_limit_scalar<scalar_t>(1, 1, F), F, pool);
	auto mrk_Et = Et;
	nmhv_limit_add_scaled(mrk_Et, Et_cyclic2,
		nmhv_limit_scalar<scalar_t>(-1, 1, F), F, pool);
	mrk_E.append(std::move(mrk_Et));

	auto collinear = sparse_mat_mul(legal, collinear_E, F, pool);
	auto mrk = sparse_mat_mul(legal, mrk_E, F, pool);
	collinear.compress(pool);
	mrk.compress(pool);
	return {std::move(collinear), std::move(mrk)};
}

template <typename scalar_t, typename index_t>
std::vector<std::vector<scalar_t>> nmhv_limit_project(
	const sparse_mat<scalar_t, index_t>& amplitude_map,
	const nmhv_limit_chain_values_t<scalar_t>& chain, size_t degree,
	const field_t& F, thread_pool* pool) {
	if (amplitude_map.ncol != chain.basis_count || degree >= chain.stride)
		throw std::runtime_error("Projection dimensions are inconsistent.");
	std::vector<std::vector<scalar_t>> result(amplitude_map.nrow,
		std::vector<scalar_t>(chain.sample_count));
	pool->detach_loop(0, amplitude_map.nrow * chain.sample_count, [&](size_t task) {
		const size_t row = task / chain.sample_count;
		const size_t sample = task % chain.sample_count;
		scalar_t value = 0;
		for (size_t pos = 0; pos < amplitude_map[row].nnz(); ++pos)
			value = scalar_add(value, scalar_mul(amplitude_map[row][pos],
				chain(sample, amplitude_map[row](pos), degree), F), F);
		result[row][sample] = value;
	});
	pool->wait();
	return result;
}

template <typename scalar_t>
void nmhv_limit_append_columns(std::vector<std::vector<scalar_t>>& target,
	const std::vector<std::vector<scalar_t>>& source) {
	if (target.size() != source.size())
		throw std::runtime_error("Constraint matrices have different row counts.");
	for (size_t row = 0; row < target.size(); ++row)
		target[row].insert(target[row].end(), source[row].begin(), source[row].end());
}
