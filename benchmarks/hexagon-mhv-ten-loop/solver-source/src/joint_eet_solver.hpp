#pragma once

#include "joint_eet.hpp"
#include "nmhv/fec_provider_d3.hpp"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename index_t>
std::pair<sparse_mat<T, index_t>, sparse_mat<T, index_t>>
joint_eet_expand_last0_components(
	const sparse_mat<T, index_t>& sewing,
	const sparse_mat<T, index_t>& E_generators,
	const sparse_mat<T, index_t>& Et_generators,
	const field_t& F, thread_pool* pool) {
	if (E_generators.ncol != Et_generators.ncol
		|| sewing.ncol != E_generators.nrow + Et_generators.nrow)
		throw std::runtime_error(
			"Last-0 sewing and E/Et generator maps are incompatible.");
	sparse_mat<T, index_t> E_coefficients(sewing.nrow, E_generators.nrow);
	sparse_mat<T, index_t> Et_coefficients(sewing.nrow, Et_generators.nrow);
	pool->detach_loop(0, sewing.nrow, [&](size_t row) {
		for (size_t pos = 0; pos < sewing[row].nnz(); ++pos) {
			const size_t column = sewing[row](pos);
			if (column < E_generators.nrow)
				E_coefficients[row].push_back(
					static_cast<index_t>(column), sewing[row][pos]);
			else
				Et_coefficients[row].push_back(static_cast<index_t>(
					column - E_generators.nrow), sewing[row][pos]);
		}
	});
	pool->wait();
	return {
		sparse_mat_mul(E_coefficients, E_generators, F, pool),
		sparse_mat_mul(Et_coefficients, Et_generators, F, pool)};
}

template <typename index_t>
struct joint_eet_last0_state_t {
	nmhv_sew_pivot_stage_t<index_t> last_entry;
	nmhv_sew_pivot_stage_t<index_t> collinear;
	nmhv_sew_pivot_stage_t<index_t> first_kernel;
	nmhv_sew_pivot_stage_t<index_t> second_kernel;
};

template <typename index_t>
struct joint_eet_last0_prime_result_t {
	sparse_mat<ulong, index_t> basis;
	size_t candidate_dimension = 0;
	size_t last_entry_rank = 0;
	size_t last_entry_nnz = 0;
	size_t line_image_columns = 0;
	size_t line_image_nnz = 0;
	size_t collinear_rank = 0;
	size_t collinear_kernel_nnz = 0;
};

template <typename index_t>
joint_eet_collinear_prime_result_t<index_t> joint_eet_last0_restrict_basis(
	sparse_mat<ulong, index_t> basis,
	const std::filesystem::path& FEC_root, size_t weight, ulong prime,
	joint_eet_last0_state_t<index_t>& state, rref_option_t& opt) {
	field_t F(FIELD_Fp, prime);
	auto* pool = &(opt->pool);
	auto data = nmhv_l0c_d3_final_entry_data<ulong, index_t>(
		FEC_root, weight, F, pool);
	auto generator_map = nmhv_l0c_collinear_generator_map(
		data.generators.E, data.generators.Et, data.cyclic2, F, pool);
	auto projection = nmhv_l0c_v2_d3_exact_collinear<index_t>(
		FEC_root, weight, F, opt);
	auto head_projection = sparse_mat_mul(
		generator_map, projection.coordinates, F, pool);
	sparse_mat<ulong, index_t> tail_projection(1, 1);
	tail_projection[0].push_back(0, 1);
	return joint_eet_restrict_basis(std::move(basis),
		head_projection, tail_projection, state.collinear,
		state.first_kernel, state.second_kernel, 2, F, opt);
}

template <typename index_t>
joint_eet_last0_prime_result_t<index_t> joint_eet_last0_prime(
	const std::filesystem::path& FEC_root, size_t weight,
	const sparse_mat<rat_t, index_t>& exact_relations,
	bool impose_collinear, ulong prime,
	joint_eet_last0_state_t<index_t>& state, rref_option_t& opt) {
	field_t F(FIELD_Fp, prime);
	auto* pool = &(opt->pool);
	auto data = nmhv_l0c_d3_final_entry_data<ulong, index_t>(
		FEC_root, weight, F, pool);
	auto relations = nmhv_sew_modular_reduce(exact_relations, F, pool);
	auto constraints = nmhv_l0c_final_entry_constraints(
		data.letters, data.cyclic, data.cyclic2,
		data.generators.E, data.generators.Et, relations, F, pool);
	joint_eet_last0_prime_result_t<index_t> result;
	result.candidate_dimension = data.generators.E.nrow
		+ data.generators.Et.nrow;
	result.last_entry_nnz = constraints.nnz();
	auto basis = nmhv_sew_modular_kernel_rows(
		std::move(constraints), state.last_entry, F, opt);
	result.last_entry_rank = nmhv_sew_modular_rank(state.last_entry.pivots);
	if (!impose_collinear) {
		result.basis = std::move(basis);
		return result;
	}

	auto restricted = joint_eet_last0_restrict_basis(
		std::move(basis), FEC_root, weight, prime, state, opt);
	result.basis = std::move(restricted.basis);
	result.line_image_columns = restricted.image_ncol;
	result.line_image_nnz = restricted.image_nnz;
	result.collinear_rank = restricted.collinear_rank;
	result.collinear_kernel_nnz = restricted.kernel_nnz;
	return result;
}

template <typename index_t>
nmhv_sew_prime_result_t<index_t> joint_eet_modular_prime_natural(
	const sparse_mat<rat_t, index_t>& exact_candidates,
	std::vector<sparse_mat<ulong, index_t>> FEC_blocks,
	const tail_state_t<rat_t, index_t>& tail,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	ulong prime, nmhv_sew_modular_cache_t<index_t>& cache,
	rref_option_t& opt) {
	field_t F(FIELD_Fp, prime);
	auto* pool = &(opt->pool);
	auto basis = nmhv_sew_modular_reduce(exact_candidates, F, pool);
	auto tail_blocks = nmhv_sew_tail_letter_blocks(tail, F, pool);
	const auto order = nmhv_sew_condition_order(condition);
	if (cache.directions.empty()) {
		cache.directions = order;
		cache.seams.resize(order.size());
	}
	else if (cache.directions != order || cache.seams.size() != order.size()) {
		throw std::runtime_error(
			"Condition direction order changed between primes.");
	}

	nmhv_sew_prime_result_t<index_t> result;
	result.dimensions.push_back(basis.nrow);
	for (size_t stage = 0; stage < order.size(); ++stage) {
		auto image = nmhv_sew_factorized_image(
			basis, FEC_blocks, tail_blocks, condition, order[stage], F, pool);
		result.seam_nnz.push_back(image.nnz());
		if (image.nnz() == 0) {
			cache.seams[stage].ready = true;
			cache.seams[stage].zero = true;
			result.dimensions.push_back(basis.nrow);
			continue;
		}
		auto kernel = nmhv_sew_modular_kernel_rows(
			image.transpose(), cache.seams[stage], F, opt);
		basis = sparse_mat_mul(kernel, basis, F, pool);
		nmhv_sew_modular_canonicalize(basis, pool);
		result.dimensions.push_back(basis.nrow);
	}
	result.basis = std::move(basis);
	return result;
}

inline std::filesystem::path joint_eet_cache_path(
	const std::filesystem::path& directory, ulong prime) {
	return directory / ("p" + std::to_string(prime) + ".bin");
}

inline nmhv_sew_cache_header_t joint_eet_read_cache_header(
	const std::filesystem::path& path) {
	std::ifstream input(path, std::ios::binary);
	nmhv_sew_cache_header_t header{};
	input.read(reinterpret_cast<char*>(&header), sizeof(header));
	if (!input || std::memcmp(header.magic, "NMHVSEW", 7) != 0
		|| header.version != 1)
		throw std::runtime_error("Invalid modular cache header: " + path.string());
	return header;
}

inline void joint_eet_validate_cache_header(
	const nmhv_sew_cache_header_t& header,
	size_t loop, size_t head_weight, size_t tail_weight,
	ulong prime, uint64_t fingerprint,
	const std::filesystem::path& path) {
	if (header.loop != loop || header.head_weight != head_weight
		|| header.tail_weight != tail_weight || header.prime != prime
		|| header.fingerprint != fingerprint)
		throw std::runtime_error(
			"Modular cache metadata mismatch: " + path.string());
	if (header.nrow > static_cast<uint64_t>(
		std::numeric_limits<int32_t>::max())
		|| header.ncol > static_cast<uint64_t>(
			std::numeric_limits<int32_t>::max()))
		throw std::runtime_error(
			"Modular cache exceeds int32 sparse indices: " + path.string());
	const uintmax_t cells = static_cast<uintmax_t>(header.nrow)
		* static_cast<uintmax_t>(header.ncol);
	if (header.nrow != 0
		&& cells / static_cast<uintmax_t>(header.nrow) != header.ncol)
		throw std::runtime_error("Modular cache size overflows: " + path.string());
	if (cells > (std::numeric_limits<uintmax_t>::max() - sizeof(header))
		/ sizeof(ulong))
		throw std::runtime_error("Modular cache byte size overflows: " + path.string());
	const uintmax_t expected = sizeof(header) + cells * sizeof(ulong);
	if (!std::filesystem::is_regular_file(path)
		|| std::filesystem::file_size(path) != expected)
		throw std::runtime_error("Modular cache has the wrong size: " + path.string());
}

template <typename index_t>
sparse_mat<ulong, index_t> joint_eet_read_cache(
	const std::filesystem::path& path,
	size_t loop, size_t head_weight, size_t tail_weight,
	ulong prime, uint64_t fingerprint) {
	const auto header = joint_eet_read_cache_header(path);
	joint_eet_validate_cache_header(header, loop, head_weight,
		tail_weight, prime, fingerprint, path);
	std::ifstream input(path, std::ios::binary);
	input.seekg(sizeof(header), std::ios::beg);
	sparse_mat<ulong, index_t> basis(header.nrow, header.ncol);
	std::vector<ulong> dense(header.ncol);
	for (size_t row = 0; row < header.nrow; ++row) {
		input.read(reinterpret_cast<char*>(dense.data()),
			static_cast<std::streamsize>(dense.size() * sizeof(ulong)));
		if (!input)
			throw std::runtime_error("Truncated modular cache: " + path.string());
		for (size_t column = 0; column < dense.size(); ++column)
			if (dense[column] != 0)
				basis[row].push_back(static_cast<index_t>(column), dense[column]);
	}
	return basis;
}

inline std::vector<ulong> joint_eet_primes(size_t reconstruction_count) {
	if (reconstruction_count == 0)
		throw std::runtime_error("At least one reconstruction prime is required.");
	std::vector<ulong> result;
	result.reserve(reconstruction_count + 2);
	ulong prime = n_nextprime(1ULL << 60, 0);
	for (size_t i = 0; i < reconstruction_count + 2; ++i) {
		result.push_back(prime);
		prime = n_nextprime(prime, 0);
	}
	return result;
}
