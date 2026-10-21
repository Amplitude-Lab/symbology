#pragma once

#include "nmhv_sew.hpp"
#include "reduced_d3/hex_hash.hpp"

#include <atomic>
#include <cstring>
#include <optional>
#include <ranges>

template <typename index_t>
struct nmhv_sew_pivot_stage_t {
	bool ready = false;
	bool zero = false;
	std::vector<std::vector<pivot_t<index_t>>> pivots;
};

template <typename index_t>
struct nmhv_sew_modular_cache_t {
	std::vector<size_t> directions;
	std::vector<nmhv_sew_pivot_stage_t<index_t>> seams;
	nmhv_sew_pivot_stage_t<index_t> final_basis;
};

template <typename index_t>
struct nmhv_sew_prime_result_t {
	sparse_mat<ulong, index_t> basis;
	std::vector<size_t> dimensions;
	std::vector<size_t> seam_nnz;
};

template <typename index_t>
struct nmhv_sew_reconstruction_result_t {
	sparse_mat<rat_t, index_t> basis;
	std::vector<size_t> dimensions;
	std::vector<size_t> seam_nnz;
	size_t modulus_bits = 0;
};

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_modular_reduce(
	const sparse_mat<rat_t, index_t>& source,
	const field_t& F, thread_pool* pool) {
	sparse_mat<ulong, index_t> result(source.nrow, source.ncol);
	pool->detach_loop(0, source.nrow, [&](size_t row) {
		result[row].reserve(source[row].nnz());
		for (size_t pos = 0; pos < source[row].nnz(); ++pos) {
			const ulong value = source[row][pos] % F.mod;
			if (value != 0)
				result[row].push_back(source[row](pos), value);
		}
	});
	pool->wait();
	return result;
}

template <typename index_t>
void nmhv_sew_modular_canonicalize(
	sparse_mat<ulong, index_t>& matrix, thread_pool* pool) {
	pool->detach_loop(0, matrix.nrow, [&](size_t row) {
		matrix[row].canonicalize();
		matrix[row].sort_indices();
	});
	pool->wait();
}

template <typename index_t>
size_t nmhv_sew_modular_rank(
	const std::vector<std::vector<pivot_t<index_t>>>& pivots) {
	size_t rank = 0;
	for (const auto& group : pivots)
		rank += group.size();
	return rank;
}

template <typename index_t>
void nmhv_sew_modular_rref(
	sparse_mat<ulong, index_t>& matrix,
	nmhv_sew_pivot_stage_t<index_t>& stage,
	const field_t& F, rref_option_t& opt) {
	matrix.compress(&(opt->pool));
	matrix.clear_zero_row();
	if (!stage.ready) {
		stage.pivots = sparse_mat_rref(matrix, F, opt);
		stage.ready = true;
		stage.zero = stage.pivots.empty();
	}
	else if (!stage.zero) {
		sparse_mat_direct_rref(matrix, stage.pivots, F, opt);
	}
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_modular_kernel_rows(
	sparse_mat<ulong, index_t>&& matrix,
	nmhv_sew_pivot_stage_t<index_t>& stage,
	const field_t& F, rref_option_t& opt) {
	const size_t variables = matrix.ncol;
	nmhv_sew_modular_rref(matrix, stage, F, opt);
	if (stage.zero) {
		sparse_mat<ulong, index_t> identity(variables, variables);
		for (size_t row = 0; row < variables; ++row)
			identity[row].push_back(static_cast<index_t>(row), 1);
		return identity;
	}
	auto kernel = sparse_mat_rref_kernel(matrix, stage.pivots, F, opt);
	matrix.clear();
	kernel.transpose_replace();
	nmhv_sew_modular_canonicalize(kernel, &(opt->pool));
	return kernel;
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_modular_pivot_rows(
	sparse_mat<ulong, index_t>&& matrix,
	const std::vector<std::vector<pivot_t<index_t>>>& pivots) {
	std::vector<pivot_t<index_t>> flat;
	for (const auto& group : pivots)
		flat.insert(flat.end(), group.begin(), group.end());
	std::ranges::sort(flat, std::less{}, &pivot_t<index_t>::c);
	sparse_mat<ulong, index_t> result(flat.size(), matrix.ncol);
	for (size_t row = 0; row < flat.size(); ++row)
		result[row] = std::move(matrix[flat[row].r]);
	return result;
}

template <typename index_t>
std::vector<sparse_mat<ulong, index_t>> nmhv_sew_tail_letter_blocks(
	const tail_state_t<rat_t, index_t>& tail,
	const field_t& F, thread_pool* pool) {
	const size_t rest_count = nmhv_sew_checked_power(
		tail.letter_count, tail.weight - 1, "Tail rest dimension");
	const size_t block_columns = nmhv_sew_checked_index_product<index_t>(
		{tail.r_count, rest_count}, "Tail-letter block columns");
	std::vector<sparse_mat<ulong, index_t>> blocks;
	blocks.reserve(tail.letter_count);
	for (size_t letter = 0; letter < tail.letter_count; ++letter)
		blocks.emplace_back(tail.basis.nrow, block_columns);

	pool->detach_loop(0, tail.basis.nrow, [&](size_t row) {
		for (size_t pos = 0; pos < tail.basis[row].nnz(); ++pos) {
			const size_t column = static_cast<size_t>(tail.basis[row](pos));
			const size_t rest = column % rest_count;
			const size_t prefix = column / rest_count;
			const size_t letter = prefix % tail.letter_count;
			const size_t r = prefix / tail.letter_count;
			const ulong value = tail.basis[row][pos] % F.mod;
			if (value != 0) {
				blocks[letter][row].push_back(
					static_cast<index_t>(r * rest_count + rest), value);
			}
		}
	});
	pool->wait();
	return blocks;
}

template <typename index_t>
std::vector<sparse_mat<ulong, index_t>> nmhv_sew_fec_letter_blocks(
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& FEC,
	const field_t& F) {
	if (FEC.rank() != 3)
		throw std::runtime_error("FEC must have layout {head,previous,letter}.");
	std::vector<sparse_mat<ulong, index_t>> blocks;
	blocks.reserve(FEC.dim(2));
	for (size_t letter = 0; letter < FEC.dim(2); ++letter)
		blocks.emplace_back(FEC.dim(0), FEC.dim(1));
	for (size_t pos = 0; pos < FEC.nnz(); ++pos) {
		const ulong value = FEC.val(pos) % F.mod;
		if (value != 0) {
			blocks[FEC.index(pos)[2]][FEC.index(pos)[0]].push_back(
				FEC.index(pos)[1], value);
		}
	}
	return blocks;
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_head_condition_block(
	const std::vector<sparse_mat<ulong, index_t>>& FEC_blocks,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	size_t tail_letter, size_t direction,
	const field_t& F, thread_pool* pool) {
	if (FEC_blocks.empty())
		throw std::runtime_error("FEC letter blocks are empty.");
	sparse_mat<ulong, index_t> result(
		FEC_blocks.front().nrow, FEC_blocks.front().ncol);
	for (size_t pos = 0; pos < condition.nnz(); ++pos) {
		if (static_cast<size_t>(condition.index(pos)[1]) != tail_letter
			|| static_cast<size_t>(condition.index(pos)[2]) != direction)
			continue;
		const size_t head_letter = condition.index(pos)[0];
		const ulong scale = condition.val(pos) % F.mod;
		pool->detach_loop(0, result.nrow, [&](size_t row) {
			sparse_vec_add_mul(
				result[row], FEC_blocks[head_letter][row], scale, F);
		});
		pool->wait();
	}
	nmhv_sew_modular_canonicalize(result, pool);
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_replicate_head_map(
	const sparse_mat<ulong, index_t>& map, size_t tail_count) {
	const size_t rows = nmhv_sew_checked_index_product<index_t>(
		{map.nrow, tail_count}, "Replicated head-map rows");
	const size_t columns = nmhv_sew_checked_index_product<index_t>(
		{map.ncol, tail_count}, "Replicated head-map columns");
	sparse_mat<ulong, index_t> result(rows, columns);
	for (size_t head = 0; head < map.nrow; ++head) {
		for (size_t tail = 0; tail < tail_count; ++tail) {
			auto& row = result[head * tail_count + tail];
			row.reserve(map[head].nnz());
			for (size_t pos = 0; pos < map[head].nnz(); ++pos) {
				row.push_back(static_cast<index_t>(
					static_cast<size_t>(map[head](pos)) * tail_count + tail),
					map[head][pos]);
			}
		}
	}
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_replicate_tail_map(
	const sparse_mat<ulong, index_t>& map, size_t head_count) {
	const size_t rows = nmhv_sew_checked_index_product<index_t>(
		{head_count, map.nrow}, "Replicated tail-map rows");
	const size_t columns = nmhv_sew_checked_index_product<index_t>(
		{head_count, map.ncol}, "Replicated tail-map columns");
	sparse_mat<ulong, index_t> result(rows, columns);
	for (size_t head = 0; head < head_count; ++head) {
		for (size_t tail = 0; tail < map.nrow; ++tail) {
			auto& row = result[head * map.nrow + tail];
			row.reserve(map[tail].nnz());
			for (size_t pos = 0; pos < map[tail].nnz(); ++pos) {
				row.push_back(static_cast<index_t>(
					head * map.ncol + static_cast<size_t>(map[tail](pos))),
					map[tail][pos]);
			}
		}
	}
	return result;
}

template <typename index_t>
void nmhv_sew_modular_add(
	sparse_mat<ulong, index_t>& target,
	const sparse_mat<ulong, index_t>& source,
	const field_t& F, thread_pool* pool) {
	if (target.nrow != source.nrow || target.ncol != source.ncol)
		throw std::runtime_error("Cannot add seam images with different dimensions.");
	pool->detach_loop(0, target.nrow, [&](size_t row) {
		sparse_vec_add_mul(target[row], source[row], 1, F);
	});
	pool->wait();
}

template <typename index_t>
sparse_mat<ulong, index_t> nmhv_sew_factorized_image(
	const sparse_mat<ulong, index_t>& basis,
	const std::vector<sparse_mat<ulong, index_t>>& FEC_blocks,
	const std::vector<sparse_mat<ulong, index_t>>& tail_blocks,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	size_t direction, const field_t& F, thread_pool* pool) {
	if (tail_blocks.empty() || FEC_blocks.size() != tail_blocks.size())
		throw std::runtime_error("Letter blocks are incompatible.");
	const size_t head_count = FEC_blocks.front().nrow;
	const size_t previous_count = FEC_blocks.front().ncol;
	const size_t tail_count = tail_blocks.front().nrow;
	const size_t rest_count = tail_blocks.front().ncol;
	if (basis.ncol != nmhv_sew_checked_index_product<index_t>(
		{head_count, tail_count}, "Natural sewing columns"))
		throw std::runtime_error("Natural sewing basis has incompatible columns.");
	const size_t image_columns = nmhv_sew_checked_index_product<index_t>(
		{previous_count, rest_count}, "Seam image columns");
	sparse_mat<ulong, index_t> image(basis.nrow, image_columns);

	for (size_t letter = 0; letter < tail_blocks.size(); ++letter) {
		auto head = nmhv_sew_head_condition_block(
			FEC_blocks, condition, letter, direction, F, pool);
		if (head.nnz() == 0 || tail_blocks[letter].nnz() == 0)
			continue;
		auto head_map = nmhv_sew_replicate_head_map(head, tail_count);
		auto after_head = sparse_mat_mul(basis, head_map, F, pool);
		head_map.clear();
		auto tail_map = nmhv_sew_replicate_tail_map(
			tail_blocks[letter], previous_count);
		auto term = sparse_mat_mul(after_head, tail_map, F, pool);
		after_head.clear();
		tail_map.clear();
		nmhv_sew_modular_add(image, term, F, pool);
	}
	nmhv_sew_modular_canonicalize(image, pool);
	return image;
}

template <typename index_t>
std::vector<size_t> nmhv_sew_condition_order(
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition) {
	std::vector<size_t> counts(condition.dim(2), 0);
	for (size_t pos = 0; pos < condition.nnz(); ++pos)
		++counts[condition.index(pos)[2]];
	std::vector<size_t> order(condition.dim(2));
	std::iota(order.begin(), order.end(), 0);
	std::ranges::stable_sort(order, [&](size_t left, size_t right) {
		return counts[left] < counts[right];
	});
	return order;
}

template <typename index_t>
nmhv_sew_prime_result_t<index_t> nmhv_sew_modular_prime_from_blocks(
	const sparse_mat<rat_t, index_t>& exact_candidates,
	std::vector<sparse_mat<ulong, index_t>> FEC_blocks,
	const tail_state_t<rat_t, index_t>& tail,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	ulong prime, nmhv_sew_modular_cache_t<index_t>& cache,
	rref_option_t& opt) {
	field_t F(FIELD_Fp, prime);
	thread_pool* pool = &(opt->pool);
	auto basis = nmhv_sew_modular_reduce(exact_candidates, F, pool);
	auto tail_blocks = nmhv_sew_tail_letter_blocks(tail, F, pool);
	const auto order = nmhv_sew_condition_order(condition);
	if (cache.directions.empty()) {
		cache.directions = order;
		cache.seams.resize(order.size());
	}
	else if (cache.directions != order || cache.seams.size() != order.size()) {
		throw std::runtime_error("Condition direction order changed between primes.");
	}

	nmhv_sew_prime_result_t<index_t> result;
	result.dimensions.push_back(basis.nrow);
	for (size_t stage = 0; stage < order.size(); ++stage) {
		auto image = nmhv_sew_factorized_image(
			basis, FEC_blocks, tail_blocks, condition,
			order[stage], F, pool);
		result.seam_nnz.push_back(image.nnz());
		if (image.nnz() == 0) {
			cache.seams[stage].ready = true;
			cache.seams[stage].zero = true;
			result.dimensions.push_back(basis.nrow);
			continue;
		}
		auto image_transpose = image.transpose();
		image.clear();
		auto kernel = nmhv_sew_modular_kernel_rows(
			std::move(image_transpose), cache.seams[stage], F, opt);
		basis = sparse_mat_mul(kernel, basis, F, pool);
		nmhv_sew_modular_canonicalize(basis, pool);
		result.dimensions.push_back(basis.nrow);
	}

	nmhv_sew_modular_rref(basis, cache.final_basis, F, opt);
	result.basis = nmhv_sew_modular_pivot_rows(
		std::move(basis), cache.final_basis.pivots);
	nmhv_sew_modular_canonicalize(result.basis, pool);
	return result;
}

template <typename index_t>
nmhv_sew_prime_result_t<index_t> nmhv_sew_modular_prime(
	const sparse_mat<rat_t, index_t>& exact_candidates,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& FEC,
	const tail_state_t<rat_t, index_t>& tail,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	ulong prime, nmhv_sew_modular_cache_t<index_t>& cache,
	rref_option_t& opt) {
	field_t F(FIELD_Fp, prime);
	auto FEC_blocks = nmhv_sew_fec_letter_blocks(FEC, F);
	return nmhv_sew_modular_prime_from_blocks(
		exact_candidates, std::move(FEC_blocks), tail, condition,
		prime, cache, opt);
}

struct nmhv_sew_cache_header_t {
	char magic[8];
	uint64_t version;
	uint64_t loop;
	uint64_t head_weight;
	uint64_t tail_weight;
	uint64_t nrow;
	uint64_t ncol;
	uint64_t prime;
	uint64_t fingerprint;
};

inline nmhv_sew_cache_header_t nmhv_sew_cache_header(
	size_t loop, size_t head_weight, size_t tail_weight,
	size_t nrow, size_t ncol, ulong prime, uint64_t fingerprint) {
	nmhv_sew_cache_header_t header{};
	std::memcpy(header.magic, "NMHVSEW", 7);
	header.version = 1;
	header.loop = loop;
	header.head_weight = head_weight;
	header.tail_weight = tail_weight;
	header.nrow = nrow;
	header.ncol = ncol;
	header.prime = prime;
	header.fingerprint = fingerprint;
	return header;
}

inline bool nmhv_sew_cache_header_matches(
	const nmhv_sew_cache_header_t& header,
	size_t loop, size_t head_weight, size_t tail_weight,
	size_t nrow, size_t ncol, ulong prime, uint64_t fingerprint) {
	return std::memcmp(header.magic, "NMHVSEW", 7) == 0
		&& header.version == 1 && header.loop == loop
		&& header.head_weight == head_weight
		&& header.tail_weight == tail_weight
		&& header.nrow == nrow && header.ncol == ncol
		&& header.prime == prime && header.fingerprint == fingerprint;
}

inline std::optional<std::filesystem::path> nmhv_sew_find_modular_cache(
	const std::filesystem::path& cache_dir,
	size_t loop, size_t head_weight, size_t tail_weight,
	size_t nrow, size_t ncol, ulong prime, uint64_t fingerprint) {
	for (const auto& entry : std::filesystem::directory_iterator(cache_dir)) {
		if (!entry.is_regular_file())
			continue;
		std::ifstream input(entry.path(), std::ios::binary);
		nmhv_sew_cache_header_t header{};
		input.read(reinterpret_cast<char*>(&header), sizeof(header));
		if (input && nmhv_sew_cache_header_matches(
			header, loop, head_weight, tail_weight,
			nrow, ncol, prime, fingerprint))
			return entry.path();
	}
	return std::nullopt;
}

template <typename index_t>
void nmhv_sew_write_modular_basis(
	const std::filesystem::path& path,
	const sparse_mat<ulong, index_t>& basis, ulong prime,
	size_t loop, size_t head_weight, size_t tail_weight,
	uint64_t fingerprint) {
	if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
		throw std::runtime_error("Modular cache already exists: " + path.string());
	const auto staging =
		nmhv_sew_create_staging_directory(path.parent_path());
	const auto temporary = staging / "payload";
	std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
	if (!output)
		throw std::runtime_error("Cannot write modular cache: " + temporary.string());
	const auto header = nmhv_sew_cache_header(
		loop, head_weight, tail_weight, basis.nrow, basis.ncol,
		prime, fingerprint);
	output.write(reinterpret_cast<const char*>(&header), sizeof(header));
	std::vector<ulong> dense(basis.ncol);
	for (size_t row = 0; row < basis.nrow; ++row) {
		std::fill(dense.begin(), dense.end(), 0);
		for (size_t pos = 0; pos < basis[row].nnz(); ++pos)
			dense[basis[row](pos)] = basis[row][pos];
		output.write(reinterpret_cast<const char*>(dense.data()),
			static_cast<std::streamsize>(dense.size() * sizeof(ulong)));
	}
	output.flush();
	if (!output)
		throw std::runtime_error("Cannot finish modular cache: " + temporary.string());
	output.close();
	nmhv_sew_sync_file(temporary);
	std::error_code error;
	std::filesystem::create_hard_link(temporary, path, error);
	if (error) {
		if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
			throw std::runtime_error(
				"Modular cache already exists: " + path.string());
		throw std::runtime_error(
			"Cannot publish modular cache " + path.string()
				+ ": " + error.message());
	}
	error.clear();
	std::filesystem::remove(temporary, error);
	if (error)
		throw std::runtime_error(
			"Published modular cache but cannot remove staged payload: "
				+ error.message());
	error.clear();
	std::filesystem::remove(staging, error);
	if (error)
		throw std::runtime_error(
			"Published modular cache but cannot remove staging directory: "
				+ error.message());
	std::ofstream checksum(path.string()+".sha256",std::ios::binary);checksum<<hex_sha256_file(path)<<'\n';checksum.close();
	if(!checksum)throw std::runtime_error("Cannot publish modular cache checksum");
}

template <typename index_t>
bool nmhv_sew_reconstruct_basis(
	sparse_mat<rat_t, index_t>& result,
	const std::vector<std::filesystem::path>& paths,
	const std::vector<ulong>& primes,
	const std::vector<std::filesystem::path>& validation_paths,
	const std::vector<ulong>& validation_primes,
	size_t loop, size_t head_weight, size_t tail_weight,
	size_t nrow, size_t ncol, uint64_t fingerprint,
	thread_pool* pool) {
	std::vector<std::ifstream> inputs;
	for (size_t i = 0; i < paths.size(); ++i) {
		inputs.emplace_back(paths[i], std::ios::binary);
		nmhv_sew_cache_header_t header{};
		inputs.back().read(reinterpret_cast<char*>(&header), sizeof(header));
		if (!inputs.back() || !nmhv_sew_cache_header_matches(
			header, loop, head_weight, tail_weight,
			nrow, ncol, primes[i], fingerprint)) {
			throw std::runtime_error("Invalid modular cache header: " + paths[i].string());
		}
	}
	std::vector<std::ifstream> validation_inputs;
	std::vector<field_t> validation_fields;
	for (size_t i = 0; i < validation_paths.size(); ++i) {
		validation_inputs.emplace_back(validation_paths[i], std::ios::binary);
		nmhv_sew_cache_header_t header{};
		validation_inputs.back().read(
			reinterpret_cast<char*>(&header), sizeof(header));
		if (!validation_inputs.back() || !nmhv_sew_cache_header_matches(
			header, loop, head_weight, tail_weight,
			nrow, ncol, validation_primes[i], fingerprint)) {
			throw std::runtime_error(
				"Invalid validation cache header: " + validation_paths[i].string());
		}
		validation_fields.emplace_back(FIELD_Fp, validation_primes[i]);
	}

	std::vector<int_t> moduli(primes.size());
	moduli[0] = primes[0];
	for (size_t i = 1; i < primes.size(); ++i)
		moduli[i] = moduli[i - 1] * primes[i];
	std::vector<std::vector<ulong>> residues(
		primes.size(), std::vector<ulong>(ncol));
	std::vector<std::vector<ulong>> checks(
		validation_primes.size(), std::vector<ulong>(ncol));
	std::vector<rat_t> rational(ncol);
	result = sparse_mat<rat_t, index_t>(nrow, ncol);

	for (size_t row = 0; row < nrow; ++row) {
		for (size_t i = 0; i < inputs.size(); ++i) {
			inputs[i].read(reinterpret_cast<char*>(residues[i].data()),
				static_cast<std::streamsize>(ncol * sizeof(ulong)));
			if (!inputs[i])
				throw std::runtime_error("Truncated modular cache: " + paths[i].string());
		}
		for (size_t i = 0; i < validation_inputs.size(); ++i) {
			validation_inputs[i].read(reinterpret_cast<char*>(checks[i].data()),
				static_cast<std::streamsize>(ncol * sizeof(ulong)));
			if (!validation_inputs[i])
				throw std::runtime_error(
					"Truncated validation cache: " + validation_paths[i].string());
		}
		std::atomic<bool> success = true;
		pool->detach_loop(0, ncol, [&](size_t column) {
			int_t residue = residues[0][column];
			for (size_t i = 1; i < primes.size(); ++i)
				residue = CRT(
					residue, moduli[i - 1], residues[i][column], primes[i]);
			if (!rational_reconstruct(rational[column], residue, moduli.back())) {
				success = false;
				return;
			}
			for (size_t i = 0; i < validation_fields.size(); ++i) {
				if (rational[column] % validation_fields[i].mod != checks[i][column]) {
					success = false;
					return;
				}
			}
		});
		pool->wait();
		if (!success)
			return false;
		for (size_t column = 0; column < ncol; ++column) {
			if (rational[column] != 0)
				result[row].push_back(
					static_cast<index_t>(column), std::move(rational[column]));
		}
		std::cout << "Rational reconstruction rows: " << row + 1
			<< "/" << nrow << std::endl;
	}
	return true;
}

template <typename index_t, typename prime_evaluator_t>
nmhv_sew_reconstruction_result_t<index_t> nmhv_sew_reconstruct_with_evaluator(
	size_t loop, size_t head_weight, size_t tail_weight, size_t prime_count,
	const std::filesystem::path& cache_dir,
	uint64_t fingerprint, rref_option_t& opt,
	prime_evaluator_t&& evaluate_prime, bool reuse_initial = false) {
	if (prime_count == 0)
		throw std::runtime_error("At least one reconstruction prime is required.");
	if (std::filesystem::exists(cache_dir)
		&& !std::filesystem::is_directory(cache_dir))
		throw std::runtime_error(
			"Modular cache path is not a directory: " + cache_dir.string());
	std::filesystem::create_directories(cache_dir);

	nmhv_sew_modular_cache_t<index_t> cache;
	std::vector<std::filesystem::path> paths;
	std::vector<std::filesystem::path> validation_paths;
	std::vector<ulong> primes;
	std::vector<ulong> validation_primes;
	int_t modulus = 1;
	size_t nrow = 0;
	size_t ncol = 0;
	nmhv_sew_reconstruction_result_t<index_t> result;
	ulong prime = n_nextprime(1ULL << 60, 0);
	const bool verbose = opt->verbose;
	opt->verbose = verbose;
	// Stateful evaluators still initialize their pivot charts on the first prime.
	// The direct multiplicity evaluator is stateless and can resume from this receipt.
	nmhv_sew_prime_result_t<index_t> current;
	bool evaluated_initial = false;
	const auto receipt = cache_dir / "initial_prime.tsv";
	if (reuse_initial && std::filesystem::exists(receipt)) {
		std::ifstream digest(receipt.string()+".sha256");std::string expected;digest>>expected;
		if(!digest||expected!=hex_sha256_file(receipt))throw std::runtime_error("Initial-prime receipt checksum mismatch");
		std::ifstream in(receipt);uint64_t id;size_t l,h,k,count;ulong p;
		in>>id>>l>>h>>k>>p>>nrow>>ncol;
		if(!in||id!=fingerprint||l!=loop||h!=head_weight||k!=tail_weight||p!=prime)
			throw std::runtime_error("Initial-prime receipt identity mismatch");
		for(auto* values:{&result.dimensions,&result.seam_nnz}) {
			in>>count;if(!in||count>1000000)throw std::runtime_error("Invalid initial-prime receipt");
			values->resize(count);for(auto& x:*values)in>>x;
		}
		if(!in)throw std::runtime_error("Truncated initial-prime receipt");
	}
	else {
		current = evaluate_prime(prime, cache, opt);evaluated_initial = true;
		nrow = current.basis.nrow;ncol = current.basis.ncol;
		result.dimensions = current.dimensions;result.seam_nnz = current.seam_nnz;
		if(reuse_initial) {
			std::ofstream out(receipt);out<<fingerprint<<' '<<loop<<' '<<head_weight<<' '<<tail_weight<<' '<<prime<<' '<<nrow<<' '<<ncol<<'\n';
			for(const auto* values:{&result.dimensions,&result.seam_nnz}){out<<values->size();for(auto x:*values)out<<' '<<x;out<<'\n';}
			out.close();if(!out)throw std::runtime_error("Cannot write initial-prime receipt");
			std::ofstream digest(receipt.string()+".sha256",std::ios::binary);digest<<hex_sha256_file(receipt)<<'\n';digest.close();
			if(!digest)throw std::runtime_error("Cannot write initial-prime receipt checksum");
		}
	}

	for (size_t iteration = 0; iteration < prime_count; ++iteration) {
		auto cached = nmhv_sew_find_modular_cache(
			cache_dir, loop, head_weight, tail_weight,
			nrow, ncol, prime, fingerprint);
		if (cached) {
			paths.push_back(*cached);
			std::cout << "Reconstruction prime " << prime
				<< ": cached" << std::endl;
		}
		else {
			if (iteration != 0 || !evaluated_initial) {
				opt->verbose = false;
				current = evaluate_prime(prime, cache, opt);
				if (current.basis.nrow != nrow || current.basis.ncol != ncol
					|| current.dimensions != result.dimensions)
					throw std::runtime_error(
						"Modular seam dimensions changed between primes.");
			}
			const auto path = cache_dir /
				("prime_" + std::to_string(iteration + 1) + ".bin");
			if (std::filesystem::exists(path))
				throw std::runtime_error(
					"Conflicting modular cache file: " + path.string());
			nmhv_sew_write_modular_basis(
				path, current.basis, prime, loop, head_weight,
				tail_weight, fingerprint);
			paths.push_back(path);
		}
		primes.push_back(prime);
		modulus *= prime;
		std::cout << "CRT modulus after prime " << iteration + 1
			<< ": 2^" << modulus.bits() << std::endl;
		prime = n_nextprime(prime, 0);
	}

	for (size_t validation = 0; validation < 2; ++validation) {
		auto cached = nmhv_sew_find_modular_cache(
			cache_dir, loop, head_weight, tail_weight,
			nrow, ncol, prime, fingerprint);
		if (cached)
			validation_paths.push_back(*cached);
		else {
			opt->verbose = false;
			current = evaluate_prime(prime, cache, opt);
			if (current.basis.nrow != nrow || current.basis.ncol != ncol
				|| current.dimensions != result.dimensions)
				throw std::runtime_error(
					"Validation-prime seam dimensions changed.");
			const auto path = cache_dir /
				("validation_p" + std::to_string(prime) + ".bin");
			nmhv_sew_write_modular_basis(
				path, current.basis, prime, loop, head_weight,
				tail_weight, fingerprint);
			validation_paths.push_back(path);
		}
		validation_primes.push_back(prime);
		std::cout << "Validation prime " << prime << ": cached" << std::endl;
		prime = n_nextprime(prime, 0);
	}

	if (!nmhv_sew_reconstruct_basis(
		result.basis, paths, primes, validation_paths, validation_primes,
		loop, head_weight, tail_weight, nrow, ncol, fingerprint,
		&(opt->pool))) {
		throw std::runtime_error(
			"Rational reconstruction failed validation; increase --modular-primes.");
	}
	result.modulus_bits = modulus.bits();
	opt->verbose = verbose;
	return result;
}

template <typename index_t>
nmhv_sew_reconstruction_result_t<index_t> nmhv_sew_reconstruct(
	const sparse_mat<rat_t, index_t>& candidates,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& FEC,
	const tail_state_t<rat_t, index_t>& tail,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	size_t loop, size_t head_weight, size_t prime_count,
	const std::filesystem::path& cache_dir,
	uint64_t fingerprint, rref_option_t& opt) {
	return nmhv_sew_reconstruct_with_evaluator<index_t>(
		loop, head_weight, tail.weight, prime_count, cache_dir,
		fingerprint, opt,
		[&](ulong prime, nmhv_sew_modular_cache_t<index_t>& cache,
			rref_option_t& current_opt) {
			return nmhv_sew_modular_prime(
				candidates, FEC, tail, condition, prime, cache, current_opt);
		});
}
