#pragma once

#include "joint_eet.hpp"
#include "reduced_d3/hex_hash.hpp"
#include "reduced_d3/platform_resources.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

struct joint_eet_rational_stats_t {
	size_t nnz = 0;
	size_t fractions = 0;
	size_t numerator_bits = 0;
	size_t denominator_bits = 0;
	size_t height_bits = 0;
};

inline size_t joint_eet_parse_size(
	std::string_view text, std::string_view label) {
	try {
		size_t parsed = 0;
		const auto value = std::stoull(std::string(text), &parsed);
		if (parsed != text.size())
			throw std::invalid_argument("trailing characters");
		return static_cast<size_t>(value);
	}
	catch (...) {
		throw std::runtime_error(
			"Invalid " + std::string(label) + ": " + std::string(text));
	}
}

inline std::string joint_eet_weight_name(size_t weight) {
	std::ostringstream output;
	output << 'w' << std::setw(2) << std::setfill('0') << weight;
	return output.str();
}

inline std::filesystem::path joint_eet_write_reconstruction_record(
	const std::filesystem::path& directory,
	const std::vector<ulong>& primes, size_t reconstruction_count) {
	if (reconstruction_count == 0 || primes.size() < reconstruction_count)
		throw std::runtime_error("Invalid reconstruction prime record.");
	std::ostringstream output;
	output << "role\tindex\tprime\n";
	for (size_t i = 0; i < primes.size(); ++i)
		output << (i < reconstruction_count ? "reconstruction" : "validation")
			<< '\t' << (i < reconstruction_count ? i + 1
				: i - reconstruction_count + 1) << '\t' << primes[i] << '\n';
	const auto path = directory / "reconstruction.tsv";
	nmhv_sew_publish_bytes(path, output.str());
	return path;
}

template <typename index_t>
joint_eet_rational_stats_t joint_eet_rational_stats(
	const sparse_mat<rat_t, index_t>& matrix) {
	joint_eet_rational_stats_t result;
	for (size_t row = 0; row < matrix.nrow; ++row) {
		for (size_t pos = 0; pos < matrix[row].nnz(); ++pos) {
			const auto& value = matrix[row][pos];
			++result.nnz;
			result.fractions += !value.is_integer();
			result.numerator_bits = std::max<size_t>(
				result.numerator_bits, value.num().abs().bits());
			result.denominator_bits = std::max<size_t>(
				result.denominator_bits, value.den().bits());
			result.height_bits = std::max<size_t>(
				result.height_bits, value.height_bits());
		}
	}
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> joint_eet_flatten_tensor(
	const sparse_tensor<T, index_t, SPARSE_CSR>& input,
	thread_pool* pool) {
	if(input.rank()==2)return input.to_sparse_mat(pool);
	sparse_tensor<T, index_t, SPARSE_COO> tensor(input);
	std::vector<size_t> trailing;
	for (size_t axis = 1; axis < tensor.rank(); ++axis)
		trailing.push_back(axis);
	tensor.flatten({{0}, trailing});
	return tensor.to_sparse_mat(pool);
}

template <typename T, typename index_t>
sparse_mat<T,index_t> joint_eet_flatten_tensor(sparse_tensor<T,index_t,SPARSE_CSR>&& input,thread_pool* pool) {
	if(input.rank()!=2)return joint_eet_flatten_tensor(static_cast<const sparse_tensor<T,index_t,SPARSE_CSR>&>(input),pool);
	sparse_mat<T,index_t> result(input.dim(0),input.dim(1));
	auto row=[&](size_t i){result[i].reserve(input.rowptr()[i+1]-input.rowptr()[i]);
		for(size_t p=input.rowptr()[i];p<input.rowptr()[i+1];++p)result[i].push_back(input.index(p)[0],std::move(input.val(p)));};
	if(pool){pool->detach_loop(0,result.nrow,row);pool->wait();}else for(size_t i=0;i<result.nrow;++i)row(i);
	return result;
}

inline void joint_eet_write_blocks(const std::filesystem::path& path,
	const std::vector<fec_irrep_axis_block_t>& blocks) {
	std::ostringstream output;
	output << "block\tparity\tirrep\tirrep_dimension\tmultiplicity\toffset\tdimension\n";
	for (const auto& block : blocks) {
		const char* irrep = block.spec.irrep == fec_irrep_t::A1 ? "A1"
			: block.spec.irrep == fec_irrep_t::A2 ? "A2" : "E";
		output << block.spec.name << '\t' << block.spec.parity << '\t'
			<< irrep << '\t' << block.spec.irrep_dimension << '\t'
			<< block.multiplicity << '\t' << block.offset << '\t'
			<< block.dimension << '\n';
	}
	nmhv_sew_publish_bytes(path, output.str());
}

template <typename index_t>
void joint_eet_write_expanded(const std::filesystem::path& path,
	const joint_eet_lec_t<index_t>& tail, thread_pool* pool) {
	if (!tail.expanded)
		return;
	std::vector<size_t> dimensions{tail.dimension()};
	dimensions.insert(dimensions.end(), tail.weight, tail.letter_count);
	dimensions.push_back(tail.component_count);
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(*tail.expanded, pool);
	tensor.reshape(dimensions);
	nmhv_sew_write_tensor(path,
		sparse_tensor<rat_t, index_t, SPARSE_CSR>(std::move(tensor), pool));
}

template <typename index_t>
void joint_eet_write_sewing(const std::filesystem::path& path,
	const sparse_mat<rat_t, index_t>& matrix,
	size_t head_dimension, size_t tail_dimension, thread_pool* pool) {
	if (matrix.ncol != head_dimension * tail_dimension)
		throw std::runtime_error("Sewing matrix width is not head x tail.");
	if (matrix.nrow == 0 || head_dimension == 0 || tail_dimension == 0)
		throw std::runtime_error("Cannot write a zero-dimensional sewing tensor.");
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(matrix, pool);
	tensor.reshape({matrix.nrow, head_dimension, tail_dimension});
	nmhv_sew_write_tensor(path,
		sparse_tensor<rat_t, index_t, SPARSE_CSR>(std::move(tensor), pool));
}

template <typename index_t>
void joint_eet_write_sewing_origin(const std::filesystem::path& path,
	const sparse_mat<rat_t, index_t>& matrix,
	size_t head_dimension, size_t tail_dimension, thread_pool* pool) {
	if (matrix.nrow != 1 || matrix.ncol != head_dimension * tail_dimension)
		throw std::runtime_error("Sewing origin must have shape 1 x (head x tail).");
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(matrix, pool);
	tensor.reshape({head_dimension, tail_dimension});
	nmhv_sew_write_tensor(path,
		sparse_tensor<rat_t, index_t, SPARSE_CSR>(std::move(tensor), pool));
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_expand_sewing_component(
	const sparse_mat<rat_t, index_t>& sewing,
	const sparse_mat<rat_t, index_t>& tail_basis,
	size_t head_dimension, size_t component,
	size_t component_count, thread_pool* pool) {
	if (sewing.ncol % head_dimension != 0
		|| sewing.ncol / head_dimension != tail_basis.nrow
		|| tail_basis.ncol % component_count != 0
		|| component >= component_count)
		throw std::runtime_error("Cannot expand incompatible sewing/tail data.");
	const size_t tail_dimension = tail_basis.nrow;
	const size_t word_count = tail_basis.ncol / component_count;
	sparse_mat<rat_t, index_t> result(
		sewing.nrow, head_dimension * word_count);
	pool->detach_loop(0, sewing.nrow, [&](size_t solution) {
		auto& target = result[solution];
		for (size_t pos = 0; pos < sewing[solution].nnz(); ++pos) {
			const size_t column = sewing[solution](pos);
			const size_t head = column / tail_dimension;
			const size_t tail = column % tail_dimension;
			for (size_t tail_pos = 0;
				tail_pos < tail_basis[tail].nnz(); ++tail_pos) {
				const size_t ambient = tail_basis[tail](tail_pos);
				if (ambient % component_count == component) {
					target.push_back(static_cast<index_t>(
						head * word_count + ambient / component_count),
						sewing[solution][pos] * tail_basis[tail][tail_pos]);
				}
			}
		}
		fec_irrep_merge_entries(target);
	});
	pool->wait();
	return result;
}

template <typename index_t>
void joint_eet_write_last_k(const std::filesystem::path& path,
	const sparse_mat<rat_t, index_t>& matrix,
	size_t head_dimension, size_t tail_weight, thread_pool* pool) {
	std::vector<size_t> dimensions{matrix.nrow, head_dimension};
	dimensions.insert(dimensions.end(), tail_weight, 9);
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(matrix, pool);
	tensor.reshape(dimensions);
	nmhv_sew_write_tensor(path,
		sparse_tensor<rat_t, index_t, SPARSE_CSR>(std::move(tensor), pool));
}

inline uint64_t joint_eet_fingerprint(
	std::string_view algorithm,
	const std::vector<std::filesystem::path>& paths,
	const std::vector<uint64_t>& parameters = {}) {
	uint64_t hash = 1469598103934665603ULL;
	auto append = [&](std::string_view text) {
		for (const unsigned char byte : text) {
			hash ^= byte;
			hash *= 1099511628211ULL;
		}
	};
	append(algorithm);
	for (const uint64_t value : parameters) {
		for (size_t byte = 0; byte < sizeof(value); ++byte) {
			hash ^= (value >> (8 * byte)) & 0xff;
			hash *= 1099511628211ULL;
		}
	}
	for (const auto& path : paths) {
		if (!std::filesystem::is_regular_file(path))
			throw std::runtime_error(
				"Fingerprint input does not exist: " + path.string());
		append(hex_sha256_file(path));
	}
	return hash;
}
