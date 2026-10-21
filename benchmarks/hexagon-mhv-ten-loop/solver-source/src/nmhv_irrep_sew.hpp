#pragma once

#include "fec_irrep_reduced.hpp"
#include "nmhv_sew_modular.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

template <typename index_t>
struct nmhv_sew_irrep_decomposition_t {
	sparse_mat<rat_t, index_t> basis;
	std::vector<fec_irrep_axis_block_t> blocks;
};

template <typename index_t>
struct nmhv_sew_irrep_fec_t {
	std::vector<fec_irrep_axis_block_t> blocks;
	std::vector<sparse_tensor<rat_t, index_t, SPARSE_CSR>> tensors;
	size_t previous_dimension = 0;
	size_t letter_count = 0;
};

inline std::vector<std::string> nmhv_sew_split_tabs(const std::string& line) {
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, '\t'))
		fields.push_back(field);
	return fields;
}

inline std::vector<fec_irrep_axis_block_t> nmhv_sew_read_irrep_axis(
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
		auto fields = nmhv_sew_split_tabs(line);
		if (fields.size() != 7 || blocks.size() >= specs.size())
			throw std::runtime_error("Malformed FEC block metadata row.");
		const auto& spec = specs[blocks.size()];
		if (fields[0] != spec.name || std::stoi(fields[1]) != spec.parity
			|| fields[2] != (spec.irrep == fec_irrep_t::A1 ? "A1"
				: spec.irrep == fec_irrep_t::A2 ? "A2" : "E")
			|| std::stoull(fields[3]) != spec.irrep_dimension) {
			throw std::runtime_error("FEC block metadata uses a noncanonical block order.");
		}
		blocks.push_back({spec, std::stoull(fields[4]),
			std::stoull(fields[5]), std::stoull(fields[6])});
	}
	if (blocks.size() != specs.size())
		throw std::runtime_error("FEC block metadata must contain six canonical blocks.");
	(void)fec_irrep_axis_dimension(blocks);
	return blocks;
}

template <typename index_t>
nmhv_sew_irrep_fec_t<index_t> nmhv_sew_read_irrep_fec(
	const std::filesystem::path& directory, size_t weight,
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, thread_pool* pool) {
	nmhv_sew_irrep_fec_t<index_t> result;
	result.blocks = blocks;
	for (const auto& block : blocks) {
		if (block.dimension == 0) {
			result.tensors.emplace_back();
			continue;
		}
		const auto path = directory /
			("FEC_" + std::to_string(weight) + "_" + block.spec.name + ".wxf");
		std::cout << "Reading file " << path.string() << " ..." << std::endl;
		auto tensor = nmhv_sew_read_tensor<rat_t, index_t>(path, F, pool);
		print_tensor_info(tensor);
		print_crc32(path.string(), path);
		if (tensor.rank() != 3 || tensor.dim(0) != block.dimension
			|| tensor.dim(2) == 0) {
			throw std::runtime_error(
				"FEC irrep block must have layout {block,previous,letter}.");
		}
		if (result.previous_dimension == 0) {
			result.previous_dimension = tensor.dim(1);
			result.letter_count = tensor.dim(2);
		}
		else if (tensor.dim(1) != result.previous_dimension
			|| tensor.dim(2) != result.letter_count) {
			throw std::runtime_error("FEC irrep blocks have inconsistent dimensions.");
		}
		result.tensors.push_back(std::move(tensor));
	}
	if (result.previous_dimension == 0 || result.letter_count == 0)
		throw std::runtime_error("FEC irrep blocks are empty.");
	return result;
}

template <typename index_t>
nmhv_sew_irrep_decomposition_t<index_t> nmhv_sew_irrep_tail_decomposition(
	const tail_state_t<rat_t, index_t>& tail,
	const field_t& F, rref_option_t opt) {
	if (tail.cyclic.nrow != tail.basis.nrow
		|| tail.cyclic.ncol != tail.basis.nrow
		|| tail.flip.nrow != tail.basis.nrow
		|| tail.flip.ncol != tail.basis.nrow) {
		throw std::runtime_error("Tail actions do not match the tail basis.");
	}

	// tail.cyclic stores the column action of C_raw, with C_raw^3=P.
	// Projectors use row actions and the order-three generator c=P C_raw.
	auto cyclic_column = tail.cyclic;
	auto reflection_column = tail.flip;
	auto cyclic = cyclic_column.transpose();
	auto reflection = reflection_column.transpose();
	auto cyclic2 = sparse_mat_mul(cyclic, cyclic, F, &(opt->pool));
	auto parity = sparse_mat_mul(cyclic2, cyclic, F, &(opt->pool));
	auto order_three = sparse_mat_mul(parity, cyclic, F, &(opt->pool));

	nmhv_sew_irrep_decomposition_t<index_t> result;
	for (const auto& spec : hhex_block_specs(fec_symmetry_t::parity_dihedral)) {
		auto basis = dlogmat_irrep_projector_basis(
			spec.parity, spec.irrep, order_three,
			reflection, parity, F, opt);
		const size_t offset = result.basis.nrow;
		if (spec.irrep == fec_irrep_t::E) {
			auto cyclic_image = sparse_mat_mul(
				basis, order_three, F, &(opt->pool));
			auto paired = dlogmat_irrep_linear_combination({
				{&cyclic_image, rat_t(2, 3)}, {&basis, rat_t(1, 3)}});
			result.blocks.push_back({spec, basis.nrow, offset, 2 * basis.nrow});
			result.basis.append(std::move(basis));
			result.basis.append(std::move(paired));
		}
		else {
			result.blocks.push_back({spec, basis.nrow, offset, basis.nrow});
			result.basis.append(std::move(basis));
		}
	}
	if (result.basis.nrow != tail.basis.nrow
		|| result.basis.ncol != tail.basis.nrow) {
		throw std::runtime_error("Tail parity-dihedral decomposition is incomplete.");
	}
	auto inverse = sparse_mat_inverse(result.basis, F, opt);
	if (inverse.nrow != result.basis.nrow || inverse.ncol != result.basis.ncol)
		throw std::runtime_error("Tail parity-dihedral basis is not invertible.");

	std::vector<dlogmat_irrep_block_t> expected_blocks;
	for (const auto& block : result.blocks)
		expected_blocks.push_back({
			block.spec, block.multiplicity, block.offset, block.dimension});
	for (const auto& [source, action] : {
		std::pair<const sparse_mat<rat_t, index_t>*, hhex_dlogmat_action_t>{
			&order_three, hhex_dlogmat_action_t::cyclic},
		std::pair<const sparse_mat<rat_t, index_t>*, hhex_dlogmat_action_t>{
			&reflection, hhex_dlogmat_action_t::reflection},
		std::pair<const sparse_mat<rat_t, index_t>*, hhex_dlogmat_action_t>{
			&parity, hhex_dlogmat_action_t::parity}}) {
		auto adapted = sparse_mat_mul(
			sparse_mat_mul(result.basis, *source, F, &(opt->pool)),
			inverse, F, &(opt->pool));
		auto expected = dlogmat_irrep_expected_action(expected_blocks, action);
		if (!dlogmat_irrep_matrices_equal(std::move(adapted), std::move(expected)))
			throw std::runtime_error("Tail parity-dihedral action convention is inconsistent.");
	}
	return result;
}

template <typename index_t>
sparse_mat<rat_t, index_t> nmhv_sew_irrep_candidates(
	const std::vector<fec_irrep_axis_block_t>& head,
	const nmhv_sew_irrep_decomposition_t<index_t>& tail,
	const field_t& F, thread_pool* pool) {
	const size_t head_dimension = fec_irrep_axis_dimension(head);
	const size_t tail_dimension = fec_irrep_axis_dimension(tail.blocks);
	if (tail.basis.nrow != tail_dimension || tail.basis.ncol != tail_dimension)
		throw std::runtime_error("Tail irrep basis has incompatible dimensions.");
	const auto& singlet = hhex_block_specs(
		fec_symmetry_t::parity_dihedral).front();
	rref_option_t opt;
	opt->verbose = false;
	opt->pool.reset(1);
	auto adapted = make_hhex_codomain_basis(head, tail.blocks, singlet, F, opt);

	sparse_mat<rat_t, index_t> result(
		adapted.embedding.nrow, head_dimension * tail_dimension);
	pool->detach_loop(0, result.nrow, [&](size_t row) {
		for (size_t pos = 0; pos < adapted.embedding[row].nnz(); ++pos) {
			const size_t coordinate = adapted.embedding[row](pos);
			const size_t head_coordinate = coordinate / tail_dimension;
			const size_t tail_coordinate = coordinate % tail_dimension;
			for (size_t tail_pos = 0;
				tail_pos < tail.basis[tail_coordinate].nnz(); ++tail_pos) {
				result[row].push_back(
					static_cast<index_t>(head_coordinate * tail_dimension
						+ tail.basis[tail_coordinate](tail_pos)),
					adapted.embedding[row][pos]
						* tail.basis[tail_coordinate][tail_pos]);
			}
		}
		result[row].sort_indices();
		result[row].compress();
	});
	pool->wait();
	return result;
}

template <typename index_t>
std::vector<sparse_mat<ulong, index_t>> nmhv_sew_irrep_fec_letter_blocks(
	const nmhv_sew_irrep_fec_t<index_t>& FEC,
	const field_t& F, thread_pool* pool) {
	const size_t head_dimension = fec_irrep_axis_dimension(FEC.blocks);
	if (FEC.tensors.size() != FEC.blocks.size())
		throw std::runtime_error("FEC tensors do not match their block metadata.");
	std::vector<sparse_mat<ulong, index_t>> letters;
	letters.reserve(FEC.letter_count);
	for (size_t letter = 0; letter < FEC.letter_count; ++letter)
		letters.emplace_back(head_dimension, FEC.previous_dimension);

	for (size_t block_index = 0; block_index < FEC.blocks.size(); ++block_index) {
		const auto& block = FEC.blocks[block_index];
		const auto& tensor = FEC.tensors[block_index];
		if (block.dimension == 0)
			continue;
		pool->detach_loop(0, block.dimension, [&](size_t row) {
			for (size_t pos = tensor.rowptr()[row];
				pos < tensor.rowptr()[row + 1]; ++pos) {
				auto index = tensor.index(pos);
				const ulong value = tensor.val(pos) % F.mod;
				if (value != 0) {
					letters[index[1]][block.offset + row].push_back(
						index[0], value);
				}
			}
		});
		pool->wait();
	}
	for (auto& letter : letters)
		nmhv_sew_modular_canonicalize(letter, pool);
	return letters;
}

template <typename index_t>
nmhv_sew_reconstruction_result_t<index_t> nmhv_sew_irrep_reconstruct(
	const sparse_mat<rat_t, index_t>& candidates,
	const nmhv_sew_irrep_fec_t<index_t>& FEC,
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
			field_t finite_field(FIELD_Fp, prime);
			auto blocks = nmhv_sew_irrep_fec_letter_blocks(
				FEC, finite_field, &(current_opt->pool));
			return nmhv_sew_modular_prime_from_blocks(
				candidates, std::move(blocks), tail, condition,
				prime, cache, current_opt);
		});
}
