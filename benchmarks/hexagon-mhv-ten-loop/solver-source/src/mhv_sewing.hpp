#pragma once

#include "mhv.hpp"

// The recurrence/action machinery is component-independent. MHV uses one
// scalar component and its own six-letter seed, never the NMHV Qbar seed.
template <typename index_t>
joint_eet_lec_t<index_t> mhv_lec_seed(const sparse_mat<rat_t, index_t>& cyclic_letter,
	const sparse_mat<rat_t, index_t>& flip_letter, rref_option_t& opt) {
	const field_t F(FIELD_QQ);
	joint_eet_lec_t<index_t> result;
	result.weight = 1;
	result.component_count = 1;
	result.cyclic_letter = cyclic_letter;
	result.flip_letter = flip_letter;
	result.cyclic_component = joint_eet_identity<rat_t, index_t>(1);
	result.flip_component = result.cyclic_component;
	sparse_mat<rat_t, index_t> basis(6, 9);
	for (size_t row = 0; row < 6; ++row) basis[row].push_back(row + 3, 1);
	auto cyclic_image = nmhv_sew_transform_product_basis(basis,
		result.cyclic_component, cyclic_letter, F, &(opt->pool));
	auto flip_image = nmhv_sew_transform_product_basis(basis,
		result.flip_component, flip_letter, F, &(opt->pool));
	result.cyclic = nmhv_sew_restrict_action(basis, cyclic_image, F, opt);
	result.flip = nmhv_sew_restrict_action(basis, flip_image, F, opt);
	joint_eet_adapt_basis(basis, result.cyclic, result.flip, result.blocks, F, opt);
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(basis, &(opt->pool));
	tensor.reshape({basis.nrow, 9, 1});
	result.LEC = sparse_tensor<rat_t, index_t, SPARSE_CSR>(std::move(tensor), &(opt->pool));
	return result;
}

inline std::vector<std::filesystem::path> mhv_lec_files(
	const std::filesystem::path& root, size_t weight) {
	const auto directory = root / joint_eet_weight_name(weight);
	const auto stem = "LEC_" + std::to_string(weight);
	return {directory / (stem + ".wxf"), directory / ("Cyclic_" + stem + ".wxf"),
		directory / ("Flip_" + stem + ".wxf"), directory / (stem + "_blocks.tsv")};
}

template <typename index_t>
std::vector<joint_eet_lec_t<index_t>> mhv_lec_chain(
	const std::filesystem::path& root, const std::filesystem::path& data,
	const std::filesystem::path& condition_path, size_t maximum, rref_option_t& opt) {
	const field_t F(FIELD_QQ);
	const std::vector<std::filesystem::path> fixed{data / "cycrepmat.wxf",
		data / "fliprepmat.wxf", condition_path};
	auto cyclic_letter = nmhv_sew_read_matrix<rat_t, index_t>(fixed[0], F);
	auto flip_letter = nmhv_sew_read_matrix<rat_t, index_t>(fixed[1], F);
	sparse_tensor<rat_t, index_t, SPARSE_COO> condition(
		sparse_tensor_read_wxf<rat_t, index_t>(condition_path, F, &(opt->pool)));
	if (condition.rank() != 3 || condition.dim(0) != 9 || condition.dim(1) != 9)
		throw std::runtime_error("MHV LEC condition must have dimensions {9,9,n}.");
	std::vector<joint_eet_lec_t<index_t>> chain;
	for (size_t weight = 1; weight <= maximum; ++weight) {
		const auto files = mhv_lec_files(root, weight);
		const auto directory = files[0].parent_path();
		auto inputs = fixed;
		if (weight > 1) {
			const auto previous = mhv_lec_files(root, weight - 1);
			inputs.insert(inputs.end(), previous.begin(), previous.end());
		}
		auto signature = [&]() {
			std::ostringstream text;
			text << "mhv-lec-v1\t" << weight << '\n';
			for (const auto& path : inputs) text << hex_sha256_file(path) << '\n';
			for (const auto& path : files) text << hex_sha256_file(path) << '\n';
			return text.str();
		};
		joint_eet_lec_t<index_t> tail;
		if (std::filesystem::exists(directory)) {
			std::ifstream stream(directory / "signature.txt");
			const std::string saved{std::istreambuf_iterator<char>(stream), {}};
			if (saved.empty() || saved != signature())
				throw std::runtime_error("Incomplete or changed MHV LEC state: " + directory.string());
			tail.weight = weight;
			tail.component_count = 1;
			tail.LEC = sparse_tensor_read_wxf<rat_t, index_t>(files[0], F, &(opt->pool));
			tail.cyclic = nmhv_sew_read_matrix<rat_t, index_t>(files[1], F);
			tail.flip = nmhv_sew_read_matrix<rat_t, index_t>(files[2], F);
			tail.blocks = nmhv_sew_read_irrep_axis(files[3]);
			tail.cyclic_letter = cyclic_letter;
			tail.flip_letter = flip_letter;
		}
		else {
			if (weight == 1) tail = mhv_lec_seed(cyclic_letter, flip_letter, opt);
			else {
				auto previous = chain.back();
				tail = extend_joint_eet_lec(std::move(previous), condition, F, opt, false);
			}
			std::filesystem::create_directories(directory);
			nmhv_sew_write_tensor(files[0], tail.LEC);
			nmhv_sew_write_matrix(files[1], tail.cyclic);
			nmhv_sew_write_matrix(files[2], tail.flip);
			joint_eet_write_blocks(files[3], tail.blocks);
			nmhv_sew_publish_bytes(directory / "signature.txt", signature());
			std::ostringstream manifest;
			manifest << "{\n  \"schema\": \"hexagon-mhv-lec-v1\",\n  \"weight\": " << weight
				<< ",\n  \"dimension\": " << tail.dimension() << ",\n  \"nnz\": " << tail.LEC.nnz()
				<< ",\n  \"condition_sha256\": \"" << hex_sha256_file(condition_path) << "\"\n}\n";
			nmhv_sew_publish_bytes(directory / "manifest.json", manifest.str());
		}
		std::cout << "MHV LEC_" << weight << ": dimension " << tail.dimension()
			<< ", nnz " << tail.LEC.nnz() << ", CRC32 " << std::hex << get_file_crc32(files[0])
			<< std::dec << std::endl;
		chain.push_back(std::move(tail));
	}
	return chain;
}

template <typename index_t>
std::vector<sparse_mat<ulong, index_t>> mhv_lec_letter_blocks(
	const sparse_tensor<rat_t, index_t, SPARSE_CSR>& LEC, const field_t& F) {
	std::vector<sparse_mat<ulong, index_t>> result;
	for (size_t letter = 0; letter < 9; ++letter) result.emplace_back(LEC.dim(0), LEC.dim(2));
	for (size_t row = 0; row < LEC.dim(0); ++row)
		for (size_t pos = LEC.rowptr()[row]; pos < LEC.rowptr()[row + 1]; ++pos) {
			const auto index = LEC.index(pos);
			const ulong value = LEC.val(pos) % F.mod;
			if (value) result[index[0]][row].push_back(index[1], value);
		}
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> mhv_sew_bulk(
	const nmhv_sew_irrep_fec_t<index_t>& head, const joint_eet_lec_t<index_t>& tail,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	nmhv_sew_modular_cache_t<index_t>& cache, const field_t& F, rref_option_t& opt) {
	const field_t QQ(FIELD_QQ);
	nmhv_sew_irrep_decomposition_t<index_t> decomposition;
	decomposition.blocks = tail.blocks;
	decomposition.basis = joint_eet_identity<rat_t, index_t>(tail.dimension());
	auto basis = nmhv_sew_modular_reduce(nmhv_sew_irrep_candidates(
		head.blocks, decomposition, QQ, &(opt->pool)), F, &(opt->pool));
	auto first = nmhv_sew_irrep_fec_letter_blocks(head, F, &(opt->pool));
	auto last = mhv_lec_letter_blocks(tail.LEC, F);
	const auto order = nmhv_sew_condition_order(condition);
	if (cache.seams.empty()) cache.seams.resize(order.size());
	std::cout << "MHV CG candidates: " << basis.nrow << ", columns " << basis.ncol << std::endl;
	for (size_t stage = 0; stage < order.size() && basis.nrow; ++stage) {
		auto image = nmhv_sew_factorized_image(basis, first, last, condition, order[stage], F, &(opt->pool));
		std::cout << "MHV seam " << order[stage] << ": " << image.nrow << 'x' << image.ncol
			<< ", nnz " << image.nnz() << std::endl;
		if (image.nnz() == 0) continue;
		auto kernel = nmhv_sew_modular_kernel_rows(image.transpose(), cache.seams[stage], F, opt);
		if (kernel.nrow == 0) return sparse_mat<ulong, index_t>(0, basis.ncol);
		basis = sparse_mat_mul(kernel, basis, F, &(opt->pool));
		nmhv_sew_modular_canonicalize(basis, &(opt->pool));
	}
	return basis;
}

// Contract the backward recurrence directly against global symbol slots.
// Neither limit evaluation nor sewing expands a 9^k letter-word tensor.
template <typename index_t>
nmhv_limit_chain_values_t<ulong> mhv_evaluate_tail(
	const std::vector<joint_eet_lec_t<index_t>>& chain, size_t head_weight,
	bool mrk, size_t samples, size_t degree, uint64_t seed, const field_t& F, thread_pool* pool) {
	const size_t stride = degree + 1, dimension = chain.back().dimension();
	nmhv_limit_chain_values_t<ulong> result{samples, dimension, stride,
		std::vector<ulong>(samples * dimension * stride)};
	pool->detach_loop(0, samples, [&](size_t sample) {
		uint64_t state = seed ^ (0xd1b54a32d192ed03ULL * (sample + 1));
		const auto kind = mrk ? nmhv_limit_kind_t::mrk : nmhv_limit_kind_t::collinear;
		for (size_t w = 1; w <= head_weight; ++w) (void)nmhv_limit_source_weights<ulong>(kind, w, state, F);
		std::vector<nmhv_limit_slot_weights_t<ulong>> slots;
		for (size_t k = 0; k < chain.size(); ++k)
			slots.push_back(nmhv_limit_source_weights<ulong>(kind, head_weight + k + 1, state, F));
		std::vector<ulong> values(stride);
		values[0] = 1;
		for (size_t k = 0; k < chain.size(); ++k) {
			const auto& LEC = chain[k].LEC;
			const auto& slot = slots[chain.size() - k - 1];
			std::vector<ulong> next(LEC.dim(0) * stride);
			for (size_t row = 0; row < LEC.dim(0); ++row)
				for (size_t pos = LEC.rowptr()[row]; pos < LEC.rowptr()[row + 1]; ++pos) {
					const auto index = LEC.index(pos);
					const ulong coefficient = LEC.val(pos) % F.mod;
					for (size_t d = 0; d <= degree; ++d) {
						const ulong value = scalar_mul(coefficient, values[index[1] * stride + d], F);
						next[row * stride + d] = scalar_add(next[row * stride + d], scalar_mul(value, slot.regular[index[0]], F), F);
						if (d < degree) next[row * stride + d + 1] = scalar_add(next[row * stride + d + 1], scalar_mul(value, slot.logarithmic[index[0]], F), F);
					}
				}
			values = std::move(next);
		}
		std::copy(values.begin(), values.end(), result.values.begin() + sample * dimension * stride);
	});
	pool->wait();
	return result;
}

template <typename index_t>
std::vector<nmhv_physical_equation_t> mhv_sew_equations(
	const std::filesystem::path& root, size_t loop,
	const std::vector<joint_eet_lec_t<index_t>>& chain,
	const nmhv_physical_affine_solution_t<index_t>& space, bool mrk,
	size_t samples, uint64_t seed, const nmhv_physical_symbol_table_t& symbols,
	const field_t& F, thread_pool* pool) {
	const size_t degree = mrk ? 2 * loop - 1 : 0;
	const size_t h = 2 * loop - chain.size();
	auto head = nmhv_limit_evaluate_fec_chain<ulong, index_t>(root, h,
		mrk ? nmhv_limit_kind_t::mrk : nmhv_limit_kind_t::collinear, samples, degree, seed, F, pool);
	auto tail = mhv_evaluate_tail(chain, h, mrk, samples, degree, seed, F, pool);
	auto values = joint_eet_physical_project_basis_all(mhv_pack(space, 0, space.origin.size()), head, tail, F, pool);
	std::vector<nmhv_physical_equation_t> result(samples * (degree + 1));
	pool->detach_loop(0, samples, [&](size_t sample) {
		const auto rhs = mrk ? nmhv_physical_mrk_rhs(loop,
			nmhv_physical_sample_mrk_slots(seed, sample, 2 * loop - 1, F), symbols, degree, F)
			: std::vector<ulong>{nmhv_physical_collinear_rhs(loop,
				nmhv_physical_sample_collinear_slots(seed, sample, 2 * loop, F), F)};
		for (size_t d = 0; d <= degree; ++d) {
			auto& equation = result[sample * (degree + 1) + d];
			equation.rhs = scalar_sub(scalar_div(rhs[d], ulong(2), F), values[d][sample][0], F);
			equation.coefficients.assign(values[d][sample].begin() + 1, values[d][sample].end());
		}
	});
	pool->wait();
	return result;
}

template <typename index_t>
sparse_mat<ulong, index_t> mhv_tail_collinear(
	const std::vector<joint_eet_lec_t<index_t>>& chain, const field_t& F, rref_option_t& opt) {
	auto previous = joint_eet_identity<ulong, index_t>(1);
	for (const auto& tail : chain) {
		const auto letters = mhv_lec_letter_blocks(tail.LEC, F);
		sparse_mat<ulong, index_t> image(tail.dimension(), 4 * previous.ncol);
		std::vector<sparse_mat<ulong, index_t>> terms;
		for (const auto& letter : letters) terms.push_back(sparse_mat_mul(letter, previous, F, &(opt->pool)));
		opt->pool.detach_loop(0, image.nrow, [&](size_t row) {
			std::vector<ulong> values(image.ncol);
			for (size_t letter = 0; letter < 9; ++letter)
				for (size_t pos = 0; pos < terms[letter][row].nnz(); ++pos)
					for (size_t target = 0; target < 4; ++target) {
						const int coefficient = nmhv_limit_collinear_map()[letter][target];
						const size_t column = target * previous.ncol + terms[letter][row](pos);
						if (coefficient) values[column] = scalar_add(values[column],
							scalar_mul(terms[letter][row][pos], nmhv_limit_scalar<ulong>(coefficient, 1, F), F), F);
					}
			for (size_t column = 0; column < values.size(); ++column)
				if (values[column]) image[row].push_back(column, values[column]);
		});
		opt->pool.wait();
		previous = nmhv_l0c_v2_reduce_line_image(std::move(image), tail.weight, F, opt).coordinates;
	}
	return previous;
}

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> mhv_sew_limit(
	const std::filesystem::path& root, size_t loop,
	const std::vector<joint_eet_lec_t<index_t>>& chain,
	const nmhv_physical_affine_solution_t<index_t>& parent, bool mrk,
	size_t samples, uint64_t seed, const nmhv_physical_symbol_table_t& symbols,
	const field_t& F, rref_option_t& opt) {
	const auto equations = mhv_sew_equations(root, loop, chain, parent, mrk, samples, seed, symbols, F, &(opt->pool));
	auto child = nmhv_physical_solve_affine_space<index_t>(equations, parent.directions.nrow, F);
	if (!mrk) {
		auto head = nmhv_l0c_v2_d3_exact_collinear<index_t>(root, 2 * loop - chain.size(), F, opt);
		auto tail = mhv_tail_collinear(chain, F, opt);
		auto image = joint_eet_project_natural(parent.directions, head.coordinates, tail, F, &(opt->pool));
		nmhv_sew_pivot_stage_t<index_t> stage;
		auto kernel = nmhv_sew_modular_kernel_rows(image.transpose(), stage, F, opt);
		if (kernel.nrow != child.directions.nrow)
			throw std::runtime_error("MHV collinear samples did not reach exact rank; increase --samples.");
		child.directions = kernel.nrow ? std::move(kernel) : sparse_mat<ulong, index_t>(0, parent.directions.nrow);
	}
	mhv_check_affine(equations, child, F, &(opt->pool));
	mhv_check_affine(mhv_sew_equations(root, loop, chain, parent, mrk, 16,
		seed ^ 0x713ab19ULL, symbols, F, &(opt->pool)), child, F, &(opt->pool));
	auto result = nmhv_physical_compose_affine(parent, child, F, &(opt->pool));
	if (!result.directions.nrow) result.directions = sparse_mat<ulong, index_t>(0, parent.origin.size());
	return result;
}
