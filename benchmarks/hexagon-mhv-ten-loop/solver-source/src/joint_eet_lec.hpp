#pragma once

#include "nmhv_irrep_sew.hpp"

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

template <typename index_t>
struct joint_eet_lec_t {
	size_t weight = 0;
	size_t letter_count = 9;
	size_t component_count = 5;
	sparse_tensor<rat_t, index_t, SPARSE_CSR> LEC;
	sparse_mat<rat_t, index_t> cyclic;
	sparse_mat<rat_t, index_t> flip;
	std::vector<fec_irrep_axis_block_t> blocks;
	std::optional<sparse_mat<rat_t, index_t>> expanded;
	sparse_mat<rat_t, index_t> cyclic_letter;
	sparse_mat<rat_t, index_t> flip_letter;
	sparse_mat<rat_t, index_t> cyclic_component;
	sparse_mat<rat_t, index_t> flip_component;

	size_t dimension() const {
		return LEC.rank() == 0 ? 0 : LEC.dim(0);
	}
};

template <typename T, typename index_t>
sparse_mat<T, index_t> joint_eet_identity(size_t dimension) {
	sparse_mat<T, index_t> result(dimension, dimension);
	for (size_t row = 0; row < dimension; ++row)
		result[row].push_back(static_cast<index_t>(row), T(1));
	return result;
}

template <typename index_t>
tail_state_t<rat_t, index_t> joint_eet_make_solver_tail(
	const sparse_mat<rat_t, index_t>& component_last,
	const sparse_mat<rat_t, index_t>& cyclic,
	const sparse_mat<rat_t, index_t>& flip,
	size_t weight, size_t letter_count, size_t component_count,
	thread_pool* pool) {
	const size_t word_count = nmhv_sew_checked_power(
		letter_count, weight, "Joint E/Et tail word count");
	const size_t ambient_columns = nmhv_sew_checked_index_product<index_t>(
		{word_count, component_count}, "Joint E/Et tail ambient columns");
	if (component_last.ncol != ambient_columns)
		throw std::runtime_error(
			"Expanded joint E/Et tail must have layout {tail,letters...,component}.");
	if (cyclic.nrow != component_last.nrow || cyclic.ncol != component_last.nrow
		|| flip.nrow != component_last.nrow || flip.ncol != component_last.nrow)
		throw std::runtime_error("Joint E/Et tail actions have incompatible dimensions.");

	tail_state_t<rat_t, index_t> result;
	result.weight = weight;
	result.r_count = component_count;
	result.letter_count = letter_count;
	result.basis = sparse_mat<rat_t, index_t>(
		component_last.nrow, ambient_columns);
	pool->detach_loop(0, component_last.nrow, [&](size_t row) {
		for (size_t pos = 0; pos < component_last[row].nnz(); ++pos) {
			const size_t column = component_last[row](pos);
			const size_t word = column / component_count;
			const size_t component = column % component_count;
			result.basis[row].push_back(static_cast<index_t>(
				component * word_count + word), component_last[row][pos]);
		}
		result.basis[row].sort_indices();
	});
	pool->wait();
	result.cyclic = cyclic;
	result.flip = flip;
	return result;
}

template <typename index_t>
struct joint_eet_tail_chart_t {
	std::vector<size_t> pivot_columns;
	sparse_mat<rat_t, index_t> inverse;
};

template <typename index_t>
joint_eet_tail_chart_t<index_t> joint_eet_make_tail_chart(
	const sparse_mat<rat_t, index_t>& tail_basis,
	const field_t& F, rref_option_t& opt) {
	auto reduced = tail_basis;
	reduced.sort_rows_by_nnz();
	auto pivots = sparse_mat_rref_reconstruct(reduced, opt);
	std::vector<size_t> pivot_columns;
	for (const auto& group : pivots)
		for (const auto& pivot : group)
			pivot_columns.push_back(static_cast<size_t>(pivot.c));
	std::ranges::sort(pivot_columns);
	pivot_columns.erase(
		std::unique(pivot_columns.begin(), pivot_columns.end()),
		pivot_columns.end());
	if (pivot_columns.size() != tail_basis.nrow)
		throw std::runtime_error("The joint E/Et tail basis is not row independent.");
	auto pivot_block = nmhv_sew_select_columns(tail_basis, pivot_columns);
	auto inverse = sparse_mat_inverse(pivot_block, F, opt);
	if (inverse.nrow != tail_basis.nrow || inverse.ncol != tail_basis.nrow)
		throw std::runtime_error("The joint E/Et tail coordinate chart is singular.");
	return {std::move(pivot_columns), std::move(inverse)};
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_apply_tail_chart(
	const sparse_mat<rat_t, index_t>& selected_ambient,
	const joint_eet_tail_chart_t<index_t>& chart,
	const field_t& F, thread_pool* pool) {
	if (selected_ambient.ncol != chart.pivot_columns.size()
		|| chart.inverse.nrow != chart.pivot_columns.size()
		|| chart.inverse.ncol != chart.pivot_columns.size())
		throw std::runtime_error(
			"Selected suffix data and the joint E/Et tail chart are incompatible.");
	return sparse_mat_mul(selected_ambient, chart.inverse, F, pool);
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_factor_tail(
	const sparse_mat<rat_t, index_t>& ambient,
	const sparse_mat<rat_t, index_t>& tail_basis,
	const field_t& F, rref_option_t& opt) {
	if (ambient.ncol != tail_basis.ncol)
		throw std::runtime_error(
			"Ambient suffix data and the joint E/Et tail basis have different widths.");
	auto chart = joint_eet_make_tail_chart(tail_basis, F, opt);
	auto selected = nmhv_sew_select_columns(ambient, chart.pivot_columns);
	auto result = joint_eet_apply_tail_chart(
		selected, chart, F, &(opt->pool));
	if (!nmhv_sew_matrices_equal(
		sparse_mat_mul(result, tail_basis, F, &(opt->pool)), ambient))
		throw std::runtime_error(
			"The suffix data do not lie in the joint E/Et tail space.");
	return result;
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_selected_suffix_data(
	const std::vector<sparse_mat<rat_t, index_t>>& components,
	const std::vector<std::vector<sparse_mat<rat_t, index_t>>>& suffix_steps,
	const std::vector<size_t>& ambient_columns,
	const field_t& F, thread_pool* pool) {
	if (components.empty() || suffix_steps.empty())
		throw std::runtime_error("Suffix assembly requires components and FEC steps.");
	const size_t component_count = components.size();
	const size_t solution_count = components.front().nrow;
	const size_t letter_count = suffix_steps.front().size();
	if (letter_count == 0)
		throw std::runtime_error("Suffix assembly requires a nonempty alphabet.");
	for (const auto& component : components)
		if (component.nrow != solution_count
			|| component.ncol != components.front().ncol)
			throw std::runtime_error("E/Et component dimensions are inconsistent.");
	size_t current_dimension = components.front().ncol;
	for (const auto& step : suffix_steps) {
		if (step.size() != letter_count)
			throw std::runtime_error("FEC steps use different alphabet sizes.");
		const size_t next_dimension = step.front().ncol;
		for (const auto& letter : step)
			if (letter.nrow != current_dimension || letter.ncol != next_dimension)
				throw std::runtime_error("FEC suffix dimensions are inconsistent.");
		current_dimension = next_dimension;
	}
	const size_t head_dimension = current_dimension;

	struct trie_node_t {
		std::vector<std::unique_ptr<trie_node_t>> children;
		std::optional<size_t> output_column;
		explicit trie_node_t(size_t size) : children(size) {}
	};
	std::vector<trie_node_t> roots;
	roots.reserve(component_count);
	for (size_t component = 0; component < component_count; ++component)
		roots.emplace_back(letter_count);
	for (size_t output_column = 0; output_column < ambient_columns.size();
		++output_column) {
		size_t word = ambient_columns[output_column] / component_count;
		const size_t component = ambient_columns[output_column] % component_count;
		auto* node = &roots[component];
		for (size_t depth = 0; depth < suffix_steps.size(); ++depth) {
			const size_t letter = word % letter_count;
			word /= letter_count;
			if (!node->children[letter])
				node->children[letter] = std::make_unique<trie_node_t>(letter_count);
			node = node->children[letter].get();
		}
		if (word != 0 || node->output_column)
			throw std::runtime_error("Invalid or repeated suffix ambient column.");
		node->output_column = output_column;
	}

	sparse_mat<rat_t, index_t> result(
		solution_count * head_dimension, ambient_columns.size());
	std::function<void(const sparse_mat<rat_t, index_t>&,
		const trie_node_t&, size_t)> visit;
	visit = [&](const sparse_mat<rat_t, index_t>& state,
		const trie_node_t& node, size_t depth) {
		if (depth == suffix_steps.size()) {
			if (!node.output_column)
				return;
			for (size_t solution = 0; solution < solution_count; ++solution)
				for (size_t pos = 0; pos < state[solution].nnz(); ++pos)
					result[solution * head_dimension + state[solution](pos)].push_back(
						static_cast<index_t>(*node.output_column), state[solution][pos]);
			return;
		}
		for (size_t letter = 0; letter < letter_count; ++letter) {
			if (!node.children[letter])
				continue;
			auto next = sparse_mat_mul(
				state, suffix_steps[depth][letter], F, pool);
			visit(next, *node.children[letter], depth + 1);
		}
	};
	for (size_t component = 0; component < component_count; ++component)
		visit(components[component], roots[component], 0);
	pool->detach_loop(0, result.nrow, [&](size_t row) {
		result[row].sort_indices();
		result[row].compress();
	});
	pool->wait();
	return result;
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_component_transform(
	const sparse_mat<rat_t, index_t>& cyclic_R,
	const field_t& F, thread_pool* pool) {
	if (cyclic_R.nrow != 5 || cyclic_R.ncol != 5)
		throw std::runtime_error("NMHV R action must have dimensions {5,5}.");
	sparse_mat<rat_t, index_t> E(5, 1);
	E[0].push_back(0, rat_t(1));
	E[3].push_back(0, rat_t(1));
	sparse_mat<rat_t, index_t> Et(5, 1);
	Et[0].push_back(0, rat_t(1));
	Et[3].push_back(0, rat_t(-1));
	auto cyclic_R_copy = cyclic_R;
	auto cyclic_R_transpose = cyclic_R_copy.transpose();
	auto CE = sparse_mat_mul(cyclic_R_transpose, E, F, pool);
	auto C2E = sparse_mat_mul(cyclic_R_transpose, CE, F, pool);
	auto CEt = sparse_mat_mul(cyclic_R_transpose, Et, F, pool);
	auto C2Et = sparse_mat_mul(cyclic_R_transpose, CEt, F, pool);
	for (size_t row = 0; row < 5; ++row) {
		sparse_vec_add_mul(Et[row], C2Et[row], rat_t(-1), F);
		sparse_vec_add_mul(CEt[row], C2Et[row], rat_t(1), F);
	}
	sparse_mat<rat_t, index_t> result(5, 5);
	for (size_t row = 0; row < 5; ++row) {
		for (const auto& [column, vector] : {
			std::pair<size_t, const sparse_mat<rat_t, index_t>*>(0, &E),
			std::pair<size_t, const sparse_mat<rat_t, index_t>*>(1, &CE),
			std::pair<size_t, const sparse_mat<rat_t, index_t>*>(2, &C2E),
			std::pair<size_t, const sparse_mat<rat_t, index_t>*>(3, &Et),
			std::pair<size_t, const sparse_mat<rat_t, index_t>*>(4, &CEt)}) {
			if ((*vector)[row].nnz() != 0)
				result[row].push_back(
					static_cast<index_t>(column), (*vector)[row][0]);
		}
	}
	return result;
}

template <typename index_t>
std::pair<sparse_mat<rat_t, index_t>, sparse_mat<rat_t, index_t>>
joint_eet_adapted_actions(
	const std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, thread_pool* pool) {
	std::vector<dlogmat_irrep_block_t> action_blocks;
	for (const auto& block : blocks)
		action_blocks.push_back({
			block.spec, block.multiplicity, block.offset, block.dimension});
	auto order_three = dlogmat_irrep_expected_action(
		action_blocks, hhex_dlogmat_action_t::cyclic);
	auto parity = dlogmat_irrep_expected_action(
		action_blocks, hhex_dlogmat_action_t::parity);
	auto reflection = dlogmat_irrep_expected_action(
		action_blocks, hhex_dlogmat_action_t::reflection);
	auto raw_cyclic_row = sparse_mat_mul(parity, order_three, F, pool);
	return {raw_cyclic_row.transpose(), reflection.transpose()};
}

template <typename index_t>
bool joint_eet_validate_action_matrices(
	const sparse_mat<rat_t, index_t>& cyclic,
	const sparse_mat<rat_t, index_t>& flip,
	const field_t& F, thread_pool* pool) {
	if (cyclic.nrow != cyclic.ncol || flip.nrow != cyclic.nrow
		|| flip.ncol != cyclic.ncol)
		return false;
	auto identity = joint_eet_identity<rat_t, index_t>(cyclic.nrow);
	auto cyclic2 = sparse_mat_mul(cyclic, cyclic, F, pool);
	auto parity = sparse_mat_mul(cyclic2, cyclic, F, pool);
	auto parity2 = sparse_mat_mul(parity, parity, F, pool);
	auto order_three = sparse_mat_mul(parity, cyclic, F, pool);
	auto order_three2 = sparse_mat_mul(order_three, order_three, F, pool);
	auto order_three3 = sparse_mat_mul(order_three2, order_three, F, pool);
	auto flip2 = sparse_mat_mul(flip, flip, F, pool);
	auto reflected = sparse_mat_mul(
		sparse_mat_mul(flip, order_three, F, pool), flip, F, pool);
	const bool parity_ok = nmhv_sew_matrices_equal(parity2, identity);
	const bool cyclic_ok = nmhv_sew_matrices_equal(order_three3, identity);
	const bool flip_ok = nmhv_sew_matrices_equal(flip2, identity);
	const bool dihedral_ok = nmhv_sew_matrices_equal(reflected, order_three2);
	if (!(parity_ok && cyclic_ok && flip_ok && dihedral_ok))
		std::cout << "action checks: parity=" << parity_ok
			<< ", cyclic=" << cyclic_ok << ", flip=" << flip_ok
			<< ", dihedral=" << dihedral_ok << std::endl;
	return parity_ok && cyclic_ok && flip_ok && dihedral_ok;
}

template <typename index_t>
void joint_eet_adapt_basis(
	sparse_mat<rat_t, index_t>& basis,
	sparse_mat<rat_t, index_t>& cyclic,
	sparse_mat<rat_t, index_t>& flip,
	std::vector<fec_irrep_axis_block_t>& blocks,
	const field_t& F, rref_option_t& opt) {
	tail_state_t<rat_t, index_t> state;
	state.weight = 1;
	state.r_count = 1;
	state.letter_count = basis.ncol;
	state.basis = basis;
	state.cyclic = cyclic;
	state.flip = flip;
	if (!joint_eet_validate_action_matrices(cyclic, flip, F, &(opt->pool)))
		throw std::runtime_error("Joint E/Et actions do not represent parity times D3.");
	auto decomposition = nmhv_sew_irrep_tail_decomposition(state, F, opt);
	basis = sparse_mat_mul(decomposition.basis, basis, F, &(opt->pool));
	blocks = std::move(decomposition.blocks);
	std::tie(cyclic, flip) = joint_eet_adapted_actions<index_t>(
		blocks, F, &(opt->pool));
}

template <typename index_t>
bool validate_joint_eet_actions(
	const joint_eet_lec_t<index_t>& tail,
	const field_t& F, rref_option_t& opt) {
	if (tail.cyclic.nrow != tail.dimension()
		|| tail.cyclic.ncol != tail.dimension()
		|| tail.flip.nrow != tail.dimension()
		|| tail.flip.ncol != tail.dimension())
		return false;
	return joint_eet_validate_action_matrices(
		tail.cyclic, tail.flip, F, &(opt->pool));
}

template <typename index_t>
joint_eet_lec_t<index_t> make_joint_eet_lec_seed(
	const sparse_mat<rat_t, index_t>& qbar_basis,
	const sparse_mat<rat_t, index_t>& cyclic_R,
	const sparse_mat<rat_t, index_t>& cyclic_letter,
	const sparse_mat<rat_t, index_t>& flip_R,
	const sparse_mat<rat_t, index_t>& flip_letter,
	const field_t& F, rref_option_t& opt) {
	if (cyclic_letter.nrow != 9 || cyclic_letter.ncol != 9
		|| flip_letter.nrow != 9 || flip_letter.ncol != 9)
		throw std::runtime_error("Hhex letter actions must have dimensions {9,9}.");
	if (qbar_basis.nrow != 18 || qbar_basis.ncol != 45
		|| cyclic_R.nrow != 5 || cyclic_R.ncol != 5
		|| flip_R.nrow != 5 || flip_R.ncol != 5)
		throw std::runtime_error(
			"NMHV Qbar seed must have layout {18,5*9} with {5,5} R actions.");
	auto basis = qbar_basis;
	nmhv_sew_normalize_rows(basis, &(opt->pool));
	auto cyclic_ambient = nmhv_sew_transform_product_basis(
		basis, cyclic_R, cyclic_letter, F, &(opt->pool));
	auto cyclic = nmhv_sew_restrict_action(
		basis, cyclic_ambient, F, opt);
	auto flip_ambient = nmhv_sew_transform_product_basis(
		basis, flip_R, flip_letter, F, &(opt->pool));
	auto flip = nmhv_sew_restrict_action(
		basis, flip_ambient, F, opt);

	std::vector<fec_irrep_axis_block_t> blocks;
	joint_eet_adapt_basis(basis, cyclic, flip, blocks, F, opt);
	auto component_transform = joint_eet_component_transform<index_t>(
		cyclic_R, F, &(opt->pool));
	auto component_inverse = sparse_mat_inverse(component_transform, F, opt);
	if (component_inverse.nrow != 5 || component_inverse.ncol != 5)
		throw std::runtime_error("The cyclic E/Et component basis is singular.");
	auto component_basis = nmhv_sew_transform_product_basis(
		basis, component_inverse.transpose(),
		joint_eet_identity<rat_t, index_t>(9), F, &(opt->pool));
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(
		component_basis, &(opt->pool));
	tensor.reshape({basis.nrow, 5, 9});
	tensor.transpose_replace({0, 2, 1}, &(opt->pool));
	tensor.sort_indices(&(opt->pool));
	auto expanded_tensor = tensor;
	expanded_tensor.flatten({{0}, {1, 2}});

	joint_eet_lec_t<index_t> result;
	result.weight = 1;
	result.LEC = sparse_tensor<rat_t, index_t, SPARSE_CSR>(
		std::move(tensor), &(opt->pool));
	result.cyclic = std::move(cyclic);
	result.flip = std::move(flip);
	result.blocks = std::move(blocks);
	result.expanded = expanded_tensor.to_sparse_mat(&(opt->pool));
	result.cyclic_letter = cyclic_letter;
	result.flip_letter = flip_letter;
	result.cyclic_component = sparse_mat_mul(
		sparse_mat_mul(component_transform.transpose(), cyclic_R,
			F, &(opt->pool)), component_inverse.transpose(), F, &(opt->pool));
	result.flip_component = sparse_mat_mul(
		sparse_mat_mul(component_transform.transpose(), flip_R,
			F, &(opt->pool)), component_inverse.transpose(), F, &(opt->pool));
	return result;
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_extend_raw(
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	const sparse_tensor<rat_t, index_t, SPARSE_CSR>& LEC,
	const field_t& F, rref_option_t& opt) {
	if (condition.rank() != 3 || condition.dim(0) != 9 || condition.dim(1) != 9)
		throw std::runtime_error("Condition tensor must have layout {9,9,condition}.");
	if (LEC.rank() != 3 || LEC.dim(1) != 9)
		throw std::runtime_error("LEC must have layout {basis,9,previous}.");
	sparse_tensor<rat_t, index_t, SPARSE_COO> old(LEC);
	auto tmp = tensor_contract(condition, old, 1, 1, F, &(opt->pool));
	tmp.flatten({{1, 3}, {0, 2}});
	auto matrix = tmp.to_sparse_mat(&(opt->pool));
	return nmhv_sew_kernel_rows(std::move(matrix), F, opt);
}

template <typename index_t>
sparse_mat<rat_t, index_t> joint_eet_expand(
	const sparse_mat<rat_t, index_t>& recurrence,
	const sparse_mat<rat_t, index_t>& old_expanded,
	size_t old_weight, size_t old_dimension,
	size_t component_count,
	const field_t& F, thread_pool* pool) {
	sparse_tensor<rat_t, index_t, SPARSE_COO> recurrence_tensor(
		recurrence, pool);
	recurrence_tensor.reshape({recurrence.nrow, 9, old_dimension});
	sparse_tensor<rat_t, index_t, SPARSE_COO> old_tensor(old_expanded, pool);
	std::vector<size_t> old_dimensions{old_dimension};
	old_dimensions.insert(old_dimensions.end(), old_weight, 9);
	old_dimensions.push_back(component_count);
	old_tensor.reshape(old_dimensions);
	auto expanded = tensor_contract(
		recurrence_tensor, old_tensor, 2, 0, F, pool);
	std::vector<size_t> ambient_axes;
	for (size_t axis = 1; axis < expanded.rank(); ++axis)
		ambient_axes.push_back(axis);
	expanded.flatten({{0}, ambient_axes});
	return expanded.to_sparse_mat(pool);
}

template <typename index_t>
joint_eet_lec_t<index_t> extend_joint_eet_lec(
	joint_eet_lec_t<index_t>&& old,
	const sparse_tensor<rat_t, index_t, SPARSE_COO>& condition,
	const field_t& F, rref_option_t& opt,
	bool keep_expanded = true) {
	const size_t old_dimension = old.dimension();
	auto basis = joint_eet_extend_raw(condition, old.LEC, F, opt);
	auto cyclic_transformed = nmhv_sew_transform_product_basis(
		basis, old.cyclic_letter, old.cyclic.transpose(), F, &(opt->pool));
	auto cyclic = nmhv_sew_restrict_action(
		basis, cyclic_transformed, F, opt);
	auto flip_transformed = nmhv_sew_transform_product_basis(
		basis, old.flip_letter, old.flip.transpose(), F, &(opt->pool));
	auto flip = nmhv_sew_restrict_action(
		basis, flip_transformed, F, opt);
	std::vector<fec_irrep_axis_block_t> blocks;
	joint_eet_adapt_basis(basis, cyclic, flip, blocks, F, opt);

	joint_eet_lec_t<index_t> result;
	result.weight = old.weight + 1;
	result.cyclic = std::move(cyclic);
	result.flip = std::move(flip);
	result.blocks = std::move(blocks);
	result.component_count = old.component_count;
	result.cyclic_letter = std::move(old.cyclic_letter);
	result.flip_letter = std::move(old.flip_letter);
	result.cyclic_component = std::move(old.cyclic_component);
	result.flip_component = std::move(old.flip_component);
	if (keep_expanded) {
		if (!old.expanded)
			throw std::runtime_error("Previous expanded LEC basis is unavailable.");
		result.expanded = joint_eet_expand(
			basis, *old.expanded, old.weight, old_dimension,
			old.component_count, F, &(opt->pool));
	}
	sparse_tensor<rat_t, index_t, SPARSE_COO> tensor(basis, &(opt->pool));
	tensor.reshape({basis.nrow, 9, old_dimension});
	result.LEC = sparse_tensor<rat_t, index_t, SPARSE_CSR>(
		std::move(tensor), &(opt->pool));
	return result;
}
