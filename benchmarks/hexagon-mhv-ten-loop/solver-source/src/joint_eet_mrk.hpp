#pragma once

#include "joint_eet_io.hpp"
#include "nmhv_limit.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct nmhv_physical_symbol_term_t {
	std::vector<uint8_t> word;
	rat_t coefficient;
};

using nmhv_physical_symbol_table_t = std::map<
	std::pair<size_t, size_t>, std::vector<nmhv_physical_symbol_term_t>>;

struct nmhv_physical_equation_t {
	std::string condition;
	size_t sample = 0;
	std::vector<ulong> coefficients;
	ulong rhs = 0;
	// OPE: endpoint power, twist, helicity, log(T), log(endpoint).
	std::array<int,5> ope_coordinate{-1,-1,0,0,0};
};

template <typename T, typename index_t>
struct nmhv_physical_fec_expansion_block_t {
	size_t offset = 0;
	sparse_mat<T, index_t> expansion;
};

enum class joint_eet_physical_limit_t { collinear, mrk };

inline uint64_t nmhv_physical_splitmix64(uint64_t& state) {
	return nmhv_limit_splitmix64(state);
}

inline std::vector<std::array<ulong, 3>> nmhv_physical_sample_mrk_slots(
	uint64_t seed, size_t sample, size_t weight, const field_t& F) {
	uint64_t state = seed ^ (0xd1b54a32d192ed03ULL * (sample + 1));
	std::vector<std::array<ulong, 3>> slots(weight);
	for (auto& slot : slots) {
		std::array<ulong, 5> target{};
		for (auto& value : target)
			value = nmhv_physical_splitmix64(state) % (F.mod.n - 1) + 1;
		slot = {target[0], target[1], target[2]};
	}
	return slots;
}

inline std::vector<std::array<ulong, 4>> nmhv_physical_sample_collinear_slots(
	uint64_t seed, size_t sample, size_t weight, const field_t& F) {
	uint64_t state = seed ^ (0xd1b54a32d192ed03ULL * (sample + 1));
	std::vector<std::array<ulong, 4>> slots(weight);
	for (auto& slot : slots)
		for (auto& value : slot)
			value = nmhv_physical_splitmix64(state) % (F.mod.n - 1) + 1;
	return slots;
}

inline ulong nmhv_physical_collinear_rhs(size_t loop,
	const std::vector<std::array<ulong, 4>>& slots, const field_t& F) {
	if (slots.size() != 2 * loop)
		throw std::runtime_error("Collinear RHS requires exactly 2*loop slots.");
	if (slots.size() >= 8 * sizeof(size_t))
		throw std::runtime_error("Collinear RHS has too many slots for a bit mask.");
	const size_t full = (size_t(1) << slots.size()) - 1;
	std::vector<std::optional<ulong>> cache(full + 1);
	cache[full] = 1;

	auto matching = [&](auto&& self, size_t mask) -> ulong {
		if (cache[mask])
			return *cache[mask];
		const size_t left = std::countr_one(mask);
		ulong result = 0;
		for (size_t right = left + 1; right < slots.size(); ++right) {
			if ((mask >> right) & 1)
				continue;
			const ulong product = scalar_add(
				scalar_mul(slots[left][0], slots[right][0], F),
				scalar_mul(slots[left][1], slots[right][1], F), F);
			const ulong one_loop = scalar_sub(ulong(0), product, F);
			result = scalar_add(result, scalar_mul(one_loop,
				self(self, mask | (size_t(1) << left)
					| (size_t(1) << right)), F), F);
		}
		cache[mask] = result;
		return result;
	};
	return scalar_mul(ulong(2), matching(matching, 0), F);
}

namespace nmhv_physical_detail {

using polynomial_t = std::vector<ulong>;

inline bool polynomial_zero(const polynomial_t& polynomial) {
	for (const auto value : polynomial)
		if (value != 0)
			return false;
	return true;
}

inline void polynomial_add(polynomial_t& target,
	const polynomial_t& source, const field_t& F) {
	for (size_t degree = 0; degree < target.size(); ++degree)
		target[degree] = scalar_add(target[degree], source[degree], F);
}

inline polynomial_t polynomial_multiply(const polynomial_t& left,
	const polynomial_t& right, size_t max_degree, const field_t& F) {
	polynomial_t result(max_degree + 1);
	for (size_t i = 0; i < left.size() && i <= max_degree; ++i) {
		if (left[i] == 0)
			continue;
		for (size_t j = 0; j < right.size() && i + j <= max_degree; ++j)
			if (right[j] != 0)
				result[i + j] = scalar_add(result[i + j],
					scalar_mul(left[i], right[j], F), F);
	}
	return result;
}

inline void polynomial_scale(polynomial_t& polynomial,
	ulong scale, const field_t& F) {
	for (auto& value : polynomial)
		value = scalar_mul(value, scale, F);
}

inline ulong power(ulong base, size_t exponent, const field_t& F) {
	ulong result = 1;
	while (exponent != 0) {
		if (exponent & 1)
			result = scalar_mul(result, base, F);
		base = scalar_mul(base, base, F);
		exponent >>= 1;
	}
	return result;
}

inline ulong factorial(size_t value, const field_t& F) {
	ulong result = 1;
	for (size_t factor = 2; factor <= value; ++factor)
		result = scalar_mul(result, static_cast<ulong>(factor), F);
	return result;
}

enum class factor_kind_t { symbol, log_tau, one_loop, derivative, phase };

struct factor_t {
	factor_kind_t kind;
	size_t loop = 0;
	size_t degree = 0;
	size_t weight = 0;
};

inline std::vector<size_t> mask_indices(size_t mask) {
	std::vector<size_t> indices;
	while (mask != 0) {
		const size_t index = std::countr_zero(mask);
		indices.push_back(index);
		mask &= mask - 1;
	}
	return indices;
}

inline polynomial_t factor_value(const factor_t& factor, size_t mask,
	const std::vector<std::array<ulong, 3>>& slots,
	const nmhv_physical_symbol_table_t& symbols,
	size_t max_degree, const field_t& F) {
	const auto indices = mask_indices(mask);
	polynomial_t result(max_degree + 1);
	const ulong inv2 = scalar_inv(ulong(2), F);
	if (factor.kind == factor_kind_t::symbol) {
		const auto found = symbols.find({factor.loop, factor.degree});
		if (found == symbols.end())
			throw std::runtime_error("Missing all-order MRK symbol input.");
		// Dense binary symbols share most word prefixes. Evaluate their
		// monomials once instead of multiplying every word from scratch.
		std::vector<ulong> products;
		if(indices.size()<8*sizeof(size_t)&&found->second.size()>=(size_t(1)<<indices.size())/4){
			products.resize(size_t(1)<<indices.size());products[0]=1;size_t count=1;
			for(auto index:indices){for(size_t j=count;j-->0;){auto v=products[j];
				products[2*j]=scalar_mul(v,slots[index][1],F);
				products[2*j+1]=scalar_mul(v,slots[index][2],F);}count*=2;}
		}
		for (const auto& term : found->second) {
			if (term.word.size() != indices.size())
				throw std::runtime_error("MRK symbol word has the wrong length.");
			ulong value = term.coefficient % F.mod;
			size_t word=0;
			for (size_t i = 0; i < indices.size(); ++i) {
				if (term.word[i] > 1)
					throw std::runtime_error("MRK symbol word is not binary.");
				if(products.empty())value = scalar_mul(value,slots[indices[i]][1 + term.word[i]], F);
				else word=2*word+term.word[i];
			}
			if(!products.empty())value=scalar_mul(value,products[word],F);
			result[0] = scalar_add(result[0], value, F);
		}
		return result;
	}
	if (indices.size() != factor.weight)
		throw std::runtime_error("MRK shuffle factor has the wrong weight.");
	if (factor.kind == factor_kind_t::log_tau) {
		result[0] = scalar_sub(scalar_mul(inv2, slots[indices[0]][1], F),
			slots[indices[0]][2], F);
		if (max_degree >= 1)
			result[1] = slots[indices[0]][0];
	}
	else if (factor.kind == factor_kind_t::one_loop) {
		polynomial_t left(max_degree + 1), right(max_degree + 1);
		left[0] = scalar_sub(scalar_mul(inv2, slots[indices[0]][1], F),
			slots[indices[0]][2], F);
		right[0] = scalar_sub(scalar_mul(inv2, slots[indices[1]][1], F),
			slots[indices[1]][2], F);
		if (max_degree >= 1) {
			left[1] = slots[indices[0]][0];
			right[1] = slots[indices[1]][0];
		}
		result = polynomial_multiply(left, right, max_degree, F);
		polynomial_scale(result, nmhv_limit_scalar<ulong>(-2, 1, F), F);
		result[0] = scalar_sub(result[0], scalar_mul(inv2,
			scalar_mul(slots[indices[0]][1], slots[indices[1]][1], F), F), F);
	}
	else if (factor.kind == factor_kind_t::derivative) {
		if (max_degree >= 1)
			result[1] = slots[indices[0]][0];
	}
	else {
		result[0] = scalar_sub(slots[indices[0]][1],
			scalar_mul(ulong(2), slots[indices[0]][2], F), F);
	}
	return result;
}

inline polynomial_t shuffled_product(const std::vector<factor_t>& factors,
	const std::vector<std::array<ulong, 3>>& slots,
	const nmhv_physical_symbol_table_t& symbols,
	size_t max_degree, const field_t& F) {
	const size_t state_count = size_t(1) << slots.size();
	const size_t full = state_count - 1;
	const size_t stride = max_degree + 1;
	std::vector<ulong> states(state_count * stride), next(states.size());
	states[0] = 1;
	for (const auto& factor : factors) {
		std::fill(next.begin(), next.end(), 0);
		std::vector<std::optional<polynomial_t>> factor_cache(state_count);
		for (size_t used = 0; used < state_count; ++used) {
			const auto left = states.data() + used * stride;
			if (std::all_of(left, left + stride, [](ulong x){return x == 0;}))
				continue;
			const size_t remaining = full ^ used;
			// Enumerate only subsets of the required weight. In particular,
			// one/two-letter factors must not scan all 2^remaining subsets.
			auto accumulate = [&](size_t chosen) {
				if (!factor_cache[chosen])
					factor_cache[chosen] = factor_value(factor, chosen,
						slots, symbols, max_degree, F);
				const auto& right = *factor_cache[chosen];
				auto target = next.data() + (used | chosen) * stride;
				for(size_t i=0;i<=max_degree;++i)if(left[i])
					for(size_t j=0;i+j<=max_degree;++j)if(right[j])
						target[i+j]=scalar_add(target[i+j],scalar_mul(left[i],right[j],F),F);
			};
			auto subsets = [&](auto&& self,size_t available,size_t count,size_t chosen)->void {
				if(!count){accumulate(chosen);return;}
				size_t left=std::popcount(available);
				if(left<count)return;
				if(left==count){accumulate(chosen|available);return;}
				while(left>=count){size_t bit=available&(~available+1);available^=bit;--left;
					self(self,available,count-1,chosen|bit);}
			};
			subsets(subsets,remaining,factor.weight,0);
		}
		states.swap(next);
	}
	return polynomial_t(states.begin() + full * stride, states.end());
}

} // namespace nmhv_physical_detail

inline std::vector<ulong> nmhv_physical_mrk_rhs(size_t loop,
	const std::vector<std::array<ulong, 3>>& slots,
	const nmhv_physical_symbol_table_t& symbols,
	size_t max_degree, const field_t& F) {
	using namespace nmhv_physical_detail;
	if (slots.size() != 2 * loop - 1)
		throw std::runtime_error("MRK RHS requires exactly 2*loop-1 slots.");
	if (slots.size() >= 8 * sizeof(size_t))
		throw std::runtime_error("MRK RHS has too many slots for a bit mask.");
	polynomial_t result(max_degree + 1);
	const ulong minus_quarter = nmhv_limit_scalar<ulong>(-1, 4, F);
	for (size_t mrk_loop = 1; mrk_loop <= loop; ++mrk_loop) {
		const size_t mhv_power = loop - mrk_loop;
		const ulong scale = scalar_mul(minus_quarter,
			scalar_mul(power(ulong(2), mrk_loop, F),
				scalar_inv(factorial(mhv_power, F), F), F), F);
		for (size_t degree = 0; degree < mrk_loop; ++degree) {
			std::vector<factor_t> factors;
			for (size_t i = 0; i < mhv_power; ++i)
				factors.push_back({factor_kind_t::one_loop, 0, 0, 2});
			factors.push_back({factor_kind_t::symbol, mrk_loop, degree,
				2 * mrk_loop - 1 - degree});
			for (size_t i = 0; i < degree; ++i)
				factors.push_back({factor_kind_t::log_tau, 0, 0, 1});
			auto contribution = shuffled_product(
				factors, slots, symbols, max_degree, F);
			polynomial_scale(contribution, scale, F);
			polynomial_add(result, contribution, F);
		}
	}
	std::vector<factor_t> factors;
	for (size_t i = 0; i + 1 < loop; ++i)
		factors.push_back({factor_kind_t::one_loop, 0, 0, 2});
	factors.push_back({factor_kind_t::derivative, 0, 0, 1});
	auto contribution = shuffled_product(factors, slots, symbols, max_degree, F);
	polynomial_scale(contribution,
		scalar_inv(factorial(loop - 1, F), F), F);
	polynomial_add(result, contribution, F);
	factors.back() = {factor_kind_t::phase, 0, 0, 1};
	contribution = shuffled_product(factors, slots, symbols, max_degree, F);
	polynomial_scale(contribution, scalar_inv(
		scalar_mul(ulong(2), factorial(loop - 1, F), F), F), F);
	polynomial_add(result, contribution, F);
	polynomial_scale(result, ulong(2), F);
	return result;
}

template <typename index_t>
struct nmhv_physical_affine_solution_t {
	std::vector<ulong> origin;
	sparse_mat<ulong, index_t> directions;
	size_t rank = 0;
	std::vector<size_t> pivot_columns;
};

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> nmhv_physical_solve_affine_space(
	const std::vector<nmhv_physical_equation_t>& equations,
	size_t variable_count, const field_t& F) {
	std::vector<std::optional<std::vector<ulong>>> pivots(variable_count);
	size_t rank = 0;
	for (const auto& equation : equations) {
		if (equation.coefficients.size() != variable_count)
			throw std::runtime_error("Physical equation has the wrong width.");
		auto row = equation.coefficients;
		row.push_back(equation.rhs);
		for (size_t pivot = 0; pivot < variable_count; ++pivot) {
			if (!pivots[pivot] || row[pivot] == 0)
				continue;
			const ulong scale = row[pivot];
			for (size_t column = pivot; column <= variable_count; ++column)
				row[column] = scalar_sub(row[column],
					scalar_mul(scale, (*pivots[pivot])[column], F), F);
		}
		size_t pivot = 0;
		while (pivot < variable_count && row[pivot] == 0)
			++pivot;
		if (pivot == variable_count) {
			if (row.back() != 0)
				throw std::runtime_error("Inconsistent " + equation.condition
					+ " sample " + std::to_string(equation.sample) + ".");
			continue;
		}
		const ulong inverse = scalar_inv(row[pivot], F);
		for (size_t column = pivot; column <= variable_count; ++column)
			row[column] = scalar_mul(row[column], inverse, F);
		for (auto& old_row : pivots) {
			if (!old_row || (*old_row)[pivot] == 0)
				continue;
			const ulong scale = (*old_row)[pivot];
			for (size_t column = pivot; column <= variable_count; ++column)
				(*old_row)[column] = scalar_sub((*old_row)[column],
					scalar_mul(scale, row[column], F), F);
		}
		pivots[pivot] = std::move(row);
		++rank;
	}

	nmhv_physical_affine_solution_t<index_t> result;
	result.origin.resize(variable_count);
	result.rank = rank;
	std::vector<size_t> free_columns;
	for (size_t column = 0; column < variable_count; ++column) {
		if (pivots[column]) {
			result.pivot_columns.push_back(column);
			result.origin[column] = (*pivots[column]).back();
		}
		else
			free_columns.push_back(column);
	}
	result.directions = sparse_mat<ulong, index_t>(
		free_columns.size(), variable_count);
	for (size_t row = 0; row < free_columns.size(); ++row) {
		const size_t free_column = free_columns[row];
		for (size_t column = 0; column < variable_count; ++column) {
			ulong value = 0;
			if (column == free_column)
				value = 1;
			else if (pivots[column] && (*pivots[column])[free_column] != 0)
				value = scalar_sub(
					ulong(0), (*pivots[column])[free_column], F);
			if (value != 0)
				result.directions[row].push_back(
					static_cast<index_t>(column), value);
		}
	}
	return result;
}

template <typename index_t>
std::vector<nmhv_physical_equation_t> nmhv_physical_restrict_equations(
	const std::vector<nmhv_physical_equation_t>& equations,
	const std::vector<ulong>& origin,
	const sparse_mat<ulong, index_t>& directions,
	const field_t& F, thread_pool* pool) {
	if (directions.ncol != origin.size())
		throw std::runtime_error("Affine origin and directions have different widths.");
	std::vector<nmhv_physical_equation_t> result(equations.size());
	pool->detach_loop(0, equations.size(), [&](size_t row) {
		const auto& equation = equations[row];
		if (equation.coefficients.size() != origin.size())
			throw std::runtime_error("Physical equation has the wrong width.");
		result[row].condition = equation.condition;
		result[row].sample = equation.sample;
		result[row].coefficients.resize(directions.nrow);
		ulong origin_value = 0;
		for (size_t column = 0; column < origin.size(); ++column)
			origin_value = scalar_add(origin_value,
				scalar_mul(equation.coefficients[column], origin[column], F), F);
		result[row].rhs = scalar_sub(equation.rhs, origin_value, F);
		for (size_t direction = 0; direction < directions.nrow; ++direction) {
			ulong value = 0;
			for (size_t pos = 0; pos < directions[direction].nnz(); ++pos)
				value = scalar_add(value, scalar_mul(directions[direction][pos],
					equation.coefficients[directions[direction](pos)], F), F);
			result[row].coefficients[direction] = value;
		}
	});
	pool->wait();
	return result;
}

template <typename index_t>
nmhv_physical_affine_solution_t<index_t> nmhv_physical_compose_affine(
	const nmhv_physical_affine_solution_t<index_t>& parent,
	const nmhv_physical_affine_solution_t<index_t>& child,
	const field_t& F, thread_pool* pool) {
	if (parent.directions.ncol != parent.origin.size()
		|| child.origin.size() != parent.directions.nrow
		|| child.directions.ncol != parent.directions.nrow)
		throw std::runtime_error("Cannot compose incompatible affine spaces.");
	nmhv_physical_affine_solution_t<index_t> result;
	result.origin = parent.origin;
	for (size_t direction = 0; direction < parent.directions.nrow; ++direction) {
		if (child.origin[direction] == 0)
			continue;
		for (size_t pos = 0; pos < parent.directions[direction].nnz(); ++pos) {
			const size_t column = parent.directions[direction](pos);
			result.origin[column] = scalar_add(result.origin[column],
				scalar_mul(child.origin[direction],
					parent.directions[direction][pos], F), F);
		}
	}
	result.directions = sparse_mat_mul(
		child.directions, parent.directions, F, pool);
	result.rank = result.origin.size() - result.directions.nrow;
	return result;
}

inline std::vector<ulong> nmhv_physical_solve_affine(
	const std::vector<nmhv_physical_equation_t>& equations,
	size_t variable_count, const field_t& F) {
	auto result = nmhv_physical_solve_affine_space<int32_t>(
		equations, variable_count, F);
	if (result.rank != variable_count)
		throw std::runtime_error("Physical affine system rank "
			+ std::to_string(result.rank) + " is below "
			+ std::to_string(variable_count) + ".");
	return std::move(result.origin);
}

template <typename T>
std::array<T, 5> joint_eet_physical_component_weights(
	joint_eet_physical_limit_t kind, const field_t& F) {
	std::array<T, 5> result{};
	result[0] = nmhv_limit_scalar<T>(1, 1, F);
	result[2] = nmhv_limit_scalar<T>(1, 1, F);
	// The stored components eliminate C^2 Et through C^2 Et = -Et + C Et.
	if (kind == joint_eet_physical_limit_t::mrk) {
		result[3] = nmhv_limit_scalar<T>(2, 1, F);
		result[4] = nmhv_limit_scalar<T>(-1, 1, F);
	}
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> joint_eet_physical_expand_component(
	const sparse_mat<T, index_t>& sewing,
	const sparse_mat<T, index_t>& tail_basis,
	size_t head_dimension, size_t component,
	size_t component_count, const field_t& F, thread_pool* pool) {
	if (sewing.ncol % head_dimension != 0
		|| sewing.ncol / head_dimension != tail_basis.nrow
		|| tail_basis.ncol % component_count != 0
		|| component >= component_count)
		throw std::runtime_error("Cannot expand incompatible sewing/tail data.");
	const size_t tail_dimension = tail_basis.nrow;
	const size_t word_count = tail_basis.ncol / component_count;
	sparse_mat<T, index_t> result(sewing.nrow, head_dimension * word_count);
	pool->detach_loop(0, sewing.nrow, [&](size_t solution) {
		std::map<index_t, T> entries;
		for (size_t pos = 0; pos < sewing[solution].nnz(); ++pos) {
			const size_t column = sewing[solution](pos);
			const size_t head = column / tail_dimension;
			const size_t tail = column % tail_dimension;
			for (size_t tail_pos = 0;
				tail_pos < tail_basis[tail].nnz(); ++tail_pos) {
				const size_t ambient = tail_basis[tail](tail_pos);
				if (ambient % component_count != component)
					continue;
				const auto output = static_cast<index_t>(
					head * word_count + ambient / component_count);
				const T value = scalar_mul(
					sewing[solution][pos], tail_basis[tail][tail_pos], F);
				auto found = entries.find(output);
				if (found == entries.end())
					entries.emplace(output, value);
				else
					found->second = scalar_add(found->second, value, F);
			}
		}
		for (const auto& [column, value] : entries)
			if (value != 0)
				result[solution].push_back(column, value);
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
nmhv_limit_chain_values_t<T> joint_eet_physical_evaluate_tail(
	const sparse_mat<T, index_t>& expanded_tail,
	size_t head_weight, size_t tail_weight,
	joint_eet_physical_limit_t kind, size_t sample_count,
	size_t max_log_degree, uint64_t seed,
	const field_t& F, thread_pool* pool) {
	constexpr size_t letter_count = 9;
	constexpr size_t component_count = 5;
	size_t word_count = 1;
	for (size_t i = 0; i < tail_weight; ++i)
		word_count *= letter_count;
	if (expanded_tail.ncol != word_count * component_count)
		throw std::runtime_error(
			"Expanded joint E/Et tail has incompatible dimensions.");
	const size_t stride = max_log_degree + 1;
	nmhv_limit_chain_values_t<T> result{sample_count, expanded_tail.nrow,
		stride, std::vector<T>(sample_count * expanded_tail.nrow * stride)};
	const auto components = joint_eet_physical_component_weights<T>(kind, F);

	pool->detach_loop(0, sample_count * expanded_tail.nrow, [&](size_t task) {
		const size_t sample = task / expanded_tail.nrow;
		const size_t basis = task % expanded_tail.nrow;
		uint64_t state = seed ^ (0xd1b54a32d192ed03ULL * (sample + 1));
		for (size_t weight = 1; weight <= head_weight; ++weight)
			(void)nmhv_limit_source_weights<T>(
				kind == joint_eet_physical_limit_t::collinear
					? nmhv_limit_kind_t::collinear : nmhv_limit_kind_t::mrk,
				weight, state, F);
		std::vector<nmhv_limit_slot_weights_t<T>> slots(tail_weight);
		for (size_t depth = 0; depth < tail_weight; ++depth)
			slots[depth] = nmhv_limit_source_weights<T>(
				kind == joint_eet_physical_limit_t::collinear
					? nmhv_limit_kind_t::collinear : nmhv_limit_kind_t::mrk,
				head_weight + depth + 1, state, F);

		for (size_t pos = 0; pos < expanded_tail[basis].nnz(); ++pos) {
			size_t word = expanded_tail[basis](pos) / component_count;
			const size_t component = expanded_tail[basis](pos) % component_count;
			if (components[component] == 0)
				continue;
			std::vector<size_t> letters(tail_weight);
			for (size_t depth = tail_weight; depth-- != 0;) {
				letters[depth] = word % letter_count;
				word /= letter_count;
			}
			std::vector<T> polynomial(stride);
			polynomial[0] = scalar_mul(
				expanded_tail[basis][pos], components[component], F);
			for (size_t depth = 0; depth < tail_weight; ++depth) {
				std::vector<T> next(stride);
				for (size_t degree = 0; degree < stride; ++degree) {
					if (polynomial[degree] == 0)
						continue;
					const size_t letter = letters[depth];
					if (slots[depth].regular[letter] != 0)
						next[degree] = scalar_add(next[degree],
							scalar_mul(polynomial[degree],
								slots[depth].regular[letter], F), F);
					if (degree + 1 < stride
						&& slots[depth].logarithmic[letter] != 0)
						next[degree + 1] = scalar_add(next[degree + 1],
							scalar_mul(polynomial[degree],
								slots[depth].logarithmic[letter], F), F);
				}
				polynomial = std::move(next);
			}
			for (size_t degree = 0; degree < stride; ++degree)
				result(sample, basis, degree) = scalar_add(
					result(sample, basis, degree), polynomial[degree], F);
		}
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
std::vector<std::vector<T>> joint_eet_physical_project_basis(
	const sparse_mat<T, index_t>& basis,
	const nmhv_limit_chain_values_t<T>& head,
	const nmhv_limit_chain_values_t<T>& tail,
	size_t degree, const field_t& F, thread_pool* pool) {
	if (head.sample_count != tail.sample_count || head.stride != tail.stride
		|| degree >= head.stride
		|| basis.ncol != head.basis_count * tail.basis_count)
		throw std::runtime_error("Physical projection dimensions are inconsistent.");
	std::vector<std::vector<T>> result(head.sample_count,
		std::vector<T>(basis.nrow));
	pool->detach_loop(0, head.sample_count * basis.nrow, [&](size_t task) {
		const size_t sample = task / basis.nrow;
		const size_t solution = task % basis.nrow;
		T value = 0;
		for (size_t pos = 0; pos < basis[solution].nnz(); ++pos) {
			const size_t column = basis[solution](pos);
			const size_t head_basis = column / tail.basis_count;
			const size_t tail_basis = column % tail.basis_count;
			T product = 0;
			for (size_t head_degree = 0;
				head_degree <= degree; ++head_degree)
				product = scalar_add(product,
					scalar_mul(head(sample, head_basis, head_degree),
						tail(sample, tail_basis, degree - head_degree), F), F);
			value = scalar_add(value,
				scalar_mul(basis[solution][pos], product, F), F);
		}
		result[sample][solution] = value;
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
std::vector<std::vector<std::vector<T>>> joint_eet_physical_project_basis_all(
	const sparse_mat<T, index_t>& basis,
	const nmhv_limit_chain_values_t<T>& head,
	const nmhv_limit_chain_values_t<T>& tail,
	const field_t& F, thread_pool* pool) {
	if (head.sample_count != tail.sample_count || head.stride != tail.stride
		|| basis.ncol != head.basis_count * tail.basis_count)
		throw std::runtime_error("Physical projection dimensions are inconsistent.");
	std::vector<std::vector<std::vector<T>>> result(head.stride,
		std::vector<std::vector<T>>(head.sample_count,
			std::vector<T>(basis.nrow)));
	pool->detach_loop(0, head.sample_count * basis.nrow, [&](size_t task) {
		const size_t sample = task / basis.nrow;
		const size_t solution = task % basis.nrow;
		for (size_t pos = 0; pos < basis[solution].nnz(); ++pos) {
			const size_t column = basis[solution](pos);
			const size_t head_basis = column / tail.basis_count;
			const size_t tail_basis = column % tail.basis_count;
			for (size_t head_degree = 0;
				head_degree < head.stride; ++head_degree) {
				const T head_value = head(sample, head_basis, head_degree);
				if (head_value == 0)
					continue;
				for (size_t tail_degree = 0;
					head_degree + tail_degree < head.stride; ++tail_degree) {
					const T tail_value = tail(sample, tail_basis, tail_degree);
					if (tail_value == 0)
						continue;
					auto& value = result[head_degree + tail_degree][sample][solution];
					value = scalar_add(value, scalar_mul(basis[solution][pos],
						scalar_mul(head_value, tail_value, F), F), F);
				}
			}
		}
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_physical_expand_fec_step(
	const sparse_mat<T, index_t>& coefficients,
	const std::vector<nmhv_physical_fec_expansion_block_t<T, index_t>>& blocks,
	size_t letter_count, const field_t& F, thread_pool* pool) {
	if (letter_count == 0 || blocks.empty())
		throw std::runtime_error("FEC expansion requires blocks and letters.");
	const size_t previous_count = blocks.front().expansion.ncol / letter_count;
	if (previous_count == 0
		|| blocks.front().expansion.ncol != previous_count * letter_count)
		throw std::runtime_error("FEC expansion block has incompatible columns.");
	sparse_mat<T, index_t> result(
		coefficients.nrow * letter_count, previous_count);
	for (const auto& block : blocks) {
		if (block.offset + block.expansion.nrow > coefficients.ncol
			|| block.expansion.ncol != previous_count * letter_count)
			throw std::runtime_error("FEC expansion block has incompatible dimensions.");
		sparse_mat<T, index_t> local(coefficients.nrow, block.expansion.nrow);
		for (size_t suffix = 0; suffix < coefficients.nrow; ++suffix)
			for (size_t pos = 0; pos < coefficients[suffix].nnz(); ++pos) {
				const size_t column = coefficients[suffix](pos);
				if (column >= block.offset
					&& column < block.offset + block.expansion.nrow)
					local[suffix].push_back(
						static_cast<index_t>(column - block.offset),
						coefficients[suffix][pos]);
			}
		if (local.nnz() == 0)
			continue;
		auto product = sparse_mat_mul(local, block.expansion, F, pool);
		sparse_mat<T, index_t> contribution(result.nrow, result.ncol);
		for (size_t suffix = 0; suffix < product.nrow; ++suffix)
			for (size_t pos = 0; pos < product[suffix].nnz(); ++pos) {
				const size_t column = product[suffix](pos);
				const size_t previous = column / letter_count;
				const size_t letter = column % letter_count;
				contribution[letter * coefficients.nrow + suffix].push_back(
					static_cast<index_t>(previous), product[suffix][pos]);
			}
		for (size_t row = 0; row < result.nrow; ++row)
			sparse_vec_add(result[row], contribution[row], F);
	}
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_physical_flatten_expanded(
	const sparse_mat<T, index_t>& coefficients, thread_pool* pool) {
	sparse_mat<T, index_t> result(1, coefficients.ncol * coefficients.nrow);
	for (size_t suffix = 0; suffix < coefficients.nrow; ++suffix)
		for (size_t pos = 0; pos < coefficients[suffix].nnz(); ++pos)
			result[0].push_back(static_cast<index_t>(
				static_cast<size_t>(coefficients[suffix](pos))
					* coefficients.nrow + suffix), coefficients[suffix][pos]);
	result.compress(pool);
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> joint_eet_physical_suffix_by_head(
	const sparse_mat<T, index_t>& flattened,
	size_t head_dimension, size_t suffix_dimension) {
	if (flattened.nrow != 1
		|| flattened.ncol != head_dimension * suffix_dimension)
		throw std::runtime_error("Expanded physical component has the wrong shape.");
	sparse_mat<T, index_t> result(suffix_dimension, head_dimension);
	for (size_t pos = 0; pos < flattened[0].nnz(); ++pos) {
		const size_t column = flattened[0](pos);
		result[column % suffix_dimension].push_back(
			static_cast<index_t>(column / suffix_dimension), flattened[0][pos]);
	}
	return result;
}

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> joint_eet_physical_last_k_tensor(
	const sparse_mat<T, index_t>& flattened, size_t head_dimension,
	size_t tail_weight, size_t letter_count, thread_pool* pool) {
	size_t suffix_dimension = 1;
	for (size_t i = 0; i < tail_weight; ++i)
		suffix_dimension *= letter_count;
	if (flattened.nrow != 1
		|| flattened.ncol != head_dimension * suffix_dimension)
		throw std::runtime_error("Physical last-k data have incompatible dimensions.");
	sparse_tensor<T, index_t, SPARSE_COO> tensor(flattened, pool);
	std::vector<size_t> dimensions{head_dimension};
	dimensions.insert(dimensions.end(), tail_weight, letter_count);
	tensor.reshape(dimensions);
	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(tensor), pool);
}
