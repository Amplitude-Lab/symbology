#pragma once

#include "bootstrap.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

template <typename T, typename index_t>
struct tail_state_t {
	size_t weight = 0;
	size_t r_count = 0;
	size_t letter_count = 0;
	sparse_mat<T, index_t> basis;
	sparse_mat<T, index_t> cyclic;
	sparse_mat<T, index_t> flip;
};

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_read_matrix(
	const std::filesystem::path& path, const field_t& F) {
	if (!std::filesystem::is_regular_file(path))
		throw std::runtime_error("Input file does not exist: " + path.string());
	auto matrix = sparse_mat_read_wxf<T, index_t>(path, F);
	if (matrix.nrow == 0 || matrix.ncol == 0)
		throw std::runtime_error("Input matrix is empty: " + path.string());
	return matrix;
}

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> nmhv_sew_read_tensor(
	const std::filesystem::path& path, const field_t& F, thread_pool* pool) {
	if (!std::filesystem::is_regular_file(path))
		throw std::runtime_error("Input file does not exist: " + path.string());
	auto tensor = sparse_tensor_read_wxf<T, index_t>(path, F, pool);
	if (tensor.rank() == 0)
		throw std::runtime_error("Input tensor is empty: " + path.string());
	return tensor;
}

template <typename T, typename index_t, SPARSE_TYPE Type>
size_t nmhv_sew_validate_condition_layout(
	const sparse_tensor<T, index_t, Type>& condition,
	size_t letter_count) {
	if (condition.rank() != 3
		|| condition.dim(0) != letter_count
		|| condition.dim(1) != letter_count) {
		throw std::runtime_error(
			"Condition tensor must have layout {letter,letter,condition}.");
	}
	return condition.dim(2);
}

inline void nmhv_sew_sync_file(const std::filesystem::path& path) {
#ifdef _WIN32
	const int descriptor = _wopen(path.c_str(), _O_RDWR | _O_BINARY);
#else
	const int descriptor = ::open(path.c_str(), O_RDONLY);
#endif
	if (descriptor == -1)
		throw std::runtime_error("Cannot open file for sync: " + path.string());
#ifdef _WIN32
	const int sync_result = _commit(descriptor);
	const int close_result = _close(descriptor);
#else
	const int sync_result = ::fsync(descriptor);
	const int close_result = ::close(descriptor);
#endif
	if (sync_result != 0 || close_result != 0)
		throw std::runtime_error("Cannot sync file: " + path.string());
}

inline std::filesystem::path nmhv_sew_create_staging_directory(
	const std::filesystem::path& requested_parent) {
	const auto parent = requested_parent.empty()
		? std::filesystem::path(".") : requested_parent;
	std::error_code error;
	std::filesystem::create_directories(parent, error);
	if (error || !std::filesystem::is_directory(parent))
		throw std::runtime_error(
			"Cannot create output parent directory: " + parent.string());

	const auto nonce = std::chrono::steady_clock::now()
		.time_since_epoch().count();
	for (size_t attempt = 0; attempt < 1000; ++attempt) {
		const auto staging = parent /
			(".nmhv_sew_stage_" + std::to_string(nonce)
				+ "_" + std::to_string(attempt));
		error.clear();
		if (std::filesystem::create_directory(staging, error))
			return staging;
		if (error
			&& error != std::make_error_code(std::errc::file_exists)) {
			throw std::runtime_error(
				"Cannot reserve output staging directory: " + error.message());
		}
	}
	throw std::runtime_error("Cannot reserve a unique output staging directory.");
}

template <typename bytes_t>
void nmhv_sew_publish_bytes(
	const std::filesystem::path& path, const bytes_t& bytes) {
	if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
		throw std::runtime_error("Output file already exists: " + path.string());
	if (bytes.size()
		> static_cast<size_t>(std::numeric_limits<std::streamsize>::max())) {
		throw std::overflow_error("Serialized WXF exceeds streamsize.");
	}

	const auto staging_directory =
		nmhv_sew_create_staging_directory(path.parent_path());
	const auto staged_path = staging_directory / "payload";
	{
		std::ofstream output(
			staged_path, std::ios::binary | std::ios::trunc);
		if (!output)
			throw std::runtime_error(
				"Cannot write staged file: " + staged_path.string());
		output.write(
			reinterpret_cast<const char*>(bytes.data()),
			static_cast<std::streamsize>(bytes.size()));
		output.flush();
		if (!output)
			throw std::runtime_error(
				"Cannot write staged file: " + staged_path.string());
	}
	 nmhv_sew_sync_file(staged_path);

	std::error_code error;
	std::filesystem::create_hard_link(staged_path, path, error);
	if (error) {
		if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
			throw std::runtime_error(
				"Output file already exists: " + path.string());
		throw std::runtime_error(
			"Cannot atomically publish output file " + path.string()
				+ ": " + error.message());
	}
	error.clear();
	std::filesystem::remove(staged_path, error);
	if (error)
		throw std::runtime_error(
			"Published output but cannot remove staged payload: "
				+ error.message());
	error.clear();
	std::filesystem::remove(staging_directory, error);
	if (error)
		throw std::runtime_error(
			"Published output but cannot remove staging directory: "
				+ error.message());
}

template <typename T, typename index_t>
void nmhv_sew_write_matrix(
	const std::filesystem::path& path,
	const sparse_mat<T, index_t>& matrix) {
	auto bytes = sparse_mat_write_wxf(matrix);
	nmhv_sew_publish_bytes(path, bytes);
}

template <typename T, typename index_t, SPARSE_TYPE Type>
void nmhv_sew_write_tensor(
	const std::filesystem::path& path,
	const sparse_tensor<T, index_t, Type>& tensor) {
	auto bytes = sparse_tensor_write_wxf(tensor);
	nmhv_sew_publish_bytes(path, bytes);
}

inline size_t nmhv_sew_checked_size_product(
	size_t left, size_t right, const std::string& label) {
	if (left != 0 && right > std::numeric_limits<size_t>::max() / left)
		throw std::overflow_error(label + " overflows size_t.");
	return left * right;
}

inline size_t nmhv_sew_checked_power(
	size_t base, size_t exponent, const std::string& label) {
	size_t result = 1;
	for (size_t i = 0; i < exponent; ++i)
		result = nmhv_sew_checked_size_product(result, base, label);
	return result;
}

template <typename index_t>
size_t nmhv_sew_checked_index_dimension(
	size_t value, const std::string& label) {
	static_assert(std::is_integral_v<index_t>);
	const auto limit = static_cast<size_t>(std::numeric_limits<index_t>::max());
	if (value > limit)
		throw std::overflow_error(
			label + " exceeds the sparse index limit (INT32_MAX for int32_t).");
	return value;
}

template <typename index_t>
size_t nmhv_sew_checked_index_product(
	std::initializer_list<size_t> factors, const std::string& label) {
	const auto limit = static_cast<size_t>(std::numeric_limits<index_t>::max());
	size_t result = 1;
	for (const size_t factor : factors) {
		nmhv_sew_checked_index_dimension<index_t>(factor, label);
		if (result != 0 && factor > limit / result)
			throw std::overflow_error(
				label + " exceeds the sparse index limit (INT32_MAX for int32_t).");
		result *= factor;
	}
	return result;
}

template <typename index_t>
size_t nmhv_sew_validate_tail_shape(
	size_t weight, size_t r_count, size_t letter_count,
	size_t condition_count) {
	if (weight == 0 || r_count == 0 || letter_count == 0
		|| condition_count == 0)
		throw std::invalid_argument("Tail shape dimensions must be positive.");

	size_t word_count = 1;
	for (size_t position = 0; position < weight; ++position) {
		word_count = nmhv_sew_checked_index_product<index_t>(
			{word_count, letter_count}, "Tail word dimension");
	}
	const size_t ambient_columns = nmhv_sew_checked_index_product<index_t>(
		{r_count, word_count}, "Tail ambient column dimension");

	size_t rest_count = 1;
	for (size_t old_weight = 1; old_weight < weight; ++old_weight) {
		(void)nmhv_sew_checked_index_product<index_t>(
			{r_count, rest_count, condition_count},
			"Tail condition row dimension");
		rest_count = nmhv_sew_checked_index_product<index_t>(
			{rest_count, letter_count}, "Tail condition rest dimension");
	}
	return ambient_columns;
}

template <typename T, typename index_t>
void nmhv_sew_normalize_rows(
	sparse_mat<T, index_t>& matrix, thread_pool* pool) {
	static_assert(std::is_same_v<T, rat_t>,
		"Exact NMHV tail construction requires rat_t.");
	pool->detach_loop(0, matrix.nrow, [&](size_t row) {
		matrix[row].canonicalize();
		if (matrix[row].nnz() != 0)
			vec_cancel_divisor(matrix[row]);
	});
	pool->wait();
	matrix.sort_rows_by_nnz();
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_kernel_rows(
	sparse_mat<T, index_t>&& matrix, const field_t& F,
	rref_option_t& opt) {
	static_assert(std::is_same_v<T, rat_t>,
		"Exact NMHV tail construction requires rat_t.");
	matrix.compress(&(opt->pool));
	matrix.clear_zero_row();
	matrix.sort_rows_by_nnz();
	auto pivots = sparse_mat_rref_reconstruct(matrix, opt);
	auto kernel = sparse_mat_rref_kernel(matrix, pivots, F, opt).transpose();
	nmhv_sew_normalize_rows(kernel, &(opt->pool));
	return kernel;
}

template <typename T, typename index_t>
bool nmhv_sew_matrices_equal(
	sparse_mat<T, index_t> left, sparse_mat<T, index_t> right) {
	if (left.nrow != right.nrow || left.ncol != right.ncol)
		return false;
	left.compress();
	right.compress();
	for (size_t row = 0; row < left.nrow; ++row) {
		if (left[row].nnz() != right[row].nnz())
			return false;
		for (size_t pos = 0; pos < left[row].nnz(); ++pos) {
			if (left[row](pos) != right[row](pos)
				|| left[row][pos] != right[row][pos])
				return false;
		}
	}
	return true;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_transform_product_basis(
	const sparse_mat<T, index_t>& basis,
	const sparse_mat<T, index_t>& left_action,
	const sparse_mat<T, index_t>& right_action,
	const field_t& F, thread_pool* pool) {
	if (left_action.nrow != left_action.ncol
		|| right_action.nrow != right_action.ncol
		|| basis.ncol != nmhv_sew_checked_index_product<index_t>(
			{left_action.nrow, right_action.nrow},
			"Product-basis column dimension")) {
		throw std::runtime_error(
			"Product basis and action dimensions are incompatible.");
	}
	nmhv_sew_checked_index_dimension<index_t>(
		basis.nrow, "Product-basis row dimension");

	sparse_tensor<T, index_t, SPARSE_COO> tensor(basis, pool);
	tensor.reshape({basis.nrow, left_action.nrow, right_action.nrow});
	sparse_tensor<T, index_t, SPARSE_COO> left(left_action, pool);
	sparse_tensor<T, index_t, SPARSE_COO> right(right_action, pool);

	auto transformed = tensor_contract(tensor, left, 1, 0, F, pool);
	transformed.transpose_replace({0, 2, 1}, pool);
	transformed.sort_indices(pool);
	auto result = tensor_contract(transformed, right, 2, 0, F, pool);
	result.flatten({{0}, {1, 2}});
	return result.to_sparse_mat(pool);
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_select_columns(
	const sparse_mat<T, index_t>& matrix,
	const std::vector<size_t>& columns) {
	const index_t missing = index_sval<index_t>();
	nmhv_sew_checked_index_dimension<index_t>(
		matrix.nrow, "Selected matrix row dimension");
	nmhv_sew_checked_index_dimension<index_t>(
		matrix.ncol, "Selected matrix source column dimension");
	nmhv_sew_checked_index_dimension<index_t>(
		columns.size(), "Selected matrix target column dimension");
	std::vector<index_t> positions(matrix.ncol, missing);
	for (size_t i = 0; i < columns.size(); ++i) {
		if (columns[i] >= matrix.ncol)
			throw std::runtime_error("Selected column is out of range.");
		positions[columns[i]] = static_cast<index_t>(i);
	}

	sparse_mat<T, index_t> result(matrix.nrow, columns.size());
	for (size_t row = 0; row < matrix.nrow; ++row) {
		for (size_t pos = 0; pos < matrix[row].nnz(); ++pos) {
			const auto selected = positions[matrix[row](pos)];
			if (selected != missing)
				result[row].push_back(selected, matrix[row][pos]);
		}
	}
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_restrict_action(
	const sparse_mat<T, index_t>& basis,
	const sparse_mat<T, index_t>& transformed,
	const field_t& F, rref_option_t& opt) {
	static_assert(std::is_same_v<T, rat_t>,
		"Exact NMHV tail construction requires rat_t.");
	if (basis.nrow != transformed.nrow || basis.ncol != transformed.ncol)
		throw std::runtime_error("Basis and transformed basis dimensions differ.");
	nmhv_sew_checked_index_dimension<index_t>(
		basis.nrow, "Restricted-action dimension");
	nmhv_sew_checked_index_dimension<index_t>(
		basis.ncol, "Restricted-action ambient dimension");

	auto reduced = basis;
	reduced.sort_rows_by_nnz();
	auto pivot_groups = sparse_mat_rref_reconstruct(reduced, opt);
	std::vector<size_t> pivot_columns;
	for (const auto& group : pivot_groups) {
		for (const auto& pivot : group)
			pivot_columns.push_back(static_cast<size_t>(pivot.c));
	}
	std::ranges::sort(pivot_columns);
	pivot_columns.erase(
		std::unique(pivot_columns.begin(), pivot_columns.end()),
		pivot_columns.end());
	if (pivot_columns.size() != basis.nrow)
		throw std::runtime_error("Tail basis rows are not linearly independent.");

	auto pivot_block = nmhv_sew_select_columns(basis, pivot_columns);
	auto transformed_block = nmhv_sew_select_columns(transformed, pivot_columns);
	auto inverse = sparse_mat_inverse(pivot_block, F, opt);
	if (inverse.nrow != basis.nrow || inverse.ncol != basis.nrow)
		throw std::runtime_error("Tail basis pivot block is not invertible.");
	auto row_action = sparse_mat_mul(
		transformed_block, inverse, F, &(opt->pool));
	auto reconstructed = sparse_mat_mul(
		row_action, basis, F, &(opt->pool));
	if (!nmhv_sew_matrices_equal(reconstructed, transformed))
		throw std::runtime_error("Ambient action does not preserve the tail basis.");
	return row_action.transpose();
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_tail_condition_matrix(
	const tail_state_t<T, index_t>& tail,
	const sparse_tensor<T, index_t, SPARSE_COO>& condition,
	const field_t& F, thread_pool* pool) {
	if (tail.weight == 0)
		throw std::runtime_error("Tail weight must be positive.");
	const size_t condition_count = nmhv_sew_validate_condition_layout(
		condition, tail.letter_count);
	const size_t expected_ambient = nmhv_sew_validate_tail_shape<index_t>(
		tail.weight, tail.r_count, tail.letter_count, condition_count);
	if (tail.basis.ncol != expected_ambient)
		throw std::runtime_error("Tail basis ambient dimension is inconsistent.");
	nmhv_sew_checked_index_dimension<index_t>(
		tail.basis.nrow, "Tail basis row dimension");
	const size_t rest_count = nmhv_sew_checked_power(
		tail.letter_count, tail.weight - 1, "Tail condition rest dimension");
	(void)nmhv_sew_checked_index_product<index_t>(
		{tail.r_count, rest_count, condition_count},
		"Tail condition row dimension");
	(void)nmhv_sew_checked_index_product<index_t>(
		{tail.letter_count, tail.basis.nrow},
		"Tail condition variable dimension");

	std::vector<size_t> tail_dims{tail.basis.nrow, tail.r_count};
	tail_dims.insert(
		tail_dims.end(), tail.weight, tail.letter_count);
	sparse_tensor<T, index_t, SPARSE_COO> tail_tensor(tail.basis, pool);
	tail_tensor.reshape(tail_dims);
	auto constrained = tensor_contract(
		tail_tensor, condition, 2, 1, F, pool);

	std::vector<size_t> constraint_axes{1};
	for (size_t axis = 2; axis <= tail.weight; ++axis)
		constraint_axes.push_back(axis);
	constraint_axes.push_back(tail.weight + 2);
	constrained.flatten({
		constraint_axes,
		{tail.weight + 1, 0}
	});
	return constrained.to_sparse_mat(pool);
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_expand_tail(
	const sparse_mat<T, index_t>& recurrence,
	const tail_state_t<T, index_t>& old_tail,
	const field_t& F, thread_pool* pool) {
	const size_t recurrence_columns =
		nmhv_sew_checked_index_product<index_t>(
			{old_tail.letter_count, old_tail.basis.nrow},
			"Tail recurrence column dimension");
	if (recurrence.ncol != recurrence_columns)
		throw std::runtime_error("Tail recurrence has incompatible columns.");
	nmhv_sew_checked_index_dimension<index_t>(
		recurrence.nrow, "Expanded tail row dimension");
	const size_t next_ambient = nmhv_sew_validate_tail_shape<index_t>(
		old_tail.weight + 1, old_tail.r_count, old_tail.letter_count, 1);

	sparse_tensor<T, index_t, SPARSE_COO> recurrence_tensor(
		recurrence, pool);
	recurrence_tensor.reshape({
		recurrence.nrow, old_tail.letter_count, old_tail.basis.nrow
	});
	std::vector<size_t> old_dims{old_tail.basis.nrow, old_tail.r_count};
	old_dims.insert(
		old_dims.end(), old_tail.weight, old_tail.letter_count);
	sparse_tensor<T, index_t, SPARSE_COO> old_tensor(old_tail.basis, pool);
	old_tensor.reshape(old_dims);

	auto expanded = tensor_contract(
		recurrence_tensor, old_tensor, 2, 0, F, pool);
	std::vector<size_t> permutation{0, 2, 1};
	for (size_t axis = 3; axis < expanded.rank(); ++axis)
		permutation.push_back(axis);
	expanded.transpose_replace(permutation, pool);
	expanded.sort_indices(pool);
	std::vector<size_t> ambient_axes;
	for (size_t axis = 1; axis < expanded.rank(); ++axis)
		ambient_axes.push_back(axis);
	expanded.flatten({{0}, ambient_axes});
	auto result = expanded.to_sparse_mat(pool);
	if (result.ncol != next_ambient)
		throw std::runtime_error("Expanded tail ambient dimension is inconsistent.");
	return result;
}

template <typename T, typename index_t>
tail_state_t<T, index_t> build_tail_state(
	size_t weight,
	const sparse_mat<T, index_t>& qbar_basis,
	const sparse_tensor<T, index_t, SPARSE_COO>& condition,
	const sparse_mat<T, index_t>& cyclic_R,
	const sparse_mat<T, index_t>& cyclic_letter,
	const sparse_mat<T, index_t>& flip_R,
	const sparse_mat<T, index_t>& flip_letter,
	const field_t& F, rref_option_t& opt) {
	static_assert(std::is_same_v<T, rat_t>,
		"Exact NMHV tail construction requires rat_t.");
	if (weight == 0)
		throw std::runtime_error("Tail weight must be positive.");
	if (cyclic_R.nrow == 0 || cyclic_R.nrow != cyclic_R.ncol
		|| flip_R.nrow != cyclic_R.nrow || flip_R.ncol != cyclic_R.ncol)
		throw std::runtime_error("R actions must be equally sized square matrices.");
	if (cyclic_letter.nrow == 0
		|| cyclic_letter.nrow != cyclic_letter.ncol
		|| flip_letter.nrow != cyclic_letter.nrow
		|| flip_letter.ncol != cyclic_letter.ncol)
		throw std::runtime_error(
			"Letter actions must be equally sized square matrices.");
	const size_t condition_count = nmhv_sew_validate_condition_layout(
		condition, cyclic_letter.nrow);
	const size_t target_ambient = nmhv_sew_validate_tail_shape<index_t>(
		weight, cyclic_R.nrow, cyclic_letter.nrow, condition_count);
	nmhv_sew_checked_index_dimension<index_t>(
		qbar_basis.nrow, "Qbar basis row dimension");
	const size_t qbar_columns = nmhv_sew_checked_index_product<index_t>(
		{cyclic_R.nrow, cyclic_letter.nrow}, "Qbar basis column dimension");
	if (qbar_basis.ncol != qbar_columns)
		throw std::runtime_error("Qbar basis must have layout {tail,R*letter}.");

	tail_state_t<T, index_t> tail;
	tail.weight = 1;
	tail.r_count = cyclic_R.nrow;
	tail.letter_count = cyclic_letter.nrow;
	tail.basis = qbar_basis;
	nmhv_sew_normalize_rows(tail.basis, &(opt->pool));
	auto cyclic_ambient = nmhv_sew_transform_product_basis(
		tail.basis, cyclic_R, cyclic_letter, F, &(opt->pool));
	tail.cyclic = nmhv_sew_restrict_action(
		tail.basis, cyclic_ambient, F, opt);
	auto flip_ambient = nmhv_sew_transform_product_basis(
		tail.basis, flip_R, flip_letter, F, &(opt->pool));
	tail.flip = nmhv_sew_restrict_action(
		tail.basis, flip_ambient, F, opt);

	while (tail.weight < weight) {
		auto condition_matrix = nmhv_sew_tail_condition_matrix(
			tail, condition, F, &(opt->pool));
		auto recurrence = nmhv_sew_kernel_rows(
			std::move(condition_matrix), F, opt);
		auto next_basis = nmhv_sew_expand_tail(
			recurrence, tail, F, &(opt->pool));

		auto cyclic_column = tail.cyclic;
		auto cyclic_row = cyclic_column.transpose();
		auto cyclic_recurrence = nmhv_sew_transform_product_basis(
			recurrence, cyclic_letter, cyclic_row, F, &(opt->pool));
		auto next_cyclic = nmhv_sew_restrict_action(
			recurrence, cyclic_recurrence, F, opt);
		auto flip_column = tail.flip;
		auto flip_row = flip_column.transpose();
		auto flip_recurrence = nmhv_sew_transform_product_basis(
			recurrence, flip_letter, flip_row, F, &(opt->pool));
		auto next_flip = nmhv_sew_restrict_action(
			recurrence, flip_recurrence, F, opt);

		tail.weight += 1;
		tail.basis = std::move(next_basis);
		tail.cyclic = std::move(next_cyclic);
		tail.flip = std::move(next_flip);
		const size_t expected_columns = nmhv_sew_validate_tail_shape<index_t>(
			tail.weight, tail.r_count, tail.letter_count, condition_count);
		if (tail.basis.ncol != expected_columns)
			throw std::runtime_error("Tail ambient basis has an invalid dimension.");
		std::cout << "Tail weight " << tail.weight
			<< " dimension: " << tail.basis.nrow << std::endl;
	}
	if (tail.basis.ncol != target_ambient)
		throw std::runtime_error("Final tail ambient dimension is inconsistent.");
	return tail;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_column_eigenspace(
	const sparse_mat<T, index_t>& action, const T& eigenvalue,
	const field_t& F, rref_option_t& opt) {
	if (action.nrow != action.ncol)
		throw std::runtime_error("Action matrix must be square.");
	nmhv_sew_checked_index_dimension<index_t>(
		action.nrow, "Action eigenspace dimension");
	auto constraint = action;
	for (size_t row = 0; row < constraint.nrow; ++row) {
		auto* diagonal = constraint[row].find(
			static_cast<index_t>(row), false);
		if (diagonal == nullptr)
			constraint[row].push_back(
				static_cast<index_t>(row), -eigenvalue);
		else
			*diagonal = *diagonal - eigenvalue;
	}
	constraint.compress(&(opt->pool));
	return nmhv_sew_kernel_rows(std::move(constraint), F, opt);
}

template <typename T, typename index_t>
std::pair<sparse_mat<T, index_t>, sparse_mat<T, index_t>>
tail_parity_sectors(
	const tail_state_t<T, index_t>& tail,
	const field_t& F, rref_option_t& opt) {
	if (tail.cyclic.nrow != tail.basis.nrow
		|| tail.cyclic.ncol != tail.basis.nrow)
		throw std::runtime_error("Tail cyclic action has incompatible dimensions.");
	auto cyclic2 = sparse_mat_mul(
		tail.cyclic, tail.cyclic, F, &(opt->pool));
	auto parity = sparse_mat_mul(
		cyclic2, tail.cyclic, F, &(opt->pool));
	auto even = nmhv_sew_column_eigenspace(parity, T(1), F, opt);
	auto odd = nmhv_sew_column_eigenspace(parity, T(-1), F, opt);
	if (even.nrow + odd.nrow != tail.basis.nrow)
		throw std::runtime_error("Tail parity sectors do not span the tail basis.");
	return {std::move(even), std::move(odd)};
}

struct d3_module_count_t {
	size_t trivial = 0;
	size_t sign = 0;
	size_t standard = 0;
	size_t generators = 0;
};

template <typename T, typename index_t>
T nmhv_sew_matrix_trace(const sparse_mat<T, index_t>& matrix) {
	if (matrix.nrow != matrix.ncol)
		throw std::runtime_error("Cannot take the trace of a nonsquare matrix.");
	T trace = 0;
	for (size_t row = 0; row < matrix.nrow; ++row) {
		const auto* diagonal = matrix[row].find(static_cast<index_t>(row), false);
		if (diagonal != nullptr)
			trace += *diagonal;
	}
	return trace;
}

inline int64_t nmhv_sew_exact_integer(
	const rat_t& value, const std::string& label) {
	if (!value.is_integer() || !value.num().fits_si())
		throw std::runtime_error(label + " is not a machine-size integer.");
	return static_cast<int64_t>(value.num().to_si());
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_restrict_coordinate_action(
	const sparse_mat<T, index_t>& basis,
	const sparse_mat<T, index_t>& action,
	const field_t& F, rref_option_t& opt) {
	if (action.nrow != action.ncol || basis.ncol != action.nrow)
		throw std::runtime_error("Coordinate action and basis dimensions differ.");
	auto action_copy = action;
	auto row_action = action_copy.transpose();
	auto transformed = sparse_mat_mul(
		basis, row_action, F, &(opt->pool));
	return nmhv_sew_restrict_action(basis, transformed, F, opt);
}

template <typename T, typename index_t>
d3_module_count_t nmhv_sew_d3_module_count(
	const sparse_mat<T, index_t>& cyclic,
	const sparse_mat<T, index_t>& flip, int parity_sign) {
	static_assert(std::is_same_v<T, rat_t>,
		"D3 character decomposition is exact over QQ.");
	if ((parity_sign != 1 && parity_sign != -1)
		|| cyclic.nrow != cyclic.ncol
		|| flip.nrow != cyclic.nrow || flip.ncol != cyclic.ncol)
		throw std::runtime_error("Invalid D3 sector actions.");

	const int64_t dimension = static_cast<int64_t>(cyclic.nrow);
	const int64_t trace_cyclic = parity_sign
		* nmhv_sew_exact_integer(nmhv_sew_matrix_trace(cyclic), "Cyclic trace");
	const int64_t trace_flip =
		nmhv_sew_exact_integer(nmhv_sew_matrix_trace(flip), "Flip trace");
	const int64_t trivial_numerator =
		dimension + 2 * trace_cyclic + 3 * trace_flip;
	const int64_t sign_numerator =
		dimension + 2 * trace_cyclic - 3 * trace_flip;
	const int64_t standard_numerator = dimension - trace_cyclic;
	if (trivial_numerator < 0 || sign_numerator < 0
		|| standard_numerator < 0 || trivial_numerator % 6 != 0
		|| sign_numerator % 6 != 0 || standard_numerator % 3 != 0) {
		throw std::runtime_error("D3 character multiplicities are not nonnegative integers.");
	}

	d3_module_count_t result;
	result.trivial = static_cast<size_t>(trivial_numerator / 6);
	result.sign = static_cast<size_t>(sign_numerator / 6);
	result.standard = static_cast<size_t>(standard_numerator / 3);
	result.generators = std::max({
		result.trivial, result.sign, (result.standard + 1) / 2
	});
	return result;
}

template <typename T, typename index_t>
std::array<sparse_mat<T, index_t>, 6> nmhv_sew_action_orbits(
	const sparse_mat<T, index_t>& vectors,
	const sparse_mat<T, index_t>& cyclic,
	const sparse_mat<T, index_t>& flip,
	const field_t& F, thread_pool* pool) {
	if (cyclic.nrow != cyclic.ncol || flip.nrow != cyclic.nrow
		|| flip.ncol != cyclic.ncol || vectors.ncol != cyclic.nrow)
		throw std::runtime_error("D3 action and vector dimensions differ.");
	auto cyclic_copy = cyclic;
	auto flip_copy = flip;
	auto cyclic_row = cyclic_copy.transpose();
	auto flip_row = flip_copy.transpose();
	std::array<sparse_mat<T, index_t>, 6> orbit;
	orbit[0] = vectors;
	orbit[1] = sparse_mat_mul(orbit[0], cyclic_row, F, pool);
	orbit[2] = sparse_mat_mul(orbit[1], cyclic_row, F, pool);
	for (size_t power = 0; power < 3; ++power)
		orbit[power + 3] = sparse_mat_mul(orbit[power], flip_row, F, pool);
	return orbit;
}

template <typename T, typename index_t>
size_t nmhv_sew_row_rank(
	sparse_mat<T, index_t> matrix, rref_option_t& opt) {
	matrix.compress(&(opt->pool));
	matrix.clear_zero_row();
	if (matrix.nrow == 0)
		return 0;
	matrix.sort_rows_by_nnz();
	auto pivots = sparse_mat_rref_reconstruct(matrix, opt);
	size_t rank = 0;
	for (const auto& group : pivots)
		rank += group.size();
	return rank;
}

template <typename T, typename index_t>
bool nmhv_sew_same_row_space(
	sparse_mat<T, index_t> left,
	sparse_mat<T, index_t> right, rref_option_t& opt) {
	if (left.ncol != right.ncol)
		return false;
	const size_t left_rank = nmhv_sew_row_rank(left, opt);
	const size_t right_rank = nmhv_sew_row_rank(right, opt);
	left.append(std::move(right));
	return left_rank == right_rank
		&& nmhv_sew_row_rank(std::move(left), opt) == left_rank;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_module_orbits(
	const sparse_mat<T, index_t>& generators,
	const sparse_mat<T, index_t>& cyclic,
	const sparse_mat<T, index_t>& flip,
	const field_t& F, thread_pool* pool) {
	auto orbit = nmhv_sew_action_orbits(
		generators, cyclic, flip, F, pool);
	auto result = std::move(orbit[0]);
	for (size_t image = 1; image < orbit.size(); ++image)
		result.append(std::move(orbit[image]));
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_module_generators(
	const sparse_mat<T, index_t>& sector_basis,
	const sparse_mat<T, index_t>& full_cyclic,
	const sparse_mat<T, index_t>& full_flip,
	int parity_sign, const field_t& F, rref_option_t& opt,
	uint64_t seed = 20260904) {
	static_assert(std::is_same_v<T, rat_t>,
		"Module generators are selected and verified over QQ.");
	if (sector_basis.nrow == 0)
		return sparse_mat<T, index_t>(0, sector_basis.ncol);
	auto cyclic = nmhv_sew_restrict_coordinate_action(
		sector_basis, full_cyclic, F, opt);
	auto flip = nmhv_sew_restrict_coordinate_action(
		sector_basis, full_flip, F, opt);
	const auto count = nmhv_sew_d3_module_count(
		cyclic, flip, parity_sign).generators;
	if (count == 0)
		return sparse_mat<T, index_t>(0, sector_basis.ncol);

	std::mt19937_64 generator(seed ^ static_cast<uint64_t>(parity_sign + 2));
	std::vector<size_t> positions(sector_basis.nrow);
	std::iota(positions.begin(), positions.end(), 0);
	for (size_t attempt = 0; attempt < 192; ++attempt) {
		const size_t support = std::min(
			sector_basis.nrow,
			attempt < 64 ? size_t(8) : attempt < 128 ? size_t(16) : sector_basis.nrow);
		sparse_mat<T, index_t> selected(count, sector_basis.nrow);
		for (size_t row = 0; row < count; ++row) {
			std::shuffle(positions.begin(), positions.end(), generator);
			for (size_t pos = 0; pos < support; ++pos) {
				int value = static_cast<int>(generator() % 5) - 2;
				if (value == 0)
					value = 1;
				selected[row].push_back(
					static_cast<index_t>(positions[pos]), T(value));
			}
			selected[row].sort_indices();
		}
		auto orbit = nmhv_sew_module_orbits(
			selected, cyclic, flip, F, &(opt->pool));
		if (nmhv_sew_row_rank(std::move(orbit), opt) == sector_basis.nrow) {
			auto result = sparse_mat_mul(
				selected, sector_basis, F, &(opt->pool));
			if (nmhv_sew_row_rank(
				nmhv_sew_module_orbits(
					result, full_cyclic, full_flip, F, &(opt->pool)), opt)
				!= sector_basis.nrow) {
				throw std::runtime_error("Selected D3 generators failed exact orbit verification.");
			}
			return result;
		}
	}
	throw std::runtime_error("Could not find the character-minimal D3 module generators.");
}

template <typename T, typename index_t>
sparse_vec<T, index_t> nmhv_sew_tensor_product_row(
	const sparse_vec<T, index_t>& left,
	const sparse_vec<T, index_t>& right,
	size_t right_dimension) {
	sparse_vec<T, index_t> result;
	result.reserve(nmhv_sew_checked_size_product(
		left.nnz(), right.nnz(), "Tensor-product row nnz"));
	for (size_t i = 0; i < left.nnz(); ++i) {
		for (size_t j = 0; j < right.nnz(); ++j) {
			const size_t column = nmhv_sew_checked_size_product(
				static_cast<size_t>(left(i)), right_dimension,
				"Tensor-product column") + static_cast<size_t>(right(j));
			nmhv_sew_checked_index_dimension<index_t>(
				column, "Tensor-product column");
			result.push_back(
				static_cast<index_t>(column), left[i] * right[j]);
		}
	}
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_group_average_candidates(
	const sparse_mat<T, index_t>& head_basis,
	const sparse_mat<T, index_t>& tail_generators,
	const sparse_mat<T, index_t>& head_cyclic,
	const sparse_mat<T, index_t>& head_flip,
	const sparse_mat<T, index_t>& tail_cyclic,
	const sparse_mat<T, index_t>& tail_flip,
	const field_t& F, thread_pool* pool) {
	const size_t columns = nmhv_sew_checked_index_product<index_t>(
		{head_basis.ncol, tail_generators.ncol},
		"Natural sewing column dimension");
	if (head_basis.nrow == 0 || tail_generators.nrow == 0)
		return sparse_mat<T, index_t>(0, columns);
	auto head_orbits = nmhv_sew_action_orbits(
		head_basis, head_cyclic, head_flip, F, pool);
	auto tail_orbits = nmhv_sew_action_orbits(
		tail_generators, tail_cyclic, tail_flip, F, pool);
	const size_t rows = nmhv_sew_checked_index_product<index_t>(
		{head_basis.nrow, tail_generators.nrow},
		"Invariant candidate row dimension");
	sparse_mat<T, index_t> result(rows, columns);
	pool->detach_loop(0, rows, [&](size_t row) {
		const size_t head = row / tail_generators.nrow;
		const size_t tail = row % tail_generators.nrow;
		for (size_t image = 0; image < 6; ++image) {
			auto term = nmhv_sew_tensor_product_row(
				head_orbits[image][head], tail_orbits[image][tail],
				tail_generators.ncol);
			sparse_vec_add_mul(result[row], term, T(1), F);
		}
		result[row].compress();
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_matched_parity_candidates(
	const sparse_mat<T, index_t>& head_even,
	const sparse_mat<T, index_t>& head_odd,
	const sparse_mat<T, index_t>& tail_even_generators,
	const sparse_mat<T, index_t>& tail_odd_generators,
	const sparse_mat<T, index_t>& head_cyclic,
	const sparse_mat<T, index_t>& head_flip,
	const sparse_mat<T, index_t>& tail_cyclic,
	const sparse_mat<T, index_t>& tail_flip,
	const field_t& F, thread_pool* pool) {
	auto even = nmhv_sew_group_average_candidates(
		head_even, tail_even_generators, head_cyclic, head_flip,
		tail_cyclic, tail_flip, F, pool);
	auto odd = nmhv_sew_group_average_candidates(
		head_odd, tail_odd_generators, head_cyclic, head_flip,
		tail_cyclic, tail_flip, F, pool);
	if (even.nrow == 0)
		return odd;
	if (odd.nrow != 0)
		even.append(std::move(odd));
	return even;
}

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> nmhv_sew_tail_basis_tensor(
	const tail_state_t<T, index_t>& tail, thread_pool* pool) {
	const size_t expected_columns = nmhv_sew_validate_tail_shape<index_t>(
		tail.weight, tail.r_count, tail.letter_count, 1);
	if (tail.basis.ncol != expected_columns)
		throw std::runtime_error("Tail basis ambient dimension is inconsistent.");
	nmhv_sew_checked_index_dimension<index_t>(
		tail.basis.nrow, "Tail tensor row dimension");
	std::vector<size_t> dimensions{tail.basis.nrow, tail.r_count};
	dimensions.insert(dimensions.end(), tail.weight, tail.letter_count);
	sparse_tensor<T, index_t, SPARSE_COO> tensor(tail.basis, pool);
	tensor.reshape(dimensions);
	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(tensor), pool);
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_exact_row_basis(
	sparse_mat<T, index_t> matrix, rref_option_t& opt) {
	static_assert(std::is_same_v<T, rat_t>,
		"Exact row-space normalization requires rat_t.");
	matrix.compress(&(opt->pool));
	matrix.clear_zero_row();
	matrix.sort_rows_by_nnz();
	sparse_mat_rref_reconstruct(matrix, opt);
	matrix.clear_zero_row();
	nmhv_sew_normalize_rows(matrix, &(opt->pool));
	std::ranges::sort(matrix.rows, [](const auto& left, const auto& right) {
		if (left.nnz() == 0)
			return false;
		if (right.nnz() == 0)
			return true;
		return left(0) < right(0);
	});
	return matrix;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_fec_expansion(
	size_t head_count, size_t letter_count,
	const std::vector<sparse_tensor<T, index_t, SPARSE_COO>>& extensions,
	const field_t& F, thread_pool* pool) {
	nmhv_sew_checked_index_dimension<index_t>(
		head_count, "FEC expansion head dimension");
	sparse_mat<T, index_t> identity(head_count, head_count);
	for (size_t row = 0; row < head_count; ++row)
		identity[row].push_back(static_cast<index_t>(row), T(1));
	sparse_tensor<T, index_t, SPARSE_COO> expansion(identity, pool);
	expansion.reshape({head_count, head_count});
	size_t current_count = head_count;
	for (const auto& FEC : extensions) {
		if (FEC.rank() != 3 || FEC.dim(1) != current_count
			|| FEC.dim(2) != letter_count) {
			throw std::runtime_error(
				"FEC extension has incompatible {new,old,letter} dimensions.");
		}
		auto next = tensor_contract(FEC, expansion, 1, 0, F, pool);
		std::vector<size_t> permutation{0, 2};
		for (size_t axis = 3; axis < next.rank(); ++axis)
			permutation.push_back(axis);
		permutation.push_back(1);
		next.transpose_replace(permutation, pool);
		next.sort_indices(pool);
		expansion = std::move(next);
		current_count = FEC.dim(0);
	}
	std::vector<size_t> columns;
	for (size_t axis = 1; axis < expansion.rank(); ++axis)
		columns.push_back(axis);
	expansion.flatten({{0}, columns});
	return expansion.to_sparse_mat(pool);
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_natural_ambient_blocks(
	const sparse_mat<T, index_t>& natural,
	const tail_state_t<T, index_t>& tail,
	size_t head_count, const field_t& F, thread_pool* pool) {
	static_assert(std::is_same_v<T, rat_t>,
		"Direct natural-ambient expansion is exact over QQ.");
	if (natural.ncol != nmhv_sew_checked_index_product<index_t>(
		{head_count, tail.basis.nrow}, "Natural sewing columns"))
		throw std::runtime_error("Natural sewing basis has incompatible columns.");
	(void)F;
	const size_t prefix_count = nmhv_sew_checked_power(
		tail.letter_count, tail.weight - 1,
		"Natural ambient tail-prefix dimension");
	const size_t ambient_columns = nmhv_sew_checked_index_product<index_t>(
		{head_count, prefix_count}, "Natural ambient columns");
	const size_t block_rows = nmhv_sew_checked_index_product<index_t>(
		{natural.nrow, tail.r_count, tail.letter_count},
		"Natural ambient block rows");
	sparse_mat<T, index_t> result(block_rows, ambient_columns);
	pool->detach_loop(0, natural.nrow, [&](size_t solution) {
		for (size_t pos = 0; pos < natural[solution].nnz(); ++pos) {
			const size_t natural_column = natural[solution](pos);
			const size_t head = natural_column / tail.basis.nrow;
			const size_t tail_index = natural_column % tail.basis.nrow;
			for (size_t tail_pos = 0;
				tail_pos < tail.basis[tail_index].nnz(); ++tail_pos) {
				const size_t tail_column = tail.basis[tail_index](tail_pos);
				const size_t final_letter = tail_column % tail.letter_count;
				const size_t prefix_and_r = tail_column / tail.letter_count;
				const size_t prefix = prefix_and_r % prefix_count;
				const size_t r = prefix_and_r / prefix_count;
				const size_t row =
					(solution * tail.r_count + r) * tail.letter_count
					+ final_letter;
				const size_t column = head * prefix_count + prefix;
				result[row].push_back(
					static_cast<index_t>(column),
					natural[solution][pos] * tail.basis[tail_index][tail_pos]);
			}
		}
	});
	pool->wait();
	pool->detach_loop(0, result.nrow, [&](size_t row_index) {
		auto& row = result[row_index];
		row.sort_indices();
		sparse_vec<T, index_t> combined;
		combined.reserve(row.nnz());
		for (size_t pos = 0; pos < row.nnz();) {
			const index_t column = row(pos);
			T value = row[pos];
			++pos;
			while (pos < row.nnz() && row(pos) == column) {
				value += row[pos];
				++pos;
			}
			if (value != 0)
				combined.push_back(column, std::move(value));
		}
		row = std::move(combined);
	});
	pool->wait();
	return result;
}

template <typename T, typename index_t>
sparse_mat<T, index_t> nmhv_sew_project_standard(
	const sparse_mat<T, index_t>& natural,
	const tail_state_t<T, index_t>& tail,
	const sparse_mat<T, index_t>& expansion,
	size_t head_count, const field_t& F, rref_option_t& opt) {
	static_assert(std::is_same_v<T, rat_t>,
		"Standard sewing projection is verified over QQ.");
	if (expansion.nrow == 0)
		throw std::runtime_error("FEC expansion is empty.");
	auto reduced = expansion;
	reduced.sort_rows_by_nnz();
	auto pivot_groups = sparse_mat_rref_reconstruct(reduced, opt);
	std::vector<size_t> pivot_columns;
	for (const auto& group : pivot_groups) {
		for (const auto& pivot : group)
			pivot_columns.push_back(static_cast<size_t>(pivot.c));
	}
	std::ranges::sort(pivot_columns);
	pivot_columns.erase(
		std::unique(pivot_columns.begin(), pivot_columns.end()),
		pivot_columns.end());
	if (pivot_columns.size() != expansion.nrow)
		throw std::runtime_error("FEC expansion rows are not independent.");
	auto pivot_block = nmhv_sew_select_columns(expansion, pivot_columns);
	auto inverse = sparse_mat_inverse(pivot_block, F, opt);
	sparse_mat<T, index_t> coordinates(0, expansion.nrow);
	for (size_t solution = 0; solution < natural.nrow; ++solution) {
		auto natural_block = natural.take(
			{solution, solution + 1}, &(opt->pool));
		auto blocks = nmhv_sew_natural_ambient_blocks(
			natural_block, tail, head_count, F, &(opt->pool));
		if (blocks.ncol != expansion.ncol)
			throw std::runtime_error(
				"Tail ambient blocks and FEC expansion differ.");
		auto selected = nmhv_sew_select_columns(blocks, pivot_columns);
		auto coordinate_block = sparse_mat_mul(
			selected, inverse, F, &(opt->pool));
		auto reconstructed = sparse_mat_mul(
			coordinate_block, expansion, F, &(opt->pool));
		if (!nmhv_sew_matrices_equal(
			std::move(reconstructed), std::move(blocks))) {
			throw std::runtime_error(
				"Sewing result is not in the target FEC span.");
		}
		coordinates.append(std::move(coordinate_block));
	}

	const size_t expected_blocks = nmhv_sew_checked_index_product<index_t>(
		{natural.nrow, tail.r_count, tail.letter_count},
		"Standard projection block rows");
	if (coordinates.nrow != expected_blocks)
		throw std::runtime_error("Standard projection has an invalid block count.");
	const size_t standard_columns = nmhv_sew_checked_index_product<index_t>(
		{tail.r_count, expansion.nrow, tail.letter_count},
		"Standard amplitude columns");
	sparse_mat<T, index_t> standard(natural.nrow, standard_columns);
	for (size_t solution = 0; solution < natural.nrow; ++solution) {
		for (size_t r = 0; r < tail.r_count; ++r) {
			for (size_t letter = 0; letter < tail.letter_count; ++letter) {
				const size_t block =
					(solution * tail.r_count + r) * tail.letter_count + letter;
				for (size_t pos = 0; pos < coordinates[block].nnz(); ++pos) {
					const size_t target = coordinates[block](pos);
					const size_t column =
						(r * expansion.nrow + target) * tail.letter_count + letter;
					standard[solution].push_back(
						static_cast<index_t>(column), coordinates[block][pos]);
				}
			}
		}
		standard[solution].sort_indices();
	}
	return nmhv_sew_exact_row_basis(std::move(standard), opt);
}
