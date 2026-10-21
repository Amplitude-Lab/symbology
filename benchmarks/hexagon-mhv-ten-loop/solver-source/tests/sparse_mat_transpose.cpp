#include "SparseRREF/sparse_type.h"
#include <iostream>
#include <memory>
#include <random>

using namespace SparseRREF;
using Flint::int_t;
using Flint::rat_t;

void check(bool ok, const char* message) {
	if (!ok) throw std::runtime_error(message);
}

template <typename T, Flint::builtin_integral index_t>
void compare(sparse_mat<T, index_t>& matrix) {
	auto expected = matrix.transpose();
	matrix.transpose_replace();
	auto& result = matrix;
	check(result.nrow == expected.nrow && result.ncol == expected.ncol, "Wrong dimensions");
	for (size_t i = 0; i < result.nrow; i++) {
		check(result[i] == expected[i], "Values/indices differ from copying transpose");
		check(result[i].alloc() == result[i].nnz(), "Target has spare capacity");
	}
	matrix.transpose_replace();
	check(matrix.nrow == expected.ncol && matrix.ncol == expected.nrow, "Wrong roundtrip dimensions");
}

template <typename T, Flint::builtin_integral index_t>
void cases() {
	for (auto [r, c] : {std::pair<size_t, size_t>{0, 0}, {0, 7}, {7, 0}, {6, 9}, {1, 1}}) {
		sparse_mat<T, index_t> matrix(r, c);
		compare(matrix);
	}
	std::mt19937 random(20260921);
	for (size_t trial = 0; trial < 100; trial++) {
		sparse_mat<T, index_t> matrix(1 + random() % 15, 1 + random() % 19);
		for (size_t i = 0; i < matrix.nrow; i++) {
			for (size_t j = 0, n = random() % 30; j < n; j++)
				matrix[i].push_back(static_cast<index_t>(random() % matrix.ncol), T(random() % 5));
		}
		compare(matrix); // Includes unsorted/duplicate indices and stored zeros.
	}
}

struct move_only {
	std::unique_ptr<int> value;
	inline static size_t moves = 0;
	move_only() = default;
	explicit move_only(int n) : value(std::make_unique<int>(n)) {}
	move_only(const move_only&) = delete;
	move_only& operator=(const move_only&) = delete;
	move_only(move_only&&) = default;
	move_only& operator=(move_only&& other) noexcept {
		value = std::move(other.value);
		moves++;
		return *this;
	}
};

int main() try {
	Flint::set_memory_functions();
	cases<rat_t, int32_t>();
	cases<int_t, int16_t>();
	cases<uint64_t, int8_t>();
	cases<uint64_t, uint16_t>();

	sparse_mat<rat_t, int32_t> large(3, 5);
	large[0].push_back(4, rat_t("123456789012345678901234567890123456789012345678901/987654321098765432109876543210987654321"));
	large[2].push_back(1, rat_t("-987654321098765432109876543210987654321098765432109/12345678901234567890123456789"));
	for (size_t i = 0; i < 100; i++) compare(large);

	sparse_mat<bool, int16_t> pattern(4, 6);
	pattern[0].push_back(5);
	pattern[0].push_back(2);
	pattern[0].push_back(2);
	pattern[3].push_back(2);
	pattern.transpose_replace();
	auto& transposed = pattern;
	check(transposed.nrow == 6 && transposed.ncol == 4 && transposed.nnz() == 4, "Bool dimensions/nnz");
	check(transposed[2].nnz() == 3 && transposed[2](0) == 0 && transposed[2](1) == 0 && transposed[2](2) == 3, "Bool indices");

	sparse_mat<move_only, int32_t> owned(3, 5);
	owned[0].reserve(1);
	owned[2].reserve(1);
	owned[0].push_back(4, move_only(7));
	owned[2].push_back(1, move_only(9));
	auto first = owned[0][0].value.get(), second = owned[2][0].value.get();
	move_only::moves = 0;
	owned.transpose_replace();
	auto& moved = owned;
	check(move_only::moves == 2, "Values not moved exactly once");
	check(moved[4][0].value.get() == first && moved[1][0].value.get() == second, "Heap ownership not transferred");

	std::cout << "PASS: copying-transpose equivalence, empty/rectangular/unsorted/duplicate/zero cases, "
		"large rationals, integer/index variants, bool, move-only ownership\n";
} catch (const std::exception& error) {
	std::cerr << error.what() << '\n';
	return 1;
}
