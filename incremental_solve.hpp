// incremental_solve.hpp — Incremental sampled solver for non-homogeneous
// linear systems arising from collinear constraints.
//
// Problem: solve  sum_i A[i, j1..jk] * c[i] = b[j1..jk]  over ALL letter
// multi-indices (j1..jk). The number of constraints (b1*...*bk, e.g. 11^6
// at weight 6) is astronomically larger than the number of unknowns s
// (the sew / basis dimension), and at most s constraints can be linearly
// independent.
//
// Algorithm (exact, over the rationals — no sampling of the *variables*,
// only of the constraint rows, with an exactness guarantee):
//
//   1. Maintain an incrementally grown RREF basis B of consumed constraint
//      rows (augmented with the b-column). B has at most s rows, one per
//      pivot column, always kept in reduced row echelon form.
//   2. Process constraints in batches of ~ sample_factor * s rows
//      (default sample_factor = 3). For each row r:
//        a. Reduce r against B (row -= r[pivot_col] * B[pivot_row], incl.
//           the augmented b entry). This is exactly "substituting the
//           current partial solution into the remaining constraints".
//        b. If the reduced row is all-zero:
//             - residual 0   -> constraint already implied, drop it;
//             - residual !=0 -> system INCONSISTENT, terminate.
//        c. Otherwise the reduced row is provably linearly independent of
//           B (its leading entry sits in a non-pivot column): normalize,
//           append to B, and eliminate its pivot column from all existing
//           basis rows to restore RREF.
//   3. When all non-trivial constraints are consumed:
//        - particular solution: pivot variables = basis rhs, free = 0;
//        - null space: read directly off the RREF basis (one vector per
//          non-pivot column) — no second, full-size RREF is ever needed.
//
// Because every dropped row was reduced to zero by exact rational
// elimination against rows spanning the same space as all previously seen
// constraints, the final result is *identical* to solving the full system,
// while the elimination work is O(rank) sparse row operations per batch
// instead of one monolithic RREF over all constraints.
//
// This file is self-contained and does not modify linear_solve.hpp; the
// result type is shared so callers can use either solver interchangeably.

#ifndef INCREMENTAL_SOLVE_HPP
#define INCREMENTAL_SOLVE_HPP

#include "linear_solve.hpp"

#include <map>

namespace incremental {

// dst += coef * src  (sparse vectors, dst kept sorted and compressed).
// Implemented through a sorted map so entries cancelled to zero disappear.
template <typename T, typename index_t>
void axpy(sparse_vec<T, index_t>& dst, const T& coef, const sparse_vec<T, index_t>& src) {
	std::map<index_t, T> acc;
	for (size_t j = 0; j < dst.nnz(); j++) {
		acc[dst(j)] += dst[j];
	}
	for (size_t j = 0; j < src.nnz(); j++) {
		acc[src(j)] += coef * src[j];
	}
	dst.clear();
	for (const auto& [idx, val] : acc) {
		if (val != (T)0) {
			dst.push_back(idx, val);
		}
	}
}

// One basis row of the incremental RREF: coefficients over the unknowns
// plus the augmented right-hand-side entry.
template <typename T, typename index_t>
struct basis_row_t {
	sparse_vec<T, index_t> coeffs;   // columns 0..n-1, sorted
	T rhs = (T)0;                    // augmented column n
	index_t pivot_col = -1;          // leading column (unique per basis)
};

// Reduce (row | rhs) against the basis. Returns the reduced augmented
// entry; `row` is updated in place to the reduced coefficient part.
// Basis rows are kept pivot-sorted; binary-search the pivot columns.
template <typename T, typename index_t>
T reduce_against_basis(sparse_vec<T, index_t>& row, T rhs,
                       const std::vector<basis_row_t<T, index_t>>& basis) {
	if (basis.empty()) {
		if (row.nnz() > 0) {
			// caller assigns pivot; keep row as is
			return rhs;
		}
		return rhs;
	}
	std::vector<index_t> pivot_cols;
	pivot_cols.reserve(basis.size());
	for (const auto& br : basis) {
		pivot_cols.push_back(br.pivot_col);
	}

	std::map<index_t, T> acc;
	for (size_t j = 0; j < row.nnz(); j++) {
		acc[row(j)] += row[j];
	}

	// Repeatedly eliminate the smallest key that is a pivot column. Each
	// elimination strictly increases the smallest non-pivot key, so the
	// loop terminates after at most nnz + total basis fill-in steps.
	while (!acc.empty()) {
		auto first = acc.begin();
		index_t c = first->first;
		auto it = std::lower_bound(pivot_cols.begin(), pivot_cols.end(), c);
		if (it == pivot_cols.end() || *it != c) {
			break;  // smallest remaining key is a non-pivot column: done
		}
		size_t r = (size_t)(it - pivot_cols.begin());
		T factor = first->second;      // pivot entry of basis row is 1 (RREF)
		acc.erase(first);
		for (size_t j = 0; j < basis[r].coeffs.nnz(); j++) {
			index_t cc = basis[r].coeffs(j);
			if (cc == c) continue;     // the pivot itself was just erased
			T nv = acc.count(cc) ? acc[cc] - factor * basis[r].coeffs[j] : -(factor * basis[r].coeffs[j]);
			if (nv != (T)0) {
				acc[cc] = nv;
			} else {
				acc.erase(cc);
			}
		}
		rhs = rhs - factor * basis[r].rhs;
	}

	row.clear();
	for (const auto& [idx, val] : acc) {
		if (val != (T)0) {
			row.push_back(idx, val);
		}
	}
	return rhs;
}

}  // namespace incremental

// Solve sum_i A[i, j...] * c[i] = b[j...] with the incremental sampled
// algorithm described at the top of this file.
//
// A: (n_unknowns, d1, ..., dk) CSR;  b: (d1, ..., dk) CSR.
// sample_factor: batch size = sample_factor * n_unknowns constraint rows.
template <typename T, typename index_t>
linear_solve_result_t<T, index_t> solve_linear_system_incremental(
	sparse_tensor<T, index_t, SPARSE_CSR>&& A_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& b_csr,
	const field_t& F, rref_option_t& opt,
	size_t sample_factor = 3) {

	(void)F;  // arithmetic is exact over T (rat_t); F kept for API parity

	thread_pool* pool = &(opt->pool);

	sparse_tensor<T, index_t, SPARSE_COO> A(std::move(A_csr));
	sparse_tensor<T, index_t, SPARSE_COO> b(std::move(b_csr));

	size_t n_unknowns = A.dim(0);
	size_t n_constraints = 1;
	for (size_t i = 1; i < A.rank(); i++) {
		n_constraints *= A.dim(i);
	}

	size_t b_size = 1;
	for (size_t i = 0; i < b.rank(); i++) {
		b_size *= b.dim(i);
	}
	if (b_size != n_constraints) {
		throw std::runtime_error("solve_linear_system_incremental: A and b dimensions mismatch (b_size="
			+ std::to_string(b_size) + ", n_constraints=" + std::to_string(n_constraints) + ")");
	}

	std::cout << "-- Incremental linear solve --" << std::endl;
	std::cout << "   Unknowns: " << n_unknowns << std::endl;
	std::cout << "   Constraints: " << n_constraints << std::endl;

	A.reshape({n_unknowns, n_constraints});
	auto A_mat = A.to_sparse_mat(pool);
	A.clear();

	auto M = A_mat.transpose();   // (n_constraints, n_unknowns)
	A_mat.clear();

	b.reshape({1, n_constraints});
	auto b_mat = b.to_sparse_mat(pool);
	b.clear();

	// Non-trivial constraint rows (M row nonzero). Zero-M rows with b != 0
	// (0 = nonzero) are detected separately below.
	std::vector<size_t> nontrivial_indices;
	{
		std::vector<char> has_a(n_constraints, 0);
		for (size_t i = 0; i < M.nrow; i++) {
			if (M[i].nnz() > 0) {
				has_a[i] = 1;
				nontrivial_indices.push_back(i);
			}
		}
		for (size_t j = 0; j < b_mat[0].nnz(); j++) {
			if (!has_a[b_mat[0](j)]) {
				std::cout << "   INCONSISTENT: constraint " << b_mat[0](j)
				          << " has A = 0 but b != 0." << std::endl;
				return {.consistent = false, .unique = false};
			}
		}
	}
	std::cout << "   Non-trivial constraints: " << nontrivial_indices.size() << std::endl;

	if (nontrivial_indices.empty()) {
		std::cout << "   No non-trivial constraints — any solution works" << std::endl;
		linear_solve_result_t<T, index_t> result;
		result.consistent = true;
		result.unique = false;
		result.n_unknowns = n_unknowns;
		result.null_space = sparse_mat<T, index_t>(n_unknowns, n_unknowns);
		for (size_t i = 0; i < n_unknowns; i++) {
			result.null_space[i].push_back((index_t)i, (T)1);
			result.null_space[i].compress();
		}
		return result;
	}

	Timer timer;
	timer.start();

	std::vector<incremental::basis_row_t<T, index_t>> basis;  // pivot-sorted
	basis.reserve(n_unknowns);

	size_t batch_size = std::max<size_t>(sample_factor * n_unknowns, 1);
	size_t n_batches = (nontrivial_indices.size() + batch_size - 1) / batch_size;
	size_t consumed = 0;

	for (size_t batch = 0; batch < n_batches && basis.size() < n_unknowns; batch++) {
		size_t lo = batch * batch_size;
		size_t hi = std::min(lo + batch_size, nontrivial_indices.size());

		for (size_t k = lo; k < hi; k++) {
			size_t idx = nontrivial_indices[k];
			sparse_vec<T, index_t> row = M[idx];   // copy: M stays intact
			T rhs = (T)0;
			auto b_ptr = b_mat[0].find((index_t)idx);
			if (b_ptr != nullptr && *b_ptr != (T)0) {
				rhs = *b_ptr;
			}

			rhs = incremental::reduce_against_basis(row, rhs, basis);

			if (row.nnz() == 0) {
				if (rhs != (T)0) {
					std::cout << "   INCONSISTENT at constraint " << idx
					          << " (reduces to 0 = " << rhs << ")" << std::endl;
					return {.consistent = false, .unique = false};
				}
				continue;  // already implied by the basis
			}

			// Independent row: normalize on its leading column and insert
			// keeping the basis pivot-sorted.
			incremental::basis_row_t<T, index_t> br;
			br.pivot_col = row(0);
			T piv = row[0];
			for (size_t j = 0; j < row.nnz(); j++) {
				br.coeffs.push_back(row(j), row[j] / piv);
			}
			br.coeffs.compress();
			br.rhs = rhs / piv;

			auto pos = std::lower_bound(basis.begin(), basis.end(), br.pivot_col,
				[](const incremental::basis_row_t<T, index_t>& r, index_t c) {
					return r.pivot_col < c;
				});
			pos = basis.insert(pos, std::move(br));

			// Restore RREF: eliminate the new pivot column from the other
			// basis rows (their rhs included).
			for (size_t r = 0; r < basis.size(); r++) {
				if ((index_t)r == (index_t)(pos - basis.begin())) continue;
				auto v = basis[r].coeffs.find(pos->pivot_col);
				if (v != nullptr && *v != (T)0) {
					T factor = *v;
				// coeffs -= factor * pos->coeffs
					std::map<index_t, T> acc;
					for (size_t j = 0; j < basis[r].coeffs.nnz(); j++) {
						acc[basis[r].coeffs(j)] += basis[r].coeffs[j];
					}
					for (size_t j = 0; j < pos->coeffs.nnz(); j++) {
						index_t cc = pos->coeffs(j);
						T nv = acc.count(cc) ? acc[cc] - factor * pos->coeffs[j]
						                     : -(factor * pos->coeffs[j]);
						if (nv != (T)0) acc[cc] = nv; else acc.erase(cc);
					}
					basis[r].coeffs.clear();
					for (const auto& [idx2, val] : acc) {
						if (val != (T)0) basis[r].coeffs.push_back(idx2, val);
					}
					basis[r].rhs = basis[r].rhs - factor * pos->rhs;
				}
			}
		}

		consumed = hi;
		std::cout << "   Batch " << (batch + 1) << "/" << n_batches
		          << ": consumed " << consumed << "/" << nontrivial_indices.size()
		          << " constraints, rank = " << basis.size() << std::endl;
	}

	// If rank saturated early, remaining constraints still need a residual
	// check? No: rank < n_unknowns and basis has rank rows; any further row
	// reduces against a full-rank-in-pivot-columns basis only if its support
	// hits pivot columns. With rank < n there ARE non-pivot columns, so a
	// later row could still be independent — hence we must continue until
	// all rows are consumed unless rank == n (then every row is dependent).
	// The batch loop above already stops early only when basis.size() ==
	// n_unknowns, in which case all remaining rows are dependent BUT their
	// residuals must still be checked for consistency.
	for (size_t k = consumed; k < nontrivial_indices.size(); k++) {
		size_t idx = nontrivial_indices[k];
		sparse_vec<T, index_t> row = M[idx];
		T rhs = (T)0;
		auto b_ptr = b_mat[0].find((index_t)idx);
		if (b_ptr != nullptr && *b_ptr != (T)0) rhs = *b_ptr;
		rhs = incremental::reduce_against_basis(row, rhs, basis);
		if (row.nnz() != 0 || rhs != (T)0) {
			// row.nnz() != 0 is impossible when rank == n; keep the check
			// defensive anyway.
			std::cout << "   INCONSISTENT at constraint " << idx
			          << " (full-rank residual check failed)" << std::endl;
			return {.consistent = false, .unique = false};
		}
	}
	if (consumed < nontrivial_indices.size()) {
		std::cout << "   Full-rank residual check passed for remaining "
		          << (nontrivial_indices.size() - consumed) << " constraints" << std::endl;
	}

	timer.stop();
	std::cout << "   Incremental elimination time: " << timer.milliseconds() << " ms" << std::endl;

	// Extract particular solution: pivot vars = rhs, free vars = 0.
	sparse_vec<T, index_t> solution;
	solution.reserve(n_unknowns);
	for (const auto& br : basis) {
		if (br.rhs != (T)0) {
			solution.push_back(br.pivot_col, br.rhs);
		}
	}
	solution.compress();

	size_t rank = basis.size();
	bool unique = (rank == n_unknowns);

	std::cout << "   Rank: " << rank << " / " << n_unknowns << std::endl;
	std::cout << "   Solution" << (unique ? " (unique):" : " (particular, system underdetermined):") << std::endl;
	for (size_t i = 0; i < solution.nnz(); i++) {
		std::cout << "      c[" << solution(i) << "] = " << solution[i] << std::endl;
	}

	// Null space from the RREF basis: one vector per non-pivot column j:
	// v[j] = 1, v[pivot_col(r)] = -B[r][j] for each basis row r.
	sparse_mat<T, index_t> null_space;
	if (!unique) {
		std::vector<char> is_pivot(n_unknowns, 0);
		for (const auto& br : basis) is_pivot[br.pivot_col] = 1;
		for (size_t j = 0; j < n_unknowns; j++) {
			if (is_pivot[j]) continue;
			sparse_vec<T, index_t> v;
			v.reserve(rank + 1);
			for (const auto& br : basis) {
				auto a = br.coeffs.find((index_t)j);
				if (a != nullptr && *a != (T)0) {
					v.push_back(br.pivot_col, -(*a));
				}
			}
			v.push_back((index_t)j, (T)1);
			v.compress();
			null_space.rows.push_back(std::move(v));
		}
		null_space.nrow = null_space.rows.size();
		null_space.ncol = n_unknowns;
		std::cout << "   Null space: " << null_space.nrow << "x" << null_space.ncol << std::endl;
	}

	M.clear();
	b_mat.clear();

	return {
		.consistent = true,
		.unique = unique,
		.n_unknowns = n_unknowns,
		.solution = std::move(solution),
		.null_space = std::move(null_space)
	};
}

#endif // INCREMENTAL_SOLVE_HPP
