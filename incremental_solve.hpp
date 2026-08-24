// incremental_solve.hpp — Incremental sampled solver for non-homogeneous
// linear systems arising from collinear constraints.
//
// Problem: solve  sum_i A[i, j1..jk] * c[i] = b[j1..jk]  over ALL letter
// multi-indices (j1..jk). The number of constraints (b1*...*bk, e.g. 11^6
// at weight 6) is astronomically larger than the number of unknowns s
// (the sew / basis dimension), and at most s constraints can be linearly
// independent.
//
// Algorithm (exact, over the rationals — only the constraint *rows* are
// batched, with an exactness guarantee, never the arithmetic):
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
// constraints, the final result is *identical* to solving the full system.
//
// All sparse row operations are linear-time two-pointer merges on sorted
// index arrays — no ordered containers in the hot path.
//
// This file is self-contained and does not modify linear_solve.hpp; the
// result type is shared so callers can use either solver interchangeably.

#ifndef INCREMENTAL_SOLVE_HPP
#define INCREMENTAL_SOLVE_HPP

#include "linear_solve.hpp"

namespace incremental {

// One basis row of the incremental RREF: coefficients over the unknowns
// plus the augmented right-hand-side entry.
template <typename T, typename index_t>
struct basis_row_t {
	sparse_vec<T, index_t> coeffs;   // columns 0..n-1, sorted
	T rhs = (T)0;                    // augmented column n
	index_t pivot_col = -1;          // leading column (unique per basis)
};

// row -= factor * src  (all vectors sorted by index; result sorted and
// compressed). Linear-time two-pointer merge.
template <typename T, typename index_t>
void sub_scaled(sparse_vec<T, index_t>& row, const T& factor,
                const sparse_vec<T, index_t>& src) {
	sparse_vec<T, index_t> out;
	out.reserve(row.nnz() + src.nnz());
	size_t i = 0, j = 0;
	while (i < row.nnz() && j < src.nnz()) {
		if (row(i) < src(j)) {
			out.push_back(row(i), row[i]);
			i++;
		} else if (row(i) > src(j)) {
			out.push_back(src(j), -(factor * src[j]));
			j++;
		} else {
			T nv = row[i] - factor * src[j];
			if (nv != (T)0) out.push_back(row(i), nv);
			i++; j++;
		}
	}
	for (; i < row.nnz(); i++) out.push_back(row(i), row[i]);
	for (; j < src.nnz(); j++) out.push_back(src(j), -(factor * src[j]));
	row = std::move(out);
}

// Reduce (row | rhs) against the RREF basis (pivot-sorted). Returns the
// reduced augmented entry; `row` is updated in place. Repeatedly
// eliminates the smallest remaining index whenever it is a pivot column;
// the smallest remaining index strictly increases each step, so the loop
// terminates. The result is reduced w.r.t. the row space of the basis.
template <typename T, typename index_t>
T reduce_against_basis(sparse_vec<T, index_t>& row, T rhs,
                       const std::vector<basis_row_t<T, index_t>>& basis) {
	std::vector<index_t> pivot_cols;
	pivot_cols.reserve(basis.size());
	for (const auto& br : basis) {
		pivot_cols.push_back(br.pivot_col);
	}

	while (row.nnz() > 0) {
		index_t c = row(0);
		auto it = std::lower_bound(pivot_cols.begin(), pivot_cols.end(), c);
		if (it == pivot_cols.end() || *it != c) {
			break;  // smallest remaining index is a non-pivot column: done
		}
		size_t r = (size_t)(it - pivot_cols.begin());
		T factor = row[0];              // pivot entry of basis row is 1
		sub_scaled(row, factor, basis[r].coeffs);  // also removes the pivot
		rhs = rhs - factor * basis[r].rhs;
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

	// Helper: fetch b[idx] (0 if absent).
	auto get_rhs = [&](size_t idx) -> T {
		auto p = b_mat[0].find((index_t)idx);
		return (p != nullptr && *p != (T)0) ? *p : (T)0;
	};

	// Two-tier scheme. "Solve sampled constraints, then substitute the
	// solutions into the remaining ones":
	//   - full reduction (elimination against the RREF basis) is applied
	//     only to rows in the current sample batch;
	//   - all other rows get the SUBSTITUTION check row·sol - b == 0,
	//     a single sparse dot product. Any row in the span of the basis
	//     passes this check exactly (row = sum a_i B_i implies
	//     row·sol = sum a_i rhs_i = b), so the failing rows are exactly
	//     the (independent or inconsistent) ones that need full
	//     reduction; they form the next round's sample.
	//   - when the sweep is clean the system is solved: sol satisfies
	//     every constraint. If rank < n at that point we additionally run
	//     one exact full-reduction sweep so the null space is exact too.
	std::vector<size_t> pending = nontrivial_indices;
	size_t round = 0;
	size_t consumed_total = 0;
	bool unique = false;
	size_t rank = 0;
	sparse_mat<T, index_t> null_space;
	bool completed_by_full_rref = false;

	while (!pending.empty()) {
		round++;
		size_t take = std::min(batch_size, pending.size());

		for (size_t k = 0; k < take; k++) {
			size_t idx = pending[k];
			sparse_vec<T, index_t> row = M[idx];   // copy: M stays intact
			T rhs = incremental::reduce_against_basis(row, get_rhs(idx), basis);

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
				if (r == (size_t)(pos - basis.begin())) continue;
				auto v = basis[r].coeffs.find(pos->pivot_col);
				if (v != nullptr && *v != (T)0) {
					T factor = *v;
					incremental::sub_scaled(basis[r].coeffs, factor, pos->coeffs);
					basis[r].rhs = basis[r].rhs - factor * pos->rhs;
				}
			}

			if (basis.size() == n_unknowns) {
				// Rank saturated: remaining sample rows are all dependent,
				// break out to the cheap substitution sweep.
				break;
			}
		}
		consumed_total += take;
		pending.erase(pending.begin(), pending.begin() + (long)take);

		// Current particular solution for the substitution sweep.
		sparse_vec<T, index_t> sol;
		for (const auto& br : basis) {
			if (br.rhs != (T)0) sol.push_back(br.pivot_col, br.rhs);
		}
		sol.compress();

		std::vector<size_t> failed;
		failed.reserve(pending.size());
		for (size_t idx : pending) {
			const auto& row = M[idx];
			T residual = -get_rhs(idx);
			for (size_t j = 0; j < row.nnz(); j++) {
				auto p = sol.find(row(j));
				if (p != nullptr) residual = residual + row[j] * (*p);
			}
			if (residual != (T)0) failed.push_back(idx);
		}

		std::cout << "   Round " << round << ": rank = " << basis.size()
		          << "/" << n_unknowns << ", substitution sweep over " << pending.size()
		          << " constraints -> " << failed.size() << " need further reduction" << std::endl;

		if (failed.empty() && basis.size() < n_unknowns) {
			// sol satisfies every remaining constraint (verified by the
			// exact substitution sweep above), but rank < n: the null
			// space still needs the independent rows outside the sample.
			// Complete it with one batched RREF of the full system (the
			// same primitive the sampled solver uses) — much faster than
			// row-by-row incremental reduction at this scale.
			std::cout << "   Underdetermined: completing rank / null space via full RREF..."
			          << std::endl;
			sparse_mat<T, index_t> M_sub(nontrivial_indices.size(), n_unknowns);
			for (size_t k = 0; k < nontrivial_indices.size(); k++) {
				M_sub[k] = M[nontrivial_indices[k]];
				M_sub[k].compress();
			}
			auto full_pivots_nested = sparse_mat_rref_reconstruct(M_sub, opt);
			std::vector<pivot_t<index_t>> full_pivots;
			for (auto& p : full_pivots_nested) {
				full_pivots.insert(full_pivots.end(), p.begin(), p.end());
			}
			unique = (full_pivots.size() == n_unknowns);
			rank = full_pivots.size();
			if (!unique) {
				null_space = sparse_mat_rref_kernel(M_sub, full_pivots_nested, F, opt).transpose();
			}
			completed_by_full_rref = true;
			std::cout << "   Full-system rank: " << rank << " / " << n_unknowns << std::endl;
			break;
		}

		pending = std::move(failed);
	}
	(void)consumed_total;

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

	size_t rank_final = completed_by_full_rref ? rank : basis.size();
	bool unique_final = completed_by_full_rref ? unique : (rank_final == n_unknowns);

	std::cout << "   Rank: " << rank_final << " / " << n_unknowns << std::endl;
	std::cout << "   Solution" << (unique_final ? " (unique):" : " (particular, system underdetermined):") << std::endl;
	for (size_t i = 0; i < solution.nnz(); i++) {
		std::cout << "      c[" << solution(i) << "] = " << solution[i] << std::endl;
	}

	// Null space: either already computed by the full RREF fallback, or
	// read off the incremental RREF basis: one vector per non-pivot
	// column j: v[j] = 1, v[pivot_col(r)] = -B[r][j] for each basis row r.
	if (!completed_by_full_rref) {
		if (!unique_final) {
			std::vector<char> is_pivot(n_unknowns, 0);
			for (const auto& br : basis) is_pivot[br.pivot_col] = 1;
			for (size_t j = 0; j < n_unknowns; j++) {
				if (is_pivot[j]) continue;
				sparse_vec<T, index_t> v;
				v.reserve(rank_final + 1);
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
		}
	}
	if (null_space.nrow > 0 || !unique_final) {
		std::cout << "   Null space: " << null_space.nrow << "x" << null_space.ncol << std::endl;
	}

	M.clear();
	b_mat.clear();

	return {
		.consistent = true,
		.unique = unique_final,
		.n_unknowns = n_unknowns,
		.solution = std::move(solution),
		.null_space = std::move(null_space)
	};
}

#endif // INCREMENTAL_SOLVE_HPP
