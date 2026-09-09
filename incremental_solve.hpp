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
//      pivot column, kept in MUTUALLY REDUCED form:
//        - before a row is inserted it is reduced against EVERY pivot
//          column it contains (reduce_fully_against_basis), not only its
//          leading ones — an unreduced tail would leave other rows' pivot
//          columns inside the basis and the "pivot var = rhs" read-off
//          would silently return a wrong solution;
//        - after insertion the new pivot column is eliminated from all
//          other basis rows.
//   2. Process constraints in batches of ~ sample_factor * s rows
//      (default: 3). For each batch row: reduce, then drop (0 = 0),
//      declare INCONSISTENT (0 = nonzero), or insert. On rank saturation
//      the batch stops early and ONLY the processed prefix is retired.
//   3. Each round ends with a substitution sweep (row·sol - b) over the
//      not-yet-consumed rows. A FAILING row needs full reduction and
//      forms the next round's batch; a PASSING row is only tentatively
//      satisfied — pass does NOT prove the row lies in the basis span,
//      and structured (repetitive) constraint systems produce false
//      passes at scale.
//   4. Termination therefore requires the outer verification loop: when
//      the rounds quiesce, EVERY constraint is re-checked exactly against
//      the final solution; violations re-enter the reduction loop (below
//      full rank they grow the basis; at full rank a violated row is a
//      genuine 0 = nonzero witness). Only a clean full verification ends
//      the solve.
//   5. Result: particular solution = basis rhs on pivots, free = 0; when
//      the system is underdetermined the rank / null space are completed
//      with one batched full RREF (reporting only — the particular
//      solution is already fully verified).
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
//
// IMPORTANT: this eliminates only LEADING pivot columns (it stops at the
// first non-pivot column). A row reduced this way may still contain pivot
// columns at LATER positions — use reduce_fully_against_basis before
// inserting a row into the basis, otherwise the basis never becomes
// mutually reduced and the RREF read-off (solution[pivot] = rhs) is wrong.
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

// Eliminate EVERY pivot column the row contains, wherever it sits. Each
// elimination subtracts a basis row whose non-pivot entries can reintroduce
// other pivot columns, so iterate until no pivot column remains (at most
// rank passes; rows are tiny). Guarantees the row is mutually reduced
// against the basis — the invariant the RREF solution extraction needs.
template <typename T, typename index_t>
T reduce_fully_against_basis(sparse_vec<T, index_t>& row, T rhs,
                             const std::vector<basis_row_t<T, index_t>>& basis) {
	if (basis.empty() || row.nnz() == 0)
		return rhs;
	std::vector<index_t> pivot_cols;
	pivot_cols.reserve(basis.size());
	for (const auto& br : basis)
		pivot_cols.push_back(br.pivot_col);

	bool changed = true;
	while (changed && row.nnz() > 0) {
		changed = false;
		for (size_t j = 0; j < row.nnz(); ) {
			index_t c = row(j);
			auto it = std::lower_bound(pivot_cols.begin(), pivot_cols.end(), c);
			if (it == pivot_cols.end())
				break;  // pivot_cols sorted: no further pivot columns exist
			if (*it != c) {
				j++;
				continue;
			}
			size_t r = (size_t)(it - pivot_cols.begin());
			T factor = row[j];
			sub_scaled(row, factor, basis[r].coeffs);
			rhs = rhs - factor * basis[r].rhs;
			changed = true;
			// sub_scaled may have shifted/removed entries: rescan from the start
			j = 0;
		}
	}
	return rhs;
}

// Row-count threshold at which the batch pre-reduction and the verification
// sweeps are spread over the thread pool; below it the overhead outweighs
// the work. Purely a performance knob — results are identical either way.
inline constexpr size_t kParallelRows = 64;

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
	//     row·sol = sum a_i rhs_i = b) — but the CONVERSE fails: rows
	//     outside the span can pass too, so a clean sweep alone does not
	//     prove the system solved.
	//   - hence the final full verification below: every constraint is
	//     re-checked exactly against the final solution, and violations
	//     re-enter the reduction loop. Termination: below full rank every
	//     violated row either grows the basis or returns INCONSISTENT; at
	//     full rank it can only be INCONSISTENT (a satisfied row reduces
	//     to 0 = 0), so each outer iteration makes progress.
	std::vector<size_t> pending = nontrivial_indices;
	size_t round = 0;
	size_t consumed_total = 0;
	sparse_mat<T, index_t> null_space;

	for (;;) {  // outer loop: rounds, then exact full verification
	while (!pending.empty()) {
		round++;
		size_t take = std::min(batch_size, pending.size());
		size_t processed = take;

		// ---- Phase A: pre-reduce the batch against the FROZEN basis ----
		// Each row's reduction is a read-only computation w.r.t. the basis as
		// of the start of the batch, so the rows are independent and reduce
		// in parallel. Phase B (insert loop below) re-reduces every row
		// against the by-then-grown basis: over exact rationals the reduction
		// is confluent, so "fully-reduce(original, grown)" and
		// "re-reduce(pre-reduced, grown)" coincide bit-for-bit, and the basis
		// evolution, printed diagnostics and solution are IDENTICAL to the
		// purely sequential interleaving. Below the row threshold the
		// pre-reduction simply runs on the calling thread.
		std::vector<sparse_vec<T, index_t>> red_rows(take);
		std::vector<T> red_rhs(take);
		{
			const std::vector<size_t>& pend = pending;
			auto preduce = [&](size_t k) {
				red_rows[k] = M[pend[k]];   // copy: M stays intact
				red_rhs[k] = incremental::reduce_against_basis(red_rows[k], get_rhs(pend[k]), basis);
				red_rhs[k] = incremental::reduce_fully_against_basis(red_rows[k], red_rhs[k], basis);
			};
			if (take >= incremental::kParallelRows) {
				pool->detach_loop(0, take, preduce);
				pool->wait();
			} else {
				for (size_t k = 0; k < take; k++) preduce(k);
			}
		}

		// ---- Phase B: sequential insert against the growing basis ----
		for (size_t k = 0; k < take; k++) {
			size_t idx = pending[k];
			sparse_vec<T, index_t> row = std::move(red_rows[k]);
			T rhs = red_rhs[k];
			rhs = incremental::reduce_against_basis(row, rhs, basis);

			if (row.nnz() == 0) {
				if (rhs != (T)0) {
					std::cout << "   INCONSISTENT at constraint " << idx
					          << " (reduces to 0 = " << rhs << ")" << std::endl;
					return {.consistent = false, .unique = false};
				}
				continue;  // already implied by the basis
			}

			// Fully reduce against ALL existing pivots (not only the leading
			// ones): otherwise the inserted row keeps other rows' pivot
			// columns, the basis never becomes mutually reduced, and the
			// RREF solution extraction (pivot var = rhs) reads off a wrong
			// solution even though every stored row is mathematically valid.
			rhs = incremental::reduce_fully_against_basis(row, rhs, basis);
			if (row.nnz() == 0) {
				if (rhs != (T)0) {
					std::cout << "   INCONSISTENT at constraint " << idx
					          << " (fully reduces to 0 = " << rhs << ")" << std::endl;
					return {.consistent = false, .unique = false};
				}
				continue;
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
				// break out to the cheap substitution sweep. Only the rows
				// processed so far may be retired below — the unprocessed
				// rest of the batch still needs the sweep.
				processed = k + 1;
				break;
			}
		}
		consumed_total += processed;
		pending.erase(pending.begin(), pending.begin() + (long)processed);

		// Current particular solution for the substitution sweep.
		sparse_vec<T, index_t> sol;
		for (const auto& br : basis) {
			if (br.rhs != (T)0) sol.push_back(br.pivot_col, br.rhs);
		}
		sol.compress();

		// Substitution sweep row·sol - b over the pending constraints: rows
		// are independent read-only dot products; flag in place, then filter
		// in order so `failed` keeps ascending indices either way.
		std::vector<size_t> failed;
		failed.reserve(pending.size());
		{
			std::vector<char> violated(pending.size(), 0);
			const std::vector<size_t>& pend = pending;
			auto sweep = [&](size_t k) {
				const auto& row = M[pend[k]];
				T residual = -get_rhs(pend[k]);
				for (size_t j = 0; j < row.nnz(); j++) {
					auto p = sol.find(row(j));
					if (p != nullptr) residual = residual + row[j] * (*p);
				}
				if (residual != (T)0) violated[k] = 1;
			};
			if (pending.size() >= incremental::kParallelRows) {
				pool->detach_loop(0, pending.size(), sweep);
				pool->wait();
			} else {
				for (size_t k = 0; k < pending.size(); k++) sweep(k);
			}
			for (size_t k = 0; k < pending.size(); k++)
				if (violated[k]) failed.push_back(pending[k]);
		}

		std::cout << "   Round " << round << ": rank = " << basis.size()
		          << "/" << n_unknowns << ", substitution sweep over " << pending.size()
		          << " constraints -> " << failed.size() << " need further reduction" << std::endl;

		pending = std::move(failed);
	}

	// ---- Final exact verification over ALL constraints ----
	// The per-round substitution sweep only covers not-yet-consumed rows,
	// and a pass (row·sol = b) does NOT prove the row lies in the span of
	// the current basis — structured systems produce false passes that
	// would otherwise be dropped forever. Re-check every constraint against
	// the final solution; violations go back through the reduction loop
	// (at full rank they expose genuine inconsistencies, below full rank
	// they grow the basis), so this terminates.
	{
		sparse_vec<T, index_t> sol;
		for (const auto& br : basis)
			if (br.rhs != (T)0) sol.push_back(br.pivot_col, br.rhs);
		sol.compress();
		// Independent read-only dot products (see the substitution sweep
		// above): flag in parallel, filter in order.
		std::vector<size_t> viol;
		{
			std::vector<char> violated(nontrivial_indices.size(), 0);
			const std::vector<size_t>& nt = nontrivial_indices;
			auto verify = [&](size_t k) {
				const auto& row = M[nt[k]];
				T residual = -get_rhs(nt[k]);
				for (size_t j = 0; j < row.nnz(); j++) {
					auto p = sol.find(row(j));
					if (p != nullptr) residual = residual + row[j] * (*p);
				}
				if (residual != (T)0) violated[k] = 1;
			};
			if (nontrivial_indices.size() >= incremental::kParallelRows) {
				pool->detach_loop(0, nontrivial_indices.size(), verify);
				pool->wait();
			} else {
				for (size_t k = 0; k < nontrivial_indices.size(); k++) verify(k);
			}
			for (size_t k = 0; k < nontrivial_indices.size(); k++)
				if (violated[k]) viol.push_back(nt[k]);
		}
		if (viol.empty())
			break;
		std::cout << "   Final verification: " << viol.size()
		          << " constraint(s) violated by the sampled solution — resuming reduction" << std::endl;
		pending = std::move(viol);
	}
	}  // end outer verification loop
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

	size_t rank_final = basis.size();
	bool unique_final = (rank_final == n_unknowns);

	// Underdetermined: the basis rank may still undercount the true rank
	// (independent rows that passed every sweep by accident are only known
	// to be CONSISTENT, not spanned). Complete rank / null space with one
	// batched full RREF — reporting only; the particular solution above is
	// already verified against every constraint by the outer loop.
	if (!unique_final) {
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
		rank_final = full_pivots.size();
		unique_final = (rank_final == n_unknowns);
		if (!unique_final) {
			null_space = sparse_mat_rref_kernel(M_sub, full_pivots_nested, F, opt).transpose();
		}
		std::cout << "   Full-system rank: " << rank_final << " / " << n_unknowns << std::endl;
	}

	std::cout << "   Rank: " << rank_final << " / " << n_unknowns << std::endl;
	std::cout << "   Solution" << (unique_final ? " (unique):" : " (particular, system underdetermined):") << std::endl;
	for (size_t i = 0; i < solution.nnz(); i++) {
		std::cout << "      c[" << solution(i) << "] = " << solution[i] << std::endl;
	}

	// Null space from the incremental RREF basis is not used any more: the
	// full-RREF completion above computes it whenever the system is
	// underdetermined.
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
