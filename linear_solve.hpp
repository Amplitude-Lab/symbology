// linear_solve.hpp — Universal non-homogeneous linear equation solver
//
// Given an expanded expression tensor A with shape (n_unknowns, d1, ..., dk)
// and a RHS tensor b with shape (d1, ..., dk), solves the system:
//   sum_i A[i, j1, ..., jk] * c[i] = b[j1, ..., jk]  for all (j1, ..., jk)
//
// Optimization: if the number of constraints >> n_unknowns, samples
// 3 * n_unknowns constraints, solves, then verifies against all constraints.
//
// Returns the particular solution (free variables = 0) and, if underdetermined,
// the null space for parameterization.

#ifndef LINEAR_SOLVE_HPP
#define LINEAR_SOLVE_HPP

#include "bootstrap.hpp"
#include "projection.hpp"

template <typename T, typename index_t>
struct linear_solve_result_t {
	bool consistent = false;
	bool unique = false;
	size_t n_unknowns = 0;                  // logical length of the solution vector
	sparse_vec<T, index_t> solution;       // particular solution (free vars = 0)
	sparse_mat<T, index_t> null_space;     // empty if unique
};

// Solve M·c = b where M = A_reshaped^T.
//
// A: expanded expression (n_unknowns, d1, ..., dk), CSR
// b: RHS (d1, ..., dk), CSR
// sample_factor: if n_constraints > sample_factor * n_unknowns, sample
//                 sample_factor * n_unknowns constraints (default: 3)
template <typename T, typename index_t>
linear_solve_result_t<T, index_t> solve_linear_system(
	sparse_tensor<T, index_t, SPARSE_CSR>&& A_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& b_csr,
	const field_t& F, rref_option_t& opt,
	size_t sample_factor = 3) {

	thread_pool* pool = &(opt->pool);

	sparse_tensor<T, index_t, SPARSE_COO> A(std::move(A_csr));
	sparse_tensor<T, index_t, SPARSE_COO> b(std::move(b_csr));

	// Compute dimensions
	size_t n_unknowns = A.dim(0);
	size_t n_constraints = 1;
	for (size_t i = 1; i < A.rank(); i++) {
		n_constraints *= A.dim(i);
	}

	// Verify b dimensions
	size_t b_size = 1;
	for (size_t i = 0; i < b.rank(); i++) {
		b_size *= b.dim(i);
	}
	if (b_size != n_constraints) {
		throw std::runtime_error("solve_linear_system: A and b dimensions mismatch (b_size="
			+ std::to_string(b_size) + ", n_constraints=" + std::to_string(n_constraints) + ")");
	}

	std::cout << "-- Linear solve --" << std::endl;
	std::cout << "   Unknowns: " << n_unknowns << std::endl;
	std::cout << "   Constraints: " << n_constraints << std::endl;

	// Reshape A to (n_unknowns, n_constraints) and convert to sparse_mat
	A.reshape({n_unknowns, n_constraints});
	auto A_mat = A.to_sparse_mat(pool);
	A.clear();

	// M = A_mat^T (n_constraints, n_unknowns)
	auto M = A_mat.transpose();
	A_mat.clear();

	// Reshape b to (1, n_constraints) and extract the single row as sparse_vec
	b.reshape({1, n_constraints});
	auto b_mat = b.to_sparse_mat(pool);
	b.clear();

	// Collect non-trivial constraint indices (where M row has non-zero entries)
	std::vector<size_t> nontrivial_indices;
	for (size_t i = 0; i < M.nrow; i++) {
		if (M[i].nnz() > 0) {
			nontrivial_indices.push_back(i);
		}
	}

	std::cout << "   Non-trivial constraints: " << nontrivial_indices.size() << std::endl;

	// Determine which constraints to use for the initial solve
	size_t max_sample = sample_factor * n_unknowns;
	std::vector<size_t> use_indices;
	if (nontrivial_indices.size() > max_sample) {
		use_indices.assign(nontrivial_indices.begin(),
		                   nontrivial_indices.begin() + max_sample);
		std::cout << "   Sampling " << max_sample << " constraints ("
		          << sample_factor << "x unknowns)" << std::endl;
	} else {
		use_indices = nontrivial_indices;
		std::cout << "   Using all " << use_indices.size() << " non-trivial constraints" << std::endl;
	}

	if (use_indices.empty()) {
		std::cout << "   No non-trivial constraints — any solution works" << std::endl;
		// Trivially consistent, underdetermined
		linear_solve_result_t<T, index_t> result;
		result.consistent = true;
		result.unique = false;
		result.n_unknowns = n_unknowns;
		// Null space = identity (all variables are free)
		result.null_space = sparse_mat<T, index_t>(n_unknowns, n_unknowns);
		for (size_t i = 0; i < n_unknowns; i++) {
			result.null_space[i].push_back((index_t)i, (T)1);
			result.null_space[i].compress();
		}
		return result;
	}

	// Build augmented matrix [M_sample | b_sample] of shape (m, n+1)
	size_t m = use_indices.size();
	sparse_mat<T, index_t> aug(m, n_unknowns + 1);
	for (size_t k = 0; k < m; k++) {
		size_t i = use_indices[k];
		aug[k] = M[i];  // copy row (column indices 0..n-1)
		// Append b[i] at column n_unknowns
		auto b_ptr = b_mat[0].find((index_t)i);
		if (b_ptr != nullptr && *b_ptr != (T)0) {
			aug[k].push_back((index_t)n_unknowns, *b_ptr);
		}
		aug[k].compress();  // sort indices and remove zeros
	}

	// RREF the augmented matrix
	std::cout << "   RREF augmented matrix (" << m << " x " << (n_unknowns + 1) << ")..." << std::endl;
	Timer timer;
	timer.start();
	auto pivots_nested = sparse_mat_rref_reconstruct(aug, opt);
	timer.stop();
	std::cout << "   RREF time: " << timer.milliseconds() << " ms" << std::endl;

	// Flatten pivots
	std::vector<pivot_t<index_t>> flat_pivots;
	for (auto& p : pivots_nested) {
		flat_pivots.insert(flat_pivots.end(), p.begin(), p.end());
	}

	// Derive the pivot structure from the RECONSTRUCTED ROWS (leading
	// entries), not from the reported pivot_t list: SparseRREF's RREF
	// eliminates columns in a free order (sparse_mat_inverse reweights
	// appended columns for exactly this reason), so a pivot reported in the
	// b column does NOT imply inconsistency, and rows are not necessarily
	// normalized to a +1 leading entry.
	// A row is a genuine 0 = 1 witness only when its ONLY nonzero sits in
	// the b column.
	bool consistent = true;
	for (size_t r = 0; r < aug.nrow && consistent; r++) {
		if (aug[r].nnz() == 0)
			continue;
		if (aug[r](0) == (index_t)n_unknowns && aug[r].nnz() == 1)
			consistent = false;
	}

	if (!consistent) {
		std::cout << "   System is INCONSISTENT — no solution" << std::endl;
		return {.consistent = false, .unique = false};
	}

	// Extract particular solution (free variables = 0) from the leading
	// entries: pivot var = leading column, value = (b entry) / (leading
	// coefficient). Rows that are not mutually reduced will be caught by
	// the full verification below.
	sparse_vec<T, index_t> solution;
	solution.reserve(n_unknowns);
	size_t n_lead = 0;
	std::vector<char> lead_seen(n_unknowns, 0);
	for (size_t r = 0; r < aug.nrow; r++) {
		if (aug[r].nnz() == 0)
			continue;
		index_t lead = aug[r](0);
		if (lead >= (index_t)n_unknowns)
			continue;  // already excluded by the consistency check
		if (lead_seen[lead])
			continue;  // duplicate leading column: keep the first
		lead_seen[lead] = 1;
		n_lead++;
		auto val_ptr = aug[r].find((index_t)n_unknowns);
		if (val_ptr != nullptr && *val_ptr != (T)0) {
			solution.push_back(lead, *val_ptr / aug[r][0]);
		}
		// else: solution[lead] = 0 (default, don't store)
	}
	solution.compress();

	bool unique = (n_lead == n_unknowns);

	std::cout << "   Solution" << (unique ? " (unique):" : " (particular, system underdetermined):") << std::endl;
	for (size_t i = 0; i < solution.nnz(); i++) {
		std::cout << "      c[" << solution(i) << "] = " << solution[i] << std::endl;
	}
	for (size_t i = 0; i < n_unknowns; i++) {
		auto ptr = solution.find((index_t)i);
		if (ptr == nullptr) {
			std::cout << "      c[" << i << "] = 0 (free or zero)" << std::endl;
		}
	}

	// Verify against ALL non-trivial constraints
	std::cout << "   Verifying against all " << nontrivial_indices.size()
	          << " non-trivial constraints..." << std::endl;
	bool verified = true;
	size_t fail_at = 0;
	for (size_t idx : nontrivial_indices) {
		// Compute M[idx] · solution
		T residual = (T)0;
		for (size_t j = 0; j < M[idx].nnz(); j++) {
			index_t col = M[idx](j);
			T entry = M[idx][j];
			auto sol_ptr = solution.find(col);
			if (sol_ptr != nullptr) {
				residual = residual + entry * (*sol_ptr);
			}
		}
		// Subtract b[idx]
		auto b_ptr = b_mat[0].find((index_t)idx);
		if (b_ptr != nullptr) {
			residual = residual - *b_ptr;
		}
		if (residual != (T)0) {
			verified = false;
			fail_at = idx;
			break;
		}
	}

	if (verified) {
		std::cout << "   All constraints verified!" << std::endl;
	} else {
		// A failed verification means the extraction above read off a wrong
		// solution (the sampled rows are rank-deficient or not mutually
		// reduced) — NOT that the system is inconsistent. Solve the full
		// system exactly instead of reporting a false verdict.
		std::cout << "   VERIFICATION FAILED at constraint " << fail_at
		          << " — falling back to the exact full-system solve" << std::endl;
		sparse_mat<T, index_t> full_aug(nontrivial_indices.size(), n_unknowns + 1);
		for (size_t k = 0; k < nontrivial_indices.size(); k++) {
			size_t i = nontrivial_indices[k];
			full_aug[k] = M[i];
			auto b_ptr = b_mat[0].find((index_t)i);
			if (b_ptr != nullptr && *b_ptr != (T)0)
				full_aug[k].push_back((index_t)n_unknowns, *b_ptr);
			full_aug[k].compress();
		}
		auto full_pivots_nested = sparse_mat_rref_reconstruct(full_aug, opt);
		(void)full_pivots_nested;
		// Re-extract from the full reconstructed system, same leading-entry
		// logic as above.
		sparse_vec<T, index_t> full_solution;
		full_solution.reserve(n_unknowns);
		size_t full_n_lead = 0;
		std::vector<char> full_seen(n_unknowns, 0);
		for (size_t r = 0; r < full_aug.nrow; r++) {
			if (full_aug[r].nnz() == 0)
				continue;
			if (full_aug[r](0) == (index_t)n_unknowns && full_aug[r].nnz() == 1) {
				std::cout << "   System is INCONSISTENT — no solution (full solve)" << std::endl;
				return {.consistent = false, .unique = false};
			}
			index_t lead = full_aug[r](0);
			if (lead >= (index_t)n_unknowns || full_seen[lead])
				continue;
			full_seen[lead] = 1;
			full_n_lead++;
			auto val_ptr = full_aug[r].find((index_t)n_unknowns);
			if (val_ptr != nullptr && *val_ptr != (T)0)
				full_solution.push_back(lead, *val_ptr / full_aug[r][0]);
		}
		full_solution.compress();
		solution = std::move(full_solution);
		unique = (full_n_lead == n_unknowns);
		// Verify the full-solve solution as well; a second failure at this
		// point means the reconstructed rows are not mutually reduced and
		// the sampled solver cannot produce a trustworthy answer.
		bool verified2 = true;
		for (size_t idx : nontrivial_indices) {
			const auto& row = M[idx];
			T residual = (T)0;
			for (size_t j = 0; j < row.nnz(); j++) {
				auto sol_ptr = solution.find(row(j));
				if (sol_ptr != nullptr)
					residual = residual + row[j] * (*sol_ptr);
			}
			auto b_ptr = b_mat[0].find((index_t)idx);
			if (b_ptr != nullptr)
				residual = residual - *b_ptr;
			if (residual != (T)0) {
				verified2 = false;
				break;
			}
		}
		if (!verified2) {
			throw std::runtime_error(
				"solve_linear_system (sampled): full-system fallback still fails verification — "
				"the reconstructed RREF rows are not mutually reduced; use --solver incremental");
		}
		std::cout << "   Full-system solution verified." << std::endl;
	}

	// Compute null space if underdetermined
	sparse_mat<T, index_t> null_space;
	if (!unique) {
		std::cout << "   Computing null space..." << std::endl;
		// RREF the full M (not augmented) to get null space
		// Use only the non-trivial rows
		sparse_mat<T, index_t> M_sub(nontrivial_indices.size(), n_unknowns);
		for (size_t k = 0; k < nontrivial_indices.size(); k++) {
			M_sub[k] = M[nontrivial_indices[k]];
			M_sub[k].compress();
		}
		auto M_pivots_nested = sparse_mat_rref_reconstruct(M_sub, opt);
		null_space = sparse_mat_rref_kernel(M_sub, M_pivots_nested, F, opt).transpose();
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

#endif // LINEAR_SOLVE_HPP
