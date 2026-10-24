// tensor_expand.hpp — Universal rank-3 tensor expansion via basis chain contraction
//
// Given a starting rank-3 tensor T^{a}_{b, c} and a list of rank-3 basis tensors
// B^{b}_{d, e}, B'^{d}_{f, g}, ..., expands T by recursively contracting axis 1
// with each basis. After each contraction, the new FEC axis goes to position 1
// and the new letter axis goes to position 2 (right after FEC), before the
// existing letters. This keeps letter slots in ascending weight order.
//
// After tensor_contract(current, basis, 1, 0), the contracted axes are:
//   [T.0, T.2, T.3, ..., T.{r-1}, B.1, B.2]
// We permute to:
//   [T.0, B.1, B.2, T.2, T.3, ..., T.{r-1}]
//
// Example: T = (2, 546, 11), bases = [(546, 210, 11), (210, 76, 11), ...]
//   Step 1: contract T.axis=1 with B1.axis=0 → (2, 210, 11, 11) → permute → (2, 210, 11, 11)
//   Step 2: contract .axis=1 with B2.axis=0 → (2, 76, 11, 11, 11) → permute → (2, 76, 11, 11, 11)
//   ...
// Final: (2, 11, 11, ..., 11) with rank = 2 + len(bases) + 1

#ifndef TENSOR_EXPAND_HPP
#define TENSOR_EXPAND_HPP

#include "bootstrap.hpp"

// Build the permutation vector for reordering after contraction.
// tensor_contract(A, B, 1, 0) convention: result = [A.remaining, B.remaining].
// For A rank-r (contracting axis 1) and B rank-3 (contracting axis 0):
//   Contracted: [A.0, A.2, A.3, ..., A.{r-1}, B.1, B.2]
//   (A.1 and B.0 are contracted; A.remaining first, then B.remaining)
// We want: [A.0, B.1, B.2, A.2, A.3, ..., A.{r-1}]
// (new FEC at position 1, new letter at position 2, then existing letters)
// This produces letter slots in ascending weight order [w1, w2, ..., wN].
inline std::vector<size_t> expansion_perm(size_t r) {
	// r = rank before contraction; result has rank r+1
	std::vector<size_t> perm;
	perm.push_back(0);        // A.0 → position 0
	perm.push_back(r - 1);    // B.1 (FEC) → position 1
	perm.push_back(r);        // B.2 (new letter) → position 2
	for (size_t i = 1; i + 1 <= r - 1; i++) {
		perm.push_back(i);    // A.2, A.3, ..., A.{r-1} → positions 3..r
	}
	return perm;
}

// Universal expansion: contract tensor with each basis in order.
// tensor: starting rank-3 CSR tensor (a, b, c)
// basis_paths: list of rank-3 CSR basis files, each (b_i, b_{i+1}, letter_dim)
//   where b_i matches axis 1 of the current tensor
// lenient: when true, a basis whose compressed side (axis 0) does not match
//   the current axis 1 is SKIPPED instead of throwing. This makes an
//   over-specified basis chain harmless: a seed that is already expanded
//   (or only needs the lower-weight bases) simply ignores the extras. A
//   genuinely wrong basis chain still fails the caller's rank check
//   ("seed has N letter slots but rhs has M").
// Returns the expanded CSR tensor.
template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> expand_tensor(
	sparse_tensor<T, index_t, SPARSE_CSR>&& tensor_csr,
	const std::vector<std::filesystem::path>& basis_paths,
	const field_t& F, thread_pool* pool, bool lenient = false) {

	sparse_tensor<T, index_t, SPARSE_COO> current(std::move(tensor_csr));

	std::cout << "-- Tensor expansion --" << std::endl;
	std::cout << "   Starting tensor: rank=" << current.rank() << " dims=";
	for (size_t i = 0; i < current.rank(); i++) {
		std::cout << current.dim(i);
		if (i + 1 < current.rank()) std::cout << "x";
	}
	std::cout << " nnz=" << current.nnz() << std::endl;

	for (size_t step = 0; step < basis_paths.size(); step++) {
		const auto& path = basis_paths[step];
		std::cout << "   Step " << step + 1 << "/" << basis_paths.size()
		          << ": contracting with " << path.filename().string() << std::endl;

		// Read basis as CSR, convert to COO for contraction
		auto basis_csr = projection_read_tensor<T, index_t>(path, F, pool);
		sparse_tensor<T, index_t, SPARSE_COO> basis(std::move(basis_csr));

		std::cout << "      Basis dims=";
		for (size_t i = 0; i < basis.rank(); i++) {
			std::cout << basis.dim(i);
			if (i + 1 < basis.rank()) std::cout << "x";
		}
		std::cout << " nnz=" << basis.nnz() << std::endl;

		// Verify dimension compatibility. In lenient mode a basis whose
		// compressed side does not match axis 1 is not needed for this tensor
		// (already expanded there, or a lower/higher weight than required):
		// skip it and let the caller's rank check surface a genuinely wrong
		// basis chain.
		if (current.dim(1) != basis.dim(0)) {
			if (!lenient) {
				throw std::runtime_error("expand_tensor: dimension mismatch at step "
					+ std::to_string(step + 1) + ": current.axis(1)=" + std::to_string(current.dim(1))
					+ " != basis.axis(0)=" + std::to_string(basis.dim(0)));
			}
			std::cout << "      Skipping " << path.filename().string()
			          << ": axis(1)=" << current.dim(1)
			          << " does not match basis.axis(0)=" << basis.dim(0)
			          << " (not needed for this tensor)" << std::endl;
			continue;
		}

		size_t r = current.rank();

		// Contract axis 1 of current with axis 0 of basis
		auto contracted = tensor_contract(current, basis, 1, 0, F, pool);

		// Permute: [T.0, B.1, T.2, ..., T.{r-1}, B.2]
		auto perm = expansion_perm(r);
		current = contracted.transpose(perm);

		std::cout << "      Result: rank=" << current.rank() << " dims=";
		for (size_t i = 0; i < current.rank(); i++) {
			std::cout << current.dim(i);
			if (i + 1 < current.rank()) std::cout << "x";
		}
		std::cout << " nnz=" << current.nnz() << std::endl;
	}

	std::cout << "   Final: rank=" << current.rank() << " nnz=" << current.nnz() << std::endl;

	// Convert to CSR for output
	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(current), pool);
}

#endif // TENSOR_EXPAND_HPP
