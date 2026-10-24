// solve_symmetry.hpp — Symmetry constraint solving module
//
// Given the projection matrix M of a target (e.g. SEW_5p1) under a symmetry
// (cyclic, flip, parity), compute the invariant subspace:
//   invariant = ker(M^T - I)
// where the transpose follows the invariant_space.cpp convention.
//
// Special case: if M is identity (M^T - I = 0), the target is already invariant
// under the symmetry. The invariant space is the full space (identity matrix).
//
// If the projection matrix is not found at output/<symmetry>/<target>.wxf,
// the projection pipeline (run_projection_pipeline) is invoked automatically.

#ifndef SOLVE_SYMMETRY_HPP
#define SOLVE_SYMMETRY_HPP

#include "projection.hpp"

// ========== Subtract identity from a sparse matrix: M -> M - I ==========
// (Adapted from invariant_space.cpp:92-109)
// Uses sparse_mat::find(row, col) instead of the archive's sparse_mat_entry.
template <typename T, typename index_t>
sparse_mat<T, index_t> subtract_identity(const sparse_mat<T, index_t>& mat, const field_t& F) {
	sparse_mat<T, index_t> result = mat;  // copy
	for (size_t i = 0; i < mat.nrow; i++) {
		T* entry = result.find(i, static_cast<index_t>(i));
		if (entry != nullptr) {
			// Diagonal element exists: subtract 1
			*entry = scalar_sub(*entry, (T)1, F);
		} else {
			// Diagonal element doesn't exist: add -1
			result[i].push_back(static_cast<index_t>(i), scalar_neg((T)1, F));
		}
	}
	return result;
}

// ========== Build an n×n identity sparse_mat ==========
template <typename T, typename index_t>
sparse_mat<T, index_t> identity_sparse_mat(size_t n) {
	sparse_mat<T, index_t> I(n, n);
	for (size_t i = 0; i < n; i++) {
		I[i].push_back(static_cast<index_t>(i), (T)1);
	}
	return I;
}

// ========== Symmetry solver entry point ==========
//
// Flow:
//   1. Ensure output/<symmetry>/<target>.wxf exists (invoke projection pipeline if not)
//   2. Load projection matrix M (must be square)
//   3. Transpose M (invariant_space.cpp convention)
//   4. Compute M^T - I
//   5. If M^T - I is zero: target is invariant, invariant space = identity (full space)
//      Otherwise: RREF + kernel = invariant space
//   6. Post-process (canonicalize, cancel divisor, sort by nnz)
//   7. Write invariant space to output/<symmetry>/<target>_invariant.wxf

template <typename T, typename index_t>
void run_symmetry_solver(
	const std::string& symmetry_name,
	const std::string& target_name,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const field_t& F,
	rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);

	auto sym = get_symmetry_info(symmetry_name);
	auto target = parse_target(target_name);

	std::cout << "========================================" << std::endl;
	std::cout << "Symmetry solving" << std::endl;
	std::cout << "  Symmetry: " << sym.name << std::endl;
	std::cout << "  Target: " << target.name << std::endl;
	std::cout << "  Data dir:   " << data_dir.string() << std::endl;
	std::cout << "  Output dir: " << output_dir.string() << std::endl;
	std::cout << "========================================" << std::endl;

	std::filesystem::path sym_dir = output_dir / sym.name;
	std::filesystem::path projection_path = sym_dir / projection_filename(target);
	std::filesystem::path invariant_path = sym_dir / (target.name + "_invariant.wxf");

	// Step 1: Ensure projection matrix exists
	if (!std::filesystem::exists(projection_path)) {
		std::cout << "Projection matrix not found: " << projection_path << std::endl;
		std::cout << "Running projection pipeline to generate it..." << std::endl;
		run_projection_pipeline<T, index_t>(symmetry_name, target_name, data_dir, output_dir, F, opt);
	}

	// Step 2: Load projection matrix
	std::cout << std::endl << "--- Loading projection matrix ---" << std::endl;
	std::cout << "   " << projection_path << std::endl;
	auto M_csr = projection_read_tensor<T, index_t>(projection_path, F, pool);
	print_tensor_info(M_csr);

	if (M_csr.rank() != 2) {
		throw std::runtime_error("Symmetry solver: projection matrix must be rank-2, got rank "
			+ std::to_string(M_csr.rank()));
	}
	if (M_csr.dim(0) != M_csr.dim(1)) {
		throw std::runtime_error("Symmetry solver: projection matrix must be square, got "
			+ std::to_string(M_csr.dim(0)) + "x" + std::to_string(M_csr.dim(1)));
	}

	size_t n = M_csr.dim(0);

	// Step 3: Transpose (following invariant_space.cpp convention)
	std::cout << "-- Transposing projection matrix..." << std::endl;
	sparse_tensor<T, index_t, SPARSE_COO> M_coo(std::move(M_csr));
	M_coo = M_coo.transpose({1, 0});
	auto M_mat = M_coo.to_sparse_mat(pool);
	M_coo.clear();

	// Step 4: Subtract identity: M^T - I
	std::cout << "-- Subtracting identity (M^T - I)..." << std::endl;
	auto M_minus_I = subtract_identity<T, index_t>(M_mat, F);
	M_mat.clear();

	// Canonicalize to remove zero entries (e.g., diagonal 1-1=0)
	pool->detach_loop(0, M_minus_I.nrow, [&](size_t i) {
		M_minus_I[i].canonicalize();
	});
	pool->wait();

	// Step 5: Check if M^T - I is zero (i.e., M was identity)
	bool is_invariant = (M_minus_I.nnz() == 0);

	sparse_mat<T, index_t> invariant_basis;

	if (is_invariant) {
		std::cout << std::endl;
		std::cout << "*** TARGET " << target.name << " IS INVARIANT UNDER "
		          << sym.name << " SYMMETRY ***" << std::endl;
		std::cout << "    (M^T - I is null, projection matrix is identity)" << std::endl;
		std::cout << "    Invariant space = full space (dimension " << n << ")" << std::endl;
		std::cout << std::endl;
		invariant_basis = identity_sparse_mat<T, index_t>(n);
	} else {
		std::cout << "-- Computing invariant space (RREF + kernel of M^T - I)..." << std::endl;
		Timer timer;
		timer.start();
		auto pivots = sparse_mat_rref_reconstruct(M_minus_I, opt);
		// sparse_mat_rref_kernel returns null vectors as columns; transpose to rows
		auto kernel = sparse_mat_rref_kernel(M_minus_I, pivots, F, opt).transpose();
		M_minus_I.clear();
		timer.stop();
		std::cout << "   Kernel: " << kernel.nrow << "x" << kernel.ncol
		          << " (" << timer.milliseconds() << " ms)" << std::endl;

		// Step 6: Post-process (canonicalize, cancel divisor, sort by nnz)
		std::cout << "-- Post-processing (cancel divisor, sort by nnz)..." << std::endl;
		timer.start();
		pool->detach_loop(0, kernel.nrow, [&](size_t i) {
			kernel[i].canonicalize();
			vec_cancel_divisor(kernel[i]);
		});
		pool->wait();
		kernel.sort_rows_by_nnz();
		timer.stop();
		std::cout << "   Post-process: " << timer.milliseconds() << " ms" << std::endl;

		invariant_basis = std::move(kernel);
		std::cout << "   Invariant space dimension: " << invariant_basis.nrow
		          << " / " << n << std::endl;
	}

	// Step 7: Write invariant space
	std::cout << std::endl << "--- Invariant space ---" << std::endl;
	sparse_tensor<T, index_t, SPARSE_COO> inv_coo(invariant_basis);
	print_tensor_info(inv_coo);
	auto inv_dims = inv_coo.dims();
	auto inv_nnz = inv_coo.nnz();
	projection_write_tensor(
		invariant_path,
		sparse_tensor<T, index_t, SPARSE_CSR>(std::move(inv_coo), pool), pool);

	// Summary
	std::cout << std::endl << "========================================" << std::endl;
	std::cout << "Symmetry solving complete." << std::endl;
	std::cout << "  Target: " << target.name << std::endl;
	std::cout << "  Symmetry: " << sym.name << std::endl;
	std::cout << "  Invariant: " << (is_invariant ? "YES (full space)" : "partial") << std::endl;
	std::cout << "  Invariant dimension: " << invariant_basis.nrow << " / " << n << std::endl;
	std::cout << "  Output: " << invariant_path << std::endl;
	std::cout << "========================================" << std::endl;
}

#endif // SOLVE_SYMMETRY_HPP
