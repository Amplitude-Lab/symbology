// solve_collinear.hpp — Collinear constraint solver
//
// Implements the collinear projection split into finite and divergent parts,
// following the algorithm in collinear_proj.cpp (archive reference).
//
// Stages:
//   1. collinear_proj_step: split a rank-3 collinear expression (basis tensor)
//      into finite/divergent parts using colprojdiv(N-1) on axis 1 and
//      colprojdiv (base seed) on axis 2.
//      Kernel of the constraint matrix = colprojfin (finite combinations).
//      Orthogonal complement = colprojdiv (divergent combinations).
//
//   2. run_collinear_proj_chain: recursively compute colprojfin/colprojdiv
//      for weights 2..N, using seeds at weight 1. The base tensor at each
//      weight is the collinear-reduced basis (first_wN_basis or SEW_FpL_basis).
//
//   3. run_collinear_solver: full orchestrator — compute projections, expand
//      the target tensor, apply the selected projection, solve the linear system.

#ifndef SOLVE_COLLINEAR_HPP
#define SOLVE_COLLINEAR_HPP

#include "bootstrap.hpp"
#include "projection.hpp"
#include "tensor_expand.hpp"
#include "linear_solve.hpp"
#include "incremental_solve.hpp"

#include <set>
#include <map>

// ========== collinear_proj_step: finite/divergent split ==========
//
// Mirrors collinear_proj.cpp (archive). Given:
//   T          — rank-3 collinear expression (basis tensor), shape (a, b, c)
//   CM_first   — map for axis 1 (colprojdiv(N-1)), shape (b, b')  [input, output]
//   CM_second  — map for axis 2 (colprojdiv base seed), shape (c, c')  [input, output]
//
// Algorithm:
//   M1 = contract(T, CM_second, axis=2, axis=0)  → (a, b, c')
//   M2 = contract(T, CM_first, axis=1, axis=0)   → (a, b', c)
//   Reshape M1 to (a, b*c'), M2 to (a, b'*c)
//   Transpose: M1_T (b*c', a), M2_T (b'*c, a)
//   M = join(M1_T, M2_T)  → (b*c' + b'*c, a)
//   RREF M → pivots
//   kernel = null(M).transpose()  → colprojfin (nullity, a), rows = finite combos
//   RREF colprojfin → pivots_k
//   kernel_temp = null(colprojfin)  → (a, rank), columns
//   colprojdiv = kernel_temp.transpose()  → (rank, a), rows = divergent combos
//
// Saved convention: colprojfin/colprojdiv are transposed to (a, nullity)/(a, rank)
// = (input, output) to match the seed file convention (colprojdiv.wxf is 11x2).
//
// Triviality checks (from archive):
//   tensor_first  = (CM_first.nnz() != 0)
//   tensor_second = (CM_second.nnz() != 1)
// If only one map is non-trivial, use just that contraction.

template <typename T, typename index_t>
struct collinear_proj_result_t {
	sparse_mat<T, index_t> colprojfin;  // (nullity, dim0) — rows = finite combinations
	sparse_mat<T, index_t> colprojdiv;  // (rank, dim0) — rows = divergent combinations
};

template <typename T, typename index_t>
collinear_proj_result_t<T, index_t> collinear_proj_step(
	sparse_tensor<T, index_t, SPARSE_CSR>&& T_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& CM_first_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& CM_second_csr,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);

	sparse_tensor<T, index_t, SPARSE_COO> T_coo(std::move(T_csr));
	sparse_tensor<T, index_t, SPARSE_COO> CM_first(std::move(CM_first_csr));
	sparse_tensor<T, index_t, SPARSE_COO> CM_second(std::move(CM_second_csr));

	if (T_coo.rank() != 3) {
		throw std::runtime_error("collinear_proj_step: T must be rank-3.");
	}

	// Triviality checks (asymmetric, matching archive)
	bool tensor_first = (CM_first.nnz() != 0);
	bool tensor_second = (CM_second.nnz() != 1);

	std::cout << "-- Collinear proj step --" << std::endl;
	std::cout << "   T dims: " << T_coo.dim(0) << "x" << T_coo.dim(1) << "x" << T_coo.dim(2) << std::endl;
	std::cout << "   CM_first: " << CM_first.dim(0) << "x" << CM_first.dim(1)
	          << " nnz=" << CM_first.nnz() << (tensor_first ? " (active)" : " (trivial)") << std::endl;
	std::cout << "   CM_second: " << CM_second.dim(0) << "x" << CM_second.dim(1)
	          << " nnz=" << CM_second.nnz() << (tensor_second ? " (active)" : " (trivial)") << std::endl;

	// Dimension check: CM.axis(0) = input must match T's axis
	if (T_coo.dim(1) != CM_first.dim(0) || T_coo.dim(2) != CM_second.dim(0)) {
		throw std::runtime_error("collinear_proj_step: dimension mismatch (T.dim(1)="
			+ std::to_string(T_coo.dim(1)) + " vs CM_first.dim(0)=" + std::to_string(CM_first.dim(0))
			+ ", T.dim(2)=" + std::to_string(T_coo.dim(2)) + " vs CM_second.dim(0)=" + std::to_string(CM_second.dim(0)) + ")");
	}

	sparse_mat<T, index_t> M;

	if (tensor_first && tensor_second) {
		// Both maps active: two contractions, reshape, transpose, join
		auto M1 = tensor_contract(T_coo, CM_second, 2, 0, F, pool);
		auto M2 = tensor_contract(T_coo, CM_first, 1, 0, F, pool);

		std::cout << "   M1: " << M1.dim(0) << "x" << M1.dim(1) << "x" << M1.dim(2) << std::endl;
		std::cout << "   M2: " << M2.dim(0) << "x" << M2.dim(1) << "x" << M2.dim(2) << std::endl;

		// Reshape to 2D and transpose
		std::vector<size_t> dims1 = {M1.dim(0), M1.dim(1) * M1.dim(2)};
		std::vector<size_t> dims2 = {M2.dim(0), M2.dim(1) * M2.dim(2)};
		M1.reshape(dims1);
		M2.reshape(dims2);
		// Convert COO -> CSR first, then to_sparse_mat (CSR version is simpler/safer)
		auto M1_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(M1), pool);
		auto M2_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(M2), pool);
		auto M1_mat = M1_csr.to_sparse_mat();
		auto M2_mat = M2_csr.to_sparse_mat();
		M1_csr.clear();
		M2_csr.clear();

		auto M1_T = M1_mat.transpose();
		auto M2_T = M2_mat.transpose();

		if (M1_T.ncol != M2_T.ncol) {
			throw std::runtime_error("collinear_proj_step: M1_T.ncol != M2_T.ncol");
		}

		M = sparse_mat_join(M1_T, M2_T, pool);
		std::cout << "   M: " << M.nrow << "x" << M.ncol << std::endl;

	} else if (tensor_second) {
		// Only second map active
		auto M1 = tensor_contract(T_coo, CM_second, 2, 0, F, pool);
		std::vector<size_t> dims1 = {M1.dim(0), M1.dim(1) * M1.dim(2)};
		M1.reshape(dims1);
		auto M1_mat = M1.to_sparse_mat(pool);
		M1.clear();
		M = M1_mat.transpose();
		std::cout << "   M (from CM_second only): " << M.nrow << "x" << M.ncol << std::endl;

	} else if (tensor_first) {
		// Only first map active
		auto M2 = tensor_contract(T_coo, CM_first, 1, 0, F, pool);
		std::vector<size_t> dims2 = {M2.dim(0), M2.dim(1) * M2.dim(2)};
		M2.reshape(dims2);
		auto M2_mat = M2.to_sparse_mat(pool);
		M2.clear();
		M = M2_mat.transpose();
		std::cout << "   M (from CM_first only): " << M.nrow << "x" << M.ncol << std::endl;

	} else {
		throw std::runtime_error("collinear_proj_step: both maps are trivial — at least one must be non-trivial");
	}

	T_coo.clear();
	CM_first.clear();
	CM_second.clear();

	// RREF M → pivots ("finite combinations")
	std::cout << "   RREF M..." << std::endl;
	Timer timer;
	timer.start();
	auto pivots_M = sparse_mat_rref_reconstruct(M, opt);
	timer.stop();
	std::cout << "   RREF time: " << timer.milliseconds() << " ms" << std::endl;

	// a = input dimension = M.ncol (T.dim(0)); capture before M is cleared.
	size_t input_dim = M.ncol;

	// kernel = null(M), columns → transpose to rows = colprojfin
	auto kernel = sparse_mat_rref_kernel(M, pivots_M, F, opt);
	M.clear();
	auto colprojfin = kernel.transpose();
	kernel.clear();

	std::cout << "   colprojfin: " << colprojfin.nrow << "x" << colprojfin.ncol << std::endl;

	sparse_mat<T, index_t> colprojdiv;

	if (colprojfin.nrow == 0) {
		// nullity = 0: no finite combinations exist — everything is divergent.
		// sparse_mat_rref_kernel returns (0, 0) when nullity=0; fix dimensions to
		// (0, input_dim) so downstream code knows the input dimension.
		// colprojdiv = identity matrix of size (input_dim, input_dim).
		std::cout << "   colprojfin is empty — all combinations are divergent!" << std::endl;
		colprojfin = sparse_mat<T, index_t>(0, input_dim);
		colprojdiv = sparse_mat<T, index_t>(input_dim, input_dim);
		for (size_t i = 0; i < input_dim; i++) {
			colprojdiv[i].push_back((index_t)i, (T)1);
			colprojdiv[i].compress();
		}
	} else {
		// RREF colprojfin → pivots ("divergent combinations")
		auto pivots_k = sparse_mat_rref_reconstruct(colprojfin, opt);

		// kernel_temp = null(colprojfin), columns
		auto kernel_temp = sparse_mat_rref_kernel(colprojfin, pivots_k, F, opt);

		// colprojdiv = kernel_temp.transpose() → rows = divergent combinations
		colprojdiv = kernel_temp.transpose();
		kernel_temp.clear();
	}

	std::cout << "   colprojdiv: " << colprojdiv.nrow << "x" << colprojdiv.ncol << std::endl;

	if (colprojdiv.nrow == 0) {
		std::cout << "   colprojdiv is empty — all results are finite!" << std::endl;
	}

	return {
		.colprojfin = std::move(colprojfin),
		.colprojdiv = std::move(colprojdiv)
	};
}

// ========== Recursive chain: build colprojfin/colprojdiv for all weights ==========
//
// Weight 1: seeds from data/colprojfin.wxf, data/colprojdiv.wxf
// Weight N (>= 2): collinear_proj_step(base_N, colprojdiv(N-1), colprojdiv_base)
//   base_N = the collinear-reduced basis tensor at weight N
//     (first_wN_basis.wxf for FEC weights, SEW_FpL_basis.wxf for the SEW level)
//   colprojdiv(N-1) from previous step (or data/ for N=2)
//   colprojdiv_base always from data/colprojdiv.wxf
//
// base_paths[i] provides the base tensor for weight (i + 2).
// So base_paths.size() = target_weight - 1 (weights 2..target_weight).
//
// For FEC targets: all weights use colprojfin_w{N}.wxf / colprojdiv_w{N}.wxf
// For SEW targets: FEC weights (2..target_weight-1) use colprojfin_w{N}.wxf,
//   the SEW level (target_weight) uses colprojfin_SEW_{sew_name}.wxf /
//   colprojdiv_SEW_{sew_name}.wxf to avoid name collision with FEC-level files.
//
// Output: output/collinear/

template <typename T, typename index_t>
void run_collinear_proj_chain(
	const std::vector<std::filesystem::path>& base_paths,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const field_t& F, rref_option_t& opt,
	const std::string& sew_name = "") {

	thread_pool* pool = &(opt->pool);
	auto collinear_dir = output_dir / "collinear";
	std::filesystem::create_directories(collinear_dir);

	size_t target_weight = base_paths.size() + 1;  // base_paths covers weights 2..target_weight
	bool is_sew_target = !sew_name.empty();

	// Weight 1: copy seeds to output
	auto seed_fin = data_dir / "colprojfin.wxf";
	auto seed_div = data_dir / "colprojdiv.wxf";

	auto w1_fin = collinear_dir / "colprojfin_w1.wxf";
	auto w1_div = collinear_dir / "colprojdiv_w1.wxf";

	if (!std::filesystem::exists(w1_fin)) {
		std::filesystem::copy_file(seed_fin, w1_fin);
	}
	if (!std::filesystem::exists(w1_div)) {
		std::filesystem::copy_file(seed_div, w1_div);
	}

	std::cout << "== Collinear projection chain (weights 1.." << target_weight << ") ==" << std::endl;
	if (is_sew_target) {
		// sew_name already includes the "SEW_" prefix (e.g. "SEW_3p1"), so the
		// SEW-level file is colprojdiv_<sew_name>.wxf (e.g. colprojdiv_SEW_3p1.wxf),
		// NOT colprojdiv_SEW_SEW_3p1.wxf.
		std::cout << "   SEW target: " << sew_name
		          << " (SEW level uses colprojdiv_" << sew_name << ".wxf naming)" << std::endl;
	}

	for (size_t i = 0; i < base_paths.size(); i++) {
		size_t N = i + 2;
		bool is_sew_level = is_sew_target && (N == target_weight);

		// Determine output file names
		std::filesystem::path fin_N, div_N;
		if (is_sew_level) {
			// sew_name already includes "SEW_" prefix (e.g. "SEW_5p1"),
			// so the file name is colprojfin_SEW_5p1.wxf, not colprojfin_SEW_SEW_5p1.wxf
			fin_N = collinear_dir / ("colprojfin_" + sew_name + ".wxf");
			div_N = collinear_dir / ("colprojdiv_" + sew_name + ".wxf");
		} else {
			fin_N = collinear_dir / ("colprojfin_w" + std::to_string(N) + ".wxf");
			div_N = collinear_dir / ("colprojdiv_w" + std::to_string(N) + ".wxf");
		}

		if (std::filesystem::exists(fin_N) && std::filesystem::exists(div_N)) {
			std::cout << "Weight " << N << (is_sew_level ? " (SEW)" : "") << ": already computed, skipping." << std::endl;
			continue;
		}

		std::cout << "Weight " << N << (is_sew_level ? " (SEW " + sew_name + ")" : "") << ":" << std::endl;

		// Load the collinear-reduced basis at weight N
		const auto& base_path = base_paths[i];
		if (!std::filesystem::exists(base_path)) {
			throw std::runtime_error("run_collinear_proj_chain: base tensor not found: " + base_path.string());
		}
		auto base_tensor = projection_read_tensor<T, index_t>(base_path, F, pool);

		// Load colprojdiv(N-1)
		auto div_prev = collinear_dir / ("colprojdiv_w" + std::to_string(N - 1) + ".wxf");
		auto CM_first = projection_read_tensor<T, index_t>(div_prev, F, pool);

		// Load colprojdiv (base seed, always the same)
		auto CM_second = projection_read_tensor<T, index_t>(seed_div, F, pool);

		// Run collinear_proj_step
		auto result = collinear_proj_step<T, index_t>(
			std::move(base_tensor), std::move(CM_first), std::move(CM_second), F, opt);

		// Write results: transpose to (input, output) convention matching seed files.
		// Create CSR tensors directly from sparse_mat to avoid the COO intermediate
		// (the COO move constructor from CSR has an early return when nnz=0, which
		// leaves dims/rank uninitialized for empty tensors like colprojfin at SEW level).
		// Skip writing if the projection is empty (0 rows) — the WXF writer cannot
		// serialize a 0-nnz tensor, and an empty projection means "no combinations
		// of this type exist" (handled by the solver).
		if (result.colprojfin.nrow > 0) {
			auto fin_mat = result.colprojfin.transpose();
			sparse_tensor<T, index_t, SPARSE_CSR> fin_csr(fin_mat);
			projection_write_tensor(fin_N, std::move(fin_csr), pool);
		} else {
			std::cout << "   colprojfin is empty — skipping write (no finite combinations)." << std::endl;
		}
		if (result.colprojdiv.nrow > 0) {
			auto div_mat = result.colprojdiv.transpose();
			sparse_tensor<T, index_t, SPARSE_CSR> div_csr(div_mat);
			projection_write_tensor(div_N, std::move(div_csr), pool);
		} else {
			std::cout << "   colprojdiv is empty — skipping write (all finite)." << std::endl;
		}
	}
}

// ========== Auto-detect chain base_paths from target name ==========
//
// For SEW_FpL (e.g. SEW_5p1, target_weight = F + L):
//   Weights 2..F: output/collinear/first_w{N}_basis.wxf
//   Weight F+1:   output/collinear/SEW_FpL_basis.wxf  (the target basis)
//
// For FEC_W (e.g. FEC_3, target_weight = W):
//   Weights 2..W: output/collinear/first_w{N}_basis.wxf
//   (last entry is the target basis)
//
// For LEC_W: similar with last_w{N}_basis (not yet fully supported)

inline std::vector<std::filesystem::path> detect_chain_base_paths(
	const target_t& target,
	const std::filesystem::path& output_dir) {

	auto collinear_dir = output_dir / "collinear";
	std::vector<std::filesystem::path> paths;

	size_t target_weight;
	std::filesystem::path target_basis_path;

	if (target.kind == target_kind_t::SEW) {
		target_weight = target.fec_weight + target.lec_weight;
		target_basis_path = collinear_dir / (target.name + "_basis.wxf");
	} else if (target.kind == target_kind_t::FEC) {
		target_weight = target.fec_weight;
		target_basis_path = collinear_dir / ("first_w" + std::to_string(target.fec_weight) + "_basis.wxf");
	} else {
		// LEC: not fully supported yet
		throw std::runtime_error("detect_chain_base_paths: LEC targets not yet supported");
	}

	for (size_t w = 2; w <= target_weight; w++) {
		if (w == target_weight && target.kind == target_kind_t::SEW) {
			paths.push_back(target_basis_path);
		} else {
			paths.push_back(collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf"));
		}
	}

	return paths;
}

// ========== Apply projection to tensor ==========
//
// Apply a projection matrix P (shape (input, output)) to axis 0 of tensor T
// (shape (input, ...)). Result: (output, ...).
// P is stored as a sparse_mat; we convert to a COO sparse_tensor for contraction.
//
// Note: the COO format prepends a row dimension, so a 2D matrix (input, output)
// becomes a rank-3 COO tensor with dims {1, input, output}. Therefore we contract
// axis 1 (input) of P_coo with axis 0 (input) of T_coo.
// Result axes: [P.0=1, P.2=output, T.1, T.2, ...] → CSR consumes the "1" → (output, ...).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> apply_projection_axis0(
	sparse_tensor<T, index_t, SPARSE_CSR>&& T_csr,
	const sparse_mat<T, index_t>& P,
	const field_t& F, thread_pool* pool) {

	sparse_tensor<T, index_t, SPARSE_COO> T_coo(std::move(T_csr));

	// Convert P to a COO tensor for contraction.
	// P_tensor will have dims {1, input, output} and rank 3 (COO prepends row dim).
	sparse_tensor<T, index_t> P_tensor(P, pool);
	sparse_tensor<T, index_t, SPARSE_COO> P_coo(std::move(P_tensor));

	std::cout << "   Applying projection: P=" << P.nrow << "x" << P.ncol
	          << " (input=" << P.nrow << ", output=" << P.ncol << ")"
	          << ", T.axis(0)=" << T_coo.dim(0) << std::endl;

	if (P.nrow != T_coo.dim(0)) {
		throw std::runtime_error("apply_projection_axis0: P.nrow(input)="
			+ std::to_string(P.nrow) + " != T.dim(0)=" + std::to_string(T_coo.dim(0)));
	}

	// Contract the INPUT axis of P_coo with axis 0 (input) of T_coo.
	// P_coo is rank-2 (input, output) — the sparse_mat -> COO conversion
	// does NOT prepend a "1" row dimension — so the input axis is axis 0.
	// (Square/1x1 projections are ambiguous positionally but symmetric
	// under this choice, so axis 0 is always safe here.)
	auto result = tensor_contract(P_coo, T_coo, 0, 0, F, pool);

	std::cout << "   Result: ";
	for (size_t i = 0; i < result.rank(); i++) {
		std::cout << result.dim(i);
		if (i + 1 < result.rank()) std::cout << "x";
	}
	std::cout << std::endl;

	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(result), pool);
}

// ========== apply_colprojdiv_slots: project each letter slot to divergent subspace ==========
//
// The collinear constraint requires matching the DIVERGENT part of the boundary.
// Both A (expanded SEW collinear basis) and the boundary live in the full letter
// space (11-dim slots). We project each slot to the 2-dim divergent subspace via
// colprojdiv_w1 (shape (11, 2), (input, output) convention) before solving c.A = b.
//
// After each contraction, the projected axis moves to the end of the tensor, so
// the next slot to project is always at `first_slot_axis`. The final axis order
// is: [non-projected axes..., projected axes...]. For A = (sew_dim, 11^k) this
// gives (sew_dim, 2^k); for boundary = (11^k) this gives (2^k).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> apply_colprojdiv_slots(
	sparse_tensor<T, index_t, SPARSE_COO>&& tensor_coo,
	const sparse_tensor<T, index_t, SPARSE_COO>& proj_coo,  // (11, m)
	size_t first_slot_axis,
	size_t n_slots,
	const field_t& F, thread_pool* pool) {

	auto result = std::move(tensor_coo);
	// SparseRREF's tensor_contract dereferences A.index(permA[0]) even when
	// A has zero nonzeros (out-of-bounds read, crashes on the all-zero RHS
	// from --rhs 0). An all-zero tensor contracts to an all-zero tensor, so
	// shortcut: keep the non-contracted axes, append the projected axes at
	// the end — the exact dimension layout tensor_contract would produce.
	if (result.nnz() == 0) {
		std::vector<size_t> dims;
		for (size_t a = 0; a < result.rank(); a++) {
			if (a < first_slot_axis || a >= first_slot_axis + n_slots) {
				dims.push_back(result.dim(a));
			}
		}
		for (size_t s = 0; s < n_slots; s++) {
			dims.push_back(proj_coo.dim(1));
		}
		result = sparse_tensor<T, index_t, SPARSE_COO>(dims);
		return result;
	}
	for (size_t s = 0; s < n_slots; s++) {
		result = tensor_contract(result, proj_coo, first_slot_axis, 0, F, pool);
	}
	return result;
}

// ========== Letter-space support filters (sentinels: divergent / finite) ==========
//
// The per-slot contraction above can only express PRODUCT-form projections:
// an entry survives iff EVERY letter slot lands in the projection's column
// space. For "finite" (colprojfin: zero rows on divergent letters) this
// coincides with "all letter indices finite" — the desired semantics.
// But "any letter index divergent" is NOT product-form: a mixed key like
// (0, 3, 7) has one divergent and two finite letters, and no per-slot matrix
// can select it. We therefore implement the divergent selection as a direct
// SUPPORT FILTER on the letter part of each key, leaving all dimensions
// unchanged:
//
//   filter=any-divergent: keep entries whose letter key contains >= 1
//                         divergent letter (exact complement of the finite
//                         part: any-div ∪ all-fin = full support, disjoint)
//   filter=all-finite:    keep entries where every letter index is finite
//                         (identical support to the colprojfin contraction,
//                         but keeps dims and letter labels — handy for
//                         diagnostics and shared plumbing)
//
// The divergent-letter set is derived from data/colprojdiv.wxf: its nonzero
// ROWS are the divergent letters (for E6: rows 0 and 1). Reading the file
// per solver run keeps the filter problem-agnostic (no hardcoded letters).
//
// Note: this is exactly the "project finite, then subtract from the whole
// nonzero support" construction, specialized: support(any-div) =
// support(all) \ support(all-fin), which a direct key scan computes in one
// pass without materializing the finite contraction first.

enum class letter_filter_t { any_divergent, all_finite };

// Read the divergent letter set from data/colprojdiv.wxf (nonzero rows).
template <typename T, typename index_t>
std::set<index_t> load_divergent_letters(
	const std::filesystem::path& data_dir,
	const field_t& F, thread_pool* pool) {

	auto seed_div = data_dir / "colprojdiv.wxf";
	if (!std::filesystem::exists(seed_div)) {
		throw std::runtime_error("load_divergent_letters: divergent-letter seed file not found: "
			+ seed_div.string() + " (data-dir must contain colprojdiv.wxf)");
	}
	// Multi-pair runs call this once per sentinel pair — cache the parsed set
	// per file so only the first call reads and banners the file.
	static std::map<std::string, std::set<index_t>> cache;
	auto cache_key = std::filesystem::weakly_canonical(seed_div).string();
	auto cached = cache.find(cache_key);
	if (cached != cache.end()) {
		std::cout << "   Divergent letters (cached from " << seed_div.filename().string()
		          << "): {";
		bool first = true;
		for (auto l : cached->second) {
			if (!first) std::cout << ", ";
			std::cout << l;
			first = false;
		}
		std::cout << "}" << std::endl;
		return cached->second;
	}
	auto seed_csr = projection_read_tensor<T, index_t>(seed_div, F, pool);
	std::set<index_t> div_letters;
	// CSR tensors expose index_vector(i) (row index at position 0) but not
	// gen_perm() — iterate by flat nonzero index. Guard on val != 0: the wxf
	// reader does not canonicalize explicit zeros on read.
	for (size_t i = 0; i < seed_csr.nnz(); i++) {
		if (seed_csr.val(i) != T(0)) {
			div_letters.insert(seed_csr.index_vector(i)[0]);
		}
	}
	if (div_letters.empty()) {
		std::cout << "   [WARNING] colprojdiv.wxf has no nonzero rows: the divergent-letter"
			<< std::endl
			<< "   set is empty. 'divergent' will keep NOTHING and 'finite' will keep"
			<< std::endl
			<< "   EVERYTHING — this is almost certainly a data mistake." << std::endl;
	}
	std::cout << "   Divergent letters (nonzero rows of " << seed_div.filename().string()
	          << "): {";
	bool first = true;
	for (auto l : div_letters) {
		if (!first) std::cout << ", ";
		std::cout << l;
		first = false;
	}
	std::cout << "}" << std::endl;
	cache[cache_key] = div_letters;
	return div_letters;
}

// Apply a letter-space support filter to `tensor_coo` (COO). `n_slots` letter
// slots start at `first_slot_axis`. Dims, rank and the non-letter axes are
// untouched; only the retained entries survive. Returns the filtered tensor.
template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> apply_letter_filter(
	sparse_tensor<T, index_t, SPARSE_COO>&& tensor_coo,
	letter_filter_t filter,
	size_t first_slot_axis,
	size_t n_slots,
	const std::set<index_t>& div_letters) {

	sparse_tensor<T, index_t, SPARSE_COO> result(tensor_coo.dims());
	size_t kept = 0, total = 0;
	for (auto i : tensor_coo.gen_perm()) {
		total++;
		auto idx = tensor_coo.index_vector(i);
		// any_divergent: keep requires PROOF of a divergent letter (init false).
		// all_finite:   keep unless DISPROVEN by a divergent letter (init true).
		bool keep = (filter == letter_filter_t::all_finite);
		for (size_t a = first_slot_axis; a < first_slot_axis + n_slots; a++) {
			if (div_letters.count((index_t)idx[a]) > 0) {
				keep = (filter == letter_filter_t::any_divergent);
				break;
			}
		}
		if (keep) {
			result.push_back(idx, tensor_coo.val(i));
			kept++;
		}
	}
	std::cout << "      letter filter: kept " << kept << " of " << total << " entries" << std::endl;
	return result;
}

// Dispatch a letter-projection spec onto the (A, b) pair of one solver run:
//   "identity"           — no-op (full letter space)
//   "divergent"          — support filter, keep entries with >= 1 divergent letter
//   "finite"             — support filter, keep entries with all letters finite
//   <file path>          — legacy per-slot contraction with the matrix in the
//                          file (product-form only; data/colprojdiv.wxf keeps
//                          only all-divergent keys, data/colprojfin.wxf only
//                          all-finite keys)
// A_coo has its letter slots starting at axis 1 (leading unknown-count axis);
// b_coo at axis 0.
template <typename T, typename index_t>
void apply_letter_projection_ab(
	const std::string& letter_projection,
	sparse_tensor<T, index_t, SPARSE_COO>& A_coo,
	sparse_tensor<T, index_t, SPARSE_COO>& b_coo,
	size_t n_letter_slots,
	const std::filesystem::path& data_dir,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);

	if (letter_projection == "identity") {
		std::cout << "== letter_projection identity: full letter space (no projection) ==" << std::endl;
		return;
	}
	if (letter_projection == "divergent" || letter_projection == "finite") {
		std::cout << "== Letter-space support filter: " << letter_projection
		          << (letter_projection == "divergent"
		              ? " (keep entries with ANY divergent letter)"
		              : " (keep entries with ALL letters finite)")
		          << " ==" << std::endl;
		auto div_letters = load_divergent_letters<T, index_t>(data_dir, F, pool);
		auto filt = (letter_projection == "divergent")
			? letter_filter_t::any_divergent
			: letter_filter_t::all_finite;
		A_coo = apply_letter_filter<T, index_t>(std::move(A_coo), filt, 1, n_letter_slots, div_letters);
		b_coo = apply_letter_filter<T, index_t>(std::move(b_coo), filt, 0, n_letter_slots, div_letters);
		std::cout << "   A_filtered: rank=" << A_coo.rank() << " dims=";
		for (size_t i = 0; i < A_coo.rank(); i++) std::cout << A_coo.dim(i) << (i + 1 < A_coo.rank() ? "x" : "");
		std::cout << " nnz=" << A_coo.nnz() << std::endl;
		std::cout << "   b_filtered: rank=" << b_coo.rank() << " dims=";
		for (size_t i = 0; i < b_coo.rank(); i++) std::cout << b_coo.dim(i) << (i + 1 < b_coo.rank() ? "x" : "");
		std::cout << " nnz=" << b_coo.nnz() << std::endl;
		return;
	}

	// Legacy: file path → per-slot contraction
	std::filesystem::path letter_proj_path(letter_projection);
	if (!std::filesystem::exists(letter_proj_path)) {
		throw std::runtime_error("apply_letter_projection_ab: letter projection file not found: "
			+ letter_proj_path.string());
	}
	auto letter_proj_csr = projection_read_tensor<T, index_t>(letter_proj_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_COO> letter_proj_coo(std::move(letter_proj_csr));
	std::cout << "== Projecting A and boundary via letter projection ==" << std::endl;
	std::cout << "   letter_projection: " << letter_proj_path.string()
	          << " (rank=" << letter_proj_coo.rank() << " dims=";
	for (size_t i = 0; i < letter_proj_coo.rank(); i++) {
		std::cout << letter_proj_coo.dim(i) << (i + 1 < letter_proj_coo.rank() ? "x" : "");
	}
	std::cout << " nnz=" << letter_proj_coo.nnz() << ")" << std::endl;
	std::cout << "   n_slots=" << n_letter_slots << std::endl;
	A_coo = apply_colprojdiv_slots<T, index_t>(std::move(A_coo), letter_proj_coo, 1, n_letter_slots, F, pool);
	b_coo = apply_colprojdiv_slots<T, index_t>(std::move(b_coo), letter_proj_coo, 0, n_letter_slots, F, pool);
	std::cout << "   A_proj: rank=" << A_coo.rank() << " dims=";
	for (size_t i = 0; i < A_coo.rank(); i++) std::cout << A_coo.dim(i) << (i + 1 < A_coo.rank() ? "x" : "");
	std::cout << " nnz=" << A_coo.nnz() << std::endl;
	std::cout << "   b_proj: rank=" << b_coo.rank() << " dims=";
	for (size_t i = 0; i < b_coo.rank(); i++) std::cout << b_coo.dim(i) << (i + 1 < b_coo.rank() ? "x" : "");
	std::cout << " nnz=" << b_coo.nnz() << std::endl;
}

// ========== Full collinear solver ==========
//
// Orchestrates: compute projections → apply projection → expand → solve.
//
// Parameters:
//   target_path     — path to the target basis tensor (e.g., SEW_5p1_basis.wxf)
//   rhs_path        — path to the RHS tensor (rank k, dims (11,...,11))
//   projection_type — "finite" (use colprojfin) or "divergent" (use colprojdiv)
//   basis_paths     — list of expansion basis files (highest weight first)
//                     e.g. [first_w5_basis, first_w4_basis, first_w3_basis, first_w2_basis]
//   chain_base_paths— list of chain base tensors for colprojfin/colprojdiv computation
//                     (lowest weight first), e.g. [first_w2_basis, ..., SEW_5p1_basis]
//   target_weight   — the weight N of the target (determines which projection to use)
//   data_dir        — directory with seed files (colprojfin.wxf, colprojdiv.wxf)
//   output_dir      — directory for computed projections
//
// The projection at the target weight is applied to the target basis (axis 0),
// then the result is expanded using the basis chain. The first axis is contracted
// with the unknown vector and compared with the RHS.

template <typename T, typename index_t>
void run_collinear_solver(
	const std::filesystem::path& target_path,
	const std::filesystem::path& rhs_path,
	const std::string& projection_type,
	const std::vector<std::filesystem::path>& basis_paths,
	const std::vector<std::filesystem::path>& chain_base_paths,
	size_t target_weight,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const field_t& F, rref_option_t& opt,
	const std::string& sew_name = "",
	const std::string& letter_projection = "identity",
	const std::string& solver = "incremental",
	const std::string& seed_name = "") {

	thread_pool* pool = &(opt->pool);
	auto collinear_dir = output_dir / "collinear";

	// Custom-seed mode (--projection none): the target tensor is used as-is,
	// no seed-space projection. seed_name names the solution output file.
	const bool custom_seed = (projection_type == "none");

	std::cout << "========================================" << std::endl;
	std::cout << "Collinear solver" << std::endl;
	std::cout << "   Target basis: " << target_path.string() << std::endl;
	std::cout << "   RHS: " << rhs_path.string() << std::endl;
	std::cout << "   Projection: " << projection_type
	          << (custom_seed ? " (custom seed: no seed-space projection)" : "")
	          << std::endl;
	std::cout << "   Expansion bases: " << basis_paths.size() << " files" << std::endl;
	if (!custom_seed) {
		std::cout << "   Chain bases: " << chain_base_paths.size() << " files" << std::endl;
	}
	std::cout << "========================================" << std::endl;

	// Step 1: Ensure colprojfin/colprojdiv are computed for the target weight.
	// For SEW targets, the SEW-level projection uses colprojdiv_SEW_{sew_name}.wxf naming.
	// For FEC targets, it uses colprojdiv_w{target_weight}.wxf.
	bool is_sew_target = !sew_name.empty();
	std::filesystem::path proj_fin, proj_div;
	if (is_sew_target) {
		// sew_name already includes "SEW_" prefix (e.g. "SEW_5p1")
		proj_fin = collinear_dir / ("colprojfin_" + sew_name + ".wxf");
		proj_div = collinear_dir / ("colprojdiv_" + sew_name + ".wxf");
	} else {
		proj_fin = collinear_dir / ("colprojfin_w" + std::to_string(target_weight) + ".wxf");
		proj_div = collinear_dir / ("colprojdiv_w" + std::to_string(target_weight) + ".wxf");
	}

	if (!custom_seed && (!std::filesystem::exists(proj_fin) || !std::filesystem::exists(proj_div))) {
		std::cout << "   Projections not found, computing chain..." << std::endl;
		run_collinear_proj_chain<T, index_t>(chain_base_paths, data_dir, output_dir, F, opt, sew_name);
	}

	// Step 2: Load the selected projection
	std::filesystem::path proj_path;
	sparse_mat<T, index_t> proj_mat;
	bool have_projection = false;
	if (custom_seed) {
		// --projection none: no seed-space projection, proj_mat stays empty.
	} else if (projection_type == "finite") {
		proj_path = proj_fin;
	} else if (projection_type == "divergent") {
		proj_path = proj_div;
	} else {
		throw std::runtime_error("Unknown projection type: " + projection_type
			+ " (expected 'finite', 'divergent' or 'none')");
	}

	// Handle empty projection: if the file doesn't exist, the projection is empty
	// (no combinations of that type exist). For "finite": no finite combinations.
	// For "divergent": no divergent combinations (all finite).
	if (!custom_seed && !std::filesystem::exists(proj_path)) {
		std::cout << "========================================" << std::endl;
		std::cout << "Empty " << projection_type << " projection (file not found)." << std::endl;
		if (projection_type == "finite") {
			std::cout << "   No finite combinations exist at weight " << target_weight
			          << " — all combinations are divergent." << std::endl;
		} else {
			std::cout << "   No divergent combinations exist at weight " << target_weight
			          << " — all combinations are finite." << std::endl;
		}
		std::cout << "   Nothing to solve." << std::endl;
		std::cout << "========================================" << std::endl;
		return;
	}

	if (!custom_seed) {
		auto proj_csr = projection_read_tensor<T, index_t>(proj_path, F, pool);
		std::cout << "--- Projection matrix ---" << std::endl;
		print_tensor_info(proj_csr);

		// Convert projection (rank-2 CSR tensor) directly to sparse_mat.
		// Avoid CSR→COO conversion: the COO format prepends a row dimension,
		// which breaks reshape and to_sparse_mat for rank-2 tensors.
		proj_mat = proj_csr.to_sparse_mat();
		proj_csr.clear();
		have_projection = true;
	}

	// Step 3: Load target basis and apply projection
	auto target = projection_read_tensor<T, index_t>(target_path, F, pool);
	std::cout << "--- Target basis ---" << std::endl;
	print_tensor_info(target);
	auto projected = have_projection
		? apply_projection_axis0<T, index_t>(std::move(target), proj_mat, F, pool)
		: std::move(target);  // custom seed: used as-is

	// Step 4: Expand the projected tensor using the basis chain
	auto expanded = expand_tensor<T, index_t>(std::move(projected), basis_paths, F, pool);

	std::cout << "   Expanded expression: rank=" << expanded.rank() << " dims=";
	for (size_t i = 0; i < expanded.rank(); i++) {
		std::cout << expanded.dim(i);
		if (i + 1 < expanded.rank()) std::cout << "x";
	}
	std::cout << std::endl;

	// Step 5: Load RHS (or construct empty RHS for "--rhs 0")
	// Q7: "--rhs 0" means the RHS is an all-zero tensor. We construct it in-memory
	// with shape = expanded.dims[1..] so that b_size = n_constraints (matches A).
	sparse_tensor<T, index_t, SPARSE_CSR> rhs;
	if (rhs_path.string() == "0") {
		std::vector<size_t> b_dims;
		for (size_t i = 1; i < expanded.rank(); i++) {
			b_dims.push_back(expanded.dim(i));
		}
		rhs = sparse_tensor<T, index_t, SPARSE_CSR>(b_dims);  // 0-nnz
		std::cout << "--- RHS: empty (all-zero, dims=";
		for (size_t i = 0; i < b_dims.size(); i++) {
			std::cout << b_dims[i] << (i + 1 < b_dims.size() ? "x" : "");
		}
		std::cout << ") ---" << std::endl;
	} else {
		rhs = projection_read_tensor<T, index_t>(rhs_path, F, pool);
		std::cout << "--- RHS ---" << std::endl;
		print_tensor_info(rhs);
	}

	// Step 5b: Project A + boundary to the letter subspace via --letter-projection.
	// The collinear constraint c.A = b is only enforced in the projected subspace.
	// Accepts "identity" (full space), the support-filter sentinels "divergent" /
	// "finite", or a projection matrix file contracted into each letter slot
	// (e.g. colprojdiv_w1, shape (11, 2)). This is required at L>=3 because E1
	// has divergent-letter entries (E1[0,0]=-2, E1[1,1]=-2), so the boundary
	// has divergent components and solving in the full 11-dim space fails.
	sparse_tensor<T, index_t, SPARSE_COO> A_coo(std::move(expanded));
	sparse_tensor<T, index_t, SPARSE_COO> b_coo(std::move(rhs));

	size_t n_letter_slots = b_coo.rank();   // = 2L
	apply_letter_projection_ab<T, index_t>(
		letter_projection, A_coo, b_coo, n_letter_slots, data_dir, F, opt);

	// Step 5c: Match positions — UNION of supports.
	// We enforce c·A = boundary at EVERY position where either A or boundary
	// is nonzero (not just the intersection). This is the "exact match"
	// semantics: after projecting both sides via --letter-projection, the
	// projected supports should coincide and c·A cancels boundary exactly.
	//
	// Three cases per letter multi-index key:
	//   - Both nonzero (intersection): c·A[key] = b[key]
	//   - A nonzero, b zero (homogeneous): c·A[key] = 0
	//   - A zero, b nonzero (b-only): 0 = b[key] → inconsistent if b[key]≠0
	//
	// The b-only case makes the system trivially inconsistent (0 = nonzero);
	// we detect it here and skip the linear solver, because the solver's
	// "non-trivial constraint" filter (M[i].nnz() > 0) would silently skip
	// these zero-A rows.
	//
	// CRITICAL: preserve the sew axis in A_match — do not collapse into std::map
	// (that overwrites duplicate sew entries and breaks multi-unknown systems).
	std::map<std::vector<index_t>, T> b_map;
	for (auto i : b_coo.gen_perm()) {
		b_map[b_coo.index_vector(i)] = b_coo.val(i);
	}

	sparse_tensor<T, index_t, SPARSE_COO> A_match(A_coo.dims());
	sparse_tensor<T, index_t, SPARSE_COO> b_match(b_coo.dims());
	std::set<std::vector<index_t>> A_keys;          // keys where A ≠ 0
	std::set<std::vector<index_t>> emitted_keys;    // keys already emitted to b_match

	// First pass: iterate A. Emit every A entry; pair with b[key] if present,
	// else with 0 (homogeneous constraint c·A[key] = 0).
	for (auto i : A_coo.gen_perm()) {
		auto full_idx = A_coo.index_vector(i);
		std::vector<index_t> key(full_idx.begin() + 1, full_idx.end());
		A_keys.insert(key);

		auto it = b_map.find(key);
		T b_val = (it != b_map.end()) ? it->second : T(0);

		A_match.push_back(full_idx, A_coo.val(i));
		if (emitted_keys.insert(key).second) {
			if (b_val != T(0)) {
				b_match.push_back(key, b_val);
			}
			// b_val == 0: don't store (canonicalize drops zeros anyway);
			// the solver treats missing b entries as 0 → enforces c·A[key] = 0.
		}
	}

	// Second pass: detect b-only positions (A=0, b≠0).
	// These make the system trivially inconsistent.
	size_t n_b_only = 0;
	for (const auto& [key, b_val] : b_map) {
		if (A_keys.find(key) == A_keys.end()) {
			n_b_only++;
		}
	}

	A_match.canonicalize();
	b_match.canonicalize();

	size_t n_intersection = 0;
	size_t n_homogeneous = 0;
	for (const auto& key : A_keys) {
		if (b_map.find(key) != b_map.end()) {
			n_intersection++;
		} else {
			n_homogeneous++;
		}
	}

	std::cout << "   Union matching:" << std::endl;
	std::cout << "      Both nonzero (intersection): " << n_intersection << std::endl;
	std::cout << "      A nonzero, b zero (homogeneous c·A=0): " << n_homogeneous << std::endl;
	std::cout << "      A zero, b nonzero (INCONSISTENT): " << n_b_only << std::endl;

	// Step 6: Solve c.A_match = b_match
	// If any b-only position exists, the system is trivially inconsistent
	// (0 = nonzero) — skip the linear solver.
	linear_solve_result_t<T, index_t> result;
	if (n_b_only > 0) {
		std::cout << "   System is INCONSISTENT: " << n_b_only
		          << " positions where boundary ≠ 0 but A = 0." << std::endl;
		std::cout << "   (Under union matching, c·A must equal boundary exactly.)" << std::endl;
		result.consistent = false;
		result.unique = false;
	} else {
		auto A_match_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(A_match), pool);
		auto b_match_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(b_match), pool);
		if (solver == "incremental") {
			result = solve_linear_system_incremental<T, index_t>(
				std::move(A_match_csr), std::move(b_match_csr), F, opt);
		} else {
			result = solve_linear_system<T, index_t>(
				std::move(A_match_csr), std::move(b_match_csr), F, opt);
		}
	}

	// Step 6b: Write solMHV_LL.wxf for SEW targets
	if (!sew_name.empty() && result.consistent) {
		size_t L = target_weight / 2;
		auto L_loop_dir = output_dir / (std::to_string(L) + "loop");
		std::filesystem::create_directories(L_loop_dir);
		auto solMHV_path = L_loop_dir / ("solMHV_" + std::to_string(L) + "L.wxf");

		sparse_mat<T, index_t> sol_mat(1, result.n_unknowns);
		for (size_t i = 0; i < result.solution.nnz(); i++) {
			sol_mat[0].push_back(result.solution(i), result.solution[i]);
		}
		sol_mat[0].compress();
		sparse_tensor<T, index_t, SPARSE_CSR> sol_csr(sol_mat);
		projection_write_tensor<T, index_t>(solMHV_path, std::move(sol_csr), pool);
		std::cout << "   Wrote " << solMHV_path.string() << std::endl;
	}

	// Step 6c: Write sol_<seed_name>.wxf for custom seeds (--projection none).
	// The 1 x n_unknowns coefficient vector expands the custom seed tensor to
	// the RHS: c · seed = boundary.
	if (custom_seed && !seed_name.empty() && result.consistent) {
		auto sol_path = collinear_dir / ("sol_" + seed_name + ".wxf");
		sparse_mat<T, index_t> sol_mat(1, result.n_unknowns);
		for (size_t i = 0; i < result.solution.nnz(); i++) {
			sol_mat[0].push_back(result.solution(i), result.solution[i]);
		}
		sol_mat[0].compress();
		sparse_tensor<T, index_t, SPARSE_CSR> sol_csr(sol_mat);
		projection_write_tensor<T, index_t>(sol_path, std::move(sol_csr), pool);
		std::cout << "   Wrote " << sol_path.string() << std::endl;
	}

	// Step 7: Print result
	if (result.consistent) {
		std::cout << "========================================" << std::endl;
		std::cout << "Solution" << std::endl;
		if (result.unique) {
			std::cout << "   Unique solution:" << std::endl;
		} else {
			std::cout << "   Particular solution (system underdetermined):" << std::endl;
			std::cout << "   Null space dimension: " << result.null_space.nrow << std::endl;
		}
		for (size_t i = 0; i < result.solution.nnz(); i++) {
			std::cout << "      c[" << result.solution(i) << "] = " << result.solution[i] << std::endl;
		}
		std::cout << "========================================" << std::endl;
	} else {
		std::cout << "========================================" << std::endl;
		std::cout << "No solution (system inconsistent)" << std::endl;
		std::cout << "========================================" << std::endl;
	}
}

// ========== Multi-pair non-homogeneous constraint solver ==========
//
// Several {seed, rhs, letter_projection} pairs, each producing its own set of
// non-homogeneous constraints c.A^(k) = b^(k) (after per-pair letter
// projection), plus optional pre-computed condition matrices [M | r]. All
// rows are stacked into ONE linear system and solved once:
//
//   row_j:  sum_i c[i] * M[j][i] = r[j]
//
// Each pair may use a DIFFERENT letter projection (e.g. pair 1 identity,
// pair 2 the collinear divergent projection colprojdiv_w1.wxf). Rows from
// different pairs are never merged: every (pair, letter-key) combination
// gets its own global row index, so identical letter keys from different
// pairs impose independent constraints.
//
// Seeds are custom tensors (used as-is on axis 0, like --target-basis /
// --projection none), optionally expanded with a shared --basis chain.

template <typename T, typename index_t>
struct collinear_pair_t {
	std::filesystem::path seed_path;               // custom seed tensor (axis 0 = unknowns)
	std::filesystem::path rhs_path;                // file, or the sentinel "0" (all-zero)
	std::string letter_projection;                 // "identity", "divergent", "finite" or a resolved file path
	std::string stem;                              // seed file stem (naming)
};

template <typename T, typename index_t>
struct cond_row_t {
	std::map<index_t, T> coeffs;                   // unknown index → coefficient
	T rhs = T(0);
};

// Build constraint rows for one pair: load seed → expand → load/construct rhs
// → per-pair letter projection → union matching (same semantics as the
// single-pair solver: enforce c.A = b wherever either side is nonzero).
template <typename T, typename index_t>
void build_pair_rows(
	const collinear_pair_t<T, index_t>& pair,
	const std::vector<std::filesystem::path>& basis_paths,
	std::vector<cond_row_t<T, index_t>>& rows_out,
	size_t& n_unknowns_out,
	size_t& n_b_only_out,
	const std::filesystem::path& data_dir,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);

	std::cout << "---- Pair: " << pair.stem << " ----" << std::endl;

	// Load seed tensor (custom: no seed-space projection on axis 0)
	auto seed = projection_read_tensor<T, index_t>(pair.seed_path, F, pool);
	std::cout << "--- Seed tensor ---" << std::endl;
	print_tensor_info(seed);

	// Expand with the shared basis chain (no-op when basis_paths is empty)
	auto expanded = expand_tensor<T, index_t>(std::move(seed), basis_paths, F, pool);
	std::cout << "   Expanded expression: rank=" << expanded.rank() << " dims=";
	for (size_t i = 0; i < expanded.rank(); i++) {
		std::cout << expanded.dim(i) << (i + 1 < expanded.rank() ? "x" : "");
	}
	std::cout << std::endl;
	if (expanded.rank() < 1) {
		throw std::runtime_error("build_pair_rows: expanded seed must have rank >= 1 (pair " + pair.stem + ")");
	}
	if (n_unknowns_out == 0) {
		n_unknowns_out = expanded.dim(0);
	} else if (expanded.dim(0) != n_unknowns_out) {
		throw std::runtime_error("build_pair_rows: pair '" + pair.stem + "' has n_unknowns="
			+ std::to_string(expanded.dim(0)) + ", expected " + std::to_string(n_unknowns_out)
			+ " (all pairs must share the same unknown-count on axis 0)");
	}

	// Load RHS or construct the all-zero sentinel
	sparse_tensor<T, index_t, SPARSE_CSR> rhs;
	if (pair.rhs_path.string() == "0") {
		std::vector<size_t> b_dims;
		for (size_t i = 1; i < expanded.rank(); i++) {
			b_dims.push_back(expanded.dim(i));
		}
		rhs = sparse_tensor<T, index_t, SPARSE_CSR>(b_dims);
		std::cout << "--- RHS: empty (all-zero, dims=";
		for (size_t i = 0; i < b_dims.size(); i++) {
			std::cout << b_dims[i] << (i + 1 < b_dims.size() ? "x" : "");
		}
		std::cout << ") ---" << std::endl;
	} else {
		rhs = projection_read_tensor<T, index_t>(pair.rhs_path, F, pool);
		std::cout << "--- RHS ---" << std::endl;
		print_tensor_info(rhs);
	}

	sparse_tensor<T, index_t, SPARSE_COO> A_coo(std::move(expanded));
	sparse_tensor<T, index_t, SPARSE_COO> b_coo(std::move(rhs));

	// Per-pair letter projection (the projection that differs between pairs):
	// "identity", a support-filter sentinel ("divergent"/"finite") or a
	// projection matrix file (per-slot contraction).
	size_t n_letter_slots = b_coo.rank();
	if (A_coo.rank() != n_letter_slots + 1) {
		throw std::runtime_error("build_pair_rows: pair '" + pair.stem + "' — seed has "
			+ std::to_string(A_coo.rank() - 1) + " letter slots but rhs has "
			+ std::to_string(n_letter_slots));
	}
	apply_letter_projection_ab<T, index_t>(
		pair.letter_projection, A_coo, b_coo, n_letter_slots, data_dir, F, opt);

	// Union matching → rows (map-based: deterministic key order, sew axis preserved)
	std::map<std::vector<index_t>, T> b_map;
	for (auto i : b_coo.gen_perm()) {
		b_map[b_coo.index_vector(i)] = b_coo.val(i);
	}

	std::map<std::vector<index_t>, std::map<index_t, T>> A_by_key;
	for (auto i : A_coo.gen_perm()) {
		auto full_idx = A_coo.index_vector(i);
		std::vector<index_t> key(full_idx.begin() + 1, full_idx.end());
		A_by_key[key][full_idx[0]] = A_coo.val(i);
	}

	size_t n_intersection = 0;
	size_t n_homogeneous = 0;
	for (const auto& [key, sew_map] : A_by_key) {
		auto it = b_map.find(key);
		T b_val = (it != b_map.end()) ? it->second : T(0);
		if (it != b_map.end()) {
			n_intersection++;
		} else {
			n_homogeneous++;
		}
		cond_row_t<T, index_t> row;
		row.coeffs = sew_map;
		row.rhs = b_val;
		rows_out.push_back(std::move(row));
	}
	for (const auto& [key, b_val] : b_map) {
		if (A_by_key.find(key) == A_by_key.end()) {
			n_b_only_out++;
		}
	}

	std::cout << "   Union matching:" << std::endl;
	std::cout << "      Both nonzero (intersection): " << n_intersection << std::endl;
	std::cout << "      A nonzero, b zero (homogeneous c·A=0): " << n_homogeneous << std::endl;
	std::cout << "      A zero, b nonzero (INCONSISTENT): " << n_b_only_out << std::endl;
}

// Ingest a pre-computed condition matrix [M | r]: rank-2, dims (R, n+1),
// last column = rhs. Rows are appended to the stacked system.
template <typename T, typename index_t>
void ingest_cond_file(
	const std::filesystem::path& cond_path,
	std::vector<cond_row_t<T, index_t>>& rows_out,
	size_t& n_unknowns_out,
	size_t& n_zero_rhs_only_out,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);
	std::cout << "---- Condition file: " << cond_path.filename().string() << " ----" << std::endl;

	auto cond_csr = projection_read_tensor<T, index_t>(cond_path, F, pool);
	if (cond_csr.rank() != 2) {
		throw std::runtime_error("ingest_cond_file: " + cond_path.string()
			+ " must be rank-2 [M | r] (got rank " + std::to_string(cond_csr.rank()) + ")");
	}
	if (cond_csr.dim(1) < 1) {
		throw std::runtime_error("ingest_cond_file: " + cond_path.string() + " has 0 columns");
	}
	if (n_unknowns_out == 0) {
		n_unknowns_out = cond_csr.dim(1) - 1;
	} else if (cond_csr.dim(1) - 1 != n_unknowns_out) {
		throw std::runtime_error("ingest_cond_file: " + cond_path.string() + " has n_unknowns="
			+ std::to_string(cond_csr.dim(1) - 1) + ", expected " + std::to_string(n_unknowns_out));
	}

	auto mat = cond_csr.to_sparse_mat();
	cond_csr.clear();
	for (size_t r = 0; r < mat.nrow; r++) {
		cond_row_t<T, index_t> row;
		for (size_t j = 0; j < mat[r].nnz(); j++) {
			index_t col = mat[r](j);
			T val = mat[r][j];
			if ((size_t)col < n_unknowns_out) {
				row.coeffs[col] = val;
			} else {
				row.rhs = val;
			}
		}
		if (row.coeffs.empty() && row.rhs != T(0)) {
			n_zero_rhs_only_out++;  // 0 = nonzero → inconsistent
			continue;
		}
		if (row.coeffs.empty() && row.rhs == T(0)) {
			continue;  // trivial 0 = 0
		}
		rows_out.push_back(std::move(row));
	}
	std::cout << "   Ingested " << mat.nrow << " rows (" << mat.nrow << "x" << (n_unknowns_out + 1)
	          << " [M | r])" << std::endl;
	mat.clear();
}

template <typename T, typename index_t>
void run_collinear_solver_pairs(
	const std::vector<collinear_pair_t<T, index_t>>& pairs,
	const std::vector<std::filesystem::path>& cond_paths,
	bool export_conditions,
	const std::vector<std::filesystem::path>& basis_paths,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const field_t& F, rref_option_t& opt,
	const std::string& solver = "incremental",
	const std::string& out_stem = "") {

	thread_pool* pool = &(opt->pool);
	auto collinear_dir = output_dir / "collinear";

	std::cout << "========================================" << std::endl;
	std::cout << "Non-homogeneous constraint solver (multi-pair)" << std::endl;
	std::cout << "   Pairs: " << pairs.size() << std::endl;
	std::cout << "   Condition files: " << cond_paths.size() << std::endl;
	std::cout << "   Expansion bases: " << basis_paths.size() << " files" << std::endl;
	std::cout << "   Export conditions: " << (export_conditions ? "yes" : "no") << std::endl;
	std::cout << "   Solver: " << solver << std::endl;
	std::cout << "========================================" << std::endl;

	// Build rows from every pair (per-pair letter projection), then every cond file
	std::vector<cond_row_t<T, index_t>> rows;
	size_t n_unknowns = 0;
	size_t n_b_only = 0;
	for (const auto& pair : pairs) {
		build_pair_rows<T, index_t>(pair, basis_paths, rows, n_unknowns, n_b_only, data_dir, F, opt);
	}
	for (const auto& cond_path : cond_paths) {
		if (!std::filesystem::exists(cond_path)) {
			throw std::runtime_error("run_collinear_solver_pairs: condition file not found: " + cond_path.string());
		}
		ingest_cond_file<T, index_t>(cond_path, rows, n_unknowns, n_b_only, F, opt);
	}

	if (n_unknowns == 0) {
		throw std::runtime_error("run_collinear_solver_pairs: could not determine n_unknowns (empty system?)");
	}

	size_t R = rows.size();
	std::cout << "== Stacked system: " << R << " rows x " << n_unknowns << " unknowns ==" << std::endl;

	// Output naming: --out-stem override; else single source → its stem;
	// multiple → first stem + "_x<count>"
	std::vector<std::string> names;
	for (const auto& pair : pairs) names.push_back(pair.stem);
	for (const auto& p : cond_paths) names.push_back(p.stem().string());
	std::string stem = !out_stem.empty() ? out_stem
		: (names.size() == 1) ? names[0]
		: names[0] + "_x" + std::to_string(names.size());

	// Export the combined conditions [M | r] (rank-2 CSR, rhs = last column)
	if (export_conditions) {
		auto cond_path = collinear_dir / ("cond_" + stem + ".wxf");
		if (R == 0) {
			std::cout << "   Conditions export skipped: 0 rows (nothing to write)." << std::endl;
		} else {
			sparse_mat<T, index_t> cond_mat(R, n_unknowns + 1);
			for (size_t r = 0; r < R; r++) {
				for (const auto& [i, v] : rows[r].coeffs) {
					cond_mat[r].push_back(i, v);
				}
				if (rows[r].rhs != T(0)) {
					cond_mat[r].push_back((index_t)n_unknowns, rows[r].rhs);
				}
				cond_mat[r].compress();
			}
			sparse_tensor<T, index_t, SPARSE_CSR> cond_csr(cond_mat);
			projection_write_tensor<T, index_t>(cond_path, std::move(cond_csr), pool);
			std::cout << "   Wrote " << cond_path.string() << std::endl;
		}
	}

	// Solve the stacked system (or handle degenerate cases without the solver)
	linear_solve_result_t<T, index_t> result;
	if (n_b_only > 0) {
		std::cout << "   System is INCONSISTENT: " << n_b_only
		          << " positions where boundary != 0 but A = 0." << std::endl;
		result.consistent = false;
		result.unique = false;
		result.n_unknowns = n_unknowns;
	} else if (R == 0) {
		std::cout << "   No constraints — trivially consistent (any c works)." << std::endl;
		result.consistent = true;
		result.unique = false;
		result.n_unknowns = n_unknowns;
		result.null_space = sparse_mat<T, index_t>(n_unknowns, n_unknowns);
		for (size_t i = 0; i < n_unknowns; i++) {
			result.null_space[i].push_back((index_t)i, (T)1);
			result.null_space[i].compress();
		}
	} else {
		sparse_tensor<T, index_t, SPARSE_COO> A_coo(std::vector<size_t>{n_unknowns, R});
		sparse_tensor<T, index_t, SPARSE_COO> b_coo(std::vector<size_t>{R});
		for (size_t r = 0; r < R; r++) {
			for (const auto& [i, v] : rows[r].coeffs) {
				std::vector<index_t> idx = {(index_t)i, (index_t)r};
				A_coo.push_back(idx, v);
			}
			if (rows[r].rhs != T(0)) {
				std::vector<index_t> idx = {(index_t)r};
				b_coo.push_back(idx, rows[r].rhs);
			}
		}
		A_coo.canonicalize();
		b_coo.canonicalize();

		auto A_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(A_coo), pool);
		auto b_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(b_coo), pool);
		if (solver == "incremental") {
			result = solve_linear_system_incremental<T, index_t>(
				std::move(A_csr), std::move(b_csr), F, opt);
		} else {
			result = solve_linear_system<T, index_t>(
				std::move(A_csr), std::move(b_csr), F, opt);
		}
	}

	// Write sol_<stem>.wxf (1 x n_unknowns)
	if (result.consistent) {
		auto sol_path = collinear_dir / ("sol_" + stem + ".wxf");
		sparse_mat<T, index_t> sol_mat(1, result.n_unknowns);
		for (size_t i = 0; i < result.solution.nnz(); i++) {
			sol_mat[0].push_back(result.solution(i), result.solution[i]);
		}
		sol_mat[0].compress();
		sparse_tensor<T, index_t, SPARSE_CSR> sol_csr(sol_mat);
		projection_write_tensor<T, index_t>(sol_path, std::move(sol_csr), pool);
		std::cout << "   Wrote " << sol_path.string() << std::endl;
	}

	// Print result
	if (result.consistent) {
		std::cout << "========================================" << std::endl;
		std::cout << "Solution" << std::endl;
		if (result.unique) {
			std::cout << "   Unique solution:" << std::endl;
		} else {
			std::cout << "   Particular solution (system underdetermined):" << std::endl;
			std::cout << "   Null space dimension: " << result.null_space.nrow << std::endl;
		}
		for (size_t i = 0; i < result.solution.nnz(); i++) {
			std::cout << "      c[" << result.solution(i) << "] = " << result.solution[i] << std::endl;
		}
		std::cout << "========================================" << std::endl;
	} else {
		std::cout << "========================================" << std::endl;
		std::cout << "No solution (system inconsistent)" << std::endl;
		std::cout << "========================================" << std::endl;
	}
}

#endif // SOLVE_COLLINEAR_HPP
