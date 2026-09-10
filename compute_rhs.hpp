// compute_rhs.hpp — Compute the RHS (collinear boundary) for collinear constraints
//
// This is a standalone module that calculates the right-hand side of the
// collinear constraint equation c.A = b, where:
//   - A is the expanded SEW tensor after collinear projection
//   - b (the RHS / boundary) is the expected collinear limit
//
// The collinear limit is non-zero because of the normalization difference
// between BDS ansatz and BDS-like ansatz. This difference is an exponential
// that can be expanded loop order by loop order (each loop order L has
// weight 2L, so SEW_5p1 = (5+1)/2 = 3 loop order).
//
// The one-loop seed E1 is provided in data/E1.wxf (rank 2, (11,11)).
// The boundary at L loops is computed recursively:
//   L=2: E1^2 / 2
//   L=3: E1^3 / 6 + E1 * R2
//   L=4: -E1^4 / 12 + E2^2 / 2 + E1 * R3        (not yet implemented)
//   L=5: E1^5 / 20 - E1*E2^2 / 2 + E2*R3 + E1*R4  (not yet implemented)
//
// where E_L is the expanded collinear projection of hepMHV_LL, and
// R_L = E_L - boundary_L (the "remainder" at loop L).
//
// hepMHV_LL = contract(solMHV_LL, projected_SEW_basis, axis 0)
// E_L = expand(hepMHV_LL) to rank 2L, dims (11,...,11)
//
// solMHV_LL is obtained by solving c.A = boundary_L for the SEW_{(2L-1)}p1 target.
// This module invokes the solve-collinear logic for lower loop orders recursively.
//
// Naming convention for outputs:
//   output/oneloop/E1.wxf                           (copy of data/E1.wxf)
//   output/2loop/{solMHV_2L, hepMHV_2L, E2, R2, boundary_2L}.wxf
//   output/3loop/{solMHV_3L, hepMHV_3L, E3, R3, boundary_3L}.wxf
//
// SEW-level projections use naming: colprojdiv_{sew_name}.wxf (e.g. colprojdiv_SEW_5p1.wxf)
//
// IMPORTANT (Q3): We project SEW_* to its collinear expression by using the SEW
// collinear-reduced basis (SEW_{name}_basis.wxf, produced by `--project`) directly.
// We do NOT apply colprojdiv to the SEW basis. The colprojdiv operation is only
// conceptual: it ensures all parts contain a divergent component, so the divergent
// part of the collinear limit alone can fix all the parameters. In practice, the
// SEW basis already encodes the collinear expression, and colprojdiv at the SEW
// level is identity, so applying it would be a no-op (and is conceptually wrong
// per Q3: "no need to further project to divergent part").

#ifndef COMPUTE_RHS_HPP
#define COMPUTE_RHS_HPP

#include "bootstrap.hpp"
#include "projection.hpp"
#include "solve_collinear.hpp"
#include "tensor_expand.hpp"
#include "linear_solve.hpp"
#include "tensor_shuffle.h"

#include <set>
#include <map>
#include <utility>

// ========== Utility: write a CSR tensor to a WXF file ==========

template <typename T, typename index_t>
void write_tensor_file(
	const std::filesystem::path& path,
	sparse_tensor<T, index_t, SPARSE_CSR>&& tensor_csr) {

	auto u8arr = sparse_tensor_write_wxf(tensor_csr);
	std::ofstream ofs(path, std::ios::binary);
	if (!ofs) {
		throw std::runtime_error("Cannot write file: " + path.string());
	}
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs.flush();
	if (!ofs.good()) {
		throw std::runtime_error("Failed writing file (disk full or I/O error?): " + path.string());
	}
	ofs.close();
	std::cout << "   Wrote " << path.filename().string() << std::endl;
}

// ========== Shuffle power: compute E^n via repeated shuffle product ==========
//
// E^n = shuffle(shuffle(...shuffle(E, E)..., E), E)  [n times]
// The result has rank = n * E.rank() and all dims = max dim of E.

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> shuffle_power(
	const sparse_tensor<T, index_t, SPARSE_COO>& E,
	size_t n,
	const field_t& F, thread_pool* /*pool*/) {

	if (n == 0) {
		throw std::runtime_error("shuffle_power: n=0 is not supported");
	}
	if (n == 1) {
		return E;  // copy
	}

	// Use the SEQUENTIAL version of tensor_shuffle_product (pool=nullptr).
	// The parallel version has been observed to produce incorrect results for
	// the boundary calculation (E1^2/2), so we force sequential here.
	constexpr thread_pool* seq_pool = nullptr;

	// Start with E^1 = E, then repeatedly shuffle with E
	auto result = tensor_shuffle_product_parallel(E, E, F, seq_pool);  // E^2
	std::cout << "      E^2: nnz=" << result.nnz() << std::endl;
	for (size_t i = 3; i <= n; i++) {
		auto next = tensor_shuffle_product_parallel(result, E, F, seq_pool);
		std::cout << "      E^" << i << ": nnz=" << next.nnz() << std::endl;
		result = std::move(next);
	}
	return result;
}

// ========== Expand hepMHV (rank 2) to rank 2L via basis chain ==========
//
// hepMHV has shape (FEC_basis, 11) rank 2. We prepend a dummy dimension
// to make it rank 3: (1, FEC_basis, 11), then use expand_tensor which
// contracts axis 1 with each basis. The result is (1, 11, ..., 11) rank 2L+1.
// Finally, reshape to remove the dummy dimension → (11, ..., 11) rank 2L.
//
// basis_paths: highest weight first, e.g. [first_w3_basis, first_w2_basis] for L=2.

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> expand_hepmhv(
	sparse_tensor<T, index_t, SPARSE_CSR>&& hepMHV_csr,
	const std::vector<std::filesystem::path>& basis_paths,
	const field_t& F, thread_pool* pool) {

	// Convert to COO for reshaping (CSR has no reshape method).
	sparse_tensor<T, index_t, SPARSE_COO> hep(std::move(hepMHV_csr));

	// Prepend dummy dimension: (FEC, 11) → (1, FEC, 11)
	std::vector<size_t> new_dims = {1, hep.dim(0), hep.dim(1)};
	hep.reshape(new_dims);

	std::cout << "-- Expand hepMHV --" << std::endl;
	std::cout << "   Starting: rank=" << hep.rank() << " dims=1x" << hep.dim(1) << "x" << hep.dim(2)
	          << " nnz=" << hep.nnz() << std::endl;

	// Convert back to CSR for expand_tensor (which takes CSR input).
	auto hep_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(hep), pool);

	// Use expand_tensor (contracts axis 1)
	auto expanded = expand_tensor<T, index_t>(std::move(hep_csr), basis_paths, F, pool);

	// Remove the dummy dimension: (1, 11, ..., 11) → (11, ..., 11)
	// Convert to COO for reshaping.
	sparse_tensor<T, index_t, SPARSE_COO> expanded_coo(std::move(expanded));
	std::vector<size_t> final_dims;
	for (size_t i = 1; i < expanded_coo.rank(); i++) {
		final_dims.push_back(expanded_coo.dim(i));
	}
	expanded_coo.reshape(final_dims);
	auto result = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(expanded_coo), pool);

	std::cout << "   Final: rank=" << result.rank() << " dims=";
	for (size_t i = 0; i < result.rank(); i++) {
		std::cout << result.dim(i);
		if (i + 1 < result.rank()) std::cout << "x";
	}
	std::cout << " nnz=" << result.nnz() << std::endl;

	return result;
}

// ========== Compute the boundary (RHS) for loop order L ==========
//
// The master-equation recursion (identical to the flow editor's Compute RHS
// node; proven in the shuffle word algebra at every order L <= 6, and equal
// to the historical closed forms at L = 2..4):
//
//   boundary_L = (1/L) * sum_{k=1}^{L-1} k * (R_k shuffle E_{L-k}),
//   with R_1 := E_1.
//
// It replaced the hardcoded L=2..5 closed forms on 2026-09-10: the L=5 form
// was KNOWN-WRONG (t-remainder; see skills/05_compute_rhs.md) and the closed
// forms capped the loop order at 5. The recursion is any-L; the practical
// cap remains L <= 5 only because the shuffle kernel throws above total
// weight 11 (see tensor_shuffle.h).
//
// Inputs: E[1..L-1] and R[2..L-1] as COO tensors (R_1 is E_1).
// Returns the boundary as a COO tensor (rank 2L letter slots).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> compute_boundary(
	size_t L,
	const std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>>& E_list,
	const std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>>& R_list,
	const field_t& F, thread_pool* pool) {

	std::cout << "== Computing boundary for L=" << L << " (master-equation recursion) ==" << std::endl;

	using coo_t = sparse_tensor<T, index_t, SPARSE_COO>;

	if (L < 2)
		throw std::runtime_error("compute_boundary: L=" + std::to_string(L)
			+ " < 2 — nothing to compute (E1 is the seed)");
	if (E_list.find(1) == E_list.end())
		throw std::runtime_error("compute_boundary: E1 missing from the E list");
	const size_t n_letters = E_list.at(1).dim(0);

	// rank-2L letter-slot result, all slots the alphabet size
	std::vector<size_t> bdims(2 * L, n_letters);
	coo_t boundary(bdims);

	for (size_t k = 1; k <= L - 1; k++) {
		const coo_t& Rk = (k == 1) ? E_list.at(1) : R_list.at(k);   // R_1 := E1
		const coo_t& Em = E_list.at(L - k);
		std::cout << "   term k=" << k << ": (" << k << "/" << L
		          << ") * R_" << k << " shuffle E_" << (L - k) << std::flush;
		auto term = tensor_shuffle_product_parallel(Rk, Em, F, nullptr);
		std::cout << " nnz=" << term.nnz() << std::endl;
		tensor_add_weighted(boundary, term, T(1), T((int)k, (int)L), F);
	}
	std::cout << "   boundary_L = (1/" << L << ") * sum_k k*(R_k shuffle E_{L-k}): nnz="
	          << boundary.nnz() << std::endl;
	return boundary;
}

// ========== Helper: locate the dlogmat seed in data_dir ==========
//
// The condition tensor is named dlogmat_<group>.wxf (e.g. dlogmat_E6.wxf).
// For multi-project support we scan data_dir for the first match rather
// than hardcoding the group name.

inline std::filesystem::path find_dlogmat(const std::filesystem::path& data_dir) {
	if (!std::filesystem::exists(data_dir)) {
		throw std::runtime_error("find_dlogmat: data dir not found: " + data_dir.string());
	}
	std::filesystem::path found;
	for (auto& entry : std::filesystem::directory_iterator(data_dir)) {
		auto name = entry.path().filename().string();
		if (name.rfind("dlogmat_", 0) == 0 && entry.path().extension() == ".wxf") {
			found = entry.path();
			break;
		}
	}
	if (found.empty()) {
		throw std::runtime_error("find_dlogmat: no dlogmat_*.wxf found in " + data_dir.string());
	}
	return found;
}

// ========== Helper: run a bootstrap command and log it ==========
//
// All bootstrap subprocess invocations MUST receive --data-dir and --output-dir
// so that --project / --solve-symmetry / --solve-collinear auto-invocations inside
// the subprocess read from / write to the same project directories as compute_rhs.
// The paths must be absolute because the subprocess resolves relative paths
// against its own executable directory.
//
// The bootstrap binary itself is resolved as a SIBLING of the running
// executable (set from main() via compute_rhs_exe_dir()), never the
// cwd-relative "./bootstrap" — that only worked when the caller happened to
// run from the repository root. All paths are single-quoted so directories
// containing spaces survive the shell round-trip.

inline std::filesystem::path& compute_rhs_exe_dir() {
	static std::filesystem::path d = std::filesystem::current_path();
	return d;
}

inline std::string bootstrap_exe_path() {
	auto p = compute_rhs_exe_dir() / "bootstrap";
	if (std::filesystem::exists(p)) return p.string();
	return "bootstrap";  // not next to us: fall back to a PATH lookup
}

inline std::string shell_quote(const std::string& s) {
	std::string out = "'";
	for (char c : s) out += c == '\'' ? "'\\''" : std::string(1, c);
	return out + "'";
}

inline void run_bootstrap_cmd(
	const std::string& cmd,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto abs_data = std::filesystem::absolute(data_dir);
	auto abs_output = std::filesystem::absolute(output_dir);
	std::string full = cmd
		+ " --data-dir " + shell_quote(abs_data.string())
		+ " --output-dir " + shell_quote(abs_output.string());
	std::cout << "   [bootstrap] " << full << std::endl;
	int ret = std::system(full.c_str());
	if (ret != 0) {
		throw std::runtime_error("bootstrap command failed (code " + std::to_string(ret) + "): " + full);
	}
}

// ========== Ensure FEC tensors exist (bootstrap --extend chain) ==========
//
// FEC_1 is in data/. FEC_{2..F} are generated by bootstrap --extend, stored in output/.
// This generates any missing FEC tensors in the chain.

inline void ensure_fec_tensors(
	size_t F, const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto dlogmat = find_dlogmat(data_dir);

	auto prev_fec = data_dir / "FEC_1.wxf";
	if (!std::filesystem::exists(prev_fec)) {
		throw std::runtime_error("ensure_fec_tensors: FEC_1 not found: " + prev_fec.string());
	}

	for (size_t w = 2; w <= F; w++) {
		auto curr_fec = output_dir / ("FEC_" + std::to_string(w) + ".wxf");
		native_cache::receipt cache(curr_fec.string()+".native-cache", "FEC extend", {dlogmat, prev_fec, bootstrap_exe_path()}, {curr_fec});
		if (cache.valid()) {
			prev_fec = curr_fec;
			continue;
		}
		std::cout << "   Generating FEC_" << w << " via bootstrap --extend..." << std::endl;
		std::string cmd = shell_quote(bootstrap_exe_path()) + " --extend"
			+ std::string(" -c ") + shell_quote(std::filesystem::absolute(dlogmat).string())
			+ std::string(" -f ") + shell_quote(std::filesystem::absolute(prev_fec).string())
			+ std::string(" -o ") + shell_quote(std::filesystem::absolute(curr_fec).string());
		run_bootstrap_cmd(cmd, data_dir, output_dir);
		if (!std::filesystem::exists(curr_fec)) {
			throw std::runtime_error("ensure_fec_tensors: bootstrap did not produce " + curr_fec.string());
		}
		cache.save();
		prev_fec = curr_fec;
	}
}

// ========== Ensure SEW tensor and collinear basis exist ==========
//
// The collinear basis (output/collinear/SEW_{name}_basis.wxf) is generated by
// bootstrap --project --symmetry collinear --target SEW_{name}, which also
// generates all FEC collinear bases (first_w{N}_basis.wxf).
// The --project command needs the SEW tensor (output/SEW_{name}.wxf), which
// is generated by bootstrap --sew, which needs FEC_{F} (output/FEC_{F}.wxf).

inline void ensure_sew_basis(
	const std::string& sew_name, size_t F, size_t L,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto collinear_dir = output_dir / "collinear";
	auto sew_basis_path = collinear_dir / (sew_name + "_basis.wxf");

	// Validate prerequisites even when a previous basis exists.
	ensure_fec_tensors(F, data_dir, output_dir);
	auto dlogmat = find_dlogmat(data_dir);
	auto fec_path = F == 1 ? data_dir / "FEC_1.wxf" : output_dir / ("FEC_" + std::to_string(F) + ".wxf");
	auto lec_path = data_dir / ("LEC_" + std::to_string(L) + ".wxf");

	// Step 2: Ensure SEW tensor exists (via bootstrap --sew)
	auto sew_tensor_path = output_dir / (sew_name + ".wxf");
	native_cache::receipt sew_cache(sew_tensor_path.string()+".native-cache", "SEW", {dlogmat, fec_path, lec_path, bootstrap_exe_path()}, {sew_tensor_path});
	if (!sew_cache.valid()) {
		std::cout << "   Generating SEW tensor via bootstrap --sew..." << std::endl;
		if (!std::filesystem::exists(lec_path)) {
			throw std::runtime_error("ensure_sew_basis: LEC file not found: " + lec_path.string());
		}
		std::string cmd = shell_quote(bootstrap_exe_path()) + " --sew"
			+ std::string(" -c ") + shell_quote(std::filesystem::absolute(dlogmat).string())
			+ std::string(" -f ") + shell_quote(std::filesystem::absolute(fec_path).string())
			+ std::string(" -l ") + shell_quote(std::filesystem::absolute(lec_path).string())
			+ std::string(" -o ") + shell_quote(std::filesystem::absolute(sew_tensor_path).string());
		run_bootstrap_cmd(cmd, data_dir, output_dir);
		if (!std::filesystem::exists(sew_tensor_path)) {
			throw std::runtime_error("ensure_sew_basis: bootstrap --sew did not produce " + sew_tensor_path.string());
		}
		sew_cache.save();
	}
	std::vector<std::filesystem::path> project_inputs{sew_tensor_path, data_dir/"FEC_1.wxf", data_dir/"LEC_1.wxf", data_dir/"colmat42.wxf", bootstrap_exe_path()};
	std::vector<std::filesystem::path> project_outputs{collinear_dir/"first_w1.wxf", collinear_dir/"last_w1.wxf", collinear_dir/(sew_name+".wxf"), sew_basis_path};
	for (size_t w=2; w<=F; ++w) {
		project_inputs.push_back(output_dir/("FEC_"+std::to_string(w)+".wxf"));
		project_outputs.push_back(collinear_dir/("first_w"+std::to_string(w)+".wxf"));
		project_outputs.push_back(collinear_dir/("first_w"+std::to_string(w)+"_basis.wxf"));
	}
	for (size_t w=2; w<=L; ++w) {
		project_inputs.push_back(output_dir/("LEC_"+std::to_string(w)+".wxf"));
		project_outputs.push_back(collinear_dir/("last_w"+std::to_string(w)+".wxf"));
		project_outputs.push_back(collinear_dir/("last_w"+std::to_string(w)+"_basis.wxf"));
	}
	native_cache::receipt project_cache(sew_basis_path.string()+".native-cache", "SEW collinear basis", project_inputs, project_outputs);
	if (project_cache.valid()) return;

	// Step 3: Run --project to generate collinear projections and bases
	std::cout << "   Running bootstrap --project --symmetry collinear --target " << sew_name << "..." << std::endl;
	std::string cmd = shell_quote(bootstrap_exe_path()) + " --project --symmetry collinear --target " + sew_name;
	run_bootstrap_cmd(cmd, data_dir, output_dir);

	if (!std::filesystem::exists(sew_basis_path)) {
		throw std::runtime_error("ensure_sew_basis: --project did not produce " + sew_basis_path.string());
	}
	project_cache.save();
	std::cout << "   SEW collinear basis generated successfully." << std::endl;
}

// ========== Ensure FEC collinear basis files exist ==========
//
// After ensure_sew_basis runs --project, all FEC collinear bases should exist.
// This is a verification check.

inline void ensure_fec_bases(size_t target_weight, const std::filesystem::path& output_dir) {
	auto collinear_dir = output_dir / "collinear";
	for (size_t w = 2; w < target_weight; w++) {
		auto path = collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf");
		if (!std::filesystem::exists(path)) {
			throw std::runtime_error("ensure_fec_bases: " + path.string()
				+ " not found (should have been generated by --project).");
		}
	}
}

// ========== Main: compute RHS for a given loop order L ==========
//
// This function recursively computes E[1..L] and R[1..L] (R[1] = E1 - E1 = 0).
// For each loop l from 2 to L:
//   1. Compute boundary_l using E[1..l-1] and R[2..l-1]
//   2. Ensure SEW_{(2l-1)}p1_basis exists (generate if needed)
//   3. Expand SEW_basis → A (rank 2l+1, first dim = sew_dim) [no colprojdiv, per Q3]
//   4. Solve c.A = boundary_l → solMHV_lL (vector of length sew_dim)
//   5. hepMHV_lL = contract(solMHV_lL, SEW_basis, axis 0) [no colprojdiv, per Q3]
//   6. E_l = expand(hepMHV_lL) → rank 2l
//   7. R_l = E_l - boundary_l
//
// Q3: colprojdiv is NOT applied to the SEW basis. The SEW basis IS the collinear
// expression. colprojdiv is only used conceptually (and by the leaf solver
// `bootstrap --solve-collinear` when the user provides a divergent RHS).

template <typename T, typename index_t>
void compute_rhs_for_loop(
	size_t L,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const std::string& letter_projection,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);
	auto collinear_dir = output_dir / "collinear";

	std::cout << "========================================" << std::endl;
	std::cout << "Compute RHS for L=" << L << " (SEW_" << (2*L - 1) << "p1)" << std::endl;
	std::cout << "========================================" << std::endl;

	// Ensure E1 exists
	auto oneloop_dir = output_dir / "oneloop";
	std::filesystem::create_directories(oneloop_dir);
	auto E1_path = oneloop_dir / "E1.wxf";
	std::filesystem::copy_file(data_dir / "E1.wxf", E1_path, std::filesystem::copy_options::overwrite_existing);

	// E and R storage (COO tensors, indexed by loop number)
	std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>> E_list, R_list;

	// Load E1
	{
		auto E1_csr = projection_read_tensor<T, index_t>(E1_path, F, pool);
		sparse_tensor<T, index_t, SPARSE_COO> E1_coo(std::move(E1_csr));
		E_list[1] = std::move(E1_coo);
		std::cout << "E1: rank=" << E_list[1].rank() << " dims=";
		for (size_t i = 0; i < E_list[1].rank(); i++) {
			std::cout << E_list[1].dim(i) << (i + 1 < E_list[1].rank() ? "x" : "");
		}
		std::cout << " nnz=" << E_list[1].nnz() << std::endl;
	}

	// Recursively compute E[2..L-1] and R[2..L-1]
	for (size_t l = 2; l < L; l++) {
		auto l_loop_dir = output_dir / (std::to_string(l) + "loop");
		auto E_l_path = l_loop_dir / ("E" + std::to_string(l) + ".wxf");
		auto R_l_path = l_loop_dir / ("R" + std::to_string(l) + ".wxf");

		compute_rhs_for_loop<T, index_t>(l, data_dir, output_dir, letter_projection, F, opt);
		auto E_csr = projection_read_tensor<T, index_t>(E_l_path, F, pool);
		auto R_csr = projection_read_tensor<T, index_t>(R_l_path, F, pool);
		E_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(E_csr));
		R_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(R_csr));
	}

	// Now compute for loop L
	auto L_loop_dir = output_dir / (std::to_string(L) + "loop");
	std::filesystem::create_directories(L_loop_dir);

	std::string sew_name = "SEW_" + std::to_string(2*L - 1) + "p1";
	size_t F_weight = 2*L - 1;
	size_t L_weight = 1;
	size_t target_weight = F_weight + L_weight;  // = 2L

	// Step 1: Compute boundary_L
	auto boundary = compute_boundary<T, index_t>(L, E_list, R_list, F, pool);
	auto boundary_path = L_loop_dir / ("boundary_" + std::to_string(L) + "L.wxf");
	sparse_tensor<T, index_t, SPARSE_CSR> boundary_kept;
	{
		auto boundary_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(boundary), pool);
		// the file must exist for the solve subprocess; keep an in-memory copy
		// so R_L below does not re-parse it from disk
		boundary_kept = boundary_csr;
		write_tensor_file<T, index_t>(boundary_path, std::move(boundary_csr));
	}

	// Step 2: Ensure SEW basis exists
	std::cout << "== Ensuring SEW basis ==" << std::endl;
	ensure_sew_basis(sew_name, F_weight, L_weight, data_dir, output_dir);

	// Step 3: Ensure FEC bases exist
	ensure_fec_bases(target_weight, output_dir);

	// Step 4: Invoke --solve-collinear to solve c.A = boundary.
	// --solve-collinear handles: projection chain + divergent projection + matching
	// + linear solve. It writes solMHV_LL.wxf to output/<L>loop/.
	auto solMHV_path = L_loop_dir / ("solMHV_" + std::to_string(L) + "L.wxf");
	auto hepMHV_path = L_loop_dir / ("hepMHV_" + std::to_string(L) + "L.wxf");
	auto E_L_path = L_loop_dir / ("E" + std::to_string(L) + ".wxf");
	auto R_L_path = L_loop_dir / ("R" + std::to_string(L) + ".wxf");
	auto sew_basis_path = collinear_dir / (sew_name + "_basis.wxf");

	run_collinear_proj_chain<T, index_t>(detect_chain_base_paths(parse_target(sew_name), output_dir), data_dir, output_dir, F, opt, sew_name);
	std::vector<std::filesystem::path> solve_inputs{boundary_path, sew_basis_path, data_dir/"colprojdiv.wxf", data_dir/"colprojfin.wxf", bootstrap_exe_path()};
	for (size_t w=2; w<target_weight; ++w) solve_inputs.push_back(collinear_dir/("first_w"+std::to_string(w)+"_basis.wxf"));
	if (letter_projection != "identity" && letter_projection != "finite" && letter_projection != "divergent") solve_inputs.push_back(letter_projection);
	native_cache::receipt solve_cache(E_L_path.string()+".native-cache", "RHS solve " + letter_projection,
		solve_inputs, {solMHV_path, hepMHV_path, E_L_path, R_L_path});
	if (solve_cache.valid()) { std::cout << "Verified cached loop " << L << std::endl; return; }
	std::cout << "== Invoking --solve-collinear ==" << std::endl;
	{
		auto abs_data = std::filesystem::absolute(data_dir);
		auto abs_output = std::filesystem::absolute(output_dir);
		// Resolve letter_projection: sentinels ("identity"/"divergent"/"finite")
		// pass through; a relative path must be made absolute because the
		// bootstrap subprocess resolves relative paths against its own
		// executable directory.
		std::string letter_proj_arg = letter_projection;
		if (letter_proj_arg != "identity" && letter_proj_arg != "divergent" && letter_proj_arg != "finite") {
			letter_proj_arg = shell_quote(std::filesystem::absolute(letter_proj_arg).string());
		}
		std::string cmd = shell_quote(bootstrap_exe_path()) + " --solve-collinear"
			+ std::string(" --target ") + sew_name
			+ std::string(" --rhs ") + shell_quote(std::filesystem::absolute(boundary_path).string())
			+ std::string(" --projection divergent")
			+ std::string(" --letter-projection ") + letter_proj_arg
			+ std::string(" --data-dir ") + shell_quote(abs_data.string())
			+ std::string(" --output-dir ") + shell_quote(abs_output.string());
		std::cout << "   [bootstrap] " << cmd << std::endl;
		int ret = std::system(cmd.c_str());
		if (ret != 0) {
			throw std::runtime_error("compute_rhs: --solve-collinear failed (code "
				+ std::to_string(ret) + ")");
		}
	}
	if (!std::filesystem::exists(solMHV_path)) {
		throw std::runtime_error("compute_rhs: --solve-collinear did not produce " + solMHV_path.string());
	}

	// Step 5: Read solMHV_LL.wxf and compute hepMHV = contract(sol, SEW_basis)
	// Per Q3: use the SEW basis directly (no colprojdiv). The SEW basis IS the
	// collinear expression; colprojdiv is not applied to hepMHV/E_L.
	std::cout << "== Computing hepMHV (contract solMHV with SEW basis) ==" << std::endl;
	auto solMHV_csr = projection_read_tensor<T, index_t>(solMHV_path, F, pool);
	std::cout << "   solMHV_" << L << "L: rank=" << solMHV_csr.rank() << " dims=";
	for (size_t i = 0; i < solMHV_csr.rank(); i++) {
		std::cout << solMHV_csr.dim(i) << (i + 1 < solMHV_csr.rank() ? "x" : "");
	}
	std::cout << " nnz=" << solMHV_csr.nnz() << std::endl;

	auto sew_basis_csr = projection_read_tensor<T, index_t>(sew_basis_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_COO> sol_coo(std::move(solMHV_csr));
	sparse_tensor<T, index_t, SPARSE_COO> sew_coo(std::move(sew_basis_csr));

	// Contract axis 1 of sol (n_unknowns) with axis 0 of SEW basis (sew_dim).
	auto hepMHV = tensor_contract(sol_coo, sew_coo, 1, 0, F, pool);
	std::cout << "   hepMHV (raw): rank=" << hepMHV.rank() << " dims=";
	for (size_t i = 0; i < hepMHV.rank(); i++) {
		std::cout << hepMHV.dim(i) << (i + 1 < hepMHV.rank() ? "x" : "");
	}
	std::cout << " nnz=" << hepMHV.nnz() << std::endl;

	// Remove the leading "1": (1, FEC, LEC) → (FEC, LEC)
	std::vector<size_t> hepMHV_dims;
	for (size_t i = 1; i < hepMHV.rank(); i++) {
		hepMHV_dims.push_back(hepMHV.dim(i));
	}
	hepMHV.reshape(hepMHV_dims);
	auto hepMHV_csr_out = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(hepMHV), pool);

	std::cout << "   hepMHV: rank=" << hepMHV_csr_out.rank() << " dims=";
	for (size_t i = 0; i < hepMHV_csr_out.rank(); i++) {
		std::cout << hepMHV_csr_out.dim(i) << (i + 1 < hepMHV_csr_out.rank() ? "x" : "");
	}
	std::cout << " nnz=" << hepMHV_csr_out.nnz() << std::endl;
	auto hepMHV_kept = hepMHV_csr_out;   // in-memory handoff to the expansion below
	write_tensor_file<T, index_t>(hepMHV_path, std::move(hepMHV_csr_out));

	// Expansion bases for E_L (highest weight first)
	std::vector<std::filesystem::path> expansion_bases;
	for (long w = static_cast<long>(target_weight) - 1; w >= 2; w--) {
		expansion_bases.push_back(collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf"));
	}

	// Step 11: Expand hepMHV → E_L (in-memory handoff — no disk round-trip)
	std::cout << "== Expanding hepMHV to E" << L << " ==" << std::endl;
	auto E_L = expand_hepmhv<T, index_t>(std::move(hepMHV_kept), expansion_bases, F, pool);
	auto E_L_kept = E_L;   // COO copy for R_L below
	write_tensor_file<T, index_t>(E_L_path, std::move(E_L));

	// Step 12: Compute R_L = E_L - boundary (all in memory)
	std::cout << "== Computing R" << L << " ==" << std::endl;
	auto E_L_csr2 = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(E_L_kept), pool);
	auto boundary_csr2 = std::move(boundary_kept);
	sparse_tensor<T, index_t, SPARSE_COO> E_L_coo(std::move(E_L_csr2));
	sparse_tensor<T, index_t, SPARSE_COO> boundary_coo(std::move(boundary_csr2));
	tensor_add_weighted(E_L_coo, boundary_coo, T(1), T(-1), F);  // E_L = E_L - boundary
	auto R_L_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(E_L_coo), pool);
	auto R_L_kept = R_L_csr;   // in-memory handoff to the verification below
	write_tensor_file<T, index_t>(R_L_path, std::move(R_L_csr));

	// Step 13: Verify R* is free of divergent letters (indicator-vector method).
	// R_L = c·A - boundary. The collinear constraint requires the divergent part to
	// match; R_L is the finite remainder R* if and only if no divergent letter appears
	// in its support. Method: collect all distinct letter indices from R_L's nonzero
	// entries, form an 11-dim indicator vector v (v[k]=1 if letter k appears), project
	// v with the --letter-projection matrix (11, m_div). If the projection is zero →
	// R_L = R*.
	// If --letter-projection identity: this check is skipped (no divergent subspace
	// is defined — the solve was done in the full letter space and the constraint is
	// c.A = boundary at all matching positions, not just the divergent subspace).
	std::cout << "== Verifying R* is divergent-free (indicator-vector method) ==" << std::endl;

	if (letter_projection == "identity") {
		auto R_L_verify_csr_id = std::move(R_L_kept);
		std::cout << "   [SKIP] --letter-projection identity: no divergent subspace defined." << std::endl;
		std::cout << "   R" << L << " nnz=" << R_L_verify_csr_id.nnz() << " (constraint enforced at "
		          << "matching positions only)." << std::endl;
	} else if (letter_projection == "divergent" || letter_projection == "finite") {
		// Sentinel support filters: derive the divergent-letter set from
		// data/colprojdiv.wxf and check R_L's support directly (any-semantics).
		std::cout << "   [sentinel " << letter_projection << "] checking R* letter support directly" << std::endl;
		auto div_letters = load_divergent_letters<T, index_t>(data_dir, F, pool);
		sparse_tensor<T, index_t, SPARSE_COO> R_L_verify(std::move(R_L_kept));

		std::set<size_t> letter_indices;
		size_t r_rank = R_L_verify.rank();
		for (size_t i = 0; i < R_L_verify.nnz(); i++) {
			auto coords = R_L_verify.index(i);
			for (size_t d = 0; d < r_rank; d++) {
				letter_indices.insert(static_cast<size_t>(coords[d]));
			}
		}
		std::cout << "   R_L nnz=" << R_L_verify.nnz()
		          << ", distinct letter indices: " << letter_indices.size() << std::endl;
		bool has_div = false;
		for (auto l : letter_indices) {
			if (div_letters.count((index_t)l) > 0) {
				has_div = true;
				break;
			}
		}
		if (has_div) {
			std::cout << "========================================" << std::endl;
			std::cout << "FAILED: R* contains divergent letters." << std::endl;
			std::cout << "   R_L's support contains a divergent letter — the" << std::endl;
			std::cout << "   constraint was not fully matched in the selected subspace." << std::endl;
			std::cout << "========================================" << std::endl;
			throw std::runtime_error("compute_rhs: R* divergent-letter check failed at L="
				+ std::to_string(L));
		}
		std::cout << "   [OK] R* is free of divergent letters — R_L = R*" << std::endl;
	} else {
		// Load letter-projection matrix for the indicator-vector projection
		std::filesystem::path letter_proj_path(letter_projection);
		if (!std::filesystem::exists(letter_proj_path)) {
			throw std::runtime_error("compute_rhs: --letter-projection file not found: "
				+ letter_proj_path.string());
		}
		auto letter_proj_csr = projection_read_tensor<T, index_t>(letter_proj_path, F, pool);
		sparse_tensor<T, index_t, SPARSE_COO> letter_proj_coo(std::move(letter_proj_csr));
		sparse_tensor<T, index_t, SPARSE_COO> R_L_verify(std::move(R_L_kept));

		// Collect distinct letter indices from R_L's nonzero entries
		std::set<size_t> letter_indices;
		size_t r_rank = R_L_verify.rank();
		for (size_t i = 0; i < R_L_verify.nnz(); i++) {
			auto coords = R_L_verify.index(i);  // pointer to r_rank index_t values
			for (size_t d = 0; d < r_rank; d++) {
				letter_indices.insert(static_cast<size_t>(coords[d]));
			}
		}
		std::cout << "   R_L nnz=" << R_L_verify.nnz()
		          << ", distinct letter indices: " << letter_indices.size() << std::endl;
		std::cout << "   Letters appearing: {";
		{
			bool first = true;
			for (auto idx : letter_indices) {
				if (!first) std::cout << ",";
				std::cout << idx;
				first = false;
			}
		}
		std::cout << "}" << std::endl;

		// Build the indicator over the FULL alphabet, sized from the projection
		// matrix (rank-2: n_letters x m_div) — never a hardcoded 11, which
		// silently wrote out-of-range indices for any other alphabet.
		const size_t n_letters = letter_proj_coo.dim(0);
		for (auto idx : letter_indices) {
			if (idx >= n_letters)
				throw std::runtime_error("compute_rhs: R_L contains letter index "
					+ std::to_string(idx) + " outside the alphabet size "
					+ std::to_string(n_letters) + " of " + letter_proj_path.string()
					+ " — the RHS and the letter projection use different alphabets");
		}
		sparse_tensor<T, index_t, SPARSE_COO> indicator({n_letters});
		indicator.resize(letter_indices.size());
		size_t pos = 0;
		for (auto idx : letter_indices) {
			indicator.index(pos)[0] = static_cast<index_t>(idx);
			indicator.val(pos) = T(1);
			pos++;
		}
		indicator.canonicalize();
		indicator.sort_indices();
		std::cout << "   Indicator vector (" << n_letters << "-dim): nnz=" << indicator.nnz() << std::endl;

		// Project: indicator (n_letters) · letter_proj (n_letters, m_div) → (m_div,)
		auto div_check = tensor_contract(indicator, letter_proj_coo, 0, 0, F, pool);
		std::cout << "   Projection to divergent subspace: rank=" << div_check.rank()
		          << " dims=";
		for (size_t i = 0; i < div_check.rank(); i++) {
			std::cout << div_check.dim(i) << (i + 1 < div_check.rank() ? "x" : "");
		}
		std::cout << " nnz=" << div_check.nnz() << std::endl;

		if (div_check.nnz() != 0) {
			std::cout << "========================================" << std::endl;
			std::cout << "FAILED: R* contains divergent letters." << std::endl;
			std::cout << "   The indicator vector projected to a non-zero vector in the" << std::endl;
			std::cout << "   divergent subspace — R_L is NOT free of divergent letters." << std::endl;
			std::cout << "========================================" << std::endl;
			throw std::runtime_error("compute_rhs: R* divergent-letter check failed at L="
				+ std::to_string(L));
		}
		std::cout << "   [OK] R* is free of divergent letters — R_L = R*" << std::endl;
	}

	std::cout << "========================================" << std::endl;
	std::cout << "Done: L=" << L << std::endl;
	std::cout << "   E" << L << ": " << (L_loop_dir / ("E" + std::to_string(L) + ".wxf")).string() << std::endl;
	std::cout << "   R" << L << ": " << (L_loop_dir / ("R" + std::to_string(L) + ".wxf")).string() << std::endl;
	std::cout << "   boundary_" << L << "L: " << boundary_path.string() << std::endl;
	std::cout << "   solMHV_" << L << "L: " << (L_loop_dir / ("solMHV_" + std::to_string(L) + "L.wxf")).string() << std::endl;
	std::cout << "   hepMHV_" << L << "L: " << (L_loop_dir / ("hepMHV_" + std::to_string(L) + "L.wxf")).string() << std::endl;
	std::cout << "========================================" << std::endl;
	solve_cache.save();
}

#endif // COMPUTE_RHS_HPP
