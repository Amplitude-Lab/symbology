#pragma once

// Universal projection module: computes projection matrices recursively for
// collinear and symmetry constraints on sparse tensors (symbols).
//
// Stage 1 (universal, common to collinear AND symmetry):
//   Given T^{a}_{b,c} and per-axis maps M_{b,d}, M_{c,e}:
//   1. Contract T with M_{b,d} and M_{c,e} -> T'^{a}_{d,e}
//   2. RREF T' -> reduced basis P^{a'}_{d,e}  (collinear only)
//   3. Join P with T' (collinear) or T (symmetry), extract projection P^{a'}_{a}
//
// Weight 1 is special: no Stage 1, direct contraction with the seed map.
//   Collinear:  *first@w1 = reshape(FEC_1, {7,42}) . colmat42
//   Symmetry:   *first@w1 = reshape(FEC_1, {7,42}) . S . reshape(FEC_1, {7,42})^T
//
// The projection matrix is stored as (input_dim, output_dim) so it can be
// used directly as the per-axis map in the next weight's contraction.

#include "bootstrap.hpp"
#include <algorithm>
#include <functional>
#include <regex>

// ========== File I/O helpers ==========

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> projection_read_tensor(
	const std::filesystem::path& path, const field_t& F, thread_pool* pool) {
	std::cout << "Reading file " << path.string() << " ..." << std::endl;
	auto tensor = sparse_tensor_read_wxf<T, index_t>(path, F, pool);
	print_crc32(path.string(), path);
	return tensor;
}

template <typename T, typename index_t>
void projection_write_tensor(
	const std::filesystem::path& path,
	sparse_tensor<T, index_t, SPARSE_CSR>&& tensor_csr,
	thread_pool* /*pool*/) {
	if (!path.parent_path().empty()) {
		std::filesystem::create_directories(path.parent_path());
	}
	std::cout << "Writing file " << path.string() << " ..." << std::endl;
	Timer timer;
	timer.start();
	auto u8arr = sparse_tensor_write_wxf(tensor_csr);
	std::ofstream ofs(path, std::ios::binary);
	if (!ofs) {
		throw std::runtime_error("Cannot write file: " + path.string());
	}
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs.close();
	u8arr.clear();
	u8arr.shrink_to_fit();
	timer.stop();
	std::cout << "** Write time: " << timer.milliseconds() << " ms" << std::endl;
	print_crc32(path.string(), path);
}

// ========== Ternary contraction ==========
// Contract A.axis1 with B.axis0, then result.axis1 with C.axis0.
// Result dims: (A.dim(0), B.dim(1), C.dim(1)).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> tenary_contract(
	const sparse_tensor<T, index_t, SPARSE_COO>& A,
	const sparse_tensor<T, index_t, SPARSE_COO>& B,
	const sparse_tensor<T, index_t, SPARSE_COO>& C,
	const field_t& F, thread_pool* pool) {

	std::cout << "-- Ternary contraction --" << std::endl;
	std::cout << "   A dims: ";
	for (size_t i = 0; i < A.rank(); i++) { std::cout << A.dim(i); if (i < A.rank()-1) std::cout << "x"; }
	std::cout << std::endl;
	std::cout << "   B dims: ";
	for (size_t i = 0; i < B.rank(); i++) { std::cout << B.dim(i); if (i < B.rank()-1) std::cout << "x"; }
	std::cout << std::endl;
	std::cout << "   C dims: ";
	for (size_t i = 0; i < C.rank(); i++) { std::cout << C.dim(i); if (i < C.rank()-1) std::cout << "x"; }
	std::cout << std::endl;

	Timer timer;
	timer.start();
	auto E = tensor_contract(A, B, 1, 0, F, pool);
	auto D = tensor_contract(E, C, 1, 0, F, pool);
	timer.stop();

	std::cout << "   Result dims: ";
	for (size_t i = 0; i < D.rank(); i++) { std::cout << D.dim(i); if (i < D.rank()-1) std::cout << "x"; }
	std::cout << ", nnz: " << D.nnz() << " (" << timer.milliseconds() << " ms)" << std::endl;

	return D;
}

// ========== Weight-1 projection (special case, no Stage 1) ==========
//
// Collinear: result = reshape(seed, {d0, d1*d2}) . map42
//   -> (d0, map42.dim(1))  e.g. (7, 11)
//
// Symmetry:  result = reshape(seed, {d0, d1*d2}) . map42 . reshape(seed, {d0, d1*d2})^T
//   -> (d0, d0)  e.g. (7, 7)
//
// seed is FEC_1 (7, 1, 42) or LEC_1 (14, 42, 1).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> compute_weight1_projection(
	sparse_tensor<T, index_t, SPARSE_CSR>&& seed_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& map42_csr,
	bool is_collinear,
	const field_t& F, thread_pool* pool) {

	sparse_tensor<T, index_t, SPARSE_COO> seed(std::move(seed_csr));
	sparse_tensor<T, index_t, SPARSE_COO> map42(std::move(map42_csr));

	if (seed.rank() != 3) {
		throw std::runtime_error("compute_weight1_projection: seed must be rank-3.");
	}
	if (map42.rank() != 2) {
		throw std::runtime_error("compute_weight1_projection: map42 must be rank-2.");
	}

	// Reshape seed to 2D: {dim(0), dim(1)*dim(2)}
	size_t d0 = seed.dim(0);
	size_t d1 = seed.dim(1);
	size_t d2 = seed.dim(2);
	seed.reshape({d0, d1 * d2});

	std::cout << "-- Weight-1 projection --" << std::endl;
	std::cout << "   Seed reshaped: " << d0 << "x" << d1*d2 << std::endl;
	std::cout << "   Map42: " << map42.dim(0) << "x" << map42.dim(1) << std::endl;

	if (seed.dim(1) != map42.dim(0)) {
		throw std::runtime_error("compute_weight1_projection: seed axis 1 does not match map42 axis 0.");
	}

	// C = seed . map42
	auto C = tensor_contract(seed, map42, 1, 0, F, pool);

	sparse_tensor<T, index_t, SPARSE_COO> result_coo;
	if (is_collinear) {
		result_coo = std::move(C);
		std::cout << "   Collinear w1 result: " << result_coo.dim(0) << "x" << result_coo.dim(1) << std::endl;
	} else {
		// Symmetry: result = C . seed^T = tensor_contract(C, seed, 1, 1)
		result_coo = tensor_contract(C, seed, 1, 1, F, pool);
		std::cout << "   Symmetry w1 result: " << result_coo.dim(0) << "x" << result_coo.dim(1) << std::endl;
	}

	// Convert COO -> CSR for output, matching bootstrap.hpp convention.
	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(result_coo), pool);
}

// ========== Transformation extraction (shared col_weight + is_back_sub trick) ==========
//
// Given the kernel of [D^T | partner^T] (where partner = D1 for collinear, T for symmetry),
// extract the transformation matrix that maps the original basis to the reduced basis.
//
// The kernel has shape (kernel_original_row, kernel_original_col) where
//   kernel_original_col = D.nrow + partner.nrow
//   kernel_original_row = nullity
//
// After extraction, the result has shape (kernel_original_row, partner.nrow)
// = (input_dim, output_dim), suitable for use as a per-axis map.

template <typename T, typename index_t>
sparse_mat<T, index_t> extract_transformation(
	sparse_mat<T, index_t>&& kernel,
	rref_option_t& opt) {

	// Save old options
	auto old_col_weight = opt->col_weight;
	bool old_is_back_sub = opt->is_back_sub;

	auto kernel_original_row = kernel.nrow;
	auto kernel_original_col = kernel.ncol;

	std::cout << "   Kernel before extraction: " << kernel.nrow << "x" << kernel.ncol << std::endl;

	// Force pivots into the first kernel_original_row columns (the D-block).
	// Disallow columns >= kernel_original_row (the partner-block).
	opt->col_weight = [&old_col_weight, kernel_original_row](int64_t c) -> int64_t {
		return ((c < (int64_t)kernel_original_row) ? old_col_weight(c) : (int64_t)-1);
	};
	opt->is_back_sub = true;

	auto pivots = sparse_mat_rref_reconstruct(kernel, opt, true);

	// Flatten pivots
	std::vector<pivot_t<index_t>> flatten_pivots;
	size_t rank = 0;
	for (auto& p : pivots) {
		rank += p.size();
		flatten_pivots.insert(flatten_pivots.end(), p.begin(), p.end());
	}

	if (rank != kernel_original_row) {
		throw std::runtime_error(
			"extract_transformation: kernel vectors not independent (rank="
			+ std::to_string(rank) + ", expected="
			+ std::to_string(kernel_original_row) + ")");
	}

	// Permute rows so pivots are in ascending column order
	auto perm = perm_init(kernel_original_row);
	std::sort(perm.begin(), perm.end(),
		[&](size_t a, size_t b) {
			return flatten_pivots[a].c < flatten_pivots[b].c;
		});
	for (size_t i = 0; i < kernel_original_row; i++) {
		perm[i] = flatten_pivots[perm[i]].r;
	}
	permute(perm, kernel.rows);

	// Strip the identity (pivot = first entry in each row), flip sign,
	// shift column indices to remove the D-block offset.
	for (size_t i = 0; i < kernel.nrow; i++) {
		for (size_t j = 0; j < kernel[i].nnz() - 1; j++) {
			kernel[i][j] = -kernel[i][j + 1];
			kernel[i](j) = kernel[i](j + 1) - kernel_original_row;
		}
		kernel[i].resize(kernel[i].nnz() - 1);
	}
	kernel.ncol = kernel_original_col - kernel_original_row;
	kernel.compress();

	// Restore options
	opt->col_weight = old_col_weight;
	opt->is_back_sub = old_is_back_sub;

	std::cout << "   Transformation matrix: " << kernel.nrow << "x" << kernel.ncol << std::endl;

	return std::move(kernel);
}

// ========== Stage 1 projection (weight >= 2 and SEW) ==========

template <typename T, typename index_t>
struct projection_result_t {
	sparse_tensor<T, index_t, SPARSE_CSR> projection;
	sparse_tensor<T, index_t, SPARSE_CSR> basis;  // empty for symmetry
};

template <typename T, typename index_t>
projection_result_t<T, index_t> compute_stage1_projection(
	sparse_tensor<T, index_t, SPARSE_CSR>&& T_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& map_first_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& map_second_csr,
	bool is_collinear,
	const field_t& F,
	rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);
	sparse_tensor<T, index_t, SPARSE_COO> T_coo(std::move(T_csr));
	sparse_tensor<T, index_t, SPARSE_COO> map_first(std::move(map_first_csr));
	sparse_tensor<T, index_t, SPARSE_COO> map_second(std::move(map_second_csr));

	if (T_coo.rank() != 3) {
		throw std::runtime_error("compute_stage1_projection: T must be rank-3.");
	}

	std::cout << "-- Stage 1 projection --" << std::endl;
	std::cout << "   T dims: ";
	for (size_t i = 0; i < T_coo.rank(); i++) { std::cout << T_coo.dim(i); if (i < T_coo.rank()-1) std::cout << "x"; }
	std::cout << std::endl;

	// Step 1: Contract T with per-axis maps -> D = T'^{a}_{d,e}
	auto D = tenary_contract<T, index_t>(T_coo, map_first, map_second, F, pool);

	size_t D_dim0 = D.dim(0);
	size_t D_dim1 = D.dim(1);
	size_t D_dim2 = D.dim(2);

	// Step 2: Reshape D to matrix (rows = D_dim0, cols = D_dim1*D_dim2)
	D.reshape({D_dim0, D_dim1 * D_dim2});
	auto D_mat = D.to_sparse_mat(pool);
	D.clear();

	projection_result_t<T, index_t> result;

	if (is_collinear) {
		// Collinear: RREF D -> P_basis, then join [D | P_basis]

		// RREF D to get the reduced basis
		sparse_mat<T, index_t> D1 = D_mat;  // copy (D_mat still needed for join)
		auto pivots = sparse_mat_rref_reconstruct(D1, opt);
		D1.clear_zero_row();

		std::cout << "   RREF basis: " << D1.nrow << "x" << D1.ncol;
		if (D1.nrow < D_dim0) {
			std::cout << " (degenerate from " << D_dim0 << ")";
		}
		std::cout << std::endl;

		// Store basis: reshape D1 back to rank-3 (a', d, e)
		{
			sparse_tensor<T, index_t> D1_tensor(D1);
			D1_tensor.reshape({D1.nrow, D_dim1, D_dim2});
			std::cout << "   Basis tensor: ";
			for (size_t i = 0; i < D1_tensor.rank(); i++) { std::cout << D1_tensor.dim(i); if (i < D1_tensor.rank()-1) std::cout << "x"; }
			std::cout << std::endl;
			result.basis = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(D1_tensor), pool);
		}

		// Join [D | D1], transpose, RREF, kernel
		std::cout << "   Joining [D | D1] and extracting kernel..." << std::endl;
		Timer timer;
		timer.start();
		auto joined = sparse_mat_join(D_mat, D1, pool);
		D_mat.clear();
		D1.clear();
		auto joinedT = joined.transpose();
		joined.clear();
		auto joined_pivots = sparse_mat_rref_reconstruct(joinedT, opt);
		// sparse_mat_rref_kernel returns null vectors as columns; transpose to rows
		// so each row is a null vector (matches the extract_transformation convention).
		auto kernel = sparse_mat_rref_kernel(joinedT, joined_pivots, F, opt).transpose();
		joinedT.clear();
		timer.stop();
		std::cout << "   Kernel: " << kernel.nrow << "x" << kernel.ncol
		          << " (" << timer.milliseconds() << " ms)" << std::endl;

		// Extract transformation
		kernel = extract_transformation<T, index_t>(std::move(kernel), opt);

		sparse_tensor<T, index_t> kernel_tensor(kernel);
		result.projection = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(kernel_tensor), pool);
	} else {
		// Symmetry: join [D | T], extract transformation

		// Reshape T to matrix (same column count as D since maps are square)
		size_t T_dim0 = T_coo.dim(0);
		size_t T_dim1 = T_coo.dim(1);
		size_t T_dim2 = T_coo.dim(2);
		T_coo.reshape({T_dim0, T_dim1 * T_dim2});
		auto T_mat = T_coo.to_sparse_mat(pool);
		T_coo.clear();

		if (D_mat.ncol != T_mat.ncol) {
			throw std::runtime_error("compute_stage1_projection (symmetry): D and T column counts differ ("
				+ std::to_string(D_mat.ncol) + " vs " + std::to_string(T_mat.ncol) + ")");
		}

		// Join [D | T], transpose, RREF, kernel
		std::cout << "   Joining [D | T] and extracting kernel..." << std::endl;
		Timer timer;
		timer.start();
		auto joined = sparse_mat_join(D_mat, T_mat, pool);
		D_mat.clear();
		T_mat.clear();
		auto joinedT = joined.transpose();
		joined.clear();
		auto joined_pivots = sparse_mat_rref_reconstruct(joinedT, opt);
		auto kernel = sparse_mat_rref_kernel(joinedT, joined_pivots, F, opt).transpose();
		joinedT.clear();
		timer.stop();
		std::cout << "   Kernel: " << kernel.nrow << "x" << kernel.ncol
		          << " (" << timer.milliseconds() << " ms)" << std::endl;

		// Extract transformation
		kernel = extract_transformation<T, index_t>(std::move(kernel), opt);

		sparse_tensor<T, index_t> kernel_tensor(kernel);
		result.projection = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(kernel_tensor), pool);
		// No basis for symmetry
	}

	return result;
}

// ========== Target parsing ==========

enum class target_kind_t { SEW, FEC, LEC };

struct target_t {
	std::string name;       // "SEW_5p1", "FEC_3", "LEC_2"
	target_kind_t kind;     // SEW, FEC, or LEC
	size_t fec_weight;      // FEC chain weight (for SEW: F; for FEC_W: W; for LEC: 1)
	size_t lec_weight;      // LEC chain weight (for SEW: L; for LEC_W: W; for FEC: 1)
};

// Parse target name. Recognized formats:
//   SEW_FpL  (e.g. SEW_5p1) -> kind=SEW, fec_weight=F, lec_weight=L
//   FEC_W    (e.g. FEC_3)   -> kind=FEC, fec_weight=W, lec_weight=1
//   LEC_W    (e.g. LEC_2)   -> kind=LEC, fec_weight=1, lec_weight=W
inline target_t parse_target(const std::string& name) {
	std::smatch match;
	// SEW_FpL
	std::regex re_sew("SEW_(\\d+)p(\\d+)", std::regex::icase);
	if (std::regex_match(name, match, re_sew)) {
		return {name, target_kind_t::SEW, std::stoul(match[1].str()), std::stoul(match[2].str())};
	}
	// FEC_W
	std::regex re_fec("FEC_(\\d+)", std::regex::icase);
	if (std::regex_match(name, match, re_fec)) {
		size_t w = std::stoul(match[1].str());
		return {name, target_kind_t::FEC, w, 1};
	}
	// LEC_W
	std::regex re_lec("LEC_(\\d+)", std::regex::icase);
	if (std::regex_match(name, match, re_lec)) {
		size_t w = std::stoul(match[1].str());
		return {name, target_kind_t::LEC, 1, w};
	}
	throw std::runtime_error("Invalid target name: " + name
		+ ". Expected: SEW_FpL (e.g. SEW_5p1), FEC_W (e.g. FEC_3), or LEC_W (e.g. LEC_2)");
}

// Infer the projection filename for a target.
//   SEW_FpL -> "SEW_FpL.wxf"
//   FEC_W   -> "first_wW.wxf"
//   LEC_W   -> "last_wW.wxf"
inline std::string projection_filename(const target_t& target) {
	switch (target.kind) {
		case target_kind_t::SEW:
			return target.name + ".wxf";
		case target_kind_t::FEC:
			return "first_w" + std::to_string(target.fec_weight) + ".wxf";
		case target_kind_t::LEC:
			return "last_w" + std::to_string(target.lec_weight) + ".wxf";
	}
	throw std::runtime_error("projection_filename: unknown target kind");
}

// ========== Symmetry info ==========

struct symmetry_info_t {
	std::string name;       // "collinear", "cyclic", "flip", "parity"
	std::string seed_file;  // "colmat42.wxf", "cycrepmat.wxf", etc.
	bool is_collinear;      // true for collinear (shrinks), false for symmetry
};

inline symmetry_info_t get_symmetry_info(const std::string& symmetry) {
	if (symmetry == "collinear") {
		return {"collinear", "colmat42.wxf", true};
	} else if (symmetry == "cyclic") {
		return {"cyclic", "cycrepmat.wxf", false};
	} else if (symmetry == "flip") {
		return {"flip", "fliprepmat.wxf", false};
	} else if (symmetry == "parity") {
		return {"parity", "parityrepmat.wxf", false};
	} else if (!symmetry.empty() && symmetry != "." && symmetry != ".."
	           && symmetry.find('/') == std::string::npos
	           && symmetry.find('\\') == std::string::npos) {
		// Custom letter-symmetry matrix: data/<name>.wxf (square rep matrix).
		return {symmetry, symmetry + ".wxf", false};
	}
	throw std::runtime_error("Unknown symmetry: " + symmetry
		+ ". Valid: collinear, cyclic, flip, parity, or a custom data/<name>.wxf");
}

// ========== Pipeline orchestration ==========

// Record of one output file: name, dimensions, nnz, CRC32.
struct file_record_t {
	std::string filename;
	std::vector<size_t> dims;
	size_t nnz;
	uint32_t crc32;
};

// Helper: record a file that was just written.
inline file_record_t make_file_record(
	const std::filesystem::path& path,
	const std::vector<size_t>& dims, size_t nnz) {
	return { path.filename().string(), dims, nnz, get_file_crc32(path) };
}

// Helper: write summary.txt recording all produced files.
inline void write_summary(
	const std::filesystem::path& sym_dir,
	const std::string& symmetry_name, const std::string& target_name,
	const std::vector<file_record_t>& records) {
	std::ofstream ofs(sym_dir / "summary.txt");
	if (!ofs) return;
	ofs << "# Projection summary" << std::endl;
	ofs << "# symmetry: " << symmetry_name << std::endl;
	ofs << "# target: " << target_name << std::endl;
	ofs << "#" << std::endl;
	ofs << "# filename                       dims                nnz      CRC32" << std::endl;
	for (const auto& r : records) {
		std::string dims_str;
		for (size_t i = 0; i < r.dims.size(); i++) {
			dims_str += std::to_string(r.dims[i]);
			if (i + 1 < r.dims.size()) dims_str += "x";
		}
		ofs << std::left << std::setw(32) << r.filename
		    << std::setw(20) << dims_str
		    << std::setw(9) << r.nnz
		    << std::hex << r.crc32 << std::dec << std::endl;
	}
}

template <typename T, typename index_t>
void run_projection_pipeline(
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
	std::cout << "Projection pipeline" << std::endl;
	std::cout << "  Symmetry: " << sym.name << " (is_collinear=" << sym.is_collinear << ")" << std::endl;
	std::cout << "  Target: " << target.name << " (kind="
	          << (target.kind == target_kind_t::SEW ? "SEW" :
	              target.kind == target_kind_t::FEC ? "FEC" : "LEC")
	          << ", FEC weight=" << target.fec_weight
	          << ", LEC weight=" << target.lec_weight << ")" << std::endl;
	std::cout << "  Data dir:   " << data_dir.string() << std::endl;
	std::cout << "  Output dir: " << output_dir.string() << std::endl;
	std::cout << "========================================" << std::endl;

	std::filesystem::path sym_dir = output_dir / sym.name;
	std::filesystem::create_directories(sym_dir);

	std::vector<file_record_t> records;

	// Read seed map (reused for every weight's 42-axis contraction)
	auto seed_map = projection_read_tensor<T, index_t>(data_dir / sym.seed_file, F, pool);
	std::cout << "--- Seed map ---" << std::endl;
	print_tensor_info(seed_map);

	// ===== Weight 1 (special case) =====
	std::cout << std::endl << "=== Weight 1 (seed level) ===" << std::endl;

	auto FEC_1 = projection_read_tensor<T, index_t>(data_dir / "FEC_1.wxf", F, pool);
	auto LEC_1 = projection_read_tensor<T, index_t>(data_dir / "LEC_1.wxf", F, pool);

	// *first@w1 from FEC_1
	{
		std::cout << std::endl << "--- Computing first@w1 ---" << std::endl;
		// Copy seed_map for FEC (the original is needed for LEC too)
		auto seed_copy = projection_read_tensor<T, index_t>(data_dir / sym.seed_file, F, pool);
		auto first_w1 = compute_weight1_projection<T, index_t>(
			std::move(FEC_1), std::move(seed_copy), sym.is_collinear, F, pool);
		std::cout << "--- first@w1 result ---" << std::endl;
		print_tensor_info(first_w1);
		auto dims = first_w1.dims(); auto nnz = first_w1.nnz();
		projection_write_tensor(sym_dir / "first_w1.wxf", std::move(first_w1), pool);
		records.push_back(make_file_record(sym_dir / "first_w1.wxf", dims, nnz));
	}

	// *last@w1 from LEC_1
	{
		std::cout << std::endl << "--- Computing last@w1 ---" << std::endl;
		auto last_w1 = compute_weight1_projection<T, index_t>(
			std::move(LEC_1), std::move(seed_map), sym.is_collinear, F, pool);
		std::cout << "--- last@w1 result ---" << std::endl;
		print_tensor_info(last_w1);
		auto dims = last_w1.dims(); auto nnz = last_w1.nnz();
		projection_write_tensor(sym_dir / "last_w1.wxf", std::move(last_w1), pool);
		records.push_back(make_file_record(sym_dir / "last_w1.wxf", dims, nnz));
	}

	// ===== FEC chain: weights 2..fec_weight (skip for LEC-only targets) =====
	for (size_t w = 2; w <= target.fec_weight; w++) {
		std::cout << std::endl << "=== FEC weight " << w << " ===" << std::endl;

		auto FEC_w = projection_read_tensor<T, index_t>(
			output_dir / ("FEC_" + std::to_string(w) + ".wxf"), F, pool);
		auto first_prev = projection_read_tensor<T, index_t>(
			sym_dir / ("first_w" + std::to_string(w - 1) + ".wxf"), F, pool);
		auto seed_copy = projection_read_tensor<T, index_t>(data_dir / sym.seed_file, F, pool);

		// FEC_w dims: (basis_w, basis_{w-1}, 42)
		// map_first = *first@(w-1) for axis 1, map_second = seed for axis 2
		auto result = compute_stage1_projection<T, index_t>(
			std::move(FEC_w), std::move(first_prev), std::move(seed_copy),
			sym.is_collinear, F, opt);

		std::cout << "--- first@w" << w << " ---" << std::endl;
		print_tensor_info(result.projection);
		auto p_dims = result.projection.dims(); auto p_nnz = result.projection.nnz();
		projection_write_tensor(
			sym_dir / ("first_w" + std::to_string(w) + ".wxf"),
			std::move(result.projection), pool);
		records.push_back(make_file_record(
			sym_dir / ("first_w" + std::to_string(w) + ".wxf"), p_dims, p_nnz));

		if (sym.is_collinear) {
			std::cout << "--- first@w" << w << " basis ---" << std::endl;
			print_tensor_info(result.basis);
			auto b_dims = result.basis.dims(); auto b_nnz = result.basis.nnz();
			projection_write_tensor(
				sym_dir / ("first_w" + std::to_string(w) + "_basis.wxf"),
				std::move(result.basis), pool);
			records.push_back(make_file_record(
				sym_dir / ("first_w" + std::to_string(w) + "_basis.wxf"), b_dims, b_nnz));
		}
	}

	// ===== LEC chain: weights 2..lec_weight =====
	for (size_t w = 2; w <= target.lec_weight; w++) {
		std::cout << std::endl << "=== LEC weight " << w << " ===" << std::endl;

		auto LEC_w = projection_read_tensor<T, index_t>(
			output_dir / ("LEC_" + std::to_string(w) + ".wxf"), F, pool);
		auto seed_copy = projection_read_tensor<T, index_t>(data_dir / sym.seed_file, F, pool);
		auto last_prev = projection_read_tensor<T, index_t>(
			sym_dir / ("last_w" + std::to_string(w - 1) + ".wxf"), F, pool);

		// LEC_w dims: (basis_w, 42, basis_{w-1})
		// map_first = seed for axis 1 (42), map_second = *last@(w-1) for axis 2
		auto result = compute_stage1_projection<T, index_t>(
			std::move(LEC_w), std::move(seed_copy), std::move(last_prev),
			sym.is_collinear, F, opt);

		std::cout << "--- last@w" << w << " ---" << std::endl;
		print_tensor_info(result.projection);
		auto p_dims = result.projection.dims(); auto p_nnz = result.projection.nnz();
		projection_write_tensor(
			sym_dir / ("last_w" + std::to_string(w) + ".wxf"),
			std::move(result.projection), pool);
		records.push_back(make_file_record(
			sym_dir / ("last_w" + std::to_string(w) + ".wxf"), p_dims, p_nnz));

		if (sym.is_collinear) {
			std::cout << "--- last@w" << w << " basis ---" << std::endl;
			print_tensor_info(result.basis);
			auto b_dims = result.basis.dims(); auto b_nnz = result.basis.nnz();
			projection_write_tensor(
				sym_dir / ("last_w" + std::to_string(w) + "_basis.wxf"),
				std::move(result.basis), pool);
			records.push_back(make_file_record(
				sym_dir / ("last_w" + std::to_string(w) + "_basis.wxf"), b_dims, b_nnz));
		}
	}

	// ===== SEW projection (only for SEW targets) =====
	if (target.kind == target_kind_t::SEW) {
		std::cout << std::endl << "=== SEW projection: " << target.name << " ===" << std::endl;

		auto SEW = projection_read_tensor<T, index_t>(
			output_dir / (target.name + ".wxf"), F, pool);
		auto first_final = projection_read_tensor<T, index_t>(
			sym_dir / ("first_w" + std::to_string(target.fec_weight) + ".wxf"), F, pool);
		auto last_final = projection_read_tensor<T, index_t>(
			sym_dir / ("last_w" + std::to_string(target.lec_weight) + ".wxf"), F, pool);

		// SEW dims: (sew_basis, FEC_basis, LEC_basis)
		// map_first = *first@fec_weight for axis 1, map_second = *last@lec_weight for axis 2
		auto result = compute_stage1_projection<T, index_t>(
			std::move(SEW), std::move(first_final), std::move(last_final),
			sym.is_collinear, F, opt);

		std::cout << "--- SEW projection ---" << std::endl;
		print_tensor_info(result.projection);
		auto p_dims = result.projection.dims(); auto p_nnz = result.projection.nnz();
		projection_write_tensor(
			sym_dir / (target.name + ".wxf"),
			std::move(result.projection), pool);
		records.push_back(make_file_record(
			sym_dir / (target.name + ".wxf"), p_dims, p_nnz));

		if (sym.is_collinear) {
			std::cout << "--- SEW basis ---" << std::endl;
			print_tensor_info(result.basis);
			auto b_dims = result.basis.dims(); auto b_nnz = result.basis.nnz();
			projection_write_tensor(
				sym_dir / (target.name + "_basis.wxf"),
				std::move(result.basis), pool);
			records.push_back(make_file_record(
				sym_dir / (target.name + "_basis.wxf"), b_dims, b_nnz));
		}
	}

	// Write summary
	write_summary(sym_dir, sym.name, target.name, records);
	std::cout << std::endl << "========================================" << std::endl;
	std::cout << "Pipeline complete." << std::endl;
	std::cout << "Results in: " << sym_dir.string() << std::endl;
	std::cout << "Summary: " << (sym_dir / "summary.txt").string() << std::endl;
	std::cout << "========================================" << std::endl;
}
