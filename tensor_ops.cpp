#include "numeric_parse.hpp"
#include "native_cache.hpp"
// tensor_ops — small tensor utilities built on SparseRREF:
//   ternary  T.wxf B1 B2 out.wxf          T'[a,...] = sum_{b,c} T[a,b,c] B1[b,...] B2[c,...]
//                                          ('I' = identity of the matching axis; B1/B2 may be matrices
//                                           or any rank>=2 tensor — its first axis is contracted with
//                                           the tensor axis, the remaining axes are spliced in there)
//   power    M.wxf n out.wxf               M^n (n >= 0, exact rational arithmetic)
//   join     A.wxf B.wxf axis out.wxf      Join along 1-based axis (negative counts from the end)
//   impose   T.wxf D.wxf trans out.wxf     contract S[s,i,a].D[a,i,c], SparseRREF-solve, output kernel
//   assemble out.wxf --elems s1..sn --groups g1,..,gn --coefs c1,..,cm
//                                          aligned solution space of A = sum_j cj*Bj in the common frame
//                                          Join[c_(g1)*s1, ..., 1] projected by P_A (groups with cj == 0 zeroed)
//   squeeze  T.wxf out.wxf                 drop all size-1 axes (e.g. (a,1,b) -> (a,b) matrix)
//   project  T.wxf S.wxf P.wxf out.wxf     out[d,b',c'] = sum P[d,a] T[a,b,c] S[b',b] S[c',c]
//   symsolve T.wxf Sb.wxf Sc.wxf out.wxf   3-step composite: transform T's axes 2,3 by Sb/Sc,
//                                          derive induced R (a x a) with R.T = T', K = ker(R^T - I),
//                                          out[e,b,c] = sum K[e,a] T[a,b,c] ('-' reuses Sb)

#include <chrono>
#include <cctype>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include "./SparseRREF/sparse_mat.h"
#include "./SparseRREF/sparse_tensor.h"
#include "./SparseRREF/wxf_support.h"
#include "tensor_shuffle.h"

using namespace SparseRREF;

using index_t = int32_t;
using scalar_t = rat_t;

const size_t n_of_threads = std::thread::hardware_concurrency() > 0 ?
    std::thread::hardware_concurrency() : 8;

// Write a full WXF byte buffer and fail loudly: a failed ofstream (missing
// parent dir, full disk, permission) must never leave a silently truncated or
// empty output that a later step would read as valid data.
static void write_u8(const std::filesystem::path& path, const std::vector<uint8_t>& u8arr) {
	std::ofstream ofs(path, std::ios::binary);
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs.flush();
	if (!ofs.good()) {
		auto n = std::filesystem::weakly_canonical(path).string();
		throw std::runtime_error("failed to write output file '" + n + "' (" + std::to_string(u8arr.size())
			+ " bytes; check the parent directory exists and the disk is not full)");
	}
}

static void write_tensor(const sparse_tensor<scalar_t, index_t, SPARSE_COO>& T, const std::filesystem::path& path) {
	write_u8(path, sparse_tensor_write_wxf(sparse_tensor<scalar_t, index_t, SPARSE_CSR>(T)));
}

// sparse_mat_write_wxf / sparse_tensor_write_wxf cannot encode a matrix with 0 rows;
// build the WXF for SparseArray[Automatic, {0, ncol}, 0, {1, {{0}, {}}, {}}] directly.
static std::vector<uint8_t> write_empty_sparse_tensor_wxf(const std::vector<int64_t>& dims) {
	using namespace WXF_PARSER;
	std::string_view tpl = "SparseArray[Automatic,#dims,0,{1,{#rowptr,#colindex},#vals}]";
	std::unordered_map<std::string, std::function<void(Encoder&)>> fm;
	fm["#dims"] = [&](Encoder& enc) {
		enc.push_packed_array({ dims.size() }, dims);
	};
	fm["#rowptr"] = [&](Encoder& enc) {
		enc.push_packed_array({ 1 }, std::vector<int64_t>{ 0 });
	};
	fm["#colindex"] = [&](Encoder& enc) {
		enc.push_array_info({ 0, dims.size() - 1 }, WXF_HEAD::array, 0);
		const char* empty_data = "";
		enc.push_ustr(empty_data, 0);
	};
	fm["#vals"] = [&](Encoder& enc) {
		enc.push_function("List", 0);
	};
	Encoder enc = fullform_to_wxf(tpl, fm, true);
	return enc.buffer;
}

static std::vector<uint8_t> write_empty_sparse_mat_wxf(index_t ncol) {
	return write_empty_sparse_tensor_wxf({0, ncol});
}

static sparse_tensor<scalar_t, index_t, SPARSE_COO> mat_to_tensor2(const sparse_mat<scalar_t, index_t>& M) {
	sparse_tensor<scalar_t, index_t, SPARSE_COO> t(std::vector<size_t>{ (size_t)M.nrow, (size_t)M.ncol });
	t.reserve(M.nnz());
	for (index_t i = 0; i < M.nrow; i++)
		for (size_t j = 0; j < M[i].nnz(); j++)
			t.push_back(std::vector<index_t>{ i, M[i](j) }, M[i][j]);
	return t;
}

static int run_ternary(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops ternary <tensor.wxf> <op1.wxf|I> <op2.wxf|I> <out.wxf> "
		"(operands: 'I' = identity of the matching axis, or a matrix / rank>=2 tensor .wxf)");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	if (T.rank() != 3)
		throw std::runtime_error("ternary: the first argument must be a rank-3 tensor");

	// read one operand: 'I' -> identity tensor, else tensor first (rank >= 2),
	// matrix fallback if the file is not tensor-readable. The fallback is
	// announced on stdout and a failure of BOTH interpretations rethrows with
	// both error messages — a corrupt file must never silently switch formats.
	auto read_operand = [&](const std::string& arg, size_t axis_dim) {
		if (arg == "I") {
			sparse_tensor<scalar_t, index_t, SPARSE_COO> X(std::vector<size_t>{ axis_dim, axis_dim });
			for (size_t i = 0; i < axis_dim; i++)
				X.push_back(std::vector<index_t>{ (index_t)i, (index_t)i }, (scalar_t)1);
			return X;
		}
		try {
			auto Xcsr = sparse_tensor_read_wxf<scalar_t, index_t>(base / arg, F, pool);
			sparse_tensor<scalar_t, index_t, SPARSE_COO> X(std::move(Xcsr));
			if (X.rank() < 2)
				throw std::runtime_error("ternary: operand '" + arg + "' must be rank >= 2 (got rank "
					+ std::to_string(X.rank()) + ")");
			return X;
		} catch (const std::filesystem::filesystem_error&) {
			throw;
		} catch (const std::exception& tensor_err) {
			try {
				auto X = mat_to_tensor2(sparse_mat_read_wxf<scalar_t, index_t>(base / arg, F));
				std::cout << "ternary: operand '" << arg << "' read as a plain matrix (tensor read failed: "
				          << tensor_err.what() << ")" << std::endl;
				return X;
			} catch (const std::exception& mat_err) {
				throw std::runtime_error("ternary: operand '" + arg + "' is neither a readable tensor nor a "
					"matrix (tensor error: " + tensor_err.what() + "; matrix error: " + mat_err.what() + ")");
			}
		}
	};
	auto X1 = read_operand(argv[3], T.dim(1));
	auto X2 = read_operand(argv[4], T.dim(2));
	const size_t p = X1.rank(), q = X2.rank();
	if (T.dim(1) != X1.dim(0) || T.dim(2) != X2.dim(0))
		throw std::runtime_error("ternary: dimension mismatch — need T.dim(1) == op1 axis-1 dim and T.dim(2) == op2 axis-1 dim");

	std::vector<size_t> out_dims = { T.dim(0) };
	for (size_t k = 1; k < p; k++) out_dims.push_back(X1.dim(k));
	for (size_t k = 1; k < q; k++) out_dims.push_back(X2.dim(k));

	if (T.nnz() == 0) {
		// Empty tensor (e.g. a zero-dimensional sewing basis): tensor_contract
		// cannot handle nnz == 0, so emit the empty result directly.
		sparse_tensor<scalar_t, index_t, SPARSE_COO> out(out_dims);
		write_tensor(out, base / argv[5]);
		std::cout << "ternary: dims";
		for (auto d : out_dims) std::cout << " " << d;
		std::cout << ", nnz 0 (empty input) -> " << argv[5] << std::endl;
		return 0;
	}

	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tcoo(T);
	auto step1 = tensor_contract(Tcoo, X2, 2, 0, F, pool);   // [a, b, X2tr...]
	auto step2 = tensor_contract(step1, X1, 1, 0, F, pool);  // [a, X2tr..., X1tr...]

	// permute [a, X2tr..., X1tr...] -> [a, X1tr..., X2tr...]
	std::vector<size_t> perm{ 0 };
	for (size_t k = q; k <= q + p - 2; k++) perm.push_back(k);
	for (size_t k = 1; k <= q - 1; k++) perm.push_back(k);
	auto out = step2.transpose(perm, pool);
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "ternary: dims";
	for (size_t r = 0; r < out.rank(); r++) std::cout << " " << out.dim(r);
	std::cout << ", nnz " << out.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

static int run_matmul(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 5) throw std::runtime_error("usage: tensor_ops matmul <A.wxf> <B.wxf> <out.wxf> ('I' = identity)");
	auto read_mat = [&](const char* arg) -> sparse_mat<scalar_t, index_t> {
		if (std::string(arg) == "I") throw std::runtime_error("matmul: 'I' must be the first argument to be meaningful");
		return sparse_mat_read_wxf<scalar_t, index_t>(base / arg, F);
	};
	sparse_mat<scalar_t, index_t> A;
	if (std::string(argv[2]) == "I") {
		auto B = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[3], F);
		A = sparse_mat<scalar_t, index_t>(B.nrow, B.nrow);
		for (index_t i = 0; i < A.nrow; i++) A[i].push_back(i, (scalar_t)1);
	} else {
		A = read_mat(argv[2]);
	}
	auto B = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[3], F);
	if (A.ncol != B.nrow)
		throw std::runtime_error("matmul: dimension mismatch — need A.ncol == B.nrow (got "
			+ std::to_string(A.ncol) + " vs " + std::to_string(B.nrow) + ")");
	auto C = sparse_mat_mul(A, B, F, pool);
	write_u8(base / argv[4], sparse_mat_write_wxf(C));
	std::cout << "matmul: " << A.nrow << "x" << A.ncol << " * " << B.nrow << "x" << B.ncol
	          << " = " << C.nrow << "x" << C.ncol << ", nnz " << C.nnz() << " -> " << argv[4] << std::endl;
	return 0;
}

static int run_power(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 5) throw std::runtime_error("usage: tensor_ops power <mat.wxf> <n> <out.wxf>");
	auto M = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[2], F);
	if (M.nrow != M.ncol)
		throw std::runtime_error("power: the matrix must be square");
	long long n;
	try { n = strict_integer(argv[3]); } catch (...) { throw std::runtime_error("power: n must be a non-negative integer"); }
	if (n < 0) throw std::runtime_error("power: negative powers (matrix inverse) are not supported");

	sparse_mat<scalar_t, index_t> result(M.nrow, M.ncol);
	for (index_t i = 0; i < M.nrow; i++) {
		result[i].push_back(i, (scalar_t)1);
		result[i].compress();
	}
	auto base_mat = M;
	while (n > 0) {
		if (n & 1) result = sparse_mat_mul(result, base_mat, F, pool);
		n >>= 1;
		if (n) base_mat = sparse_mat_mul(base_mat, base_mat, F, pool);
	}
	write_u8(base / argv[4], sparse_mat_write_wxf(result));
	std::cout << "power: " << M.nrow << "x" << M.ncol << " matrix to power " << argv[3]
	          << ", nnz " << result.nnz() << " -> " << argv[4] << std::endl;
	return 0;
}

// shared: build the integrability condition matrix
//   M[(a), b*d] = sum_{i} S[(a),b,i,i'] D[i',i,d]   (rank-3 S: a -> (a, d))
// The first leading axis of S is kept intact as the matrix row; any remaining
// leading axes between the first and the contracted pair are flattened
// together with the output axis d into the column index (row-major).
static sparse_mat<scalar_t, index_t> build_condition_matrix(
    const sparse_tensor<scalar_t, index_t, SPARSE_CSR>& S,
    const sparse_tensor<scalar_t, index_t, SPARSE_CSR>& D,
    const field_t& F, thread_pool* pool) {
	if (S.rank() < 2 || D.rank() != 3 || S.dim(S.rank() - 1) != D.dim(0) || S.dim(S.rank() - 2) != D.dim(1))
		throw std::runtime_error("need S[(a),b,i,a] (rank >= 2) and D[a,i,c] (rank 3) with S's last two dims == D's first two dims");
	const index_t nlead = S.rank() - 2;         // leading axes before the contracted pair
	const index_t nrow = nlead >= 1 ? (index_t)S.dim(0) : 1;
	// middle leading axes (between first and contracted pair), flattened with
	// size_t — an index_t product overflows silently for large FEC dims
	size_t ncol_inner = 1;
	for (index_t k = 1; k < nlead; k++) ncol_inner *= S.dim(k);
	const index_t dd = D.dim(2);
	const size_t ncol = ncol_inner * (size_t)dd;
	if (ncol > (size_t)std::numeric_limits<index_t>::max())
		throw std::runtime_error("icond: flattened column count " + std::to_string(ncol)
			+ " exceeds the 32-bit index limit — the condition matrix cannot be represented");
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Scoo(S), Dcoo(D);
	auto X = tensor_contract(Scoo, Dcoo, S.rank() - 1, 0, F, pool);   // [(a), b..., i, b', c]
	const index_t xi = nlead, xb = nlead + 1, xc = nlead + 2;
	// trace over the (i, b') axes; flatten: col = (b * ncol_inner + ...) * dd + c.
	// Flat per-row accumulation (append + one sort + one merge per row) — no
	// std::map in the hot path, matching the solver's own discipline; the
	// (row-major, ascending-column) layout of the previous map is preserved.
	std::vector<std::vector<std::pair<index_t, scalar_t>>> acc((size_t)nrow);
	for (size_t k = 0; k < X.nnz(); k++) {
		const index_t* idx = X.index(k);
		if (idx[xi] != idx[xb]) continue;
		index_t col = 0;
		for (index_t a = 1; a < nlead; a++) col = col * (index_t)S.dim(a) + idx[a];
		col = col * dd + idx[xc];
		index_t row = nlead >= 1 ? idx[0] : 0;
		acc[(size_t)row].emplace_back(col, X.val(k));
	}
	sparse_mat<scalar_t, index_t> M(nrow, (index_t)ncol);
	for (index_t s = 0; s < nrow; s++) {
		auto& ents = acc[(size_t)s];
		std::sort(ents.begin(), ents.end(),
		          [](const auto& a, const auto& b) { return a.first < b.first; });
		size_t w = 0;
		for (size_t r = 0; r < ents.size(); r++) {
			if (w > 0 && ents[r].first == ents[w - 1].first) {
				ents[w - 1].second = scalar_add(ents[w - 1].second, ents[r].second, F);
			} else {
				ents[w++] = ents[r];
			}
		}
		ents.resize(w);
		for (const auto& [c, v] : ents)
			if (!(v == scalar_t(0))) M[s].push_back(c, v);
		M[s].compress();
	}
	return M;
}

static int run_icond(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 5) throw std::runtime_error("usage: tensor_ops icond <tensor.wxf> <dlogmat.wxf> <out.wxf>");
	auto S = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto D = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	auto M = build_condition_matrix(S, D, F, pool);
	std::cout << "icond: condition matrix " << M.nrow << " x " << M.ncol << ", nnz " << M.nnz() << std::endl;
	write_u8(base / argv[4], sparse_mat_write_wxf(M));
	return 0;
}

static sparse_mat<scalar_t, index_t> solve_condition_kernel(
    sparse_mat<scalar_t, index_t> M, bool trans,
    const field_t& F, rref_option_t& opt) {
	if (trans) M = M.transpose();
	auto pivots = sparse_mat_rref_reconstruct(M, opt);
	auto Kcols = sparse_mat_rref_kernel(M, pivots, F, opt);
	if (Kcols.nrow == 0 && Kcols.ncol == 0) Kcols = sparse_mat<scalar_t, index_t>(M.ncol, 0);
	return Kcols.transpose();  // rows = solution basis
}

static int run_isolve(int argc, char* argv[], const field_t& F, rref_option_t& opt, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 5) throw std::runtime_error("usage: tensor_ops isolve <condmat.wxf> <trans 0|1> <out.wxf>");
	auto M = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[2], F);
	bool trans = std::string(argv[3]) == "1";
	std::cout << "isolve: solving " << M.nrow << " x " << M.ncol << ", nnz " << M.nnz()
	          << (trans ? " (transposed before solving)" : "") << std::endl;
	auto K = solve_condition_kernel(std::move(M), trans, F, opt);
	std::vector<uint8_t> u8arr;
	if (K.nrow == 0) u8arr = write_empty_sparse_mat_wxf(K.ncol);
	else u8arr = sparse_mat_write_wxf(K);
	write_u8(base / argv[4], u8arr);
	std::cout << "isolve: solution basis " << K.nrow << " x " << K.ncol << ", nnz " << K.nnz() << std::endl;
	return 0;
}

// legacy combined mode: icond + isolve in one call
static int run_impose(int argc, char* argv[], const field_t& F, rref_option_t& opt, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops impose <tensor.wxf> <dlogmat.wxf> <transpose 0|1> <out.wxf>");
	auto S = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto D = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	auto M = build_condition_matrix(S, D, F, pool);
	bool trans = std::string(argv[4]) == "1";
	std::cout << "condition matrix: " << M.nrow << " x " << M.ncol << ", nnz " << M.nnz()
	          << (trans ? " (transposed before solving)" : "") << std::endl;
	auto K = solve_condition_kernel(std::move(M), trans, F, opt);
	std::vector<uint8_t> u8arr;
	if (K.nrow == 0) u8arr = write_empty_sparse_mat_wxf(K.ncol);
	else u8arr = sparse_mat_write_wxf(K);
	write_u8(base / argv[5], u8arr);
	std::cout << "impose: solution basis " << K.nrow << " x " << K.ncol << ", nnz " << K.nnz()
	          << " -> " << argv[5] << std::endl;
	return 0;
}

static int run_join(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops join <A.wxf> <B.wxf> <axis> <out.wxf>");
	auto A = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto B = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	if (A.rank() != B.rank())
		throw std::runtime_error("join: rank mismatch");
	long long axis_in;
	try { axis_in = strict_integer(argv[4]); } catch (...) { throw std::runtime_error("join: axis must be an integer"); }
	long long rank = (long long)A.rank();
	long long ax = axis_in > 0 ? axis_in - 1 : rank + axis_in;  // 1-based; negative counts from the end
	if (ax < 0 || ax >= rank)
		throw std::runtime_error("join: axis out of range for rank " + std::to_string(rank));
	auto dims = A.dims();
	for (size_t k = 0; k < dims.size(); k++)
		if (k != (size_t)ax && dims[k] != B.dim(k))
			throw std::runtime_error("join: all dimensions except the join axis must match");
	dims[ax] += B.dim(ax);

	sparse_tensor<scalar_t, index_t, SPARSE_COO> out(dims);
	out.reserve(A.nnz() + B.nnz());
	for (size_t i = 0; i < A.nnz(); i++)
		out.push_back(A.index_vector(i), A.val(i));
	for (size_t i = 0; i < B.nnz(); i++) {
		auto idx = B.index_vector(i);
		idx[ax] += (index_t)A.dim(ax);
		out.push_back(idx, B.val(i));
	}
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "join: axis " << axis_in << ", nnz " << out.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

static std::vector<std::string> split_commas(const std::string& s) {
	std::vector<std::string> out;
	size_t pos = 0;
	while (true) {
		size_t comma = s.find(',', pos);
		std::string tok = s.substr(pos, comma == std::string::npos ? comma : comma - pos);
		if (!tok.empty()) out.push_back(tok);
		if (comma == std::string::npos) break;
		pos = comma + 1;
	}
	return out;
}

static int run_assemble(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	std::string out_arg, groups_arg, coefs_arg;
	std::vector<std::string> elem_args;
	for (int i = 2; i < argc; i++) {
		std::string a = argv[i];
		if (a == "--elems") {
			while (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) elem_args.push_back(argv[++i]);
		} else if (a == "--groups" && i + 1 < argc) {
			groups_arg = argv[++i];
		} else if (a == "--coefs" && i + 1 < argc) {
			coefs_arg = argv[++i];
		} else if (a.rfind("--", 0) != 0 && out_arg.empty()) {
			out_arg = a;
		} else {
			throw std::runtime_error("assemble: unexpected argument '" + a + "'");
		}
	}
	std::vector<long long> group_of;
	if (!groups_arg.empty()) {
		auto group_tokens = split_commas(groups_arg);
		if (group_tokens.size() != elem_args.size())
			throw std::runtime_error("assemble: --groups must give one group index per element");
		for (auto& t : group_tokens) {
			try { group_of.push_back(strict_integer(t)); } catch (...) { throw std::runtime_error("assemble: invalid group index '" + t + "'"); }
			if (group_of.back() < 1) throw std::runtime_error("assemble: group indices are 1-based positive integers");
		}
	}

	auto coef_tokens = split_commas(coefs_arg);
	long long n_coefs_expected = groups_arg.empty() ? (long long)elem_args.size()
		: [&]{ long long n = 0; for (auto g : group_of) n = std::max(n, g); return n; }();
	if ((long long)coef_tokens.size() != n_coefs_expected)
		throw std::runtime_error("assemble: --coefs must give one rational coefficient per " +
			std::string(groups_arg.empty() ? "element" : "group") + " (" +
			(n_coefs_expected == 1 ? "1 entry" : std::to_string(n_coefs_expected) + " entries") + ")");
	std::vector<scalar_t> coefs;
	for (auto& t : coef_tokens) {
		size_t i = (t.size() > 0 && t[0] == '-') ? 1 : 0;
		size_t digits = 0;
		while (i < t.size() && std::isdigit((unsigned char)t[i])) { i++; digits++; }
		bool ok = digits > 0;
		bool zero_denominator = false;
		if (ok && i < t.size() && t[i] == '/') {
			i++;
			size_t den_start = i;
			digits = 0;
			while (i < t.size() && std::isdigit((unsigned char)t[i])) { i++; digits++; }
			ok = digits > 0;
			if (ok) zero_denominator = t.substr(den_start).find_first_not_of("0") == std::string::npos;
		}
		if (!ok || i != t.size() || zero_denominator)
			throw std::runtime_error("assemble: invalid rational coefficient '" + t + "' (examples: 1, -2, 1/2)");
		coefs.push_back(scalar_t(t));
	}

	std::vector<sparse_tensor<scalar_t, index_t, SPARSE_COO>> elems;
	elems.reserve(elem_args.size());
	for (auto& e : elem_args)
		elems.emplace_back(sparse_tensor_read_wxf<scalar_t, index_t>(base / e, F, pool));

	size_t rank = elems[0].rank();
	if (rank < 1) throw std::runtime_error("assemble: elements must have rank >= 1");
	std::vector<size_t> dims = elems[0].dims();
	size_t total_rows = 0;
	for (size_t i = 0; i < elems.size(); i++) {
		if (elems[i].rank() != rank)
			throw std::runtime_error("assemble: rank mismatch between element 1 and element " + std::to_string(i + 1));
		for (size_t r = 1; r < rank; r++)
			if (elems[i].dim(r) != dims[r])
				throw std::runtime_error("assemble: trailing dimensions of element " + std::to_string(i + 1) + " do not match element 1");
		total_rows += elems[i].dim(0);
	}
	if (total_rows > (size_t)std::numeric_limits<index_t>::max())
		throw std::runtime_error("assemble: stacked first-axis dim " + std::to_string(total_rows)
			+ " exceeds the 32-bit index limit — the assembled tensor cannot be represented");
	dims[0] = total_rows;

	sparse_tensor<scalar_t, index_t, SPARSE_COO> out(dims);
	index_t offset = 0;
	std::cout << "assemble frame: first-axis dim " << total_rows << " (order as given)" << std::endl;
	for (size_t i = 0; i < elems.size(); i++) {
		const scalar_t& c = groups_arg.empty() ? coefs[i] : coefs[group_of[i] - 1];
		std::cout << "  element " << (i + 1) << " (" << elem_args[i] << "): coef " << c
		          << ", rows [" << offset << ", " << offset + (index_t)elems[i].dim(0) << ")"
		          << (c == 0 ? "  (zero rows kept)" : "") << std::endl;
		if (!(c == 0)) {
			for (size_t k = 0; k < elems[i].nnz(); k++) {
				auto idx = elems[i].index_vector(k);
				idx[0] += offset;
				out.push_back(idx, elems[i].val(k) * c);
			}
		}
		offset += (index_t)elems[i].dim(0);
	}
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / out_arg);
	std::cout << "assemble: dims";
	for (size_t r = 0; r < rank; r++) std::cout << " " << out.dim(r);
	std::cout << ", nnz " << out.nnz() << " -> " << out_arg << std::endl;
	return 0;
}

// Apply the symmetry rep matrix S (s x s) to axes 1 and 2 of T[a,b,c], then
// contract the (possibly non-square) projection map P (d x a) with axis 0:
//   out[d,b',c'] = sum_{a,b,c} P[d,a] * T[a,b,c] * S[b',b] * S[c',c]
static int run_project(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops project <tensor.wxf> <repmat.wxf> <map.wxf> <out.wxf>");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	if (T.rank() != 3)
		throw std::runtime_error("project: the tensor must be rank 3");
	auto S = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[3], F);
	if ((index_t)T.dim(1) != S.nrow || (index_t)T.dim(2) != S.nrow)
		throw std::runtime_error("project: dimension mismatch — need T.dim(2) == T.dim(3) == S.nrow (rep applied to axes 2,3)");
	auto P = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[4], F);
	if ((index_t)T.dim(0) != P.ncol)
		throw std::runtime_error("project: dimension mismatch — need T.dim(1) == P.ncol (map is d x a)");

	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tcoo(T), Scoo = mat_to_tensor2(S);
	auto step1 = tensor_contract(Tcoo, Scoo, 2, 0, F, pool);    // [a, b, c']
	auto step2 = tensor_contract(step1, Scoo, 1, 0, F, pool);   // [a, c', b']
	std::vector<size_t> perm_dims = { step2.dim(0), step2.dim(2), step2.dim(1) }; // [a, b', c']
	sparse_tensor<scalar_t, index_t, SPARSE_COO> t2(perm_dims);
	t2.reserve(step2.nnz());
	for (size_t i = 0; i < step2.nnz(); i++) {
		auto idx = step2.index_vector(i);
		t2.push_back(std::vector<index_t>{ idx[0], idx[2], idx[1] }, step2.val(i));
	}
	t2.canonicalize();
	t2.sort_indices();
	t2.reserve(t2.nnz());

	auto out = tensor_contract(mat_to_tensor2(P), t2, 1, 0, F, pool); // [d, b', c']
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "project: dims " << out.dim(0) << " " << out.dim(1) << " " << out.dim(2)
	          << ", nnz " << out.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

// Shared core for symsolve / symderive:
//   1. Impose the symmetry transformation on the input tensor's last two
//      entries: T'[a,b',c'] = sum_{b,c} T[a,b,c] Sb[b',b] Sc[c',c]
//   2. Derive the induced (a x a) transformation matrix R on the first entry
//      by solving R . T = T' via one augmented RREF on [T^T | T'^T].
// Returns R^T (a x a). Throws if no consistent R exists.
static sparse_mat<scalar_t, index_t> sym_impose_derive(
	const sparse_tensor<scalar_t, index_t, SPARSE_COO>& Tcoo, const sparse_mat<scalar_t, index_t>& Sb,
	const sparse_mat<scalar_t, index_t>& Sc, const field_t& F, rref_option_t& opt, thread_pool* pool)
{
	if (Tcoo.rank() != 3)
		throw std::runtime_error("the tensor must be rank 3");
	if ((index_t)Tcoo.dim(1) != Sb.nrow || Sb.nrow != Sb.ncol)
		throw std::runtime_error("dimension mismatch — need T.dim(2) == Sb.nrow and Sb square");
	if ((index_t)Tcoo.dim(2) != Sc.nrow || Sc.nrow != Sc.ncol)
		throw std::runtime_error("dimension mismatch — need T.dim(3) == Sc.nrow and Sc square");

	auto Tp = tensor_contract(Tcoo, mat_to_tensor2(Sc), 2, 0, F, pool);   // [a, b, c']
	Tp = tensor_contract(Tp, mat_to_tensor2(Sb), 1, 0, F, pool);         // [a, c', b']
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tpc({ Tp.dim(0), Tp.dim(2), Tp.dim(1) });
	Tpc.reserve(Tp.nnz());
	for (size_t i = 0; i < Tp.nnz(); i++) {
		auto idx = Tp.index_vector(i);
		Tpc.push_back(std::vector<index_t>{ idx[0], idx[2], idx[1] }, Tp.val(i));
	}
	Tpc.canonicalize();
	Tpc.sort_indices();
	Tpc.reserve(Tpc.nnz());

	const index_t a = (index_t)Tcoo.dim(0);
	const size_t N = Tcoo.dim(1) * Tcoo.dim(2);
	sparse_mat<scalar_t, index_t> TT((index_t)N, a);
	sparse_mat<scalar_t, index_t> TpT((index_t)N, a);
	for (size_t i = 0; i < Tcoo.nnz(); i++) {
		auto idx = Tcoo.index_vector(i);
		TT[idx[1] * Tcoo.dim(2) + idx[2]].push_back(idx[0], Tcoo.val(i));
	}
	for (size_t i = 0; i < Tpc.nnz(); i++) {
		auto idx = Tpc.index_vector(i);
		TpT[idx[1] * Tpc.dim(2) + idx[2]].push_back(idx[0], Tpc.val(i));
	}
	pool->detach_loop(0, (size_t)N, [&](size_t n) { TT[n].canonicalize(); TpT[n].canonicalize(); });
	pool->wait();
	sparse_mat<scalar_t, index_t> aug((index_t)N, 2 * a);
	for (index_t n = 0; n < (index_t)N; n++) {
		aug[n] = TT[n];
		for (size_t j = 0; j < TpT[n].nnz(); j++)
			aug[n].push_back(a + TpT[n](j), TpT[n][j]);
		aug[n].canonicalize();
	}
	// [T^T | T'^T] solves for the induced map. SparseRREF chooses sparse
	// pivots in a free column order; an appended-column pivot is not a
	// consistency witness. Restrict pivots to the unknowns, as its inverse
	// routine does, and check the reduced residual rows instead.
	auto old_weight = opt->col_weight;
	opt->col_weight = [old_weight, a](int64_t c) { return c < a ? old_weight(c) : -1; };
	std::vector<std::vector<pivot_t<index_t>>> pivots;
	try { pivots = sparse_mat_rref_reconstruct(aug, opt); }
	catch (...) { opt->col_weight = old_weight; throw; }
	opt->col_weight = old_weight;
	for (index_t r = 0; r < aug.nrow; ++r)
		if (aug[r].nnz() && aug[r](0) >= a)
			throw std::runtime_error("no transformation matrix R satisfies R.T = T' — "
			                         "the symmetry does not act on the tensor's first axis space");

	sparse_mat<scalar_t, index_t> Rt(a, a);
	for (auto& pv : pivots) for (auto& p : pv) {
		if (p.c >= a) continue;
		for (index_t e = 0; e < a; e++) {
			auto v = aug[p.r].find(a + e);
			if (v != nullptr && *v != (scalar_t)0) Rt[p.c].push_back(e, *v);
		}
	}
	for (index_t i = 0; i < a; i++) Rt[i].canonicalize();
	return Rt;  // holds R^T
}

static int run_symsolve(int argc, char* argv[], const field_t& F, rref_option_t& opt, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops symsolve <tensor.wxf> <Sb.wxf> <Sc.wxf> <out.wxf> ('-' reuses Sb for Sc, 'I' = identity of the matching axis)");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	if (T.rank() != 3) throw std::runtime_error("symmetry operations require a rank-3 tensor");
	std::string sb_arg = argv[3];
	sparse_mat<scalar_t, index_t> Sb((index_t)T.dim(1), (index_t)T.dim(1));
	if (sb_arg == "I") { for (index_t i = 0; i < Sb.nrow; i++) Sb[i].push_back(i, (scalar_t)1); }
	else Sb = sparse_mat_read_wxf<scalar_t, index_t>(base / sb_arg, F);
	std::string sc_arg = argv[4];
	sparse_mat<scalar_t, index_t> Sc;
	if (sc_arg == "-") Sc = Sb;
	else if (sc_arg == "I") { Sc = sparse_mat<scalar_t, index_t>((index_t)T.dim(2), (index_t)T.dim(2)); for (index_t i = 0; i < Sc.nrow; i++) Sc[i].push_back(i, (scalar_t)1); }
	else Sc = sparse_mat_read_wxf<scalar_t, index_t>(base / sc_arg, F);
	if (T.rank() != 3 || Sb.nrow != Sb.ncol || Sc.nrow != Sc.ncol ||
	    Sb.nrow != (index_t)T.dim(1) || Sc.nrow != (index_t)T.dim(2))
		throw std::runtime_error("symsolve: require a rank-3 tensor and square maps matching its letter axes");
	auto write_empty = [&] {
		write_u8(base / argv[5], write_empty_sparse_tensor_wxf({0, (int64_t)T.dim(1), (int64_t)T.dim(2)}));
		std::cout << "symsolve: empty invariant space, dims 0 " << T.dim(1) << " " << T.dim(2) << std::endl;
	};
	if (T.dim(0) == 0) { write_empty(); return 0; }
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tcoo(T);
	auto Rt = sym_impose_derive(Tcoo, Sb, Sc, F, opt, pool);  // holds R^T

	const index_t a = (index_t)T.dim(0);
	auto& RtI = Rt;
	for (index_t i = 0; i < RtI.nrow; i++) {
		scalar_t* entry = RtI.find(i, i);
		if (entry != nullptr) *entry = scalar_sub(*entry, (scalar_t)1, F);
		else RtI[i].push_back(i, scalar_neg((scalar_t)1, F));
	}
	pool->detach_loop(0, RtI.nrow, [&](size_t i) { RtI[i].canonicalize(); });
	pool->wait();
	auto piv2 = sparse_mat_rref_reconstruct(RtI, opt);
	auto Kcols = sparse_mat_rref_kernel(RtI, piv2, F, opt);
	auto K = Kcols.transpose();  // e x a
	if (K.nrow == 0) { write_empty(); return 0; }
	std::cout << "symsolve: derived R (" << a << "x" << a << "), invariant space " << K.nrow
	          << " x " << K.ncol << ", nnz " << K.nnz() << std::endl;

	auto out = tensor_contract(mat_to_tensor2(K), Tcoo, 1, 0, F, pool);  // [e, b, c]
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "symsolve: dims " << out.dim(0) << " " << out.dim(1) << " " << out.dim(2)
	          << ", nnz " << out.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

// Symmetry Derive — like symsolve, but writes the derived (a x a)
// transformation matrix R itself (stored transposed, R^T) instead of the
// invariant-space projection. Used to recursively derive the symmetry
// transformation of an extended solution space: from T (s_n, s_{n-1}, b) with
// maps (s_{n-1}, s_{n-1}) on entry 2 and (b, b) on entry 3, obtain R (s_n, s_n).
static int run_transpose(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 4) throw std::runtime_error("usage: tensor_ops transpose <mat.wxf> <out.wxf>");
	auto M = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[2], F);
	sparse_mat<scalar_t, index_t> Mt(M.ncol, M.nrow);
	for (index_t i = 0; i < M.nrow; i++)
		for (size_t j = 0; j < M[i].nnz(); j++)
			Mt[M[i](j)].push_back(i, M[i][j]);
	pool->detach_loop(0, Mt.nrow, [&](index_t r) { Mt[r].canonicalize(); });
	pool->wait();
	write_u8(base / argv[3], sparse_mat_write_wxf(Mt));
	std::cout << "transpose: (" << M.nrow << "x" << M.ncol << ") -> (" << Mt.nrow << "x" << Mt.ncol
	          << "), nnz " << Mt.nnz() << " -> " << argv[3] << std::endl;
	return 0;
}

static int run_symderive(int argc, char* argv[], const field_t& F, rref_option_t& opt, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops symderive <tensor.wxf> <Sb.wxf> <Sc.wxf> <out.wxf> ('-' reuses Sb for Sc, 'I' = identity of the matching axis)");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	if (T.rank() != 3) throw std::runtime_error("symmetry operations require a rank-3 tensor");
	std::string sb_arg = argv[3];
	sparse_mat<scalar_t, index_t> Sb((index_t)T.dim(1), (index_t)T.dim(1));
	if (sb_arg == "I") { for (index_t i = 0; i < Sb.nrow; i++) Sb[i].push_back(i, (scalar_t)1); }
	else Sb = sparse_mat_read_wxf<scalar_t, index_t>(base / sb_arg, F);
	std::string sc_arg = argv[4];
	sparse_mat<scalar_t, index_t> Sc;
	if (sc_arg == "-") Sc = Sb;
	else if (sc_arg == "I") { Sc = sparse_mat<scalar_t, index_t>((index_t)T.dim(2), (index_t)T.dim(2)); for (index_t i = 0; i < Sc.nrow; i++) Sc[i].push_back(i, (scalar_t)1); }
	else Sc = sparse_mat_read_wxf<scalar_t, index_t>(base / sc_arg, F);
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tcoo(T);
	auto Rt = sym_impose_derive(Tcoo, Sb, Sc, F, opt, pool);
	write_u8(base / argv[5], sparse_mat_write_wxf(Rt));
	std::cout << "symderive: wrote R^T (" << Rt.nrow << "x" << Rt.ncol << "), nnz "
	          << Rt.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

static int run_squeeze(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 4) throw std::runtime_error("usage: tensor_ops squeeze <tensor.wxf> <out.wxf>");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	std::vector<size_t> dims;
	std::vector<bool> keep(T.rank());
	for (size_t r = 0; r < T.rank(); r++) {
		keep[r] = T.dim(r) != 1;
		if (keep[r]) dims.push_back(T.dim(r));
	}
	if (dims.size() < 2)
		throw std::runtime_error("squeeze: removing the size-1 axes must leave at least rank 2");

	auto squeeze_idx = [&](const std::vector<index_t>& idx) {
		std::vector<index_t> s;
		s.reserve(dims.size());
		for (size_t r = 0; r < idx.size(); r++)
			if (keep[r]) s.push_back(idx[r]);
		return s;
	};

	if (dims.size() == 2) {
		sparse_mat<scalar_t, index_t> M(dims[0], dims[1]);
		for (size_t k = 0; k < T.nnz(); k++) {
			auto idx = squeeze_idx(T.index_vector(k));
			M[idx[0]].push_back(idx[1], T.val(k));
		}
		for (index_t i = 0; i < M.nrow; i++) M[i].compress();
		write_u8(base / argv[3], sparse_mat_write_wxf(M));
		std::cout << "squeeze: rank " << T.rank() << " -> matrix " << M.nrow << "x" << M.ncol
		          << ", nnz " << M.nnz() << " -> " << argv[3] << std::endl;
		return 0;
	}

	sparse_tensor<scalar_t, index_t, SPARSE_COO> out(dims);
	out.reserve(T.nnz());
	for (size_t k = 0; k < T.nnz(); k++)
		out.push_back(squeeze_idx(T.index_vector(k)), T.val(k));
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[3]);
	std::cout << "squeeze: rank " << T.rank() << " -> rank " << dims.size() << ", nnz " << out.nnz()
	          << " -> " << argv[3] << std::endl;
	return 0;
}

// Shuffle product: A ⊗ B with all entries scaled by a rational weight.
// Runs the SEQUENTIAL shuffle (pool = nullptr) — the parallel version has
// been observed to produce incorrect results (see compute_rhs.hpp).
// Usage: tensor_ops shuf <A.wxf> <B.wxf> <w_num/w_den> <out.wxf>
static int run_shuf(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops shuf <A.wxf> <B.wxf> <num/den (rational weight)> <out.wxf>");
	(void)pool;
	auto A = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto B = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	validate_rational(argv[4]);
	rat_t w(argv[4]);
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Acoo(A), Bcoo(B);
	auto out = tensor_shuffle_product_parallel(Acoo, Bcoo, F, nullptr);
	if (!(w == rat_t(1))) out = tensor_scalar_mul(out, w, F);
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "shuf: rank " << A.rank() << " ⊗ rank " << B.rank() << " -> rank " << out.rank()
	          << ", nnz " << out.nnz() << ", weight " << argv[4] << " -> " << argv[5] << std::endl;
	return 0;
}

// Expansion: contract the FEC axis (axis 2) of the hepMHV-style input with a
// chain of rank-3 bases, exactly mirroring expand_hepmhv / expand_tensor in
// the bootstrap: after each contraction the new FEC axis goes to position 1
// and the new letter axis to position 2 (letters in ascending weight order).
// Usage: tensor_ops expand <T.wxf> <out.wxf> <basis1.wxf> [basis2.wxf ...]
//   (bases listed highest weight first, e.g. first_w3_basis first_w2_basis)
static std::vector<size_t> expansion_perm(size_t r) {
	std::vector<size_t> perm;
	perm.push_back(0);
	perm.push_back(r - 1);
	perm.push_back(r);
	for (size_t i = 1; i + 1 <= r - 1; i++) perm.push_back(i);
	return perm;
}

static int run_expand(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc < 5) throw std::runtime_error("usage: tensor_ops expand <T.wxf> <out.wxf> <basis1.wxf> [basis2.wxf ...]");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	std::vector<std::filesystem::path> basis_paths;
	for (int i = 4; i < argc; i++) basis_paths.emplace_back(base / argv[i]);

	sparse_tensor<scalar_t, index_t, SPARSE_COO> current(T);

	// Normalize the input to the (1, FEC, letter) convention:
	//   (FEC, letter)          -> prepend dummy (1, FEC, letter)
	//   (1, FEC, letter)       -> as-is (already has the dummy solution axis)
	if (current.rank() == 2) {
		std::vector<size_t> new_dims = {1, current.dim(0), current.dim(1)};
		current.reshape(new_dims);
	} else if (current.rank() < 2 || current.dim(0) != 1) {
		throw std::runtime_error("expand: input must be rank 2 (FEC, letter) or rank 3 with a "
			"leading size-1 axis (1, FEC, letter) — the hepMHV / tdot convention");
	}

	for (size_t step = 0; step < basis_paths.size(); step++) {
		auto basis_csr = sparse_tensor_read_wxf<scalar_t, index_t>(basis_paths[step], F, pool);
		sparse_tensor<scalar_t, index_t, SPARSE_COO> basis(std::move(basis_csr));
		if (basis.rank() != 3)
			throw std::runtime_error("expand: basis " + basis_paths[step].filename().string()
				+ " must be a rank-3 tensor (FEC, FEC', letter)");
		if (current.dim(1) != basis.dim(0))
			throw std::runtime_error("expand: dimension mismatch at basis " + std::to_string(step + 1)
				+ " (" + basis_paths[step].filename().string() + "): current axis 2 has dim "
				+ std::to_string(current.dim(1)) + " != basis axis 1 dim " + std::to_string(basis.dim(0)));
		size_t r = current.rank();
		auto contracted = tensor_contract(current, basis, 1, 0, F, pool);
		current = contracted.transpose(expansion_perm(r));
	}

	// Remove the leading dummy axis: (1, 11, ..., 11) -> (11, ..., 11).
	// Fires for both input conventions — a prepended dummy (rank-2 input) and
	// a tdot-produced leading solution axis (rank-3 input (1, FEC, letter)).
	if (current.rank() >= 2 && current.dim(0) == 1) {
		std::vector<size_t> final_dims;
		for (size_t i = 1; i < current.rank(); i++) final_dims.push_back(current.dim(i));
		current.reshape(final_dims);
	}

	current.canonicalize();
	current.sort_indices();
	current.reserve(current.nnz());
	write_tensor(current, base / argv[3]);
	std::cout << "expand: rank " << T.rank() << " + " << basis_paths.size() << " bases -> rank "
	          << current.rank() << ", nnz " << current.nnz() << " -> " << argv[3] << std::endl;
	return 0;
}

// Tensor dot: contract axis i of A with axis j of B (1-based; negative counts
// from the end). Result axes = A's remaining axes (in order) followed by B's
// remaining axes (in order).
static int run_tdot(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 7) throw std::runtime_error("usage: tensor_ops tdot <A.wxf> <B.wxf> <axisA> <axisB> <out.wxf>");
	auto A = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto B = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	auto axis_of = [](const std::string& s, long long rank, const char* who) -> long long {
		long long v;
		try { v = strict_integer(s); } catch (...) { throw std::runtime_error(std::string("tdot: axis of ") + who + " must be an integer"); }
		long long ax = v > 0 ? v - 1 : rank + v;
		if (ax < 0 || ax >= rank)
			throw std::runtime_error(std::string("tdot: axis of ") + who + " out of range for rank " + std::to_string(rank));
		return ax;
	};
	long long ia = axis_of(argv[4], (long long)A.rank(), "A");
	long long ib = axis_of(argv[5], (long long)B.rank(), "B");
	if (A.dim(ia) != B.dim(ib))
		throw std::runtime_error("tdot: dimension mismatch — A.dim(axis " + std::to_string(ia + 1) + ") = "
			+ std::to_string(A.dim(ia)) + " != B.dim(axis " + std::to_string(ib + 1) + ") = " + std::to_string(B.dim(ib)));
	sparse_tensor<scalar_t, index_t, SPARSE_COO> Acoo(A), Bcoo(B);
	auto out = tensor_contract(Acoo, Bcoo, (size_t)ia, (size_t)ib, F, pool);
	if (out.rank() == 0)
		throw std::runtime_error("tdot: full contraction to a scalar (rank-0) is not supported by the WXF tensor format — "
			"keep at least one free axis (e.g. dot with a matrix instead of two vectors)");
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[6]);
	std::cout << "tdot: contract A[" << (ia + 1) << "] x B[" << (ib + 1) << "] (" << A.dim(ia) << "), result rank "
		<< out.rank() << ", nnz " << out.nnz() << " -> " << argv[6] << std::endl;
	return 0;
}

static int run_dims(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 3) throw std::runtime_error("usage: tensor_ops dims <tensor.wxf>");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	std::cout << "rank " << T.rank();
	for (size_t i = 0; i < T.rank(); i++) std::cout << " " << T.dim(i);
	std::cout << " nnz " << T.nnz() << std::endl;
	return 0;
}

int main(int argc, char* argv[]) {
	try {
		if (argc < 2) throw std::runtime_error("usage: tensor_ops <ternary|matmul|power|join|impose|icond|isolve|assemble|squeeze|project|tdot|shuf|expand|symsolve|symderive|transpose|dims> ...");
		field_t F(FIELD_QQ);
		rref_option_t opt;
		opt->pool.reset(n_of_threads);
		thread_pool* pool = &(opt->pool);
		std::filesystem::path base = native_cache::executable_path(argv[0]).parent_path();

		std::string mode = argv[1];
		auto start = std::chrono::steady_clock::now();
		int rc;
		if (mode == "ternary") rc = run_ternary(argc, argv, F, pool, base);
		else if (mode == "matmul") rc = run_matmul(argc, argv, F, pool, base);
		else if (mode == "power") rc = run_power(argc, argv, F, pool, base);
		else if (mode == "join") rc = run_join(argc, argv, F, pool, base);
		else if (mode == "impose") rc = run_impose(argc, argv, F, opt, pool, base);
		else if (mode == "icond") rc = run_icond(argc, argv, F, pool, base);
		else if (mode == "isolve") rc = run_isolve(argc, argv, F, opt, pool, base);
		else if (mode == "assemble") rc = run_assemble(argc, argv, F, pool, base);
		else if (mode == "squeeze") rc = run_squeeze(argc, argv, F, pool, base);
		else if (mode == "project") rc = run_project(argc, argv, F, pool, base);
		else if (mode == "tdot") rc = run_tdot(argc, argv, F, pool, base);
		else if (mode == "shuf") rc = run_shuf(argc, argv, F, pool, base);
		else if (mode == "expand") rc = run_expand(argc, argv, F, pool, base);
		else if (mode == "symsolve") rc = run_symsolve(argc, argv, F, opt, pool, base);
		else if (mode == "symderive") rc = run_symderive(argc, argv, F, opt, pool, base);
	else if (mode == "transpose") rc = run_transpose(argc, argv, F, pool, base);
	else if (mode == "dims") rc = run_dims(argc, argv, F, pool, base);
		else throw std::runtime_error("unknown mode '" + mode + "' (expected ternary|matmul|power|join|impose|icond|isolve|assemble|squeeze|project|tdot|shuf|expand|symsolve|symderive|transpose|dims)");
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
		std::cout << "** tensor_ops " << mode << " finished in " << ms << " ms **" << std::endl;
		return rc;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
