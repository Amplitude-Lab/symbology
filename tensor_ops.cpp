// tensor_ops — small tensor utilities built on SparseRREF:
//   ternary  T.wxf M1.wxf M2.wxf out.wxf   T'[a,b',c'] = sum_{b,c} T[a,b,c] M1[b,b'] M2[c,c']
//   power    M.wxf n out.wxf               M^n (n >= 0, exact rational arithmetic)
//   join     A.wxf B.wxf axis out.wxf      Join along 1-based axis (negative counts from the end)
//   impose   T.wxf D.wxf trans out.wxf     contract S[s,i,a].D[a,i,c], SparseRREF-solve, output kernel
//   assemble out.wxf --elems s1..sn --groups g1,..,gn --coefs c1,..,cm
//                                          aligned solution space of A = sum_j cj*Bj in the common frame
//                                          Join[c_(g1)*s1, ..., 1] projected by P_A (groups with cj == 0 zeroed)
//   squeeze  T.wxf out.wxf                 drop all size-1 axes (e.g. (a,1,b) -> (a,b) matrix)

#include <chrono>
#include <cctype>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include "./SparseRREF/sparse_mat.h"
#include "./SparseRREF/sparse_tensor.h"
#include "./SparseRREF/wxf_support.h"

using namespace SparseRREF;

using index_t = int32_t;
using scalar_t = rat_t;

const size_t n_of_threads = std::thread::hardware_concurrency() > 0 ?
    std::thread::hardware_concurrency() : 8;

static void write_tensor(const sparse_tensor<scalar_t, index_t, SPARSE_COO>& T, const std::filesystem::path& path) {
	auto u8arr = sparse_tensor_write_wxf(sparse_tensor<scalar_t, index_t, SPARSE_CSR>(T));
	std::ofstream ofs(path, std::ios::binary);
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
}

// sparse_mat_write_wxf / sparse_tensor_write_wxf cannot encode a matrix with 0 rows;
// build the WXF for SparseArray[Automatic, {0, ncol}, 0, {1, {{0}, {}}, {}}] directly.
static std::vector<uint8_t> write_empty_sparse_mat_wxf(index_t ncol) {
	using namespace WXF_PARSER;
	std::string_view tpl = "SparseArray[Automatic,#dims,0,{1,{#rowptr,#colindex},#vals}]";
	std::unordered_map<std::string, std::function<void(Encoder&)>> fm;
	fm["#dims"] = [&](Encoder& enc) {
		enc.push_packed_array({ 2 }, std::vector<int64_t>{ 0, (int64_t)ncol });
	};
	fm["#rowptr"] = [&](Encoder& enc) {
		enc.push_packed_array({ 1 }, std::vector<int64_t>{ 0 });
	};
	fm["#colindex"] = [&](Encoder& enc) {
		enc.push_array_info({ 0, 1 }, WXF_HEAD::array, 0);
		const char* empty_data = "";
		enc.push_ustr(empty_data, 0);
	};
	fm["#vals"] = [&](Encoder& enc) {
		enc.push_function("List", 0);
	};
	Encoder enc = fullform_to_wxf(tpl, fm, true);
	return enc.buffer;
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
	if (argc != 6) throw std::runtime_error("usage: tensor_ops ternary <tensor.wxf> <mat1.wxf> <mat2.wxf> <out.wxf>");
	auto T = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	if (T.rank() != 3)
		throw std::runtime_error("ternary: the first argument must be a rank-3 tensor");
	auto M1 = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[3], F);
	auto M2 = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[4], F);
	if ((index_t)T.dim(1) != M1.nrow || (index_t)T.dim(2) != M2.nrow)
		throw std::runtime_error("ternary: dimension mismatch — need T.dim(2) == M2.nrow and T.dim(1) == M1.nrow");

	sparse_tensor<scalar_t, index_t, SPARSE_COO> Tcoo(T);
	auto step1 = tensor_contract(Tcoo, mat_to_tensor2(M2), 2, 0, F, pool);   // axes [a, b, c']
	auto step2 = tensor_contract(step1, mat_to_tensor2(M1), 1, 0, F, pool);  // axes [a, c', b']

	std::vector<size_t> dims = { step2.dim(0), step2.dim(2), step2.dim(1) }; // [a, b', c']
	sparse_tensor<scalar_t, index_t, SPARSE_COO> out(dims);
	out.reserve(step2.nnz());
	for (size_t i = 0; i < step2.nnz(); i++) {
		auto idx = step2.index_vector(i);
		out.push_back(std::vector<index_t>{ idx[0], idx[2], idx[1] }, step2.val(i));
	}
	out.canonicalize();
	out.sort_indices();
	out.reserve(out.nnz());
	write_tensor(out, base / argv[5]);
	std::cout << "ternary: dims " << out.dim(0) << " " << out.dim(1) << " " << out.dim(2)
	          << ", nnz " << out.nnz() << " -> " << argv[5] << std::endl;
	return 0;
}

static int run_power(int argc, char* argv[], const field_t& F, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 5) throw std::runtime_error("usage: tensor_ops power <mat.wxf> <n> <out.wxf>");
	auto M = sparse_mat_read_wxf<scalar_t, index_t>(base / argv[2], F);
	if (M.nrow != M.ncol)
		throw std::runtime_error("power: the matrix must be square");
	long long n;
	try { n = std::stoll(argv[3]); } catch (...) { throw std::runtime_error("power: n must be a non-negative integer"); }
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
	auto u8arr = sparse_mat_write_wxf(result);
	std::ofstream ofs(base / argv[4], std::ios::binary);
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	std::cout << "power: " << M.nrow << "x" << M.ncol << " matrix to power " << argv[3]
	          << ", nnz " << result.nnz() << " -> " << argv[4] << std::endl;
	return 0;
}

static int run_impose(int argc, char* argv[], const field_t& F, rref_option_t& opt, thread_pool* pool, const std::filesystem::path& base) {
	if (argc != 6) throw std::runtime_error("usage: tensor_ops impose <tensor.wxf> <dlogmat.wxf> <transpose 0|1> <out.wxf>");
	auto S = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[2], F, pool);
	auto D = sparse_tensor_read_wxf<scalar_t, index_t>(base / argv[3], F, pool);
	bool trans = std::string(argv[4]) == "1";
	if (S.rank() != 3 || D.rank() != 3 || S.dim(2) != D.dim(0) || S.dim(1) != D.dim(1))
		throw std::runtime_error("impose: need S[s,i,a] and D[a,i,c] (both rank 3) with S.dim(2)==D.dim(0) and S.dim(1)==D.dim(1)");

	sparse_tensor<scalar_t, index_t, SPARSE_COO> Scoo(S), Dcoo(D);
	auto X = tensor_contract(Scoo, Dcoo, 2, 0, F, pool);   // [s, i, b, c]

	// trace over the (i, b) axes: M[s,c] = sum_i X[s,i,i,c]
	std::vector<std::map<index_t, scalar_t>> acc(S.dim(0));
	for (size_t k = 0; k < X.nnz(); k++) {
		auto idx = X.index_vector(k);
		if (idx[1] != idx[2]) continue;
		acc[idx[0]][idx[3]] += X.val(k);
	}
	sparse_mat<scalar_t, index_t> M(S.dim(0), D.dim(2));
	for (index_t s = 0; s < (index_t)S.dim(0); s++) {
		for (auto& [c, v] : acc[s])
			if (!(v == 0)) M[s].push_back(c, v);
		M[s].compress();
	}
	std::cout << "condition matrix: " << M.nrow << " x " << M.ncol << ", nnz " << M.nnz()
	          << (trans ? " (transposed before solving)" : "") << std::endl;

	if (trans) M = M.transpose();
	auto pivots = sparse_mat_rref_reconstruct(M, opt);
	// kernel vectors as columns; transpose to rows = solution basis
	auto Kcols = sparse_mat_rref_kernel(M, pivots, F, opt);
	if (Kcols.nrow == 0 && Kcols.ncol == 0) Kcols = sparse_mat<scalar_t, index_t>(M.ncol, 0);  // nullity 0 keeps column count
	auto K = Kcols.transpose();
	std::vector<uint8_t> u8arr;
	if (K.nrow == 0) {
		u8arr = write_empty_sparse_mat_wxf(K.ncol);
	} else {
		u8arr = sparse_mat_write_wxf(K);
	}
	std::ofstream ofs(base / argv[5], std::ios::binary);
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
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
	try { axis_in = std::stoll(argv[4]); } catch (...) { throw std::runtime_error("join: axis must be an integer"); }
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
	if (out_arg.empty() || elem_args.empty() || groups_arg.empty() || coefs_arg.empty())
		throw std::runtime_error("usage: tensor_ops assemble <out.wxf> --elems s1.wxf [s2.wxf ...] --groups g1,...,gn --coefs c1,...,cm");

	auto group_tokens = split_commas(groups_arg);
	if (group_tokens.size() != elem_args.size())
		throw std::runtime_error("assemble: --groups must give one group index per element");
	std::vector<long long> group_of;
	for (auto& t : group_tokens) {
		try { group_of.push_back(std::stoll(t)); } catch (...) { throw std::runtime_error("assemble: invalid group index '" + t + "'"); }
		if (group_of.back() < 1) throw std::runtime_error("assemble: group indices are 1-based positive integers");
	}
	long long n_groups = 0;
	for (auto g : group_of) n_groups = std::max(n_groups, g);

	auto coef_tokens = split_commas(coefs_arg);
	if ((long long)coef_tokens.size() != n_groups)
		throw std::runtime_error("assemble: --coefs must give one rational coefficient per group (groups are 1.." + std::to_string(n_groups) + ")");
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
			if (ok) zero_denominator = std::stoll(t.substr(den_start)) == 0;
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
	index_t total_rows = 0;
	for (size_t i = 0; i < elems.size(); i++) {
		if (elems[i].rank() != rank)
			throw std::runtime_error("assemble: rank mismatch between element 1 and element " + std::to_string(i + 1));
		for (size_t r = 1; r < rank; r++)
			if (elems[i].dim(r) != dims[r])
				throw std::runtime_error("assemble: trailing dimensions of element " + std::to_string(i + 1) + " do not match element 1");
		total_rows += (index_t)elems[i].dim(0);
	}
	dims[0] = (size_t)total_rows;

	sparse_tensor<scalar_t, index_t, SPARSE_COO> out(dims);
	index_t offset = 0;
	std::cout << "assemble frame: first-axis dim " << total_rows << " (order as given)" << std::endl;
	for (size_t i = 0; i < elems.size(); i++) {
		const scalar_t& c = coefs[group_of[i] - 1];
		bool included = !(c == 0);
		std::cout << "  element " << (i + 1) << " (" << elem_args[i] << "): group B" << group_of[i]
		          << ", coef " << c << ", rows [" << offset << ", " << offset + (index_t)elems[i].dim(0) << ")"
		          << (included ? "" : "  <- zeroed by P_A") << std::endl;
		if (included) {
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
		auto u8arr = sparse_mat_write_wxf(M);
		std::ofstream ofs(base / argv[3], std::ios::binary);
		ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
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

int main(int argc, char* argv[]) {
	try {
		if (argc < 2) throw std::runtime_error("usage: tensor_ops <ternary|power|join|impose|assemble|squeeze> ...");
		field_t F(FIELD_QQ);
		rref_option_t opt;
		opt->pool.reset(n_of_threads);
		thread_pool* pool = &(opt->pool);
		std::filesystem::path base = std::filesystem::path(argv[0]).parent_path();

		std::string mode = argv[1];
		auto start = std::chrono::steady_clock::now();
		int rc;
		if (mode == "ternary") rc = run_ternary(argc, argv, F, pool, base);
		else if (mode == "power") rc = run_power(argc, argv, F, pool, base);
		else if (mode == "join") rc = run_join(argc, argv, F, pool, base);
		else if (mode == "impose") rc = run_impose(argc, argv, F, opt, pool, base);
		else if (mode == "assemble") rc = run_assemble(argc, argv, F, pool, base);
		else if (mode == "squeeze") rc = run_squeeze(argc, argv, F, pool, base);
		else throw std::runtime_error("unknown mode '" + mode + "' (expected ternary|power|join|impose|assemble|squeeze)");
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
		std::cout << "** tensor_ops " << mode << " finished in " << ms << " ms **" << std::endl;
		return rc;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
