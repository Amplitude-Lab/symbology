// tensor_ops — small tensor utilities built on SparseRREF:
//   ternary  T.wxf M1.wxf M2.wxf out.wxf   T'[a,b',c'] = sum_{b,c} T[a,b,c] M1[b,b'] M2[c,c']
//   power    M.wxf n out.wxf               M^n (n >= 0, exact rational arithmetic)
//   join     A.wxf B.wxf axis out.wxf      Join along 1-based axis (negative counts from the end)

#include <chrono>
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

int main(int argc, char* argv[]) {
	try {
		if (argc < 2) throw std::runtime_error("usage: tensor_ops <ternary|power|join> ...");
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
		else throw std::runtime_error("unknown mode '" + mode + "' (expected ternary|power|join)");
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
		std::cout << "** tensor_ops " << mode << " finished in " << ms << " ms **" << std::endl;
		return rc;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
