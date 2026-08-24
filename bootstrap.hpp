#pragma once

#include <array>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#define USE_MIMALLOC 1
#include "SparseRREF/sparse_mat.h"
#include "SparseRREF/sparse_tensor.h"
#include "SparseRREF/wxf_support.h"

using namespace SparseRREF;

template <typename T, typename index_t, SPARSE_TYPE Type>
void print_tensor_info(const sparse_tensor<T, index_t, Type>& S) {
	std::cout << "rank: " << S.rank() << std::endl;
	std::cout << "dims: ";
	for (size_t i = 0; i < S.rank(); i++)
		std::cout << S.dim(i) << " ";
	std::cout << std::endl;
	std::cout << "nnz: " << S.nnz() << std::endl;
	std::cout << "alloc: " << S.alloc() << std::endl;
}

template <typename T, typename index_t>
void print_tensor_info(const sparse_mat<T, index_t>& S) {
	std::cout << "rank: " << 2 << std::endl;
	std::cout << "dims: ";
	std::cout << S.nrow << " ";
	std::cout << S.ncol << " ";
	std::cout << std::endl;
	std::cout << "nnz: " << S.nnz() << std::endl;
	std::cout << "alloc: " << S.alloc() << std::endl;
}

inline uint32_t get_file_crc32(const std::filesystem::path& filename) {
	std::ifstream file(filename, std::ios::binary);
	if (!file) throw std::runtime_error("Cannot open file: " + filename.string());

	constexpr auto crc_table = [] {
		std::array<uint32_t, 256> table{};
		for (uint32_t i = 0; i < 256; ++i) {
			uint32_t c = i;
			for (size_t j = 0; j < 8; ++j) {
				c = (c & 1) ? (0xEDB88320 ^ (c >> 1)) : (c >> 1);
			}
			table[i] = c;
		}
		return table;
	}();

	return std::accumulate(
		std::istreambuf_iterator<char>(file),
		std::istreambuf_iterator<char>(),
		uint32_t(0xFFFFFFFF),
		[&](uint32_t crc, char c) {
			return crc_table[(crc ^ c) & 0xFF] ^ (crc >> 8);
		}
	) ^ 0xFFFFFFFF;
}

inline void print_crc32(const std::string& label, const std::filesystem::path& path) {
	auto crc = get_file_crc32(path);
	std::cout << "CRC32 of " << label << " : " << std::hex << crc << std::dec << std::endl;
}

template <typename index_t>
void vec_cancel_divisor(sparse_vec<rat_t, index_t>& vec) {
	// RREF over QQ leaves rational rows; archive bootstrap output stores
	// primitive integer rows after clearing common denominators.
	int_t den = 1;
	for (size_t i = 0; i < vec.nnz(); i++) {
		den = LCM(den, vec.entries[i].den());
	}
	sparse_vec_rescale(vec, (rat_t)den, field_t(FIELD_QQ));
	int_t gcd = 1;
	for (size_t i = 0; i < vec.nnz(); i++) {
		auto num = vec.entries[i].num();
		if (num != 0) {
			gcd = GCD(gcd, num);
		}
	}
	sparse_vec_rescale(vec, rat_t(1, gcd), field_t(FIELD_QQ));
}

// Small reusable stopwatch. milliseconds()/seconds() also work while running,
// which is useful for long RREF stages if progress reporting is added later.
class Timer {
public:
	void start() {
		start_time = std::chrono::steady_clock::now();
		stop_time = start_time;
		stopped = false;
	}

	void stop() {
		stop_time = std::chrono::steady_clock::now();
		stopped = true;
	}

	int64_t milliseconds() const {
		auto end_time = stopped ? stop_time : std::chrono::steady_clock::now();
		return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	}

	double seconds() const {
		auto end_time = stopped ? stop_time : std::chrono::steady_clock::now();
		return std::chrono::duration<double>(end_time - start_time).count();
	}

private:
	std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point stop_time = start_time;
	bool stopped = true;
};

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> extend_forward(
	sparse_tensor<T, index_t, SPARSE_CSR>&& dlogmat_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& FEC_csr,
	const field_t& F,
	rref_option_t& opt) {

	static_assert(std::is_same_v<T, rat_t>, "v1 bootstrap only supports rat_t.");

	// FEC layout is {new_basis, old_basis, letter}. Contracting the letter
	// with dlogmat gives the integrability equations for the next weight.
	thread_pool* pool = &(opt->pool);
	sparse_tensor<T, index_t, SPARSE_COO> dlogmat(std::move(dlogmat_csr));
	sparse_tensor<T, index_t, SPARSE_COO> FEC(std::move(FEC_csr));

	if (dlogmat.rank() != 3) {
		throw std::runtime_error("dlogmat must be a rank-3 sparse tensor.");
	}
	if (FEC.rank() != 3) {
		throw std::runtime_error("FEC must be a rank-3 sparse tensor.");
	}
	if (FEC.dim(2) != dlogmat.dim(0)) {
		throw std::runtime_error("FEC letter dimension does not match dlogmat first letter dimension.");
	}

	std::cout << "--- FEC ---" << std::endl;
	print_tensor_info(FEC);
	std::cout << "--- dlogmat ---" << std::endl;
	print_tensor_info(dlogmat);

	size_t alphabet_size = dlogmat.dim(1);
	size_t old_size = FEC.dim(0);

	Timer timer;
	timer.start();
	auto tmp = tensor_contract(FEC, dlogmat, 2, 0, F, pool);
	tmp.flatten({{1, 3}, {0, 2}});

	std::cout << "--- tmp ---" << std::endl;
	print_tensor_info(tmp);

	sparse_mat<T, index_t> FEC_dlogmat = tmp.to_sparse_mat(pool);
	tmp.clear();

	std::cout << "--- FEC.dlogmat ---" << std::endl;
	print_tensor_info(FEC_dlogmat);
	timer.stop();
	std::cout << "** FEC.dlogmat time: " << timer.milliseconds() << " ms" << std::endl;

	// Null vectors of the flattened constraint matrix are the next FEC rows.
	timer.start();
	FEC_dlogmat.sort_rows_by_nnz();
	sparse_mat<T, index_t> kernel;
	auto piv = sparse_mat_rref_reconstruct(FEC_dlogmat, opt);
	std::cout << "--- FEC.dlogmat RREF ---" << std::endl;
	print_tensor_info(FEC_dlogmat);
	kernel = sparse_mat_rref_kernel(FEC_dlogmat, piv, F, opt).transpose();
	timer.stop();
	std::cout << "** RREF time: " << timer.milliseconds() << " ms" << std::endl;

	timer.start();
	std::cout << "-- Post-processing (cancel common divisor and sort) ..." << std::endl;
	pool->detach_loop(0, kernel.nrow, [&](size_t i) {
		kernel[i].canonicalize();
		vec_cancel_divisor(kernel[i]);
	});
	pool->wait();
	kernel.sort_rows_by_nnz();
	timer.stop();
	std::cout << "** Post-process time: " << timer.milliseconds() << " ms" << std::endl;

	sparse_tensor<T, index_t, SPARSE_COO> kernel_tensor(kernel, pool);
	kernel_tensor.reshape({kernel.nrow, old_size, alphabet_size});

	std::cout << "--- kernel ---" << std::endl;
	print_tensor_info(kernel_tensor);

	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(kernel_tensor), pool);
}

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> extend_backward(
	sparse_tensor<T, index_t, SPARSE_CSR>&& dlogmat_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& LEC_csr,
	const field_t& F,
	rref_option_t& opt) {

	static_assert(std::is_same_v<T, rat_t>, "v1 bootstrap only supports rat_t.");

	// LEC layout is {new_basis, letter, old_basis}; this is the natural
	// backward-growth convention, not the old archive LCC axis order.
	thread_pool* pool = &(opt->pool);
	sparse_tensor<T, index_t, SPARSE_COO> dlogmat(std::move(dlogmat_csr));
	sparse_tensor<T, index_t, SPARSE_COO> LEC(std::move(LEC_csr));

	if (dlogmat.rank() != 3) {
		throw std::runtime_error("dlogmat must be a rank-3 sparse tensor.");
	}
	if (LEC.rank() != 3) {
		throw std::runtime_error("LEC must be a rank-3 sparse tensor.");
	}
	if (LEC.dim(1) != dlogmat.dim(1)) {
		throw std::runtime_error("LEC letter dimension does not match dlogmat second letter dimension.");
	}

	std::cout << "--- LEC ---" << std::endl;
	print_tensor_info(LEC);
	std::cout << "--- dlogmat ---" << std::endl;
	print_tensor_info(dlogmat);

	size_t alphabet_size = dlogmat.dim(0);
	size_t old_size = LEC.dim(0);

	Timer timer;
	timer.start();
	auto tmp = tensor_contract(dlogmat, LEC, 1, 1, F, pool);
	tmp.flatten({{1, 3}, {0, 2}});

	std::cout << "--- tmp ---" << std::endl;
	print_tensor_info(tmp);

	sparse_mat<T, index_t> LEC_dlogmat = tmp.to_sparse_mat(pool);
	tmp.clear();

	std::cout << "--- LEC.dlogmat ---" << std::endl;
	print_tensor_info(LEC_dlogmat);
	timer.stop();
	std::cout << "** LEC.dlogmat time: " << timer.milliseconds() << " ms" << std::endl;

	// The kernel is reshaped back into the same LEC axis order.
	timer.start();
	LEC_dlogmat.sort_rows_by_nnz();
	sparse_mat<T, index_t> kernel;
	auto piv = sparse_mat_rref_reconstruct(LEC_dlogmat, opt);
	std::cout << "--- LEC.dlogmat RREF ---" << std::endl;
	print_tensor_info(LEC_dlogmat);
	kernel = sparse_mat_rref_kernel(LEC_dlogmat, piv, F, opt).transpose();
	timer.stop();
	std::cout << "** RREF time: " << timer.milliseconds() << " ms" << std::endl;

	timer.start();
	std::cout << "-- Post-processing (cancel common divisor and sort) ..." << std::endl;
	pool->detach_loop(0, kernel.nrow, [&](size_t i) {
		kernel[i].canonicalize();
		vec_cancel_divisor(kernel[i]);
	});
	pool->wait();
	kernel.sort_rows_by_nnz();
	timer.stop();
	std::cout << "** Post-process time: " << timer.milliseconds() << " ms" << std::endl;

	sparse_tensor<T, index_t, SPARSE_COO> kernel_tensor(kernel, pool);
	kernel_tensor.reshape({kernel.nrow, alphabet_size, old_size});

	std::cout << "--- kernel ---" << std::endl;
	print_tensor_info(kernel_tensor);

	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(kernel_tensor), pool);
}

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> sew_first_last(
	sparse_tensor<T, index_t, SPARSE_CSR>&& dlogmat_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& FEC_csr,
	sparse_tensor<T, index_t, SPARSE_CSR>&& LEC_csr,
	const field_t& F,
	rref_option_t& opt) {

	static_assert(std::is_same_v<T, rat_t>, "v1 bootstrap only supports rat_t.");

	// Sewing first reduces FEC.dlogmat to an independent row basis, then
	// contracts that basis with LEC to impose first/last consistency.
	thread_pool* pool = &(opt->pool);
	sparse_tensor<T, index_t, SPARSE_COO> dlogmat(std::move(dlogmat_csr));
	sparse_tensor<T, index_t, SPARSE_COO> FEC(std::move(FEC_csr));
	sparse_tensor<T, index_t, SPARSE_COO> LEC(std::move(LEC_csr));

	if (dlogmat.rank() != 3) {
		throw std::runtime_error("dlogmat must be a rank-3 sparse tensor.");
	}
	if (FEC.rank() != 3) {
		throw std::runtime_error("FEC must be a rank-3 sparse tensor.");
	}
	if (LEC.rank() != 3) {
		throw std::runtime_error("LEC must be a rank-3 sparse tensor.");
	}
	if (FEC.dim(2) != dlogmat.dim(0)) {
		throw std::runtime_error("FEC letter dimension does not match dlogmat first letter dimension.");
	}
	if (LEC.dim(1) != dlogmat.dim(1)) {
		throw std::runtime_error("LEC letter dimension does not match dlogmat second letter dimension.");
	}

	std::cout << "--- FEC ---" << std::endl;
	print_tensor_info(FEC);
	std::cout << "--- LEC ---" << std::endl;
	print_tensor_info(LEC);
	std::cout << "--- dlogmat ---" << std::endl;
	print_tensor_info(dlogmat);

	size_t alphabet_size = dlogmat.dim(1);
	size_t FEC_size = FEC.dim(0);
	size_t LEC_size = LEC.dim(0);

	Timer timer;
	timer.start();
	std::cout << "-- Computing FEC.dlogmat ..." << std::endl;
	auto tmp = tensor_contract(FEC, dlogmat, 2, 0, F, pool);
	tmp.flatten({{1, 3}, {0, 2}});
	std::cout << "-- Converting to matrix ..." << std::endl;
	auto FEC_dlogmat = tmp.to_sparse_mat(pool);
	tmp.clear();
	std::cout << "--- FEC.dlogmat ---" << std::endl;
	print_tensor_info(FEC_dlogmat);
	timer.stop();
	std::cout << "** FEC.dlogmat time: " << timer.milliseconds() << " ms" << std::endl;

	timer.start();
	std::cout << "-- Computing FEC.dlogmat RREF ..." << std::endl;
	FEC_dlogmat.sort_rows_by_nnz();
	sparse_mat_rref_reconstruct(FEC_dlogmat, opt);
	std::cout << "--- FEC.dlogmat RREF ---" << std::endl;
	print_tensor_info(FEC_dlogmat);
	timer.stop();
	std::cout << "** FEC.dlogmat RREF time: " << timer.milliseconds() << " ms" << std::endl;

	// Keep the row space but discard zero rows before converting the matrix
	// back to a tensor for the second contraction.
	timer.start();
	std::cout << "-- Post-processing (cancel common divisor, clear zero and sort) ..." << std::endl;
	pool->detach_loop(0, FEC_dlogmat.nrow, [&](size_t i) {
		FEC_dlogmat[i].compress();
		vec_cancel_divisor(FEC_dlogmat[i]);
	});
	pool->wait();
	FEC_dlogmat.clear_zero_row();
	FEC_dlogmat.sort_rows_by_nnz();
	std::cout << "--- FEC.dlogmat RREF simplified ---" << std::endl;
	print_tensor_info(FEC_dlogmat);
	timer.stop();
	std::cout << "** Post-process time: " << timer.milliseconds() << " ms" << std::endl;

	std::cout << "-- Converting to tensor ..." << std::endl;
	sparse_tensor<T, index_t, SPARSE_COO> FEC_dlogmat_tensor(FEC_dlogmat, pool);
	const size_t FEC_dlogmat_nrow = FEC_dlogmat.nrow;
	FEC_dlogmat.clear();
	FEC_dlogmat_tensor.reshape({FEC_dlogmat_nrow, FEC_size, alphabet_size});

	std::cout << "-- Computing FEC.dlogmat.LEC ..." << std::endl;
	timer.start();
	tmp = tensor_contract(FEC_dlogmat_tensor, LEC, 2, 1, F, pool);
	tmp.flatten({{0, 3}, {1, 2}});
	std::cout << "-- Converting to matrix ..." << std::endl;
	auto FEC_dlogmat_LEC = tmp.to_sparse_mat(pool);
	tmp.clear();
	std::cout << "--- FEC.dlogmat.LEC ---" << std::endl;
	print_tensor_info(FEC_dlogmat_LEC);
	timer.stop();
	std::cout << "** FEC.dlogmat.LEC time: " << timer.milliseconds() << " ms" << std::endl;

	// The final kernel is the sewing matrix, with layout
	// {sew_basis, FEC_basis, LEC_basis}.
	timer.start();
	std::cout << "-- Computing FEC.dlogmat.LEC RREF ..." << std::endl;
	auto pivots = sparse_mat_rref_reconstruct(FEC_dlogmat_LEC, opt);
	std::cout << "--- FEC.dlogmat.LEC RREF ---" << std::endl;
	print_tensor_info(FEC_dlogmat_LEC);
	auto kernel = sparse_mat_rref_kernel(FEC_dlogmat_LEC, pivots, F, opt).transpose();
	timer.stop();
	std::cout << "** FEC.dlogmat.LEC RREF time: " << timer.milliseconds() << " ms" << std::endl;

	timer.start();
	std::cout << "-- Post-processing (cancel common divisor, clear zero and sort) ..." << std::endl;
	pool->detach_loop(0, kernel.nrow, [&](size_t i) {
		kernel[i].canonicalize();
		vec_cancel_divisor(kernel[i]);
	});
	pool->wait();
	kernel.sort_rows_by_nnz();
	timer.stop();
	std::cout << "** Post-process time: " << timer.milliseconds() << " ms" << std::endl;

	// A weight-1 side can yield an empty kernel (zero sewing basis); the
	// sparse_tensor constructor degenerates on 0 rows, so build the empty
	// tensor explicitly.
	sparse_tensor<T, index_t, SPARSE_COO> kernel_tensor(
		kernel.nrow > 0 ? sparse_tensor<T, index_t, SPARSE_COO>(kernel, pool)
		                : sparse_tensor<T, index_t, SPARSE_COO>({kernel.nrow, FEC_size, LEC_size}));
	kernel_tensor.reshape({kernel.nrow, FEC_size, LEC_size});

	std::cout << "--- kernel ---" << std::endl;
	print_tensor_info(kernel_tensor);

	return sparse_tensor<T, index_t, SPARSE_CSR>(std::move(kernel_tensor), pool);
}
