#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <array>
#include <numeric>
#include <stdexcept>
#include "./SparseRREF/sparse_mat.h"
#include "./SparseRREF/sparse_tensor.h"
#include "./SparseRREF/wxf_support.h"
#include "tensor_shuffle.h"

using namespace SparseRREF;

using index_t = int32_t;

// Use hardware concurrency for optimal performance
const size_t n_of_threads = std::thread::hardware_concurrency() > 0 ? 
    std::thread::hardware_concurrency() : 8;  // Fallback to 8 if detection fails

template <typename T>
void print_tensor_info(const sparse_tensor<T, index_t>& S) {
	std::cout << "dims: ";
	for (size_t i = 0; i < S.rank(); i++)
		std::cout << S.dim(i) << " ";
	std::cout << std::endl;
	std::cout << "nnz: " << S.nnz() << " alloc: " << S.alloc() << std::endl;
}

template <typename T>
void print_tensor_info(const sparse_mat<T, index_t>& S) {
	std::cout << "dims: " << S.nrow << " " << S.ncol << std::endl;
	std::cout << "nnz: " << S.nnz() << " alloc: " << S.alloc() << std::endl;
}

int main(int argc, char* argv[]) {
    try {

    field_t F(FIELD_QQ);
    using scalar_t = rat_t;

    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <tensor1.wxf> <tensor2.wxf> <weight1> <weight2> <output.wxf>" << std::endl;
        std::cerr << "  weight1: rational number for first tensor (e.g., '1/2' or '3')" << std::endl;
        std::cerr << "  weight2: rational number for second tensor (e.g., '2/3' or '5')" << std::endl;
        return 1;
    }

    // Parse rational weights from command line arguments
    rat_t weight1, weight2;
    try {
        weight1 = rat_t(argv[3]);
        weight2 = rat_t(argv[4]);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing rational numbers: " << e.what() << std::endl;
        std::cerr << "Please use format like '1/2' or '3' for rational numbers" << std::endl;
        return 1;
    }

    std::cout << "Weighted tensor addition: " << weight1 << " * " << argv[1] << " + " << weight2 << " * " << argv[2] << " will be stored in " << argv[5] << std::endl;

    auto start_time_chrono = std::chrono::system_clock::now();
    std::time_t start_time_c = std::chrono::system_clock::to_time_t(start_time_chrono);
    std::cout << "Task begin at: " << std::ctime(&start_time_c);

    rref_option_t opt;
	// opt->method = 0;
	// opt->verbose = true;
	opt->pool.reset(n_of_threads);  // number of threads
	thread_pool* pool = &(opt->pool);

    std::filesystem::path base = std::filesystem::path(argv[0]).parent_path();

    std::cout << "Reading file " << argv[1] << " and " << argv[2] << " ..." << std::endl;
    std::filesystem::path path = base / argv[1];
	sparse_tensor<scalar_t, index_t, SPARSE_COO> tensor_first(sparse_tensor_read_wxf<scalar_t, index_t>(path, F, pool));
    print_tensor_info(tensor_first);
    std::cout << std::endl;

    path = base / argv[2];
    sparse_tensor<scalar_t, index_t, SPARSE_COO> tensor_second(sparse_tensor_read_wxf<scalar_t, index_t>(path, F, pool));
    print_tensor_info(tensor_second);
    std::cout << "Reading file " << argv[1] << " and " << argv[2] << " finished" << std::endl;
    std::cout << std::endl;
    
    auto start_time = std::chrono::steady_clock::now();
    tensor_add_weighted(tensor_first, tensor_second, weight1, weight2, F);
    tensor_first.canonicalize();
    tensor_first.sort_indices();
    tensor_first.reserve(tensor_first.nnz());

    auto end_time = std::chrono::steady_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << "** weighted tensor addition finished in: " << duration.count() << " ms **" << std::endl;
    print_tensor_info(tensor_first);
    std::cout << std::endl;

    std::cout << "Writing results to " << argv[5] << " ..." << std::endl;
    path = base / argv[5];
	auto u8arr = sparse_tensor_write_wxf(sparse_tensor<scalar_t, index_t, SPARSE_CSR>(tensor_first));
	std::ofstream ofs_transform(path, std::ios::binary);
	if (!ofs_transform) {
		throw std::runtime_error("Cannot write file: " + path.string());
	}
	ofs_transform.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs_transform.flush();
	if (!ofs_transform.good()) {
		throw std::runtime_error("Failed writing file (disk full or I/O error?): " + path.string());
	}
	ofs_transform.close();
	u8arr.clear();
	u8arr.shrink_to_fit();
	end_time = std::chrono::steady_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << "--- Writing finished in: " << duration.count() << " ms ---" << std::endl;

    auto end_time_chrono = std::chrono::system_clock::now();
	std::time_t end_time_c = std::chrono::system_clock::to_time_t(end_time_chrono);
	std::cout << "Task end at: " << std::ctime(&end_time_c) << std::endl;

    return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}