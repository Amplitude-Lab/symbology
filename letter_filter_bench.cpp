// letter_filter_bench.cpp — Semantics verification + efficiency benchmark for
// the letter-space support filters (solve_collinear.hpp sentinels).
//
// Part A (semantics): hand-crafted (5,11,11,11) seed with divergent letters
//   {0,1}. Checks the production dispatcher apply_letter_projection_ab for
//   "divergent" / "finite" / "identity" against exact expected (key, value)
//   sets, and verifies the subtract construction
//     support(any-div) = support(whole) \ support(finite contraction)
//   where the finite contraction uses a label-preserving diagonal selector
//   (11x11 identity restricted to finite letters), contracted per letter slot
//   via the production apply_colprojdiv_slots.
//
// Part B (benchmark): D=1009-letter alphabet, 3 letter slots, 8 unknowns,
//   N=1,000,000 unique entries (identical across all runs). For each
//   divergent-letter count n_div in {1,2,3,5,8,16,64,256,504}:
//     method 1 (direct):     apply_letter_filter(any_divergent) — one key scan
//     method 2 (subtract):   diagonal finite selector contracted into each
//                            letter slot, key-set build, then whole-scan
//                            keep-not-in-set
//   Both must produce identical supports and values (verified). Median of R
//   runs, with subtract broken into contract / set / scan phases.
//
// Also writes /tmp/lfbench/{data/colprojdiv.wxf, seed.wxf} for the bootstrap
// end-to-end check, and provides --info <file.wxf> mode that prints
// rank/dims/nnz of a WXF tensor (used to count exported condition rows).
//
// Usage:
//   ./letter_filter_bench              # Part A + Part B
//   ./letter_filter_bench --info f.wxf
#include "solve_collinear.hpp"

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>

using index_t = int32_t;
using scalar_t = rat_t;

using key_val_map = std::map<std::vector<index_t>, scalar_t>;

static double ms_since(std::chrono::steady_clock::time_point t0) {
	return std::chrono::duration<double, std::milli>(
		std::chrono::steady_clock::now() - t0).count();
}

static key_val_map collect_entries(
	const sparse_tensor<scalar_t, index_t, SPARSE_COO>& t) {
	key_val_map m;
	for (auto i : t.gen_perm()) {
		m[t.index_vector(i)] = t.val(i);
	}
	return m;
}

static bool maps_equal(const key_val_map& a, const key_val_map& b, const char* what) {
	if (a.size() != b.size()) {
		std::cout << "   [FAIL] " << what << ": size " << a.size() << " != " << b.size() << std::endl;
		return false;
	}
	for (const auto& [k, v] : a) {
		auto it = b.find(k);
		if (it == b.end() || !(it->second == v)) {
			std::cout << "   [FAIL] " << what << ": key/value mismatch at [";
			for (auto x : k) std::cout << x << ",";
			std::cout << "]" << std::endl;
			return false;
		}
	}
	std::cout << "   [OK] " << what << " (" << a.size() << " entries match)" << std::endl;
	return true;
}

// Diagonal label-preserving finite selector: (D, D) with 1 at (i, i) for
// finite letters i, zero rows on divergent letters.
static sparse_tensor<scalar_t, index_t, SPARSE_COO> diagonal_finite_selector(
	size_t D, const std::set<index_t>& div_letters) {
	sparse_tensor<scalar_t, index_t, SPARSE_COO> sel(std::vector<size_t>{D, D});
	for (index_t i = 0; i < (index_t)D; i++) {
		if (div_letters.count(i) == 0) {
			sel.push_back(std::vector<index_t>{i, i}, scalar_t(1));
		}
	}
	return sel;  // pushed in ascending i: already lexicographically sorted
}

// Subtract variant: contract a diagonal finite selector into each letter slot
// of a COPY of the whole tensor (labels preserved by diagonality), build the
// finite key set, then scan the WHOLE tensor keeping entries not in the set.
struct subtract_timing_t {
	double contract_ms = 0, set_ms = 0, scan_ms = 0;
	double total_ms() const { return contract_ms + set_ms + scan_ms; }
};

static uint64_t pack_key(const std::vector<index_t>& k, size_t D) {
	uint64_t p = 0;
	for (auto x : k) p = p * D + (uint64_t)(x < 0 ? 0 : x);
	return p;
}

static sparse_tensor<scalar_t, index_t, SPARSE_COO> filter_any_div_subtract(
	const sparse_tensor<scalar_t, index_t, SPARSE_COO>& whole,
	size_t first_slot_axis, size_t n_slots, size_t D,
	const std::set<index_t>& div_letters,
	const field_t& F, thread_pool* pool,
	subtract_timing_t& tm) {

	auto sel = diagonal_finite_selector(D, div_letters);

	// Phase 1: finite contraction on a copy of the whole tensor.
	auto t0 = std::chrono::steady_clock::now();
	sparse_tensor<scalar_t, index_t, SPARSE_COO> finite_view(whole);
	auto finite = apply_colprojdiv_slots<scalar_t, index_t>(
		std::move(finite_view), sel, first_slot_axis, n_slots, F, pool);
	auto t1 = std::chrono::steady_clock::now();
	tm.contract_ms = ms_since(t0);

	// Phase 2: finite key set (packed into uint64 for fast hashing).
	std::unordered_set<uint64_t> finite_keys;
	finite_keys.reserve(2 * finite.nnz() + 16);
	for (auto i : finite.gen_perm()) {
		finite_keys.insert(pack_key(finite.index_vector(i), D));
	}
	auto t2 = std::chrono::steady_clock::now();
	tm.set_ms = ms_since(t1);

	// Phase 3: whole-scan, keep entries whose key is not in the finite set.
	sparse_tensor<scalar_t, index_t, SPARSE_COO> result(whole.dims());
	for (auto i : whole.gen_perm()) {
		auto k = whole.index_vector(i);
		if (finite_keys.count(pack_key(k, D)) == 0) {
			result.push_back(k, whole.val(i));
		}
	}
	auto t3 = std::chrono::steady_clock::now();
	tm.scan_ms = ms_since(t2);
	return result;
}

static int run_part_a(const std::filesystem::path& bench_dir, field_t& F, rref_option_t& opt) {
	thread_pool* pool = &(opt->pool);
	std::cout << "================ PART A: semantics (5,11,11,11), div={0,1} ================" << std::endl;

	// Crafted letter keys over slots (l1,l2,l3):
	//   (2,3,4) all-finite      (0,3,4) mixed   (2,0,5) mixed
	//   (5,6,1) mixed           (0,1,2) all-divergent
	const std::vector<std::vector<index_t>> keys = {
		{2, 3, 4}, {0, 3, 4}, {2, 0, 5}, {5, 6, 1}, {0, 1, 2}};

	sparse_tensor<scalar_t, index_t, SPARSE_COO> A(std::vector<size_t>{5, 11, 11, 11});
	sparse_tensor<scalar_t, index_t, SPARSE_COO> b(std::vector<size_t>{11, 11, 11});
	for (size_t j = 0; j < keys.size(); j++) {
		for (index_t u = 0; u < 5; u++) {
			std::vector<index_t> idx{u, keys[j][0], keys[j][1], keys[j][2]};
			A.push_back(idx, scalar_t((long)(j + 1), (long)(u + 2)));
		}
	}
	b.push_back(std::vector<index_t>{2, 3, 4}, scalar_t(7));
	b.push_back(std::vector<index_t>{0, 1, 2}, scalar_t(-3));

	// Expected supports
	key_val_map exp_any_A, exp_fin_A, exp_any_b, exp_fin_b;
	for (size_t j = 0; j < keys.size(); j++) {
		bool has_div = false, all_fin = true;
		for (auto l : keys[j]) {
			if (l == 0 || l == 1) { has_div = true; all_fin = false; }
		}
		(void)all_fin;  // has_div == !all_fin for 3 slots
		for (index_t u = 0; u < 5; u++) {
			std::vector<index_t> idx{u, keys[j][0], keys[j][1], keys[j][2]};
			auto v = scalar_t((long)(j + 1), (long)(u + 2));
			if (has_div) exp_any_A[idx] = v; else exp_fin_A[idx] = v;
		}
	}
	exp_any_b[{0, 1, 2}] = scalar_t(-3);
	exp_fin_b[{2, 3, 4}] = scalar_t(7);

	// Write bench data dir: colprojdiv.wxf with nonzero rows {0,1} (11x2)
	auto data_dir = bench_dir / "data";
	std::filesystem::create_directories(data_dir);
	{
		sparse_tensor<scalar_t, index_t, SPARSE_COO> div(std::vector<size_t>{11, 2});
		div.push_back(std::vector<index_t>{0, 0}, scalar_t(1));
		div.push_back(std::vector<index_t>{1, 1}, scalar_t(1));
		div.canonicalize();
		div.sort_indices();
		sparse_tensor<scalar_t, index_t, SPARSE_CSR> div_csr(std::move(div));
		projection_write_tensor<scalar_t, index_t>(data_dir / "colprojdiv.wxf", std::move(div_csr), pool);
	}

	int failures = 0;

	// dispatcher: divergent
	{
		auto Aw = A, bw = b;
		apply_letter_projection_ab<scalar_t, index_t>("divergent", Aw, bw, b.rank(), data_dir, F, opt);
		failures += !maps_equal(collect_entries(Aw), exp_any_A, "divergent: A entries");
		failures += !maps_equal(collect_entries(bw), exp_any_b, "divergent: b entries");
	}
	// dispatcher: finite
	{
		auto Aw = A, bw = b;
		apply_letter_projection_ab<scalar_t, index_t>("finite", Aw, bw, b.rank(), data_dir, F, opt);
		failures += !maps_equal(collect_entries(Aw), exp_fin_A, "finite: A entries");
		failures += !maps_equal(collect_entries(bw), exp_fin_b, "finite: b entries");
	}
	// dispatcher: identity
	{
		auto Aw = A, bw = b;
		apply_letter_projection_ab<scalar_t, index_t>("identity", Aw, bw, b.rank(), data_dir, F, opt);
		failures += !maps_equal(collect_entries(Aw), collect_entries(A), "identity: A unchanged");
		failures += !maps_equal(collect_entries(bw), collect_entries(b), "identity: b unchanged");
	}

	// Partition: any-div ∪ all-fin = whole (sizes add up, disjoint)
	{
		size_t whole = A.nnz();
		size_t parts = exp_any_A.size() + exp_fin_A.size();
		bool ok = (whole == parts);
		std::cout << "   [" << (ok ? "OK" : "FAIL") << "] partition: |any-div| " << exp_any_A.size()
		          << " + |all-fin| " << exp_fin_A.size() << " = " << parts
		          << " == |whole| " << whole << std::endl;
		failures += !ok;
	}

	// Subtract construction: diagonal finite contraction support == finite
	// filter support; whole \ that == any-div filter support.
	{
		std::set<index_t> div_letters{0, 1};
		auto sel = diagonal_finite_selector(11, div_letters);
		auto Af = A;
		auto contracted = apply_colprojdiv_slots<scalar_t, index_t>(
			std::move(Af), sel, 1, 3, F, pool);
		failures += !maps_equal(collect_entries(contracted), exp_fin_A,
			"subtract: diagonal finite contraction == finite filter");

		std::unordered_set<uint64_t> fin_keys;
		for (const auto& [k, v] : exp_fin_A) {
			fin_keys.insert(pack_key(k, 11));
		}
		key_val_map subtract_any;
		for (const auto& [k, v] : collect_entries(A)) {
			if (fin_keys.count(pack_key(k, 11)) == 0) subtract_any[k] = v;
		}
		failures += !maps_equal(subtract_any, exp_any_A,
			"subtract: whole \\ finite-contraction == any-divergent filter");
	}

	// Write the seed for the bootstrap end-to-end run
	{
		auto seed = A;
		seed.canonicalize();
		seed.sort_indices();
		sparse_tensor<scalar_t, index_t, SPARSE_CSR> seed_csr(std::move(seed));
		projection_write_tensor<scalar_t, index_t>(bench_dir / "seed.wxf", std::move(seed_csr), pool);
	}

	std::cout << (failures == 0 ? "PART A: ALL PASS" : "PART A: FAILURES!") << std::endl << std::endl;
	return failures;
}

static int run_part_b(field_t& F, rref_option_t& opt) {
	thread_pool* pool = &(opt->pool);
	constexpr size_t D = 1009;      // letter alphabet size
	constexpr size_t n_slots = 3;   // letter slots (A also has the unknown axis)
	constexpr size_t n_unk = 8;
	constexpr size_t N = 1000000;   // unique entries
	constexpr int R = 5;            // repetitions (median)

	std::cout << "================ PART B: benchmark (D=" << D << ", slots=" << n_slots
	          << ", N=" << N << ") ================" << std::endl;

	// Deterministic unique random entries (LCG), shared across all n_div
	sparse_tensor<scalar_t, index_t, SPARSE_COO> A(
		std::vector<size_t>{n_unk, D, D, D});
	{
		std::unordered_set<uint64_t> seen;
		seen.reserve(2 * N);
		uint64_t x = 0x9E3779B97F4A7C15ULL;
		while (seen.size() < N) {
			x = x * 6364136223846793005ULL + 1442695040888963407ULL;
			uint64_t s = x ^ (x >> 29);
			uint64_t u = s % n_unk;
			uint64_t l1 = (s / n_unk) % D;
			uint64_t l2 = (s / (n_unk * D)) % D;
			uint64_t l3 = (s / (n_unk * D * D)) % D;
			uint64_t p = ((u * D + l1) * D + l2) * D + l3;
			if (seen.insert(p).second) {
				scalar_t val((long)(1 + p % 7), (long)(1 + p % 5));
				A.push_back(std::vector<index_t>{
					(index_t)u, (index_t)l1, (index_t)l2, (index_t)l3}, val);
			}
		}
	}
	std::cout << "Generated " << A.nnz() << " unique entries" << std::endl;

	std::cout << std::left << std::setw(8) << "n_div"
	          << std::setw(10) << "any-kept" << std::setw(10) << "fin-kept"
	          << std::setw(14) << "direct[ms]" << std::setw(14) << "subtract[ms]"
	          << std::setw(12) << "contract" << std::setw(10) << "set" << std::setw(10) << "scan"
	          << std::setw(8) << "ratio" << "equal?" << std::endl;
	std::cout << std::string(100, '-') << std::endl;

	int failures = 0;
	for (size_t n_div : {1u, 2u, 3u, 5u, 8u, 16u, 64u, 256u, 504u}) {
		std::set<index_t> div_letters;
		for (index_t l = 0; l < (index_t)n_div; l++) div_letters.insert(l);

		// sanity/equality: one direct run vs one subtract run
		auto A_dir = A;
		auto direct = apply_letter_filter<scalar_t, index_t>(
			std::move(A_dir), letter_filter_t::any_divergent, 1, n_slots, div_letters);
		size_t kept_any = direct.nnz();

		subtract_timing_t tm_chk;
		auto sub = filter_any_div_subtract(A, 1, n_slots, D, div_letters, F, pool, tm_chk);
		bool eq = maps_equal(collect_entries(direct), collect_entries(sub),
			("n_div=" + std::to_string(n_div) + " equality").c_str());
		if (!eq) failures++;

		// timings (median of R)
		std::vector<double> t_direct, t_sub;
		std::vector<subtract_timing_t> tms;
		for (int r = 0; r < R; r++) {
			{
				auto Aw = A;
				auto t0 = std::chrono::steady_clock::now();
				auto out = apply_letter_filter<scalar_t, index_t>(
					std::move(Aw), letter_filter_t::any_divergent, 1, n_slots, div_letters);
				t_direct.push_back(ms_since(t0));
				(void)out;
			}
			{
				subtract_timing_t tm;
				auto t0 = std::chrono::steady_clock::now();
				auto out = filter_any_div_subtract(A, 1, n_slots, D, div_letters, F, pool, tm);
				t_sub.push_back(ms_since(t0));
				tms.push_back(tm);
				(void)out;
			}
		}
		auto median = [](std::vector<double> v) {
			std::sort(v.begin(), v.end());
			return v[v.size() / 2];
		};
		double md = median(t_direct);
		double ms_ = median(t_sub);
		size_t mi = 0;
		for (size_t i = 1; i < tms.size(); i++) {
			if (tms[i].total_ms() < tms[mi].total_ms()) mi = i;
		}

		// finite kept count (for the report)
		auto Af = A;
		auto fin = apply_letter_filter<scalar_t, index_t>(
			std::move(Af), letter_filter_t::all_finite, 1, n_slots, div_letters);

		std::cout << std::left << std::setw(8) << n_div
		          << std::setw(10) << kept_any << std::setw(10) << fin.nnz()
		          << std::setw(14) << std::fixed << std::setprecision(1) << md
		          << std::setw(14) << ms_
		          << std::setw(12) << tms[mi].contract_ms
		          << std::setw(10) << tms[mi].set_ms
		          << std::setw(10) << tms[mi].scan_ms
		          << std::setw(8) << std::setprecision(2) << (ms_ / md)
		          << (eq ? "yes" : "NO") << std::endl;
	}
	std::cout << "\nPART B: " << (failures == 0 ? "ALL EQUAL — benchmark valid" : "EQUALITY FAILURES!") << std::endl;
	return failures;
}

int main(int argc, char* argv[]) {
	std::filesystem::path bench_dir = "/tmp/lfbench";

	if (argc == 3 && std::string(argv[1]) == "--info") {
		field_t F(FIELD_QQ);
		rref_option_t opt;
		opt->method = 0;
		opt->verbose = false;
		opt->pool.reset();
		thread_pool* pool = &(opt->pool);
		auto t = projection_read_tensor<scalar_t, index_t>(argv[2], F, pool);
		std::cout << "rank=" << t.rank() << " dims=";
		for (size_t i = 0; i < t.rank(); i++) std::cout << t.dim(i) << (i + 1 < t.rank() ? "x" : "");
		std::cout << " nnz=" << t.nnz() << std::endl;
		return 0;
	}

	field_t F(FIELD_QQ);
	rref_option_t opt;
	opt->method = 0;
	opt->verbose = false;
	opt->pool.reset();

	int failures = run_part_a(bench_dir, F, opt);
	failures += run_part_b(F, opt);

	std::cout << "\nTOTAL: " << (failures == 0 ? "ALL PASS" : std::to_string(failures) + " FAILURES") << std::endl;
	return failures == 0 ? 0 : 1;
}
