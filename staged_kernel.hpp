#pragma once
#include "structured_kernel.hpp"
#include <tuple>

// Exact kernel intersection for replayable operators. No symmetry, tensor,
// alphabet, expected physical dimension, or polygon is part of this layer.
namespace staged_kernel {
using namespace symrep;
using MV = sparse_vec<ulong, I>;
using MM = sparse_mat<ulong, I>;

struct Options {
    size_t seed_terms = 64;
    size_t target_nullity = 128;
    size_t complete_through_terms = 1024;
    bool remember_equations = true;
    bool sample_equations = true;
    size_t batch_rows = 4096;
    size_t batch_nonzeros = 4 * 1024 * 1024;
    // Bound the initial global front too. Retained elimination fill and the
    // required output basis are separate costs, not bounded by these options.
    size_t seed_nonzeros = 64 * 1024 * 1024;
    // A discovery plateau may justify intersection before the strict target.
    // This is bounded by candidate dimension and never changes acceptance.
    size_t stagnation_batches = 4;
    size_t stagnation_rank_gain = 2;
    size_t max_residual_nullity = 256;
};

struct Selection {
    size_t max_terms = SIZE_MAX;
    const streaming_echelon::Substitution* substitution = nullptr;
    unsigned sample_shift = 0;
};

// A non-owning conjunction of operators. Each operator provides exact short
// rows, modular discovery rows, complete residual blocks, and a height bound.
// Discovery masks are local to one prime and never used for residual checks.
template<class... Operators> class Stack {
    std::tuple<const Operators&...> operators;
public:
    explicit Stack(const Operators&... values) : operators(values...) {
        static_assert(sizeof...(Operators) > 0);
        std::apply([&](const auto&... op) {
            (require(op.ncol() == ncol(), "stacked operator dimensions differ"), ...);
        }, operators);
    }
    size_t ncol() const { return std::get<0>(operators).ncol(); }
    auto masks(bool remember) const {
        return std::apply([&](const auto&... op) {
            std::vector<std::vector<unsigned char>> result;
            result.reserve(sizeof...(Operators));
            (result.emplace_back(remember ? op.nrow() : 0), ...);
            return result;
        }, operators);
    }
    template<class Emit> void operator()(structured_kernel::Relations& r, size_t limit, Emit&& emit) const {
        std::apply([&](const auto&... op) { (op(r, limit, emit), ...); }, operators);
    }
    template<class Emit> void discover(structured_kernel::Relations& r, const std::vector<I>& columns,
        const field_t& fp, const Selection& select, std::vector<std::vector<unsigned char>>& seen,
        Emit&& emit, thread_pool* pool) const {
        size_t i = 0;
        std::apply([&](const auto&... op) {
            auto run = [&](const auto& source) {
                auto& mask = seen.at(i++);
                source.modular(r, columns, fp, emit, pool, select.max_terms, 0,
                    select.substitution, mask.empty() ? nullptr : &mask, select.sample_shift);
            };
            (run(op), ...);
        }, operators);
    }
    size_t residual_blocks() const {
        return std::apply([](const auto&... op) { return (op.residual_blocks() + ...); }, operators);
    }
    template<class Emit> void residual(const MM& k, const field_t& fp, size_t block,
        Emit&& emit, thread_pool* pool) const {
        require(block < residual_blocks(), "invalid residual block");
        bool done = false;
        std::apply([&](const auto&... op) {
            auto run = [&](const auto& source) {
                if (done) return;
                if (block < source.residual_blocks()) {
                    source.residual(k, fp, block, emit, pool);
                    done = true;
                } else block -= source.residual_blocks();
            };
            (run(op), ...);
        }, operators);
    }
    int_t integral_residual_bound(const Mat& k) const {
        int_t bound = 0;
        std::apply([&](const auto&... op) {
            auto include = [&](const auto& source) {
                auto h = source.integral_residual_bound(k);
                if (h > bound) bound = std::move(h);
            };
            (include(op), ...);
        }, operators);
        return bound;
    }
};
template<class... Operators> auto stack(Operators&... op) { return Stack<Operators...>(op...); }

inline void normalize_mod(MV& row, const field_t& fp) {
    row.sort_indices();
    size_t out = 0;
    for (size_t first = 0; first < row.nnz();) {
        size_t end = first + 1;
        I c = row(first);
        ulong x = row[first];
        while (end < row.nnz() && row(end) == c) x = nmod_add(x, row[end++], fp.mod);
        if (x) { row(out) = c; row[out] = x; ++out; }
        first = end;
    }
    row.resize(out); row.compress();
}

// Rows span the right kernel. The private free coordinate is -1, matching
// certified_kernel's rational-reconstruction convention.
inline MM modular_kernel(const MM& a, const std::vector<std::vector<pivot_t<I>>>& pivots, const field_t& fp) {
    std::vector<I> pr(a.ncol, -1), free(a.ncol, -1);
    for (const auto& batch : pivots) for (auto [r,c] : batch) pr[c] = r;
    size_t n = 0;
    for (size_t c = 0; c < a.ncol; ++c) if (pr[c] < 0) free[c] = index(n++);
    MM k(n, a.ncol); std::vector<size_t> count(n, 1);
    for (size_t c = 0; c < a.ncol; ++c) if (pr[c] >= 0)
        for (auto [j,x] : a[pr[c]]) if (size_t(j) != c && x) {
            require(free[j] >= 0, "staged kernel needs RREF"); ++count[free[j]];
        }
    for (size_t i = 0; i < n; ++i) k[i].reserve(count[i]);
    for (size_t c = 0; c < a.ncol; ++c) {
        if (pr[c] < 0) k[free[c]].push_back(index(c), fp.mod.n - 1);
        else for (auto [j,x] : a[pr[c]]) if (size_t(j) != c && x) k[free[j]].push_back(index(c), x);
    }
    return k;
}

// Once K spans the current candidate space, imposing B replaces it by
// ker(B K^T) K. Every restriction keeps all previously satisfied equations.
template<class Source> struct RefinedRows {
    const Source& source;
    rref_option* opt;
    Options options;
    template<class Emit> void operator()(structured_kernel::Relations& r, size_t limit, Emit&& emit) const {
        source(r, limit, emit);
    }
    template<class Emit> void modular(structured_kernel::Relations& relations, const std::vector<I>& columns,
        const field_t& fp, Emit&& emit, thread_pool* pool) const;
};

template<class Source> template<class Emit>
void RefinedRows<Source>::modular(structured_kernel::Relations& relations, const std::vector<I>& columns,
    const field_t& fp, Emit&& emit, thread_pool* pool) const {
    require(options.seed_terms && options.seed_nonzeros, "staged front limits must be positive");
    size_t n = 0;
    for (auto c : columns) if (c >= 0) n = std::max(n, size_t(c) + 1);
    auto start = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); };
    auto seen = source.masks(options.remember_equations || options.sample_equations);
    Selection selected{options.seed_terms};
    streaming_echelon::Reducer reducer(n, fp, opt, options.batch_rows, options.batch_nonzeros);
    MM a(0, n); size_t seed_nnz = 0, seed_rows = 0; bool seeded = false;
    source.discover(relations, columns, fp, selected, seen, [&](MV&& row) {
        ++seed_rows; seed_nnz += row.nnz();
        if (seeded) reducer.consume(std::move(row));
        else {
            a.rows.push_back(std::move(row));
            if (seed_nnz >= options.seed_nonzeros) {
                a.nrow = a.rows.size(); reducer.seed(std::move(a)); seeded = true;
            }
        }
    }, pool);
    if (!seeded) { a.nrow = a.rows.size(); reducer.seed(std::move(a)); }
    else reducer.flush_pending();
    std::cout << "REFINE_SEED rows=" << seed_rows << " cols=" << n << " nnz=" << seed_nnz
              << " bounded=" << seeded << " elapsed=" << elapsed() << std::endl;
    struct Enough {};
    size_t threshold = options.seed_terms;
    bool complete_discovery = threshold == SIZE_MAX;
    bool plateau = false;
    while (n - reducer.rank() > options.target_nullity && !complete_discovery && !plateau) {
        size_t zeros = reducer.prune_zeros();
        auto active_columns = columns;
        for (auto& c : active_columns) if (c >= 0 && reducer.is_zero(c)) c = -1;
        selected.substitution = &reducer.substitution();
        if (threshold == SIZE_MAX) selected.sample_shift = selected.sample_shift > 2 ? selected.sample_shift - 2 : 0;
        else {
            threshold = threshold >= 16384 ? SIZE_MAX : threshold * 4;
            selected.sample_shift = options.sample_equations && threshold > options.complete_through_terms ? 4 : 0;
        }
        selected.max_terms = threshold;
        std::cout << "REFINE_ZERO_SUBSTITUTION columns=" << zeros << " limit=" << threshold
                  << " sample_denominator=" << (UINT64_C(1) << selected.sample_shift) << std::endl;
        bool stopped = false;
        size_t last_batch = reducer.completed_batches(), last_rank = reducer.rank(), stagnant = 0;
        try {
            source.discover(relations, active_columns, fp, selected, seen, [&](MV&& row) {
                reducer.consume(std::move(row));
                if (n - reducer.rank() <= options.target_nullity) throw Enough{};
                if (reducer.completed_batches() != last_batch) {
                    size_t gain = reducer.rank() - last_rank;
                    stagnant = gain <= options.stagnation_rank_gain ? stagnant + 1 : 0;
                    last_batch = reducer.completed_batches(); last_rank = reducer.rank();
                    if (options.stagnation_batches && stagnant >= options.stagnation_batches
                        && n - reducer.rank() <= options.max_residual_nullity) {
                        plateau = true;
                        std::cout << "REFINE_STAGNATION nullity=" << n - reducer.rank()
                                  << " batches=" << stagnant << " elapsed=" << elapsed() << std::endl;
                        throw Enough{};
                    }
                }
            }, pool);
        } catch (const Enough&) { stopped = true; }
        reducer.flush_pending();
        complete_discovery = !stopped && threshold == SIZE_MAX && !selected.sample_shift;
        std::cout << "REFINE_DISCOVERY limit=" << threshold << " rank=" << reducer.rank()
                  << " nullity=" << n - reducer.rank() << " elapsed=" << elapsed() << std::endl;
    }
    auto reduced = reducer.finish();
    if (complete_discovery) {
        // Every original equation is already represented. Reuse that row
        // basis directly: building, canonicalizing and then annihilating a
        // large kernel would undo the advantage of ordinary elimination.
        size_t rank = 0;
        for (const auto& batch : reduced.second) for (auto [r,c] : batch) {
            ++rank; emit(std::move(reduced.first[r]));
        }
        std::cout << "REFINE_COMPLETE nullity=" << n - rank
                  << " complete_discovery=1 full_operator=1 elapsed=" << elapsed() << std::endl;
        return;
    }
    auto k = modular_kernel(reduced.first, reduced.second, fp);
    reduced.first.clear();
    size_t visited = 0;
    // Discovery was partial. Every omitted equation is imposed here, without
    // sampling, by contracting the original operators into the candidate.
    if (!complete_discovery && k.nrow) {
        std::map<I, MV> residual;
        auto accept = [&](MV&& r) {
            ++visited;
            while (r.nnz()) {
                auto it = residual.find(r(0));
                if (it == residual.end()) {
                    sparse_vec_rescale(r, nmod_inv(r[0], fp.mod), fp);
                    I c = r(0); residual.emplace(c, std::move(r)); return;
                }
                auto x = r[0]; sparse_vec_sub_mul(r, it->second, x, fp);
            }
        };
        std::vector<I> mapped(columns.size(), -1);
        std::vector<ulong> scale(columns.size());
        for (size_t c = 0; c < columns.size(); ++c) {
            auto [r,x] = relations.image(c);
            if (x == 0) continue;
            if (x.den() % fp.mod == 0) throw certified_kernel::BadPrime("residual lifting denominator");
            mapped[c] = columns[r]; scale[c] = x % fp.mod;
        }
        auto lift = [&]() {
            auto kt = k.transpose(); MM full(k.nrow, columns.size());
            for (size_t c = 0; c < columns.size(); ++c) if (mapped[c] >= 0)
                for (auto [i,x] : kt[mapped[c]]) {
                    auto v = nmod_mul(x, scale[c], fp.mod);
                    if (v) full[i].push_back(index(c), v);
                }
            return full;
        };
        auto full = lift();
        for (size_t block = 0; block < source.residual_blocks() && k.nrow; ++block) {
            source.residual(full, fp, block, accept, pool);
            if (!residual.empty()) {
                MM r(0, k.nrow);
                for (auto& [c,row] : residual) r.rows.push_back(std::move(row));
                residual.clear(); r.nrow = r.rows.size();
                auto pivots = sparse_mat_rref(r, fp, opt);
                auto small = modular_kernel(r, pivots, fp);
                auto before = k.nrow;
                k = sparse_mat_mul(small, k, fp, pool); full = lift();
                std::cout << "REFINE_INTERSECT dimension=" << before << "->" << k.nrow
                          << " elapsed=" << elapsed() << std::endl;
            }
        }
    }
    // Canonicalize via a small coordinate chart, then expose its annihilator
    // to the shared rank/CRT reconstruction backend.
    auto kp = sparse_mat_rref(k, fp, opt);
    std::vector<I> free(k.nrow, -1), owner(n, -1);
    for (const auto& batch : kp) for (auto [r,c] : batch) { free[r] = c; owner[c] = r; }
    for (size_t r = 0; r < k.nrow; ++r) require(free[r] >= 0, "dependent staged kernel rows");
    auto kt = k.transpose();
    for (size_t c = 0; c < n; ++c) if (owner[c] < 0) {
        MV row; row.push_back(index(c), ulong(1));
        for (auto [r,x] : kt[c]) row.push_back(free[r], nmod_neg(x, fp.mod));
        normalize_mod(row, fp); emit(std::move(row));
    }
    std::cout << "REFINE_COMPLETE nullity=" << k.nrow << " residual_nonzero_rows=" << visited
              << " complete_discovery=" << complete_discovery << " full_operator=1 elapsed=" << elapsed() << std::endl;
}

// The caller holds COMPLETE modular kernel congruences for crt_modulus.
// Additional primes certify the small candidate without repeating discovery.
template<class Source> bool verify(const Mat& candidate, const Source& source, rref_option_t opt,
    ulong last_prime, const int_t& crt_modulus) {
    auto bound = source.integral_residual_bound(candidate);
    int_t modulus = crt_modulus;
    std::cout << "combined_certificate_bound_bits=" << bound.bits() << std::endl;
    struct Nonzero {};
    for (size_t attempt = 0; !(bound < modulus) && attempt < 16; ++attempt) {
        last_prime = n_nextprime(last_prime, 0); field_t fp(FIELD_Fp, last_prime);
        try {
            MM k(candidate.nrow, candidate.ncol);
            for (size_t r = 0; r < candidate.nrow; ++r) k[r] = candidate[r] % fp.mod;
            auto reject = [](MV&& row) { if (row.nnz()) throw Nonzero{}; };
            for (size_t block = 0; block < source.residual_blocks(); ++block)
                source.residual(k, fp, block, reject, &opt->pool);
            modulus *= last_prime;
            std::cout << "extra_residual_prime=" << last_prime << " certified_modulus_bits=" << modulus.bits() << std::endl;
        } catch (const certified_kernel::BadPrime&) {
            std::cout << "skipped_residual_prime=" << last_prime << std::endl;
        } catch (const Nonzero&) { return false; }
    }
    if (!(bound < modulus)) return false;
    std::cout << "exact_certificate=complete_factored_residues_plus_height modulus_bits=" << modulus.bits() << std::endl;
    return true;
}

template<class Source> Mat solve(const Source& source, rref_option_t opt, const Options& options = {},
    structured_kernel::Statistics* stats = nullptr) {
    RefinedRows<Source> refined{source, opt, options};
    return structured_kernel::solve_stream(source.ncol(), refined, opt, stats,
        [&](const Mat& k, ulong prime, const int_t& modulus) { return verify(k, source, opt, prime, modulus); });
}
} // namespace staged_kernel
