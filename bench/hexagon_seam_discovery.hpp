#pragma once
#include "staged_kernel.hpp"
namespace staged_kernel {
template<class Source> MM discover_candidate(const Source& source, const field_t& fp, rref_option_t& options_rref, Options options) {
    auto opt=options_rref; auto pool=&opt->pool;
    structured_kernel::Relations relations(source.ncol());
    std::vector<I> columns(source.ncol());std::iota(columns.begin(),columns.end(),0);
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
    auto k = modular_kernel(reduced.first, reduced.second, fp);
    std::cout << "JOINT_DISCOVERY nullity=" << k.nrow << " complete=" << complete_discovery << " nnz=" << k.nnz() << std::endl;
    return k;
}
}
