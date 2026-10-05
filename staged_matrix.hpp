#pragma once
#include "staged_kernel.hpp"

namespace staged_kernel {
// Adapter for ordinary rational matrices. Owns no second copy of the input.
class MatrixRows {
    const Mat& matrix;
    int_t height = 0;
public:
    explicit MatrixRows(const Mat& a) : matrix(a) {
        for (const auto& row : a.rows) {
            int_t denominator = 1, norm = 0;
            for (auto [c,x] : row) denominator = LCM(denominator, x.den());
            for (auto [c,x] : row) { int_t term = denominator; term /= x.den(); term *= x.num().abs(); norm += term; }
            if (norm > height) height = std::move(norm);
        }
    }
    size_t ncol() const { return matrix.ncol; }
    size_t nrow() const { return matrix.nrow; }
    template<class Emit> void operator()(structured_kernel::Relations& relations, size_t limit, Emit&& emit) const {
        for (const auto& original : matrix.rows) {
            if (limit != SIZE_MAX) {
                size_t terms = 0;
                for (auto c : original.index_span()) if (!relations.is_zero(c) && ++terms > limit) break;
                if (terms > limit) continue;
            }
            auto row = original; relations.compress(row);
            if (row.nnz()) emit(std::move(row));
        }
    }
    template<class Emit> void modular(structured_kernel::Relations& relations, const std::vector<I>& columns,
        const field_t& fp, Emit&& emit, thread_pool* pool = nullptr, size_t max_terms = SIZE_MAX,
        size_t min_terms = 0, const streaming_echelon::Substitution* substitution = nullptr,
        std::vector<unsigned char>* seen = nullptr, unsigned sample_shift = 0) const {
        require(!seen || seen->size() == nrow(), "matrix discovery mask dimensions");
        require(sample_shift<64,"invalid matrix sampling shift");
        std::vector<I> mapped(ncol(), -1); std::vector<ulong> scales(ncol());
        for (size_t c = 0; c < ncol(); ++c) {
            auto [r,x] = relations.image(c);
            if (x == 0) continue;
            if (x.den() % fp.mod == 0) throw certified_kernel::BadPrime("matrix relation denominator");
            mapped[c] = columns[r]; scales[c] = x % fp.mod;
            if (substitution && mapped[c] >= 0) {
                I j = mapped[c]; mapped[c] = substitution->image[j];
                scales[c] = nmod_mul(scales[c], substitution->scale[j], fp.mod);
            }
        }
        for (size_t r = 0; r < nrow(); ++r) {
            if ((seen && (*seen)[r]) || !streaming_echelon::sampled_equation(r, sample_shift, UINT64_C(0x9e3779b97f4a7c15))) continue;
            MV row; bool skipped = false;
            for (auto [c,x] : matrix[r]) {
                if (x.den() % fp.mod == 0) throw certified_kernel::BadPrime("matrix denominator");
                if (mapped[c] < 0) continue;
                ulong value = nmod_mul(x % fp.mod, scales[c], fp.mod);
                if (value) row.push_back(mapped[c], value);
                if (row.nnz() > max_terms) { skipped = true; break; }
            }
            if (skipped || (row.nnz() && row.nnz() <= min_terms)) continue;
            normalize_mod(row, fp);
            if (seen) (*seen)[r] = 1;
            if (row.nnz()) emit(std::move(row));
        }
    }
    size_t residual_blocks() const { return (nrow() + 127) / 128; }
    template<class Emit> void residual(const MM& k, const field_t& fp, size_t block, Emit&& emit, thread_pool*) const {
        require(k.ncol == ncol() && block < residual_blocks(), "matrix residual dimensions");
        MM kt(k.ncol,k.nrow);
        for(size_t r=0;r<k.nrow;++r)for(auto [c,x]:k[r])kt[c].push_back(index(r),x);
        for (size_t r = block * 128; r < std::min(nrow(), block * 128 + 128); ++r) {
            std::vector<ulong> values(k.nrow);
            for (auto [c,x] : matrix[r]) {
                if (x.den() % fp.mod == 0) throw certified_kernel::BadPrime("matrix residual denominator");
                ulong a = x % fp.mod;
                for (auto [j,y] : kt[c]) values[j]=nmod_add(values[j],nmod_mul(a,y,fp.mod),fp.mod);
            }
            MV row;for(size_t j=0;j<values.size();++j)if(values[j])row.push_back(index(j),values[j]);
            if (row.nnz()) emit(std::move(row));
        }
    }
    int_t integral_residual_bound(const Mat& k) const {
        int_t maximum = 0;
        for (const auto& row : k.rows) for (auto [c,x] : row) {
            require(x.is_integer(), "matrix certificate requires primitive integral rows");
            auto h = x.num().abs(); if (h > maximum) maximum = std::move(h);
        }
        return height * maximum;
    }
};

inline Mat solve(const Mat& a, rref_option_t opt, const Options& options = {}) {
    MatrixRows rows(a); return solve(stack(rows), opt, options);
}

struct AffineSolution {
    bool consistent;
    Mat particular; // one row, or zero rows if inconsistent
    Mat directions; // all homogeneous directions
};

// Homogenize A x = b with one extra coordinate t, then intersect with t=1.
// There is no such affine section exactly when the system is inconsistent.
inline AffineSolution solve_affine(const Mat& a, const std::vector<Q>& b, rref_option_t opt, const Options& options = {}) {
    require(b.size() == a.nrow, "affine right-hand side dimensions");
    Mat augmented(a.nrow, a.ncol + 1);
    for (size_t r = 0; r < a.nrow; ++r) {
        augmented[r] = a[r]; if (b[r] != 0) augmented[r].push_back(index(a.ncol), -b[r]);
    }
    auto k = solve(augmented, opt, options);
    size_t pivot = k.nrow;
    for (size_t r = 0; r < k.nrow; ++r) if (k[r].find(index(a.ncol))) { pivot = r; break; }
    if (pivot == k.nrow) {
        k.ncol = a.ncol;
        return {false, Mat(0, a.ncol), std::move(k)};
    }
    Q scale = *k[pivot].find(index(a.ncol)); Mat point(1, a.ncol), directions(0, a.ncol);
    for (auto [c,x] : k[pivot]) if (size_t(c) != a.ncol) point[0].push_back(c, x / scale);
    for (size_t r = 0; r < k.nrow; ++r) if (r != pivot) {
        Vec row;
        for (auto [c,x] : k[r]) if (size_t(c) != a.ncol) row.push_back(c,x);
        if (auto t = k[r].find(index(a.ncol))) for (auto [c,x] : point[0]) row.push_back(c, -(*t) * x);
        normalize(row); directions.rows.push_back(std::move(row));
    }
    directions.nrow = directions.rows.size();
    return {true, std::move(point), std::move(directions)};
}
} // namespace staged_kernel
