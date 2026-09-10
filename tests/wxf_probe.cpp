#include "bootstrap.hpp"
// Tiny in-process reader for sanitizer/fuzz tests. Every file must either
// construct a valid sparse tensor or throw a normal exception.
int main(int argc, char** argv) {
    field_t field(FIELD_QQ);
    for (int i = 1; i < argc; ++i) {
        try {
            auto t = sparse_tensor_read_wxf<rat_t, int32_t>(std::filesystem::path(argv[i]), field);
            if (t.rank() == 2) {
                auto m = sparse_mat_read_wxf<rat_t, int32_t>(std::filesystem::path(argv[i]), field);
                if (m.nrow != (int32_t)t.dim(0) || m.ncol != (int32_t)t.dim(1) || m.nnz() != t.nnz())
                    throw std::runtime_error("matrix/tensor reader disagreement");
            }
            std::cout << "OK " << i << " " << t.rank() << " " << t.nnz() << "\n";
        } catch (const std::exception& e) {
            std::cout << "REJECT " << i << " " << e.what() << "\n";
        }
    }
}
