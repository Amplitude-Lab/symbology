#pragma once
#include <cstddef>

namespace hexagon_adaptive_physics {
// This state survives a callback deadline. Only completed work advances it.
struct WorkCursor {
    std::size_t collinear_batches = 0;
    std::size_t completed_cycles = 0;
    bool ope_done = false;
};

template<class Needed, class WithinBudget, class ColDone, class MrkDone,
         class Collinear, class Mrk, class Ope>
void run_budgeted_physics(WorkCursor& cursor, std::size_t collinear_quota,
                         std::size_t maximum_cycles, std::size_t cycles_before_ope,
                         bool enable_ope, Needed needed, WithinBudget within_budget,
                         ColDone col_done, MrkDone mrk_done,
                         Collinear collinear, Mrk mrk, Ope ope) {
    std::size_t callback_cycles = 0;
    while (needed() && within_budget()) {
        // A due OPE block runs before starting another cycle, including when
        // the preceding MRK block used the remainder of the previous budget.
        if (enable_ope && !cursor.ope_done &&
            (cursor.completed_cycles >= cycles_before_ope || (col_done() && mrk_done()))) {
            ope();
            cursor.ope_done = true;
            continue;
        }
        if (callback_cycles >= maximum_cycles || (col_done() && mrk_done())) break;
        if (!col_done() && cursor.collinear_batches < collinear_quota) {
            collinear();
            ++cursor.collinear_batches;
            continue;
        }
        if (!mrk_done()) mrk();
        cursor.collinear_batches = 0;
        ++cursor.completed_cycles;
        ++callback_cycles;
    }
}
} // namespace hexagon_adaptive_physics
