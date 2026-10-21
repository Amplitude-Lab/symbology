#include "hexagon_physics_work_cursor.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using hexagon_adaptive_physics::WorkCursor;
using hexagon_adaptive_physics::run_budgeted_physics;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct Result { std::string trace; std::size_t col=0, mrk=0, ope=0, calls=0; };
static Result simulate(double budget, double col_cost, double mrk_cost,
                       std::size_t col_max, std::size_t mrk_max, bool enable_ope,
                       std::size_t stop_after=10000, std::size_t maximum_cycles=1) {
    WorkCursor cursor;
    Result result;
    while (result.trace.size()<stop_after &&
           (result.col<col_max || result.mrk<mrk_max || (enable_ope && !cursor.ope_done))) {
        require(++result.calls<1000, "Scheduler failed to make bounded progress");
        double elapsed=0;
        auto record=[&](char kind, double cost) { result.trace+=kind; elapsed+=cost; };
        run_budgeted_physics(cursor, 8, maximum_cycles, 2, enable_ope,
            [&]{return result.trace.size()<stop_after;}, [&]{return elapsed<budget;},
            [&]{return result.col>=col_max;}, [&]{return result.mrk>=mrk_max;},
            [&]{++result.col;record('C',col_cost);},
            [&]{++result.mrk;record('M',mrk_cost);},
            [&]{++result.ope;record('O',3*mrk_cost);});
    }
    return result;
}

int main() try {
    std::size_t cases=0;
    // Deadline placement must not change the ordered stream of admitted
    // operators, including stream exhaustion and early target attainment.
    for (auto col_max : {std::size_t(0),std::size_t(3),std::size_t(16),std::size_t(40)})
    for (auto mrk_max : {std::size_t(0),std::size_t(1),std::size_t(5)})
    for (bool ope : {false,true})
    for (auto stop : {std::size_t(1),std::size_t(9),std::size_t(10000)}) {
        auto reference=simulate(100000,1,1,col_max,mrk_max,ope,stop);
        for (double budget : {0.001,1.0,7.0,120.0})
        for (double col_cost : {1.0,21.0,125.0})
        for (double mrk_cost : {1.0,80.0,250.0})
        for (std::size_t max_cycles : {1,2}) {
            auto value=simulate(budget,col_cost,mrk_cost,col_max,mrk_max,ope,stop,max_cycles);
            require(value.trace==reference.trace, "Deadline changed operator order or replayed work");
            require(value.ope<=1, "OPE block replayed");
            ++cases;
        }
    }
    auto measured=simulate(120,21,80,40,5,true);
    require(measured.trace.substr(0,19)==std::string(8,'C')+"M"+std::string(8,'C')+"MO",
            "MRK or OPE starved at measured ten-loop collinear cost");
    std::cout << "PHYSICS_CURSOR_PROBE passed=1 cases=" << cases
              << " measured_trace=" << measured.trace << " callbacks=" << measured.calls << '\n';
} catch (const std::exception& error) {
    std::cerr << "PHYSICS_CURSOR_PROBE_ERROR " << error.what() << '\n';
    return 1;
}
