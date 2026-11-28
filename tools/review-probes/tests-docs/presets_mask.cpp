// Reproduces the mask / connectivity assertions of test_layout_presets (tests/test_qpuforge.cpp lines 1915-1923) for
// the three committed chips, without the (expensive) cell extraction, which does not depend on the control lines or
// the junction regions. Compiled against an include directory given with -I, so it can be run on mutated copies.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;

int main(int argc, char** argv) {
    const std::string src = argc > 1 ? argv[1] : "/home/user/wt-fab";
    const int only = argc > 2 ? std::atoi(argv[2]) : -1;
    struct Chip { const char* name; int which; const char* cells; };
    const std::vector<Chip> chips{{"AF-17", 0, "designs/alpha1-v11/layout_cells.json"},
                                  {"BF-49", 1, "designs/beta1/layout_cells.json"},
                                  {"GF-97", 2, "designs/beta1/layout_cells.json"}};
    int fails = 0;
    for (int ci = 0; ci < int(chips.size()); ++ci) {
        if (only >= 0 && ci != only) continue;
        const auto& c0 = chips[std::size_t(ci)];
        struct { const char* name; ChipDesign D; const char* cells; } c{c0.name, c0.which == 0 ? alpha1_v11(3) : c0.which == 1 ? beta1(5) : gamma1(7), c0.cells};
        std::ifstream f(src + "/" + c.cells);
        std::stringstream ss; ss << f.rdbuf();
        const CellSet C = cells_from_json(Json::parse(ss.str()));
        const auto t0 = std::chrono::steady_clock::now();
        const auto L = assemble_chip(c.D, C);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const bool placement = L.problems.empty();
        const bool drc = L.drc_count("width") == 0 && L.drc_count("space") == L.conn.fills_reverted;
        const bool conn = L.conn.clean();
        const bool counts = L.conn.junctions == int(c.D.qubits.size() + c.D.couplings.size()) && L.conn.flux_lines == int(c.D.qubits.size());
        std::printf("%s (%.1f s): placement %s; DRC %d+%d (%s); conn.clean %s; counts %s (jj %d, flux %d)\n  %s\n", c.name, sec,
                    placement ? "ok" : L.problems.front().c_str(), L.drc_count("width"), L.drc_count("space"), drc ? "ok" : "FAIL",
                    conn ? "yes" : "NO", counts ? "ok" : "FAIL", L.conn.junctions, L.conn.flux_lines, connectivity_summary(L).c_str());
        const bool pass = placement && drc && conn && counts;
        std::printf("  => layout_presets mask/connectivity assertions for %s: %s\n", c.name, pass ? "PASS" : "FAIL");
        fails += !pass;
    }
    return fails;
}
