// Probe 2 on the straps of a committed chip: any strap span over qubit metal, any foot on island metal, the metal
// area each strap spans per crossed island (coupler, bus, line), parallel-plate capacitance at 3 um bridge height,
// and the summed strap capacitance per coupler/bus island.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
using namespace qpuforge;
using namespace qpuforge::layout;
using detail::kStrapW;

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "AF-17";
    ChipDesign D; std::string cells;
    if (which == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else if (which == "BF-49") { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    else { D = gamma1(7); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    const double eps0 = 8.8541878128e-12, h = 3e-6;
    int over_q = 0, foot_on_island = 0;
    std::map<std::string, double> area_by_island;     // um^2 of island metal under strap spans (between the feet)
    double max_single = 0; std::string max_single_name;
    for (const auto& s : L.straps) {
        const bool horiz = s.span.w() > s.span.h();
        const Rect mid = horiz ? Rect{std::min(s.foot_a.x1, s.foot_b.x1) , s.span.y0, std::max(s.foot_a.x0, s.foot_b.x0), s.span.y1}
                               : Rect{s.span.x0, std::min(s.foot_a.y1, s.foot_b.y1), s.span.x1, std::max(s.foot_a.y0, s.foot_b.y0)};
        for (const auto& [name, rs] : L.islands) {
            double a = 0;
            for (const auto& r : rs) {
                if (r.intersects(s.span) && name[0] == 'Q') ++over_q;
                if (r.intersects(s.foot_a) || r.intersects(s.foot_b)) ++foot_on_island;
                if (r.intersects(mid)) a += (std::min(r.x1, mid.x1) - std::max(r.x0, mid.x0)) * (std::min(r.y1, mid.y1) - std::max(r.y0, mid.y0));
            }
            if (a > 0) { area_by_island[name] += a; if (a > max_single) { max_single = a; max_single_name = name; } }
        }
    }
    std::printf("%s: %zu straps; spans over qubit metal: %d; feet touching island metal: %d\n", which.c_str(), L.straps.size(), over_q, foot_on_island);
    // per island type
    double worst_c = 0, worst_b = 0; std::string wc, wb; int nC = 0, nB = 0;
    for (const auto& [n, a] : area_by_island) {
        const double cf = eps0 * a * 1e-12 / h * 1e15;
        if (n[0] == 'C') { ++nC; if (cf > worst_c) { worst_c = cf; wc = n; } }
        if (n[0] == 'B') { ++nB; if (cf > worst_b) { worst_b = cf; wb = n; } }
    }
    std::printf("  largest single strap-over-island overlap: %.0f um^2 over %s (%.2f fF parallel-plate at 3 um)\n", max_single, max_single_name.c_str(),
                eps0 * max_single * 1e-12 / h * 1e15);
    std::printf("  couplers under straps: %d, worst summed overlap %s %.2f fF; buses under straps: %d, worst %s %.2f fF (parallel plate, no fringe)\n",
                nC, wc.c_str(), worst_c, nB, wb.c_str(), worst_b);
    int shown = 0;
    for (const auto& [n, a] : area_by_island) if ((n[0] == 'C' || n[0] == 'B') && shown < 6 && eps0 * a * 1e-12 / h * 1e15 > 0.5) {
        std::printf("    %s: %.0f um^2 -> %.2f fF\n", n.c_str(), a, eps0 * a * 1e-12 / h * 1e15); ++shown; }
    return 0;
}
