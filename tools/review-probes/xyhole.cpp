#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main() {
    const auto D = alpha1_v11(3);
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    const auto M = detail::metal_map(L.outline, L.etch);
    for (int id : {1, 6}) {
        double cx = 0, cy = 0; for (const auto& q : D.qubits) if (q.id == id) { cx = q.x_um; cy = q.y_um; }
        std::printf("Q%d: metal component along the line y = -280 (rel), x = 250..335:", id);
        for (double x = 250; x <= 335; x += 5) std::printf(" %g:%d", x, M.at(cx + x, cy - 280));
        std::printf("\n     along x = +280 (rel), y = -250..-335:");
        for (double y = -250; y >= -335; y -= 5) std::printf(" %g:%d", y, M.at(cx + 280, cy + y));
        std::printf("\n     strip between XY pad hole (x1 = +305) and coupler pad hole (x0 = +307) at y = -280: at(+306) = %d\n", M.at(cx + 306, cy - 280));
    }
}
