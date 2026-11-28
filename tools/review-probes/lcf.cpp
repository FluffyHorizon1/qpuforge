#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main() {
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet S = cells_from_json(Json::parse(ss.str()));
    const MeshOptions mo{3.0, 50.0, 2.0};
    const auto c55 = line_coupling(S.meas, S.cg, S.pitch, 5.5, mo), c599 = line_coupling(S.meas, S.cg, S.pitch, 5.99, mo);
    std::printf("C_flux_eff weighted at 5.50 GHz: %.5f fF; at 5.99 GHz: %.5f fF (ratio %.3f); open %.3f fF; C_xy %.4f fF\n", c55.C_flux_fF, c599.C_flux_fF,
                c599.C_flux_fF / c55.C_flux_fF, c55.C_flux_open_fF, c55.C_xy_fF);
    const double w = 2 * PI * 5.99e9;
    auto t1 = [&](double C, double Cs) { return Cs * 1e-15 / (w * w * C * 1e-15 * C * 1e-15 * 50) * 1e6; };
    const double a = t1(c55.C_flux_fF, c55.C_sigma_fF), b = t1(c599.C_flux_fF, c599.C_sigma_fF);
    std::printf("Q16 (5.99 GHz, design T1 65.5 us) flux T1 limit: as line_budget computes it %.0f us = %.1fx; with the 5.99 GHz weighting %.0f us = %.1fx\n", a, a / 65.5, b, b / 65.5);
}
