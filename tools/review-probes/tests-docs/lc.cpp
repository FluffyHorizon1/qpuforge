// line_coupling of the measure / data cells of a cell file, and line_budget of a chip, at a given mesh
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main(int argc, char** argv) {
    const std::string cells = argv[1];
    const double h = std::atof(argv[2]), r = std::atof(argv[3]);   // MeshOptions{h, 50, r}
    const int chip = argc > 4 ? std::atoi(argv[4]) : -1;            // -1: cells only; 0 AF-17, 1 BF-49, 2 GF-97: line_budget
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet S = cells_from_json(Json::parse(ss.str()));
    const MeshOptions mo{h, 50.0, r};
    const auto t0 = std::chrono::steady_clock::now();
    if (chip < 0) {
        const auto c = line_coupling(S.meas, S.cg, S.pitch, 5.5, mo);
        std::printf("measure cell @5.5 GHz: C_xy %.4f fF, C_flux_eff %.4f fF, C_flux_open %.3f fF, C_sigma %.2f fF, T1_xy %.3g us, T1_flux %.3g us (%.0f s)\n",
                    c.C_xy_fF, c.C_flux_fF, c.C_flux_open_fF, c.C_sigma_fF, c.T1_xy_us, c.T1_flux_us,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        return 0;
    }
    const ChipDesign D = chip == 0 ? alpha1_v11(3) : chip == 1 ? beta1(5) : gamma1(7);
    const auto R = analyze(D);
    const auto B = line_budget(D, R, S, mo);
    std::printf("%s\n(%.0f s)\n", line_budget_summary(B).c_str(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    // per-qubit: f01, design T1, and the flux T1 limit with C_flux scaled to the qubit frequency (C_eff ~ sin(beta z) ~ f)
    double worst_fixed = 1e300, worst_scaled = 1e300; int wq_fixed = -1, wq_scaled = -1;
    double fmin = 1e9, fmax = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const LineCoupling& c = D.qubits[i].role != "data" ? B.meas : B.data;
        const double f = R.q[i].f01, w = 2 * 3.14159265358979 * f * 1e9;
        fmin = std::min(fmin, f); fmax = std::max(fmax, f);
        auto t1 = [&](double cc) { return c.C_sigma_fF * 1e-15 / (w * w * cc * 1e-15 * cc * 1e-15 * 50) * 1e6; };
        const double fixed = t1(c.C_flux_fF) / R.q[i].T1_us, scaled = t1(c.C_flux_fF * f / 5.5) / R.q[i].T1_us;
        if (fixed < worst_fixed) { worst_fixed = fixed; wq_fixed = D.qubits[i].id; }
        if (scaled < worst_scaled) { worst_scaled = scaled; wq_scaled = D.qubits[i].id; }
        if (argc > 5) std::printf("  Q%d %s f01 %.3f GHz T1 %.1f us: flux ratio %.2fx (C at 5.5 GHz) vs %.2fx (C scaled to f01)\n", D.qubits[i].id,
                                  D.qubits[i].role.c_str(), f, R.q[i].T1_us, fixed, scaled);
    }
    std::printf("f01 range %.3f-%.3f GHz; worst flux ratio %.2fx (Q%d) as coded, %.2fx (Q%d) with C_flux scaled to f01\n", fmin, fmax,
                worst_fixed, wq_fixed, worst_scaled, wq_scaled);
}
