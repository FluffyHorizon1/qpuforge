// layout_gates: the mask gates of one chip from a cell file, with the time each step takes. A development and review
// driver for docs/design/mask-connectivity.md: it prints what the review programs print about the mask (DRC, metal
// connectivity and ground bumps, the extracted cells against their targets, the control-line budget) without the rest
// of the review, and the flux coupling of any one line to any one SQUID at a chosen mesh.
//
//   layout_gates gates <chip> <cells.json> [--flux-mesh H] [--cross-mesh H] [--window-mesh H] [--no-crosstalk] [--bump-sensitivity] [--gds <out.gds>]
//   layout_gates flux  <chip> <cells.json> <owner> <line> [--flux-mesh H]
//
// <chip> is a preset as the reviews use it: alpha1-v11 (AF-17, distance 3), beta1 (BF-49, distance 5) or gamma1
// (GF-97, distance 7). A mesh is given by its edge panel size H in um (largest panel 50 um, growth 1.6 below 2 um and
// 2.0 from 2 um). Defaults are the reviews' settings: cells and control windows 1.5 um, flux solves 1.5 um
// (kFluxMesh), crosstalk solves 2.5 um (kFluxCrossMesh). --bump-sensitivity evaluates the budget again (no new solve)
// with the bump inductance halved and doubled, and with all of a qubit's current and half of a neighbour's taken through
// the bumps (quadrant fractions 1 and 0.5 instead of 0.75 and 0.25).
// Exit status: 0 when every gate passes (gates) or the solve ran (flux); 1 when a gate fails or the solve did not
// run; 2 on a usage error or a file that cannot be read.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace qpuforge;
using namespace qpuforge::layout;

namespace {
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }
MeshOptions mesh_of(double h) { return MeshOptions{h, 50.0, h < 2 ? 1.6 : 2.0}; }

int usage(const char* why) {
    std::fprintf(stderr,
                 "layout_gates: %s\n"
                 "usage: layout_gates gates <chip> <cells.json> [--flux-mesh H] [--cross-mesh H] [--window-mesh H] [--no-crosstalk] [--bump-sensitivity] [--gds <out.gds>]\n"
                 "       layout_gates flux  <chip> <cells.json> <owner> <line> [--flux-mesh H]\n"
                 "  <chip>: alpha1-v11 | beta1 | gamma1 (the presets of the reviews); H: edge panel size in um\n",
                 why);
    return 2;
}
bool chip_of(const std::string& s, ChipDesign& D) {
    if (s == "alpha1-v11") D = alpha1_v11(3);
    else if (s == "beta1") D = beta1(5);
    else if (s == "gamma1") D = gamma1(7);
    else return false;
    return true;
}
void print_flux(const char* what, const std::string& owner, const std::string& line, const FluxCoupling& F, double h) {
    std::printf("  %s %s <- %s (flux mesh %.2f um): %s M %.4f pH, M_cm %.4f pH (signed %.4f), I(Phi0/2) %.4g mA; %d apertures, %d panels, %.1f s%s%s\n", what,
                owner.c_str(), line.c_str(), h, F.ok ? "" : "NOT EVALUATED", F.M_pH, F.M_cm_pH, F.M_cm_signed_pH, std::min(F.I_half_mA, 1e9), F.apertures, F.panels,
                F.seconds, F.note.empty() ? "" : "; ", F.note.c_str());
}
}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const auto t_start = Clock::now();
    std::vector<std::string> a(argv + 1, argv + argc);
    if (a.size() < 3) return usage("too few arguments");
    const std::string cmd = a[0];
    if (cmd != "gates" && cmd != "flux") return usage(("unknown command '" + cmd + "'").c_str());
    ChipDesign D;
    if (!chip_of(a[1], D)) return usage(("unknown chip '" + a[1] + "'").c_str());
    std::ifstream f(a[2]);
    if (!f) { std::fprintf(stderr, "layout_gates: cannot open %s\n", a[2].c_str()); return 2; }
    std::stringstream ss; ss << f.rdbuf();
    CellSet C;
    try { C = cells_from_json(Json::parse(ss.str())); } catch (const std::exception& e) { std::fprintf(stderr, "layout_gates: %s: %s\n", a[2].c_str(), e.what()); return 2; }
    std::size_t first_opt = cmd == "flux" ? 5 : 3;
    if (cmd == "flux" && a.size() < 5) return usage("flux needs <owner> <line>");
    double h_flux = kFluxMesh.h_min, h_cross = kFluxCrossMesh.h_min, h_window = 1.5;
    bool crosstalk = true, bump_sensitivity = false;
    std::string gds;
    for (std::size_t i = first_opt; i < a.size(); ++i) {
        auto num = [&](double& v) {
            if (i + 1 >= a.size()) return false;
            char* end = nullptr; v = std::strtod(a[i + 1].c_str(), &end);
            ++i;
            return end && *end == 0 && v > 0.2 && v < 20;
        };
        if (a[i] == "--flux-mesh") { if (!num(h_flux)) return usage("--flux-mesh needs a size in um (0.2 .. 20)"); }
        else if (a[i] == "--cross-mesh" && cmd == "gates") { if (!num(h_cross)) return usage("--cross-mesh needs a size in um (0.2 .. 20)"); }
        else if (a[i] == "--window-mesh" && cmd == "gates") { if (!num(h_window)) return usage("--window-mesh needs a size in um (0.2 .. 20)"); }
        else if (a[i] == "--no-crosstalk" && cmd == "gates") crosstalk = false;
        else if (a[i] == "--bump-sensitivity" && cmd == "gates") bump_sensitivity = true;
        else if (a[i] == "--gds" && cmd == "gates" && i + 1 < a.size()) gds = a[++i];
        else return usage(("unknown option '" + a[i] + "'").c_str());
    }
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif

    auto t0 = Clock::now();
    const ChipLayout L = assemble_chip(D, C);
    const double t_assemble = since(t0);
    std::printf("%s from %s: control lines %s, %zu qubits, %zu couplings\n", D.name.c_str(), a[2].c_str(), C.separate_xy ? "separate (XY lines)" : "combined (flux and drive on one line)",
                D.qubits.size(), D.couplings.size());
    std::printf("  %s\n  %s\n", drc_summary(L).c_str(), connectivity_summary(L).c_str());

    if (cmd == "flux") {
        t0 = Clock::now();
        const FluxCoupling F = flux_coupling(L, a[3], a[4], mesh_of(h_flux));
        print_flux("flux coupling", a[3], a[4], F, h_flux);
        std::printf("  time (%d threads): assemble %.1f s, flux solve %.1f s, total %.1f s\n", threads, t_assemble, since(t0), since(t_start));
        return F.ok ? 0 : 1;
    }

    if (!gds.empty()) write_gds(L, gds, D.name);
    t0 = Clock::now();
    const ChipReport R = analyze(D);
    const LayoutMetrics m = evaluate_cells(C, mesh_of(1.5));
    const LayoutTargets T = targets_from_design(D, R);
    const PostLayout P = post_layout(D, R, m);
    const double t_cells = since(t0);
    auto rel = [](double x, double t) { return std::abs(x / t - 1); };
    const double ec = std::max({rel(m.EC_m, T.EC_m), rel(m.EC_d, T.EC_d), rel(m.EC_c, T.EC_c)});
    const double g = std::max({rel(m.G_mc, T.G_mc), rel(m.G_cd, T.G_cd), rel(m.G_md, T.G_md)});
    std::printf("  cells (mesh 1.5 um) against the design targets: E_C within %.2f%%, couplings within %.2f%%, readout C within %.2f%%; parked exchange re-nulled to %.1e MHz "
                "(couplers parked at %.3f GHz)\n", 100 * ec, 100 * g, 100 * rel(m.C_qr_fF, T.C_qr_fF), P.park_gap_MHz, P.park_fc_GHz);

    t0 = Clock::now();
    LineBudget B = line_budget(D, R, L, m, P.park_fc_GHz, mesh_of(h_window), 50, mesh_of(h_flux));
    const double t_budget = since(t0);
    double t_cross = 0;
    if (crosstalk) { t0 = Clock::now(); line_budget_crosstalk(B, D, R, L, m, mesh_of(h_cross)); t_cross = since(t0); }
    print_flux("measure qubit", "Q" + std::to_string(B.meas.qubit), "FQ" + std::to_string(B.meas.qubit), B.meas.flux, h_flux);
    print_flux("data qubit   ", "Q" + std::to_string(B.data.qubit), "FQ" + std::to_string(B.data.qubit), B.data.flux, h_flux);
    print_flux("coupler      ", "C" + B.coupler, "FC" + B.coupler, B.coupler_flux, h_flux);
    std::printf("  control windows (mesh %.2f um): Q%d %d panels %.1f s, Q%d %d panels %.1f s\n", h_window, B.meas.qubit, B.meas.win.C.panels, B.meas.win.C.seconds,
                B.data.qubit, B.data.win.C.panels, B.data.win.C.seconds);
    std::printf("  %s\n", line_budget_summary(B).c_str());

    if (bump_sensitivity) {
        struct Case { const char* what; double bump_pH, own, neighbour; };
        for (const Case& k : {Case{"bump inductance halved", 0.5 * kBumpInductance_pH, kBumpOwnFraction, kBumpNeighbourFraction},
                              Case{"bump inductance doubled", 2 * kBumpInductance_pH, kBumpOwnFraction, kBumpNeighbourFraction},
                              Case{"quadrant fractions 1 / 0.5", kBumpInductance_pH, 1.0, 0.5}}) {
            LineBudget S = B;
            S.bump_pH = k.bump_pH; S.bump_own_fraction = k.own; S.bump_neighbour_fraction = k.neighbour;
            line_budget_evaluate(S, D, R, L, m);
            std::printf("  sensitivity, %s (%.1f pH per bump, fractions %.2f / %.2f): bump term up to %.3f pH (coupler %.3f pH); common-mode T1 %.1fx (Q%d), drive T1 %.1fx "
                        "(Q%d), coupler T1 %.1fx; drive on a neighbour %.1f%% (%s); 25 MHz Rabi at %.1f dBm -> drive %s, flux %s\n", k.what, S.bump_pH, S.bump_own_fraction,
                        S.bump_neighbour_fraction, S.max_M_bump_pH, S.M_bump_coupler_pH, std::min(S.min_ratio_flux_cm, 9999.0), S.worst_flux_cm_qubit,
                        std::min(S.min_ratio_drive, 9999.0), S.worst_drive_qubit, std::min(S.ratio_coupler, 9999.0), 100 * std::min(S.max_drive_crosstalk, 99.0),
                        S.worst_crosstalk.c_str(), S.drive_power_dBm, S.drive_ok() ? "ok" : "FAIL", S.flux_ok() ? "ok" : "FAIL");
        }
    }
    const bool drc = L.drc_clean(), conn = L.conn.nets_clean(), ground = L.conn.ground_ok(), bumps = L.conn.bumps_ok();
    const bool pass = drc && conn && ground && B.ok();
    std::printf("  gates: mask DRC and placement %s, nets %s, ground quality %s (ground bumps %s), control-line drive %s, flux lines %s%s -> %s\n", drc ? "ok" : "FAIL",
                conn ? "ok" : "FAIL", ground ? "ok" : "FAIL", bumps ? "ok" : "FAIL", B.drive_ok() ? "ok" : "FAIL", B.flux_ok() ? "ok" : "FAIL",
                crosstalk ? "" : " (crosstalk not evaluated: the drive gate cannot pass)", pass ? "PASS" : "FAIL");
    std::printf("  time (%d threads): assemble and check the mask %.1f s, cells and targets %.1f s, line budget %.1f s (windows %.1f s, flux solves %.1f s), crosstalk %.1f s, "
                "total %.1f s\n", threads, t_assemble, t_cells, t_budget, B.meas.win.C.seconds + B.data.win.C.seconds,
                B.meas.flux.seconds + B.data.flux.seconds + B.coupler_flux.seconds, t_cross, since(t_start));
    return pass ? 0 : 1;
}
