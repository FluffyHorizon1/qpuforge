// Probe 3: the SQUID pocket as its own conductor in the PR's line_coupling cell; galvanic T1 through the pocket strap;
// Biot-Savart mutual inductance of the flux-line current loop to the SQUID region.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;
using cd = std::complex<double>;

int main(int argc, char** argv) {
    const double hmin = argc > 1 ? std::atof(argv[1]) : 3.0;
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet S = cells_from_json(Json::parse(ss.str()));
    const XmonParams& q = S.meas; const CouplerGeom& c = S.cg; const double p = S.pitch;
    // --- the line_coupling cell, verbatim
    std::vector<Conductor> cs;
    std::vector<Rect> holes;
    Conductor Q; Q.name = "qubit"; xmon(q, 0, 0, Q.rects, holes); cs.push_back(Q);
    for (int k = 0; k < 4; ++k) {
        auto E = coupler_and_bus(q, c, p, k, 0, 0, false);
        Conductor C; C.name = "coupler" + std::to_string(k); C.rects = E.coupler; cs.push_back(C);
        if (!E.bus.empty()) { Conductor B; B.name = "bus" + std::to_string(k); B.rects = E.bus; cs.push_back(B); }
        holes.insert(holes.end(), E.holes.begin(), E.holes.end());
    }
    Conductor R; R.name = "readout"; R.rects = readout_claw(q, 0, 0, 170, &holes); cs.push_back(R);
    std::vector<Rect> fm, xm;
    qubit_control_lines(0, 0, fm, holes, xm, holes, false);
    const int n = int(std::lround(95 / 5.0));
    Conductor F; F.name = "flux";
    for (int k = 0; k < int(fm.size()); ++k) if (k < n - 1 || k > n) F.rects.push_back(fm[std::size_t(k)]);
    cs.push_back(F);
    Conductor X; X.name = "xy"; X.rects = xm; cs.push_back(X);
    Conductor G; G.name = "ground"; G.ground = true;
    const double ext = q.L + c.g_c + c.t_w + 200 + 100;
    G.outline = {-ext, -ext, ext, ext}; G.holes = cell_ground_holes(holes, cs);
    // --- find the pocket that holds the SQUID ground lead (the tongue) in this cell's metal
    std::vector<Rect> metal; for (const auto& cc : cs) metal.insert(metal.end(), cc.rects.begin(), cc.rects.end());
    const auto etch = detail::etch_of(G.holes, metal);
    const auto M = detail::metal_map(G.outline, etch);
    const double s0 = -q.w_arm / 2 - q.squid_stub;
    const double tx = s0 - kSquidGap - kTongueW / 2, ty = -q.w_arm / 2 - kSquidGap - 5;   // inside the tongue
    const int pk = M.at(tx, ty), plane = M.at(G.outline.x0 + 0.5, G.outline.y0 + 0.5);
    Conductor P; P.name = "pocket";
    double parea = 0;
    for (std::size_t i = 0; i < M.comp.size(); ++i) if (M.comp[i] == pk) { P.rects.push_back({M.ix0[i], M.iy0[i], M.ix1[i], M.iy1[i]}); parea += P.rects.back().area(); }
    const Rect pb = bbox_of(P.rects);
    std::printf("cell: tongue at (%.1f, %.1f) is metal component %d (cell plane %d); pocket area %.0f um^2, bbox x[%.1f,%.1f] y[%.1f,%.1f], %zu rects\n", tx, ty, pk, plane, parea,
                pb.x0, pb.x1, pb.y0, pb.y1, P.rects.size());
    // the SQUID region must bridge stub -> pocket
    Rect jjr; { std::vector<Rect> m2, h2, j2; xmon(q, 0, 0, m2, h2, &j2); jjr = j2.front(); }
    std::printf("SQUID region x[%.1f,%.1f] y[%.1f,%.1f]: left end comp %d (pocket %d), right end comp %d\n", jjr.x0, jjr.x1, jjr.y0, jjr.y1, M.at(jjr.x0 - 0.5, jjr.cy()), pk, M.at(jjr.x1 + 0.5, jjr.cy()));
    G.holes.insert(G.holes.end(), P.rects.begin(), P.rects.end());
    cs.push_back(P);
    cs.push_back(G);
    const MeshOptions mo{hmin, 50.0, 2.0};
    const auto Cm = extract_capacitance(cs, 11.45, mo);
    std::printf("extraction: %d panels, %.0f s, asymmetry %.2e\n", Cm.panels, Cm.seconds, Cm.asymmetry);
    const double Cqq = Cm.self("qubit"), Cqp = Cm.mutual("qubit", "pocket"), Cpp = Cm.self("pocket"), Cqg = Cm.mutual("qubit", "ground");
    const double Cpf = Cm.mutual("pocket", "flux"), Cqf = Cm.mutual("qubit", "flux"), Cqx = Cm.mutual("qubit", "xy");
    std::printf("C_sigma %.2f fF; qubit-pocket %.2f fF (%.1f%% of C_sigma); qubit-plane %.2f fF; pocket self %.2f fF; pocket-flux %.2f fF; qubit-flux (open, truncated) %.3f fF; qubit-xy %.4f fF\n",
                Cqq, Cqp, 100 * Cqp / Cqq, Cqg, Cpp, Cpf, Cqf, Cqx);
    // --- galvanic T1: island q -- L_J -- pocket p -- (L_s || Z0) -- plane; capacitances from the extraction
    const double Phi0 = 2.067833848e-15, hP = 6.62607015e-34, Z0 = 50;
    const std::vector<std::pair<double, double>> qubits{{5.110, 23.73}, {5.990, 23.73}};   // placeholders replaced below
    (void)qubits;
    auto T1_galv = [&](double EJ_GHz, double Ls) {
        const double LJ = std::pow(Phi0 / (2 * PI), 2) / (hP * EJ_GHz * 1e9);
        const double Cq = Cqq * 1e-15, Cp = Cpp * 1e-15, Cx = Cqp * 1e-15;
        auto det = [&](cd w) {
            const cd j(0, 1);
            const cd yL = 1.0 / (j * w * LJ);
            const cd a11 = j * w * Cq + yL, a12 = -j * w * Cx - yL;
            const cd a22 = j * w * Cp + yL + 1.0 / (j * w * Ls) + 1.0 / Z0;
            return a11 * a22 - a12 * a12;
        };
        cd w = 1.0 / std::sqrt(LJ * (Cq - Cx * Cx / Cp));
        for (int it = 0; it < 60; ++it) { const cd h = w * 1e-7; const cd d = det(w), dp = (det(w + h) - det(w - h)) / (2.0 * h); w -= d / dp; }
        return std::make_pair(1.0 / (2 * w.imag()) * 1e6, w.real() / (2 * PI) / 1e9);   // energy T1 (us), f (GHz)
    };
    for (double EJ : {17.72, 23.73}) {
        const double LJ = std::pow(Phi0 / (2 * PI), 2) / (hP * EJ * 1e9);
        std::printf("E_J %.2f GHz (L_J %.2f nH):", EJ, LJ * 1e9);
        for (double Ls : {1e-12, 2e-12, 5e-12, 20e-12, 50e-12}) { const auto [t1, fp] = T1_galv(EJ, Ls); std::printf("  L_s %2.0f pH -> T1 %.3g us (f_p %.2f GHz)", Ls * 1e12, t1, fp); }
        std::printf("\n");
    }
    // partial self-inductance of a flat strip (Grover): l = bridge length, w = 10 um, t = 0.3 um
    auto Lstrip = [](double l, double w, double t) { return 2e-7 * l * (std::log(2 * l / (w + t)) + 0.5 + 0.2235 * (w + t) / l); };
    for (double l : {36e-6, 98e-6, 208e-6}) std::printf("strip partial self-inductance, %.0f um x 10 um: %.0f pH\n", l * 1e6, Lstrip(l, 10e-6, 0.3e-6) * 1e12);
    // --- Biot-Savart: flux-line loop (stair centre line -> short -> pocket strap foot -> bridge -> plane foot -> back to the
    // stair start) through the SQUID region (all in-plane filaments; B_z integrated over the SQUID rectangle)
    const std::vector<std::pair<double, double>> loop{{-156, -156}, {-61, -61}, {-51, -76.5}, {-159, -76.5}, {-156, -156}};
    auto Bz = [&](double x, double y) {
        double b = 0; const int ns = 4000;
        for (std::size_t s = 0; s + 1 < loop.size(); ++s) {
            const double ax = loop[s].first, ay = loop[s].second, bx = loop[s + 1].first, by = loop[s + 1].second;
            for (int k = 0; k < ns; ++k) {
                const double t0 = (k + 0.5) / ns, px = ax + (bx - ax) * t0, py = ay + (by - ay) * t0;
                const double dlx = (bx - ax) / ns, dly = (by - ay) / ns, rx = x - px, ry = y - py, r2 = rx * rx + ry * ry;
                b += (dlx * ry - dly * rx) / std::pow(r2, 1.5);   // (dl x r)_z / r^3, um^-1
            }
        }
        return 1e-7 * b * 1e6;   // T per A (mu0/4pi = 1e-7; um^-1 -> m^-1)
    };
    double phi = 0; const int g = 40;
    for (int i = 0; i < g; ++i) for (int j = 0; j < g; ++j) {
        const double x = jjr.x0 + (i + 0.5) * jjr.w() / g, y = jjr.y0 + (j + 0.5) * jjr.h() / g;
        phi += Bz(x, y) * (jjr.w() / g) * (jjr.h() / g) * 1e-12;
    }
    std::printf("mutual inductance flux-line loop -> SQUID region (%.0f um^2): %.3f pH; current for Phi0/2: %.1f mA\n", jjr.area(), std::abs(phi) * 1e12,
                Phi0 / 2 / std::abs(phi) * 1e3);
    return 0;
}
