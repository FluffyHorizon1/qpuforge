// Probe: control-line budget. (1) layout::line_budget as the reviews call it (reproduce 23.8x / 22.4x and 11.4x / 9.7x);
// (2) per-qubit frequency of the worst qubits vs the 5.5 GHz at which the flux-line standing-wave weights are computed;
// (3) the same cell as line_coupling with the (-x,-y) pocket split out as its own conductor: the qubit-pocket
//     capacitance, and the flux-line effective coupling when the short is not ideal but closes through an inductance
//     L_s (the pocket's single ground strap): a lumped L_s at the short is an extra line length dz = L_s v / Z0, so every
//     panel weight becomes sin(beta (z + dz)) and the pocket itself (charge C_qp) gets sin(beta dz).
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;

static CellSet load(const std::string& p) { std::ifstream f(p); std::stringstream ss; ss << f.rdbuf(); return cells_from_json(Json::parse(ss.str())); }

struct Split { double C_open, C_eff0, C_qp, C_sigma, C_xy; std::vector<std::pair<double, double>> zq; };   // (z, q) per flux panel
static Split split_cell(const XmonParams& q, const CouplerGeom& c, double p, double f_GHz, const MeshOptions& mo) {
    std::vector<Conductor> cs; std::vector<Rect> holes;
    Conductor Q; Q.name = "qubit"; xmon(q, 0, 0, Q.rects, holes); cs.push_back(Q);
    for (int k = 0; k < 4; ++k) {
        auto E = coupler_and_bus(q, c, p, k, 0, 0, false);
        Conductor Cc; Cc.name = "coupler" + std::to_string(k); Cc.rects = E.coupler; cs.push_back(Cc);
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
    // pocket = the ground component holding the SQUID tongue
    std::vector<Rect> metal; for (const auto& cc : cs) metal.insert(metal.end(), cc.rects.begin(), cc.rects.end());
    const auto etch = detail::etch_of(G.holes, metal);
    const auto Mp = detail::metal_map(G.outline, etch);
    const double hw = q.w_arm / 2, s0 = -hw - q.squid_stub;
    const int pc = Mp.at(s0 - kSquidGap - kTongueW / 2, -hw - kSquidGap - 2);
    Conductor P; P.name = "pocket";
    double parea = 0;
    for (std::size_t i = 0; i < Mp.comp.size(); ++i) if (Mp.comp[i] == pc) { P.rects.push_back({Mp.ix0[i], Mp.iy0[i], Mp.ix1[i], Mp.iy1[i]}); parea += P.rects.back().area(); }
    const Rect pb = bbox_of(P.rects);
    std::printf("    pocket: %zu rects, %.0f um^2, bbox (%.1f,%.1f)-(%.1f,%.1f); plane comp %d, pocket comp %d\n", P.rects.size(), parea, pb.x0, pb.y0, pb.x1, pb.y1,
                Mp.at(G.outline.x0 + 1, G.outline.y0 + 1), pc);
    G.holes.insert(G.holes.end(), P.rects.begin(), P.rects.end());
    cs.push_back(P);
    cs.push_back(G);
    std::vector<PanelCharge> pq;
    const auto Mx = extract_capacitance(cs, 11.45, mo, 0, &pq);
    Split s{};
    const double w = 2 * 3.14159265358979 * f_GHz * 1e9, beta = w * std::sqrt(0.5 * (11.45 + 1)) / 299792458.0 * 1e-6, ls = 95 * std::sqrt(2.0);
    const int fi = Mx.index("flux");
    for (const auto& c2 : pq) {
        if (c2.cond != fi) continue;
        const double x = c2.r.cx(), y = c2.r.cy();
        const double z = (y < -161 && x < -140) ? ls + (-156 - y) : std::max(0.0, ((-61 - x) + (-61 - y)) / std::sqrt(2.0));
        s.C_open -= c2.q_fF; s.C_eff0 -= c2.q_fF * std::sin(beta * z); s.zq.push_back({z, -c2.q_fF});
    }
    s.C_qp = Mx.mutual("qubit", "pocket"); s.C_sigma = Mx.self("qubit"); s.C_xy = Mx.mutual("qubit", "xy");
    std::printf("    split extraction: %d panels, %.1f s\n", Mx.panels, Mx.seconds);
    return s;
}

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "AF-17";
    const bool af = which == "AF-17";
    const ChipDesign D = af ? alpha1_v11(3) : beta1(5);
    const CellSet S = load(af ? "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json" : "/home/user/wt-fab/designs/beta1/layout_cells.json");
    const MeshOptions mo{1.5, 50.0, 1.6};
    const auto R = analyze(D);
    const auto B = line_budget(D, R, S, mo);
    std::printf("%s: %s\n", which.c_str(), line_budget_summary(B).c_str());
    std::printf("  meas: C_xy %.4f C_flux %.4f (open %.3f) C_sigma %.2f | data: C_xy %.4f C_flux %.4f (open %.3f) C_sigma %.2f fF\n",
                B.meas.C_xy_fF, B.meas.C_flux_fF, B.meas.C_flux_open_fF, B.meas.C_sigma_fF, B.data.C_xy_fF, B.data.C_flux_fF, B.data.C_flux_open_fF, B.data.C_sigma_fF);
    std::printf("  cell file C_sigma: meas %.2f data %.2f fF\n", 0.0, 0.0);
    // (2) per-qubit f01 and the flux ratio if the sin(beta z) weights are taken at the qubit's own f01 (C_eff ~ f)
    double worst_f = 1e300; int wq = -1; double fmin = 1e9, fmax = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const auto& c = D.qubits[i].role != "data" ? B.meas : B.data;
        const double f = R.q[i].f01, wv = 2 * 3.14159265358979 * f * 1e9;
        const double Ceff = c.C_flux_fF * (std::sin(1.0) ? f / 5.5 : 1.0);     // small-angle: weights scale with f
        const double t1 = c.C_sigma_fF * 1e-15 / (wv * wv * Ceff * 1e-15 * Ceff * 1e-15 * 50) * 1e6;
        fmin = std::min(fmin, f); fmax = std::max(fmax, f);
        if (t1 / R.q[i].T1_us < worst_f) { worst_f = t1 / R.q[i].T1_us; wq = D.qubits[i].id; }
        if (D.qubits[i].id == B.worst_flux_qubit || D.qubits[i].id == wq)
            std::printf("  Q%d (%s): f01 %.3f GHz, design T1 %.1f us\n", D.qubits[i].id, D.qubits[i].role.c_str(), f, R.q[i].T1_us);
    }
    std::printf("  qubit f01 range %.3f-%.3f GHz; flux ratio with the weights at each qubit's f01: %.2fx (Q%d) vs reported %.2fx (Q%d)\n",
                fmin, fmax, worst_f, wq, B.min_ratio_flux, B.worst_flux_qubit);
    // (3) split pocket, inductive short
    const double v = 299792458.0 / std::sqrt(0.5 * (11.45 + 1)), Z0 = 50;
    for (int cell = 0; cell < 2; ++cell) {
        const XmonParams& q = cell ? S.data : S.meas;
        std::printf("  %s cell:\n", cell ? "data" : "measure");
        const auto s = split_cell(q, S.cg, S.pitch, 5.5, mo);
        std::printf("    C_open %.3f fF, C_eff(ideal short) %.4f fF (line_coupling: %.4f), C_qubit-pocket %.3f fF, C_sigma %.2f fF, C_xy %.4f fF\n",
                    s.C_open, s.C_eff0, cell ? B.data.C_flux_fF : B.meas.C_flux_fF, s.C_qp, s.C_sigma, s.C_xy);
        // T1 ratio at the worst flux qubit of this cell type (by its f01 and design T1)
        for (double Ls_pH : {0.0, 10.0, 25.0, 50.0, 75.0}) {
            const double dz = Ls_pH * 1e-12 * v / Z0 * 1e6;    // um
            double worst = 1e300; int wqq = -1;
            for (std::size_t i = 0; i < D.qubits.size(); ++i) {
                if ((D.qubits[i].role != "data") != (cell == 0)) continue;
                const double f = R.q[i].f01, wv = 2 * 3.14159265358979 * f * 1e9;
                const double beta = 2 * 3.14159265358979 * 5.5e9 / v * 1e-6;    // as line_coupling (5.5 GHz)
                double ce = s.C_qp * std::sin(beta * dz);
                for (const auto& [z, qf] : s.zq) ce += qf * std::sin(beta * (z + dz));
                const double t1 = s.C_sigma * 1e-15 / (wv * wv * ce * 1e-15 * ce * 1e-15 * Z0) * 1e6;
                if (t1 / R.q[i].T1_us < worst) { worst = t1 / R.q[i].T1_us; wqq = D.qubits[i].id; }
            }
            double ce_line = 0, ce_all = 0; const double beta = 2 * 3.14159265358979 * 5.5e9 / v * 1e-6;
            for (const auto& [z, qf] : s.zq) ce_line += qf * std::sin(beta * (z + dz));
            ce_all = ce_line + s.C_qp * std::sin(beta * dz);
            std::printf("    L_s %4.0f pH -> dz %5.0f um: C_eff line only %.4f fF, with pocket %.4f fF; worst flux T1 limit %.2fx design T1 (Q%d)\n",
                        Ls_pH, dz, ce_line, ce_all, worst, wqq);
        }
    }
    return 0;
}
