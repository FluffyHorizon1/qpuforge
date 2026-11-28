// layout_gates: the mask gates of one chip from a cell file, with the time each step takes. A development and review
// driver for docs/design/mask-connectivity.md: it prints what the review programs print about the mask (DRC, metal
// connectivity and ground bumps, the extracted cells against their targets, the control-line budget) without the rest
// of the review, and the flux coupling of any one line to any one SQUID with a chosen solver. With combined lines the
// gates also report, as a check and not a gate, M_cm by reciprocity with the island's own charge (flux_cm_by_charge, 1 um
// cells) for the two qubits' own lines and the worst drive crosstalk; flux reports it for a qubit and a grid solve.
//
//   layout_gates gates <chip> <cells.json> [--cell UM] [--cross-mesh H] [--grid-margin F] [--window-mesh H] [--no-crosstalk] [--bump-sensitivity] [--gds <out.gds>]
//   layout_gates flux  <chip> <cells.json> <owner> <line> [--cell UM | --moments H] [--grid-margin F] [--zone UM] [--margin UM] [--local UM] [--dump <window.json>]
//
// <chip> is a preset as the reviews use it: alpha1-v11 (AF-17, distance 3), beta1 (BF-49, distance 5) or gamma1
// (GF-97, distance 7). Flux couplings are solved on a grid of square cells (--cell, default 0.5 um, kFluxGrid), on a
// periodic grid that extends the window with film by --grid-margin of its larger side (default 0.1); the crosstalk
// solves by the method of moments (--cross-mesh H, default 2.5 um, kFluxCrossMesh); flux --moments H solves by the
// method of moments instead, on edge panels of H um. A capacitance mesh is given by its edge panel size H in um (largest panel 50 um, growth 1.6
// below 2 um, 2.0 from 2 um): cells and control windows use 1.5 um, as the reviews do. --bump-sensitivity evaluates the budget again (no new solve)
// with the bump inductance halved and doubled, and with all of a qubit's current and half of a neighbour's taken through
// the bumps (quadrant fractions 1 and 0.5 instead of 0.75 and 0.25). flux: --margin sets the solve window's margin
// around the line, the loop, the owner and the owner's whole etch (default 300 um, kFluxMargin_um, as line_budget),
// --zone the branch rule's zone around the owner (default 40 um, kFluxZone_um; a grid solve
// only), --local keeps only the etch within UM of the line and of the owner; --dump writes the solve's window and
// result for tools/flux_check.py.
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
                 "usage: layout_gates gates <chip> <cells.json> [--cell UM] [--cross-mesh H] [--grid-margin F] [--window-mesh H] [--no-crosstalk] [--bump-sensitivity] [--gds <out.gds>]\n"
                 "       layout_gates flux  <chip> <cells.json> <owner> <line> [--cell UM | --moments H] [--grid-margin F] [--zone UM] [--margin UM] [--local UM] [--dump <window.json>]\n"
                 "  <chip>: alpha1-v11 | beta1 | gamma1 (the presets of the reviews); UM: grid cell in um; H: edge panel size in um\n",
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
// The charge of a qubit's island on its panels for 1 V on the island, from a control window (1.5 um panels): for the
// reciprocity check of M_cm (flux_cm_by_charge) and for tools/flux_check.py.
struct IslandCharge { std::vector<double> x, y, q; double owner{0}, ground{0}, other{0}, seconds{0}; bool ok{false}; };
IslandCharge island_charge(const ChipLayout& L, int qubit) {
    IslandCharge c;
    const auto t0 = Clock::now();
    PanelSolution ps;
    const ControlWindow cw = control_window(L, qubit, mesh_of(1.5), 200, &ps);
    c.seconds = since(t0);
    const int qi = cw.C.index("Q" + std::to_string(qubit)), gi = cw.C.index("ground");
    if (qi < 0 || std::size_t(qi) >= ps.q_fF.size()) return c;
    for (std::size_t p = 0; p < ps.panels.size(); ++p) {
        const double q = ps.q_fF[std::size_t(qi)][p];
        if (ps.panels[p].cond == qi) { c.owner += q; c.x.push_back(ps.panels[p].r.cx()); c.y.push_back(ps.panels[p].r.cy()); c.q.push_back(q); }
        else if (ps.panels[p].cond == gi) c.ground += q;
        else c.other += q;
    }
    c.ok = !c.x.empty();
    return c;
}
// M_cm of `owner` (a qubit) from `line` with the even spread and with the island's own charge, on `grid`
void print_reciprocity(const ChipLayout& L, const std::string& owner, const std::string& line, const FluxGrid& grid) {
    if (owner.size() < 2 || owner[0] != 'Q') return;
    FluxWindow fw;
    const FluxCoupling F = flux_coupling(L, owner, line, grid, kFluxMargin_um, 0, 0, &fw);
    const IslandCharge c = island_charge(L, std::stoi(owner.substr(1)));
    if (!F.ok || !c.ok) { std::printf("    reciprocity %s <- %s: not evaluated (%s)\n", owner.c_str(), line.c_str(), F.ok ? "no island charge" : F.note.c_str()); return; }
    const double r = flux_cm_by_charge(fw, c.x, c.y, c.q);
    std::printf("    reciprocity %s <- %s (%.2f um cells; the island's charge from 1.5 um panels, %.1f fF of it to ground, %.1f fF to other conductors): M_cm %+.4f pH with the "
                "even spread, %+.4f pH with the charge (%+.1f%%); %.1f s + %.1f s\n", owner.c_str(), line.c_str(), grid.cell_um, c.ground, c.other, F.M_cm_signed_pH, r,
                F.M_cm_signed_pH != 0 ? 100 * (r / F.M_cm_signed_pH - 1) : 0.0, F.seconds, c.seconds);
}
void print_flux(const char* what, const std::string& owner, const std::string& line, const FluxCoupling& F, const char* solver) {
    std::printf("  %s %s <- %s (%s): %s M %.4f pH, M_cm %.4f pH (signed %.4f), I(Phi0/2) %.4g mA; %d apertures, %d %s, %.1f s%s%s\n", what, owner.c_str(),
                line.c_str(), solver, F.ok ? "" : "NOT EVALUATED", F.M_pH, F.M_cm_pH, F.M_cm_signed_pH, std::min(F.I_half_mA, 1e9), F.apertures, F.panels,
                F.iterations > 0 ? "cells" : "panels", F.seconds, F.note.empty() ? "" : "; ", F.note.c_str());
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
    double cell = kFluxGrid.cell_um, h_cross = kFluxCrossMesh.h_min, grid_margin = kFluxGrid.margin_frac, zone = kFluxZone_um, h_window = 1.5, margin = kFluxMargin_um,
           local = 0, moments = 0;
    bool crosstalk = true, bump_sensitivity = false;
    std::string gds, dump;
    for (std::size_t i = first_opt; i < a.size(); ++i) {
        auto num = [&](double& v, double lo = 0.2, double hi = 20) {     // the next argument, a number in [lo, hi]
            if (i + 1 >= a.size()) return false;
            char* end = nullptr; v = std::strtod(a[i + 1].c_str(), &end);
            ++i;
            return end && end != a[i].c_str() && *end == 0 && v >= lo && v <= hi;
        };
        if (a[i] == "--cell") { if (!num(cell, 0.1, 5)) return usage("--cell needs a size in um (0.1 .. 5)"); }
        else if (a[i] == "--zone" && cmd == "flux") { if (!num(zone, 1, 1e6)) return usage("--zone needs a size in um (1 .. 1e6)"); }
        else if (a[i] == "--grid-margin") { if (!num(grid_margin, 0.05, 1)) return usage("--grid-margin needs a fraction of the window (0.05 .. 1)"); }
        else if (a[i] == "--cross-mesh" && cmd == "gates") { if (!num(h_cross)) return usage("--cross-mesh needs a size in um (0.2 .. 20)"); }
        else if (a[i] == "--moments" && cmd == "flux") { if (!num(moments)) return usage("--moments needs an edge panel size in um (0.2 .. 20)"); }
        else if (a[i] == "--window-mesh" && cmd == "gates") { if (!num(h_window)) return usage("--window-mesh needs a size in um (0.2 .. 20)"); }
        else if (a[i] == "--no-crosstalk" && cmd == "gates") crosstalk = false;
        else if (a[i] == "--margin" && cmd == "flux") { if (!num(margin, 20, 1000)) return usage("--margin needs a size in um (20 .. 1000)"); }
        else if (a[i] == "--local" && cmd == "flux") { if (!num(local, 20, 1000)) return usage("--local needs a size in um (20 .. 1000)"); }
        else if (a[i] == "--bump-sensitivity" && cmd == "gates") bump_sensitivity = true;
        else if (a[i] == "--gds" && cmd == "gates" && i + 1 < a.size()) gds = a[++i];
        else if (a[i] == "--dump" && cmd == "flux" && i + 1 < a.size()) dump = a[++i];
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
        FluxWindow fw;
        char solver[64];
        std::snprintf(solver, sizeof solver, moments > 0 ? "method of moments, edge panels %.2f um" : "grid, %.2f um cells", moments > 0 ? moments : cell);
        FluxGrid fg = kFluxGrid;
        fg.cell_um = cell; fg.margin_frac = grid_margin; fg.zone_um = zone;
        const FluxCoupling F = moments > 0 ? flux_coupling(L, a[3], a[4], mesh_of(moments), margin, 0, local, &fw) : flux_coupling(L, a[3], a[4], fg, margin, 0, local, &fw);
        print_flux("flux coupling", a[3], a[4], F, solver);
        // the reciprocity check of M_cm (a grid solve of a qubit), and the island's charge for the dump
        IslandCharge charge;
        if (a[3].size() > 1 && a[3][0] == 'Q' && F.ok && (moments == 0 || !dump.empty())) charge = island_charge(L, std::stoi(a[3].substr(1)));
        if (charge.ok && moments == 0) {
            const double r = flux_cm_by_charge(fw, charge.x, charge.y, charge.q);
            std::printf("  reciprocity: M_cm %+.4f pH with the island's own charge (1.5 um panels; %.1f fF to ground, %.1f fF to other conductors) against %+.4f pH with the "
                        "even spread (%+.1f%%); %.1f s\n", r, charge.ground, charge.other, F.M_cm_signed_pH, F.M_cm_signed_pH != 0 ? 100 * (r / F.M_cm_signed_pH - 1) : 0.0,
                        charge.seconds);
        }
        if (!dump.empty()) {
            auto rect = [](const Rect& r) { Json j = Json::array(); for (double v : {r.x0, r.y0, r.x1, r.y1}) j.push(Json(v)); return j; };
            auto rects = [&](const std::vector<Rect>& v) { Json j = Json::array(); for (const auto& r : v) j.push(rect(r)); return j; };
            auto ints = [](const std::vector<int>& v) { Json j = Json::array(); for (int x : v) j.push(Json(double(x))); return j; };
            Json j = Json::object();
            j.set("chip", D.name); j.set("cells", a[2]); j.set("owner", a[3]); j.set("line", a[4]);
            j.set("window", rect(fw.window)); j.set("loop", rect(fw.loop)); j.set("pieces", rects(fw.pieces)); j.set("aperture", ints(fw.aperture));
            j.set("kind", ints(fw.kind)); j.set("gap_a", fw.gap_a); j.set("gap_b", fw.gap_b); j.set("own_moat", fw.own_moat);
            Json centre = Json::array(); centre.push(Json(fw.qx)); centre.push(Json(fw.qy));
            j.set("island_centre", centre);
            j.set("theta0", fw.theta0); j.set("etch", rects(fw.etch)); j.set("island", rects(fw.island));
            if (charge.ok) {                                     // the owner's charge distribution, for the reciprocity check of M_cm
                Json ch = Json::array();
                for (std::size_t p = 0; p < charge.x.size(); ++p) {
                    Json e = Json::array();
                    for (double v : {charge.x[p], charge.y[p], charge.q[p]}) e.push(Json(v));
                    ch.push(e);
                }
                Json c = Json::object();
                c.set("panels", ch); c.set("owner_fF", charge.owner); c.set("ground_fF", charge.ground); c.set("other_conductors_fF", charge.other);
                j.set("island_charge", c);
            }
            Json r = Json::object();
            r.set("mesh_um", moments > 0 ? moments : cell); r.set("panels", F.panels); r.set("method", moments > 0 ? "method of moments" : "grid"); r.set("M_pH", F.M_pH); r.set("M_cm_pH", F.M_cm_pH); r.set("M_cm_signed_pH", F.M_cm_signed_pH);
            j.set("method_of_moments", r);
            std::ofstream o(dump);
            j.dump(o);
            o << "\n";
            if (!o) { std::fprintf(stderr, "layout_gates: cannot write %s\n", dump.c_str()); return 2; }
        }
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
    FluxGrid fg = kFluxGrid, fc{1.0};                       // fc: the reciprocity check
    fg.cell_um = cell; fg.margin_frac = fc.margin_frac = grid_margin;
    LineBudget B = line_budget(D, R, L, m, P.park_fc_GHz, mesh_of(h_window), 50, fg);
    const double t_budget = since(t0);
    double t_cross = 0;
    if (crosstalk) { t0 = Clock::now(); line_budget_crosstalk(B, D, R, L, m, mesh_of(h_cross)); t_cross = since(t0); }
    char solver[64];
    std::snprintf(solver, sizeof solver, "grid, %.2f um cells", cell);
    print_flux("measure qubit", "Q" + std::to_string(B.meas.qubit), "FQ" + std::to_string(B.meas.qubit), B.meas.flux, solver);
    print_flux("data qubit   ", "Q" + std::to_string(B.data.qubit), "FQ" + std::to_string(B.data.qubit), B.data.flux, solver);
    print_flux("coupler      ", "C" + B.coupler, "FC" + B.coupler, B.coupler_flux, solver);
    std::printf("  control windows (mesh %.2f um): Q%d %d panels %.1f s, Q%d %d panels %.1f s\n", h_window, B.meas.qubit, B.meas.win.C.panels, B.meas.win.C.seconds,
                B.data.qubit, B.data.win.C.panels, B.data.win.C.seconds);
    std::printf("  %s\n", line_budget_summary(B).c_str());
    // a check of M_cm's even spread, not part of the gates: M_cm by reciprocity with the island's own charge, on the
    // crosstalk grid, for the two qubits' own lines and for the worst drive crosstalk
    double t_recip = 0;
    if (!B.separate_xy && B.solved) {
        t0 = Clock::now();
        std::printf("  reciprocity check of M_cm (not a gate):\n");
        print_reciprocity(L, "Q" + std::to_string(B.meas.qubit), "FQ" + std::to_string(B.meas.qubit), fc);
        print_reciprocity(L, "Q" + std::to_string(B.data.qubit), "FQ" + std::to_string(B.data.qubit), fc);
        const std::size_t arrow = B.worst_crosstalk.find(" -> ");
        if (arrow != std::string::npos) print_reciprocity(L, B.worst_crosstalk.substr(arrow + 4), B.worst_crosstalk.substr(0, arrow), fc);
        t_recip = since(t0);
    }

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
    std::printf("  time (%d threads): assemble and check the mask %.1f s, cells and targets %.1f s, line budget %.1f s (windows %.1f s, flux solves at the gate's cell %.1f s), "
                "crosstalk %.1f s, reciprocity check %.1f s, total %.1f s\n", threads, t_assemble, t_cells, t_budget, B.meas.win.C.seconds + B.data.win.C.seconds,
                B.meas.flux.seconds + B.data.flux.seconds + B.coupler_flux.seconds, t_cross, t_recip, since(t_start));
    return pass ? 0 : 1;
}
