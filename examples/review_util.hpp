// Shared helpers of the Beta-1 and Gamma-1 design reviews. Command line and provenance: review_io.hpp.
#pragma once
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include "review_io.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace review {
using namespace qpuforge;

inline Options g_opt;                           // command line of the running review (set in main)
using Tweak = std::function<void(NoiseModel&)>;

inline ChipReport prepared_report(const ChipDesign& D) {
    auto R = analyze(D);
    if (D.coupler.physical) calibrate_chip_cz(D, R, false, g_opt.cache);
    return R;
}

struct QecPoint {
    double eps{0}, eps_err{0}, det{0}, round_ns{0}, leak_d{0}, leak_m{0}, p_last{0};
    double mean_defects{0}, frac_over_22{0}, decode_us_per_round{0};
    int max_defects{0}, rounds_last{0};
    long failures_last{0}, shots{0};
    std::vector<double> detection_by_round;
};

// eps from a fit over the round counts (several) or from the one run (single round count)
inline QecPoint qec(const ChipDesign& D, const std::vector<int>& rounds, long shots, unsigned seed, char basis = 'Z', const Tweak& tw = {},
                    DecodeMethod method = DecodeMethod::Auto) {
    const auto R = prepared_report(D);
    auto nm = noise_model(D, R);
    if (tw) tw(nm);
    FrameRunOptions opt; opt.method = method;
    std::vector<MemoryResult> rs;
    for (int r : rounds) rs.push_back(run_memory_frame(prepare_memory(D, R, nm, r, basis), shots, seed + unsigned(r) + (basis == 'X' ? 7919u : 0u), opt));
    QecPoint q;
    const auto& last = rs.back();
    q.eps = rs.size() > 1 ? fit_eps_per_round(rs) : last.eps_per_round;
    q.eps_err = fit_eps_per_round_err(rs);
    q.det = last.detection_rate; q.round_ns = last.round_ns; q.leak_d = last.leak_frac_data; q.leak_m = last.leak_frac_meas; q.p_last = last.p_logical;
    q.mean_defects = last.mean_defects; q.max_defects = last.max_defects; q.frac_over_22 = last.frac_over_22;
    q.decode_us_per_round = last.decode_us_per_shot / last.rounds;
    q.rounds_last = last.rounds; q.failures_last = last.logical_errors; q.shots = last.shots;
    q.detection_by_round = last.detection_by_round;
    return q;
}

// relative change of the detection fraction between the first and the last `window` syndrome rounds
inline double detection_drift(const std::vector<double>& by_round, int window) { return chz::detection_drift(by_round, window); }

struct WiringCount { int xy{0}, flux_qubit{0}, flux_coupler{0}, feedlines{0}, total{0}; double per_feedline{0}; };
inline WiringCount wiring(const ChipDesign& D) {
    std::set<int> feeds; for (const auto& q : D.qubits) feeds.insert(q.feedline);
    WiringCount w;
    w.xy = w.flux_qubit = int(D.qubits.size()); w.flux_coupler = int(D.couplings.size()); w.feedlines = int(feeds.size());
    w.total = w.xy + w.flux_qubit + w.flux_coupler + 2 * w.feedlines;
    w.per_feedline = double(D.qubits.size()) / std::max(w.feedlines, 1);
    return w;
}
inline Json wiring_json(const WiringCount& w) {
    Json j = Json::object();
    j.set("xy", w.xy); j.set("flux_qubit", w.flux_qubit); j.set("flux_coupler", w.flux_coupler); j.set("readout_in", w.feedlines); j.set("readout_out", w.feedlines);
    j.set("total_lines", w.total); j.set("feedlines", w.feedlines); j.set("qubits_per_feedline", w.per_feedline);
    return j;
}

// Chip assembly from a synthesised cell set, GDSII, the placement checks and the mask-level width/space DRC, and the
// extracted circuit against the design's targets. `gds_path` is the committed location (redirected by --out).
struct LayoutCheck {
    bool found{false}, placement_clean{false}, mask_drc_clean{false}, connectivity_clean{false}, lines_ok{false}, within_tolerance{false}, renull_ok{false};
    double max_ec_err{0}, max_g_err{0}, readout_c_err{0};
    Json json;
};
inline LayoutCheck layout_section(const ChipDesign& D, const std::string& cells_path, const std::string& gds_path, const std::string& cell_name) {
    using namespace layout;
    LayoutCheck out; out.json = Json::object();
    std::ifstream f(cells_path);
    if (!f) { std::printf("   (cell geometry %s not found — synthesise it first: qpuforge gds <chip.json> <out.gds> --cells %s --synthesize)\n", cells_path.c_str(), cells_path.c_str()); return out; }
    out.found = true;
    std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto R = analyze(D);
    const auto T = targets_from_design(D, R);
    const auto m = evaluate_cells(C, g_opt.smoke ? MeshOptions{3.0, 50.0, 2.0} : MeshOptions{1.5, 50.0, 1.6});   // smoke: coarse mesh
    const auto L = assemble_chip(D, C);
    const std::string gds_out = g_opt.path(gds_path);
    write_gds(L, gds_out, cell_name);
    const auto P = post_layout(D, R, m);
    auto rel = [](double a, double t) { return a / t - 1; };
    out.max_ec_err = std::max({std::abs(rel(m.EC_m, T.EC_m)), std::abs(rel(m.EC_d, T.EC_d)), std::abs(rel(m.EC_c, T.EC_c))});
    out.max_g_err = std::max({std::abs(rel(m.G_mc, T.G_mc)), std::abs(rel(m.G_cd, T.G_cd)), std::abs(rel(m.G_md, T.G_md))});
    out.readout_c_err = rel(m.C_qr_fF, T.C_qr_fF);
    out.placement_clean = L.problems.empty();
    out.mask_drc_clean = L.drc_clean();
    out.connectivity_clean = L.conn.clean();
    const auto LB = line_budget(D, R, C, g_opt.smoke ? MeshOptions{3.0, 50.0, 2.0} : MeshOptions{1.5, 50.0, 1.6});
    out.lines_ok = LB.ok();
    out.within_tolerance = out.max_ec_err <= 0.03 && out.max_g_err <= 0.06 && std::abs(out.readout_c_err) <= 0.06;
    out.renull_ok = std::abs(P.park_gap_MHz) < 0.05;
    std::printf("   chip %.1f x %.1f mm, %zu etch polygons, %zu junction regions, %zu bump pads -> %s\n   DRC: %s\n",
                L.outline.w() / 1e3, L.outline.h() / 1e3, L.etch.size(), L.jj.size(), L.bumps.size(), gds_path.c_str(), drc_summary(L).c_str());
    for (const auto& p : L.problems) std::printf("     %s\n", p.c_str());
    std::printf("   %s\n", connectivity_summary(L).c_str());
    std::printf("   %s\n", line_budget_summary(LB).c_str());
    std::printf("   extracted vs target: E_C within %.1f%%, couplings within %.1f%%, qubit-readout capacitance %+.1f%% (%.2f fF, target %.2f)\n",
                100 * out.max_ec_err, 100 * out.max_g_err, 100 * out.readout_c_err, m.C_qr_fF, T.C_qr_fF);
    std::printf("   post-layout: g_qc %.1f MHz, readout g %.1f MHz (design %.1f), parked exchange re-nulled at %.3f GHz (residual gap %.3f MHz, ZZ %+.2f kHz)\n",
                P.g_qc_MHz, P.g_r_MHz, P.g_r_target_MHz, P.park_fc_GHz, P.park_gap_MHz, P.park_zz_kHz);
    Json& j = out.json;
    j.set("width_mm", L.outline.w() / 1e3); j.set("height_mm", L.outline.h() / 1e3);
    j.set("drc_clean", out.mask_drc_clean); j.set("placement_clean", out.placement_clean); j.set("drc", drc_json(L)); j.set("connectivity", connectivity_json(L)); j.set("control_lines", line_budget_json(LB));
    j.set("etch_polygons", int(L.etch.size())); j.set("junction_regions", int(L.jj.size())); j.set("bumps", int(L.bumps.size()));
    j.set("max_EC_err", out.max_ec_err); j.set("max_coupling_err", out.max_g_err); j.set("readout_C_err", out.readout_c_err);
    j.set("readout_C_fF", m.C_qr_fF); j.set("readout_C_target_fF", T.C_qr_fF); j.set("readout_g_MHz", P.g_r_MHz); j.set("readout_g_target_MHz", P.g_r_target_MHz);
    j.set("g_qc_MHz", P.g_qc_MHz); j.set("park_fc_GHz", P.park_fc_GHz); j.set("park_gap_MHz", P.park_gap_MHz); j.set("park_zz_kHz", P.park_zz_kHz);
    j.set("within_tolerance", out.within_tolerance); j.set("renull_ok", out.renull_ok); j.set("cells", cells_path);
    return out;
}

} // namespace review
