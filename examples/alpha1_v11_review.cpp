// Alpha-1 v1.1 design review — the Alpha-1 chip re-evaluated with the new physics layers
//   1. pulse-level CZ calibration of all 24 pairs (physical tunable couplers)
//   2. headline: logical error per round, v1.0 model vs v1.1 model
//   3. leakage: no removal / leakage-removing reset / + data-qubit leakage removal, short and long experiments
//   4. leave-one-out error budget
//   5. sensitivity (what the next revision should buy)
//   6. distance-5 projection, suppression factor Lambda
//   7. physical layout: capacitance synthesis vs targets, post-layout circuit check
// Writes designs/alpha1-v11/alpha1_v11_review.json and prints a report.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include "review_io.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>

using namespace qpuforge;

static review::Options g_opt;   // command line (review_io.hpp)

struct QecPoint { double eps{0}, eps_err{0}, det{0}, round_ns{0}, leak_d{0}, leak_m{0}, p_last{0}; int r_last{0}; };
using Tweak = std::function<void(NoiseModel&)>;

static ChipReport prepared_report(const ChipDesign& D) {
    auto R = analyze(D);
    if (D.coupler.physical) calibrate_chip_cz(D, R, false, g_opt.cache);
    return R;
}
static QecPoint qec(const ChipDesign& D, const std::vector<int>& rounds, long shots, unsigned seed, const Tweak& tw = {}, char basis = 'Z') {
    const auto R = prepared_report(D);
    auto nm = noise_model(D, R);
    if (tw) tw(nm);
    std::vector<MemoryResult> rs;
    for (int r : rounds) rs.push_back(run_memory_frame(prepare_memory(D, R, nm, r, basis), shots, seed + unsigned(r)));
    QecPoint q;
    q.eps = fit_eps_per_round(rs);
    const auto& last = rs.back();
    q.eps_err = fit_eps_per_round_err(rs);
    q.det = last.detection_rate; q.round_ns = last.round_ns; q.leak_d = last.leak_frac_data; q.leak_m = last.leak_frac_meas;
    q.p_last = last.p_logical; q.r_last = last.rounds;
    return q;
}
static Json jq(const QecPoint& q) {
    Json o = Json::object();
    o.set("eps_per_round", q.eps); o.set("eps_err", q.eps_err); o.set("detection_rate", q.det); o.set("round_ns", q.round_ns);
    o.set("leaked_data", q.leak_d); o.set("leaked_measure", q.leak_m); o.set("P_L_last", q.p_last); o.set("rounds_last", q.r_last);
    return o;
}

int main(int argc, char** argv) {
    g_opt = review::parse(argc, argv, 100000);
    const long shots = g_opt.shots;
    const auto t0 = std::chrono::steady_clock::now();
    Json out = Json::object();
    const std::vector<int> R5{5, 10, 15, 20, 25};

    const ChipDesign V10 = alpha1(3);
    const ChipDesign V11 = alpha1_v11(3);
    save_design(V11, g_opt.path("designs/alpha1-v11/alpha1_v11.json"));

    // ------------------------------------------------------------------ 1. CZ calibration
    std::puts("1. Pulse-level CZ calibration (physical tunable couplers, 3-knob calibration per pair)");
    auto R11 = analyze(V11);
    const auto cal = calibrate_chip_cz(V11, R11, true, g_opt.cache);
    {
        double tmin = 1e9, tmax = 0, tm = 0, coh = 0, lk = 0, lkmax = 0, ctrl = 0, tot = 0, zz = 0, fl = 0;
        for (const auto& c : cal) {
            tmin = std::min(tmin, c.t_cz_ns); tmax = std::max(tmax, c.t_cz_ns); tm += c.t_cz_ns; coh += c.err_coherent; lk += c.leakage;
            lkmax = std::max(lkmax, c.leakage); ctrl += c.err_ctrl; tot += c.err_total; zz += std::abs(c.zz_idle_kHz); fl += c.err_flux_local + c.err_flux_quasistatic;
        }
        const double n = double(cal.size());
        std::printf("   %zu pairs: t_cz %.1f-%.1f ns (mean %.1f; v1.0 abstract model 23.3 ns)\n", cal.size(), tmin, tmax, tm / n);
        std::printf("   coherent pulse error %.1e, leakage %.1e per CZ (max %.1e), flux-noise error %.1e, residual idle ZZ %.1f kHz\n", coh / n, lk / n, lkmax, fl / n, zz / n);
        std::printf("   control error handed to the emulator %.2e (incl. %.0e residual allowance); total CZ error incl. decoherence %.2e\n",
                    ctrl / n, V11.coupler.residual_ctrl_err, tot / n);
        Json c = Json::object();
        c.set("pairs", int(cal.size())); c.set("t_cz_min_ns", tmin); c.set("t_cz_max_ns", tmax); c.set("t_cz_mean_ns", tm / n);
        c.set("err_coherent_mean", coh / n); c.set("leakage_mean", lk / n); c.set("leakage_max", lkmax); c.set("err_flux_mean", fl / n);
        c.set("zz_idle_abs_mean_kHz", zz / n); c.set("err_ctrl_mean", ctrl / n); c.set("err_total_mean", tot / n);
        c.set("per_pair", cz_calibration_json(cal));
        out.set("cz_calibration", c);
    }

    // ------------------------------------------------------------------ 2. headline
    std::puts("\n2. Headline: distance-3 memory, logical error per round (fit over 5-25 rounds)");
    const auto h10 = qec(V10, R5, shots, 101);
    const auto h11 = qec(V11, R5, 2 * shots, 101);
    double best_T1 = 0; for (const auto& q : R11.q) best_T1 = std::max(best_T1, q.T1_us);
    std::printf("   v1.0 model (abstract coupler, no leakage):     eps %.4f%% +- %.4f%%, round %.0f ns\n", 100 * h10.eps, 100 * h10.eps_err, h10.round_ns);
    std::printf("   v1.1 model (pulse-level CZ, leakage, LRR):     eps %.4f%% +- %.4f%%, round %.0f ns, detection %.4f\n", 100 * h11.eps, 100 * h11.eps_err, h11.round_ns, h11.det);
    std::printf("   logical lifetime %.0f us vs best physical T1 %.1f us\n", h11.round_ns / 1e3 / h11.eps, best_T1);
    const auto h11x = qec(V11, R5, shots, 111, {}, 'X');
    std::printf("   v1.1 model, memory-X (|+>_L, phase flips):     eps %.4f%% +- %.4f%%\n", 100 * h11x.eps, 100 * h11x.eps_err);
    { Json h = Json::object(); h.set("v1_0", jq(h10)); h.set("v1_1", jq(h11)); h.set("v1_1_memory_X", jq(h11x)); h.set("logical_lifetime_us", h11.round_ns / 1e3 / h11.eps); h.set("best_T1_us", best_T1); out.set("headline", h); }

    // ------------------------------------------------------------------ 3. leakage
    std::puts("\n3. Leakage handling (eps fit over 5-25 rounds; 50-round experiment for the long-run behaviour)");
    Json lj = Json::object();
    struct LV { const char* key; const char* name; bool model, lrr, dqlr; };
    for (const LV& v : {LV{"off", "leakage not modelled", false, false, false}, LV{"none", "no leakage removal", true, false, false},
                        LV{"lrr", "leakage-removing reset (v1.1)", true, true, false}, LV{"lrr_dqlr", "LRR + data-qubit leakage removal", true, true, true}}) {
        ChipDesign D = V11; D.leakage.model = v.model; D.leakage.reset_removes_leakage = v.lrr; D.leakage.dqlr = v.dqlr;
        const auto q = qec(D, R5, shots, 202);
        const auto q50 = qec(D, {50}, shots / 2, 303);
        std::printf("   %-34s eps %.4f%% | 50 rounds: P_L %.4f, leaked data %.1e meas %.1e, round %.0f ns\n", v.name, 100 * q.eps, q50.p_last, q50.leak_d, q50.leak_m, q50.round_ns);
        Json o = jq(q); o.set("P_L_50_rounds", q50.p_last); o.set("leaked_data_50", q50.leak_d); o.set("leaked_meas_50", q50.leak_m);
        lj.set(v.key, o);
    }
    { // stress: 10x leakage per CZ
        for (bool dq : {false, true}) {
            ChipDesign D = V11; D.leakage.dqlr = dq;
            const auto q50 = qec(D, {50}, shots / 2, 404, [](NoiseModel& nm) { for (auto& kv : nm.leak_cz) kv.second.second *= 10; });
            std::printf("   stress: 10x CZ leakage, %-20s 50 rounds: P_L %.4f, leaked data %.1e\n", dq ? "LRR + DQLR" : "LRR only", q50.p_last, q50.leak_d);
            lj.set(dq ? "stress10_lrr_dqlr_P_L_50" : "stress10_lrr_P_L_50", q50.p_last);
        }
    }
    out.set("leakage", lj);

    // ------------------------------------------------------------------ 4. leave-one-out budget
    std::puts("\n4. Error budget (leave one source out; change in eps)");
    Json bj = Json::object();
    for (const char* src : {"T1", "Tphi", "control", "readout", "reset", "zz", "leakage", "thermal"}) {
        const std::string s = src;
        const auto q = qec(V11, R5, shots, 505, [s](NoiseModel& nm) { nm.disable(s); });
        const double d = q.eps / h11.eps - 1;
        std::printf("   without %-8s eps %.4f%%  (%+.0f%%)\n", src, 100 * q.eps, 100 * d);
        bj.set(src, d);
    }
    out.set("budget_leave_one_out", bj);

    // ------------------------------------------------------------------ 5. sensitivity
    std::puts("\n5. Sensitivity (single changes to v1.1)");
    Json sj = Json::object();
    struct SV { const char* key; const char* name; std::function<void(ChipDesign&)> f; };
    const std::vector<SV> sens = {
        {"resid_0p5e-3", "CZ residual control error 1e-3 -> 5e-4", [](ChipDesign& D) { D.coupler.residual_ctrl_err = 5e-4; }},
        {"resid_2e-3", "CZ residual control error 1e-3 -> 2e-3", [](ChipDesign& D) { D.coupler.residual_ctrl_err = 2e-3; }},
        {"Q_x2", "dielectric Q x2 (T1 x~2)", [](ChipDesign& D) { D.Q_dielectric *= 2; }},
        {"readout_300", "readout 500 -> 300 ns", [](ChipDesign& D) { D.t_readout_ns = 300; }},
        {"Teff_35", "effective temperature 45 -> 35 mK", [](ChipDesign& D) { D.T_eff_mK = 35; }},
        {"flux_x3", "flux noise x3 (6 uPhi0)", [](ChipDesign& D) { D.flux_noise_uPhi0 *= 3; }},
        {"sep_0p25", "readout separation error 0.5% -> 0.25%", [](ChipDesign& D) { D.readout_sep_err = 0.0025; }},
    };
    for (const auto& v : sens) {
        ChipDesign D = V11; v.f(D);
        const auto q = qec(D, R5, shots, 606);
        std::printf("   %-42s eps %.4f%%  (%+.0f%%)\n", v.name, 100 * q.eps, 100 * (q.eps / h11.eps - 1));
        sj.set(v.key, q.eps / h11.eps - 1);
    }
    out.set("sensitivity", sj);

    // ------------------------------------------------------------------ 6. distance 5
    std::puts("\n6. Distance-5 projection with the v1.1 physics (49 qubits)");
    {
        const ChipDesign D5 = alpha1_v11(5);
        const auto q5 = qec(D5, {5, 10, 15}, shots / 2, 707);
        const auto q3 = qec(V11, {5, 10, 15}, shots, 707);
        const double lam = q3.eps / q5.eps, lam_err = lam * std::sqrt(std::pow(q3.eps_err / q3.eps, 2) + std::pow(q5.eps_err / q5.eps, 2));
        std::printf("   fit over 5-15 rounds: eps d=3 %.4f%% +- %.4f%%, d=5 %.4f%% +- %.4f%%  ->  Lambda = %.2f +- %.2f\n", 100 * q3.eps, 100 * q3.eps_err, 100 * q5.eps,
                    100 * q5.eps_err, lam, lam_err);
        Json d = Json::object(); d.set("eps_d3", q3.eps); d.set("eps_d5", q5.eps); d.set("Lambda", lam); d.set("round_ns_d5", q5.round_ns);
        d.set("eps_err_d3", q3.eps_err); d.set("eps_err_d5", q5.eps_err); d.set("Lambda_err", lam_err); d.set("fit_rounds", "5, 10, 15");
        out.set("distance5", d);
    }

    // ------------------------------------------------------------------ 7. layout
    std::puts("\n7. Physical layout (designs/alpha1-v11/layout_cells.json, AF17_v11.gds)");
    {
        std::ifstream f("designs/alpha1-v11/layout_cells.json");
        if (f) {
            std::stringstream ss; ss << f.rdbuf();
            const Json j = Json::parse(ss.str());
            const auto C = layout::cells_from_json(j);
            const auto m = layout::evaluate_cells(C, layout::MeshOptions{1.5, 50.0, 1.6});
            const auto P = layout::post_layout(V11, R11, m);
            const auto L = layout::assemble_chip(V11, C);
            std::printf("   cells: measure arm %.1f um, data arm %.1f um, coupler pad %.0f um, coupler T-bar %.1f um, readout pad %.1f um\n",
                        C.meas.L, C.data.L, C.cg.pad, C.cg.l_c, C.meas.r_pad);
            std::printf("   extracted: C_sigma %.1f/%.1f fF, C_coupler %.0f fF, C_qc %.2f fF, C_qq %.3f fF, C_qr %.2f fF\n", m.C_sigma_m, m.C_sigma_d, m.C_c, m.C_qc,
                        m.C_qq, m.C_qr_fF);
            std::printf("   post-layout: alpha %.0f/%.0f MHz, coupler alpha %.0f MHz, g_qc %.1f MHz, g12 %.2f MHz (needs %.2f), residual exchange gap %.3f MHz,\n"
                        "                spectator transfer %.1e per crossing, ZZ at park %+.2f kHz, readout g %.1f MHz (design %.1f)\n",
                        P.alpha_m * 1e3, P.alpha_d * 1e3, P.alpha_c * 1e3, P.g_qc_MHz, P.g12_MHz, P.g12_null_MHz, P.residual_gap_MHz, P.spectator_LZ, P.zz_park_kHz,
                        P.g_r_MHz, P.g_r_target_MHz);
            std::printf("   re-null: couplers parked at %.3f GHz (%.0f MHz below sweet spot, df/dPhi %.2f GHz/Phi0) -> residual gap %.3f MHz, ZZ %+.2f kHz\n",
                        P.park_fc_GHz, (V11.coupler.f_max_GHz - P.park_fc_GHz) * 1e3, P.park_dfdphi_GHz, P.park_gap_MHz, P.park_zz_kHz);
            std::printf("   chip %.1f x %.1f mm, %zu etch polygons\n   DRC: %s\n", L.outline.w() / 1e3, L.outline.h() / 1e3, L.etch.size(), layout::drc_summary(L).c_str());
            const auto T = layout::targets_from_design(V11, R11);
            auto rel = [](double a, double t) { return a / t - 1; };
            const double ec_err = std::max({std::abs(rel(m.EC_m, T.EC_m)), std::abs(rel(m.EC_d, T.EC_d)), std::abs(rel(m.EC_c, T.EC_c))});
            const double g_err = std::max({std::abs(rel(m.G_mc, T.G_mc)), std::abs(rel(m.G_cd, T.G_cd)), std::abs(rel(m.G_md, T.G_md))});
            const double cr_err = rel(m.C_qr_fF, T.C_qr_fF);
            std::printf("   extracted vs target: E_C within %.1f%%, couplings within %.1f%%, qubit-readout capacitance %+.1f%% (%.2f fF, target %.2f)\n",
                        100 * ec_err, 100 * g_err, 100 * cr_err, m.C_qr_fF, T.C_qr_fF);
            Json lj2 = Json::object();
            lj2.set("alpha_measure_MHz", P.alpha_m * 1e3); lj2.set("alpha_data_MHz", P.alpha_d * 1e3); lj2.set("alpha_coupler_MHz", P.alpha_c * 1e3);
            lj2.set("g_qc_MHz", P.g_qc_MHz); lj2.set("g12_MHz", P.g12_MHz); lj2.set("g12_null_MHz", P.g12_null_MHz);
            lj2.set("residual_gap_MHz", P.residual_gap_MHz); lj2.set("spectator_transfer", P.spectator_LZ); lj2.set("zz_park_kHz", P.zz_park_kHz);
            lj2.set("renull_park_fc_GHz", P.park_fc_GHz); lj2.set("renull_gap_MHz", P.park_gap_MHz); lj2.set("renull_zz_kHz", P.park_zz_kHz);
            lj2.set("readout_g_MHz", P.g_r_MHz); lj2.set("readout_g_design_MHz", P.g_r_target_MHz); lj2.set("drc_clean", L.drc_clean()); lj2.set("drc", layout::drc_json(L));
            lj2.set("max_EC_err", ec_err); lj2.set("max_coupling_err", g_err); lj2.set("readout_C_err", cr_err);
            lj2.set("within_tolerance", ec_err <= 0.03 && g_err <= 0.06 && std::abs(cr_err) <= 0.06); lj2.set("renull_ok", std::abs(P.park_gap_MHz) < 0.05);
            out.set("layout", lj2);
        } else std::puts("   (no layout_cells.json yet: run `qpuforge gds ... --cells designs/alpha1-v11/layout_cells.json`)");
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.set("shots", double(shots)); out.set("runtime_s", secs);
    out.set("provenance", review::provenance(g_opt, "frame Monte Carlo seeds 101 + rounds (v1.0 vs v1.1), 111 + rounds (memory-X), 202/303/404 + rounds (leakage schemes), "
                                                    "505 + rounds (budget), 606 + rounds (sensitivity), 707 + rounds (distance 5)"));
    const std::string jpath = g_opt.path("designs/alpha1-v11/alpha1_v11_review.json");
    std::ofstream f(jpath); out.dump(f); f << "\n";
    std::printf("\nreview complete in %.0f s -> %s\n", secs, jpath.c_str());
    return 0;
}
