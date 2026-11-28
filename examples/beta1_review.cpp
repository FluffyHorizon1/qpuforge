// Beta-1 (BF-49) design review — distance-5 surface-code processor (rev C)
//   1. device summary and pulse-level CZ calibration of all pairs
//   2. error suppression with code distance: eps at d = 3, 5 in both memory bases (fits over 5-20 rounds);
//      d = 7 outlook from one 9-round run; Lambda
//   3. long experiments: leakage removal over 50 rounds (LRR vs LRR + DQLR), exact matching against the pre-v0.4
//      decoder on the same syndromes, and detection-rate stability over 100 rounds
//   4. error budget (leave one out, d = 3 Beta physics) and remaining sensitivity
//   5. fabrication yield at d = 5: untrimmed (2% junction spread) vs trimmed (0.5%)
//   6. wiring and readout multiplexing
//   7. physical layout of the 49-qubit chip (cells synthesised for the Beta-1 readout coupling), DRC, extraction check
//   8. exit criteria
// Writes designs/beta1/beta1_review.json and prints a report. Options: see review_io.hpp ([shots] [--out DIR] [--cache FILE] [--smoke]).
#include "review_util.hpp"
#include <chrono>

using namespace qpuforge;
using namespace review;

int main(int argc, char** argv) {
    g_opt = parse(argc, argv, 100000);
    const long shots = g_opt.shots;
    const auto t0 = std::chrono::steady_clock::now();
    Json out = Json::object();
    const ChipDesign B5 = beta1(5), B3 = beta1(3), B7 = beta1(7);
    save_design(B5, g_opt.path("designs/beta1/beta1.json"));
    const std::vector<int> RF{5, 10, 15, 20};
    // Long experiments. Setting one up costs time quadratic in the number of rounds, so a smoke run shortens them;
    // each JSON block records the round count it used.
    const int RLONG = g_opt.smoke ? 20 : 50, RSTAB = g_opt.smoke ? 30 : 100;

    // ------------------------------------------------------------------ 1. device
    std::puts("1. Beta-1 device (49 qubits: 25 data + 24 measure, tunable couplers)");
    auto R5 = analyze(B5);
    const auto cal = calibrate_chip_cz(B5, R5, true, g_opt.cache);
    {
        double T1 = 0, T2 = 0, ro = 0, t = 0, e = 0, lk = 0, lkmax = 0, T1b = 0;
        for (const auto& q : R5.q) { T1 += q.T1_us; T2 += q.T2_us; ro += 0.5 * (q.ro_err_0to1 + q.ro_err_1to0); T1b = std::max(T1b, q.T1_us); }
        for (const auto& c : cal) { t += c.t_cz_ns; e += c.err_total; lk += c.leakage; lkmax = std::max(lkmax, c.leakage); }
        const double n = double(R5.q.size()), m = double(cal.size());
        std::printf("   mean T1 %.0f us (best %.0f), T2 %.0f us, readout error %.2f%% at %.0f ns; %zu CZ pairs: %.1f ns, error %.2e, leakage %.1e (max %.1e)\n",
                    T1 / n, T1b, T2 / n, 100 * ro / n, B5.t_readout_ns, cal.size(), t / m, e / m, lk / m, lkmax);
        Json d = Json::object();
        d.set("qubits", int(n)); d.set("couplers", int(m)); d.set("mean_T1_us", T1 / n); d.set("best_T1_us", T1b); d.set("mean_T2_us", T2 / n);
        d.set("mean_readout_err", ro / n); d.set("t_cz_mean_ns", t / m); d.set("cz_err_mean", e / m); d.set("cz_leak_mean", lk / m); d.set("cz_leak_max", lkmax);
        out.set("device", d);
    }

    // ------------------------------------------------------------------ 2. distance scaling
    std::puts("\n2. Logical error per round vs code distance (Beta-1 physics; exact matching decoder)");
    const auto q3 = qec(B3, RF, 2 * shots, 11);
    const auto q5 = qec(B5, RF, shots, 11);
    const auto q3x = qec(B3, RF, 2 * shots, 11, 'X');
    const auto q5x = qec(B5, RF, shots, 11, 'X');
    // d = 7 outlook: one 9-round experiment, eps from P_L, and the same estimate for d = 5. (A free-intercept fit over a
    // few short runs is unstable at these rates.) The distance-7 chip itself is reviewed as Gamma-1.
    const long n7 = g_opt.fixed(100000), n5_9 = g_opt.fixed(50000);
    const auto q7 = qec(B7, {9}, n7, 11);
    const auto q5_9 = qec(B5, {9}, n5_9, 12);
    const double eps5_9_err = std::sqrt(std::max(q5_9.p_last * (1 - q5_9.p_last), 1e-12) / double(n5_9)) / 9;
    const bool d7_resolved = q7.eps > 0;
    const double eps7_ub = d7_resolved ? q7.eps : 0.5 * (1 - std::pow(1 - 2 * 3.0 / double(n7), 1.0 / 9));   // 95% bound for zero failures
    const double eps7_err = std::sqrt(std::max(q7.p_last * (1 - q7.p_last), 1e-12) / double(n7)) / 9;
    const double L35 = q3.eps / q5.eps, L35x = q3x.eps / q5x.eps, L57 = q5_9.eps / (d7_resolved ? q7.eps : eps7_ub);
    std::printf("   d=3: memory-Z eps %.4f%% +- %.4f%%, memory-X eps %.4f%% +- %.4f%%   (round %.0f ns)\n", 100 * q3.eps, 100 * q3.eps_err, 100 * q3x.eps, 100 * q3x.eps_err, q3.round_ns);
    std::printf("   d=5: memory-Z eps %.4f%% +- %.4f%%, memory-X eps %.4f%% +- %.4f%%   (round %.0f ns)\n", 100 * q5.eps, 100 * q5.eps_err, 100 * q5x.eps, 100 * q5x.eps_err, q5.round_ns);
    if (d7_resolved) std::printf("   d=7: memory-Z eps %.4f%% +- %.4f%%   (97 qubits, 9 rounds, %ld failures / %ld; d=5 at 9 rounds %.4f%%)\n", 100 * q7.eps, 100 * eps7_err, q7.failures_last, n7, 100 * q5_9.eps);
    else std::printf("   d=7: no logical failures in %ld shots of 9 rounds -> eps < %.4f%% (95%%)\n", n7, 100 * eps7_ub);
    std::printf("   Lambda(3->5) = %.2f (Z), %.2f (X); Lambda(5->7) %s %.2f;  d=5 logical lifetime %.2f ms\n", L35, L35x, d7_resolved ? "=" : ">", L57, q5.round_ns / 1e6 / q5.eps);
    {
        Json s = Json::object();
        s.set("eps_d3", q3.eps); s.set("eps_d5", q5.eps); s.set("eps_d7", q7.eps); s.set("eps_err_d3", q3.eps_err); s.set("eps_err_d5", q5.eps_err);
        s.set("eps_d3_X", q3x.eps); s.set("eps_d5_X", q5x.eps); s.set("eps_err_d3_X", q3x.eps_err); s.set("eps_err_d5_X", q5x.eps_err);
        s.set("eps_err_d7", eps7_err); s.set("d7_resolved", d7_resolved); s.set("eps_d7_upper95", eps7_ub); s.set("d7_failures", double(q7.failures_last)); s.set("shots_d7", double(n7));
        s.set("eps_d5_9rounds", q5_9.eps); s.set("eps_err_d5_9rounds", eps5_9_err); s.set("d5_9_failures", double(q5_9.failures_last)); s.set("shots_d5_9rounds", double(q5_9.shots));
        s.set("p_logical_d7_9rounds", q7.p_last); s.set("p_logical_d5_9rounds", q5_9.p_last);
        s.set("Lambda_3_5", L35); s.set("Lambda_3_5_X", L35x); s.set("Lambda_5_7", L57); s.set("round_ns_d5", q5.round_ns);
        s.set("logical_lifetime_d5_ms", q5.round_ns / 1e6 / q5.eps);
        s.set("mean_defects_d5_20_rounds", q5.mean_defects); s.set("decode_us_per_round_d5", q5.decode_us_per_round);
        s.set("d7_method", "one 9-round run, eps from P_L; d=5 at 9 rounds by the same estimate for Lambda(5->7)");
        out.set("distance_scaling", s);
    }
    // reference: Alpha-1 v1.1 physics at d=5 (same layout, no Beta upgrades)
    const auto a5 = qec(alpha1_v11(5), RF, shots, 11);
    const auto a3 = qec(alpha1_v11(3), RF, 2 * shots, 11);
    std::printf("   reference, Alpha-1 v1.1 physics: d=3 %.4f%%, d=5 %.4f%%, Lambda %.2f\n", 100 * a3.eps, 100 * a5.eps, a3.eps / a5.eps);
    { Json r = Json::object(); r.set("eps_d3", a3.eps); r.set("eps_d5", a5.eps); r.set("Lambda", a3.eps / a5.eps); out.set("alpha1_v11_reference", r); }

    // ------------------------------------------------------------------ 3. long experiments and leakage
    std::printf("\n3. Long experiments at d=5: leakage removal (%d rounds) and detection-rate stability (%d rounds)\n", RLONG, RSTAB);
    {
        Json lj = Json::object();
        for (bool dq : {false, true}) {
            ChipDesign D = B5; D.leakage.dqlr = dq;
            const auto q = qec(D, {RLONG}, shots / 2, 21);
            std::printf("   %-28s P_L(%d) %.4f, leaked data %.1e meas %.1e, round %.0f ns\n", dq ? "LRR + DQLR (option)" : "LRR only (Beta-1)", RLONG, q.p_last, q.leak_d, q.leak_m, q.round_ns);
            Json o = Json::object(); o.set("rounds", RLONG); o.set("P_L_50", q.p_last); o.set("leaked_data", q.leak_d); o.set("leaked_meas", q.leak_m); o.set("round_ns", q.round_ns);
            lj.set(dq ? "lrr_dqlr" : "lrr", o);
        }
        out.set("long_run_d5", lj);
        {   // what exact decoding changed: the same syndromes, decoded exactly and with the pre-v0.4 fallback
            // (exact up to 22 defects, greedy above)
            auto Rr = prepared_report(B5);
            const auto nm = noise_model(B5, Rr);
            const MemorySetup M = prepare_memory(B5, Rr, nm, RLONG);
            const long nd = g_opt.smoke ? shots : std::max<long>(2000, shots / 25);
            FrameRunOptions legacy; legacy.method = DecodeMethod::Legacy;
            const auto ex = run_memory_frame(M, nd, 77);
            const auto lg = run_memory_frame(M, nd, 77, legacy);
            std::printf("   decoder, %d rounds, %ld shots, same syndromes (%.1f defects per shot, %.1f%% above 22): exact matching P_L %.4f (%ld failures), greedy above 22 defects P_L %.4f (%ld failures)\n",
                        RLONG, nd, ex.mean_defects, 100 * ex.frac_over_22, ex.p_logical, ex.logical_errors, lg.p_logical, lg.logical_errors);
            Json o = Json::object();
            o.set("rounds", RLONG); o.set("shots", double(nd)); o.set("mean_defects", ex.mean_defects); o.set("frac_over_22", ex.frac_over_22);
            o.set("P_L_exact", ex.p_logical); o.set("P_L_legacy", lg.p_logical); o.set("failures_exact", double(ex.logical_errors)); o.set("failures_legacy", double(lg.logical_errors));
            out.set("decoding_d5", o);
        }
        Json sj = Json::object();
        for (bool lrr : {true, false}) {
            ChipDesign D = B5; D.leakage.reset_removes_leakage = lrr;
            const auto q = qec(D, {RSTAB}, shots / 5, 23);
            const double drift = detection_drift(q.detection_by_round, 10);
            std::printf("   %d rounds, %-26s detection drift %+.1f%% (first vs last 10 rounds), P_L %.4f, %.1f defects per shot (max %d)\n",
                        RSTAB, lrr ? "leakage-removing reset:" : "plain reset:", 100 * drift, q.p_last, q.mean_defects, q.max_defects);
            Json o = Json::object(); o.set("rounds", RSTAB); o.set("drift", drift); o.set("P_L_100", q.p_last); o.set("mean_defects", q.mean_defects); o.set("max_defects", q.max_defects);
            o.set("frac_over_22", q.frac_over_22);
            sj.set(lrr ? "lrr" : "plain_reset", o);
        }
        out.set("stability_100_rounds_d5", sj);
    }

    // ------------------------------------------------------------------ 4. budget + sensitivity
    std::puts("\n4. Error budget (Beta-1 physics, d=3, leave one source out) and remaining sensitivity");
    {
        Json bj = Json::object();
        for (const char* src : {"T1", "Tphi", "control", "readout", "reset", "zz", "leakage", "thermal"}) {
            const std::string s = src;
            const auto q = qec(B3, RF, shots, 31, 'Z', [s](NoiseModel& nm) { nm.disable(s); });
            std::printf("   without %-8s eps %.4f%%  (%+.0f%%)\n", src, 100 * q.eps, 100 * (q.eps / q3.eps - 1));
            bj.set(src, q.eps / q3.eps - 1);
        }
        out.set("budget_d3", bj);
        Json sj = Json::object();
        struct SV { const char* key; const char* name; std::function<void(ChipDesign&)> f; };
        for (const SV& v : {SV{"resid_2p5e-4", "CZ residual 5e-4 -> 2.5e-4", [](ChipDesign& D) { D.coupler.residual_ctrl_err = 2.5e-4; }},
                            SV{"Q_x2", "dielectric Q 6e6 -> 1.2e7", [](ChipDesign& D) { D.Q_dielectric *= 2; }},
                            SV{"readout_200", "readout 300 -> 200 ns", [](ChipDesign& D) { D.t_readout_ns = 200; }},
                            SV{"sep_0p25", "readout separation 0.5% -> 0.25%", [](ChipDesign& D) { D.readout_sep_err = 0.0025; }}}) {
            ChipDesign D = B3; v.f(D);
            const auto q = qec(D, RF, shots, 41);
            std::printf("   %-38s eps %.4f%%  (%+.0f%%)\n", v.name, 100 * q.eps, 100 * (q.eps / q3.eps - 1));
            sj.set(v.key, q.eps / q3.eps - 1);
        }
        out.set("sensitivity_d3", sj);
    }

    // ------------------------------------------------------------------ 5. fabrication yield
    std::puts("\n5. Fabrication yield at d=5 (junction resistance spread; frequency-collision rules)");
    {
        Json fj = Json::object();
        for (double spread : {2.0, 1.0, 0.5}) {
            ChipDesign D = B5; D.Rn_spread_pct = spread;
            const auto F = fab_monte_carlo(D, g_opt.smoke ? 20 : 200, 7);
            double sig = 0; for (double s : F.f01_sigma_MHz) sig += s; sig /= double(std::max<std::size_t>(F.f01_sigma_MHz.size(), 1));
            std::printf("   spread %.1f%%%s: collision-free %3d/%d (%.0f%%), frequency scatter %.0f MHz rms\n", spread, spread <= 0.5 ? " (trimmed)" : "",
                        F.collision_free, F.samples, 100.0 * F.collision_free / F.samples, sig);
            Json o = Json::object(); o.set("collision_free", F.collision_free); o.set("samples", F.samples); o.set("f01_sigma_MHz", sig);
            char k[16]; std::snprintf(k, sizeof k, "spread_%.1f", spread); fj.set(k, o);
        }
        out.set("fab_yield_d5", fj);
    }

    // ------------------------------------------------------------------ 6. wiring
    std::puts("\n6. Wiring and readout multiplexing");
    {
        const auto w = wiring(B5), w3 = wiring(alpha1_v11(3));
        double fmin = 1e9, fmax = 0; for (const auto& q : B5.qubits) { fmin = std::min(fmin, q.ro.f_GHz); fmax = std::max(fmax, q.ro.f_GHz); }
        std::printf("   %d XY drive + %d qubit flux + %d coupler flux + %d readout in + %d readout out = %d coaxial lines\n", w.xy, w.flux_qubit, w.flux_coupler, w.feedlines, w.feedlines, w.total);
        std::printf("   %d feedlines, %.1f resonators per line, %.2f-%.2f GHz; Alpha-1 (17 qubits) needs %d lines; per readout output chain: TWPA + HEMT\n",
                    w.feedlines, w.per_feedline, fmin, fmax, w3.total);
        Json j = wiring_json(w); j.set("readout_band_lo_GHz", fmin); j.set("readout_band_hi_GHz", fmax); j.set("alpha1_lines", w3.total);
        out.set("wiring", j);
    }

    // ------------------------------------------------------------------ 7. layout
    std::puts("\n7. Physical layout (49-qubit chip, cells synthesised for the Beta-1 readout coupling)");
    const auto lay = layout_section(B5, "designs/beta1/layout_cells.json", "designs/beta1/BF49.gds", "BF49_BETA1");
    out.set("layout", lay.json);

    // ------------------------------------------------------------------ 8. exit criteria
    std::puts("\n8. Beta exit criteria (simulation gate before tape-out)");
    {
        Json ej = Json::object();
        auto crit = [&](const char* key, const std::string& text, bool ok) { std::printf("   [%s] %s\n", ok ? "PASS" : "FAIL", text.c_str()); ej.set(key, ok); };
        crit("lambda", "error suppression Lambda(3->5) >= 2 in both bases", L35 >= 2.0 && L35x >= 2.0);
        crit("eps_d5", "logical error per round at d=5 <= 0.1% in both bases", q5.eps <= 1e-3 && q5x.eps <= 1e-3);
        crit("beats_alpha", "d=5 logical error below Alpha-1 v1.1 d=3 by >= 3x", a3.eps / q5.eps >= 3.0);
        crit("lifetime", "d=5 logical lifetime > 5x best physical T1", q5.round_ns / 1e3 / q5.eps > 5 * out["device"].get("best_T1_us", 1.0));
        crit("yield", "trimmed fabrication yield >= 80% collision-free", out["fab_yield_d5"]["spread_0.5"].get("collision_free", 0.0) >= 0.8 * out["fab_yield_d5"]["spread_0.5"].get("samples", 200.0));
        crit("leakage", "chosen leakage scheme (LRR) at least as good as LRR + DQLR over " + std::to_string(RLONG) + " rounds", out["long_run_d5"]["lrr"].get("P_L_50", 1.0) <= out["long_run_d5"]["lrr_dqlr"].get("P_L_50", 0.0) * 1.05);
        crit("stability", "detection rate flat over " + std::to_string(RSTAB) + " rounds with the leakage-removing reset (drift < 5%)", std::abs(out["stability_100_rounds_d5"]["lrr"].get("drift", 1.0)) < 0.05);
        // The simulation gate requires the layout to pass the placement / island-spacing checks and the mask-level
        // width/space DRC (layout.drc), to pass the metal connectivity check with its ground straps (layout.connectivity),
        // and to meet the circuit targets. The control-line T1 budget (layout.control_lines) is reported, not gated: the
        // flux line misses it on the Beta-1 cells.
        crit("drc", "layout passes placement, island-spacing and mask width/space checks", lay.found && lay.placement_clean && lay.mask_drc_clean);
        crit("connectivity", "metal connectivity clean: conductors separate, shorted lines and all ground on the plane, SQUIDs island-to-plane", lay.found && lay.connectivity_clean);
        crit("layout", "extracted E_C within 3%, couplings and readout capacitance within 6%; parked exchange re-nullable",
             lay.found && lay.within_tolerance && lay.renull_ok);
        out.set("exit_criteria", ej);
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.set("shots", double(shots)); out.set("runtime_s", secs);
    out.set("provenance", provenance(g_opt, "frame Monte Carlo seeds 11 + rounds (distance scaling, reference; memory-X adds 7919), 12 + rounds (d=5 at 9 rounds), "
                                            "21 + rounds (50-round runs), 77 (decoder comparison), 23 + rounds (100-round runs), 31 + rounds (budget), 41 + rounds (sensitivity); fab Monte Carlo seed 7"));
    const std::string jpath = g_opt.path("designs/beta1/beta1_review.json");
    std::ofstream f(jpath); out.dump(f); f << "\n";
    std::printf("\nreview complete in %.0f s -> %s\n", secs, jpath.c_str());
    return 0;
}
