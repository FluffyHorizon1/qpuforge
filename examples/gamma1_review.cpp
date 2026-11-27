// Gamma-1 (GF-97) design review — distance-7 surface-code memory, 97 qubits
//   1. device summary and pulse-level CZ calibration of all pairs
//   2. error suppression with code distance in both memory bases: eps at d = 3, 5, 7 (same estimator), Lambda
//   3. decoding: defects per shot, exact matching at any defect count, software decode time against the round time,
//      and what the pre-v0.4 fallback (greedy above 22 defects) costs on a long experiment
//   4. high-energy impacts: logical-error floor, unprotected time, visibility in the syndrome
//   5. long experiment at d = 7: detection rate per round with and without the leakage-removing reset
//   6. fabrication yield at 97 qubits
//   7. wiring and readout multiplexing
//   8. physical layout (Beta-1 rev C cells), DRC and the post-layout circuit check
//   9. exit criteria (Gamma memory gate; lattice surgery is reviewed separately)
// Writes designs/gamma1/gamma1_review.json and prints a report. Options: see review_io.hpp ([shots] [--out DIR] [--cache FILE] [--smoke]).
#include "review_util.hpp"
#include <chrono>

using namespace qpuforge;
using namespace review;

int main(int argc, char** argv) {
    g_opt = parse(argc, argv, 100000);                              // d = 7 runs use 2x the shots, burst ladders 0.2x
    const long shots = g_opt.shots;
    const auto t0 = std::chrono::steady_clock::now();
    Json out = Json::object();
    const ChipDesign G7 = gamma1(7), G5 = gamma1(5), G3 = gamma1(3);
    save_design(G7, g_opt.path("designs/gamma1/gamma1.json"));
    // fixed-size runs: at least `n` shots in a full review, the smoke-run shot count otherwise
    auto at_least = [&](long n, long scaled) { return g_opt.smoke ? shots : std::max(n, scaled); };
    const int R = 14;                                              // rounds of every memory run in section 2 (= 2 d at d = 7)

    // ------------------------------------------------------------------ 1. device
    std::puts("1. Gamma-1 device (97 qubits: 49 data + 48 measure, tunable couplers)");
    auto R7 = analyze(G7);
    const auto cal = calibrate_chip_cz(G7, R7, true, g_opt.cache);
    double best_T1 = 0;
    {
        double T1 = 0, T2 = 0, ro = 0, t = 0, e = 0, lk = 0, lkmax = 0;
        for (const auto& q : R7.q) { T1 += q.T1_us; T2 += q.T2_us; ro += 0.5 * (q.ro_err_0to1 + q.ro_err_1to0); best_T1 = std::max(best_T1, q.T1_us); }
        for (const auto& c : cal) { t += c.t_cz_ns; e += c.err_total; lk += c.leakage; lkmax = std::max(lkmax, c.leakage); }
        const double n = double(R7.q.size()), m = double(cal.size());
        std::printf("   mean T1 %.0f us (best %.0f), T2 %.0f us, readout error %.2f%% at %.0f ns; %zu CZ pairs: %.1f ns, error %.2e, leakage %.1e (max %.1e); collisions %d\n",
                    T1 / n, best_T1, T2 / n, 100 * ro / n, G7.t_readout_ns, cal.size(), t / m, e / m, lk / m, lkmax, R7.n_collisions);
        Json d = Json::object();
        d.set("qubits", int(n)); d.set("couplers", int(m)); d.set("mean_T1_us", T1 / n); d.set("best_T1_us", best_T1); d.set("mean_T2_us", T2 / n);
        d.set("mean_readout_err", ro / n); d.set("t_cz_mean_ns", t / m); d.set("cz_err_mean", e / m); d.set("cz_leak_mean", lk / m); d.set("cz_leak_max", lkmax); d.set("collisions", R7.n_collisions);
        out.set("device", d);
    }

    // ------------------------------------------------------------------ 2. distance scaling, both bases
    std::printf("\n2. Logical error per round vs code distance (%d-round memory, both bases, exact matching decoder)\n", R);
    struct Row { int d; QecPoint z, x; };
    std::vector<Row> rows;
    {
        Json sj = Json::object();
        for (int d : {3, 5, 7}) {
            const ChipDesign& D = d == 3 ? G3 : d == 5 ? G5 : G7;
            const long n = d == 7 ? 2 * shots : shots;
            Row r{d, qec(D, {R}, n, 11, 'Z'), qec(D, {R}, n, 11, 'X')};
            rows.push_back(r);
            std::printf("   d=%d: memory-Z eps %.5f%% +- %.5f%% (%ld failures / %ld), memory-X eps %.5f%% +- %.5f%% (%ld / %ld), round %.0f ns\n", d,
                        100 * r.z.eps, 100 * r.z.eps_err, r.z.failures_last, r.z.shots, 100 * r.x.eps, 100 * r.x.eps_err, r.x.failures_last, r.x.shots, r.z.round_ns);
            Json o = Json::object();
            o.set("eps_Z", r.z.eps); o.set("eps_Z_err", r.z.eps_err); o.set("failures_Z", double(r.z.failures_last)); o.set("shots_Z", double(r.z.shots));
            o.set("eps_X", r.x.eps); o.set("eps_X_err", r.x.eps_err); o.set("failures_X", double(r.x.failures_last)); o.set("shots_X", double(r.x.shots));
            o.set("round_ns", r.z.round_ns); o.set("detection_rate", r.z.det);
            sj.set("d" + std::to_string(d), o);
        }
        auto lam = [](double a, double b) { return b > 0 ? a / b : 0.0; };
        const double L35z = lam(rows[0].z.eps, rows[1].z.eps), L57z = lam(rows[1].z.eps, rows[2].z.eps);
        const double L35x = lam(rows[0].x.eps, rows[1].x.eps), L57x = lam(rows[1].x.eps, rows[2].x.eps);
        const double eps7 = 0.5 * (rows[2].z.eps + rows[2].x.eps);
        const double life_ms = rows[2].z.round_ns / 1e6 / std::max(eps7, 1e-12);
        std::printf("   Lambda(3->5): Z %.2f, X %.2f;  Lambda(5->7): Z %.2f, X %.2f;  d=7 logical lifetime %.1f ms (%.0fx best physical T1)\n",
                    L35z, L35x, L57z, L57x, life_ms, life_ms * 1e3 / best_T1);
        sj.set("rounds", R); sj.set("Lambda_3_5_Z", L35z); sj.set("Lambda_3_5_X", L35x); sj.set("Lambda_5_7_Z", L57z); sj.set("Lambda_5_7_X", L57x);
        sj.set("eps_d7_mean", eps7); sj.set("logical_lifetime_d7_ms", life_ms); sj.set("lifetime_over_best_T1", life_ms * 1e3 / best_T1);
        out.set("distance_scaling", sj);
    }

    // ------------------------------------------------------------------ 3. decoding
    std::puts("\n3. Decoding (exact minimum-weight matching: DP up to 8 defects, blossom beyond)");
    // Long experiment, shared with section 5. Setting it up costs time quadratic in the number of rounds, so a smoke
    // run shortens it; the JSON records the round count used.
    const int RL = g_opt.smoke ? 16 : 40;
    const auto nm7 = noise_model(G7, R7);
    const MemorySetup ML = prepare_memory(G7, R7, nm7, RL);
    {
        Json dj = Json::object();
        const auto& q7 = rows[2].z;
        std::printf("   d=7, %d rounds: %.1f defects per shot (max %d), %.3f%% of shots above 22; software decode %.2f us per round on one core (round %.2f us)\n",
                    R, q7.mean_defects, q7.max_defects, 100 * q7.frac_over_22, q7.decode_us_per_round, q7.round_ns / 1e3);
        dj.set("mean_defects_d7", q7.mean_defects); dj.set("max_defects_d7", q7.max_defects); dj.set("frac_over_22_d7", q7.frac_over_22);
        dj.set("decode_us_per_round_d7", q7.decode_us_per_round); dj.set("round_us_d7", q7.round_ns / 1e3);
        // long experiment: the same syndromes decoded exactly and with the pre-v0.4 fallback
        const long nl = at_least(2000, shots / 25);
        FrameRunOptions legacy; legacy.method = DecodeMethod::Legacy;
        const auto ex = run_memory_frame(ML, nl, 77);
        const auto lg = run_memory_frame(ML, nl, 77, legacy);
        std::printf("   d=7, %d rounds (%ld shots, same syndromes): %.1f defects per shot, %.1f%% above 22\n", RL, nl, ex.mean_defects, 100 * ex.frac_over_22);
        std::printf("     exact matching          P_L %.5f (%ld failures), decode %.1f us per round\n", ex.p_logical, ex.logical_errors, ex.decode_us_per_shot / RL);
        std::printf("     greedy above 22 defects P_L %.5f (%ld failures), decode %.1f us per round\n", lg.p_logical, lg.logical_errors, lg.decode_us_per_shot / RL);
        Json lj = Json::object();
        lj.set("rounds", RL); lj.set("shots", double(nl)); lj.set("mean_defects", ex.mean_defects); lj.set("frac_over_22", ex.frac_over_22);
        lj.set("P_L_exact", ex.p_logical); lj.set("P_L_legacy", lg.p_logical); lj.set("failures_exact", double(ex.logical_errors)); lj.set("failures_legacy", double(lg.logical_errors));
        lj.set("decode_us_per_round_exact", ex.decode_us_per_shot / RL); lj.set("decode_us_per_round_legacy", lg.decode_us_per_shot / RL);
        dj.set("long_run", lj);
        out.set("decoding", dj);
    }

    // ------------------------------------------------------------------ 4. high-energy impacts
    std::puts("\n4. High-energy impacts (correlated error bursts; the decoder keeps its quiet-chip priors)");
    {
        Json bj = Json::object();
        const auto& nm = nm7;
        const double eps7 = out["distance_scaling"].get("eps_d7_mean", 0.0);
        const double lam = 0.5 * (out["distance_scaling"].get("Lambda_5_7_Z", 0.0) + out["distance_scaling"].get("Lambda_5_7_X", 0.0));
        struct Sc { const char* key; const char* name; BurstSpec B; char basis; };
        BurstSpec raw; raw.rate_hz = 0.1; raw.tau_ms = 25; raw.t1_peak_us = 1.0; raw.tphi_peak_us = 0;
        BurstSpec gap; gap.rate_hz = 0.1; gap.tau_ms = 1; gap.t1_peak_us = 0; gap.tphi_peak_us = 0.1;
        for (const Sc& sc : {Sc{"unmitigated", "no mitigation: T1 burst, 25 ms recovery (memory-Z)", raw, 'Z'},
                             Sc{"gap_engineered", "gap-engineered junctions: dephasing burst, 1 ms (memory-X)", gap, 'X'}}) {
            const auto b = burst_floor(G7, R7, nm, sc.B, R, at_least(4000, shots / 5), 505, sc.basis);
            const double d_star = (lam > 1 && b.floor_per_round > 0 && eps7 > b.floor_per_round) ? 7 + 2 * std::log(eps7 / b.floor_per_round) / std::log(lam) : 0.0;
            std::printf("   %s\n", sc.name);
            std::printf("     quiet chip: eps %.6f per round, %.1f defects per shot\n", b.eps0, b.mean_defects0);
            for (const auto& pt : b.curve) {
                if (pt.saturated) std::printf("     severity %.1e: logical qubit randomised within the run (eps >= %.3f), %.1f defects per shot\n", pt.s, eps_from_PL(0.45, R), pt.mean_defects);
                else std::printf("     severity %.1e: eps %.6f per round, %.1f defects per shot\n", pt.s, pt.eps, pt.mean_defects);
            }
            std::printf("     one impact flips the logical qubit with probability %.2f; unprotected for %.1f ms; syndrome above 1.5x quiet for %.1f ms\n",
                        b.p_fail, b.dead_time_ms, b.detectable_ms);
            std::printf("     floor %.2e per round = %.2f%% of the d=7 error; impacts would dominate from about d = %.0f at this Lambda\n",
                        b.floor_per_round, 100 * b.floor_per_round / std::max(eps7, 1e-30), d_star);
            Json o = Json::object();
            o.set("rate_hz", sc.B.rate_hz); o.set("tau_ms", sc.B.tau_ms); o.set("t1_peak_us", sc.B.t1_peak_us); o.set("tphi_peak_us", sc.B.tphi_peak_us);
            o.set("basis", std::string(1, sc.basis)); o.set("p_fail", b.p_fail); o.set("floor_per_round", b.floor_per_round);
            o.set("floor_over_eps_d7", b.floor_per_round / std::max(eps7, 1e-30)); o.set("dead_time_ms", b.dead_time_ms); o.set("detectable_ms", b.detectable_ms);
            o.set("lost_time_fraction", sc.B.rate_hz * b.dead_time_ms * 1e-3); o.set("distance_where_impacts_dominate", d_star);
            Json cv = Json::array();
            for (const auto& pt : b.curve) { Json p = Json::object(); p.set("severity", pt.s); p.set("eps", pt.eps); p.set("saturated", pt.saturated); p.set("defects", pt.mean_defects); cv.push(p); }
            o.set("eps_quiet", b.eps0); o.set("defects_quiet", b.mean_defects0);
            o.set("curve", cv);
            bj.set(sc.key, o);
        }
        out.set("impacts", bj);
    }

    // ------------------------------------------------------------------ 5. long experiment and leakage
    std::printf("\n5. Long experiment at d=7 (%d rounds): detection rate per round\n", RL);
    {
        Json lj = Json::object();
        const long nl = at_least(4000, shots / 2);
        for (bool lrr : {true, false}) {
            // The reset type changes what happens to a leaked qubit, not the Pauli faults: both cases share one
            // detector error model and one decoder.
            ChipDesign D = G7; D.leakage.reset_removes_leakage = lrr;
            const auto nm = noise_model(D, R7);
            const MemorySetup M = prepare_memory(D, R7, nm, RL, 'Z', false);
            FrameRunOptions opt; opt.decoder = ML.decoder.get();
            const auto q = run_memory_frame(M, nl, 21 + std::uint64_t(RL), opt);
            const double drift = detection_drift(q.detection_by_round, 5);
            const auto& dr = q.detection_by_round;
            std::printf("   %-34s P_L(%d) %.5f (%ld failures / %ld), detection %.4f -> %.4f (drift %+.1f%%), leaked data %.1e meas %.1e\n", lrr ? "leakage-removing reset (Gamma-1)" : "plain reset (no leakage removal)",
                        RL, q.p_logical, q.logical_errors, q.shots, dr[1], dr[dr.size() - 2], 100 * drift, q.leak_frac_data, q.leak_frac_meas);
            Json o = Json::object(); o.set("rounds", RL); o.set("P_L_40", q.p_logical); o.set("failures", double(q.logical_errors)); o.set("shots", double(q.shots)); o.set("drift", drift);
            o.set("leaked_data", q.leak_frac_data); o.set("leaked_meas", q.leak_frac_meas);
            Json a = Json::array(); for (double v : dr) a.push(v);
            o.set("detection_by_round", a);
            lj.set(lrr ? "lrr" : "plain_reset", o);
        }
        out.set("long_run_d7", lj);
    }

    // ------------------------------------------------------------------ 6. fabrication yield
    std::puts("\n6. Fabrication yield at 97 qubits (junction resistance spread; frequency-collision rules)");
    {
        Json fj = Json::object();
        for (double spread : {1.0, 0.5}) {
            ChipDesign D = G7; D.Rn_spread_pct = spread;
            const auto F = fab_monte_carlo(D, g_opt.smoke ? 20 : 100, 7);
            double sig = 0; for (double s : F.f01_sigma_MHz) sig += s; sig /= double(std::max<std::size_t>(F.f01_sigma_MHz.size(), 1));
            std::printf("   spread %.1f%%%s: collision-free %3d/%d (%.0f%%), frequency scatter %.0f MHz rms\n", spread, spread <= 0.5 ? " (trimmed)" : "",
                        F.collision_free, F.samples, 100.0 * F.collision_free / F.samples, sig);
            Json o = Json::object(); o.set("collision_free", F.collision_free); o.set("samples", F.samples); o.set("f01_sigma_MHz", sig);
            char k[16]; std::snprintf(k, sizeof k, "spread_%.1f", spread); fj.set(k, o);
        }
        out.set("fab_yield_d7", fj);
    }

    // ------------------------------------------------------------------ 7. wiring
    std::puts("\n7. Wiring and readout multiplexing");
    {
        const auto w = wiring(G7), w5 = wiring(gamma1(5)), w3 = wiring(alpha1_v11(3));
        double fmin = 1e9, fmax = 0; for (const auto& q : G7.qubits) { fmin = std::min(fmin, q.ro.f_GHz); fmax = std::max(fmax, q.ro.f_GHz); }
        std::printf("   %d XY drive + %d qubit flux + %d coupler flux + %d readout in + %d readout out = %d coaxial lines\n", w.xy, w.flux_qubit, w.flux_coupler, w.feedlines, w.feedlines, w.total);
        std::printf("   %d feedlines, %.1f resonators per line, %.2f-%.2f GHz; Beta-1 needs %d lines, Alpha-1 %d\n", w.feedlines, w.per_feedline, fmin, fmax, w5.total, w3.total);
        Json j = wiring_json(w); j.set("readout_band_lo_GHz", fmin); j.set("readout_band_hi_GHz", fmax); j.set("beta1_lines", w5.total); j.set("alpha1_lines", w3.total);
        out.set("wiring", j);
    }

    // ------------------------------------------------------------------ 8. layout
    std::puts("\n8. Physical layout (97-qubit chip from the Beta-1 rev C cells)");
    const auto lay = layout_section(G7, "designs/beta1/layout_cells.json", "designs/gamma1/GF97.gds", "GF97_GAMMA1");
    out.set("layout", lay.json);

    // ------------------------------------------------------------------ 9. exit criteria
    std::puts("\n9. Gamma memory gate (simulation, before tape-out)");
    {
        Json ej = Json::object();
        auto crit = [&](const char* key, const std::string& text, bool ok) { std::printf("   [%s] %s\n", ok ? "PASS" : "FAIL", text.c_str()); ej.set(key, ok); };
        const Json& s = out["distance_scaling"];
        crit("lambda", "Lambda(5->7) >= 2 in both bases", s.get("Lambda_5_7_Z", 0.0) >= 2 && s.get("Lambda_5_7_X", 0.0) >= 2);
        crit("eps_d7", "distance-7 logical error per round <= 0.01% in both bases", rows[2].z.eps <= 1e-4 && rows[2].x.eps <= 1e-4);
        crit("lifetime", "distance-7 logical lifetime > 50x best physical T1", s.get("lifetime_over_best_T1", 0.0) > 50);
        crit("impacts", "logical-error floor from impacts below 10% of the distance-7 error (both scenarios)",
             out["impacts"]["unmitigated"].get("floor_over_eps_d7", 1.0) < 0.1 && out["impacts"]["gap_engineered"].get("floor_over_eps_d7", 1.0) < 0.1);
        crit("stability", "detection rate flat over " + std::to_string(RL) + " rounds with the leakage-removing reset (drift < 5%)", std::abs(out["long_run_d7"]["lrr"].get("drift", 1.0)) < 0.05);
        crit("yield", "trimmed fabrication yield >= 80% collision-free", out["fab_yield_d7"]["spread_0.5"].get("collision_free", 0.0) >= 0.8 * out["fab_yield_d7"]["spread_0.5"].get("samples", 100.0));
        // As for Beta-1: placement / island-spacing checks, mask-level width/space DRC and circuit targets are gated.
        crit("drc", "layout passes placement, island-spacing and mask width/space DRC (4 um / 4 um, merged etch layer)", lay.found && lay.mask_drc_clean);
        crit("layout", "extracted E_C within 3%, couplings and readout capacitance within 6%; parked exchange re-nullable",
             lay.found && lay.within_tolerance && lay.renull_ok);
        out.set("exit_criteria", ej);
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.set("shots", double(shots)); out.set("runtime_s", secs);
    out.set("provenance", provenance(g_opt, "frame Monte Carlo seeds 11 + rounds (distance scaling; memory-X adds 7919), 77 (decoder comparison), 505 (impacts), "
                                            "21 + rounds (long runs); fab Monte Carlo seed 7"));
    const std::string jpath = g_opt.path("designs/gamma1/gamma1_review.json");
    std::ofstream f(jpath); out.dump(f); f << "\n";
    std::printf("\nreview complete in %.0f s -> %s\n", secs, jpath.c_str());
    return 0;
}
