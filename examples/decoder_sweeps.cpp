// Parameter sweeps behind the design choices of the Phase Gamma decoding work, on the d = 5 test chip of
// tests/test_qpuforge.cpp (surface_chip):
//   1. herald rule and factor (docs/design/heralded-leakage.md): 15 rounds, measure-qubit CZ leakage 2e-3 per CZ,
//      transport 0.05, leakage-removing reset; the same syndromes decoded with each rule
//   2. veto window shape (docs/design/impact-veto.md): 10 rounds, chip-wide T1 impacts in 5% of shots
//   3. latched-veto hold time: burst_floor_veto over the T1 burst profile (tau 25 ms)
// Writes results/decoder_sweeps.{txt,json}. Options as for the reviews ([shots] [--out DIR] [--smoke]); the shot
// counts below are fixed (a smoke run divides them by 20).
#include "review_util.hpp"
#include <cstdarg>

using namespace qpuforge;
using namespace review;

static std::string g_txt;
static void say(const char* fmt, ...) {
    char buf[4096];
    va_list ap; va_start(ap, fmt); std::vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    std::fputs(buf, stdout); std::fflush(stdout);
    g_txt += buf;
}

static ChipDesign surface_chip(int d) {
    GenSpec s; s.topology = "surface"; s.distance = d; s.architecture = "tunable-coupler-cz"; s.f_hi = 5.9; s.f_lo = 5.2; s.g = 0.004;
    s.feedlines = 3; s.kappa_MHz = 4; s.chi_MHz = 4; s.purcell_filter = 100;
    auto D = generate(s);
    D.t_readout_ns = 500; D.active_reset = true; D.reset_err = 0.003; D.t_reset_ns = 200;
    return D;
}

int main(int argc, char** argv) {
    g_opt = parse(argc, argv, 0);
    const long div = g_opt.smoke ? 20 : 1;
    const auto t_all = std::chrono::steady_clock::now();
    Json out = Json::object();
    const auto D = surface_chip(5); const auto R = analyze(D);

    // ------------------------------------------------------------------ 1. herald rule and factor
    {
        auto nm = noise_model(D, R);
        nm.leakage = true; nm.leak_1q = 0; nm.leak_transport = 0.05; nm.reset_removes_leakage = true;
        for (const auto& p : R.pairs) if (p.designed) {
            const int a = D.index_of(p.a), b = D.index_of(p.b);
            const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
            nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
        }
        const auto M = prepare_memory(D, R, nm, 15);
        const long n = 40000 / div;
        say("1. Herald rules (d = 5 test chip, 15 rounds, CZ leakage 2e-3, leakage-removing reset, %ld shots, exact matching, same syndromes)\n", n);
        struct Case { const char* name; bool use, footprint; double erase, foot, touch; };
        const Case cases[] = {
            {"no heralds", false, true, 0, 0.4, 1},
            {"detector rule: every touching edge erased", true, false, 0, 1, 0},
            {"detector rule: measurement edge erased only", true, false, 0, 1, 1},
            {"footprint rule, factor 1 (measurement edge erased only)", true, true, 0, 1.0, 1},
            {"footprint rule, factor 0.5", true, true, 0, 0.5, 1},
            {"footprint rule, factor 0.4 (default)", true, true, 0, 0.4, 1},
            {"footprint rule, factor 0.3", true, true, 0, 0.3, 1},
            {"footprint rule, factor 0.15", true, true, 0, 0.15, 1},
            {"footprint rule, factor 0", true, true, 0, 0.0, 1},
            {"footprint rule, factor 0.4, measurement edge x0.4", true, true, 0.4, 0.4, 1},
        };
        Json hj = Json::array();
        for (const auto& c : cases) {
            FrameRunOptions o; o.herald.on = true; o.herald.use_in_decoder = c.use; o.herald.footprint = c.footprint;
            o.herald.erase_factor = c.erase; o.herald.footprint_factor = c.foot; o.herald.touch_factor = c.touch;
            const auto r = run_memory_frame(M, n, 77, o);
            say("   %-56s eps %.5f (%ld failures)\n", c.name, r.eps_per_round, r.logical_errors);
            Json j = Json::object(); j.set("case", std::string(c.name)); j.set("eps", r.eps_per_round); j.set("failures", double(r.logical_errors)); hj.push(j);
        }
        out.set("herald_rules", hj);
    }

    const auto nm = noise_model(D, R);
    const auto M = prepare_memory(D, R, nm, 10);
    // ------------------------------------------------------------------ 2. veto window shape
    {
        const long nc = 50000 / div, ni = 20000 / div, ns = 10000 / div;
        say("\n2. Veto window shape (d = 5 test chip, 10 rounds; threshold for a 1e-3 false-veto rate on %ld burst-free shots; T1 impacts, peak 1 us, "
            "in 5%% of %ld shots; discarded fraction at chip-wide severity s over %ld shots)\n", nc, ni, ns);
        Json wj = Json::array();
        for (auto [W, rad] : std::vector<std::pair<int, double>>{{1, 1}, {2, 1}, {3, 2}, {3, 10}, {11, 10}}) {
            VetoSpec V; V.window_rounds = W; V.radius = rad;
            const auto cal = calibrate_veto(M, nc, 3, 1e-3, V);
            V.on = true; V.threshold = cal.threshold;
            FrameRunOptions o; o.impact.p_shot = 0.05; o.veto = V;
            const auto r = run_memory_frame(M, ni, 9, o);
            const double caught = double(r.impacts_vetoed) / double(std::max(r.impacts, 1L));
            say("   %2d rounds x radius %4.1f: threshold %2d (false %.1e), caught %.1f%% of %ld impacts;", W, rad, V.threshold, cal.false_rate, 100 * caught, r.impacts);
            Json j = Json::object(); j.set("window_rounds", W); j.set("radius", rad); j.set("threshold", V.threshold); j.set("caught", caught);
            BurstSpec B;
            for (double sv : {0.1, 0.03, 0.01}) {
                const auto Ms = prepare_memory(D, R, burst_noise(nm, B, sv), 10, 'Z', false);
                FrameRunOptions so; so.decoder = M.decoder.get(); so.veto = V;
                const auto rs = run_memory_frame(Ms, ns, 5, so);
                const double v = double(rs.discarded) / double(rs.shots);
                say(" v(%.2f) %.3f", sv, v);
                j.set("veto_at_" + std::to_string(sv).substr(0, 4), v);
            }
            say("\n");
            wj.push(j);
        }
        out.set("veto_windows", wj);
    }

    // ------------------------------------------------------------------ 3. latched hold
    {
        const long n = 20000 / div;
        VetoSpec V;
        V.threshold = calibrate_veto(M, 50000 / div, 3, 1e-3, V).threshold;
        BurstSpec B;   // T1 burst: peak excess T1 1 us, tau 25 ms, 0.1 per s
        const auto bv = burst_floor_veto(D, R, nm, B, 10, n, 505, V);
        say("\n3. Latched veto hold (d = 5 test chip, 10 rounds, T1 burst tau 25 ms, %ld shots per severity, threshold %d)\n", n, V.threshold);
        say("   floor %.3e per round without the veto, %.3e with shot-by-shot discard (%.2f ms discarded per impact)\n",
            bv.floor_per_round, bv.floor_per_round_veto, bv.discarded_ms);
        say("   ladder:");
        for (const auto& p : bv.curve) say(" s %.1e eps %.5f veto %.4f;", p.s, p.eps, p.veto_frac);
        say("\n");
        Json lj = Json::array();
        for (const auto& L : bv.latched) {
            const bool below = std::exp(-L.hold_tau) < bv.curve.back().s;
            say("   hold %2.0f tau: rounds lost %.3e, floor %.3e per round, dead time %.2e%s\n", L.hold_tau, L.rounds_lost, L.floor_per_round, L.dead_time_frac,
                below ? " (below the ladder: linear-tail extrapolation)" : "");
            Json j = Json::object(); j.set("hold_tau", L.hold_tau); j.set("floor_per_round", L.floor_per_round); j.set("dead_time_frac", L.dead_time_frac);
            j.set("below_ladder", below); lj.push(j);
        }
        out.set("latched_hold", lj);
        out.set("floor_no_veto", bv.floor_per_round); out.set("floor_shot_discard", bv.floor_per_round_veto);
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_all).count();
    say("\nsweeps complete in %.0f s\n", secs);
    out.set("runtime_s", secs);
    out.set("provenance", provenance(g_opt, "heralds: 77; veto windows: 3 (calibration), 9 (injection), 5 (severity runs); hold: 3 (calibration), 505 (ladder)"));
    std::ofstream f(g_opt.path("results/decoder_sweeps.json")); out.dump(f); f << "\n";
    std::ofstream t(g_opt.path("results/decoder_sweeps.txt")); t << g_txt;
    return 0;
}
