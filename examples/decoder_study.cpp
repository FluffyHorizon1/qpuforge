// Decoder study for Phase Gamma (Gamma-1, GF-97) — the decoding items of ROADMAP "Phase Gamma":
//   1. linear-time fault enumeration: setup time at d = 7 for 25 and 50 rounds, sparse (quadratic) against
//      linear-time, signatures compared fault by fault
//   2. union-find decoder against exact matching: logical error per round at d = 3, 5, 7 in both bases on the same
//      syndromes, and the software decode time
//   3. leakage-aware decoding from heralded leakage (3-state readout of the measure qubits, leakage-removing reset):
//      Gamma-1's own leakage and a 10x stressed case, decoded with and without the heralds
//   4. impact veto: per-shot impact injection at d = 7, the veto's false-veto rate, catch rate and discard fraction,
//      logical error with and without it; the floor over the burst profile with the veto, against burst_floor()
// Writes results/decoder_study.{txt,json}. Options as for the reviews ([shots] [--out DIR] [--cache FILE] [--smoke]).
#include "review_util.hpp"
#include <chrono>
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
static double since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

int main(int argc, char** argv) {
    g_opt = parse(argc, argv, 100000);
    const long shots = g_opt.shots;
    const auto t_all = std::chrono::steady_clock::now();
    Json out = Json::object();
    const ChipDesign G7 = gamma1(7), G5 = gamma1(5), G3 = gamma1(3);
    auto R7 = prepared_report(G7);
    const auto nm7 = noise_model(G7, R7);
    const int R = 14;                                              // rounds of the memory runs (as in the Gamma-1 review)
    auto at_least = [&](long n, long scaled) { return g_opt.smoke ? shots : std::max(n, scaled); };
    say("Decoder study, Gamma-1 (GF-97): %ld shots per point, %d-round memories unless stated\n", shots, R);

    // ------------------------------------------------------------------ 1. fault enumeration
    say("\n1. Single-fault enumeration at d = 7 (setup of the detector error model)\n");
    {
        Json ej = Json::array();
        const std::vector<int> rr = g_opt.smoke ? std::vector<int>{6, 12} : std::vector<int>{25, 50};
        for (int r : rr) {
            const MemorySetup M = prepare_memory(G7, R7, nm7, r, 'Z', false);
            auto t0 = std::chrono::steady_clock::now();
            const auto a = enumerate_faults_direct(M.frame, M.exp);
            const double t_direct = since(t0);
            t0 = std::chrono::steady_clock::now();
            const auto b = enumerate_faults(M.frame, M.exp);
            const double t_linear = since(t0);
            long differ = a.size() == b.size() ? 0 : long(std::max(a.size(), b.size()));
            for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
                if (a[i].dets != b[i].dets || a[i].obs != b[i].obs || a[i].p != b[i].p || a[i].moment != b[i].moment || a[i].what != b[i].what) ++differ;
            t0 = std::chrono::steady_clock::now();
            const MemorySetup Mf = prepare_memory(G7, R7, nm7, r, 'Z', true);
            const double t_setup = since(t0);
            say("   %2d rounds: %zu faults, %zu time slices; sparse %.2f s, linear-time %.2f s (%.1fx), %ld signatures differ; full setup %.2f s\n",
                r, b.size(), M.frame.moments.size(), t_direct, t_linear, t_direct / std::max(t_linear, 1e-9), differ, t_setup);
            Json o = Json::object();
            o.set("rounds", r); o.set("faults", double(b.size())); o.set("time_slices", double(M.frame.moments.size()));
            o.set("t_sparse_s", t_direct); o.set("t_linear_s", t_linear); o.set("signatures_differ", double(differ)); o.set("t_setup_s", t_setup);
            o.set("dem_edges", double(Mf.dem.edges.size()));
            ej.push(o);
        }
        out.set("fault_enumeration_d7", ej);
    }

    // ------------------------------------------------------------------ 2. union-find against exact matching
    say("\n2. Union-find decoder against exact matching (same syndromes; seeds as in the Gamma-1 review)\n");
    {
        Json uj = Json::object();
        for (int d : {3, 5, 7}) {
            const ChipDesign& D = d == 3 ? G3 : d == 5 ? G5 : G7;
            const auto Rd = prepared_report(D);
            const auto nm = noise_model(D, Rd);
            const long n = d == 7 ? 2 * shots : shots;
            Json dj = Json::object();
            for (char b : {'Z', 'X'}) {
                const auto M = prepare_memory(D, Rd, nm, R, b);
                const std::uint64_t seed = 11 + std::uint64_t(R) + (b == 'X' ? 7919u : 0u);
                FrameRunOptions uf; uf.method = DecodeMethod::UnionFind;
                const auto ex = run_memory_frame(M, n, seed), u = run_memory_frame(M, n, seed, uf);
                const double ratio = u.eps_per_round / std::max(ex.eps_per_round, 1e-30);
                say("   d=%d memory-%c: exact eps %.3e (%ld failures / %ld), union-find eps %.3e (%ld), ratio %.2f; decode %.2f vs %.2f us per round (%.1f defects per shot)\n",
                    d, b, ex.eps_per_round, ex.logical_errors, n, u.eps_per_round, u.logical_errors, ratio, ex.decode_us_per_shot / R, u.decode_us_per_shot / R, ex.mean_defects);
                Json o = Json::object();
                o.set("shots", double(n)); o.set("eps_exact", ex.eps_per_round); o.set("eps_uf", u.eps_per_round);
                o.set("failures_exact", double(ex.logical_errors)); o.set("failures_uf", double(u.logical_errors)); o.set("ratio", ratio);
                o.set("decode_us_per_round_exact", ex.decode_us_per_shot / R); o.set("decode_us_per_round_uf", u.decode_us_per_shot / R);
                o.set("mean_defects", ex.mean_defects);
                dj.set(std::string(1, b), o);
            }
            uj.set("d" + std::to_string(d), dj);
        }
        out.set("union_find", uj);
    }

    // ------------------------------------------------------------------ 3. heralded leakage
    say("\n3. Heralded leakage (3-state readout of the measure qubits; leakage-removing reset; footprint rule, factor 0.4)\n");
    {
        Json hj = Json::object();
        for (int d : {5, 7}) {
            const ChipDesign& D = d == 5 ? G5 : G7;
            const auto Rd = prepared_report(D);
            for (double scale : {1.0, 10.0}) {
                auto nm = noise_model(D, Rd);
                for (auto& kv : nm.leak_cz) kv.second.second *= scale;
                const auto M = prepare_memory(D, Rd, nm, R, 'Z');
                const long n = d == 7 ? 2 * shots : shots;
                Json o = Json::object();
                for (DecodeMethod meth : {DecodeMethod::Auto, DecodeMethod::UnionFind}) {
                    FrameRunOptions on; on.method = meth; on.herald.on = true;
                    FrameRunOptions off = on; off.herald.use_in_decoder = false;
                    const std::uint64_t seed = 300 + std::uint64_t(d);
                    const auto rb = run_memory_frame(M, n, seed, off), rh = run_memory_frame(M, n, seed, on);
                    const char* nmth = meth == DecodeMethod::Auto ? "exact" : "union-find";
                    say("   d=%d, CZ leakage x%-2.0f (%.2f events, %.2f heralds per shot, %.0f%% leaked), %-10s: eps %.3e without heralds (%ld failures), %.3e with (%ld): %+.0f%%\n",
                        d, scale, rh.leak_events_per_shot, rh.heralds_per_shot, 100 * rh.herald_true_frac, nmth, rb.eps_per_round, rb.logical_errors,
                        rh.eps_per_round, rh.logical_errors, 100 * (rh.eps_per_round / std::max(rb.eps_per_round, 1e-30) - 1));
                    Json m = Json::object();
                    m.set("eps_without", rb.eps_per_round); m.set("eps_with", rh.eps_per_round);
                    m.set("failures_without", double(rb.logical_errors)); m.set("failures_with", double(rh.logical_errors));
                    o.set(meth == DecodeMethod::Auto ? "exact" : "union_find", m);
                    o.set("leak_events_per_shot", rh.leak_events_per_shot); o.set("heralds_per_shot", rh.heralds_per_shot); o.set("herald_true_frac", rh.herald_true_frac);
                    o.set("shots", double(n));
                }
                hj.set("d" + std::to_string(d) + "_x" + std::to_string(int(scale)), o);
            }
        }
        out.set("heralds", hj);
    }

    // ------------------------------------------------------------------ 4. impact veto
    say("\n4. Impact veto at d = 7 (window 3 rounds x radius 2 lattice units, threshold for a 1e-3 false-veto rate)\n");
    {
        Json vj = Json::object();
        struct Sc { const char* key; const char* name; BurstSpec B; char basis; };
        BurstSpec raw; raw.rate_hz = 0.1; raw.tau_ms = 25; raw.t1_peak_us = 1.0; raw.tphi_peak_us = 0;
        BurstSpec gap; gap.rate_hz = 0.1; gap.tau_ms = 1; gap.t1_peak_us = 0; gap.tphi_peak_us = 0.1;
        for (const Sc& sc : {Sc{"unmitigated", "no mitigation: T1 burst, 25 ms recovery (memory-Z)", raw, 'Z'},
                             Sc{"gap_engineered", "gap-engineered junctions: dephasing burst, 1 ms (memory-X)", gap, 'X'}}) {
            const auto M = prepare_memory(G7, R7, nm7, R, sc.basis);
            VetoSpec V;
            const auto cal = calibrate_veto(M, shots, 401, 1e-3, V);
            V.on = true; V.threshold = cal.threshold;
            FrameRunOptions clean; clean.veto = V;
            FrameRunOptions inj; inj.impact.burst = sc.B; inj.impact.p_shot = 0.05;
            FrameRunOptions vet = inj; vet.veto = V;
            const auto rc = run_memory_frame(M, shots, 403, clean), ri = run_memory_frame(M, shots, 403, inj), rv = run_memory_frame(M, shots, 403, vet);
            const double caught = double(rv.impacts_vetoed) / double(std::max(rv.impacts, 1L));
            say("   %s\n", sc.name);
            say("     threshold %d events: false-veto rate %.1e on %ld burst-free calibration shots, %.1e on fresh ones\n",
                V.threshold, cal.false_rate, cal.shots, double(rc.discarded) / double(rc.shots));
            say("     impacts in 5%% of shots (chip-wide, onset uniform in the shot): %ld impacts, %.0f%% caught; %.2f%% of all shots discarded\n",
                rv.impacts, 100 * caught, 100.0 * double(rv.discarded) / double(rv.shots));
            say("     eps per round: %.3e burst-free, %.3e with impacts and no veto, %.3e with the veto; failures of impacted shots %ld -> %ld\n",
                rc.eps_kept, ri.eps_per_round, rv.eps_kept, ri.errors_impacted, rv.errors_impacted_kept);
            const long nl = at_least(4000, shots / 5);
            const auto bf = burst_floor(G7, R7, nm7, sc.B, R, nl, 505, sc.basis);
            const auto bv = burst_floor_veto(G7, R7, nm7, sc.B, R, nl, 505, V, sc.basis);
            say("     over the burst profile (%ld shots per severity): burst_floor() %.2e per round; ladder with the veto: %.2e without it, %.2e with shot-by-shot discard (%.1f ms discarded per impact)\n",
                nl, bf.floor_per_round, bv.floor_per_round, bv.floor_per_round_veto, bv.discarded_ms);
            say("     latched veto, hold after the onset:");
            for (const auto& L : bv.latched) say(" %.0f tau %.1e (dead %.1e)", L.hold_tau, L.floor_per_round, L.dead_time_frac);
            say("\n     ladder:");
            for (const auto& p : bv.curve) say(" s %.1e eps %.2e veto %.2f;", p.s, p.eps, p.veto_frac);
            say("\n");
            Json o = Json::object();
            o.set("threshold", V.threshold); o.set("false_veto_calibration", cal.false_rate); o.set("false_veto_fresh", double(rc.discarded) / double(rc.shots));
            o.set("impacts", double(rv.impacts)); o.set("caught", caught); o.set("discard_fraction", double(rv.discarded) / double(rv.shots));
            o.set("eps_burst_free", rc.eps_kept); o.set("eps_impacts_no_veto", ri.eps_per_round); o.set("eps_impacts_veto", rv.eps_kept);
            o.set("failures_impacted_no_veto", double(ri.errors_impacted)); o.set("failures_impacted_kept", double(rv.errors_impacted_kept));
            o.set("burst_floor", bf.floor_per_round); o.set("floor_ladder_no_veto", bv.floor_per_round); o.set("floor_shot_discard", bv.floor_per_round_veto);
            o.set("discarded_ms_per_impact", bv.discarded_ms); o.set("ladder_shots", double(nl));
            Json lj = Json::array();
            for (const auto& L : bv.latched) { Json l = Json::object(); l.set("hold_tau", L.hold_tau); l.set("floor_per_round", L.floor_per_round); l.set("dead_time_frac", L.dead_time_frac); l.set("rounds_lost", L.rounds_lost); lj.push(l); }
            o.set("latched", lj);
            Json cj = Json::array();
            for (const auto& p : bv.curve) { Json c = Json::object(); c.set("severity", p.s); c.set("eps", p.eps); c.set("eps_kept", p.eps_kept); c.set("veto_frac", p.veto_frac); c.set("saturated", p.saturated); cj.push(c); }
            o.set("curve", cj);
            vj.set(sc.key, o);
        }
        out.set("impact_veto_d7", vj);
    }

    const double secs = since(t_all);
    say("\nstudy complete in %.0f s\n", secs);
    out.set("shots", double(shots)); out.set("runtime_s", secs);
    out.set("provenance", provenance(g_opt, "union-find: 11 + rounds (memory-X adds 7919, as the Gamma-1 review); heralds: 300 + d; veto: 401 (calibration), "
                                            "403 (injection), 505 (burst ladders)"));
    {
        std::ofstream f(g_opt.path("results/decoder_study.json")); out.dump(f); f << "\n";
        std::ofstream t(g_opt.path("results/decoder_study.txt")); t << g_txt;
    }
    return 0;
}
