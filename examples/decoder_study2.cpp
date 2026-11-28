// Decoder study, part 2, for Phase Gamma (Gamma-1, GF-97) -- the decoding items that followed decoder_study:
//   1. herald weights derived from the leakage model: logical error with no heralds, with the fixed-factor footprint
//      rule and with the derived weights, at d = 3, 5, 7 and at Gamma-1's CZ leakage and 10x, for exact matching and
//      union-find, on the same syndromes (docs/design/derived-herald-weights.md)
//   2. data-qubit leakage heralding: DQLR that leaves the leakage on the measure qubit, whose next 3-state readout
//      heralds it, against the leakage-removing reset alone and against DQLR that dumps the leakage; on Gamma-1 as
//      designed and with the data qubits as the leaking partner of every CZ (docs/design/data-leakage-herald.md)
//   3. operation counts of the streaming union-find decoder at d = 3, 5, 7 for the window / commit settings of
//      decoder_study section 5, and the cycle counts of two hardware models against the QEC round, with and without an
//      isolated-pair pre-decoder (docs/design/streaming-union-find-fpga.md). A model: nothing here ran on hardware.
// Writes results/decoder_study2.{txt,json}. Options as for the reviews ([shots] [--out DIR] [--cache FILE] [--smoke]).
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
static double true_heralds_model(const MemorySetup& M, double fidelity) { double s = 0; for (double p : M.herald_model.prior) s += p; return s * fidelity; }
static const char* method_name(DecodeMethod m) { return m == DecodeMethod::Auto ? "exact" : "union-find"; }
static const char* method_key(DecodeMethod m) { return m == DecodeMethod::Auto ? "exact" : "union_find"; }
static Json run_json(const MemoryResult& r) {
    Json o = Json::object();
    o.set("failures", double(r.logical_errors)); o.set("eps", r.eps_per_round);
    return o;
}

int main(int argc, char** argv) {
    g_opt = parse(argc, argv, 100000);
    const long shots = g_opt.shots;
    const auto t_all = std::chrono::steady_clock::now();
    Json out = Json::object();
    const int R = 14;                                              // rounds of the memory runs (as in decoder_study)
    const std::vector<int> dists = g_opt.smoke ? std::vector<int>{3, 5} : std::vector<int>{3, 5, 7};
    say("Decoder study, part 2, Gamma-1 (GF-97): %ld shots per point, %d-round memory-Z runs unless stated\n", shots, R);

    // ------------------------------------------------------------------ 1. derived herald weights
    say("\n1. Herald weights derived from the leakage model (3-state readout of the measure qubits, fidelity 0.95, false-herald rate 1e-3;\n"
        "   leakage-removing reset; same syndromes; fixed = footprint rule, factor 0.4)\n");
    {
        Json hj = Json::object();
        for (int d : dists) {
            const ChipDesign D = gamma1(d);
            const auto Rd = prepared_report(D);
            for (double scale : {1.0, 10.0}) {
                auto nm = noise_model(D, Rd);
                for (auto& kv : nm.leak_cz) kv.second.second *= scale;
                auto M = prepare_memory(D, Rd, nm, R, 'Z');
                const auto t0 = std::chrono::steady_clock::now();
                prepare_herald_model(M);
                const double t_model = since(t0);
                const auto& H = M.herald_model;
                const long n = d == 7 ? 2 * shots : shots;
                const std::uint64_t seed = 300 + std::uint64_t(d);
                Json o = Json::object();
                for (DecodeMethod meth : {DecodeMethod::Auto, DecodeMethod::UnionFind}) {
                    FrameRunOptions fx; fx.method = meth; fx.herald.on = true;
                    FrameRunOptions off = fx; off.herald.use_in_decoder = false;
                    FrameRunOptions dv = fx; dv.herald.derived = true;
                    const auto rb = run_memory_frame(M, n, seed, off), rf = run_memory_frame(M, n, seed, fx), rd = run_memory_frame(M, n, seed, dv);
                    if (meth == DecodeMethod::Auto) {
                        std::size_t entries = 0; for (const auto& v : H.cond) entries += v.size();
                        say("   d=%d, CZ leakage x%-2.0f: %.3f leakage events and %.2f heralds per shot (%.0f%% leaked); true heralds per shot %.4f sampled, %.4f from the model;\n"
                            "       model: %zu (bit, edge) entries, at most %.3f of a herald unexplained by the %d-round look-back, built in %.2f s\n",
                            d, scale, rb.leak_events_per_shot, rb.heralds_per_shot, 100 * rb.herald_true_frac, rb.heralds_per_shot * rb.herald_true_frac,
                            true_heralds_model(M, fx.herald.fidelity), entries, H.worst_unexplained, H.max_rounds, t_model);
                        o.set("shots", double(n)); o.set("leak_events_per_shot", rb.leak_events_per_shot); o.set("heralds_per_shot", rb.heralds_per_shot);
                        o.set("herald_true_frac", rb.herald_true_frac); o.set("true_heralds_per_shot_model", true_heralds_model(M, fx.herald.fidelity));
                        o.set("model_entries", double(entries)); o.set("model_worst_unexplained", H.worst_unexplained);
                        o.set("model_dropped_mass", H.dropped_mass); o.set("model_pruned_mass", H.pruned_mass); o.set("model_kept_mass", H.kept_mass);
                    }
                    auto pct = [&](const MemoryResult& r) { return rb.logical_errors > 0 ? 100.0 * (double(r.logical_errors) / double(rb.logical_errors) - 1) : 0.0; };
                    say("       %-10s: failures of %ld: %ld without heralds, %ld fixed (%+.0f%%), %ld derived (%+.0f%%); eps %.3e / %.3e / %.3e\n", method_name(meth), n,
                        rb.logical_errors, rf.logical_errors, pct(rf), rd.logical_errors, pct(rd), rb.eps_per_round, rf.eps_per_round, rd.eps_per_round);
                    Json m = Json::object();
                    m.set("none", run_json(rb)); m.set("fixed", run_json(rf)); m.set("derived", run_json(rd));
                    o.set(method_key(meth), m);
                }
                hj.set("d" + std::to_string(d) + "_x" + std::to_string(int(scale)), o);
            }
        }
        out.set("derived_herald_weights", hj);
    }

    // ------------------------------------------------------------------ 2. data-qubit leakage heralding
    say("\n2. Data-qubit leakage heralding (derived weights; same syndromes within a scheme, different syndromes between schemes)\n"
        "   schemes: LRR = leakage-removing reset alone (Gamma-1); +DQLR = data-qubit leakage removal that dumps the leakage; +DQLR-herald = the\n"
        "   removal pulse leaves the leakage on the measure qubit, heralded by that qubit's next readout\n");
    {
        Json dj = Json::object();
        for (int d : dists) {
            if (d == 3) continue;
            const long n = shots;
            const std::uint64_t seed = 320 + std::uint64_t(d);
            const std::vector<DecodeMethod> meths = d == 7 ? std::vector<DecodeMethod>{DecodeMethod::Auto} : std::vector<DecodeMethod>{DecodeMethod::Auto, DecodeMethod::UnionFind};
            for (int variant : {0, 1}) {
                say("   d=%d, %s\n", d, variant == 0 ? "Gamma-1 as designed (the measure qubit is the leaking partner of every CZ; data qubits leak by transport only)"
                                                     : "stress: the data qubit is the leaking partner of every CZ");
                for (double scale : {1.0, 10.0}) {
                    Json sj = Json::object();
                    for (int scheme : {0, 1, 2}) {
                        ChipDesign D = gamma1(d);
                        if (scheme >= 1) D.leakage.dqlr = true;                // the removal pulse lengthens the reset (and the round)
                        if (scheme == 2) D.leakage.dqlr_transfer = true;
                        const auto Rd = prepared_report(D);
                        auto nm = noise_model(D, Rd);
                        for (auto& kv : nm.leak_cz) {
                            kv.second.second *= scale;
                            if (variant == 1) { const int a = kv.first.first, b = kv.first.second; kv.second.first = D.qubits[std::size_t(a)].role == "data" ? a : b; }
                        }
                        auto M = prepare_memory(D, Rd, nm, R, 'Z');
                        prepare_herald_model(M);
                        const char* sname = scheme == 0 ? "LRR" : scheme == 1 ? "+DQLR" : "+DQLR-herald";
                        Json o = Json::object();
                        for (DecodeMethod meth : meths) {
                            FrameRunOptions dv; dv.method = meth; dv.herald.on = true; dv.herald.derived = true;
                            FrameRunOptions off = dv; off.herald.use_in_decoder = false;
                            const auto rb = run_memory_frame(M, n, seed, off), rd = run_memory_frame(M, n, seed, dv);
                            if (meth == meths.front()) {
                                say("     x%-2.0f %-12s: round %.0f ns, leaked data %.1e meas %.1e, %.2f heralds per shot (%.0f%% leaked; true %.4f sampled, %.4f model)\n", scale, sname,
                                    rb.round_ns, rb.leak_frac_data, rb.leak_frac_meas, rb.heralds_per_shot, 100 * rb.herald_true_frac, rb.heralds_per_shot * rb.herald_true_frac,
                                    true_heralds_model(M, dv.herald.fidelity));
                                o.set("round_ns", rb.round_ns); o.set("leak_frac_data", rb.leak_frac_data); o.set("leak_frac_meas", rb.leak_frac_meas);
                                o.set("leak_events_per_shot", rb.leak_events_per_shot);
                                o.set("heralds_per_shot", rb.heralds_per_shot); o.set("herald_true_frac", rb.herald_true_frac);
                                o.set("true_heralds_per_shot_model", true_heralds_model(M, dv.herald.fidelity)); o.set("model_worst_unexplained", M.herald_model.worst_unexplained);
                            }
                            say("         %-10s: failures of %ld: %ld without heralds (eps %.3e), %ld with (eps %.3e): %+.0f%%\n", method_name(meth), n, rb.logical_errors,
                                rb.eps_per_round, rd.logical_errors, rd.eps_per_round, rb.logical_errors > 0 ? 100.0 * (double(rd.logical_errors) / double(rb.logical_errors) - 1) : 0.0);
                            Json m = Json::object();
                            m.set("none", run_json(rb)); m.set("derived", run_json(rd));
                            o.set(method_key(meth), m);
                        }
                        sj.set(scheme == 0 ? "lrr" : scheme == 1 ? "lrr_dqlr" : "lrr_dqlr_herald", o);
                    }
                    dj.set("d" + std::to_string(d) + (variant == 0 ? "_designed" : "_data_leak") + "_x" + std::to_string(int(scale)), sj);
                }
            }
        }
        dj.set("shots", double(shots));
        out.set("data_leakage_herald", dj);
    }

    // ------------------------------------------------------------------ 3. operation counts and the hardware models
    {
        const int RS = g_opt.smoke ? 12 : 50;
        const long n = g_opt.smoke ? shots : std::max(10000L, shots / 2);
        const std::vector<UFCostModel> models = {UFCostModel::serial_bram(), UFCostModel::parallel_growth()};
        say("\n3. Operation counts of the streaming union-find decoder, %d rounds, memory-Z, %ld shots (seed 601: the syndromes of decoder_study section 5)\n", RS, n);
        say("   per window: mean / 99th percentile / 99.9th percentile / maximum. Hardware models (stated in docs/design/streaming-union-find-fpga.md), %.0f MHz:\n",
            models[0].clock_mhz);
        say("     S  serial, block RAM: one table access per cycle, 2-cycle reads, pipeline depth 4\n");
        say("     P  parallel growth: one cell per vertex and per edge, 5 cycles per growth step, 2 per union, serial peeling\n");
        say("     +pre: an isolated-pair pre-decoder runs first (2 cycles when it clears the window)\n");
        Json pj = Json::object();
        for (int d : dists) {
            const ChipDesign D = gamma1(d);
            const auto Rd = prepared_report(D);
            const auto nm = noise_model(D, Rd);
            const auto M = prepare_memory(D, Rd, nm, RS, 'Z');
            const double round_ns = M.prog.duration_ns / double(RS), budget = round_ns * models[0].clock_mhz * 1e-3;
            say("   d=%d: round %.1f ns = %.1f cycles\n", d, round_ns, budget);
            Json dj = Json::array();
            for (auto [W, C] : std::vector<std::pair<int, int>>{{2, 1}, {4, 2}, {6, 3}, {10, 5}, {20, 10}}) {
                if (W > RS) continue;
                const StreamingUnionFind S(int(M.exp.detectors.size()), M.dem.edges, M.exp.detector_round, W, C);
                const auto P = run_stream_profile(M, S, n, 601, false, models), Q = run_stream_profile(M, S, n, 601, true, models);
                const long with_def = Q.windows - Q.windows_empty;
                say("     window %2d / commit %2d (%d vertices, %d edges per window; %ld windows): %ld failures, %ld with the pre-decoder; windows without defects %.1f%%;\n"
                    "         the pre-decoder clears %.1f%% of all windows (%.1f%% of those with defects)\n",
                    W, C, S.window_vertices(S.n_windows() > 2 ? 1 : 0), S.window_edges(S.n_windows() > 2 ? 1 : 0), P.windows, P.logical_errors, Q.logical_errors,
                    100.0 * double(P.windows_empty) / double(P.windows), 100.0 * double(Q.windows_cleared) / double(Q.windows),
                    with_def > 0 ? 100.0 * double(Q.windows_cleared) / double(with_def) : 0.0);
                Json o = Json::object();
                o.set("window", W); o.set("commit", C); o.set("windows", double(P.windows)); o.set("failures", double(P.logical_errors));
                o.set("failures_predecode", double(Q.logical_errors)); o.set("frac_windows_empty", double(P.windows_empty) / double(P.windows));
                o.set("frac_windows_cleared_by_predecoder", double(Q.windows_cleared) / double(Q.windows));
                o.set("frac_defect_windows_cleared_by_predecoder", with_def > 0 ? double(Q.windows_cleared) / double(with_def) : 0.0);
                o.set("window_vertices", S.window_vertices(S.n_windows() > 2 ? 1 : 0)); o.set("window_edges", S.window_edges(S.n_windows() > 2 ? 1 : 0));
                auto stats = [](const StreamProfile::Hist& h, double scale) {
                    Json s = Json::object();
                    s.set("mean", StreamProfile::mean(h) * scale); s.set("p99", double(StreamProfile::quantile(h, 0.99)) * scale);
                    s.set("p999", double(StreamProfile::quantile(h, 0.999)) * scale); s.set("max", double(StreamProfile::max(h)) * scale);
                    return s;
                };
                Json cj = Json::object();
                say("         counts per window:");
                for (int k = 0; k < UFOpCounts::kFields - 1; ++k) {        // pre_pairs is zero without the pre-decoder
                    const auto& h = P.count_hist[std::size_t(k)];
                    say("%s %s %.2f / %ld / %ld / %ld;", k == 5 ? "\n                           " : "", UFOpCounts::name(k), StreamProfile::mean(h),
                        StreamProfile::quantile(h, 0.99), StreamProfile::quantile(h, 0.999), StreamProfile::max(h));
                    cj.set(UFOpCounts::name(k), stats(h, 1.0));
                }
                say("\n");
                o.set("counts_per_window", cj);
                Json mj = Json::object();
                for (int pre = 0; pre < 2; ++pre)
                    for (std::size_t mi = 0; mi < models.size(); ++mi) {
                        const StreamProfile& X = pre ? Q : P;
                        const auto& rh = X.round_hist[mi]; const auto& lh = X.latency_hist[mi];
                        const double mean = StreamProfile::mean(rh) / 1e3;
                        const double over = StreamProfile::above(rh, std::lround(1e3 * budget));
                        say("         %s%-4s cycles per committed round: %.1f / %.1f / %.1f / %.1f (mean = %.0f%% of the round; %.3f%% of windows over their budget);"
                            " latency in rounds: mean %.2f, 99.9th %.1f, max %.1f\n", mi == 0 ? "S" : "P", pre ? "+pre" : "", mean,
                            double(StreamProfile::quantile(rh, 0.99)) / 1e3, double(StreamProfile::quantile(rh, 0.999)) / 1e3, double(StreamProfile::max(rh)) / 1e3,
                            100 * mean / budget, 100 * over, StreamProfile::mean(lh) / budget, double(StreamProfile::quantile(lh, 0.999)) / budget,
                            double(StreamProfile::max(lh)) / budget);
                        Json m = Json::object();
                        m.set("cycles_per_committed_round", stats(rh, 1e-3)); m.set("cycles_per_window", stats(X.cycles_hist[mi], 1.0));
                        m.set("mean_fraction_of_round", mean / budget); m.set("frac_windows_over_budget", over);
                        m.set("latency_rounds", stats(lh, 1.0 / budget));
                        mj.set(std::string(mi == 0 ? "serial" : "parallel") + (pre ? "_predecode" : ""), m);
                    }
                o.set("models", mj);
                dj.push(o);
            }
            Json o = Json::object();
            o.set("round_ns", round_ns); o.set("cycles_per_round", budget); o.set("settings", dj);
            pj.set("d" + std::to_string(d), o);
        }
        pj.set("rounds", RS); pj.set("shots", double(n)); pj.set("clock_mhz", models[0].clock_mhz);
        Json cm = Json::object();
        for (const auto& m : models) {
            Json o = Json::object();
            o.set("empty", double(m.empty)); o.set("cleared", double(m.cleared)); o.set("fixed", double(m.fixed));
            for (int k = 0; k < UFOpCounts::kFields; ++k) o.set(std::string("per_") + UFOpCounts::name(k), double(m.per[k]));
            cm.set(m.name, o);
        }
        pj.set("cost_models", cm);
        out.set("stream_operation_model", pj);
    }

    const double secs = since(t_all);
    say("\nstudy complete in %.0f s\n", secs);
    out.set("shots", double(shots)); out.set("runtime_s", secs);
    out.set("provenance", provenance(g_opt, "derived herald weights: 300 + d (as the heralds of decoder_study); data-qubit heralding: 320 + d; operation counts: 601 "
                                            "(as the streaming section of decoder_study)"));
    {
        std::ofstream f(g_opt.path("results/decoder_study2.json")); out.dump(f); f << "\n";
        std::ofstream t(g_opt.path("results/decoder_study2.txt")); t << g_txt;
    }
    return 0;
}
