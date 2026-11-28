// Flux-line distortion study: what the CZ gates of Alpha-1 v1.1 and Beta-1 need from their flux lines
// (fluxline.hpp, docs/design/flux-line-predistortion.md).
//
//   flux_distortion_study [--out DIR] [--cache FILE] [--smoke]     (review_io.hpp)
//
// Writes results/flux_distortion_study.txt and results/flux_distortion_study.json (or into DIR). No Monte Carlo: every
// number is a deterministic pulse simulation of a calibrated gate from the CZ cache.
//   1. static flux offsets on the calibrated gates
//   2. one exponential tail on both lines, gate on its own: programmed as calibrated, and with the amplitudes re-calibrated
//   3. the same tail in the memory cycle of each chip: gates calibrated alone, then run back to back, round after round;
//      per gate, and per measure qubit (the phase shifts of its gates add up between its two Hadamards)
//   4. predistortion of an example line with fit errors
//   5. AWG sample rate and bandwidth (with a local re-calibration), DAC resolution, timing skew between the two lines
//   6. the pair class that the calibration's leakage fallback handled, re-calibrated for 1 - F
#include "qpuforge/qpuforge.hpp"
#include "review_io.hpp"
#include <cstdarg>
#include <exception>
#include <fstream>
#include <map>

using namespace qpuforge;

static std::FILE* g_txt = nullptr;
static void out(const char* fmt, ...) {
    char buf[4096];
    va_list ap; va_start(ap, fmt); std::vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    std::fputs(buf, stdout); std::fflush(stdout);
    if (g_txt) std::fputs(buf, g_txt);
}

namespace {

struct PairClass { double fm{0}, ft{0}, curv{0}; CZSystem S; CZPulse P; };        // one calibrated frequency combination
struct Slot { int a{0}, b{0}, mover{0}, cls{0}, line{0}; double offset{0}; };       // one CZ of the cycle
struct Line { int mover{0}; std::vector<int> gates; };                              // the flux line of one measure qubit
struct Chip {
    std::string name; double allowance{0};
    double period{0};                  // round period, ns: the same for every line (the slowest line of the compiled schedule)
    double period_min{0};              // shortest line period of that schedule (its lines are not synchronised)
    std::vector<PairClass> cls; std::vector<Slot> gates; std::vector<Line> lines;
};

// An exception must not leave an OpenMP loop: the first one is kept and thrown again after the loop.
struct LoopGuard {
    std::exception_ptr first;
    template <class F> void run(F&& body) {
        try { body(); }
        catch (...) {
#pragma omp critical(flux_study_guard)
            if (!first) first = std::current_exception();
        }
    }
    void rethrow() { if (first) std::rethrow_exception(first); }
};

// Calibrated pair classes of a chip and the pulse pattern of its memory experiment: for every measure qubit (the mover
// of all its CZ gates) the start of each CZ within a round, taken from a round of the compiled schedule in which the
// pattern of every line has stopped changing (it changes during the first ten rounds). The compiler starts every gate
// as soon as its qubits are free, so its lines keep slightly different periods; a machine runs all of them on one
// clock, and the study gives every line the period of the slowest.
Chip load_chip(const ChipDesign& D, const std::string& cache, int max_classes) {
    Chip C; C.name = D.name; C.allowance = D.coupler.residual_ctrl_err;
    auto R = analyze(D);
    const auto cal = calibrate_chip_cz(D, R, false, cache);
    const double sig2 = flux_noise_variance(D.flux_noise_uPhi0);
    const int rounds = 30, ref = 24;
    CompileOptions co; co.layout = "trivial";
    const auto prog = compile(memory_experiment(rotated_surface_code(D.surface_d), rounds, 'Z').circuit, chip_graph(D, R), co);
    std::map<std::pair<int, int>, std::vector<double>> starts;
    for (const auto& s : cz_starts(prog)) {
        const int a = D.qubits[std::size_t(s.a)].id, b = D.qubits[std::size_t(s.b)].id;
        starts[{std::min(a, b), std::max(a, b)}].push_back(s.t_ns);
    }
    std::map<int, int> line_of;
    for (const auto& c : cal) {
        int ci = -1;
        for (std::size_t i = 0; i < C.cls.size(); ++i) if (C.cls[i].fm == c.f_mover && C.cls[i].ft == c.f_target) ci = int(i);
        if (ci < 0) {
            if (int(C.cls.size()) >= max_classes) continue;
            PairClass k; k.fm = c.f_mover; k.ft = c.f_target; k.P = c.pulse; k.curv = sig2 > 0 ? c.err_flux_quasistatic / (0.5 * sig2) : 0.0;
            k.S = cz_system_of(c, R.q[std::size_t(D.index_of(c.mover))].alpha, R.q[std::size_t(D.index_of(c.target))].alpha, D.coupler);
            C.cls.push_back(k); ci = int(C.cls.size()) - 1;
        }
        if (!line_of.count(c.mover)) { line_of[c.mover] = int(C.lines.size()); Line L; L.mover = c.mover; C.lines.push_back(L); }
        Slot g; g.a = std::min(c.a, c.b); g.b = std::max(c.a, c.b); g.mover = c.mover; g.cls = ci; g.line = line_of[c.mover];
        const auto& v = starts.at({g.a, g.b});
        if (int(v.size()) != rounds) throw std::runtime_error("flux_distortion_study: a pair does not have one CZ per round");
        g.offset = v[std::size_t(ref)];
        // settled: the pair repeats with one period over the rounds around the reference round
        const double per = v[std::size_t(ref) + 1] - v[std::size_t(ref)];
        for (int r = ref - 4; r + 1 < rounds; ++r)
            if (std::abs(v[std::size_t(r) + 1] - v[std::size_t(r)] - per) > 1e-6) throw std::runtime_error("flux_distortion_study: the schedule has not settled by the reference round");
        C.period = std::max(C.period, per);
        C.period_min = C.period_min > 0 ? std::min(C.period_min, per) : per;
        C.lines[std::size_t(g.line)].gates.push_back(int(C.gates.size()));
        C.gates.push_back(g);
    }
    if (C.gates.empty()) throw std::runtime_error("flux_distortion_study: no calibrated CZ gate");
    if (C.period - C.period_min > 0.02 * C.period) throw std::runtime_error("flux_distortion_study: the line periods of the schedule differ by more than 2 %");
    for (auto& L : C.lines) {
        std::sort(L.gates.begin(), L.gates.end(), [&](int x, int y) { return C.gates[std::size_t(x)].offset < C.gates[std::size_t(y)].offset; });
        const double t0 = C.gates[std::size_t(L.gates.front())].offset;
        for (int g : L.gates) C.gates[std::size_t(g)].offset -= t0;
    }
    return C;
}

// A gate on its own, after a long idle: programmed steps, virtual Z gates and result.
struct Alone { CZCommand cmd; CZFrame frame; double err{0}, leak{0}, cphase{0}; };
std::vector<Alone> alone(const Chip& C, const FluxSetup& F, bool recalibrate_amplitude) {
    std::vector<Alone> v(C.cls.size());
    LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < int(C.cls.size()); ++i) guard.run([&] {
        const auto& c = C.cls[std::size_t(i)];
        Alone A; A.cmd = recalibrate_amplitude ? cz_matched_command(c.S, c.P, F) : cz_nominal_command(c.S, c.P);
        const auto r = simulate_cz_drive(c.S, cz_gate_drive(c.S, c.P, F, A.cmd));
        A.frame = cz_frame(r); A.err = 1 - r.fidelity; A.leak = r.leakage; A.cphase = r.cphase;
        v[std::size_t(i)] = A;
    });
    guard.rethrow();
    return v;
}

// Every gate of the chip in round `round` of the cycle (0: the first round after a long idle), relative to the same
// gate on its own: with the virtual Z gates kept from the calibration alone (fixed) and re-calibrated in place (free).
// est_*: the quasi-static estimate of the same numbers, from the flux offsets averaged over the flat top and the
// gate's response to constant offsets (cz_offset_response), with no pulse simulation. offset: the mover line's flux
// averaged over the flat top, against the gate on its own, as a fraction of the mover's flux step.
struct InCycle { double fixed{0}, free{0}, dleak{0}, dphase{0}, est_fixed{0}, est_free{0}, offset{0}; };   // offset: mover flat top, fraction of its step
std::vector<InCycle> in_cycle(const Chip& C, const FluxSetup& F, const std::vector<Alone>& al, const std::vector<CZOffsetResponse>& resp, int round) {
    const double Ts = F.awg.period(), settle = F.settle_ns();
    std::vector<double> start(C.gates.size());
    double busy = 0;                   // the longest time a line needs for the pulses of one round
    for (std::size_t l = 0; l < C.lines.size(); ++l) {
        // pulses start on the AWG sample grid, and never before the previous pulse on the line has ended and settled
        double free_at = -1e300;
        for (int j : C.lines[l].gates) {
            const auto& g = C.gates[std::size_t(j)];
            double t = std::round(g.offset / Ts) * Ts;
            if (t < free_at) t = std::ceil(free_at / Ts - 1e-9) * Ts;
            start[std::size_t(j)] = t;
            free_at = t + C.cls[std::size_t(g.cls)].P.duration() + settle;
        }
        busy = std::max(busy, free_at);
    }
    const double T = std::ceil(std::max(C.period, busy) / Ts - 1e-9) * Ts;   // the round period, on the sample grid
    std::vector<InCycle> v(C.gates.size());
    LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
    for (int gi = 0; gi < int(C.gates.size()); ++gi) guard.run([&] {
        const auto& g = C.gates[std::size_t(gi)];
        const auto& c = C.cls[std::size_t(g.cls)];
        const auto& A = al[std::size_t(g.cls)];
        std::vector<FluxPulse> mover, coupler;
        for (int r = 0; r <= round; ++r) {
            for (int j : C.lines[std::size_t(g.line)].gates) {
                const auto& h = C.gates[std::size_t(j)];
                const auto& hc = C.cls[std::size_t(h.cls)];
                mover.push_back({r * T + start[std::size_t(j)], al[std::size_t(h.cls)].cmd.amp_mover, hc.P.t_rise, hc.P.t_flat});
            }
            coupler.push_back({r * T + start[std::size_t(gi)], A.cmd.amp_coupler, c.P.t_rise, c.P.t_flat});
        }
        const auto W = cz_train_drive(c.S, c.P, F, round * T + start[std::size_t(gi)], mover, coupler), Wa = cz_gate_drive(c.S, c.P, F, A.cmd);
        const auto r = simulate_cz_drive(c.S, W);
        InCycle x;
        double dm = 0, dc = 0, hold = 0, t = 0;
        for (std::size_t s = 0; s < W.h.size(); ++s) {
            const double mid = t + 0.5 * W.h[s];
            if (mid > c.P.t_rise && mid < c.P.t_rise + c.P.t_flat) { dm += (W.f1[s] - Wa.f1[s]) * W.h[s]; dc += (W.fc[s] - Wa.fc[s]) * W.h[s]; hold += W.h[s]; }
            t += W.h[s];
        }
        if (hold > 0) { dm /= hold; dc /= hold; }
        x.offset = dm / (c.P.flux1_int - c.S.flux1_idle);
        x.est_fixed = resp[std::size_t(g.cls)].added_fixed(dm, dc);
        x.est_free = resp[std::size_t(g.cls)].added_free(dm, dc);
        x.fixed = (1 - cz_fidelity_in_frame(r, A.frame)) - A.err;
        x.free = (1 - r.fidelity) - A.err;
        x.dleak = r.leakage - A.leak;
        x.dphase = cz_frame(r).th_mover - A.frame.th_mover;
        x.dphase = std::remainder(x.dphase, 2 * PI);
        v[std::size_t(gi)] = x;
    });
    guard.rethrow();
    return v;
}
// The measure qubit of a line is in a superposition between its two Hadamards, and the phase shifts of its gates add
// up there. With the virtual Z gates kept from the calibration alone, that sum is a Z rotation nobody corrects, and
// the qubit's outcome flips with probability sin^2(sum / 2). `independent` counts the same phases gate by gate, the
// way a sum of per-gate errors does. An estimate from the mover phases of the pair simulations: the data qubits'
// phases and the conditional-phase errors are not in it, and nothing here simulates more than one pair.
struct LinePhase { double sum{0}, flip{0}, independent{0}; int gates{0}; };
std::vector<LinePhase> line_phases(const Chip& C, const std::vector<InCycle>& v) {
    std::vector<LinePhase> out(C.lines.size());
    for (std::size_t l = 0; l < C.lines.size(); ++l) {
        LinePhase& p = out[l];
        for (int g : C.lines[l].gates) { const double d = v[std::size_t(g)].dphase; p.sum += d; p.independent += std::pow(std::sin(0.5 * d), 2); ++p.gates; }
        p.flip = std::pow(std::sin(0.5 * p.sum), 2);
    }
    return out;
}
// Probability that depolarising CZ gates with average infidelity r flip the outcome of a measure qubit with n gates:
// a two-qubit depolarising channel with average gate infidelity r applies each of the 15 Paulis with probability
// r / 12, and 8 of them have Z or Y on the measure qubit. To first order in r.
double depolarising_flip(double r, double gates) { return gates * r * 8.0 / 12.0; }
// Response of every class to constant flux offsets, on ideal lines.
std::vector<CZOffsetResponse> responses(const Chip& C) {
    std::vector<CZOffsetResponse> v(C.cls.size());
    LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < int(C.cls.size()); ++i)
        guard.run([&] { v[std::size_t(i)] = cz_offset_response(C.cls[std::size_t(i)].S, cz_ideal_drive(C.cls[std::size_t(i)].S, C.cls[std::size_t(i)].P)); });
    guard.rethrow();
    return v;
}
// A round late enough for the slowest tail to have reached its periodic state (exp(-14) of the start remains).
int steady_round(const Chip& C, double tau_max) { return std::min(600, int(std::ceil(14 * tau_max / C.period)) + 2); }

// Local three-knob calibration on a given setup, as on hardware: programmed mover step, flat-top length and programmed
// coupler step, minimising 1 - F from the amplitude-matched start.
struct Recal { CZCommand cmd; CZPulse P; double err{0}, leak{0}; };
Recal recalibrate(const PairClass& c, const FluxSetup& F, int max_eval) {
    const CZCommand c0 = cz_matched_command(c.S, c.P, F);
    auto sim = [&](const std::vector<double>& x) {
        CZPulse P = c.P; P.t_flat = std::max(0.0, x[1]);
        return simulate_cz_drive(c.S, cz_gate_drive(c.S, P, F, {x[0], x[2]}));
    };
    auto cost = [&](const std::vector<double>& x) { return 1 - sim(x).fidelity; };
    const auto best = nelder_mead(cost, {c0.amp_mover, c.P.t_flat, c0.amp_coupler}, {2e-4, 0.25, 5e-4}, max_eval, 1e-9);
    Recal R; R.cmd = {best[0], best[2]}; R.P = c.P; R.P.t_flat = std::max(0.0, best[1]);
    const auto r = sim(best);
    R.err = 1 - r.fidelity; R.leak = r.leakage;
    return R;
}

struct Stat { double mean{0}, max{-1e300}; };
template <class V, class F> Stat stat(const V& v, F f) {
    Stat s;
    for (const auto& x : v) { const double y = f(x); s.mean += y / double(v.size()); s.max = std::max(s.max, y); }
    return s;
}
double median(std::vector<double> v) {
    if (v.empty()) throw std::invalid_argument("median of no values");
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
}
Json jstat(const Stat& s) { Json o = Json::object(); o.set("mean", s.mean); o.set("max", s.max); return o; }
FluxSetup with_tail(double a, double tau) {
    FluxSetup F; F.pad_ns = 0;
    F.mover.line.sections = {{a, tau}}; F.coupler.line.sections = {{a, tau}};
    return F;
}

int run(int argc, char** argv) {
    const auto t_start = std::chrono::steady_clock::now();
    const auto opt = review::parse(argc, argv, 0);
    const std::string txt_path = opt.path("results/flux_distortion_study.txt"), json_path = opt.path("results/flux_distortion_study.json");
    g_txt = std::fopen(txt_path.c_str(), "w");
    if (!g_txt) { std::fprintf(stderr, "cannot write %s\n", txt_path.c_str()); return 1; }
    const int max_cls = opt.smoke ? 3 : 1000;
    const Chip A = load_chip(alpha1_v11(), opt.cache, max_cls);
    Json J = Json::object();
    J.set("schema", "qpuforge-flux-distortion-study");
    J.set("provenance", review::provenance(opt, "none (deterministic pulse simulations)"));

    out("qpuforge flux-line distortion study (fluxline.hpp): CZ gates of %s from the CZ cache, model version %d\n", A.name.c_str(), CZModel::kVersion);
    {
        double m0 = 1e300, m1 = 0, c0 = 1e300, c1 = 0, t0 = 1e300, t1 = 0;
        for (const auto& c : A.cls) {
            const double sm = c.P.flux1_int - c.S.flux1_idle, sc = c.P.fluxc_on - c.S.fluxc_idle;
            m0 = std::min(m0, sm); m1 = std::max(m1, sm); c0 = std::min(c0, sc); c1 = std::max(c1, sc);
            t0 = std::min(t0, c.P.duration()); t1 = std::max(t1, c.P.duration());
        }
        out("%zu pair classes, %.1f-%.1f ns; flux step of the mover (the measure qubit) %.3f-%.3f Phi0, of the coupler %.3f-%.3f Phi0\n", A.cls.size(), t0, t1, m0, m1, c0, c1);
    }
    out("added error = change of 1 - F (average gate fidelity) against the same gate on ideal lines, or against the gate on its own\n");
    out("  free: the virtual Z gates after the CZ are re-calibrated for the case shown;  fixed: they stay as calibrated\n\n");

    FluxSetup ideal; ideal.pad_ns = 0;
    const auto ref = alone(A, ideal, false);
    const auto respA = responses(A);
    std::vector<double> slope(A.cls.size());     // the mover's frequency per flux at the interaction point (dressed |10>, coupler on), GHz / Phi0
    for (std::size_t i = 0; i < A.cls.size(); ++i) {
        const auto& c = A.cls[i];
        auto f10 = [&](double f1) { const auto d = c.S.dressed(f1, c.P.fluxc_on, c.S.flux2_idle); return d.E[2] - d.E[0]; };
        slope[i] = (f10(c.P.flux1_int + 1e-4) - f10(c.P.flux1_int - 1e-4)) / 2e-4;
    }
    {
        Json jc = Json::array();
        for (std::size_t i = 0; i < A.cls.size(); ++i) {
            const auto& c = A.cls[i];
            const auto& R = respA[i];
            Json o = Json::object();
            o.set("f_mover_GHz", c.fm); o.set("f_target_GHz", c.ft); o.set("t_cz_ns", c.P.duration());
            o.set("step_mover_Phi0", c.P.flux1_int - c.S.flux1_idle); o.set("step_coupler_Phi0", c.P.fluxc_on - c.S.fluxc_idle);
            o.set("err_ideal_lines", ref[i].err); o.set("leakage_ideal_lines", ref[i].leak);
            o.set("mover_slope_GHz_per_Phi0", slope[i]);
            Json jr = Json::object();           // added 1 - F for constant offsets dm, dc (Phi0): m dm + c dc + mm dm^2 + cc dc^2 + mc dm dc
            jr.set("free_m", R.free_m); jr.set("free_c", R.free_c); jr.set("free_mm", R.free_mm); jr.set("free_cc", R.free_cc); jr.set("free_mc", R.free_mc);
            jr.set("fixed_m", R.fixed_m); jr.set("fixed_c", R.fixed_c); jr.set("fixed_mm", R.fixed_mm); jr.set("fixed_cc", R.fixed_cc); jr.set("fixed_mc", R.fixed_mc);
            jr.set("mover_phase_rad_per_Phi0", R.phase_mover);
            o.set("offset_response", jr);
            jc.push(o);
        }
        J.set("pair_classes", jc);
    }

    // ---------------------------------------------------------------- 1. static offsets
    out("1. Static flux offset during the gate, as a fraction of the line's flux step (%zu classes, median / worst)\n", A.cls.size());
    out("   offset   line       free             fixed\n");
    Json j1 = Json::array();
    std::vector<double> curv_ratio(A.cls.size(), 0.0);
    for (double eps : opt.smoke ? std::vector<double>{3e-4} : std::vector<double>{1e-4, 3e-4, 1e-3}) {
        for (int line = 0; line < 2; ++line) {
            std::vector<double> dfree(A.cls.size()), dfixed(A.cls.size());
            LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < int(A.cls.size()); ++i) guard.run([&] {
                const auto& c = A.cls[std::size_t(i)];
                double fr = 0, fx = 0;
                for (int sgn : {+1, -1}) {
                    CZCommand cmd = cz_nominal_command(c.S, c.P);
                    (line == 0 ? cmd.amp_mover : cmd.amp_coupler) *= 1 + sgn * eps;
                    const auto r = simulate_cz_drive(c.S, cz_gate_drive(c.S, c.P, ideal, cmd));
                    fr += 0.5 * ((1 - r.fidelity) - ref[std::size_t(i)].err);
                    fx += 0.5 * ((1 - cz_fidelity_in_frame(r, ref[std::size_t(i)].frame)) - ref[std::size_t(i)].err);
                }
                dfree[std::size_t(i)] = fr; dfixed[std::size_t(i)] = fx;
            });
            guard.rethrow();
            const Stat sf = stat(dfree, [](double x) { return x; }), sx = stat(dfixed, [](double x) { return x; });
            out("   %.0e    %-8s   %.2e / %.2e   %.2e / %.2e\n", eps, line == 0 ? "mover" : "coupler", median(dfree), sf.max, median(dfixed), sx.max);
            Json o = Json::object();
            o.set("offset", eps); o.set("line", line == 0 ? "mover" : "coupler");
            o.set("free_median", median(dfree)); o.set("free_worst", sf.max); o.set("fixed_median", median(dfixed)); o.set("fixed_worst", sx.max);
            j1.push(o);
            if (eps == (opt.smoke ? 3e-4 : 1e-4))            // curvature of 1 - F in flux, from the symmetric difference
                for (std::size_t i = 0; i < A.cls.size(); ++i) {
                    const auto& c = A.cls[i];
                    const double step = line == 0 ? c.P.flux1_int - c.S.flux1_idle : c.P.fluxc_on - c.S.fluxc_idle;
                    curv_ratio[i] += 2 * dfree[i] / (eps * step * eps * step) / c.curv;
                }
        }
    }
    {
        const double lo = -stat(curv_ratio, [](double x) { return -x; }).max, hi = stat(curv_ratio, [](double x) { return x; }).max;
        out("   curvature of 1 - F in flux from these runs / the curvature stored with the calibrations: %.2f-%.2f\n", lo, hi);
        std::vector<double> ratio, phase;
        for (const auto& R : respA) { ratio.push_back(R.fixed_mm / R.free_mm); phase.push_back(std::abs(R.phase_mover) * 1e-3); }
        out("   constant offset on the mover line (cz_offset_response): fixed costs %.1f times free (median; %.1f-%.1f); the mover's phase moves by %.2f rad per mPhi0 (median; %.2f-%.2f)\n",
            median(ratio), -stat(ratio, [](double x) { return -x; }).max, stat(ratio, [](double x) { return x; }).max, median(phase), -stat(phase, [](double x) { return -x; }).max,
            stat(phase, [](double x) { return x; }).max);
        std::vector<double> sl;
        for (double x : slope) sl.push_back(std::abs(x));
        out("   the mover's frequency at its interaction flux moves by %.2f-%.2f GHz per Phi0 (median %.2f)\n\n", -stat(sl, [](double x) { return -x; }).max, stat(sl, [](double x) { return x; }).max, median(sl));
        J.set("static_offsets", j1);
        Json o = Json::object(); o.set("min", lo); o.set("max", hi);
        J.set("curvature_ratio_to_calibration", o);
    }

    // ---------------------------------------------------------------- 2. one tail, gate on its own
    const std::vector<double> taus = opt.smoke ? std::vector<double>{30, 1000} : std::vector<double>{2, 5, 10, 30, 100, 300, 1000, 3000, 10000};
    const double a_ref = -1e-3;
    out("2. One tail on both lines (step response 1 + a exp(-t / tau)), gate on its own after a long idle (mean / max over classes)\n");
    out("       a       tau ns   programmed as calibrated   amplitudes re-calibrated   of which leakage\n");
    Json j2 = Json::array();
    struct TailCase { double a, tau; };
    std::vector<TailCase> tail_cases;
    for (double tau : taus) tail_cases.push_back({a_ref, tau});
    if (!opt.smoke) { tail_cases.push_back({-a_ref, 30}); tail_cases.push_back({3 * a_ref, 30}); }
    for (const auto& tc : tail_cases) {
        const FluxSetup F = with_tail(tc.a, tc.tau);
        const auto nom = alone(A, F, false), mat = alone(A, F, true);
        std::vector<double> dn(A.cls.size()), dm(A.cls.size()), dl(A.cls.size());
        for (std::size_t i = 0; i < A.cls.size(); ++i) { dn[i] = nom[i].err - ref[i].err; dm[i] = mat[i].err - ref[i].err; dl[i] = mat[i].leak - ref[i].leak; }
        const Stat sn = stat(dn, [](double x) { return x; }), sm = stat(dm, [](double x) { return x; }), sl = stat(dl, [](double x) { return x; });
        out("   %+.1e  %7.0f   %+.2e / %+.2e      %+.2e / %+.2e      %+.2e\n", tc.a, tc.tau, sn.mean, sn.max, sm.mean, sm.max, sl.mean);
        Json o = Json::object();
        o.set("a", tc.a); o.set("tau_ns", tc.tau); o.set("as_calibrated", jstat(sn)); o.set("amplitude_recalibrated", jstat(sm)); o.set("leakage_change_mean", sl.mean);
        j2.push(o);
    }
    out("\n");
    J.set("tail_alone", j2);

    // ---------------------------------------------------------------- 3. one tail in the memory cycle
    std::vector<Chip> chips{A};
    if (!opt.smoke) chips.push_back(load_chip(beta1(), opt.cache, max_cls));
    Json j3 = Json::array(), jtol = Json::array();
    out("3. The same tail in the memory cycle. Each gate keeps the amplitudes and virtual Z gates calibrated on its own;\n");
    out("   a) per gate: added error against the gate on its own (mean / max over the CZ gates of a round)\n");
    out("   b) per measure qubit, virtual Z gates fixed: the phase shifts of its gates add up between its two Hadamards.\n");
    out("      flip = sin^2(sum / 2), the probability that this rotation flips the outcome (mean / max over the measure qubits);\n");
    out("      gate by gate = the same phases counted independently (mean). An estimate from the mover phases of the pair simulations.\n");
    out("      first round, Z gates of the steady state = the flip of the first round when the virtual Z gates are calibrated in the running cycle\n");
    struct Tol { std::string chip; double tau, a_fixed, a_summed, a_free; };
    std::vector<Tol> tols;
    struct Est { std::string chip; double tau, sim_fixed, est_fixed, sim_free, est_free; };
    std::vector<Est> ests;
    struct FlipRow { double a, tau; Stat first, steady, settling; double first_ind, steady_ind, sum_max, offset; };
    double spread_lo = 1e300, spread_hi = 0;      // worst gate / mean, fixed, steady state, for tau >= 30 ns
    for (const auto& C : chips) {
        const auto resp = &C == &chips.front() ? respA : responses(C);
        int n4 = 0; for (const auto& L : C.lines) if (L.gates.size() == 4) ++n4;
        out("   %s: %zu CZ per round on %zu measure-qubit lines (%d with four pulses), round period %.1f ns for every line\n", C.name.c_str(), C.gates.size(), C.lines.size(), n4, C.period);
        out("   (the slowest line of the compiled schedule, whose lines run at %.1f-%.1f ns; pattern of a round after the schedule has settled)\n", C.period_min, C.period);
        out("   a)  a       tau ns   first round: fixed              free                   steady state: fixed             free                   mover phase, rad rms\n");
        std::vector<TailCase> cases;
        for (double tau : taus) cases.push_back({a_ref, tau});
        if (!opt.smoke && &C == &chips.front()) for (double tau : {10.0, 30.0, 1000.0}) cases.push_back({3 * a_ref, tau});
        if (!opt.smoke && &C == &chips.front()) cases.push_back({-a_ref, 30});
        std::vector<FlipRow> flips;
        const double ref_flip = depolarising_flip(0.1 * C.allowance, double(C.gates.size()) / double(C.lines.size()));
        for (const auto& tc : cases) {
            const FluxSetup F = with_tail(tc.a, tc.tau);
            const auto al = alone(C, F, true);
            const int rs = steady_round(C, tc.tau);
            const auto c0 = in_cycle(C, F, al, resp, 0), c1 = in_cycle(C, F, al, resp, rs);
            const Stat f0 = stat(c0, [](const InCycle& x) { return x.fixed; }), g0 = stat(c0, [](const InCycle& x) { return x.free; });
            const Stat f1 = stat(c1, [](const InCycle& x) { return x.fixed; }), g1 = stat(c1, [](const InCycle& x) { return x.free; });
            const double ph = std::sqrt(stat(c1, [](const InCycle& x) { return x.dphase * x.dphase; }).mean);
            const double ef = stat(c1, [](const InCycle& x) { return x.est_fixed; }).mean, eg = stat(c1, [](const InCycle& x) { return x.est_free; }).mean;
            const auto l0 = line_phases(C, c0), l1 = line_phases(C, c1);
            std::vector<double> settle(l0.size());     // first round with the virtual Z gates of the steady state: what is left while the lines settle
            for (std::size_t l = 0; l < l0.size(); ++l) settle[l] = std::pow(std::sin(0.5 * (l0[l].sum - l1[l].sum)), 2);
            FlipRow fr{tc.a, tc.tau, stat(l0, [](const LinePhase& x) { return x.flip; }), stat(l1, [](const LinePhase& x) { return x.flip; }), stat(settle, [](double x) { return x; }),
                       stat(l0, [](const LinePhase& x) { return x.independent; }).mean, stat(l1, [](const LinePhase& x) { return x.independent; }).mean,
                       stat(l1, [](const LinePhase& x) { return std::abs(x.sum); }).max, stat(c1, [](const InCycle& x) { return x.offset; }).mean};
            flips.push_back(fr);
            out("   %+.1e  %7.0f   %+.2e / %+.2e   %+.2e / %+.2e    %+.2e / %+.2e   %+.2e / %+.2e    %.4f\n", tc.a, tc.tau, f0.mean, f0.max, g0.mean, g0.max, f1.mean, f1.max, g1.mean, g1.max, ph);
            Json o = Json::object();
            o.set("chip", C.name); o.set("a", tc.a); o.set("tau_ns", tc.tau); o.set("steady_round", rs);
            o.set("first_round_fixed", jstat(f0)); o.set("first_round_free", jstat(g0)); o.set("steady_fixed", jstat(f1)); o.set("steady_free", jstat(g1));
            o.set("steady_mover_phase_rms_rad", ph);
            o.set("steady_fixed_estimate_mean", ef); o.set("steady_free_estimate_mean", eg);
            o.set("first_round_flip", jstat(fr.first)); o.set("steady_flip", jstat(fr.steady)); o.set("first_round_flip_with_steady_state_frames", jstat(fr.settling));
            o.set("first_round_flip_gate_by_gate_mean", fr.first_ind); o.set("steady_flip_gate_by_gate_mean", fr.steady_ind);
            o.set("steady_summed_phase_max_rad", fr.sum_max); o.set("steady_mover_flat_top_offset_mean", fr.offset);
            j3.push(o);
            if (tc.a == a_ref) ests.push_back({C.name, tc.tau, f1.mean, ef, g1.mean, eg});
            if (tc.a == a_ref) {
                const double target = 0.1 * C.allowance;
                auto tol = [&](double e, double limit) { return e > 0 ? std::abs(a_ref) * std::sqrt(limit / e) : 1e9; };
                tols.push_back({C.name, tc.tau, tol(std::max(f0.mean, f1.mean), target), tol(std::max(fr.first.mean, fr.steady.mean), ref_flip), tol(std::max(g0.mean, g1.mean), target)});
                if (tc.tau >= 30 && f1.mean > 0) { spread_lo = std::min(spread_lo, f1.max / f1.mean); spread_hi = std::max(spread_hi, f1.max / f1.mean); }
            }
        }
        out("   b)  a       tau ns   first round: flip       gate by gate   steady state: flip      gate by gate   first round, Z gates of the steady state   summed phase, rad (max)   mover flat top / step (mean)\n");
        for (const auto& fr : flips)
            out("   %+.1e  %7.0f   %.2e / %.2e     %.2e       %.2e / %.2e     %.2e       %.2e / %.2e                        %.4f                    %+.2e\n", fr.a, fr.tau, fr.first.mean, fr.first.max, fr.first_ind,
                fr.steady.mean, fr.steady.max, fr.steady_ind, fr.settling.mean, fr.settling.max, fr.sum_max, fr.offset);
        out("      for comparison, depolarising gates at a tenth of the allowance flip a measure qubit with probability %.2e per round (%.2f gates on average)\n",
            ref_flip, double(C.gates.size()) / double(C.lines.size()));
    }
    out("   largest tail amplitude |a| that keeps the added error below a tenth of the chip's residual control-error allowance, from the a^2\n");
    out("   growth of the added error%s; the larger of first round and steady state.\n", opt.smoke ? "" : " (checked above at 10, 30 and 1000 ns between 0.1 % and 0.3 %)");
    out("   fixed, per gate: mean over the gates;  fixed, summed: mean flip of a measure qubit against that of depolarising gates at the\n");
    out("   same tenth (3b);  free: mean over the gates\n");
    out("   chip                     allowance   tau ns    fixed, per gate   fixed, summed   free\n");
    const double tol_cap = 0.01;                 // no statement beyond 1 %
    auto pct = [&](double a) { char b[32]; if (a > tol_cap) std::snprintf(b, sizeof b, "> 1 %%"); else std::snprintf(b, sizeof b, "%.3f %%", 100 * a); return std::string(b); };
    for (const auto& t : tols) {
        const double allow = t.chip == A.name ? A.allowance : chips.back().allowance;
        out("   %-24s %.0e    %7.0f   %-15s   %-13s   %s\n", t.chip.c_str(), allow, t.tau, pct(t.a_fixed).c_str(), pct(t.a_summed).c_str(), pct(t.a_free).c_str());
        Json o = Json::object();
        o.set("chip", t.chip); o.set("allowance", allow); o.set("target_added_error", 0.1 * allow); o.set("tau_ns", t.tau);
        o.set("a_max_fixed", std::min(t.a_fixed, tol_cap)); o.set("a_max_fixed_summed", std::min(t.a_summed, tol_cap)); o.set("a_max_free", std::min(t.a_free, tol_cap));
        o.set("fixed_above_cap", t.a_fixed > tol_cap); o.set("fixed_summed_above_cap", t.a_summed > tol_cap); o.set("free_above_cap", t.a_free > tol_cap);
        jtol.push(o);
    }
    if (spread_hi > 0) out("   the worst gate of a round has %.1f-%.1f times the mean added error (fixed, steady state, tau >= 30 ns)\n", spread_lo, spread_hi);
    out("   quasi-static estimate / simulation, steady state (mean added error; the estimate takes the flux offsets averaged over the\n");
    out("   flat top and each gate's response to constant offsets, without a pulse simulation; shown where the simulation gives more than 1e-6):\n");
    out("   chip                     tau ns    fixed    free\n");
    auto ratio = [](double est, double sim) { char b[16]; if (sim > 1e-6) std::snprintf(b, sizeof b, "%.2f", est / sim); else std::snprintf(b, sizeof b, "-"); return std::string(b); };
    for (const auto& e : ests) out("   %-24s %7.0f    %-6s   %s\n", e.chip.c_str(), e.tau, ratio(e.est_fixed, e.sim_fixed).c_str(), ratio(e.est_free, e.sim_free).c_str());
    out("\n");
    {
        Json jc = Json::array();
        for (const auto& C : chips) {
            Json o = Json::object();
            o.set("chip", C.name); o.set("allowance", C.allowance); o.set("cz_per_round", int(C.gates.size())); o.set("measure_qubit_lines", int(C.lines.size()));
            o.set("round_period_ns", C.period); o.set("schedule_line_period_min_ns", C.period_min);
            o.set("depolarising_flip_at_a_tenth_of_the_allowance", depolarising_flip(0.1 * C.allowance, double(C.gates.size()) / double(C.lines.size())));
            jc.push(o);
        }
        J.set("cycles", jc);
        if (spread_hi > 0) { Json o = Json::object(); o.set("min", spread_lo); o.set("max", spread_hi); J.set("worst_gate_over_mean_fixed", o); }
    }
    J.set("tail_in_cycle", j3); J.set("tail_tolerance", jtol);

    // ---------------------------------------------------------------- 4. predistortion of an example line
    const std::vector<FluxSection> ex = {{-0.05, 15}, {-0.025, 150}, {-0.012, 1500}};
    {
        FluxLine L; L.sections = ex;
        out("4. Predistortion of an example line in the %s cycle. Sections (a, tau ns):", A.name.c_str());
        for (const auto& t : ex) out(" (%+.3f, %.0f)", t.a, t.tau_ns);
        out("\n   uncorrected step response: %.3f after 20 ns, %.3f after 60 ns, %.3f after 650 ns, %.4f after 2.5 us (an assumption, not a measured line)\n", L.step(20), L.step(60), L.step(650), L.step(2500));
        out("   the filter inverts a fit of the sections; no sample-and-hold, no bandwidth limit (those are in section 5)\n");
        out("   case                                   gate alone   first round: fixed              free                   steady state: fixed             free\n");
    }
    struct PreCase { std::string label; bool on; double da, dtau; int drop; };
    std::vector<PreCase> pre_cases = {{"no predistortion", false, 0, 0, -1}, {"exact fit", true, 0, 0, -1}};
    if (!opt.smoke) {
        for (double d : {0.01, 0.03, 0.10}) pre_cases.push_back({"amplitudes fitted " + std::to_string(int(std::lround(100 * d))) + " % too large", true, d, 0, -1});
        pre_cases.push_back({"amplitudes fitted 3 % too small", true, -0.03, 0, -1});
        for (double d : {0.03, 0.10}) pre_cases.push_back({"time constants fitted " + std::to_string(int(std::lround(100 * d))) + " % too long", true, 0, d, -1});
        for (int k = 0; k < 3; ++k) pre_cases.push_back({"section with tau " + std::to_string(int(ex[std::size_t(k)].tau_ns)) + " ns not corrected", true, 0, 0, k});
    } else pre_cases.push_back({"amplitudes fitted 3 % too large", true, 0.03, 0, -1});
    Json j4 = Json::array();
    double amp3_fixed = 0, amp3_free = 0, tau10_fixed = 0, tau10_free = 0;
    for (const auto& pc : pre_cases) {
        FluxSetup F; F.pad_ns = 0;
        F.mover.line.sections = ex; F.coupler.line.sections = ex;
        if (pc.on) {
            std::vector<FluxSection> fit;
            for (std::size_t k = 0; k < ex.size(); ++k) if (int(k) != pc.drop) fit.push_back({ex[k].a * (1 + pc.da), ex[k].tau_ns * (1 + pc.dtau)});
            F.mover.correction = fit; F.coupler.correction = fit;
        }
        const auto al = alone(A, F, true);
        std::vector<double> da(A.cls.size());
        for (std::size_t i = 0; i < A.cls.size(); ++i) da[i] = al[i].err - ref[i].err;
        const Stat sa = stat(da, [](double x) { return x; });
        const auto c0 = in_cycle(A, F, al, respA, 0), c1 = in_cycle(A, F, al, respA, steady_round(A, 1500 * 1.1));
        const Stat f0 = stat(c0, [](const InCycle& x) { return x.fixed; }), g0 = stat(c0, [](const InCycle& x) { return x.free; });
        const Stat f1 = stat(c1, [](const InCycle& x) { return x.fixed; }), g1 = stat(c1, [](const InCycle& x) { return x.free; });
        out("   %-38s %+.1e     %+.2e / %+.2e   %+.2e / %+.2e    %+.2e / %+.2e   %+.2e / %+.2e\n", pc.label.c_str(), sa.mean, f0.mean, f0.max, g0.mean, g0.max, f1.mean, f1.max, g1.mean, g1.max);
        Json o = Json::object();
        o.set("case", pc.label); o.set("alone_added", jstat(sa));
        o.set("first_round_fixed", jstat(f0)); o.set("first_round_free", jstat(g0)); o.set("steady_fixed", jstat(f1)); o.set("steady_free", jstat(g1));
        j4.push(o);
        if (pc.on && pc.da == 0.03 && pc.dtau == 0) { amp3_fixed = std::max(f0.mean, f1.mean); amp3_free = std::max(g0.mean, g1.mean); }
        if (pc.on && pc.dtau == 0.10) { tau10_fixed = std::max(f0.mean, f1.mean); tau10_free = std::max(g0.mean, g1.mean); }
    }
    out("\n");
    {
        Json o = Json::object();
        Json jt = Json::array();
        for (const auto& t : ex) { Json e = Json::object(); e.set("a", t.a); e.set("tau_ns", t.tau_ns); jt.push(e); }
        o.set("sections", jt); o.set("cases", j4);
        J.set("predistortion_example", o);
    }

    // ---------------------------------------------------------------- 5. AWG and timing
    out("5. AWG and timing\n");
    out("   a) sample-and-hold and bandwidth (single-pole low-pass, 10-90 %% rise time), no tails; three pair classes, gate on its own:\n");
    out("      amplitudes re-calibrated only, and after a local re-calibration of the three knobs (mover step, flat-top length, coupler step)\n");
    out("      sample rate   rise ns    class          amplitudes only   three knobs   leakage     flat top\n");
    // three classes that the plain 1 - F calibration produced (not the leakage fallback); the smoke run takes the first class
    auto find_class = [&](double fm, double ft) {
        for (std::size_t i = 0; i < A.cls.size(); ++i) if (std::abs(A.cls[i].fm - fm) < 1e-6 && std::abs(A.cls[i].ft - ft) < 1e-6) return int(i);
        throw std::runtime_error("flux_distortion_study: the CZ cache has no " + std::to_string(fm) + " / " + std::to_string(ft) + " GHz pair class");
    };
    std::vector<int> pick;
    if (opt.smoke) pick.push_back(0);
    else pick = {find_class(5.99, 5.29), find_class(5.87, 5.29), find_class(5.87, 5.11)};
    struct AwgCase { double gs, rise; };
    const std::vector<AwgCase> awg_cases = opt.smoke ? std::vector<AwgCase>{{2.4, 1.0}} : std::vector<AwgCase>{{0, 1.0}, {0, 2.0}, {2.4, 0.5}, {2.4, 1.0}, {2.4, 2.0}, {1.0, 1.0}, {1.0, 2.0}};
    const int max_eval = opt.smoke ? 12 : 80;
    Json j5a = Json::array();
    std::vector<Recal> at24(pick.size());
    std::string awg_ok;                          // combinations whose re-calibrated gates stay within a tenth of the allowance
    for (const auto& ac : awg_cases) {
        FluxSetup F; F.awg.sample_ns = ac.gs > 0 ? 1.0 / ac.gs : 0.0;
        F.mover.line.rise_ns = ac.rise; F.coupler.line.rise_ns = ac.rise;
        std::vector<Recal> rc(pick.size()); std::vector<double> am(pick.size());
        LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < int(pick.size()); ++k) guard.run([&] {
            const auto& c = A.cls[std::size_t(pick[std::size_t(k)])];
            am[std::size_t(k)] = 1 - simulate_cz_drive(c.S, cz_gate_drive(c.S, c.P, F, cz_matched_command(c.S, c.P, F))).fidelity;
            rc[std::size_t(k)] = recalibrate(c, F, max_eval);
        });
        guard.rethrow();
        for (std::size_t k = 0; k < pick.size(); ++k) {
            const auto& c = A.cls[std::size_t(pick[k])];
            const double e0 = ref[std::size_t(pick[k])].err;
            char rate[32];
            if (ac.gs > 0) std::snprintf(rate, sizeof rate, "%.1f GS/s", ac.gs); else std::snprintf(rate, sizeof rate, "none");
            out("      %-12s  %.1f        %.2f/%.2f      %+.2e         %+.2e     %.2e    %+.2f ns\n", rate, ac.rise, c.fm, c.ft, am[k] - e0, rc[k].err - e0, rc[k].leak, rc[k].P.t_flat - c.P.t_flat);
            Json o = Json::object();
            o.set("sample_rate_GSps", ac.gs); o.set("rise_ns", ac.rise); o.set("f_mover_GHz", c.fm); o.set("f_target_GHz", c.ft);
            o.set("added_amplitudes_only", am[k] - e0); o.set("added_three_knobs", rc[k].err - e0); o.set("leakage_three_knobs", rc[k].leak);
            o.set("added_leakage_three_knobs", rc[k].leak - ref[std::size_t(pick[k])].leak);
            o.set("flat_top_change_ns", rc[k].P.t_flat - c.P.t_flat);
            j5a.push(o);
        }
        if (ac.gs == 2.4 && ac.rise == 1.0) at24 = rc;
        double worst = -1e300;
        for (std::size_t k = 0; k < pick.size(); ++k) worst = std::max(worst, rc[k].err - ref[std::size_t(pick[k])].err);
        if (worst <= 0.1 * A.allowance) {
            char b[64];
            if (ac.gs > 0) std::snprintf(b, sizeof b, "%.1f GS/s with %.1f ns", ac.gs, ac.rise); else std::snprintf(b, sizeof b, "no sampling with %.1f ns", ac.rise);
            awg_ok += (awg_ok.empty() ? "" : "; ") + std::string(b);
        }
    }
    J.set("awg_sampling", j5a);

    out("   b) DAC resolution, full scale 1 Phi0 peak to peak: the re-calibrated gates of 2.4 GS/s and 1.0 ns played through the DAC (max over the classes),\n");
    out("      and the bound for a flat top that is off by half a step on both lines, in the worse combination of signs (response to constant offsets, all classes)\n");
    out("      bits   step, Phi0    added error   bound\n");
    Json j5b = Json::array();
    int bits_needed = 0;
    for (int bits : opt.smoke ? std::vector<int>{14} : std::vector<int>{10, 12, 14, 16}) {
        FluxSetup F; F.awg.sample_ns = 1.0 / 2.4; F.mover.line.rise_ns = 1.0; F.coupler.line.rise_ns = 1.0;
        FluxSetup Fq = F; Fq.awg.bits = bits; Fq.awg.full_scale = 1.0;
        std::vector<double> d(pick.size());
        LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < int(pick.size()); ++k) guard.run([&] {
            const auto& c = A.cls[std::size_t(pick[std::size_t(k)])];
            const auto& rc = at24[std::size_t(k)];
            d[std::size_t(k)] = (1 - simulate_cz_drive(c.S, cz_gate_drive(c.S, rc.P, Fq, rc.cmd)).fidelity) - rc.err;
        });
        guard.rethrow();
        const double lsb = 1.0 / std::ldexp(1.0, bits);
        double bound = 0;                       // half a step on each line, either sign: the worst class and combination (free virtual Z gates)
        for (const auto& R : respA) for (int sm : {+1, -1}) for (int sc : {+1, -1}) bound = std::max(bound, R.added_free(sm * 0.5 * lsb, sc * 0.5 * lsb));
        const Stat s = stat(d, [](double x) { return x; });
        out("      %2d     %.1e      %+.1e      %.1e\n", bits, lsb, s.max, bound);
        if (!bits_needed && bound <= 0.1 * A.allowance) bits_needed = bits;
        Json o = Json::object(); o.set("bits", bits); o.set("step_Phi0", lsb); o.set("added_max", s.max); o.set("bound_half_step", bound);
        j5b.push(o);
    }
    J.set("dac_resolution", j5b);

    out("   c) timing skew between the mover and the coupler line (free). All %zu classes as calibrated (median / worst); the three classes of 5a\n", A.cls.size());
    out("      as calibrated and after the local three-knob re-calibration with the skew present (worst of the three)\n");
    out("      late line   by ns    as calibrated            three classes: as calibrated   three knobs\n");
    Json j5c = Json::array();
    double skew_fit = 0, recal_keeps_lo = 1e300, recal_keeps_hi = 0;   // worst added error / skew^2; share of the error the re-calibration leaves
    for (double d : opt.smoke ? std::vector<double>{0.2} : std::vector<double>{0.05, 0.1, 0.2, 0.5}) for (int which = 0; which < 2; ++which) {
        FluxSetup F; F.pad_ns = d;
        (which == 0 ? F.coupler.line.delay_ns : F.mover.line.delay_ns) = d;
        FluxSetup F0 = F; F0.coupler.line.delay_ns = 0; F0.mover.line.delay_ns = 0;     // same time grid without the skew
        const auto base = alone(A, F0, false), sk = alone(A, F, false);
        std::vector<double> v(A.cls.size());
        for (std::size_t i = 0; i < A.cls.size(); ++i) v[i] = sk[i].err - base[i].err;
        const Stat s = stat(v, [](double x) { return x; });
        std::vector<double> rk(pick.size());
        LoopGuard guard;
#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < int(pick.size()); ++k)
            guard.run([&] { rk[std::size_t(k)] = recalibrate(A.cls[std::size_t(pick[std::size_t(k)])], F, max_eval).err - base[std::size_t(pick[std::size_t(k)])].err; });
        guard.rethrow();
        double as_cal = -1e300, three = -1e300;
        for (std::size_t k = 0; k < pick.size(); ++k) { as_cal = std::max(as_cal, v[std::size_t(pick[k])]); three = std::max(three, rk[k]); }
        out("      %-8s    %.2f     %.2e / %.2e      %.2e                       %.2e\n", which == 0 ? "coupler" : "mover", d, median(v), s.max, as_cal, three);
        Json o = Json::object(); o.set("late_line", which == 0 ? "coupler" : "mover"); o.set("skew_ns", d); o.set("added_median", median(v)); o.set("added_worst", s.max);
        o.set("three_classes_as_calibrated_worst", as_cal); o.set("three_classes_three_knobs_worst", three);
        j5c.push(o);
        skew_fit = std::max(skew_fit, s.max / (d * d));
        if (as_cal > 1e-6) { recal_keeps_lo = std::min(recal_keeps_lo, three / as_cal); recal_keeps_hi = std::max(recal_keeps_hi, three / as_cal); }
    }
    if (recal_keeps_hi > 0) out("      the re-calibration leaves %.0f-%.0f %% of the added error of the three classes\n", 100 * recal_keeps_lo, 100 * recal_keeps_hi);
    out("\n");
    J.set("timing_skew", j5c);

    // ---------------------------------------------------------------- 6. the class of the leakage fallback
    if (!opt.smoke) {
        const int k = find_class(5.81, 5.29);
        const auto& c = A.cls[std::size_t(k)];
        const Recal rp = recalibrate(c, ideal, 400);
        const CZCommand n0 = cz_nominal_command(c.S, c.P);
        out("6. The %.2f/%.2f GHz class: the cache entry that the leakage fallback of CZ model version 4 changed (results/PROVENANCE.md).\n", c.fm, c.ft);
        out("   Stored pulse: 1 - F %.2e, leakage %.2e. Re-calibrated from it for 1 - F alone (three knobs, ideal lines, 400 evaluations):\n", ref[std::size_t(k)].err, ref[std::size_t(k)].leak);
        out("   1 - F %.2e, leakage %.2e (flat top %+.2f ns, mover step %+.1e Phi0, coupler step %+.1e Phi0); limits of the Alpha gate: %.0e and %.0e.\n",
            rp.err, rp.leak, rp.P.t_flat - c.P.t_flat, rp.cmd.amp_mover - n0.amp_mover, rp.cmd.amp_coupler - n0.amp_coupler, CZModel::coh_limit, CZModel::leak_limit);
        out("   Nothing uses this pulse: the cache and the reviews keep the stored one.\n\n");
        Json o = Json::object();
        o.set("f_mover_GHz", c.fm); o.set("f_target_GHz", c.ft);
        o.set("stored_err", ref[std::size_t(k)].err); o.set("stored_leakage", ref[std::size_t(k)].leak);
        o.set("repolished_err", rp.err); o.set("repolished_leakage", rp.leak);
        o.set("flat_top_change_ns", rp.P.t_flat - c.P.t_flat); o.set("mover_step_change_Phi0", rp.cmd.amp_mover - n0.amp_mover); o.set("coupler_step_change_Phi0", rp.cmd.amp_coupler - n0.amp_coupler);
        o.set("coherent_limit", CZModel::coh_limit); o.set("leakage_limit", CZModel::leak_limit);
        J.set("leakage_fallback_class_repolished", o);
    }

    // ---------------------------------------------------------------- summary
    out("Summary: what keeps each contribution below a tenth of the residual control-error allowance\n");
    Json js = Json::object();
    {
        out("   residual tail of the corrected step response, at the worst time constant:\n");
        Json jt = Json::array();
        for (const auto& C : chips) {
            double tf = 1, ts = 1, tr = 1;
            for (const auto& t : tols) if (t.chip == C.name) { tf = std::min(tf, t.a_fixed); ts = std::min(ts, t.a_summed); tr = std::min(tr, t.a_free); }
            out("      %s (allowance %.0e): %.3f %% (fixed, per gate), %.3f %% (fixed, phases summed per measure qubit), %.3f %% (free)\n", C.name.c_str(), C.allowance, 100 * tf, 100 * ts, 100 * tr);
            Json o = Json::object(); o.set("chip", C.name); o.set("tail_amplitude_fixed", tf); o.set("tail_amplitude_fixed_summed", ts); o.set("tail_amplitude_free", tr);
            jt.push(o);
        }
        js.set("tail_tolerance", jt);
        const double target = 0.1 * A.allowance;
        out("   %s:\n", A.name.c_str());
        if (amp3_fixed > 0 && amp3_free > 0) {
            const double df = 0.03 * std::sqrt(target / amp3_fixed), dr = 0.03 * std::sqrt(target / amp3_free);
            out("      example line, fit of the section amplitudes: within %.1f %% (fixed), %.1f %% (free)\n", 100 * df, 100 * dr);
            js.set("example_amplitude_fit_fixed", df); js.set("example_amplitude_fit_free", dr);
        }
        if (tau10_fixed > 0 && tau10_free > 0) {
            const double df = 0.10 * std::sqrt(target / tau10_fixed), dr = 0.10 * std::sqrt(target / tau10_free);
            out("      example line, fit of the time constants: within %.1f %% (fixed), %.1f %% (free)\n", 100 * df, 100 * dr);
            js.set("example_tau_fit_fixed", df); js.set("example_tau_fit_free", dr);
        }
        if (!awg_ok.empty()) out("      sample rate with rise time, after re-calibration: %s\n", awg_ok.c_str());
        js.set("awg_within_target", awg_ok);
        if (bits_needed) out("      DAC resolution at 1 Phi0 full scale: %d bits\n", bits_needed);
        js.set("dac_bits", bits_needed);
        if (skew_fit > 0) { const double sk = std::sqrt(target / skew_fit); out("      timing skew between mover and coupler line: %.2f ns (gates as calibrated; the three-knob re-calibration does not remove it, 5c)\n", sk); js.set("skew_ns", sk); }
    }
    J.set("summary", js);

    const double runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    out("\nruntime %.0f s\n", runtime);
    std::fclose(g_txt); g_txt = nullptr;
    J.set("runtime_s", runtime);
    std::ofstream f(json_path);
    J.dump(f); f << "\n";
    std::printf("wrote %s and %s\n", txt_path.c_str(), json_path.c_str());
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try { return run(argc, argv); }
    catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 2; }
}
