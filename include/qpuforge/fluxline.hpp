// qpuforge — flux-line distortion and predistortion for the pulse-level CZ
//
// simulate_cz (pulse.hpp) assumes that the flux at each SQUID follows the programmed sin^2 pulse exactly. On hardware
// the pulse first passes through a linear system: the AWG (sample-and-hold, finite resolution), amplifier, bias tee,
// filters, cable and the on-chip line. This header models that system and its digital correction (predistortion) and
// runs the CZ of pulse.hpp on the flux that actually arrives:
//   FluxLine           response of one line: first-order sections in cascade (each with the step response
//                      1 + a exp(-t / tau)), an optional single-pole low-pass and a delay
//   FluxLineFilter     that response for a drive held constant over each interval (exact, no time-stepping error)
//   FluxChannel        a line and the sections its predistortion filter was built from (a fit of the line)
//   flux_line_output   flux delivered by a channel for a train of pulses, averaged over given time intervals
//   CZDrive            piecewise-constant mover and coupler flux; simulate_cz_drive propagates it with simulate_cz's
//                      propagator and fidelity metric (cz_ideal_drive reproduces simulate_cz)
//   CZFrame            the virtual Z gates of a CZ; cz_fidelity_in_frame scores a gate with corrections fixed elsewhere
//   cz_offset_response the added error per unit of constant flux offset, with free and with fixed virtual Z gates
//   cz_gate_drive      the drive of one CZ on two channels; cz_train_drive the same gate inside pulse trains (history)
// Nothing here is used by calibrate_chip_cz, so no review output or cache depends on it.
// docs/design/flux-line-predistortion.md describes the model and its limits; examples/flux_distortion_study.cpp uses it.
#pragma once
#include "compiler.hpp"
#include "pulse.hpp"
#include <numeric>
#include <stdexcept>

namespace qpuforge {

// One first-order section of a flux line. Step response: 1 + a exp(-t / tau).
//   a > 0: an overshoot that decays (high-pass character, e.g. a bias tee);
//   a < 0: slow settling towards the final value (skin effect, filters).
struct FluxSection { double a{0}, tau_ns{1}; };

struct FluxLine {
    std::vector<FluxSection> sections;   // in cascade; for small amplitudes the step response is ~ 1 + sum a_k exp(-t / tau_k)
    double rise_ns{0};             // 10-90 % rise time of a single-pole low-pass (0: none)
    double delay_ns{0};            // pure delay

    bool ideal() const {
        for (const auto& t : sections) if (t.a != 0) return false;
        return rise_ns <= 0 && delay_ns == 0;
    }
    double lp_tau_ns() const { return rise_ns > 0 ? rise_ns / std::log(9.0) : 0.0; }

    // Partial fractions: H(s) = c0 + sum_j r_j p_j / (s + p_j), poles p_j in 1/ns, DC gain c0 + sum_j r_j = 1.
    // Each term is a one-pole low-pass, so the response to a held drive has a closed form (FluxLineFilter).
    struct Parallel { double c0{1}; std::vector<double> p, r; };
    Parallel parallel() const {
        std::vector<double> a, p;
        for (const auto& t : sections) {
            if (t.a == 0) continue;
            if (!(t.tau_ns > 0) || !(t.a > -1)) throw std::invalid_argument("FluxLine: a section needs tau > 0 and a > -1");
            a.push_back(t.a); p.push_back(1.0 / t.tau_ns);
        }
        if (rise_ns < 0 || delay_ns < 0) throw std::invalid_argument("FluxLine: rise time and delay must not be negative");
        const double pl = rise_ns > 0 ? std::log(9.0) / rise_ns : 0.0;
        auto same = [](double x, double y) { return std::abs(x - y) <= 1e-9 * std::max(x, y); };
        for (std::size_t i = 0; i < p.size(); ++i) {
            if (pl > 0 && same(p[i], pl)) throw std::invalid_argument("FluxLine: a section has the time constant of the low-pass");
            for (std::size_t j = 0; j < i; ++j) if (same(p[i], p[j])) throw std::invalid_argument("FluxLine: two sections share a time constant");
        }
        Parallel P;
        P.c0 = pl > 0 ? 0.0 : 1.0;
        if (pl <= 0) for (double ak : a) P.c0 *= 1 + ak;
        for (std::size_t k = 0; k < p.size(); ++k) {
            double r = -a[k];
            for (std::size_t j = 0; j < p.size(); ++j) if (j != k) r *= 1 - a[j] * p[k] / (p[j] - p[k]);
            if (pl > 0) r *= pl / (pl - p[k]);
            P.p.push_back(p[k]); P.r.push_back(r);
        }
        if (pl > 0) {
            double r = 1;
            for (std::size_t k = 0; k < p.size(); ++k) r *= 1 - a[k] * pl / (p[k] - pl);
            P.p.push_back(pl); P.r.push_back(r);
        }
        return P;
    }
    // Step response t ns after a unit step of the drive.
    double step(double t) const {
        t -= delay_ns;
        if (t < 0) return 0;
        const Parallel P = parallel();
        double y = P.c0;
        for (std::size_t j = 0; j < P.p.size(); ++j) y -= P.r[j] * std::expm1(-P.p[j] * t);
        return y;
    }
};

// Running state of a line: one low-pass per pole. Exact for a drive that is constant over each interval.
class FluxLineFilter {
public:
    explicit FluxLineFilter(const FluxLine& L, double x_rest = 0.0) : P_(L.parallel()), w_(P_.p.size(), x_rest) {}
    // Hold the drive at x for h ns. Returns the integral of the output over the interval and advances the state.
    double advance(double x, double h) {
        double acc = P_.c0 * x * h;
        for (std::size_t j = 0; j < w_.size(); ++j) {
            const double q = std::exp(-P_.p[j] * h), e = -std::expm1(-P_.p[j] * h);   // exp(-p h) and 1 - exp(-p h)
            const double d = w_[j] - x;
            acc += P_.r[j] * (x * h + d * e / P_.p[j]);
            w_[j] = x + d * q;
        }
        return acc;
    }
    double value(double x) const {             // output now, for the drive value x
        double y = P_.c0 * x;
        for (std::size_t j = 0; j < w_.size(); ++j) y += P_.r[j] * w_[j];
        return y;
    }
private:
    FluxLine::Parallel P_;
    std::vector<double> w_;
};

// The pulse shape of simulate_cz: sin^2 ramp 0 -> 1 over t_rise, 1 for t_flat, sin^2 ramp back; 0 outside.
inline double cz_pulse_shape(double t, double t_rise, double t_flat) {
    const double T = 2 * t_rise + t_flat;
    if (t <= 0 || t >= T) return 0;
    if (t >= t_rise && t <= t_rise + t_flat) return 1;
    const double s = std::sin(0.5 * PI * (t < t_rise ? t : T - t) / t_rise);
    return s * s;
}

struct FluxPulse {                 // one pulse on a line
    double t0{0};                  // start, ns
    double amp{0};                 // programmed flux step, Phi0
    double t_rise{5}, t_flat{20};
    double end() const { return t0 + 2 * t_rise + t_flat; }
};

struct FluxAwg {
    double sample_ns{0};           // sample period; 0: no sample-and-hold (the waveform is then held on a 0.05 ns grid)
    int bits{0};                   // DAC resolution; 0: not quantised
    double full_scale{1.0};        // DAC range in flux, Phi0 peak to peak
    double period() const { return sample_ns > 0 ? sample_ns : 0.05; }
};

struct FluxChannel {
    FluxLine line;                       // what the hardware does
    std::vector<FluxSection> correction; // the sections the predistortion inverts: a fit of line.sections (empty: none)
};

// Digital predistortion at the AWG sample rate: one first-order inverse per fitted section, in cascade. For a single
// section with an exact fit, the line output averaged over each sample equals the programmed sample exactly. The
// low-pass and the delay of a line are not corrected.
class Predistorter {
public:
    Predistorter(const std::vector<FluxSection>& fit, double sample_ns) {
        for (const auto& t : fit) {
            if (t.a == 0) continue;
            if (!(t.tau_ns > 0) || !(t.a > -1)) throw std::invalid_argument("Predistorter: a section needs tau > 0 and a > -1");
            const double al = std::exp(-sample_ns / t.tau_ns);
            s_.push_back({-t.a * std::expm1(-sample_ns / t.tau_ns) * t.tau_ns / sample_ns, al, 0.0});
        }
    }
    bool active() const { return !s_.empty(); }
    double operator()(double x) {
        for (auto& s : s_) {
            x = (x + s.ak * s.w) / (1 + s.ak);
            s.w = s.al * s.w + (1 - s.al) * x;
        }
        return x;
    }
private:
    struct Section { double ak, al, w; };   // a * (mean weight of the tail over one sample), exp(-Ts / tau), low-pass state
    std::vector<Section> s_;
};

// Flux delivered by a channel (relative to its idle flux) for a train of pulses, averaged over the intervals between
// `edges` (ascending, ns). The programmed waveform is sampled at the middle of each AWG sample, predistorted, quantised
// and held; the line is at rest before the first pulse. Pulses must be in time order and must not overlap.
inline std::vector<double> flux_line_output(const FluxChannel& ch, const FluxAwg& awg, const std::vector<FluxPulse>& train,
                                            const std::vector<double>& edges) {
    if (edges.size() < 2) return {};
    for (std::size_t i = 1; i < edges.size(); ++i) if (edges[i] < edges[i - 1]) throw std::invalid_argument("flux_line_output: edges must ascend");
    for (std::size_t i = 1; i < train.size(); ++i) if (train[i].t0 < train[i - 1].end()) throw std::invalid_argument("flux_line_output: pulses overlap or are out of order");
    std::vector<double> out(edges.size() - 1, 0.0);
    if (train.empty()) return out;
    const double Ts = awg.period(), delay = ch.line.delay_ns;
    const double lsb = awg.bits > 0 ? awg.full_scale / std::ldexp(1.0, awg.bits) : 0.0;
    FluxLineFilter filt(ch.line);
    Predistorter pre(ch.correction, Ts);
    const long n0 = long(std::floor(train.front().t0 / Ts));
    const long n1 = std::max(n0, long(std::ceil((edges.back() - delay) / Ts)));
    double t = std::min(edges.front(), double(n0) * Ts + delay);
    std::size_t ie = 0;                                    // first edge after t
    auto hold = [&](double x, double t_to) {               // drive x from t to t_to, split at the edges
        while (t < t_to) {
            while (ie < edges.size() && edges[ie] <= t) ++ie;
            const double t_next = (ie < edges.size() && edges[ie] < t_to) ? edges[ie] : t_to;
            const double integral = filt.advance(x, t_next - t);
            if (ie >= 1 && ie < edges.size()) out[ie - 1] += integral;
            t = t_next;
        }
    };
    hold(0.0, double(n0) * Ts + delay);
    std::size_t ip = 0;
    for (long n = n0; n < n1; ++n) {
        const double tm = (double(n) + 0.5) * Ts;
        while (ip < train.size() && tm >= train[ip].end()) ++ip;
        double x = 0;
        if (ip < train.size() && tm > train[ip].t0) x = train[ip].amp * cz_pulse_shape(tm - train[ip].t0, train[ip].t_rise, train[ip].t_flat);
        if (x == 0 && !pre.active()) {                     // nothing is played until the next pulse: one closed-form step
            long nn = n1;
            if (ip < train.size()) nn = std::min(n1, std::max(n + 1, long(std::floor(train[ip].t0 / Ts + 0.5))));
            hold(0.0, double(nn) * Ts + delay);
            n = nn - 1;
            continue;
        }
        x = pre(x);
        if (lsb > 0) x = lsb * std::round(x / lsb);
        hold(x, double(n + 1) * Ts + delay);
    }
    for (std::size_t i = 0; i + 1 < edges.size(); ++i) { const double h = edges[i + 1] - edges[i]; if (h > 0) out[i] /= h; }
    return out;
}

// ---------------------------------------------------------------- CZ on a given flux drive

struct CZDrive {
    std::vector<double> h, f1, fc;     // step lengths (ns), mover flux and coupler flux (Phi0) on each step
    double duration() const { return std::accumulate(h.begin(), h.end(), 0.0); }
};

// simulate_cz's own discretisation of an undistorted pulse: midpoint rule on the ramps, one step on the flat top.
inline CZDrive cz_ideal_drive(const CZSystem& S, const CZPulse& P, double dt = CZModel::dt_ns) {
    CZDrive W;
    auto add = [&](double s, double tau) {
        W.h.push_back(tau);
        W.f1.push_back(S.flux1_idle + s * (P.flux1_int - S.flux1_idle));
        W.fc.push_back(S.fluxc_idle + s * (P.fluxc_on - S.fluxc_idle));
    };
    const int nr = std::max(1, int(std::ceil(P.t_rise / dt)));
    const double h = P.t_rise / nr;
    for (int k = 0; k < nr; ++k) { const double x = (k + 0.5) * h / P.t_rise; add(std::pow(std::sin(PI * x / 2), 2), h); }
    if (P.t_flat > 0) add(1.0, P.t_flat);
    for (int k = 0; k < nr; ++k) { const double x = 1 - (k + 0.5) * h / P.t_rise; add(std::pow(std::sin(PI * x / 2), 2), h); }
    return W;
}

// simulate_cz for an arbitrary drive: the same propagator, rotating frame and fidelity metric (free single-qubit Z
// corrections). The time after the pulse counts as part of the gate.
inline CZResult simulate_cz_drive(const CZSystem& S, const CZDrive& W) {
    if (W.h.size() != W.f1.size() || W.h.size() != W.fc.size()) throw std::invalid_argument("simulate_cz_drive: drive arrays differ in length");
    const int D = S.dim();
    const auto d0 = S.dressed(S.flux1_idle, S.fluxc_idle, S.flux2_idle);
    std::array<std::vector<cplx>, 4> psi;
    for (int j = 0; j < 4; ++j) { psi[std::size_t(j)].resize(std::size_t(D)); for (int r = 0; r < D; ++r) psi[std::size_t(j)][std::size_t(r)] = d0.vec[std::size_t(j)][std::size_t(r)]; }
    std::vector<cplx> tmp(static_cast<std::size_t>(D)), work(static_cast<std::size_t>(D));
    std::vector<double> V, ev;
    double f1_prev = 0, fc_prev = 0;
    bool have = false;
    for (std::size_t s = 0; s < W.h.size(); ++s) {
        const double tau = W.h[s];
        if (tau <= 0) continue;
        if (!have || W.f1[s] != f1_prev || W.fc[s] != fc_prev) {       // same flux as the last step: same eigenbasis
            ev = sym_eigen_tri(S.hamiltonian(W.f1[s], W.fc[s], S.flux2_idle), D, V);
            f1_prev = W.f1[s]; fc_prev = W.fc[s]; have = true;
        }
        for (auto& p : psi) {
            for (int m = 0; m < D; ++m) {          // tmp = V^T p
                cplx acc = 0;
                for (int r = 0; r < D; ++r) acc += V[std::size_t(r) * D + m] * p[std::size_t(r)];
                tmp[std::size_t(m)] = acc * std::exp(cplx(0, -2 * PI * ev[std::size_t(m)] * tau));
            }
            for (int r = 0; r < D; ++r) {          // p = V tmp
                cplx acc = 0;
                for (int m = 0; m < D; ++m) acc += V[std::size_t(r) * D + m] * tmp[std::size_t(m)];
                work[std::size_t(r)] = acc;
            }
            p.swap(work);
        }
    }
    CZResult R;
    R.duration_ns = W.duration();
    const std::array<double, 4> Ef = {d0.E[0], d0.E[1], d0.E[2], d0.E[1] + d0.E[2] - d0.E[0]};
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
        cplx acc = 0;
        for (int r = 0; r < D; ++r) acc += d0.vec[std::size_t(i)][std::size_t(r)] * psi[std::size_t(j)][std::size_t(r)];
        R.M[std::size_t(i)][std::size_t(j)] = acc * std::exp(cplx(0, 2 * PI * Ef[std::size_t(i)] * R.duration_ns));
    }
    const double th0 = std::arg(R.M[0][0]), th2 = std::arg(R.M[1][1]) - th0, th1 = std::arg(R.M[2][2]) - th0;
    const int x1[4] = {0, 0, 1, 1}, x2[4] = {0, 1, 0, 1};
    std::array<std::array<cplx, 4>, 4> Mc = R.M;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) Mc[std::size_t(i)][std::size_t(j)] *= std::exp(cplx(0, -(th0 + x1[i] * th1 + x2[i] * th2)));
    double tr_mm = 0;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) tr_mm += std::norm(Mc[std::size_t(i)][std::size_t(j)]);
    const cplx tr = Mc[0][0] + Mc[1][1] + Mc[2][2] - Mc[3][3];
    R.fidelity = (std::norm(tr) + tr_mm) / 20.0;
    R.leakage = 1 - tr_mm / 4;
    double ph = std::arg(Mc[3][3]);
    if (ph < 0) ph += 2 * PI;
    R.cphase = ph;
    return R;
}

// Single-qubit Z corrections that follow a CZ (virtual Z gates), rad. simulate_cz and simulate_cz_drive choose them
// freely for every result; on hardware they are calibrated once, in one context, and then kept.
struct CZFrame { double th_mover{0}, th_target{0}; };
inline CZFrame cz_frame(const CZResult& R) {
    const double th0 = std::arg(R.M[0][0]);
    return {std::arg(R.M[2][2]) - th0, std::arg(R.M[1][1]) - th0};
}
// Average gate fidelity of a simulated gate when its virtual Z gates are those of `frame` (calibrated elsewhere), so
// that a shift of the single-qubit phases counts as an error.
inline double cz_fidelity_in_frame(const CZResult& R, const CZFrame& frame) {
    const int x1[4] = {0, 0, 1, 1}, x2[4] = {0, 1, 0, 1};
    const double sign[4] = {1, 1, 1, -1};
    cplx tr = 0;
    double tr_mm = 0;
    for (int i = 0; i < 4; ++i) {
        tr += sign[i] * R.M[std::size_t(i)][std::size_t(i)] * std::exp(cplx(0, -(x1[i] * frame.th_mover + x2[i] * frame.th_target)));
        for (int j = 0; j < 4; ++j) tr_mm += std::norm(R.M[std::size_t(i)][std::size_t(j)]);
    }
    return (std::norm(tr) + tr_mm) / 20.0;
}

// Response of a gate to small constant flux offsets during its drive, to second order:
//   added 1 - F = m dm + c dc + mm dm^2 + cc dc^2 + mc dm dc      (dm on the mover line, dc on the coupler line, Phi0),
// with the virtual Z gates re-calibrated for the offset (free) and with those of the undisturbed gate (fixed). The
// fixed form includes the shift of the single-qubit phases, which is first order in the offset and usually the
// largest term. The linear coefficients vanish for a gate at its optimum; a residual conditional-phase error makes
// them finite. This is what a quasi-static estimate of pulse history, drift or DAC steps needs; the curvature stored
// with a calibration (CZCalibration::err_flux_quasistatic) is 2 (free_mm + free_cc) before its noise scaling.
// Seven simulations.
struct CZOffsetResponse {
    double free_m{0}, free_c{0}, free_mm{0}, free_cc{0}, free_mc{0};
    double fixed_m{0}, fixed_c{0}, fixed_mm{0}, fixed_cc{0}, fixed_mc{0};
    double phase_mover{0};          // shift of the mover's virtual Z per Phi0 on the mover line, rad
    double added_free(double dm, double dc) const { return free_m * dm + free_c * dc + free_mm * dm * dm + free_cc * dc * dc + free_mc * dm * dc; }
    double added_fixed(double dm, double dc) const { return fixed_m * dm + fixed_c * dc + fixed_mm * dm * dm + fixed_cc * dc * dc + fixed_mc * dm * dc; }
};
inline CZOffsetResponse cz_offset_response(const CZSystem& S, const CZDrive& W, double h = 2e-5) {
    const auto r0 = simulate_cz_drive(S, W);
    const CZFrame f0 = cz_frame(r0);
    struct Run { double free, fixed, th; };
    auto run = [&](double dm, double dc) {
        CZDrive V = W;
        for (auto& f : V.f1) f += dm;
        for (auto& f : V.fc) f += dc;
        const auto r = simulate_cz_drive(S, V);
        return Run{r0.fidelity - r.fidelity, r0.fidelity - cz_fidelity_in_frame(r, f0), cz_frame(r).th_mover};
    };
    const Run mp = run(h, 0), mn = run(-h, 0), cp = run(0, h), cn = run(0, -h), bp = run(h, h), bn = run(-h, -h);
    CZOffsetResponse R;
    const double q = 1 / (2 * h * h), l = 1 / (2 * h);   // sums give the quadratic terms, differences the linear ones
    R.free_m = (mp.free - mn.free) * l; R.free_c = (cp.free - cn.free) * l;
    R.free_mm = (mp.free + mn.free) * q; R.free_cc = (cp.free + cn.free) * q; R.free_mc = (bp.free + bn.free) * q - R.free_mm - R.free_cc;
    R.fixed_m = (mp.fixed - mn.fixed) * l; R.fixed_c = (cp.fixed - cn.fixed) * l;
    R.fixed_mm = (mp.fixed + mn.fixed) * q; R.fixed_cc = (cp.fixed + cn.fixed) * q; R.fixed_mc = (bp.fixed + bn.fixed) * q - R.fixed_mm - R.fixed_cc;
    R.phase_mover = std::remainder(mp.th - mn.th, 2 * PI) / (2 * h);
    return R;
}

// ---------------------------------------------------------------- a CZ on two channels

struct FluxSetup {
    FluxChannel mover, coupler;
    FluxAwg awg;
    double dt_ns{CZModel::dt_ns};   // step on the ramps and in the settling time after the pulse
    double flat_dt_ns{0.5};         // step on the flat top
    double pad_ns{-1};              // time simulated after the nominal end of the pulse; < 0: settle_ns()
    // enough for a delayed, low-passed falling edge to finish (the tails that remain only shift single-qubit phases)
    double settle_ns() const {
        if (pad_ns >= 0) return pad_ns;
        auto need = [](const FluxLine& L) { return L.delay_ns + 8 * L.lp_tau_ns(); };
        return std::max(need(mover.line), need(coupler.line));
    }
};

// Time grid of a gate that starts at t0: ramp steps of dt, flat-top steps of flat_dt, ramp, then the settling time.
inline std::vector<double> cz_gate_edges(const CZPulse& P, const FluxSetup& F, double t0 = 0) {
    std::vector<double> e{t0};
    auto span = [&](double len, double step) {
        if (len <= 0) return;
        const int n = std::max(1, int(std::ceil(len / step - 1e-9)));
        const double start = e.back();
        for (int k = 1; k <= n; ++k) e.push_back(start + len * k / n);
    };
    span(P.t_rise, F.dt_ns); span(P.t_flat, F.flat_dt_ns); span(P.t_rise, F.dt_ns); span(F.settle_ns(), F.dt_ns);
    return e;
}

// Programmed flux steps of the two lines for one gate.
struct CZCommand { double amp_mover{0}, amp_coupler{0}; };

inline CZCommand cz_nominal_command(const CZSystem& S, const CZPulse& P) { return {P.flux1_int - S.flux1_idle, P.fluxc_on - S.fluxc_idle}; }

// Amplitude calibration on a distorted line: the programmed step is scaled so that an isolated pulse delivers the
// calibrated flux on average over its flat top (the window is shifted by the delay and the low-pass time constant).
// This is the part of the distortion that a hardware amplitude calibration removes. The DAC resolution is ignored
// here, so that quantisation shows up as an error and is not calibrated away.
inline double flux_amplitude_scale(const FluxChannel& ch, FluxAwg awg, double t_rise, double t_flat) {
    awg.bits = 0;
    const double tg = ch.line.delay_ns + ch.line.lp_tau_ns();
    const std::vector<FluxPulse> unit{{0.0, 1.0, t_rise, t_flat}};
    const double m = t_flat > 0 ? flux_line_output(ch, awg, unit, {t_rise + tg, t_rise + t_flat + tg})[0]
                                : 2 * flux_line_output(ch, awg, unit, {tg, 2 * t_rise + tg})[0];
    return 1.0 / m;
}
inline CZCommand cz_matched_command(const CZSystem& S, const CZPulse& P, const FluxSetup& F) {
    const CZCommand n = cz_nominal_command(S, P);
    return {n.amp_mover * flux_amplitude_scale(F.mover, F.awg, P.t_rise, P.t_flat), n.amp_coupler * flux_amplitude_scale(F.coupler, F.awg, P.t_rise, P.t_flat)};
}

// Drive of the gate with pulse P that starts at t0, when the two lines carry the given pulse trains (which include the
// gate's own pulses). Earlier pulses reach into the gate through the tails of the lines.
inline CZDrive cz_train_drive(const CZSystem& S, const CZPulse& P, const FluxSetup& F, double t0, const std::vector<FluxPulse>& mover_train,
                              const std::vector<FluxPulse>& coupler_train) {
    const auto e = cz_gate_edges(P, F, t0);
    const auto y1 = flux_line_output(F.mover, F.awg, mover_train, e), yc = flux_line_output(F.coupler, F.awg, coupler_train, e);
    CZDrive W;
    for (std::size_t i = 0; i + 1 < e.size(); ++i) {
        W.h.push_back(e[i + 1] - e[i]);
        W.f1.push_back(S.flux1_idle + y1[i]);
        W.fc.push_back(S.fluxc_idle + yc[i]);
    }
    return W;
}
// The same gate on its own, after a long idle.
inline CZDrive cz_gate_drive(const CZSystem& S, const CZPulse& P, const FluxSetup& F, const CZCommand& cmd) {
    return cz_train_drive(S, P, F, 0.0, {{0.0, cmd.amp_mover, P.t_rise, P.t_flat}}, {{0.0, cmd.amp_coupler, P.t_rise, P.t_flat}});
}

// The pair system a stored calibration was made for (design_pair_system's system, with the direct coupling it found).
inline CZSystem cz_system_of(const CZCalibration& c, double alpha_mover, double alpha_target, const CouplerElement& ce) {
    return pair_system(c.f_mover, alpha_mover, c.f_target, alpha_target, ce, c.g12_MHz * 1e-3);
}

// Start times of the CZ gates of a compiled program, in program order: (ns, physical qubits a < b).
struct CZStart { double t_ns{0}; int a{0}, b{0}; };
inline std::vector<CZStart> cz_starts(const Compiled& prog) {
    std::vector<CZStart> out;
    double t = 0;
    for (const auto& m : prog.moments) {
        for (const auto& g : m.ops) if (g.op == Op::CZ) out.push_back({t, std::min(g.q[0], g.q[1]), std::max(g.q[0], g.q[1])});
        t += m.duration_ns;
    }
    return out;
}

} // namespace qpuforge
