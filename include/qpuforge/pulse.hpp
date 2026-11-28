// qpuforge — pulse-level simulation of a tunable-coupler CZ gate
//
// System: qubit 1 (moves) — tunable coupler (transmon) — qubit 2 (parked), all flux-tunable dc-SQUID transmons.
//   H(t) = sum_k E_k(Phi_k(t)) + G_1c n_1 n_c + G_2c n_2 n_c + G_12 n_1 n_2        (charge coupling, lab frame)
// Each element is represented in its own instantaneous eigenbasis (exact charge-basis spectrum, sign-fixed gauge),
// truncated to a few levels. The gate: smooth flux pulses bring qubit 1 to the |11> <-> |20> interaction point while
// the coupler swings from its ZZ-free idle bias towards the qubits, then everything returns. The propagator is built
// from exact exponentials of piecewise-constant H (fine steps on the ramps, one step on the flat top).
// Outputs: average gate fidelity to CZ (with free single-qubit Z corrections), leakage, conditional phase.
#pragma once
#include "device.hpp"
#include <fstream>
#include <mutex>
#include <sstream>

namespace qpuforge {

struct FluxTransmon {
    double EC{0.22}, EJsum{20}, d{0};
    int ncut{12};
    double EJ(double flux) const { return EJ_at_flux(EJsum, d, flux); }
    double f01(double flux) const { return transmon_fa(EJ(flux), EC, 0.0, ncut).f01; }
    double alpha(double flux) const { return transmon_fa(EJ(flux), EC, 0.0, ncut).alpha; }
    double flux_for_f01(double f) const {             // symmetric/asymmetric SQUID: f01 decreases on [0, 0.5)
        if (f >= f01(0)) return 0;
        double lo = 0, hi = 0.4999;
        for (int it = 0; it < 70; ++it) { const double mid = 0.5 * (lo + hi); (f01(mid) > f ? lo : hi) = mid; }
        return 0.5 * (lo + hi);
    }
    static FluxTransmon from_target(double fmax, double alpha, double asym = 0) {
        auto [ej, ec] = solve_EJ_EC(fmax, alpha);
        FluxTransmon t; t.EC = ec; t.EJsum = ej; t.d = asym;
        return t;
    }
};

struct ElementSpectrum { std::vector<double> E, N; int L{0}; double n(int i, int j) const { return N[std::size_t(i) * L + j]; } };

// Levels and charge matrix elements in a smooth sign gauge: psi_0 > 0 at n = 0, then <i|n|i+1> > 0 for every i
// (the harmonic-oscillator convention; these elements never pass through zero, so the gauge is continuous in flux).
inline ElementSpectrum element_spectrum(const FluxTransmon& t, double flux, int L) {
    const int ncut = t.ncut, dim = 2 * ncut + 1;
    std::vector<double> d(dim), e(dim - 1, -t.EJ(flux) / 2), z;
    for (int k = 0; k < dim; ++k) { const double nn = double(k - ncut); d[k] = 4 * t.EC * nn * nn; }
    tridiag_eigen(d, e, z);
    auto Z = [&](int k, int i) -> double& { return z[std::size_t(k) * dim + i]; };
    auto nme = [&](int i, int j) { double v = 0; for (int k = 0; k < dim; ++k) v += Z(k, i) * double(k - ncut) * Z(k, j); return v; };
    if (Z(ncut, 0) < 0) for (int k = 0; k < dim; ++k) Z(k, 0) = -Z(k, 0);
    for (int i = 0; i + 1 < L; ++i)
        if (nme(i, i + 1) < 0) for (int k = 0; k < dim; ++k) Z(k, i + 1) = -Z(k, i + 1);
    ElementSpectrum s; s.L = L;
    for (int i = 0; i < L; ++i) s.E.push_back(d[i] - d[0]);
    s.N.assign(std::size_t(L) * L, 0.0);
    for (int i = 0; i < L; ++i) for (int j = 0; j < L; ++j) s.N[std::size_t(i) * L + j] = nme(i, j);
    return s;
}

struct CZSystem {
    FluxTransmon q1, c, q2;
    double G1c{0}, G2c{0}, G12{0};       // GHz per (n n)
    int L1{4}, Lc{4}, L2{4};             // levels kept per element; Lc = 3 leaves a ~4e-3 rad CZ phase error (see tests)
    double flux1_idle{0}, flux2_idle{0}, fluxc_idle{0.25};
    int dim() const { return L1 * Lc * L2; }
    int index(int i, int k, int j) const { return (i * Lc + k) * L2 + j; }

    std::vector<double> hamiltonian(double f1, double fc, double f2) const {
        const auto s1 = element_spectrum(q1, f1, L1), sc = element_spectrum(c, fc, Lc), s2 = element_spectrum(q2, f2, L2);
        const int D = dim();
        std::vector<double> H(std::size_t(D) * D, 0.0);
        for (int i = 0; i < L1; ++i) for (int k = 0; k < Lc; ++k) for (int j = 0; j < L2; ++j) {
            const int a = index(i, k, j);
            H[std::size_t(a) * D + a] = s1.E[i] + sc.E[k] + s2.E[j];
            for (int i2 = 0; i2 < L1; ++i2) for (int k2 = 0; k2 < Lc; ++k2) for (int j2 = 0; j2 < L2; ++j2) {
                double v = 0;
                if (j == j2) v += G1c * s1.n(i, i2) * sc.n(k, k2);
                if (i == i2) v += G2c * sc.n(k, k2) * s2.n(j, j2);
                if (k == k2) v += G12 * s1.n(i, i2) * s2.n(j, j2);
                if (v != 0) H[std::size_t(a) * D + index(i2, k2, j2)] += v;
            }
        }
        return H;
    }

    // dressed computational states |00>,|01>,|10>,|11> (coupler in 0): energies + eigenvector indices
    struct Dressed { std::array<double, 4> E; std::array<std::vector<double>, 4> vec; double zz; };
    Dressed dressed(double f1, double fc, double f2) const {
        const int D = dim();
        std::vector<double> V;
        const auto ev = sym_eigen_tri(hamiltonian(f1, fc, f2), D, V);
        Dressed out;
        const int bare[4] = {index(0, 0, 0), index(0, 0, 1), index(1, 0, 0), index(1, 0, 1)};
        for (int s = 0; s < 4; ++s) {
            int best = 0; double bw = -1;
            for (int m = 0; m < D; ++m) { const double w = V[std::size_t(bare[s]) * D + m] * V[std::size_t(bare[s]) * D + m]; if (w > bw) { bw = w; best = m; } }
            out.E[std::size_t(s)] = ev[std::size_t(best)];
            out.vec[std::size_t(s)].resize(std::size_t(D));
            for (int r = 0; r < D; ++r) out.vec[std::size_t(s)][std::size_t(r)] = V[std::size_t(r) * D + best];
        }
        out.zz = out.E[3] - out.E[2] - out.E[1] + out.E[0];
        return out;
    }
    double zz(double fc) const { return dressed(flux1_idle, fc, flux2_idle).zz; }

    // Build from target frequencies (GHz) and couplings (GHz) defined at the idle point.
    static CZSystem make(double f1, double a1, double f2, double a2, double fc_max, double ac, double g1c, double g2c, double g12,
                         double fc_idle_guess = 7.3) {
        CZSystem S;
        S.q1 = FluxTransmon::from_target(f1, a1);
        S.q2 = FluxTransmon::from_target(f2, a2);
        S.c = FluxTransmon::from_target(fc_max, ac);
        S.fluxc_idle = S.c.flux_for_f01(fc_idle_guess);
        const double n1 = std::abs(element_spectrum(S.q1, 0, 2).n(0, 1)), n2 = std::abs(element_spectrum(S.q2, 0, 2).n(0, 1));
        const double nc = std::abs(element_spectrum(S.c, S.fluxc_idle, 2).n(0, 1));
        S.G1c = g1c / (n1 * nc); S.G2c = g2c / (n2 * nc); S.G12 = g12 / (n1 * n2);
        return S;
    }

    // Idle coupler bias with the smallest |ZZ|: the zero crossing if one exists between the coupler sweet spot and
    // fc_min_GHz, otherwise the minimum of |ZZ| (qubits detuned by more than their anharmonicity often have no
    // exact zero; the residual minimum is what the chip lives with).
    double find_zz_free_bias(double fc_min_GHz = 6.3) const {
        const double hi_flux = c.flux_for_f01(fc_min_GHz);
        const int Nscan = 160;
        std::vector<double> xs, vs;
        for (int k = 0; k <= Nscan; ++k) { xs.push_back(hi_flux * k / Nscan); vs.push_back(zz(xs.back())); }
        for (int k = 1; k <= Nscan; ++k)
            if ((vs[std::size_t(k)] > 0) != (vs[std::size_t(k - 1)] > 0)) {
                double lo = xs[std::size_t(k - 1)], hi = xs[std::size_t(k)], vlo = vs[std::size_t(k - 1)];
                for (int it = 0; it < 50; ++it) {
                    const double mid = 0.5 * (lo + hi), vm = zz(mid);
                    if ((vm > 0) == (vlo > 0)) { lo = mid; vlo = vm; } else hi = mid;
                }
                return 0.5 * (lo + hi);
            }
        std::size_t kb = 0;
        for (std::size_t k = 1; k < vs.size(); ++k) if (std::abs(vs[k]) < std::abs(vs[kb])) kb = k;
        double a = xs[kb > 0 ? kb - 1 : 0], b = xs[std::min(kb + 1, vs.size() - 1)];
        const double gr = (std::sqrt(5.0) - 1) / 2;
        double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1v = std::abs(zz(x1)), f2v = std::abs(zz(x2));
        for (int it = 0; it < 50; ++it) {
            if (f1v < f2v) { b = x2; x2 = x1; f2v = f1v; x1 = b - gr * (b - a); f1v = std::abs(zz(x1)); }
            else { a = x1; x1 = x2; f1v = f2v; x2 = a + gr * (b - a); f2v = std::abs(zz(x2)); }
        }
        return 0.5 * (a + b);
    }
};

struct CZPulse {
    double flux1_int{0}, fluxc_on{0};   // qubit-1 flux at the interaction point; coupler flux while "on"
    double t_rise{5}, t_flat{20};       // ns (sin^2 ramps, flat top)
    double duration() const { return 2 * t_rise + t_flat; }
};

struct CZResult {
    double fidelity{0}, leakage{0}, cphase{0}, duration_ns{0};
    std::array<std::array<cplx, 4>, 4> M{};   // computational block in the rotating frame (rows: output)
};

inline CZResult simulate_cz(const CZSystem& S, const CZPulse& P, double dt = 0.1) {
    const int D = S.dim();
    const auto d0 = S.dressed(S.flux1_idle, S.fluxc_idle, S.flux2_idle);
    std::array<std::vector<cplx>, 4> psi;
    for (int j = 0; j < 4; ++j) { psi[std::size_t(j)].resize(std::size_t(D)); for (int r = 0; r < D; ++r) psi[std::size_t(j)][std::size_t(r)] = d0.vec[std::size_t(j)][std::size_t(r)]; }
    std::vector<cplx> tmp(static_cast<std::size_t>(D)), work(static_cast<std::size_t>(D));
    auto step = [&](double s, double tau) {
        const double f1 = S.flux1_idle + s * (P.flux1_int - S.flux1_idle);
        const double fc = S.fluxc_idle + s * (P.fluxc_on - S.fluxc_idle);
        std::vector<double> V;
        const auto ev = sym_eigen_tri(S.hamiltonian(f1, fc, S.flux2_idle), D, V);
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
    };
    const int nr = std::max(1, int(std::ceil(P.t_rise / dt)));
    const double h = P.t_rise / nr;
    for (int k = 0; k < nr; ++k) { const double x = (k + 0.5) * h / P.t_rise; step(std::pow(std::sin(PI * x / 2), 2), h); }
    if (P.t_flat > 0) step(1.0, P.t_flat);
    for (int k = 0; k < nr; ++k) { const double x = 1 - (k + 0.5) * h / P.t_rise; step(std::pow(std::sin(PI * x / 2), 2), h); }

    CZResult R;
    R.duration_ns = P.duration();
    // single-qubit rotating frame (dressed qubit frequencies): |11> rotates at E10 + E01 - E00, so the idle ZZ is NOT
    // absorbed into the frame and the conditional phase is the physical one
    const std::array<double, 4> Ef = {d0.E[0], d0.E[1], d0.E[2], d0.E[1] + d0.E[2] - d0.E[0]};
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
        cplx acc = 0;
        for (int r = 0; r < D; ++r) acc += d0.vec[std::size_t(i)][std::size_t(r)] * psi[std::size_t(j)][std::size_t(r)];
        R.M[std::size_t(i)][std::size_t(j)] = acc * std::exp(cplx(0, 2 * PI * Ef[std::size_t(i)] * R.duration_ns));
    }
    // free single-qubit Z corrections: null the phases of |00>, |01>, |10>
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

// The |11> <-> |20> avoided crossing while the coupler is "on": qubit-1 flux at the minimum gap and the gap (GHz).
inline std::pair<double, double> cz_avoided_crossing(const CZSystem& S, double fluxc_on) {
    const int D = S.dim();
    const int b11 = S.index(1, 0, 1), b20 = S.index(2, 0, 0);
    auto gap = [&](double f1) {
        std::vector<double> V;
        const auto ev = sym_eigen_tri(S.hamiltonian(f1, fluxc_on, S.flux2_idle), D, V);
        // the two eigenstates with the largest weight in span{|11>, |20>}
        std::vector<std::pair<double, int>> w;
        for (int m = 0; m < D; ++m) w.push_back({V[std::size_t(b11) * D + m] * V[std::size_t(b11) * D + m] + V[std::size_t(b20) * D + m] * V[std::size_t(b20) * D + m], m});
        std::sort(w.rbegin(), w.rend());
        return std::abs(ev[std::size_t(w[0].second)] - ev[std::size_t(w[1].second)]);
    };
    // bare resonance guess: f1(flux) + alpha1(flux) = f2
    const double f2 = S.q2.f01(S.flux2_idle);
    double lo = 0, hi = 0.45;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (S.q1.f01(mid) + S.q1.alpha(mid) > f2 ? lo : hi) = mid;
    }
    double a = std::max(0.0, 0.5 * (lo + hi) - 0.03), b = std::min(0.49, 0.5 * (lo + hi) + 0.03);
    const double gr = (std::sqrt(5.0) - 1) / 2;
    double x1 = b - gr * (b - a), x2 = a + gr * (b - a), g1 = gap(x1), g2 = gap(x2);
    for (int it = 0; it < 60; ++it) {
        if (g1 < g2) { b = x2; x2 = x1; g2 = g1; x1 = b - gr * (b - a); g1 = gap(x1); }
        else { a = x1; x1 = x2; g1 = g2; x2 = a + gr * (b - a); g2 = gap(x2); }
    }
    const double xm = 0.5 * (a + b);
    return {xm, gap(xm)};
}

// Generic Nelder-Mead minimiser (standard coefficients); returns the best point.
template <class F>
inline std::vector<double> nelder_mead(F&& f, std::vector<double> x0, const std::vector<double>& step, int max_eval, double ftol = 1e-10,
                                       double* fbest = nullptr) {
    const int n = int(x0.size());
    std::vector<std::vector<double>> s(std::size_t(n) + 1, x0);
    std::vector<double> fv(static_cast<std::size_t>(n) + 1);
    for (int k = 0; k < n; ++k) s[std::size_t(k) + 1][std::size_t(k)] += step[std::size_t(k)];
    int evals = 0;
    for (int k = 0; k <= n; ++k) { fv[std::size_t(k)] = f(s[std::size_t(k)]); ++evals; }
    std::vector<int> o(static_cast<std::size_t>(n) + 1);
    while (evals < max_eval) {
        std::iota(o.begin(), o.end(), 0);
        std::sort(o.begin(), o.end(), [&](int a, int b) { return fv[std::size_t(a)] < fv[std::size_t(b)]; });
        const int ib = o.front(), iw = o.back(), is = o[std::size_t(n) - 1];
        if (fv[std::size_t(iw)] - fv[std::size_t(ib)] < ftol) break;
        std::vector<double> cen(static_cast<std::size_t>(n), 0.0);
        for (int k = 0; k <= n; ++k) if (k != iw) for (int d = 0; d < n; ++d) cen[std::size_t(d)] += s[std::size_t(k)][std::size_t(d)] / n;
        auto at = [&](double c) { std::vector<double> x(static_cast<std::size_t>(n)); for (int d = 0; d < n; ++d) x[std::size_t(d)] = cen[std::size_t(d)] + c * (s[std::size_t(iw)][std::size_t(d)] - cen[std::size_t(d)]); return x; };
        auto xr = at(-1); const double fr = f(xr); ++evals;
        if (fr < fv[std::size_t(ib)]) {
            auto xe = at(-2); const double fe = f(xe); ++evals;
            if (fe < fr) { s[std::size_t(iw)] = xe; fv[std::size_t(iw)] = fe; } else { s[std::size_t(iw)] = xr; fv[std::size_t(iw)] = fr; }
        } else if (fr < fv[std::size_t(is)]) { s[std::size_t(iw)] = xr; fv[std::size_t(iw)] = fr; }
        else {
            auto xc = (fr < fv[std::size_t(iw)]) ? at(-0.5) : at(0.5);
            const double fc = f(xc); ++evals;
            if (fc < std::min(fr, fv[std::size_t(iw)])) { s[std::size_t(iw)] = xc; fv[std::size_t(iw)] = fc; }
            else {
                for (int k = 0; k <= n; ++k) if (k != ib) {
                    for (int d = 0; d < n; ++d) s[std::size_t(k)][std::size_t(d)] = 0.5 * (s[std::size_t(k)][std::size_t(d)] + s[std::size_t(ib)][std::size_t(d)]);
                    fv[std::size_t(k)] = f(s[std::size_t(k)]); ++evals;
                }
            }
        }
    }
    int ib = 0;
    for (int k = 1; k <= n; ++k) if (fv[std::size_t(k)] < fv[std::size_t(ib)]) ib = k;
    if (fbest) *fbest = fv[std::size_t(ib)];
    return s[std::size_t(ib)];
}

// Calibrate a CZ: three knobs, as on hardware: qubit-1 interaction flux (|11>-|20> detuning), flat-top length (one full
// |11> -> |20> -> |11> cycle) and coupler "on" flux (coupling strength, which also sets the off-resonant |01>-|10> swap
// so that it closes at the end of the gate). Minimises 1 - F_avg (leakage and conditional-phase error included).
// Leakage fallback: F_avg weighs a leaked population no more than a phase error, and the coarse scan can pick a
// short-flat-top branch that leaks several times more than the neighbouring branch (docs/design/cz-leakage.md). Only if
// the plain calibration fails the leakage limit (leakage > leak_limit) is it repeated with the objective
// 1 - F_avg + leak_weight * leakage. The repeat is kept only if it passes the limit, its 1 - F_avg is no worse than the
// plain one's and stays within coh_limit; otherwise the plain calibration stands. Pairs within the limit are
// calibrated exactly as before; leak_weight = 0 disables the fallback.
inline std::pair<CZPulse, CZResult> optimize_cz(const CZSystem& S, double fluxc_on, double t_rise, double dt = 0.1, int max_eval = 260,
                                                double leak_weight = 0, double leak_limit = 1.0, double coh_limit = 1.0) {
    // named variables, not a structured binding: x0 is captured by the lambda below (capturing bindings is C++20)
    const auto crossing = cz_avoided_crossing(S, fluxc_on);
    const double x0 = crossing.first, g0 = crossing.second;
    CZPulse P; P.fluxc_on = fluxc_on; P.t_rise = t_rise;
    const double t_full = 1.0 / std::max(g0, 1e-6);
    auto run = [&](double fx, double tf, double fc) { CZPulse q = P; q.flux1_int = fx; q.t_flat = std::max(0.0, tf); q.fluxc_on = fc; return simulate_cz(S, q, dt); };
    auto calibrate = [&](double w) {
        auto objective = [&](const CZResult& r) { return 1 - r.fidelity + w * r.leakage; };
        double bx = x0, bt = std::max(0.0, t_full - t_rise), bc = 1e9;
        for (int b = -5; b <= 5; ++b) {
            const double tf = std::max(0.0, t_full - t_rise + b * 0.1 * t_full);
            for (int a = -2; a <= 2; ++a) {
                const double fx = x0 + a * 0.0010;
                const double c = objective(run(fx, tf, fluxc_on));
                if (c < bc) { bc = c; bx = fx; bt = tf; }
            }
        }
        auto cost = [&](const std::vector<double>& x) { return objective(run(x[0], x[1], x[2])); };
        const auto best = nelder_mead(cost, {bx, bt, fluxc_on}, {0.0008, 0.05 * t_full + 0.3, 0.002}, max_eval, 1e-10);
        CZPulse q = P; q.flux1_int = best[0]; q.t_flat = std::max(0.0, best[1]); q.fluxc_on = best[2];
        return std::make_pair(q, simulate_cz(S, q, dt));
    };
    auto out = calibrate(0.0);
    if (leak_weight > 0 && out.second.leakage > leak_limit) {
        auto alt = calibrate(leak_weight);
        const double e_alt = 1 - alt.second.fidelity, e_out = 1 - out.second.fidelity;
        if (alt.second.leakage <= leak_limit && e_alt <= e_out && e_alt <= coh_limit) out = alt;
    }
    return out;
}

// ---------------------------------------------------------------- chip-level calibration

// The pair system for a given direct coupling g12 (GHz), coupler parked at its sweet spot.
inline CZSystem pair_system(double f_hi, double a_hi, double f_lo, double a_lo, const CouplerElement& ce, double g12) {
    const double g = ce.g_qc_MHz * 1e-3;
    auto S = CZSystem::make(f_hi, a_hi, f_lo, a_lo, ce.f_max_GHz, ce.alpha_GHz, g, g, g12, ce.f_max_GHz);
    S.fluxc_idle = 0.0;
    return S;
}
// Pair system for a physical coupler parked at its sweet spot. The direct coupling g12 (set per pair by layout) is
// chosen so that the residual |11>-|20> coupling with the coupler parked vanishes: the qubits are then exchange-free
// at idle, and a measure qubit sweeping past a spectator's |11>-|20> point while its other couplers are parked sees
// (almost) no avoided crossing.
inline CZSystem design_pair_system(double f_hi, double a_hi, double f_lo, double a_lo, const CouplerElement& ce, double* g12_out = nullptr,
                                   double* gap_out = nullptr) {
    auto build = [&](double g12) { return pair_system(f_hi, a_hi, f_lo, a_lo, ce, g12); };
    auto gap = [&](double g12) { return cz_avoided_crossing(build(g12), 0.0).second; };
    // coarse scan then golden section on the V-shaped |gap|(g12)
    double best = 0, bg = 1e9;
    for (int k = 0; k <= 24; ++k) { const double x = 0.0015 * k, v = gap(x); if (v < bg) { bg = v; best = x; } }
    double a = std::max(0.0, best - 0.0015), b = best + 0.0015;
    const double gr = (std::sqrt(5.0) - 1) / 2;
    double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = gap(x1), f2 = gap(x2);
    for (int it = 0; it < 30; ++it) {
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = gap(x1); }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = gap(x2); }
    }
    const double g12 = 0.5 * (a + b);
    if (g12_out) *g12_out = g12;
    if (gap_out) *gap_out = gap(g12);
    return build(g12);
}

// Numerical model behind a coherent CZ calibration. Everything here changes the cached result, so all of it is part of
// the cache key; bump kVersion whenever the CZ Hamiltonian, pulse shape, calibration procedure or fidelity metric
// changes, so stale entries miss instead of being reused silently.
struct CZModel {
    static constexpr int kVersion = 4;        // v1 = qpuforge 0.3 pre-release (key had no model fields); v4 adds the leakage fallback
                                              // (v3, never released: triggered at 1e-4 and kept the repeat on leakage alone)
    static constexpr double dt_ns = 0.1;      // ramp time step of simulate_cz during calibration
    static constexpr int max_eval = 260;      // Nelder-Mead budget of optimize_cz
    static constexpr double leak_limit = 2e-4; // Alpha gate: leakage per CZ; optimize_cz re-calibrates a pair above it
    static constexpr double coh_limit = 5e-4;  // Alpha gate: coherent error per CZ; the repeat must stay within it
    static constexpr double leak_weight = 10;  // objective of the repeat: 1 - F_avg + leak_weight * leakage
};
struct CZCalibration {
    int a{0}, b{0}, mover{0}, target{0};
    double f_mover{0}, f_target{0}, g12_MHz{0}, park_fc_GHz{0}, zz_idle_kHz{0}, residual_gap_MHz{0};
    double f_int_GHz{0}, excursion_MHz{0}, fc_on_GHz{0}, t_rise_ns{0}, t_flat_ns{0}, t_cz_ns{0};
    double err_coherent{0}, leakage{0}, cphase_err{0};
    double err_flux_local{0}, err_flux_quasistatic{0}, leak_spectators{0}, err_residual{0};
    bool distortion_modeled{false};   // ChipDesign::flux_line.enabled (q18): err_distortion is part of err_ctrl
    double flux_eps_eff{0};           // effective fractional flux-step error over the hold (|static + averaged tails|)
    double err_distortion{0};         // flux-line distortion error, on top of err_residual (0 unless modeled)
    double err_ctrl{0};           // Pauli-type error handed to the emulators (T1/Tphi during the gate are simulated there)
    double err_total{0};          // average gate infidelity incl. decoherence during the gate, for reporting
    std::vector<std::pair<int, double>> spectators;   // (spectator qubit id, |11>-|20> transfer probability per CZ)
    CZPulse pulse;
};

// Coherent part of one pair's calibration (pure function of frequencies + coupler design).
inline CZCalibration calibrate_cz_pair(double f_hi, double a_hi, double f_lo, double a_lo, const CouplerElement& ce) {
    CZCalibration c;
    double g12 = 0, gap = 0;
    const CZSystem S = design_pair_system(f_hi, a_hi, f_lo, a_lo, ce, &g12, &gap);
    c.f_mover = f_hi; c.f_target = f_lo;
    c.g12_MHz = g12 * 1e3; c.residual_gap_MHz = gap * 1e3;
    c.park_fc_GHz = S.c.f01(S.fluxc_idle);
    c.zz_idle_kHz = S.zz(S.fluxc_idle) * 1e6;
    // named variables, not a structured binding: P is captured by the lambda below (capturing bindings is C++20)
    const auto opt = optimize_cz(S, S.c.flux_for_f01(ce.f_on_GHz), ce.t_rise_ns, CZModel::dt_ns, CZModel::max_eval, CZModel::leak_weight, CZModel::leak_limit, CZModel::coh_limit);
    const CZPulse P = opt.first;
    const CZResult R = opt.second;
    c.pulse = P;
    c.fc_on_GHz = S.c.f01(P.fluxc_on);
    c.f_int_GHz = S.q1.f01(P.flux1_int);
    c.excursion_MHz = (f_hi - c.f_int_GHz) * 1e3;
    c.t_rise_ns = P.t_rise; c.t_flat_ns = P.t_flat; c.t_cz_ns = P.duration();
    c.err_coherent = 1 - R.fidelity; c.leakage = R.leakage; c.cphase_err = R.cphase - PI;
    // quasi-static flux noise on the mover and the coupler: curvature of 1 - F (local Z phases are corrected away in F,
    // so this captures conditional-phase and leakage sensitivity)
    const double h = 2e-4;
    auto inf = [&](double d1, double dc) { CZPulse q = P; q.flux1_int += d1; q.fluxc_on += dc; return 1 - simulate_cz(S, q, CZModel::dt_ns).fidelity; };
    const double e0 = c.err_coherent;
    const double c1 = (inf(h, 0) - 2 * e0 + inf(-h, 0)) / (h * h), cc = (inf(0, h) - 2 * e0 + inf(0, -h)) / (h * h);
    c.err_flux_quasistatic = std::max(0.0, c1) + std::max(0.0, cc);   // times sigma^2 / 2 later (needs the chip's A)
    return c;
}

// 1/f flux-noise variance (Phi0^2) seen by a gate, with the same infrared cut-off as the coherence model.
inline double flux_noise_variance(double A_uPhi0) {
    const double A = A_uPhi0 * 1e-6;
    const double ln_ir = std::abs(std::log(2 * PI * 1.0 * 10e-6));
    return 2 * A * A * ln_ir;
}

// q18 fast path (docs/design/flux-line-distortion.md): the residual step error of the flux lines, averaged over the
// hold of the CZ. A tail a*exp(-t/tau) averages to a*(tau/T)*(1 - exp(-T/tau)) over a hold of T ns; the static offset
// is taken as is. Both lines (mover and coupler) are given the same fractional error.
inline double flux_line_eps_eff(const FluxLineSpec& L, double t_hold_ns) {
    double e = L.static_offset;
    const double T = std::max(t_hold_ns, 1e-9);
    for (const auto& t : L.tails) e += t.amplitude * (t.tau_ns / T) * (1 - std::exp(-T / t.tau_ns));
    return std::abs(e);
}
// Infidelity of a fractional error eps on both flux steps, from the stored quasi-static curvature of 1 - F in flux
// (curv = c_mover + c_coupler, per Phi0^2): 1/2 (c_m (eps dPhi_m)^2 + c_c (eps dPhi_c)^2) <= 1/2 curv (eps max|dPhi|)^2.
// The bound is used because the cache stores only the sum; it is the "worst pair" column of the q18 sensitivity table.
// `curv` is the raw curvature (CZCalibration::err_flux_quasistatic before calibrate_chip_cz scales it by sigma^2 / 2).
// Both idle fluxes are 0 (mover at its sweet spot, coupler parked at f_max), so the steps are flux1_int and fluxc_on.
inline double flux_distortion_err(double curv, const CZPulse& P, double eps) {
    const double dphi = std::max(std::abs(P.flux1_int), std::abs(P.fluxc_on));
    return 0.5 * std::max(0.0, curv) * (eps * dphi) * (eps * dphi);
}

// Coherent calibrations depend only on the pair's frequencies/anharmonicities and the coupler design: memoised in
// process and optionally in a JSON file so reviews and parameter sweeps do not re-run the pulse optimisation.
inline std::string cz_cache_key(double fh, double ah, double fl, double al, const CouplerElement& ce) {
    const CZSystem ref;                       // level truncations / charge cutoff actually used by design_pair_system
    char b[320];
    std::snprintf(b, sizeof b, "v%d|dt%.4f_ev%d_lw%.3g_ll%.3g_cl%.3g_L%d.%d.%d_nc%d.%d.%d|%.5f_%.5f_%.5f_%.5f|%.4f_%.4f_%.3f_%.4f_%.3f", CZModel::kVersion,
                  CZModel::dt_ns, CZModel::max_eval, CZModel::leak_weight, CZModel::leak_limit, CZModel::coh_limit, ref.L1, ref.Lc, ref.L2, ref.q1.ncut, ref.c.ncut, ref.q2.ncut, fh, ah, fl, al,
                  ce.f_max_GHz, ce.alpha_GHz, ce.g_qc_MHz, ce.f_on_GHz, ce.t_rise_ns);
    return b;
}
inline Json cz_coherent_json(const CZCalibration& c) {
    Json o = Json::object();
    o.set("f_mover", c.f_mover); o.set("f_target", c.f_target); o.set("g12_MHz", c.g12_MHz); o.set("residual_gap_MHz", c.residual_gap_MHz);
    o.set("park_fc_GHz", c.park_fc_GHz); o.set("zz_idle_kHz", c.zz_idle_kHz); o.set("fc_on_GHz", c.fc_on_GHz); o.set("f_int_GHz", c.f_int_GHz);
    o.set("excursion_MHz", c.excursion_MHz); o.set("t_rise_ns", c.t_rise_ns); o.set("t_flat_ns", c.t_flat_ns); o.set("t_cz_ns", c.t_cz_ns);
    o.set("err_coherent", c.err_coherent); o.set("leakage", c.leakage); o.set("cphase_err", c.cphase_err);
    o.set("curv_flux", c.err_flux_quasistatic);
    o.set("flux1_int", c.pulse.flux1_int); o.set("fluxc_on", c.pulse.fluxc_on); o.set("p_t_rise", c.pulse.t_rise); o.set("p_t_flat", c.pulse.t_flat);
    return o;
}
inline CZCalibration cz_coherent_from_json(const Json& o) {
    CZCalibration c;
    c.f_mover = o.get("f_mover", 0.0); c.f_target = o.get("f_target", 0.0); c.g12_MHz = o.get("g12_MHz", 0.0);
    c.residual_gap_MHz = o.get("residual_gap_MHz", 0.0); c.park_fc_GHz = o.get("park_fc_GHz", 0.0); c.zz_idle_kHz = o.get("zz_idle_kHz", 0.0);
    c.fc_on_GHz = o.get("fc_on_GHz", 0.0); c.f_int_GHz = o.get("f_int_GHz", 0.0); c.excursion_MHz = o.get("excursion_MHz", 0.0);
    c.t_rise_ns = o.get("t_rise_ns", 0.0); c.t_flat_ns = o.get("t_flat_ns", 0.0); c.t_cz_ns = o.get("t_cz_ns", 0.0);
    c.err_coherent = o.get("err_coherent", 0.0); c.leakage = o.get("leakage", 0.0); c.cphase_err = o.get("cphase_err", 0.0);
    c.err_flux_quasistatic = o.get("curv_flux", 0.0);
    c.pulse.flux1_int = o.get("flux1_int", 0.0); c.pulse.fluxc_on = o.get("fluxc_on", 0.0);
    c.pulse.t_rise = o.get("p_t_rise", 0.0); c.pulse.t_flat = o.get("p_t_flat", 0.0);
    return c;
}
struct CZCache {
    std::map<std::string, CZCalibration> m;
    std::mutex mu;
    static CZCache& instance() { static CZCache c; return c; }
    void load(const std::string& path) {
        if (path.empty()) return;
        std::ifstream f(path);
        if (!f) return;
        std::stringstream ss; ss << f.rdbuf();
        const Json j = Json::parse(ss.str());
        // versioned file: {"schema": "qpuforge-cz-cache", "model_version": N, "entries": {key: calibration}}.
        // Anything else (including the unversioned v0.3 pre-release format) is rejected and regenerated.
        const int ver = int(j.get("model_version", -1.0));
        if (j.get("schema", "") != std::string("qpuforge-cz-cache") || ver != CZModel::kVersion || !j["entries"].is_obj()) {
            std::fprintf(stderr, "CZ cache %s: model version %d does not match %d, ignored (calibrations will be re-simulated)\n",
                         path.c_str(), ver, CZModel::kVersion);
            return;
        }
        std::lock_guard<std::mutex> lk(mu);
        for (const auto& [k, v] : j["entries"].items()) if (!m.count(k)) m[k] = cz_coherent_from_json(v);
    }
    void save(const std::string& path) {
        if (path.empty()) return;
        std::lock_guard<std::mutex> lk(mu);
        Json e = Json::object();
        for (const auto& [k, v] : m) e.set(k, cz_coherent_json(v));
        Json j = Json::object();
        j.set("schema", "qpuforge-cz-cache"); j.set("model_version", CZModel::kVersion); j.set("entries", e);
        std::ofstream f(path); j.dump(f); f << "\n";
    }
};

// Calibrate every coupled pair of a chip with a physical coupler. Unique frequency combinations are simulated once
// (in parallel); results are written into the report (t_cz, ZZ at park, control error, leakage) and returned.
inline std::vector<CZCalibration> calibrate_chip_cz(const ChipDesign& D, ChipReport& R, bool verbose = false, const std::string& cache_path = "") {
    auto& cache = CZCache::instance();
    cache.load(cache_path);
    std::vector<int> pidx;
    for (std::size_t k = 0; k < R.pairs.size(); ++k) if (R.pairs[k].designed) pidx.push_back(int(k));
    struct Key { long fh, ah, fl, al; bool operator<(const Key& o) const { return std::tie(fh, ah, fl, al) < std::tie(o.fh, o.ah, o.fl, o.al); } };
    auto key_of = [&](const PairReport& p) {
        const auto& qa = R.q[std::size_t(D.index_of(p.a))]; const auto& qb = R.q[std::size_t(D.index_of(p.b))];
        const auto& hi = qa.f01 >= qb.f01 ? qa : qb; const auto& lo = qa.f01 >= qb.f01 ? qb : qa;
        return Key{std::lround(hi.f01 * 1e5), std::lround(hi.alpha * 1e5), std::lround(lo.f01 * 1e5), std::lround(lo.alpha * 1e5)};
    };
    std::map<Key, int> uniq; std::vector<Key> keys;
    for (int k : pidx) { const Key kk = key_of(R.pairs[std::size_t(k)]); if (!uniq.count(kk)) { uniq[kk] = int(keys.size()); keys.push_back(kk); } }
    std::vector<CZCalibration> base(keys.size());
    std::vector<int> todo;
    for (int u = 0; u < int(keys.size()); ++u) {
        const Key& k = keys[std::size_t(u)];
        const auto ck = cz_cache_key(k.fh * 1e-5, k.ah * 1e-5, k.fl * 1e-5, k.al * 1e-5, D.coupler);
        std::lock_guard<std::mutex> lk(cache.mu);
        auto it = cache.m.find(ck);
        if (it != cache.m.end()) base[std::size_t(u)] = it->second; else todo.push_back(u);
    }
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < int(todo.size()); ++i) {
        const int u = todo[std::size_t(i)];
        const Key& k = keys[std::size_t(u)];
        base[std::size_t(u)] = calibrate_cz_pair(k.fh * 1e-5, k.ah * 1e-5, k.fl * 1e-5, k.al * 1e-5, D.coupler);
        const auto ck = cz_cache_key(k.fh * 1e-5, k.ah * 1e-5, k.fl * 1e-5, k.al * 1e-5, D.coupler);
        std::lock_guard<std::mutex> lk(cache.mu);
        cache.m[ck] = base[std::size_t(u)];
    }
    if (!todo.empty()) cache.save(cache_path);
    if (verbose) std::fprintf(stderr, "CZ calibration: %zu unique frequency combinations for %zu pairs (%zu simulated, %zu cached)\n", keys.size(),
                              pidx.size(), todo.size(), keys.size() - todo.size());
    const double sig2 = flux_noise_variance(D.flux_noise_uPhi0);
    std::vector<CZCalibration> out;
    std::map<std::pair<int, int>, const CZCalibration*> by_pair;
    std::vector<CZCalibration> tmp;
    tmp.reserve(pidx.size());
    for (int k : pidx) {
        auto& p = R.pairs[std::size_t(k)];
        CZCalibration c = base[std::size_t(uniq[key_of(p)])];
        const auto& qa = R.q[std::size_t(D.index_of(p.a))]; const auto& qb = R.q[std::size_t(D.index_of(p.b))];
        c.a = p.a; c.b = p.b;
        c.mover = qa.f01 >= qb.f01 ? p.a : p.b; c.target = c.mover == p.a ? p.b : p.a;
        tmp.push_back(c);
    }
    for (auto& c : tmp) by_pair[{std::min(c.a, c.b), std::max(c.a, c.b)}] = &c;
    for (auto& c : tmp) {
        const auto& qm = R.q[std::size_t(D.index_of(c.mover))];
        const auto& qt = R.q[std::size_t(D.index_of(c.target))];
        // mover's local dephasing at the excursion (Gaussian 1/f decay over the time spent away from the sweet spot)
        const FluxTransmon tm = FluxTransmon::from_target(c.f_mover, qm.alpha);
        const double hh = 1e-4;
        const double dfdphi = (tm.f01(c.pulse.flux1_int + hh) - tm.f01(c.pulse.flux1_int - hh)) / (2 * hh);   // GHz/Phi0
        const double t_eff = c.t_flat_ns + c.t_rise_ns;
        const double var_theta = std::pow(2 * PI * dfdphi * t_eff, 2) * sig2;     // rad^2 (GHz*ns = cycles)
        const double pz = 0.5 * (1 - std::exp(-0.5 * var_theta));
        c.err_flux_local = 0.8 * pz;
        if (D.flux_line.enabled) {   // q18, opt-in: from the raw curvature, before it is scaled to the noise variance
            c.distortion_modeled = true;
            c.flux_eps_eff = flux_line_eps_eff(D.flux_line, c.t_flat_ns);
            c.err_distortion = flux_distortion_err(c.err_flux_quasistatic, c.pulse, c.flux_eps_eff);
        }
        c.err_flux_quasistatic *= 0.5 * sig2;
        // spectator crossings: the mover sweeps through the |11>-|20> point of each other neighbour s with
        // f_s - alpha_m between f_int and f_idle; Landau-Zener transfer with that pair's parked residual gap, twice per CZ
        c.leak_spectators = 0; c.spectators.clear();
        for (const auto& q2 : R.pairs) {
            if (!q2.designed) continue;
            int s = -1;
            if (q2.a == c.mover && q2.b != c.target) s = q2.b;
            if (q2.b == c.mover && q2.a != c.target) s = q2.a;
            if (s < 0) continue;
            const auto& qs = R.q[std::size_t(D.index_of(s))];
            const double fx = qs.f01 - qm.alpha;
            if (!(fx > c.f_int_GHz && fx < c.f_mover)) continue;
            auto it = by_pair.find({std::min(c.mover, s), std::max(c.mover, s)});
            if (it == by_pair.end()) continue;
            const double gapx = it->second->residual_gap_MHz * 1e-3;
            const double phx = tm.flux_for_f01(fx);
            const double sx = std::clamp(phx / std::max(c.pulse.flux1_int, 1e-9), 1e-6, 1 - 1e-6);
            const double dsdt = PI / c.t_rise_ns * std::sqrt(sx * (1 - sx));
            const double dfdphx = std::abs(tm.f01(phx + hh) - tm.f01(phx - hh)) / (2 * hh);
            const double v = dfdphx * c.pulse.flux1_int * dsdt;    // GHz/ns
            const double P = std::min(1.0, 2 * PI * PI * gapx * gapx / std::max(v, 1e-9));
            c.spectators.push_back({s, P});
            c.leak_spectators += P;
        }
        c.err_residual = D.coupler.residual_ctrl_err;
        c.err_ctrl = std::max(0.0, c.err_coherent - c.leakage) + c.err_flux_local + c.err_flux_quasistatic + c.err_residual;
        if (c.distortion_modeled) c.err_ctrl += c.err_distortion;   // additive, after the unchanged sum (off: same bits)
        auto rates = [&](const QubitReport& q) { return 1e6 / q.T1_us + (std::isfinite(q.Tphi_us) ? 1e6 / q.Tphi_us : 0.0); };
        c.err_total = c.err_ctrl + c.leakage + c.leak_spectators + 0.8 * (c.t_cz_ns * 1e-9 / 2) * (rates(qm) + rates(qt));
    }
    // write back
    double sum = 0; int n = 0;
    for (const auto& c : tmp) {
        for (auto& p : R.pairs) {
            if (!((p.a == c.a && p.b == c.b))) continue;
            p.t_cz_ns = c.t_cz_ns;
            p.zz_idle_kHz = c.zz_idle_kHz;
            p.err_2q_ctrl = c.err_ctrl;
            p.leak_2q = c.leakage + c.leak_spectators;
            p.leak_qubit = c.mover;
            p.err_2q = c.err_total;
            p.cz_interaction_GHz = c.f_int_GHz;
            p.cz_excursion_MHz = c.excursion_MHz;
            p.cz_calibrated = true;
            sum += c.err_total; ++n;
        }
    }
    if (n) R.mean_err_2q = sum / n;
    return tmp;
}

inline Json cz_calibration_json(const std::vector<CZCalibration>& cal) {
    Json arr = Json::array();
    for (const auto& c : cal) {
        Json o = Json::object();
        o.set("a", c.a); o.set("b", c.b); o.set("mover", c.mover);
        o.set("f_mover_GHz", c.f_mover); o.set("f_target_GHz", c.f_target); o.set("g12_MHz", c.g12_MHz);
        o.set("park_fc_GHz", c.park_fc_GHz); o.set("zz_idle_kHz", c.zz_idle_kHz); o.set("residual_gap_MHz", c.residual_gap_MHz);
        o.set("f_int_GHz", c.f_int_GHz); o.set("excursion_MHz", c.excursion_MHz); o.set("fc_on_GHz", c.fc_on_GHz);
        o.set("t_rise_ns", c.t_rise_ns); o.set("t_flat_ns", c.t_flat_ns); o.set("t_cz_ns", c.t_cz_ns);
        o.set("flux_mover_int", c.pulse.flux1_int); o.set("flux_coupler_on", c.pulse.fluxc_on);
        o.set("err_coherent", c.err_coherent); o.set("leakage", c.leakage); o.set("cphase_err_rad", c.cphase_err);
        o.set("err_flux_local", c.err_flux_local); o.set("err_flux_quasistatic", c.err_flux_quasistatic);
        o.set("leak_spectators", c.leak_spectators); o.set("err_residual", c.err_residual);
        if (c.distortion_modeled) { o.set("flux_eps_eff", c.flux_eps_eff); o.set("err_distortion", c.err_distortion); }   // q18, only when on
        o.set("err_ctrl", c.err_ctrl); o.set("err_total", c.err_total);
        arr.push(o);
    }
    return arr;
}

} // namespace qpuforge
