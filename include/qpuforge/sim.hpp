// qpuforge — simulation backends
//   StateVector   : pure state, noise applied by quantum-trajectory (Monte Carlo wave-function) unravelling
//   DensityMatrix : exact mixed-state evolution, used to validate the trajectory backend (small n)
// Qubit k <-> bit k of the basis index (little-endian).
#pragma once
#include "core.hpp"

namespace qpuforge {

class StateVector {
public:
    explicit StateVector(int n) : n_(n), a_(std::size_t(1) << n) { reset(); }
    int n() const { return n_; }
    std::size_t dim() const { return a_.size(); }
    const std::vector<cplx>& amps() const { return a_; }
    std::vector<cplx>& amps() { return a_; }

    void reset() { std::fill(a_.begin(), a_.end(), cplx(0)); a_[0] = 1; }
    void set_basis(std::uint64_t k) { std::fill(a_.begin(), a_.end(), cplx(0)); a_[k] = 1; }

    void gate1(int q, const M2& u) {
        const std::size_t m = std::size_t(1) << q, D = a_.size();
        for (std::size_t hi = 0; hi < D; hi += 2 * m)
            for (std::size_t i = hi; i < hi + m; ++i) {
                const cplx x = a_[i], y = a_[i | m];
                a_[i] = u.a * x + u.b * y;
                a_[i | m] = u.c * x + u.d * y;
            }
    }
    void cz(int p, int q) {
        const std::size_t mask = (std::size_t(1) << p) | (std::size_t(1) << q);
        for (std::size_t i = 0; i < a_.size(); ++i) if ((i & mask) == mask) a_[i] = -a_[i];
    }
    void cphase(int p, int q, double phi) {
        const std::size_t mask = (std::size_t(1) << p) | (std::size_t(1) << q);
        const cplx ph = std::exp(cplx(0, phi));
        for (std::size_t i = 0; i < a_.size(); ++i) if ((i & mask) == mask) a_[i] *= ph;
    }
    void cx(int c, int t) {
        const std::size_t mc = std::size_t(1) << c, mt = std::size_t(1) << t;
        for (std::size_t i = 0; i < a_.size(); ++i) if ((i & mc) && !(i & mt)) std::swap(a_[i], a_[i | mt]);
    }
    void swapq(int p, int q) {
        const std::size_t mp = std::size_t(1) << p, mq = std::size_t(1) << q;
        for (std::size_t i = 0; i < a_.size(); ++i) if ((i & mp) && !(i & mq)) std::swap(a_[i], a_[(i & ~mp) | mq]);
    }
    // U acts on (|a=1,b=0>, |a=0,b=1>); identity on |00>, |11>
    void exch(int a, int b, const M2& u) {
        const std::size_t ma = std::size_t(1) << a, mb = std::size_t(1) << b;
        for (std::size_t i = 0; i < a_.size(); ++i) if ((i & ma) && !(i & mb)) {
            const std::size_t j = i ^ (ma | mb);
            const cplx x = a_[i], y = a_[j];
            a_[i] = u.a * x + u.b * y;
            a_[j] = u.c * x + u.d * y;
        }
    }
    double prob1(int q) const {
        const std::size_t m = std::size_t(1) << q; double s = 0;
        for (std::size_t i = 0; i < a_.size(); ++i) if (i & m) s += std::norm(a_[i]);
        return s;
    }
    void scale_bit(int q, double s0, double s1) {
        const std::size_t m = std::size_t(1) << q;
        for (std::size_t i = 0; i < a_.size(); ++i) a_[i] *= (i & m) ? s1 : s0;
    }
    void lower(int q) {   // sigma^- : |1> -> |0>
        const std::size_t m = std::size_t(1) << q, D = a_.size();
        for (std::size_t hi = 0; hi < D; hi += 2 * m) for (std::size_t i = hi; i < hi + m; ++i) { a_[i] = a_[i | m]; a_[i | m] = 0; }
    }
    void raise(int q) {   // sigma^+ : |0> -> |1>
        const std::size_t m = std::size_t(1) << q, D = a_.size();
        for (std::size_t hi = 0; hi < D; hi += 2 * m) for (std::size_t i = hi; i < hi + m; ++i) { a_[i | m] = a_[i]; a_[i] = 0; }
    }
    void normalize() {
        double s = 0; for (auto& x : a_) s += std::norm(x);
        const double f = 1 / std::sqrt(s); for (auto& x : a_) x *= f;
    }
    cplx inner(const StateVector& o) const {   // <this|o>
        cplx s = 0; for (std::size_t i = 0; i < a_.size(); ++i) s += std::conj(a_[i]) * o.a_[i];
        return s;
    }

    // projective Z measurement: sample, collapse, return outcome
    template <class RNG> int measure(int q, RNG& rng) {
        const double P1 = prob1(q);
        const int out = std::uniform_real_distribution<double>(0, 1)(rng) < P1 ? 1 : 0;
        const std::size_t m = std::size_t(1) << q;
        for (std::size_t i = 0; i < a_.size(); ++i) if ((((i & m) != 0) ? 1 : 0) != out) a_[i] = 0;
        normalize();
        return out;
    }
    template <class RNG> void reset(int q, RNG& rng) { if (measure(q, rng)) gate1(q, gates::X()); }
    template <class RNG> void bitflip(int q, double p, RNG& rng) {
        if (p > 0 && std::uniform_real_distribution<double>(0, 1)(rng) < p) gate1(q, gates::X());
    }
    // many conditional phases in one pass: amplitude *= prod over (mask fully set) of phase
    void multi_cphase(const std::vector<std::pair<std::size_t, cplx>>& terms) {
        if (terms.empty()) return;
        for (std::size_t i = 0; i < a_.size(); ++i) {
            cplx f = 1; bool hit = false;
            for (const auto& [mask, ph] : terms) if ((i & mask) == mask) { f *= ph; hit = true; }
            if (hit) a_[i] *= f;
        }
    }

    // ---------------- stochastic channels (trajectory unravelling; exact in expectation)
    // Generalised amplitude damping: relaxation towards a thermal state.
    // gamma = 1 - exp(-t/T1);  p_down = Gamma_down / (Gamma_down + Gamma_up)
    template <class RNG> void gad(int q, double gamma, double p_down, RNG& rng) {
        if (gamma <= 0) return;
        const double P1 = prob1(q), P0 = 1 - P1;
        const double w_dn = p_down * gamma * P1, w_up = (1 - p_down) * gamma * P0;
        const double w_k0 = p_down * (P0 + (1 - gamma) * P1);
        const double r = std::uniform_real_distribution<double>(0, 1)(rng);
        if (r < w_dn) { lower(q); normalize(); }
        else if (r < w_dn + w_up) { raise(q); normalize(); }
        else if (r < w_dn + w_up + w_k0) {     // no-jump branches: the new norm is known analytically, one pass
            const double nn = std::sqrt(P0 + (1 - gamma) * P1);
            scale_bit(q, 1.0 / nn, std::sqrt(1 - gamma) / nn);
        } else {
            const double nn = std::sqrt((1 - gamma) * P0 + P1);
            scale_bit(q, std::sqrt(1 - gamma) / nn, 1.0 / nn);
        }
    }
    template <class RNG> void dephase(int q, double pz, RNG& rng) {
        if (pz > 0 && std::uniform_real_distribution<double>(0, 1)(rng) < pz) gate1(q, gates::Z());
    }
    template <class RNG> void pauli1(int q, double p, RNG& rng) {
        std::uniform_real_distribution<double> U(0, 1);
        if (p <= 0 || U(rng) >= p) return;
        const int k = int(U(rng) * 3);
        gate1(q, k == 0 ? gates::X() : (k == 1 ? gates::Y() : gates::Z()));
    }
    template <class RNG> void pauli2(int a, int b, double p, RNG& rng) {
        std::uniform_real_distribution<double> U(0, 1);
        if (p <= 0 || U(rng) >= p) return;
        const int k = 1 + int(U(rng) * 15);        // 1..15, never II
        const M2 P[4] = {gates::I(), gates::X(), gates::Y(), gates::Z()};
        if (k % 4) gate1(a, P[k % 4]);
        if (k / 4) gate1(b, P[k / 4]);
    }

private:
    int n_;
    std::vector<cplx> a_;
};

class DensityMatrix {
public:
    explicit DensityMatrix(int n) : n_(n), D_(std::size_t(1) << n), r_(D_ * D_, 0.0) { r_[0] = 1; }
    int n() const { return n_; }
    std::size_t dim() const { return D_; }
    cplx& at(std::size_t i, std::size_t j) { return r_[i * D_ + j]; }
    cplx at(std::size_t i, std::size_t j) const { return r_[i * D_ + j]; }

    void set_product_thermal(const std::vector<double>& p1) {   // diagonal product state
        std::fill(r_.begin(), r_.end(), cplx(0));
        for (std::size_t i = 0; i < D_; ++i) {
            double p = 1; for (int q = 0; q < n_; ++q) p *= ((i >> q) & 1) ? p1[q] : 1 - p1[q];
            r_[i * D_ + i] = p;
        }
    }
    void gate1(int q, const M2& u) { left(q, u); right_dagger(q, u); }
    void cz(int p, int q) { diag2(p, q, cplx(-1)); }
    void exch(int a, int b, const M2& u) {
        const std::size_t ma = std::size_t(1) << a, mb = std::size_t(1) << b, mm = ma | mb;
        for (std::size_t i = 0; i < D_; ++i) if ((i & ma) && !(i & mb)) {          // rows
            const std::size_t k = i ^ mm;
            for (std::size_t j = 0; j < D_; ++j) {
                const cplx x = r_[i * D_ + j], y = r_[k * D_ + j];
                r_[i * D_ + j] = u.a * x + u.b * y; r_[k * D_ + j] = u.c * x + u.d * y;
            }
        }
        for (std::size_t j = 0; j < D_; ++j) if ((j & ma) && !(j & mb)) {          // columns (rho U^dagger)
            const std::size_t k = j ^ mm;
            for (std::size_t i = 0; i < D_; ++i) {
                const cplx x = r_[i * D_ + j], y = r_[i * D_ + k];
                r_[i * D_ + j] = x * std::conj(u.a) + y * std::conj(u.b); r_[i * D_ + k] = x * std::conj(u.c) + y * std::conj(u.d);
            }
        }
    }
    void cphase(int p, int q, double phi) { diag2(p, q, std::exp(cplx(0, phi))); }

    template <class RNG> int measure(int q, RNG&) {   // non-selective (no record possible)
        const std::size_t m = std::size_t(1) << q;
        for (std::size_t i = 0; i < D_; ++i) for (std::size_t j = 0; j < D_; ++j) if ((i ^ j) & m) r_[i * D_ + j] = 0;
        return -1;
    }
    template <class RNG> void reset(int q, RNG&) { kraus1(q, {M2{1, 0, 0, 0}, M2{0, 1, 0, 0}}); }
    template <class RNG> void bitflip(int q, double p, RNG&) {
        if (p > 0) kraus1(q, {scaled(gates::I(), std::sqrt(1 - p)), scaled(gates::X(), std::sqrt(p))});
    }
    void multi_cphase(const std::vector<std::pair<std::size_t, cplx>>& terms) {
        for (const auto& [mask, ph] : terms) {
            for (std::size_t i = 0; i < D_; ++i) for (std::size_t j = 0; j < D_; ++j) {
                cplx f = 1;
                if ((i & mask) == mask) f *= ph;
                if ((j & mask) == mask) f *= std::conj(ph);
                r_[i * D_ + j] *= f;
            }
        }
    }
    template <class RNG> void gad(int q, double gamma, double p_down, RNG&) {
        if (gamma <= 0) return;
        const double s = std::sqrt(1 - gamma), g = std::sqrt(gamma);
        const double a = std::sqrt(p_down), b = std::sqrt(1 - p_down);
        kraus1(q, {M2{a, 0, 0, a * s}, M2{0, a * g, 0, 0}, M2{b * s, 0, 0, b}, M2{0, 0, b * g, 0}});
    }
    template <class RNG> void dephase(int q, double pz, RNG&) {
        const std::size_t m = std::size_t(1) << q; const double f = 1 - 2 * pz;
        for (std::size_t i = 0; i < D_; ++i) for (std::size_t j = 0; j < D_; ++j) if (((i ^ j) & m)) r_[i * D_ + j] *= f;
    }
    template <class RNG> void pauli1(int q, double p, RNG&) {
        if (p <= 0) return;
        const double a = std::sqrt(1 - p), b = std::sqrt(p / 3);
        kraus1(q, {scaled(gates::I(), a), scaled(gates::X(), b), scaled(gates::Y(), b), scaled(gates::Z(), b)});
    }
    template <class RNG> void pauli2(int a, int b, double p, RNG&) {
        if (p <= 0) return;
        const M2 P[4] = {gates::I(), gates::X(), gates::Y(), gates::Z()};
        std::vector<cplx> acc(r_.size(), 0.0);
        const auto orig = r_;
        for (int k = 0; k < 16; ++k) {
            r_ = orig;
            if (k % 4) gate1(a, P[k % 4]);
            if (k / 4) gate1(b, P[k / 4]);
            const double w = (k == 0) ? 1 - p : p / 15;
            for (std::size_t i = 0; i < r_.size(); ++i) acc[i] += w * r_[i];
        }
        r_.swap(acc);
    }
    double fidelity(const StateVector& psi) const {
        cplx s = 0;
        const auto& v = psi.amps();
        for (std::size_t i = 0; i < D_; ++i) for (std::size_t j = 0; j < D_; ++j) s += std::conj(v[i]) * r_[i * D_ + j] * v[j];
        return s.real();
    }
    double prob1(int q) const { double s = 0; for (std::size_t i = 0; i < D_; ++i) if ((i >> q) & 1) s += r_[i * D_ + i].real(); return s; }
    double trace() const { double s = 0; for (std::size_t i = 0; i < D_; ++i) s += r_[i * D_ + i].real(); return s; }
    cplx coherence(int q) const {   // <0|rho_q|1> of the reduced single-qubit state
        const std::size_t m = std::size_t(1) << q; cplx s = 0;
        for (std::size_t i = 0; i < D_; ++i) if (!(i & m)) s += r_[i * D_ + (i | m)];
        return s;
    }

private:
    int n_;
    std::size_t D_;
    std::vector<cplx> r_;

    void left(int q, const M2& u) {      // rho <- U rho
        const std::size_t m = std::size_t(1) << q;
        for (std::size_t i = 0; i < D_; ++i) if (!(i & m))
            for (std::size_t j = 0; j < D_; ++j) {
                const cplx x = r_[i * D_ + j], y = r_[(i | m) * D_ + j];
                r_[i * D_ + j] = u.a * x + u.b * y;
                r_[(i | m) * D_ + j] = u.c * x + u.d * y;
            }
    }
    void right_dagger(int q, const M2& u) {   // rho <- rho U^dagger
        const std::size_t m = std::size_t(1) << q;
        for (std::size_t i = 0; i < D_; ++i)
            for (std::size_t j = 0; j < D_; ++j) if (!(j & m)) {
                const cplx x = r_[i * D_ + j], y = r_[i * D_ + (j | m)];
                r_[i * D_ + j] = x * std::conj(u.a) + y * std::conj(u.b);
                r_[i * D_ + (j | m)] = x * std::conj(u.c) + y * std::conj(u.d);
            }
    }
    void diag2(int p, int q, cplx ph) {
        const std::size_t mask = (std::size_t(1) << p) | (std::size_t(1) << q);
        for (std::size_t i = 0; i < D_; ++i) for (std::size_t j = 0; j < D_; ++j) {
            cplx f = 1;
            if ((i & mask) == mask) f *= ph;
            if ((j & mask) == mask) f *= std::conj(ph);
            r_[i * D_ + j] *= f;
        }
    }
    void kraus1(int q, const std::vector<M2>& K) {
        const auto orig = r_;
        std::vector<cplx> acc(r_.size(), 0.0);
        for (const auto& k : K) {
            r_ = orig; left(q, k); right_dagger(q, k);
            for (std::size_t i = 0; i < r_.size(); ++i) acc[i] += r_[i];
        }
        r_.swap(acc);
    }
};

} // namespace qpuforge
