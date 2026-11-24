// qpuforge — numerical core: constants, 2x2 gates, eigen-solvers, dense inverse (header-only, C++17)
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace qpuforge {

using cplx = std::complex<double>;
inline constexpr double PI = 3.14159265358979323846;

namespace phys {
inline constexpr double h = 6.62607015e-34;          // J s
inline constexpr double hbar = h / (2 * PI);
inline constexpr double e = 1.602176634e-19;         // C
inline constexpr double kB = 1.380649e-23;           // J/K
inline constexpr double Phi0 = h / (2 * e);          // Wb
inline constexpr double ueV = 1e-6 * e;              // J
} // namespace phys

// ------------------------------------------------------------------ 2x2 complex matrices
struct M2 {
    cplx a{1}, b{0}, c{0}, d{1};   // [[a, b], [c, d]]
};
inline M2 operator*(const M2& x, const M2& y) {
    return {x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d};
}
inline M2 dagger(const M2& x) { return {std::conj(x.a), std::conj(x.c), std::conj(x.b), std::conj(x.d)}; }
inline M2 scaled(const M2& x, cplx s) { return {x.a * s, x.b * s, x.c * s, x.d * s}; }

// Distance between unitaries up to a global phase: min_phi ||A - e^{i phi} B||_max
inline double phase_distance(const M2& A, const M2& B) {
    const cplx tr = std::conj(B.a) * A.a + std::conj(B.b) * A.b + std::conj(B.c) * A.c + std::conj(B.d) * A.d;
    const cplx ph = std::abs(tr) > 1e-300 ? tr / std::abs(tr) : cplx(1);
    const M2 Bp = scaled(B, ph);
    return std::max({std::abs(A.a - Bp.a), std::abs(A.b - Bp.b), std::abs(A.c - Bp.c), std::abs(A.d - Bp.d)});
}

namespace gates {
inline const cplx I1{0, 1};
inline M2 I() { return {1, 0, 0, 1}; }
inline M2 X() { return {0, 1, 1, 0}; }
inline M2 Y() { return {0, -I1, I1, 0}; }
inline M2 Z() { return {1, 0, 0, -1}; }
inline M2 H() { const double s = 1 / std::sqrt(2.0); return {s, s, s, -s}; }
inline M2 S() { return {1, 0, 0, I1}; }
inline M2 Sdg() { return {1, 0, 0, -I1}; }
inline M2 T() { return {1, 0, 0, std::exp(I1 * (PI / 4))}; }
inline M2 Tdg() { return {1, 0, 0, std::exp(-I1 * (PI / 4))}; }
inline M2 SX() { return {cplx(0.5, 0.5), cplx(0.5, -0.5), cplx(0.5, -0.5), cplx(0.5, 0.5)}; }
inline M2 SXdg() { return dagger(SX()); }
inline M2 RX(double t) { return {std::cos(t / 2), -I1 * std::sin(t / 2), -I1 * std::sin(t / 2), std::cos(t / 2)}; }
inline M2 RY(double t) { return {std::cos(t / 2), -std::sin(t / 2), std::sin(t / 2), std::cos(t / 2)}; }
inline M2 RZ(double t) { return {std::exp(-I1 * (t / 2)), 0, 0, std::exp(I1 * (t / 2))}; }
inline M2 P(double l) { return {1, 0, 0, std::exp(I1 * l)}; }
inline M2 U3(double t, double p, double l) {
    return {std::cos(t / 2), -std::exp(I1 * l) * std::sin(t / 2), std::exp(I1 * p) * std::sin(t / 2), std::exp(I1 * (p + l)) * std::cos(t / 2)};
}
// sqrt(Y) and sqrt(W), W = (X+Y)/sqrt2 — the Google XEB single-qubit set
inline M2 SY() { return RY(PI / 2); }
inline M2 SW() {
    const double s = 1 / std::sqrt(2.0);
    const M2 W{0, cplx(s, -s), cplx(s, s), 0};
    // sqrt(W) = (I - i W)/sqrt2 up to phase (W^2 = I)
    return {cplx(s, 0), -I1 * W.b * s, -I1 * W.c * s, cplx(s, 0)};
}
} // namespace gates

// ZYZ Euler angles: U = e^{i g} RZ(phi) RY(theta) RZ(lam)
struct Euler { double theta, phi, lam, gphase; };
inline Euler zyz(const M2& U) {
    const cplx det = U.a * U.d - U.b * U.c;
    const cplx s = std::sqrt(det);                      // U / s is in SU(2)
    const cplx a = U.a / s, c = U.c / s;                // SU(2): [[a, -conj(c)], [c, conj(a)]]
    const double theta = 2 * std::atan2(std::abs(c), std::abs(a));
    const double sum = 2 * std::arg(std::conj(a));       // phi + lam  (a = e^{-i(phi+lam)/2} cos)
    const double diff = 2 * std::arg(c);                 // phi - lam  (c = e^{i(phi-lam)/2} sin)
    Euler E;
    E.theta = theta;
    if (std::abs(c) < 1e-12) { E.phi = sum; E.lam = 0; }
    else if (std::abs(a) < 1e-12) { E.phi = diff; E.lam = 0; }
    else { E.phi = 0.5 * (sum + diff); E.lam = 0.5 * (sum - diff); }
    E.gphase = std::arg(s);
    return E;
}

// ------------------------------------------------------------------ real symmetric eigen-solvers
// Implicit QL for a symmetric tridiagonal matrix (diag d, off-diag e[i] between i and i+1).
// On return d holds eigenvalues (ascending), z (n x n, row-major) holds eigenvectors as columns.
inline void tridiag_eigen(std::vector<double>& d, std::vector<double> e_off, std::vector<double>& z) {
    const int n = int(d.size());
    std::vector<double> e(n, 0.0);
    for (int i = 0; i + 1 < n; ++i) e[i] = e_off[i];
    z.assign(std::size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) z[std::size_t(i) * n + i] = 1.0;
    for (int l = 0; l < n; ++l) {
        int iter = 0, m;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) <= 1e-15 * dd) break;
            }
            if (m != l) {
                if (iter++ == 100) throw std::runtime_error("tridiag_eigen: no convergence");
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? std::abs(r) : -std::abs(r)));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; --i) {
                    double f = s * e[i];
                    const double b = c * e[i];
                    e[i + 1] = (r = std::hypot(f, g));
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    for (int k = 0; k < n; ++k) {
                        f = z[std::size_t(k) * n + i + 1];
                        z[std::size_t(k) * n + i + 1] = s * z[std::size_t(k) * n + i] + c * f;
                        z[std::size_t(k) * n + i] = c * z[std::size_t(k) * n + i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }
    // sort ascending
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int x, int y) { return d[x] < d[y]; });
    std::vector<double> d2(n), z2(z.size());
    for (int j = 0; j < n; ++j) {
        d2[j] = d[idx[j]];
        for (int k = 0; k < n; ++k) z2[std::size_t(k) * n + j] = z[std::size_t(k) * n + idx[j]];
    }
    d.swap(d2); z.swap(z2);
}

// Cyclic Jacobi for a dense real symmetric matrix A (n x n row-major). Returns ascending eigenvalues; V columns = eigenvectors.
inline std::vector<double> sym_eigen(std::vector<double> A, int n, std::vector<double>& V) {
    V.assign(std::size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) V[std::size_t(i) * n + i] = 1.0;
    auto at = [&](int i, int j) -> double& { return A[std::size_t(i) * n + j]; };
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0, diag = 0;
        for (int i = 0; i < n; ++i) { diag += at(i, i) * at(i, i); for (int j = i + 1; j < n; ++j) off += at(i, j) * at(i, j); }
        if (off <= 1e-30 * (diag + 1e-300)) break;
        for (int p = 0; p < n; ++p) for (int q = p + 1; q < n; ++q) {
            const double apq = at(p, q);
            if (std::abs(apq) < 1e-300) continue;
            const double theta = (at(q, q) - at(p, p)) / (2 * apq);
            const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
            const double c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < n; ++k) {
                const double akp = at(k, p), akq = at(k, q);
                at(k, p) = c * akp - s * akq; at(k, q) = s * akp + c * akq;
            }
            for (int k = 0; k < n; ++k) {
                const double apk = at(p, k), aqk = at(q, k);
                at(p, k) = c * apk - s * aqk; at(q, k) = s * apk + c * aqk;
            }
            for (int k = 0; k < n; ++k) {
                const double vkp = V[std::size_t(k) * n + p], vkq = V[std::size_t(k) * n + q];
                V[std::size_t(k) * n + p] = c * vkp - s * vkq; V[std::size_t(k) * n + q] = s * vkp + c * vkq;
            }
        }
    }
    std::vector<double> ev(n);
    for (int i = 0; i < n; ++i) ev[i] = at(i, i);
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int x, int y) { return ev[x] < ev[y]; });
    std::vector<double> ev2(n), V2(V.size());
    for (int j = 0; j < n; ++j) {
        ev2[j] = ev[idx[j]];
        for (int k = 0; k < n; ++k) V2[std::size_t(k) * n + j] = V[std::size_t(k) * n + idx[j]];
    }
    V.swap(V2);
    return ev2;
}

// Dense inverse by Gauss-Jordan with partial pivoting.
inline std::vector<double> inverse(std::vector<double> A, int n) {
    std::vector<double> I(std::size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) I[std::size_t(i) * n + i] = 1;
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::abs(A[std::size_t(r) * n + c]) > std::abs(A[std::size_t(piv) * n + c])) piv = r;
        if (std::abs(A[std::size_t(piv) * n + c]) < 1e-300) throw std::runtime_error("inverse: singular");
        if (piv != c) for (int k = 0; k < n; ++k) { std::swap(A[std::size_t(c) * n + k], A[std::size_t(piv) * n + k]); std::swap(I[std::size_t(c) * n + k], I[std::size_t(piv) * n + k]); }
        const double inv = 1.0 / A[std::size_t(c) * n + c];
        for (int k = 0; k < n; ++k) { A[std::size_t(c) * n + k] *= inv; I[std::size_t(c) * n + k] *= inv; }
        for (int r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = A[std::size_t(r) * n + c];
            if (f == 0) continue;
            for (int k = 0; k < n; ++k) { A[std::size_t(r) * n + k] -= f * A[std::size_t(c) * n + k]; I[std::size_t(r) * n + k] -= f * I[std::size_t(c) * n + k]; }
        }
    }
    return I;
}

} // namespace qpuforge
