// qpuforge — physical layout layer
//
//  * quasi-static capacitance extraction by the method of moments: zero-thickness conductors on the surface of a thick
//    substrate (permittivity eps_r below, vacuum above). For charges on the interface the Green's function is exactly
//    1 / (4 pi eps0 eps_eff r) with eps_eff = (eps_r + 1) / 2, so the problem is a free-space one with eps_eff.
//    Conductors are unions of axis-aligned rectangles (the ground plane = an outline minus cut-outs); they are meshed by
//    an adaptive quadtree graded towards edges (where the surface charge diverges as d^-1/2) and clipped exactly to
//    the geometry. Uniform charge per panel, collocation at panel centres with the analytic potential of a uniformly
//    charged rectangle, dense LU.
//  * parametric cells (Xmon-style qubit, tunable coupler with claw ends, readout claw) and synthesis of their
//    dimensions against circuit targets (qubit E_C, qubit-coupler g, direct g12, coupler E_C);
//  * GDSII stream writer for the assembled chip.
#pragma once
#include "device.hpp"
#include "gds.hpp"
#include "pulse.hpp"
#include <chrono>
#include <fstream>

namespace qpuforge {
namespace layout {

struct Rect {
    double x0{0}, y0{0}, x1{0}, y1{0};
    double w() const { return x1 - x0; }
    double h() const { return y1 - y0; }
    double area() const { return w() * h(); }
    double cx() const { return 0.5 * (x0 + x1); }
    double cy() const { return 0.5 * (y0 + y1); }
    bool intersects(const Rect& o) const { return x0 < o.x1 && o.x0 < x1 && y0 < o.y1 && o.y0 < y1; }
    bool covers(const Rect& o) const { return x0 <= o.x0 && x1 >= o.x1 && y0 <= o.y0 && y1 >= o.y1; }
    bool contains(double x, double y) const { return x > x0 && x < x1 && y > y0 && y < y1; }
    Rect grow(double g) const { return {x0 - g, y0 - g, x1 + g, y1 + g}; }
    Rect shift(double dx, double dy) const { return {x0 + dx, y0 + dy, x1 + dx, y1 + dy}; }
};
inline Rect rect_c(double cx, double cy, double w, double h) { return {cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2}; }
inline Rect bbox_of(const std::vector<Rect>& rs) {
    Rect b{1e300, 1e300, -1e300, -1e300};
    for (const auto& r : rs) { b.x0 = std::min(b.x0, r.x0); b.y0 = std::min(b.y0, r.y0); b.x1 = std::max(b.x1, r.x1); b.y1 = std::max(b.y1, r.y1); }
    return b;
}

// Metal conductor: union of rects. Ground: outline minus the union of cut-outs.
struct Conductor {
    std::string name;
    std::vector<Rect> rects;
    bool ground{false};
    Rect outline;
    std::vector<Rect> holes;
    bool inside(double x, double y) const {
        if (!ground) { for (const auto& r : rects) if (r.contains(x, y)) return true; return false; }
        if (!outline.contains(x, y)) return false;
        for (const auto& h : holes) if (h.contains(x, y)) return false;
        return true;
    }
    Rect extent() const { return ground ? outline : bbox_of(rects); }
};

struct MeshOptions {
    double h_min{1.0};       // um, panel size at an exposed edge (edge resolution)
    double h_max{40.0};      // um, largest panel
    double growth{1.6};      // geometric growth of panel size away from exposed edges
};

struct Panel { Rect r; int cond; };

// Disjoint rectangular pieces of a conductor, each with flags telling which sides are exposed edges
// (bottom, right, top, left) — the only places the surface charge needs fine resolution.
struct Piece { Rect r; std::array<bool, 4> exposed{{true, true, true, true}}; };

inline std::vector<Piece> decompose(const Conductor& c) {
    std::vector<double> xs, ys;
    auto addr = [&](const Rect& r) { xs.push_back(r.x0); xs.push_back(r.x1); ys.push_back(r.y0); ys.push_back(r.y1); };
    if (c.ground) { addr(c.outline); for (const auto& h : c.holes) addr(h); }
    else for (const auto& r : c.rects) addr(r);
    std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    const int nx = int(xs.size()) - 1, ny = int(ys.size()) - 1;
    std::vector<std::uint8_t> in(std::size_t(std::max(nx, 0)) * std::size_t(std::max(ny, 0)), 0), used(in.size(), 0);
    auto id = [&](int i, int j) { return std::size_t(j) * std::size_t(nx) + std::size_t(i); };
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
        in[id(i, j)] = c.inside(0.5 * (xs[std::size_t(i)] + xs[std::size_t(i) + 1]), 0.5 * (ys[std::size_t(j)] + ys[std::size_t(j) + 1]));
    auto isin = [&](int i, int j) { return i >= 0 && j >= 0 && i < nx && j < ny && in[id(i, j)]; };
    std::vector<Piece> out;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        if (!in[id(i, j)] || used[id(i, j)]) continue;
        int i1 = i;
        while (i1 + 1 < nx && in[id(i1 + 1, j)] && !used[id(i1 + 1, j)]) ++i1;
        int j1 = j;
        for (;;) {
            if (j1 + 1 >= ny) break;
            bool ok = true;
            for (int k = i; k <= i1; ++k) if (!in[id(k, j1 + 1)] || used[id(k, j1 + 1)]) { ok = false; break; }
            if (!ok) break;
            ++j1;
        }
        for (int jj = j; jj <= j1; ++jj) for (int k = i; k <= i1; ++k) used[id(k, jj)] = 1;
        Piece p; p.r = {xs[std::size_t(i)], ys[std::size_t(j)], xs[std::size_t(i1) + 1], ys[std::size_t(j1) + 1]};
        bool b = false, r = false, tp = false, l = false;
        for (int k = i; k <= i1; ++k) { b |= !isin(k, j - 1); tp |= !isin(k, j1 + 1); }
        for (int jj = j; jj <= j1; ++jj) { l |= !isin(i - 1, jj); r |= !isin(i1 + 1, jj); }
        p.exposed = {{b, r, tp, l}};
        out.push_back(p);
    }
    return out;
}

// 1-D breakpoints on [a, b], graded geometrically from exposed ends
inline std::vector<double> graded_points(double a, double b, bool fine_a, bool fine_b, const MeshOptions& o) {
    const double L = b - a;
    std::vector<double> left{0.0}, right{0.0};
    auto ramp = [&](std::vector<double>& v) { double h = o.h_min, x = 0; while (x + h < 0.5 * L && h < o.h_max) { x += h; v.push_back(x); h *= o.growth; } };
    if (fine_a) ramp(left);
    if (fine_b) ramp(right);
    const double la = left.back(), lb = right.back();
    const double mid = L - la - lb;
    std::vector<double> pts;
    for (double x : left) pts.push_back(a + x);
    const int n = std::max(1, int(std::ceil(mid / o.h_max - 1e-9)));
    for (int k = 1; k < n; ++k) pts.push_back(a + la + mid * k / n);
    for (auto it = right.rbegin(); it != right.rend(); ++it) pts.push_back(b - *it);
    std::sort(pts.begin(), pts.end());
    std::vector<double> u;
    for (double x : pts) if (u.empty() || x - u.back() > 1e-9) u.push_back(x);
    if (u.back() < b - 1e-9) u.push_back(b); else u.back() = b;
    return u;
}

inline std::vector<Panel> mesh(const std::vector<Conductor>& conds, const MeshOptions& o = {}) {
    std::vector<Panel> out;
    for (int k = 0; k < int(conds.size()); ++k)
        for (const auto& p : decompose(conds[std::size_t(k)])) {
            const auto X = graded_points(p.r.x0, p.r.x1, p.exposed[3], p.exposed[1], o);
            const auto Y = graded_points(p.r.y0, p.r.y1, p.exposed[0], p.exposed[2], o);
            for (std::size_t i = 0; i + 1 < X.size(); ++i)
                for (std::size_t j = 0; j + 1 < Y.size(); ++j) out.push_back({{X[i], Y[j], X[i + 1], Y[j + 1]}, k});
        }
    return out;
}

// potential integral of a uniformly charged rectangle (unit density) at an in-plane point: int int dA / r
inline double rect_potential(const Rect& r, double x, double y) {
    auto F = [](double u, double v) {
        double s = 0;
        if (u != 0) s += u * std::asinh(v / std::abs(u));
        if (v != 0) s += v * std::asinh(u / std::abs(v));
        return s;
    };
    const double u1 = r.x0 - x, u2 = r.x1 - x, v1 = r.y0 - y, v2 = r.y1 - y;
    return F(u2, v2) - F(u2, v1) - F(u1, v2) + F(u1, v1);
}

struct CapacitanceMatrix {
    std::vector<std::string> names;
    int n{0};
    std::vector<double> C_fF;      // Maxwell matrix (row-major), symmetrised
    double asymmetry{0};           // max |C - C^T| / max |C| before symmetrisation (discretisation indicator)
    int panels{0};
    double seconds{0};
    double at(int i, int j) const { return C_fF[std::size_t(i) * n + j]; }
    int index(const std::string& s) const { for (int i = 0; i < n; ++i) if (names[std::size_t(i)] == s) return i; throw std::runtime_error("no conductor " + s); }
    double mutual(const std::string& a, const std::string& b) const { return -at(index(a), index(b)); }
    double self(const std::string& a) const { return at(index(a), index(a)); }
};

// Geometry sanity: distinct metal conductors must not overlap, and metal must sit inside the ground cut-outs.
inline std::string check_geometry(const std::vector<Conductor>& conds) {
    for (std::size_t i = 0; i < conds.size(); ++i) {
        if (conds[i].ground) continue;
        for (std::size_t j = i + 1; j < conds.size(); ++j) {
            if (conds[j].ground) continue;
            for (const auto& a : conds[i].rects) for (const auto& b : conds[j].rects)
                if (a.intersects(b)) return conds[i].name + " overlaps " + conds[j].name;
        }
        for (const auto& g : conds) if (g.ground)
            for (const auto& a : conds[i].rects)
                for (double fx : {0.03, 0.47, 0.97}) for (double fy : {0.03, 0.47, 0.97})
                    if (g.inside(a.x0 + fx * a.w(), a.y0 + fy * a.h())) return conds[i].name + " touches the ground plane";
    }
    return "";
}

inline CapacitanceMatrix extract_capacitance(const std::vector<Conductor>& conds, double eps_r = 11.45, const MeshOptions& mo = {}) {
    const auto t0 = std::chrono::steady_clock::now();
    if (const auto err = check_geometry(conds); !err.empty()) throw std::runtime_error("extract_capacitance: " + err);
    const auto P = mesh(conds, mo);
    const int N = int(P.size()), K = int(conds.size());
    std::vector<double> A(std::size_t(N) * N);
#pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < N; ++i) {
        const double x = P[std::size_t(i)].r.cx(), y = P[std::size_t(i)].r.cy();
        for (int j = 0; j < N; ++j) A[std::size_t(i) * N + j] = rect_potential(P[std::size_t(j)].r, x, y) / P[std::size_t(j)].r.area();
    }
    // LU with partial pivoting (in place)
    std::vector<int> piv(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) {
        int p = k; double best = std::abs(A[std::size_t(k) * N + k]);
        for (int i = k + 1; i < N; ++i) { const double v = std::abs(A[std::size_t(i) * N + k]); if (v > best) { best = v; p = i; } }
        piv[std::size_t(k)] = p;
        if (p != k) for (int j = 0; j < N; ++j) std::swap(A[std::size_t(k) * N + j], A[std::size_t(p) * N + j]);
        const double d = A[std::size_t(k) * N + k];
        if (d == 0) throw std::runtime_error("extract_capacitance: singular system");
        double* rowk = &A[std::size_t(k) * N];
#pragma omp parallel for schedule(static) if (N - k > 256)
        for (int i = k + 1; i < N; ++i) {
            double* rowi = &A[std::size_t(i) * N];
            const double f = rowi[k] / d;
            rowi[k] = f;
            if (f != 0) for (int j = k + 1; j < N; ++j) rowi[j] -= f * rowk[j];
        }
    }
    const double eps_eff = 0.5 * (eps_r + 1);
    const double Kc = 1e6 / (4 * PI * phys::eps0 * eps_eff);     // phi = Kc * sum_j (I/A)_um * q_j
    CapacitanceMatrix C;
    C.n = K; C.panels = N;
    for (const auto& c : conds) C.names.push_back(c.name);
    std::vector<double> Cr(std::size_t(K) * K, 0.0);
    std::vector<double> b(static_cast<std::size_t>(N));
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < N; ++i) b[std::size_t(i)] = P[std::size_t(i)].cond == k ? 1.0 : 0.0;
        for (int i = 0; i < N; ++i) std::swap(b[std::size_t(i)], b[std::size_t(piv[std::size_t(i)])]);
        for (int i = 0; i < N; ++i) { double s = b[std::size_t(i)]; for (int j = 0; j < i; ++j) s -= A[std::size_t(i) * N + j] * b[std::size_t(j)]; b[std::size_t(i)] = s; }
        for (int i = N - 1; i >= 0; --i) { double s = b[std::size_t(i)]; for (int j = i + 1; j < N; ++j) s -= A[std::size_t(i) * N + j] * b[std::size_t(j)]; b[std::size_t(i)] = s / A[std::size_t(i) * N + i]; }
        for (int i = 0; i < N; ++i) Cr[std::size_t(k) * K + std::size_t(P[std::size_t(i)].cond)] += b[std::size_t(i)] / Kc * 1e15;
    }
    double mx = 0, asym = 0;
    for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) { mx = std::max(mx, std::abs(Cr[std::size_t(i) * K + j])); asym = std::max(asym, std::abs(Cr[std::size_t(i) * K + j] - Cr[std::size_t(j) * K + i])); }
    C.asymmetry = mx > 0 ? asym / mx : 0;
    C.C_fF.assign(std::size_t(K) * K, 0.0);
    for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) C.C_fF[std::size_t(i) * K + j] = 0.5 * (Cr[std::size_t(i) * K + j] + Cr[std::size_t(j) * K + i]);
    C.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return C;
}

// ---------------------------------------------------------------- analytic references (validation)
inline double ellipK(double k) {           // complete elliptic integral of the first kind K(k), via AGM
    double a = 1, b = std::sqrt(1 - k * k);
    for (int i = 0; i < 40; ++i) { const double an = 0.5 * (a + b); b = std::sqrt(a * b); a = an; }
    return PI / (2 * a);
}
// CPW on a thick substrate with finite-width grounds (conformal mapping): capacitance per unit length, F/m
inline double cpw_capacitance_per_m(double w, double s, double g, double eps_r) {
    const double a = w / 2, b = w / 2 + s, c = w / 2 + s + g;
    const double k = (a / b) * std::sqrt((1 - b * b / (c * c)) / (1 - a * a / (c * c)));
    const double kp = std::sqrt(1 - k * k);
    return 2 * phys::eps0 * (eps_r + 1) * ellipK(k) / ellipK(kp);
}

// ---------------------------------------------------------------- parametric cells
// Local frames: an element attached to a qubit's +x arm is drawn for the +x direction and rotated by k * 90 deg.
inline Rect rot90(const Rect& r, int k, double cx = 0, double cy = 0) {
    Rect q = r.shift(-cx, -cy);
    for (int i = 0; i < ((k % 4) + 4) % 4; ++i) q = {-q.y1, q.x0, -q.y0, q.x1};   // (x, y) -> (-y, x)
    return q.shift(cx, cy);
}

struct XmonParams {
    double L{130};          // arm length, centre to the outer edge of the T-bar (um)
    double w_arm{16}, gap{30};
    double t_w{12}, T{120}; // arm-end T-bar: thickness along the arm, length across it
    double r_pad{40};       // readout coupling pad in the (+x,+y) inner corner
    double claw_gap{5}, claw_w{10};
    double squid_stub{16};
};
struct CouplerGeom {
    double g_c{5}, t_w{12};
    double y1{32.5};        // coupler at -y1, direct-coupling bus at +y1 (local +x frame)
    double l_c{52};         // coupler T-bar length (sets C_qc)
    double w_body{10}, gap{20};
    double pad{150};        // coupler mid pad (sets the coupler E_C)
    double l_b{30};         // bus T-bar length (sets the direct coupling)
    double w_bus{6}, gap_bus{15};
};

struct Drawn {
    std::vector<Rect> metal;     // islands inside cut-outs (qubit, coupler, bus, claws, centre conductors)
    std::vector<Rect> holes;     // ground cut-outs
    std::vector<Rect> jj;        // Josephson-junction / SQUID regions (GDS layer 2)
    std::vector<Rect> bumps;     // flip-chip bump pads (GDS layer 3)
    std::vector<std::pair<std::string, std::pair<double, double>>> labels;
};

// Qubit metal (cross + T-bars + readout pad + SQUID stub) and its cut-out, centred at (cx, cy).
inline void xmon(const XmonParams& q, double cx, double cy, std::vector<Rect>& metal, std::vector<Rect>& holes, std::vector<Rect>* jj = nullptr) {
    const double hw = q.w_arm / 2;
    for (int k = 0; k < 4; ++k) {
        metal.push_back(rot90({0, -hw, q.L - q.t_w, hw}, k).shift(cx, cy));
        metal.push_back(rot90({q.L - q.t_w, -q.T / 2, q.L, q.T / 2}, k).shift(cx, cy));
        holes.push_back(rot90({-hw - q.gap, -hw - q.gap, q.L - q.t_w + q.gap, hw + q.gap}, k).shift(cx, cy));
        holes.push_back(rot90({q.L - q.t_w - q.gap, -q.T / 2 - q.gap, q.L + q.gap, q.T / 2 + q.gap}, k).shift(cx, cy));
    }
    metal.push_back(rect_c(cx, cy, q.w_arm, q.w_arm));
    metal.push_back(Rect{hw, hw, hw + q.r_pad, hw + q.r_pad}.shift(cx, cy));
    const double cl = hw + q.r_pad + q.claw_gap + q.claw_w + q.gap;
    holes.push_back(Rect{hw, hw, cl, cl}.shift(cx, cy));
    metal.push_back(Rect{-hw - q.squid_stub, -hw - q.squid_stub, -hw, -hw}.shift(cx, cy));
    holes.push_back(Rect{-hw - q.squid_stub - q.gap, -hw - q.squid_stub - q.gap, -hw, -hw}.shift(cx, cy));
    if (jj) jj->push_back(Rect{-hw - q.squid_stub - q.gap, -hw - q.squid_stub + 2, -hw - q.squid_stub, -hw - 2}.shift(cx, cy));
}
// Readout claw (L-shape around the corner pad) + the start of the resonator centre conductor, leaving along the
// diagonal between the two T-bars as a staircase of 10 um squares (CPW: 10 um centre, 6 um gaps).
inline std::vector<Rect> readout_claw(const XmonParams& q, double cx, double cy, double reach = 170, std::vector<Rect>* holes = nullptr) {
    const double hw = q.w_arm / 2, a = hw + q.r_pad + q.claw_gap, b = a + q.claw_w, lo = hw + q.claw_gap + 10;
    std::vector<Rect> m{Rect{lo, a, b, b}.shift(cx, cy), Rect{a, lo, b, a}.shift(cx, cy)};
    for (double s = b - 10; s + 10 <= reach; s += 5) {
        m.push_back(Rect{s, s, s + 10, s + 10}.shift(cx, cy));
        if (holes) holes->push_back(Rect{s - 6, s - 6, s + 16, s + 16}.shift(cx, cy));
    }
    return m;
}
// Coupler (grounded tunable transmon: T-bars at both ends, body, mid pad) and bus (floating direct-coupling strip)
// between a qubit at (cx, cy) and its neighbour at distance p in direction k.
struct PairElements { std::vector<Rect> coupler, bus, holes, jj; };
inline PairElements coupler_and_bus(const XmonParams& q, const CouplerGeom& c, double p, int k, double cx, double cy, bool far_end = true,
                                    double stub = 200, double L_far = -1) {
    PairElements E;
    const double xa = q.L + c.g_c, xb = p - (L_far > 0 ? L_far : q.L) - c.g_c;     // faces of the T-bars opposite each qubit
    const double xend = far_end ? xb : xa + c.t_w + stub;    // truncated near end (extraction of the qubit's loading)
    auto add = [&](std::vector<Rect>& v, Rect r) { v.push_back(rot90(r, k).shift(cx, cy)); };
    // coupler
    add(E.coupler, {xa, -c.y1 - c.l_c / 2, xa + c.t_w, -c.y1 + c.l_c / 2});
    add(E.coupler, {xa + c.t_w, -c.y1 - c.w_body / 2, std::min(xend, xb - c.t_w), -c.y1 + c.w_body / 2});
    if (far_end) add(E.coupler, {xb - c.t_w, -c.y1 - c.l_c / 2, xb, -c.y1 + c.l_c / 2});
    if (far_end) add(E.coupler, {p / 2 - c.pad / 2, -c.y1 - c.w_body / 2 - c.pad, p / 2 + c.pad / 2, -c.y1 - c.w_body / 2});
    add(E.holes, {xa - c.gap, -c.y1 - c.l_c / 2 - c.gap, xa + c.t_w + c.gap, -c.y1 + c.l_c / 2 + c.gap});
    add(E.holes, {xa, -c.y1 - c.w_body / 2 - c.gap, std::min(xend, xb - c.t_w) + (far_end ? 0 : c.gap), -c.y1 + c.w_body / 2 + c.gap});
    if (far_end) add(E.holes, {xb - c.t_w - c.gap, -c.y1 - c.l_c / 2 - c.gap, xb + c.gap, -c.y1 + c.l_c / 2 + c.gap});
    if (far_end) {
        add(E.holes, {p / 2 - c.pad / 2 - c.gap, -c.y1 - c.w_body / 2 - c.pad - c.gap, p / 2 + c.pad / 2 + c.gap, -c.y1 - c.w_body / 2});
        add(E.jj, {p / 2 - 4, -c.y1 - c.w_body / 2 - c.pad - c.gap, p / 2 + 4, -c.y1 - c.w_body / 2 - c.pad});
    }
    // bus
    if (c.l_b > 0) {
        add(E.bus, {xa, c.y1 - c.l_b / 2, xa + c.t_w, c.y1 + c.l_b / 2});
        add(E.bus, {xa + c.t_w, c.y1 - c.w_bus / 2, far_end ? xb - c.t_w : xend, c.y1 + c.w_bus / 2});
        if (far_end) add(E.bus, {xb - c.t_w, c.y1 - c.l_b / 2, xb, c.y1 + c.l_b / 2});
        add(E.holes, {xa - c.gap_bus, c.y1 - c.l_b / 2 - c.gap_bus, xa + c.t_w + c.gap_bus, c.y1 + c.l_b / 2 + c.gap_bus});
        add(E.holes, {xa, c.y1 - c.w_bus / 2 - c.gap_bus, (far_end ? xb - c.t_w : xend) + (far_end ? 0 : c.gap_bus), c.y1 + c.w_bus / 2 + c.gap_bus});
        if (far_end) add(E.holes, {xb - c.t_w - c.gap_bus, c.y1 - c.l_b / 2 - c.gap_bus, xb + c.gap_bus, c.y1 + c.l_b / 2 + c.gap_bus});
    }
    return E;
}

// ---------------------------------------------------------------- extraction cells and synthesis
struct CellSet { XmonParams meas, data; CouplerGeom cg; double pitch{900}; };

// Star cell: one qubit with its 4 couplers + buses (truncated beyond the mid pad) and the readout claw.
inline CapacitanceMatrix star_cell(const XmonParams& q, const CouplerGeom& c, double p, const MeshOptions& mo) {
    std::vector<Conductor> cs;
    Conductor Q; Q.name = "qubit";
    std::vector<Rect> holes;
    xmon(q, 0, 0, Q.rects, holes);
    cs.push_back(Q);
    for (int k = 0; k < 4; ++k) {
        auto E = coupler_and_bus(q, c, p, k, 0, 0, false);
        Conductor C; C.name = "coupler" + std::to_string(k); C.rects = E.coupler; cs.push_back(C);
        if (!E.bus.empty()) { Conductor B; B.name = "bus" + std::to_string(k); B.rects = E.bus; cs.push_back(B); }
        holes.insert(holes.end(), E.holes.begin(), E.holes.end());
    }
    Conductor R; R.name = "readout"; R.rects = readout_claw(q, 0, 0, 170, &holes);
    cs.push_back(R);
    Conductor G; G.name = "ground"; G.ground = true;
    const double ext = q.L + c.g_c + c.t_w + 200 + 100;
    G.outline = {-ext, -ext, ext, ext}; G.holes = holes;
    cs.push_back(G);
    return extract_capacitance(cs, 11.45, mo);
}
// Pair cell: measure qubit at 0, data qubit at +p, coupler + bus between (other arms bare).
inline CapacitanceMatrix pair_cell(const XmonParams& qa, const XmonParams& qb, const CouplerGeom& c, double p, const MeshOptions& mo) {
    std::vector<Conductor> cs;
    std::vector<Rect> holes;
    Conductor A; A.name = "qa"; xmon(qa, 0, 0, A.rects, holes);
    Conductor B; B.name = "qb"; xmon(qb, p, 0, B.rects, holes);
    auto E = coupler_and_bus(qa, c, p, 0, 0, 0, true, 200, qb.L);
    Conductor C; C.name = "coupler"; C.rects = E.coupler;
    Conductor U; U.name = "bus"; U.rects = E.bus;
    holes.insert(holes.end(), E.holes.begin(), E.holes.end());
    cs.push_back(A); cs.push_back(C);
    if (!U.rects.empty()) cs.push_back(U);
    cs.push_back(B);
    Conductor G; G.name = "ground"; G.ground = true;
    const double ex = std::max(qa.L, qb.L) + 150;
    G.outline = {-ex, -ex, p + ex, ex}; G.holes = holes;
    cs.push_back(G);
    return extract_capacitance(cs, 11.45, mo);
}

// Lumped circuit of a (measure, coupler, data) triple from the two cells: qubit self-capacitances from the star cells,
// everything else from the pair cell; the floating bus is eliminated (Schur complement). Returns the 3x3 Maxwell
// matrix (fF) in the order (measure, coupler, data).
inline std::array<double, 9> triple_capacitance(const CapacitanceMatrix& star_m, const CapacitanceMatrix& star_d, const CapacitanceMatrix& pc) {
    const bool has_bus = std::find(pc.names.begin(), pc.names.end(), "bus") != pc.names.end();
    const std::vector<std::string> nodes = has_bus ? std::vector<std::string>{"qa", "coupler", "qb", "bus"} : std::vector<std::string>{"qa", "coupler", "qb"};
    const int n = int(nodes.size());
    std::vector<double> M(std::size_t(n) * n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) M[std::size_t(i) * n + j] = pc.at(pc.index(nodes[std::size_t(i)]), pc.index(nodes[std::size_t(j)]));
    // replace qubit self terms by the star-cell totals (all four couplers/buses + readout attached)
    M[0] = star_m.self("qubit");
    M[std::size_t(2) * n + 2] = star_d.self("qubit");
    std::array<double, 9> out{};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = M[std::size_t(i) * n + j];
        if (has_bus) v -= M[std::size_t(i) * n + 3] * M[std::size_t(3) * n + j] / M[std::size_t(3) * n + 3];
        out[std::size_t(i) * 3 + j] = v;
    }
    return out;
}

// Circuit parameters of the triple: charging energies (GHz) and charge-charge couplings G_ij (GHz per n_i n_j),
// H = sum 4 E_C,i n_i^2 + sum_{i<j} G_ij n_i n_j with E_C,i = e^2 (C^-1)_ii / 2h, G_ij = 4 e^2 (C^-1)_ij / h.
struct TripleCircuit { double EC_m{0}, EC_c{0}, EC_d{0}, G_mc{0}, G_cd{0}, G_md{0}; };
inline TripleCircuit triple_circuit(const std::array<double, 9>& C_fF) {
    std::vector<double> C(9);
    for (int k = 0; k < 9; ++k) C[std::size_t(k)] = C_fF[std::size_t(k)] * 1e-15;
    const auto Ci = inverse(C, 3);
    const double k1 = phys::e * phys::e / (2 * phys::h) / 1e9, k2 = 4 * phys::e * phys::e / phys::h / 1e9;
    TripleCircuit t;
    t.EC_m = k1 * Ci[0]; t.EC_c = k1 * Ci[4]; t.EC_d = k1 * Ci[8];
    t.G_mc = k2 * Ci[1]; t.G_cd = k2 * Ci[5]; t.G_md = k2 * Ci[2];
    return t;
}

// Synthesis targets for one (measure, coupler, data) triple + readout.
struct LayoutTargets {
    double EC_m{0.2}, EC_c{0.126}, EC_d{0.2};          // GHz
    double G_mc{0.0575}, G_cd{0.061}, G_md{0.0056};     // GHz per n n
    double C_qr_fF{6.6};                                 // qubit-readout claw capacitance
};
struct LayoutMetrics {
    double EC_m{0}, EC_c{0}, EC_d{0}, G_mc{0}, G_cd{0}, G_md{0}, C_qr_fF{0};
    double C_sigma_m{0}, C_sigma_d{0}, C_c{0}, C_qc{0}, C_qb{0}, C_bus{0}, C_qq{0}, C_qground{0};
    int panels{0}; double seconds{0};
};
inline LayoutMetrics evaluate_cells(const CellSet& S, const MeshOptions& mo) {
    const auto sm = star_cell(S.meas, S.cg, S.pitch, mo), sd = star_cell(S.data, S.cg, S.pitch, mo);
    const auto pc = pair_cell(S.meas, S.data, S.cg, S.pitch, mo);
    const auto tc = triple_circuit(triple_capacitance(sm, sd, pc));
    LayoutMetrics m;
    m.EC_m = tc.EC_m; m.EC_c = tc.EC_c; m.EC_d = tc.EC_d; m.G_mc = tc.G_mc; m.G_cd = tc.G_cd; m.G_md = tc.G_md;
    m.C_qr_fF = 0.5 * (sm.mutual("qubit", "readout") + sd.mutual("qubit", "readout"));
    m.C_sigma_m = sm.self("qubit"); m.C_sigma_d = sd.self("qubit"); m.C_c = pc.self("coupler");
    m.C_qc = 0.5 * (pc.mutual("qa", "coupler") + pc.mutual("qb", "coupler"));
    const bool bus = std::find(pc.names.begin(), pc.names.end(), "bus") != pc.names.end();
    m.C_qb = bus ? 0.5 * (pc.mutual("qa", "bus") + pc.mutual("qb", "bus")) : 0.0;
    m.C_bus = bus ? pc.self("bus") : 0.0;
    m.C_qq = pc.mutual("qa", "qb");
    m.C_qground = sm.mutual("qubit", "ground");
    m.panels = sm.panels + sd.panels + pc.panels; m.seconds = sm.seconds + sd.seconds + pc.seconds;
    return m;
}

// Bounded Gauss-Newton (finite-difference Jacobian each step, backtracking line search) on 6 geometric knobs:
//   measure arm length, data arm length, coupler pad, coupler T-bar, bus T-bar, readout pad
// against 6 targets: E_C,m, E_C,d, E_C,c, mean(G_mc, G_cd), G_md, C_qr (relative residuals). Runs on a coarse mesh,
// then re-targets with the measured fine/coarse ratio and converges again (mesh-error correction).
struct SynthesisLog { std::vector<std::string> lines; LayoutMetrics final_fine; int evaluations{0}; double max_rel_err{0}; std::array<double, 6> final_rel{}; };
inline CellSet synthesize_cells(CellSet S, const LayoutTargets& T, SynthesisLog* log = nullptr, int newton_iters = 5,
                                MeshOptions coarse = {3.0, 50.0, 2.0}, MeshOptions fine = {1.5, 50.0, 1.6}) {
    constexpr int K = 6;
    using Vec = std::array<double, K>;
    auto get = [](const CellSet& s) { return Vec{s.meas.L, s.data.L, s.cg.pad, s.cg.l_c, s.cg.l_b, s.meas.r_pad}; };
    auto put = [](CellSet& s, const Vec& x) {
        s.meas.L = x[0]; s.data.L = x[1]; s.cg.pad = x[2]; s.cg.l_c = x[3]; s.cg.l_b = x[4]; s.meas.r_pad = s.data.r_pad = x[5];
    };
    const double Tb = S.meas.T;
    auto bounds = [&](const Vec& x, Vec& lo, Vec& hi) {
        lo = {70, 70, 20, 8, 0, 12};
        hi = {300, 300, 420, Tb - 26, std::max(0.0, Tb - 26 - x[3]), 80};
        for (int q = 0; q < 2; ++q) {          // claw clear of the arm T-bars; resonator exit clear of the T-bar corners
            const XmonParams& Q = q ? S.data : S.meas;
            const double claw_outer = Q.w_arm / 2 + x[5] + Q.claw_gap + Q.claw_w;
            lo[std::size_t(q)] = std::max({lo[std::size_t(q)], Tb / 2 + 20, claw_outer + Q.t_w + 12, 80.0});
        }
    };
    auto clampx = [&](Vec x) { Vec lo, hi; bounds(x, lo, hi); for (int i = 0; i < K; ++i) x[std::size_t(i)] = std::clamp(x[std::size_t(i)], lo[std::size_t(i)], std::max(lo[std::size_t(i)], hi[std::size_t(i)])); return x; };
    const Vec w{1, 1, 1, 1, 0.7, 0.5};    // residual weights (G_md can be trimmed by the coupler park point; C_qr by the resonator)
    auto metric = [&](const LayoutMetrics& m) { return Vec{m.EC_m, m.EC_d, m.EC_c, 0.5 * (m.G_mc + m.G_cd), m.G_md, m.C_qr_fF}; };
    Vec target{T.EC_m, T.EC_d, T.EC_c, 0.5 * (T.G_mc + T.G_cd), T.G_md, T.C_qr_fF};
    const Vec target0 = target;
    auto say = [&](const std::string& s) { if (log) log->lines.push_back(s); };
    auto line = [&](const char* tag, const Vec& x, const Vec& f) {
        char b[400];
        std::snprintf(b, sizeof b, "%-18s L_m %6.1f L_d %6.1f pad %6.1f l_c %5.1f l_b %5.1f r %5.1f | EC_m %.4f EC_d %.4f EC_c %.4f G %.4f G_md %.5f C_qr %.2f",
                      tag, x[0], x[1], x[2], x[3], x[4], x[5], f[0], f[1], f[2], f[3], f[4], f[5]);
        say(b);
    };
    auto eval = [&](const Vec& x, const MeshOptions& mo) { CellSet s = S; put(s, x); if (log) log->evaluations++; return metric(evaluate_cells(s, mo)); };
    auto res = [&](const Vec& f) { Vec r; for (int i = 0; i < K; ++i) r[std::size_t(i)] = w[std::size_t(i)] * (f[std::size_t(i)] / target[std::size_t(i)] - 1); return r; };
    auto norm = [](const Vec& r) { double s = 0; for (double v : r) s += v * v; return std::sqrt(s); };
    Vec x = clampx(get(S));
    for (int pass = 0; pass < 2; ++pass) {
        Vec f = eval(x, coarse);
        line(pass ? "pass 2 start" : "start", x, f);
        for (int it = 0; it < newton_iters; ++it) {
            Vec r = res(f);
            double worst = 0; for (double v : r) worst = std::max(worst, std::abs(v));
            if (worst < 0.004) break;
            std::array<Vec, K> J{};
            for (int k = 0; k < K; ++k) {
                Vec xp = x; const double h = std::max(1.5, 0.04 * std::abs(x[std::size_t(k)]));
                xp[std::size_t(k)] += h;
                const Vec rp = res(eval(xp, coarse));
                for (int i = 0; i < K; ++i) J[std::size_t(i)][std::size_t(k)] = (rp[std::size_t(i)] - r[std::size_t(i)]) / h;
            }
            Vec lo, hi; bounds(x, lo, hi);
            std::array<bool, K> fixed{};
            Vec dx{};
            for (int round = 0; round < K; ++round) {
                // least squares on the free knobs: (J_f^T J_f + mu I) dx_f = -J_f^T r
                std::vector<int> fr; for (int k = 0; k < K; ++k) if (!fixed[std::size_t(k)]) fr.push_back(k);
                const int m = int(fr.size());
                std::vector<double> A(std::size_t(m) * m, 0.0), b(std::size_t(m), 0.0);
                for (int a = 0; a < m; ++a) {
                    for (int c = 0; c < m; ++c) for (int i = 0; i < K; ++i) A[std::size_t(a * m + c)] += J[std::size_t(i)][std::size_t(fr[std::size_t(a)])] * J[std::size_t(i)][std::size_t(fr[std::size_t(c)])];
                    for (int i = 0; i < K; ++i) b[std::size_t(a)] -= J[std::size_t(i)][std::size_t(fr[std::size_t(a)])] * r[std::size_t(i)];
                    A[std::size_t(a * m + a)] *= 1.0 + 1e-6;
                }
                std::vector<double> sol(static_cast<std::size_t>(m), 0.0);
                if (m > 0) {
                    const auto Ai = inverse(A, m);
                    for (int a = 0; a < m; ++a) for (int c = 0; c < m; ++c) sol[std::size_t(a)] += Ai[std::size_t(a * m + c)] * b[std::size_t(c)];
                }
                dx = {};
                for (int a = 0; a < m; ++a) dx[std::size_t(fr[std::size_t(a)])] = sol[std::size_t(a)];
                bool changed = false;
                for (int k = 0; k < K; ++k) if (!fixed[std::size_t(k)]) {
                    if ((x[std::size_t(k)] <= lo[std::size_t(k)] + 1e-9 && dx[std::size_t(k)] < 0) || (x[std::size_t(k)] >= hi[std::size_t(k)] - 1e-9 && dx[std::size_t(k)] > 0)) { fixed[std::size_t(k)] = true; changed = true; }
                }
                if (!changed) break;
            }
            for (int k = 0; k < K; ++k) { const double lim = 0.3 * std::max(std::abs(x[std::size_t(k)]), 15.0); dx[std::size_t(k)] = std::clamp(dx[std::size_t(k)], -lim, lim); }
            const double n0 = norm(r);
            Vec xb = x, fb = f; double nb = n0;
            for (double alpha : {1.0, 0.5, 0.25}) {
                Vec xn = x; for (int k = 0; k < K; ++k) xn[std::size_t(k)] += alpha * dx[std::size_t(k)];
                xn = clampx(xn);
                const Vec fn = eval(xn, coarse);
                const double nn = norm(res(fn));
                if (nn < nb) { xb = xn; fb = fn; nb = nn; break; }
            }
            if (nb >= n0) break;
            x = xb; f = fb;
            line(("  gauss-newton " + std::to_string(it + 1)).c_str(), x, f);
        }
        if (pass == 0) {
            const Vec ff = eval(x, fine);
            line("fine-mesh check", x, ff);
            for (int i = 0; i < K; ++i) target[std::size_t(i)] *= f[std::size_t(i)] / ff[std::size_t(i)];
        }
    }
    put(S, x);
    const LayoutMetrics mf = evaluate_cells(S, fine);
    if (log) {
        log->final_fine = mf;
        const Vec ff = metric(mf);
        for (int i = 0; i < K; ++i) {
            log->final_rel[std::size_t(i)] = ff[std::size_t(i)] / target0[std::size_t(i)] - 1;
            log->max_rel_err = std::max(log->max_rel_err, std::abs(log->final_rel[std::size_t(i)]));
        }
        line("final (fine mesh)", x, ff);
    }
    return S;
}

// ---------------------------------------------------------------- chip assembly
namespace detail {
// axis-aligned CPW polyline: centre width w, gaps s; the far end is open (trench continues) or shorted to ground
inline void cpw_path(const std::vector<std::pair<double, double>>& P, double w, double s, std::vector<Rect>& metal, std::vector<Rect>& holes, bool short_end) {
    for (std::size_t i = 0; i + 1 < P.size(); ++i) {
        const auto [x0, y0] = P[i]; const auto [x1, y1] = P[i + 1];
        const bool last = i + 2 == P.size();
        const double a = w / 2, g = w / 2 + s;
        if (std::abs(y1 - y0) < 1e-9) {           // horizontal
            const double lo = std::min(x0, x1), hi = std::max(x0, x1);
            const bool end_hi = x1 > x0;
            double mlo = lo - a, mhi = hi + a, hlo = lo - g, hhi = hi + g;
            if (last && short_end) { if (end_hi) { hhi = hi; mhi = hi + s + 2; } else { hlo = lo; mlo = lo - s - 2; } }
            metal.push_back({mlo, y0 - a, mhi, y0 + a}); holes.push_back({hlo, y0 - g, hhi, y0 + g});
        } else {                                   // vertical
            const double lo = std::min(y0, y1), hi = std::max(y0, y1);
            const bool end_hi = y1 > y0;
            double mlo = lo - a, mhi = hi + a, hlo = lo - g, hhi = hi + g;
            if (last && short_end) { if (end_hi) { hhi = hi; mhi = hi + s + 2; } else { hlo = lo; mlo = lo - s - 2; } }
            metal.push_back({x0 - a, mlo, x0 + a, mhi}); holes.push_back({x0 - g, hlo, x0 + g, hhi});
        }
    }
}
// diagonal staircase CPW of 10 um squares from (x0,y0) to (x1,y1) (|dx| = |dy|)
inline void stair(double x0, double y0, double x1, double y1, std::vector<Rect>& metal, std::vector<Rect>& holes, bool short_end) {
    const int n = int(std::lround(std::abs(x1 - x0) / 5.0));
    const double sx = x1 > x0 ? 5 : -5, sy = y1 > y0 ? 5 : -5;
    for (int k = 0; k <= n; ++k) {
        const double cx = x0 + k * sx, cy = y0 + k * sy;
        metal.push_back(rect_c(cx, cy, 10, 10));
        if (!(short_end && k == n)) holes.push_back(rect_c(cx, cy, 22, 22));
    }
}
inline std::vector<Rect> subtract(const Rect& a, const Rect& b) {
    if (!a.intersects(b)) return {a};
    std::vector<Rect> out;
    if (b.y0 > a.y0) out.push_back({a.x0, a.y0, a.x1, b.y0});
    if (b.y1 < a.y1) out.push_back({a.x0, b.y1, a.x1, a.y1});
    const double y0 = std::max(a.y0, b.y0), y1 = std::min(a.y1, b.y1);
    if (b.x0 > a.x0) out.push_back({a.x0, y0, b.x0, y1});
    if (b.x1 < a.x1) out.push_back({b.x1, y0, a.x1, y1});
    return out;
}
} // namespace detail

struct ChipLayout {
    std::string name;
    Rect outline;
    std::vector<std::pair<std::string, std::vector<Rect>>> islands;   // named metal islands (qubits, couplers, buses, lines)
    std::vector<Rect> holes, jj, bumps, etch;
    std::vector<std::pair<std::string, std::pair<double, double>>> labels;
    struct Resonator { int qubit; double f_GHz, target_um, drawn_um; int rows; };
    std::vector<Resonator> resonators;
    std::vector<std::string> problems;
    double etch_area_mm2{0};
};

inline ChipLayout assemble_chip(const ChipDesign& D, const CellSet& cs, double eps_r = 11.45, double min_space_um = 10) {
    using detail::cpw_path; using detail::stair;
    ChipLayout L; L.name = D.name;
    const double p = cs.pitch, eps_eff = 0.5 * (eps_r + 1);
    std::map<int, int> idx;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) idx[D.qubits[i].id] = int(i);
    Rect bb{1e300, 1e300, -1e300, -1e300};
    for (const auto& q : D.qubits) { bb.x0 = std::min(bb.x0, q.x_um); bb.y0 = std::min(bb.y0, q.y_um); bb.x1 = std::max(bb.x1, q.x_um); bb.y1 = std::max(bb.y1, q.y_um); }
    L.outline = bb.grow(1100);
    for (const auto& q : D.qubits) {
        const bool meas = q.role != "data";
        const XmonParams& X = meas ? cs.meas : cs.data;
        const double cx = q.x_um, cy = q.y_um;
        std::vector<Rect> qm, qh;
        xmon(X, cx, cy, qm, qh, &L.jj);
        L.islands.push_back({"Q" + std::to_string(q.id), qm});
        L.holes.insert(L.holes.end(), qh.begin(), qh.end());
        L.labels.push_back({"Q" + std::to_string(q.id), {cx - 60, cy - 90}});
        // readout: claw + diagonal exit + lead + serpentine, shorted end with a flip-chip coupling pad
        std::vector<Rect> rm, rh;
        const double hw = X.w_arm / 2, b = hw + X.r_pad + X.claw_gap + X.claw_w, lo = hw + X.claw_gap + 10;
        rm.push_back(Rect{lo, X.r_pad + hw + X.claw_gap, b, b}.shift(cx, cy));
        rm.push_back(Rect{X.r_pad + hw + X.claw_gap, lo, b, X.r_pad + hw + X.claw_gap}.shift(cx, cy));
        const double s_end = 170;
        stair(cx + b - 5, cy + b - 5, cx + s_end, cy + s_end, rm, rh, false);
        const double f_r = q.ro.f_GHz;
        const double target = 299792458.0 / (4 * f_r * 1e9 * std::sqrt(eps_eff)) * 1e6;
        double used = (s_end - (b - 5)) * std::sqrt(2.0);
        std::vector<std::pair<double, double>> path{{cx + s_end, cy + s_end}};
        const double bx0 = cx + 330, bx1 = cx + 810, pitch_r = 36;
        path.push_back({bx0 - 20, cy + s_end});
        used += (bx0 - 20) - (cx + s_end);
        double y = cy + s_end, x = bx0 - 20;
        int rows = 0;
        bool right = true;
        const double pad_reserve = 60;
        while (used < target - pad_reserve) {
            const double xt = right ? bx1 : bx0 - 20;
            const double run = std::abs(xt - x);
            if (used + run >= target - pad_reserve) { const double rem = target - pad_reserve - used; x += right ? rem : -rem; path.push_back({x, y}); used = target - pad_reserve; break; }
            x = xt; path.push_back({x, y}); used += run;
            if (used + pitch_r >= target - pad_reserve) { const double rem = target - pad_reserve - used; y += rem; path.push_back({x, y}); used = target - pad_reserve; break; }
            y += pitch_r; path.push_back({x, y}); used += pitch_r; ++rows; right = !right;
        }
        // coupling pad then a short to ground, both going up from the last point: the serpentine ends on its top row,
        // so the space above is free (going sideways could run into the next plaquette's coupler bus)
        const auto [ex, ey] = path.back();
        cpw_path(path, 10, 6, rm, rh, false);
        const double pcx = ex, pcy = ey + 20;
        rm.push_back(rect_c(pcx, pcy, 36, 36)); rh.push_back(rect_c(pcx, pcy, 52, 52));
        std::vector<std::pair<double, double>> tail{{pcx, pcy}, {pcx, pcy + 40}};
        cpw_path(tail, 10, 6, rm, rh, true);
        used += 60;
        L.bumps.push_back(rect_c(pcx, pcy, 30, 30));
        L.resonators.push_back({q.id, f_r, target, used, rows});
        L.islands.push_back({"R" + std::to_string(q.id), rm});
        L.holes.insert(L.holes.end(), rh.begin(), rh.end());
        if (y > cy + 480) L.problems.push_back("resonator R" + std::to_string(q.id) + " exceeds its serpentine box");
        // flux line (inductive, shorted next to the SQUID) and XY line (open end at the (+x,-y) corner), both to bumps
        std::vector<Rect> fm, fh;
        stair(cx - 61, cy - 61, cx - 156, cy - 156, fm, fh, false);
        cpw_path({{cx - 156, cy - 156}, {cx - 156, cy - 250}}, 10, 6, fm, fh, false);
        std::reverse(fm.begin(), fm.end());
        // mark the line end next to the SQUID as shorted: drop the trench square closest to the qubit
        for (auto it = fh.begin(); it != fh.end(); ++it) if (std::abs(it->cx() - (cx - 61)) < 1e-6 && std::abs(it->cy() - (cy - 61)) < 1e-6) { fh.erase(it); break; }
        fm.push_back(rect_c(cx - 156, cy - 262, 36, 36)); fh.push_back(rect_c(cx - 156, cy - 262, 52, 52));
        L.bumps.push_back(rect_c(cx - 156, cy - 262, 30, 30));
        L.islands.push_back({"FQ" + std::to_string(q.id), fm});
        L.holes.insert(L.holes.end(), fh.begin(), fh.end());
        std::vector<Rect> xm, xh;
        stair(cx + 52, cy - 52, cx + 156, cy - 156, xm, xh, false);
        cpw_path({{cx + 156, cy - 156}, {cx + 250, cy - 156}}, 10, 6, xm, xh, false);
        xm.push_back(rect_c(cx + 262, cy - 156, 36, 36)); xh.push_back(rect_c(cx + 262, cy - 156, 52, 52));
        L.bumps.push_back(rect_c(cx + 262, cy - 156, 30, 30));
        L.islands.push_back({"XY" + std::to_string(q.id), xm});
        L.holes.insert(L.holes.end(), xh.begin(), xh.end());
    }
    int nc = 0;
    for (const auto& c : D.couplings) {
        const auto& A = D.qubits[std::size_t(idx[c.a])]; const auto& B = D.qubits[std::size_t(idx[c.b])];
        const auto* lo_q = &A; const auto* hi_q = &B;
        int k = 0;
        if (std::abs(A.y_um - B.y_um) < 1e-6) { k = 0; if (B.x_um < A.x_um) std::swap(lo_q, hi_q); }
        else { k = 1; if (B.y_um < A.y_um) std::swap(lo_q, hi_q); }
        if (std::abs(std::hypot(A.x_um - B.x_um, A.y_um - B.y_um) - p) > 1) { L.problems.push_back("coupler " + std::to_string(c.a) + "-" + std::to_string(c.b) + " is not a nearest-neighbour pair"); continue; }
        const XmonParams& Xl = lo_q->role != "data" ? cs.meas : cs.data;
        const XmonParams& Xh = hi_q->role != "data" ? cs.meas : cs.data;
        auto E = coupler_and_bus(Xl, cs.cg, p, k, lo_q->x_um, lo_q->y_um, true, 200, Xh.L);
        const std::string tag = std::to_string(c.a) + "-" + std::to_string(c.b);
        L.islands.push_back({"C" + tag, E.coupler});
        if (!E.bus.empty()) L.islands.push_back({"B" + tag, E.bus});
        L.holes.insert(L.holes.end(), E.holes.begin(), E.holes.end());
        // coupler SQUID on the pad's far side (local +x of the pad), flux line to a bump
        const double ycen = -cs.cg.y1 - cs.cg.w_body / 2 - cs.cg.pad / 2;
        const double xj0 = p / 2 + cs.cg.pad / 2, xj1 = xj0 + cs.cg.gap;
        L.jj.push_back(rot90({xj0, ycen - 6, xj1, ycen + 6}, k).shift(lo_q->x_um, lo_q->y_um));
        std::vector<Rect> fm, fh, fm_l, fh_l;
        cpw_path({{xj1 + 70, ycen}, {xj1 + 8, ycen}}, 10, 6, fm_l, fh_l, true);
        fm_l.push_back(rect_c(xj1 + 82, ycen, 36, 36)); fh_l.push_back(rect_c(xj1 + 82, ycen, 52, 52));
        for (const auto& r : fm_l) fm.push_back(rot90(r, k).shift(lo_q->x_um, lo_q->y_um));
        for (const auto& r : fh_l) fh.push_back(rot90(r, k).shift(lo_q->x_um, lo_q->y_um));
        L.bumps.push_back(rot90(rect_c(xj1 + 82, ycen, 30, 30), k).shift(lo_q->x_um, lo_q->y_um));
        L.islands.push_back({"FC" + tag, fm});
        L.holes.insert(L.holes.end(), fh.begin(), fh.end());
        const auto mid = rot90(rect_c(p / 2, ycen, 1, 1), k).shift(lo_q->x_um, lo_q->y_um);
        L.labels.push_back({"C" + tag, {mid.cx() - 40, mid.cy()}});
        ++nc;
    }
    // design-rule checks: islands must not overlap each other; all metal must lie inside the outline
    std::vector<std::pair<Rect, int>> all;
    for (int i = 0; i < int(L.islands.size()); ++i) for (const auto& r : L.islands[std::size_t(i)].second) all.push_back({r, i});
    std::sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.first.x0 < b.first.x0; });
    int clashes = 0;
    // designed coupling gaps (qubit <-> its couplers, buses and readout claw) only need the 4 um process minimum
    auto designed_gap = [&](const std::string& a, const std::string& b) {
        auto qid = [](const std::string& s) { return s.size() > 1 && s[0] == 'Q' ? std::stoi(s.substr(1)) : -1; };
        auto ends = [](const std::string& s) -> std::pair<int, int> {
            if (s.size() < 2 || (s[0] != 'C' && s[0] != 'B') || s[1] == 'Q') return {-1, -1};
            const auto d = s.find('-'); if (d == std::string::npos) return {-1, -1};
            return {std::stoi(s.substr(1, d - 1)), std::stoi(s.substr(d + 1))};
        };
        for (int pass = 0; pass < 2; ++pass) {
            const std::string& x = pass ? b : a; const std::string& y = pass ? a : b;
            const int q = qid(x);
            if (q < 0) continue;
            if (y == "R" + std::to_string(q)) return true;
            const auto [e0, e1] = ends(y);
            if (e0 == q || e1 == q) return true;
        }
        return false;
    };
    for (std::size_t i = 0; i < all.size(); ++i)
        for (std::size_t j = i + 1; j < all.size() && all[j].first.x0 < all[i].first.x1 + min_space_um; ++j) {
            if (all[i].second == all[j].second) continue;
            const double need = designed_gap(L.islands[std::size_t(all[i].second)].first, L.islands[std::size_t(all[j].second)].first) ? 4.0 : min_space_um;
            if (all[i].first.grow(need - 1e-6).intersects(all[j].first)) {
                if (clashes < 12) L.problems.push_back("metal spacing < " + std::to_string(int(need)) + " um between " + L.islands[std::size_t(all[i].second)].first + " and " + L.islands[std::size_t(all[j].second)].first);
                ++clashes;
            }
        }
    for (const auto& r : all) if (!L.outline.covers(r.first)) { L.problems.push_back("metal outside the chip outline"); break; }
    // etch = cut-outs minus metal islands
    std::vector<Rect> metal; for (const auto& r : all) metal.push_back(r.first);
    std::sort(metal.begin(), metal.end(), [](const Rect& a, const Rect& b) { return a.x0 < b.x0; });
    double area = 0;
    for (const auto& h : L.holes) {
        std::vector<Rect> pieces{h};
        auto it = std::lower_bound(metal.begin(), metal.end(), h.x0 - 2000.0, [](const Rect& a, double v) { return a.x0 < v; });
        for (; it != metal.end() && it->x0 < h.x1; ++it) {
            if (!it->intersects(h)) continue;
            std::vector<Rect> next;
            for (const auto& pc : pieces) { auto s = detail::subtract(pc, *it); next.insert(next.end(), s.begin(), s.end()); }
            pieces.swap(next);
        }
        for (const auto& pc : pieces) if (pc.w() > 1e-6 && pc.h() > 1e-6) { L.etch.push_back(pc); area += pc.area(); }
    }
    L.etch_area_mm2 = area * 1e-6;   // with overlaps between cut-outs counted twice (upper bound)
    (void)nc;
    return L;
}

inline void write_gds(const ChipLayout& L, const std::string& path, const std::string& cell = "TOP") {
    gds::Structure s; s.name = cell;
    auto poly = [](const Rect& r, int layer) { return gds::Polygon{layer, 0, {{r.x0, r.y0}, {r.x1, r.y0}, {r.x1, r.y1}, {r.x0, r.y1}}}; };
    for (const auto& r : L.etch) s.polys.push_back(poly(r, 1));
    for (const auto& r : L.jj) s.polys.push_back(poly(r, 2));
    for (const auto& r : L.bumps) s.polys.push_back(poly(r, 3));
    s.polys.push_back(poly(L.outline, 10));
    for (const auto& [t, xy] : L.labels) s.texts.push_back({20, xy.first, xy.second, t});
    gds::write(path, "QPUFORGE", {s});
}

// ---------------------------------------------------------------- targets from a chip design, cell I/O, post-layout check
inline LayoutTargets targets_from_design(const ChipDesign& D, const ChipReport& R) {
    double fm = 0, fd = 0, am = 0, ad = 0, fr = 0, gr = 0; int nm = 0, nd = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const bool meas = D.qubits[i].role != "data";
        (meas ? fm : fd) += R.q[i].f01; (meas ? am : ad) += R.q[i].alpha; (meas ? nm : nd)++;
        fr += D.qubits[i].ro.f_GHz; gr += D.qubits[i].ro.g_MHz * 1e-3;
    }
    fm /= std::max(nm, 1); fd /= std::max(nd, 1); am /= std::max(nm, 1); ad /= std::max(nd, 1);
    fr /= double(D.qubits.size()); gr /= double(D.qubits.size());
    double g12 = 0, gap = 0;
    const auto S = design_pair_system(fm, am, fd, ad, D.coupler, &g12, &gap);
    LayoutTargets T;
    T.EC_m = S.q1.EC; T.EC_c = S.c.EC; T.EC_d = S.q2.EC; T.G_mc = S.G1c; T.G_cd = S.G2c; T.G_md = S.G12;
    const double eps_eff = 0.5 * (11.45 + 1), ell = 299792458.0 / (4 * fr * 1e9 * std::sqrt(eps_eff)) * 1e6;
    const double c_l = cpw_capacitance_per_m(10e-6, 6e-6, 200e-6, 11.45) * 1e15 * 1e-6;       // fF/um
    const double Cr = ell * c_l / 2, Cs = 0.5 * (C_from_EC_F(T.EC_m) + C_from_EC_F(T.EC_d)) * 1e15;
    T.C_qr_fF = 2 * gr * std::sqrt(Cs * Cr) / std::sqrt(0.5 * (fm + fd) * fr);
    return T;
}

inline Json cells_json(const CellSet& c, const LayoutMetrics& m, const LayoutTargets& T, const SynthesisLog* log = nullptr) {
    Json j = Json::object();
    auto xm = [](const XmonParams& q) {
        Json o = Json::object();
        o.set("arm_length_um", q.L); o.set("arm_width_um", q.w_arm); o.set("gap_um", q.gap); o.set("tbar_thickness_um", q.t_w);
        o.set("tbar_length_um", q.T); o.set("readout_pad_um", q.r_pad); o.set("claw_gap_um", q.claw_gap); o.set("claw_width_um", q.claw_w);
        o.set("squid_stub_um", q.squid_stub);
        return o;
    };
    j.set("pitch_um", c.pitch); j.set("measure_qubit", xm(c.meas)); j.set("data_qubit", xm(c.data));
    Json g = Json::object();
    g.set("tbar_gap_um", c.cg.g_c); g.set("tbar_thickness_um", c.cg.t_w); g.set("offset_um", c.cg.y1); g.set("coupler_tbar_um", c.cg.l_c);
    g.set("body_width_um", c.cg.w_body); g.set("gap_um", c.cg.gap); g.set("pad_um", c.cg.pad); g.set("bus_tbar_um", c.cg.l_b);
    g.set("bus_width_um", c.cg.w_bus); g.set("bus_gap_um", c.cg.gap_bus);
    j.set("coupler", g);
    Json a = Json::object();
    a.set("EC_measure_GHz", m.EC_m); a.set("EC_data_GHz", m.EC_d); a.set("EC_coupler_GHz", m.EC_c);
    a.set("G_measure_coupler_GHz", m.G_mc); a.set("G_coupler_data_GHz", m.G_cd); a.set("G_measure_data_GHz", m.G_md);
    a.set("C_qubit_readout_fF", m.C_qr_fF); a.set("C_sigma_measure_fF", m.C_sigma_m); a.set("C_sigma_data_fF", m.C_sigma_d);
    a.set("C_coupler_fF", m.C_c); a.set("C_qubit_coupler_fF", m.C_qc); a.set("C_qubit_bus_fF", m.C_qb); a.set("C_bus_fF", m.C_bus);
    a.set("C_qubit_qubit_fF", m.C_qq); a.set("C_qubit_ground_fF", m.C_qground);
    j.set("achieved_fine_mesh", a);
    Json tt = Json::object();
    tt.set("EC_measure_GHz", T.EC_m); tt.set("EC_data_GHz", T.EC_d); tt.set("EC_coupler_GHz", T.EC_c); tt.set("G_measure_coupler_GHz", T.G_mc);
    tt.set("G_coupler_data_GHz", T.G_cd); tt.set("G_measure_data_GHz", T.G_md); tt.set("C_qubit_readout_fF", T.C_qr_fF);
    j.set("targets", tt);
    if (log) { Json l = Json::array(); for (const auto& s : log->lines) l.push(Json(s)); j.set("synthesis_log", l); j.set("max_rel_err", log->max_rel_err); }
    return j;
}
inline CellSet cells_from_json(const Json& j) {
    CellSet c;
    auto xm = [](const Json& o, XmonParams& q) {
        q.L = o.get("arm_length_um", q.L); q.w_arm = o.get("arm_width_um", q.w_arm); q.gap = o.get("gap_um", q.gap); q.t_w = o.get("tbar_thickness_um", q.t_w);
        q.T = o.get("tbar_length_um", q.T); q.r_pad = o.get("readout_pad_um", q.r_pad); q.claw_gap = o.get("claw_gap_um", q.claw_gap);
        q.claw_w = o.get("claw_width_um", q.claw_w); q.squid_stub = o.get("squid_stub_um", q.squid_stub);
    };
    c.pitch = j.get("pitch_um", c.pitch);
    xm(j["measure_qubit"], c.meas); xm(j["data_qubit"], c.data);
    const Json& g = j["coupler"];
    c.cg.g_c = g.get("tbar_gap_um", c.cg.g_c); c.cg.t_w = g.get("tbar_thickness_um", c.cg.t_w); c.cg.y1 = g.get("offset_um", c.cg.y1);
    c.cg.l_c = g.get("coupler_tbar_um", c.cg.l_c); c.cg.w_body = g.get("body_width_um", c.cg.w_body); c.cg.gap = g.get("gap_um", c.cg.gap);
    c.cg.pad = g.get("pad_um", c.cg.pad); c.cg.l_b = g.get("bus_tbar_um", c.cg.l_b); c.cg.w_bus = g.get("bus_width_um", c.cg.w_bus);
    c.cg.gap_bus = g.get("bus_gap_um", c.cg.gap_bus);
    return c;
}

// Post-layout check: rebuild the representative CZ pair from the extracted circuit (charging energies and charge
// couplings from the layout, junctions retargeted to the design frequencies) and report what changed.
struct PostLayout {
    double alpha_m{0}, alpha_d{0}, alpha_c{0};          // GHz
    double g_qc_MHz{0}, g12_MHz{0}, g12_null_MHz{0};    // at park
    double residual_gap_MHz{0}, zz_park_kHz{0};
    double spectator_LZ{0};                              // |11>-|20> transfer per spectator crossing (both passages)
    double g_r_MHz{0}, g_r_target_MHz{0};
    // re-nulling by parking the coupler below its sweet spot (possible when the layout's g12 exceeds the null value)
    double park_fc_GHz{0}, park_gap_MHz{0}, park_zz_kHz{0}, park_dfdphi_GHz{0};
};
inline PostLayout post_layout(const ChipDesign& D, const ChipReport& R, const LayoutMetrics& m) {
    double fm = 0, fd = 0, am = 0, ad = 0, fr = 0, grt = 0; int nm = 0, nd = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const bool meas = D.qubits[i].role != "data";
        (meas ? fm : fd) += R.q[i].f01; (meas ? am : ad) += R.q[i].alpha; (meas ? nm : nd)++;
        fr += D.qubits[i].ro.f_GHz; grt += D.qubits[i].ro.g_MHz;
    }
    fm /= std::max(nm, 1); fd /= std::max(nd, 1); am /= std::max(nm, 1); ad /= std::max(nd, 1);
    fr /= double(D.qubits.size()); grt /= double(D.qubits.size());
    PostLayout P;
    // design pair (null-designed g12) for reference
    double g12n = 0, gapn = 0;
    const auto Sd = design_pair_system(fm, am, fd, ad, D.coupler, &g12n, &gapn);
    P.g12_null_MHz = g12n * 1e3;
    // layout pair: same junction targets (f01), charging energies and couplings from the extraction
    auto retarget = [](double f01, double EC) { FluxTransmon t; t.EC = EC; t.EJsum = solve_EJ_for_f01(f01, EC); return t; };
    CZSystem S = Sd;
    S.q1 = retarget(fm, m.EC_m); S.q2 = retarget(fd, m.EC_d); S.c = retarget(D.coupler.f_max_GHz, m.EC_c);
    S.G1c = m.G_mc; S.G2c = m.G_cd; S.G12 = m.G_md; S.fluxc_idle = 0;
    P.alpha_m = S.q1.alpha(0); P.alpha_d = S.q2.alpha(0); P.alpha_c = S.c.alpha(0);
    const auto s1 = element_spectrum(S.q1, 0, 2), s2 = element_spectrum(S.q2, 0, 2), sc = element_spectrum(S.c, 0, 2);
    P.g_qc_MHz = S.G1c * s1.n(0, 1) * sc.n(0, 1) * 1e3;
    P.g12_MHz = S.G12 * s1.n(0, 1) * s2.n(0, 1) * 1e3;
    P.zz_park_kHz = S.zz(0.0) * 1e6;
    P.residual_gap_MHz = cz_avoided_crossing(S, 0.0).second * 1e3;
    const double v = 0.08;   // GHz/ns: typical mover sweep rate at a spectator crossing (7 ns ramps, ~400 MHz excursion)
    P.spectator_LZ = std::min(1.0, 2 * PI * PI * std::pow(P.residual_gap_MHz * 1e-3, 2) / v);
    // readout coupling
    const double eps_eff = 0.5 * (11.45 + 1), ell = 299792458.0 / (4 * fr * 1e9 * std::sqrt(eps_eff)) * 1e6;
    const double c_l = cpw_capacitance_per_m(10e-6, 6e-6, 200e-6, 11.45) * 1e15 * 1e-6;
    const double Cr = ell * c_l / 2, Cs = 0.5 * (m.C_sigma_m + m.C_sigma_d);
    P.g_r_MHz = 0.5 * m.C_qr_fF / std::sqrt(Cs * Cr) * std::sqrt(0.5 * (fm + fd) * fr) * 1e3;
    P.g_r_target_MHz = grt;
    // best park: minimise the residual |11>-|20> gap over coupler flux in [0, 0.12] Phi0 (golden section after a scan)
    auto gap_at = [&](double x) { return cz_avoided_crossing(S, x).second; };
    double bx = 0, bg = gap_at(0.0);
    for (int k = 1; k <= 24; ++k) { const double x = 0.005 * k, g = gap_at(x); if (g < bg) { bg = g; bx = x; } }
    double a = std::max(0.0, bx - 0.005), b = bx + 0.005;
    const double gr = (std::sqrt(5.0) - 1) / 2;
    double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = gap_at(x1), f2 = gap_at(x2);
    for (int it = 0; it < 30; ++it) {
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = gap_at(x1); }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = gap_at(x2); }
    }
    const double xp = 0.5 * (a + b);
    P.park_fc_GHz = S.c.f01(xp); P.park_gap_MHz = gap_at(xp) * 1e3; P.park_zz_kHz = S.zz(xp) * 1e6;
    P.park_dfdphi_GHz = (S.c.f01(xp + 1e-4) - S.c.f01(std::max(0.0, xp - 1e-4))) / (xp + 1e-4 - std::max(0.0, xp - 1e-4));
    return P;
}

} // namespace layout

} // namespace qpuforge
