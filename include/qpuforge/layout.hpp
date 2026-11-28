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
#include <queue>
#include <set>
#include <tuple>
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

// panel_charge (optional): the charge on every panel (fF x 1 V) with conductor `excite` at 1 V and the others at 0.
struct PanelCharge { Rect r; int cond; double q_fF; };
inline CapacitanceMatrix extract_capacitance(const std::vector<Conductor>& conds, double eps_r = 11.45, const MeshOptions& mo = {},
                                             int excite = -1, std::vector<PanelCharge>* panel_charge = nullptr) {
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
        if (panel_charge && k == excite)
            for (int i = 0; i < N; ++i) panel_charge->push_back({P[std::size_t(i)].r, P[std::size_t(i)].cond, b[std::size_t(i)] / Kc * 1e15});
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
    double l_b{30};         // bus T-bar length (sets the direct coupling); drawn no shorter than w_bus
    double w_bus{6}, gap_bus{15};
    double g_b{-1};         // gap between a bus end and the qubit T-bar face (< 0: g_c); retracted bus for a weak direct coupling
    double bus_gap() const { return g_b < 0 ? g_c : g_b; }
    double bus_tbar() const { return std::max(l_b, w_bus); }
};
// The synthesis drives the direct coupling with one continuous knob u: for u >= w_bus it is the bus T-bar length; below
// that the T-bar is gone (a bus no wider than its strip) and the bus end is pulled back from the qubit by
// kBusRetract * (w_bus - u). The two branches meet at u = w_bus, so no metal narrower than the bus strip is ever drawn.
inline constexpr double kBusRetract = 3.0;
inline double bus_knob(const CouplerGeom& c) { return c.l_b > c.w_bus ? c.l_b : c.w_bus - std::max(0.0, c.bus_gap() - c.g_c) / kBusRetract; }
inline void set_bus_knob(CouplerGeom& c, double u) {
    c.l_b = std::max(u, c.w_bus);
    c.g_b = c.g_c + kBusRetract * std::max(0.0, c.w_bus - u);
}

struct Drawn {
    std::vector<Rect> metal;     // islands inside cut-outs (qubit, coupler, bus, claws, centre conductors)
    std::vector<Rect> holes;     // ground cut-outs
    std::vector<Rect> jj;        // Josephson-junction / SQUID regions (GDS layer 2)
    std::vector<Rect> bumps;     // flip-chip bump pads (GDS layer 3)
    std::vector<std::pair<std::string, std::pair<double, double>>> labels;
};

// Qubit metal (cross + T-bars + readout pad + SQUID stub) and its cut-out, centred at (cx, cy).
namespace detail { inline std::vector<Rect> subtract(const Rect& a, const Rect& b); }
inline constexpr double kSquidGap = 10.0, kTongueW = 10.0;   // qubit SQUID span (stub to ground tongue) and tongue width (um)
inline void xmon(const XmonParams& q, double cx, double cy, std::vector<Rect>& metal, std::vector<Rect>& holes, std::vector<Rect>* jj = nullptr) {
    const double hw = q.w_arm / 2;
    const std::size_t first_hole = holes.size();
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
    // SQUID ground lead: a tongue of ground plane from the (-x,-y) corner pocket up beside the stub, kSquidGap from its
    // -x face, so the SQUID bridges stub -> tongue (island to ground). The tongue is ground: it is cut out of every
    // cut-out drawn so far. It stops kSquidGap below the -x arm.
    const double s0 = -hw - q.squid_stub;                                 // stub's outer corner (local)
    const Rect tongue = Rect{s0 - kSquidGap - kTongueW, -hw - q.gap - q.squid_stub, s0 - kSquidGap, -hw - kSquidGap}.shift(cx, cy);
    std::vector<Rect> cut;
    for (std::size_t i = first_hole; i < holes.size(); ++i) { auto pcs = detail::subtract(holes[i], tongue); cut.insert(cut.end(), pcs.begin(), pcs.end()); }
    holes.resize(first_hole); holes.insert(holes.end(), cut.begin(), cut.end());
    if (jj) jj->push_back(Rect{s0 - kSquidGap, s0 + 1, s0, s0 + 5}.shift(cx, cy));
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
    // bus (T-bar no shorter than the strip is wide; its ends may be retracted from the qubits, see set_bus_knob)
    if (c.l_b > 0) {
        const double lb = c.bus_tbar(), ua = q.L + c.bus_gap(), ub = p - (L_far > 0 ? L_far : q.L) - c.bus_gap();
        add(E.bus, {ua, c.y1 - lb / 2, ua + c.t_w, c.y1 + lb / 2});
        add(E.bus, {ua + c.t_w, c.y1 - c.w_bus / 2, far_end ? ub - c.t_w : xend, c.y1 + c.w_bus / 2});
        if (far_end) add(E.bus, {ub - c.t_w, c.y1 - lb / 2, ub, c.y1 + lb / 2});
        add(E.holes, {ua - c.gap_bus, c.y1 - lb / 2 - c.gap_bus, ua + c.t_w + c.gap_bus, c.y1 + lb / 2 + c.gap_bus});
        add(E.holes, {ua, c.y1 - c.w_bus / 2 - c.gap_bus, (far_end ? ub - c.t_w : xend) + (far_end ? 0 : c.gap_bus), c.y1 + c.w_bus / 2 + c.gap_bus});
        if (far_end) add(E.holes, {ub - c.t_w - c.gap_bus, c.y1 - lb / 2 - c.gap_bus, ub + c.gap_bus, c.y1 + lb / 2 + c.gap_bus});
    }
    return E;
}

// ---------------------------------------------------------------- extraction cells and synthesis
struct CellSet { XmonParams meas, data; CouplerGeom cg; double pitch{900}; };

namespace detail { inline int fill_ground_slivers(std::vector<Rect>& holes, const std::vector<Rect>& metal, double min_ground, int max_passes = 6); }
// Ground cut-outs of an extraction cell after the same mask clean-up as the assembled chip (ground slivers below the
// 4 um minimum merged away), so that the cells are extracted as they will be fabricated.
inline std::vector<Rect> cell_ground_holes(std::vector<Rect> holes, const std::vector<Conductor>& metal, double min_ground = 4.0) {
    std::vector<Rect> m;
    for (const auto& c : metal) if (!c.ground) m.insert(m.end(), c.rects.begin(), c.rects.end());
    detail::fill_ground_slivers(holes, m, min_ground);
    return holes;
}
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
    G.outline = {-ext, -ext, ext, ext}; G.holes = cell_ground_holes(holes, cs);
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
    G.outline = {-ex, -ex, p + ex, ex}; G.holes = cell_ground_holes(holes, cs);
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
    // The synthesis aims the direct coupling this fraction above G_md. Parking the coupler below its sweet spot can
    // only add coupler-mediated (negative) exchange, so a layout g12 slightly above the null value can be re-nulled and
    // one slightly below cannot; the margin keeps mesh error and the G_mc/G_cd asymmetry on the re-nullable side.
    double G_md_renull_margin{0.015};
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
    auto get = [](const CellSet& s) { return Vec{s.meas.L, s.data.L, s.cg.pad, s.cg.l_c, bus_knob(s.cg), s.meas.r_pad}; };
    auto put = [](CellSet& s, const Vec& x) {
        s.meas.L = x[0]; s.data.L = x[1]; s.cg.pad = x[2]; s.cg.l_c = x[3]; set_bus_knob(s.cg, x[4]); s.meas.r_pad = s.data.r_pad = x[5];
    };
    const double Tb = S.meas.T;
    auto bounds = [&](const Vec& x, Vec& lo, Vec& hi) {
        // bus knob: down to a 60 um retraction of the bus ends (see set_bus_knob)
        lo = {70, 70, 20, 8, S.cg.w_bus - 60 / kBusRetract, 12};
        hi = {300, 300, 420, Tb - 26, std::max(S.cg.w_bus, Tb - 26 - x[3]), 80};
        for (int q = 0; q < 2; ++q) {          // claw clear of the arm T-bars; resonator exit clear of the T-bar corners
            const XmonParams& Q = q ? S.data : S.meas;
            const double claw_outer = Q.w_arm / 2 + x[5] + Q.claw_gap + Q.claw_w;
            lo[std::size_t(q)] = std::max({lo[std::size_t(q)], Tb / 2 + 20, claw_outer + Q.t_w + 12, 80.0});
        }
    };
    auto clampx = [&](Vec x) { Vec lo, hi; bounds(x, lo, hi); for (int i = 0; i < K; ++i) x[std::size_t(i)] = std::clamp(x[std::size_t(i)], lo[std::size_t(i)], std::max(lo[std::size_t(i)], hi[std::size_t(i)])); return x; };
    const Vec w{1, 1, 1, 1, 0.7, 0.5};    // residual weights (G_md can be trimmed by the coupler park point; C_qr by the resonator)
    auto metric = [&](const LayoutMetrics& m) { return Vec{m.EC_m, m.EC_d, m.EC_c, 0.5 * (m.G_mc + m.G_cd), m.G_md, m.C_qr_fF}; };
    Vec target{T.EC_m, T.EC_d, T.EC_c, 0.5 * (T.G_mc + T.G_cd), T.G_md * (1 + T.G_md_renull_margin), T.C_qr_fF};
    const Vec target0{T.EC_m, T.EC_d, T.EC_c, 0.5 * (T.G_mc + T.G_cd), T.G_md, T.C_qr_fF};   // reported errors: against the design targets
    auto say = [&](const std::string& s) { if (log) log->lines.push_back(s); };
    auto line = [&](const char* tag, const Vec& x, const Vec& f) {
        char b[400];
        std::snprintf(b, sizeof b, "%-18s L_m %6.1f L_d %6.1f pad %6.1f l_c %5.1f bus %5.1f r %5.1f | EC_m %.4f EC_d %.4f EC_c %.4f G %.4f G_md %.5f C_qr %.2f",
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
// diagonal staircase CPW of 10 um squares from (x0,y0) to (x1,y1) (|dx| = |dy|). A shorted far end drops the trench of
// the last square and trims the one before it to the last square's outer edges, so the short meets the ground plane
// without leaving a 1 um etch slit beside it.
inline void stair(double x0, double y0, double x1, double y1, std::vector<Rect>& metal, std::vector<Rect>& holes, bool short_end) {
    const int n = int(std::lround(std::abs(x1 - x0) / 5.0));
    const double sx = x1 > x0 ? 5 : -5, sy = y1 > y0 ? 5 : -5;
    const double xn = x0 + n * sx, yn = y0 + n * sy;
    for (int k = 0; k <= n; ++k) {
        const double cx = x0 + k * sx, cy = y0 + k * sy;
        metal.push_back(rect_c(cx, cy, 10, 10));
        if (short_end && k == n) continue;
        Rect h = rect_c(cx, cy, 22, 22);
        if (short_end && k == n - 1) {
            if (sx > 0) h.x1 = std::min(h.x1, xn + 5); else h.x0 = std::max(h.x0, xn - 5);
            if (sy > 0) h.y1 = std::min(h.y1, yn + 5); else h.y0 = std::max(h.y0, yn - 5);
        }
        holes.push_back(h);
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

// Exact minimum-width / minimum-space check of a merged Manhattan region (union of rects), projection metric: in every
// horizontal strip between consecutive edge coordinates the region is a set of disjoint x-intervals; an interval shorter
// than min_width is a width violation and a gap between neighbouring intervals shorter than min_space is a space
// violation (likewise in vertical strips). Corner-to-corner (Euclidean) spacing is not checked. Runs in tiles with a halo;
// markers that continue across strips or tiles are merged, so one violation = one marker.
struct DrcMarker { std::string rule; Rect box; double value; };
inline std::vector<DrcMarker> width_space_check(const std::vector<Rect>& region, double min_width, double min_space,
                                                const std::string& name, double tile = 250) {
    std::vector<DrcMarker> raw;
    if (region.empty()) return raw;
    const Rect ext = bbox_of(region);
    const double halo = std::max(min_width, min_space) + 1, eps = 1e-7;
    for (int dir = 0; dir < 2; ++dir) {                           // dir 0: x-intervals in horizontal strips; 1: transposed
        std::vector<Rect> rs;
        for (const auto& r : region) rs.push_back(dir ? Rect{r.y0, r.x0, r.y1, r.x1} : r);
        const Rect E = dir ? Rect{ext.y0, ext.x0, ext.y1, ext.x1} : ext;
        std::sort(rs.begin(), rs.end(), [](const Rect& a, const Rect& b) { return a.y0 < b.y0; });
        const int ntx = std::max(1, int(std::ceil(E.w() / tile))), nty = std::max(1, int(std::ceil(E.h() / tile)));
        for (int ty = 0; ty < nty; ++ty) for (int tx = 0; tx < ntx; ++tx) {
            const Rect core{E.x0 + tx * tile, E.y0 + ty * tile, E.x0 + (tx + 1) * tile, E.y0 + (ty + 1) * tile};
            const Rect win = core.grow(halo);
            std::vector<Rect> w;
            for (const auto& r : rs) {
                if (r.y0 >= win.y1) break;
                if (!r.intersects(win)) continue;
                w.push_back({std::max(r.x0, win.x0), std::max(r.y0, win.y0), std::min(r.x1, win.x1), std::min(r.y1, win.y1)});
            }
            std::vector<double> ys;
            for (const auto& r : w) { ys.push_back(r.y0); ys.push_back(r.y1); }
            std::sort(ys.begin(), ys.end());
            ys.erase(std::unique(ys.begin(), ys.end(), [&](double a, double b) { return b - a < eps; }), ys.end());
            for (std::size_t k = 0; k + 1 < ys.size(); ++k) {
                const double ya = ys[k], yb = ys[k + 1], ym = 0.5 * (ya + yb);
                if (yb - ya < eps || ym < core.y0 || ym >= core.y1) continue;
                std::vector<std::pair<double, double>> iv;
                for (const auto& r : w) if (r.y0 < ym && r.y1 > ym) iv.push_back({r.x0, r.x1});
                std::sort(iv.begin(), iv.end());
                std::vector<std::pair<double, double>> m;
                for (const auto& v : iv) { if (!m.empty() && v.first <= m.back().second + eps) m.back().second = std::max(m.back().second, v.second); else m.push_back(v); }
                auto emit = [&](const char* rule, double a, double b) {
                    const double xm = 0.5 * (a + b);
                    if (xm < core.x0 || xm >= core.x1) return;
                    const Rect bx = dir ? Rect{ya, a, yb, b} : Rect{a, ya, b, yb};
                    raw.push_back({name + (dir ? ".v." : ".h.") + rule, bx, b - a});
                };
                for (std::size_t i = 0; i < m.size(); ++i) {
                    const bool clipped = m[i].first <= win.x0 + eps || m[i].second >= win.x1 - eps;
                    if (!clipped && m[i].second - m[i].first < min_width - eps) emit("width", m[i].first, m[i].second);
                    if (i + 1 < m.size() && m[i + 1].first - m[i].second < min_space - eps) emit("space", m[i].second, m[i + 1].first);
                }
            }
        }
    }
    // merge markers that continue into the next strip / tile (same rule and interval, touching boxes)
    std::sort(raw.begin(), raw.end(), [](const DrcMarker& a, const DrcMarker& b) {
        const bool va = a.rule.find(".v.") != std::string::npos;
        const double a0 = va ? a.box.y0 : a.box.x0, b0 = va ? b.box.y0 : b.box.x0;
        if (a.rule != b.rule) return a.rule < b.rule;
        if (std::abs(a0 - b0) > 1e-7) return a0 < b0;
        if (std::abs(a.value - b.value) > 1e-7) return a.value < b.value;
        return (va ? a.box.x0 : a.box.y0) < (va ? b.box.x0 : b.box.y0);
    });
    std::vector<DrcMarker> out;
    for (const auto& r : raw) {
        if (!out.empty()) {
            auto& o = out.back();
            const bool v = r.rule.find(".v.") != std::string::npos;
            const bool same = o.rule == r.rule && std::abs(o.value - r.value) < 1e-7 &&
                              std::abs((v ? o.box.y0 : o.box.x0) - (v ? r.box.y0 : r.box.x0)) < 1e-7;
            const bool touch = v ? r.box.x0 <= o.box.x1 + 1e-7 : r.box.y0 <= o.box.y1 + 1e-7;
            if (same && touch) { o.box.x1 = std::max(o.box.x1, r.box.x1); o.box.y1 = std::max(o.box.y1, r.box.y1); continue; }
        }
        out.push_back(r);
    }
    return out;
}
// Etch layer: the cut-outs minus the metal islands, as disjoint-per-cut-out rectangles. `area_um2` (if given) receives
// the summed area, with overlaps between cut-outs counted twice (an upper bound).
inline std::vector<Rect> etch_of(const std::vector<Rect>& holes, std::vector<Rect> metal, double* area_um2 = nullptr) {
    std::sort(metal.begin(), metal.end(), [](const Rect& a, const Rect& b) { return a.x0 < b.x0; });
    double wmax = 0; for (const auto& m : metal) wmax = std::max(wmax, m.w());
    std::vector<Rect> etch;
    double area = 0;
    for (const auto& h : holes) {
        std::vector<Rect> pieces{h};
        auto it = std::lower_bound(metal.begin(), metal.end(), h.x0 - wmax - 1.0, [](const Rect& a, double v) { return a.x0 < v; });
        for (; it != metal.end() && it->x0 < h.x1; ++it) {
            if (!it->intersects(h)) continue;
            std::vector<Rect> next;
            for (const auto& pc : pieces) { auto s = subtract(pc, *it); next.insert(next.end(), s.begin(), s.end()); }
            pieces.swap(next);
        }
        for (const auto& pc : pieces) if (pc.w() > 1e-6 && pc.h() > 1e-6) { etch.push_back(pc); area += pc.area(); }
    }
    if (area_um2) *area_um2 = area;
    return etch;
}

struct UnionFindLite {
    std::vector<int> p;
    explicit UnionFindLite(int n) : p(std::size_t(n)) { for (int i = 0; i < n; ++i) p[std::size_t(i)] = i; }
    int find(int i) { while (p[std::size_t(i)] != i) { p[std::size_t(i)] = p[std::size_t(p[std::size_t(i)])]; i = p[std::size_t(i)]; } return i; }
    bool unite(int a, int b) { a = find(a); b = find(b); if (a == b) return false; p[std::size_t(std::max(a, b))] = std::min(a, b); return true; }
};
// Connected components of the metal (the outline minus the etch), exact for Manhattan geometry. The chip is cut into
// vertical columns; in each column, every horizontal strip between consecutive edge coordinates holds a set of disjoint
// metal x-intervals. Intervals that overlap (by more than a rounding tolerance) in adjacent strips, or along a column
// boundary, are joined by union-find. Corner-to-corner contact does not connect. at(x, y) returns the component of the
// metal at a point, or -1 on etch or outside the outline.
struct MetalMap {
    Rect outline;
    double col_w{400};
    struct Strip { double y0, y1; int first, count; };
    struct Col { double x0, x1; std::vector<Strip> strips; };
    std::vector<Col> cols;
    std::vector<double> ix0, ix1, iy0, iy1;
    std::vector<int> comp;
    int n_comp{0};
    int at(double x, double y) const {
        if (!(x > outline.x0 && x < outline.x1 && y > outline.y0 && y < outline.y1)) return -1;
        const int c = std::clamp(int((x - outline.x0) / col_w), 0, int(cols.size()) - 1);
        const auto& S = cols[std::size_t(c)].strips;
        auto it = std::upper_bound(S.begin(), S.end(), y, [](double v, const Strip& st) { return v < st.y0; });
        if (it == S.begin()) return -1;
        const Strip& st = *(it - 1);
        if (y >= st.y1) return -1;
        int lo = st.first, hi = st.first + st.count;         // last interval with ix0 <= x
        while (hi - lo > 1) { const int m = (lo + hi) / 2; (ix0[std::size_t(m)] <= x ? lo : hi) = m; }
        if (st.count == 0 || ix0[std::size_t(lo)] > x || x >= ix1[std::size_t(lo)]) return -1;
        return comp[std::size_t(lo)];
    }
};
inline MetalMap metal_map(const Rect& outline, const std::vector<Rect>& etch, double col_w = 400) {
    MetalMap M; M.outline = outline; M.col_w = col_w;
    const double eps = 1e-7;
    const int nc = std::max(1, int(std::ceil(outline.w() / col_w)));
    std::vector<std::vector<Rect>> bucket(static_cast<std::size_t>(nc));
    for (const auto& r : etch) {
        const int a = std::clamp(int((r.x0 - outline.x0) / col_w), 0, nc - 1), b = std::clamp(int((r.x1 - outline.x0) / col_w), 0, nc - 1);
        for (int c = a; c <= b; ++c) bucket[std::size_t(c)].push_back(r);
    }
    std::vector<int> parent;
    auto find = [&](int i) { while (parent[std::size_t(i)] != i) { parent[std::size_t(i)] = parent[std::size_t(parent[std::size_t(i)])]; i = parent[std::size_t(i)]; } return i; };
    auto unite = [&](int a, int b) { a = find(a); b = find(b); if (a != b) parent[std::size_t(std::max(a, b))] = std::min(a, b); };
    M.cols.resize(std::size_t(nc));
    for (int c = 0; c < nc; ++c) {
        MetalMap::Col& C = M.cols[std::size_t(c)];
        C.x0 = outline.x0 + c * col_w; C.x1 = c + 1 == nc ? outline.x1 : std::min(outline.x1, C.x0 + col_w);
        auto& rs = bucket[std::size_t(c)];
        for (auto& r : rs) { r.x0 = std::max(r.x0, C.x0); r.x1 = std::min(r.x1, C.x1); r.y0 = std::max(r.y0, outline.y0); r.y1 = std::min(r.y1, outline.y1); }
        rs.erase(std::remove_if(rs.begin(), rs.end(), [&](const Rect& r) { return r.w() <= eps || r.h() <= eps; }), rs.end());
        std::sort(rs.begin(), rs.end(), [](const Rect& a, const Rect& b) { return a.y0 < b.y0; });
        std::vector<double> ys{outline.y0, outline.y1};
        for (const auto& r : rs) { ys.push_back(r.y0); ys.push_back(r.y1); }
        std::sort(ys.begin(), ys.end());
        ys.erase(std::unique(ys.begin(), ys.end(), [&](double a, double b) { return b - a < eps; }), ys.end());
        std::vector<const Rect*> active;
        std::size_t next = 0;
        int prev_first = -1, prev_count = 0;
        for (std::size_t k = 0; k + 1 < ys.size(); ++k) {
            const double ya = ys[k], yb = ys[k + 1], ym = 0.5 * (ya + yb);
            while (next < rs.size() && rs[next].y0 < ym) active.push_back(&rs[next++]);
            active.erase(std::remove_if(active.begin(), active.end(), [&](const Rect* r) { return r->y1 <= ym; }), active.end());
            std::vector<std::pair<double, double>> iv;
            for (const Rect* r : active) if (r->y0 < ym) iv.push_back({r->x0, r->x1});
            std::sort(iv.begin(), iv.end());
            const int first = int(M.ix0.size());
            double x = C.x0;
            auto emit = [&](double a, double b) {
                if (b - a <= eps) return;
                M.ix0.push_back(a); M.ix1.push_back(b); M.iy0.push_back(ya); M.iy1.push_back(yb);
                parent.push_back(int(parent.size()));
            };
            for (const auto& v : iv) { if (v.first > x + eps) emit(x, v.first); x = std::max(x, v.second); }
            emit(x, C.x1);
            const int count = int(M.ix0.size()) - first;
            // join with the strip below (they share the edge y = ya when ya is the previous strip's top)
            if (prev_first >= 0) {
                int i = prev_first, j = first;
                while (i < prev_first + prev_count && j < first + count) {
                    const double o = std::min(M.ix1[std::size_t(i)], M.ix1[std::size_t(j)]) - std::max(M.ix0[std::size_t(i)], M.ix0[std::size_t(j)]);
                    if (o > eps) unite(i, j);
                    (M.ix1[std::size_t(i)] < M.ix1[std::size_t(j)] ? i : j)++;
                }
            }
            C.strips.push_back({ya, yb, first, count});
            prev_first = first; prev_count = count;
        }
    }
    // join across column boundaries: intervals ending at the boundary in column c and starting there in column c + 1
    for (int c = 0; c + 1 < nc; ++c) {
        const MetalMap::Col& A = M.cols[std::size_t(c)]; const MetalMap::Col& B = M.cols[std::size_t(c + 1)];
        std::vector<int> ra, lb;
        for (const auto& st : A.strips) if (st.count && M.ix1[std::size_t(st.first + st.count - 1)] >= A.x1 - eps) ra.push_back(st.first + st.count - 1);
        for (const auto& st : B.strips) if (st.count && M.ix0[std::size_t(st.first)] <= B.x0 + eps) lb.push_back(st.first);
        std::size_t i = 0, j = 0;
        while (i < ra.size() && j < lb.size()) {
            const int a = ra[i], b = lb[j];
            const double o = std::min(M.iy1[std::size_t(a)], M.iy1[std::size_t(b)]) - std::max(M.iy0[std::size_t(a)], M.iy0[std::size_t(b)]);
            if (o > eps) unite(a, b);
            (M.iy1[std::size_t(a)] < M.iy1[std::size_t(b)] ? i : j)++;
        }
    }
    M.comp.resize(parent.size());
    std::map<int, int> ids;
    for (std::size_t i = 0; i < parent.size(); ++i) {
        const int r = find(int(i));
        auto it = ids.find(r);
        if (it == ids.end()) it = ids.emplace(r, int(ids.size())).first;
        M.comp[i] = it->second;
    }
    M.n_comp = int(ids.size());
    return M;
}

// Mask clean-up: a strip of ground plane narrower than `min_ground` between two cut-outs (projection metric, as in
// width_space_check) is etched away, i.e. the two cut-outs are merged. This happens where a diagonal staircase line
// passes the corner of a qubit cut-out. Only ground is removed: a gap that touches island metal is left alone (and
// still reported by the DRC), so no island changes shape and no net is shorted or split. Repeats until nothing changes
// (a fill can expose a new notch); returns the number of rectangles added to `holes`.
inline int fill_ground_slivers(std::vector<Rect>& holes, const std::vector<Rect>& metal, double min_ground, int max_passes) {
    std::vector<Rect> ms = metal;
    std::sort(ms.begin(), ms.end(), [](const Rect& a, const Rect& b) { return a.x0 < b.x0; });
    double wmax = 0; for (const auto& m : ms) wmax = std::max(wmax, m.w());
    auto touches_metal = [&](const Rect& b) {
        const Rect g = b.grow(1e-6);
        auto it = std::lower_bound(ms.begin(), ms.end(), g.x0 - wmax - 1.0, [](const Rect& a, double v) { return a.x0 < v; });
        for (; it != ms.end() && it->x0 < g.x1; ++it) if (it->intersects(g)) return true;
        return false;
    };
    int added = 0;
    for (int pass = 0; pass < max_passes; ++pass) {
        const auto etch = etch_of(holes, ms);
        int now = 0;
        for (const auto& mk : width_space_check(etch, 0.0, min_ground, "fill")) {
            if (mk.rule.find("space") == std::string::npos || touches_metal(mk.box)) continue;
            holes.push_back(mk.box); ++now;
        }
        added += now;
        if (!now) break;
    }
    return added;
}
// fill_ground_slivers with a connectivity guard: a fill may only remove ground, never cut ground off the plane (the
// metal component of the outline border). Every fill that borders metal which was part of the plane before the clean-up
// and is cut off from it after is undone, until nothing is cut off; the DRC then reports the narrow neck it leaves.
// Returns the fills kept and undone, and the metal map of the final cut-outs.
struct GuardedFill { int kept{0}, undone{0}; MetalMap map; };
inline GuardedFill guarded_fill(const Rect& outline, std::vector<Rect>& holes, const std::vector<Rect>& metal, double min_ground) {
    GuardedFill G;
    const std::size_t n_holes = holes.size();
    const auto map0 = metal_map(outline, etch_of(holes, metal));
    const int plane0 = map0.at(outline.x0 + 0.5, outline.y0 + 0.5);
    fill_ground_slivers(holes, metal, min_ground);
    std::vector<bool> keep(holes.size() - n_holes, true);
    for (int pass = 0; pass < 8; ++pass) {
        std::vector<Rect> hs(holes.begin(), holes.begin() + std::ptrdiff_t(n_holes));
        for (std::size_t k = 0; k < keep.size(); ++k) if (keep[k]) hs.push_back(holes[n_holes + k]);
        G.map = metal_map(outline, etch_of(hs, metal));
        const int plane1 = G.map.at(outline.x0 + 0.5, outline.y0 + 0.5);
        std::vector<char> detached(std::size_t(G.map.n_comp), 0), seen(std::size_t(G.map.n_comp), 0);
        for (std::size_t i = 0; i < G.map.comp.size(); ++i) {
            const int c = G.map.comp[i];
            if (seen[std::size_t(c)] || c == plane1) continue;
            seen[std::size_t(c)] = 1;
            if (map0.at(0.5 * (G.map.ix0[i] + G.map.ix1[i]), 0.5 * (G.map.iy0[i] + G.map.iy1[i])) == plane0) detached[std::size_t(c)] = 1;
        }
        int undone = 0;
        for (std::size_t k = 0; k < keep.size(); ++k) {
            if (!keep[k]) continue;
            const Rect& f = holes[n_holes + k];
            const double d = 0.25;
            const std::array<std::pair<double, double>, 4> probe{{{f.x0 - d, f.cy()}, {f.x1 + d, f.cy()}, {f.cx(), f.y0 - d}, {f.cx(), f.y1 + d}}};
            for (const auto& pr : probe) {
                const int c = G.map.at(pr.first, pr.second);
                if (c >= 0 && detached[std::size_t(c)]) { keep[k] = false; ++undone; break; }
            }
        }
        G.undone += undone;
        if (!undone) break;
    }
    std::vector<Rect> hs(holes.begin(), holes.begin() + std::ptrdiff_t(n_holes));
    for (std::size_t k = 0; k < keep.size(); ++k) if (keep[k]) hs.push_back(holes[n_holes + k]);
    holes.swap(hs);
    G.kept = int(holes.size() - n_holes);
    return G;
}

// Ground straps (modelled airbridges): a strap is a kStrapW-wide bridge whose two kStrapW x kStrapW feet sit on ground
// metal of two different metal components; it ties them together electrically without touching anything it spans.
// place_ground_straps ties every component that may be grounded -- floating ground fragments, and components holding
// only shorted lines (resonators, flux lines) with their pocket of ground -- to the ground plane, with the fewest and
// shortest straps (Kruskal over candidate straps, cost = span + a penalty per island crossed). A strap may span etch,
// ground, line centre conductors (R*, FQ*, FC*, XY*; penalty kStrapLinePenalty) and coupler bodies or bus strips (C*,
// B*; penalty kStrapCouplerPenalty), never a qubit. Couplers and buses must be crossable: each plaquette's ground and
// each strip between a coupler and its bus is enclosed by the coupler/bus bands between the qubits. A bridge over a
// coupler or bus adds roughly eps0 * kStrapW * width / (bridge height) of capacitance to it, which the extraction cells
// do not contain (about 0.2 fF over a 6 um bus at 3 um height). Feet never touch island metal; spans are at most
// kStrapMaxSpan.
inline constexpr double kStrapW = 10.0, kStrapMaxSpan = 150.0, kStrapLinePenalty = 40.0, kStrapCouplerPenalty = 80.0;
struct Strap { Rect span, foot_a, foot_b; int lines_crossed{0}, couplers_crossed{0}; };
// Partial self-inductance (pH) of a straight rectangular bar of length l, width w and thickness t (um): Grover's
// formula, L = 0.2 l [ln(2l/(w+t)) + 1/2 + 0.2235 (w+t)/l] nH with l in mm. A strap is modelled as such a bar from foot
// centre to foot centre; kStrapT is an assumed airbridge thickness (the result depends on it only logarithmically).
inline constexpr double kStrapT = 0.5;
inline double grover_bar_pH(double l_um, double w_um, double t_um = kStrapT) {
    if (l_um <= 0) return 0;
    const double s = w_um + t_um;
    return 0.2 * l_um * (std::log(2 * l_um / s) + 0.5 + 0.2235 * s / l_um);
}
inline double strap_inductance_pH(const Strap& st) {
    return grover_bar_pH(std::hypot(st.foot_a.cx() - st.foot_b.cx(), st.foot_a.cy() - st.foot_b.cy()), kStrapW);
}
struct IslandIndex {
    std::vector<std::pair<Rect, int>> rs;    // island rects sorted by x0, with the island index
    std::vector<char> kind;                  // 0: qubit (never bridged), 1: line, 2: coupler or bus
    double wmax{0};
    template <class F> void each(const Rect& q, F&& f) const {
        auto it = std::lower_bound(rs.begin(), rs.end(), q.x0 - wmax - 1.0, [](const std::pair<Rect, int>& a, double v) { return a.first.x0 < v; });
        for (; it != rs.end() && it->first.x0 < q.x1; ++it) if (it->first.intersects(q)) f(it->second);
    }
};
inline IslandIndex island_index(const std::vector<std::pair<std::string, std::vector<Rect>>>& islands) {
    IslandIndex I;
    for (int i = 0; i < int(islands.size()); ++i) {
        const std::string& n = islands[std::size_t(i)].first;
        I.kind.push_back(n[0] == 'Q' ? 0 : (n[0] == 'R' || n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n.rfind("XY", 0) == 0) ? 1 : 2);
        for (const auto& r : islands[std::size_t(i)].second) { I.rs.push_back({r, i}); I.wmax = std::max(I.wmax, r.w()); }
    }
    std::sort(I.rs.begin(), I.rs.end(), [](const auto& a, const auto& b) { return a.first.x0 < b.first.x0; });
    return I;
}
// component of the ground metal under a foot square, or -1 if the square is not wholly ground of one component
inline int foot_comp(const MetalMap& M, const IslandIndex& I, const Rect& f) {
    bool hit = false;
    I.each(f.grow(-1e-6), [&](int) { hit = true; });
    if (hit) return -1;
    const double e = 0.25;
    const std::array<double, 3> xs{f.x0 + e, f.cx(), f.x1 - e}, ys{f.y0 + e, f.cy(), f.y1 - e};
    int c = -2;
    for (double x : xs) for (double y : ys) {
        const int k = M.at(x, y);
        if (k < 0 || (c != -2 && k != c)) return -1;
        c = k;
    }
    return c;
}
// eligible[c]: component c may carry a strap foot (plane, floating ground, shorted-line pockets)
inline std::vector<Strap> place_ground_straps(const MetalMap& M, const IslandIndex& I, const std::vector<char>& eligible, int plane) {
    struct Cand { double cost; int a, b; Strap s; };
    std::map<std::pair<int, int>, Cand> best;
    std::vector<std::vector<std::size_t>> ivs(std::size_t(M.n_comp));
    for (std::size_t i = 0; i < M.comp.size(); ++i) ivs[std::size_t(M.comp[i])].push_back(i);
    const double F = kStrapW;
    const std::array<std::pair<int, int>, 4> dirs{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
    for (int T = 0; T < M.n_comp; ++T) {
        if (T == plane || !eligible[std::size_t(T)]) continue;
        // sample foot positions on the component's intervals (every 5 um along each, at most ~200 per component)
        std::vector<std::pair<double, double>> pts;
        for (std::size_t i : ivs[std::size_t(T)]) {
            if (M.ix1[i] - M.ix0[i] < F || M.iy1[i] - M.iy0[i] < 1.0) continue;
            for (double x = M.ix0[i] + F / 2; x <= M.ix1[i] - F / 2 + 1e-9; x += 5) pts.push_back({x, 0.5 * (M.iy0[i] + M.iy1[i])});
        }
        const std::size_t stride = std::max<std::size_t>(1, pts.size() / 200);
        for (std::size_t k = 0; k < pts.size(); k += stride) {
            // slide the foot vertically onto the component (the interval may be thinner than the foot)
            for (double dy : {0.0, -F / 2 + 0.5, F / 2 - 0.5}) {
                const double px = pts[k].first, py = pts[k].second + dy;
                const Rect fa = rect_c(px, py, F, F);
                if (foot_comp(M, I, fa) != T) continue;
                for (const auto& [dx, dyy] : dirs) {
                    int step = 0;
                    for (double t = F; t <= kStrapMaxSpan + F; t += 2.0, ++step) {
                        const double qx = px + dx * t, qy = py + dyy * t;
                        const int c0 = M.at(qx, qy);
                        if (c0 == T) continue;                   // still on the starting component (or back on it)
                        const Rect fb = rect_c(qx, qy, F, F);
                        if (c0 < 0 || !eligible[std::size_t(c0)] || foot_comp(M, I, fb) != c0) {
                            // spanning: every 10 um, give up this direction once the span would cross a qubit
                            if (step % 5 == 0) {
                                const Rect sp{std::min(fa.x0, fb.x0), std::min(fa.y0, fb.y0), std::max(fa.x1, fb.x1), std::max(fa.y1, fb.y1)};
                                bool bad = false; I.each(sp.grow(-1e-6), [&](int i) { if (!I.kind[std::size_t(i)]) bad = true; });
                                if (bad) break;
                            }
                            continue;
                        }
                        const int U = c0;
                        Strap s; s.foot_a = fa; s.foot_b = fb;
                        s.span = {std::min(fa.x0, fb.x0), std::min(fa.y0, fb.y0), std::max(fa.x1, fb.x1), std::max(fa.y1, fb.y1)};
                        std::set<int> crossed; I.each(s.span.grow(-1e-6), [&](int i) { crossed.insert(i); });
                        bool bad = false;
                        for (int i : crossed) { const int kd = I.kind[std::size_t(i)]; if (!kd) bad = true; else (kd == 1 ? s.lines_crossed : s.couplers_crossed)++; }
                        if (bad) break;
                        const double cost = (t - F) + kStrapLinePenalty * s.lines_crossed + kStrapCouplerPenalty * s.couplers_crossed;
                        const auto key = std::make_pair(std::min(T, U), std::max(T, U));
                        auto it = best.find(key);
                        if (it == best.end() || cost < it->second.cost) best[key] = {cost, T, U, s};
                        break;
                    }
                }
            }
        }
    }
    std::vector<Cand> cs; for (auto& [k, c] : best) cs.push_back(c);
    std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.cost < b.cost || (a.cost == b.cost && std::tie(a.a, a.b) < std::tie(b.a, b.b)); });
    UnionFindLite uf(M.n_comp);
    std::vector<Strap> out;
    for (const auto& c : cs) if (uf.unite(c.a, c.b)) out.push_back(c.s);
    return out;
}
} // namespace detail

// Flux line (FQ) and XY drive line (XY) of the qubit at (cx, cy), each from a bump pad. The flux line is shorted
// (short_flux) into the (-x,-y) corner pocket next to the SQUID; the pocket is enclosed by the arm and T-bar cut-outs and
// the ground straps (place_ground_straps) tie it to the plane. The XY line ends open below the qubit's (+x,-y) corner,
// outside its cut-out, at x = cx + xy_end, y = cy - 156: its coupling capacitance is set by xy_end (see line_coupling).
inline constexpr double kXYEnd = 247.0, kXYLevel = 279.0;
// centre line of the qubit flux line as qubit_control_lines draws it: bump pad -> straight run -> staircase -> short
inline std::vector<std::pair<double, double>> flux_line_centre(double cx, double cy) {
    return {{cx - 156, cy - 262}, {cx - 156, cy - 156}, {cx - 61, cy - 61}};
}
inline void qubit_control_lines(double cx, double cy, std::vector<Rect>& fm, std::vector<Rect>& fh, std::vector<Rect>& xm, std::vector<Rect>& xh,
                                bool short_flux, double xy_end = kXYEnd, double xy_level = kXYLevel) {
    using detail::cpw_path; using detail::stair;
    stair(cx - 156, cy - 156, cx - 61, cy - 61, fm, fh, short_flux);
    cpw_path({{cx - 156, cy - 156}, {cx - 156, cy - 250}}, 10, 6, fm, fh, false);
    fm.push_back(rect_c(cx - 156, cy - 262, 36, 36)); fh.push_back(rect_c(cx - 156, cy - 262, 52, 52));
    cpw_path({{cx + xy_end, cy - xy_level}, {cx + xy_end + 20, cy - xy_level}}, 10, 6, xm, xh, false);
    xm.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 36, 36)); xh.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 52, 52));
}
// Control-line coupling of one qubit: star-cell geometry (qubit, couplers and buses truncated, readout claw) plus its
// flux and XY lines up to their bump pads. The T1 limit of a coupling C_c into a line matched (Z0) at its far end:
// 1/T1 = omega^2 C_c^2 Z0 / C_sigma, with C_c the standing-wave weighted capacitance for the shorted flux line.
struct LineCoupling {
    double f_ref_GHz{0};        // frequency of C_flux_fF, T1_flux_us and T1_xy_us
    double C_flux_fF{0};        // effective (standing-wave weighted) qubit-flux-line capacitance at f_ref_GHz
    double C_flux_open_fF{0};   // plain capacitance to the same pieces (what an open line would couple)
    double C_xy_fF{0}, C_sigma_fF{0}, T1_flux_us{0}, T1_xy_us{0};
    double C_xy_coupler_fF{0};  // largest XY-line-to-coupler capacitance in the cell (couplers truncated: no far pads)
    // (path distance from the short in um, qubit-induced charge in fF) of every flux-line panel; the panel charges do not
    // depend on frequency, the sin(beta z) weight does, so the effective capacitance can be had at any frequency
    std::vector<std::pair<double, double>> flux_panels;
    double C_flux_at(double f_GHz, double eps_r = 11.45) const {
        const double beta = 2 * 3.14159265358979 * f_GHz * 1e9 * std::sqrt(0.5 * (eps_r + 1)) / 299792458.0 * 1e-6;   // rad / um
        double c = 0;
        for (const auto& [z, q] : flux_panels) c -= q * std::sin(beta * z);
        return c;
    }
    // T1 limit (us) of a line coupling C_c (fF) at f_GHz into Z0: 1/T1 = omega^2 C_c^2 Z0 / C_sigma
    double T1_us(double C_c_fF, double f_GHz, double Z0 = 50) const {
        const double w = 2 * 3.14159265358979 * f_GHz * 1e9;
        return C_c_fF > 0 ? C_sigma_fF * 1e-15 / (w * w * C_c_fF * 1e-15 * C_c_fF * 1e-15 * Z0) * 1e6 : 1e300;
    }
};
inline LineCoupling line_coupling(const XmonParams& q, const CouplerGeom& c, double p, double f_GHz, const MeshOptions& mo, double Z0 = 50,
                                  double xy_end = kXYEnd, double xy_level = kXYLevel) {
    std::vector<Conductor> cs;
    std::vector<Rect> holes;
    Conductor Q; Q.name = "qubit"; xmon(q, 0, 0, Q.rects, holes); cs.push_back(Q);
    for (int k = 0; k < 4; ++k) {
        auto E = coupler_and_bus(q, c, p, k, 0, 0, false);
        Conductor C; C.name = "coupler" + std::to_string(k); C.rects = E.coupler; cs.push_back(C);
        if (!E.bus.empty()) { Conductor B; B.name = "bus" + std::to_string(k); B.rects = E.bus; cs.push_back(B); }
        holes.insert(holes.end(), E.holes.begin(), E.holes.end());
    }
    Conductor R; R.name = "readout"; R.rects = readout_claw(q, 0, 0, 170, &holes); cs.push_back(R);
    // flux line drawn open and cut short by its last two squares (the short and the square joining it, both at ground
    // potential in the real line), so that it is a separate conductor
    std::vector<Rect> fm, xm;
    qubit_control_lines(0, 0, fm, holes, xm, holes, false, xy_end, xy_level);
    const int n = int(std::lround(95 / 5.0));                       // stair squares 0..n (see qubit_control_lines)
    Conductor F; F.name = "flux";
    for (int k = 0; k < int(fm.size()); ++k) if (k < n - 1 || k > n) F.rects.push_back(fm[std::size_t(k)]);
    cs.push_back(F);
    Conductor X; X.name = "xy"; X.rects = xm; cs.push_back(X);
    Conductor G; G.name = "ground"; G.ground = true;
    const double ext = q.L + c.g_c + c.t_w + 200 + 100;
    G.outline = {-ext, -ext, ext, ext}; G.holes = cell_ground_holes(holes, cs);
    cs.push_back(G);
    std::vector<PanelCharge> pq;
    const auto M = extract_capacitance(cs, 11.45, mo, 0, &pq);
    LineCoupling o;
    // a point at path distance z from the short of a line matched at its far end sees Re Z = Z0 sin^2(beta z): the
    // qubit's coupling to the shorted flux line is its capacitance to each piece of the line weighted by sin(beta z)
    // (panel charges with the qubit at 1 V; the weight is applied by C_flux_at)
    const double ls = 95 * std::sqrt(2.0);
    const int fi = M.index("flux");
    for (const auto& pc : pq) {
        if (pc.cond != fi) continue;
        const double x = pc.r.cx(), y = pc.r.cy();
        const double z = (y < -161 && x < -140) ? ls + (-156 - y) : std::max(0.0, ((-61 - x) + (-61 - y)) / std::sqrt(2.0));
        o.C_flux_open_fF -= pc.q_fF; o.flux_panels.push_back({z, pc.q_fF});
    }
    o.f_ref_GHz = f_GHz;
    o.C_flux_fF = o.C_flux_at(f_GHz);
    o.C_xy_fF = M.mutual("qubit", "xy"); o.C_sigma_fF = M.self("qubit");
    for (int k = 0; k < 4; ++k) o.C_xy_coupler_fF = std::max(o.C_xy_coupler_fF, M.mutual("coupler" + std::to_string(k), "xy"));
    o.T1_flux_us = o.T1_us(o.C_flux_fF, f_GHz, Z0); o.T1_xy_us = o.T1_us(o.C_xy_fF, f_GHz, Z0);
    return o;
}

// Control-line budget of a chip: the line couplings of its measure and data cells, and for every qubit the T1 limit
// of each line at the qubit's idle frequency against the qubit's design T1 (analyze()). The flux line's effective
// capacitance is re-weighted at each qubit's f01 (LineCoupling::C_flux_at); the cells are extracted once, at 5.5 GHz. The budget asks each line's
// limit to be at least kLineT1Margin times the design T1, so that the lines take at most ~10 % of the T1 budget each.
inline constexpr double kLineT1Margin = 10.0;
struct LineBudget {
    LineCoupling meas, data;
    double min_ratio_flux{0}, min_ratio_xy{0}, min_T1_flux_us{0}, min_T1_xy_us{0};
    int worst_flux_qubit{-1}, worst_xy_qubit{-1};
    bool ok() const { return min_ratio_flux >= kLineT1Margin && min_ratio_xy >= kLineT1Margin; }
};
inline LineBudget line_budget(const ChipDesign& D, const ChipReport& R, const CellSet& S, const MeshOptions& mo, double Z0 = 50) {
    LineBudget B;
    B.meas = line_coupling(S.meas, S.cg, S.pitch, 5.5, mo, Z0);
    B.data = line_coupling(S.data, S.cg, S.pitch, 5.5, mo, Z0);
    B.min_ratio_flux = B.min_ratio_xy = B.min_T1_flux_us = B.min_T1_xy_us = 1e300;
    for (std::size_t i = 0; i < D.qubits.size() && i < R.q.size(); ++i) {
        const LineCoupling& c = D.qubits[i].role != "data" ? B.meas : B.data;
        const double f = R.q[i].f01;
        const double tf = c.T1_us(c.C_flux_at(f), f, Z0), tx = c.T1_us(c.C_xy_fF, f, Z0);
        B.min_T1_flux_us = std::min(B.min_T1_flux_us, tf); B.min_T1_xy_us = std::min(B.min_T1_xy_us, tx);
        if (tf / R.q[i].T1_us < B.min_ratio_flux) { B.min_ratio_flux = tf / R.q[i].T1_us; B.worst_flux_qubit = D.qubits[i].id; }
        if (tx / R.q[i].T1_us < B.min_ratio_xy) { B.min_ratio_xy = tx / R.q[i].T1_us; B.worst_xy_qubit = D.qubits[i].id; }
    }
    return B;
}
inline std::string line_budget_summary(const LineBudget& B) {
    char b[600];
    std::snprintf(b, sizeof b, "control lines: XY %.3f/%.3f fF, flux %.3f/%.3f fF effective (%.2f/%.2f fF if open; measure/data); "
                  "T1 limits >= %.2f ms (XY), %.2f ms (flux): %.1fx (Q%d) and %.1fx (Q%d) the design T1 (budget %.0fx: %s)",
                  B.meas.C_xy_fF, B.data.C_xy_fF, B.meas.C_flux_fF, B.data.C_flux_fF, B.meas.C_flux_open_fF, B.data.C_flux_open_fF,
                  B.min_T1_xy_us / 1e3, B.min_T1_flux_us / 1e3, B.min_ratio_xy, B.worst_xy_qubit, B.min_ratio_flux, B.worst_flux_qubit,
                  kLineT1Margin, B.ok() ? "met" : "NOT met");
    return b;
}
inline Json line_budget_json(const LineBudget& B) {
    Json j = Json::object();
    j.set("C_xy_measure_fF", B.meas.C_xy_fF); j.set("C_xy_data_fF", B.data.C_xy_fF);
    j.set("C_flux_measure_fF", B.meas.C_flux_fF); j.set("C_flux_data_fF", B.data.C_flux_fF);
    j.set("C_flux_open_measure_fF", B.meas.C_flux_open_fF); j.set("C_flux_open_data_fF", B.data.C_flux_open_fF);
    j.set("min_T1_xy_us", B.min_T1_xy_us); j.set("min_T1_flux_us", B.min_T1_flux_us);
    j.set("min_ratio_xy", B.min_ratio_xy); j.set("min_ratio_flux", B.min_ratio_flux);
    j.set("worst_xy_qubit", B.worst_xy_qubit); j.set("worst_flux_qubit", B.worst_flux_qubit);
    j.set("margin", kLineT1Margin); j.set("ok", B.ok());
    return j;
}

// Flux-line to SQUID mutual inductance (model, see docs/design/mask-connectivity.md). M = Phi / I with Phi the flux of
// B_z through the SQUID loop from a current I along an in-plane filament path (Biot-Savart, closed form per straight
// segment, integrated on a grid over the loop). Assumptions, stated plainly: the SQUID loop is the drawn junction
// region (layer-2 rect of xmon); the current follows the flux line's centre line (flux_line_centre) to its short; the
// return current in the ground plane is not included (on a wide plane it spreads out; a return that hugs the loop
// would change M). Points within 0.05 um of a filament are skipped (a real line has width).
namespace detail {
// B_z (T per A) at (px, py) from a straight filament a -> b, coordinates in um
inline double segment_Bz_per_A(double ax, double ay, double bx, double by, double px, double py) {
    const double L = std::hypot(bx - ax, by - ay);
    if (L <= 0) return 0;
    const double ux = (bx - ax) / L, uy = (by - ay) / L;
    const double d = ux * (py - ay) - uy * (px - ax);          // signed perpendicular distance (um)
    if (std::abs(d) < 0.05) return 0;
    const double t = ux * (px - ax) + uy * (py - ay);
    return 1e-7 / (d * 1e-6) * ((L - t) / std::hypot(L - t, d) + t / std::hypot(t, d));
}
// flux (Wb per A) through rect `loop` (um) from a current along the polyline `path` (um)
inline double filament_flux_per_A(const std::vector<std::pair<double, double>>& path, const Rect& loop, int nx = 40, int ny = 40) {
    double phi = 0;
    const double dx = loop.w() / nx, dy = loop.h() / ny;
    for (int i = 0; i < nx; ++i) for (int j = 0; j < ny; ++j) {
        const double x = loop.x0 + (i + 0.5) * dx, y = loop.y0 + (j + 0.5) * dy;
        double b = 0;
        for (std::size_t k = 0; k + 1 < path.size(); ++k) b += segment_Bz_per_A(path[k].first, path[k].second, path[k + 1].first, path[k + 1].second, x, y);
        phi += b * dx * dy * 1e-12;
    }
    return phi;
}
} // namespace detail
inline constexpr double kPhi0 = 2.067833848e-15;   // Wb
struct FluxMutual { double M_pH{0}, I_half_mA{0}; Rect loop; };
inline FluxMutual flux_line_mutual(const XmonParams& q) {
    std::vector<Rect> m, h, jj;
    xmon(q, 0, 0, m, h, &jj);
    FluxMutual o; o.loop = jj.front();
    o.M_pH = std::abs(detail::filament_flux_per_A(flux_line_centre(0, 0), o.loop)) * 1e12;
    o.I_half_mA = o.M_pH > 0 ? kPhi0 / 2 / (o.M_pH * 1e-12) * 1e3 : 1e300;
    return o;
}

// Width / space rules on the fabricated etch layer (layer 1). Metal is the outline minus etch, so the etch-space rule is
// also the minimum metal (island and ground-plane) width, and the etch-width rule the minimum metal-to-metal gap.
// Default: the 4 um process minimum that assemble_chip already uses for designed coupling gaps.
struct DrcRules { double etch_min_width_um{4.0}, etch_min_space_um{4.0}; };
inline constexpr double kRemoveFragmentUm2 = 2500;   // unreachable ground fragments below this area are etched away

// Connectivity of the fabricated metal (LVS-lite), from detail::metal_map on the final etch layer:
//   - every island (qubit, coupler, bus, line) is its own conductor: no two islands share a metal component, and only
//     the lines meant to end in a short (resonators R*, qubit flux lines FQ*, coupler flux lines FC*) touch the plane;
//   - each shorted line reaches the ground plane (the component of the chip border), so it has a DC return;
//   - no ground metal is cut off from the plane (a floating ground fragment is an unintended island);
//   - each junction (layer 2) lands on metal at both ends of its long axis, on two different components.
struct ConnectivityReport {
    int components{0}, islands{0}, islands_merged{0}, islands_on_plane{0};
    int shorted_lines{0}, shorted_lines_on_plane{0}, flux_lines{0}, flux_lines_on_plane{0};
    int floating_ground{0}, junctions{0}, junctions_landed{0};
    int fills_reverted{0};                               // sliver fills undone because they cut ground off the plane (0: fills are no longer guarded)
    int straps{0}, straps_landed{0}, strap_lines_crossed{0}, strap_couplers_crossed{0};   // ground straps (modelled airbridges), see place_ground_straps
    int fragments_removed{0};                            // small unreachable ground fragments etched away
    // junction ownership: the island end of each junction is a qubit (Q*) or coupler (C*) island -- not a line, a bus
    // or a pad -- and every qubit and coupler island carries exactly one junction
    int junctions_on_device{0}, devices{0}, devices_one_junction{0};
    // RF ground (a strap is a ~70 pH airbridge, not an ideal ground): the ground end of a junction on the plane's own
    // metal (no strap in the junction's current path); shorted lines ending on the plane's own metal; components that
    // carry a junction's ground end or a line's short and reach the plane through straps, and how many of those hang on
    // a single strap (removing that one strap cuts them off); junction ground components that also terminate a flux line
    int junction_grounds_galvanic{0}, shorted_lines_galvanic{0}, flux_shorts_galvanic{0};
    int strapped_returns{0}, returns_single_strap{0}, shared_returns{0};
    // per junction (L.jj order): inductance (pH, strap chain to the plane, Grover) shared between the junction's current
    // and a flux line's return, i.e. when the junction's ground component also terminates a qubit or coupler flux line;
    // 0 when the junction grounds on the plane's own metal or its component terminates no flux line
    std::vector<double> junction_shared_pH;
    bool clean() const {
        return islands_merged == 0 && islands_on_plane == 0 && shorted_lines_on_plane == shorted_lines && floating_ground == 0 &&
               junctions_landed == junctions && straps_landed == straps && junctions_on_device == junctions && devices_one_junction == devices;
    }
    // RF grounding: every junction grounds on plane metal, and no line return or junction ground hangs on one strap.
    // Reported, not part of clean() yet: the committed masks fail it (docs/design/mask-connectivity.md, open items).
    bool rf_ground_ok() const { return junction_grounds_galvanic == junctions && returns_single_strap == 0; }
};
using detail::DrcMarker;

struct ChipLayout {
    std::string name;
    Rect outline;
    std::vector<std::pair<std::string, std::vector<Rect>>> islands;   // named metal islands (qubits, couplers, buses, lines)
    std::vector<Rect> holes, jj, bumps, etch;
    std::vector<detail::Strap> straps;                 // ground straps (GDS layer 4: span, layer 5: feet)
    std::vector<std::pair<std::string, std::pair<double, double>>> labels;
    struct Resonator { int qubit; double f_GHz, target_um, drawn_um; int rows; };
    std::vector<Resonator> resonators;
    std::vector<std::string> problems;                 // placement / spacing / outline problems found during assembly
    std::vector<DrcMarker> drc;                         // width / space violations on the merged etch layer
    DrcRules rules;
    int slivers_filled{0};                              // ground-sliver fills added by the mask clean-up (fill_ground_slivers)
    ConnectivityReport conn;                            // metal connectivity of the final mask (see ConnectivityReport)
    double etch_area_mm2{0};
    int drc_count(const std::string& rule_suffix) const {
        int n = 0; for (const auto& m : drc) if (m.rule.size() >= rule_suffix.size() && m.rule.compare(m.rule.size() - rule_suffix.size(), rule_suffix.size(), rule_suffix) == 0) ++n;
        return n;
    }
    bool drc_clean() const { return problems.empty() && drc.empty(); }
};

namespace detail {
inline bool is_shorted_line(const std::string& n) { return n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n[0] == 'R'; }
inline bool is_flux_line(const std::string& n) { return n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0; }
// metal component of each island (at the centre of its largest rect)
inline std::vector<int> island_components(const ChipLayout& L, const MetalMap& M) {
    std::vector<int> v;
    for (const auto& [name, rs] : L.islands) {
        const Rect* big = &rs.front();
        for (const auto& r : rs) if (r.area() > big->area()) big = &r;
        v.push_back(M.at(big->cx(), big->cy()));
    }
    return v;
}
} // namespace detail

// Connectivity of the final mask (ConnectivityReport) from its metal map, the islands, the junctions and the straps.
// Keeps fills_reverted and fragments_removed (set by the clean-up in assemble_chip). Can be rerun on a modified layout.
inline void compute_connectivity(ChipLayout& L, const detail::MetalMap& map1) {
    ConnectivityReport C;
    C.fills_reverted = L.conn.fills_reverted; C.fragments_removed = L.conn.fragments_removed;
    C.components = map1.n_comp;
    // nets: metal components joined by the ground straps
    detail::UnionFindLite uf(map1.n_comp);
    struct Edge { int a, b; double pH; };
    std::vector<Edge> edges;
    for (const auto& st : L.straps) {
        const int a = map1.at(st.foot_a.cx(), st.foot_a.cy()), b = map1.at(st.foot_b.cx(), st.foot_b.cy());
        C.straps++; C.strap_lines_crossed += st.lines_crossed; C.strap_couplers_crossed += st.couplers_crossed;
        if (a >= 0 && b >= 0 && a != b) { C.straps_landed++; uf.unite(a, b); edges.push_back({a, b, detail::strap_inductance_pH(st)}); }
    }
    const int plane0 = map1.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    auto net = [&](int c) { return c < 0 ? -1 : uf.find(c); };
    const int plane = net(plane0);
    std::map<int, int> island_of;                  // net -> number of islands in it
    const auto raw = detail::island_components(L, map1);
    std::vector<int> isl_comp;
    for (int c : raw) { isl_comp.push_back(net(c)); island_of[net(c)]++; }
    std::map<int, std::vector<int>> island_at_comp;   // raw component -> islands whose largest rect is on it
    for (std::size_t i = 0; i < raw.size(); ++i) if (raw[i] >= 0) island_at_comp[raw[i]].push_back(int(i));
    std::set<int> line_comps, flux_comps;           // raw components holding a shorted line / a flux line
    C.islands = int(L.islands.size());
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const std::string& n = L.islands[i].first;
        const bool fq = n.rfind("FQ", 0) == 0, line = detail::is_shorted_line(n);
        const int c = isl_comp[i];
        if (line) {
            C.shorted_lines++; if (c == plane) C.shorted_lines_on_plane++;
            if (raw[i] == plane0) C.shorted_lines_galvanic++;
            if (raw[i] >= 0) line_comps.insert(raw[i]);
            if (detail::is_flux_line(n) && raw[i] >= 0) flux_comps.insert(raw[i]);
        }
        else if (c == plane) C.islands_on_plane++;
        if (fq) { C.flux_lines++; if (c == plane) C.flux_lines_on_plane++; if (raw[i] == plane0) C.flux_shorts_galvanic++; }
        // a conductor (qubit, coupler, bus, XY line) sharing its metal with any other island is a short; two shorted
        // lines ending in the same enclosed pocket count under shorted lines that miss the plane instead
        if (!line && island_of[c] > 1) C.islands_merged++;
    }
    std::set<int> nets; for (int c = 0; c < map1.n_comp; ++c) nets.insert(net(c));
    for (int n : nets) if (n != plane && !island_of.count(n)) C.floating_ground++;
    // a junction (qubit or coupler SQUID: grounded transmons) lands on an island at one end and on the plane at the
    // other; the island end must be a qubit or coupler island
    std::map<int, int> jj_per_island;
    std::vector<int> jj_ground;                     // raw component of each junction's ground end (-1 if not landed)
    for (const auto& j : L.jj) {
        C.junctions++;
        const bool horiz = j.w() >= j.h();
        const double d = 0.5;
        const int ra = horiz ? map1.at(j.x0 - d, j.cy()) : map1.at(j.cx(), j.y0 - d);
        const int rb = horiz ? map1.at(j.x1 + d, j.cy()) : map1.at(j.cx(), j.y1 + d);
        const int a = net(ra), b = net(rb);
        int g = -1;
        if (a >= 0 && b >= 0 && a != b && (a == plane || b == plane)) {
            C.junctions_landed++;
            const int dev = a == plane ? rb : ra;
            g = a == plane ? ra : rb;
            bool own = false;
            if (auto it = island_at_comp.find(dev); it != island_at_comp.end() && it->second.size() == 1) {
                const std::string& n = L.islands[std::size_t(it->second.front())].first;
                if (n[0] == 'Q' || n[0] == 'C') { own = true; jj_per_island[it->second.front()]++; }
            }
            if (own) C.junctions_on_device++;
            if (g == plane0) C.junction_grounds_galvanic++;
        }
        jj_ground.push_back(g);
    }
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const char k = L.islands[i].first[0];
        if (k != 'Q' && k != 'C') continue;
        C.devices++;
        if (jj_per_island[int(i)] == 1) C.devices_one_junction++;
    }
    // strapped returns: components (not the plane's own metal) holding a junction's ground end or a line's short
    std::set<int> need;
    for (int g : jj_ground) if (g >= 0 && g != plane0) need.insert(g);
    for (int c : line_comps) if (c != plane0) need.insert(c);
    C.strapped_returns = int(need.size());
    std::set<int> fragile;
    for (std::size_t e = 0; e < edges.size() && !need.empty(); ++e) {
        detail::UnionFindLite u2(map1.n_comp);
        for (std::size_t f = 0; f < edges.size(); ++f) if (f != e) u2.unite(edges[f].a, edges[f].b);
        for (int c : need) if (u2.find(c) != u2.find(plane0)) fragile.insert(c);
    }
    C.returns_single_strap = int(fragile.size());
    // shared inductance: shortest strap chain (sum of strap inductances) from each component to the plane's own metal
    std::vector<double> dist(std::size_t(map1.n_comp), 1e300);
    if (plane0 >= 0) {
        std::vector<std::vector<std::pair<int, double>>> adj(std::size_t(map1.n_comp));
        for (const auto& e : edges) { adj[std::size_t(e.a)].push_back({e.b, e.pH}); adj[std::size_t(e.b)].push_back({e.a, e.pH}); }
        using QE = std::pair<double, int>;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        dist[std::size_t(plane0)] = 0; pq.push({0, plane0});
        while (!pq.empty()) {
            const auto [dd, u] = pq.top(); pq.pop();
            if (dd > dist[std::size_t(u)]) continue;
            for (const auto& [v, w] : adj[std::size_t(u)]) if (dd + w < dist[std::size_t(v)]) { dist[std::size_t(v)] = dd + w; pq.push({dd + w, v}); }
        }
    }
    std::set<int> shared;
    for (int g : jj_ground) {
        double ls = 0;
        if (g >= 0 && g != plane0 && flux_comps.count(g)) { ls = dist[std::size_t(g)]; shared.insert(g); }
        C.junction_shared_pH.push_back(ls);
    }
    C.shared_returns = int(shared.size());
    L.conn = C;
}
// One line on the RF ground (see ConnectivityReport::rf_ground_ok)
inline std::string rf_ground_summary(const ChipLayout& L) {
    const auto& C = L.conn;
    double worst = 0; for (double v : C.junction_shared_pH) worst = std::max(worst, v);
    char b[600];
    std::snprintf(b, sizeof b, "RF ground %s: junctions grounding on plane metal %d/%d; flux shorts on plane metal %d/%d; strapped returns %d "
                  "(%d on a single strap); junction grounds shared with a flux return %d (worst shared strap chain %.0f pH); junctions on "
                  "their own qubit/coupler %d/%d",
                  C.rf_ground_ok() ? "ok" : "NOT ok", C.junction_grounds_galvanic, C.junctions, C.flux_shorts_galvanic, C.flux_lines,
                  C.strapped_returns, C.returns_single_strap, C.shared_returns, worst, C.junctions_on_device, C.junctions);
    return b;
}

// XY-line crosstalk from a window of the assembled mask around one qubit and its XY line: the chip's own cut-outs
// (after the sliver fill) and islands, shorted lines treated as ground, extracted as is. ratio = the largest
// capacitance from the XY line to any other conductor (couplers, buses, other qubits, other lines) over its
// capacitance to its own qubit.
struct XYCrosstalk { double C_own_fF{0}, C_worst_fF{0}, ratio{0}; std::string worst; int panels{0}; double seconds{0}; };
inline XYCrosstalk xy_window_crosstalk(const ChipLayout& L, int qubit_id, const MeshOptions& mo = MeshOptions{3.0, 50.0, 2.0}) {
    const std::string Q = "Q" + std::to_string(qubit_id), XY = "XY" + std::to_string(qubit_id);
    Rect win{1e300, 1e300, -1e300, -1e300};
    for (const auto& [n, rs] : L.islands) if (n == Q || n == XY) for (const auto& r : rs) win = {std::min(win.x0, r.x0), std::min(win.y0, r.y0), std::max(win.x1, r.x1), std::max(win.y1, r.y1)};
    if (win.x0 > win.x1) throw std::runtime_error("xy_window_crosstalk: no " + Q + " / " + XY);
    win = win.grow(150);
    auto clip = [&](const Rect& r, Rect& o) {
        o = {std::max(r.x0, win.x0), std::max(r.y0, win.y0), std::min(r.x1, win.x1), std::min(r.y1, win.y1)};
        return o.x1 > o.x0 + 1e-9 && o.y1 > o.y0 + 1e-9;
    };
    std::vector<Conductor> cs;
    std::vector<Rect> holes;
    for (const auto& hr : L.holes) { Rect o; if (clip(hr, o)) holes.push_back(o); }
    for (const auto& [name, rs] : L.islands) {
        if (detail::is_shorted_line(name)) {               // ground metal: carve it out of the cut-outs
            for (const auto& r : rs) {
                Rect o; if (!clip(r, o)) continue;
                std::vector<Rect> nh;
                for (const auto& hr : holes) { auto pcs = detail::subtract(hr, o); nh.insert(nh.end(), pcs.begin(), pcs.end()); }
                holes.swap(nh);
            }
            continue;
        }
        Conductor c; c.name = name;
        for (const auto& r : rs) { Rect o; if (clip(r, o)) c.rects.push_back(o); }
        if (!c.rects.empty()) cs.push_back(c);
    }
    Conductor G; G.name = "ground"; G.ground = true; G.outline = win; G.holes = holes;
    cs.push_back(G);
    const auto X = extract_capacitance(cs, 11.45, mo);
    XYCrosstalk o; o.panels = X.panels; o.seconds = X.seconds;
    o.C_own_fF = X.mutual(XY, Q);
    for (const auto& n : X.names) {
        if (n == XY || n == Q || n == "ground") continue;
        const double c = X.mutual(XY, n);
        if (c > o.C_worst_fF) { o.C_worst_fF = c; o.worst = n; }
    }
    o.ratio = o.C_own_fF > 0 ? o.C_worst_fF / o.C_own_fF : 1e300;
    return o;
}

inline ChipLayout assemble_chip(const ChipDesign& D, const CellSet& cs, double eps_r = 11.45, double min_space_um = 10, const DrcRules& rules = {}) {
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
        std::vector<Rect> fm, fh, xm, xh;
        qubit_control_lines(cx, cy, fm, fh, xm, xh, true);
        L.bumps.push_back(rect_c(cx - 156, cy - 262, 30, 30));
        L.islands.push_back({"FQ" + std::to_string(q.id), fm});
        L.holes.insert(L.holes.end(), fh.begin(), fh.end());
        L.bumps.push_back(rect_c(cx + kXYEnd + 32, cy - kXYLevel, 30, 30));
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
    // etch = cut-outs minus metal islands, after merging cut-outs separated by ground slivers below the minimum width
    std::vector<Rect> metal; for (const auto& r : all) metal.push_back(r.first);
    // Sliver fills are not guarded here: ground they cut off the plane is strapped back (or removed) below.
    L.slivers_filled = detail::fill_ground_slivers(L.holes, metal, rules.etch_min_space_um);
    double area = 0;
    L.etch = detail::etch_of(L.holes, metal, &area);
    const auto I = detail::island_index(L.islands);
    auto shorted = [](const std::string& n) { return n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n[0] == 'R'; };
    auto island_comps = [&](const detail::MetalMap& M) {
        std::vector<int> v;
        for (const auto& [name, rs] : L.islands) {
            const Rect* big = &rs.front();
            for (const auto& r : rs) if (r.area() > big->area()) big = &r;
            v.push_back(M.at(big->cx(), big->cy()));
        }
        return v;
    };
    // ground straps; ground fragments no strap can reach and smaller than kRemoveFragmentUm2 are etched away instead
    detail::MetalMap map1 = detail::metal_map(L.outline, L.etch);
    for (int round = 0; round < 3; ++round) {
        const int plane = map1.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
        const auto ic = island_comps(map1);
        std::vector<char> eligible(std::size_t(map1.n_comp), 1), has_island(std::size_t(map1.n_comp), 0);
        for (std::size_t i = 0; i < ic.size(); ++i) if (ic[i] >= 0) {
            has_island[std::size_t(ic[i])] = 1;
            if (!shorted(L.islands[i].first)) eligible[std::size_t(ic[i])] = 0;
        }
        L.straps = detail::place_ground_straps(map1, I, eligible, plane);
        detail::UnionFindLite uf(map1.n_comp);
        for (const auto& st : L.straps) uf.unite(map1.at(st.foot_a.cx(), st.foot_a.cy()), map1.at(st.foot_b.cx(), st.foot_b.cy()));
        std::vector<double> carea(std::size_t(map1.n_comp), 0);
        for (std::size_t i = 0; i < map1.comp.size(); ++i) carea[std::size_t(map1.comp[i])] += (map1.ix1[i] - map1.ix0[i]) * (map1.iy1[i] - map1.iy0[i]);
        std::vector<Rect> removed;
        for (std::size_t i = 0; i < map1.comp.size(); ++i) {
            const int c = map1.comp[i];
            if (has_island[std::size_t(c)] || uf.find(c) == uf.find(plane) || carea[std::size_t(c)] >= kRemoveFragmentUm2) continue;
            removed.push_back({map1.ix0[i], map1.iy0[i], map1.ix1[i], map1.iy1[i]});
        }
        if (removed.empty()) break;
        std::set<int> rc; for (const auto& r : removed) rc.insert(map1.at(r.cx(), r.cy()));
        L.conn.fragments_removed += int(rc.size());
        L.holes.insert(L.holes.end(), removed.begin(), removed.end());
        L.etch = detail::etch_of(L.holes, metal, &area);
        map1 = detail::metal_map(L.outline, L.etch);
    }
    // connectivity of the final mask
    compute_connectivity(L, map1);
    L.etch_area_mm2 = area * 1e-6;   // with overlaps between cut-outs counted twice (upper bound)
    // minimum width / space on the merged etch layer (what the mask actually contains)
    L.rules = rules;
    L.drc = detail::width_space_check(L.etch, rules.etch_min_width_um, rules.etch_min_space_um, "etch");
    (void)nc;
    return L;
}

inline std::string connectivity_summary(const ChipLayout& L) {
    const auto& C = L.conn;
    char b[800];
    std::snprintf(b, sizeof b, "metal connectivity %s: %d components; islands merged %d, on the plane %d; shorted lines reaching the plane %d/%d "
                  "(qubit flux lines %d/%d); floating ground fragments %d; junctions landing island-to-plane %d/%d; ground straps %d "
                  "(landed %d; crossing %d lines, %d coupler/bus strips); small fragments removed %d",
                  C.clean() ? "clean" : "NOT clean", C.components, C.islands_merged, C.islands_on_plane, C.shorted_lines_on_plane, C.shorted_lines,
                  C.flux_lines_on_plane, C.flux_lines, C.floating_ground, C.junctions_landed, C.junctions, C.straps, C.straps_landed,
                  C.strap_lines_crossed, C.strap_couplers_crossed, C.fragments_removed);
    return b;
}
inline Json connectivity_json(const ChipLayout& L) {
    const auto& C = L.conn;
    Json j = Json::object();
    j.set("clean", C.clean()); j.set("components", C.components); j.set("islands", C.islands); j.set("islands_merged", C.islands_merged);
    j.set("islands_on_plane", C.islands_on_plane); j.set("shorted_lines", C.shorted_lines); j.set("shorted_lines_on_plane", C.shorted_lines_on_plane);
    j.set("flux_lines", C.flux_lines); j.set("flux_lines_on_plane", C.flux_lines_on_plane); j.set("floating_ground", C.floating_ground);
    j.set("junctions", C.junctions); j.set("junctions_landed", C.junctions_landed); j.set("sliver_fills_undone", C.fills_reverted);
    j.set("ground_straps", C.straps); j.set("ground_straps_landed", C.straps_landed); j.set("strap_line_crossings", C.strap_lines_crossed); j.set("strap_coupler_crossings", C.strap_couplers_crossed);
    j.set("fragments_removed", C.fragments_removed);
    j.set("junctions_on_device", C.junctions_on_device); j.set("devices", C.devices); j.set("devices_one_junction", C.devices_one_junction);
    Json rf = Json::object();
    rf.set("ok", C.rf_ground_ok()); rf.set("junction_grounds_on_plane_metal", C.junction_grounds_galvanic);
    rf.set("shorted_lines_on_plane_metal", C.shorted_lines_galvanic); rf.set("flux_shorts_on_plane_metal", C.flux_shorts_galvanic);
    rf.set("strapped_returns", C.strapped_returns); rf.set("returns_on_single_strap", C.returns_single_strap);
    rf.set("junction_grounds_shared_with_flux_return", C.shared_returns);
    double worst = 0; for (double v : C.junction_shared_pH) worst = std::max(worst, v);
    rf.set("worst_shared_strap_chain_pH", worst);
    j.set("rf_ground", rf);
    return j;
}

// One-line DRC summary and the JSON block the review programs write.
inline std::string drc_summary(const ChipLayout& L) {
    char b[400];
    std::snprintf(b, sizeof b, "placement/island-spacing checks %s (%zu problems); merged-etch width/space (%.0f/%.0f um, projection metric): "
                  "%d width + %d space violations", L.problems.empty() ? "clean" : "FAIL", L.problems.size(), L.rules.etch_min_width_um,
                  L.rules.etch_min_space_um, L.drc_count("width"), L.drc_count("space"));
    return b;
}
inline Json drc_json(const ChipLayout& L) {
    Json j = Json::object();
    j.set("placement_clean", L.problems.empty()); j.set("placement_problems", int(L.problems.size()));
    j.set("etch_min_width_um", L.rules.etch_min_width_um); j.set("etch_min_space_um", L.rules.etch_min_space_um);
    j.set("etch_width_violations", L.drc_count("width")); j.set("etch_space_violations", L.drc_count("space"));
    j.set("metric", "projection"); j.set("ground_sliver_fills", L.slivers_filled); j.set("clean", L.drc_clean());
    return j;
}

inline void write_gds(const ChipLayout& L, const std::string& path, const std::string& cell = "TOP") {
    gds::Structure s; s.name = cell;
    auto poly = [](const Rect& r, int layer) { return gds::Polygon{layer, 0, {{r.x0, r.y0}, {r.x1, r.y0}, {r.x1, r.y1}, {r.x0, r.y1}}}; };
    for (const auto& r : L.etch) s.polys.push_back(poly(r, 1));
    for (const auto& r : L.jj) s.polys.push_back(poly(r, 2));
    for (const auto& r : L.bumps) s.polys.push_back(poly(r, 3));
    for (const auto& st : L.straps) { s.polys.push_back(poly(st.span, 4)); s.polys.push_back(poly(st.foot_a, 5)); s.polys.push_back(poly(st.foot_b, 5)); }
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
    g.set("bus_width_um", c.cg.w_bus); g.set("bus_gap_um", c.cg.gap_bus); g.set("bus_end_gap_um", c.cg.bus_gap());
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
    c.cg.gap_bus = g.get("bus_gap_um", c.cg.gap_bus); c.cg.g_b = g.get("bus_end_gap_um", c.cg.g_b);
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
