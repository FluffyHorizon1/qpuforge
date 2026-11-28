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
#include <set>
#include <tuple>
#include <fstream>
#include <functional>

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

// PanelSolution (optional output): the panels, and for every conductor k the charge on every panel (fF x 1 V) with
// conductor k at 1 V and the others at 0. The charges of row k summed over the panels of conductor j give C(k, j)
// before symmetrisation.
struct PanelSolution {
    std::vector<Panel> panels;
    std::vector<std::vector<double>> q_fF;          // [excited conductor][panel]
    // charge on every panel for the conductor potentials V (volts): the superposition of the rows
    std::vector<double> charge(const std::vector<double>& V) const {
        std::vector<double> q(panels.size(), 0.0);
        for (std::size_t k = 0; k < q_fF.size() && k < V.size(); ++k) if (V[k] != 0) for (std::size_t i = 0; i < q.size(); ++i) q[i] += V[k] * q_fF[k][i];
        return q;
    }
};
inline CapacitanceMatrix extract_capacitance(const std::vector<Conductor>& conds, double eps_r = 11.45, const MeshOptions& mo = {},
                                             PanelSolution* sol = nullptr) {
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
        if (sol) { sol->q_fF.emplace_back(std::size_t(N)); for (int i = 0; i < N; ++i) sol->q_fF.back()[std::size_t(i)] = b[std::size_t(i)] / Kc * 1e15; }
    }
    if (sol) sol->panels = P;
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
    double relief_gap{10};  // reduced T-bar gap at the three free inner corners (< gap: opens their pockets, see xmon); <= 0: none
    double squid_w{20};     // SQUID loop width along its ground edge
    double squid_gap{10};   // SQUID loop depth: the T-bar end carries a finger that comes this close to the ground edge
    double slot_up{14}, slot_down{30};  // arms of the flux launcher's collar slot towards corner 2 and away from it (see kCollarW)
    double slot_feed{-15};              // where the feed slot joins the collar slot, along the edge from the SQUID's axis (< 0: away from corner 2)
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
    double slot_lo{38}, slot_hi{12};   // arms of the coupler flux launcher's collar slot, on the side where the line's gap stops short and on the other
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

// Qubit metal (cross + T-bars + readout pad + SQUID finger) and its cut-out, centred at (cx, cy).
//
// Inner corners. Corner k lies between arms k and k + 1 (0: (+x,+y), the readout corner; 1: (-x,+y); 2: (-x,-y);
// 3: (+x,-y)). The cut-outs of the two T-bars that meet at a corner leave the corner's ground (the pocket) joined to
// the ground outside by a neck between the two T-bar ends.
//   - Corners 1 and 3: the T-bar cut-out is drawn with the smaller gap `relief_gap` at the T-bar's end and at its back
//     face, on the half of the T-bar that faces the corner. The neck is sqrt(2) (L - t_w - T/2 - 2 relief_gap) wide.
//   - Corner 2 gets the relief only at the back face of the -x arm's T-bar. The -y arm's T-bar keeps the full gap
//     all round: the SQUID spans its end gap, and the straight ground edge beside that gap, x = -(T/2 + gap), then
//     runs the same distance (gap + t_w / 2) to either side of the SQUID. The pocket joins the ground outside
//     through a straight neck along the top of that edge, L - t_w - relief_gap - T/2 - gap wide (assemble_chip
//     reports a neck below the minimum width as a problem). The qubit's return current reaches the SQUID along
//     that edge from both sides.
//
// SQUID. It sits at the -x end of the -y arm's T-bar, in the (-x,-y) plaquette: a finger squid_w wide runs from the
// T-bar end to squid_gap from the straight ground edge x = -(T/2 + gap), and the junction region (layer 2) covers the
// loop between the finger and that edge, overlapping both metals by 1 um.
namespace detail { inline std::vector<Rect> subtract(const Rect& a, const Rect& b); }
inline constexpr double kSquidW = 20.0;      // coupler SQUID loop width along its ground edge (um); a qubit's is XmonParams::squid_w
inline double corner_neck(const XmonParams& q, int corner) {   // width of a corner's neck (<= 0: sealed)
    const bool relief_on = q.relief_gap > 0 && q.relief_gap < q.gap;
    if (corner == 2) return q.L - q.t_w - q.T / 2 - q.gap - (relief_on ? q.relief_gap : q.gap);        // straight, beside the -y T-bar's cut-out
    const double g = relief_on && (corner == 1 || corner == 3) ? q.relief_gap : q.gap;
    return std::sqrt(2.0) * (q.L - q.t_w - q.T / 2 - 2 * g);                                           // diagonal, between the two T-bar ends
}
inline void xmon(const XmonParams& q, double cx, double cy, std::vector<Rect>& metal, std::vector<Rect>& holes, std::vector<Rect>* jj = nullptr) {
    const double hw = q.w_arm / 2;
    const bool relief_on = q.relief_gap > 0 && q.relief_gap < q.gap;
    for (int k = 0; k < 4; ++k) {
        metal.push_back(rot90({0, -hw, q.L - q.t_w, hw}, k).shift(cx, cy));
        metal.push_back(rot90({q.L - q.t_w, -q.T / 2, q.L, q.T / 2}, k).shift(cx, cy));
        holes.push_back(rot90({-hw - q.gap, -hw - q.gap, q.L - q.t_w + q.gap, hw + q.gap}, k).shift(cx, cy));
        // T-bar cut-out: the band of the arm's own cut-out, then one piece per half. The +y half of arm k faces
        // corner k, its -y half corner k - 1.
        for (int half = 0; half < 2; ++half) {
            const int corner = half == 0 ? k : (k + 3) % 4;
            const bool c13 = corner == 1 || corner == 3, c2x = corner == 2 && k == 2;         // corner 2: the -x arm's T-bar only, back face only
            const double back = relief_on && (c13 || c2x) ? q.relief_gap : q.gap, end = relief_on && c13 ? q.relief_gap : q.gap;
            Rect h{q.L - q.t_w - back, hw + q.gap, q.L + q.gap, q.T / 2 + end};
            if (h.y1 <= h.y0) continue;
            if (half == 1) h = {h.x0, -h.y1, h.x1, -h.y0};
            holes.push_back(rot90(h, k).shift(cx, cy));
        }
    }
    metal.push_back(rect_c(cx, cy, q.w_arm, q.w_arm));
    metal.push_back(Rect{hw, hw, hw + q.r_pad, hw + q.r_pad}.shift(cx, cy));
    const double cl = hw + q.r_pad + q.claw_gap + q.claw_w + q.gap;
    holes.push_back(Rect{hw, hw, cl, cl}.shift(cx, cy));
    const double sy = -q.L + q.t_w / 2, xf = -q.T / 2 - (q.gap - q.squid_gap);       // finger tip
    metal.push_back(Rect{xf, sy - q.squid_w / 2, -q.T / 2, sy + q.squid_w / 2}.shift(cx, cy));
    if (jj) jj->push_back(Rect{-q.T / 2 - q.gap - 1, sy - q.squid_w / 2, xf + 1, sy + q.squid_w / 2}.shift(cx, cy));
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
// separate_xy selects the control-line mode of the mask: false (default), one line per qubit that carries both the
// flux bias and the microwave drive (the drive couples through the flux line's common-mode mutual, see flux_coupling);
// true, a separate XY line per qubit next to its flux line (qubit_control_lines; xy_tip sets its coupling).
struct CellSet { XmonParams meas, data; CouplerGeom cg; double pitch{900}; bool separate_xy{false}; double xy_tip{100.0}; };

namespace detail { inline int fill_ground_slivers(std::vector<Rect>& holes, const std::vector<Rect>& metal, double min_ground, int max_passes = 6); }
// Ground cut-outs of an extraction cell after the same mask clean-up as the assembled chip (ground slivers below the
// 4 um minimum merged away), so that the cells are extracted as they will be fabricated.
inline std::vector<Rect> cell_ground_holes(std::vector<Rect> holes, const std::vector<Conductor>& metal, double min_ground = 4.0) {
    std::vector<Rect> m;
    for (const auto& c : metal) if (!c.ground) m.insert(m.end(), c.rects.begin(), c.rects.end());
    detail::fill_ground_slivers(holes, m, min_ground);
    return holes;
}
// The near end of the coupler and bus on arm k of a qubit at the origin, as on the chip: on the +-x arms the coupler
// runs below the axis and the bus above, on the +-y arms the coupler to the right and the bus to the left (the -x and
// -y arms hold the far ends of the neighbours' pairs: the +x and +y drawings mirrored).
inline PairElements arm_elements(const XmonParams& q, const CouplerGeom& c, double p, int k) {
    PairElements E = coupler_and_bus(q, c, p, k % 2, 0, 0, false);
    if (k >= 2)
        for (auto* v : {&E.coupler, &E.bus, &E.holes, &E.jj})
            for (auto& r : *v) r = k == 2 ? Rect{-r.x1, r.y0, -r.x0, r.y1} : Rect{r.x0, -r.y1, r.x1, -r.y0};
    return E;
}
// Star cell: one qubit with its 4 couplers + buses (truncated beyond the mid pad) and the readout claw.
inline CapacitanceMatrix star_cell(const XmonParams& q, const CouplerGeom& c, double p, const MeshOptions& mo) {
    std::vector<Conductor> cs;
    Conductor Q; Q.name = "qubit";
    std::vector<Rect> holes;
    xmon(q, 0, 0, Q.rects, holes);
    cs.push_back(Q);
    for (int k = 0; k < 4; ++k) {
        auto E = arm_elements(q, c, p, k);
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
            // ... and the neck of the SQUID corner (full gap) no narrower than 4 um across (5.7 um diagonal)
            lo[std::size_t(q)] = std::max({lo[std::size_t(q)], Tb / 2 + 20, claw_outer + Q.t_w + 12, 80.0, Q.t_w + Q.T / 2 + Q.gap + std::min(Q.gap, std::max(Q.relief_gap, 0.0) > 0 ? Q.relief_gap : Q.gap) + 4});
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
inline constexpr double kStrapFootEnclosure = 2.0;      // ground metal required around every foot (um)
inline constexpr double kStrapQubitClearance = 10.0;    // a strap keeps this far from qubit metal (um)
inline constexpr double kStrapSeparation = 25.0;        // two straps between the same two fragments are at least this far apart (um)
struct Strap { Rect span, foot_a, foot_b; int lines_crossed{0}, couplers_crossed{0}; };
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
        I.kind.push_back(n[0] == 'Q' ? 0 : (n[0] == 'R' || n[0] == 'A' || n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n.rfind("XY", 0) == 0) ? 1 : 2);
        for (const auto& r : islands[std::size_t(i)].second) { I.rs.push_back({r, i}); I.wmax = std::max(I.wmax, r.w()); }
    }
    std::sort(I.rs.begin(), I.rs.end(), [](const auto& a, const auto& b) { return a.first.x0 < b.first.x0; });
    return I;
}
// component of the ground metal under a foot square, or -1 if the square with kStrapFootEnclosure around it is not
// wholly ground of one component
inline int foot_comp(const MetalMap& M, const IslandIndex& I, const Rect& foot) {
    const Rect f = foot.grow(kStrapFootEnclosure);
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
// eligible[c]: component c may carry a strap foot (plane, floating ground, shorted-line pockets).
// joined (optional): components already tied together (plane metal, straps placed before); the placer adds the
// straps that tie the remaining eligible components to them (a spanning forest, cheapest first), and then a second
// strap, at least kStrapSeparation away, between each pair of components it tied, so that no fragment hangs on one strap.
inline std::vector<Strap> place_ground_straps(const MetalMap& M, const IslandIndex& I, const std::vector<char>& eligible, int plane,
                                              UnionFindLite* joined = nullptr) {
    struct Cand { double cost; int a, b; Strap s; };
    std::map<std::pair<int, int>, std::vector<Cand>> cands;
    std::vector<std::vector<std::size_t>> ivs(std::size_t(M.n_comp));
    for (std::size_t i = 0; i < M.comp.size(); ++i) ivs[std::size_t(M.comp[i])].push_back(i);
    const double F = kStrapW;
    const std::array<std::pair<int, int>, 4> dirs{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
    auto near_qubit = [&](const Rect& sp) { bool bad = false; I.each(sp.grow(kStrapQubitClearance - 1e-6), [&](int i) { if (!I.kind[std::size_t(i)]) bad = true; }); return bad; };
    for (int T = 0; T < M.n_comp; ++T) {
        if (T == plane || !eligible[std::size_t(T)]) continue;
        // sample foot positions on the component's intervals (every 5 um along each, at most ~200 per component)
        std::vector<std::pair<double, double>> pts;
        for (std::size_t i : ivs[std::size_t(T)]) {
            if (M.ix1[i] - M.ix0[i] < F || M.iy1[i] - M.iy0[i] < 1.0) continue;
            const double ym = 0.5 * (M.iy0[i] + M.iy1[i]);
            for (double x = M.ix0[i] + F / 2; x <= M.ix1[i] - F / 2 + 1e-9; x += 5) {
                pts.push_back({x, ym});
                // a tall interval: rows of positions every 10 um as well (the middle may be taken by a line's metal)
                for (double y = M.iy0[i] + F / 2 + kStrapFootEnclosure; y <= M.iy1[i] - F / 2 - kStrapFootEnclosure + 1e-9; y += 10) if (std::abs(y - ym) > 2.5) pts.push_back({x, y});
            }
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
                            // spanning: every 10 um, give up this direction once the span would come too near a qubit
                            if (step % 5 == 0 && near_qubit({std::min(fa.x0, fb.x0), std::min(fa.y0, fb.y0), std::max(fa.x1, fb.x1), std::max(fa.y1, fb.y1)})) break;
                            continue;
                        }
                        const int U = c0;
                        Strap s; s.foot_a = fa; s.foot_b = fb;
                        s.span = {std::min(fa.x0, fb.x0), std::min(fa.y0, fb.y0), std::max(fa.x1, fb.x1), std::max(fa.y1, fb.y1)};
                        if (near_qubit(s.span)) break;
                        std::set<int> crossed; I.each(s.span.grow(-1e-6), [&](int i) { crossed.insert(i); });
                        for (int i : crossed) (I.kind[std::size_t(i)] == 1 ? s.lines_crossed : s.couplers_crossed)++;
                        const double cost = (t - F) + kStrapLinePenalty * s.lines_crossed + kStrapCouplerPenalty * s.couplers_crossed;
                        cands[{std::min(T, U), std::max(T, U)}].push_back({cost, T, U, s});
                        break;
                    }
                }
            }
        }
    }
    // per pair: the cheapest candidate, and the cheapest one at least kStrapSeparation from it
    struct Pair { Cand first; bool has_second{false}; Cand second; };
    std::vector<Pair> ps;
    for (auto& [key, v] : cands) {
        std::sort(v.begin(), v.end(), [](const Cand& a, const Cand& b) { return a.cost < b.cost || (a.cost == b.cost && std::tie(a.s.span.x0, a.s.span.y0) < std::tie(b.s.span.x0, b.s.span.y0)); });
        Pair p; p.first = v.front();
        for (const auto& c : v)
            if (std::hypot(c.s.span.cx() - p.first.s.span.cx(), c.s.span.cy() - p.first.s.span.cy()) >= kStrapSeparation + kStrapW) { p.second = c; p.has_second = true; break; }
        ps.push_back(p);
    }
    std::sort(ps.begin(), ps.end(), [](const Pair& x, const Pair& y) { const Cand &a = x.first, &b = y.first; return a.cost < b.cost || (a.cost == b.cost && std::tie(a.a, a.b) < std::tie(b.a, b.b)); });
    UnionFindLite local(M.n_comp);
    UnionFindLite& uf = joined ? *joined : local;
    std::vector<Strap> out;
    for (const auto& p : ps) if (uf.unite(p.first.a, p.first.b)) { out.push_back(p.first.s); if (p.has_second) out.push_back(p.second.s); }
    return out;
}
} // namespace detail

// Flux launcher: how a flux line ends at a SQUID that sits on a straight ground edge (frame: the edge is v = 0, ground
// at v > 0, the SQUID centred at u = 0).
//   - collar: the strip of ground 0 < v < kCollarW along the edge; the SQUID's ground side lands on it;
//   - collar slot: etch kSlotW wide behind the collar, running both ways from the SQUID's axis;
//   - the line's centre conductor ends kTipGround behind the slot and is shorted to the ground on ONE side (that gap
//     stops kShortLen before the tip); the gap on the other side runs on into the collar slot.
// The line's return current on the slot side has to pass between the slot and the SQUID's cut-out, i.e. along the
// collar, past the loop, in one direction and symmetrically about it: half the line current couples to the loop. The
// qubit's own current leaves the collar at both ends into wide ground; both of the line's gaps lie behind the slot,
// on the same side of those two streams, so the streams do not drive the line.
inline constexpr double kCollarW = 4.0, kSlotW = 5.0, kTipGround = 7.0, kShortLen = 15.0;
// The qubit's launcher has an unequal slot (XmonParams::slot_up towards the (-x,-y) corner, slot_down the other way).
// The line's two gaps are at opposite potentials (they carry its current), and each couples a little to the edge of
// the qubit's cut-out on its own side; with equal arms that leaves a common-mode mutual M_cm of about 1.5 pH. A longer
// arm on the side away from the corner puts more of the collar's flux there and cancels most of it. The balance is
// set per cell file (it depends on the arm length): M_cm changes by about 0.06 pH per um of the long arm and
// 0.025 pH per um of the short one.

// Flux line (FQ) of the qubit at (cx, cy), from a bump pad in the (-x,-y) plaquette to the flux launcher at the SQUID.
// The ground edge is x = cx - (T/2 + gap) and the SQUID sits at the -x end of the -y arm's T-bar. The line leaves its
// pad as a diagonal staircase CPW of kFluxStair squares; kFluxSlotSquares squares before the launcher the gap on the
// lower right (the side of the qubit's cut-out) stops and the centre conductor joins the ground there, so that over
// the last stretch only ONE slot, on the side away from the qubit, runs on to the collar slot (the "feed slot"). The
// feed slot meets the collar slot XmonParams::slot_feed from the SQUID's axis. The pad has a ground bump on its two
// far sides (flux_ground_bumps), which is where the line's return current leaves the chip.
// The line is kept away from the edges of the qubit's and couplers' cut-outs: the qubit's return current flows along
// those edges to the collar, and a line lying in that stream picks it up on both of its grounds
// (FluxCoupling::M_cm_pH). With short_flux false the line is drawn without the launcher, ending open (tests).
// With separate_xy (CellSet::separate_xy) a separate XY line is drawn as well: a staircase CPW that ends open in
// front of the neck of the (-x,+y) pocket, at (cx - xy_tip, cy + xy_tip), fed through a kXYGap break from its pad.
struct ControlLineGeom {
    std::pair<double, double> flux_bump, xy_bump;
    std::vector<std::pair<double, double>> flux_ground_bumps, xy_ground_bumps;      // ground bumps beside each pad
};
inline constexpr double kXYRun = 125.0, kXYAntennaRun = 20.0, kXYGap = 4.0, kXYFeedRun = 20.0, kFluxSlotRun = 10.0, kGroundBumpPitch = 50.0;
inline constexpr int kFluxStair = 8, kFluxSlotSquares = 10;
inline ControlLineGeom qubit_control_lines(const XmonParams& q, double cx, double cy, std::vector<Rect>& fm, std::vector<Rect>& fh, std::vector<Rect>& xm,
                                           std::vector<Rect>& xh, std::vector<Rect>& am, bool short_flux, bool separate_xy = false, double xy_tip = 100.0) {
    using detail::cpw_path; using detail::stair;
    ControlLineGeom G;
    const double edge = -q.T / 2 - q.gap, yc = -q.L + q.t_w / 2;                    // ground edge and SQUID centre (local)
    const double xs1 = edge - kCollarW, xs0 = xs1 - kSlotW;
    auto add = [&](std::vector<Rect>& v, const Rect& r) { v.push_back(r.shift(cx, cy)); };
    // A staircase of n + 1 squares from the pad end (k = 0) up to the launcher (k = n). Squares 0 .. n_cpw are a CPW
    // (centre conductor, both gaps); square n_cpw is the conductor's tip, shorted to the ground on its upper-right
    // faces, where the lower-right gap ends. From there on only the upper-left gap continues, as a single slot (the
    // feed slot): its last stretch runs straight along the SQUID's axis into the collar slot.
    const int n_slot = kFluxSlotSquares, n = kFluxStair + n_slot, n_cpw = n - n_slot;
    const double s0 = kFluxSlotRun;
    const double feed = q.slot_feed;
    const double xn = xs0 - s0 - 5, yn = yc + feed - 8;                                // centre of square n: its top strip lies on the axis
    const double xe = xn - 5 * n, ye = yn - 5 * n;
    for (int k = 0; k <= n; ++k) {
        const double px = xe + 5 * k, py = ye + 5 * k;
        if (!short_flux) { add(fm, rect_c(px, py, 10, 10)); add(fh, rect_c(px, py, 22, 22)); continue; }
        if (k <= n_cpw) add(fm, rect_c(px, py, 10, 10));
        if (k < n_cpw - 1) add(fh, rect_c(px, py, 22, 22));
        else if (k == n_cpw - 1) add(fh, {px - 11, py - 11, px + 10, py + 10});         // trimmed to the tip's outer faces
        else if (k == n_cpw) {
            add(fh, {px - 11, py - 11, px + 11, py - 5}); add(fh, {px - 11, py - 5, px - 5, py + 11});       // bottom and left strips
            add(fh, {px - 11, py + 5, px, py + 11}); add(fh, {px + 5, py - 5, px + 11, py});                  // top and right strips, stopping at the tip's corner
            add(fh, {px - 5, py - 5, px + 5, py + 5});                                                        // under the tip (no etch)
        } else { add(fh, {px - 11, py - 5, px - 5, py + 11}); add(fh, {px - 11, py + 5, px + 5, py + 11}); } // feed slot: the upper-left gap alone
    }
    if (short_flux) {
        add(fh, {xn + 5, yc + feed - 3, xs0, yc + feed + 3});                           // the feed slot's straight end
        add(fh, {xs0, yc - q.slot_down, xs1, yc + q.slot_up});                          // collar slot
    }
    G.flux_bump = {cx + xe - 18, cy + ye - 18};
    fm.push_back(rect_c(G.flux_bump.first, G.flux_bump.second, 36, 36)); fh.push_back(rect_c(G.flux_bump.first, G.flux_bump.second, 52, 52));
    G.flux_ground_bumps = {{G.flux_bump.first - kGroundBumpPitch, G.flux_bump.second}, {G.flux_bump.first, G.flux_bump.second - kGroundBumpPitch}};
    if (!separate_xy) return G;
    // XY: antenna (floating island: open tip at (-xy_tip, +xy_tip), staircase along the anti-diagonal to the level
    // kXYRun, kXYAntennaRun of straight CPW), a kXYGap gap in the centre conductor, then the feed to the pad
    const double tipv = xy_tip, yrun = kXYRun, alen = kXYAntennaRun, flen = kXYFeedRun, gapc = kXYGap;
    const int m = int(std::lround((yrun - tipv) / 5.0));                           // squares up (m > 0) or down to the run
    const double ex = cx - tipv - 5 * std::abs(m), ey = cy + tipv + 5 * m;
    std::vector<Rect> none;
    stair(ex, ey, cx - tipv, cy + tipv, am, xh, false);
    cpw_path({{ex, ey}, {ex - alen - gapc - flen, ey}}, 10, 6, none, xh, false);      // the trench runs through
    am.push_back({ex - alen, ey - 5, ex + 5, ey + 5});
    xm.push_back({ex - alen - gapc - flen - 5, ey - 5, ex - alen - gapc, ey + 5});
    G.xy_bump = {ex - alen - gapc - flen - 12, ey};
    xm.push_back(rect_c(G.xy_bump.first, G.xy_bump.second, 36, 36)); xh.push_back(rect_c(G.xy_bump.first, G.xy_bump.second, 52, 52));
    G.xy_ground_bumps = {{G.xy_bump.first - kGroundBumpPitch, G.xy_bump.second}};      // the resonator and the bus leave no room above or below
    return G;
}

// Width / space rules on the fabricated etch layer (layer 1). Metal is the outline minus etch, so the etch-space rule is
// also the minimum metal (island and ground-plane) width, and the etch-width rule the minimum metal-to-metal gap.
// Default: the 4 um process minimum that assemble_chip already uses for designed coupling gaps.
struct DrcRules { double etch_min_width_um{4.0}, etch_min_space_um{4.0}; };
inline constexpr double kRemoveFragmentUm2 = 2500;   // unreachable ground fragments below this area are etched away

// Connectivity of the fabricated metal (LVS-lite), from detail::metal_map on the final etch layer.
// Nets (metal components joined by the ground straps):
//   - every island (qubit, coupler, bus, line) is its own conductor: no two islands share a metal component, and only
//     the lines meant to end in a short (resonators R*, qubit flux lines FQ*, coupler flux lines FC*) touch the ground;
//   - each shorted line reaches the ground net, so it has a DC return; no ground metal is cut off from it;
//   - each junction region (layer 2) has an owner island (its qubit or coupler). It lands if the metal around its
//     boundary is the owner's and ground, and nothing else; every qubit and coupler island owns exactly one junction.
// Ground quality (a strap is a modelled airbridge with inductance, so a net is not yet a ground):
//   - plane metal is the metal component of the chip border, and every component that carries a ground bump
//     (ChipLayout::ground_bumps: flip-chip bumps to the carrier's ground; none unless the layout draws them);
//   - for the ground side of every SQUID and for every shorted line end: is it plane metal itself, and if not, the
//     fewest straps in series to plane metal and their summed bridge length (ChipLayout::grounds);
//   - per strap: how many of those circuits return through it (along their fewest-strap paths);
//   - strapped fragments that hang on a single strap (a strap whose removal cuts them off: a bridge of the strap graph).
// ground_ok() is the gate: the ground side of every qubit SQUID and every qubit flux-line short is plane metal with no
// strap in the path; the ground on both sides of each qubit flux line's bump pad is that same metal component (the
// line's return is continuous metal); every other strapped fragment keeps its connection when any one strap is
// removed, and lies at most kMaxStrapDepth straps from plane metal.
inline constexpr int kMaxStrapDepth = 2;
struct GroundRef {
    std::string name;         // SQ<id> (qubit SQUID), SC<a-b> (coupler SQUID), FQ<id>, FC<a-b>, R<id> (shorted line ends)
    int comp{-1};             // metal component of its ground side (-1: no ground there)
    bool on_plane{false};     // that component is plane metal (no strap in the path)
    int straps{-1};           // fewest straps in series to plane metal (0 on it; -1: not reachable)
    double span_um{0};        // summed bridge length (foot centre to foot centre) along that path
    bool on_border{false};    // plane metal joined in 2D to the chip border (no bump needed for its reference)
    int bumps{0};             // ground bumps on its component
    double bump_dist_um{-1};  // SQUIDs: distance to the nearest ground bump on its component (-1: none)
};
// Ground bumps (flip-chip bumps from ground metal to the top chip's ground plane, GDS layer 6). Rules checked by
// analyze_connectivity:
//   - a plaquette of plane metal that is not joined in 2D to the chip border takes its reference through its bumps
//     alone: at least kMinPlaquetteBumps of them;
//   - every SQUID ground on such a plaquette has a ground bump within kMaxSquidBumpDist_um, and every signal bump has
//     a ground bump within kSignalGroundBump_um (its return);
//   - bump geometry: bumps are kBumpSize squares, at least kBumpPitch apart centre to centre; a ground bump has
//     kBumpEnclosure of ground metal around it and keeps kBumpClearance from every island (qubit, coupler, bus,
//     resonator, line).
// kBumpInductance_pH is the inductance assumed per bump in the line budgets: the partial self-inductance of a post
// 30 um across and 10 um tall is about 1.1 pH; 2 pH allows for the current spreading into the bump. An assumption.
inline constexpr int kMinPlaquetteBumps = 4;
inline constexpr double kMaxSquidBumpDist_um = 400.0, kSignalGroundBump_um = 75.0;
inline constexpr double kBumpSize = 30.0, kBumpPitch = 50.0, kBumpEnclosure = 3.0, kBumpClearance = 10.0;
inline constexpr double kBumpInductance_pH = 2.0;
inline constexpr double kBumpGrid = 50.0, kBumpGridKeepout = 40.0, kBumpGridLineKeepout = 20.0;   // the ground-bump lattice of assemble_chip: pitch, distance kept from qubits, couplers and buses, and from lines
struct ConnectivityReport {
    int components{0}, islands{0}, islands_merged{0}, islands_on_plane{0};
    int shorted_lines{0}, shorted_lines_on_plane{0}, flux_lines{0}, flux_lines_on_plane{0};
    int floating_ground{0}, junctions{0}, junctions_landed{0};
    int junction_owner_errors{0};                        // qubit/coupler islands without exactly one junction, junctions without such an owner
    int straps{0}, straps_landed{0}, strap_lines_crossed{0}, strap_couplers_crossed{0};   // ground straps (modelled airbridges), see place_ground_straps
    int fragments_removed{0};                            // small unreachable ground fragments etched away
    // ground quality
    int ground_bumps{0}, plane_components{0};
    int qubit_squids{0}, qubit_squids_on_plane{0}, qubit_flux_on_plane{0}, flux_returns_continuous{0};
    int coupler_squids{0}, coupler_squids_on_plane{0}, other_shorts{0}, other_shorts_on_plane{0};
    int strapped_fragments{0}, fragments_on_one_strap{0}, cut_straps{0};
    int max_strap_depth{0};                              // over all strapped fragments
    int max_circuit_depth{0}; double max_circuit_span_um{0};   // over the SQUID grounds and shorted line ends
    int max_returns_per_strap{0};
    // ground bumps: which plane metal is joined to the border in 2D, and which is grounded by bumps alone
    int border_bumps{0}, bump_plaquettes{0}, min_plaquette_bumps{0};
    int qubit_squids_on_border{0}, qubit_squids_on_bumps{0}, coupler_squids_on_border{0}, coupler_squids_on_bumps{0};
    double max_squid_bump_dist_um{0};                    // over the SQUID grounds on bump-grounded plaquettes
    int squids_far_from_bump{0}, signal_bumps{0}, signal_bumps_without_ground{0}, bump_rule_violations{0};
    bool bumps_ok() const {
        return (bump_plaquettes == 0 || min_plaquette_bumps >= kMinPlaquetteBumps) && squids_far_from_bump == 0 &&
               (ground_bumps == 0 || signal_bumps_without_ground == 0) && bump_rule_violations == 0;
    }
    bool nets_clean() const {
        return islands_merged == 0 && islands_on_plane == 0 && shorted_lines_on_plane == shorted_lines && floating_ground == 0 &&
               junctions_landed == junctions && junction_owner_errors == 0 && straps_landed == straps;
    }
    bool ground_ok() const {
        return qubit_squids_on_plane == qubit_squids && coupler_squids_on_plane == coupler_squids && qubit_flux_on_plane == flux_lines &&
               flux_returns_continuous == flux_lines && fragments_on_one_strap == 0 && max_strap_depth <= kMaxStrapDepth && bumps_ok();
    }
    bool clean() const { return nets_clean() && ground_ok(); }
};
using detail::DrcMarker;

struct ChipLayout {
    std::string name;
    Rect outline;
    std::vector<std::pair<std::string, std::vector<Rect>>> islands;   // named metal islands (qubits, couplers, buses, lines)
    std::vector<Rect> holes, jj, bumps, etch;
    std::vector<int> jj_owner;                         // per junction region: index of its qubit or coupler island (-1: none)
    std::vector<int> bump_island;                      // per signal bump: index of the line island it feeds (-1: none)
    std::vector<Rect> ground_bumps;                    // flip-chip bumps on ground metal, to the carrier's ground (GDS layer 6)
    std::vector<int> ground_bump_line;                 // per ground bump: the line island whose pad it flanks (-1: none)
    int lattice_bumps{0};                              // ground bumps placed on the kBumpGrid lattice (the last ones in ground_bumps)
    std::vector<detail::Strap> straps;                 // ground straps (GDS layer 4: span, layer 5: feet)
    std::vector<std::pair<std::string, std::pair<double, double>>> labels;
    struct Resonator { int qubit; double f_GHz, target_um, drawn_um; int rows; };
    std::vector<Resonator> resonators;
    std::vector<std::string> problems;                 // placement / spacing / outline problems found during assembly
    std::vector<DrcMarker> drc;                         // width / space violations on the merged etch layer
    DrcRules rules;
    int slivers_filled{0};                              // ground-sliver fills added by the mask clean-up (fill_ground_slivers)
    ConnectivityReport conn;                            // metal connectivity of the final mask (see ConnectivityReport)
    std::vector<GroundRef> grounds;                     // ground reference of every SQUID and shorted line end
    std::vector<int> strap_returns;                     // per strap: circuits returning through it
    double etch_area_mm2{0};
    int island(const std::string& n) const { for (std::size_t i = 0; i < islands.size(); ++i) if (islands[i].first == n) return int(i); return -1; }
    int drc_count(const std::string& rule_suffix) const {
        int n = 0; for (const auto& m : drc) if (m.rule.size() >= rule_suffix.size() && m.rule.compare(m.rule.size() - rule_suffix.size(), rule_suffix.size(), rule_suffix) == 0) ++n;
        return n;
    }
    bool drc_clean() const { return problems.empty() && drc.empty(); }
};

// Component of every island's metal (the component under the centre of its largest rectangle).
inline std::vector<int> island_components(const ChipLayout& L, const detail::MetalMap& M) {
    std::vector<int> v;
    for (const auto& [name, rs] : L.islands) {
        const Rect* big = &rs.front();
        for (const auto& r : rs) if (r.area() > big->area()) big = &r;
        v.push_back(M.at(big->cx(), big->cy()));
    }
    return v;
}
inline bool is_shorted_line(const std::string& n) { return n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n[0] == 'R'; }

// Fills L.conn, L.grounds and L.strap_returns from L.outline, L.etch, L.islands, L.jj (+ jj_owner), L.bumps
// (+ bump_island), L.ground_bumps and L.straps. conn.fragments_removed is left as it is.
inline void analyze_connectivity(ChipLayout& L) {
    const detail::MetalMap M = detail::metal_map(L.outline, L.etch);
    ConnectivityReport C; C.fragments_removed = L.conn.fragments_removed;
    L.grounds.clear(); L.strap_returns.assign(L.straps.size(), 0);
    C.components = M.n_comp;
    // nets: metal components joined by the ground straps
    detail::UnionFindLite uf(M.n_comp);
    struct Edge { int a, b; double len; int strap; };
    std::vector<Edge> edges;
    for (std::size_t k = 0; k < L.straps.size(); ++k) {
        const auto& st = L.straps[k];
        const int a = M.at(st.foot_a.cx(), st.foot_a.cy()), b = M.at(st.foot_b.cx(), st.foot_b.cy());
        C.straps++; C.strap_lines_crossed += st.lines_crossed; C.strap_couplers_crossed += st.couplers_crossed;
        if (a >= 0 && b >= 0) {                               // a strap between two points of one component is a redundant stitch
            C.straps_landed++; uf.unite(a, b);
            if (a != b) edges.push_back({a, b, std::hypot(st.foot_a.cx() - st.foot_b.cx(), st.foot_a.cy() - st.foot_b.cy()), int(k)});
        }
    }
    const int plane0 = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    auto net = [&](int c) { return c < 0 ? -1 : uf.find(c); };
    const int plane = net(plane0);
    const std::vector<int> icomp = island_components(L, M);
    std::map<int, int> island_of;                  // net -> number of islands in it
    for (int c : icomp) island_of[net(c)]++;
    C.islands = int(L.islands.size());
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const std::string& n = L.islands[i].first;
        const bool fq = n.rfind("FQ", 0) == 0, line = is_shorted_line(n);
        const int c = net(icomp[i]);
        if (line) { C.shorted_lines++; if (c == plane) C.shorted_lines_on_plane++; }
        else if (c == plane) C.islands_on_plane++;
        if (fq) { C.flux_lines++; if (c == plane) C.flux_lines_on_plane++; }
        // a conductor (qubit, coupler, bus, XY line) sharing its metal with any other island is a short; two shorted
        // lines ending in the same enclosed pocket count under shorted lines that miss the plane instead
        if (!line && island_of[c] > 1) C.islands_merged++;
    }
    std::set<int> nets; for (int c = 0; c < M.n_comp; ++c) nets.insert(net(c));
    for (int n : nets) if (n != plane && !island_of.count(n)) C.floating_ground++;
    // plane metal: the border component and every component under a ground bump
    std::vector<char> is_plane(std::size_t(M.n_comp), 0);
    if (plane0 >= 0) is_plane[std::size_t(plane0)] = 1;
    C.ground_bumps = int(L.ground_bumps.size());
    for (const auto& gb : L.ground_bumps) { const int c = M.at(gb.cx(), gb.cy()); if (c >= 0) is_plane[std::size_t(c)] = 1; }
    for (char p : is_plane) C.plane_components += p;
    // bumps per plane component; plaquettes grounded by bumps alone (not the border component)
    std::vector<int> comp_bumps(std::size_t(M.n_comp), 0);
    for (const auto& gb : L.ground_bumps) { const int c = M.at(gb.cx(), gb.cy()); if (c >= 0) comp_bumps[std::size_t(c)]++; }
    for (int c = 0; c < M.n_comp; ++c) {
        if (!is_plane[std::size_t(c)]) continue;
        if (c == plane0) { C.border_bumps = comp_bumps[std::size_t(c)]; continue; }
        C.min_plaquette_bumps = C.bump_plaquettes == 0 ? comp_bumps[std::size_t(c)] : std::min(C.min_plaquette_bumps, comp_bumps[std::size_t(c)]);
        C.bump_plaquettes++;
    }
    // signal bumps and their ground bump; bump geometry
    C.signal_bumps = int(L.bumps.size());
    {
        auto dist = [](const Rect& a, const Rect& b) { return std::hypot(a.cx() - b.cx(), a.cy() - b.cy()); };
        for (const auto& sb : L.bumps) {
            bool has = false; for (const auto& gb : L.ground_bumps) has = has || dist(sb, gb) <= kSignalGroundBump_um;
            if (!has) C.signal_bumps_without_ground++;
        }
        std::vector<Rect> all = L.bumps; all.insert(all.end(), L.ground_bumps.begin(), L.ground_bumps.end());
        for (std::size_t a = 0; a < all.size(); ++a) {
            if (std::abs(all[a].w() - kBumpSize) > 1e-6 || std::abs(all[a].h() - kBumpSize) > 1e-6) C.bump_rule_violations++;
            for (std::size_t b = a + 1; b < all.size(); ++b)
                if (std::abs(all[a].cx() - all[b].cx()) < kBumpPitch - 1e-6 && std::abs(all[a].cy() - all[b].cy()) < kBumpPitch - 1e-6) C.bump_rule_violations++;
        }
        for (const auto& gb : L.ground_bumps) {
            const Rect e = gb.grow(kBumpEnclosure);
            bool bad = false;
            for (const double x : {e.x0 + 0.25, e.cx(), e.x1 - 0.25}) for (const double y : {e.y0 + 0.25, e.cy(), e.y1 - 0.25}) bad = bad || M.at(x, y) < 0;
            const Rect k = gb.grow(kBumpClearance - 1e-6);
            for (const auto& isl : L.islands) for (const auto& r : isl.second) bad = bad || k.intersects(r);
            for (const auto& st : L.straps) bad = bad || k.intersects(st.span);
            if (bad) C.bump_rule_violations++;
        }
    }
    // fewest straps from plane metal to every component (ties: shorter summed bridge length)
    std::vector<std::vector<int>> adj(std::size_t(M.n_comp));
    for (std::size_t e = 0; e < edges.size(); ++e) { adj[std::size_t(edges[e].a)].push_back(int(e)); adj[std::size_t(edges[e].b)].push_back(int(e)); }
    std::vector<int> depth(std::size_t(M.n_comp), -1), via(std::size_t(M.n_comp), -1);
    std::vector<double> span(std::size_t(M.n_comp), 0.0);
    {
        std::set<std::tuple<int, double, int>> pq;
        for (int c = 0; c < M.n_comp; ++c) if (is_plane[std::size_t(c)]) { depth[std::size_t(c)] = 0; pq.insert({0, 0.0, c}); }
        while (!pq.empty()) {
            const auto [d, sp, c] = *pq.begin(); pq.erase(pq.begin());
            if (d != depth[std::size_t(c)] || sp != span[std::size_t(c)]) continue;
            for (int e : adj[std::size_t(c)]) {
                const int o = edges[std::size_t(e)].a == c ? edges[std::size_t(e)].b : edges[std::size_t(e)].a;
                const int nd = d + 1; const double ns = sp + edges[std::size_t(e)].len;
                if (depth[std::size_t(o)] < 0 || nd < depth[std::size_t(o)] || (nd == depth[std::size_t(o)] && ns < span[std::size_t(o)] - 1e-9)) {
                    depth[std::size_t(o)] = nd; span[std::size_t(o)] = ns; via[std::size_t(o)] = e; pq.insert({nd, ns, o});
                }
            }
        }
    }
    // straps whose removal cuts a fragment off plane metal: bridges of the strap graph with plane metal as one node
    {
        const int root = M.n_comp;                                 // all plane components contracted
        auto node = [&](int c) { return is_plane[std::size_t(c)] ? root : c; };
        std::vector<std::vector<std::pair<int, int>>> g(std::size_t(M.n_comp) + 1);
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const int a = node(edges[e].a), b = node(edges[e].b);
            if (a == b) continue;
            g[std::size_t(a)].push_back({b, int(e)}); g[std::size_t(b)].push_back({a, int(e)});
        }
        std::vector<int> tin(g.size(), -1), low(g.size(), 0);
        std::vector<char> is_bridge(edges.size(), 0);
        int timer = 0;
        struct Frame { int v, pe; std::size_t i; };
        std::vector<Frame> st{{root, -1, 0}};
        tin[std::size_t(root)] = low[std::size_t(root)] = timer++;
        while (!st.empty()) {
            Frame& f = st.back();
            if (f.i < g[std::size_t(f.v)].size()) {
                const auto [to, e] = g[std::size_t(f.v)][f.i++];
                if (e == f.pe) continue;
                if (tin[std::size_t(to)] >= 0) low[std::size_t(f.v)] = std::min(low[std::size_t(f.v)], tin[std::size_t(to)]);
                else { tin[std::size_t(to)] = low[std::size_t(to)] = timer++; st.push_back({to, e, 0}); }
            } else {
                const Frame done = f; st.pop_back();
                if (!st.empty()) {
                    const int pv = st.back().v;
                    low[std::size_t(pv)] = std::min(low[std::size_t(pv)], low[std::size_t(done.v)]);
                    if (low[std::size_t(done.v)] > tin[std::size_t(pv)]) is_bridge[std::size_t(done.pe)] = 1;
                }
            }
        }
        for (char bbr : is_bridge) C.cut_straps += bbr;
        // fragments that depend on a single strap: not reachable from plane metal once the bridges are removed
        std::vector<char> safe(g.size(), 0); std::vector<int> q{root}; safe[std::size_t(root)] = 1;
        while (!q.empty()) { const int v = q.back(); q.pop_back(); for (const auto& [to, e] : g[std::size_t(v)]) if (!is_bridge[std::size_t(e)] && !safe[std::size_t(to)]) { safe[std::size_t(to)] = 1; q.push_back(to); } }
        for (int c = 0; c < M.n_comp; ++c) {
            if (is_plane[std::size_t(c)] || depth[std::size_t(c)] < 0) continue;
            C.strapped_fragments++; C.max_strap_depth = std::max(C.max_strap_depth, depth[std::size_t(c)]);
            if (!safe[std::size_t(c)]) C.fragments_on_one_strap++;
        }
    }
    auto ground_ref = [&](const std::string& name, int comp, const Rect* at = nullptr) {
        GroundRef g; g.name = name; g.comp = comp;
        if (comp >= 0) {
            g.on_plane = is_plane[std::size_t(comp)] != 0; g.straps = depth[std::size_t(comp)]; g.span_um = span[std::size_t(comp)];
            g.on_border = comp == plane0; g.bumps = comp_bumps[std::size_t(comp)];
            if (at) for (const auto& gb : L.ground_bumps) {
                if (M.at(gb.cx(), gb.cy()) != comp) continue;
                const double d = std::hypot(gb.cx() - at->cx(), gb.cy() - at->cy());
                if (g.bump_dist_um < 0 || d < g.bump_dist_um) g.bump_dist_um = d;
            }
        }
        if (g.straps > 0) {
            C.max_circuit_depth = std::max(C.max_circuit_depth, g.straps); C.max_circuit_span_um = std::max(C.max_circuit_span_um, g.span_um);
            for (int c = comp; via[std::size_t(c)] >= 0 && !is_plane[std::size_t(c)];) {
                const Edge& e = edges[std::size_t(via[std::size_t(c)])];
                L.strap_returns[std::size_t(e.strap)]++;
                c = e.a == c ? e.b : e.a;
            }
        }
        L.grounds.push_back(g);
        return g;
    };
    // junctions: owner island at one end, the ground net at the other; exactly one per qubit and coupler island
    std::vector<int> owned(L.islands.size(), 0);
    for (std::size_t k = 0; k < L.jj.size(); ++k) {
        const Rect& j = L.jj[k];
        C.junctions++;
        const int o = k < L.jj_owner.size() ? L.jj_owner[k] : -1;
        const bool owner_ok = o >= 0 && o < int(L.islands.size()) && (L.islands[std::size_t(o)].first[0] == 'Q' || L.islands[std::size_t(o)].first[0] == 'C');
        if (!owner_ok) { C.junction_owner_errors++; continue; }
        owned[std::size_t(o)]++;
        // metal just outside the region's boundary: the owner's and ground, nothing else
        const int oc = icomp[std::size_t(o)];
        int gside = -1; bool on_owner = false, stray = false;
        for (int side = 0; side < 4; ++side)
            for (double t : {0.1, 0.3, 0.5, 0.7, 0.9}) {
                const double d = 0.5;
                const double x = side == 0 ? j.x0 - d : side == 1 ? j.x1 + d : j.x0 + t * j.w(), y = side == 2 ? j.y0 - d : side == 3 ? j.y1 + d : j.y0 + t * j.h();
                const int c = M.at(x, y);
                if (c < 0) continue;
                if (c == oc) on_owner = true;
                else if (net(c) == plane) { if (gside < 0 || is_plane[std::size_t(c)]) gside = c; }
                else stray = true;
            }
        const bool landed = oc >= 0 && on_owner && gside >= 0 && !stray && net(oc) != plane;
        if (landed) C.junctions_landed++;
        const std::string& on = L.islands[std::size_t(o)].first;
        const bool qubit = on[0] == 'Q';
        const GroundRef g = ground_ref((qubit ? "SQ" : "SC") + on.substr(1), landed ? gside : -1, &j);
        if (qubit) { C.qubit_squids++; C.qubit_squids_on_plane += g.on_plane; } else { C.coupler_squids++; C.coupler_squids_on_plane += g.on_plane; }
        if (g.on_plane) {
            (qubit ? (g.on_border ? C.qubit_squids_on_border : C.qubit_squids_on_bumps) : (g.on_border ? C.coupler_squids_on_border : C.coupler_squids_on_bumps))++;
            if (!g.on_border) {
                if (g.bump_dist_um < 0 || g.bump_dist_um > kMaxSquidBumpDist_um) C.squids_far_from_bump++;
                C.max_squid_bump_dist_um = std::max(C.max_squid_bump_dist_um, g.bump_dist_um);
            }
        }
    }
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const char t = L.islands[i].first[0];
        if ((t == 'Q' || t == 'C') && owned[i] != 1) C.junction_owner_errors++;
    }
    // shorted line ends; for qubit flux lines also the ground around the bump pad (the line's return)
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const std::string& n = L.islands[i].first;
        if (!is_shorted_line(n)) continue;
        const bool fq = n.rfind("FQ", 0) == 0;
        // a line that ends open is a component of its own with no ground in it: it has no ground reference
        const bool grounded = net(icomp[i]) == plane;
        const GroundRef g = ground_ref(n, grounded ? icomp[i] : -1);
        if (fq) C.qubit_flux_on_plane += g.on_plane; else { C.other_shorts++; C.other_shorts_on_plane += g.on_plane; }
        if (!fq) continue;
        bool cont = false;
        for (std::size_t b = 0; b < L.bumps.size() && b < L.bump_island.size(); ++b) {
            if (L.bump_island[b] != int(i)) continue;
            const Rect ring = L.bumps[b].grow(11 + 3);             // 3 um outside the pad's 52 um cut-out
            int metal = 0; cont = true;
            for (const auto& pt : {std::pair{ring.x0, ring.y0}, std::pair{ring.cx(), ring.y0}, std::pair{ring.x1, ring.y0}, std::pair{ring.x1, ring.cy()},
                                   std::pair{ring.x1, ring.y1}, std::pair{ring.cx(), ring.y1}, std::pair{ring.x0, ring.y1}, std::pair{ring.x0, ring.cy()}}) {
                const int c = M.at(pt.first, pt.second);
                if (c < 0) continue;
                ++metal; if (c != icomp[i]) cont = false;
            }
            if (metal < 2) cont = false;
        }
        C.flux_returns_continuous += cont && grounded;
    }
    for (int r : L.strap_returns) C.max_returns_per_strap = std::max(C.max_returns_per_strap, r);
    L.conn = C;
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
        L.jj_owner.resize(L.jj.size(), int(L.islands.size()) - 1);
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
        L.bump_island.push_back(int(L.islands.size()) - 1);
        // the readout bump's return: a ground bump on either side, above the serpentine's top row
        for (double dx : {-kGroundBumpPitch, kGroundBumpPitch}) { L.ground_bumps.push_back(rect_c(pcx + dx, pcy + 25, 30, 30)); L.ground_bump_line.push_back(int(L.islands.size()) - 1); }
        L.holes.insert(L.holes.end(), rh.begin(), rh.end());
        if (y > cy + 480) L.problems.push_back("resonator R" + std::to_string(q.id) + " exceeds its serpentine box");
        // flux line (inductive, shorted next to the SQUID) and XY line (open end at the (+x,-y) corner), both to bumps
        std::vector<Rect> fm, fh, xm, xh, am;
        const ControlLineGeom CL = qubit_control_lines(X, cx, cy, fm, fh, xm, xh, am, true, cs.separate_xy, cs.xy_tip);
        L.bumps.push_back(rect_c(CL.flux_bump.first, CL.flux_bump.second, 30, 30));
        L.islands.push_back({"FQ" + std::to_string(q.id), fm});
        L.bump_island.push_back(int(L.islands.size()) - 1);
        for (const auto& [gx, gy] : CL.flux_ground_bumps) { L.ground_bumps.push_back(rect_c(gx, gy, 30, 30)); L.ground_bump_line.push_back(int(L.islands.size()) - 1); }
        L.holes.insert(L.holes.end(), fh.begin(), fh.end());
        if (!xm.empty()) {
        L.bumps.push_back(rect_c(CL.xy_bump.first, CL.xy_bump.second, 30, 30));
        L.islands.push_back({"XY" + std::to_string(q.id), xm});
        L.bump_island.push_back(int(L.islands.size()) - 1);
        for (const auto& [gx, gy] : CL.xy_ground_bumps) { L.ground_bumps.push_back(rect_c(gx, gy, 30, 30)); L.ground_bump_line.push_back(int(L.islands.size()) - 1); }
        if (!am.empty()) L.islands.push_back({"A" + std::to_string(q.id), am});
        }
        L.holes.insert(L.holes.end(), xh.begin(), xh.end());
    }
    int nc = 0;
    std::vector<detail::Strap> stitch;
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
        const int coupler_island = int(L.islands.size()) - 1;
        if (!E.bus.empty()) L.islands.push_back({"B" + tag, E.bus});
        // stitching straps across the band: over the coupler body just beyond each T-bar and just before the pad, over
        // the bus at the same four places and in the middle. They tie the ground on both sides of the band to the strip
        // between coupler and bus, give the qubit's return current a path from one quadrant of its ground to the next
        // close to the qubit, and cut the long slots beside the bus into short ones (a closed slot conserves flux, and a
        // long one carries a flux line's field from one qubit to the next).
        {
            const CouplerGeom& g = cs.cg;
            const double yc0 = -g.y1 - g.w_body / 2 - g.gap, yc1 = -g.y1 + g.w_body / 2 + g.gap;          // coupler body cut-out
            const double yb0 = g.y1 - g.w_bus / 2 - g.gap_bus, yb1 = g.y1 + g.w_bus / 2 + g.gap_bus;       // bus cut-out
            const double e = detail::kStrapFootEnclosure + 1, F = detail::kStrapW;
            const bool bus_ok = !E.bus.empty() && yb0 - yc1 >= F + 2 * e;
            auto strap = [&](double x, double y0, double y1) {
                detail::Strap st;
                st.foot_a = rot90(rect_c(x, y0 - e - F / 2, F, F), k).shift(lo_q->x_um, lo_q->y_um);
                st.foot_b = rot90(rect_c(x, y1 + e + F / 2, F, F), k).shift(lo_q->x_um, lo_q->y_um);
                st.span = rot90(Rect{x - F / 2, y0 - e - F, x + F / 2, y1 + e + F}, k).shift(lo_q->x_um, lo_q->y_um);
                st.couplers_crossed = 1;
                stitch.push_back(st);
            };
            for (int end = 0; end < 2; ++end) {
                const double Lq = end ? Xh.L : Xl.L;
                auto at = [&](double x) { return end ? p - x : x; };
                const double xc = Lq + g.g_c + g.t_w + g.gap + e + F / 2 + 8;                              // beyond the coupler T-bar cut-out
                const double xp = p / 2 - g.pad / 2 - g.gap - e - F / 2 - 8;                               // before the pad's cut-out
                const double xbus = std::max(xc + F + 2 * e + 4, Lq + g.bus_gap() + g.t_w + g.gap_bus + e + F / 2 + 8);
                strap(at(xc), yc0, yc1);
                if (xp - xc >= 60) strap(at(xp), yc0, yc1);
                if (bus_ok) { strap(at(xbus), yb0, yb1); if (xp - F - 2 * e - 4 - xbus >= 60) strap(at(xp - F - 2 * e - 4), yb0, yb1); }
            }
            if (bus_ok) strap(p / 2, yb0, yb1);
        }
        L.holes.insert(L.holes.end(), E.holes.begin(), E.holes.end());
        // coupler SQUID on the pad's near side (towards the lower qubit: the other side belongs to the upper qubit's
        // flux line), flux line to a bump. Drawn for the far side and mirrored about the pad centre x = p / 2.
        const double ycen = -cs.cg.y1 - cs.cg.w_body / 2 - cs.cg.pad / 2;
        const double xj0 = p / 2 + cs.cg.pad / 2, xj1 = xj0 + cs.cg.gap;
        auto place = [&](const Rect& r) { return rot90(Rect{p - r.x1, r.y0, p - r.x0, r.y1}, k).shift(lo_q->x_um, lo_q->y_um); };
        L.jj.push_back(place({xj0 - 1, ycen - kSquidW / 2, xj1 + 1, ycen + kSquidW / 2}));
        L.jj_owner.push_back(coupler_island);
        // flux launcher (see kCollarW): the pad cut-out's edge x = xj1 is the ground edge; the line comes straight in
        std::vector<Rect> fm, fh;
        const double xs0 = xj1 + kCollarW, xs1 = xs0 + kSlotW, xtip = xs1 + kTipGround, xpad = xtip + kShortLen + 36;
        fm.push_back(place({xtip, ycen - 5, xpad, ycen + 5}));
        fh.push_back(place({xs0, ycen - cs.cg.slot_lo, xs1, ycen + cs.cg.slot_hi}));          // collar slot
        fh.push_back(place({xs1, ycen + 5, xpad, ycen + 11}));                               // the gap that runs into the slot
        fh.push_back(place({xtip, ycen - 5, xpad, ycen + 5}));                               // under the centre conductor (no etch)
        fh.push_back(place({xtip + kShortLen, ycen - 11, xpad, ycen - 5}));                  // the gap that stops short of the tip
        fm.push_back(place(rect_c(xpad, ycen, 36, 36))); fh.push_back(place(rect_c(xpad, ycen, 52, 52)));
        L.bumps.push_back(place(rect_c(xpad, ycen, 30, 30)));
        L.islands.push_back({"FC" + tag, fm});
        L.bump_island.push_back(int(L.islands.size()) - 1);
        for (double dy : {-kGroundBumpPitch, kGroundBumpPitch}) { L.ground_bumps.push_back(place(rect_c(xpad, ycen + dy, 30, 30))); L.ground_bump_line.push_back(int(L.islands.size()) - 1); }
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
            if (x == "XY" + std::to_string(q) && y == "A" + std::to_string(q)) return true;
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
    // ground straps: the stitching straps drawn with the couplers (kept where both feet stand on ground metal), then
    // whatever else is needed to tie the remaining ground to plane metal (the border component and every component
    // under a ground bump). Ground fragments no strap can reach and smaller than kRemoveFragmentUm2 are etched away.
    detail::MetalMap map1 = detail::metal_map(L.outline, L.etch);
    // ground-bump lattice: besides the bumps that flank the signal pads, a ground bump on every point of a
    // kBumpGrid lattice that keeps kBumpGridKeepout from every qubit, coupler and bus, kBumpGridLineKeepout from every
    // line, and has ground metal all round. They set how
    // many bumps share a plaquette's return currents (kBumpInductance_pH each, see line_budget).
    if (!L.ground_bumps.empty()) {
        const std::size_t n0 = L.ground_bumps.size();
        const int border = map1.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
        std::vector<Rect> taken = L.bumps; taken.insert(taken.end(), L.ground_bumps.begin(), L.ground_bumps.end());
        std::sort(taken.begin(), taken.end(), [](const Rect& a, const Rect& b) { return a.x0 < b.x0; });
        for (double y = L.outline.y0 + kBumpGrid; y <= L.outline.y1 - kBumpGrid + 1e-9; y += kBumpGrid)
            for (double x = L.outline.x0 + kBumpGrid; x <= L.outline.x1 - kBumpGrid + 1e-9; x += kBumpGrid) {
                const Rect b = rect_c(x, y, kBumpSize, kBumpSize);
                bool ok = true;
                I.each(b.grow(kBumpGridKeepout), [&](int i) { if (I.kind[std::size_t(i)] != 1) ok = false; });
                I.each(b.grow(kBumpGridLineKeepout), [&](int) { ok = false; });
                if (!ok) continue;
                for (auto it = std::lower_bound(taken.begin(), taken.end(), b.x0 - 2 * kBumpGrid, [](const Rect& a, double v) { return a.x0 < v; });
                     ok && it != taken.end() && it->x0 < b.x1 + 2 * kBumpGrid; ++it)
                    ok = std::abs(it->cx() - x) >= kBumpPitch + 10 || std::abs(it->cy() - y) >= kBumpPitch + 10;
                const Rect e = b.grow(kBumpEnclosure + 1);
                int comp = -2;
                for (int i = 0; ok && i <= 10; ++i) for (int j = 0; ok && j <= 10; ++j) {
                    const int c = map1.at(e.x0 + 0.05 + (e.w() - 0.1) * i / 10.0, e.y0 + 0.05 + (e.h() - 0.1) * j / 10.0);
                    ok = c >= 0 && (comp == -2 || c == comp); comp = c;
                }
                // the metal joined to the chip border has its reference in 2D: every third lattice point is enough there
                if (ok && comp == border && (std::lround((x - L.outline.x0) / kBumpGrid) % 3 != 0 || std::lround((y - L.outline.y0) / kBumpGrid) % 3 != 0)) ok = false;
                if (ok) { L.ground_bumps.push_back(b); L.ground_bump_line.push_back(-1); }
            }
        L.lattice_bumps = int(L.ground_bumps.size() - n0);
    }
    for (int round = 0; round < 3; ++round) {
        const int plane = map1.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
        const auto ic = island_components(L, map1);
        std::vector<char> eligible(std::size_t(map1.n_comp), 1), has_island(std::size_t(map1.n_comp), 0);
        for (std::size_t i = 0; i < ic.size(); ++i) if (ic[i] >= 0) {
            has_island[std::size_t(ic[i])] = 1;
            if (!is_shorted_line(L.islands[i].first)) eligible[std::size_t(ic[i])] = 0;
        }
        detail::UnionFindLite uf(map1.n_comp);
        for (const auto& gb : L.ground_bumps) { const int c = detail::foot_comp(map1, I, gb); if (c >= 0 && eligible[std::size_t(c)]) uf.unite(c, plane); }
        L.straps.clear();
        for (const auto& st : stitch) {
            const int a = detail::foot_comp(map1, I, st.foot_a), b = detail::foot_comp(map1, I, st.foot_b);
            if (a < 0 || b < 0 || !eligible[std::size_t(a)] || !eligible[std::size_t(b)]) continue;
            L.straps.push_back(st); uf.unite(a, b);
        }
        const auto more = detail::place_ground_straps(map1, I, eligible, plane, &uf);
        L.straps.insert(L.straps.end(), more.begin(), more.end());
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
    // ground bumps must stand on ground metal; a qubit's corner neck must not be narrower than the minimum width
    {
        int bad = 0;
        for (const auto& gb : L.ground_bumps) bad += detail::foot_comp(map1, I, gb.grow(3)) < 0;
        if (bad) L.problems.push_back(std::to_string(bad) + " ground bump(s) not on ground metal with 5 um around");
        for (const XmonParams* X : {&cs.meas, &cs.data})
            for (int corner = 1; corner <= 3; ++corner) {
                const double w = corner_neck(*X, corner);
                if (w > 0 && w < rules.etch_min_space_um) L.problems.push_back("corner " + std::to_string(corner) + " neck of a qubit is " + std::to_string(w) + " um wide");
            }
    }
    analyze_connectivity(L);          // nets and ground quality of the final mask
    L.etch_area_mm2 = area * 1e-6;   // with overlaps between cut-outs counted twice (upper bound)
    // minimum width / space on the merged etch layer (what the mask actually contains)
    L.rules = rules;
    L.drc = detail::width_space_check(L.etch, rules.etch_min_width_um, rules.etch_min_space_um, "etch");
    (void)nc;
    return L;
}

inline std::string connectivity_summary(const ChipLayout& L) {
    const auto& C = L.conn;
    char b[2400];
    std::snprintf(b, sizeof b, "metal connectivity %s (nets %s, ground %s): %d components; islands merged %d, on the plane %d; shorted lines on the "
                  "ground net %d/%d (qubit flux lines %d/%d); floating ground fragments %d; junctions landing owner-to-ground %d/%d (owner errors %d); "
                  "ground straps %d (landed %d; crossing %d lines, %d coupler/bus strips); small fragments removed %d\n"
                  "    ground: plane metal %d component(s), %d ground bumps; on plane metal with no strap in the path: qubit SQUID grounds %d/%d, "
                  "qubit flux-line shorts %d/%d (return continuous to the bump pad %d/%d), coupler SQUID grounds %d/%d, other line shorts %d/%d; "
                  "strapped fragments %d, on a single strap %d (cut straps %d), deepest %d straps (limit %d); circuits: deepest %d straps, "
                  "longest summed bridge %.0f um, most returns through one strap %d\n"
                  "    ground bumps %s: SQUID grounds on metal joined to the chip border in 2D: qubits %d, couplers %d; on plaquette metal grounded by "
                  "bumps alone: qubits %d, couplers %d (%d plaquettes, fewest bumps on one %d, rule %d; farthest SQUID ground from a bump %.0f um, "
                  "rule %.0f um); signal bumps without a ground bump within %.0f um: %d/%d; bump rule violations %d",
                  C.clean() ? "clean" : "NOT clean", C.nets_clean() ? "clean" : "NOT clean", C.ground_ok() ? "ok" : "NOT ok", C.components,
                  C.islands_merged, C.islands_on_plane, C.shorted_lines_on_plane, C.shorted_lines, C.flux_lines_on_plane, C.flux_lines,
                  C.floating_ground, C.junctions_landed, C.junctions, C.junction_owner_errors, C.straps, C.straps_landed, C.strap_lines_crossed,
                  C.strap_couplers_crossed, C.fragments_removed, C.plane_components, C.ground_bumps, C.qubit_squids_on_plane, C.qubit_squids,
                  C.qubit_flux_on_plane, C.flux_lines, C.flux_returns_continuous, C.flux_lines, C.coupler_squids_on_plane, C.coupler_squids,
                  C.other_shorts_on_plane, C.other_shorts, C.strapped_fragments, C.fragments_on_one_strap, C.cut_straps, C.max_strap_depth,
                  kMaxStrapDepth, C.max_circuit_depth, C.max_circuit_span_um, C.max_returns_per_strap, C.bumps_ok() ? "ok" : "NOT ok",
                  C.qubit_squids_on_border, C.coupler_squids_on_border, C.qubit_squids_on_bumps, C.coupler_squids_on_bumps, C.bump_plaquettes,
                  C.min_plaquette_bumps, kMinPlaquetteBumps, C.max_squid_bump_dist_um, kMaxSquidBumpDist_um, kSignalGroundBump_um,
                  C.signal_bumps_without_ground, C.signal_bumps, C.bump_rule_violations);
    return b;
}
inline Json connectivity_json(const ChipLayout& L) {
    const auto& C = L.conn;
    Json j = Json::object();
    j.set("clean", C.clean()); j.set("nets_clean", C.nets_clean()); j.set("ground_ok", C.ground_ok());
    j.set("components", C.components); j.set("islands", C.islands); j.set("islands_merged", C.islands_merged);
    j.set("islands_on_plane", C.islands_on_plane); j.set("shorted_lines", C.shorted_lines); j.set("shorted_lines_on_plane", C.shorted_lines_on_plane);
    j.set("flux_lines", C.flux_lines); j.set("flux_lines_on_plane", C.flux_lines_on_plane); j.set("floating_ground", C.floating_ground);
    j.set("junctions", C.junctions); j.set("junctions_landed", C.junctions_landed); j.set("junction_owner_errors", C.junction_owner_errors);
    j.set("ground_straps", C.straps); j.set("ground_straps_landed", C.straps_landed); j.set("strap_line_crossings", C.strap_lines_crossed); j.set("strap_coupler_crossings", C.strap_couplers_crossed);
    j.set("fragments_removed", C.fragments_removed);
    Json g = Json::object();
    g.set("plane_components", C.plane_components); g.set("ground_bumps", C.ground_bumps);
    g.set("qubit_squids", C.qubit_squids); g.set("qubit_squids_on_plane_metal", C.qubit_squids_on_plane);
    g.set("qubit_flux_shorts_on_plane_metal", C.qubit_flux_on_plane); g.set("qubit_flux_returns_continuous", C.flux_returns_continuous);
    g.set("coupler_squids", C.coupler_squids); g.set("coupler_squids_on_plane_metal", C.coupler_squids_on_plane);
    g.set("other_line_shorts", C.other_shorts); g.set("other_line_shorts_on_plane_metal", C.other_shorts_on_plane);
    g.set("strapped_fragments", C.strapped_fragments); g.set("fragments_on_one_strap", C.fragments_on_one_strap); g.set("cut_straps", C.cut_straps);
    g.set("max_strap_depth", C.max_strap_depth); g.set("max_strap_depth_allowed", kMaxStrapDepth);
    g.set("max_circuit_strap_depth", C.max_circuit_depth); g.set("max_circuit_bridge_um", C.max_circuit_span_um);
    g.set("max_returns_per_strap", C.max_returns_per_strap);
    Json b = Json::object();
    b.set("ok", C.bumps_ok()); b.set("on_border_metal", C.border_bumps); b.set("bump_grounded_plaquettes", C.bump_plaquettes);
    b.set("min_bumps_per_plaquette", C.min_plaquette_bumps); b.set("min_bumps_per_plaquette_required", kMinPlaquetteBumps);
    b.set("qubit_squids_on_border_metal", C.qubit_squids_on_border); b.set("qubit_squids_on_bump_grounded_metal", C.qubit_squids_on_bumps);
    b.set("coupler_squids_on_border_metal", C.coupler_squids_on_border); b.set("coupler_squids_on_bump_grounded_metal", C.coupler_squids_on_bumps);
    b.set("max_squid_to_bump_um", C.max_squid_bump_dist_um); b.set("max_squid_to_bump_allowed_um", kMaxSquidBumpDist_um);
    b.set("squids_far_from_bump", C.squids_far_from_bump); b.set("signal_bumps", C.signal_bumps);
    b.set("signal_bumps_without_ground_bump", C.signal_bumps_without_ground); b.set("signal_ground_bump_distance_um", kSignalGroundBump_um);
    b.set("bump_rule_violations", C.bump_rule_violations); b.set("assumed_bump_inductance_pH", kBumpInductance_pH);
    g.set("bumps", b);
    j.set("ground", g);
    return j;
}

// ---------------------------------------------------------------- flux coupling: thin-film magnetostatics by duality
// Model. The metal is a perfectly conducting film of zero thickness (penetration depth and kinetic inductance
// neglected), carrying in-plane sheet currents; the etched regions are apertures. The field of such currents derives
// from a magnetic scalar potential that is antisymmetric about the film plane: it is constant on every connected
// aperture (no current there), and its normal derivative vanishes on the metal (no flux through a superconductor). That
// is the electrostatic problem of conductors shaped like the apertures, in vacuum (Babinet duality), which
// extract_capacitance solves (eps_r = 1). With psi the "potential" of an aperture:
//   - the current flowing in the metal between two apertures is 2 (psi_a - psi_b);
//   - the flux through an aperture is mu0 / (2 eps0) times its dual "charge", and the flux through a part of an
//     aperture (a SQUID loop) is the same factor times the dual charge on that part (panel charges).
// A shorted CPW has two gap apertures; a current I in its centre conductor is psi_1 - psi_2 = I / 2, and the split of
// the return between the two grounds follows from flux conservation, Q_1 + Q_2 = 0. A ground strap is a
// superconducting bridge across the etch it spans; an aperture closed inside the window (by metal and straps) is
// flux-conserving (Q = 0, "floating"), SQUID loops included; an aperture that leaves the window is held at psi = 0.
// Check (unit test): the inductance per length of a long CPW comes out as mu0 K(k') / (4 K(k)).
// Not in the model: the junction leads (thin; open at zero circulating current), the bump and the carrier chip (the
// bump pad's cut-out is taken out of the window, so the line is fed at its pad), finite penetration depth, and
// everything outside the window.
namespace detail {
inline bool clip_rect(const Rect& r, const Rect& w, Rect& o) {
    o = {std::max(r.x0, w.x0), std::max(r.y0, w.y0), std::min(r.x1, w.x1), std::min(r.y1, w.y1)};
    return o.x1 > o.x0 + 1e-9 && o.y1 > o.y0 + 1e-9;
}
// Connected components of a set of rectangles: overlap or a shared edge of positive length joins; corner contact does not.
inline std::vector<int> rect_components(const std::vector<Rect>& rs, int* n_comp = nullptr) {
    std::vector<std::size_t> order(rs.size());
    for (std::size_t i = 0; i < rs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rs[a].x0 < rs[b].x0; });
    UnionFindLite uf(int(rs.size()));
    const double eps = 1e-7;
    for (std::size_t a = 0; a < order.size(); ++a) {
        const Rect& A = rs[order[a]];
        for (std::size_t b = a + 1; b < order.size() && rs[order[b]].x0 <= A.x1 + eps; ++b) {
            const Rect& B = rs[order[b]];
            const double ox = std::min(A.x1, B.x1) - std::max(A.x0, B.x0), oy = std::min(A.y1, B.y1) - std::max(A.y0, B.y0);
            if (ox > -eps && oy > -eps && (ox > eps || oy > eps)) uf.unite(int(order[a]), int(order[b]));
        }
    }
    std::vector<int> comp(rs.size());
    std::map<int, int> ids;
    for (std::size_t i = 0; i < rs.size(); ++i) { const int r = uf.find(int(i)); auto it = ids.find(r); if (it == ids.end()) it = ids.emplace(r, int(ids.size())).first; comp[i] = it->second; }
    if (n_comp) *n_comp = int(ids.size());
    return comp;
}
// dense linear solve (Gaussian elimination with partial pivoting); A row-major n x n
inline std::vector<double> solve_dense(std::vector<double> A, std::vector<double> b, int n) {
    for (int k = 0; k < n; ++k) {
        int p = k; for (int i = k + 1; i < n; ++i) if (std::abs(A[std::size_t(i) * n + k]) > std::abs(A[std::size_t(p) * n + k])) p = i;
        if (A[std::size_t(p) * n + k] == 0) throw std::runtime_error("solve_dense: singular system");
        if (p != k) { for (int j = 0; j < n; ++j) std::swap(A[std::size_t(k) * n + j], A[std::size_t(p) * n + j]); std::swap(b[std::size_t(k)], b[std::size_t(p)]); }
        for (int i = k + 1; i < n; ++i) {
            const double f = A[std::size_t(i) * n + k] / A[std::size_t(k) * n + k];
            if (f == 0) continue;
            for (int j = k; j < n; ++j) A[std::size_t(i) * n + j] -= f * A[std::size_t(k) * n + j];
            b[std::size_t(i)] -= f * b[std::size_t(k)];
        }
    }
    for (int i = n - 1; i >= 0; --i) { double v = b[std::size_t(i)]; for (int j = i + 1; j < n; ++j) v -= A[std::size_t(i) * n + j] * b[std::size_t(j)]; b[std::size_t(i)] = v / A[std::size_t(i) * n + i]; }
    return b;
}
} // namespace detail

// flux through an aperture per unit dual charge: mu0 / (2 eps0) = 1 / (2 eps0^2 c^2), in pH per fF (x 1 A)
inline double dual_flux_pH_per_fF() { return 1.0 / (2 * phys::eps0 * phys::eps0 * 299792458.0 * 299792458.0) * 1e-3; }

// Mutual inductance between the flux line of qubit `qubit_id` and SQUID loops, from the assembled mask.
// The window is the bounding box of the flux line (without its bump pad) and the qubit's SQUID loop, grown by `margin`;
// with reach > 0 it also takes in every junction region within `reach` of the loop, to report the flux those loops
// receive (crosstalk). The owner may be any island with a junction region, the line any flux line.
// M_cm: the flux line also links the path of the owner's own oscillating current (island -> SQUID -> ground, closed by
// the displacement current across the island's gap). With that current's return spread evenly in angle around the
// island and shared equally between the two ways round, the linked flux is the dual charge on the etch around the
// island weighted by (theta / 2 pi - 1/2), theta the position around the island measured from the SQUID (the polar
// angle, about the island's centre, of the ground metal nearest to each point of the etch). "The
// etch around the island" is the aperture that holds the SQUID loop, with every ground strap taken as a metal bridge
// across the etch it spans: beyond the straps nearest the island the return is taken to have left the etch edge.
// M_cm vanishes when the line's field is symmetric about the SQUID, and is large when the line lies in the return's
// stream. A series inductance M_cm between the junction current and a line of impedance Z0 limits T1 to
// Z0 L_J / (omega M_cm)^2.
struct FluxCoupling {
    bool ok{false};
    std::string note;
    double M_pH{0}, M_cm_pH{0}, M_cm_signed_pH{0};
    double I_half_mA{1e300};                          // current for half a flux quantum, Phi0 / (2 M)
    struct Other { std::string owner; double M_pH; };
    std::vector<Other> others;                        // other SQUID loops in the window
    double worst_other_ratio{0}; std::string worst_other;
    int apertures{0}, panels{0}; double seconds{0};
    Rect window;
};
inline FluxCoupling flux_coupling(const ChipLayout& L, const std::string& owner, const std::string& flux_line, const MeshOptions& mo, double margin = 150, double reach = 0,
                                  double local_um = 0) {
    FluxCoupling F;
    const int iq = L.island(owner), line_island = L.island(flux_line);
    if (iq < 0 || line_island < 0) { F.note = "no island " + owner + " or flux line " + flux_line; return F; }
    int own = -1;
    for (std::size_t k = 0; k < L.jj.size() && k < L.jj_owner.size(); ++k) if (L.jj_owner[k] == iq) { own = int(k); break; }
    if (own < 0) { F.note = owner + " has no junction region"; return F; }
    const Rect loop = L.jj[std::size_t(own)];
    // the line without its bump pad. The line is fed at the pad by its bump, and its return leaves the plane through
    // the ground bumps beside the pad: the cut-out ring around the pad is cut towards each of them (a 4 um strip from
    // the pad centre), so that each gap keeps the part of the ring on its side.
    std::vector<Rect> line, pad_rects;
    for (const auto& r : L.islands[std::size_t(line_island)].second) {
        bool pad = false;
        for (std::size_t b = 0; b < L.bumps.size() && b < L.bump_island.size(); ++b) if (L.bump_island[b] == line_island && r.covers(L.bumps[b])) pad = true;
        (pad ? pad_rects : line).push_back(r);
    }
    if (line.empty()) { F.note = "the flux line has no metal besides its pad"; return F; }
    std::vector<Rect> ext = line; ext.push_back(loop);
    ext.push_back(bbox_of(L.islands[std::size_t(iq)].second).grow(100));        // the owner's whole cut-out, up to the straps across its arms
    if (reach > 0) for (const auto& j : L.jj) if (std::hypot(j.cx() - loop.cx(), j.cy() - loop.cy()) <= reach) ext.push_back(j);
    std::vector<Rect> pads;                                  // the cuts
    if (!pad_rects.empty()) {
        const Rect pb = bbox_of(pad_rects), ring = pb.grow(8 + 0.5);           // the pad's cut-out (8 um around the 36 um pad)
        for (std::size_t g = 0; g < L.ground_bumps.size() && g < L.ground_bump_line.size(); ++g) {
            if (L.ground_bump_line[g] != line_island) continue;
            const double dx = L.ground_bumps[g].cx() - pb.cx(), dy = L.ground_bumps[g].cy() - pb.cy();
            if (std::abs(dx) >= std::abs(dy)) pads.push_back(dx > 0 ? Rect{pb.cx(), pb.cy() - 2, ring.x1, pb.cy() + 2} : Rect{ring.x0, pb.cy() - 2, pb.cx(), pb.cy() + 2});
            else pads.push_back(dy > 0 ? Rect{pb.cx() - 2, pb.cy(), pb.cx() + 2, ring.y1} : Rect{pb.cx() - 2, ring.y0, pb.cx() + 2, pb.cy()});
        }
        if (pads.empty()) {                                  // no ground bumps drawn: cut across the line, through the pad centre
            const Rect lb = bbox_of(line);
            const bool along_x = std::abs(lb.cx() - pb.cx()) >= std::abs(lb.cy() - pb.cy());
            pads.push_back(along_x ? Rect{pb.cx() - 2, ring.y0, pb.cx() + 2, ring.y1} : Rect{ring.x0, pb.cy() - 2, ring.x1, pb.cy() + 2});
        }
        ext.insert(ext.end(), pad_rects.begin(), pad_rects.end());
    }
    Rect W = bbox_of(ext).grow(margin);
    W = {std::max(W.x0, L.outline.x0), std::max(W.y0, L.outline.y0), std::min(W.x1, L.outline.x1), std::min(W.y1, L.outline.y1)};
    F.window = W;
    // apertures: the etch in the window
    // (with local_um > 0: only the etch within local_um of the line or of the owner, for a line and an owner far apart)
    const Rect near_line = bbox_of(line).grow(local_um), near_owner = bbox_of(L.islands[std::size_t(iq)].second).grow(local_um);
    std::vector<Rect> E;
    for (const auto& r : L.etch) {
        Rect o; if (!detail::clip_rect(r, W, o)) continue;
        if (local_um > 0 && !o.intersects(near_line) && !o.intersects(near_owner)) continue;
        std::vector<Rect> pcs{o};
        for (const auto& st : L.straps) {                    // a ground strap bridges the etch it spans
            if (!st.span.intersects(o)) continue;
            std::vector<Rect> nx; for (const auto& pc : pcs) { auto sb = detail::subtract(pc, st.span); nx.insert(nx.end(), sb.begin(), sb.end()); }
            pcs.swap(nx);
        }
        for (const auto& pd : pads) {
            if (!pd.intersects(o)) continue;
            std::vector<Rect> nx; for (const auto& pc : pcs) { auto sb = detail::subtract(pc, pd); nx.insert(nx.end(), sb.begin(), sb.end()); }
            pcs.swap(nx);
        }
        for (const auto& pc : pcs) if (pc.w() > 1e-6 && pc.h() > 1e-6) E.push_back(pc);
    }
    if (E.empty()) { F.note = "no etch in the window"; return F; }
    int n = 0;
    const std::vector<int> comp = detail::rect_components(E, &n);
    F.apertures = n;
    auto aperture_at = [&](double x, double y) { for (std::size_t i = 0; i < E.size(); ++i) if (E[i].contains(x, y)) return comp[i]; return -1; };
    enum { Floating = 0, Held = 1, Gap = 2 };
    std::vector<int> type(std::size_t(n), Floating);
    const double eb = 1e-6;
    for (std::size_t i = 0; i < E.size(); ++i)
        if (E[i].x0 <= W.x0 + eb || E[i].y0 <= W.y0 + eb || E[i].x1 >= W.x1 - eb || E[i].y1 >= W.y1 - eb) type[std::size_t(comp[i])] = Held;
    std::vector<char> has_loop(std::size_t(n), 0);
    for (const auto& j : L.jj) { const int a = aperture_at(j.cx(), j.cy()); if (a >= 0) has_loop[std::size_t(a)] = 1; }
    // the line's gap apertures: the etch bordering its metal
    std::map<int, int> touch;
    auto on_line = [&](double x, double y) { for (const auto& r : line) if (r.contains(x, y)) return true; return false; };
    for (const auto& r0 : line) {
        Rect r; if (!detail::clip_rect(r0, W, r)) continue;
        const double d = 0.5;
        auto probe = [&](double x, double y) { if (on_line(x, y)) return; const int a = aperture_at(x, y); if (a >= 0) touch[a]++; };
        for (double y = r.y0 + 1; y < r.y1; y += 2) { probe(r.x0 - d, y); probe(r.x1 + d, y); }      // every 2 um along each side
        for (double x = r.x0 + 1; x < r.x1; x += 2) { probe(x, r.y0 - d); probe(x, r.y1 + d); }
    }
    std::vector<std::pair<int, int>> gaps;                          // (contact samples, aperture), most contact first
    for (const auto& [a, c] : touch) gaps.push_back({c, a});
    std::sort(gaps.begin(), gaps.end(), [](const auto& x, const auto& y) { return x.first > y.first || (x.first == y.first && x.second < y.second); });
    if (gaps.empty()) { F.note = "no etch borders the flux line"; return F; }
    if (gaps.size() == 1) {
        F.ok = true; F.note = "the etch around the flux line is one aperture: the line carries no current past its end (open, or its gaps merge)";
        return F;
    }
    if (gaps.size() > 2) F.note = "more than two apertures border the flux line; the two with most contact are taken as its gaps";
    const int ga = gaps[0].second, gb = gaps[1].second;
    for (int g : {ga, gb}) if (!has_loop[std::size_t(g)]) type[std::size_t(g)] = Gap;
    std::vector<Conductor> cs(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) cs[std::size_t(k)].name = "a" + std::to_string(k);
    for (std::size_t i = 0; i < E.size(); ++i) cs[std::size_t(comp[i])].rects.push_back(E[i]);
    PanelSolution sol;
    const CapacitanceMatrix C = extract_capacitance(cs, 1.0, mo, &sol);
    F.panels = C.panels; F.seconds = C.seconds;
    // potentials for 1 A in the centre conductor
    std::vector<double> A(std::size_t(n) * n, 0.0), b(std::size_t(n), 0.0);
    auto row_charge = [&](int row, int k) { for (int j = 0; j < n; ++j) A[std::size_t(row) * n + j] += C.at(k, j); };
    for (int k = 0; k < n; ++k) {
        if (k == ga || k == gb) continue;
        if (type[std::size_t(k)] == Held) A[std::size_t(k) * n + k] = 1; else row_charge(k, k);
    }
    const bool fa = type[std::size_t(ga)] == Gap, fb = type[std::size_t(gb)] == Gap;
    if (fa && fb) {
        A[std::size_t(ga) * n + ga] = 1; A[std::size_t(ga) * n + gb] = -1; b[std::size_t(ga)] = 0.5;
        row_charge(gb, ga); row_charge(gb, gb);
    } else {                                           // a gap that is part of the held etch network: psi = 0 there
        A[std::size_t(ga) * n + ga] = 1; b[std::size_t(ga)] = fa ? 0.5 : 0.0;
        A[std::size_t(gb) * n + gb] = 1; b[std::size_t(gb)] = fb ? -0.5 : 0.0;
        if (!fa && !fb) { F.ok = true; F.note = "both gaps of the flux line open into the etch around the qubit: no current path"; return F; }
    }
    const std::vector<double> V = detail::solve_dense(A, b, n);
    const std::vector<double> q = sol.charge(V);
    const double K = dual_flux_pH_per_fF();
    auto loop_flux = [&](const Rect& lp) {
        double s = 0;
        for (std::size_t p = 0; p < sol.panels.size(); ++p) {
            Rect o; if (!detail::clip_rect(sol.panels[p].r, lp, o)) continue;
            s += q[p] * o.area() / sol.panels[p].r.area();
        }
        return K * s;
    };
    F.M_pH = std::abs(loop_flux(loop));
    const Rect qb = bbox_of(L.islands[std::size_t(iq)].second);
    const double qx = qb.cx(), qy = qb.cy();
    const int own_moat = aperture_at(loop.cx(), loop.cy());      // the etch around the owner, up to the nearest straps across it
    // position of a point of that etch around the island: the polar angle, about the island's centre, of the nearest
    // ground metal (the outer edge of the island's cut-out is star-shaped about the centre; the island's own outline,
    // with its T-bars and fingers, is not)
    std::vector<Rect> etch_w;
    for (const auto& r : L.etch) { Rect o; if (detail::clip_rect(r, W, o)) etch_w.push_back(o); }
    const auto& own_rects = L.islands[std::size_t(iq)].second;
    auto ground_at = [&](double x, double y) {
        if (!W.contains(x, y)) return false;
        for (const auto& r : etch_w) if (r.contains(x, y)) return false;
        for (const auto& r : own_rects) if (r.grow(1e-9).contains(x, y)) return false;
        return true;
    };
    auto ring_angle = [&](double x, double y) {
        double best = 1e300, bx = x, by = y;
        for (int d = 0; d < 16; ++d) {
            const double ux = std::cos(2 * PI * d / 16), uy = std::sin(2 * PI * d / 16);
            for (double t = 0.5; t < std::min(best, 80.0); t += 0.5)
                if (ground_at(x + t * ux, y + t * uy)) { best = t; bx = x + t * ux; by = y + t * uy; break; }
        }
        return std::atan2(by - qy, bx - qx);
    };
    const double th0 = ring_angle(loop.cx(), loop.cy());
    double cm = 0;
    for (std::size_t p = 0; p < sol.panels.size(); ++p) {
        if (sol.panels[p].cond != own_moat || type[std::size_t(own_moat)] == Gap) continue;
        double th = ring_angle(sol.panels[p].r.cx(), sol.panels[p].r.cy()) - th0;
        while (th < 0) th += 2 * PI;
        while (th >= 2 * PI) th -= 2 * PI;
        // the loop itself lies astride the SQUID: its flux counts half to each side, i.e. not at all
        Rect o; const double outside = detail::clip_rect(sol.panels[p].r, loop, o) ? 1 - o.area() / sol.panels[p].r.area() : 1.0;
        cm += q[p] * outside * (th / (2 * PI) - 0.5);
    }
    F.M_cm_pH = std::abs(K * cm); F.M_cm_signed_pH = K * cm * (loop_flux(loop) < 0 ? -1 : 1);   // sign relative to the loop's flux
    F.I_half_mA = F.M_pH > 0 ? phys::Phi0 / 2 / (F.M_pH * 1e-12) * 1e3 : 1e300;
    for (std::size_t k = 0; k < L.jj.size(); ++k) {
        if (int(k) == own || !W.covers(L.jj[k])) continue;
        const int o = k < L.jj_owner.size() ? L.jj_owner[k] : -1;
        FluxCoupling::Other x; x.owner = o >= 0 ? L.islands[std::size_t(o)].first : "?"; x.M_pH = std::abs(loop_flux(L.jj[k]));
        if (F.M_pH > 0 && x.M_pH / F.M_pH > F.worst_other_ratio) { F.worst_other_ratio = x.M_pH / F.M_pH; F.worst_other = x.owner; }
        F.others.push_back(x);
    }
    F.ok = true;
    return F;
}
inline FluxCoupling flux_coupling(const ChipLayout& L, int qubit_id, const MeshOptions& mo, double margin = 150, double reach = 0) {
    return flux_coupling(L, "Q" + std::to_string(qubit_id), "FQ" + std::to_string(qubit_id), mo, margin, reach);
}

// ---------------------------------------------------------------- control lines, from a window of the assembled mask
// Window extraction around one qubit: the qubit, its XY line, its flux line (cut free of its short), and every other
// qubit, coupler, bus and XY island in the window are conductors, clipped to the window; resonators and the other flux
// lines are ground metal, and all ground (fragments and straps included) is one ideal conductor. Coupler pads are in
// the window as drawn.
struct ControlWindow {
    bool ok{false};
    std::string note;
    int qubit{-1};
    Rect window;
    CapacitanceMatrix C;                                  // names: island names, "flux" (own flux line), "ground"
    std::vector<std::pair<double, double>> flux_zq;       // own flux line: (path distance from the short in um, charge in fF with the qubit at 1 V) per panel
    std::vector<std::string> clipped;                     // conductors cut by the window (their self-capacitance is not complete)
    bool has(const std::string& n) const { return std::find(C.names.begin(), C.names.end(), n) != C.names.end(); }
    double mutual(const std::string& a, const std::string& b) const { return has(a) && has(b) ? C.mutual(a, b) : 0.0; }
};
inline ControlWindow control_window(const ChipLayout& L, int qubit_id, const MeshOptions& mo, double margin = 200) {
    ControlWindow O; O.qubit = qubit_id;
    const std::string qn = "Q" + std::to_string(qubit_id), xn = "XY" + std::to_string(qubit_id), fn = "FQ" + std::to_string(qubit_id);
    const int iq = L.island(qn), ix = L.island(xn), ifq = L.island(fn);
    if (iq < 0) { O.note = "no qubit " + qn; return O; }
    std::vector<Rect> ext = L.islands[std::size_t(iq)].second;
    if (ix >= 0) ext.insert(ext.end(), L.islands[std::size_t(ix)].second.begin(), L.islands[std::size_t(ix)].second.end());
    if (ifq >= 0) ext.insert(ext.end(), L.islands[std::size_t(ifq)].second.begin(), L.islands[std::size_t(ifq)].second.end());
    Rect W = bbox_of(ext).grow(margin);
    W = {std::max(W.x0, L.outline.x0), std::max(W.y0, L.outline.y0), std::min(W.x1, L.outline.x1), std::min(W.y1, L.outline.y1)};
    O.window = W;
    std::vector<Rect> holes;
    for (const auto& h : L.holes) { Rect o; if (detail::clip_rect(h, W, o)) holes.push_back(o); }
    std::vector<Rect> etch;
    for (const auto& e : L.etch) { Rect o; if (detail::clip_rect(e, W, o)) etch.push_back(o); }
    auto in_etch = [&](double x, double y) { for (const auto& e : etch) if (e.contains(x, y)) return true; return false; };
    auto carve = [&](const Rect& m) {                    // metal that is ground: taken out of the cut-outs
        std::vector<Rect> nh;
        for (const auto& h : holes) { if (!h.intersects(m)) { nh.push_back(h); continue; } auto pcs = detail::subtract(h, m); nh.insert(nh.end(), pcs.begin(), pcs.end()); }
        holes.swap(nh);
    };
    // own flux line: the rectangles that touch ground metal (the short) are ground; the rest, kept 3 um clear of them, is "flux"
    std::vector<Rect> flux, shorts;
    if (ifq >= 0) {
        const auto& line = L.islands[std::size_t(ifq)].second;
        auto ground_at = [&](double x, double y) {
            if (!W.contains(x, y) || in_etch(x, y)) return false;
            for (const auto& r : line) if (r.grow(1e-9).contains(x, y)) return false;
            return true;
        };
        for (const auto& r0 : line) {
            Rect r; if (!detail::clip_rect(r0, W, r)) continue;
            bool touches = false;
            const double d = 0.25;
            for (double t : {0.1, 0.5, 0.9})
                touches = touches || ground_at(r.x0 - d, r.y0 + t * r.h()) || ground_at(r.x1 + d, r.y0 + t * r.h()) ||
                          ground_at(r.x0 + t * r.w(), r.y0 - d) || ground_at(r.x0 + t * r.w(), r.y1 + d);
            (touches ? shorts : flux).push_back(r);
        }
        for (const auto& sh : shorts) {
            std::vector<Rect> nx;
            for (const auto& r : flux) { auto pcs = detail::subtract(r, sh.grow(3.0)); for (const auto& pc : pcs) if (pc.w() > 0.5 && pc.h() > 0.5) nx.push_back(pc); }
            flux.swap(nx);
        }
    }
    std::vector<Conductor> cs;
    for (std::size_t i = 0; i < L.islands.size(); ++i) {
        const std::string& n = L.islands[i].first;
        if (int(i) == ifq) continue;
        std::vector<Rect> rs; bool cut = false;
        for (const auto& r0 : L.islands[i].second) { Rect r; if (detail::clip_rect(r0, W, r)) { rs.push_back(r); if (!W.covers(r0)) cut = true; } else cut = true; }
        if (rs.empty()) continue;
        if (is_shorted_line(n)) { for (const auto& r : rs) carve(r); continue; }
        Conductor c; c.name = n; c.rects = rs; cs.push_back(c);
        if (cut) O.clipped.push_back(n);
    }
    for (const auto& sh : shorts) carve(sh);
    if (!flux.empty()) { Conductor c; c.name = "flux"; c.rects = flux; cs.push_back(c); }
    Conductor G; G.name = "ground"; G.ground = true; G.outline = W; G.holes = holes;
    cs.push_back(G);
    PanelSolution sol;
    try { O.C = extract_capacitance(cs, 11.45, mo, &sol); }
    catch (const std::exception& e) { O.note = e.what(); return O; }
    // path distance from the short along the flux line's metal (shortest path over touching panels)
    if (!flux.empty() && !shorts.empty()) {
        const int fi = O.C.index("flux"), qi = O.C.index(qn);
        std::vector<std::size_t> fp;
        for (std::size_t p = 0; p < sol.panels.size(); ++p) if (sol.panels[p].cond == fi) fp.push_back(p);
        const std::size_t m = fp.size();
        std::vector<double> z(m, 1e300);
        std::set<std::pair<double, std::size_t>> pq;
        for (std::size_t a = 0; a < m; ++a) {
            const Rect& r = sol.panels[fp[a]].r;
            for (const auto& sh : shorts) if (r.grow(1e-6).intersects(sh.grow(3.0))) z[a] = std::min(z[a], std::hypot(r.cx() - sh.cx(), r.cy() - sh.cy()));
            if (z[a] < 1e299) pq.insert({z[a], a});
        }
        while (!pq.empty()) {
            const auto [d, a] = *pq.begin(); pq.erase(pq.begin());
            if (d > z[a]) continue;
            const Rect ra = sol.panels[fp[a]].r.grow(1e-6);
            for (std::size_t b2 = 0; b2 < m; ++b2) {
                const Rect& rb = sol.panels[fp[b2]].r;
                if (!ra.intersects(rb)) continue;
                const double nd = d + std::hypot(ra.cx() - rb.cx(), ra.cy() - rb.cy());
                if (nd < z[b2] - 1e-12) { z[b2] = nd; pq.insert({nd, b2}); }
            }
        }
        for (std::size_t a = 0; a < m; ++a) O.flux_zq.push_back({z[a] < 1e299 ? z[a] : 0.0, -sol.q_fF[std::size_t(qi)][fp[a]]});
    }
    O.ok = true;
    return O;
}

// Control-line budget of a chip. One window (control_window) and one flux-coupling solve (flux_coupling) for a
// representative measure qubit and data qubit (the qubit of that role with the most couplers, nearest the chip
// centre), and one flux-coupling solve for a coupler next to the measure qubit. For every qubit, at its own f01 and
// against its design T1 from analyze():
//   flux line, capacitive: 1/T1 = omega^2 C_eff^2 Z0 / C_sigma with C_eff the sum over the line's panels of q sin(beta z),
//     z the path distance from the short (a point of a shorted, matched line sees Re Z = Z0 sin^2(beta z)), beta at the
//     qubit's f01;
//   flux line, through the SQUID loop: Fermi's golden rule with the current noise 2 hbar omega / Z0 of the matched line
//     and <1|dH/dPhi|0> = pi E_JSigma d phi_01 / (Phi0 r), r = sqrt(cos^2(pi f) + d^2 sin^2(pi f)), phi_01 =
//     (2 E_C / E_J)^(1/4), at the design's flux bias f and SQUID asymmetry d. The device model has d = 0 and f = 0 on
//     the committed chips (no loss through the loop at first order), so the limit is quoted for an assumed junction
//     asymmetry of kAssumedSquidAsymmetry where the design's is smaller;
//   flux line, common mode: T1 = Z0 L_J / (omega M_cm)^2 (FluxCoupling::M_cm_pH), L_J = (Phi0 / 2 pi)^2 / E_J.
// The microwave drive.
//   Without separate XY lines (the committed chips) a qubit is driven through its flux line: a current I at the short
//     puts the emf omega M_cm I in series with the junction, so the common-mode mutual is the drive coupling and its T1
//     limit is the drive port's. Reported: the power on the chip for a kRabiMHz Rabi rate (I = hbar Omega Phi0 /
//     (2 pi E_J phi_01 M_cm) at the short, twice the incident wave's), and the drive a neighbouring qubit receives from
//     the line, M_cm(line, neighbour) / M_cm(neighbour's own line) (line_budget_crosstalk).
//   With XY islands in the layout ("XY<id>", optionally a floating antenna "A<id>"): 1/T1 = omega^2 C_eff^2 Z0 /
//     C_sigma, where C_eff comes from the window's Maxwell matrix with the line at 1 V, the target at 0 V and every
//     coupler, bus and antenna in the window as a passive node: C_nn (1 - f_n^2 / f^2) V_n + sum_m C_nm V_m = 0 (a
//     coupler resonates at its parked frequency f_n, a bus or antenna floats, f_n = 0), C_eff = C(XY, t) + sum_n
//     C(n, t) V_n. The same for each neighbouring qubit at its f01 (drive crosstalk, relative to its own drive), the
//     couplers' T1 limit through the line, C_c / (omega_c^2 C(XY, coupler)^2 Z0), and the drive power.
// Flux coupling: I(Phi0/2) of the qubits' and the coupler's flux lines; the coupler's T1 limit through the common mode
// of its own flux line; and the flux other SQUID loops receive (line_budget_crosstalk).
// Gates: every T1 limit at least kLineT1Margin times the design T1 (couplers: times the longest qubit design T1, the
// device model has no coupler T1); drive crosstalk at most kMaxDriveCrosstalk; drive power at most kMaxDrivePower_dBm;
// I(Phi0/2) at most kMaxHalfFluxQuantum_mA.
inline constexpr double kLineT1Margin = 10.0, kAssumedSquidAsymmetry = 0.1, kMaxDriveCrosstalk = 0.1, kMaxDrivePower_dBm = -55.0,
                        kMaxHalfFluxQuantum_mA = 2.0, kRabiMHz = 25.0;
// The flux-coupling solve is small (a few thousand panels): it always runs on this mesh; the solves between a line
// and a neighbour's SQUID (two clusters of etch, far apart) on the coarser one.
inline const MeshOptions kFluxMesh{1.5, 50.0, 1.6}, kFluxCrossMesh{2.5, 50.0, 2.0};
inline double line_T1_us(double C_sigma_fF, double C_c_fF, double f_GHz, double Z0 = 50) {
    const double w = 2 * PI * f_GHz * 1e9;
    return C_c_fF != 0 ? C_sigma_fF * 1e-15 / (w * w * C_c_fF * 1e-15 * C_c_fF * 1e-15 * Z0) * 1e6 : 1e300;
}
// T1 limit through the SQUID loop (see above): M in pH, E_JSigma and E_C in GHz, flux in Phi0
inline double flux_loop_T1_us(double M_pH, double EJsum_GHz, double EC_GHz, double d, double flux, double f_GHz, double Z0 = 50) {
    const double r = std::sqrt(std::pow(std::cos(PI * flux), 2) + d * d * std::pow(std::sin(PI * flux), 2));
    const double EJ = EJsum_GHz * r;
    const double m01 = PI * phys::h * EJsum_GHz * 1e9 * d * std::pow(2 * EC_GHz / EJ, 0.25) / (phys::Phi0 * r);     // J / Wb
    const double w = 2 * PI * f_GHz * 1e9, rate = std::pow(m01 * M_pH * 1e-12, 2) * 2 * w / (phys::hbar * Z0);
    return rate > 0 ? 1e6 / rate : 1e300;
}
// Shared return through the ground bumps. A line's return current leaves the chip through the ground bumps of the
// plane metal its short stands on (N bumps, kBumpInductance_pH each, taken in parallel: the metal is one node); that
// metal then moves against the top chip's ground by j omega (L_b / N) I. A qubit whose SQUID ground is on that metal
// has the other three quadrants of its cut-out on metal that does not move: kBumpOwnFraction of its mode current
// crosses the bumps, a mutual of kBumpOwnFraction L_b / N with the line, added to M_cm (worst sign). A neighbour has
// one quadrant on that metal: kBumpNeighbourFraction L_b / N. The fractions are the island's charge by quadrant, taken
// as equal; the ground straps between plaquettes are ignored (they lower the term). Metal with no ground bump: 0.
inline constexpr double kBumpOwnFraction = 0.75, kBumpNeighbourFraction = 0.25;
inline double bump_return_pH(const ChipLayout& L, const std::string& ground_ref) {
    for (const auto& g : L.grounds) if (g.name == ground_ref) return g.bumps > 0 ? kBumpInductance_pH / g.bumps : 0.0;
    return 0.0;
}
inline double flux_cm_T1_us(double M_cm_pH, double EJ_GHz, double f_GHz, double Z0 = 50) {
    const double LJ = std::pow(phys::Phi0 / (2 * PI), 2) / (phys::h * EJ_GHz * 1e9), w = 2 * PI * f_GHz * 1e9;
    return M_cm_pH > 0 ? Z0 * LJ / std::pow(w * M_cm_pH * 1e-12, 2) * 1e6 : 1e300;
}
// power (dBm) on an open-ended line for a Rabi rate f_R through C_eff: Omega = (2e / hbar) (C_eff / C_sigma) V n01, V = 2 V_in
inline double drive_power_dBm(double C_eff_fF, double C_sigma_fF, double n01, double rabi_MHz = kRabiMHz, double Z0 = 50) {
    if (C_eff_fF == 0 || n01 == 0) return 1e300;
    const double V = 2 * PI * rabi_MHz * 1e6 * phys::hbar * C_sigma_fF / (2 * phys::e * std::abs(C_eff_fF) * n01);
    return 10 * std::log10(V * V / (8 * Z0) / 1e-3);
}
// power (dBm) on a shorted line for a Rabi rate f_R through the common-mode mutual: hbar Omega = E_J phi_01 2 pi M_cm I / Phi0, I = 2 I_in
inline double flux_drive_power_dBm(double M_cm_pH, double EJ_GHz, double EC_GHz, double rabi_MHz = kRabiMHz, double Z0 = 50) {
    if (M_cm_pH <= 0 || EJ_GHz <= 0) return 1e300;
    const double I = phys::hbar * 2 * PI * rabi_MHz * 1e6 * phys::Phi0 / (2 * PI * phys::h * EJ_GHz * 1e9 * std::pow(2 * EC_GHz / EJ_GHz, 0.25) * M_cm_pH * 1e-12);
    return 10 * std::log10(I * I * Z0 / 8 / 1e-3);
}
struct LineCell {
    int qubit{-1};
    bool separate_xy{false};
    ControlWindow win;
    FluxCoupling flux;
    double C_sigma_fF{0}, C_xy_fF{0}, C_flux_open_fF{0};
    struct Node { std::string name; double C_self_fF, f_GHz; };          // passive nodes of the drive: couplers, buses, antenna
    std::vector<Node> nodes;
    double max_C_xy_coupler_fF{0}; std::string max_C_xy_coupler;
    // effective capacitance from the XY line to `target` at f; C_nt(n, t) may replace the window's C(n, t) (a
    // coupler or bus to one of its own qubits, which the window clips), or return a negative number to keep it
    double C_xy_eff(const std::string& target, double f_GHz, const std::function<double(const std::string&, const std::string&)>& C_nt = {}) const {
        const std::string xn = "XY" + std::to_string(qubit);
        const int n = int(nodes.size());
        std::vector<double> A(std::size_t(n) * n, 0.0), b(std::size_t(n), 0.0);
        for (int i = 0; i < n; ++i) {
            const Node& a = nodes[std::size_t(i)];
            A[std::size_t(i) * n + i] = a.C_self_fF * (1 - (a.f_GHz > 0 ? a.f_GHz * a.f_GHz / (f_GHz * f_GHz) : 0.0));
            for (int j = 0; j < n; ++j) if (j != i) A[std::size_t(i) * n + j] = -win.mutual(a.name, nodes[std::size_t(j)].name);
            b[std::size_t(i)] = win.mutual(a.name, xn);
        }
        const std::vector<double> V = n ? detail::solve_dense(A, b, n) : std::vector<double>{};
        double c = win.mutual(xn, target);
        for (int i = 0; i < n; ++i) {
            double cnt = C_nt ? C_nt(nodes[std::size_t(i)].name, target) : -1.0;
            if (cnt < 0) cnt = win.mutual(nodes[std::size_t(i)].name, target);
            c += cnt * V[std::size_t(i)];
        }
        return c;
    }
    double C_flux_eff(double f_GHz) const {
        const double beta = 2 * PI * f_GHz * 1e9 * std::sqrt(0.5 * (11.45 + 1)) / 299792458.0 * 1e-6;   // rad / um
        double c = 0; for (const auto& [z, q] : win.flux_zq) c += q * std::sin(beta * z);
        return c;
    }
};
struct LineBudget {
    bool evaluated{false}, crosstalk_evaluated{false};
    std::string note;
    bool separate_xy{false};                             // the layout has XY lines; otherwise the drive goes through the flux lines
    LineCell meas, data;
    std::string coupler; FluxCoupling coupler_flux;      // representative coupler and its flux line
    double f_coupler_GHz{0};
    double min_ratio_drive{1e300}, min_ratio_flux_cap{1e300}, min_ratio_flux_loop{1e300}, min_ratio_flux_cm{1e300};
    double min_T1_drive_us{1e300}, min_T1_flux_cap_us{1e300}, min_T1_flux_loop_us{1e300}, min_T1_flux_cm_us{1e300};
    int worst_drive_qubit{-1}, worst_flux_cap_qubit{-1}, worst_flux_loop_qubit{-1}, worst_flux_cm_qubit{-1};
    double min_T1_coupler_us{1e300}, ratio_coupler{1e300}; std::string worst_coupler;      // coupler T1 limit through an XY line, or through its own flux line
    double max_drive_crosstalk{0}; std::string worst_crosstalk;
    double drive_power_dBm{-1e300};                      // worst (highest) power for the Rabi rate
    double min_mediated_fraction{1e300}, max_mediated_fraction{-1e300};   // separate XY: (C_eff - C_direct) / C_direct over the qubits
    double max_I_half_mA{0}, min_M_pH{1e300}, max_flux_crosstalk{0}; std::string worst_flux_crosstalk;
    // shared return through the ground bumps (see bump_return_pH): the largest term over the qubits, the fewest bumps
    // on a qubit's plaquette, and the coupler's term
    double max_M_bump_pH{0}, M_bump_coupler_pH{0}; int min_bumps{-1};
    bool drive_ok() const { return min_ratio_drive >= kLineT1Margin && (!separate_xy || ratio_coupler >= kLineT1Margin) && crosstalk_evaluated && max_drive_crosstalk <= kMaxDriveCrosstalk && drive_power_dBm <= kMaxDrivePower_dBm; }
    bool flux_ok() const { return min_ratio_flux_cap >= kLineT1Margin && min_ratio_flux_loop >= kLineT1Margin && min_ratio_flux_cm >= kLineT1Margin && (separate_xy || ratio_coupler >= kLineT1Margin) && max_I_half_mA <= kMaxHalfFluxQuantum_mA; }
    bool ok() const { return evaluated && drive_ok() && flux_ok(); }
};
namespace detail {
inline std::pair<int, int> pair_ends(const std::string& s) {      // "C5-8" -> (5, 8)
    const auto d = s.find('-'); if (d == std::string::npos || s.size() < 2) return {-1, -1};
    return {std::stoi(s.substr(1, d - 1)), std::stoi(s.substr(d + 1))};
}
}
inline LineBudget line_budget(const ChipDesign& D, const ChipReport& R, const ChipLayout& L, const LayoutMetrics& m, double f_coupler_GHz,
                              const MeshOptions& mo, double Z0 = 50) {
    LineBudget B; B.f_coupler_GHz = f_coupler_GHz;
    // representative qubits
    Rect bb{1e300, 1e300, -1e300, -1e300};
    for (const auto& q : D.qubits) { bb.x0 = std::min(bb.x0, q.x_um); bb.y0 = std::min(bb.y0, q.y_um); bb.x1 = std::max(bb.x1, q.x_um); bb.y1 = std::max(bb.y1, q.y_um); }
    std::map<int, int> deg;
    for (const auto& c : D.couplings) { deg[c.a]++; deg[c.b]++; }
    auto pick = [&](bool meas) {
        int best = -1; double bd = 0; int bdeg = -1;
        for (const auto& q : D.qubits) {
            if ((q.role != "data") != meas) continue;
            const double d = std::hypot(q.x_um - bb.cx(), q.y_um - bb.cy());
            if (deg[q.id] > bdeg || (deg[q.id] == bdeg && d < bd - 1e-9)) { best = q.id; bd = d; bdeg = deg[q.id]; }
        }
        return best;
    };
    auto cell = [&](LineCell& c, bool meas) {
        c.qubit = pick(meas);
        if (c.qubit < 0) { B.note = "no qubit of one role"; return false; }
        c.win = control_window(L, c.qubit, mo);
        if (!c.win.ok) { B.note = c.win.note; return false; }
        c.flux = flux_coupling(L, c.qubit, kFluxMesh);
        if (!c.flux.ok) { B.note = c.flux.note; return false; }
        const std::string qn = "Q" + std::to_string(c.qubit), xn = "XY" + std::to_string(c.qubit);
        c.separate_xy = c.win.has(xn);
        c.C_sigma_fF = meas ? m.C_sigma_m : m.C_sigma_d;
        c.C_xy_fF = c.win.mutual(xn, qn);
        for (const auto& [z, q] : c.win.flux_zq) c.C_flux_open_fF += q;
        for (const auto& n : c.win.C.names) {
            if (n[0] != 'C' && n[0] != 'B' && n[0] != 'A') continue;
            const bool cut = std::find(c.win.clipped.begin(), c.win.clipped.end(), n) != c.win.clipped.end();
            c.nodes.push_back({n, cut && n[0] != 'A' ? (n[0] == 'C' ? m.C_c : m.C_bus) : c.win.C.self(n), n[0] == 'C' ? f_coupler_GHz : 0.0});
            if (n[0] == 'C' && std::abs(c.win.mutual(xn, n)) > c.max_C_xy_coupler_fF) { c.max_C_xy_coupler_fF = std::abs(c.win.mutual(xn, n)); c.max_C_xy_coupler = n; }
        }
        return true;
    };
    if (!cell(B.meas, true) || !cell(B.data, false)) return B;
    B.separate_xy = B.meas.separate_xy || B.data.separate_xy;
    // representative coupler: one of the measure qubit's
    for (const auto& c : D.couplings) if (c.a == B.meas.qubit || c.b == B.meas.qubit) { B.coupler = std::to_string(c.a) + "-" + std::to_string(c.b); break; }
    if (!B.coupler.empty()) B.coupler_flux = flux_coupling(L, "C" + B.coupler, "FC" + B.coupler, kFluxMesh);
    // C(n, target): a coupler or bus to one of its own two qubits from the cell metrics (the window may clip it), else from the window
    auto C_nt = [&](const std::string& n, const std::string& t) {
        if (n[0] != 'C' && n[0] != 'B') return -1.0;
        const auto [a, b] = detail::pair_ends(n);
        const int tq = std::stoi(t.substr(1));
        return a == tq || b == tq ? (n[0] == 'C' ? m.C_qc : m.C_qb) : -1.0;
    };
    double T1max = 0; for (const auto& q : R.q) T1max = std::max(T1max, q.T1_us);
    for (std::size_t i = 0; i < D.qubits.size() && i < R.q.size(); ++i) {
        const bool meas = D.qubits[i].role != "data";
        const LineCell& c = meas ? B.meas : B.data;       // the cell's own-qubit couplings stand for this qubit's
        const auto& rq = R.q[i];
        const double f = rq.f01;
        const double tcap = line_T1_us(c.C_sigma_fF, c.C_flux_eff(f), f, Z0);
        const double d = std::max(std::abs(D.qubits[i].squid_d), kAssumedSquidAsymmetry);
        const double EJsum = EJ_from_Rn_GHz(D.qubits[i].Rn_kOhm, D.gap_ueV);
        const double mb = bump_return_pH(L, "SQ" + std::to_string(D.qubits[i].id));
        B.max_M_bump_pH = std::max(B.max_M_bump_pH, kBumpOwnFraction * mb);
        for (const auto& g : L.grounds) if (g.name == "SQ" + std::to_string(D.qubits[i].id)) B.min_bumps = B.min_bumps < 0 ? g.bumps : std::min(B.min_bumps, g.bumps);
        const double tloop = flux_loop_T1_us(c.flux.M_pH, EJsum, rq.EC, d, D.qubits[i].flux, f, Z0),
                     tcm = flux_cm_T1_us(c.flux.M_cm_pH + kBumpOwnFraction * mb, rq.EJ, f, Z0);
        auto upd = [&](double t, double& mt, double& mr, int& wq) { mt = std::min(mt, t); if (t / rq.T1_us < mr) { mr = t / rq.T1_us; wq = D.qubits[i].id; } };
        upd(tcap, B.min_T1_flux_cap_us, B.min_ratio_flux_cap, B.worst_flux_cap_qubit);
        upd(tloop, B.min_T1_flux_loop_us, B.min_ratio_flux_loop, B.worst_flux_loop_qubit);
        upd(tcm, B.min_T1_flux_cm_us, B.min_ratio_flux_cm, B.worst_flux_cm_qubit);
        if (B.separate_xy) {
            const double ceff = c.C_xy_eff("Q" + std::to_string(c.qubit), f, C_nt);
            upd(line_T1_us(c.C_sigma_fF, ceff, f, Z0), B.min_T1_drive_us, B.min_ratio_drive, B.worst_drive_qubit);
            B.drive_power_dBm = std::max(B.drive_power_dBm, drive_power_dBm(ceff, c.C_sigma_fF, rq.n01));
            if (c.C_xy_fF != 0) {
                const double mf = (ceff - c.C_xy_fF) / c.C_xy_fF;
                B.min_mediated_fraction = std::min(B.min_mediated_fraction, mf); B.max_mediated_fraction = std::max(B.max_mediated_fraction, mf);
            }
        } else {                                           // driven through the flux line: the common mode is the drive port
            upd(tcm, B.min_T1_drive_us, B.min_ratio_drive, B.worst_drive_qubit);
            B.drive_power_dBm = std::max(B.drive_power_dBm, flux_drive_power_dBm(c.flux.M_cm_pH - kBumpOwnFraction * mb, rq.EJ, rq.EC));
        }
    }
    for (const LineCell* c : {&B.meas, &B.data}) {
        if (B.separate_xy)
            for (const auto& n : c->win.C.names) if (n[0] == 'C') {              // coupler decay into the XY line
                const double t = line_T1_us(m.C_c, c->win.mutual("XY" + std::to_string(c->qubit), n), f_coupler_GHz, Z0);
                if (t < B.min_T1_coupler_us) { B.min_T1_coupler_us = t; B.worst_coupler = n; }
            }
        B.max_I_half_mA = std::max(B.max_I_half_mA, c->flux.I_half_mA); B.min_M_pH = std::min(B.min_M_pH, c->flux.M_pH);
    }
    if (B.coupler_flux.ok) {
        B.max_I_half_mA = std::max(B.max_I_half_mA, B.coupler_flux.I_half_mA);
        if (!B.separate_xy) {                              // coupler decay into its own flux line (common mode)
            const FluxTransmon ct = FluxTransmon::from_target(D.coupler.f_max_GHz, D.coupler.alpha_GHz);
            B.M_bump_coupler_pH = kBumpOwnFraction * bump_return_pH(L, "SC" + B.coupler);
            B.min_T1_coupler_us = flux_cm_T1_us(B.coupler_flux.M_cm_pH + B.M_bump_coupler_pH, ct.EJ(ct.flux_for_f01(f_coupler_GHz)), f_coupler_GHz, Z0);
            B.worst_coupler = "C" + B.coupler;
        }
    } else if (!B.coupler.empty()) { B.note = "coupler flux line: " + B.coupler_flux.note; B.max_I_half_mA = 1e300; }
    B.ratio_coupler = T1max > 0 ? B.min_T1_coupler_us / T1max : 1e300;
    B.evaluated = true;
    return B;
}
// Crosstalk of the representative qubits' flux lines, kept apart from line_budget because of its cost: one
// flux-coupling solve per victim, on the etch within kCrossLocal of the line and of the victim (two clusters). The
// victims are the two coupled qubits nearest the line's pad and the qubit's couplers towards them.
//   Flux: the flux in the victim's SQUID loop as a fraction of the own loop's.
//   Drive (through the flux lines): the common-mode mutual between the line and a neighbouring qubit, relative to that
//     of the neighbour's own line. With XY lines: C_eff to each coupled neighbour at its f01 relative to its own.
inline constexpr double kCrossLocal = 150.0;
inline void line_budget_crosstalk(LineBudget& B, const ChipDesign& D, const ChipReport& R, const ChipLayout& L, const LayoutMetrics& m) {
    if (!B.evaluated) return;
    std::map<int, int> index;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) index[D.qubits[i].id] = int(i);
    auto C_nt = [&](const std::string& n, const std::string& t) {
        if (n[0] != 'C' && n[0] != 'B') return -1.0;
        const auto [a, b] = detail::pair_ends(n);
        const int tq = std::stoi(t.substr(1));
        return a == tq || b == tq ? (n[0] == 'C' ? m.C_qc : m.C_qb) : -1.0;
    };
    for (LineCell* c : {&B.meas, &B.data}) {
        if (c->qubit < 0) continue;
        const std::string line = "FQ" + std::to_string(c->qubit);
        const int il = L.island(line);
        std::pair<double, double> pad{D.qubits[std::size_t(index[c->qubit])].x_um, D.qubits[std::size_t(index[c->qubit])].y_um};
        for (std::size_t b = 0; b < L.bumps.size() && b < L.bump_island.size(); ++b) if (L.bump_island[b] == il) pad = {L.bumps[b].cx(), L.bumps[b].cy()};
        std::vector<std::pair<double, int>> nb;            // coupled neighbours, nearest to the pad first
        for (const auto& cp : D.couplings) {
            const int j = cp.a == c->qubit ? cp.b : cp.b == c->qubit ? cp.a : -1;
            if (j < 0 || !index.count(j)) continue;
            nb.push_back({std::hypot(D.qubits[std::size_t(index[j])].x_um - pad.first, D.qubits[std::size_t(index[j])].y_um - pad.second), j});
        }
        std::sort(nb.begin(), nb.end());
        c->flux.others.clear(); c->flux.worst_other.clear(); c->flux.worst_other_ratio = 0;
        auto flux_victim = [&](const std::string& owner, const FluxCoupling& X) {
            if (!X.ok) return;
            c->flux.others.push_back({owner, X.M_pH});
            const double r = c->flux.M_pH > 0 ? X.M_pH / c->flux.M_pH : 0.0;
            if (r > c->flux.worst_other_ratio) { c->flux.worst_other_ratio = r; c->flux.worst_other = owner; }
        };
        for (std::size_t k = 0; k < nb.size(); ++k) {
            const int j = nb[k].second;
            const bool jmeas = D.qubits[std::size_t(index[j])].role != "data";
            const LineCell& own = jmeas ? B.meas : B.data;         // the neighbour's own drive: its role's cell
            const std::string qj = "Q" + std::to_string(j);
            double x = 0;
            if (B.separate_xy) {
                const double fj = R.q[std::size_t(index[j])].f01;
                const double there = c->C_xy_eff(qj, fj, C_nt), self = own.C_xy_eff("Q" + std::to_string(own.qubit), fj, C_nt);
                x = self != 0 ? std::abs(there / self) : 1e300;
            }
            if (k < 2) {                                           // the two nearest: flux and (through the flux lines) drive
                const FluxCoupling X = flux_coupling(L, qj, line, kFluxCrossMesh, 60, 0, kCrossLocal);
                flux_victim(qj, X);
                if (!B.separate_xy) {                              // with the shared bump return on both sides, worst signs
                    const double mb = bump_return_pH(L, "SQ" + std::to_string(c->qubit)), mbj = bump_return_pH(L, "S" + qj);
                    const double self = own.flux.M_cm_pH - kBumpOwnFraction * mbj;
                    x = !X.ok || self <= 0 ? 1e300 : (X.M_cm_pH + kBumpNeighbourFraction * mb) / self;
                }
                for (const std::string& tag : {std::to_string(c->qubit) + "-" + std::to_string(j), std::to_string(j) + "-" + std::to_string(c->qubit)})
                    if (L.island("C" + tag) >= 0) flux_victim("C" + tag, flux_coupling(L, "C" + tag, line, kFluxCrossMesh, 60, 0, kCrossLocal));
            }
            if (x > B.max_drive_crosstalk) { B.max_drive_crosstalk = x; B.worst_crosstalk = (B.separate_xy ? "XY" : "FQ") + std::to_string(c->qubit) + " -> " + qj; }
        }
        if (c->flux.worst_other_ratio > B.max_flux_crosstalk) { B.max_flux_crosstalk = c->flux.worst_other_ratio; B.worst_flux_crosstalk = line + " -> " + c->flux.worst_other; }
    }
    B.crosstalk_evaluated = true;
}
inline std::string line_budget_summary(const LineBudget& B) {
    if (!B.evaluated) return "control lines: not evaluated (" + B.note + ")";
    auto x = [](double v) { return std::min(v, 9999.0); };
    char b[3200], d[1200];
    if (B.separate_xy)
        std::snprintf(d, sizeof d, "    drive (XY lines) %s: C(XY, qubit) %.4f / %.4f fF; largest C(XY, coupler) %.4f / %.4f fF (%s / %s); coupler- and bus-mediated "
                      "part of the drive %+.0f%% .. %+.0f%%; T1 limit >= %.3g ms = %.1fx the design T1 (Q%d); coupler T1 limit %.3g ms = %.1fx (%s); drive crosstalk "
                      "%s%.0f%% (%s); %.0f MHz Rabi at <= %.1f dBm",
                      B.drive_ok() ? "ok" : "NOT ok", B.meas.C_xy_fF, B.data.C_xy_fF, B.meas.max_C_xy_coupler_fF, B.data.max_C_xy_coupler_fF,
                      B.meas.max_C_xy_coupler.c_str(), B.data.max_C_xy_coupler.c_str(), 100 * B.min_mediated_fraction, 100 * B.max_mediated_fraction,
                      B.min_T1_drive_us / 1e3, x(B.min_ratio_drive), B.worst_drive_qubit, B.min_T1_coupler_us / 1e3, x(B.ratio_coupler), B.worst_coupler.c_str(),
                      B.crosstalk_evaluated ? "<= " : "not evaluated ", 100 * x(B.max_drive_crosstalk), B.worst_crosstalk.c_str(), kRabiMHz, B.drive_power_dBm);
    else
        std::snprintf(d, sizeof d, "    drive (through the flux lines, no XY lines) %s: common-mode mutual %.3f / %.3f pH; T1 limit >= %.3g ms = %.1fx the design T1 (Q%d); "
                      "%.0f MHz Rabi at <= %.1f dBm (limit %.0f); drive on a neighbouring qubit %s%.1f%% of its own (%s; limit %.0f%%)",
                      B.drive_ok() ? "ok" : "NOT ok", B.meas.flux.M_cm_pH, B.data.flux.M_cm_pH, B.min_T1_drive_us / 1e3, x(B.min_ratio_drive), B.worst_drive_qubit,
                      kRabiMHz, B.drive_power_dBm, kMaxDrivePower_dBm, B.crosstalk_evaluated ? "<= " : "not evaluated ", 100 * x(B.max_drive_crosstalk),
                      B.worst_crosstalk.c_str(), 100 * kMaxDriveCrosstalk);
    std::snprintf(b, sizeof b,
                  "control lines %s (Q%d measure / Q%d data, coupler C%s parked at %.3f GHz):\n%s\n"
                  "    flux %s: M %.3f / %.3f pH, coupler %.3f pH; I(Phi0/2) %.3g / %.3g mA, coupler %.3g mA (limit %.0f); flux in other SQUID loops %s%.1f%% "
                  "of the own (%s); T1 limits: capacitive (%.2f / %.2f fF open) >= %.3g ms = %.1fx (Q%d); through the loop at asymmetry >= %.2f %.3g ms = %.1fx "
                  "(Q%d); common mode %.3g ms = %.1fx (Q%d)%s; budget %.0fx",
                  B.ok() ? "ok" : "NOT ok", B.meas.qubit, B.data.qubit, B.coupler.c_str(), B.f_coupler_GHz, d,
                  B.flux_ok() ? "ok" : "NOT ok", B.meas.flux.M_pH, B.data.flux.M_pH, B.coupler_flux.M_pH, x(B.meas.flux.I_half_mA), x(B.data.flux.I_half_mA),
                  x(B.coupler_flux.I_half_mA), kMaxHalfFluxQuantum_mA, B.crosstalk_evaluated ? "<= " : "not evaluated ", 100 * B.max_flux_crosstalk,
                  B.worst_flux_crosstalk.c_str(), B.meas.C_flux_open_fF, B.data.C_flux_open_fF, B.min_T1_flux_cap_us / 1e3, x(B.min_ratio_flux_cap),
                  B.worst_flux_cap_qubit, kAssumedSquidAsymmetry, std::min(B.min_T1_flux_loop_us / 1e3, 1e9), x(B.min_ratio_flux_loop), B.worst_flux_loop_qubit,
                  std::min(B.min_T1_flux_cm_us / 1e3, 1e9), x(B.min_ratio_flux_cm), B.worst_flux_cm_qubit,
                  B.separate_xy ? "" : (std::string("; coupler through its own flux line (M_cm ") + std::to_string(B.coupler_flux.M_cm_pH).substr(0, 5) + " pH) " +
                                        std::to_string(std::min(B.min_T1_coupler_us / 1e3, 1e9)).substr(0, 6) + " ms = " + std::to_string(x(B.ratio_coupler)).substr(0, 5) + "x the longest qubit T1").c_str(),
                  kLineT1Margin);
    char e[400];
    std::snprintf(e, sizeof e, "\n    shared return through the ground bumps (%.1f pH per bump assumed, fewest on a qubit's metal %d): up to %.3f pH added to M_cm "
                  "of a qubit, %.3f pH for the coupler; included in the common-mode T1 limits, the drive power and the drive crosstalk (worst signs)",
                  kBumpInductance_pH, B.min_bumps, B.max_M_bump_pH, B.M_bump_coupler_pH);
    return std::string(b) + e;
}
inline Json line_budget_json(const LineBudget& B) {
    Json j = Json::object();
    j.set("bump_inductance_assumed_pH", kBumpInductance_pH); j.set("min_ground_bumps_on_qubit_metal", B.min_bumps);
    j.set("max_M_bump_return_pH", B.max_M_bump_pH); j.set("M_bump_return_coupler_pH", B.M_bump_coupler_pH);
    j.set("evaluated", B.evaluated); j.set("ok", B.ok());
    if (!B.evaluated) { j.set("note", B.note); return j; }
    auto fin = [](double v) { return std::min(v, 1e12); };
    j.set("drive_ok", B.drive_ok()); j.set("flux_ok", B.flux_ok()); j.set("margin", kLineT1Margin);
    j.set("separate_xy_lines", B.separate_xy); j.set("crosstalk_evaluated", B.crosstalk_evaluated);
    j.set("measure_qubit", B.meas.qubit); j.set("data_qubit", B.data.qubit); j.set("coupler", "C" + B.coupler); j.set("coupler_f_GHz", B.f_coupler_GHz);
    if (B.separate_xy) {
        j.set("C_xy_measure_fF", B.meas.C_xy_fF); j.set("C_xy_data_fF", B.data.C_xy_fF);
        j.set("max_C_xy_coupler_measure_fF", B.meas.max_C_xy_coupler_fF); j.set("max_C_xy_coupler_data_fF", B.data.max_C_xy_coupler_fF);
        j.set("mediated_drive_fraction_min", B.min_mediated_fraction); j.set("mediated_drive_fraction_max", B.max_mediated_fraction);
    }
    j.set("min_T1_drive_us", fin(B.min_T1_drive_us)); j.set("min_ratio_drive", fin(B.min_ratio_drive)); j.set("worst_drive_qubit", B.worst_drive_qubit);
    j.set("min_T1_coupler_us", fin(B.min_T1_coupler_us)); j.set("ratio_coupler", fin(B.ratio_coupler)); j.set("worst_coupler", B.worst_coupler);
    j.set("max_drive_crosstalk", fin(B.max_drive_crosstalk)); j.set("worst_drive_crosstalk", B.worst_crosstalk); j.set("max_drive_crosstalk_allowed", kMaxDriveCrosstalk);
    j.set("rabi_MHz", kRabiMHz); j.set("drive_power_dBm", B.drive_power_dBm); j.set("max_drive_power_dBm", kMaxDrivePower_dBm);
    j.set("M_measure_pH", B.meas.flux.M_pH); j.set("M_data_pH", B.data.flux.M_pH); j.set("M_coupler_pH", B.coupler_flux.M_pH);
    j.set("M_cm_measure_pH", B.meas.flux.M_cm_pH); j.set("M_cm_data_pH", B.data.flux.M_cm_pH); j.set("M_cm_coupler_pH", B.coupler_flux.M_cm_pH);
    j.set("I_half_flux_quantum_measure_mA", fin(B.meas.flux.I_half_mA)); j.set("I_half_flux_quantum_data_mA", fin(B.data.flux.I_half_mA));
    j.set("I_half_flux_quantum_coupler_mA", fin(B.coupler_flux.I_half_mA)); j.set("max_I_half_flux_quantum_mA", kMaxHalfFluxQuantum_mA);
    j.set("max_flux_crosstalk", B.max_flux_crosstalk); j.set("worst_flux_crosstalk", B.worst_flux_crosstalk);
    j.set("C_flux_open_measure_fF", B.meas.C_flux_open_fF); j.set("C_flux_open_data_fF", B.data.C_flux_open_fF);
    j.set("min_T1_flux_capacitive_us", fin(B.min_T1_flux_cap_us)); j.set("min_ratio_flux_capacitive", fin(B.min_ratio_flux_cap)); j.set("worst_flux_capacitive_qubit", B.worst_flux_cap_qubit);
    j.set("assumed_squid_asymmetry", kAssumedSquidAsymmetry);
    j.set("min_T1_flux_loop_us", fin(B.min_T1_flux_loop_us)); j.set("min_ratio_flux_loop", fin(B.min_ratio_flux_loop)); j.set("worst_flux_loop_qubit", B.worst_flux_loop_qubit);
    j.set("min_T1_flux_common_mode_us", fin(B.min_T1_flux_cm_us)); j.set("min_ratio_flux_common_mode", fin(B.min_ratio_flux_cm)); j.set("worst_flux_common_mode_qubit", B.worst_flux_cm_qubit);
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
    for (const auto& r : L.ground_bumps) s.polys.push_back(poly(r, 6));
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
        o.set("corner_relief_gap_um", q.relief_gap); o.set("squid_width_um", q.squid_w); o.set("squid_gap_um", q.squid_gap);
        o.set("flux_slot_up_um", q.slot_up); o.set("flux_slot_down_um", q.slot_down); o.set("flux_slot_feed_um", q.slot_feed);
        return o;
    };
    j.set("pitch_um", c.pitch); j.set("control_lines", c.separate_xy ? "separate" : "combined"); if (c.separate_xy) j.set("xy_tip_um", c.xy_tip);
    j.set("measure_qubit", xm(c.meas)); j.set("data_qubit", xm(c.data));
    Json g = Json::object();
    g.set("tbar_gap_um", c.cg.g_c); g.set("tbar_thickness_um", c.cg.t_w); g.set("offset_um", c.cg.y1); g.set("coupler_tbar_um", c.cg.l_c);
    g.set("body_width_um", c.cg.w_body); g.set("gap_um", c.cg.gap); g.set("pad_um", c.cg.pad); g.set("bus_tbar_um", c.cg.l_b);
    g.set("bus_width_um", c.cg.w_bus); g.set("bus_gap_um", c.cg.gap_bus); g.set("bus_end_gap_um", c.cg.bus_gap());
    g.set("flux_slot_short_side_um", c.cg.slot_lo); g.set("flux_slot_other_side_um", c.cg.slot_hi);
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
        q.claw_w = o.get("claw_width_um", q.claw_w); q.relief_gap = o.get("corner_relief_gap_um", q.relief_gap); q.squid_w = o.get("squid_width_um", q.squid_w); q.squid_gap = o.get("squid_gap_um", q.squid_gap);
        q.slot_up = o.get("flux_slot_up_um", q.slot_up); q.slot_down = o.get("flux_slot_down_um", q.slot_down); q.slot_feed = o.get("flux_slot_feed_um", q.slot_feed);
    };
    c.pitch = j.get("pitch_um", c.pitch);
    {   // "combined": flux bias and drive on one line per qubit; "separate": an XY line as well
        const std::string mode = j.get("control_lines", "combined");
        if (mode != "combined" && mode != "separate") throw std::runtime_error("cell file: control_lines must be \"combined\" or \"separate\"");
        c.separate_xy = mode == "separate"; c.xy_tip = j.get("xy_tip_um", c.xy_tip);
    }
    xm(j["measure_qubit"], c.meas); xm(j["data_qubit"], c.data);
    const Json& g = j["coupler"];
    c.cg.g_c = g.get("tbar_gap_um", c.cg.g_c); c.cg.t_w = g.get("tbar_thickness_um", c.cg.t_w); c.cg.y1 = g.get("offset_um", c.cg.y1);
    c.cg.l_c = g.get("coupler_tbar_um", c.cg.l_c); c.cg.w_body = g.get("body_width_um", c.cg.w_body); c.cg.gap = g.get("gap_um", c.cg.gap);
    c.cg.pad = g.get("pad_um", c.cg.pad); c.cg.l_b = g.get("bus_tbar_um", c.cg.l_b); c.cg.w_bus = g.get("bus_width_um", c.cg.w_bus);
    c.cg.gap_bus = g.get("bus_gap_um", c.cg.gap_bus); c.cg.g_b = g.get("bus_end_gap_um", c.cg.g_b);
    c.cg.slot_lo = g.get("flux_slot_short_side_um", c.cg.slot_lo); c.cg.slot_hi = g.get("flux_slot_other_side_um", c.cg.slot_hi);
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
