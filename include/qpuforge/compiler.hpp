// qpuforge — compiler: logical circuit -> scheduled native operations on the chip
//
//   lower (CCX, CP, SWAP, CX)  ->  layout (interaction-aware placement)  ->  routing (SWAP insertion
//   along shortest paths of the coupler graph)  ->  native gates {RZ (virtual, zero time), SX, X, CZ}
//   ->  ASAP schedule into time "moments" whose duration is the slowest gate in the moment.
#pragma once
#include "circuit.hpp"
#include <deque>
#include <map>

namespace qpuforge {

struct ChipGraph {
    int n{0};
    std::vector<std::vector<int>> adj;
    std::map<std::pair<int, int>, double> t_cz_ns;   // key (min,max)
    double t_1q_ns{25}, t_ro_ns{800};
    std::vector<std::vector<int>> dist;

    void finalize() {
        dist.assign(n, std::vector<int>(n, 1 << 28));
        for (int s = 0; s < n; ++s) {
            std::deque<int> Q{s}; dist[s][s] = 0;
            while (!Q.empty()) { int u = Q.front(); Q.pop_front(); for (int v : adj[u]) if (dist[s][v] > dist[s][u] + 1) { dist[s][v] = dist[s][u] + 1; Q.push_back(v); } }
        }
    }
    bool adjacent(int a, int b) const { for (int v : adj[a]) if (v == b) return true; return false; }
    std::vector<int> path(int a, int b) const {      // shortest path a -> b
        std::vector<int> p{a};
        while (a != b) {
            int nxt = -1;
            for (int v : adj[a]) if (dist[v][b] == dist[a][b] - 1) { nxt = v; break; }
            if (nxt < 0) throw std::runtime_error("chip graph is disconnected");
            p.push_back(a = nxt);
        }
        return p;
    }
    double cz_time(int a, int b) const {
        auto it = t_cz_ns.find({std::min(a, b), std::max(a, b)});
        if (it == t_cz_ns.end()) throw std::runtime_error("no coupler between physical qubits " + std::to_string(a) + " and " + std::to_string(b));
        return it->second;
    }
};

struct Moment {
    std::vector<Gate> ops;      // virtual RZ frame updates are listed before the physical op on the same qubit
    double duration_ns{0};
};

struct Compiled {
    int n_phys{0}, nclbits{0};
    std::vector<int> layout_initial, layout_final;   // logical -> physical
    std::vector<Moment> moments;
    std::vector<std::pair<int, int>> measures;       // (physical qubit, classical bit)
    std::vector<int> active;                          // physical qubits that must be simulated (sorted)
    int n_swaps{0}, n_cz{0}, n_pulses{0}, n_virtual_rz{0};
    double duration_ns{0};
};

struct CompileOptions {
    std::string layout{"auto"};          // auto | trivial | manual
    std::vector<int> manual_layout;
};

namespace compile_detail {
inline double wrap_angle(double a) { a = std::fmod(a, 2 * PI); if (a > PI) a -= 2 * PI; if (a < -PI) a += 2 * PI; return a; }
// Decompose a 1q unitary into time-ordered native ops (at most two physical pulses)
inline std::vector<Gate> native_1q(int q, const M2& U) {
    std::vector<Gate> out;
    auto rz = [&](double a) { a = wrap_angle(a); if (std::abs(a) > 1e-12) out.push_back({Op::RZ, {q}, gates::RZ(a), a, -1, "rz"}); };
    const Euler E = zyz(U);
    if (std::abs(E.theta) < 1e-10) { rz(E.phi + E.lam); return out; }
    if (std::abs(E.theta - PI / 2) < 1e-10) { rz(E.lam - PI / 2); out.push_back({Op::SX, {q}, gates::SX(), 0, -1, "sx"}); rz(E.phi + PI / 2); return out; }
    if (std::abs(E.theta - PI) < 1e-10) { rz(E.lam - PI / 2); out.push_back({Op::X, {q}, gates::X(), 0, -1, "x"}); rz(E.phi + PI / 2); return out; }
    // U3(theta, phi, lam) = RZ(phi + pi) SX RZ(theta + pi) SX RZ(lam)   (up to global phase)
    rz(E.lam);
    out.push_back({Op::SX, {q}, gates::SX(), 0, -1, "sx"});
    rz(E.theta + PI);
    out.push_back({Op::SX, {q}, gates::SX(), 0, -1, "sx"});
    rz(E.phi + PI);
    return out;
}
inline bool is_identity(const M2& U) { return phase_distance(U, gates::I()) < 1e-12; }

// Lower to {U, CX, CZ}; collect terminal measurements.
inline std::vector<Gate> lower(const Circuit& C, std::vector<std::pair<int, int>>& meas) {
    std::vector<Gate> out;
    std::vector<int> measured(C.n, 0);
    auto U = [&](int q, const M2& m) { out.push_back({Op::U, {q}, m, 0, -1, "u"}); };
    auto CX = [&](int c, int t) { out.push_back({Op::CX, {c, t}, {}, 0, -1, "cx"}); };
    for (const auto& g : C.gates) {
        for (int q : g.q) if (q < 0 || q >= C.n) throw std::runtime_error("gate on qubit out of range");
        if (g.op != Op::MEASURE && g.op != Op::BARRIER)
            for (int q : g.q) if (measured[q]) throw std::runtime_error("mid-circuit measurement is not supported yet (gate after measure on q" + std::to_string(q) + ")");
        switch (g.op) {
        case Op::U: U(g.q[0], g.u); break;
        case Op::CX: CX(g.q[0], g.q[1]); break;
        case Op::CZ: out.push_back({Op::CZ, {g.q[0], g.q[1]}, {}, 0, -1, "cz"}); break;
        case Op::SWAP: CX(g.q[0], g.q[1]); CX(g.q[1], g.q[0]); CX(g.q[0], g.q[1]); break;
        case Op::CP: {
            const int a = g.q[0], b = g.q[1]; const double l = g.theta;
            U(a, gates::P(l / 2)); CX(a, b); U(b, gates::P(-l / 2)); CX(a, b); U(b, gates::P(l / 2));
            break;
        }
        case Op::CCX: {
            const int a = g.q[0], b = g.q[1], c = g.q[2];
            U(c, gates::H()); CX(b, c); U(c, gates::Tdg()); CX(a, c); U(c, gates::T()); CX(b, c); U(c, gates::Tdg()); CX(a, c);
            U(b, gates::T()); U(c, gates::T()); U(c, gates::H()); CX(a, b); U(a, gates::T()); U(b, gates::Tdg()); CX(a, b);
            break;
        }
        case Op::MEASURE: meas.push_back({g.q[0], g.cbit}); measured[g.q[0]] = 1; break;
        case Op::BARRIER: break;
        default: throw std::runtime_error("lower: unexpected native op in logical circuit");
        }
    }
    return out;
}

inline std::vector<int> auto_layout(int nl, const std::vector<Gate>& gs, const ChipGraph& G) {
    std::vector<std::vector<double>> w(nl, std::vector<double>(nl, 0));
    for (const auto& g : gs) if (g.q.size() == 2) { w[g.q[0]][g.q[1]] += 1; w[g.q[1]][g.q[0]] += 1; }
    std::vector<double> tot(nl, 0);
    for (int i = 0; i < nl; ++i) for (int j = 0; j < nl; ++j) tot[i] += w[i][j];
    std::vector<int> l2p(nl, -1), used(G.n, 0);
    // centre of the chip graph: min eccentricity, then max degree
    int centre = 0; int best_ecc = 1 << 30;
    for (int p = 0; p < G.n; ++p) {
        int ecc = 0; for (int q = 0; q < G.n; ++q) ecc = std::max(ecc, G.dist[p][q]);
        if (ecc < best_ecc || (ecc == best_ecc && G.adj[p].size() > G.adj[centre].size())) { best_ecc = ecc; centre = p; }
    }
    std::vector<int> placed;
    for (int step = 0; step < nl; ++step) {
        int l = -1; double bw = -1;
        for (int c = 0; c < nl; ++c) {
            if (l2p[c] >= 0) continue;
            double s = 0; for (int m : placed) s += w[c][m];
            const double key = s * 1e6 + tot[c] - c * 1e-9;
            if (key > bw) { bw = key; l = c; }
        }
        int bp = -1; double bc = 1e300;
        for (int p = 0; p < G.n; ++p) {
            if (used[p]) continue;
            double cost = 0;
            if (placed.empty()) cost = G.dist[centre][p];
            else {
                double ws = 0;
                for (int m : placed) { cost += w[l][m] * G.dist[p][l2p[m]]; ws += w[l][m]; }
                if (ws == 0) { cost = 0; for (int m : placed) cost += G.dist[p][l2p[m]]; }
                else cost += 1e-3 * G.dist[centre][p];
            }
            if (cost < bc) { bc = cost; bp = p; }
        }
        l2p[l] = bp; used[bp] = 1; placed.push_back(l);
    }
    return l2p;
}
} // namespace compile_detail

inline Compiled compile(const Circuit& C, const ChipGraph& G, const CompileOptions& opt = {}) {
    using namespace compile_detail;
    if (C.n > G.n) throw std::runtime_error("circuit needs " + std::to_string(C.n) + " qubits, chip has " + std::to_string(G.n));
    Compiled out; out.n_phys = G.n; out.nclbits = C.nclbits;
    std::vector<std::pair<int, int>> meas_logical;
    const auto lowered = lower(C, meas_logical);

    // ---- layout
    std::vector<int> l2p(C.n);
    if (opt.layout == "manual") {
        if (int(opt.manual_layout.size()) != C.n) throw std::runtime_error("manual layout size mismatch");
        l2p = opt.manual_layout;
    } else if (opt.layout == "trivial") {
        std::iota(l2p.begin(), l2p.end(), 0);
    } else {
        l2p = auto_layout(C.n, lowered, G);
    }
    std::vector<int> p2l(G.n, -1);
    for (int l = 0; l < C.n; ++l) {
        if (l2p[l] < 0 || l2p[l] >= G.n || p2l[l2p[l]] >= 0) throw std::runtime_error("invalid layout");
        p2l[l2p[l]] = l;
    }
    out.layout_initial = l2p;

    // ---- routing (physical gate list in {U, CX, CZ})
    std::vector<Gate> phys;
    std::vector<int> touched(G.n, 0);
    for (int l = 0; l < C.n; ++l) touched[l2p[l]] = 1;
    auto do_swap = [&](int pa, int pb) {
        phys.push_back({Op::CX, {pa, pb}, {}, 0, -1, "cx"});
        phys.push_back({Op::CX, {pb, pa}, {}, 0, -1, "cx"});
        phys.push_back({Op::CX, {pa, pb}, {}, 0, -1, "cx"});
        const int la = p2l[pa], lb = p2l[pb];
        p2l[pa] = lb; p2l[pb] = la;
        if (la >= 0) l2p[la] = pb;
        if (lb >= 0) l2p[lb] = pa;
        touched[pa] = touched[pb] = 1;
        ++out.n_swaps;
    };
    for (const auto& g : lowered) {
        if (g.q.size() == 1) { phys.push_back({Op::U, {l2p[g.q[0]]}, g.u, 0, -1, "u"}); continue; }
        const int la = g.q[0], lb = g.q[1];
        while (G.dist[l2p[la]][l2p[lb]] > 1) {
            const auto p = G.path(l2p[la], l2p[lb]);
            do_swap(p[0], p[1]);
        }
        phys.push_back({g.op, {l2p[la], l2p[lb]}, {}, 0, -1, g.name});
    }
    out.layout_final = l2p;
    for (auto [lq, cb] : meas_logical) out.measures.push_back({l2p[lq], cb});

    // ---- native lowering with 1q merging
    std::vector<Gate> native;
    std::vector<M2> pending(G.n, gates::I());
    std::vector<int> has_pending(G.n, 0);
    auto flush = [&](int q) {
        if (!has_pending[q]) return;
        if (!is_identity(pending[q])) for (auto& ng : native_1q(q, pending[q])) native.push_back(ng);
        pending[q] = gates::I(); has_pending[q] = 0;
    };
    auto add1 = [&](int q, const M2& m) { pending[q] = m * pending[q]; has_pending[q] = 1; };
    for (const auto& g : phys) {
        if (g.op == Op::U) { add1(g.q[0], g.u); continue; }
        int a = g.q[0], b = g.q[1];
        if (g.op == Op::CX) {           // CX(c,t) = H_t CZ H_t
            add1(b, gates::H());
            flush(a); flush(b);
            native.push_back({Op::CZ, {a, b}, {}, 0, -1, "cz"});
            add1(b, gates::H());
        } else {
            flush(a); flush(b);
            native.push_back({Op::CZ, {a, b}, {}, 0, -1, "cz"});
        }
    }
    for (int q = 0; q < G.n; ++q) flush(q);

    // ---- ASAP schedule into moments; virtual RZ rides along with the next physical op on that qubit
    std::vector<int> next_free(G.n, 0);
    std::vector<std::vector<Gate>> pending_rz(G.n);
    auto ensure = [&](int m) { while (int(out.moments.size()) <= m) out.moments.emplace_back(); };
    for (const auto& g : native) {
        if (g.op == Op::RZ) { pending_rz[g.q[0]].push_back(g); ++out.n_virtual_rz; continue; }
        int m = 0;
        for (int q : g.q) m = std::max(m, next_free[q]);
        ensure(m);
        for (int q : g.q) {
            for (auto& rz : pending_rz[q]) out.moments[m].ops.push_back(rz);
            pending_rz[q].clear();
        }
        out.moments[m].ops.push_back(g);
        double dur = (g.op == Op::CZ) ? G.cz_time(g.q[0], g.q[1]) : G.t_1q_ns;
        out.moments[m].duration_ns = std::max(out.moments[m].duration_ns, dur);
        for (int q : g.q) next_free[q] = m + 1;
        if (g.op == Op::CZ) ++out.n_cz; else ++out.n_pulses;
    }
    bool tail = false;
    for (int q = 0; q < G.n; ++q) if (!pending_rz[q].empty()) tail = true;
    if (tail) {
        out.moments.emplace_back();
        for (int q = 0; q < G.n; ++q) { for (auto& rz : pending_rz[q]) out.moments.back().ops.push_back(rz); pending_rz[q].clear(); }
    }
    for (const auto& m : out.moments) out.duration_ns += m.duration_ns;
    for (int q = 0; q < G.n; ++q) if (touched[q]) out.active.push_back(q);
    return out;
}

} // namespace qpuforge
