// qpuforge — quantum error correction on an emulated chip
//
//   memory_z_experiment : rotated surface-code memory (prepare |0>_L, R rounds of stabiliser measurement with
//                         mid-circuit measurement + reset of the measure qubits, final data readout)
//   FrameProgram        : Pauli-frame model of the compiled (Clifford) schedule with physics-derived Pauli noise
//                         (Pauli-twirled T1/T_phi, control errors, ZZ crosstalk, readout flips, reset errors)
//   enumerate_faults    : every single fault location -> which detectors / the logical it flips
//   build_dem           : detector error model (graph with fault-probability weights) from those signatures
//   MatchingDecoder     : exact minimum-weight perfect matching (bitmask DP over defects, all-pairs paths)
//   run_memory_frame    : fast Monte Carlo of logical error rate (Pauli frames)
//   run_memory_statevector : the same experiment with the exact non-Pauli noise (amplitude damping, coherent ZZ),
//                         one statevector trajectory per shot — used to validate the Pauli-frame approximation
#pragma once
#include "emulator.hpp"
#include "surface.hpp"
#include <unordered_map>

namespace qpuforge {

// ================================================================== memory experiment
struct MemoryExperiment {
    Circuit circuit;
    std::vector<std::vector<int>> detectors;   // decoded (Z-type) detectors: XOR of these classical bits
    std::vector<int> observable;               // logical Z = XOR of these classical bits
    std::vector<std::string> detector_label;
    std::vector<int> reference;                // noiseless value of each classical bit: 0/1 if deterministic, -1 if random
    int rounds{0}, d{0};
};

inline MemoryExperiment memory_z_experiment(const SurfaceCode& S, int rounds) {
    if (rounds < 1) throw std::runtime_error("memory experiment needs >= 1 round");
    MemoryExperiment E; E.rounds = rounds; E.d = S.d;
    Circuit& C = E.circuit;
    C.n = S.n_qubits();
    C.name = "memory-Z d=" + std::to_string(S.d) + " R=" + std::to_string(rounds);
    for (int q = 0; q < C.n; ++q) C.reset(q);
    std::vector<std::vector<int>> mb(S.stabs.size());
    int cb = 0;
    for (int r = 0; r < rounds; ++r) {
        for (const auto& st : S.stabs) C.h(st.qubit);
        for (int step = 0; step < 4; ++step) {
            for (const auto& st : S.stabs) if (st.type == 'X' && st.order[step] >= 0) C.h(S.data[std::size_t(st.order[step])].qubit);
            for (const auto& st : S.stabs) if (st.order[step] >= 0) C.cz(st.qubit, S.data[std::size_t(st.order[step])].qubit);
            for (const auto& st : S.stabs) if (st.type == 'X' && st.order[step] >= 0) C.h(S.data[std::size_t(st.order[step])].qubit);
        }
        for (const auto& st : S.stabs) C.h(st.qubit);
        for (std::size_t s = 0; s < S.stabs.size(); ++s) { C.measure(S.stabs[s].qubit, cb); mb[s].push_back(cb++); }
        if (r + 1 < rounds) for (const auto& st : S.stabs) C.reset(st.qubit);
    }
    std::vector<int> db(S.data.size());
    for (std::size_t k = 0; k < S.data.size(); ++k) { C.measure(S.data[k].qubit, cb); db[k] = cb++; }
    // Z-stabiliser outcomes are deterministic (+1 -> bit 0) for |0>_L; X-stabiliser and data outcomes are random
    E.reference.assign(std::size_t(cb), -1);
    for (std::size_t s = 0; s < S.stabs.size(); ++s) if (S.stabs[s].type == 'Z') for (int c : mb[s]) E.reference[std::size_t(c)] = 0;
    for (std::size_t s = 0; s < S.stabs.size(); ++s) {
        if (S.stabs[s].type != 'Z') continue;
        const std::string nm = "Z@(" + std::to_string(S.stabs[s].x) + "," + std::to_string(S.stabs[s].y) + ")";
        E.detectors.push_back({mb[s][0]}); E.detector_label.push_back(nm + " r0");
        for (int r = 1; r < rounds; ++r) { E.detectors.push_back({mb[s][std::size_t(r)], mb[s][std::size_t(r - 1)]}); E.detector_label.push_back(nm + " r" + std::to_string(r)); }
        std::vector<int> fin{mb[s][std::size_t(rounds - 1)]};
        for (int k : S.stabs[s].support) fin.push_back(db[std::size_t(k)]);
        E.detectors.push_back(fin); E.detector_label.push_back(nm + " final");
    }
    for (int k : S.logical_z) E.observable.push_back(db[std::size_t(k)]);
    return E;
}

// ================================================================== fast RNG for frame sampling
struct Xoshiro256 {
    std::uint64_t s[4];
    explicit Xoshiro256(std::uint64_t seed) {
        for (auto& v : s) { seed += 0x9E3779B97F4A7C15ULL; std::uint64_t z = seed; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; v = z ^ (z >> 31); }
    }
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    std::uint64_t next() {
        const std::uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        return r;
    }
    double uniform() { return double(next() >> 11) * 0x1.0p-53; }
};

// ================================================================== Pauli-frame program
struct FrameProgram {
    // p: readout flip / reset error; ref: noiseless outcome (-1 random); mover/pleak: CZ leakage injection
    struct FOp { Op op; int a{-1}, b{-1}; int quarter{0}; int cbit{-1}; double p{0}; int ref{-1}; int mover{-1}; double pleak{0}; };
    struct Ctrl2 { int a, b; double p; };
    struct ZZTwirl { int a, b; double pa, pb, pab; };                                    // Z_a, Z_b, Z_aZ_b probabilities
    struct Pauli1 { int q; double px, py, pz; };
    struct FMoment {
        std::vector<FOp> ops;
        std::vector<int> ctrl1;
        std::vector<Ctrl2> ctrl2;
        std::vector<int> dqlr;            // data qubits whose leakage-removal pulse follows this slice
        std::vector<ZZTwirl> zz;
        std::vector<Pauli1> deco;
        double duration_ns{0};
    };
    int n{0}, nclbits{0};
    double p1{0}, p2{0};
    std::vector<FMoment> moments;
    int ignored_exchange_pairs{0};
    double duration_ns{0};
    // leakage model
    bool leakage{false}, reset_removes_leakage{false};
    double leak_1q{0}, transport{0}, dqlr_eff{0}, dqlr_err{0};
    std::vector<double> seep_rate;        // per sim qubit, 1/ns (|2> -> |1> decay, ~2/T1 for a transmon)
};

// reference[cbit] (optional): noiseless outcome if deterministic, used to pick the asymmetric readout flip
// probability P(1|0) or P(0|1); random outcomes use the average.
inline FrameProgram build_frame_program(const Compiled& prog, const NoiseModel& nm, bool noisy = true,
                                        const std::vector<int>& reference = {},
                                        const std::map<int, std::vector<int>>& dqlr_partners = {}) {
    // Noise is converted to Pauli channels with two refinements that matter for QEC circuits:
    //  * state-aware decoherence: while a qubit is in a known computational basis state of the noiseless circuit
    //    (after RESET, before its next non-diagonal gate), amplitude damping is exactly a classical bit flip
    //    (decay from |1>, thermal excitation from |0>) and dephasing does nothing;
    //  * coherent ZZ windows: conditional phases commute with CZ and RZ, so each pair's phase is accumulated until a
    //    non-commuting op (SX, X, MEASURE, RESET) on either qubit and only then Pauli-twirled. A spectator known to be
    //    in |0> makes the phase vanish; one known in |1> turns it into a Z rotation of its partner.
    FrameProgram F;
    const auto sim = sim_index(prog);
    F.n = int(prog.active.size());
    F.nclbits = prog.nclbits;
    F.p1 = (noisy && nm.control) ? 1.5 * nm.err_1q_ctrl : 0.0;
    F.p2 = (noisy && nm.control) ? 1.25 * nm.err_2q_ctrl : 0.0;
    F.leakage = noisy && nm.leakage;
    if (F.leakage) {
        F.leak_1q = nm.leak_1q; F.transport = nm.leak_transport; F.reset_removes_leakage = nm.reset_removes_leakage;
        F.dqlr_eff = nm.dqlr ? nm.dqlr_eff : 0.0; F.dqlr_err = nm.dqlr ? nm.dqlr_err : 0.0;
        F.seep_rate.assign(std::size_t(F.n), 0.0);
        for (int k = 0; k < F.n; ++k) {
            const double T1 = nm.T1_ns[std::size_t(prog.active[std::size_t(k)])];
            F.seep_rate[std::size_t(k)] = std::isfinite(T1) ? 2.0 / T1 : 0.0;
        }
    }
    struct ZZPair { int a, b; double rate, phi; };
    std::vector<ZZPair> zz;
    if (noisy && nm.zz)
        for (const auto& [pr, rate] : nm.zz_GHz) {
            const int a = sim[pr.first], b = sim[pr.second];
            if (a >= 0 && b >= 0 && rate != 0) zz.push_back({std::min(a, b), std::max(a, b), rate, 0.0});
        }
    if (noisy && nm.exchange_on)
        for (const auto& [pr, v] : nm.exchange) if (sim[pr.first] >= 0 && sim[pr.second] >= 0) ++F.ignored_exchange_pairs;
    std::vector<int> kn(std::size_t(F.n), 0);   // 0: |0>, 1: |1>, 2: superposition / unknown
    auto twirl_zz = [&](const ZZPair& z, std::vector<FrameProgram::ZZTwirl>& out) {
        if (z.phi == 0) return;
        const double s2 = std::pow(std::sin(z.phi / 2), 2);
        const int ka = kn[std::size_t(z.a)], kb = kn[std::size_t(z.b)];
        if (ka == 0 || kb == 0) return;                                        // conditional phase is the identity
        if (ka == 1 && kb == 1) return;                                        // global phase
        if (ka == 1) { out.push_back({z.a, z.b, 0.0, s2, 0.0}); return; }      // Z rotation of b
        if (kb == 1) { out.push_back({z.a, z.b, s2, 0.0, 0.0}); return; }      // Z rotation of a
        out.push_back({z.a, z.b, s2 / 4, s2 / 4, s2 / 4});
    };
    auto non_diagonal = [](Op op) { return op == Op::SX || op == Op::X || op == Op::MEASURE || op == Op::RESET; };
    for (const auto& m : prog.moments) {
        FrameProgram::FMoment fm;
        fm.duration_ns = m.duration_ns;
        F.duration_ns += m.duration_ns;
        // close ZZ windows of pairs touched by a non-commuting op at the start of this slice
        std::vector<char> touched(std::size_t(F.n), 0);
        for (const auto& g : m.ops) if (non_diagonal(g.op)) touched[std::size_t(sim[g.q[0]])] = 1;
        for (auto& z : zz)
            if (touched[std::size_t(z.a)] || touched[std::size_t(z.b)]) {
                if (!F.moments.empty()) twirl_zz(z, F.moments.back().zz);
                z.phi = 0;
            }
        std::vector<std::pair<int, int>> busy;
        for (auto [a, b] : m.busy) busy.push_back({std::min(sim[a], sim[b]), std::max(sim[a], sim[b])});
        for (const auto& g : m.ops) {
            FrameProgram::FOp o; o.op = g.op; o.a = sim[g.q[0]];
            switch (g.op) {
            case Op::RZ: {
                const double qt = g.theta / (PI / 2);
                const long k = std::lround(qt);
                if (std::abs(qt - double(k)) > 1e-6) throw std::runtime_error("Pauli-frame simulation needs a Clifford circuit (RZ angle not a multiple of pi/2)");
                o.quarter = int(((k % 4) + 4) % 4);
                break;
            }
            case Op::SX: if (noisy) fm.ctrl1.push_back(o.a); kn[std::size_t(o.a)] = 2; break;
            case Op::X: if (noisy) fm.ctrl1.push_back(o.a); if (kn[std::size_t(o.a)] < 2) kn[std::size_t(o.a)] ^= 1; break;
            case Op::CZ:
                o.b = sim[g.q[1]];
                if (noisy) fm.ctrl2.push_back({o.a, o.b, nm.control ? 1.25 * nm.err_2q(g.q[0], g.q[1]) : 0.0});
                if (F.leakage) {
                    auto it = nm.leak_cz.find({std::min(g.q[0], g.q[1]), std::max(g.q[0], g.q[1])});
                    if (it != nm.leak_cz.end()) { o.mover = sim[it->second.first]; o.pleak = it->second.second; }
                }
                if (std::find(busy.begin(), busy.end(), std::make_pair(std::min(o.a, o.b), std::max(o.a, o.b))) == busy.end())
                    busy.push_back({std::min(o.a, o.b), std::max(o.a, o.b)});
                break;
            case Op::MEASURE:
                o.cbit = g.cbit;
                if (noisy && nm.readout) {
                    const int ref = (g.cbit >= 0 && std::size_t(g.cbit) < reference.size()) ? reference[std::size_t(g.cbit)] : -1;
                    const double p01 = nm.ro01[std::size_t(g.q[0])], p10 = nm.ro10[std::size_t(g.q[0])];
                    o.p = ref == 0 ? p01 : (ref == 1 ? p10 : 0.5 * (p01 + p10));
                }
                o.ref = (g.cbit >= 0 && std::size_t(g.cbit) < reference.size()) ? reference[std::size_t(g.cbit)] : -1;
                break;
            case Op::RESET: {
                o.p = (noisy && nm.reset) ? nm.reset_err : 0.0; kn[std::size_t(o.a)] = 0;
                auto it = dqlr_partners.find(g.q[0]);
                if (F.leakage && nm.dqlr && it != dqlr_partners.end()) for (int dq : it->second) if (sim[dq] >= 0) fm.dqlr.push_back(sim[dq]);
                break;
            }
            default: throw std::runtime_error("frame program: unexpected op");
            }
            fm.ops.push_back(o);
        }
        if (noisy && m.duration_ns > 0) {
            const double tau = m.duration_ns;
            for (auto& z : zz) {
                if (std::find(busy.begin(), busy.end(), std::make_pair(z.a, z.b)) != busy.end()) continue;
                z.phi += 2 * PI * z.rate * tau;
            }
            for (int k = 0; k < F.n; ++k) {
                const int p = prog.active[std::size_t(k)];
                const bool hasT1 = nm.T1 && std::isfinite(nm.T1_ns[std::size_t(p)]);
                const double gamma = hasT1 ? 1 - std::exp(-tau / nm.T1_ns[std::size_t(p)]) : 0.0;
                const double pdown = nm.thermal ? 1 - nm.p_thermal[std::size_t(p)] : 1.0;
                if (kn[std::size_t(k)] == 0) { if (gamma * (1 - pdown) > 0) fm.deco.push_back({k, gamma * (1 - pdown), 0.0, 0.0}); continue; }
                if (kn[std::size_t(k)] == 1) { if (gamma * pdown > 0) fm.deco.push_back({k, gamma * pdown, 0.0, 0.0}); continue; }
                double rate2 = hasT1 ? 0.5 / nm.T1_ns[std::size_t(p)] : 0.0;
                if (nm.Tphi && std::isfinite(nm.Tphi_ns[std::size_t(p)])) rate2 += 1.0 / nm.Tphi_ns[std::size_t(p)];
                const double px = gamma / 4;
                const double pz = std::max(0.0, (1 - std::exp(-tau * rate2)) / 2 - px);
                if (px > 0 || pz > 0) fm.deco.push_back({k, px, px, pz});
            }
        }
        F.moments.push_back(std::move(fm));
    }
    if (!F.moments.empty()) for (auto& z : zz) { twirl_zz(z, F.moments.back().zz); z.phi = 0; }
    return F;
}

namespace frame_detail {
inline void apply_op(const FrameProgram::FOp& o, std::uint8_t* x, std::uint8_t* z, std::uint8_t* flips) {
    switch (o.op) {
    case Op::RZ: if (o.quarter & 1) z[o.a] ^= x[o.a]; break;
    case Op::SX: x[o.a] ^= z[o.a]; break;
    case Op::CZ: { const std::uint8_t xa = x[o.a], xb = x[o.b]; z[o.a] ^= xb; z[o.b] ^= xa; break; }
    case Op::MEASURE: flips[o.cbit] ^= x[o.a]; z[o.a] = 0; break;
    case Op::RESET: x[o.a] = 0; z[o.a] = 0; break;
    default: break;   // X: Pauli frame unchanged
    }
}
inline void apply_pauli(int P, int q, std::uint8_t* x, std::uint8_t* z) {   // 1=X 2=Y 3=Z
    if (P == 1 || P == 2) x[q] ^= 1;
    if (P == 2 || P == 3) z[q] ^= 1;
}
} // namespace frame_detail

// ================================================================== single-fault signatures and detector error model
struct FaultSignature {
    std::vector<int> dets;
    bool obs{false};
    double p{0};
    int moment{0};
    std::string what;
};

inline std::vector<FaultSignature> enumerate_faults(const FrameProgram& F, const MemoryExperiment& E) {
    std::vector<FaultSignature> out;
    std::vector<std::uint8_t> x(std::size_t(F.n)), z(std::size_t(F.n)), flips(std::size_t(std::max(F.nclbits, 1)));
    auto evaluate = [&](FaultSignature& s) {
        for (std::size_t d = 0; d < E.detectors.size(); ++d) {
            std::uint8_t v = 0; for (int c : E.detectors[d]) v ^= flips[std::size_t(c)];
            if (v) s.dets.push_back(int(d));
        }
        std::uint8_t o = 0; for (int c : E.observable) o ^= flips[std::size_t(c)];
        s.obs = o;
    };
    auto propagate = [&](std::size_t from_moment) {   // inject already placed in x/z after the ops of from_moment
        for (std::size_t m = from_moment + 1; m < F.moments.size(); ++m)
            for (const auto& o : F.moments[m].ops) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
    };
    auto inject = [&](std::size_t m, const std::vector<std::pair<int, int>>& paulis, double p, const std::string& what) {
        if (p <= 0) return;
        std::fill(x.begin(), x.end(), 0); std::fill(z.begin(), z.end(), 0); std::fill(flips.begin(), flips.end(), 0);
        for (auto [P, q] : paulis) frame_detail::apply_pauli(P, q, x.data(), z.data());
        propagate(m);
        FaultSignature s; s.p = p; s.moment = int(m); s.what = what;
        evaluate(s);
        out.push_back(std::move(s));
    };
    for (std::size_t m = 0; m < F.moments.size(); ++m) {
        const auto& fm = F.moments[m];
        for (const auto& o : fm.ops) {
            if (o.op == Op::MEASURE && o.p > 0) {
                std::fill(flips.begin(), flips.end(), 0);
                flips[std::size_t(o.cbit)] = 1;
                FaultSignature s; s.p = o.p; s.moment = int(m); s.what = "readout flip q" + std::to_string(o.a);
                evaluate(s);
                out.push_back(std::move(s));
            }
            if (o.op == Op::RESET) inject(m, {{1, o.a}}, o.p, "reset error q" + std::to_string(o.a));
        }
        for (int q : fm.ctrl1) for (int P = 1; P <= 3; ++P) inject(m, {{P, q}}, F.p1 / 3, "1q control q" + std::to_string(q));
        for (const auto& c2 : fm.ctrl2)
            for (int k = 1; k < 16; ++k) {
                std::vector<std::pair<int, int>> ps;
                if (k % 4) ps.push_back({k % 4, c2.a});
                if (k / 4) ps.push_back({k / 4, c2.b});
                inject(m, ps, c2.p / 15, "CZ control q" + std::to_string(c2.a) + ",q" + std::to_string(c2.b));
            }
        if (F.dqlr_err > 0) for (int q : fm.dqlr) for (int P = 1; P <= 3; ++P) inject(m, {{P, q}}, F.dqlr_err / 3, "DQLR q" + std::to_string(q));
        for (const auto& t : fm.zz) {
            inject(m, {{3, t.a}}, t.pa, "ZZ");
            inject(m, {{3, t.b}}, t.pb, "ZZ");
            inject(m, {{3, t.a}, {3, t.b}}, t.pab, "ZZ");
        }
        for (const auto& dq : fm.deco) {
            inject(m, {{1, dq.q}}, dq.px, "decoherence q" + std::to_string(dq.q));
            inject(m, {{2, dq.q}}, dq.py, "decoherence q" + std::to_string(dq.q));
            inject(m, {{3, dq.q}}, dq.pz, "decoherence q" + std::to_string(dq.q));
        }
    }
    return out;
}

struct DEMEdge { int a, b; bool obs; double p; };   // b = -1: boundary
struct DetectorErrorModel {
    int n_det{0};
    std::vector<DEMEdge> edges;
    int hyperedge_faults{0};          // faults flipping > 2 decoded detectors (dropped from the graph)
    int undetectable_logical{0};      // faults flipping the logical but no detector (distance-killing; must be 0)
};

inline DetectorErrorModel build_dem(const std::vector<FaultSignature>& faults, int n_det) {
    DetectorErrorModel M; M.n_det = n_det;
    std::map<std::tuple<int, int, int>, double> acc;
    for (const auto& f : faults) {
        if (f.dets.empty()) { if (f.obs) ++M.undetectable_logical; continue; }
        if (f.dets.size() > 2) { ++M.hyperedge_faults; continue; }
        const int a = f.dets[0], b = f.dets.size() == 2 ? f.dets[1] : -1;
        auto key = std::make_tuple(a, b, int(f.obs));
        double& p = acc[key];
        p = p * (1 - f.p) + f.p * (1 - p);
    }
    for (const auto& [k, p] : acc) M.edges.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k) != 0, p});
    return M;
}

// ================================================================== minimum-weight perfect matching decoder
class MatchingDecoder {
public:
    explicit MatchingDecoder(const DetectorErrorModel& M) : N_(M.n_det), V_(M.n_det + 1) {
        const double INF = 1e300;
        dist_.assign(std::size_t(V_) * V_, INF);
        par_.assign(std::size_t(V_) * V_, 0);
        for (int i = 0; i < V_; ++i) dist_[idx(i, i)] = 0;
        for (const auto& e : M.edges) {
            const double p = std::clamp(e.p, 1e-15, 0.5 - 1e-12);
            const double w = std::log((1 - p) / p);
            const int u = e.a, v = e.b < 0 ? N_ : e.b;
            if (w < dist_[idx(u, v)]) { dist_[idx(u, v)] = dist_[idx(v, u)] = w; par_[idx(u, v)] = par_[idx(v, u)] = e.obs; }
        }
        for (int k = 0; k < V_; ++k) for (int i = 0; i < V_; ++i) {
            const double dik = dist_[idx(i, k)];
            if (dik >= INF) continue;
            for (int j = 0; j < V_; ++j) {
                const double c = dik + dist_[idx(k, j)];
                if (c < dist_[idx(i, j)]) { dist_[idx(i, j)] = c; par_[idx(i, j)] = par_[idx(i, k)] ^ par_[idx(k, j)]; }
            }
        }
    }
    // returns the predicted logical flip for a set of fired detectors
    int decode(const std::vector<int>& def, std::vector<double>& cost, std::vector<int>& choice) const {
        const int k = int(def.size());
        if (k == 0) return 0;
        if (k > 22) return greedy(def);
        const std::size_t full = (std::size_t(1) << k) - 1;
        cost.assign(full + 1, 0.0);
        choice.assign(full + 1, -1);
        for (std::size_t mask = 1; mask <= full; ++mask) {
            int i = 0; while (!((mask >> i) & 1)) ++i;
            const std::size_t rest = mask & ~(std::size_t(1) << i);
            double best = dist_[idx(def[std::size_t(i)], N_)] + cost[rest];
            int ch = -1;
            for (int j = i + 1; j < k; ++j) {
                if (!((rest >> j) & 1)) continue;
                const double c = dist_[idx(def[std::size_t(i)], def[std::size_t(j)])] + cost[rest & ~(std::size_t(1) << j)];
                if (c < best) { best = c; ch = j; }
            }
            cost[mask] = best; choice[mask] = ch;
        }
        int parity = 0;
        std::size_t mask = full;
        while (mask) {
            int i = 0; while (!((mask >> i) & 1)) ++i;
            const int ch = choice[mask];
            if (ch < 0) { parity ^= par_[idx(def[std::size_t(i)], N_)]; mask &= ~(std::size_t(1) << i); }
            else { parity ^= par_[idx(def[std::size_t(i)], def[std::size_t(ch)])]; mask &= ~((std::size_t(1) << i) | (std::size_t(1) << ch)); }
        }
        return parity;
    }

private:
    int N_, V_;
    std::vector<double> dist_;
    std::vector<std::uint8_t> par_;
    std::size_t idx(int i, int j) const { return std::size_t(i) * std::size_t(V_) + std::size_t(j); }
    int greedy(std::vector<int> def) const {
        int parity = 0;
        while (!def.empty()) {
            double best = 1e300; int bi = 0, bj = -1;
            for (std::size_t i = 0; i < def.size(); ++i) {
                if (dist_[idx(def[i], N_)] < best) { best = dist_[idx(def[i], N_)]; bi = int(i); bj = -1; }
                for (std::size_t j = i + 1; j < def.size(); ++j) if (dist_[idx(def[i], def[j])] < best) { best = dist_[idx(def[i], def[j])]; bi = int(i); bj = int(j); }
            }
            if (bj < 0) { parity ^= par_[idx(def[std::size_t(bi)], N_)]; def.erase(def.begin() + bi); }
            else { parity ^= par_[idx(def[std::size_t(bi)], def[std::size_t(bj)])]; def.erase(def.begin() + bj); def.erase(def.begin() + bi); }
        }
        return parity;
    }
};

// ================================================================== experiment runners
struct MemoryResult {
    int d{0}, rounds{0};
    long shots{0}, logical_errors{0};
    double p_logical{0}, p_logical_sem{0};
    double eps_per_round{0};          // from P_L(R) = (1 - (1 - 2 eps)^R) / 2
    double detection_rate{0};         // mean fraction of detectors that fire
    double round_ns{0}, total_ns{0};
    int hyperedge_faults{0}, undetectable_logical{0}, ignored_exchange_pairs{0};
    double leak_frac_data{0}, leak_frac_meas{0};   // time-averaged fraction of leaked qubits (frame model with leakage)
    double leak_events_per_shot{0};                // leakage injections per shot
};

inline double eps_from_PL(double PL, int R) { return PL < 0.5 ? 0.5 * (1 - std::pow(1 - 2 * PL, 1.0 / R)) : 0.5; }

struct MemorySetup {
    SurfaceCode code;
    MemoryExperiment exp;
    Compiled prog;
    FrameProgram frame;
    DetectorErrorModel dem;
    std::unique_ptr<MatchingDecoder> decoder;
    std::vector<char> is_data;                 // by sim index
};

inline MemorySetup prepare_memory(const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, int rounds) {
    if (D.surface_d <= 0) throw std::runtime_error("chip is not a surface-code layout (generate it with topology 'surface')");
    MemorySetup M;
    M.code = rotated_surface_code(D.surface_d);
    M.exp = memory_z_experiment(M.code, rounds);
    CompileOptions co; co.layout = "trivial";
    M.prog = compile(M.exp.circuit, chip_graph(D, R), co);
    if (M.prog.n_swaps) throw std::runtime_error("surface-code circuit needed SWAPs: couplers do not match the code layout");
    // data-qubit leakage removal: each data qubit is paired with a neighbouring measure qubit (least-loaded first)
    std::map<int, std::vector<int>> partners;
    if (nm.leakage && nm.dqlr) {
        std::map<int, int> load;
        for (std::size_t k = 0; k < M.code.data.size(); ++k) {
            int best = -1;
            for (const auto& st : M.code.stabs)
                if (std::find(st.support.begin(), st.support.end(), int(k)) != st.support.end())
                    if (best < 0 || load[st.qubit] < load[best]) best = st.qubit;
            if (best >= 0) { ++load[best]; partners[D.index_of(best)].push_back(D.index_of(M.code.data[k].qubit)); }
        }
    }
    M.frame = build_frame_program(M.prog, nm, true, M.exp.reference, partners);
    {
        const auto sim = sim_index(M.prog);
        M.is_data.assign(std::size_t(M.frame.n), 0);
        for (const auto& dq : M.code.data) { const int s = sim[std::size_t(D.index_of(dq.qubit))]; if (s >= 0) M.is_data[std::size_t(s)] = 1; }
    }
    M.dem = build_dem(enumerate_faults(M.frame, M.exp), int(M.exp.detectors.size()));
    M.decoder = std::make_unique<MatchingDecoder>(M.dem);
    return M;
}

inline MemoryResult summarize(const MemorySetup& M, long shots, long errors, double fired_sum) {
    MemoryResult r;
    r.d = M.code.d; r.rounds = M.exp.rounds; r.shots = shots; r.logical_errors = errors;
    r.p_logical = double(errors) / double(shots);
    r.p_logical_sem = std::sqrt(std::max(r.p_logical * (1 - r.p_logical), 1e-12) / double(shots));
    r.eps_per_round = eps_from_PL(r.p_logical, r.rounds);
    r.detection_rate = fired_sum / double(shots) / double(M.exp.detectors.size());
    r.total_ns = M.prog.duration_ns;
    r.round_ns = M.prog.duration_ns / double(r.rounds);
    r.hyperedge_faults = M.dem.hyperedge_faults;
    r.undetectable_logical = M.dem.undetectable_logical;
    r.ignored_exchange_pairs = M.frame.ignored_exchange_pairs;
    return r;
}

// Fast Monte Carlo with Pauli frames (+ leakage: a per-qubit flag outside the frame).
//  * CZ leakage: with the pair's probability the moving qubit leaves the subspace (|11> -> |20>); its partner is
//    projected, modelled as random X and Z.
//  * a leaked qubit ignores single-qubit gates; a CZ with a leaked partner applies an uncontrolled conditional phase
//    (random Z on the partner) and, with probability `transport`, moves the leakage to the partner;
//  * a leaked qubit reads out as |1>; RESET removes leakage only if the reset is leakage-removing;
//  * seepage |2> -> |1> at ~2/T1 returns the qubit in a random computational state; DQLR does so on demand.
inline MemoryResult run_memory_frame(const MemorySetup& M, long shots, std::uint64_t seed) {
    const auto& F = M.frame;
    long errors = 0; double fired_sum = 0, leak_d = 0, leak_m = 0, leak_ev = 0;
    int n_data = 0; for (char c : M.is_data) n_data += c;
    const int n_meas = F.n - n_data;
#pragma omp parallel reduction(+ : errors, fired_sum, leak_d, leak_m, leak_ev)
    {
        int tid = 0, nth = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nth = omp_get_num_threads();
#endif
        Xoshiro256 rng(seed * 0x9E3779B97F4A7C15ULL + std::uint64_t(tid) * 7919ULL + 1);
        std::vector<std::uint8_t> x(std::size_t(F.n)), z(std::size_t(F.n)), flips(std::size_t(std::max(F.nclbits, 1))), L(std::size_t(F.n));
        std::vector<int> fired; std::vector<double> cost; std::vector<int> choice;
        auto randomize = [&](int q) { x[std::size_t(q)] = rng.uniform() < 0.5; z[std::size_t(q)] = rng.uniform() < 0.5; };
        for (long s = tid; s < shots; s += nth) {
            std::fill(x.begin(), x.end(), 0); std::fill(z.begin(), z.end(), 0); std::fill(flips.begin(), flips.end(), 0);
            std::fill(L.begin(), L.end(), 0);
            int nleaked = 0;
            for (const auto& fm : F.moments) {
                for (const auto& o : fm.ops) {
                    if (nleaked == 0 || !F.leakage) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                    else {
                        const bool la = L[std::size_t(o.a)] != 0;
                        if (o.op == Op::CZ) {
                            const bool lb = L[std::size_t(o.b)] != 0;
                            if (!la && !lb) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                            else if (la != lb) {
                                const int leaked = la ? o.a : o.b, other = la ? o.b : o.a;
                                if (rng.uniform() < 0.5) z[std::size_t(other)] ^= 1;
                                if (F.transport > 0 && rng.uniform() < F.transport) { L[std::size_t(leaked)] = 0; randomize(leaked); L[std::size_t(other)] = 1; }
                            }
                        } else if (o.op == Op::MEASURE) {
                            if (la) flips[std::size_t(o.cbit)] = o.ref >= 0 ? std::uint8_t(1 ^ o.ref) : std::uint8_t(rng.uniform() < 0.5);
                            else frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                        } else if (o.op == Op::RESET) {
                            if (!la) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                            else if (F.reset_removes_leakage) { L[std::size_t(o.a)] = 0; --nleaked; frame_detail::apply_op(o, x.data(), z.data(), flips.data()); }
                        } else if (!la) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                    }
                    if (o.p > 0 && rng.uniform() < o.p) {
                        if (o.op == Op::MEASURE) flips[std::size_t(o.cbit)] ^= 1;
                        else if (o.op == Op::RESET) x[std::size_t(o.a)] ^= 1;
                    }
                    if (F.leakage) {
                        if (o.op == Op::CZ && o.pleak > 0 && !L[std::size_t(o.mover)] && rng.uniform() < o.pleak) {
                            const int other = o.mover == o.a ? o.b : o.a;
                            L[std::size_t(o.mover)] = 1; ++nleaked; leak_ev += 1;
                            if (!L[std::size_t(other)]) randomize(other);
                        } else if ((o.op == Op::SX || o.op == Op::X) && F.leak_1q > 0 && !L[std::size_t(o.a)] && rng.uniform() < F.leak_1q) {
                            L[std::size_t(o.a)] = 1; ++nleaked; leak_ev += 1;
                        }
                    }
                }
                if (F.p1 > 0) for (int q : fm.ctrl1) if (rng.uniform() < F.p1) frame_detail::apply_pauli(1 + int(rng.uniform() * 3), q, x.data(), z.data());
                for (const auto& c2 : fm.ctrl2) if (c2.p > 0 && rng.uniform() < c2.p) {
                    const int k = 1 + int(rng.uniform() * 15);
                    if (k % 4) frame_detail::apply_pauli(k % 4, c2.a, x.data(), z.data());
                    if (k / 4) frame_detail::apply_pauli(k / 4, c2.b, x.data(), z.data());
                }
                for (const auto& t : fm.zz) {
                    const double r = rng.uniform();
                    if (r < t.pa) z[std::size_t(t.a)] ^= 1;
                    else if (r < t.pa + t.pb) z[std::size_t(t.b)] ^= 1;
                    else if (r < t.pa + t.pb + t.pab) { z[std::size_t(t.a)] ^= 1; z[std::size_t(t.b)] ^= 1; }
                }
                for (const auto& dq : fm.deco) {
                    const double r = rng.uniform();
                    if (r < dq.px) x[std::size_t(dq.q)] ^= 1;
                    else if (r < dq.px + dq.py) { x[std::size_t(dq.q)] ^= 1; z[std::size_t(dq.q)] ^= 1; }
                    else if (r < dq.px + dq.py + dq.pz) z[std::size_t(dq.q)] ^= 1;
                }
                if (F.leakage) {
                    for (int q : fm.dqlr) {
                        if (L[std::size_t(q)] && rng.uniform() < F.dqlr_eff) { L[std::size_t(q)] = 0; --nleaked; randomize(q); }
                        if (F.dqlr_err > 0 && rng.uniform() < F.dqlr_err) frame_detail::apply_pauli(1 + int(rng.uniform() * 3), q, x.data(), z.data());
                    }
                    if (nleaked > 0) {
                        int nd = 0, nmq = 0;
                        for (int q = 0; q < F.n; ++q) if (L[std::size_t(q)]) { if (M.is_data[std::size_t(q)]) ++nd; else ++nmq; }
                        leak_d += fm.duration_ns * nd; leak_m += fm.duration_ns * nmq;
                        for (int q = 0; q < F.n; ++q)
                            if (L[std::size_t(q)] && rng.uniform() < 1 - std::exp(-fm.duration_ns * F.seep_rate[std::size_t(q)])) {
                                L[std::size_t(q)] = 0; --nleaked; randomize(q);
                            }
                    }
                }
            }
            fired.clear();
            for (std::size_t d = 0; d < M.exp.detectors.size(); ++d) {
                std::uint8_t v = 0; for (int c : M.exp.detectors[d]) v ^= flips[std::size_t(c)];
                if (v) fired.push_back(int(d));
            }
            std::uint8_t obs = 0; for (int c : M.exp.observable) obs ^= flips[std::size_t(c)];
            fired_sum += double(fired.size());
            if (M.decoder->decode(fired, cost, choice) != obs) ++errors;
        }
    }
    auto r = summarize(M, shots, errors, fired_sum);
    const double T = std::max(M.frame.duration_ns, 1e-12) * double(shots);
    r.leak_frac_data = n_data ? leak_d / T / n_data : 0.0;
    r.leak_frac_meas = n_meas ? leak_m / T / n_meas : 0.0;
    r.leak_events_per_shot = leak_ev / double(shots);
    return r;
}

// Exact-noise cross-check: statevector trajectories with mid-circuit measurement (one per shot).
inline MemoryResult run_memory_statevector(const MemorySetup& M, const NoiseModel& nm, long shots, unsigned seed, bool noisy = true) {
    const auto recs = run_records(M.prog, nm, int(shots), seed, noisy);
    long errors = 0; double fired_sum = 0;
    std::vector<int> fired; std::vector<double> cost; std::vector<int> choice;
    for (const auto& rc : recs) {
        fired.clear();
        for (std::size_t d = 0; d < M.exp.detectors.size(); ++d) {
            std::uint8_t v = 0; for (int c : M.exp.detectors[d]) v ^= rc[std::size_t(c)];
            if (v) fired.push_back(int(d));
        }
        std::uint8_t obs = 0; for (int c : M.exp.observable) obs ^= rc[std::size_t(c)];
        fired_sum += double(fired.size());
        if (M.decoder->decode(fired, cost, choice) != obs) ++errors;
    }
    return summarize(M, shots, errors, fired_sum);
}

// Logical error per round from several round counts: 1 - 2 P_L(R) = A (1 - 2 eps)^R  (least squares in log space)
inline double fit_eps_per_round(const std::vector<MemoryResult>& rs) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
    for (const auto& r : rs) {
        if (r.p_logical >= 0.5 || r.p_logical <= 0) continue;
        const double x = r.rounds, y = std::log(1 - 2 * r.p_logical);
        sx += x; sy += y; sxx += x * x; sxy += x * y; ++n;
    }
    if (n < 2) return rs.empty() ? 0.0 : rs.back().eps_per_round;
    const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    return 0.5 * (1 - std::exp(slope));
}

} // namespace qpuforge
