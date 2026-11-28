// qpuforge — quantum error correction on an emulated chip
//
//   memory_experiment   : rotated surface-code memory in the Z or the X basis (prepare |0>_L or |+>_L, R rounds of
//                         stabiliser measurement with mid-circuit measurement + reset of the measure qubits, final
//                         data readout in the same basis)
//   FrameProgram        : Pauli-frame model of the compiled (Clifford) schedule with physics-derived Pauli noise
//                         (Pauli-twirled T1/T_phi, control errors, ZZ crosstalk, readout flips, reset errors)
//   enumerate_faults    : every single fault location -> which detectors / the logical it flips
//   build_dem           : detector error model (graph with fault-probability weights) from those signatures
//   MatchingDecoder     : exact minimum-weight perfect matching for any number of defects (bitmask DP for a few
//                         defects, Edmonds' blossom algorithm beyond; all-pairs shortest paths in the detector graph)
//   run_memory_frame    : fast Monte Carlo of logical error rate (Pauli frames)
//   run_memory_statevector : the same experiment with the exact non-Pauli noise (amplitude damping, coherent ZZ),
//                         one statevector trajectory per shot — used to validate the Pauli-frame approximation
#pragma once
#include "emulator.hpp"
#include "matching.hpp"
#include "surface.hpp"
#include "unionfind.hpp"
#include <chrono>
#include <cstdio>
#include <limits>
#include <unordered_map>

namespace qpuforge {

// ================================================================== memory experiment
struct MemoryExperiment {
    Circuit circuit;
    std::vector<std::vector<int>> detectors;   // decoded detectors (stabilisers of the memory basis): XOR of these bits
    std::vector<int> observable;               // the logical operator of the memory basis = XOR of these classical bits
    std::vector<std::string> detector_label;
    std::vector<int> detector_round;           // syndrome round of each detector (rounds = the final data readout)
    std::vector<int> reference;                // noiseless value of each classical bit: 0/1 if deterministic, -1 if random
    int rounds{0}, d{0};
    char basis{'Z'};
};

// basis 'Z': |0>_L, Z-type detectors, logical Z (sensitive to bit flips).
// basis 'X': |+>_L, X-type detectors, logical X (sensitive to phase flips).
inline MemoryExperiment memory_experiment(const SurfaceCode& S, int rounds, char basis = 'Z') {
    if (rounds < 1) throw std::runtime_error("memory experiment needs >= 1 round");
    if (basis != 'Z' && basis != 'X') throw std::runtime_error("memory basis must be 'Z' or 'X'");
    MemoryExperiment E; E.rounds = rounds; E.d = S.d; E.basis = basis;
    Circuit& C = E.circuit;
    C.n = S.n_qubits();
    C.name = std::string("memory-") + basis + " d=" + std::to_string(S.d) + " R=" + std::to_string(rounds);
    for (int q = 0; q < C.n; ++q) C.reset(q);
    if (basis == 'X') for (const auto& dq : S.data) C.h(dq.qubit);
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
    if (basis == 'X') for (const auto& dq : S.data) C.h(dq.qubit);
    for (std::size_t k = 0; k < S.data.size(); ++k) { C.measure(S.data[k].qubit, cb); db[k] = cb++; }
    // stabilisers of the memory basis are deterministic (+1 -> bit 0); the other type and the data outcomes are random
    E.reference.assign(std::size_t(cb), -1);
    for (std::size_t s = 0; s < S.stabs.size(); ++s) if (S.stabs[s].type == basis) for (int c : mb[s]) E.reference[std::size_t(c)] = 0;
    for (std::size_t s = 0; s < S.stabs.size(); ++s) {
        if (S.stabs[s].type != basis) continue;
        const std::string nm = std::string(1, basis) + "@(" + std::to_string(S.stabs[s].x) + "," + std::to_string(S.stabs[s].y) + ")";
        E.detectors.push_back({mb[s][0]}); E.detector_label.push_back(nm + " r0"); E.detector_round.push_back(0);
        for (int r = 1; r < rounds; ++r) {
            E.detectors.push_back({mb[s][std::size_t(r)], mb[s][std::size_t(r - 1)]});
            E.detector_label.push_back(nm + " r" + std::to_string(r)); E.detector_round.push_back(r);
        }
        std::vector<int> fin{mb[s][std::size_t(rounds - 1)]};
        for (int k : S.stabs[s].support) fin.push_back(db[std::size_t(k)]);
        E.detectors.push_back(fin); E.detector_label.push_back(nm + " final"); E.detector_round.push_back(rounds);
    }
    for (int k : (basis == 'Z' ? S.logical_z : S.logical_x)) E.observable.push_back(db[std::size_t(k)]);
    return E;
}
inline MemoryExperiment memory_z_experiment(const SurfaceCode& S, int rounds) { return memory_experiment(S, rounds, 'Z'); }

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

// Counter-based seeding for the frame Monte Carlo (thread-count independent).
//   splitmix64_mix: the SplitMix64 finaliser, a bijection on 64-bit words.
//   frame_shot_seed(seed, shot): the Xoshiro256 seed for one shot. A frame's stream is keyed by mix(seed), so adjacent
//   frame seeds (the review programs use seed + rounds) get unrelated keys; shot k then takes SplitMix64 outputs
//   4k+1 .. 4k+4 of that key's sequence (Xoshiro256's constructor consumes four, one golden-ratio step each). Shots of
//   one frame therefore never share state, and every shot's draws depend only on (seed, shot), not on which thread
//   ran it or in what order. See docs/design/deterministic-mc.md.
inline std::uint64_t splitmix64_mix(std::uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31);
}
inline std::uint64_t frame_shot_seed(std::uint64_t seed, std::uint64_t shot) {
    return splitmix64_mix(seed ^ 0x6A09E667F3BCC909ULL) + shot * (4 * 0x9E3779B97F4A7C15ULL);
}
// The impact-injection stream of one shot (ImpactSpec): keyed like frame_shot_seed but on a different constant, so it
// never overlaps the shot's noise stream.
inline std::uint64_t impact_shot_seed(std::uint64_t seed, std::uint64_t shot) {
    return splitmix64_mix(seed ^ 0xBB67AE8584CAA73BULL) + shot * (4 * 0x9E3779B97F4A7C15ULL);
}

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
        std::vector<std::uint8_t> known;  // per qubit during this slice: 0 |0>, 1 |1>, 2 unknown (impact injection)
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
    std::vector<std::uint8_t> kn(std::size_t(F.n), 0);   // 0: |0>, 1: |1>, 2: superposition / unknown
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
        fm.known.assign(kn.begin(), kn.end());
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

// Straightforward version: every fault is propagated through every later operation of the whole chip. Kept as the
// reference that the fast enumeration below is tested against.
inline std::vector<FaultSignature> enumerate_faults_reference(const FrameProgram& F, const MemoryExperiment& E) {
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

namespace frame_detail {
// One single-fault mechanism: Paulis placed after the operations of `moment` (kind 0) or a readout flip (kind 1).
// owner: the qubits of the operation the fault belongs to (both qubits of a CZ, even for a Pauli on one of them).
struct FaultDesc { int moment; int kind; int cbit; int np; int P[2]; int q[2]; double p; std::string what; int owner[2]{-1, -1}; };

// Every fault mechanism of the program, in the order of enumerate_faults_reference.
inline std::vector<FaultDesc> list_faults(const FrameProgram& F) {
    std::vector<FaultDesc> list;
    auto pauli = [&](std::size_t m, std::initializer_list<std::pair<int, int>> ps, double p, std::string what, int oa = -1, int ob = -1) {
        if (p <= 0) return;
        FaultDesc d{int(m), 0, -1, 0, {0, 0}, {0, 0}, p, std::move(what)};
        for (auto [P, q] : ps) { d.P[d.np] = P; d.q[d.np] = q; ++d.np; }
        d.owner[0] = oa >= 0 ? oa : d.q[0]; d.owner[1] = oa >= 0 ? ob : (d.np > 1 ? d.q[1] : -1);
        list.push_back(std::move(d));
    };
    for (std::size_t m = 0; m < F.moments.size(); ++m) {
        const auto& fm = F.moments[m];
        for (const auto& o : fm.ops) {
            if (o.op == Op::MEASURE && o.p > 0) list.push_back({int(m), 1, o.cbit, 0, {0, 0}, {0, 0}, o.p, "readout flip q" + std::to_string(o.a), {o.a, -1}});
            if (o.op == Op::RESET) pauli(m, {{1, o.a}}, o.p, "reset error q" + std::to_string(o.a));
        }
        for (int q : fm.ctrl1) for (int P = 1; P <= 3; ++P) pauli(m, {{P, q}}, F.p1 / 3, "1q control q" + std::to_string(q));
        for (const auto& c2 : fm.ctrl2)
            for (int k = 1; k < 16; ++k) {
                const std::string what = "CZ control q" + std::to_string(c2.a) + ",q" + std::to_string(c2.b);
                if (k % 4 && k / 4) pauli(m, {{k % 4, c2.a}, {k / 4, c2.b}}, c2.p / 15, what, c2.a, c2.b);
                else if (k % 4) pauli(m, {{k % 4, c2.a}}, c2.p / 15, what, c2.a, c2.b);
                else pauli(m, {{k / 4, c2.b}}, c2.p / 15, what, c2.a, c2.b);
            }
        if (F.dqlr_err > 0) for (int q : fm.dqlr) for (int P = 1; P <= 3; ++P) pauli(m, {{P, q}}, F.dqlr_err / 3, "DQLR q" + std::to_string(q));
        for (const auto& t : fm.zz) {
            pauli(m, {{3, t.a}}, t.pa, "ZZ");
            pauli(m, {{3, t.b}}, t.pb, "ZZ");
            pauli(m, {{3, t.a}, {3, t.b}}, t.pab, "ZZ");
        }
        for (const auto& dq : fm.deco) {
            pauli(m, {{1, dq.q}}, dq.px, "decoherence q" + std::to_string(dq.q));
            pauli(m, {{2, dq.q}}, dq.py, "decoherence q" + std::to_string(dq.q));
            pauli(m, {{3, dq.q}}, dq.pz, "decoherence q" + std::to_string(dq.q));
        }
    }
    return list;
}

// Index of a program for sparse propagation: operations of each moment by qubit (in program order), the moments in
// which each qubit has an operation, detectors of each classical bit, observable membership.
struct SparseIndex {
    std::size_t nq{0}, nm{0}, ncl{0}, ndet{0};
    std::vector<std::vector<std::vector<int>>> ops_of;
    std::vector<std::vector<int>> active, dets_of;
    std::vector<char> in_obs;
    SparseIndex(const FrameProgram& F, const MemoryExperiment& E)
        : nq(std::size_t(F.n)), nm(F.moments.size()), ncl(std::size_t(std::max(F.nclbits, 1))), ndet(E.detectors.size()) {
        ops_of.assign(nm, std::vector<std::vector<int>>(nq));
        for (std::size_t m = 0; m < nm; ++m)
            for (std::size_t k = 0; k < F.moments[m].ops.size(); ++k) {
                const auto& o = F.moments[m].ops[k];
                ops_of[m][std::size_t(o.a)].push_back(int(k));
                if (o.op == Op::CZ) ops_of[m][std::size_t(o.b)].push_back(int(k));
            }
        dets_of.assign(ncl, {});
        for (std::size_t d = 0; d < E.detectors.size(); ++d) for (int c : E.detectors[d]) dets_of[std::size_t(c)].push_back(int(d));
        in_obs.assign(ncl, 0);
        for (int c : E.observable) in_obs[std::size_t(c)] ^= 1;
        active.assign(nq, {});
        for (std::size_t m = 0; m < nm; ++m) for (std::size_t q = 0; q < nq; ++q) if (!ops_of[m][q].empty()) active[q].push_back(int(m));
    }
};

// Per-thread workspace: a sparse Pauli frame. Only the operations on qubits that carry a Pauli are applied, and
// moments without an operation on those qubits are skipped altogether.
struct SparseFrame {
    const FrameProgram* F;
    const SparseIndex* I;
    std::vector<std::uint8_t> x, z, flips, dpar;
    std::vector<char> in_sup;
    std::vector<std::size_t> next_of;                        // per qubit: index into active[q] of its next moment
    std::vector<int> sup, cb_touched, det_touched, pending, extra;
    int cur{0};                                              // the last moment processed
    std::uint8_t obs{0};
    long steps{0};                                           // moments processed (a work count for tests)
    SparseFrame(const FrameProgram& F_, const SparseIndex& I_)
        : F(&F_), I(&I_), x(I_.nq, 0), z(I_.nq, 0), flips(I_.ncl, 0), dpar(I_.ndet, 0), in_sup(I_.nq, 0), next_of(I_.nq, 0) {}
    void add(int q) {
        if (in_sup[std::size_t(q)]) return;
        in_sup[std::size_t(q)] = 1; sup.push_back(q);
        const auto& a = I->active[std::size_t(q)];
        next_of[std::size_t(q)] = std::size_t(std::upper_bound(a.begin(), a.end(), cur) - a.begin());
    }
    void place(int P, int q) { apply_pauli(P, q, x.data(), z.data()); add(q); }
    // apply the moments cur+1 .. stop (inclusive)
    void run(int stop) {
        const int none = std::numeric_limits<int>::max();
        for (;;) {
            // the next moment with an operation on a qubit that carries a Pauli
            int nxt = none;
            for (int q : sup) {
                if (!(x[std::size_t(q)] | z[std::size_t(q)])) continue;
                const auto& a = I->active[std::size_t(q)];
                std::size_t& p = next_of[std::size_t(q)];
                while (p < a.size() && a[p] <= cur) ++p;
                if (p < a.size()) nxt = std::min(nxt, a[p]);
            }
            if (nxt == none || nxt > stop) { cur = std::max(cur, stop); return; }
            cur = nxt;
            ++steps;
            const std::size_t m = std::size_t(cur);
            // operations of this moment on the qubits that carry a Pauli, in program order; a CZ can bring its
            // partner in, whose later operations in the same moment then join the list
            pending.clear();
            for (int q : sup) if (x[std::size_t(q)] | z[std::size_t(q)]) for (int k : I->ops_of[m][std::size_t(q)]) pending.push_back(k);
            std::sort(pending.begin(), pending.end());
            pending.erase(std::unique(pending.begin(), pending.end()), pending.end());
            for (std::size_t n = 0; n < pending.size(); ++n) {
                const int k = pending[n];
                const auto& o = F->moments[m].ops[std::size_t(k)];
                if (o.op == Op::CZ) {
                    for (int q : {o.a, o.b}) {
                        add(q);
                        extra.clear();
                        for (int k2 : I->ops_of[m][std::size_t(q)])
                            if (k2 > k && !std::binary_search(pending.begin() + std::ptrdiff_t(n) + 1, pending.end(), k2)) extra.push_back(k2);
                        if (!extra.empty()) {
                            pending.insert(pending.end(), extra.begin(), extra.end());
                            std::sort(pending.begin() + std::ptrdiff_t(n) + 1, pending.end());
                        }
                    }
                } else if (o.op == Op::MEASURE) cb_touched.push_back(o.cbit);
                apply_op(o, x.data(), z.data(), flips.data());
            }
        }
    }
    void flip_bit(int c) { flips[std::size_t(c)] ^= 1; cb_touched.push_back(c); }
    void toggle_det(int d) { if (!dpar[std::size_t(d)]) det_touched.push_back(d); dpar[std::size_t(d)] ^= 1; }
    // fold the flipped classical bits into detector parities and the observable
    void fold_bits() {
        for (int c : cb_touched) if (flips[std::size_t(c)]) {
            flips[std::size_t(c)] = 0;                                   // also clears duplicates in cb_touched
            obs ^= std::uint8_t(I->in_obs[std::size_t(c)]);
            for (int dt : I->dets_of[std::size_t(c)]) toggle_det(dt);
        }
        cb_touched.clear();
    }
    // the fired detectors (ascending); clears the detector parities
    void take_dets(std::vector<int>& out) {
        std::sort(det_touched.begin(), det_touched.end());
        det_touched.erase(std::unique(det_touched.begin(), det_touched.end()), det_touched.end());
        for (int dt : det_touched) { if (dpar[std::size_t(dt)]) out.push_back(dt); dpar[std::size_t(dt)] = 0; }
        det_touched.clear();
    }
    void clear_frame() {
        for (int q : sup) { x[std::size_t(q)] = 0; z[std::size_t(q)] = 0; in_sup[std::size_t(q)] = 0; }
        sup.clear();
    }
};
} // namespace frame_detail

// Sparse version, same output in the same order. A single fault touches a handful of qubits, so only the operations
// on the qubits that currently carry a Pauli are applied, and only the detectors of the classical bits it flipped are
// evaluated. Moments without an operation on those qubits are skipped altogether: a chip with per-pair CZ durations
// is cut into several times more time slices than it has gate layers, and most of them only carry idle noise.
// A data-qubit Pauli survives every later round, so each fault still costs time proportional to the rounds after it
// and the whole enumeration is quadratic in the number of rounds. Kept for timing comparisons (enumerate_faults is
// the linear-time version).
// steps (optional): total time slices processed, a deterministic work count.
inline std::vector<FaultSignature> enumerate_faults_direct(const FrameProgram& F, const MemoryExperiment& E, long* steps = nullptr) {
    const auto list = frame_detail::list_faults(F);
    const frame_detail::SparseIndex I(F, E);
    std::vector<FaultSignature> out(list.size());
    long total = 0;
#pragma omp parallel reduction(+ : total)
    {
        frame_detail::SparseFrame S(F, I);
#pragma omp for schedule(dynamic, 256)
        for (long i = 0; i < long(list.size()); ++i) {
            const auto& d = list[std::size_t(i)];
            S.cur = d.moment;
            if (d.kind == 1) S.flip_bit(d.cbit);
            else {
                for (int k = 0; k < d.np; ++k) S.place(d.P[k], d.q[k]);
                S.run(int(I.nm) - 1);
            }
            S.fold_bits();
            FaultSignature s; s.p = d.p; s.moment = d.moment; s.what = d.what;
            S.take_dets(s.dets);
            s.obs = S.obs != 0;
            S.obs = 0;
            S.clear_frame();
            out[std::size_t(i)] = std::move(s);
        }
        total += S.steps;
    }
    if (steps) *steps = total;
    return out;
}

// Linear-time version, same output in the same order (docs/design/linear-fault-enumeration.md).
// Pauli-frame propagation is linear over GF(2), so the effect of a frame on the rest of the circuit is the XOR of the
// effects of its single-qubit X and Z components. The program is cut after every moment that contains a measurement.
// For each cut c and each qubit q, tail[c][q][X or Z] holds the detectors and the logical flipped by that one Pauli,
// placed after cut c, from there to the end of the program. The tails are built backwards: propagate the Pauli to the
// next cut, collect the bits it flips there, and XOR in the tails of whatever frame is left. A fault is then
// propagated only to the first cut at or after it and finished with tails. Each propagation spans one segment (about
// one QEC round), so the enumeration costs O(rounds x faults per round) instead of O(rounds^2 x faults per round).
// steps (optional): total time slices processed (tails and faults), a deterministic work count.
inline std::vector<FaultSignature> enumerate_faults(const FrameProgram& F, const MemoryExperiment& E, long* steps = nullptr) {
    const auto list = frame_detail::list_faults(F);
    const frame_detail::SparseIndex I(F, E);
    std::vector<int> cuts;                                       // moments that contain a measurement
    for (std::size_t m = 0; m < I.nm; ++m)
        for (const auto& o : F.moments[m].ops) if (o.op == Op::MEASURE) { cuts.push_back(int(m)); break; }
    const std::size_t nc = cuts.size();
    std::vector<int> cut_at(I.nm, -1);                           // index of the first cut at or after moment m
    for (std::size_t m = 0, k = 0; m < I.nm; ++m) { while (k < nc && cuts[k] < int(m)) ++k; cut_at[m] = k < nc ? int(k) : -1; }
    struct Tail { std::vector<int> dets; std::uint8_t obs{0}; };
    std::vector<Tail> tail(nc * I.nq * 2);
    auto tail_of = [&](std::size_t c, int q, int t) -> Tail& { return tail[(c * I.nq + std::size_t(q)) * 2 + std::size_t(t)]; };
    // finish a frame that has been propagated up to cut c: XOR in the tails of its remaining Paulis
    auto finish = [&](frame_detail::SparseFrame& S, int c) {
        if (c < 0) return;
        for (int q : S.sup)
            for (int t = 0; t < 2; ++t) {
                if (!(t == 0 ? S.x[std::size_t(q)] : S.z[std::size_t(q)])) continue;
                const Tail& T = tail_of(std::size_t(c), q, t);
                S.obs ^= T.obs;
                for (int dt : T.dets) S.toggle_det(dt);
            }
    };
    // propagate the frame in S from moment `from` to the first cut after it (or the end), fold, finish with tails
    auto segment = [&](frame_detail::SparseFrame& S, int next_cut) {
        S.run(next_cut >= 0 ? cuts[std::size_t(next_cut)] : int(I.nm) - 1);
        S.fold_bits();
        finish(S, next_cut);
    };
    // one workspace per thread for all levels (allocating it per level would cost O(rounds^2)); the implicit barrier
    // at the end of each omp for orders the levels
    long total = 0;
#pragma omp parallel reduction(+ : total)
    {
        frame_detail::SparseFrame S(F, I);
        for (std::size_t c = nc; c-- > 0;) {
            const int next = c + 1 < nc ? int(c + 1) : -1;
#pragma omp for schedule(dynamic, 8)
            for (long j = 0; j < long(I.nq * 2); ++j) {
                const int q = int(j / 2), t = int(j % 2);
                S.cur = cuts[c];
                S.place(t == 0 ? 1 : 3, q);
                segment(S, next);
                Tail& T = tail_of(c, q, t);
                S.take_dets(T.dets);
                T.obs = S.obs; S.obs = 0;
                S.clear_frame();
            }
        }
        total += S.steps;
    }
    std::vector<FaultSignature> out(list.size());
#pragma omp parallel reduction(+ : total)
    {
        frame_detail::SparseFrame S(F, I);
#pragma omp for schedule(dynamic, 256)
        for (long i = 0; i < long(list.size()); ++i) {
            const auto& d = list[std::size_t(i)];
            S.cur = d.moment;
            if (d.kind == 1) { S.flip_bit(d.cbit); S.fold_bits(); }
            else {
                for (int k = 0; k < d.np; ++k) S.place(d.P[k], d.q[k]);
                const int c = cut_at[std::size_t(d.moment)];
                if (c >= 0 && cuts[std::size_t(c)] == d.moment) finish(S, c);   // placed right after a cut
                else segment(S, c);
            }
            FaultSignature s; s.p = d.p; s.moment = d.moment; s.what = d.what;
            S.take_dets(s.dets);
            s.obs = S.obs != 0;
            S.obs = 0;
            S.clear_frame();
            out[std::size_t(i)] = std::move(s);
        }
        total += S.steps;
    }
    if (steps) *steps = total;
    return out;
}

struct DEMEdge { int a, b; bool obs; double p; };   // b = -1: boundary
struct DetectorErrorModel {
    int n_det{0};
    std::vector<DEMEdge> edges;
    int hyperedge_faults{0};          // faults flipping > 2 decoded detectors (dropped from the graph)
    double hyperedge_p{0};            // their summed probability (error mass the decoder cannot see)
    int undetectable_logical{0};      // faults flipping the logical but no detector (distance-killing; must be 0)
};

inline DetectorErrorModel build_dem(const std::vector<FaultSignature>& faults, int n_det) {
    DetectorErrorModel M; M.n_det = n_det;
    std::map<std::tuple<int, int, int>, double> acc;
    for (const auto& f : faults) {
        if (f.dets.empty()) { if (f.obs) ++M.undetectable_logical; continue; }
        if (f.dets.size() > 2) { ++M.hyperedge_faults; M.hyperedge_p += f.p; continue; }
        const int a = f.dets[0], b = f.dets.size() == 2 ? f.dets[1] : -1;
        auto key = std::make_tuple(a, b, int(f.obs));
        double& p = acc[key];
        p = p * (1 - f.p) + f.p * (1 - p);
    }
    for (const auto& [k, p] : acc) M.edges.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k) != 0, p});
    return M;
}

// ================================================================== minimum-weight perfect matching decoder
// Defects are matched to each other or to the boundary along all-pairs shortest paths of the detector graph (the
// boundary is a node of that graph, so "both to the boundary" is one of the paths between two defects).
//   DP       exact bitmask dynamic programme, O(2^k k): fastest for a few defects
//   Blossom  exact for any number of defects (matching.hpp), O(k^3)
//   Auto     DP up to dp_max defects, blossom beyond
//   Greedy   closest-pair heuristic (not exact)
//   Legacy   the v0.3 behaviour, DP up to 22 defects and greedy beyond; kept for regression comparisons
//   UnionFind  weighted union-find (unionfind.hpp): almost-linear time, not exact; the real-time candidate
enum class DecodeMethod { Auto, DP, Blossom, Greedy, Legacy, UnionFind };

struct DecodeScratch {
    std::vector<double> cost;
    std::vector<int> choice, nodes;
    BlossomMatcher bm;
    std::vector<BlossomMatcher::W> w;
    UFScratch uf;
    struct Herald { std::vector<std::pair<int, double>> scaled; std::vector<int> term; std::vector<double> d; std::vector<std::uint8_t> p; } herald;
    double weight{0};                          // total weight of the last matching
    DecodeMethod used{DecodeMethod::Auto};     // what the last call ran
};

class MatchingDecoder {
public:
    explicit MatchingDecoder(const DetectorErrorModel& M) : N_(M.n_det), V_(M.n_det + 1), uf_(M.n_det, M.edges) {
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
        double mx = 1;
        for (double v : dist_) if (v < INF) mx = std::max(mx, v);
        cap_ = 4 * mx + 16;                    // stands in for "unreachable" in the integer matching
    }
    int n_detectors() const { return N_; }
    double dist(int a, int b) const { return dist_[idx(a, b < 0 ? N_ : b)]; }   // b = -1: the boundary
    void set_dp_max(int k) { dp_max_ = std::clamp(k, 0, 24); }
    const UnionFindDecoder& union_find() const { return uf_; }
    int dp_max() const { return dp_max_; }

    // predicted logical flip for a set of fired detectors
    int decode(const std::vector<int>& def, DecodeScratch& ws, DecodeMethod method = DecodeMethod::Auto) const {
        const int k = int(def.size());
        ws.weight = 0; ws.used = method;
        if (method == DecodeMethod::UnionFind) {
            const int parity = uf_.decode(def, ws.uf);
            for (int e : ws.uf.correction) ws.weight += uf_.edge[std::size_t(e)].w;
            return parity;
        }
        if (k == 0) return 0;
        auto D = [&](int i, int j) { return dist_[idx(i < k ? def[std::size_t(i)] : N_, j < k ? def[std::size_t(j)] : N_)]; };
        auto P = [&](int i, int j) { return int(par_[idx(i < k ? def[std::size_t(i)] : N_, j < k ? def[std::size_t(j)] : N_)]); };
        return match(k, D, P, ws, method);
    }

    // Edges touched by heralds in the "detector" rule: for each heralded bit, given as the detectors built from it (one
    // or two; -1 for none), the edge joining them (that bit's measurement error) gets `erase_factor`, every other edge
    // touching one of them `touch_factor`. Appends (edge, factor) to out.
    void herald_touch_edges(const std::vector<std::pair<int, int>>& heralded, double erase_factor, double touch_factor,
                            std::vector<std::pair<int, double>>& out) const {
        for (auto [a, b] : heralded) {
            const int ro = b >= 0 ? uf_.edge_index(a, b) : uf_.edge_index(a, -1);
            for (int h : {a, b}) {
                if (h < 0) continue;
                for (int e : uf_.incident(h)) out.push_back({e, e == ro ? erase_factor : touch_factor});
            }
        }
    }

    // Decoding with modified edge weights, for heralds (leakage seen by 3-state readout, docs/design/heralded-leakage.md):
    // scaled lists (edge index of union_find().edge, factor); each edge's weight is multiplied by its smallest listed
    // factor (0: erased, the edge is free). Shortest paths of the modified graph between the defects are the closure of
    // the static all-pairs distances and the modified edges, taken over the defects, the modified edges' ends and the
    // boundary (Floyd-Warshall on that small set), so the matching stays exact on the modified graph. The union-find
    // decoder takes the modified lengths directly.
    int decode_scaled(const std::vector<int>& def, const std::vector<std::pair<int, double>>& scaled, DecodeScratch& ws,
                      DecodeMethod method = DecodeMethod::Auto) const {
        if (scaled.empty()) return decode(def, ws, method);
        ws.weight = 0; ws.used = method;
        auto& H = ws.herald;
        H.scaled = scaled;
        std::sort(H.scaled.begin(), H.scaled.end());
        {
            std::size_t w = 0;
            for (std::size_t i = 0; i < H.scaled.size(); ++i)
                if (w > 0 && H.scaled[w - 1].first == H.scaled[i].first) H.scaled[w - 1].second = std::min(H.scaled[w - 1].second, H.scaled[i].second);
                else H.scaled[w++] = H.scaled[i];
            H.scaled.resize(w);
            // edges left at their own weight need no work
            w = 0;
            for (std::size_t i = 0; i < H.scaled.size(); ++i) if (H.scaled[i].second < 1.0) H.scaled[w++] = H.scaled[i];
            H.scaled.resize(w);
        }
        if (H.scaled.empty()) return decode(def, ws, method);
        if (method == DecodeMethod::UnionFind) {
            const int parity = uf_.decode(def, ws.uf, &H.scaled);
            for (int e : ws.uf.correction) ws.weight += uf_.edge[std::size_t(e)].w;
            return parity;
        }
        const int k = int(def.size());
        if (k == 0) return 0;
        // terminals: the defects (0..k-1), the boundary (k), then the ends of the modified edges
        H.term.assign(def.begin(), def.end());
        H.term.push_back(N_);
        for (const auto& ef : H.scaled) for (int v : {uf_.edge[std::size_t(ef.first)].u, uf_.edge[std::size_t(ef.first)].v}) H.term.push_back(v);
        std::sort(H.term.begin() + k + 1, H.term.end());
        H.term.erase(std::unique(H.term.begin() + k + 1, H.term.end()), H.term.end());
        // drop extra terminals that are defects or the boundary already
        {
            std::size_t w = std::size_t(k) + 1;
            for (std::size_t i = w; i < H.term.size(); ++i) {
                const int v = H.term[i];
                if (v == N_ || std::find(def.begin(), def.end(), v) != def.end()) continue;
                H.term[w++] = v;
            }
            H.term.resize(w);
        }
        const std::size_t T = H.term.size();
        H.d.assign(T * T, 0.0); H.p.assign(T * T, 0);
        for (std::size_t i = 0; i < T; ++i) for (std::size_t j = 0; j < T; ++j) { H.d[i * T + j] = dist_[idx(H.term[i], H.term[j])]; H.p[i * T + j] = par_[idx(H.term[i], H.term[j])]; }
        auto pos = [&](int v) { return std::size_t(std::find(H.term.begin(), H.term.end(), v) - H.term.begin()); };
        for (const auto& ef : H.scaled) {
            const auto& E = uf_.edge[std::size_t(ef.first)];
            const std::size_t i = pos(E.u), j = pos(E.v);
            const double w = ef.second * E.w;
            if (w < H.d[i * T + j]) { H.d[i * T + j] = H.d[j * T + i] = w; H.p[i * T + j] = H.p[j * T + i] = E.obs; }
        }
        for (std::size_t m = 0; m < T; ++m) for (std::size_t i = 0; i < T; ++i) {
            const double dim = H.d[i * T + m];
            for (std::size_t j = 0; j < T; ++j) {
                const double c = dim + H.d[m * T + j];
                if (c < H.d[i * T + j]) { H.d[i * T + j] = c; H.p[i * T + j] = H.p[i * T + m] ^ H.p[m * T + j]; }
            }
        }
        auto D = [&](int i, int j) { return H.d[std::size_t(i) * T + std::size_t(j)]; };
        auto P = [&](int i, int j) { return int(H.p[std::size_t(i) * T + std::size_t(j)]); };
        return match(k, D, P, ws, method);
    }

private:
    int N_, V_, dp_max_{8};        // measured crossover: the DP wins only up to about 8 defects
    double cap_{0};
    UnionFindDecoder uf_;
    std::vector<double> dist_;
    std::vector<std::uint8_t> par_;
    std::size_t idx(int i, int j) const { return std::size_t(i) * std::size_t(V_) + std::size_t(j); }

    // Matching of k defects on a distance oracle D(i, j) and parity oracle P(i, j) over local indices 0..k-1, with k
    // the boundary.
    template <class Dist, class Par>
    int match(int k, const Dist& D, const Par& P, DecodeScratch& ws, DecodeMethod method) const {
        if (method == DecodeMethod::Legacy) method = k > 22 ? DecodeMethod::Greedy : DecodeMethod::DP;
        if (method == DecodeMethod::Auto) method = k <= dp_max_ ? DecodeMethod::DP : DecodeMethod::Blossom;
        if (method == DecodeMethod::DP && k > 24) method = DecodeMethod::Blossom;   // the DP table would not fit
        ws.used = method;
        if (method == DecodeMethod::Greedy) return greedy(k, D, P, ws.weight);
        if (method == DecodeMethod::Blossom) return blossom(k, D, P, ws);
        return dp(k, D, P, ws);
    }
    template <class Dist, class Par>
    int dp(int k, const Dist& D, const Par& P, DecodeScratch& ws) const {
        auto& cost = ws.cost; auto& choice = ws.choice;
        const std::size_t full = (std::size_t(1) << k) - 1;
        cost.assign(full + 1, 0.0);
        choice.assign(full + 1, -1);
        for (std::size_t mask = 1; mask <= full; ++mask) {
            int i = 0; while (!((mask >> i) & 1)) ++i;
            const std::size_t rest = mask & ~(std::size_t(1) << i);
            double best = D(i, k) + cost[rest];
            int ch = -1;
            for (int j = i + 1; j < k; ++j) {
                if (!((rest >> j) & 1)) continue;
                const double c = D(i, j) + cost[rest & ~(std::size_t(1) << j)];
                if (c < best) { best = c; ch = j; }
            }
            cost[mask] = best; choice[mask] = ch;
        }
        ws.weight = cost[full];
        int parity = 0;
        std::size_t mask = full;
        while (mask) {
            int i = 0; while (!((mask >> i) & 1)) ++i;
            const int ch = choice[mask];
            if (ch < 0) { parity ^= P(i, k); mask &= ~(std::size_t(1) << i); }
            else { parity ^= P(i, ch); mask &= ~((std::size_t(1) << i) | (std::size_t(1) << ch)); }
        }
        return parity;
    }
    // Perfect matching of the defects (plus the boundary node when their number is odd) on the complete graph of
    // shortest-path distances. Distances are scaled to integers; subtracting them from a constant larger than any
    // matching's total turns "minimum-weight perfect" into "maximum-weight".
    template <class Dist, class Par>
    int blossom(int k, const Dist& D, const Par& P, DecodeScratch& ws) const {
        using W = BlossomMatcher::W;
        ws.nodes.resize(std::size_t(k));
        for (int i = 0; i < k; ++i) ws.nodes[std::size_t(i)] = i;
        if (k % 2) ws.nodes.push_back(k);
        const int n = int(ws.nodes.size());
        const std::size_t sn = std::size_t(n);
        const double scale = 1048576.0;
        ws.w.assign(sn * sn, 0);
        W wmax = 0;
        for (std::size_t i = 0; i < sn; ++i) for (std::size_t j = i + 1; j < sn; ++j) {
            const W wi = static_cast<W>(std::llround(std::min(D(ws.nodes[i], ws.nodes[j]), cap_) * scale));
            ws.w[i * sn + j] = ws.w[j * sn + i] = wi;
            wmax = std::max(wmax, wi);
        }
        const W big = (static_cast<W>(n) / 2 + 1) * wmax + 1;
        for (std::size_t i = 0; i < sn; ++i) for (std::size_t j = 0; j < sn; ++j) if (i != j) ws.w[i * sn + j] = big - ws.w[i * sn + j];
        const auto& mate = ws.bm.solve(n, ws.w);
        int parity = 0;
        for (std::size_t i = 0; i < sn; ++i) {
            const int j = mate[i];
            if (j < 0) { ws.used = DecodeMethod::Greedy; return greedy(k, D, P, ws.weight); }   // cannot happen on a complete graph
            if (int(i) < j) { parity ^= P(ws.nodes[i], ws.nodes[std::size_t(j)]); ws.weight += D(ws.nodes[i], ws.nodes[std::size_t(j)]); }
        }
        return parity;
    }
    template <class Dist, class Par>
    int greedy(int k, const Dist& D, const Par& P, double& weight) const {
        std::vector<int> def(static_cast<std::size_t>(k));
        for (int i = 0; i < k; ++i) def[std::size_t(i)] = i;
        int parity = 0;
        weight = 0;
        while (!def.empty()) {
            double best = 1e300; int bi = 0, bj = -1;
            for (std::size_t i = 0; i < def.size(); ++i) {
                if (D(def[i], k) < best) { best = D(def[i], k); bi = int(i); bj = -1; }
                for (std::size_t j = i + 1; j < def.size(); ++j) if (D(def[i], def[j]) < best) { best = D(def[i], def[j]); bi = int(i); bj = int(j); }
            }
            weight += best;
            if (bj < 0) { parity ^= P(def[std::size_t(bi)], k); def.erase(def.begin() + bi); }
            else { parity ^= P(def[std::size_t(bi)], def[std::size_t(bj)]); def.erase(def.begin() + bj); def.erase(def.begin() + bi); }
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
    double hyperedge_p{0};
    double leak_frac_data{0}, leak_frac_meas{0};   // time-averaged fraction of leaked qubits (frame model with leakage)
    double leak_events_per_shot{0};                // leakage injections per shot
    std::vector<double> detection_by_round;        // detection fraction per syndrome round (last entry: final data readout)
    double mean_defects{0};                        // fired detectors per shot
    int max_defects{0};
    double frac_over_22{0};                        // shots with more than 22 defects (decoded greedily before v0.4)
    double decode_us_per_shot{0};                  // software decoding time per shot on one core
    double heralds_per_shot{0};                    // heralded measure-qubit readouts per shot (HeraldSpec)
    double herald_true_frac{0};                    // fraction of those that were leaked
    // impact injection and veto (ImpactSpec, VetoSpec); logical_errors counts every shot, kept or not
    long impacts{0}, impacts_vetoed{0}, errors_impacted{0}; // shots with an impact; those discarded; their failures
    long errors_impacted_kept{0};                  // failures of impacted shots that the veto kept (silent impact failures)
    long discarded{0}, errors_kept{0};             // shots discarded by the veto; failures among the kept shots
    double p_logical_kept{0}, eps_kept{0};         // over the kept shots
    std::vector<long> veto_hist;                   // shots by burst statistic (VetoSpec), index = statistic
};

// ================================================================== correlated error bursts (high-energy impacts)
// A cosmic-ray or gamma impact fills the chip with quasiparticles: every qubit gains an excess decay (and/or
// dephasing) rate that then relaxes exponentially. The burst lasts for thousands of QEC rounds, so it is treated as a
// slowly varying severity s(t) = exp(-t / tau) on top of the normal noise:
//     1/T1(t) = 1/T1 + s(t) / t1_peak,      1/T_phi(t) = 1/T_phi + s(t) / tphi_peak.
// The decoder is not told: it keeps the priors of the quiet chip. The logical error per round is sampled at a ladder
// of severities and integrated over the burst.
struct BurstSpec {
    double rate_hz{0.1};        // impacts per second on the chip
    double tau_ms{25};          // recovery time constant of the excess error rates
    double t1_peak_us{1.0};     // excess energy-decay time at the peak of the burst (<= 0: T1 unaffected)
    double tphi_peak_us{0.0};   // excess pure-dephasing time at the peak (<= 0: dephasing unaffected)
    double radius{0.0};         // per-shot injection only: spatial extent, Gaussian sigma in lattice units (<= 0: whole chip)
};

// Per-shot impact injection (docs/design/impact-veto.md). With probability p_shot a shot contains an impact: its onset
// time is uniform over the shot and its centre is a uniformly chosen qubit. From the onset, every qubit q gains the
// excess rates of BurstSpec at severity s = exp(-(t - t0) / tau) * g(q), with g = exp(-r^2 / (2 radius^2)) of its
// lattice distance r from the centre (g = 1 for radius <= 0). The excess decay and dephasing are applied as Pauli
// channels with the same state-aware rule as the quiet noise (a qubit known to be in |1> decays as a bit flip, one in
// |0> is untouched). Every draw of the impact uses its own generator keyed on (seed, shot) (impact_shot_seed), so a
// shot without an impact is bit for bit the same shot as in a run without injection.
struct ImpactSpec {
    BurstSpec burst;
    double p_shot{0};
};

// Impact veto: the largest number of detection events in a sliding space-time window (detectors whose stabilisers lie
// within `radius` lattice units of a centre stabiliser, Chebyshev distance, over `window_rounds` consecutive rounds) is
// the shot's burst statistic. A shot whose statistic reaches `threshold` is discarded (post-selection).
// calibrate_veto picks the threshold for a target false-veto rate on burst-free shots.
struct VetoSpec {
    bool on{false};
    int window_rounds{3};
    double radius{2.0};
    int threshold{std::numeric_limits<int>::max()};   // max(): record the statistic only
};

// Heralded leakage (docs/design/heralded-leakage.md): 3-state readout of the measure qubits. A leaked measure qubit is
// reported as |2> with probability `fidelity`; an unleaked one is misreported as |2> with probability `false_rate`.
// The bit value passed to the detectors is unchanged (a leaked qubit still reads as |1>). The decoder lowers the
// weights of the graph edges a heralded bit points at (factors in [0, 1]; 0 erases the edge), by one of two rules:
//   footprint (default)  the edges of every single fault on that measure qubit's operations in the round that ended
//                        with the heralded measurement (its CZs, both qubits; its single-qubit gates and idling) get
//                        `footprint_factor`; that bit's own measurement-error edge gets `erase_factor`. This is where
//                        a leaked measure qubit puts its errors, including the data errors that a leaked X-type
//                        measure qubit causes in a Z memory, whose bits feed no detector there.
//   detector             the edge between the two detectors built from the heralded bit gets `erase_factor`, every
//                        other edge touching them `touch_factor`.
// use_in_decoder = false samples the same heralds (same random draws, so the same syndromes) and ignores them.
struct HeraldSpec {
    bool on{false};
    double fidelity{0.95};
    double false_rate{1e-3};
    bool footprint{true};
    double erase_factor{0.0};
    double footprint_factor{0.4};
    double touch_factor{1.0};
    bool use_in_decoder{true};
};

struct FrameRunOptions {
    const MatchingDecoder* decoder{nullptr};         // decode with this instead of the setup's own (model-mismatch studies)
    DecodeMethod method{DecodeMethod::Auto};
    HeraldSpec herald;
    ImpactSpec impact;
    VetoSpec veto;
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
    std::vector<std::pair<double, double>> sim_xy;   // lattice position of each sim qubit (unit: the data-qubit pitch)
    std::vector<std::pair<double, double>> det_xy;   // lattice position of each detector's stabiliser
    // herald footprints (HeraldSpec::footprint), by classical bit of a measure-qubit measurement: the decoder-graph
    // edges (union_find().edge index) of the single faults in that qubit's round, flagged 1 for the bit's own
    // measurement-error edge. Built with the decoder.
    std::vector<std::vector<std::pair<int, char>>> herald_edges;
};

// Footprint of each measure-qubit measurement in the decoder graph (see HeraldSpec): every single fault owned by the
// measure qubit (a fault on one of its operations, or on either qubit of one of its CZs) after its previous
// measurement and up to this one, plus this bit's readout flip.
inline std::vector<std::vector<std::pair<int, char>>> herald_footprints(const FrameProgram& F, const std::vector<char>& is_data,
                                                                        const std::vector<FaultSignature>& faults, const MatchingDecoder& dec) {
    std::vector<std::vector<std::pair<int, char>>> out(std::size_t(std::max(F.nclbits, 1)));
    const auto list = frame_detail::list_faults(F);
    if (list.size() != faults.size()) throw std::runtime_error("herald_footprints: fault list does not match the program");
    // measurements of each measure qubit: (moment, cbit), ascending
    std::vector<std::vector<std::pair<int, int>>> meas(std::size_t(F.n));
    for (std::size_t m = 0; m < F.moments.size(); ++m)
        for (const auto& o : F.moments[m].ops)
            if (o.op == Op::MEASURE && !is_data[std::size_t(o.a)]) meas[std::size_t(o.a)].push_back({int(m), o.cbit});
    const auto& uf = dec.union_find();
    for (std::size_t i = 0; i < list.size(); ++i) {
        const auto& f = faults[i];
        if (f.dets.empty() || f.dets.size() > 2) continue;
        const int e = uf.edge_index(f.dets[0], f.dets.size() > 1 ? f.dets[1] : -1);
        if (e < 0) continue;
        const auto& d = list[i];
        if (d.kind == 1) { out[std::size_t(d.cbit)].push_back({e, 1}); continue; }
        for (int q : d.owner) {
            if (q < 0 || is_data[std::size_t(q)]) continue;
            const auto& mq = meas[std::size_t(q)];
            // the first measurement of q after the fault (a fault placed after a measurement belongs to the next round)
            auto it = std::upper_bound(mq.begin(), mq.end(), std::make_pair(d.moment, std::numeric_limits<int>::max()));
            if (it != mq.end()) out[std::size_t(it->second)].push_back({e, 0});
        }
    }
    for (auto& v : out) {
        std::sort(v.begin(), v.end());
        // keep each edge once; the readout flag wins
        std::vector<std::pair<int, char>> u;
        for (const auto& x : v) { if (!u.empty() && u.back().first == x.first) u.back().second |= x.second; else u.push_back(x); }
        v.swap(u);
    }
    return out;
}

// build_decoder = false skips the fault enumeration, the detector error model and the decoder (the expensive part):
// for runs that are decoded with another setup's decoder.
inline MemorySetup prepare_memory(const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, int rounds, char basis = 'Z', bool build_decoder = true) {
    if (D.surface_d <= 0) throw std::runtime_error("chip is not a surface-code layout (generate it with topology 'surface')");
    MemorySetup M;
    M.code = rotated_surface_code(D.surface_d);
    M.exp = memory_experiment(M.code, rounds, basis);
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
        M.sim_xy.assign(std::size_t(M.frame.n), {0.0, 0.0});
        for (const auto& dq : M.code.data) { const int s = sim[std::size_t(D.index_of(dq.qubit))]; if (s >= 0) M.sim_xy[std::size_t(s)] = {0.5 * dq.x, 0.5 * dq.y}; }
        for (const auto& st : M.code.stabs) { const int s = sim[std::size_t(D.index_of(st.qubit))]; if (s >= 0) M.sim_xy[std::size_t(s)] = {0.5 * st.x, 0.5 * st.y}; }
        // detectors come in blocks of rounds + 1 per stabiliser of the memory basis (memory_experiment)
        for (const auto& st : M.code.stabs)
            if (st.type == basis) for (int r = 0; r <= rounds; ++r) M.det_xy.push_back({0.5 * st.x, 0.5 * st.y});
        if (M.det_xy.size() != M.exp.detectors.size()) throw std::runtime_error("prepare_memory: detector layout mismatch");
    }
    if (build_decoder) {
        const auto faults = enumerate_faults(M.frame, M.exp);
        M.dem = build_dem(faults, int(M.exp.detectors.size()));
        M.decoder = std::make_unique<MatchingDecoder>(M.dem);
        M.herald_edges = herald_footprints(M.frame, M.is_data, faults, *M.decoder);
    }
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
    r.hyperedge_p = M.dem.hyperedge_p;
    r.undetectable_logical = M.dem.undetectable_logical;
    r.ignored_exchange_pairs = M.frame.ignored_exchange_pairs;
    // The matching graph only has edges for faults that flip one or two detectors. Say so when the circuit has others:
    // their errors are still simulated but the decoder has no prior for them, so the logical error is pessimistic for
    // this decoder and not comparable to a hypergraph decoder.
    if (r.hyperedge_faults > 0)
        std::fprintf(stderr, "qec warning: d=%d, %d rounds: %d fault mechanisms flip >2 detectors (total p %.2e) and are missing from the "
                     "matching graph\n", r.d, r.rounds, r.hyperedge_faults, r.hyperedge_p);
    if (r.undetectable_logical > 0)
        std::fprintf(stderr, "qec warning: d=%d, %d rounds: %d fault mechanisms flip the logical without any detector (distance-reducing)\n",
                     r.d, r.rounds, r.undetectable_logical);
    return r;
}

// Fast Monte Carlo with Pauli frames (+ leakage: a per-qubit flag outside the frame).
//  * CZ leakage: with the pair's probability the moving qubit leaves the subspace (|11> -> |20>); its partner is
//    projected, modelled as random X and Z.
//  * a leaked qubit ignores single-qubit gates; a CZ with a leaked partner applies an uncontrolled conditional phase
//    (random Z on the partner) and, with probability `transport`, moves the leakage to the partner;
//  * a leaked qubit reads out as |1>; RESET removes leakage only if the reset is leakage-removing;
//  * seepage |2> -> |1> at ~2/T1 returns the qubit in a random computational state; DQLR does so on demand.
//
// Determinism: the result is a function of (M, shots, seed, opt) only, bit for bit, whatever OMP_NUM_THREADS or the
// schedule (the decode time decode_us_per_shot is a wall-clock measurement and is the one exception).
//  * every shot reseeds its generator from (seed, shot index) -- frame_shot_seed above;
//  * shots are grouped in fixed blocks of kFrameBlock; each block sums its own tallies in shot order, and the block
//    partials are added in block order after the parallel loop. The integer tallies (errors, detections per round,
//    defect counts) are exact in any order; the leaked-time sums (leak_d, leak_m) are not, so this fixed order is what
//    makes them reproducible.
inline constexpr long kFrameBlock = 1024;
inline MemoryResult run_memory_frame(const MemorySetup& M, long shots, std::uint64_t seed, const FrameRunOptions& opt = {}) {
    const auto& F = M.frame;
    if (!opt.decoder && !M.decoder) throw std::runtime_error("run_memory_frame: the setup has no decoder and none was given");
    const MatchingDecoder& decoder = opt.decoder ? *opt.decoder : *M.decoder;
    const DecodeMethod method = opt.method;
    int n_data = 0; for (char c : M.is_data) n_data += c;
    const int n_meas = F.n - n_data;
    const std::size_t n_layers = std::size_t(M.exp.rounds) + 1;
    std::vector<double> per_layer(n_layers, 0.0);
    for (int r : M.exp.detector_round) per_layer[std::size_t(r)] += 1;
    struct Tally {
        long errors{0}, fired{0}, leak_ev{0}, over22{0}, heralds{0}, heralds_true{0};
        long impacts{0}, impacts_vetoed{0}, errors_impacted{0}, errors_impacted_kept{0}, discarded{0}, errors_kept{0};
        int max_def{0}; double leak_d{0}, leak_m{0}, decode_s{0};
        std::vector<long> hist;
    };
    // impact injection: start time of each slice
    const ImpactSpec& IS = opt.impact;
    const bool inject = IS.p_shot > 0;
    std::vector<double> t_start(F.moments.size(), 0.0);
    for (std::size_t m = 1; m < F.moments.size(); ++m) t_start[m] = t_start[m - 1] + F.moments[m - 1].duration_ns;
    if (inject && M.sim_xy.size() != std::size_t(F.n)) throw std::runtime_error("run_memory_frame: impact injection needs qubit positions");
    // veto: the stabilisers within the window radius of each stabiliser of the memory basis
    const VetoSpec& VS = opt.veto;
    const int Rp1 = M.exp.rounds + 1;
    const std::size_t n_bs = M.exp.detectors.size() / std::size_t(Rp1);
    std::vector<std::vector<int>> vnbr(VS.on ? n_bs : 0);
    if (VS.on) {
        if (M.det_xy.size() != M.exp.detectors.size()) throw std::runtime_error("run_memory_frame: the veto needs detector positions");
        for (std::size_t c = 0; c < n_bs; ++c)
            for (std::size_t o = 0; o < n_bs; ++o) {
                const auto& a = M.det_xy[c * std::size_t(Rp1)]; const auto& b = M.det_xy[o * std::size_t(Rp1)];
                if (std::max(std::abs(a.first - b.first), std::abs(a.second - b.second)) <= VS.radius + 1e-9) vnbr[c].push_back(int(o));
            }
    }
    const HeraldSpec& HS = opt.herald;
    if (HS.on && HS.use_in_decoder && HS.footprint && (M.herald_edges.empty() || &decoder != M.decoder.get()))
        throw std::runtime_error("run_memory_frame: herald footprints need the setup's own decoder");
    std::vector<std::vector<int>> dets_of;                       // classical bit -> detectors (heralds)
    if (HS.on) {
        dets_of.assign(std::size_t(std::max(F.nclbits, 1)), {});
        for (std::size_t d = 0; d < M.exp.detectors.size(); ++d) for (int c : M.exp.detectors[d]) dets_of[std::size_t(c)].push_back(int(d));
    }
    const long nblocks = shots > 0 ? (shots + kFrameBlock - 1) / kFrameBlock : 0;
    std::vector<Tally> part(static_cast<std::size_t>(nblocks));
    std::vector<long> by_round(n_layers, 0);
#pragma omp parallel
    {
        Xoshiro256 rng(0);
        std::vector<std::uint8_t> x(std::size_t(F.n)), z(std::size_t(F.n)), flips(std::size_t(std::max(F.nclbits, 1))), L(std::size_t(F.n));
        std::vector<int> fired, heralded; std::vector<std::pair<int, int>> heralded_dets; std::vector<std::pair<int, double>> scaled; DecodeScratch ws;
        std::vector<long> my_round(n_layers, 0);
        Xoshiro256 irng(0);
        std::vector<double> g(inject ? std::size_t(F.n) : 0);
        std::vector<int> vcount(VS.on ? n_bs * std::size_t(Rp1) : 0, 0), vtouched;
        auto randomize = [&](int q) { x[std::size_t(q)] = rng.uniform() < 0.5; z[std::size_t(q)] = rng.uniform() < 0.5; };
#pragma omp for schedule(dynamic, 1)
        for (long b = 0; b < nblocks; ++b) {
            Tally& T_ = part[std::size_t(b)];
            const long s_end = std::min(shots, (b + 1) * kFrameBlock);
            for (long s = b * kFrameBlock; s < s_end; ++s) {
                rng = Xoshiro256(frame_shot_seed(seed, std::uint64_t(s)));
                std::fill(x.begin(), x.end(), 0); std::fill(z.begin(), z.end(), 0); std::fill(flips.begin(), flips.end(), 0);
                std::fill(L.begin(), L.end(), 0);
                heralded.clear();
                int nleaked = 0;
                bool hit = false; double t_hit = 0;
                if (inject) {
                    irng = Xoshiro256(impact_shot_seed(seed, std::uint64_t(s)));
                    hit = irng.uniform() < IS.p_shot;
                    if (hit) {
                        T_.impacts += 1;
                        t_hit = irng.uniform() * F.duration_ns;
                        const int c = std::min(F.n - 1, int(irng.uniform() * F.n));
                        const double r = IS.burst.radius;
                        for (int q = 0; q < F.n; ++q) {
                            const double dx = M.sim_xy[std::size_t(q)].first - M.sim_xy[std::size_t(c)].first;
                            const double dy = M.sim_xy[std::size_t(q)].second - M.sim_xy[std::size_t(c)].second;
                            g[std::size_t(q)] = r > 0 ? std::exp(-(dx * dx + dy * dy) / (2 * r * r)) : 1.0;
                        }
                    }
                }
                for (std::size_t mi = 0; mi < F.moments.size(); ++mi) {
                    const auto& fm = F.moments[mi];
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
                        if (HS.on && o.op == Op::MEASURE && !M.is_data[std::size_t(o.a)]) {
                            const bool leaked = L[std::size_t(o.a)] != 0;
                            if (rng.uniform() < (leaked ? HS.fidelity : HS.false_rate)) {
                                T_.heralds += 1; T_.heralds_true += leaked;
                                heralded.push_back(o.cbit);
                            }
                        }
                        if (F.leakage) {
                            if (o.op == Op::CZ && o.pleak > 0 && !L[std::size_t(o.mover)] && rng.uniform() < o.pleak) {
                                const int other = o.mover == o.a ? o.b : o.a;
                                L[std::size_t(o.mover)] = 1; ++nleaked; T_.leak_ev += 1;
                                if (!L[std::size_t(other)]) randomize(other);
                            } else if ((o.op == Op::SX || o.op == Op::X) && F.leak_1q > 0 && !L[std::size_t(o.a)] && rng.uniform() < F.leak_1q) {
                                L[std::size_t(o.a)] = 1; ++nleaked; T_.leak_ev += 1;
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
                    if (hit && fm.duration_ns > 0 && t_start[mi] + 0.5 * fm.duration_ns >= t_hit) {
                        // excess decay / dephasing of the impact, state-aware as in build_frame_program
                        const double tmid = t_start[mi] + 0.5 * fm.duration_ns, dur = fm.duration_ns;
                        const double decay = std::exp(-(tmid - t_hit) / (IS.burst.tau_ms * 1e6));
                        for (int q = 0; q < F.n; ++q) {
                            const double sv = decay * g[std::size_t(q)];
                            if (sv < 1e-9) continue;
                            const int kq = fm.known[std::size_t(q)];
                            if (kq == 0) continue;
                            const double r1 = IS.burst.t1_peak_us > 0 ? sv / (IS.burst.t1_peak_us * 1e3) : 0.0;
                            const double rphi = IS.burst.tphi_peak_us > 0 ? sv / (IS.burst.tphi_peak_us * 1e3) : 0.0;
                            const double gamma = 1 - std::exp(-dur * r1);
                            if (kq == 1) { if (irng.uniform() < gamma) x[std::size_t(q)] ^= 1; continue; }
                            const double px = gamma / 4, pz = std::max(0.0, (1 - std::exp(-dur * (0.5 * r1 + rphi))) / 2 - px);
                            const double u = irng.uniform();
                            if (u < px) x[std::size_t(q)] ^= 1;
                            else if (u < 2 * px) { x[std::size_t(q)] ^= 1; z[std::size_t(q)] ^= 1; }
                            else if (u < 2 * px + pz) z[std::size_t(q)] ^= 1;
                        }
                    }
                    if (F.leakage) {
                        for (int q : fm.dqlr) {
                            if (L[std::size_t(q)] && rng.uniform() < F.dqlr_eff) { L[std::size_t(q)] = 0; --nleaked; randomize(q); }
                            if (F.dqlr_err > 0 && rng.uniform() < F.dqlr_err) frame_detail::apply_pauli(1 + int(rng.uniform() * 3), q, x.data(), z.data());
                        }
                        if (nleaked > 0) {
                            int nd = 0, nmq = 0;
                            for (int q = 0; q < F.n; ++q) if (L[std::size_t(q)]) { if (M.is_data[std::size_t(q)]) ++nd; else ++nmq; }
                            T_.leak_d += fm.duration_ns * nd; T_.leak_m += fm.duration_ns * nmq;
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
                    if (v) { fired.push_back(int(d)); my_round[std::size_t(M.exp.detector_round[d])] += 1; }
                }
                std::uint8_t obs = 0; for (int c : M.exp.observable) obs ^= flips[std::size_t(c)];
                T_.fired += long(fired.size());
                T_.max_def = std::max(T_.max_def, int(fired.size()));
                if (fired.size() > 22) ++T_.over22;
                const auto t0 = std::chrono::steady_clock::now();
                int predicted;
                if (HS.on && HS.use_in_decoder && !heralded.empty()) {
                    scaled.clear();
                    if (HS.footprint) {
                        for (int c : heralded)
                            for (const auto& [e, ro] : M.herald_edges[std::size_t(c)]) scaled.push_back({e, ro ? HS.erase_factor : HS.footprint_factor});
                    } else {
                        heralded_dets.clear();
                        for (int c : heralded) {
                            const auto& dd = dets_of[std::size_t(c)];
                            if (!dd.empty()) heralded_dets.push_back({dd[0], dd.size() > 1 ? dd[1] : -1});
                        }
                        decoder.herald_touch_edges(heralded_dets, HS.erase_factor, HS.touch_factor, scaled);
                    }
                    predicted = decoder.decode_scaled(fired, scaled, ws, method);
                } else predicted = decoder.decode(fired, ws, method);
                T_.decode_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                bool vetoed = false;
                if (VS.on) {
                    // the burst statistic: the most detection events in one space-time window
                    int stat = 0;
                    for (int f : fired) {
                        const int sf = f / Rp1, rf = f % Rp1;
                        for (int c : vnbr[std::size_t(sf)])
                            for (int r = std::max(0, rf - VS.window_rounds + 1); r <= rf; ++r) {
                                int& v = vcount[std::size_t(c) * std::size_t(Rp1) + std::size_t(r)];
                                if (v == 0) vtouched.push_back(c * Rp1 + r);
                                stat = std::max(stat, ++v);
                            }
                    }
                    for (int k : vtouched) vcount[std::size_t(k)] = 0;
                    vtouched.clear();
                    if (T_.hist.size() <= std::size_t(stat)) T_.hist.resize(std::size_t(stat) + 1, 0);
                    T_.hist[std::size_t(stat)] += 1;
                    vetoed = stat >= VS.threshold;
                    if (vetoed) { T_.discarded += 1; if (hit) T_.impacts_vetoed += 1; }
                }
                if (predicted != obs) {
                    ++T_.errors;
                    if (!vetoed) ++T_.errors_kept;
                    if (hit) { ++T_.errors_impacted; if (!vetoed) ++T_.errors_impacted_kept; }
                }
            }
        }
#pragma omp critical(qpuforge_frame_rounds)
        for (std::size_t r = 0; r < n_layers; ++r) by_round[r] += my_round[r];   // integer counts: exact in any order
    }
    long errors = 0, fired_tot = 0, leak_ev = 0, over22 = 0, heralds = 0, heralds_true = 0; int max_def = 0; double leak_d = 0, leak_m = 0, decode_s = 0;
    for (const Tally& t : part) {
        errors += t.errors; fired_tot += t.fired; leak_ev += t.leak_ev; over22 += t.over22; max_def = std::max(max_def, t.max_def);
        heralds += t.heralds; heralds_true += t.heralds_true;
        leak_d += t.leak_d; leak_m += t.leak_m; decode_s += t.decode_s;
    }
    long impacts = 0, impacts_vetoed = 0, errors_impacted = 0, errors_impacted_kept = 0, discarded = 0, errors_kept = 0;
    std::vector<long> hist;
    for (const Tally& t : part) {
        impacts += t.impacts; impacts_vetoed += t.impacts_vetoed; errors_impacted += t.errors_impacted; errors_impacted_kept += t.errors_impacted_kept;
        discarded += t.discarded; errors_kept += t.errors_kept;
        if (hist.size() < t.hist.size()) hist.resize(t.hist.size(), 0);
        for (std::size_t k = 0; k < t.hist.size(); ++k) hist[k] += t.hist[k];
    }
    auto r = summarize(M, shots, errors, double(fired_tot));
    r.detection_by_round.resize(n_layers);
    for (std::size_t k = 0; k < n_layers; ++k)
        r.detection_by_round[k] = per_layer[k] > 0 ? double(by_round[k]) / (per_layer[k] * double(shots)) : 0.0;
    r.mean_defects = double(fired_tot) / double(shots);
    r.max_defects = max_def;
    r.frac_over_22 = double(over22) / double(shots);
    r.decode_us_per_shot = 1e6 * decode_s / double(shots);
    const double T = std::max(M.frame.duration_ns, 1e-12) * double(shots);
    r.leak_frac_data = n_data ? leak_d / T / n_data : 0.0;
    r.leak_frac_meas = n_meas ? leak_m / T / n_meas : 0.0;
    r.leak_events_per_shot = double(leak_ev) / double(shots);
    r.heralds_per_shot = double(heralds) / double(shots);
    r.herald_true_frac = heralds ? double(heralds_true) / double(heralds) : 0.0;
    r.impacts = impacts; r.impacts_vetoed = impacts_vetoed; r.errors_impacted = errors_impacted; r.errors_impacted_kept = errors_impacted_kept;
    if (!VS.on) errors_kept = errors;
    r.discarded = discarded; r.errors_kept = errors_kept;
    r.p_logical_kept = shots > discarded ? double(errors_kept) / double(shots - discarded) : 0.0;
    r.eps_kept = eps_from_PL(r.p_logical_kept, r.rounds);
    r.veto_hist = std::move(hist);
    return r;
}

// Exact-noise cross-check: statevector trajectories with mid-circuit measurement (one per shot).
inline MemoryResult run_memory_statevector(const MemorySetup& M, const NoiseModel& nm, long shots, unsigned seed, bool noisy = true) {
    const auto recs = run_records(M.prog, nm, int(shots), seed, noisy);
    long errors = 0; double fired_sum = 0;
    std::vector<int> fired; DecodeScratch ws;
    for (const auto& rc : recs) {
        fired.clear();
        for (std::size_t d = 0; d < M.exp.detectors.size(); ++d) {
            std::uint8_t v = 0; for (int c : M.exp.detectors[d]) v ^= rc[std::size_t(c)];
            if (v) fired.push_back(int(d));
        }
        std::uint8_t obs = 0; for (int c : M.exp.observable) obs ^= rc[std::size_t(c)];
        fired_sum += double(fired.size());
        if (M.decoder->decode(fired, ws) != obs) ++errors;
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

// One-sigma Monte Carlo error of fit_eps_per_round: the binomial errors of the points, propagated through the same
// least-squares fit. (The slope of a free-intercept fit is less certain than the last point alone would suggest: by
// about 20% for five round counts and by 40% for four.) With fewer than two usable points it is the error of the
// single-run estimate eps_from_PL.
inline double fit_eps_per_round_err(const std::vector<MemoryResult>& rs) {
    std::vector<double> x, y, var;
    for (const auto& r : rs) {
        if (r.p_logical >= 0.5 || r.p_logical <= 0) continue;
        const double sy = 2 * r.p_logical_sem / (1 - 2 * r.p_logical);
        x.push_back(r.rounds); y.push_back(std::log(1 - 2 * r.p_logical)); var.push_back(sy * sy);
    }
    if (x.size() < 2) {
        if (rs.empty()) return 0.0;
        const auto& r = rs.back();
        return r.p_logical_sem / r.rounds / std::max(1e-9, std::pow(std::max(1e-9, 1 - 2 * r.p_logical), double(r.rounds - 1) / r.rounds));
    }
    double xm = 0; for (double v : x) xm += v;
    xm /= double(x.size());
    double sxx = 0, sxy = 0, num = 0;
    for (std::size_t i = 0; i < x.size(); ++i) { const double d = x[i] - xm; sxx += d * d; sxy += d * y[i]; num += d * d * var[i]; }
    return 0.5 * std::exp(sxy / sxx) * std::sqrt(num) / sxx;
}

// ================================================================== correlated error bursts (high-energy impacts)
// (BurstSpec is defined above, with the frame-run options.)
inline NoiseModel burst_noise(const NoiseModel& nm, const BurstSpec& B, double s) {
    NoiseModel N = nm;
    if (B.t1_peak_us > 0)
        for (auto& t : N.T1_ns) { const double base = std::isfinite(t) ? 1.0 / t : 0.0; t = 1.0 / (base + s / (B.t1_peak_us * 1e3)); }
    if (B.tphi_peak_us > 0)
        for (auto& t : N.Tphi_ns) { const double base = std::isfinite(t) ? 1.0 / t : 0.0; t = 1.0 / (base + s / (B.tphi_peak_us * 1e3)); }
    return N;
}

struct BurstPoint { double s{0}, eps{0}, eps_sem{0}, mean_defects{0}; bool saturated{false}; };   // saturated: P_L near 1/2, eps set to 1/2
struct BurstResult {
    std::vector<BurstPoint> curve;     // logical error per round against severity
    double eps0{0}, eps0_sem{0}, round_ns{0}, mean_defects0{0};
    double rounds_lost{0};             // integral of the excess logical error over one burst (expected logical flips)
    double p_fail{0};                  // probability that one impact flips the logical qubit
    double floor_per_round{0};         // long-run logical error per round from impacts alone
    double dead_time_ms{0};            // time after an impact for which the excess logical error exceeds the quiet value
    double detectable_ms{0};           // time for which the detection-event count is more than 1.5x the quiet value
};

inline BurstResult burst_floor(const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, const BurstSpec& B,
                               int rounds, long shots, std::uint64_t seed, char basis = 'Z') {
    BurstResult out;
    const MemorySetup M0 = prepare_memory(D, R, nm, rounds, basis);
    const MemoryResult r0 = run_memory_frame(M0, shots, seed);
    out.eps0 = r0.eps_per_round; out.round_ns = r0.round_ns; out.mean_defects0 = r0.mean_defects;
    // standard error of eps; a saturated run (P_L near 1/2) only says "eps is large", which is all the integral needs
    auto sem_of = [](const MemoryResult& r) {
        if (r.p_logical > 0.45) return 0.05;
        return r.p_logical_sem / r.rounds / std::pow(1 - 2 * r.p_logical, double(r.rounds - 1) / r.rounds);
    };
    out.eps0_sem = sem_of(r0);
    FrameRunOptions opt; opt.decoder = M0.decoder.get();
    double s = 1.0, last_eps = 0.5;
    for (int k = 0; k < 16; ++k, s /= std::sqrt(10.0)) {
        // near the peak nearly every detector fires: decoding is slow and the answer ("eps is large") needs few shots
        const long n = last_eps > 0.02 ? std::max<long>(1000, shots / 5) : shots;
        const MemorySetup Ms = prepare_memory(D, R, burst_noise(nm, B, s), rounds, basis, false);
        const MemoryResult rs = run_memory_frame(Ms, n, seed + 17 * std::uint64_t(k + 1), opt);
        // A run whose failure probability is close to 1/2 cannot resolve eps: it only says the logical qubit is
        // randomised within the run. Such a point is recorded as eps = 1/2 and never ends the ladder.
        const bool saturated = rs.p_logical > 0.45;
        const double eps = saturated ? 0.5 : rs.eps_per_round;
        out.curve.push_back({s, eps, sem_of(rs), rs.mean_defects, saturated});
        last_eps = eps;
        if (!saturated && eps - out.eps0 < 2 * std::hypot(sem_of(rs), out.eps0_sem)) break;   // lost in the noise from here on
    }
    // excess(s) = eps(s) - eps0 is integrated over the burst: dt = -tau ds / s
    auto excess = [&](std::size_t i) { return std::max(0.0, out.curve[i].eps - out.eps0); };
    double integral = excess(out.curve.size() - 1);          // tail below the last point: excess ~ linear in s
    for (std::size_t i = 0; i + 1 < out.curve.size(); ++i)
        integral += 0.5 * (excess(i) + excess(i + 1)) * std::log(out.curve[i].s / out.curve[i + 1].s);
    const double tau_ns = B.tau_ms * 1e6;
    out.rounds_lost = integral * tau_ns / out.round_ns;
    out.p_fail = 0.5 * (1 - std::exp(-2 * out.rounds_lost));
    out.floor_per_round = B.rate_hz * out.round_ns * 1e-9 * out.p_fail;
    // crossing times: log-log interpolation between ladder points (extrapolated linearly in s below the ladder)
    auto crossing = [&](const std::vector<double>& y, double level) {
        if (y.empty() || y[0] <= level) return 0.0;
        for (std::size_t i = 0; i + 1 < y.size(); ++i)
            if (y[i + 1] <= level) {
                const double a = std::log(std::max(y[i], 1e-300)), b = std::log(std::max(y[i + 1], 1e-300)), l = std::log(level);
                const double f = (a - l) / std::max(a - b, 1e-12);
                return std::log(1.0 / out.curve[i].s) + f * std::log(out.curve[i].s / out.curve[i + 1].s);
            }
        return std::log(1.0 / out.curve.back().s) + std::log(std::max(y.back(), 1e-300) / level);
    };
    std::vector<double> ex, dx;
    for (std::size_t i = 0; i < out.curve.size(); ++i) { ex.push_back(excess(i)); dx.push_back(std::max(0.0, out.curve[i].mean_defects - out.mean_defects0)); }
    out.dead_time_ms = B.tau_ms * crossing(ex, std::max(out.eps0, 1e-12));
    out.detectable_ms = B.tau_ms * crossing(dx, 0.5 * std::max(out.mean_defects0, 1e-9));
    return out;
}

// ================================================================== impact veto
// Threshold for a target false-veto rate: the burst statistic of `shots` burst-free shots (same options otherwise),
// and the smallest threshold whose tail fraction is at most `target`.
struct VetoCalibration { int threshold{0}; double false_rate{0}; long shots{0}; std::vector<long> hist; };
inline VetoCalibration calibrate_veto(const MemorySetup& M, long shots, std::uint64_t seed, double target, const VetoSpec& spec,
                                      const FrameRunOptions& base = {}) {
    FrameRunOptions o = base;
    o.impact.p_shot = 0;
    o.veto = spec; o.veto.on = true; o.veto.threshold = std::numeric_limits<int>::max();
    const auto r = run_memory_frame(M, shots, seed, o);
    VetoCalibration c; c.shots = shots; c.hist = r.veto_hist;
    long tail = 0;
    int T = int(c.hist.size());                                   // above the largest value seen: nothing is vetoed
    for (int k = int(c.hist.size()) - 1; k >= 0; --k) {
        if (double(tail + c.hist[std::size_t(k)]) > target * double(shots)) break;
        tail += c.hist[std::size_t(k)]; T = k;
    }
    c.threshold = std::max(T, 1);
    long t2 = 0; for (std::size_t k = std::size_t(c.threshold); k < c.hist.size(); ++k) t2 += c.hist[k];
    c.false_rate = double(t2) / double(std::max(shots, 1L));
    return c;
}

// The impact floor with a veto (docs/design/impact-veto.md): burst_floor's ladder of chip-wide severities, each run
// with the veto on. A vetoed shot is discarded, so a burst costs logical errors only through the shots the veto keeps:
//   rounds_lost_veto = integral over the burst of (1 - v(s)) * max(0, eps_kept(s) - eps_kept(0)) dt / t_round
// and the veto costs the time it discards, integral of (v(s) - v(0)) dt (v: discarded fraction). The ladder stops
// where neither the excess error nor the excess discard rate is resolved any more.
// A latched veto holds the discard for `hold` after it first fires at the onset; only the tail after the hold costs
// logical errors: rounds_lost(hold) = integral over t > hold of max(0, eps(s) - eps0) dt / t_round (no-veto curve,
// interpolated linearly in log s; below the ladder the excess is taken as linear in s, as in burst_floor), and the
// time lost is rate x hold. Failures before the veto first fires are measured separately with impact injection.
struct BurstVetoPoint { double s{0}, eps{0}, eps_kept{0}, veto_frac{0}, mean_defects{0}; long kept{0}; bool saturated{false}; };
struct BurstVetoResult {
    std::vector<BurstVetoPoint> curve;
    double eps0{0}, eps0_kept{0}, veto0{0}, round_ns{0};
    double rounds_lost{0}, p_fail{0}, floor_per_round{0};                // without the veto (as burst_floor)
    double rounds_lost_veto{0}, p_fail_veto{0}, floor_per_round_veto{0}; // with it
    double discarded_ms{0};            // burst time discarded per impact, beyond the quiet false-veto rate
    double discarded_time_frac{0};     // of the chip's time (rate x discarded_ms)
    struct Latched { double hold_tau{0}, rounds_lost{0}, p_fail{0}, floor_per_round{0}, dead_time_frac{0}; };
    std::vector<Latched> latched;      // by hold time (units of tau)
};
inline BurstVetoResult burst_floor_veto(const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, const BurstSpec& B, int rounds,
                                        long shots, std::uint64_t seed, const VetoSpec& veto, char basis = 'Z', DecodeMethod method = DecodeMethod::Auto,
                                        const std::vector<double>& holds_tau = {1, 2, 3, 5, 7, 10, 15, 20}) {
    BurstVetoResult out;
    const MemorySetup M0 = prepare_memory(D, R, nm, rounds, basis);
    FrameRunOptions opt; opt.decoder = M0.decoder.get(); opt.veto = veto; opt.veto.on = true; opt.method = method;
    const MemoryResult r0 = run_memory_frame(M0, shots, seed, opt);
    out.eps0 = r0.eps_per_round; out.eps0_kept = r0.eps_kept; out.round_ns = r0.round_ns;
    out.veto0 = double(r0.discarded) / double(shots);
    auto sem_of = [](double p, long n, int rr) {
        if (p > 0.45) return 0.05;
        const double sem = std::sqrt(std::max(p * (1 - p), 1e-12) / double(std::max(n, 1L)));
        return sem / rr / std::pow(1 - 2 * p, double(rr - 1) / rr);
    };
    const double sem0 = sem_of(r0.p_logical, shots, rounds), sem0k = sem_of(r0.p_logical_kept, shots - r0.discarded, rounds);
    double s = 1.0, last_eps = 0.5;
    for (int k = 0; k < 16; ++k, s /= std::sqrt(10.0)) {
        // as in burst_floor: near the peak decoding is slow and a few shots resolve "eps is large"
        const long n = last_eps > 0.02 ? std::max<long>(1000, shots / 5) : shots;
        const MemorySetup Ms = prepare_memory(D, R, burst_noise(nm, B, s), rounds, basis, false);
        const MemoryResult rs = run_memory_frame(Ms, n, seed + 17 * std::uint64_t(k + 1), opt);
        BurstVetoPoint pt;
        pt.s = s; pt.saturated = rs.p_logical > 0.45; pt.eps = pt.saturated ? 0.5 : rs.eps_per_round;
        pt.kept = rs.shots - rs.discarded; pt.veto_frac = double(rs.discarded) / double(rs.shots);
        pt.eps_kept = pt.kept > 0 ? (rs.p_logical_kept > 0.45 ? 0.5 : rs.eps_kept) : 0.0;
        pt.mean_defects = rs.mean_defects;
        out.curve.push_back(pt);
        last_eps = pt.eps;
        const double sv = std::sqrt(std::max(pt.veto_frac * (1 - pt.veto_frac), 1e-12) / double(n));
        const bool err_lost = !pt.saturated && pt.eps - out.eps0 < 2 * std::hypot(sem_of(rs.p_logical, n, rounds), sem0);
        const bool kept_lost = pt.kept == 0 || (1 - pt.veto_frac) * (pt.eps_kept - out.eps0_kept) < 2 * std::hypot(sem_of(rs.p_logical_kept, pt.kept, rounds), sem0k);
        const bool veto_lost = pt.veto_frac - out.veto0 < 2 * sv;
        if (err_lost && kept_lost && veto_lost) break;
    }
    auto integrate = [&](auto f) {
        double I = f(out.curve.size() - 1);                       // tail below the last point: linear in s
        for (std::size_t i = 0; i + 1 < out.curve.size(); ++i) I += 0.5 * (f(i) + f(i + 1)) * std::log(out.curve[i].s / out.curve[i + 1].s);
        return I;
    };
    const double tau_ns = B.tau_ms * 1e6;
    out.rounds_lost = integrate([&](std::size_t i) { return std::max(0.0, out.curve[i].eps - out.eps0); }) * tau_ns / out.round_ns;
    out.rounds_lost_veto = integrate([&](std::size_t i) {
        const auto& c = out.curve[i];
        return c.kept > 0 ? (1 - c.veto_frac) * std::max(0.0, c.eps_kept - out.eps0_kept) : 0.0;
    }) * tau_ns / out.round_ns;
    out.p_fail = 0.5 * (1 - std::exp(-2 * out.rounds_lost));
    out.p_fail_veto = 0.5 * (1 - std::exp(-2 * out.rounds_lost_veto));
    out.floor_per_round = B.rate_hz * out.round_ns * 1e-9 * out.p_fail;
    out.floor_per_round_veto = B.rate_hz * out.round_ns * 1e-9 * out.p_fail_veto;
    out.discarded_ms = integrate([&](std::size_t i) { return std::max(0.0, out.curve[i].veto_frac - out.veto0); }) * B.tau_ms;
    out.discarded_time_frac = B.rate_hz * out.discarded_ms * 1e-3;
    // latched hold: integral of the excess below s_h = exp(-hold / tau)
    auto ex = [&](std::size_t i) { return std::max(0.0, out.curve[i].eps - out.eps0); };
    auto below = [&](double sh) {
        const auto& c = out.curve;
        const std::size_t n = c.size();
        if (sh <= c[n - 1].s) return ex(n - 1) * sh / c[n - 1].s;   // linear tail
        double I = ex(n - 1);
        for (std::size_t i = n - 1; i-- > 0;) {
            const double L0 = std::log(c[i].s), L1 = std::log(c[i + 1].s);
            if (sh >= c[i].s) { I += 0.5 * (ex(i) + ex(i + 1)) * (L0 - L1); continue; }
            const double Lh = std::log(sh), f = (Lh - L1) / (L0 - L1);
            I += 0.5 * (ex(i + 1) + (ex(i + 1) + f * (ex(i) - ex(i + 1)))) * (Lh - L1);
            break;
        }
        return I;
    };
    for (double h : holds_tau) {
        BurstVetoResult::Latched L;
        L.hold_tau = h;
        L.rounds_lost = below(std::exp(-h)) * tau_ns / out.round_ns;
        L.p_fail = 0.5 * (1 - std::exp(-2 * L.rounds_lost));
        L.floor_per_round = B.rate_hz * out.round_ns * 1e-9 * L.p_fail;
        L.dead_time_frac = B.rate_hz * h * B.tau_ms * 1e-3;
        out.latched.push_back(L);
    }
    return out;
}

} // namespace qpuforge
