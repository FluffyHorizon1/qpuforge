// qpuforge — whole-processor emulation of a chip design
//
// Every moment of the compiled schedule:
//   1. ideal native ops (virtual RZ, SX/X pulses, CZ)
//   2. control errors (Pauli channels sized to the calibration-level infidelity)
//   3. always-on ZZ crosstalk on every coupled pair that is not executing a CZ (conditional phase)
//   4. relaxation towards the thermal state (T1, thermal excitation) and pure dephasing (T_phi)
//      on every simulated qubit for the duration of the moment — idling qubits decohere too
// Initial state: thermal excited-state population (or reset error). Readout: assignment errors.
#pragma once
#include "compiler.hpp"
#include "device.hpp"
#include "sim.hpp"
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace qpuforge {

struct NoiseModel {
    std::vector<double> T1_ns, Tphi_ns, p_thermal, P1_init, ro01, ro10;
    double err_1q_ctrl{0}, err_2q_ctrl{0};
    double reset_err{0};                              // probability a RESET leaves the qubit in |1>
    std::map<std::pair<int, int>, double> zz_GHz;     // physical pair (min,max) -> conditional-phase rate
    // near-degenerate pairs (|detuning| < 10 g): resonant flip-flop exchange is simulated explicitly
    std::map<std::pair<int, int>, std::array<double, 3>> exchange;   // (min,max) -> {g, f_min, f_max} GHz
    std::map<std::pair<int, int>, double> err_2q_pair;               // (min,max) -> per-pair control error (overrides err_2q_ctrl)
    std::map<std::pair<int, int>, std::pair<int, double>> leak_cz;    // (min,max) -> (qubit that leaks, probability per CZ)
    double err_2q(int a, int b) const {
        auto it = err_2q_pair.find({std::min(a, b), std::max(a, b)});
        return it == err_2q_pair.end() ? err_2q_ctrl : it->second;
    }
    bool T1{true}, Tphi{true}, thermal{true}, zz{true}, exchange_on{true}, control{true}, readout{true}, reset{true};
    // leakage (Pauli-frame QEC model only; the state-vector emulator is qubit-only)
    bool leakage{false};
    double leak_1q{0}, leak_transport{0}, dqlr_eff{0}, dqlr_err{0};
    bool reset_removes_leakage{false}, dqlr{false};

    // switch a named source off ("T1", "Tphi", "thermal", "zz", "control", "readout", "all")
    void disable(const std::string& s) {
        if (s == "T1" || s == "all") T1 = false;
        if (s == "Tphi" || s == "all") Tphi = false;
        if (s == "thermal" || s == "all") thermal = false;
        if (s == "zz" || s == "all") zz = false;
        if (s == "exchange" || s == "all") exchange_on = false;
        if (s == "control" || s == "all") control = false;
        if (s == "readout" || s == "all") readout = false;
        if (s == "reset" || s == "all") reset = false;
        if (s == "leakage" || s == "all") leakage = false;
    }
    void only(const std::string& s) { disable("all"); if (s == "T1") T1 = true; else if (s == "Tphi") Tphi = true; else if (s == "thermal") thermal = true; else if (s == "zz") zz = true; else if (s == "exchange") exchange_on = true; else if (s == "control") control = true; else if (s == "readout") readout = true; else if (s == "reset") reset = true; else if (s == "leakage") leakage = true; }
};

inline NoiseModel noise_model(const ChipDesign& D, const ChipReport& R) {
    NoiseModel N;
    for (const auto& q : R.q) {
        N.T1_ns.push_back(q.T1_us * 1e3);
        N.Tphi_ns.push_back(q.Tphi_us * 1e3);
        N.p_thermal.push_back(q.p_thermal);
        N.P1_init.push_back(D.active_reset ? D.reset_err : q.p_thermal);
        N.ro01.push_back(q.ro_err_0to1);
        N.ro10.push_back(q.ro_err_1to0);
    }
    N.err_1q_ctrl = D.err_1q_control;
    N.err_2q_ctrl = D.err_2q_control;
    N.reset_err = D.reset_err;
    N.leakage = D.leakage.model;
    N.leak_1q = D.leakage.leak_1q; N.leak_transport = D.leakage.transport;
    N.reset_removes_leakage = D.leakage.reset_removes_leakage;
    N.dqlr = D.leakage.dqlr; N.dqlr_eff = D.leakage.dqlr_eff; N.dqlr_err = D.leakage.dqlr_err;
    for (const auto& p : R.pairs) {
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        N.zz_GHz[{std::min(a, b), std::max(a, b)}] = p.zz_idle_kHz * 1e-6;
        if (p.designed && p.err_2q_ctrl >= 0) N.err_2q_pair[{std::min(a, b), std::max(a, b)}] = p.err_2q_ctrl;
        if (p.designed && p.leak_2q > 0 && p.leak_qubit >= 0) N.leak_cz[{std::min(a, b), std::max(a, b)}] = {D.index_of(p.leak_qubit), p.leak_2q};
        const double g = p.g_MHz * 1e-3;
        const bool coupler_off = p.designed && D.tunable_coupler();   // tunable couplers idle at zero coupling
        if (!coupler_off && std::abs(p.delta_MHz) * 1e-3 < 10 * g) {
            const int lo = std::min(a, b), hi = std::max(a, b);
            N.exchange[{lo, hi}] = {g, R.q[lo].f01, R.q[hi].f01};
        }
    }
    return N;
}

inline ChipGraph chip_graph(const ChipDesign& D, const ChipReport& R) {
    ChipGraph G;
    G.n = int(D.qubits.size());
    G.adj.assign(G.n, {});
    G.t_1q_ns = D.t_1q_ns; G.t_ro_ns = D.t_readout_ns; G.t_reset_ns = D.t_reset_ns + (D.leakage.model && D.leakage.dqlr ? D.leakage.dqlr_ns : 0.0);
    for (const auto& p : R.pairs) {
        if (!p.designed) continue;
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        G.adj[a].push_back(b); G.adj[b].push_back(a);
        G.t_cz_ns[{std::min(a, b), std::max(a, b)}] = p.t_cz_ns;
    }
    G.finalize();
    return G;
}

// Execute the schedule on a backend (StateVector with sampling, or DensityMatrix exactly).
// sim[p] = simulation index of physical qubit p, or -1 if it is not simulated.
template <class B, class RNG>
void execute(B& S, const Compiled& prog, const std::vector<int>& sim, const NoiseModel& nm, bool noisy, RNG& rng,
             std::vector<std::uint8_t>* record = nullptr) {
    std::uniform_real_distribution<double> U01(0, 1);
    const double p1 = 1.5 * nm.err_1q_ctrl, p2 = 1.25 * nm.err_2q_ctrl;   // Pauli prob from average infidelity
    std::vector<std::pair<std::pair<int, int>, double>> zz;                 // (sim pair, rate)
    for (const auto& [pr, rate] : nm.zz_GHz) {
        const int a = sim[pr.first], b = sim[pr.second];
        if (a >= 0 && b >= 0 && rate != 0) zz.push_back({{a, b}, rate});
    }
    struct Ex { int a, b; double g, fa, fb; };
    std::vector<Ex> ex;
    for (const auto& [pr, v] : nm.exchange) {
        const int a = sim[pr.first], b = sim[pr.second];
        if (a >= 0 && b >= 0) ex.push_back({a, b, v[0], v[1], v[2]});
    }
    std::vector<int> act;
    for (std::size_t p = 0; p < sim.size(); ++p) if (sim[p] >= 0) act.push_back(int(p));
    double t_ns = 0;
    for (const auto& m : prog.moments) {
        const double t_start = t_ns;
        t_ns += m.duration_ns;
        std::vector<std::pair<int, int>> busy;
        for (auto [a, b] : m.busy) busy.push_back({std::min(sim[a], sim[b]), std::max(sim[a], sim[b])});
        for (const auto& g : m.ops) {
            switch (g.op) {
            case Op::RZ: case Op::SX: case Op::X: S.gate1(sim[g.q[0]], g.u); break;
            case Op::CZ: {
                S.cz(sim[g.q[0]], sim[g.q[1]]);
                const std::pair<int, int> key{std::min(sim[g.q[0]], sim[g.q[1]]), std::max(sim[g.q[0]], sim[g.q[1]])};
                if (std::find(busy.begin(), busy.end(), key) == busy.end()) busy.push_back(key);
                break;
            }
            case Op::MEASURE: {
                int bit = S.measure(sim[g.q[0]], rng);
                if (record) {
                    if (bit < 0) throw std::runtime_error("mid-circuit measurement records need the statevector backend");
                    if (noisy && nm.readout) {
                        const int p = g.q[0];
                        if (bit == 0 && U01(rng) < nm.ro01[p]) bit = 1;
                        else if (bit == 1 && U01(rng) < nm.ro10[p]) bit = 0;
                    }
                    (*record)[std::size_t(g.cbit)] = std::uint8_t(bit);
                }
                break;
            }
            case Op::RESET:
                S.reset(sim[g.q[0]], rng);
                if (noisy && nm.reset) S.bitflip(sim[g.q[0]], nm.reset_err, rng);
                break;
            default: throw std::runtime_error("execute: non-native op in schedule");
            }
        }
        if (!noisy || m.duration_ns <= 0) continue;
        const double tau = m.duration_ns;
        if (nm.control)
            for (const auto& g : m.ops) {
                if (g.op == Op::SX || g.op == Op::X) S.pauli1(sim[g.q[0]], p1, rng);
                else if (g.op == Op::CZ) S.pauli2(sim[g.q[0]], sim[g.q[1]], nm.err_2q_pair.empty() ? p2 : 1.25 * nm.err_2q(g.q[0], g.q[1]), rng);
            }
        if (nm.zz) {
            std::vector<std::pair<std::size_t, cplx>> terms;
            for (const auto& [pr, rate] : zz) {
                const std::pair<int, int> key{std::min(pr.first, pr.second), std::max(pr.first, pr.second)};
                if (std::find(busy.begin(), busy.end(), key) != busy.end()) continue;   // ZZ is part of the calibrated CZ
                terms.push_back({(std::size_t(1) << pr.first) | (std::size_t(1) << pr.second), std::exp(cplx(0, -2 * PI * rate * tau))});
            }
            S.multi_cphase(terms);
        }
        if (nm.exchange_on)
            for (const auto& e : ex) {
                const std::pair<int, int> key{std::min(e.a, e.b), std::max(e.a, e.b)};
                if (std::find(busy.begin(), busy.end(), key) != busy.end()) continue;
                // rotating-frame propagator of H = f_a n_a + f_b n_b + g (s+ s- + h.c.) on {|10>,|01>} over [t, t+tau]
                const double d = 0.5 * (e.fa - e.fb), Om = std::sqrt(d * d + e.g * e.g);
                const double c = std::cos(2 * PI * Om * tau), s = std::sin(2 * PI * Om * tau);
                const cplx I(0, 1);
                M2 U;
                U.a = std::exp(I * (2 * PI * d * tau)) * (c - I * (d / Om) * s);
                U.d = std::exp(-I * (2 * PI * d * tau)) * (c + I * (d / Om) * s);
                U.b = std::exp(I * (2 * PI * (2 * d * t_start + d * tau))) * (-I * (e.g / Om) * s);
                U.c = std::exp(-I * (2 * PI * (2 * d * t_start + d * tau))) * (-I * (e.g / Om) * s);
                S.exch(e.a, e.b, U);
            }
        for (int p : act) {
            const int q = sim[p];
            if (nm.T1 && std::isfinite(nm.T1_ns[p])) {
                const double gamma = 1 - std::exp(-tau / nm.T1_ns[p]);
                const double pdown = nm.thermal ? 1 - nm.p_thermal[p] : 1.0;
                S.gad(q, gamma, pdown, rng);
            }
            if (nm.Tphi && std::isfinite(nm.Tphi_ns[p])) S.dephase(q, 0.5 * (1 - std::exp(-tau / nm.Tphi_ns[p])), rng);
        }
    }
}

struct RunOptions {
    int trajectories{200};
    int shots_per_traj{20};
    unsigned seed{1234};
    bool noisy{true};
    bool xeb{false};
};

struct RunResult {
    int n_sim{0}, nclbits{0};
    double fidelity{1}, fidelity_sem{0};    // <psi_ideal| rho |psi_ideal> before measurement
    double xeb{std::numeric_limits<double>::quiet_NaN()};       // linear XEB: 2^m <p_ideal(x)> - 1
    double xeb_ideal{std::numeric_limits<double>::quiet_NaN()}; // same quantity for a perfect device: 2^m sum p^2 - 1
    double xeb_norm{std::numeric_limits<double>::quiet_NaN()};  // xeb / xeb_ideal (fidelity estimator for small n)
    double p_ideal_top{0};                  // probability of the most likely ideal outcome, measured noisily
    std::string ideal_top;
    long shots{0};
    std::map<std::string, long> counts;
    double duration_ns{0};
    int n_cz{0}, n_pulses{0}, n_swaps{0}, depth{0};
};

namespace emu_detail {
inline std::string key_to_bits(std::uint64_t key, int nc) {   // Qiskit order: c[nc-1] ... c[0]
    std::string s(std::size_t(nc), '0');
    for (int c = 0; c < nc; ++c) if ((key >> c) & 1ULL) s[std::size_t(nc - 1 - c)] = '1';
    return s;
}
} // namespace emu_detail

inline std::vector<int> sim_index(const Compiled& prog) {
    std::vector<int> sim(prog.n_phys, -1);
    for (std::size_t k = 0; k < prog.active.size(); ++k) sim[prog.active[k]] = int(k);
    return sim;
}

inline RunResult run(const Compiled& prog, const NoiseModel& nm, const RunOptions& opt = {}) {
    RunResult R;
    const auto sim = sim_index(prog);
    const int n = int(prog.active.size());
    if (n > 28) throw std::runtime_error("too many active qubits for the statevector backend (" + std::to_string(n) + ")");
    R.n_sim = n; R.nclbits = prog.nclbits; R.duration_ns = prog.duration_ns;
    R.n_cz = prog.n_cz; R.n_pulses = prog.n_pulses; R.n_swaps = prog.n_swaps; R.depth = int(prog.moments.size());

    std::mt19937_64 rng0(opt.seed);
    StateVector ideal(n);
    execute(ideal, prog, sim, nm, false, rng0);

    // measured bits: (sim qubit, cbit)
    std::vector<std::pair<int, int>> meas;
    for (auto [p, c] : prog.measures) meas.push_back({sim[p], c});
    const int nc = prog.nclbits;
    auto cbits_of = [&](std::size_t idx) { std::uint64_t k = 0; for (auto [q, c] : meas) if ((idx >> q) & 1) k |= (1ULL << c); return k; };

    // ideal distribution over classical bits
    std::unordered_map<std::uint64_t, double> pideal;
    for (std::size_t i = 0; i < ideal.dim(); ++i) { const double p = std::norm(ideal.amps()[i]); if (p > 1e-15) pideal[cbits_of(i)] += p; }
    double ptop = -1; std::uint64_t ktop = 0;
    for (auto& [k, p] : pideal) if (p > ptop) { ptop = p; ktop = k; }
    R.ideal_top = emu_detail::key_to_bits(ktop, nc);

    if (prog.dynamic) {
        // each trajectory is one shot; MEASURE ops write the classical register during execution
        const int T = std::max(1, opt.trajectories);
        std::vector<std::vector<std::uint8_t>> recs(T, std::vector<std::uint8_t>(std::size_t(std::max(nc, 1)), 0));
#pragma omp parallel for schedule(dynamic, 1)
        for (int t = 0; t < T; ++t) {
            std::mt19937_64 rng(opt.seed * 1000003ULL + std::uint64_t(t) * 7919ULL + 17ULL);
            std::uniform_real_distribution<double> U(0, 1);
            StateVector psi(n);
            std::uint64_t init = 0;
            if (opt.noisy && nm.thermal) for (int k = 0; k < n; ++k) if (U(rng) < nm.P1_init[prog.active[k]]) init |= (1ULL << k);
            psi.set_basis(init);
            execute(psi, prog, sim, nm, opt.noisy, rng, &recs[t]);
        }
        for (const auto& rc : recs) {
            std::uint64_t key = 0;
            for (int c = 0; c < nc && c < 64; ++c) if (rc[std::size_t(c)]) key |= (1ULL << c);
            ++R.counts[emu_detail::key_to_bits(key, nc)];
        }
        R.shots = T;
        R.fidelity = std::numeric_limits<double>::quiet_NaN();
        return R;
    }
    const int T = opt.noisy ? std::max(1, opt.trajectories) : 1;
    std::vector<double> F(T, 0.0);
    std::vector<std::unordered_map<std::uint64_t, long>> local_counts(64);
    std::vector<double> xeb_sum(64, 0.0);
    std::vector<long> shots_n(64, 0);

#pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < T; ++t) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num() % 64;
#endif
        std::mt19937_64 rng(opt.seed * 1000003ULL + std::uint64_t(t) * 7919ULL + 17ULL);
        std::uniform_real_distribution<double> U(0, 1);
        StateVector psi(n);
        std::uint64_t init = 0;
        if (opt.noisy && nm.thermal)
            for (int k = 0; k < n; ++k) if (U(rng) < nm.P1_init[prog.active[k]]) init |= (1ULL << k);
        psi.set_basis(init);
        execute(psi, prog, sim, nm, opt.noisy, rng);
        F[t] = std::norm(ideal.inner(psi));
        if (meas.empty() || opt.shots_per_traj <= 0) continue;
        // sample shots from this trajectory
        std::vector<double> cdf(psi.dim());
        double acc = 0;
        for (std::size_t i = 0; i < psi.dim(); ++i) { acc += std::norm(psi.amps()[i]); cdf[i] = acc; }
        for (int s = 0; s < opt.shots_per_traj; ++s) {
            const double r = U(rng) * acc;
            const std::size_t idx = std::size_t(std::lower_bound(cdf.begin(), cdf.end(), r) - cdf.begin());
            std::uint64_t key = 0;
            for (auto [q, c] : meas) {
                int bit = int((std::min(idx, psi.dim() - 1) >> q) & 1);
                if (opt.noisy && nm.readout) {
                    const int p = prog.active[q];
                    if (bit == 0 && U(rng) < nm.ro01[p]) bit = 1;
                    else if (bit == 1 && U(rng) < nm.ro10[p]) bit = 0;
                }
                if (bit) key |= (1ULL << c);
            }
            ++local_counts[tid][key];
            ++shots_n[tid];
            auto it = pideal.find(key);
            xeb_sum[tid] += (it == pideal.end()) ? 0.0 : it->second;
        }
    }
    double mean = 0; for (double f : F) mean += f; mean /= T;
    double var = 0; for (double f : F) var += (f - mean) * (f - mean);
    R.fidelity = mean;
    R.fidelity_sem = T > 1 ? std::sqrt(var / (T - 1) / T) : 0.0;
    double xs = 0; long ns = 0;
    std::unordered_map<std::uint64_t, long> merged;
    for (int k = 0; k < 64; ++k) { xs += xeb_sum[k]; ns += shots_n[k]; for (auto& [key, c] : local_counts[k]) merged[key] += c; }
    R.shots = ns;
    for (auto& [key, c] : merged) R.counts[emu_detail::key_to_bits(key, nc)] += c;
    if (ns > 0) {
        double p2 = 0; for (auto& kv : pideal) p2 += kv.second * kv.second;
        R.xeb = std::pow(2.0, double(meas.size())) * (xs / ns) - 1.0;
        R.xeb_ideal = std::pow(2.0, double(meas.size())) * p2 - 1.0;
        if (std::abs(R.xeb_ideal) > 1e-9) R.xeb_norm = R.xeb / R.xeb_ideal;
        auto it = merged.find(ktop);
        R.p_ideal_top = it == merged.end() ? 0.0 : double(it->second) / ns;
    }
    return R;
}

// Per-shot classical records of a dynamic program (one statevector trajectory per shot).
inline std::vector<std::vector<std::uint8_t>> run_records(const Compiled& prog, const NoiseModel& nm, int shots, unsigned seed, bool noisy = true) {
    const auto sim = sim_index(prog);
    const int n = int(prog.active.size());
    if (n > 28) throw std::runtime_error("too many active qubits for the statevector backend");
    std::vector<std::vector<std::uint8_t>> recs(std::size_t(shots), std::vector<std::uint8_t>(std::size_t(std::max(prog.nclbits, 1)), 0));
#pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < shots; ++t) {
        std::mt19937_64 rng(seed * 1000003ULL + std::uint64_t(t) * 7919ULL + 101ULL);
        std::uniform_real_distribution<double> U(0, 1);
        StateVector psi(n);
        std::uint64_t init = 0;
        if (noisy && nm.thermal) for (int k = 0; k < n; ++k) if (U(rng) < nm.P1_init[prog.active[k]]) init |= (1ULL << k);
        psi.set_basis(init);
        execute(psi, prog, sim, nm, noisy, rng, &recs[std::size_t(t)]);
    }
    return recs;
}

// GHZ fidelity as measured experimentally (populations + parity-oscillation amplitude):
//   F = (P_0...0 + P_1...1)/2 + |rho_{0...0, 1...1}|
// i.e. fidelity to the GHZ state with the best relative phase — insensitive to a deterministic phase offset.
struct GHZResult { double fidelity{0}, populations{0}, coherence{0}; };
inline GHZResult run_ghz_fidelity(const Compiled& prog, const NoiseModel& nm, const std::vector<int>& logical_qubits, const RunOptions& opt = {}) {
    const auto sim = sim_index(prog);
    const int n = int(prog.active.size());
    std::size_t ones = 0;
    for (int l : logical_qubits) ones |= std::size_t(1) << sim[prog.layout_final[l]];
    const int T = std::max(1, opt.trajectories);
    std::vector<double> P0(T), P1(T); std::vector<cplx> C(T);
#pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < T; ++t) {
        std::mt19937_64 rng(opt.seed * 1000003ULL + std::uint64_t(t) * 7919ULL + 29ULL);
        std::uniform_real_distribution<double> U(0, 1);
        StateVector psi(n);
        std::uint64_t init = 0;
        if (nm.thermal) for (int k = 0; k < n; ++k) if (U(rng) < nm.P1_init[prog.active[k]]) init |= (1ULL << k);
        psi.set_basis(init);
        execute(psi, prog, sim, nm, true, rng);
        P0[t] = std::norm(psi.amps()[0]); P1[t] = std::norm(psi.amps()[ones]);
        C[t] = psi.amps()[0] * std::conj(psi.amps()[ones]);
    }
    GHZResult G; cplx c = 0;
    for (int t = 0; t < T; ++t) { G.populations += (P0[t] + P1[t]) / T; c += C[t] / double(T); }
    G.coherence = 2 * std::abs(c);
    G.fidelity = 0.5 * G.populations + std::abs(c);
    return G;
}

// Exact density-matrix fidelity for small programs (validation of the trajectory backend).
inline double run_density_fidelity(const Compiled& prog, const NoiseModel& nm) {
    const auto sim = sim_index(prog);
    const int n = int(prog.active.size());
    if (n > 10) throw std::runtime_error("density-matrix backend limited to 10 qubits");
    std::mt19937_64 rng(1);
    StateVector ideal(n);
    execute(ideal, prog, sim, nm, false, rng);
    DensityMatrix rho(n);
    std::vector<double> p1(n, 0.0);
    if (nm.thermal) for (int k = 0; k < n; ++k) p1[k] = nm.P1_init[prog.active[k]];
    rho.set_product_thermal(p1);
    execute(rho, prog, sim, nm, true, rng);
    return rho.fidelity(ideal);
}

// Direct (uncompiled) ideal simulation of a logical circuit — reference for compiler tests.
inline StateVector simulate_logical(const Circuit& C) {
    StateVector s(C.n);
    for (const auto& g : C.gates) {
        switch (g.op) {
        case Op::U: s.gate1(g.q[0], g.u); break;
        case Op::CX: s.cx(g.q[0], g.q[1]); break;
        case Op::CZ: s.cz(g.q[0], g.q[1]); break;
        case Op::SWAP: s.swapq(g.q[0], g.q[1]); break;
        case Op::CP: s.cphase(g.q[0], g.q[1], g.theta); break;
        case Op::CCX: {
            const std::size_t ma = std::size_t(1) << g.q[0], mb = std::size_t(1) << g.q[1], mt = std::size_t(1) << g.q[2];
            auto& a = s.amps();
            for (std::size_t i = 0; i < a.size(); ++i) if ((i & ma) && (i & mb) && !(i & mt)) std::swap(a[i], a[i | mt]);
            break;
        }
        default: break;
        }
    }
    return s;
}

} // namespace qpuforge
