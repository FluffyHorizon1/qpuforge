// qpuforge — hardware characterisation: measured chip data -> noise model -> the roadmap's hardware gates
//
//   Characterisation          : one cooldown's measured parameters (file format "qpuforge.characterisation.v0")
//   model_from_measurement    : the design model with every measured quantity put in place of the predicted one
//   evaluate_hardware_gate    : the hardware-gate tables of ROADMAP.md, checked against a characterisation file,
//                               including "model vs measurement" (the memory experiment re-simulated with the measured
//                               parameters and compared with the measured logical error)
//   synthesize_characterisation : a virtual chip (design values with fabrication-like scatter, optional hidden error
//                               source) and the data file it would produce; used to test this chain before hardware
//
// A field that was not measured is left out of the file (or null); the model value is then kept and the gate report
// says so. Times are ns, coherence times us, frequencies GHz, ZZ kHz, errors are probabilities.
#pragma once
#include "qec.hpp"
#include <limits>
#include <random>

namespace qpuforge {

namespace chz {
inline double nan() { return std::numeric_limits<double>::quiet_NaN(); }
inline bool ok(double v) { return std::isfinite(v); }
inline double median(std::vector<double> v) {
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
    if (v.empty()) return nan();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
inline double worst_high(const std::vector<double>& v) { double m = nan(); for (double x : v) if (std::isfinite(x) && !(x <= m)) m = x; return m; }
inline double worst_low(const std::vector<double>& v) { double m = nan(); for (double x : v) if (std::isfinite(x) && !(x >= m)) m = x; return m; }
} // namespace chz

struct MeasuredQubit {
    int id{-1};
    double f01_GHz{chz::nan()};
    double T1_us{chz::nan()}, T2e_us{chz::nan()};          // energy relaxation; Hahn-echo coherence
    double err_1q{chz::nan()};                              // average error per single-qubit gate (simultaneous RB)
    double ro_p01{chz::nan()}, ro_p10{chz::nan()};          // P(read 1 | prepared 0), P(read 0 | prepared 1)
    double p_thermal{chz::nan()};                           // residual excited-state population
    double reset_err{chz::nan()};                           // P(|1> after active reset)
};
struct MeasuredPair {
    int a{-1}, b{-1};
    double err_cz{chz::nan()};      // average error per CZ (interleaved RB or XEB): includes decoherence and leakage
    double leak_cz{chz::nan()};     // leakage per CZ (leakage RB)
    double zz_kHz{chz::nan()};      // idle ZZ with the coupler parked (conditional Ramsey)
    double t_cz_ns{chz::nan()};     // calibrated CZ duration
};
struct MeasuredMemory {
    int d{0}, rounds{0};
    char basis{'Z'};
    long shots{0}, logical_errors{0};
    std::vector<double> detection_by_round;     // detection fraction per syndrome round (optional)
    double p_logical() const { return shots > 0 ? double(logical_errors) / double(shots) : chz::nan(); }
};
struct Characterisation {
    std::string schema{"qpuforge.characterisation.v0"};
    std::string chip, design, date, source{"measured"};      // source: "measured" or "synthetic"
    double t_1q_ns{chz::nan()}, t_readout_ns{chz::nan()}, t_reset_ns{chz::nan()};
    int leakage_removing_reset{-1};                           // 1 / 0, or -1 if not stated
    std::vector<MeasuredQubit> qubits;
    std::vector<MeasuredPair> pairs;
    std::vector<MeasuredMemory> memory;
    std::vector<std::string> notes;
};

// ------------------------------------------------------------------ file format
inline Json characterisation_to_json(const Characterisation& C) {
    auto put = [](Json& o, const char* k, double v) { if (chz::ok(v)) o.set(k, v); };
    Json j = Json::object();
    j.set("schema", C.schema); j.set("chip", C.chip); j.set("design", C.design); j.set("date", C.date); j.set("source", C.source);
    Json t = Json::object();
    put(t, "t_1q_ns", C.t_1q_ns); put(t, "t_readout_ns", C.t_readout_ns); put(t, "t_reset_ns", C.t_reset_ns);
    if (C.leakage_removing_reset >= 0) t.set("leakage_removing_reset", C.leakage_removing_reset != 0);
    j.set("operation", t);
    Json qs = Json::array();
    for (const auto& q : C.qubits) {
        Json o = Json::object();
        o.set("id", q.id);
        put(o, "f01_GHz", q.f01_GHz); put(o, "T1_us", q.T1_us); put(o, "T2e_us", q.T2e_us); put(o, "err_1q", q.err_1q);
        put(o, "ro_p01", q.ro_p01); put(o, "ro_p10", q.ro_p10); put(o, "p_thermal", q.p_thermal); put(o, "reset_err", q.reset_err);
        qs.push(o);
    }
    j.set("qubits", qs);
    Json ps = Json::array();
    for (const auto& p : C.pairs) {
        Json o = Json::object();
        o.set("a", p.a); o.set("b", p.b);
        put(o, "err_cz", p.err_cz); put(o, "leak_cz", p.leak_cz); put(o, "zz_kHz", p.zz_kHz); put(o, "t_cz_ns", p.t_cz_ns);
        ps.push(o);
    }
    j.set("pairs", ps);
    Json ms = Json::array();
    for (const auto& m : C.memory) {
        Json o = Json::object();
        o.set("d", m.d); o.set("rounds", m.rounds); o.set("basis", std::string(1, m.basis));
        o.set("shots", double(m.shots)); o.set("logical_errors", double(m.logical_errors));
        if (!m.detection_by_round.empty()) { Json a = Json::array(); for (double v : m.detection_by_round) a.push(v); o.set("detection_by_round", a); }
        ms.push(o);
    }
    j.set("memory", ms);
    Json ns = Json::array(); for (const auto& s : C.notes) ns.push(s);
    j.set("notes", ns);
    return j;
}

inline Characterisation characterisation_from_json(const Json& j) {
    Characterisation C;
    C.schema = j.get("schema", "");
    if (C.schema != "qpuforge.characterisation.v0") throw std::runtime_error("characterisation file: schema must be qpuforge.characterisation.v0");
    C.chip = j.get("chip", ""); C.design = j.get("design", ""); C.date = j.get("date", ""); C.source = j.get("source", "measured");
    if (j.has("operation")) {
        const Json& t = j["operation"];
        C.t_1q_ns = t.get("t_1q_ns", chz::nan()); C.t_readout_ns = t.get("t_readout_ns", chz::nan()); C.t_reset_ns = t.get("t_reset_ns", chz::nan());
        if (t.has("leakage_removing_reset")) C.leakage_removing_reset = t.get("leakage_removing_reset", false) ? 1 : 0;
    }
    if (j.has("qubits")) {
        const Json& qs = j["qubits"];
        for (std::size_t i = 0; i < qs.size(); ++i) {
            const Json& o = qs[i];
            MeasuredQubit q;
            q.id = int(o.get("id", -1.0));
            q.f01_GHz = o.get("f01_GHz", chz::nan()); q.T1_us = o.get("T1_us", chz::nan()); q.T2e_us = o.get("T2e_us", chz::nan());
            q.err_1q = o.get("err_1q", chz::nan()); q.ro_p01 = o.get("ro_p01", chz::nan()); q.ro_p10 = o.get("ro_p10", chz::nan());
            q.p_thermal = o.get("p_thermal", chz::nan()); q.reset_err = o.get("reset_err", chz::nan());
            C.qubits.push_back(q);
        }
    }
    if (j.has("pairs")) {
        const Json& ps = j["pairs"];
        for (std::size_t i = 0; i < ps.size(); ++i) {
            const Json& o = ps[i];
            MeasuredPair p;
            p.a = int(o.get("a", -1.0)); p.b = int(o.get("b", -1.0));
            p.err_cz = o.get("err_cz", chz::nan()); p.leak_cz = o.get("leak_cz", chz::nan());
            p.zz_kHz = o.get("zz_kHz", chz::nan()); p.t_cz_ns = o.get("t_cz_ns", chz::nan());
            C.pairs.push_back(p);
        }
    }
    if (j.has("memory")) {
        const Json& ms = j["memory"];
        for (std::size_t i = 0; i < ms.size(); ++i) {
            const Json& o = ms[i];
            MeasuredMemory m;
            m.d = int(o.get("d", 0.0)); m.rounds = int(o.get("rounds", 0.0));
            const std::string b = o.get("basis", "Z");
            m.basis = (b == "X" || b == "x") ? 'X' : 'Z';
            m.shots = long(o.get("shots", 0.0)); m.logical_errors = long(o.get("logical_errors", 0.0));
            if (o.has("detection_by_round")) { const Json& a = o["detection_by_round"]; for (std::size_t k = 0; k < a.size(); ++k) m.detection_by_round.push_back(a[k].num()); }
            if (m.d < 3 || m.rounds < 1 || m.shots < 1) throw std::runtime_error("characterisation file: memory entry needs d >= 3, rounds >= 1, shots >= 1");
            C.memory.push_back(m);
        }
    }
    if (j.has("notes")) { const Json& ns = j["notes"]; for (std::size_t i = 0; i < ns.size(); ++i) C.notes.push_back(ns[i].str()); }
    return C;
}

// ------------------------------------------------------------------ measured parameters -> model
// Gate errors are measured as average infidelities that already contain decoherence (and, for the CZ, leakage). The
// emulators add decoherence and leakage themselves, so only the remainder is passed on as control error:
//   single qubit : err_1q - (t/3) (1/T1 + 1/T_phi)
//   CZ           : err_cz - (4/5)(t/2) sum_q (1/T1 + 1/T_phi) - leak_cz
// Both are floored at zero.
struct MeasuredModel {
    ChipDesign D;
    ChipReport R;
    NoiseModel nm;
    int qubits_measured{0}, pairs_measured{0}, pairs_designed{0};
    std::vector<std::string> notes;
};

inline MeasuredModel model_from_measurement(const ChipDesign& D0, const ChipReport& R0, const Characterisation& C) {
    MeasuredModel M;
    M.D = D0; M.R = R0;
    ChipDesign& D = M.D; ChipReport& R = M.R;
    if (chz::ok(C.t_1q_ns)) D.t_1q_ns = C.t_1q_ns;
    if (chz::ok(C.t_readout_ns)) D.t_readout_ns = C.t_readout_ns;
    if (chz::ok(C.t_reset_ns)) D.t_reset_ns = C.t_reset_ns;
    if (C.leakage_removing_reset >= 0) D.leakage.reset_removes_leakage = C.leakage_removing_reset != 0;
    auto rate = [](const QubitReport& q) { return 1.0 / q.T1_us + (std::isfinite(q.Tphi_us) ? 1.0 / q.Tphi_us : 0.0); };   // 1/us
    std::vector<double> ctrl1, resets;
    for (const auto& mq : C.qubits) {
        int i = -1;
        for (std::size_t k = 0; k < D.qubits.size(); ++k) if (D.qubits[k].id == mq.id) i = int(k);
        if (i < 0) { M.notes.push_back("qubit id " + std::to_string(mq.id) + " is not in the design; ignored"); continue; }
        QubitReport& q = R.q[std::size_t(i)];
        ++M.qubits_measured;
        if (chz::ok(mq.T1_us) && mq.T1_us > 0) q.T1_us = mq.T1_us;
        if (chz::ok(mq.T2e_us) && mq.T2e_us > 0) {
            const double g = 1.0 / mq.T2e_us - 0.5 / q.T1_us;
            q.Tphi_us = g > 1e-9 ? 1.0 / g : std::numeric_limits<double>::infinity();
            q.T2_us = mq.T2e_us;
        }
        if (chz::ok(mq.p_thermal)) q.p_thermal = std::clamp(mq.p_thermal, 0.0, 0.5);
        if (chz::ok(mq.ro_p01)) q.ro_err_0to1 = std::clamp(mq.ro_p01, 0.0, 0.5);
        if (chz::ok(mq.ro_p10)) q.ro_err_1to0 = std::clamp(mq.ro_p10, 0.0, 0.5);
        if (chz::ok(mq.err_1q)) { q.err_1q = mq.err_1q; ctrl1.push_back(std::max(0.0, mq.err_1q - D.t_1q_ns * 1e-3 / 3.0 * rate(q))); }
        if (chz::ok(mq.reset_err)) resets.push_back(mq.reset_err);
    }
    if (!ctrl1.empty()) { double s = 0; for (double v : ctrl1) s += v; D.err_1q_control = s / double(ctrl1.size()); }
    if (!resets.empty()) { double s = 0; for (double v : resets) s += v; D.reset_err = s / double(resets.size()); }
    for (const auto& p : R.pairs) if (p.designed) ++M.pairs_designed;
    for (const auto& mp : C.pairs) {
        PairReport* p = nullptr;
        for (auto& x : R.pairs) if (x.designed && ((x.a == mp.a && x.b == mp.b) || (x.a == mp.b && x.b == mp.a))) p = &x;
        if (!p) { M.notes.push_back("pair " + std::to_string(mp.a) + "-" + std::to_string(mp.b) + " has no coupler in the design; ignored"); continue; }
        ++M.pairs_measured;
        if (chz::ok(mp.t_cz_ns) && mp.t_cz_ns > 0) p->t_cz_ns = mp.t_cz_ns;
        if (chz::ok(mp.zz_kHz)) p->zz_idle_kHz = mp.zz_kHz;
        if (chz::ok(mp.leak_cz)) {
            p->leak_2q = std::max(0.0, mp.leak_cz);
            if (p->leak_qubit < 0) {      // no calibration in the model: the higher-frequency qubit is the one that moves
                const auto& qa = R.q[std::size_t(D.index_of(p->a))]; const auto& qb = R.q[std::size_t(D.index_of(p->b))];
                p->leak_qubit = qa.f01 >= qb.f01 ? p->a : p->b;
            }
        }
        if (chz::ok(mp.err_cz)) {
            const auto& qa = R.q[std::size_t(D.index_of(p->a))]; const auto& qb = R.q[std::size_t(D.index_of(p->b))];
            const double deco = 0.8 * (p->t_cz_ns * 1e-3 / 2.0) * (rate(qa) + rate(qb));
            p->err_2q = mp.err_cz;
            p->err_2q_ctrl = std::max(0.0, mp.err_cz - deco - p->leak_2q);
        }
    }
    M.nm = noise_model(D, R);
    return M;
}

// ------------------------------------------------------------------ hardware gates
struct GateCheck {
    std::string name, criterion, note;
    double value{chz::nan()}, threshold{0}, worst{chz::nan()};
    bool at_least{false};        // true: value >= threshold passes; false: value <= threshold passes
    bool evaluated{false}, pass{false};
};
struct MemoryComparison { int d{0}; char basis{'Z'}; double eps_measured{0}, eps_model{0}, ratio{0}; int points{0}; };
struct HardwareGate {
    std::string phase;
    std::vector<GateCheck> checks;
    std::vector<MemoryComparison> comparisons;
    int evaluated{0}, passed{0}, missing{0};
    bool pass{false};            // every criterion evaluated and met
    MeasuredModel model;
};

namespace chz {
// logical error per round of one (d, basis) group: fit over the round counts when there are several
inline double eps_of(const std::vector<std::pair<int, double>>& rounds_pl) {
    if (rounds_pl.empty()) return nan();
    if (rounds_pl.size() == 1) return eps_from_PL(rounds_pl[0].second, rounds_pl[0].first);
    std::vector<MemoryResult> rs;
    for (auto [r, pl] : rounds_pl) { MemoryResult m; m.rounds = r; m.p_logical = pl; m.eps_per_round = eps_from_PL(pl, r); rs.push_back(m); }
    return fit_eps_per_round(rs);
}
// relative change of the detection fraction between the early and the late syndrome rounds (first round and the final
// data-readout layer excluded: both are a different kind of detector)
inline double detection_drift(const std::vector<double>& by_round, int window) {
    const int n = int(by_round.size()) - 1;                   // syndrome rounds 0 .. n-1; entry n is the final layer
    if (n < 2 * window + 1) return nan();
    double early = 0, late = 0;
    for (int k = 0; k < window; ++k) { early += by_round[std::size_t(1 + k)]; late += by_round[std::size_t(n - 1 - k)]; }
    return early > 0 ? late / early - 1.0 : nan();
}
} // namespace chz

// phase: "alpha" or "beta" (thresholds of ROADMAP.md). R0 must already carry the CZ calibration for designs with
// physical couplers, so that pairs without a measurement fall back to calibrated values.
inline HardwareGate evaluate_hardware_gate(const std::string& phase, const ChipDesign& D0, const ChipReport& R0,
                                           const Characterisation& C, long shots, std::uint64_t seed) {
    if (phase != "alpha" && phase != "beta") throw std::runtime_error("hardware gate: phase must be alpha or beta");
    const bool alpha = phase == "alpha";
    HardwareGate G; G.phase = phase;
    G.model = model_from_measurement(D0, R0, C);
    std::vector<double> T1, T2, e1, ecz, lcz, ro, zz;
    for (const auto& q : C.qubits) {
        T1.push_back(q.T1_us); T2.push_back(q.T2e_us); e1.push_back(q.err_1q);
        ro.push_back((chz::ok(q.ro_p01) && chz::ok(q.ro_p10)) ? 0.5 * (q.ro_p01 + q.ro_p10) : chz::nan());
    }
    for (const auto& p : C.pairs) { ecz.push_back(p.err_cz); lcz.push_back(p.leak_cz); zz.push_back(chz::ok(p.zz_kHz) ? std::abs(p.zz_kHz) : chz::nan()); }
    auto add = [&](const std::string& name, const std::string& crit, double value, double thr, bool at_least, double worst = chz::nan(), const std::string& note = "") {
        GateCheck c; c.name = name; c.criterion = crit; c.value = value; c.threshold = thr; c.at_least = at_least; c.worst = worst; c.note = note;
        c.evaluated = chz::ok(value);
        c.pass = c.evaluated && (at_least ? value >= thr : value <= thr);
        if (!c.evaluated && c.note.empty()) c.note = "not in the data file";
        G.checks.push_back(c);
    };
    add("Median T1 (us)", alpha ? ">= 50" : ">= 100", chz::median(T1), alpha ? 50 : 100, true, chz::worst_low(T1));
    if (alpha) add("Median T2 echo (us)", ">= 40", chz::median(T2), 40, true, chz::worst_low(T2));
    if (alpha) add("Median single-qubit gate error", "<= 1e-3", chz::median(e1), 1e-3, false, chz::worst_high(e1));
    add("Median CZ error", alpha ? "<= 5e-3" : "<= 3e-3", chz::median(ecz), alpha ? 5e-3 : 3e-3, false, chz::worst_high(ecz));
    if (alpha) add("Median CZ leakage per gate", "<= 1e-3", chz::median(lcz), 1e-3, false, chz::worst_high(lcz));
    add(std::string("Median readout assignment error"), alpha ? "<= 2%" : "<= 1%", chz::median(ro), alpha ? 0.02 : 0.01, false, chz::worst_high(ro));
    if (alpha) add("Median idle |ZZ| (kHz)", "<= 20", chz::median(zz), 20, false, chz::worst_high(zz));

    // memory experiments, grouped by (distance, basis)
    std::map<std::pair<int, char>, std::vector<const MeasuredMemory*>> groups;
    for (const auto& m : C.memory) groups[{m.d, m.basis}].push_back(&m);
    auto measured_eps = [&](int d) {              // mean over the bases that were measured
        double s = 0; int n = 0;
        for (const auto& [key, v] : groups) if (key.first == d) {
            std::vector<std::pair<int, double>> pts; for (const auto* m : v) pts.push_back({m->rounds, m->p_logical()});
            s += chz::eps_of(pts); ++n;
        }
        return n ? s / n : chz::nan();
    };
    const int d_chip = D0.surface_d;
    if (alpha) add("Distance-3 logical error per round", "<= 1.5%", measured_eps(3), 0.015, false);
    else {
        const double e3 = measured_eps(3), e5 = measured_eps(5);
        add("Measured Lambda(3->5)", ">= 1.5", (chz::ok(e3) && chz::ok(e5) && e5 > 0) ? e3 / e5 : chz::nan(), 1.5, true, chz::nan(),
            (chz::ok(e3) && chz::ok(e5) && e5 > 0 && e3 / e5 >= 2) ? "stretch target (>= 2) met" : "");
        add("Distance-5 logical error per round", "<= 0.3%", e5, 0.003, false);
    }
    // detection-event stability on the longest run of the full chip
    {
        const MeasuredMemory* longest = nullptr;
        for (const auto& m : C.memory) if (m.d == d_chip && !m.detection_by_round.empty() && (!longest || m.rounds > longest->rounds)) longest = &m;
        const int need = alpha ? 25 : 100, window = alpha ? 5 : 10;
        double drift = chz::nan(); std::string note;
        if (!longest) note = "no per-round detection data";
        else if (longest->rounds < need) note = "longest run has " + std::to_string(longest->rounds) + " rounds; the gate needs " + std::to_string(need);
        else { drift = chz::detection_drift(longest->detection_by_round, window); note = std::to_string(longest->rounds) + " rounds, first vs last " + std::to_string(window); }
        add(alpha ? "Detection-rate drift over 25 rounds" : "Detection-rate growth over >= 100 rounds", alpha ? "< 10%" : "< 5%", drift, alpha ? 0.10 : 0.05, false, chz::nan(), note);
    }
    // model vs measurement: the same experiments simulated with the measured parameters
    {
        double worst_ratio = chz::nan();
        for (const auto& [key, v] : groups) {
            if (key.first != d_chip) continue;
            std::vector<std::pair<int, double>> meas, pred;
            int k = 0;
            for (const auto* m : v) {
                meas.push_back({m->rounds, m->p_logical()});
                const MemorySetup S = prepare_memory(G.model.D, G.model.R, G.model.nm, m->rounds, key.second);
                const MemoryResult r = run_memory_frame(S, shots, seed + 101 * std::uint64_t(++k) + (key.second == 'X' ? 7 : 0));
                pred.push_back({m->rounds, r.p_logical});
            }
            MemoryComparison c; c.d = key.first; c.basis = key.second; c.points = int(v.size());
            c.eps_measured = chz::eps_of(meas); c.eps_model = chz::eps_of(pred);
            c.ratio = c.eps_model > 0 ? c.eps_measured / c.eps_model : chz::nan();
            G.comparisons.push_back(c);
            const double dev = chz::ok(c.ratio) && c.ratio > 0 ? std::max(c.ratio, 1.0 / c.ratio) : chz::nan();
            if (chz::ok(dev) && !(dev <= worst_ratio)) worst_ratio = dev;
        }
        add("Model vs measurement (larger of measured/model and its inverse)", "<= 2", worst_ratio, 2.0, false, chz::nan(),
            G.comparisons.empty() ? "no memory data at the chip's distance" : "");
    }
    for (const auto& c : G.checks) { if (c.evaluated) { ++G.evaluated; if (c.pass) ++G.passed; } else ++G.missing; }
    G.pass = G.missing == 0 && G.passed == G.evaluated;
    return G;
}

inline Json hardware_gate_json(const HardwareGate& G) {
    Json j = Json::object();
    j.set("phase", G.phase); j.set("pass", G.pass);
    j.set("evaluated", G.evaluated); j.set("passed", G.passed); j.set("missing", G.missing);
    Json cs = Json::array();
    for (const auto& c : G.checks) {
        Json o = Json::object();
        o.set("name", c.name); o.set("criterion", c.criterion); o.set("value", c.value); o.set("threshold", c.threshold);
        o.set("worst", c.worst); o.set("evaluated", c.evaluated); o.set("pass", c.pass); o.set("note", c.note);
        cs.push(o);
    }
    j.set("checks", cs);
    Json ms = Json::array();
    for (const auto& c : G.comparisons) {
        Json o = Json::object();
        o.set("d", c.d); o.set("basis", std::string(1, c.basis)); o.set("eps_measured", c.eps_measured); o.set("eps_model", c.eps_model);
        o.set("ratio", c.ratio); o.set("round_counts", c.points);
        ms.push(o);
    }
    j.set("model_vs_measurement", ms);
    j.set("qubits_measured", G.model.qubits_measured); j.set("pairs_measured", G.model.pairs_measured); j.set("pairs_designed", G.model.pairs_designed);
    Json ns = Json::array(); for (const auto& s : G.model.notes) ns.push(s);
    j.set("notes", ns);
    return j;
}

// ------------------------------------------------------------------ virtual chip (for testing the chain)
struct SynthOptions {
    unsigned seed{1};
    double t1_sigma{0.25};                 // log-normal spread of T1 around the design value
    double tls_fraction{0.06}, tls_factor{0.35};   // qubits sitting on a two-level-system defect, and their T1 factor
    double tphi_sigma{0.30};
    double ro_sigma{0.30};
    double cz_sigma{0.40}, leak_sigma{0.50};
    double zz_sigma_kHz{2.0};
    double thermal_sigma{0.30};
    double est_noise{0.03};                // relative statistical error of each reported parameter
    double zz_est_kHz{0.2};                // absolute error of each reported idle ZZ
    double hidden_dephasing_us{0};         // > 0: a dephasing time that acts during QEC but is absent from the reported T2
    std::vector<int> rounds{5, 15, 25};
    long shots{20000};
    bool both_bases{true};
};

// The "true" chip is the design with scattered parameters; the file reports those parameters with estimation noise and
// the memory experiments simulated on the true chip. The truth report is returned through `truth` when wanted.
inline Characterisation synthesize_characterisation(const ChipDesign& D, const ChipReport& R, const SynthOptions& o, ChipReport* truth = nullptr) {
    std::mt19937_64 rng(o.seed);
    std::normal_distribution<double> nd(0, 1);
    std::uniform_real_distribution<double> U(0, 1);
    auto logn = [&](double sigma) { return std::exp(sigma * nd(rng) - 0.5 * sigma * sigma); };
    auto est = [&](double v) { return v * (1 + o.est_noise * nd(rng)); };
    ChipReport T = R;
    auto rate = [](const QubitReport& q) { return 1.0 / q.T1_us + (std::isfinite(q.Tphi_us) ? 1.0 / q.Tphi_us : 0.0); };
    for (auto& q : T.q) {
        q.T1_us *= logn(o.t1_sigma);
        if (U(rng) < o.tls_fraction) q.T1_us *= o.tls_factor;
        if (std::isfinite(q.Tphi_us)) q.Tphi_us *= logn(o.tphi_sigma);
        q.T2_us = 1.0 / (0.5 / q.T1_us + (std::isfinite(q.Tphi_us) ? 1.0 / q.Tphi_us : 0.0));
        q.ro_err_0to1 = std::min(0.4, q.ro_err_0to1 * logn(o.ro_sigma));
        q.ro_err_1to0 = std::min(0.4, q.ro_err_1to0 * logn(o.ro_sigma));
        q.p_thermal = std::min(0.4, q.p_thermal * logn(o.thermal_sigma));
    }
    for (auto& p : T.pairs) {
        if (!p.designed) continue;
        const double base = p.err_2q_ctrl >= 0 ? p.err_2q_ctrl : D.err_2q_control;
        p.err_2q_ctrl = base * logn(o.cz_sigma);
        if (p.leak_2q > 0) p.leak_2q *= logn(o.leak_sigma);
        p.zz_idle_kHz += o.zz_sigma_kHz * nd(rng);
    }
    Characterisation C;
    C.source = "synthetic"; C.chip = D.name + " (virtual chip, seed " + std::to_string(o.seed) + ")"; C.design = D.name;
    C.t_1q_ns = D.t_1q_ns; C.t_readout_ns = D.t_readout_ns; C.t_reset_ns = D.t_reset_ns;
    C.leakage_removing_reset = D.leakage.reset_removes_leakage ? 1 : 0;
    for (std::size_t i = 0; i < T.q.size(); ++i) {
        const auto& q = T.q[i];
        MeasuredQubit m; m.id = D.qubits[i].id;
        m.f01_GHz = q.f01; m.T1_us = est(q.T1_us); m.T2e_us = est(q.T2_us);
        m.err_1q = est(D.err_1q_control + D.t_1q_ns * 1e-3 / 3.0 * rate(q));
        m.ro_p01 = est(q.ro_err_0to1); m.ro_p10 = est(q.ro_err_1to0); m.p_thermal = est(q.p_thermal); m.reset_err = est(D.reset_err);
        C.qubits.push_back(m);
    }
    for (const auto& p : T.pairs) {
        if (!p.designed) continue;
        const auto& qa = T.q[std::size_t(D.index_of(p.a))]; const auto& qb = T.q[std::size_t(D.index_of(p.b))];
        MeasuredPair m; m.a = p.a; m.b = p.b; m.t_cz_ns = p.t_cz_ns;
        m.leak_cz = p.leak_2q > 0 ? est(p.leak_2q) : 0.0;
        m.err_cz = est(p.err_2q_ctrl + 0.8 * (p.t_cz_ns * 1e-3 / 2.0) * (rate(qa) + rate(qb)) + p.leak_2q);
        m.zz_kHz = p.zz_idle_kHz + o.zz_est_kHz * nd(rng);
        C.pairs.push_back(m);
    }
    // the experiment runs on the true chip (with the hidden error source, if any)
    ChipReport Tq = T;
    if (o.hidden_dephasing_us > 0)
        for (auto& q : Tq.q) q.Tphi_us = 1.0 / ((std::isfinite(q.Tphi_us) ? 1.0 / q.Tphi_us : 0.0) + 1.0 / o.hidden_dephasing_us);
    const NoiseModel nm = noise_model(D, Tq);
    std::uint64_t s = 1000 + o.seed;
    for (char basis : {'Z', 'X'}) {
        if (basis == 'X' && !o.both_bases) continue;
        for (int r : o.rounds) {
            const MemorySetup S = prepare_memory(D, Tq, nm, r, basis);
            const MemoryResult res = run_memory_frame(S, o.shots, ++s);
            MeasuredMemory m; m.d = D.surface_d; m.rounds = r; m.basis = basis; m.shots = res.shots; m.logical_errors = res.logical_errors;
            m.detection_by_round = res.detection_by_round;
            C.memory.push_back(m);
        }
    }
    C.notes.push_back("Synthetic data from a virtual chip: design values with fabrication-like scatter. Not a measurement.");
    if (o.hidden_dephasing_us > 0) C.notes.push_back("A hidden dephasing source acts during the memory experiments and is absent from the reported T2.");
    if (truth) *truth = Tq;
    return C;
}

} // namespace qpuforge
