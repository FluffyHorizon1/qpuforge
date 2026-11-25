// qpuforge — superconducting transmon chip physics
//
// Lumped-element model of a chip of flux-tunable transmons:
//   * Maxwell capacitance matrix C (shunts + coupling caps) -> E_C,i and couplings via C^-1 (exact in the lumped model)
//   * Josephson energy from junction normal resistance (Ambegaokar–Baratoff), SQUID flux tuning
//   * exact transmon spectra in the charge basis (frequencies, anharmonicity, charge dispersion, n-matrix)
//   * exact two-transmon diagonalisation -> dressed frequencies and static ZZ (always-on crosstalk)
//   * optional three-transmon coupler: park below both qubits so the exchange cancels, and read the ZZ that remains
//   * coherence: dielectric loss (Q), Purcell decay, thermal excitation, flux / charge / thermal-photon dephasing
//   * CZ gate model (|11> <-> |20> interaction), flux-noise dephasing while the measure qubit is off its sweet spot,
//     feasibility, frequency-collision rules
//   * inverse design (targets -> capacitances & junction resistances) and fabrication Monte Carlo
// All frequencies inside are GHz (E/h) unless the name says otherwise.
#pragma once
#include "core.hpp"
#include "json.hpp"
#include "surface.hpp"
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace qpuforge {

// ================================================================== transmon
struct TransmonLevels {
    double EJ{0}, EC{0}, ng{0};
    int L{0};
    std::vector<double> E;     // GHz, E[0] = 0
    std::vector<double> nme;   // <i|n|j>, L x L
    double f01() const { return E[1]; }
    double alpha() const { return E[2] - 2 * E[1]; }
    double n(int i, int j) const { return nme[std::size_t(i) * L + j]; }
};

// H = 4 E_C (n - n_g)^2 - E_J cos(phi) in the charge basis |n|, n = -ncut..ncut
inline TransmonLevels transmon_levels(double EJ, double EC, double ng = 0.0, int L = 6, int ncut = 30) {
    const int dim = 2 * ncut + 1;
    std::vector<double> d(dim), e(dim - 1, -EJ / 2), z;
    for (int k = 0; k < dim; ++k) { const double nn = double(k - ncut) - ng; d[k] = 4 * EC * nn * nn; }
    tridiag_eigen(d, e, z);
    TransmonLevels T; T.EJ = EJ; T.EC = EC; T.ng = ng; T.L = L;
    for (int i = 0; i < L; ++i) T.E.push_back(d[i] - d[0]);
    T.nme.assign(std::size_t(L) * L, 0.0);
    for (int i = 0; i < L; ++i) for (int j = 0; j < L; ++j) {
        double s = 0;
        for (int k = 0; k < dim; ++k) s += z[std::size_t(k) * dim + i] * double(k - ncut) * z[std::size_t(k) * dim + j];
        T.nme[std::size_t(i) * L + j] = s;
    }
    return T;
}

// Energies only (no eigenvectors): f01 and anharmonicity, for the many repeated solves in design loops.
struct TransmonFA { double f01, alpha; };
inline TransmonFA transmon_fa(double EJ, double EC, double ng = 0.0, int ncut = 30) {
    const int dim = 2 * ncut + 1;
    std::vector<double> d(dim), e(dim - 1, -EJ / 2), z;
    for (int k = 0; k < dim; ++k) { const double nn = double(k - ncut) - ng; d[k] = 4 * EC * nn * nn; }
    tridiag_eigen(d, e, z, false);
    return {d[1] - d[0], (d[2] - d[1]) - (d[1] - d[0])};
}

inline double EJ_from_Rn_GHz(double Rn_kOhm, double gap_ueV) {            // Ambegaokar–Baratoff, T = 0
    return gap_ueV * phys::ueV / (8 * phys::e * phys::e * Rn_kOhm * 1e3) / 1e9;
}
inline double Rn_from_EJ_kOhm(double EJ_GHz, double gap_ueV) {
    return gap_ueV * phys::ueV / (8 * phys::e * phys::e * EJ_GHz * 1e9) / 1e3;
}
inline double EC_from_C_GHz(double C_F) { return phys::e * phys::e / (2 * C_F * phys::h) / 1e9; }
inline double C_from_EC_F(double EC_GHz) { return phys::e * phys::e / (2 * EC_GHz * 1e9 * phys::h); }
// Asymmetric dc-SQUID: E_J(Phi) = E_JSigma sqrt(cos^2(pi Phi) + d^2 sin^2(pi Phi)), Phi in units of Phi0
inline double EJ_at_flux(double EJsum, double d, double flux) {
    const double c = std::cos(PI * flux), s = std::sin(PI * flux);
    return EJsum * std::sqrt(c * c + d * d * s * s);
}

// Solve (E_J, E_C) for a target f01 and anharmonicity (both GHz, alpha < 0).
inline double solve_EJ_for_f01(double f01, double EC) {
    // f01 is monotonic in E_J: bracketed secant/bisection hybrid on log(E_J)
    double lo = std::log(0.5), hi = std::log(2000.0);
    for (int b = 0; b < 200 && hi - lo > 1e-13; ++b) {
        const double mid = 0.5 * (lo + hi);
        (transmon_fa(std::exp(mid), EC).f01 < f01 ? lo : hi) = mid;
    }
    return std::exp(0.5 * (lo + hi));
}
inline std::pair<double, double> solve_EJ_EC(double f01, double alpha) {
    double EC = -alpha, EJ = 10;
    for (int it = 0; it < 60; ++it) {
        EJ = solve_EJ_for_f01(f01, EC);
        const double ratio = alpha / transmon_fa(EJ, EC).alpha;
        EC *= ratio;
        if (std::abs(ratio - 1) < 1e-11) break;
    }
    return {solve_EJ_for_f01(f01, EC), EC};
}

// f01 of a transmon with this EC such that f01 + alpha(EJ) = f_partner.
// That is the |11> = |20> point. Idle alpha is not enough: |alpha| grows as EJ falls.
inline double f01_for_11_20(double EC, double f_partner) {
    double lo = std::log(0.2), hi = std::log(4000.0);
    for (int b = 0; b < 80 && hi - lo > 1e-13; ++b) {
        const double mid = 0.5 * (lo + hi);
        const auto fa = transmon_fa(std::exp(mid), EC);
        ((fa.f01 + fa.alpha) < f_partner ? lo : hi) = mid;
    }
    return transmon_fa(std::exp(0.5 * (lo + hi)), EC).f01;
}

// 1/f flux noise, S_Phi(f) = A^2 / |f|, A in Phi0. The same rate the idle T_phi uses,
// with the experiment duration made explicit (idle Ramsey reference is 10 us).
// Angular-to-rate conversion matches the historical idle formula exactly at t = 10 us.
inline double flux_dephasing_rate(double A_Phi0, double df_dflux_GHz, double d2f_dflux2_GHz, double t_s) {
    const double ln_ir = std::abs(std::log(2 * PI * std::max(t_s, 1e-30)));
    return 2 * PI * (A_Phi0 * std::abs(df_dflux_GHz) * 1e9 * std::sqrt(2 * ln_ir)
                     + A_Phi0 * A_Phi0 * std::abs(d2f_dflux2_GHz) * 1e9 * ln_ir);
}

// Short-time gate infidelity from that rate. 1/f phase variance grows as t^2, so a
// Markovian t * Gamma overstates a 20 ns pulse by ~1/(Gamma t). One Z on a two-qubit
// gate has average infidelity (4/5) * p_Z, and p_Z = (1 - exp(-(Gamma t)^2)) / 2.
inline double flux_pulse_infidelity(double gamma_per_s, double t_s) {
    const double x = gamma_per_s * t_s;
    return 0.4 * (1.0 - std::exp(-x * x));
}

// Capacitive coupler parked below both qubits. g_ab and g_qc are |0>-|1> exchange
// couplings in GHz (same sign: the mediated term cancels g_ab only from below).
// The exchange null is not a ZZ null. This returns the ZZ the three-body spectrum
// actually has once the perturbative exchange has been tuned through zero.
struct CouplerNull { double f_GHz{0}, J_MHz{0}, zz_kHz{0}; bool ok{false}; };
inline CouplerNull coupler_exchange_null(const TransmonLevels& A, const TransmonLevels& B,
                                         double g_ab_GHz, double g_qc_GHz, double EC_c_GHz) {
    CouplerNull R;
    if (!(g_qc_GHz > 0) || !(EC_c_GHz > 0) || A.L < 3 || B.L < 3) return R;
    const double fa = A.f01(), fb = B.f01();
    auto Jp = [&](double fc) {
        return g_ab_GHz + 0.5 * g_qc_GHz * g_qc_GHz * (1.0 / (fc - fa) + 1.0 / (fc - fb));
    };
    double lo = 2.0, hi = std::min(fa, fb) - 0.30;
    if (!(hi > lo) || Jp(lo) * Jp(hi) > 0) return R;
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        (Jp(mid) > 0 ? lo : hi) = mid;
    }
    const double fc = 0.5 * (lo + hi);
    const auto C = transmon_levels(solve_EJ_for_f01(fc, EC_c_GHz), EC_c_GHz, 0.0, 3);
    const int L = 3, dim = L * L * L;
    const double nA = std::abs(A.n(0, 1)), nB = std::abs(B.n(0, 1)), nC = std::abs(C.n(0, 1));
    if (nA < 1e-8 || nB < 1e-8 || nC < 1e-8) return R;
    const double Gab = g_ab_GHz / (nA * nB);
    const double Gac = g_qc_GHz / (nA * nC);
    const double Gbc = g_qc_GHz / (nB * nC);
    std::vector<double> H(std::size_t(dim) * dim, 0.0), V;
    auto idx = [&](int a, int b, int c) { return (a * L + b) * L + c; };
    for (int a = 0; a < L; ++a) for (int b = 0; b < L; ++b) for (int c = 0; c < L; ++c) {
        const int s = idx(a, b, c);
        H[std::size_t(s) * dim + s] = A.E[a] + B.E[b] + C.E[c];
        for (int a2 = 0; a2 < L; ++a2) for (int b2 = 0; b2 < L; ++b2) for (int c2 = 0; c2 < L; ++c2) {
            double v = 0;
            if (c == c2) v += Gab * A.n(a, a2) * B.n(b, b2);
            if (b == b2) v += Gac * A.n(a, a2) * C.n(c, c2);
            if (a == a2) v += Gbc * B.n(b, b2) * C.n(c, c2);
            H[std::size_t(s) * dim + idx(a2, b2, c2)] += v;
        }
    }
    const auto ev = sym_eigen(H, dim, V);
    auto dressed = [&](int a, int b, int c) {
        const int s = idx(a, b, c);
        int best = 0; double bw = -1;
        for (int m = 0; m < dim; ++m) {
            const double w = V[std::size_t(s) * dim + m] * V[std::size_t(s) * dim + m];
            if (w > bw) { bw = w; best = m; }
        }
        return ev[best];
    };
    const double E00 = dressed(0, 0, 0), E10 = dressed(1, 0, 0), E01 = dressed(0, 1, 0), E11 = dressed(1, 1, 0);
    R.f_GHz = C.f01();
    R.J_MHz = Jp(R.f_GHz) * 1e3;
    R.zz_kHz = (E11 - E10 - E01 + E00) * 1e6;
    R.ok = true;
    return R;
}

// ================================================================== design description
struct ReadoutDesign { double f_GHz{7.0}, g_MHz{90}, kappa_MHz{1.0}, purcell_filter{30}; };
struct QubitDesign {
    int id{0};
    double x_um{0}, y_um{0};
    double C_fF{80};          // shunt capacitance to ground
    double Rn_kOhm{10};       // total normal-state resistance of the junction / SQUID
    bool tunable{true};
    double squid_d{0};        // junction asymmetry of the SQUID
    double flux{0};           // idle flux bias, Phi/Phi0
    ReadoutDesign ro;
    std::string role{"qubit"};   // data | measure-X | measure-Z | qubit
    int feedline{0};             // readout feedline the resonator sits on
};
struct CouplingDesign { int a{0}, b{0}; double C_fF{1.0}; };

struct ChipDesign {
    std::string name{"chip"};
    std::string architecture{"fixed-coupling-cz"};    // or "tunable-coupler-cz"
    // environment & materials
    double T_eff_mK{45}, flux_noise_uPhi0{2.0}, charge_noise_e{1e-4}, Tphi_other_us{1e4};
    double gap_ueV{180}, Rn_spread_pct{2.0}, Q_dielectric{3e6};
    // gate & readout calibration-level parameters
    double t_1q_ns{25}, t_readout_ns{800}, t_reset_ns{200}, cz_shape_factor{1.5}, cz_ramp_ns{5};
    double err_1q_control{1e-4}, err_2q_control{2e-3}, readout_sep_err{0.005};
    double tc_g_max_MHz{40}, tc_residual_zz_kHz{5};
    double tc_g_qc_MHz{40}, tc_coupler_EC{0.20};   // coupler–qubit exchange, coupler charging energy (GHz)
    bool tc_model{false};                          // compute idle ZZ from the three-transmon null; else use tc_residual_zz_kHz
    bool active_reset{false};
    double reset_err{0.002};
    std::vector<QubitDesign> qubits;
    std::vector<CouplingDesign> couplings;
    int surface_d{0};             // > 0: chip is laid out as a rotated surface code of this distance (qubit id = code index)

    int index_of(int id) const {
        for (std::size_t i = 0; i < qubits.size(); ++i) if (qubits[i].id == id) return int(i);
        throw std::runtime_error("unknown qubit id " + std::to_string(id));
    }
    bool tunable_coupler() const { return architecture == "tunable-coupler-cz"; }
};

inline ChipDesign design_from_json(const Json& j) {
    ChipDesign D;
    D.name = j.get("name", "chip");
    D.architecture = j.get("architecture", "fixed-coupling-cz");
    const Json& env = j["environment"];
    D.T_eff_mK = env.get("T_eff_mK", D.T_eff_mK);
    D.flux_noise_uPhi0 = env.get("flux_noise_uPhi0", D.flux_noise_uPhi0);
    D.charge_noise_e = env.get("charge_noise_e", D.charge_noise_e);
    D.Tphi_other_us = env.get("Tphi_other_us", D.Tphi_other_us);
    const Json& mat = j["materials"];
    D.gap_ueV = mat.get("junction_gap_ueV", D.gap_ueV);
    D.Rn_spread_pct = mat.get("Rn_spread_pct", D.Rn_spread_pct);
    D.Q_dielectric = mat.get("Q_dielectric", D.Q_dielectric);
    const Json& g = j["gates"];
    D.t_1q_ns = g.get("t_1q_ns", D.t_1q_ns);
    D.t_readout_ns = g.get("t_readout_ns", D.t_readout_ns);
    D.t_reset_ns = g.get("t_reset_ns", D.t_reset_ns);
    D.cz_shape_factor = g.get("cz_shape_factor", D.cz_shape_factor);
    D.cz_ramp_ns = g.get("cz_ramp_ns", D.cz_ramp_ns);
    D.err_1q_control = g.get("err_1q_control", D.err_1q_control);
    D.err_2q_control = g.get("err_2q_control", D.err_2q_control);
    D.readout_sep_err = g.get("readout_separation_err", D.readout_sep_err);
    D.active_reset = g.get("active_reset", D.active_reset);
    D.reset_err = g.get("reset_err", D.reset_err);
    const Json& tc = j["tunable_coupler"];
    D.tc_g_max_MHz = tc.get("g_max_MHz", D.tc_g_max_MHz);
    D.tc_residual_zz_kHz = tc.get("residual_zz_kHz", D.tc_residual_zz_kHz);
    D.tc_g_qc_MHz = tc.get("g_qc_MHz", D.tc_g_qc_MHz);
    D.tc_coupler_EC = tc.get("coupler_EC_GHz", D.tc_coupler_EC);
    D.tc_model = tc.get("model_coupler", D.tc_model);
    const Json& qs = j["qubits"];
    for (std::size_t i = 0; i < qs.size(); ++i) {
        const Json& q = qs[i];
        QubitDesign Q;
        Q.id = int(q.get("id", double(i)));
        Q.x_um = q.get("x_um", 0.0); Q.y_um = q.get("y_um", 0.0);
        Q.C_fF = q.get("C_fF", Q.C_fF);
        Q.Rn_kOhm = q.get("Rn_kOhm", Q.Rn_kOhm);
        Q.tunable = q.get("tunable", Q.tunable);
        Q.squid_d = q.get("squid_asymmetry", Q.squid_d);
        Q.flux = q.get("flux_bias", Q.flux);
        const Json& r = q["readout"];
        Q.ro.f_GHz = r.get("f_GHz", Q.ro.f_GHz);
        Q.ro.g_MHz = r.get("g_MHz", Q.ro.g_MHz);
        Q.ro.kappa_MHz = r.get("kappa_MHz", Q.ro.kappa_MHz);
        Q.ro.purcell_filter = r.get("purcell_filter", Q.ro.purcell_filter);
        Q.role = q.get("role", "qubit");
        Q.feedline = int(q.get("feedline", 0.0));
        D.qubits.push_back(Q);
    }
    D.surface_d = int(j["surface_code"].get("distance", 0.0));
    const Json& cs = j["couplings"];
    for (std::size_t i = 0; i < cs.size(); ++i) {
        CouplingDesign c;
        c.a = int(cs[i].get("a", 0.0)); c.b = int(cs[i].get("b", 0.0));
        c.C_fF = cs[i].get("C_fF", 1.0);
        D.couplings.push_back(c);
    }
    return D;
}

inline Json design_to_json(const ChipDesign& D) {
    Json j = Json::object();
    j.set("name", D.name);
    j.set("architecture", D.architecture);
    Json env = Json::object();
    env.set("T_eff_mK", D.T_eff_mK); env.set("flux_noise_uPhi0", D.flux_noise_uPhi0);
    env.set("charge_noise_e", D.charge_noise_e); env.set("Tphi_other_us", D.Tphi_other_us);
    j.set("environment", env);
    Json mat = Json::object();
    mat.set("junction_gap_ueV", D.gap_ueV); mat.set("Rn_spread_pct", D.Rn_spread_pct); mat.set("Q_dielectric", D.Q_dielectric);
    j.set("materials", mat);
    Json g = Json::object();
    g.set("t_1q_ns", D.t_1q_ns); g.set("t_readout_ns", D.t_readout_ns); g.set("t_reset_ns", D.t_reset_ns); g.set("cz_shape_factor", D.cz_shape_factor);
    g.set("cz_ramp_ns", D.cz_ramp_ns); g.set("err_1q_control", D.err_1q_control); g.set("err_2q_control", D.err_2q_control);
    g.set("readout_separation_err", D.readout_sep_err); g.set("active_reset", D.active_reset); g.set("reset_err", D.reset_err);
    j.set("gates", g);
    Json tc = Json::object();
    tc.set("g_max_MHz", D.tc_g_max_MHz); tc.set("residual_zz_kHz", D.tc_residual_zz_kHz);
    tc.set("g_qc_MHz", D.tc_g_qc_MHz); tc.set("coupler_EC_GHz", D.tc_coupler_EC); tc.set("model_coupler", D.tc_model);
    j.set("tunable_coupler", tc);
    Json qs = Json::array();
    for (const auto& q : D.qubits) {
        Json o = Json::object();
        o.set("id", q.id); o.set("role", q.role); o.set("x_um", q.x_um); o.set("y_um", q.y_um); o.set("feedline", q.feedline);
        o.set("C_fF", q.C_fF); o.set("Rn_kOhm", q.Rn_kOhm); o.set("tunable", q.tunable);
        o.set("squid_asymmetry", q.squid_d); o.set("flux_bias", q.flux);
        Json r = Json::object();
        r.set("f_GHz", q.ro.f_GHz); r.set("g_MHz", q.ro.g_MHz); r.set("kappa_MHz", q.ro.kappa_MHz); r.set("purcell_filter", q.ro.purcell_filter);
        o.set("readout", r);
        qs.push(o);
    }
    if (D.surface_d > 0) { Json sc = Json::object(); sc.set("distance", D.surface_d); j.set("surface_code", sc); }
    j.set("qubits", qs);
    Json cs = Json::array();
    for (const auto& c : D.couplings) { Json o = Json::object(); o.set("a", c.a); o.set("b", c.b); o.set("C_fF", c.C_fF); cs.push(o); }
    j.set("couplings", cs);
    return j;
}

inline ChipDesign load_design(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    return design_from_json(Json::parse(ss.str()));
}
inline void save_design(const ChipDesign& D, const std::string& path) {
    std::ofstream f(path);
    design_to_json(D).dump(f);
    f << "\n";
}

// ================================================================== analysis
struct QubitReport {
    int id{0};
    double EJ{0}, EC{0}, EJ_EC{0}, f01{0}, alpha{0}, f01_dressed{0}, charge_disp_kHz{0};
    double df_dflux{0}, d2f_dflux2{0};               // GHz/Phi0, GHz/Phi0^2
    double f_max{0};                                   // GHz at zero flux (top of tuning range)
    double T1_us{0}, T1_diel_us{0}, T1_purcell_us{0};
    double Tphi_us{0}, Tphi_flux_us{0}, Tphi_photon_us{0}, Tphi_charge_us{0}, T2_us{0};
    double p_thermal{0}, chi_MHz{0}, ro_err_0to1{0}, ro_err_1to0{0}, err_1q{0};
    double n01{0};
    TransmonLevels levels;
};
struct PairReport {
    int a{0}, b{0};                 // qubit ids
    bool designed{false};           // has a coupler (2q gates allowed)
    double g_MHz{0}, delta_MHz{0}, zz_kHz{0}, zz_pert_kHz{0}, zz_idle_kHz{0};
    double t_cz_ns{0}, err_2q{0}, cz_interaction_GHz{0}, cz_excursion_MHz{0}, cz_alpha_GHz{0};
    double cz_flux{0}, cz_df_dflux{0}, cz_flux_err{0};   // Phi/Phi0 at the CZ point, GHz/Phi0, short-time infidelity
    double coupler_f_GHz{0}, coupler_J_MHz{0};          // set when the three-transmon idle model runs
    std::string flags;
};
struct ChipReport {
    std::vector<QubitReport> q;
    std::vector<PairReport> pairs;
    std::vector<std::string> warnings;
    int n_collisions{0};
    double mean_err_1q{0}, mean_err_2q{0}, mean_T1_us{0}, mean_T2_us{0}, max_abs_zz_kHz{0};
};

struct AnalyzeOptions { int levels{5}; double parasitic_g_MHz{0.1}; };

inline ChipReport analyze(const ChipDesign& D, const AnalyzeOptions& opt = {}) {
    const int n = int(D.qubits.size());
    ChipReport R;
    if (n == 0) return R;

    // --- Maxwell capacitance matrix and its inverse
    std::vector<double> C(std::size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) C[std::size_t(i) * n + i] += D.qubits[i].C_fF * 1e-15;
    std::set<std::pair<int, int>> designed;
    for (const auto& c : D.couplings) {
        const int i = D.index_of(c.a), j = D.index_of(c.b);
        const double cc = c.C_fF * 1e-15;
        C[std::size_t(i) * n + i] += cc; C[std::size_t(j) * n + j] += cc;
        C[std::size_t(i) * n + j] -= cc; C[std::size_t(j) * n + i] -= cc;
        designed.insert({std::min(i, j), std::max(i, j)});
    }
    const std::vector<double> Ci = inverse(C, n);
    auto ci = [&](int i, int j) { return Ci[std::size_t(i) * n + j]; };

    const double T = D.T_eff_mK * 1e-3;
    auto bose = [&](double f_Hz) { const double x = phys::h * f_Hz / (phys::kB * T); return x > 700 ? 0.0 : 1.0 / std::expm1(x); };

    // --- per-qubit physics
    R.q.resize(n);
    for (int i = 0; i < n; ++i) {
        const auto& Q = D.qubits[i];
        auto& r = R.q[i];
        r.id = Q.id;
        r.EC = phys::e * phys::e * ci(i, i) / (2 * phys::h) / 1e9;
        const double EJsum = EJ_from_Rn_GHz(Q.Rn_kOhm, D.gap_ueV);
        r.EJ = Q.tunable ? EJ_at_flux(EJsum, Q.squid_d, Q.flux) : EJsum;
        r.EJ_EC = r.EJ / r.EC;
        r.levels = transmon_levels(r.EJ, r.EC, 0.0, std::max(opt.levels, 3) + 1);
        r.f01 = r.levels.f01();
        r.alpha = r.levels.alpha();
        r.n01 = std::abs(r.levels.n(0, 1));
        r.f_max = Q.tunable ? transmon_fa(EJsum, r.EC).f01 : r.f01;
        r.charge_disp_kHz = std::abs(transmon_fa(r.EJ, r.EC, 0.5).f01 - r.f01) * 1e6;
        if (Q.tunable) {
            const double h1 = 1e-4;
            auto f_at = [&](double fl) { return transmon_fa(EJ_at_flux(EJsum, Q.squid_d, fl), r.EC).f01; };
            const double fp = f_at(Q.flux + h1), fm = f_at(Q.flux - h1);
            r.df_dflux = (fp - fm) / (2 * h1);
            r.d2f_dflux2 = (fp - 2 * r.f01 + fm) / (h1 * h1);
        }
        const double f = r.f01 * 1e9;
        // relaxation
        const double G_diel = 2 * PI * f / D.Q_dielectric;
        const double dr = (Q.ro.f_GHz - r.f01) * 1e9, gr = Q.ro.g_MHz * 1e6, kap = 2 * PI * Q.ro.kappa_MHz * 1e6;
        const double G_purcell = kap * (gr / dr) * (gr / dr) / std::max(Q.ro.purcell_filter, 1.0);
        const double G1 = G_diel + G_purcell;
        // dephasing
        const double A = D.flux_noise_uPhi0 * 1e-6;
        const double G_flux = Q.tunable ? flux_dephasing_rate(A, r.df_dflux, r.d2f_dflux2, 10e-6) : 0.0;
        const double G_charge = 2 * PI * PI * D.charge_noise_e * r.charge_disp_kHz * 1e3;
        const double delta_r = r.f01 - Q.ro.f_GHz;   // GHz
        const double chi = (Q.ro.g_MHz * 1e-3) * (Q.ro.g_MHz * 1e-3) * r.alpha / (delta_r * (delta_r + r.alpha));   // GHz (Koch)
        r.chi_MHz = chi * 1e3;
        const double chi2 = 2 * PI * std::abs(2 * chi) * 1e9;   // full dispersive shift, rad/s
        const double nr = bose(Q.ro.f_GHz * 1e9);
        const double G_photon = nr * kap * chi2 * chi2 / (kap * kap + chi2 * chi2);
        const double G_other = 1.0 / (D.Tphi_other_us * 1e-6);
        const double Gphi = G_flux + G_charge + G_photon + G_other;
        auto us = [](double rate) { return rate > 0 ? 1e6 / rate : std::numeric_limits<double>::infinity(); };
        r.T1_us = us(G1); r.T1_diel_us = us(G_diel); r.T1_purcell_us = us(G_purcell);
        r.Tphi_us = us(Gphi); r.Tphi_flux_us = us(G_flux); r.Tphi_charge_us = us(G_charge); r.Tphi_photon_us = us(G_photon);
        r.T2_us = us(0.5 * G1 + Gphi);
        r.p_thermal = 1.0 / (1.0 + std::exp(phys::h * f / (phys::kB * T)));
        (void)bose(f);
        r.ro_err_1to0 = (1 - std::exp(-D.t_readout_ns * 1e-9 * G1 / 2)) + D.readout_sep_err;
        r.ro_err_0to1 = D.readout_sep_err;
        r.err_1q = D.t_1q_ns * 1e-9 / 3 * (G1 + Gphi) + D.err_1q_control;
        r.f01_dressed = r.f01;

        if (r.EJ_EC < 30) R.warnings.push_back("Q" + std::to_string(Q.id) + ": E_J/E_C = " + std::to_string(r.EJ_EC).substr(0, 5) + " (< 30) -> charge-noise sensitive");
        if (std::abs(delta_r) < 0.75) R.warnings.push_back("Q" + std::to_string(Q.id) + ": readout resonator only " + std::to_string(int(std::abs(delta_r) * 1000)) + " MHz from qubit");
        if (G_purcell > G_diel) R.warnings.push_back("Q" + std::to_string(Q.id) + ": T1 is Purcell-limited (add filtering / detune resonator)");
    }

    // --- bare exchange couplings g_ij = 4 e^2 (C^-1)_ij / h * |n01_i| |n01_j|  (all positive for capacitive coupling)
    std::vector<double> gb(std::size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        if (i != j) gb[std::size_t(i) * n + j] = 4 * phys::e * phys::e * ci(i, j) / phys::h / 1e9 * R.q[i].n01 * R.q[j].n01;
    // effective coupling including exchange mediated by every other qubit (2nd-order Schrieffer-Wolff)
    auto J_eff = [&](int a, int b) {
        double J = gb[std::size_t(a) * n + b];
        for (int k = 0; k < n; ++k) {
            if (k == a || k == b) continue;
            if (D.tunable_coupler() && (designed.count({std::min(a, k), std::max(a, k)}) || designed.count({std::min(k, b), std::max(k, b)})))
                continue;   // a parked tunable coupler cancels the a-k (or k-b) exchange, so nothing is mediated through it
            const double da = R.q[a].f01 - R.q[k].f01, db = R.q[b].f01 - R.q[k].f01;
            if (std::abs(da) < 1e-6 || std::abs(db) < 1e-6) continue;
            J += 0.5 * gb[std::size_t(a) * n + k] * gb[std::size_t(k) * n + b] * (1 / da + 1 / db);
        }
        return J;
    };

    // --- pairs: designed couplers + parasitic couplings above threshold + near-degenerate pairs
    const double tc_g = D.tc_g_max_MHz * 1e-3;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
        const bool des = designed.count({i, j}) > 0;
        const double G = 4 * phys::e * phys::e * ci(i, j) / phys::h / 1e9;   // GHz per (n_i n_j)
        const auto& qa = R.q[i]; const auto& qb = R.q[j];
        const double g_direct = std::abs(G * qa.levels.n(0, 1) * qb.levels.n(0, 1));
        const double g = des ? g_direct : std::abs(J_eff(i, j));
        // keep designed couplers, parasitic couplings above threshold, and ANY weak coupling between near-degenerate qubits
        const bool near_degenerate = g > 1e-6 && std::abs(qa.f01 - qb.f01) < 10 * g;
        if (!des && g * 1e3 < opt.parasitic_g_MHz && !near_degenerate) continue;
        PairReport p;
        p.a = qa.id; p.b = qb.id; p.designed = des;
        p.g_MHz = g * 1e3;
        p.delta_MHz = (qa.f01 - qb.f01) * 1e3;
        // exact diagonalisation of the coupled pair
        const int L = opt.levels, dim = L * L;
        std::vector<double> H(std::size_t(dim) * dim, 0.0), V;
        for (int x = 0; x < L; ++x) for (int y = 0; y < L; ++y) {
            const int s = x * L + y;
            H[std::size_t(s) * dim + s] = qa.levels.E[x] + qb.levels.E[y];
            for (int x2 = 0; x2 < L; ++x2) for (int y2 = 0; y2 < L; ++y2)
                H[std::size_t(s) * dim + (x2 * L + y2)] += G * qa.levels.n(x, x2) * qb.levels.n(y, y2);
        }
        const auto ev = sym_eigen(H, dim, V);   // direct capacitive coupling only; mediated exchange enters via J_eff
        auto dressed = [&](int x, int y) {
            const int s = x * L + y; int best = 0; double bw = -1;
            for (int m = 0; m < dim; ++m) { const double w = V[std::size_t(s) * dim + m] * V[std::size_t(s) * dim + m]; if (w > bw) { bw = w; best = m; } }
            return ev[best];
        };
        const double E00 = dressed(0, 0), E01 = dressed(0, 1), E10 = dressed(1, 0), E11 = dressed(1, 1);
        p.zz_kHz = (E11 - E10 - E01 + E00) * 1e6;
        const double Dl = qa.f01 - qb.f01;
        p.zz_pert_kHz = 2 * g * g * (qa.alpha + qb.alpha) / ((Dl + qa.alpha) * (Dl - qb.alpha)) * 1e6;
        R.q[i].f01_dressed += (E10 - E00) - qa.f01;
        R.q[j].f01_dressed += (E01 - E00) - qb.f01;
        p.zz_idle_kHz = (des && D.tunable_coupler()) ? D.tc_residual_zz_kHz : p.zz_kHz;
        if (des && D.tunable_coupler() && D.tc_model) {
            const auto cn = coupler_exchange_null(qa.levels, qb.levels, g, D.tc_g_qc_MHz * 1e-3, D.tc_coupler_EC);
            if (cn.ok) {
                p.coupler_f_GHz = cn.f_GHz;
                p.coupler_J_MHz = cn.J_MHz;
                p.zz_idle_kHz = cn.zz_kHz;
            } else {
                R.warnings.push_back("pair Q" + std::to_string(qa.id) + "-Q" + std::to_string(qb.id) +
                                     ": no coupler exchange-null below the qubits; idle ZZ left at the typed residual");
            }
        }

        if (des) {
            const double geff = D.tunable_coupler() ? tc_g : g;
            p.t_cz_ns = D.cz_shape_factor / (2 * std::sqrt(2.0) * geff) + 2 * D.cz_ramp_ns;
            const auto& hi = (qa.f01 >= qb.f01) ? qa : qb;
            const auto& lo = (qa.f01 >= qb.f01) ? qb : qa;
            const auto& hiD = D.qubits[(qa.f01 >= qb.f01) ? i : j];
            p.cz_interaction_GHz = f01_for_11_20(hi.EC, lo.f01);   // |20> of the tuned qubit meets |11>
            {
                const double EJ = solve_EJ_for_f01(p.cz_interaction_GHz, hi.EC);
                p.cz_alpha_GHz = transmon_fa(EJ, hi.EC).alpha;
            }
            p.cz_excursion_MHz = (hi.f01 - p.cz_interaction_GHz) * 1e3;
            const double tcz = p.t_cz_ns * 1e-9;
            if (hiD.tunable) {
                const double EJsum = EJ_from_Rn_GHz(hiD.Rn_kOhm, D.gap_ueV);
                auto f_at = [&](double fl) { return transmon_fa(EJ_at_flux(EJsum, hiD.squid_d, fl), hi.EC).f01; };
                const double f_top = f_at(hiD.flux), f_bot = f_at(0.49);
                if (p.cz_interaction_GHz < f_top && p.cz_interaction_GHz > f_bot) {
                    double blo = std::abs(hiD.flux), bhi = 0.49;
                    for (int b = 0; b < 50; ++b) {
                        const double mid = 0.5 * (blo + bhi);
                        (f_at(mid) > p.cz_interaction_GHz ? blo : bhi) = mid;
                    }
                    const double flux = 0.5 * (blo + bhi);
                    const double h1 = 1e-4;
                    if (flux > 5 * h1) {
                        const double fp = f_at(flux + h1), fm = f_at(flux - h1), f0 = f_at(flux);
                        p.cz_flux = flux;
                        p.cz_df_dflux = (fp - fm) / (2 * h1);
                        const double d2 = (fp - 2 * f0 + fm) / (h1 * h1);
                        const double An = D.flux_noise_uPhi0 * 1e-6;
                        p.cz_flux_err = flux_pulse_infidelity(flux_dephasing_rate(An, p.cz_df_dflux, d2, tcz), tcz);
                    }
                }
            }
            auto rates = [&](const QubitReport& q) { return 1e6 / q.T1_us + (std::isfinite(q.Tphi_us) ? 1e6 / q.Tphi_us : 0.0); };
            p.err_2q = 0.8 * (tcz / 2) * (rates(qa) + rates(qb)) + D.err_2q_control + p.cz_flux_err;
            std::string fl;
            const double win = std::max(0.03, 5 * g);
            if (std::abs(qa.f01 - qb.f01) < win) fl += "DEGENERATE ";
            if (std::abs(qa.f01 + qa.alpha - qb.f01) < std::max(0.02, 4 * g) || std::abs(qb.f01 + qb.alpha - qa.f01) < std::max(0.02, 4 * g)) fl += "11-20-RESONANT ";
            if (!D.tunable_coupler() && std::abs(p.zz_kHz) > 1000) fl += "ZZ>1MHz ";
            if (!hiD.tunable) fl += "CZ-NEEDS-TUNABLE ";
            else if (p.cz_interaction_GHz > hi.f_max + 1e-9) fl += "CZ-INFEASIBLE ";
            p.flags = fl;
            if (!fl.empty()) {
                ++R.n_collisions;
                R.warnings.push_back("pair Q" + std::to_string(p.a) + "-Q" + std::to_string(p.b) + ": " + fl);
            }
        } else if (std::abs(qa.f01 - qb.f01) < 10 * g) {
            // unintended coupling between (near-)degenerate qubits: excitations swap resonantly at idle
            p.flags = "PARASITIC-RESONANT ";
            ++R.n_collisions;
            R.warnings.push_back("pair Q" + std::to_string(p.a) + "-Q" + std::to_string(p.b) + ": PARASITIC-RESONANT (g=" +
                                 std::to_string(p.g_MHz * 1e3).substr(0, 6) + " kHz, detuning " + std::to_string(std::abs(p.delta_MHz)).substr(0, 6) +
                                 " MHz; swap time " + std::to_string(1.0 / (4 * g) / 1e3).substr(0, 5) + " us)");
        }
        R.pairs.push_back(p);
    }

    // --- aggregates
    int n2 = 0;
    for (const auto& q : R.q) { R.mean_err_1q += q.err_1q / n; R.mean_T1_us += q.T1_us / n; R.mean_T2_us += q.T2_us / n; }
    for (const auto& p : R.pairs) {
        if (p.designed) { R.mean_err_2q += p.err_2q; ++n2; }
        R.max_abs_zz_kHz = std::max(R.max_abs_zz_kHz, std::abs(p.zz_idle_kHz));
    }
    if (n2) R.mean_err_2q /= n2;
    return R;
}

// ================================================================== design generator (inverse design)
struct GenSpec {
    std::string topology{"grid"};     // grid | line | ring | surface
    int rows{3}, cols{3};
    int distance{3};                  // surface-code distance (topology "surface")
    double f_hi{5.4}, f_lo{4.8};      // GHz idle frequencies (checkerboard; for "surface": measure / data qubits)
    double stagger{0.06};             // GHz spacing between frequency sub-groups
    int plan{8};                      // 2 = plain checkerboard, 4 = row-parity stagger, 8 = king-graph colouring (every
                                      // same-type qubit within distance 2 detuned: kills parasitic & mediated resonances)
    double alpha{-0.22};              // GHz
    double g{0.012};                  // GHz designed nearest-neighbour coupling
    std::string architecture{"fixed-coupling-cz"};
    double Q{3e6};
    double pitch_um{800};
    double ro_f0{6.9}, ro_step{0.06};
    double flux_noise_uPhi0{2.0};
    // readout: if chi_MHz > 0 the resonator coupling is set per qubit for a full dispersive shift 2chi = chi_MHz
    double kappa_MHz{1.0}, chi_MHz{0.0}, purcell_filter{30};
    int feedlines{1};
};

// Choose shunt capacitances, junction resistances and coupling capacitances so that the full
// capacitance-matrix model realises the target frequencies, anharmonicity and nearest-neighbour coupling.
inline void inverse_design(ChipDesign& D, const std::vector<double>& ftarget, double alpha, double g) {
    const int n = int(D.qubits.size());
    std::vector<double> EC(n), EJ(n), CS(n);
    for (int i = 0; i < n; ++i) {
        auto [ej, ec] = solve_EJ_EC(ftarget[i], alpha);
        EJ[i] = ej; EC[i] = ec; CS[i] = C_from_EC_F(ec);
        D.qubits[i].Rn_kOhm = Rn_from_EJ_kOhm(ej, D.gap_ueV);
    }
    for (auto& c : D.couplings) {
        const int a = D.index_of(c.a), b = D.index_of(c.b);
        const double na = std::abs(transmon_levels(EJ[a], EC[a], 0, 3).n(0, 1)), nb = std::abs(transmon_levels(EJ[b], EC[b], 0, 3).n(0, 1));
        c.C_fF = g * 1e9 * phys::h * CS[a] * CS[b] / (4 * phys::e * phys::e * na * nb) * 1e15;
    }
    for (int i = 0; i < n; ++i) {
        double sumc = 0;
        for (auto& c : D.couplings) if (D.index_of(c.a) == i || D.index_of(c.b) == i) sumc += c.C_fF;
        D.qubits[i].C_fF = CS[i] * 1e15 - sumc;
    }
    for (int it = 0; it < 6; ++it) {
        const auto R = analyze(D);
        for (int i = 0; i < n; ++i) {
            D.qubits[i].C_fF *= R.q[i].EC / EC[i];                 // >1 -> too little capacitance
            D.qubits[i].Rn_kOhm = Rn_from_EJ_kOhm(solve_EJ_for_f01(ftarget[i], R.q[i].EC), D.gap_ueV);
        }
        for (auto& c : D.couplings)
            for (const auto& p : R.pairs)
                if ((p.a == c.a && p.b == c.b) || (p.a == c.b && p.b == c.a)) { c.C_fF *= (g * 1e3) / p.g_MHz; break; }
    }
}

// Resonator coupling for a target full dispersive shift 2chi (transmon formula chi = g^2 alpha / (Delta (Delta + alpha))).
inline void tune_readout(ChipDesign& D, double kappa_MHz, double chi_full_MHz, double purcell_filter) {
    const auto R = analyze(D);
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        auto& ro = D.qubits[i].ro;
        ro.kappa_MHz = kappa_MHz;
        ro.purcell_filter = purcell_filter;
        if (chi_full_MHz <= 0) continue;
        const double dl = R.q[i].f01 - ro.f_GHz, a = R.q[i].alpha, chi = 0.5 * chi_full_MHz * 1e-3;
        ro.g_MHz = std::sqrt(chi * std::abs(dl * (dl + a) / a)) * 1e3;
    }
}

inline ChipDesign generate(const GenSpec& s) {
    ChipDesign D;
    D.architecture = s.architecture;
    D.Q_dielectric = s.Q;
    D.flux_noise_uPhi0 = s.flux_noise_uPhi0;
    std::vector<double> ftarget;
    auto offset4 = [&](int color) { return s.plan == 2 ? 0.0 : (color - 1.5) * s.stagger; };
    if (s.topology == "surface") {
        const auto S = rotated_surface_code(s.distance);
        D.surface_d = s.distance;
        D.name = "surface-code d=" + std::to_string(s.distance) + " (" + std::to_string(S.n_qubits()) + " qubits)";
        D.qubits.resize(std::size_t(S.n_qubits()));
        ftarget.assign(std::size_t(S.n_qubits()), 0.0);
        auto place = [&](int k, int x, int y, const std::string& role, double f) {
            auto& q = D.qubits[std::size_t(k)];
            q.id = k; q.role = role;
            q.x_um = (x + y) * s.pitch_um / 2;                  // rotate 45 deg: couplers become horizontal/vertical
            q.y_um = (y - x + 2 * s.distance) * s.pitch_um / 2;
            ftarget[std::size_t(k)] = f;
        };
        for (const auto& dq : S.data) place(dq.qubit, dq.x, dq.y, "data", s.f_lo + offset4(dq.color));
        for (const auto& st : S.stabs) place(st.qubit, st.x, st.y, st.type == 'X' ? "measure-X" : "measure-Z", s.f_hi + offset4(st.color));
        for (const auto& st : S.stabs) for (int k : st.support) D.couplings.push_back({st.qubit, S.data[std::size_t(k)].qubit, 1.0});
    } else if (s.topology == "grid") {
        D.name = std::to_string(s.rows) + "x" + std::to_string(s.cols) + " grid";
        for (int r = 0; r < s.rows; ++r) for (int c = 0; c < s.cols; ++c) {
            QubitDesign q; q.id = r * s.cols + c; q.x_um = c * s.pitch_um; q.y_um = r * s.pitch_um;
            D.qubits.push_back(q);
            double off = 0;
            if (s.plan == 4) off = ((r % 2) ? 0.5 : -0.5) * s.stagger;
            if (s.plan == 8) {
                const int u = (r + c) >> 1, v = (r - c + 1000) >> 1;      // same-type distance<=2 neighbours = king moves in (u,v)
                off = ((u % 2) + 2 * (v % 2) - 1.5) * s.stagger;
            }
            ftarget.push_back((((r + c) % 2 == 0) ? s.f_hi : s.f_lo) + off);
            if (c + 1 < s.cols) D.couplings.push_back({r * s.cols + c, r * s.cols + c + 1, 1.0});
            if (r + 1 < s.rows) D.couplings.push_back({r * s.cols + c, (r + 1) * s.cols + c, 1.0});
        }
    } else {
        const int N = s.rows * s.cols;
        D.name = std::to_string(N) + "-qubit " + s.topology;
        for (int i = 0; i < N; ++i) {
            QubitDesign q; q.id = i; q.x_um = i * s.pitch_um; D.qubits.push_back(q);
            ftarget.push_back((i % 2 == 0 ? s.f_hi : s.f_lo) + (s.plan == 2 ? 0.0 : (((i / 2) % 2) ? 0.5 : -0.5) * s.stagger));
            if (i + 1 < N) D.couplings.push_back({i, i + 1, 1.0});
        }
        if (s.topology == "ring") {
            if (N % 2) throw std::runtime_error("ring needs an even qubit count for a 2-colour frequency pattern");
            D.couplings.push_back({N - 1, 0, 1.0});
        }
    }
    // readout resonators: feedlines split the chip into horizontal bands; resonators spaced along each feedline
    const int n = int(D.qubits.size()), F = std::max(1, s.feedlines);
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return D.qubits[a].y_um < D.qubits[b].y_um || (D.qubits[a].y_um == D.qubits[b].y_um && D.qubits[a].x_um < D.qubits[b].x_um); });
    const int per = (n + F - 1) / F;
    for (int k = 0; k < n; ++k) {
        auto& q = D.qubits[std::size_t(order[k])];
        q.feedline = k / per;
        q.ro.f_GHz = s.ro_f0 + s.ro_step * (k % per);
    }
    inverse_design(D, ftarget, s.alpha, s.g);
    tune_readout(D, s.kappa_MHz, s.chi_MHz, s.purcell_filter);
    return D;
}

// ================================================================== chip floor-plan (SVG)
inline void write_layout_svg(const ChipDesign& D, const ChipReport& R, const std::string& path, const std::string& subtitle = "") {
    double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
    for (const auto& q : D.qubits) { x0 = std::min(x0, q.x_um); x1 = std::max(x1, q.x_um); y0 = std::min(y0, q.y_um); y1 = std::max(y1, q.y_um); }
    const double span = std::max({x1 - x0, y1 - y0, 1.0});
    const double S = 620.0 / span, pad = 70, top = 90;
    const double legend_w = 60 * 4 + 7.5 * (10 + 22 + 22 + 5) + (D.tunable_coupler() ? 150 : 0);
    const double W = std::max((x1 - x0) * S + 2 * pad, legend_w + 2 * pad), H = (y1 - y0) * S + 2 * pad + top + 70;
    auto X = [&](double x) { return pad + (x - x0) * S; };
    auto Y = [&](double y) { return top + pad + (y - y0) * S; };
    auto color = [](const std::string& r) {
        if (r == "data") return "#2563eb";
        if (r == "measure-X") return "#ea580c";
        if (r == "measure-Z") return "#16a34a";
        return "#6b7280";
    };
    std::ofstream f(path);
    char buf[512];
    std::snprintf(buf, sizeof buf, "<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' viewBox='0 0 %.0f %.0f' font-family='Helvetica,Arial,sans-serif'>\n", W, H, W, H);
    f << buf << "<rect width='100%' height='100%' fill='#ffffff'/>\n";
    f << "<text x='" << pad << "' y='36' font-size='22' font-weight='bold' fill='#111827'>" << D.name << "</text>\n";
    f << "<text x='" << pad << "' y='62' font-size='14' fill='#4b5563'>" << (subtitle.empty() ? D.architecture : subtitle) << "</text>\n";
    for (const auto& c : D.couplings) {
        const auto& a = D.qubits[std::size_t(D.index_of(c.a))];
        const auto& b = D.qubits[std::size_t(D.index_of(c.b))];
        std::snprintf(buf, sizeof buf, "<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' stroke='#9ca3af' stroke-width='%d' stroke-linecap='round'/>\n",
                      X(a.x_um), Y(a.y_um), X(b.x_um), Y(b.y_um), D.tunable_coupler() ? 9 : 4);
        f << buf;
        if (D.tunable_coupler()) {   // coupler element at the midpoint
            std::snprintf(buf, sizeof buf, "<rect x='%.1f' y='%.1f' width='14' height='14' rx='3' fill='#374151'/>\n",
                          0.5 * (X(a.x_um) + X(b.x_um)) - 7, 0.5 * (Y(a.y_um) + Y(b.y_um)) - 7);
            f << buf;
        }
    }
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const auto& q = D.qubits[i];
        std::snprintf(buf, sizeof buf, "<circle cx='%.1f' cy='%.1f' r='27' fill='%s' stroke='#111827' stroke-width='1.5'/>\n", X(q.x_um), Y(q.y_um), color(q.role));
        f << buf;
        std::snprintf(buf, sizeof buf, "<text x='%.1f' y='%.1f' font-size='15' font-weight='bold' fill='#ffffff' text-anchor='middle'>Q%d</text>\n", X(q.x_um), Y(q.y_um) + 5, q.id);
        f << buf;
        if (i < R.q.size()) {
            std::snprintf(buf, sizeof buf, "<text x='%.1f' y='%.1f' font-size='12' fill='#111827' text-anchor='middle'>%.3f GHz</text>\n", X(q.x_um), Y(q.y_um) + 44, R.q[i].f01);
            f << buf;
        }
    }
    const double ly = H - 34;
    const char* roles[4] = {"data", "measure-X", "measure-Z", "qubit"};
    const char* labels[4] = {"data qubit", "X-stabiliser (measure)", "Z-stabiliser (measure)", "qubit"};
    double lx = pad;
    for (int k = 0; k < 4; ++k) {
        bool used = false;
        for (const auto& q : D.qubits) if (q.role == roles[k]) used = true;
        if (!used) continue;
        std::snprintf(buf, sizeof buf, "<circle cx='%.1f' cy='%.1f' r='9' fill='%s'/><text x='%.1f' y='%.1f' font-size='13' fill='#374151'>%s</text>\n",
                      lx + 9, ly, color(roles[k]), lx + 24, ly + 5, labels[k]);
        f << buf;
        lx += 60 + 7.5 * std::strlen(labels[k]);
    }
    if (D.tunable_coupler()) {
        std::snprintf(buf, sizeof buf, "<rect x='%.1f' y='%.1f' width='14' height='14' rx='3' fill='#374151'/><text x='%.1f' y='%.1f' font-size='13' fill='#374151'>tunable coupler</text>\n", lx, ly - 7, lx + 22, ly + 5);
        f << buf;
    }
    f << "</svg>\n";
}

// ================================================================== fabrication Monte Carlo
struct FabResult {
    int samples{0}, collision_free{0};
    std::vector<double> err2q_mean, f01_sigma_MHz, max_zz_kHz;
    std::map<std::string, int> failure_modes;
};
inline FabResult fab_monte_carlo(const ChipDesign& D0, int N, unsigned seed = 1) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0, 1);
    FabResult F; F.samples = N;
    const auto R0 = analyze(D0);
    const double sig = D0.Rn_spread_pct / 100.0;
    for (int k = 0; k < N; ++k) {
        ChipDesign D = D0;
        for (auto& q : D.qubits) {
            q.Rn_kOhm *= (1 + sig * nd(rng));
            if (q.tunable) q.squid_d = std::abs(q.squid_d + sig / std::sqrt(2.0) * nd(rng));
        }
        const auto R = analyze(D);
        if (R.n_collisions == 0) ++F.collision_free;
        for (const auto& p : R.pairs) {
            std::istringstream ss(p.flags); std::string w;
            while (ss >> w) ++F.failure_modes[w];
        }
        double s2 = 0;
        for (std::size_t i = 0; i < R.q.size(); ++i) { const double d = (R.q[i].f01 - R0.q[i].f01) * 1e3; s2 += d * d; }
        F.f01_sigma_MHz.push_back(std::sqrt(s2 / R.q.size()));
        F.err2q_mean.push_back(R.mean_err_2q);
        F.max_zz_kHz.push_back(R.max_abs_zz_kHz);
    }
    return F;
}

} // namespace qpuforge
