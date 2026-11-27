// qpuforge — named chip designs
//
//   gen0_surface17 : the best fixed-coupling design from the 3x3 study, laid out as a distance-3 surface code
//   alpha1         : Alpha-1, the first error-correction-ready design (decisions and rationale in ALPHA-1.md)
//   alpha1_v11     : Alpha-1 v1.1 (AF-17), physical couplers and leakage model
//   beta1          : Beta-1 (BF-49), distance 5
//   gamma1         : Gamma-1 (GF-97), distance 7
#pragma once
#include "device.hpp"
#include <numeric>

namespace qpuforge {

inline GenSpec surface_spec(int d) {
    GenSpec s;
    s.topology = "surface"; s.distance = d;
    s.f_hi = 5.90;            // measure qubits (they move down to the CZ interaction point)
    s.f_lo = 5.20;            // data qubits: lower frequency -> longer Q-limited T1 while idling through readout
    s.stagger = 0.06; s.plan = 8;
    s.alpha = -0.22;
    s.pitch_um = 900;
    s.feedlines = 3; s.ro_f0 = 6.90; s.ro_step = 0.06;
    return s;
}

// Generation 0: fixed capacitive coupling (g = 3 MHz, the optimum of the 3x3 coupling sweep), slow readout,
// measurement-based reset.
inline ChipDesign gen0_surface17(int d = 3) {
    GenSpec s = surface_spec(d);
    s.architecture = "fixed-coupling-cz"; s.g = 0.003;
    s.kappa_MHz = 1.0; s.chi_MHz = 1.0; s.purcell_filter = 30;
    ChipDesign D = generate(s);
    D.name = "Gen-0 surface-" + std::to_string(2 * d * d - 1) + " (fixed coupling)";
    D.t_readout_ns = 800; D.readout_sep_err = 0.01;
    D.active_reset = false; D.t_reset_ns = 300; D.reset_err = 0.01;
    return D;
}

// Alpha-1 chip with a given bare direct coupling g (GHz) per data-measure link. alpha1() uses 1 MHz (PR #2: the coupler
// transmon supplies the gate coupling); alpha1_v11() and beta1() keep the 4 MHz (0.14 fF) their v0.3 outputs were built with.
inline ChipDesign alpha1_with_coupling(int d, double g_GHz) {
    GenSpec s = surface_spec(d);
    s.architecture = "tunable-coupler-cz";
    s.g = g_GHz;
    s.kappa_MHz = 4.0; s.chi_MHz = 4.0;            // 2chi = kappa: optimal readout SNR
    s.purcell_filter = 100;                        // dedicated band-pass Purcell filter per feedline
    ChipDesign D = generate(s);
    D.name = d == 3 ? "Alpha-1 (AF-17)" : "Alpha-1 extended to d=" + std::to_string(d) + " (" + std::to_string(2 * d * d - 1) + " qubits)";
    D.tc_g_max_MHz = 40; D.tc_residual_zz_kHz = 5;
    D.coupler.g_qc_MHz = 40; D.tc_coupler_EC = 0.20; D.tc_model = true; D.cz_flux_noise = true;
    D.t_readout_ns = 500; D.readout_sep_err = 0.005;
    D.active_reset = true; D.t_reset_ns = 200; D.reset_err = 0.003;
    D.err_1q_control = 1e-4; D.err_2q_control = 2e-3;
    return D;
}

// Alpha-1: tunable couplers, fast Purcell-filtered readout, active reset.
inline ChipDesign alpha1(int d = 3) {
    return alpha1_with_coupling(d, 0.001);         // direct capacitance only; the coupler supplies the gate coupling
}

// Alpha-1 v1.1: the same chip with the abstract couplers replaced by physical coupler transmons (calibrated pulse-level
// CZ, per-pair direct coupling for an exchange-free park point at the coupler sweet spot) and a leakage model.
inline ChipDesign alpha1_v11(int d = 3) {
    ChipDesign D = alpha1_with_coupling(d, 0.004); // bare g = 4 MHz (0.14 fF), as before PR #2; pulse.hpp sets g12 per pair
    D.name = d == 3 ? "Alpha-1 v1.1 (AF-17)" : "Alpha-1 v1.1 extended to d=" + std::to_string(d);
    D.tc_model = false; D.cz_flux_noise = false;  // pulse.hpp calibration has its own flux-noise terms
    D.coupler.physical = true;
    D.coupler.f_max_GHz = 8.5; D.coupler.alpha_GHz = -0.13; D.coupler.g_qc_MHz = 160;
    D.coupler.f_on_GHz = 6.1; D.coupler.t_rise_ns = 7; D.coupler.residual_ctrl_err = 1e-3;
    D.leakage.model = true;
    D.leakage.reset_removes_leakage = true;    // multi-level (|2> -> |0>) measure-qubit reset
    D.leakage.dqlr = false;
    return D;
}

// Readout plan: at most `per_line` resonators on a feedline, `step_GHz` apart from `f0_GHz`. Feedlines take the qubits
// in horizontal bands, as generate() does. Call tune_readout() afterwards: the couplings depend on the detunings.
inline void plan_readout(ChipDesign& D, int per_line, double f0_GHz = 6.90, double step_GHz = 0.06) {
    const int n = int(D.qubits.size());
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const auto& qa = D.qubits[std::size_t(a)]; const auto& qb = D.qubits[std::size_t(b)];
        return qa.y_um < qb.y_um || (qa.y_um == qb.y_um && qa.x_um < qb.x_um);
    });
    const int lines = (n + per_line - 1) / std::max(per_line, 1), per = (n + lines - 1) / std::max(lines, 1);
    for (int k = 0; k < n; ++k) {
        auto& q = D.qubits[std::size_t(order[std::size_t(k)])];
        q.feedline = k / per;
        q.ro.f_GHz = f0_GHz + step_GHz * (k % per);
    }
}

// Beta-1 (BF-49): distance-5 surface-code processor. Alpha-1 v1.1 architecture plus the upgrades the Alpha-1 reviews
// ranked highest:
//   materials        dielectric loss halved (Q 3e6 -> 6e6: tantalum base layer, improved substrate cleaning)
//   thermalisation   effective qubit temperature 45 -> 35 mK (IR-tight shielding, better line filtering)
//   readout          300 ns (kappa = 2 chi = 8 MHz) behind individual Purcell filters (suppression 1000): the rev-A review
//                    showed that the stronger readout coupling behind the per-feedline filters (100) made Purcell decay
//                    the T1 limit and cancelled the materials upgrade
//   CZ               residual control error 1e-3 -> 5e-4 (flux-pulse predistortion, closed-loop recalibration)
//   leakage          multi-level (leakage-removing) reset; DQLR is kept as an option only: the rev-A review showed it
//                    costs more than it saves at d = 5 over 50 rounds
//   fabrication      junction trimming (laser anneal) to 0.5% resistance spread
//   readout plan     (rev C) at most 9 resonators per feedline, 6.90-7.38 GHz. Three feedlines for 49 qubits meant 17
//                    resonators per line across 1 GHz: beyond the 6-8 per line that readout chains are run at today
inline ChipDesign beta1(int d = 5) {
    ChipDesign D = alpha1_v11(d);
    D.name = d == 5 ? "Beta-1 (BF-49)" : "Beta-1 physics at d=" + std::to_string(d);
    D.Q_dielectric = 6e6;
    D.T_eff_mK = 35;
    D.t_readout_ns = 300;
    plan_readout(D, 9);
    tune_readout(D, 8.0, 8.0, 1000);
    D.coupler.residual_ctrl_err = 5e-4;
    D.leakage.reset_removes_leakage = true;
    D.leakage.dqlr = false;
    D.Rn_spread_pct = 0.5;
    return D;
}

// Gamma-1 (GF-97): distance-7 surface-code memory, 97 qubits. Beta-1 physics and cells unchanged; the readout plan
// scales with the chip (11 feedlines of up to 9 resonators).
inline ChipDesign gamma1(int d = 7) {
    ChipDesign D = beta1(d);
    D.name = d == 7 ? "Gamma-1 (GF-97)" : "Gamma-1 physics at d=" + std::to_string(d);
    return D;
}

} // namespace qpuforge
