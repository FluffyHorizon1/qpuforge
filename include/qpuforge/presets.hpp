// qpuforge — named chip designs
//
//   gen0_surface17 : the best fixed-coupling design from the 3x3 study, laid out as a distance-3 surface code
//   alpha1         : Alpha-1, the first error-correction-ready design (decisions and rationale in ALPHA-1.md)
#pragma once
#include "device.hpp"

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

// Alpha-1: tunable couplers, fast Purcell-filtered readout, active reset.
inline ChipDesign alpha1(int d = 3) {
    GenSpec s = surface_spec(d);
    s.architecture = "tunable-coupler-cz";
    s.g = 0.001;                                   // direct capacitance only. The coupler, not this cap, does the gate.
    s.kappa_MHz = 4.0; s.chi_MHz = 4.0;            // 2chi = kappa: optimal readout SNR
    s.purcell_filter = 100;                        // dedicated band-pass Purcell filter per feedline
    ChipDesign D = generate(s);
    D.name = d == 3 ? "Alpha-1 (AF-17)" : "Alpha-1 extended to d=" + std::to_string(d) + " (" + std::to_string(2 * d * d - 1) + " qubits)";
    D.tc_g_max_MHz = 40; D.tc_residual_zz_kHz = 5;
    D.tc_g_qc_MHz = 40; D.tc_coupler_EC = 0.20; D.tc_model = true;
    D.t_readout_ns = 500; D.readout_sep_err = 0.005;
    D.active_reset = true; D.t_reset_ns = 200; D.reset_err = 0.003;
    D.err_1q_control = 1e-4; D.err_2q_control = 2e-3;
    return D;
}

} // namespace qpuforge
