#include "qpuforge/qpuforge.hpp"
#include <cstdio>
using namespace qpuforge;
int main() {
    const auto D = alpha1_v11(3); const auto R = analyze(D);
    for (std::size_t i = 0; i < D.qubits.size(); ++i) std::printf("Q%d %s f01 %.3f GHz EJ %.2f GHz T1 %.1f us\n", D.qubits[i].id, D.qubits[i].role.c_str(), R.q[i].f01, R.q[i].EJ, R.q[i].T1_us);
    // XY-line T1 limit of Q6 from the window extraction (C_sigma 93.8 fF): own XY6 (direct 0.0342 + via C6-3 -0.0497 + via C6-9 -0.0462 fF)
    const double w = 2 * PI * 5.930e9, Cs = 93.8e-15;
    for (double c : {0.0342, 0.0342 - 0.0497 - 0.0462}) std::printf("C_eff %.4f fF -> T1 limit %.0f us\n", c, Cs / (w * w * c * 1e-15 * c * 1e-15 * 50) * 1e6);
    // Rabi drive power at the chip for a 25 MHz Rabi rate: P = Omega^2 hbar w / (4 Gamma1)
    for (double c : {0.03, 0.0102}) { const double w1 = 2 * PI * 5.11e9, G = w1 * w1 * c * 1e-15 * c * 1e-15 * 50 / 91.5e-15, Om = 2 * PI * 25e6;
        const double P = Om * Om * 1.054571817e-34 * w1 / (4 * G); std::printf("C_eff %.4f fF at 5.11 GHz: T1 %.2f ms; 25 MHz Rabi needs %.1f dBm at the chip\n", c, 1e3 / G, 10 * std::log10(P / 1e-3)); }
}
