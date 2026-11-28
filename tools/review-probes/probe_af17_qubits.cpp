#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
using namespace qpuforge;
int main() {
    const ChipDesign D = alpha1_v11(3);
    const auto R = analyze(D);
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const auto& q = D.qubits[i];
        std::printf("Q%-2d %-10s (%5.0f,%5.0f) f01 %.3f GHz T1 %.1f us; couplings:", q.id, q.role.c_str(), q.x_um, q.y_um, R.q[i].f01, R.q[i].T1_us);
        for (const auto& c : D.couplings) if (c.a == q.id || c.b == q.id) {
            const int o = c.a == q.id ? c.b : c.a;
            for (const auto& p : D.qubits) if (p.id == o) std::printf(" Q%d(%+.0f,%+.0f)", o, p.x_um - q.x_um, p.y_um - q.y_um);
        }
        std::printf("\n");
    }
    std::printf("coupler f_max %.3f GHz\n", D.coupler.f_max_GHz);
}
