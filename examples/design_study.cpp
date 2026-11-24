// Design study: how do chip-design choices change whole-processor performance?
//   A. coupling strength g on a fixed-coupling CZ chip (faster gates vs. always-on ZZ crosstalk)
//   B. tunable-coupler architecture vs. the best fixed coupling
//   C. frequency plan: naive 2-frequency checkerboard vs. staggered 4-frequency plan (parasitic resonances)
//   D. materials: dielectric quality factor (T1) -> algorithm fidelity
// Output: results/design_study.csv + console summary.
#include "qpuforge/qpuforge.hpp"
#include <cstdio>
#include <fstream>

using namespace qpuforge;

struct Metrics { double t_cz, zz, e2q, T1, ghz, ghz_sem, xeb_f, qft, n_coll; };

static Metrics evaluate(const ChipDesign& D, int traj) {
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    const auto nm = noise_model(D, R);
    Metrics M{};
    double tcz = 0, zz = 0; int np = 0;
    for (const auto& p : R.pairs) if (p.designed) { tcz += p.t_cz_ns; zz += std::abs(p.zz_idle_kHz); ++np; }
    M.t_cz = tcz / np; M.zz = zz / np; M.e2q = R.mean_err_2q; M.T1 = R.mean_T1_us; M.n_coll = R.n_collisions;
    RunOptions o; o.trajectories = traj; o.shots_per_traj = 0; o.seed = 2026;
    const int n = int(D.qubits.size());
    {   // GHZ across the whole chip — experimental definition (populations + parity-oscillation amplitude)
        const auto P = compile(ghz(n), G);
        std::vector<int> lq(n); std::iota(lq.begin(), lq.end(), 0);
        M.ghz = run_ghz_fidelity(P, nm, lq, o).fidelity;
        M.ghz_sem = 0;
    }
    {   // random circuits native to the coupler graph, 10 cycles, averaged over 3 instances
        std::vector<std::pair<int, int>> edges;
        for (const auto& p : R.pairs) if (p.designed) edges.push_back({D.index_of(p.a), D.index_of(p.b)});
        CompileOptions co; co.layout = "trivial";
        double f = 0;
        for (int s = 0; s < 3; ++s) f += run(compile(random_xeb(n, edges, 10, 500 + s), G, co), nm, o).fidelity / 3;
        M.xeb_f = f;
    }
    {   // QFT on 6 qubits: routing-heavy algorithm (long circuit -> idle errors & crosstalk matter)
        M.qft = run(compile(qft(6), G), nm, o).fidelity;
    }
    return M;
}

int main(int argc, char** argv) {
    const int traj = argc > 1 ? std::atoi(argv[1]) : 300;
    std::ofstream csv("results/design_study.csv");
    csv << "study,variant,param,t_cz_ns,mean_abs_zz_kHz,mean_e2q,mean_T1_us,collisions,ghz9_fidelity,ghz9_sem,rcs10_fidelity,qft6_fidelity\n";
    auto emit = [&](const char* study, const std::string& variant, double param, const Metrics& m) {
        csv << study << "," << variant << "," << param << "," << m.t_cz << "," << m.zz << "," << m.e2q << "," << m.T1 << "," << m.n_coll << ","
            << m.ghz << "," << m.ghz_sem << "," << m.xeb_f << "," << m.qft << "\n";
        std::printf("  %-26s t_CZ %6.1f ns  |ZZ| %7.1f kHz  e2q %.2e  T1 %6.1f us  coll %2.0f |  GHZ-9 %.3f  RCS-10 %.3f  QFT-6 %.3f\n",
                    variant.c_str(), m.t_cz, m.zz, m.e2q, m.T1, m.n_coll, m.ghz, m.xeb_f, m.qft);
        std::fflush(stdout);
    };

    GenSpec base; base.rows = 3; base.cols = 3;     // default: 8-frequency plan

    std::puts("A. fixed-coupling CZ: coupling strength sweep (3x3, 8-frequency plan, Q=3e6)");
    for (double g : {2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0}) {
        GenSpec s = base; s.g = g * 1e-3;
        emit("A", "fixed g=" + std::to_string(int(g)) + " MHz", g, evaluate(generate(s), traj));
    }

    std::puts("B. tunable-coupler CZ (g_max 40 MHz, residual ZZ 5 kHz)");
    {
        GenSpec s = base; s.architecture = "tunable-coupler-cz"; s.g = 0.004;
        auto D = generate(s);
        emit("B", "tunable coupler", 40, evaluate(D, traj));
        D.tc_residual_zz_kHz = 50;
        emit("B", "tunable coupler, ZZ 50kHz", 40, evaluate(D, traj));
    }

    std::puts("C. frequency plan (fixed coupling g=4 MHz)");
    for (int plan : {2, 4, 8}) {
        GenSpec s = base; s.g = 0.004; s.plan = plan;
        emit("C", std::to_string(plan) + "-frequency plan", plan, evaluate(generate(s), traj));
    }

    std::puts("D. materials: dielectric Q (tunable coupler)");
    for (double Q : {5e5, 1e6, 3e6, 1e7}) {
        GenSpec s = base; s.architecture = "tunable-coupler-cz"; s.g = 0.004; s.Q = Q;
        char buf[64]; std::snprintf(buf, sizeof buf, "Q=%.0e", Q);
        emit("D", buf, Q, evaluate(generate(s), traj));
    }
    std::puts("wrote results/design_study.csv");
    return 0;
}
