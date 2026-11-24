// qpuforge validation suite — every layer checked against independent physics or an independent computation
#include "qpuforge/qpuforge.hpp"
#include <cstdio>

using namespace qpuforge;
static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static M2 random_unitary(std::mt19937& rng) {
    std::uniform_real_distribution<double> U(0, 2 * PI);
    const double th = std::acos(1 - 2 * std::uniform_real_distribution<double>(0, 1)(rng));
    return scaled(gates::U3(th, U(rng), U(rng)), std::exp(cplx(0, U(rng))));
}

static void test_decomposition() {
    std::puts("[decomposition] arbitrary 1q unitary -> RZ/SX/X native sequence");
    std::mt19937 rng(3);
    std::vector<M2> cases = {gates::I(), gates::X(), gates::Y(), gates::Z(), gates::H(), gates::S(), gates::T(), gates::SX(), gates::SXdg(),
                             gates::RY(PI / 2), gates::RY(PI), gates::RX(PI / 2), gates::RX(-PI / 2), gates::RZ(1.234), gates::SW(), gates::SY()};
    for (int k = 0; k < 3000; ++k) cases.push_back(random_unitary(rng));
    double worst = 0; int pulses_max = 0;
    for (const auto& U : cases) {
        const auto seq = compile_detail::native_1q(0, U);
        M2 P = gates::I(); int pulses = 0;
        for (const auto& g : seq) { P = g.u * P; if (g.op != Op::RZ) ++pulses; }
        if (seq.empty()) P = gates::I();
        worst = std::max(worst, phase_distance(P, U));
        pulses_max = std::max(pulses_max, pulses);
    }
    std::printf("  %zu unitaries, worst deviation %.2e, max physical pulses %d\n", cases.size(), worst, pulses_max);
    CHECK(worst < 1e-9, "decomposition error %.3e", worst);
    CHECK(pulses_max <= 2, "more than 2 pulses");
}

static void test_transmon() {
    std::puts("[transmon] charge-basis spectrum vs Koch et al. (2007) asymptotics");
    const double EC = 0.25;
    for (double ratio : {30.0, 50.0, 80.0}) {
        const double EJ = ratio * EC;
        const auto T = transmon_levels(EJ, EC, 0.0, 4);
        const double f_asym = std::sqrt(8 * EJ * EC) - EC;
        // charge dispersion of the 0-1 transition: eps_m = (-1)^m EC 2^(4m+5)/m! sqrt(2/pi) (EJ/2EC)^(m/2+3/4) exp(-sqrt(8EJ/EC))
        auto eps = [&](int m) { double f = 1; for (int k = 2; k <= m; ++k) f *= k;
            return std::pow(-1.0, m) * EC * std::pow(2.0, 4 * m + 5) / f * std::sqrt(2 / PI) * std::pow(EJ / (2 * EC), m / 2.0 + 0.75) * std::exp(-std::sqrt(8 * EJ / EC)); };
        const double disp_asym = std::abs(eps(1) - eps(0));
        const double disp = std::abs(transmon_levels(EJ, EC, 0.5, 3).f01() - T.f01());
        const double n01_asym = std::pow(EJ / (8 * EC), 0.25) / std::sqrt(2.0);
        std::printf("  EJ/EC=%4.0f  f01 %.4f (asym %.4f)  alpha %.4f (-EC=%.3f)  disp %.3e (asym %.3e)  |n01| %.4f (asym %.4f)\n",
                    ratio, T.f01(), f_asym, T.alpha(), -EC, disp, disp_asym, std::abs(T.n(0, 1)), n01_asym);
        CHECK(std::abs(T.f01() - f_asym) / f_asym < 0.01, "f01 vs asymptotic");
        // independent reference: numpy diagonalisation on a 2000-point phase grid (see README)
        const double alpha_ref = ratio == 30 ? -0.3091 : (ratio == 50 ? -0.2873 : -0.2774);
        const double f01_ref = ratio == 30 ? 3.6034 : (ratio == 50 ? 4.7355 : 6.0634);
        CHECK(std::abs(T.alpha() - alpha_ref) < 2e-4, "alpha vs phase-grid reference");
        CHECK(std::abs(T.f01() - f01_ref) < 2e-4, "f01 vs phase-grid reference");
        CHECK(std::abs(disp - disp_asym) / disp_asym < 0.25, "charge dispersion vs asymptotic (%.3e vs %.3e)", disp, disp_asym);
        CHECK(std::abs(std::abs(T.n(0, 1)) - n01_asym) / n01_asym < 0.05, "n01");
    }
    const double EJ = EJ_from_Rn_GHz(10.0, 180.0);
    std::printf("  E_J(Rn=10 kOhm, Delta=180 ueV) = %.3f GHz (hand calc 14.04)\n", EJ);
    CHECK(std::abs(EJ - 14.04) < 0.05, "Ambegaokar-Baratoff E_J");
    CHECK(std::abs(EC_from_C_GHz(80e-15) - 0.2421) < 0.001, "E_C(80 fF)");
}

static void test_zz_and_inverse_design() {
    std::puts("[pair physics] exact two-transmon ZZ vs 2nd-order perturbation theory; inverse design hits targets");
    for (double g : {0.003, 0.006, 0.012}) {
        GenSpec s; s.topology = "line"; s.rows = 1; s.cols = 2; s.g = g; s.f_hi = 5.4; s.f_lo = 4.8; s.stagger = 0;
        const auto D = generate(s);
        const auto R = analyze(D);
        const auto& p = R.pairs[0];
        std::printf("  g=%5.1f MHz: realised g %.3f MHz, ZZ exact %8.2f kHz, perturbative %8.2f kHz, t_CZ %.1f ns\n", g * 1e3, p.g_MHz, p.zz_kHz, p.zz_pert_kHz, p.t_cz_ns);
        CHECK(std::abs(p.g_MHz - g * 1e3) < 0.05, "inverse design coupling");
        CHECK(std::abs(p.zz_kHz - p.zz_pert_kHz) / std::abs(p.zz_pert_kHz) < 0.15, "ZZ exact vs perturbative");
        CHECK(std::abs(R.q[0].f01 - 5.4) < 1e-3 && std::abs(R.q[1].f01 - 4.8) < 1e-3, "inverse design f01 (%.4f, %.4f)", R.q[0].f01, R.q[1].f01);
        CHECK(std::abs(R.q[0].alpha + 0.22) < 2e-3, "inverse design alpha %.4f", R.q[0].alpha);
    }
    // ZZ scales as g^2
    GenSpec a; a.topology = "line"; a.rows = 1; a.cols = 2; a.g = 0.004; a.stagger = 0;
    GenSpec b = a; b.g = 0.008;
    const double r = analyze(generate(b)).pairs[0].zz_kHz / analyze(generate(a)).pairs[0].zz_kHz;
    std::printf("  ZZ(2g)/ZZ(g) = %.3f (expect ~4)\n", r);
    CHECK(std::abs(r - 4) < 0.2, "ZZ ~ g^2");
}

static ChipDesign small_chip(int rows, int cols) {
    GenSpec s; s.rows = rows; s.cols = cols; s.g = 0.012;
    return generate(s);
}

static void test_compiler() {
    std::puts("[compiler] routed + decomposed + scheduled circuits reproduce the logical state exactly");
    std::mt19937 rng(11);
    for (auto [rows, cols, nl] : std::vector<std::tuple<int, int, int>>{{2, 3, 6}, {3, 3, 5}, {3, 3, 9}, {1, 5, 5}}) {
        const auto D = small_chip(rows, cols);
        const auto R = analyze(D);
        const auto G = chip_graph(D, R);
        double worst = 0; int swaps = 0;
        for (int trial = 0; trial < 12; ++trial) {
            Circuit C; C.n = nl;
            for (int k = 0; k < 40; ++k) {
                const int t = int(rng() % 6), a = int(rng() % nl); int b = int(rng() % nl); while (b == a) b = int(rng() % nl);
                int c = int(rng() % nl); while (c == a || c == b) c = int(rng() % nl);
                if (t == 0) C.u(a, random_unitary(rng));
                else if (t == 1) C.cx(a, b);
                else if (t == 2) C.cp(a, b, 0.3 + k * 0.1);
                else if (t == 3) C.swap(a, b);
                else if (t == 4 && nl >= 3) C.gates.push_back({Op::CCX, {a, b, c}, {}, 0, -1, "ccx"});
                else C.cz(a, b);
            }
            C.measure_all();
            const auto ref = simulate_logical(C);
            const auto P = compile(C, G);
            swaps += P.n_swaps;
            const auto sim = sim_index(P);
            StateVector out(int(P.active.size()));
            std::mt19937_64 r2(1);
            execute(out, P, sim, NoiseModel{}, false, r2);
            // map logical basis index -> simulated basis index via the final layout
            cplx ov = 0;
            for (std::size_t x = 0; x < ref.dim(); ++x) {
                std::size_t y = 0;
                for (int l = 0; l < nl; ++l) if ((x >> l) & 1) y |= std::size_t(1) << sim[P.layout_final[l]];
                ov += std::conj(ref.amps()[x]) * out.amps()[y];
            }
            worst = std::max(worst, 1 - std::norm(ov));
        }
        std::printf("  %dx%d chip, %d logical qubits: worst infidelity %.2e over 12 random circuits (%d SWAPs inserted)\n", rows, cols, nl, worst, swaps);
        CHECK(worst < 1e-10, "compiled circuit infidelity %.3e", worst);
    }
}

static void test_qasm() {
    std::puts("[qasm] OpenQASM 2.0 reader");
    const char* src = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[3]; creg c[3];
h q[0];          // comment
cx q[0],q[1];
cp(pi/2) q[1],q[2];
u3(pi/2, 0, pi) q[2];
rz(-pi/4 + 2*pi) q[0];
barrier q;
measure q -> c;
)";
    const auto C = parse_qasm(src);
    CHECK(C.n == 3 && C.nclbits == 3, "registers");
    int nm = 0; for (auto& g : C.gates) nm += g.op == Op::MEASURE;
    CHECK(nm == 3, "broadcast measure");
    const auto s = simulate_logical(C);
    // reference: built by hand
    Circuit R; R.n = 3; R.h(0); R.cx(0, 1); R.cp(1, 2, PI / 2); R.u(2, gates::U3(PI / 2, 0, PI)); R.u(0, gates::RZ(-PI / 4 + 2 * PI));
    const auto r = simulate_logical(R);
    CHECK(std::abs(std::norm(s.inner(r)) - 1) < 1e-12, "qasm semantics");
}

// Hand-built schedule: one gate followed by idle moments.
static Compiled single_qubit_program(Op op, const M2& u, int idle_moments, double idle_ns, double t1q) {
    Compiled P; P.n_phys = 1; P.active = {0}; P.layout_initial = P.layout_final = {0};
    Moment m; m.ops.push_back({op, {0}, u, 0, -1, "g"}); m.duration_ns = t1q; P.moments.push_back(m);
    for (int k = 0; k < idle_moments; ++k) { Moment e; e.duration_ns = idle_ns; P.moments.push_back(e); }
    return P;
}

static void test_decoherence() {
    std::puts("[noise] T1 / T2 decay, thermal equilibrium: density matrix exact, trajectories statistical");
    NoiseModel nm;
    nm.T1_ns = {40000}; nm.Tphi_ns = {60000}; nm.p_thermal = {0.012}; nm.P1_init = {0.012}; nm.ro01 = {0}; nm.ro10 = {0};
    nm.control = false; nm.zz = false;
    const double t1q = 25, idle = 1000; const int K = 60;
    // --- T1: X then idle
    {
        auto P = single_qubit_program(Op::X, gates::X(), K, idle, t1q);
        std::mt19937_64 rng(1);
        DensityMatrix rho(1); rho.set_product_thermal(nm.P1_init);
        const auto sim = sim_index(P);
        execute(rho, P, sim, nm, true, rng);
        const double t = t1q + K * idle;
        const double p_exact = nm.p_thermal[0] + (1 - nm.P1_init[0] - nm.p_thermal[0]) * std::exp(-t / nm.T1_ns[0]);
        double acc = 0; const int T = 4000;
        for (int k = 0; k < T; ++k) {
            StateVector s(1); if (std::uniform_real_distribution<double>(0, 1)(rng) < nm.P1_init[0]) s.set_basis(1);
            execute(s, P, sim, nm, true, rng); acc += s.prob1(0);
        }
        const double ptraj = acc / T, sem = std::sqrt(ptraj * (1 - ptraj) / T);
        std::printf("  T1:   P1(t=%.0f us) analytic %.5f | density %.5f | trajectories %.4f +- %.4f\n", t / 1e3, p_exact, rho.prob1(0), ptraj, sem);
        CHECK(std::abs(rho.prob1(0) - p_exact) < 1e-12, "density T1");
        CHECK(std::abs(ptraj - p_exact) < 4 * sem + 1e-3, "trajectory T1");
    }
    // --- T2 (Ramsey): SX then idle; coherence decays with 1/T2 = 1/(2 T1) + 1/Tphi
    {
        auto P = single_qubit_program(Op::SX, gates::SX(), K, idle, t1q);
        std::mt19937_64 rng(2);
        DensityMatrix rho(1); rho.set_product_thermal(nm.P1_init);
        execute(rho, P, sim_index(P), nm, true, rng);
        const double t = t1q + K * idle;
        const double T2 = 1 / (1 / (2 * nm.T1_ns[0]) + 1 / nm.Tphi_ns[0]);
        const double c_exact = 0.5 * (1 - 2 * nm.P1_init[0]) * std::exp(-t / T2);
        std::printf("  T2:   |rho01|(t=%.0f us) analytic %.6f | density %.6f   (T2 = %.1f us)\n", t / 1e3, c_exact, std::abs(rho.coherence(0)), T2 / 1e3);
        CHECK(std::abs(std::abs(rho.coherence(0)) - c_exact) < 1e-12, "density T2");
    }
    // --- long time: relaxes to the thermal population
    {
        auto P = single_qubit_program(Op::X, gates::X(), 50, 100000, t1q);
        std::mt19937_64 rng(3);
        DensityMatrix rho(1); rho.set_product_thermal(nm.P1_init);
        execute(rho, P, sim_index(P), nm, true, rng);
        std::printf("  equilibrium P1 %.5f (thermal %.5f)\n", rho.prob1(0), nm.p_thermal[0]);
        CHECK(std::abs(rho.prob1(0) - nm.p_thermal[0]) < 1e-9, "thermal equilibrium");
    }
}

static void test_exchange() {
    std::puts("[noise] resonant / detuned exchange: moment-by-moment propagators compose to the exact Rabi formula");
    for (double det_MHz : {0.0, 0.8, 3.0}) {
        NoiseModel nm;
        nm.T1_ns = {1e300, 1e300}; nm.Tphi_ns = {1e300, 1e300}; nm.p_thermal = {0, 0}; nm.P1_init = {0, 0}; nm.ro01 = {0, 0}; nm.ro10 = {0, 0};
        nm.T1 = nm.Tphi = nm.thermal = nm.zz = nm.control = false;
        const double g = 0.5e-3;   // GHz
        nm.exchange[{0, 1}] = {g, 5.0 + det_MHz * 1e-3, 5.0};
        Compiled P; P.n_phys = 2; P.active = {0, 1};
        Moment m0; m0.ops.push_back({Op::X, {0}, gates::X(), 0, -1, "x"}); m0.duration_ns = 25; P.moments.push_back(m0);
        std::mt19937 rng(5); std::uniform_real_distribution<double> U(10, 400);
        double t = 25;
        for (int k = 0; k < 40; ++k) { Moment e; e.duration_ns = U(rng); t += e.duration_ns; P.moments.push_back(e); }
        StateVector s(2); s.reset(); std::mt19937_64 r2(1);
        execute(s, P, sim_index(P), nm, true, r2);
        const double d = 0.5 * det_MHz * 1e-3, Om = std::sqrt(d * d + g * g);
        const double expect = (g * g / (Om * Om)) * std::pow(std::sin(2 * PI * Om * t), 2);
        DensityMatrix rho(2); execute(rho, P, sim_index(P), nm, true, r2);
        std::printf("  detuning %.1f MHz, g 0.5 MHz, t=%.0f ns: P(transfer) emulator %.6f | density %.6f | analytic %.6f\n", det_MHz, t, s.prob1(1), rho.prob1(1), expect);
        CHECK(std::abs(s.prob1(1) - expect) < 1e-9, "exchange vs Rabi formula");
        CHECK(std::abs(rho.prob1(1) - expect) < 1e-9, "density exchange vs Rabi formula");
    }
}

static void test_trajectories_vs_exact() {
    std::puts("[noise] full noise model: trajectory fidelity vs exact density matrix");
    auto D = small_chip(2, 3);
    D.Q_dielectric = 3e5;             // exaggerate errors so the comparison is sensitive
    D.err_2q_control = 0.01; D.err_1q_control = 0.002;
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    auto nm = noise_model(D, R);
    for (auto& [k, v] : nm.zz_GHz) v *= 20;   // strong crosstalk
    for (const auto& C : {ghz(5), qft(4)}) {
        const auto P = compile(C, G);
        const double Fdm = run_density_fidelity(P, nm);
        RunOptions o; o.trajectories = 3000; o.shots_per_traj = 0;
        const auto r = run(P, nm, o);
        std::printf("  %-6s on %zu sim qubits, %d moments: density %.5f | trajectories %.5f +- %.5f\n", C.name.c_str(), P.active.size(), int(P.moments.size()), Fdm, r.fidelity, r.fidelity_sem);
        CHECK(std::abs(Fdm - r.fidelity) < 4 * r.fidelity_sem + 2e-3, "trajectory vs density fidelity");
    }
}

static void test_ghz_metric() {
    std::puts("[benchmark] experimental GHZ fidelity = max over relative phase of <GHZ_phi|rho|GHZ_phi>");
    auto D = small_chip(2, 2);
    D.Q_dielectric = 5e5;
    const auto R = analyze(D);
    auto nm = noise_model(D, R);
    for (auto& [k, v] : nm.zz_GHz) v *= 30;   // big coherent phase offset
    const auto P = compile(ghz(4), chip_graph(D, R));
    RunOptions o; o.trajectories = 4000;
    const auto g = run_ghz_fidelity(P, nm, {0, 1, 2, 3}, o);
    // exact reference from the density matrix
    const auto sim = sim_index(P);
    std::mt19937_64 rng(1);
    DensityMatrix rho(int(P.active.size()));
    std::vector<double> p1; for (int q : P.active) p1.push_back(nm.P1_init[q]);
    rho.set_product_thermal(p1);
    execute(rho, P, sim, nm, true, rng);
    std::size_t ones = 0; for (int l = 0; l < 4; ++l) ones |= std::size_t(1) << sim[P.layout_final[l]];
    double best = 0;
    for (int k = 0; k < 3600; ++k) {
        const double phi = 2 * PI * k / 3600;
        StateVector s(int(P.active.size())); auto& a = s.amps(); std::fill(a.begin(), a.end(), cplx(0));
        a[0] = 1 / std::sqrt(2.0); a[ones] = std::exp(cplx(0, phi)) / std::sqrt(2.0);
        best = std::max(best, rho.fidelity(s));
    }
    std::printf("  GHZ-4 with strong ZZ: trajectories %.4f | exact max_phi %.4f\n", g.fidelity, best);
    CHECK(std::abs(g.fidelity - best) < 0.01, "GHZ fidelity metric");
}

static void test_xeb() {
    std::puts("[benchmark] XEB estimator vs exact density-matrix XEB; XEB-vs-fidelity relation (reported)");
    auto D = small_chip(2, 3);
    D.Q_dielectric = 1e6;
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    auto nm = noise_model(D, R);
    nm.readout = false;
    std::vector<std::pair<int, int>> edges;
    for (const auto& p : R.pairs) if (p.designed) edges.push_back({D.index_of(p.a), D.index_of(p.b)});
    CompileOptions co; co.layout = "trivial";
    // (a) estimator correctness: sampled linear XEB must converge to sum_x p_ideal(x) rho_xx computed exactly
    {
        const auto C = random_xeb(6, edges, 12, 42);
        const auto P = compile(C, G, co);
        CHECK(P.n_swaps == 0, "XEB circuits must be native to the coupler graph");
        const auto sim = sim_index(P);
        std::mt19937_64 rng(3);
        StateVector ideal(int(P.active.size())); execute(ideal, P, sim, nm, false, rng);
        DensityMatrix rho(int(P.active.size()));
        std::vector<double> p1; for (int q : P.active) p1.push_back(nm.P1_init[q]);
        rho.set_product_thermal(p1);
        execute(rho, P, sim, nm, true, rng);
        double s = 0; for (std::size_t x = 0; x < ideal.dim(); ++x) s += std::norm(ideal.amps()[x]) * rho.at(x, x).real();
        const double xeb_exact = double(ideal.dim()) * s - 1;
        RunOptions o; o.trajectories = 3000; o.shots_per_traj = 20; o.seed = 9;
        const auto r = run(P, nm, o);
        std::printf("  12-cycle circuit: XEB exact (density matrix) %.4f | sampled %.4f   [state fidelity %.4f]\n", xeb_exact, r.xeb, rho.fidelity(ideal));
        CHECK(std::abs(r.xeb - xeb_exact) < 0.03, "XEB estimator");
    }
    // (b) physics: how well does normalised XEB track fidelity on this noise model (average over circuits)
    for (int cycles : {6, 14, 20}) {
        double F = 0, X = 0; const int inst = 8;
        for (int seed = 0; seed < inst; ++seed) {
            const auto P = compile(random_xeb(6, edges, cycles, 100 + seed), G, co);
            RunOptions o; o.trajectories = 300; o.shots_per_traj = 40; o.seed = 77 + seed;
            const auto r = run(P, nm, o);
            F += r.fidelity / inst; X += r.xeb_norm / inst;
        }
        std::printf("  %2d cycles (avg of %d circuits): state fidelity %.4f   normalised XEB %.4f\n", cycles, inst, F, X);
    }
}

int main() {
    test_decomposition();
    test_transmon();
    test_zz_and_inverse_design();
    test_qasm();
    test_compiler();
    test_decoherence();
    test_exchange();
    test_trajectories_vs_exact();
    test_ghz_metric();
    test_xeb();
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
