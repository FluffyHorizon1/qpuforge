// qpuforge validation suite — every layer checked against independent physics or an independent computation
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>

using namespace qpuforge;
static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// unique scratch path in the system temp directory (no fixed /tmp names: parallel ctest runs and non-POSIX systems)
static std::string temp_path(const std::string& name) {
    static const std::string tag = std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    return (std::filesystem::temp_directory_path() / ("qpuforge_" + tag + "_" + name)).string();
}

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

static void test_coupler_and_flux() {
    std::puts("[coupler] exchange null is not a ZZ null; CZ flux noise uses the short-time 1/f error");
    auto levels = [](double f, double alpha) {
        auto [ej, ec] = solve_EJ_EC(f, alpha);
        return transmon_levels(ej, ec, 0.0, 4);
    };
    const auto A = levels(5.20, -0.22), B = levels(5.90, -0.22);
    const auto good = coupler_exchange_null(A, B, 0.001, 0.040, 0.20);
    std::printf("  g_qq=1 MHz g_qc=40 MHz: null %.3f GHz, J %.3f MHz, ZZ %.2f kHz\n", good.f_GHz, good.J_MHz, good.zz_kHz);
    CHECK(good.ok, "exchange null exists");
    CHECK(std::abs(good.J_MHz) < 0.05, "exchange cancelled");
    CHECK(good.f_GHz > 3.6 && good.f_GHz < 4.1, "null below the data band");
    CHECK(good.zz_kHz < -2 && good.zz_kHz > -20, "ZZ at the null is a few kHz (%.2f)", good.zz_kHz);
    const auto old = coupler_exchange_null(A, B, 0.004, 0.070, 0.20);
    std::printf("  g_qq=4 MHz g_qc=70 MHz: null %.3f GHz, ZZ %.1f kHz (typed residual was 5 kHz)\n", old.f_GHz, old.zz_kHz);
    CHECK(old.ok && old.zz_kHz < -80, "4 MHz direct capacitance does not leave a 5 kHz residual");

    GenSpec s; s.topology = "line"; s.rows = 1; s.cols = 2; s.g = 0.004; s.f_hi = 5.90; s.f_lo = 5.20; s.stagger = 0;
    s.architecture = "tunable-coupler-cz";
    const auto R = analyze(generate(s));
    const auto& p = R.pairs[0];
    std::printf("  CZ flux %.3f Phi0, df/dPhi %.2f GHz/Phi0, flux infidelity %.3e (e_2q %.3e)\n",
                p.cz_flux, p.cz_df_dflux, p.cz_flux_err, p.err_2q);
    CHECK(p.cz_flux > 0.05 && p.cz_flux < 0.45, "CZ sits off the sweet spot");
    CHECK(std::abs(p.cz_df_dflux) > 1.0, "first-order flux sensitivity at the gate");
    CHECK(p.cz_flux_err > 1e-6 && p.cz_flux_err < 5e-4, "short-time flux error");
    const int hi = R.q[0].f01 >= R.q[1].f01 ? 0 : 1;
    CHECK(std::abs(p.cz_df_dflux) > 100 * std::abs(R.q[hi].df_dflux), "idle sweet spot is first-order insensitive");

    const auto RA = analyze(alpha1(3));
    double zlo = 1e99, zhi = -1e99, flo = 1e99, fhi = 0, gdev = 0;
    int nd = 0;
    for (const auto& q : RA.pairs) if (q.designed) {
        ++nd;
        zlo = std::min(zlo, q.zz_idle_kHz); zhi = std::max(zhi, q.zz_idle_kHz);
        flo = std::min(flo, q.coupler_f_GHz); fhi = std::max(fhi, q.coupler_f_GHz);
        gdev = std::max(gdev, std::abs(q.g_MHz - 1.0));
        CHECK(q.coupler_f_GHz > 3.5 && q.coupler_f_GHz < 4.3, "coupler parking %.3f", q.coupler_f_GHz);
        CHECK(std::abs(q.coupler_J_MHz) < 0.1, "exchange cancelled");
        CHECK(std::abs(q.zz_idle_kHz) < 25, "idle ZZ %.2f kHz", q.zz_idle_kHz);
    }
    std::printf("  Alpha-1: %d links, |g-1 MHz| max %.3f, couplers %.3f–%.3f GHz, ZZ %.2f to %.2f kHz, collisions %d\n",
                nd, gdev, flo, fhi, zlo, zhi, RA.n_collisions);
    CHECK(nd == 24 && RA.n_collisions == 0, "Alpha-1 still collision-free");
    CHECK(gdev < 0.05, "direct coupling realised at 1 MHz");
    const auto roundtrip = design_from_json(design_to_json(alpha1(3)));
    CHECK(roundtrip.tc_model && std::abs(roundtrip.tc_coupler_EC - 0.20) < 1e-12 && std::abs(roundtrip.coupler.g_qc_MHz - 40) < 1e-12,
          "three-transmon parameters survive design JSON round-trip");
    CHECK(!alpha1_v11(3).tc_model, "v1.1 pulse coupler retains its existing idle model");
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

static void test_dynamic_circuits() {
    std::puts("[dynamic] mid-circuit measurement, reset, classical record");
    GenSpec s; s.topology = "line"; s.rows = 1; s.cols = 3; s.architecture = "tunable-coupler-cz";
    auto D = generate(s);
    D.reset_err = 0.05;
    const auto R = analyze(D); const auto G = chip_graph(D, R); auto nm = noise_model(D, R);
    const auto C = parse_qasm("OPENQASM 2.0; qreg q[2]; creg c[4]; x q[0]; measure q[0] -> c[0]; reset q[0]; measure q[0] -> c[1];"
                              " h q[1]; cx q[1],q[0]; measure q[1] -> c[2]; measure q[0] -> c[3];");
    const auto P = compile(C, G);
    CHECK(P.dynamic && P.n_measure == 4 && P.n_reset == 1, "dynamic compile");
    RunOptions o; o.trajectories = 2000; o.noisy = false;
    const auto r = run(P, nm, o);
    long good = 0; for (auto& [k, c] : r.counts) if (k == "0001" || k == "1101") good += c;
    std::printf("  noiseless: outcomes {0001, 1101} in %ld / %ld shots, split %ld / %ld\n", good, r.shots, r.counts.count("0001") ? r.counts.at("0001") : 0L, r.counts.count("1101") ? r.counts.at("1101") : 0L);
    CHECK(good == r.shots, "measure/reset/Bell semantics");
    // reset error statistics: only reset noise on; P(read 1 after reset of |1>) = reset_err
    auto nm2 = nm; nm2.only("reset");
    const auto P2 = compile(parse_qasm("OPENQASM 2.0; qreg q[1]; creg c[1]; x q[0]; reset q[0]; measure q[0] -> c[0];"), G);
    RunOptions o2; o2.trajectories = 20000;
    const auto r2 = run(P2, nm2, o2);
    const double p1 = double(r2.counts.count("1") ? r2.counts.at("1") : 0) / r2.shots, sem = std::sqrt(0.05 * 0.95 / r2.shots);
    std::printf("  reset error: measured %.4f, expected %.4f +- %.4f\n", p1, 0.05, sem);
    CHECK(std::abs(p1 - 0.05) < 4 * sem, "reset error rate");
}

static void test_surface_geometry() {
    std::puts("[qec] rotated surface code geometry: stabilisers commute, logicals anticommute, CZ schedule conflict-free");
    for (int d : {3, 5, 7}) {
        const auto S = rotated_surface_code(d);
        int bad = 0, nx = 0;
        for (const auto& a : S.stabs) { if (a.type == 'X') ++nx; for (const auto& b : S.stabs) if (a.type == 'X' && b.type == 'Z') {
            int c = 0; for (int x : a.support) for (int y : b.support) c += (x == y); bad += c % 2; } }
        int lz = 0, lx = 0, ov = 0;
        for (const auto& st : S.stabs) {
            int c = 0;
            for (int x : st.support) for (int y : (st.type == 'X' ? S.logical_z : S.logical_x)) c += (x == y);
            (st.type == 'X' ? lz : lx) += c % 2;
        }
        for (int x : S.logical_z) for (int y : S.logical_x) ov += (x == y);
        int conflicts = 0;
        for (int step = 0; step < 4; ++step) {
            std::vector<int> use(std::size_t(S.n_qubits()), 0);
            for (const auto& st : S.stabs) if (st.order[std::size_t(step)] >= 0) {
                conflicts += use[std::size_t(st.qubit)]++ > 0;
                conflicts += use[std::size_t(S.data[std::size_t(st.order[std::size_t(step)])].qubit)]++ > 0;
            }
        }
        std::printf("  d=%d: %d qubits (%zu data, %zu stabilisers), anticommuting stab pairs %d, logical violations %d/%d, Z_L.X_L overlap %d, schedule conflicts %d\n",
                    d, S.n_qubits(), S.data.size(), S.stabs.size(), bad, lz, lx, ov, conflicts);
        CHECK(S.n_qubits() == 2 * d * d - 1 && nx == (d * d - 1) / 2, "qubit / stabiliser counts");
        CHECK(bad == 0 && lz == 0 && lx == 0 && ov % 2 == 1 && conflicts == 0, "surface code algebra");
    }
}

static ChipDesign surface_chip(int d) {
    GenSpec s; s.topology = "surface"; s.distance = d; s.architecture = "tunable-coupler-cz"; s.f_hi = 5.9; s.f_lo = 5.2; s.g = 0.004;
    s.feedlines = 3; s.kappa_MHz = 4; s.chi_MHz = 4; s.purcell_filter = 100;
    auto D = generate(s);
    D.t_readout_ns = 500; D.active_reset = true; D.reset_err = 0.003; D.t_reset_ns = 200;
    return D;
}

static void test_qec_circuits() {
    std::puts("[qec] memory experiment: noiseless determinism, detector error model, circuit distance");
    {
        const auto D = surface_chip(3); const auto R = analyze(D); const auto nm = noise_model(D, R);
        const auto M = prepare_memory(D, R, nm, 3);
        const auto recs = run_records(M.prog, nm, 30, 5, false);
        int nonzero = 0;
        for (const auto& rc : recs) {
            for (const auto& det : M.exp.detectors) { int v = 0; for (int c : det) v ^= rc[std::size_t(c)]; nonzero += v; }
            int o = 0; for (int c : M.exp.observable) o ^= rc[std::size_t(c)]; nonzero += o;
        }
        std::printf("  d=3 R=3 noiseless statevector: %d detection events / logical flips over 30 shots (%.0f ns per round)\n", nonzero, M.prog.duration_ns / 3);
        CHECK(nonzero == 0, "noiseless memory experiment is deterministic");
        CHECK(M.dem.hyperedge_faults == 0 && M.dem.undetectable_logical == 0, "DEM is graphlike with no undetectable logicals");
        DetectorErrorModel U = M.dem; for (auto& e : U.edges) e.p = 0.01;
        MatchingDecoder dec(U);
        const auto faults = enumerate_faults(M.frame, M.exp);
        int wrong = 0; std::vector<double> c; std::vector<int> ch;
        for (const auto& f : faults) wrong += dec.decode(f.dets, c, ch) != int(f.obs);
        std::printf("  d=3: %zu single faults, %d uncorrected by uniform-weight matching (circuit distance 3 needs 0)\n", faults.size(), wrong);
        CHECK(wrong == 0, "every single fault correctable at d=3");
    }
    {
        const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
        const auto M = prepare_memory(D, R, nm, 2);
        DetectorErrorModel U = M.dem; for (auto& e : U.edges) e.p = 0.01;
        MatchingDecoder dec(U);
        const auto faults = enumerate_faults(M.frame, M.exp);
        std::mt19937 rng(3); std::uniform_int_distribution<std::size_t> pick(0, faults.size() - 1);
        int wrong1 = 0, wrong2 = 0; std::vector<double> c; std::vector<int> ch;
        for (const auto& f : faults) wrong1 += dec.decode(f.dets, c, ch) != int(f.obs);
        for (int k = 0; k < 20000; ++k) {
            const auto& f = faults[pick(rng)]; const auto& g = faults[pick(rng)];
            std::vector<int> dets;
            std::set_symmetric_difference(f.dets.begin(), f.dets.end(), g.dets.begin(), g.dets.end(), std::back_inserter(dets));
            wrong2 += dec.decode(dets, c, ch) != int(f.obs ^ g.obs);
        }
        std::printf("  d=5: %zu single faults (%d uncorrected), 20000 random fault pairs (%d uncorrected); hyperedges %d\n", faults.size(), wrong1, wrong2, M.dem.hyperedge_faults);
        CHECK(wrong1 == 0 && wrong2 == 0 && M.dem.undetectable_logical == 0, "every pair of faults correctable at d=5");
    }
}

static void test_frame_vs_statevector() {
    std::puts("[qec] Pauli-frame noise model vs exact statevector trajectories (amplitude damping only, exaggerated)");
    GenSpec s; s.topology = "surface"; s.distance = 3; s.architecture = "tunable-coupler-cz"; s.f_hi = 5.9; s.f_lo = 5.2; s.g = 0.004; s.Q = 3e5;
    auto D = generate(s); D.t_readout_ns = 500; D.active_reset = true; D.t_reset_ns = 200;
    const auto R = analyze(D); auto nm = noise_model(D, R); nm.only("T1");
    const auto M = prepare_memory(D, R, nm, 2);
    const auto fr = run_memory_frame(M, 200000, 5);
    // 1000 trajectories: P_L standard error ~0.005 (was 240 shots, +-0.012, which could not resolve a 2x model error).
    // The frame model Pauli-twirls amplitude damping, so a small systematic offset (~0.003 in detection rate) is allowed.
    const long nsv = 1000;
    const auto sv = run_memory_statevector(M, nm, nsv, 9);
    const double sem = std::sqrt(sv.detection_rate * (1 - sv.detection_rate) / (double(nsv) * double(M.exp.detectors.size())));
    const double sem_pl = std::hypot(sv.p_logical_sem, fr.p_logical_sem);
    std::printf("  detection rate: statevector %.4f +- %.4f | Pauli frame %.4f ; P_L statevector %.4f +- %.4f | frame %.4f +- %.4f\n",
                sv.detection_rate, sem, fr.detection_rate, sv.p_logical, sv.p_logical_sem, fr.p_logical, fr.p_logical_sem);
    CHECK(std::abs(sv.detection_rate - fr.detection_rate) < 3 * sem + 0.004, "frame model detection rate");
    CHECK(std::abs(sv.p_logical - fr.p_logical) < 3 * sem_pl + 0.003, "frame model logical error |%.4f - %.4f| vs %.4f", sv.p_logical, fr.p_logical,
          3 * sem_pl + 0.003);
}

static void test_eigensolver() {
    std::puts("[eigensolver] Householder + implicit QL vs cyclic Jacobi (random symmetric matrices)");
    std::mt19937 rng(41); std::normal_distribution<double> nd(0, 1);
    double worst_e = 0, worst_r = 0, worst_o = 0;
    for (int n : {1, 2, 3, 7, 16, 36, 48}) {
        std::vector<double> A(std::size_t(n) * n);
        for (int i = 0; i < n; ++i) for (int j = 0; j <= i; ++j) A[std::size_t(i) * n + j] = A[std::size_t(j) * n + i] = nd(rng);
        std::vector<double> V, W;
        const auto e1 = sym_eigen_tri(A, n, V); const auto e2 = sym_eigen(A, n, W);
        for (int i = 0; i < n; ++i) worst_e = std::max(worst_e, std::abs(e1[std::size_t(i)] - e2[std::size_t(i)]));
        for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
            double s = 0; for (int k = 0; k < n; ++k) s += A[std::size_t(i) * n + k] * V[std::size_t(k) * n + j];
            worst_r = std::max(worst_r, std::abs(s - e1[std::size_t(j)] * V[std::size_t(i) * n + j]));
        }
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) {
            double s = 0; for (int k = 0; k < n; ++k) s += V[std::size_t(k) * n + a] * V[std::size_t(k) * n + b];
            worst_o = std::max(worst_o, std::abs(s - (a == b ? 1.0 : 0.0)));
        }
    }
    std::printf("  max |eig diff| %.1e, residual %.1e, orthogonality %.1e\n", worst_e, worst_r, worst_o);
    CHECK(worst_e < 1e-11 && worst_r < 1e-11 && worst_o < 1e-12, "eigensolver accuracy");
}

static void test_pulse_cz() {
    std::puts("[pulse] tunable-coupler CZ: gauge continuity, idle evolution = ZZ, |11>-|20> gap vs perturbation theory, calibration");
    // smooth sign gauge of the charge matrix elements across a flux sweep
    {
        const auto T = FluxTransmon::from_target(5.99, -0.22);
        double jump = 0;
        auto prev = element_spectrum(T, 0.0, 4);
        for (int k = 1; k <= 160; ++k) {
            const auto cur = element_spectrum(T, 0.002 * k, 4);
            for (std::size_t i = 0; i < cur.N.size(); ++i) jump = std::max(jump, std::abs(cur.N[i] - prev.N[i]));
            prev = cur;
        }
        std::printf("  charge matrix elements over a 0..0.32 Phi0 sweep: max step change %.4f (no sign flips)\n", jump);
        CHECK(jump < 0.02, "gauge discontinuity %.3f", jump);
    }
    CouplerElement ce; ce.physical = true;
    double g12 = 0, gap0 = 0;
    const auto S = design_pair_system(5.99, -0.22, 5.29, -0.22, ce, &g12, &gap0);
    std::printf("  pair 5.99/5.29 GHz, g_qc %.0f MHz: designed g12 %.2f MHz, residual parked |11>-|20> gap %.4f MHz, ZZ at park %+.2f kHz\n",
                ce.g_qc_MHz, g12 * 1e3, gap0 * 1e3, S.zz(S.fluxc_idle) * 1e6);
    CHECK(gap0 < 2e-5, "parked residual exchange not nulled: %.3e GHz", gap0);
    CHECK(std::abs(S.zz(S.fluxc_idle)) < 1e-5, "ZZ at park above 10 kHz");
    // idle evolution reproduces the static ZZ phase exactly, with no leakage
    {
        CZPulse P; P.flux1_int = S.flux1_idle; P.fluxc_on = S.fluxc_idle; P.t_rise = 5; P.t_flat = 40;
        const auto r = simulate_cz(S, P);
        double ph = r.cphase; if (ph > PI) ph -= 2 * PI;
        const double expect = -2 * PI * S.zz(S.fluxc_idle) * P.duration();   // e^{-i 2 pi E t}: conditional phase -2 pi ZZ t
        std::printf("  idle %.0f ns: leakage %.1e, conditional phase %.3e rad vs -2 pi ZZ t = %.3e rad\n", P.duration(), r.leakage, ph, expect);
        CHECK(r.leakage < 1e-10, "idle leakage");
        CHECK(std::abs(ph - expect) < 2e-6, "idle conditional phase vs static ZZ");
    }
    // avoided crossing with the coupler on vs second-order perturbation theory: gap = 2 sqrt2 |g12 + g1c g2c / (f2 - fc)|
    {
        const double xc = S.c.flux_for_f01(6.1);
        const auto [x1, gap] = cz_avoided_crossing(S, xc);
        const auto s1 = element_spectrum(S.q1, x1, 2), sc = element_spectrum(S.c, xc, 2), s2 = element_spectrum(S.q2, 0.0, 2);
        const double g1c = S.G1c * s1.n(0, 1) * sc.n(0, 1), g2c = S.G2c * s2.n(0, 1) * sc.n(0, 1), g12e = S.G12 * s1.n(0, 1) * s2.n(0, 1);
        const double fc = S.c.f01(xc), f2 = S.q2.f01(0.0);
        const double pert = 2 * std::sqrt(2.0) * std::abs(g12e + g1c * g2c / (f2 - fc));
        // the neglected next order is O((g/Delta)^2) relative; allow 1.5x the larger of the two dispersive ratios (was a flat 20 %)
        const double f1 = S.q1.f01(x1), nlo = std::max(std::pow(g1c / (f1 - fc), 2), std::pow(g2c / (f2 - fc), 2));
        std::printf("  coupler on at %.2f GHz: |11>-|20> gap %.2f MHz, perturbative %.2f MHz (%+.1f%%; next order ~%.1f%%)\n", fc, gap * 1e3, pert * 1e3,
                    100 * (gap / pert - 1), 100 * nlo);
        CHECK(std::abs(gap - pert) < 1.5 * nlo * pert, "gap vs perturbation theory");
    }
    // calibration of the 3 knobs
    const auto [P, R] = optimize_cz(S, S.c.flux_for_f01(ce.f_on_GHz), ce.t_rise_ns);
    const auto Rf = simulate_cz(S, P, 0.05);
    std::printf("  calibrated CZ: %.1f ns, 1-F %.2e, leakage %.2e, phase error %+.1e rad; dt 0.1 -> 0.05 ns changes F by %.1e\n", R.duration_ns,
                1 - R.fidelity, R.leakage, R.cphase - PI, std::abs(Rf.fidelity - R.fidelity));
    CHECK(1 - R.fidelity < 3e-4, "CZ coherent error %.2e", 1 - R.fidelity);
    CHECK(R.leakage < 2e-4, "CZ leakage %.2e", R.leakage);
    CHECK(std::abs(R.cphase - PI) < 5e-3, "CZ conditional phase");
    CHECK(R.duration_ns > 30 && R.duration_ns < 60, "CZ duration %.1f", R.duration_ns);
    CHECK(std::abs(Rf.fidelity - R.fidelity) < 2e-6, "time-step convergence");
    // level-truncation convergence: the calibrated pulse re-simulated with one more coupler / qubit level
    // (Lc = 3 was not converged: it left ~4e-3 rad of conditional-phase error, hence the default Lc = 4)
    {
        CZSystem Tc = S; Tc.Lc = S.Lc + 1;
        CZSystem Tq = S; Tq.L1 = S.L1 + 1; Tq.L2 = S.L2 + 1;
        const auto rc = simulate_cz(Tc, P), rq = simulate_cz(Tq, P);
        std::printf("  truncation L%d/%d/%d: Lc+1 changes F by %.1e, phase by %+.1e rad; L1,L2+1 changes F by %.1e, phase by %+.1e rad\n",
                    S.L1, S.Lc, S.L2, rc.fidelity - R.fidelity, rc.cphase - R.cphase, rq.fidelity - R.fidelity, rq.cphase - R.cphase);
        CHECK(S.Lc >= 4, "default coupler truncation Lc = %d", S.Lc);
        CHECK(std::abs(rc.fidelity - R.fidelity) < 2e-6 && std::abs(rc.cphase - R.cphase) < 5e-4, "coupler-level truncation not converged");
        CHECK(std::abs(rq.fidelity - R.fidelity) < 2e-6 && std::abs(rq.cphase - R.cphase) < 5e-4, "qubit-level truncation not converged");
    }
}

static void test_leakage_frame() {
    std::puts("[leakage] frame model: inert when off, Markov steady state, leakage-removing reset, DQLR");
    const auto D = surface_chip(3); const auto R = analyze(D);
    // (a) leakage enabled with zero rates is bit-identical to no leakage
    {
        auto nm0 = noise_model(D, R); nm0.leakage = false;
        auto nm1 = noise_model(D, R); nm1.leakage = true; nm1.leak_1q = 0; nm1.leak_cz.clear();
        const auto M0 = prepare_memory(D, R, nm0, 4), M1 = prepare_memory(D, R, nm1, 4);
        const auto r0 = run_memory_frame(M0, 20000, 7), r1 = run_memory_frame(M1, 20000, 7);
        CHECK(r0.logical_errors == r1.logical_errors && r0.detection_rate == r1.detection_rate, "zero-rate leakage changed results");
    }
    // (b) single-qubit-gate leakage + seepage only: time-averaged leaked fraction vs the two-state Markov solution
    {
        auto nm = noise_model(D, R); nm.only("leakage"); nm.leak_1q = 2e-3; nm.leak_transport = 0; nm.leak_cz.clear();
        const auto M = prepare_memory(D, R, nm, 25);
        const auto r = run_memory_frame(M, 20000, 11);
        const auto& F = M.frame;
        double pred_d = 0, pred_m = 0; int nd = 0, nmq = 0;
        for (int q = 0; q < F.n; ++q) {
            int gates = 0;
            for (const auto& fm : F.moments) for (const auto& o : fm.ops) if ((o.op == Op::SX || o.op == Op::X) && o.a == q) ++gates;
            const double lam = gates * nm.leak_1q / F.duration_ns, gam = F.seep_rate[std::size_t(q)], k = lam + gam, T = F.duration_ns;
            const double avg = lam / k * (1 - (1 - std::exp(-k * T)) / (k * T));
            if (M.is_data[std::size_t(q)]) { pred_d += avg; ++nd; } else { pred_m += avg; ++nmq; }
        }
        pred_d /= nd; pred_m /= nmq;
        std::printf("  gate leakage 2e-3 + seepage 2/T1, 25 rounds: leaked fraction data %.4f (Markov %.4f), measure %.4f (Markov %.4f)\n",
                    r.leak_frac_data, pred_d, r.leak_frac_meas, pred_m);
        CHECK(std::abs(r.leak_frac_data - pred_d) < 0.05 * pred_d, "data leakage vs Markov");
        CHECK(std::abs(r.leak_frac_meas - pred_m) < 0.05 * pred_m, "measure leakage vs Markov");
    }
    // (c) CZ leakage of the measure qubits: a leakage-removing reset clears it every round; (d) DQLR clears data leakage
    {
        auto base = noise_model(D, R); base.leak_1q = 0; base.leak_transport = 0.2;
        for (const auto& p : R.pairs) if (p.designed) {
            const int a = D.index_of(p.a), b = D.index_of(p.b);
            const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
            base.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
        }
        auto run = [&](bool lrr, bool dqlr) {
            auto nm = base; nm.reset_removes_leakage = lrr; nm.dqlr = dqlr; nm.dqlr_eff = 0.98; nm.dqlr_err = 0;
            return run_memory_frame(prepare_memory(D, R, nm, 15), 20000, 13);
        };
        const auto none = run(false, false), lrr = run(true, false), both = run(true, true);
        std::printf("  CZ leakage 2e-3 (15 rounds): leaked measure fraction %.1e -> %.1e with leakage-removing reset (%.1e with DQLR too);"
                    " data %.1e -> %.1e with DQLR; detection %.4f / %.4f / %.4f\n",
                    none.leak_frac_meas, lrr.leak_frac_meas, both.leak_frac_meas, lrr.leak_frac_data, both.leak_frac_data, none.detection_rate,
                    lrr.detection_rate, both.detection_rate);
        // observed ratios: 0.19 (LRR, measure), 0.29 (DQLR, data), 0.63 (detection); bounds leave ~30 % margin
        CHECK(lrr.leak_frac_meas < 0.25 * none.leak_frac_meas, "LRR does not remove measure-qubit leakage");
        CHECK(both.leak_frac_data < 0.4 * lrr.leak_frac_data, "DQLR does not remove data leakage");
        CHECK(lrr.detection_rate < 0.8 * none.detection_rate, "leakage removal should reduce detection events");
        CHECK(both.leak_frac_meas <= 1.2 * lrr.leak_frac_meas, "DQLR should not add measure-qubit leakage");
    }
}

static void test_calibrated_noise_plumbing() {
    std::puts("[plumbing] per-pair CZ calibration reaches the noise model and the Pauli-frame program");
    const auto D = surface_chip(3); auto R = analyze(D);
    int k0 = -1;
    for (std::size_t k = 0; k < R.pairs.size(); ++k) if (R.pairs[k].designed) { k0 = int(k); break; }
    auto& p = R.pairs[std::size_t(k0)];
    p.err_2q_ctrl = 0.0; p.leak_2q = 1e-3; p.leak_qubit = p.a;
    const auto nm = noise_model(D, R);
    const int a = D.index_of(p.a), b = D.index_of(p.b);
    CHECK(nm.err_2q_pair.count({std::min(a, b), std::max(a, b)}) == 1 && nm.err_2q(a, b) == 0.0, "per-pair error missing");
    CHECK(nm.leak_cz.count({std::min(a, b), std::max(a, b)}) == 1 && nm.leak_cz.at({std::min(a, b), std::max(a, b)}).second == 1e-3, "leak map");
    const auto M = prepare_memory(D, R, nm, 2);
    const auto sim = sim_index(M.prog);
    int zero = 0, other = 0, leaky = 0;
    for (const auto& fm : M.frame.moments) {
        for (const auto& c : fm.ctrl2) {
            const bool mine = (c.a == sim[std::size_t(a)] && c.b == sim[std::size_t(b)]) || (c.a == sim[std::size_t(b)] && c.b == sim[std::size_t(a)]);
            if (mine) zero += c.p == 0.0; else other += std::abs(c.p - 1.25 * D.err_2q_control) < 1e-15;
        }
        for (const auto& o : fm.ops) if (o.op == Op::CZ && o.pleak > 0) leaky += o.mover == sim[std::size_t(a)];
    }
    std::printf("  overridden pair: %d CZs with zero control error, %d other CZs at the default, %d CZs carrying leakage\n", zero, other, leaky);
    CHECK(zero == 2 && other == M.prog.n_cz - 2 && leaky == 2, "frame program per-pair parameters (2 rounds: 2 CZs per pair)");
}

static void test_layout() {
    std::puts("[layout] MoM capacitance vs exact results, geometry checks, width/space DRC, GDSII round trip, chip assembly");
    using namespace layout;
    {   // unit square plate: C = 0.3667874 * 4 pi eps0 a (Read 1997)
        Conductor sq; sq.name = "plate"; sq.rects = {{0, 0, 1000, 1000}};
        const auto C = extract_capacitance({sq}, 1.0, MeshOptions{1.0, 100, 1.6});
        const double ref = 0.3667874 * 4 * PI * phys::eps0 * 1e-3 * 1e15;
        std::printf("  square plate 1 mm: %.3f fF vs %.3f fF (%+.2f%%)\n", C.at(0, 0), ref, 100 * (C.at(0, 0) / ref - 1));
        CHECK(std::abs(C.at(0, 0) / ref - 1) < 0.003, "square plate capacitance");
    }
    {   // CPW on silicon: per-unit-length capacitance vs conformal mapping (finite grounds)
        const double w = 10, s = 6, g = 150;
        auto cpw = [&](double L) {
            Conductor c; c.name = "c"; c.rects = {{-L / 2, -w / 2, L / 2, w / 2}};
            Conductor G; G.name = "g"; G.ground = true; G.outline = {-L / 2 - s - g, -w / 2 - s - g, L / 2 + s + g, w / 2 + s + g}; G.holes = {{-L / 2 - s, -w / 2 - s, L / 2 + s, w / 2 + s}};
            return extract_capacitance({c, G}, 11.45, MeshOptions{0.5, 40, 1.6}).at(0, 0);
        };
        const double per_um = (cpw(800) - cpw(400)) / 400, ref = cpw_capacitance_per_m(w * 1e-6, s * 1e-6, g * 1e-6, 11.45) * 1e9;
        std::printf("  CPW 10/6 um on Si: %.4f fF/um vs conformal mapping %.4f fF/um (%+.2f%%)\n", per_um, ref, 100 * (per_um / ref - 1));
        CHECK(std::abs(per_um / ref - 1) < 0.025, "CPW capacitance per length");
    }
    {   // touching conductors are rejected
        Conductor a; a.name = "a"; a.rects = {{0, 0, 10, 10}};
        Conductor b; b.name = "b"; b.rects = {{5, 5, 20, 20}};
        CHECK(!check_geometry({a, b}).empty(), "overlap not detected");
    }
    {   // merged-region width/space check on a synthetic layer with known violations
        const std::vector<Rect> rs{{0, 0, 10, 10}, {10.5, 0, 20, 10},     // 0.5 um gap: 1 space marker
                                   {30, 0, 30.3, 10},                     // 0.3 um sliver: 1 width marker
                                   {40, 0, 50, 10}, {45, 0, 55, 10},      // overlapping rects merge to 15 um: clean
                                   {60, 0, 70, 10}, {70, 0, 80, 10},      // abutting rects merge: no false 0-width gap
                                   {100, 20, 100.5, 1000}};               // 0.5 um line through 4 tiles: 1 marker
        const auto m = detail::width_space_check(rs, 1.0, 1.0, "t");
        int nw = 0, ns = 0; for (const auto& x : m) (x.rule.find("width") != std::string::npos ? nw : ns)++;
        std::printf("  width/space check on a synthetic layer: %d width, %d space markers (expected 2, 1)\n", nw, ns);
        CHECK(nw == 2 && ns == 1, "width/space markers %d/%d", nw, ns);
    }
    {   // chip assembly: placement checks clean, resonators at quarter-wave length, GDSII parses back
        const auto D = alpha1_v11(3);
        CellSet cs; cs.meas.L = 129.3; cs.data.L = 131.9; cs.cg.pad = 252.6; cs.cg.l_c = 43.5; cs.cg.l_b = 0.1; cs.meas.r_pad = cs.data.r_pad = 40.9;
        cs.meas.T = cs.data.T = 120; cs.cg.y1 = 32.5;
        const auto L = assemble_chip(D, cs);
        double worst = 0; for (const auto& r : L.resonators) worst = std::max(worst, std::abs(r.drawn_um - r.target_um));
        CHECK(L.problems.empty(), "chip placement/spacing: %s", L.problems.empty() ? "" : L.problems.front().c_str());
        std::printf("  AF-17 mask DRC (not asserted; first-cut layout): %s\n", drc_summary(L).c_str());
        CHECK(worst < 1.0 && L.resonators.size() == 17, "resonator lengths");
        const std::string gds_path = temp_path("af17.gds");
        write_gds(L, gds_path, "T");
        std::ifstream f(gds_path, std::ios::binary);
        std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        std::size_t i = 0; int boundaries = 0, records = 0; bool ok = true; int last = -1;
        while (i + 4 <= b.size()) {
            const int n = (b[i] << 8) | b[i + 1];
            if (n < 4 || i + std::size_t(n) > b.size()) { ok = false; break; }
            last = b[i + 2]; if (last == 0x08) ++boundaries;
            ++records; i += std::size_t(n);
        }
        const std::size_t expect = L.etch.size() + L.jj.size() + L.bumps.size() + 1;
        std::printf("  AF-17 layout: %zu etch polygons, placement %s; GDSII: %d records, %d boundaries (expected %zu), ends with ENDLIB: %s\n", L.etch.size(),
                    L.problems.empty() ? "clean" : "violations", records, boundaries, expect, last == 0x04 ? "yes" : "no");
        CHECK(ok && last == 0x04 && std::size_t(boundaries) == expect, "GDSII stream structure");
    }
}


// ---------------------------------------------------------------- json_rfc: RFC 8259 parser/writer (q10; named json_rfc so it does not collide with the json group on test/qpuforge-json-gds-roundtrip)
static bool json_throws(const std::string& text) {
    try { Json::parse(text); } catch (const std::runtime_error&) { return true; } return false;
}
static std::string json_dump(const Json& j, int indent = 0) {
    std::ostringstream o; j.dump(o, indent); return o.str();
}
static void test_json_rfc() {
    std::puts("[json_rfc] RFC 8259: well-formed numbers/strings/escapes; reject nan/inf/hex/leading + or dot; dump escapes controls");
    const Json v = Json::parse("{\"n\": -1.5e2, \"i\": 42, \"b\": true, \"s\": \"a\\\"b\\\\c\\nd\\te\", \"a\": [1, {\"k\": []}]}");
    CHECK(v["n"].num() == -150.0 && v["i"].num() == 42 && v["b"].boolean() && v["s"].str() == std::string("a\"b\\c\nd\te"), "well-formed scalars");
    CHECK(v["a"].size() == 2 && v["a"][1]["k"].size() == 0, "nested containers");
    CHECK(Json::parse("\"\\u0041\"").str() == "A", "\\u0041 decodes to A");
    CHECK(Json::parse("\"\\u00e9\"").str() == "\xc3\xa9", "\\u00e9 decodes to UTF-8");
    CHECK(Json::parse("[0, 0.5, 1e-3, -2]").size() == 4, "canonical numbers");
    const char* bad[] = {
        "{\"a\": nan}", "{\"a\": inf}", "[0x10]", "[+5]", "[.5]", "[1e999]",
        "\"\\u0041", "\"\\q\"", "\"a\\qb\"", "\"\t\"",
        "[1.]", "[01]", "[-]", "[1e]", "[+1]", "NaN", "Infinity", "-Infinity",
        "\"\\u12\"", "\"\\uGGGG\""
    };
    int accepted = 0;
    for (const char* b : bad) if (!json_throws(b)) { ++accepted; std::printf("  accepted non-JSON: %s\n", b); }
    CHECK(accepted == 0, "%d non-RFC documents were accepted", accepted);

    Json w = Json::array();
    w.push(std::string("a\tb\rc") + char(1) + "d");
    const std::string dumped = json_dump(w, 0);
    CHECK(dumped.find('\t') == std::string::npos && dumped.find('\r') == std::string::npos && dumped.find('\x01') == std::string::npos,
          "dump must escape controls, got %s", dumped.c_str());
    CHECK(dumped.find("\\t") != std::string::npos && dumped.find("\\r") != std::string::npos && dumped.find("\\u0001") != std::string::npos,
          "dump should use \\t \\r \\u0001: %s", dumped.c_str());
    const Json back = Json::parse(dumped);
    CHECK(back[0].str() == std::string("a\tb\rc") + char(1) + "d", "round-trip of escaped controls");
    // documented writer behaviour: non-finite numbers become null (and that is legal JSON)
    Json nf = Json::array(); nf.push(std::nan("")); nf.push(HUGE_VAL);
    CHECK(json_dump(nf, 0) == "[null, null]", "non-finite dump: %s", json_dump(nf, 0).c_str());
}

static const std::vector<std::pair<const char*, void (*)()>> kTests = {
    {"decomposition", test_decomposition},
    {"transmon", test_transmon},
    {"coupler_and_flux", test_coupler_and_flux},
    {"zz_and_inverse_design", test_zz_and_inverse_design},
    {"qasm", test_qasm},
    {"compiler", test_compiler},
    {"decoherence", test_decoherence},
    {"exchange", test_exchange},
    {"trajectories_vs_exact", test_trajectories_vs_exact},
    {"ghz_metric", test_ghz_metric},
    {"xeb", test_xeb},
    {"dynamic_circuits", test_dynamic_circuits},
    {"surface_geometry", test_surface_geometry},
    {"qec_circuits", test_qec_circuits},
    {"frame_vs_statevector", test_frame_vs_statevector},
    {"eigensolver", test_eigensolver},
    {"pulse_cz", test_pulse_cz},
    {"leakage_frame", test_leakage_frame},
    {"calibrated_noise_plumbing", test_calibrated_noise_plumbing},
    {"layout", test_layout},
    {"json_rfc", test_json_rfc},
};

// test_qpuforge            run every test
// test_qpuforge NAME...    run the named tests (CTest runs one test per entry, see CMakeLists.txt)
// test_qpuforge --list     print the test names
int main(int argc, char** argv) {
    std::vector<std::string> want(argv + 1, argv + argc);
    if (want.size() == 1 && want[0] == "--list") { for (const auto& t : kTests) std::puts(t.first); return 0; }
    int ran = 0;
    for (const auto& [name, fn] : kTests) {
        if (!want.empty() && std::find(want.begin(), want.end(), name) == want.end()) continue;
        const auto t0 = std::chrono::steady_clock::now();
        fn(); ++ran;
        std::printf("  (%s: %.1f s)\n", name, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    if (ran == 0 || (!want.empty() && ran != int(want.size()))) { std::printf("unknown test name(s)\n"); return 2; }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
