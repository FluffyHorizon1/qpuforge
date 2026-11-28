// qpuforge validation suite — every layer checked against independent physics or an independent computation
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include "qpuforge/gds.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numeric>
#include <set>
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
    // |11> = |20> must use alpha at the tuned EJ, not the idle alpha.
    {
        const double EC = 0.20, f_partner = 5.20;
        const auto idle = transmon_fa(solve_EJ_for_f01(5.90, EC), EC);
        const double naive = f_partner - idle.alpha;
        const double tuned = f01_for_11_20(EC, f_partner, idle.alpha);
        const auto at = transmon_fa(solve_EJ_for_f01(tuned, EC), EC);
        std::printf("  11-20: idle alpha %.1f MHz, naive f %.4f GHz, tuned f %.4f GHz, alpha there %.1f MHz\n",
                    idle.alpha * 1e3, naive, tuned, at.alpha * 1e3);
        CHECK(std::abs((at.f01 + at.alpha) - f_partner) < 1e-6, "tuned 11-20 resonance");
        CHECK(tuned > naive + 1e-3, "idle alpha places the interaction too low");
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
    // the flux term is opt-in: reported for every design, counted in e_2q only with cz_flux_noise (alpha1() and its JSON)
    {
        auto on = generate(s); on.cz_flux_noise = true;
        const auto po = analyze(on).pairs[0];
        std::printf("  flux term opt-in: e_2q %.4e without, %.4e with (difference %.3e)\n", p.err_2q, po.err_2q, po.err_2q - p.err_2q);
        CHECK(!generate(s).cz_flux_noise && std::abs((po.err_2q - p.err_2q) - p.cz_flux_err) < 1e-12, "flux term counted only when opted in");
        CHECK(!gen0_surface17(3).cz_flux_noise && !alpha1_v11(3).cz_flux_noise && !beta1(5).cz_flux_noise && alpha1(3).cz_flux_noise,
              "presets: flux term on for Alpha-1 only");
        const bool old_on = design_from_json(Json::parse(R"({"tunable_coupler": {"model_coupler": true}})")).cz_flux_noise;
        const bool old_off = design_from_json(Json::parse(R"({"tunable_coupler": {"g_max_MHz": 40}})")).cz_flux_noise;
        CHECK(design_from_json(design_to_json(alpha1(3))).cz_flux_noise && old_on && !old_off,
              "design JSON: flux term round-trips, and files without the key follow model_coupler");
    }

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
    // Alpha-1 v1.1 and Beta-1 keep the 4 MHz direct coupling their v0.3 outputs were generated with (alpha1() is 1 MHz)
    for (const auto& D : {alpha1_v11(3), beta1(5)}) {
        const auto R = analyze(D);
        double glo = 1e99, ghi = -1e99; int n = 0;
        for (const auto& q : R.pairs) if (q.designed) { ++n; glo = std::min(glo, q.g_MHz); ghi = std::max(ghi, q.g_MHz); }
        std::printf("  %s: %d links, direct g %.3f-%.3f MHz\n", D.name.c_str(), n, glo, ghi);
        CHECK(n > 0 && std::abs(glo - 4.0) < 0.05 && std::abs(ghi - 4.0) < 0.05, "%s direct coupling is 4 MHz", D.name.c_str());
    }
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
        int wrong = 0; DecodeScratch ws;
        for (const auto& f : faults) wrong += dec.decode(f.dets, ws) != int(f.obs);
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
        int wrong1 = 0, wrong2 = 0; DecodeScratch ws;
        for (const auto& f : faults) wrong1 += dec.decode(f.dets, ws) != int(f.obs);
        for (int k = 0; k < 20000; ++k) {
            const auto& f = faults[pick(rng)]; const auto& g = faults[pick(rng)];
            std::vector<int> dets;
            std::set_symmetric_difference(f.dets.begin(), f.dets.end(), g.dets.begin(), g.dets.end(), std::back_inserter(dets));
            wrong2 += dec.decode(dets, ws) != int(f.obs ^ g.obs);
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
    const auto [P, R] = optimize_cz(S, S.c.flux_for_f01(ce.f_on_GHz), ce.t_rise_ns, CZModel::dt_ns, CZModel::max_eval, CZModel::leak_weight, CZModel::leak_limit, CZModel::coh_limit);
    const auto Rf = simulate_cz(S, P, 0.05);
    std::printf("  calibrated CZ: %.1f ns, 1-F %.2e, leakage %.2e, phase error %+.1e rad; dt 0.1 -> 0.05 ns changes F by %.1e\n", R.duration_ns,
                1 - R.fidelity, R.leakage, R.cphase - PI, std::abs(Rf.fidelity - R.fidelity));
    CHECK(1 - R.fidelity < 3e-4, "CZ coherent error %.2e", 1 - R.fidelity);
    CHECK(R.leakage < 2e-4, "CZ leakage %.2e", R.leakage);
    CHECK(std::abs(R.cphase - PI) < 5e-3, "CZ conditional phase");
    CHECK(R.duration_ns > 30 && R.duration_ns < 60, "CZ duration %.1f", R.duration_ns);
    CHECK(std::abs(Rf.fidelity - R.fidelity) < 2e-6, "time-step convergence");
    {   // leakage fallback on the two pair classes with the most leakage under the plain 1 - F calibration:
        //   5.81/5.29 GHz sits on a leaky short-flat-top branch (2.3e-4 > the 2e-4 limit): the repeat must move it to the
        //     low-leakage branch, with 1 - F no worse and within 5e-4;
        //   5.87/5.29 GHz leaks 1.2e-4, within the limit: the fallback must leave its plain calibration untouched.
        for (double fl : {5.81, 5.87}) {
            const auto Sw = design_pair_system(fl, -0.22, 5.29, -0.22, ce);
            const double xw = Sw.c.flux_for_f01(ce.f_on_GHz);
            const auto plain = optimize_cz(Sw, xw, ce.t_rise_ns, CZModel::dt_ns, CZModel::max_eval);
            const auto cal = optimize_cz(Sw, xw, ce.t_rise_ns, CZModel::dt_ns, CZModel::max_eval, CZModel::leak_weight, CZModel::leak_limit, CZModel::coh_limit);
            std::printf("  pair %.2f/5.29 GHz: plain calibration 1-F %.2e, leakage %.2e (%.1f ns); with the leakage fallback 1-F %.2e, leakage %.2e (%.1f ns)\n",
                        fl, 1 - plain.second.fidelity, plain.second.leakage, plain.second.duration_ns, 1 - cal.second.fidelity, cal.second.leakage, cal.second.duration_ns);
            if (fl < 5.84) {
                CHECK(plain.second.leakage > CZModel::leak_limit, "5.81/5.29 plain calibration no longer above the limit (%.2e): revisit the fallback", plain.second.leakage);
                CHECK(cal.second.leakage <= CZModel::leak_limit && 1 - cal.second.fidelity <= 1 - plain.second.fidelity && 1 - cal.second.fidelity <= CZModel::coh_limit,
                      "5.81/5.29 fallback: leakage %.2e, 1-F %.2e (plain %.2e)", cal.second.leakage, 1 - cal.second.fidelity, 1 - plain.second.fidelity);
            } else {
                CHECK(plain.second.leakage <= CZModel::leak_limit, "5.87/5.29 plain leakage %.2e", plain.second.leakage);
                CHECK(cal.second.fidelity == plain.second.fidelity && cal.first.t_flat == plain.first.t_flat, "5.87/5.29 within the limit must keep its plain calibration");
            }
        }
    }
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
        std::printf("  AF-17 mask DRC: %s\n", drc_summary(L).c_str());
        CHECK(L.drc_count("width") == 0 && L.drc_count("space") == 0, "AF-17 mask DRC: %d width + %d space violations",
              L.drc_count("width"), L.drc_count("space"));
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
        const std::size_t expect = L.etch.size() + L.jj.size() + L.bumps.size() + L.ground_bumps.size() + 3 * L.straps.size() + 1;
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

// ---------------------------------------------------------------- json_depth: nesting limit so a hostile file cannot overflow the stack (q11)
static void test_json_depth() {
    std::puts("[json_depth] a 200-level nest parses; a 10,000-level nest throws instead of crashing");
    auto nest = [](int n) { return std::string(std::size_t(n), '[') + std::string(std::size_t(n), ']'); };
    bool deep_ok = false, over_threw = false;
    try {
        const Json j = Json::parse(nest(200));
        // 200 nested one-element arrays end in a single empty array at the centre
        const Json* p = &j;
        for (int k = 0; k < 199; ++k) p = &(*p)[0];
        deep_ok = p->is_arr() && p->size() == 0;
    } catch (const std::exception& e) { std::printf("  200-level nest threw: %s\n", e.what()); }
    CHECK(deep_ok, "200-level nest must parse");
    try {
        Json::parse(nest(10000));
    } catch (const std::runtime_error& e) {
        over_threw = std::string(e.what()).find("nesting too deep") != std::string::npos;
        if (!over_threw) std::printf("  10000-level nest threw a different error: %s\n", e.what());
    }
    CHECK(over_threw, "10000-level nest must throw json: nesting too deep (not crash)");
}

// ---------------------------------------------------------------- gds_write: ofstream errors must be reported (q12)
static void test_gds_write() {
    std::puts("[gds_write] gds::write throws when the output path cannot be created");
    gds::Structure c; c.name = "TOP";
    gds::Polygon p; p.pts = {{0, 0}, {1, 0}, {0, 1}};
    c.polys = {p};
    bool threw = false; std::string msg;
    try { gds::write("/nonexistent-dir-qpuforge-q12/x.gds", "LIB", {c}); }
    catch (const std::runtime_error& e) { threw = true; msg = e.what(); }
    CHECK(threw && msg.find("gds:") == 0, "missing directory: threw=%d msg='%s'", threw, msg.c_str());
    CHECK(!std::filesystem::exists("/nonexistent-dir-qpuforge-q12/x.gds"), "no file should have been created");
    // a successful write still works
    const std::string okp = (std::filesystem::temp_directory_path() / "qpuforge-q12-ok.gds").string();
    gds::write(okp, "LIB", {c});
    CHECK(std::filesystem::file_size(okp) > 0, "successful write left a file");
    std::filesystem::remove(okp);
}

// ---------------------------------------------------------------- json: the parser/writer behind every design and cache file
static std::string json_text(const Json& j, int indent = 2) { std::ostringstream o; j.dump(o, indent); return o.str(); }
static void test_json() {
    std::puts("[json] parser/writer round trips, accessors, error paths, and the design-file round trip for every preset");
    // scalars and containers
    const Json v = Json::parse("{\"n\": -1.5e2, \"i\": 42, \"b\": true, \"f\": false, \"z\": null, \"s\": \"a\\\"b\\\\c\\nd\\te\", \"a\": [1, [2, 3], {\"k\": []}], \"o\": {}}");
    CHECK(v.is_obj() && v.size() == 8 && v["n"].num() == -150.0 && v["i"].num() == 42 && v["b"].boolean() && !v["f"].boolean() && v["z"].is_null(), "scalar members");
    CHECK(v["s"].str() == std::string("a\"b\\c\nd\te"), "string escapes decode to '%s'", v["s"].str().c_str());
    CHECK(v["a"].is_arr() && v["a"].size() == 3 && v["a"][1].size() == 2 && v["a"][1][1].num() == 3 && v["a"][2]["k"].is_arr() && v["a"][2]["k"].size() == 0 && v["o"].is_obj() && v["o"].size() == 0, "nested containers");
    // lookups: has / missing keys give null, get() defaults, typed accessors throw on the wrong type
    CHECK(v.has("n") && !v.has("nope") && v["nope"].is_null() && v["nope"]["deeper"].is_null(), "has / missing key");
    CHECK(v.get("i", 7.0) == 42 && v.get("nope", 7.0) == 7.0 && v.get("s", 7.0) == 7.0, "get(double) falls back on missing or mistyped");
    CHECK(v.get("b", false) && !v.get("nope", false) && v.get("nope", true) && v.get("i", true), "get(bool) falls back on missing or mistyped");
    CHECK(v.get("s", "dflt") == std::string("a\"b\\c\nd\te") && v.get("nope", "dflt") == "dflt" && v.get("i", "dflt") == "dflt", "get(string) falls back on missing or mistyped");
    bool t1 = false, t2 = false, t3 = false, t4 = false;
    try { v["s"].num(); } catch (const std::runtime_error&) { t1 = true; }
    try { v["i"].str(); } catch (const std::runtime_error&) { t2 = true; }
    try { v["i"].boolean(); } catch (const std::runtime_error&) { t3 = true; }
    try { v["a"][9]; } catch (const std::out_of_range&) { t4 = true; }
    CHECK(t1 && t2 && t3 && t4, "wrong-type accessors and out-of-range index must throw (%d %d %d %d)", t1, t2, t3, t4);
    // objects keep insertion order; set() replaces in place
    Json o = Json::object(); o.set("z", 1); o.set("a", 2); o.set("m", 3); o.set("a", 20);
    CHECK(o.size() == 3 && o.items()[0].first == "z" && o.items()[1].first == "a" && o.items()[1].second.num() == 20 && o.items()[2].first == "m", "insertion order / overwrite");
    // comments and whitespace
    CHECK(Json::parse("  // header\n [1, // one\n 2]\n").size() == 2 && Json::parse("\t\n{ }\n").is_obj(), "// comments and whitespace");
    // dump -> parse is the identity (structure and values), at several indents, and dump is a fixed point
    Json w = Json::object();
    w.set("name", "q\"uote\\back\nnl"); w.set("pi", 3.14159265358979); w.set("tiny", 1.25e-17); w.set("big", 6.02e23); w.set("neg", -0.5); w.set("yes", true); w.set("none", Json());
    Json arr = Json::array(); for (int k = 0; k < 5; ++k) arr.push(k * 0.25); arr.push(Json::array()); arr.push(Json::object());
    w.set("list", arr);
    Json inner = Json::object(); inner.set("x", 1); inner.set("y", Json::array()); w.set("inner", inner);
    for (int indent : {0, 1, 2, 4}) {
        const std::string t = json_text(w, indent);
        const Json back = Json::parse(t);
        CHECK(json_text(back, indent) == t, "dump(parse(dump(x))) differs at indent %d", indent);
        CHECK(back["name"].str() == w["name"].str() && std::abs(back["pi"].num() - 3.14159265358979) < 1e-11 && std::abs(back["tiny"].num() / 1.25e-17 - 1) < 1e-11 &&
              std::abs(back["big"].num() / 6.02e23 - 1) < 1e-11 && back["yes"].boolean() && back["none"].is_null() && back["list"].size() == 7 && back["list"][4].num() == 1.0,
              "values after round trip at indent %d", indent);
    }
    // JSON has no inf/nan: the writer emits null so a file always parses
    Json nf = Json::array(); nf.push(std::nan("")); nf.push(HUGE_VAL); nf.push(-HUGE_VAL); nf.push(1.0);
    const std::string nft = json_text(nf, 0);
    const Json nb = Json::parse(nft);
    CHECK(nft == "[null, null, null, 1]" && nb.size() == 4 && nb[0].is_null() && nb[3].num() == 1.0, "non-finite numbers are written as null: %s", nft.c_str());
    // 12 significant digits survive (the precision the design and cache files rely on)
    for (double x : {0.0352055795967, 123456.789012, 1.0 / 3.0, -2.718281828459, 5.81, 1e-9}) {
        const double y = Json::parse(json_text(Json(x), 0)).num();
        CHECK(std::abs(y - x) <= 1e-11 * std::abs(x), "number %.15g became %.15g", x, y);
    }
    // malformed input must throw, never crash or accept
    const char* bad[] = {"", "   ", "{", "[", "[1", "[1,", "[1,]", "{\"a\"}", "{\"a\":}", "{\"a\" 1}", "{a: 1}", "{\"a\": 1,}", "{\"a\": 1 \"b\": 2}", "[1 2]", "\"abc",
                         "tru", "nul", "[truex]", "1 2", "{} {}", "[1]]", "}", "]", ":", "--1", "[,1]"};
    int accepted = 0;
    for (const char* b : bad) if (!json_throws(b)) { ++accepted; std::printf("  accepted malformed JSON: '%s'\n", b); }
    CHECK(accepted == 0, "%d malformed documents were accepted", accepted);
    // a plain top-level scalar is fine
    CHECK(Json::parse("12").num() == 12 && Json::parse("\"x\"").str() == "x" && Json::parse("true").boolean() && Json::parse("null").is_null(), "top-level scalars");
    // duplicate keys: last one wins, no duplicate member is created
    const Json dup = Json::parse("{\"a\": 1, \"a\": 2}");
    CHECK(dup.size() == 1 && dup["a"].num() == 2, "duplicate key handling");

    // design files: JSON round trip must be lossless for every preset (field written but not read back would show up here)
    struct P { const char* name; ChipDesign D; };
    const P presets[] = {{"gen0_surface17", gen0_surface17(3)}, {"alpha1", alpha1(3)}, {"alpha1_v11", alpha1_v11(3)}, {"beta1", beta1(5)}};
    for (const auto& p : presets) {
        const std::string a = json_text(design_to_json(p.D));
        const ChipDesign D2 = design_from_json(Json::parse(a));
        const std::string b = json_text(design_to_json(D2));
        CHECK(a == b, "%s: design JSON changes after a read/write cycle", p.name);
        CHECK(D2.qubits.size() == p.D.qubits.size() && D2.couplings.size() == p.D.couplings.size() && D2.surface_d == p.D.surface_d && D2.name == p.D.name &&
              D2.leakage.model == p.D.leakage.model && D2.coupler.physical == p.D.coupler.physical && D2.active_reset == p.D.active_reset,
              "%s: structure after round trip (qubits %zu/%zu, couplings %zu/%zu)", p.name, D2.qubits.size(), p.D.qubits.size(), D2.couplings.size(), p.D.couplings.size());
    }
    {   // a design saved to disk and loaded back analyses identically
        const ChipDesign D = gen0_surface17(3);
        const std::string path = temp_path("design_rt.json");
        save_design(D, path);
        const ChipDesign L = load_design(path);
        std::remove(path.c_str());
        const auto R0 = analyze(D), R1 = analyze(L);
        double df = 0; for (std::size_t i = 0; i < R0.q.size(); ++i) df = std::max(df, std::abs(R0.q[i].f01 - R1.q[i].f01));
        CHECK(R0.q.size() == R1.q.size() && df < 1e-9 && R0.n_collisions == R1.n_collisions && std::abs(R0.mean_err_2q - R1.mean_err_2q) <= 1e-8 * R0.mean_err_2q,
              "analysis after save/load: max f01 diff %.2e GHz, collisions %d/%d", df, R0.n_collisions, R1.n_collisions);
        bool threw = false;
        try { load_design(temp_path("does_not_exist.json")); } catch (const std::runtime_error& e) { threw = std::string(e.what()).find("cannot open") != std::string::npos; }
        CHECK(threw, "load_design of a missing file must throw 'cannot open'");
    }
    // sparse design file: defaults apply for everything omitted
    const ChipDesign S = design_from_json(Json::parse("{\"qubits\": [{\"C_fF\": 90}, {\"Rn_kOhm\": 12}], \"couplings\": [{\"a\": 0, \"b\": 1}]}"));
    CHECK(S.qubits.size() == 2 && S.qubits[0].id == 0 && S.qubits[1].id == 1 && S.qubits[0].C_fF == 90 && S.qubits[1].Rn_kOhm == 12 && S.couplings.size() == 1 &&
          S.couplings[0].C_fF == 1.0 && S.name == "chip" && S.surface_d == 0, "sparse design uses documented defaults");
}

// ---------------------------------------------------------------- gds: byte-level checks of the stream writer
struct GdsRec { int type{0}, dtype{0}; std::vector<unsigned char> data; };
static bool gds_read(const std::string& path, std::vector<GdsRec>& out) {
    std::ifstream f(path, std::ios::binary);
    std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    out.clear();
    std::size_t i = 0;
    while (i + 4 <= b.size()) {
        const std::size_t len = (std::size_t(b[i]) << 8) | b[i + 1];
        if (len < 4 || len % 2 || i + len > b.size()) return false;
        GdsRec r; r.type = b[i + 2]; r.dtype = b[i + 3]; r.data.assign(b.begin() + long(i) + 4, b.begin() + long(i + len));
        out.push_back(std::move(r));
        i += len;
    }
    return i == b.size();
}
static double gds_real8(const unsigned char* p) {
    const int e = (p[0] & 0x7f) - 64;
    std::uint64_t m = 0; for (int k = 1; k < 8; ++k) m = (m << 8) | p[k];
    return (p[0] & 0x80 ? -1.0 : 1.0) * double(m) / std::pow(2.0, 56) * std::pow(16.0, e);
}
static int gds_i32(const std::vector<unsigned char>& d, std::size_t at) { return int((std::uint32_t(d[at]) << 24) | (std::uint32_t(d[at + 1]) << 16) | (std::uint32_t(d[at + 2]) << 8) | d[at + 3]); }
static int gds_u16(const std::vector<unsigned char>& d, std::size_t at) { return (d[at] << 8) | d[at + 1]; }

static void test_gds() {
    std::puts("[gds] GDSII stream writer: record framing, UNITS encoding, XY coordinates, guards for oversize records and coordinates, reproducible dates");
    gds::Structure c; c.name = "TOP";
    gds::Polygon tri; tri.layer = 1; tri.datatype = 0; tri.pts = {{0, 0}, {10, 0}, {0.5, 7.25}};
    gds::Polygon neg; neg.layer = 2; neg.datatype = 5; neg.pts = {{-1.5, -2.0}, {1.5, -2.0}, {1.5, 2.0}, {-1.5, 2.0}};
    gds::Polygon degenerate; degenerate.layer = 3; degenerate.pts = {{0, 0}, {1, 1}};   // fewer than 3 points: skipped
    c.polys = {tri, degenerate, neg};
    gds::Text t; t.layer = 20; t.x = 2.5; t.y = -3.25; t.s = "ODD";   // odd length: padded to even with a NUL
    c.texts = {t};
    gds::Structure c2; c2.name = "B";
    const std::string path = temp_path("rt.gds");
    gds::write(path, "LIBX", {c, c2});
    std::vector<GdsRec> R;
    CHECK(gds_read(path, R), "stream is not a sequence of well-formed records");
    if (R.size() < 10) { std::remove(path.c_str()); CHECK(false, "only %zu records", R.size()); return; }
    CHECK(R.front().type == 0x00 && R.front().dtype == 0x02 && gds_u16(R.front().data, 0) == 600, "HEADER version 600");
    CHECK(R.back().type == 0x04 && R.back().data.empty(), "stream ends with ENDLIB");
    CHECK(R[1].type == 0x01 && R[1].data.size() == 24, "BGNLIB carries two 12-byte timestamps (%zu bytes)", R[1].data.size());
    CHECK(R[2].type == 0x02 && std::string(R[2].data.begin(), R[2].data.end()) == "LIBX", "LIBNAME");
    CHECK(R[3].type == 0x03 && R[3].data.size() == 16, "UNITS is two 8-byte reals");
    if (R[3].data.size() == 16) {
        const double uu = gds_real8(R[3].data.data()), dbu = gds_real8(R[3].data.data() + 8);
        CHECK(std::abs(uu / 1e-3 - 1) < 1e-14 && std::abs(dbu / 1e-9 - 1) < 1e-14, "UNITS decode to %.17g / %.17g", uu, dbu);
        // GDSII spec example bytes: 1e-3 user unit per dbu = 3E 41 89 37 4B C6 A7 F0, 1e-9 m = 39 44 B8 2F A0 9B 5A 54
        const unsigned char e1[8] = {0x3E, 0x41, 0x89, 0x37, 0x4B, 0xC6, 0xA7, 0xF0}, e2[8] = {0x39, 0x44, 0xB8, 0x2F, 0xA0, 0x9B, 0x5A, 0x54};
        CHECK(std::memcmp(R[3].data.data(), e1, 8) == 0 && std::memcmp(R[3].data.data() + 8, e2, 8) == 0, "UNITS bytes differ from the GDSII reference encoding");
    }
    int boundaries = 0, texts = 0, strs = 0, endstr = 0;
    std::vector<std::vector<int>> xys; std::vector<std::pair<int, int>> layer_dt; std::string text_string; int text_layer = -1;
    int cur_layer = -1, cur_dt = -1; bool in_text = false;
    std::vector<std::string> names;
    for (const auto& r : R) {
        switch (r.type) {
        case 0x05: ++strs; break;
        case 0x06: { std::string n(r.data.begin(), r.data.end()); while (!n.empty() && n.back() == '\0') n.pop_back(); names.push_back(n); break; }   // odd-length names are NUL-padded
        case 0x07: ++endstr; break;
        case 0x08: ++boundaries; in_text = false; break;
        case 0x0C: ++texts; in_text = true; break;
        case 0x0D: cur_layer = gds_u16(r.data, 0); if (in_text) text_layer = cur_layer; break;
        case 0x0E: cur_dt = gds_u16(r.data, 0); layer_dt.emplace_back(cur_layer, cur_dt); break;
        case 0x10: { std::vector<int> v; for (std::size_t k = 0; k + 4 <= r.data.size(); k += 4) v.push_back(gds_i32(r.data, k)); xys.push_back(v); break; }
        case 0x19: text_string.assign(r.data.begin(), r.data.end()); break;
        default: break;
        }
    }
    CHECK(strs == 2 && endstr == 2 && names.size() == 2 && names[0] == "TOP" && names[1] == "B", "two structures with names (%d, %d)", strs, endstr);
    CHECK(boundaries == 2 && texts == 1, "degenerate polygon skipped: %d boundaries, %d texts", boundaries, texts);
    CHECK(layer_dt.size() == 2 && layer_dt[0] == std::make_pair(1, 0) && layer_dt[1] == std::make_pair(2, 5), "layer/datatype pairs");
    CHECK(xys.size() == 3, "XY records: %zu", xys.size());
    if (xys.size() == 3) {
        const std::vector<int> e0 = {0, 0, 10000, 0, 500, 7250, 0, 0};                                         // closed ring, nm
        const std::vector<int> e1 = {-1500, -2000, 1500, -2000, 1500, 2000, -1500, 2000, -1500, -2000};
        const std::vector<int> e2 = {2500, -3250};
        CHECK(xys[0] == e0 && xys[1] == e1 && xys[2] == e2, "XY coordinates in database units (um * 1000, ring closed, signed)");
    }
    CHECK(text_layer == 20 && text_string.size() == 4 && text_string.substr(0, 3) == "ODD" && text_string[3] == '\0', "TEXT layer and NUL-padded string (len %zu)", text_string.size());

    // guards: the failing write must throw and leave no file behind
    auto write_throws = [&](const gds::Structure& s, const char* what) {
        const std::string p = temp_path("bad.gds");
        std::remove(p.c_str());
        bool threw = false; std::string msg;
        try { gds::write(p, "L", {s}); } catch (const std::runtime_error& e) { threw = true; msg = e.what(); }
        const bool left = std::filesystem::exists(p);
        std::remove(p.c_str());
        CHECK(threw && !left && msg.find("gds:") == 0, "%s: threw=%d file left=%d msg='%s'", what, threw, left, msg.c_str());
    };
    auto ring = [](int n) { gds::Structure s; s.name = "R"; gds::Polygon p; for (int k = 0; k < n; ++k) p.pts.push_back({double(k % 100), double(k / 100)}); s.polys = {p}; return s; };
    {   // exactly at the limit: 8190 vertices + closing point = 8191 points = 65528 data bytes
        const std::string p = temp_path("max.gds");
        bool threw = false; try { gds::write(p, "L", {ring(8190)}); } catch (const std::runtime_error&) { threw = true; }
        std::vector<GdsRec> Rm; const bool parsed = !threw && gds_read(p, Rm);
        std::size_t xy = 0; for (const auto& r : Rm) if (r.type == 0x10) xy = r.data.size();
        std::remove(p.c_str());
        CHECK(!threw && parsed && xy == 65528, "8190-vertex polygon fits one record (threw=%d, XY bytes %zu)", threw, xy);
    }
    write_throws(ring(8191), "8191-vertex polygon exceeds the 65530-byte record limit");
    { gds::Structure s = ring(4); s.polys[0].pts[2].first = 3.0e6; write_throws(s, "coordinate 3e6 um = 3e9 dbu exceeds int32"); }
    { gds::Structure s = ring(4); s.polys[0].pts[1].second = -3.0e6; write_throws(s, "coordinate -3e6 um exceeds int32"); }
    { gds::Structure s = ring(4); s.polys[0].pts[0].first = std::nan(""); write_throws(s, "NaN coordinate"); }
    { gds::Structure s = ring(4); s.polys[0].pts[0].second = HUGE_VAL; write_throws(s, "infinite coordinate"); }
    { gds::Structure s; s.name = "T"; gds::Text tx; tx.x = 5e6; tx.s = "x"; s.texts = {tx}; write_throws(s, "text position out of range"); }
    {   // the largest representable coordinate still works
        gds::Structure s = ring(4); s.polys[0].pts[0].first = 2147483.0;
        const std::string p = temp_path("edge.gds"); bool threw = false;
        try { gds::write(p, "L", {s}); } catch (const std::runtime_error&) { threw = true; }
        std::vector<GdsRec> Re; const bool ok = !threw && gds_read(p, Re); std::remove(p.c_str());
        CHECK(ok, "coordinate 2147483 um (2147483000 dbu) is accepted");
    }
    // dates: SOURCE_DATE_EPOCH makes the stream reproducible byte for byte
#if !defined(_WIN32)
    {
        const char* old = std::getenv("SOURCE_DATE_EPOCH"); const std::string saved = old ? old : "";
        setenv("SOURCE_DATE_EPOCH", "1700000000", 1);   // 2023-11-14 22:13:20 UTC
        const std::string p1 = temp_path("d1.gds"), p2 = temp_path("d2.gds");
        gds::write(p1, "L", {c}); gds::write(p2, "L", {c});
        std::ifstream f1(p1, std::ios::binary), f2(p2, std::ios::binary);
        const std::string b1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>()), b2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
        std::vector<GdsRec> Rd; gds_read(p1, Rd);
        std::remove(p1.c_str()); std::remove(p2.c_str());
        CHECK(!b1.empty() && b1 == b2, "two writes with a fixed SOURCE_DATE_EPOCH differ");
        if (Rd.size() > 1 && Rd[1].data.size() == 24) {
            const int y = gds_u16(Rd[1].data, 0), mo = gds_u16(Rd[1].data, 2), d = gds_u16(Rd[1].data, 4), h = gds_u16(Rd[1].data, 6), mi = gds_u16(Rd[1].data, 8), s = gds_u16(Rd[1].data, 10);
            CHECK(y == 2023 && mo == 11 && d == 14 && h == 22 && mi == 13 && s == 20, "BGNLIB date %04d-%02d-%02d %02d:%02d:%02d, expected 2023-11-14 22:13:20", y, mo, d, h, mi, s);
        } else CHECK(false, "BGNLIB record missing");
        unsetenv("SOURCE_DATE_EPOCH");
        const std::string p3 = temp_path("d3.gds"); gds::write(p3, "L", {c});
        std::vector<GdsRec> Rn; gds_read(p3, Rn); std::remove(p3.c_str());
        const int yy = Rn.size() > 1 && Rn[1].data.size() >= 2 ? gds_u16(Rn[1].data, 0) : 0;
        CHECK(yy >= 2024 && yy < 2200, "without SOURCE_DATE_EPOCH the stamp is the current year (got %d)", yy);
        if (old) setenv("SOURCE_DATE_EPOCH", saved.c_str(), 1);
    }
#endif
    std::remove(path.c_str());
}

// ---------------------------------------------------------------- v0.4: exact matching for any number of defects
static void test_matching() {
    std::puts("[matching] blossom matcher vs brute force; decoder exact above the old 22-defect limit");
    // (a) the matcher alone: maximum-weight matching on random graphs (complete and sparse, even and odd sizes)
    //     against a subset dynamic programme that may leave vertices unmatched
    {
        std::mt19937_64 rng(42);
        BlossomMatcher bm;
        int bad = 0, cases = 0;
        for (int rep = 0; rep < 600; ++rep) {
            const int n = 2 + int(rng() % 11);                       // 2 .. 12 vertices
            const bool sparse = rep % 3 == 0;
            std::vector<BlossomMatcher::W> w(std::size_t(n) * std::size_t(n), 0);
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
                BlossomMatcher::W v = 1 + BlossomMatcher::W(rng() % 1000);
                if (sparse && rng() % 3 == 0) v = 0;                 // no edge
                w[std::size_t(i) * std::size_t(n) + std::size_t(j)] = w[std::size_t(j) * std::size_t(n) + std::size_t(i)] = v;
            }
            std::vector<BlossomMatcher::W> best(std::size_t(1) << n, 0);
            for (std::size_t mask = 1; mask < best.size(); ++mask) {
                int i = 0; while (!((mask >> i) & 1)) ++i;
                const std::size_t rest = mask & ~(std::size_t(1) << i);
                BlossomMatcher::W b = best[rest];                    // i stays unmatched
                for (int j = i + 1; j < n; ++j) if ((rest >> j) & 1) {
                    const BlossomMatcher::W e = w[std::size_t(i) * std::size_t(n) + std::size_t(j)];
                    if (e > 0) b = std::max(b, e + best[rest & ~(std::size_t(1) << j)]);
                }
                best[mask] = b;
            }
            const auto& mate = bm.solve(n, w);
            BlossomMatcher::W tot = 0; bool consistent = true;
            for (int i = 0; i < n; ++i) {
                const int j = mate[std::size_t(i)];
                if (j < 0) continue;
                if (mate[std::size_t(j)] != i || w[std::size_t(i) * std::size_t(n) + std::size_t(j)] <= 0) consistent = false;
                if (i < j) tot += w[std::size_t(i) * std::size_t(n) + std::size_t(j)];
            }
            ++cases;
            if (!consistent || tot != best.back()) ++bad;
        }
        std::printf("  matcher: %d random graphs (2-12 vertices, a third of them sparse), %d differ from brute force\n", cases, bad);
        CHECK(bad == 0, "blossom matcher is not optimal on %d of %d graphs", bad, cases);
    }
    // (b) the decoder: blossom against the exact DP on random defect sets, then large sets against greedy
    const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
    const auto M = prepare_memory(D, R, nm, 12);
    {
        std::mt19937_64 rng(7);
        const int N = int(M.exp.detectors.size());
        DecodeScratch a, b;
        int bad_w = 0, parity_diff = 0, n = 0;
        double worst = 0;
        for (int rep = 0; rep < 2500; ++rep) {
            const int k = 1 + rep % 18;
            std::vector<int> def;
            if (rep % 2) {                                           // clustered: sum of random error-model edges
                std::vector<char> on(std::size_t(N), 0);
                for (int e = 0; e < k; ++e) { const auto& ed = M.dem.edges[rng() % M.dem.edges.size()]; on[std::size_t(ed.a)] ^= 1; if (ed.b >= 0) on[std::size_t(ed.b)] ^= 1; }
                for (int i = 0; i < N; ++i) if (on[std::size_t(i)]) def.push_back(i);
                if (def.empty() || def.size() > 20) continue;
            } else {                                                 // uniformly random detectors
                std::vector<int> all(static_cast<std::size_t>(N)); std::iota(all.begin(), all.end(), 0); std::shuffle(all.begin(), all.end(), rng);
                def.assign(all.begin(), all.begin() + k); std::sort(def.begin(), def.end());
            }
            const int pa = M.decoder->decode(def, a, DecodeMethod::DP), pb = M.decoder->decode(def, b, DecodeMethod::Blossom);
            ++n; worst = std::max(worst, std::abs(a.weight - b.weight));
            if (std::abs(a.weight - b.weight) > 1e-5) ++bad_w;
            if (pa != pb) ++parity_diff;
        }
        std::printf("  decoder: %d defect sets (1-20 defects): %d weight mismatches (worst %.1e), %d parity differences\n", n, bad_w, worst, parity_diff);
        CHECK(bad_w == 0, "blossom and DP matchings differ in weight");
        CHECK(parity_diff == 0, "blossom and DP predict different logical flips");
        int heavier = 0, big = 0;
        for (int k : {30, 60, 100, 140}) for (int rep = 0; rep < 15; ++rep) {
            std::vector<int> all(static_cast<std::size_t>(N)); std::iota(all.begin(), all.end(), 0); std::shuffle(all.begin(), all.end(), rng);
            std::vector<int> def(all.begin(), all.begin() + std::min(k, N)); std::sort(def.begin(), def.end());
            M.decoder->decode(def, b, DecodeMethod::Blossom); M.decoder->decode(def, a, DecodeMethod::Greedy);
            ++big; if (b.weight > a.weight + 1e-6) ++heavier;
        }
        CHECK(heavier == 0, "blossom heavier than greedy on %d of %d large defect sets", heavier, big);
    }
    // (c) a long experiment, same syndromes: exact matching must not lose to the old fallback, and must gain where the
    //     fallback was in use
    {
        const auto ML = prepare_memory(D, R, nm, 30);
        FrameRunOptions legacy; legacy.method = DecodeMethod::Legacy;
        const auto ex = run_memory_frame(ML, 1200, 5), lg = run_memory_frame(ML, 1200, 5, legacy);
        std::printf("  d=5, 30 rounds, 1200 shots: %.1f defects per shot, %.0f%% above 22; logical failures exact %ld, greedy fallback %ld; decode %.0f vs %.0f us per shot\n",
                    ex.mean_defects, 100 * ex.frac_over_22, ex.logical_errors, lg.logical_errors, ex.decode_us_per_shot, lg.decode_us_per_shot);
        CHECK(ex.frac_over_22 > 0.05, "test needs shots above 22 defects");
        CHECK(ex.logical_errors < lg.logical_errors, "exact matching should beat the greedy fallback on the same syndromes");
        CHECK(int(ex.detection_by_round.size()) == 31 && ex.max_defects > 22, "per-round detection and defect statistics");
    }
}

// ---------------------------------------------------------------- v0.4: fast fault enumeration
static void test_fault_enumeration() {
    std::puts("[faults] sparse and linear-time single-fault enumeration against the reference, fault by fault");
    long total = 0, differ = 0, differ_direct = 0; double t_ref = 0, t_direct = 0, t_fast = 0;
    auto same = [](const FaultSignature& a, const FaultSignature& b) {
        return a.dets == b.dets && a.obs == b.obs && a.p == b.p && a.moment == b.moment && a.what == b.what;
    };
    auto compare_frame = [&](const FrameProgram& F, const MemoryExperiment& E) {
        const auto t0 = std::chrono::steady_clock::now();
        const auto a = enumerate_faults_reference(F, E);
        const auto t1 = std::chrono::steady_clock::now();
        const auto b = enumerate_faults_direct(F, E);
        const auto t2 = std::chrono::steady_clock::now();
        const auto c = enumerate_faults(F, E);
        const auto t3 = std::chrono::steady_clock::now();
        t_ref += std::chrono::duration<double>(t1 - t0).count(); t_direct += std::chrono::duration<double>(t2 - t1).count();
        t_fast += std::chrono::duration<double>(t3 - t2).count();
        if (a.size() != b.size() || a.size() != c.size()) { differ += long(std::max(a.size(), c.size())); return; }
        for (std::size_t i = 0; i < a.size(); ++i) {
            ++total;
            if (!same(a[i], c[i])) ++differ;
            if (!same(a[i], b[i])) ++differ_direct;
        }
    };
    auto compare = [&](const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, int rounds, char basis) {
        const auto M = prepare_memory(D, R, nm, rounds, basis, false);
        compare_frame(M.frame, M.exp);
    };
    {
        const auto D = surface_chip(3); const auto R = analyze(D); const auto nm = noise_model(D, R);
        compare(D, R, nm, 5, 'Z'); compare(D, R, nm, 5, 'X'); compare(D, R, nm, 1, 'Z'); compare(D, R, nm, 2, 'X');
        auto lk = nm; lk.leakage = true; lk.dqlr = true; lk.dqlr_err = 5e-4; lk.dqlr_eff = 0.98;      // DQLR faults
        const auto ML = prepare_memory(D, R, lk, 3, 'Z', false);                                      // (also builds the DQLR partner map)
        compare_frame(ML.frame, ML.exp);
    }
    {
        const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
        compare(D, R, nm, 8, 'Z'); compare(D, R, nm, 4, 'X');
    }
    // Unequal CZ durations, as per-pair calibration gives them: the program is cut into several times more time
    // slices, most of them idle for any given qubit. The fast enumeration skips those; the output must not change.
    std::size_t slices_equal = 0, slices_unequal = 0;
    {
        const auto D = surface_chip(3); auto R = analyze(D);
        slices_equal = prepare_memory(D, R, noise_model(D, R), 6, 'Z', false).frame.moments.size();
        int k = 0; for (auto& p : R.pairs) if (p.designed) p.t_cz_ns += 1.7 * (k++ % 11);
        const auto nm = noise_model(D, R);
        slices_unequal = prepare_memory(D, R, nm, 6, 'Z', false).frame.moments.size();
        compare(D, R, nm, 6, 'Z'); compare(D, R, nm, 6, 'X');
    }
    std::printf("  %ld faults compared over 10 experiments (d = 3 and 5, both bases, 1 to 8 rounds, with DQLR, equal and unequal CZ durations):\n"
                "  linear-time %ld differ, sparse %ld differ; reference %.2f s, sparse %.2f s, linear-time %.2f s\n",
                total, differ, differ_direct, t_ref, t_direct, t_fast);
    std::printf("  time slices of a 6-round d=3 memory: %zu with equal CZ durations, %zu with unequal ones\n", slices_equal, slices_unequal);
    CHECK(total > 100000 && differ == 0, "linear-time fault enumeration differs from the reference on %ld of %ld faults", differ, total);
    CHECK(differ_direct == 0, "sparse fault enumeration differs from the reference on %ld of %ld faults", differ_direct, total);
    CHECK(2 * slices_unequal > 3 * slices_equal, "the unequal-duration case does not exercise the many-slice regime (%zu vs %zu slices)", slices_unequal, slices_equal);
    // scaling: the work per fault (time slices processed, a deterministic count) must not grow with the number of
    // rounds (d = 5, 10 against 40 rounds); for the sparse version it grows about 4x. Wall-clock time is printed only.
    {
        const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
        double work[2] = {0, 0}, work_direct[2] = {0, 0}, secs[2] = {0, 0}; std::size_t nf[2] = {0, 0};
        const int rr[2] = {10, 40};
        for (int i = 0; i < 2; ++i) {
            const auto M = prepare_memory(D, R, nm, rr[i], 'Z', false);
            long st = 0, sd = 0;
            const auto t0 = std::chrono::steady_clock::now();
            nf[i] = enumerate_faults(M.frame, M.exp, &st).size();
            secs[i] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            enumerate_faults_direct(M.frame, M.exp, &sd);
            work[i] = double(st) / double(nf[i]); work_direct[i] = double(sd) / double(nf[i]);
        }
        std::printf("  time slices processed per fault, d=5: linear-time %.2f at 10 rounds, %.2f at 40 (sparse %.2f, %.2f); %zu and %zu faults, %.3f and %.3f us per fault\n",
                    work[0], work[1], work_direct[0], work_direct[1], nf[0], nf[1], 1e6 * secs[0] / double(nf[0]), 1e6 * secs[1] / double(nf[1]));
        CHECK(work[1] < 1.25 * work[0] && work_direct[1] > 2.5 * work_direct[0],
              "work per fault: linear-time %.2f -> %.2f, sparse %.2f -> %.2f", work[0], work[1], work_direct[0], work_direct[1]);
    }
}

// ---------------------------------------------------------------- union-find decoder
static void test_union_find() {
    std::puts("[union-find] weighted union-find decoder: clears every syndrome, decodes single faults, logical error near the exact decoder");
    long syndromes = 0, uncleared = 0, singles = 0, single_wrong = 0, single_wrong_exact = 0, outweighed = 0;
    auto check_setup = [&](const MemorySetup& M, std::uint64_t seed) {
        const MatchingDecoder& dec = *M.decoder;
        const auto& uf = dec.union_find();
        const int nd = dec.n_detectors();
        DecodeScratch ws;
        // syndromes from the detector error model (every mechanism at 3x its probability) and uniformly random ones
        Xoshiro256 rng(seed);
        std::vector<std::uint8_t> on(std::size_t(nd), 0), got(std::size_t(nd), 0);
        std::vector<int> def;
        for (int trial = 0; trial < 3000; ++trial) {
            std::fill(on.begin(), on.end(), 0);
            if (trial % 3 == 2) { for (int d = 0; d < nd; ++d) on[std::size_t(d)] = rng.uniform() < 0.1; }
            else for (const auto& e : M.dem.edges) if (rng.uniform() < 3 * e.p) { on[std::size_t(e.a)] ^= 1; if (e.b >= 0) on[std::size_t(e.b)] ^= 1; }
            def.clear(); for (int d = 0; d < nd; ++d) if (on[std::size_t(d)]) def.push_back(d);
            dec.decode(def, ws, DecodeMethod::UnionFind);
            std::fill(got.begin(), got.end(), 0);
            for (int k : ws.uf.correction) { const auto& e = uf.edge[std::size_t(k)]; for (int v : {e.u, e.v}) if (v < nd) got[std::size_t(v)] ^= 1; }
            ++syndromes; if (got != on) ++uncleared;
        }
        // every single fault mechanism of the circuit that the graph represents
        // (a mechanism whose detectors are joined by a more likely path with the other logical parity is decoded
        // "wrongly" by any decoder that uses these priors: counted separately, and its weight must exceed that path's)
        for (const auto& f : enumerate_faults(M.frame, M.exp)) {
            if (f.dets.empty() || f.dets.size() > 2) continue;
            const double w = std::log((1 - f.p) / f.p);
            const int b = f.dets.size() == 2 ? f.dets[1] : -1;
            const int ex = dec.decode(f.dets, ws);
            if (ex != int(f.obs)) {
                ++outweighed;
                if (!(dec.dist(f.dets[0], b) < w - 1e-9)) ++single_wrong_exact;
                continue;
            }
            ++singles;
            if (dec.decode(f.dets, ws, DecodeMethod::UnionFind) != int(f.obs)) ++single_wrong;
        }
    };
    // peeling on a hand-made graph: 0 -1- 1 -2- 2 -3- boundary, plus a heavy 0 - boundary edge; edge bits 0/1/0/1
    {
        struct E { int a, b; bool obs; double p; };
        const std::vector<E> g{{0, 1, false, 0.1}, {1, 2, true, 0.1}, {2, -1, false, 0.1}, {0, -1, true, 1e-4}};
        const UnionFindDecoder U(3, g);
        UFScratch us;
        const int p02 = U.decode({0, 2}, us);            // the pair is joined through 1: edges 0-1 and 1-2, parity 1
        std::vector<int> c02 = us.correction; std::sort(c02.begin(), c02.end());
        const int p2 = U.decode({2}, us);                // one defect next to the boundary: edge 2-B, parity 0
        const std::size_t n2 = us.correction.size();
        const int p0 = U.decode({0}, us);                // one defect: 0-1-2-B (weight 3 x 2.2) beats 0-B (9.2): parity 1
        std::printf("  hand-made graph: {0,2} -> parity %d via %zu edges; {2} -> parity %d via %zu edge; {0} -> parity %d via %zu edges\n",
                    p02, c02.size(), p2, n2, p0, us.correction.size());
        CHECK(p02 == 1 && c02.size() == 2 && p2 == 0 && n2 == 1 && p0 == 1 && us.correction.size() == 3, "peeling on the hand-made graph");
    }
    struct Row { int d; char b; int rounds; long shots; MemoryResult ex, uf; };
    std::vector<Row> rows;
    for (int d : {3, 5}) {
        const auto D = surface_chip(d); const auto R = analyze(D); const auto nm = noise_model(D, R);
        for (char b : {'Z', 'X'}) {
            const int rounds = d == 3 ? 6 : 5;
            const auto M = prepare_memory(D, R, nm, rounds, b);
            check_setup(M, 31 + std::uint64_t(d) + (b == 'X' ? 100u : 0u));
            const long shots = d == 3 ? 40000 : 60000;
            FrameRunOptions o; o.method = DecodeMethod::UnionFind;
            rows.push_back({d, b, rounds, shots, run_memory_frame(M, shots, 900 + std::uint64_t(d)), run_memory_frame(M, shots, 900 + std::uint64_t(d), o)});
        }
    }
    std::printf("  %ld syndromes (sampled from the error model and uniformly random): %ld not cleared by the correction\n", syndromes, uncleared);
    std::printf("  %ld single-fault mechanisms: union-find wrong on %ld; %ld more are outweighed by a likelier path of the other parity "
                "(%ld of those not explained by the weights)\n", singles, single_wrong, outweighed, single_wrong_exact);
    double worst = 0;
    for (const auto& r : rows) {
        const double ratio = r.uf.p_logical / std::max(r.ex.p_logical, 1e-12);
        worst = std::max(worst, ratio);
        std::printf("  d=%d memory-%c, %d rounds, %ld shots (same syndromes): exact P_L %.5f (%ld), union-find P_L %.5f (%ld), ratio %.2f; decode %.2f vs %.2f us per shot\n",
                    r.d, r.b, r.rounds, r.shots, r.ex.p_logical, r.ex.logical_errors, r.uf.p_logical, r.uf.logical_errors, ratio,
                    r.ex.decode_us_per_shot, r.uf.decode_us_per_shot);
        CHECK(r.uf.logical_errors >= r.ex.logical_errors * 9 / 10, "union-find beats the exact decoder by more than noise (d=%d %c: %ld vs %ld)", r.d, r.b, r.uf.logical_errors, r.ex.logical_errors);
    }
    CHECK(uncleared == 0, "union-find correction leaves %ld of %ld syndromes uncleared", uncleared, syndromes);
    CHECK(single_wrong == 0 && single_wrong_exact == 0 && outweighed < singles / 1000, "single faults decoded wrongly: union-find %ld, exact %ld (outweighed %ld)",
          single_wrong, single_wrong_exact, outweighed);
    CHECK(worst < 1.5, "union-find logical error more than 1.5x the exact decoder's (worst ratio %.2f)", worst);
}

// ---------------------------------------------------------------- heralded leakage
static NoiseModel stressed_leakage(const ChipDesign& D, const ChipReport& R, double p_cz) {
    auto nm = noise_model(D, R);
    nm.leakage = true; nm.leak_1q = 0; nm.leak_transport = 0.05; nm.reset_removes_leakage = true;
    for (const auto& p : R.pairs) if (p.designed) {
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;   // the measure qubit leaks
        nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, p_cz};
    }
    return nm;
}

static void test_heralds() {
    std::puts("[heralds] heralded leakage: modified-graph matching is exact, herald statistics, logical gain with the leakage-removing reset");
    // 1. decode_scaled against brute force: Floyd-Warshall on the whole modified graph, then the DP over the defects
    {
        const auto D = surface_chip(3); const auto R = analyze(D); const auto nm = noise_model(D, R);
        const auto M = prepare_memory(D, R, nm, 3);
        const MatchingDecoder& dec = *M.decoder;
        const auto& uf = dec.union_find();
        const int nd = dec.n_detectors(), V = nd + 1;
        Xoshiro256 rng(4242);
        DecodeScratch ws;
        int bad_w = 0, bad_uf = 0, trials = 300;
        double worst = 0;
        std::vector<double> dist(std::size_t(V) * std::size_t(V));
        for (int t = 0; t < trials; ++t) {
            std::vector<std::pair<int, double>> sc;
            const int ne = 1 + int(rng.uniform() * 20);
            for (int k = 0; k < ne; ++k) sc.push_back({int(rng.uniform() * double(uf.edge.size())), rng.uniform() < 0.3 ? 0.0 : rng.uniform()});
            std::vector<int> def;
            for (int d = 0; d < nd; ++d) if (rng.uniform() < 0.15 && def.size() < 10) def.push_back(d);
            const double INF = 1e300;
            std::fill(dist.begin(), dist.end(), INF);
            for (int i = 0; i < V; ++i) dist[std::size_t(i) * std::size_t(V) + std::size_t(i)] = 0;
            std::vector<double> f(uf.edge.size(), 1.0);
            for (auto [e, x] : sc) f[std::size_t(e)] = std::min(f[std::size_t(e)], x);
            for (std::size_t e = 0; e < uf.edge.size(); ++e) {
                const auto& E = uf.edge[e];
                double& a = dist[std::size_t(E.u) * std::size_t(V) + std::size_t(E.v)];
                a = std::min(a, f[e] * E.w); dist[std::size_t(E.v) * std::size_t(V) + std::size_t(E.u)] = a;
            }
            for (int k = 0; k < V; ++k) for (int i = 0; i < V; ++i) for (int j = 0; j < V; ++j) {
                const double c = dist[std::size_t(i) * std::size_t(V) + std::size_t(k)] + dist[std::size_t(k) * std::size_t(V) + std::size_t(j)];
                double& x = dist[std::size_t(i) * std::size_t(V) + std::size_t(j)];
                if (c < x) x = c;
            }
            // minimum-weight matching of def (or to the boundary) by DP over subsets
            const int k = int(def.size());
            std::vector<double> cost(std::size_t(1) << k, 0.0);
            for (std::size_t mask = 1; mask < cost.size(); ++mask) {
                int i = 0; while (!((mask >> i) & 1)) ++i;
                const std::size_t rest = mask & ~(std::size_t(1) << i);
                double best = dist[std::size_t(def[std::size_t(i)]) * std::size_t(V) + std::size_t(nd)] + cost[rest];
                for (int j = i + 1; j < k; ++j) if ((rest >> j) & 1)
                    best = std::min(best, dist[std::size_t(def[std::size_t(i)]) * std::size_t(V) + std::size_t(def[std::size_t(j)])] + cost[rest & ~(std::size_t(1) << j)]);
                cost[mask] = best;
            }
            dec.decode_scaled(def, sc, ws);
            const double err = std::abs(ws.weight - cost.back());
            worst = std::max(worst, err);
            if (err > 1e-9 * std::max(1.0, cost.back())) ++bad_w;
            // union-find with the same modified lengths clears the syndrome
            dec.decode_scaled(def, sc, ws, DecodeMethod::UnionFind);
            std::vector<std::uint8_t> on(std::size_t(nd), 0), got(std::size_t(nd), 0);
            for (int d : def) on[std::size_t(d)] = 1;
            for (int e : ws.uf.correction) for (int v : {uf.edge[std::size_t(e)].u, uf.edge[std::size_t(e)].v}) if (v < nd) got[std::size_t(v)] ^= 1;
            if (got != on) ++bad_uf;
        }
        std::printf("  %d random defect sets on randomly re-weighted graphs (d=3): matching weight differs from brute force on %d (worst %.1e); union-find leaves %d uncleared\n",
                    trials, bad_w, worst, bad_uf);
        CHECK(bad_w == 0 && bad_uf == 0, "re-weighted decoding: %d weight mismatches, %d uncleared", bad_w, bad_uf);
    }
    // 2. herald statistics and the logical gain (d = 5, CZ leakage of the measure qubits 2e-3, leakage-removing reset)
    {
        const auto D = surface_chip(5); const auto R = analyze(D);
        const auto nm = stressed_leakage(D, R, 2e-3);
        const auto M = prepare_memory(D, R, nm, 10);
        const long shots = 20000;
        FrameRunOptions on; on.herald.on = true;
        FrameRunOptions blind = on; blind.herald.use_in_decoder = false;
        FrameRunOptions perfect = on; perfect.herald.fidelity = 1; perfect.herald.false_rate = 0;
        FrameRunOptions none = on; none.herald.fidelity = 0; none.herald.false_rate = 0;
        FrameRunOptions det = on; det.herald.footprint = false; det.herald.touch_factor = 0.4;
        const auto rb = run_memory_frame(M, shots, 515, blind), ru = run_memory_frame(M, shots, 515, on);
        const auto rp = run_memory_frame(M, shots, 515, perfect), rn = run_memory_frame(M, shots, 515, none);
        const auto rd = run_memory_frame(M, shots, 515, det);
        FrameRunOptions ufo = on; ufo.method = DecodeMethod::UnionFind;
        FrameRunOptions ufb = blind; ufb.method = DecodeMethod::UnionFind;
        const auto rub = run_memory_frame(M, shots, 515, ufb), ruu = run_memory_frame(M, shots, 515, ufo);
        std::printf("  d=5, 10 rounds, %ld shots: %.2f leakage events, %.2f heralds per shot (%.0f%% leaked); perfect heralds %.2f (%.0f%% leaked), none %.2f\n",
                    shots, ru.leak_events_per_shot, ru.heralds_per_shot, 100 * ru.herald_true_frac, rp.heralds_per_shot, 100 * rp.herald_true_frac, rn.heralds_per_shot);
        std::printf("  same syndromes: exact matching %ld failures without heralds, %ld with (footprint rule), %ld with the detector rule; union-find %ld -> %ld\n",
                    rb.logical_errors, ru.logical_errors, rd.logical_errors, rub.logical_errors, ruu.logical_errors);
        CHECK(rp.herald_true_frac == 1.0 && rn.heralds_per_shot == 0.0 && ru.herald_true_frac > 0.7 && ru.herald_true_frac < 1.0,
              "herald sampling (perfect: %.3f leaked, none: %.3f per shot, default: %.3f leaked)", rp.herald_true_frac, rn.heralds_per_shot, ru.herald_true_frac);
        CHECK(10 * ru.logical_errors < 9 * rb.logical_errors, "heralds do not lower the logical error by 10%% (%ld vs %ld failures)", ru.logical_errors, rb.logical_errors);
        CHECK(10 * ruu.logical_errors < 9 * rub.logical_errors, "heralds do not help the union-find decoder by 10%% (%ld vs %ld failures)", ruu.logical_errors, rub.logical_errors);
    }
}

// ---------------------------------------------------------------- impact injection and veto
static void test_impact_veto() {
    std::puts("[impact veto] per-shot impact injection (deterministic, separate stream), veto calibration, residual floor");
    const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
    const auto M = prepare_memory(D, R, nm, 8);
    const long shots = 20000;
    // 1. an impact with no excess noise changes nothing: the injection draws come from their own stream
    FrameRunOptions clean;
    FrameRunOptions null_hit; null_hit.impact.p_shot = 0.5; null_hit.impact.burst.t1_peak_us = 0; null_hit.impact.burst.tphi_peak_us = 0;
    const auto r0 = run_memory_frame(M, shots, 61, clean), rn = run_memory_frame(M, shots, 61, null_hit);
    const double fh = double(rn.impacts) / double(shots);
    std::printf("  impacts without excess noise in %.3f of shots: %ld vs %ld failures, detection %.6f vs %.6f\n", fh, rn.logical_errors, r0.logical_errors,
                rn.detection_rate, r0.detection_rate);
    CHECK(rn.logical_errors == r0.logical_errors && rn.detection_rate == r0.detection_rate && rn.mean_defects == r0.mean_defects,
          "impact injection disturbs the shot's own noise stream");
    CHECK(std::abs(fh - 0.5) < 4 * std::sqrt(0.25 / double(shots)), "impact rate %.4f, expected 0.5", fh);
    // 2. veto calibrated on burst-free shots holds its false-veto rate on fresh burst-free shots
    VetoSpec V;
    const auto cal = calibrate_veto(M, 40000, 63, 2e-3, V);
    V.on = true; V.threshold = cal.threshold;
    FrameRunOptions vclean; vclean.veto = V;
    const auto rf = run_memory_frame(M, 40000, 65, vclean);
    const double fv = double(rf.discarded) / double(rf.shots);
    std::printf("  veto (%d rounds x radius %.0f): threshold %d, false-veto rate %.2e on the calibration shots, %.2e on fresh ones (target 2e-3)\n",
                V.window_rounds, V.radius, V.threshold, cal.false_rate, fv);
    CHECK(cal.false_rate <= 2e-3 && fv < 2e-3 + 4 * std::sqrt(2e-3 / 40000.0), "false-veto rate %.2e (calibration %.2e)", fv, cal.false_rate);
    // 3. chip-wide T1 impacts in 5% of shots: the veto catches most and brings the kept error back near the clean one
    FrameRunOptions inj; inj.impact.p_shot = 0.05;
    FrameRunOptions vet = inj; vet.veto = V;
    const auto ri = run_memory_frame(M, shots, 67, inj), rv = run_memory_frame(M, shots, 67, vet), rc = run_memory_frame(M, shots, 67, vclean);
    const double caught = double(rv.impacts_vetoed) / double(std::max(rv.impacts, 1L));
    std::printf("  %ld chip-wide impacts: %.0f%% caught; eps %.5f (no impacts), %.5f (impacts, no veto), %.5f (impacts, veto; %.1f%% discarded)\n",
                rv.impacts, 100 * caught, rc.eps_kept, ri.eps_per_round, rv.eps_kept, 100.0 * double(rv.discarded) / double(shots));
    std::printf("  failures in impacted shots: %ld without the veto, %ld kept by it\n", ri.errors_impacted, rv.errors_impacted_kept);
    CHECK(caught > 0.6, "veto catches only %.0f%% of chip-wide impacts", 100 * caught);
    // (impacts late in a shot leave too few rounds to be seen: the residual excess comes from those)
    CHECK(ri.eps_per_round > 1.5 * rc.eps_kept && rv.eps_kept - rc.eps_kept < 0.3 * (ri.eps_per_round - rc.eps_kept),
          "veto does not remove most of the impact excess (%.5f, %.5f, %.5f)", rc.eps_kept, ri.eps_per_round, rv.eps_kept);
    // 4. a local impact (sigma 1 lattice unit) does less damage than a chip-wide one
    FrameRunOptions loc = inj; loc.impact.burst.radius = 1.0;
    const auto rl = run_memory_frame(M, shots, 67, loc);
    std::printf("  failures per impact: chip-wide %.3f, local (sigma 1) %.3f\n", double(ri.errors_impacted) / double(ri.impacts), double(rl.errors_impacted) / double(rl.impacts));
    CHECK(rl.errors_impacted * ri.impacts < ri.errors_impacted * rl.impacts, "a local impact is not milder than a chip-wide one");
    // 5. the floor over the burst profile with the veto (d = 3, gap-engineered dephasing burst, 1 ms)
    {
        const auto D3 = surface_chip(3); const auto R3 = analyze(D3); const auto nm3 = noise_model(D3, R3);
        BurstSpec B; B.tau_ms = 1; B.t1_peak_us = 0; B.tphi_peak_us = 0.1;
        VetoSpec V3; V3.threshold = calibrate_veto(prepare_memory(D3, R3, nm3, 6), 20000, 69, 2e-3, V3).threshold;
        const auto bv = burst_floor_veto(D3, R3, nm3, B, 6, 6000, 71, V3, 'X');
        std::printf("  d=3 dephasing burst: floor %.2e without the veto, %.2e with shot-by-shot discard (%.3f ms discarded per impact); latched:",
                    bv.floor_per_round, bv.floor_per_round_veto, bv.discarded_ms);
        bool mono = true;
        for (std::size_t i = 0; i < bv.latched.size(); ++i) {
            std::printf(" %.0ftau %.1e", bv.latched[i].hold_tau, bv.latched[i].floor_per_round);
            if (i && bv.latched[i].floor_per_round > bv.latched[i - 1].floor_per_round) mono = false;
        }
        std::printf("\n");
        CHECK(bv.floor_per_round_veto <= bv.floor_per_round && bv.discarded_ms > 0 && mono && bv.latched.back().floor_per_round < 0.01 * bv.floor_per_round,
              "veto floor %.2e vs %.2e, discarded %.3f ms, latched monotone %d", bv.floor_per_round_veto, bv.floor_per_round, bv.discarded_ms, int(mono));
    }
}

// ---------------------------------------------------------------- streaming union-find
static void test_streaming_uf() {
    std::puts("[streaming] sliding-window union-find: clears the syndrome, one window = batch, logical error near batch, latency per round");
    long syndromes = 0, uncleared = 0, singles = 0, single_wrong = 0;
    double p_wrong = 0, p_all = 0;
    // 1. corrections clear the syndrome; single faults the batch decoder gets right are decoded right
    for (int d : {3, 5}) {
        const auto D = surface_chip(d); const auto R = analyze(D); const auto nm = noise_model(D, R);
        const auto M = prepare_memory(D, R, nm, 9);
        const int nd = int(M.exp.detectors.size());
        for (auto [W, C] : std::vector<std::pair<int, int>>{{2, 1}, {4, 2}, {6, 3}}) {
            const StreamingUnionFind S(nd, M.dem.edges, M.exp.detector_round, W, C);
            StreamScratch ss; DecodeScratch ws;
            Xoshiro256 rng(77 + std::uint64_t(d * 10 + W));
            std::vector<std::uint8_t> on(static_cast<std::size_t>(nd)), got(static_cast<std::size_t>(nd));
            std::vector<int> def;
            for (int trial = 0; trial < 1000; ++trial) {
                std::fill(on.begin(), on.end(), 0);
                if (trial % 2) { for (int k = 0; k < nd; ++k) on[std::size_t(k)] = rng.uniform() < 0.08; }
                else for (const auto& e : M.dem.edges) if (rng.uniform() < 3 * e.p) { on[std::size_t(e.a)] ^= 1; if (e.b >= 0) on[std::size_t(e.b)] ^= 1; }
                def.clear(); for (int k = 0; k < nd; ++k) if (on[std::size_t(k)]) def.push_back(k);
                S.decode(def, ss);
                std::fill(got.begin(), got.end(), 0);
                for (auto [a, b] : ss.committed) { got[std::size_t(a)] ^= 1; if (b >= 0) got[std::size_t(b)] ^= 1; }
                ++syndromes; if (got != on) ++uncleared;
            }
            for (const auto& f : enumerate_faults(M.frame, M.exp)) {
                if (f.dets.empty() || f.dets.size() > 2) continue;
                if (M.decoder->decode(f.dets, ws, DecodeMethod::UnionFind) != int(f.obs)) continue;
                ++singles; p_all += f.p;
                if (S.decode(f.dets, ss) != int(f.obs)) { ++single_wrong; p_wrong += f.p; }
            }
        }
    }
    std::printf("  %ld syndromes (d = 3, 5; windows 2/1, 4/2, 6/3 rounds): %ld not cleared; %ld single faults: %ld decoded wrongly "
                "(%.1e of their total probability)\n", syndromes, uncleared, singles, single_wrong, p_wrong / std::max(p_all, 1e-300));
    CHECK(uncleared == 0 && single_wrong == 0, "streaming union-find: %ld uncleared syndromes, %ld single faults wrong (probability share %.1e)",
          uncleared, single_wrong, p_wrong / std::max(p_all, 1e-300));
    // 2. logical error against batch union-find on the same syndromes (d = 5, 15 rounds)
    const auto D = surface_chip(5); const auto R = analyze(D); const auto nm = noise_model(D, R);
    const auto M = prepare_memory(D, R, nm, 15);
    const int nd = int(M.exp.detectors.size());
    const long shots = 30000;
    FrameRunOptions bo; bo.method = DecodeMethod::UnionFind;
    const auto rb = run_memory_frame(M, shots, 811, bo);
    const StreamingUnionFind whole(nd, M.dem.edges, M.exp.detector_round, 17, 16);   // one window: the batch decoder
    FrameRunOptions wo; wo.stream = &whole;
    const auto rw = run_memory_frame(M, shots, 811, wo);
    std::printf("  d=5, 15 rounds, %ld shots: batch union-find %ld failures; one window over the whole shot %ld\n", shots, rb.logical_errors, rw.logical_errors);
    CHECK(whole.n_windows() == 1 && rw.logical_errors == rb.logical_errors, "a single window does not reproduce batch union-find (%ld vs %ld)", rw.logical_errors, rb.logical_errors);
    double worst = 0;
    for (auto [W, C] : std::vector<std::pair<int, int>>{{2, 1}, {3, 1}, {4, 2}, {6, 3}, {10, 5}}) {
        const StreamingUnionFind S(nd, M.dem.edges, M.exp.detector_round, W, C);
        FrameRunOptions so; so.stream = &S;
        const auto rs = run_memory_frame(M, shots, 811, so);
        const double ratio = double(rs.logical_errors) / double(std::max(rb.logical_errors, 1L));
        if (W >= 4 && W >= 2 * C) worst = std::max(worst, ratio);
        std::printf("  window %2d / commit %d rounds (%d windows): %ld failures (x%.2f of batch); per window %.2f us mean, %.1f us max -> %.2f us per committed round\n",
                    W, C, S.n_windows(), rs.logical_errors, ratio, rs.stream_window_us_mean, rs.stream_window_us_max, rs.stream_window_us_mean / C);
    }
    CHECK(worst < 1.25, "streaming union-find with window >= 4 rounds and >= 2 x commit is %.2fx batch union-find's logical error", worst);
}

// ---------------------------------------------------------------- v0.4: memory in both bases
static void test_memory_bases() {
    std::puts("[memory-X] X-basis memory: determinism, distance, logical operators (exact emulator), agreement with the Z basis");
    for (int d : {3, 5}) {
        const auto D = surface_chip(d); const auto R = analyze(D); const auto nm = noise_model(D, R);
        auto quiet = nm; quiet.disable("all");
        const auto Q = prepare_memory(D, R, quiet, 3, 'X');
        const auto r0 = run_memory_frame(Q, 500, 1);
        const auto M = prepare_memory(D, R, nm, d == 3 ? 6 : 4, 'X');
        DetectorErrorModel U = M.dem; for (auto& e : U.edges) e.p = 0.01;
        MatchingDecoder dec(U);
        const auto faults = enumerate_faults(M.frame, M.exp);
        int wrong = 0; DecodeScratch ws;
        for (const auto& f : faults) wrong += dec.decode(f.dets, ws) != int(f.obs);
        std::printf("  d=%d memory-X: noiseless failures %ld, detection %.3g; %zu single faults, %d uncorrected; hyperedges %d, undetectable logicals %d\n",
                    d, r0.logical_errors, r0.detection_rate, faults.size(), wrong, M.dem.hyperedge_faults, M.dem.undetectable_logical);
        CHECK(r0.logical_errors == 0 && r0.detection_rate == 0, "noiseless memory-X is not deterministic");
        CHECK(M.dem.hyperedge_faults == 0 && M.dem.undetectable_logical == 0, "memory-X DEM is not graph-like");
        CHECK(wrong == 0, "memory-X: %d single faults uncorrected (full circuit distance needs 0)", wrong);
    }
    // Logical operators, checked with the exact state-vector emulator: a chain of Z along the logical-Z row (or of X
    // along the logical-X column), applied between two rounds, fires no detector; it flips the observable of the
    // other basis and leaves its own alone.
    {
        const auto D = surface_chip(3); const auto R = analyze(D);
        auto nm = noise_model(D, R); nm.disable("all");
        auto chain = [&](char basis, char pauli) {
            const SurfaceCode S = rotated_surface_code(3);
            MemoryExperiment E = memory_experiment(S, 3, basis);
            auto& gs = E.circuit.gates;
            std::size_t i = 0;
            while (gs[i].op != Op::MEASURE) ++i;                     // the measure-qubit readout that ends round 0
            while (gs[i].op != Op::RESET) ++i;                       // their resets
            while (i < gs.size() && gs[i].op == Op::RESET) ++i;      // first gate of round 1
            std::vector<Gate> ins;
            for (int k : (pauli == 'Z' ? S.logical_z : S.logical_x)) {
                Gate g; g.op = Op::U; g.q = {S.data[std::size_t(k)].qubit}; g.u = pauli == 'Z' ? gates::Z() : gates::X(); g.name = pauli == 'Z' ? "z" : "x";
                ins.push_back(g);
            }
            gs.insert(gs.begin() + std::ptrdiff_t(i), ins.begin(), ins.end());
            CompileOptions co; co.layout = "trivial";
            const auto prog = compile(E.circuit, chip_graph(D, R), co);
            const auto recs = run_records(prog, nm, 6, 5, false);
            int fired = 0, flipped = 0;
            for (const auto& rc : recs) {
                for (const auto& det : E.detectors) { std::uint8_t v = 0; for (int c : det) v ^= rc[std::size_t(c)]; fired += v; }
                std::uint8_t o = 0; for (int c : E.observable) o ^= rc[std::size_t(c)];
                flipped += o;
            }
            return std::make_pair(fired, flipped);
        };
        const auto xz = chain('X', 'Z'), xx = chain('X', 'X'), zx = chain('Z', 'X'), zz = chain('Z', 'Z');
        std::printf("  logical chains on 6 exact shots (detectors fired, observable flips): memory-X + Z_L %d,%d; memory-X + X_L %d,%d; memory-Z + X_L %d,%d; memory-Z + Z_L %d,%d\n",
                    xz.first, xz.second, xx.first, xx.second, zx.first, zx.second, zz.first, zz.second);
        CHECK(xz.first == 0 && xz.second == 6 && zx.first == 0 && zx.second == 6, "a logical operator must flip the conjugate memory without firing detectors");
        CHECK(xx.first == 0 && xx.second == 0 && zz.first == 0 && zz.second == 0, "a logical operator must act trivially on its own memory");
    }
    const auto D = surface_chip(3); const auto R = analyze(D);
    // with the full noise model the two bases agree to within statistics and a modest asymmetry
    const auto nf = noise_model(D, R);
    const auto fz = run_memory_frame(prepare_memory(D, R, nf, 10, 'Z'), 40000, 5);
    const auto fx = run_memory_frame(prepare_memory(D, R, nf, 10, 'X'), 40000, 5);
    std::printf("  full noise: memory-Z eps %.4f%%, memory-X eps %.4f%%\n", 100 * fz.eps_per_round, 100 * fx.eps_per_round);
    CHECK(fx.eps_per_round > 0.6 * fz.eps_per_round && fx.eps_per_round < 1.6 * fz.eps_per_round, "X and Z memories disagree by more than 60%%");
    // Error bar of the fitted logical error per round: against the scatter of the fit over synthetic experiments drawn
    // from binomial statistics, and against linear propagation through the fit by finite differences.
    {
        const std::vector<int> rounds{5, 10, 15, 20};
        const double eps = 2e-3, A = 0.98; const long N = 100000;
        auto make = [&](const std::vector<double>& p) {
            std::vector<MemoryResult> rs;
            for (std::size_t i = 0; i < rounds.size(); ++i) {
                MemoryResult r; r.rounds = rounds[i]; r.shots = N; r.p_logical = p[i];
                r.p_logical_sem = std::sqrt(p[i] * (1 - p[i]) / double(N)); r.eps_per_round = eps_from_PL(p[i], rounds[i]);
                rs.push_back(r);
            }
            return rs;
        };
        std::vector<double> p0;
        for (int r : rounds) p0.push_back(0.5 * (1 - A * std::pow(1 - 2 * eps, r)));
        const auto base = make(p0);
        const double err = fit_eps_per_round_err(base);
        double lin = 0;
        for (std::size_t i = 0; i < p0.size(); ++i) {
            auto up = p0, dn = p0; const double h = 1e-3 * base[i].p_logical_sem;
            up[i] += h; dn[i] -= h;
            const double dedp = (fit_eps_per_round(make(up)) - fit_eps_per_round(make(dn))) / (2 * h);
            lin += std::pow(dedp * base[i].p_logical_sem, 2);
        }
        lin = std::sqrt(lin);
        std::mt19937_64 rng(17);
        const int reps = 4000; double m1 = 0, m2 = 0;
        for (int k = 0; k < reps; ++k) {
            std::vector<double> p;
            for (double q : p0) p.push_back(double(std::binomial_distribution<long>(N, q)(rng)) / double(N));
            const double e = fit_eps_per_round(make(p)); m1 += e; m2 += e * e;
        }
        m1 /= reps; const double scatter = std::sqrt(m2 / reps - m1 * m1);
        const auto& last = base.back();
        const double last_only = last.p_logical_sem / last.rounds / std::pow(1 - 2 * last.p_logical, double(last.rounds - 1) / last.rounds);
        std::printf("  fit error: fitted eps %.5f%% (true %.5f%%); error bar %.3e, scatter of %d synthetic fits %.3e, linear propagation %.3e; last point alone %.3e\n",
                    100 * fit_eps_per_round(base), 100 * eps, err, reps, scatter, lin, last_only);
        CHECK(std::abs(fit_eps_per_round(base) / eps - 1) < 1e-9, "the fit does not recover the logical error per round of exact data");
        CHECK(std::abs(err / lin - 1) < 1e-4, "fit error bar %.4e vs linear propagation %.4e", err, lin);
        CHECK(std::abs(err / scatter - 1) < 0.06, "fit error bar %.4e vs scatter of synthetic fits %.4e", err, scatter);
        const std::vector<MemoryResult> one{base.back()};
        CHECK(std::abs(fit_eps_per_round_err(one) / last_only - 1) < 1e-12, "single-run error bar");
    }
}

// ---------------------------------------------------------------- v0.4: correlated bursts
static void test_bursts() {
    std::puts("[bursts] correlated error bursts: excess rates, severity ladder, logical-error floor");
    const auto D = surface_chip(3); const auto R = analyze(D); const auto nm = noise_model(D, R);
    BurstSpec B;                                                    // 0.1 per s, 25 ms, peak excess T1 1 us
    const auto n0 = burst_noise(nm, B, 0.0), n1 = burst_noise(nm, B, 1.0), nh = burst_noise(nm, B, 0.5);
    double e0 = 0, e1 = 0, eh = 0;
    for (std::size_t i = 0; i < nm.T1_ns.size(); ++i) {
        e0 = std::max(e0, std::abs(n0.T1_ns[i] / nm.T1_ns[i] - 1));
        e1 = std::max(e1, std::abs(1.0 / n1.T1_ns[i] - (1.0 / nm.T1_ns[i] + 1.0 / 1000.0)) * nm.T1_ns[i]);
        eh = std::max(eh, std::abs(1.0 / nh.T1_ns[i] - (1.0 / nm.T1_ns[i] + 0.5 / 1000.0)) * nm.T1_ns[i]);
    }
    CHECK(e0 < 1e-12 && e1 < 1e-9 && eh < 1e-9, "burst_noise excess decay rate (errors %.1e %.1e %.1e)", e0, e1, eh);
    CHECK(n1.Tphi_ns == nm.Tphi_ns, "a T1 burst must leave pure dephasing alone");
    const auto b = burst_floor(D, R, nm, B, 6, 20000, 9);
    bool mono = true;
    for (std::size_t i = 0; i + 1 < b.curve.size(); ++i) if (b.curve[i + 1].eps > b.curve[i].eps + 3 * (b.curve[i].eps_sem + b.curve[i + 1].eps_sem)) mono = false;
    std::printf("  ladder of %zu severities: eps %.4f at the peak -> %.5f at the last point (quiet %.5f); P(fail) %.3f, floor %.2e per round, unprotected %.0f ms\n",
                b.curve.size(), b.curve.front().eps, b.curve.back().eps, b.eps0, b.p_fail, b.floor_per_round, b.dead_time_ms);
    CHECK(mono && b.curve.front().eps > 20 * b.eps0, "logical error must rise with burst severity");
    CHECK(std::abs(b.p_fail - 0.5) < 1e-6, "a 25 ms burst randomises the logical qubit");
    CHECK(std::abs(b.floor_per_round - B.rate_hz * b.round_ns * 1e-9 * b.p_fail) < 1e-18, "floor = rate x round time x P(fail)");
    CHECK(b.dead_time_ms > B.tau_ms && b.dead_time_ms < 12 * B.tau_ms, "unprotected time should be a few recovery times (%.1f ms)", b.dead_time_ms);
    // the same burst a million times shorter costs a millionth of the flips: the integral is linear in the duration
    BurstSpec S = B; S.tau_ms = B.tau_ms * 1e-6;
    const auto bs = burst_floor(D, R, nm, S, 6, 20000, 9);
    std::printf("  25 ns burst: expected logical flips %.2e (25 ms burst: %.2e), P(fail) %.2e\n", bs.rounds_lost, b.rounds_lost, bs.p_fail);
    CHECK(std::abs(bs.rounds_lost / b.rounds_lost - 1e-6) < 1e-9 && bs.p_fail < 0.05, "burst integral scales with the recovery time");
}

// ---------------------------------------------------------------- v0.4: measured data -> model -> hardware gate
static void test_characterisation() {
    std::puts("[hardware gate] characterisation file round trip, measured-parameter model, gate evaluation on a virtual chip");
    const auto D = surface_chip(3); auto R = analyze(D);
    {   // give the pairs distinct control errors and leakage so every path is exercised
        int k = 0;
        for (auto& p : R.pairs) if (p.designed) { p.err_2q_ctrl = 1e-3 + 1e-4 * (k % 5); p.leak_2q = 5e-5 * (1 + k % 3); p.leak_qubit = (k % 2) ? p.a : p.b; ++k; }
    }
    // (a) a noiseless, scatter-free "measurement" of the design must give the design's own noise model back
    SynthOptions exact; exact.t1_sigma = exact.tphi_sigma = exact.ro_sigma = exact.cz_sigma = exact.leak_sigma = exact.thermal_sigma = 0;
    exact.tls_fraction = 0; exact.zz_sigma_kHz = 0; exact.est_noise = 0; exact.zz_est_kHz = 0; exact.rounds = {4}; exact.shots = 2000; exact.both_bases = false;
    const auto C0 = synthesize_characterisation(D, R, exact);
    const auto text = [&] { std::ostringstream os; characterisation_to_json(C0).dump(os); return os.str(); }();
    const auto C1 = characterisation_from_json(Json::parse(text));
    CHECK(C1.qubits.size() == D.qubits.size() && C1.pairs.size() == C0.pairs.size() && C1.memory.size() == 1 && C1.source == "synthetic", "file round trip: sizes");
    const auto MM = model_from_measurement(D, R, C1);
    const auto ref = noise_model(D, R);
    double eT1 = 0, eTp = 0, eRo = 0, eTh = 0, e2 = 0, eLk = 0, eZZ = 0;
    for (std::size_t i = 0; i < ref.T1_ns.size(); ++i) {
        eT1 = std::max(eT1, std::abs(MM.nm.T1_ns[i] / ref.T1_ns[i] - 1));
        eTp = std::max(eTp, std::isfinite(ref.Tphi_ns[i]) ? std::abs(MM.nm.Tphi_ns[i] / ref.Tphi_ns[i] - 1) : (std::isfinite(MM.nm.Tphi_ns[i]) ? 1.0 : 0.0));
        eRo = std::max({eRo, std::abs(MM.nm.ro01[i] - ref.ro01[i]), std::abs(MM.nm.ro10[i] - ref.ro10[i])});
        eTh = std::max(eTh, std::abs(MM.nm.p_thermal[i] - ref.p_thermal[i]));
    }
    for (const auto& [key, v] : ref.err_2q_pair) e2 = std::max(e2, std::abs(MM.nm.err_2q_pair.at(key) - v));
    for (const auto& [key, v] : ref.leak_cz) { const auto& m = MM.nm.leak_cz.at(key); eLk = std::max(eLk, std::abs(m.second - v.second) + (m.first == v.first ? 0.0 : 1.0)); }
    for (const auto& [key, v] : ref.zz_GHz) eZZ = std::max(eZZ, std::abs(MM.nm.zz_GHz.at(key) - v));
    std::printf("  measured-parameter model vs design: T1 %.1e, T_phi %.1e, readout %.1e, thermal %.1e, CZ control %.1e, leakage %.1e, ZZ %.1e, 1q control %.1e\n",
                eT1, eTp, eRo, eTh, e2, eLk, eZZ, std::abs(MM.nm.err_1q_ctrl - ref.err_1q_ctrl));
    CHECK(eT1 < 1e-9 && eTp < 1e-6 && eRo < 1e-12 && eTh < 1e-12, "coherence and readout parameters do not round-trip");
    CHECK(e2 < 1e-12 && eLk < 1e-15 && eZZ < 1e-15, "pair parameters do not round-trip");
    CHECK(std::abs(MM.nm.err_1q_ctrl - ref.err_1q_ctrl) < 1e-12 && std::abs(MM.nm.reset_err - ref.reset_err) < 1e-12, "single-qubit control / reset error do not round-trip");
    CHECK(MM.qubits_measured == int(D.qubits.size()) && MM.pairs_measured == MM.pairs_designed, "measured counts");
    // (b) a virtual chip with fabrication-like scatter passes the model-vs-measurement criterion ...
    SynthOptions so; so.seed = 3; so.shots = 20000;
    const auto Cv = synthesize_characterisation(D, R, so);
    const auto G = evaluate_hardware_gate("alpha", D, R, Cv, 20000, 77);
    double ratio = 0; int n_eval = 0;
    for (const auto& c : G.checks) { n_eval += c.evaluated; if (c.name.rfind("Model vs measurement", 0) == 0) ratio = c.value; }
    std::printf("  virtual chip: %d of %zu criteria evaluated, %d met; model vs measurement %.3f (Z and X, 3 round counts each)\n", n_eval, G.checks.size(), G.passed, ratio);
    CHECK(G.missing == 0 && n_eval == int(G.checks.size()) && G.comparisons.size() == 2, "every Alpha criterion should be evaluated from a full data file");
    CHECK(ratio >= 1.0 && ratio < 1.25, "model with measured parameters should reproduce the measured logical error (ratio %.3f)", ratio);
    // ... and fails it when an error source acts during the experiment that the parameter list does not contain
    so.hidden_dephasing_us = 8;
    const auto Gh = evaluate_hardware_gate("alpha", D, R, synthesize_characterisation(D, R, so), 20000, 77);
    double ratio_h = 0; for (const auto& c : Gh.checks) if (c.name.rfind("Model vs measurement", 0) == 0) ratio_h = c.value;
    std::printf("  hidden dephasing (8 us): model vs measurement %.2f -> gate %s\n", ratio_h, Gh.pass ? "passed" : "not passed");
    CHECK(ratio_h > 2.0 && !Gh.pass, "a hidden error source must fail the model-vs-measurement criterion");
    // (c) thresholds and missing data
    Characterisation Cb = Cv;
    for (auto& q : Cb.qubits) q.T1_us *= 0.4;                        // median T1 well below 50 us
    for (auto& p : Cb.pairs) p.leak_cz = chz::nan();                 // leakage RB not done
    const auto Gb = evaluate_hardware_gate("alpha", D, R, Cb, 5000, 77);
    bool t1_fail = false, leak_missing = false;
    for (const auto& c : Gb.checks) {
        if (c.name.rfind("Median T1", 0) == 0) t1_fail = c.evaluated && !c.pass;
        if (c.name.rfind("Median CZ leakage", 0) == 0) leak_missing = !c.evaluated;
    }
    CHECK(t1_fail && leak_missing && Gb.missing == 1 && !Gb.pass, "threshold failure and missing data are reported");
    bool threw = false;
    try { characterisation_from_json(Json::parse("{\"schema\": \"something.else\"}")); } catch (const std::exception&) { threw = true; }
    CHECK(threw, "a file with the wrong schema must be rejected");
}

// ---------------------------------------------------------------- v0.4: designs
static void test_designs_v04() {
    std::puts("[designs] readout plan, Beta-1 rev C, Gamma-1, CZ interaction point at the tuned anharmonicity");
    for (int d : {5, 7}) {
        const auto D = d == 5 ? beta1(5) : gamma1(7);
        std::map<int, int> per; double lo = 1e9, hi = 0;
        for (const auto& q : D.qubits) { ++per[q.feedline]; lo = std::min(lo, q.ro.f_GHz); hi = std::max(hi, q.ro.f_GHz); }
        int mx = 0; for (const auto& kv : per) mx = std::max(mx, kv.second);
        const auto R = analyze(D);
        std::printf("  %s: %zu qubits, %zu couplers, %zu feedlines (at most %d resonators), readout band %.2f-%.2f GHz, collisions %d\n",
                    D.name.c_str(), D.qubits.size(), D.couplings.size(), per.size(), mx, lo, hi, R.n_collisions);
        CHECK(int(D.qubits.size()) == 2 * d * d - 1 && int(D.couplings.size()) == 4 * d * (d - 1), "qubit / coupler count at d=%d", d);
        CHECK(mx <= 9 && int(per.size()) == (d == 5 ? 6 : 11), "feedline plan at d=%d", d);
        CHECK(lo > 6.89 && hi < 7.40 && hi < D.coupler.f_max_GHz - 0.5, "readout band at d=%d", d);
        CHECK(R.n_collisions == 0, "frequency collisions at d=%d", d);
    }
    // |11> = |20> solved at the tuned E_J, against the definition and against the idle-anharmonicity estimate
    const auto A = alpha1(3); const auto R = analyze(A);
    double worst_def = 0, shift_min = 1e9, shift_max = -1e9;
    for (const auto& p : R.pairs) {
        if (!p.designed) continue;
        const auto& qa = R.q[std::size_t(A.index_of(p.a))]; const auto& qb = R.q[std::size_t(A.index_of(p.b))];
        const auto& hi = qa.f01 >= qb.f01 ? qa : qb; const auto& lo = qa.f01 >= qb.f01 ? qb : qa;
        const auto fa = transmon_fa(solve_EJ_for_f01(p.cz_interaction_GHz, hi.EC), hi.EC);
        worst_def = std::max(worst_def, std::abs(fa.f01 + fa.alpha - lo.f01));
        const double shift = (p.cz_interaction_GHz - (lo.f01 - hi.alpha)) * 1e3;      // MHz, vs the idle-alpha estimate
        shift_min = std::min(shift_min, shift); shift_max = std::max(shift_max, shift);
    }
    std::printf("  Alpha-1 CZ interaction points: |f01 + alpha - f_partner| <= %.1e GHz; %.2f to %.2f MHz above the idle-alpha estimate\n", worst_def, shift_min, shift_max);
    CHECK(worst_def < 1e-7, "interaction point does not satisfy f01 + alpha = f_partner");
    CHECK(shift_min > 0.5 && shift_max < 6.0, "tuned-anharmonicity correction outside the expected few MHz");
}

// ---------------------------------------------------------------- layout_geometry: DRC-aware cell and mask clean-up logic
static void test_layout_geometry() {
    std::puts("[layout_geometry] bus knob continuity, shorted staircase end, ground-sliver fill (merges ground only, never island metal), metal connectivity, guarded fill");
    using namespace layout;
    auto count = [](const std::vector<DrcMarker>& m, const char* rule) { int n = 0; for (const auto& x : m) n += x.rule.find(rule) != std::string::npos; return n; };
    {   // bus knob: T-bar length above w_bus, retraction below it, continuous at u = w_bus, round trip, never narrower than the strip
        CouplerGeom c;
        for (double u : {40.0, 12.0, 6.0, 3.0, 0.0, -10.0}) {
            set_bus_knob(c, u);
            CHECK(std::abs(bus_knob(c) - u) < 1e-12, "bus knob round trip at u=%g: %g", u, bus_knob(c));
            CHECK(c.bus_tbar() >= c.w_bus - 1e-12, "bus T-bar %g narrower than the strip", c.bus_tbar());
        }
        set_bus_knob(c, c.w_bus); CHECK(std::abs(c.bus_gap() - c.g_c) < 1e-12 && std::abs(c.bus_tbar() - c.w_bus) < 1e-12, "branches meet at u = w_bus");
        set_bus_knob(c, 0.0); CHECK(std::abs(c.bus_gap() - (c.g_c + kBusRetract * c.w_bus)) < 1e-12, "retraction below w_bus");
        c.l_b = 0.075; c.g_b = -1;                                     // a pre-0.4.1 cell file with a sliver T-bar
        XmonParams q;
        const auto E = coupler_and_bus(q, c, 900, 0, 0, 0, true, 200, q.L);
        double narrow = 1e9; for (const auto& r : E.bus) narrow = std::min({narrow, r.w(), r.h()});
        CHECK(narrow >= c.w_bus - 1e-9, "bus drawn with a %.3f um feature", narrow);
    }
    {   // shorted staircase end: no trench under the last square, the one before trimmed to its edges -> no etch slit
        std::vector<Rect> m, h;
        detail::stair(0, 0, 50, 50, m, h, true);
        const auto etch = detail::etch_of(h, m);
        const auto mk = detail::width_space_check(etch, 4, 4, "s");
        std::printf("  shorted staircase: %zu squares, %zu trench squares, %d width + %d space markers\n", m.size(), h.size(), count(mk, "width"), count(mk, "space"));
        CHECK(h.size() + 1 == m.size() && count(mk, "width") == 0 && count(mk, "space") == 0, "shorted staircase end leaves DRC markers");
        // the last square's outer edges (x = 55, y = 55) must border ground, i.e. the short really reaches the plane
        bool etched = false; for (const auto& r : h) etched |= r.contains(55.5, 50) || r.contains(50, 55.5);
        CHECK(!etched, "the shorted square does not touch the ground plane");
    }
    {   // ground sliver between two cut-outs is merged; a metal island in the same gap is left alone; a wide strip is kept
        std::vector<Rect> holes{{0, 0, 20, 20}, {22, 0, 42, 20},           // 2 um ground strip: filled
                                {100, 0, 120, 20}, {130, 0, 150, 20}};     // 10 um ground strip: kept
        const std::vector<Rect> none;
        const int added = detail::fill_ground_slivers(holes, none, 4.0);
        const auto mk = detail::width_space_check(detail::etch_of(holes, none), 4, 4, "f");
        std::printf("  ground-sliver fill: %d rectangle(s) added, %d space markers after\n", added, count(mk, "space"));
        CHECK(added == 1 && mk.empty(), "sliver fill: %d added, %zu markers", added, mk.size());
        bool wide_kept = true; for (const auto& r : holes) if (r.intersects({121, 1, 129, 19})) wide_kept = false;
        CHECK(wide_kept, "a 10 um ground strip was etched");
        std::vector<Rect> h2{{0, 0, 20, 20}, {22, 0, 42, 20}};
        const std::vector<Rect> island{{20, 5, 22, 15}};                   // a 2 um island touching the gap: not ground
        const int added2 = detail::fill_ground_slivers(h2, island, 4.0);
        CHECK(added2 == 0 && h2.size() == 2, "fill touched a gap containing island metal (%d added)", added2);
    }
    {   // metal connectivity: a closed ring of cut-outs encloses a floating ground pocket; an island inside a cut-out is
        // its own component; etch and points outside the outline have none
        const Rect outline{0, 0, 100, 100};
        const std::vector<Rect> ring{{20, 20, 80, 30}, {20, 70, 80, 80}, {20, 30, 30, 70}, {70, 30, 80, 70}};
        auto M = detail::metal_map(outline, ring, 37);       // columns narrower than the ring: joins across column edges
        std::printf("  metal map of a closed ring: %d components (plane %d, pocket %d, etch %d)\n", M.n_comp, M.at(5, 5), M.at(50, 50), M.at(25, 50));
        CHECK(M.n_comp == 2 && M.at(5, 5) != M.at(50, 50) && M.at(25, 50) == -1 && M.at(95, 95) == M.at(5, 5) && M.at(-1, 5) == -1,
              "closed ring: %d components", M.n_comp);
        std::vector<Rect> open = ring; open[3] = {70, 30, 80, 45}; open.push_back({70, 55, 80, 70});      // 10 um opening
        CHECK(detail::metal_map(outline, open, 37).n_comp == 1, "ring with an opening must leave one component");
    }
    {   // guarded sliver fill: a pocket joined to the plane only by a 3 um neck keeps the neck (the DRC reports it); an
        // ordinary 2 um sliver between two cut-outs elsewhere is still merged
        const Rect outline{0, 0, 200, 100};
        std::vector<Rect> holes{{20, 20, 80, 30}, {20, 70, 80, 80}, {20, 30, 30, 70}, {70, 30, 80, 48}, {70, 51, 80, 70},   // neck y 48..51
                                {120, 20, 150, 80}, {152, 20, 180, 80}};                                                 // 2 um sliver
        const auto G = detail::guarded_fill(outline, holes, {}, 4.0);
        const auto mk = detail::width_space_check(detail::etch_of(holes, {}), 4, 4, "g");
        std::printf("  guarded fill: %d kept, %d undone; pocket on the plane: %s; %zu DRC marker(s) left\n", G.kept, G.undone,
                    G.map.at(50, 50) == G.map.at(5, 5) ? "yes" : "no", mk.size());
        CHECK(G.map.at(50, 50) == G.map.at(5, 5), "guarded fill cut the pocket off the plane");
        CHECK(G.undone >= 1 && G.kept >= 1 && mk.size() == 1, "guarded fill: kept %d, undone %d, markers %zu", G.kept, G.undone, mk.size());
        bool merged = false; for (const auto& h : holes) merged |= h.contains(151, 50);
        CHECK(merged, "the 2 um sliver away from the pocket was not merged");
    }
    {   // ground straps: a pocket inside a closed 10 um ring gets one strap across the ring; a strap never spans a qubit
        // island (here a qubit drawn in the ring's left side forces the strap to another side), and a pocket ringed
        // only by a qubit gets none
        const Rect outline{0, 0, 300, 100};
        const std::vector<Rect> ring{{20, 20, 80, 30}, {20, 70, 80, 80}, {20, 30, 30, 70}, {70, 30, 80, 70}};
        std::vector<std::pair<std::string, std::vector<Rect>>> isl{{"Q1", {{22, 32, 28, 68}}}};
        std::vector<Rect> holes = ring;
        for (const auto& r : ring) holes.push_back(r.shift(100, 0));                                         // second ring
        isl.push_back({"Q2", {{121, 77, 179, 79}, {121, 21, 179, 23}, {121, 21, 123, 79}, {177, 21, 179, 79}}});  // qubit ring 1 um inside the outer edge
        std::vector<Rect> metal; for (const auto& [n, rs] : isl) metal.insert(metal.end(), rs.begin(), rs.end());
        const auto e = detail::etch_of(holes, metal);
        const auto M = detail::metal_map(outline, e);
        const auto I = detail::island_index(isl);
        std::vector<char> eligible(std::size_t(M.n_comp), 1);
        for (const auto& [n, rs] : isl) eligible[std::size_t(M.at(rs.front().cx(), rs.front().cy()))] = 0;
        const int plane = M.at(1, 1);
        const auto S = detail::place_ground_straps(M, I, eligible, plane);
        bool over_q = false; for (const auto& st : S) for (const auto& [n, rs] : isl) for (const auto& r : rs) over_q |= st.span.intersects(r);
        const bool p1 = std::any_of(S.begin(), S.end(), [&](const detail::Strap& st) { return M.at(st.foot_a.cx(), st.foot_a.cy()) == M.at(50, 50) || M.at(st.foot_b.cx(), st.foot_b.cy()) == M.at(50, 50); });
        std::printf("  ground straps: %zu placed; first pocket strapped: %s; any strap over a qubit: %s\n", S.size(), p1 ? "yes" : "no", over_q ? "yes" : "no");
        CHECK(S.size() >= 1 && S.size() <= 2 && p1 && !over_q, "ground straps: %zu, pocket %d, over qubit %d", S.size(), int(p1), int(over_q));
    }
    {   // strap feet stand on ground metal only, with kStrapFootEnclosure around them: a shorted line (R9) fills most of a
        // pocket and is the same metal component as the pocket, yet no foot may touch it; only the 20 um band above it
        // can take a foot. The pocket gets two straps, kStrapSeparation apart.
        const Rect outline{0, 0, 200, 140};
        const std::vector<Rect> ring{{20, 20, 100, 30}, {20, 90, 100, 100}, {20, 30, 30, 90}, {90, 30, 100, 90}};
        const std::vector<std::pair<std::string, std::vector<Rect>>> isl{{"R9", {{34, 34, 86, 70}}}};
        const auto M = detail::metal_map(outline, detail::etch_of(ring, {}));
        const auto I = detail::island_index(isl);
        const std::vector<char> eligible(std::size_t(M.n_comp), 1);
        const auto S = detail::place_ground_straps(M, I, eligible, M.at(1, 1));
        bool on_line = false, enclosed = true;
        for (const auto& st : S) for (const Rect& f : {st.foot_a, st.foot_b}) {
            on_line = on_line || f.intersects(isl[0].second[0]);
            for (const auto& r : ring) enclosed = enclosed && !f.grow(detail::kStrapFootEnclosure - 1e-6).intersects(r);
        }
        const double sep = S.size() == 2 ? std::hypot(S[0].span.cx() - S[1].span.cx(), S[0].span.cy() - S[1].span.cy()) : 0.0;
        std::printf("  strap feet beside a shorted line: %zu straps, any foot on the line: %s, feet enclosed: %s, %.0f um apart\n", S.size(), on_line ? "yes" : "no", enclosed ? "yes" : "no", sep);
        CHECK(S.size() == 2 && !on_line && enclosed && sep >= detail::kStrapSeparation, "strap feet: %zu straps, on the line %d, enclosed %d, %.0f um apart", S.size(), int(on_line), int(enclosed), sep);
    }
}

// ---------------------------------------------------------------- layout_checks: the connectivity, ground, flux and control-line checks on synthetic layouts
static void test_layout_checks() {
    std::puts("[layout_checks] junction ownership, ground quality (plane metal, straps, bumps), flux coupling by duality, control-line window and formulas, qubit cell and launcher");
    using namespace layout;
    // a synthetic mask: islands inside cut-outs, etch = cut-outs minus island metal
    auto finish = [](ChipLayout& L) {
        std::vector<Rect> metal; for (const auto& [n, rs] : L.islands) metal.insert(metal.end(), rs.begin(), rs.end());
        L.etch = detail::etch_of(L.holes, metal);
        analyze_connectivity(L);
    };
    {   // (1) junction ownership. Q1 sits in a cut-out of the plane; C2-3 sits in the etch ring around an enclosed pocket
        // of ground; XY1 is an open line in its own cut-out.
        ChipLayout B; B.outline = {0, 0, 400, 200};
        B.islands = {{"Q1", {{40, 80, 80, 120}}}, {"C2-3", {{214, 90, 226, 110}}}, {"XY1", {{330, 90, 350, 110}}}};
        B.holes = {{20, 60, 100, 140}, {210, 60, 300, 80}, {210, 120, 300, 140}, {210, 80, 240, 120}, {280, 80, 300, 120}, {320, 80, 360, 120}};
        const Rect jq{80, 96, 100, 104}, jc{226, 98, 240, 102}, jx{350, 98, 360, 102};          // Q1 -> plane, C2-3 -> pocket, XY1 -> plane
        const detail::Strap s1{{250, 60, 260, 90}, {250, 80, 260, 90}, {250, 50, 260, 60}, 0, 0};   // pocket -> plane across the ring's bottom side
        const detail::Strap s2{{250, 110, 260, 150}, {250, 110, 260, 120}, {250, 140, 260, 150}, 0, 0};  // and across its top side
        {   // one SQUID on the qubit: lands, on plane metal
            ChipLayout L = B; L.jj = {jq}; L.jj_owner = {0}; finish(L);
            CHECK(L.conn.junctions_landed == 1 && L.conn.qubit_squids == 1 && L.conn.qubit_squids_on_plane == 1 && L.conn.floating_ground == 1,
                  "junction Q1 -> plane: landed %d, on plane %d/%d, floating %d", L.conn.junctions_landed, L.conn.qubit_squids_on_plane, L.conn.qubit_squids, L.conn.floating_ground);
            CHECK(L.conn.junction_owner_errors == 1, "C2-3 has no junction: owner errors %d (expected 1)", L.conn.junction_owner_errors);
        }
        {   // the same region owned by the other island: on the wrong island
            ChipLayout L = B; L.jj = {jq}; L.jj_owner = {1}; finish(L);
            CHECK(L.conn.junctions_landed == 0 && L.conn.junction_owner_errors == 1, "junction on the wrong island: landed %d, owner errors %d", L.conn.junctions_landed, L.conn.junction_owner_errors);
        }
        {   // between a line and the plane: an XY line cannot own a junction, and a qubit does not own a junction that misses it
            ChipLayout L = B; L.jj = {jq, jx}; L.jj_owner = {0, 2}; finish(L);
            CHECK(L.conn.junctions_landed == 1 && L.conn.junction_owner_errors == 2, "junction owned by a line: landed %d of 2, owner errors %d", L.conn.junctions_landed, L.conn.junction_owner_errors);
            L.jj_owner = {0, 1}; finish(L);
            CHECK(L.conn.junctions_landed == 1 && L.conn.junction_owner_errors == 0 && !L.conn.nets_clean(), "junction between a line and the plane: landed %d of 2", L.conn.junctions_landed);
        }
        {   // two junctions on one island, none on the other
            ChipLayout L = B; L.jj = {jq, {80, 106, 100, 112}}; L.jj_owner = {0, 0}; finish(L);
            CHECK(L.conn.junction_owner_errors == 2, "two junctions on Q1, none on C2-3: owner errors %d", L.conn.junction_owner_errors);
        }
        {   // onto a floating fragment: not landed; the 0.4 rule (two different components) would have accepted it
            ChipLayout L = B; L.jj = {jq, jc}; L.jj_owner = {0, 1}; finish(L);
            CHECK(L.conn.junctions_landed == 1 && L.conn.floating_ground == 1 && !L.conn.nets_clean(), "junction onto a floating fragment: landed %d of 2", L.conn.junctions_landed);
            // (2) the pocket reached through one strap: a net, but one cut strap away from floating
            L.straps = {s1}; finish(L);
            const auto g1 = L.grounds.back();
            CHECK(L.conn.nets_clean() && L.conn.junctions_landed == 2 && L.conn.coupler_squids_on_plane == 0 && L.conn.fragments_on_one_strap == 1 &&
                  L.conn.cut_straps == 1 && L.conn.max_strap_depth == 1 && !L.conn.ground_ok() && !L.conn.clean(),
                  "pocket on one strap: nets %d, on one strap %d, cut straps %d, ground ok %d", int(L.conn.nets_clean()), L.conn.fragments_on_one_strap, L.conn.cut_straps, int(L.conn.ground_ok()));
            CHECK(L.grounds.size() == 2 && L.grounds[1].name == "SC2-3" && !L.grounds[1].on_plane && L.grounds[1].straps == 1 && std::abs(L.grounds[1].span_um - 30) < 1e-9 &&
                  L.grounds[0].on_plane && L.grounds[0].straps == 0 && L.strap_returns == std::vector<int>{1}, "ground reference of the coupler SQUID: %d strap(s), %.0f um", g1.straps, g1.span_um);
            // two independent straps: a coupler SQUID may sit there
            L.straps = {s1, s2}; finish(L);
            CHECK(L.conn.fragments_on_one_strap == 0 && L.conn.cut_straps == 0 && L.conn.ground_ok() && L.conn.clean() && L.conn.coupler_squids_on_plane == 0,
                  "pocket on two straps: on one strap %d, ground ok %d", L.conn.fragments_on_one_strap, int(L.conn.ground_ok()));
            // ... but a qubit SQUID may not: its ground must be plane metal itself
            L.islands[1].first = "Q2"; finish(L);
            CHECK(L.conn.nets_clean() && L.conn.qubit_squids_on_plane == 1 && L.conn.qubit_squids == 2 && !L.conn.ground_ok(),
                  "qubit SQUID behind straps: on plane metal %d/%d, ground ok %d", L.conn.qubit_squids_on_plane, L.conn.qubit_squids, int(L.conn.ground_ok()));
            // a ground bump makes the pocket plane metal
            L.ground_bumps = {rect_c(260, 100, 10, 10)}; finish(L);
            CHECK(L.conn.plane_components == 2 && L.conn.qubit_squids_on_plane == 2 && L.conn.strapped_fragments == 0 && L.conn.ground_ok(),
                  "ground bump on the pocket: plane components %d, on plane metal %d/2", L.conn.plane_components, L.conn.qubit_squids_on_plane);
        }
        {   // a chain: pocket A -> pocket B -> plane, two straps each (no cut strap) but depth 2; a third ring makes depth 3
            ChipLayout L; L.outline = {0, 0, 700, 200};
            auto ring = [&](double x) { for (const Rect& r : {Rect{x, 60, x + 90, 80}, Rect{x, 120, x + 90, 140}, Rect{x, 80, x + 20, 120}, Rect{x + 70, 80, x + 90, 120}}) L.holes.push_back(r); };
            ring(100); ring(180); ring(260);                  // three pockets side by side, sharing their walls' edges
            auto strap = [&](double xa, double xb, double y) { L.straps.push_back({{xa - 5, y - 5, xb + 5, y + 5}, rect_c(xa, y, 10, 10), rect_c(xb, y, 10, 10), 0, 0}); };
            strap(90, 130, 90); strap(90, 130, 110);          // plane -> pocket 1 (twice)
            strap(150, 210, 90); strap(150, 210, 110);        // pocket 1 -> pocket 2
            finish(L);
            CHECK(L.conn.strapped_fragments == 2 && L.conn.max_strap_depth == 2 && L.conn.fragments_on_one_strap == 0 && L.conn.floating_ground == 1,
                  "strap chain: %d fragments, depth %d, on one strap %d, floating %d", L.conn.strapped_fragments, L.conn.max_strap_depth, L.conn.fragments_on_one_strap, L.conn.floating_ground);
            strap(230, 290, 90); strap(230, 290, 110);        // pocket 2 -> pocket 3: three straps deep
            finish(L);
            CHECK(L.conn.max_strap_depth == 3 && L.conn.nets_clean() && !L.conn.ground_ok(), "strap chain of depth %d (limit %d): ground ok %d", L.conn.max_strap_depth, kMaxStrapDepth, int(L.conn.ground_ok()));
        }
    }
    {   // (3) a qubit flux line: shorted into the plane with its bump pad's ground continuous to the short, or shorted
        // into an enclosed pocket (strapped): the net is the same, the ground is not
        for (int pocket = 0; pocket < 2; ++pocket) {
            ChipLayout L; L.outline = {0, 0, 300, 400};
            std::vector<Rect> fm, fh;
            detail::cpw_path({{150, 100}, {150, 250}}, 10, 6, fm, fh, true);
            fm.push_back(rect_c(150, 82, 36, 36)); fh.push_back(rect_c(150, 82, 52, 52));
            L.islands = {{"FQ7", fm}}; L.holes = fh; L.bumps = {rect_c(150, 82, 30, 30)}; L.bump_island = {0};
            if (pocket) {                                     // a ring of etch around the shorted end, open only where the line passes
                for (const Rect& r : {Rect{100, 225, 139, 235}, Rect{161, 225, 200, 235}, Rect{100, 290, 200, 300}, Rect{100, 235, 110, 290}, Rect{190, 235, 200, 290}}) L.holes.push_back(r);
                L.straps = {{{175, 280, 185, 310}, rect_c(180, 285, 10, 10), rect_c(180, 305, 10, 10), 0, 0}, {{115, 280, 125, 310}, rect_c(120, 285, 10, 10), rect_c(120, 305, 10, 10), 0, 0}};
            }
            finish(L);
            CHECK(L.conn.flux_lines == 1 && L.conn.flux_lines_on_plane == 1 && L.conn.nets_clean(), "flux line (pocket %d): on the ground net %d/1: %s", pocket, L.conn.flux_lines_on_plane, connectivity_summary(L).c_str());
            CHECK(L.conn.qubit_flux_on_plane == 1 - pocket && L.conn.flux_returns_continuous == 1 - pocket && L.conn.ground_ok() == !pocket,
                  "flux line (pocket %d): short on plane metal %d, return continuous %d, ground ok %d", pocket, L.conn.qubit_flux_on_plane, L.conn.flux_returns_continuous, int(L.conn.ground_ok()));
        }
    }
    {   // (4) duality check: inductance per length of a CPW (10 um centre, 6 um gaps) from the capacitance of its two gap
        // apertures as conductors in vacuum, against mu0 K(k') / (4 K(k)). 3 mm long: end effects below 1 %.
        const double len = 3000, w = 10, sgap = 6;
        std::vector<Conductor> gaps(2);
        gaps[0].name = "g1"; gaps[0].rects = {{0, w / 2, len, w / 2 + sgap}};
        gaps[1].name = "g2"; gaps[1].rects = {{0, -w / 2 - sgap, len, -w / 2}};
        PanelSolution sol;
        const auto Cg = extract_capacitance(gaps, 1.0, MeshOptions{0.75, 60.0, 1.6}, &sol);
        // 1 A in the centre conductor: psi = +1/4, -1/4; flux through gap 1 = mu0 / (2 eps0) * its dual charge
        const double L_pH = dual_flux_pH_per_fF() * (Cg.at(0, 0) * 0.25 - Cg.at(0, 1) * 0.25);
        const double k = w / (w + 2 * sgap), ref_pH = 4e-7 * PI / 4 * ellipK(std::sqrt(1 - k * k)) / ellipK(k) * len * 1e-6 * 1e12;
        std::printf("  CPW inductance by duality: %.1f pH over %.0f um (%.4f pH/um), conformal mapping %.1f pH (%.4f pH/um): %+.2f%%\n", L_pH, len, L_pH / len, ref_pH, ref_pH / len, 100 * (L_pH / ref_pH - 1));
        CHECK(std::abs(L_pH / ref_pH - 1) < 0.03, "CPW inductance by duality %.2f pH vs %.2f pH", L_pH, ref_pH);
        // the panel solution reproduces the matrix: charges of row k on the panels of conductor j sum to C(k, j)
        double q01 = 0, q00 = 0;
        for (std::size_t p = 0; p < sol.panels.size(); ++p) (sol.panels[p].cond == 0 ? q00 : q01) += sol.q_fF[0][p];
        const std::vector<double> qs = sol.charge({0.25, -0.25});
        double q1 = 0; for (std::size_t p = 0; p < sol.panels.size(); ++p) if (sol.panels[p].cond == 0) q1 += qs[p];
        CHECK(std::abs(q00 / Cg.at(0, 0) - 1) < 1e-9 && std::abs(q01 / Cg.at(0, 1) - 1) < 2 * Cg.asymmetry + 1e-9 && sol.q_fF.size() == 2 &&
              std::abs(q1 - (Cg.at(0, 0) - Cg.at(0, 1)) * 0.25) < 1e-3 * std::abs(q1), "panel charges do not sum to the capacitance matrix (%.6g vs %.6g)", q00, Cg.at(0, 0));
    }
    {   // (5) flux coupling of a shorted CPW to a SQUID loop in the cut-out of a square island. The CPW comes from below
        // on the island's axis and is shorted `d` below the cut-out. With the loop on the axis the two halves of the
        // return current cancel in the loop (M ~ 0) but not in the common mode; a loop beside the axis sees a net
        // flux, which falls quickly as the short moves away; an open line carries no current at all.
        auto make = [&](double y_end, double loop_x, bool shorted) {
            ChipLayout L; L.outline = {0, 0, 600, 600};
            std::vector<Rect> fm, fh;
            detail::cpw_path({{300, 80}, {300, y_end}}, 10, 6, fm, fh, shorted);
            fm.push_back(rect_c(300, 60, 36, 36)); fh.push_back(rect_c(300, 60, 52, 52));
            L.islands = {{"Q1", {{250, 250, 350, 350}}}, {"FQ1", fm}};
            L.holes = {{220, 220, 380, 380}}; L.holes.insert(L.holes.end(), fh.begin(), fh.end());
            L.jj = {{loop_x - 10, 220, loop_x + 10, 250}}; L.jj_owner = {0};
            L.bumps = {rect_c(300, 60, 30, 30)}; L.bump_island = {1};
            finish(L);
            return L;
        };
        const MeshOptions mo{2.0, 40.0, 1.8};
        const auto Fc = flux_coupling(make(206, 300, true), 1, mo), Fn = flux_coupling(make(206, 330, true), 1, mo), Ff = flux_coupling(make(140, 330, true), 1, mo),
                   Fo = flux_coupling(make(206, 330, false), 1, mo);
        std::printf("  flux coupling by duality: loop on the axis M %.4f pH (M_cm %.3f pH); loop 30 um off the axis M %.3f pH, I(Phi0/2) %.2f mA (M_cm %.3f pH); short 66 um further away "
                    "M %.4f pH, I(Phi0/2) %.1f mA; open line M %.3g pH (%s)\n", Fc.M_pH, Fc.M_cm_pH, Fn.M_pH, Fn.I_half_mA, Fn.M_cm_pH, Ff.M_pH, Ff.I_half_mA, Fo.M_pH, Fo.note.c_str());
        CHECK(Fc.ok && Fn.ok && Ff.ok && Fo.ok && Fn.apertures >= 3, "flux coupling did not run (%s; %d apertures)", Fn.note.c_str(), Fn.apertures);
        CHECK(Fn.M_pH > 0.3 && Fn.I_half_mA < 5 && Fc.M_pH < 0.02 * Fn.M_pH && Fc.M_cm_pH > 0.2, "loop beside the axis %.3f pH, on the axis %.4f pH (M_cm %.3f pH)", Fn.M_pH, Fc.M_pH, Fc.M_cm_pH);
        CHECK(Ff.M_pH < 0.2 * Fn.M_pH && Ff.I_half_mA > kMaxHalfFluxQuantum_mA, "a short far from the loop: M %.4f pH, I(Phi0/2) %.2f mA (gate %.0f mA)", Ff.M_pH, Ff.I_half_mA, kMaxHalfFluxQuantum_mA);
        CHECK(Fo.M_pH == 0 && Fo.I_half_mA > 1e6, "an open flux line couples %.3g pH", Fo.M_pH);
    }
    {   // (6) control-line window. Q1 (a cross-free square island) with a coupler pad C1-2 beside it and its XY line:
        // (a) ending outside the qubit's cut-out behind 64 um of ground, (b) ending inside the cut-out, (c) as (a) but
        // with the line's pad cut-out merged with the coupler's (no ground between them).
        auto make = [&](int variant) {
            ChipLayout L; L.outline = {0, 0, 900, 600};
            std::vector<Rect> xm, xh;
            const double y_end = variant == 1 ? 215 : 150;
            detail::cpw_path({{300, 80}, {300, y_end}}, 10, 6, xm, xh, false);
            xm.push_back(rect_c(300, 60, 36, 36)); xh.push_back(rect_c(300, 60, 52, 52));
            L.islands = {{"Q1", {{250, 250, 350, 350}}}, {"XY1", xm}, {"C1-2", {{350, 30, 410, 90}}}};
            L.holes = {{220, 220, 380, 380}, {340, 20, 420, 100}};
            L.holes.insert(L.holes.end(), xh.begin(), xh.end());
            if (variant == 2) L.holes.push_back({326, 34, 340, 86});      // the 14 um of ground between the XY pad and the coupler pad etched away
            L.bumps = {rect_c(300, 60, 30, 30)}; L.bump_island = {1};
            finish(L);
            return L;
        };
        const MeshOptions mo{2.5, 40.0, 2.0};
        const auto Wa = control_window(make(0), 1, mo), Wb = control_window(make(1), 1, mo), Wc = control_window(make(2), 1, mo);
        CHECK(Wa.ok && Wb.ok && Wc.ok && Wa.has("C1-2") && Wa.has("XY1") && Wa.has("Q1") && Wa.flux_zq.empty(), "control windows did not run (%s)", Wa.note.c_str());
        const double qa = Wa.mutual("XY1", "Q1"), qb = Wb.mutual("XY1", "Q1"), ca = Wa.mutual("XY1", "C1-2"), cc = Wc.mutual("XY1", "C1-2");
        std::printf("  control window: C(XY, qubit) %.4f fF behind ground, %.3f fF with the end inside the cut-out; C(XY, coupler) %.4f fF behind ground, %.3f fF with the pad cut-outs merged\n", qa, qb, ca, cc);
        CHECK(qa > 0 && qb > 4 * qa && qb > 0.5, "XY end inside the cut-out: %.4f fF against %.4f fF outside", qb, qa);
        CHECK(ca > 0 && cc > 1.4 * ca && cc > 0.8, "XY pad beside a coupler pad: %.4f fF against %.4f fF behind ground", cc, ca);
        // the mediated drive: with the coupler as a passive node resonant at f_c, C_eff = C_xq + C_cq V_c with
        // V_c = C_xc / (C_c (1 - f_c^2 / f^2)): opposite in sign to the direct coupling below the coupler's frequency
        LineCell cell; cell.qubit = 1; cell.win = Wc;
        cell.nodes = {{"C1-2", 150.0, 8.4}};
        auto Cnt = [](const std::string&, const std::string&) { return 5.6; };
        const double f = 6.0, direct = Wc.mutual("XY1", "Q1"), eff = cell.C_xy_eff("Q1", f, Cnt);
        CHECK(std::abs(eff - (direct + cc * 5.6 / 150.0 * 36.0 / (36.0 - 8.4 * 8.4))) < 1e-12 && eff < direct, "mediated drive: C_eff %.5f fF, direct %.5f fF", eff, direct);
        cell.nodes = {{"C1-2", 150.0, 0.0}};                               // a floating node adds C_xn C_nq / C_n at every frequency
        CHECK(std::abs(cell.C_xy_eff("Q1", f, Cnt) - (direct + cc * 5.6 / 150.0)) < 1e-12 && std::abs(cell.C_xy_eff("Q1", f) - (direct + cc * Wc.mutual("C1-2", "Q1") / 150.0)) < 1e-12,
              "drive through a floating node");
    }
    {   // (8) the qubit cell: corner pockets and the SQUID. With the T-bar cut-outs drawn at the full gap all round
        // (T = 120, as in release 0.4) the three free inner corners are sealed pockets; with the corner relief
        // they join the ground outside through necks, the SQUID corner through the narrowest one.
        auto cell = [&](const XmonParams& q, ChipLayout& L) {
            L = ChipLayout{}; L.outline = {0, 0, 700, 700};
            std::vector<Rect> qm;
            xmon(q, 350, 350, qm, L.holes, &L.jj);
            L.islands = {{"Q1", qm}}; L.jj_owner = {0};
            finish(L);
            return detail::metal_map(L.outline, L.etch);
        };
        XmonParams q; q.L = 118; q.T = 100;
        ChipLayout L;
        auto M = cell(q, L);
        const int plane = M.at(1, 1);
        CHECK(M.at(350 - 60, 350 + 60) == plane && M.at(350 - 60, 350 - 60) == plane && M.at(350 + 60, 350 - 60) == plane && M.n_comp == 2,
              "corner pockets not joined to the plane (%d components)", M.n_comp);
        CHECK(std::abs(corner_neck(q, 1) - std::sqrt(2.0) * (118 - 12 - 50 - 20)) < 1e-9 && std::abs(corner_neck(q, 3) - corner_neck(q, 1)) < 1e-9 &&
              std::abs(corner_neck(q, 2) - (118 - 12 - 50 - q.gap - q.relief_gap)) < 1e-9 && corner_neck(q, 2) < corner_neck(q, 1), "corner necks %.1f / %.1f / %.1f um",
              corner_neck(q, 1), corner_neck(q, 2), corner_neck(q, 3));
        // the SQUID: one junction region, from the finger at the -x end of the -y T-bar to the plane, 10 um deep and 20 um wide
        const Rect j = L.jj.front();
        CHECK(L.jj.size() == 1 && L.conn.junctions_landed == 1 && L.conn.qubit_squids_on_plane == 1 && L.conn.junction_owner_errors == 0 &&
              std::abs(j.w() - (q.squid_gap + 2)) < 1e-9 && std::abs(j.h() - q.squid_w) < 1e-9 && std::abs(j.x0 - (350 - q.T / 2 - q.gap - 1)) < 1e-9 &&
              std::abs(j.cy() - (350 - q.L + q.t_w / 2)) < 1e-9, "qubit SQUID region (%.0f, %.0f) %.0f x %.0f um, landed %d", j.x0 - 350, j.y0 - 350, j.w(), j.h(), L.conn.junctions_landed);
        XmonParams old = q; old.T = 120; old.relief_gap = 0;
        M = cell(old, L);
        std::printf("  qubit cell: necks %.1f / %.1f / %.1f um (corners 1, 2, 3); without relief and T = 120 um: %d metal components, necks %.1f um\n", corner_neck(q, 1),
                    corner_neck(q, 2), corner_neck(q, 3), M.n_comp, corner_neck(old, 1));
        CHECK(M.n_comp >= 5 && M.at(350 - 60, 350 + 60) != M.at(1, 1) && M.at(350 - 60, 350 - 60) != M.at(1, 1) && corner_neck(old, 1) < 0 && corner_neck(old, 2) < 0,
              "release-0.4 qubit: pockets expected sealed (%d components)", M.n_comp);
        // on the chip the coupler runs below the axis on the +-x arms and to the right on the +-y arms
        const CouplerGeom cg;
        bool sides = true;
        for (int k = 0; k < 4; ++k) {
            const auto E = arm_elements(q, cg, 900, k);
            const Rect cb = bbox_of(E.coupler), bb = bbox_of(E.bus);
            sides = sides && (k % 2 == 0 ? cb.cy() < 0 && bb.cy() > 0 : cb.cx() > 0 && bb.cx() < 0) && (k == 0 ? cb.cx() > 0 : k == 1 ? cb.cy() > 0 : k == 2 ? cb.cx() < 0 : cb.cy() < 0);
        }
        CHECK(sides, "arm_elements: coupler and bus on the wrong side of an arm");
    }
    {   // (9) two coupled qubits assembled as a chip: launchers, stitching straps, ground bumps; flux coupling of the
        // qubit's and the coupler's line; the same with an equal-armed collar slot (the common mode is not balanced)
        ChipDesign D = alpha1_v11(3);
        const auto cp = D.couplings.front();
        std::vector<QubitDesign> two;
        for (const auto& qd : D.qubits) if (qd.id == cp.a || qd.id == cp.b) two.push_back(qd);
        D.qubits = two; D.couplings = {cp};
        CellSet cs; cs.meas.L = cs.data.L = 118; cs.meas.T = cs.data.T = 100; cs.cg.pad = 230; cs.cg.l_c = 52; cs.cg.l_b = 8;
        const auto L = assemble_chip(D, cs);
        const std::string tag = std::to_string(cp.a) + "-" + std::to_string(cp.b), qn = "Q" + std::to_string(cp.a);
        std::printf("  two-qubit chip: %s\n    %s\n", drc_summary(L).c_str(), connectivity_summary(L).c_str());
        CHECK(L.drc_clean() && L.conn.clean() && L.conn.qubit_squids_on_plane == 2 && L.conn.qubit_flux_on_plane == 2 && L.conn.flux_returns_continuous == 2 &&
              L.conn.coupler_squids_on_plane == 1 && L.conn.straps == 9 && L.conn.straps_landed == 9 && L.ground_bumps.size() == 6 && L.island("XY" + std::to_string(cp.a)) < 0,
              "two-qubit chip: %zu problems, %zu DRC markers, %d straps, %zu ground bumps", L.problems.size(), L.drc.size(), L.conn.straps, L.ground_bumps.size());
        const MeshOptions mo{2.5, 50.0, 2.0};
        const auto Fq = flux_coupling(L, cp.a, mo), Fc = flux_coupling(L, "C" + tag, "FC" + tag, mo);
        CellSet eq = cs; eq.meas.slot_up = eq.meas.slot_down = eq.data.slot_up = eq.data.slot_down = 18;
        const auto Fe = flux_coupling(assemble_chip(D, eq), cp.a, mo);
        std::printf("  launcher: qubit M %.2f pH, I(Phi0/2) %.2f mA, M_cm %.3f pH (%.3f pH with an 18/18 um slot); coupler M %.2f pH, I(Phi0/2) %.2f mA, M_cm %.3f pH\n", Fq.M_pH,
                    Fq.I_half_mA, Fq.M_cm_pH, Fe.M_cm_pH, Fc.M_pH, Fc.I_half_mA, Fc.M_cm_pH);
        CHECK(Fq.ok && Fc.ok && Fe.ok, "flux coupling did not run: %s %s", Fq.note.c_str(), Fc.note.c_str());
        CHECK(Fq.M_pH > 1.0 && Fq.I_half_mA < kMaxHalfFluxQuantum_mA && Fc.M_pH > 0.8 && Fc.I_half_mA < kMaxHalfFluxQuantum_mA, "launcher: qubit %.2f pH, coupler %.2f pH", Fq.M_pH, Fc.M_pH);
        CHECK(Fe.M_cm_pH > 2 * Fq.M_cm_pH && Fe.M_cm_pH > 0.5, "collar slot balance: M_cm %.3f pH with %g/%g um arms, %.3f pH with equal arms", Fq.M_cm_pH, cs.meas.slot_up, cs.meas.slot_down, Fe.M_cm_pH);
    }
    {   // (7) formulas. T1 of a coupling C into a matched line: C_sigma / (omega^2 C^2 Z0); the shorted line's weight
        // sin(beta z) makes the effective capacitance proportional to f, so that limit goes as f^-4, not f^-2.
        const double t = line_T1_us(95.0, 0.03, 6.0);
        const double w = 2 * PI * 6e9, ref = 95e-15 / (w * w * 0.03e-15 * 0.03e-15 * 50) * 1e6;
        CHECK(std::abs(t / ref - 1) < 1e-12 && std::abs(line_T1_us(95.0, -0.03, 6.0) - t) < 1e-9 * t && line_T1_us(95.0, 0.0, 6.0) > 1e200, "line_T1_us %.4g us vs %.4g us", t, ref);
        LineCell c; c.win.flux_zq = {{100.0, 1.0}, {200.0, 0.5}};          // 1 fF at 100 um and 0.5 fF at 200 um from the short
        const double beta = 2 * PI * 5.5e9 * std::sqrt(0.5 * (11.45 + 1)) / 299792458.0 * 1e-6;
        CHECK(std::abs(c.C_flux_eff(5.5) - (std::sin(beta * 100) + 0.5 * std::sin(beta * 200))) < 1e-12, "C_flux_eff %.6f", c.C_flux_eff(5.5));
        const double r_c = c.C_flux_eff(6.0) / c.C_flux_eff(5.5), r_t = line_T1_us(95, c.C_flux_eff(5.5), 5.5) / line_T1_us(95, c.C_flux_eff(6.0), 6.0);
        CHECK(std::abs(r_c - 6.0 / 5.5) < 1e-3 && std::abs(r_t - std::pow(6.0 / 5.5, 4)) < 5e-3, "flux coupling scales as f (%.4f), its T1 limit as f^-4 (%.4f)", r_c, r_t);
        // drive power: 0.030 fF into a 95 fF qubit with n01 = 1.33 needs about -74 dBm for a 25 MHz Rabi rate
        const double p = drive_power_dBm(0.030, 95.0, 1.33);
        CHECK(std::abs(p + 74.2) < 0.3 && std::abs(drive_power_dBm(0.015, 95.0, 1.33) - p - 6.02) < 0.01, "drive power %.2f dBm", p);
        // drive through the flux line: the port's rate 1/T1 = (omega M)^2 / (Z0 L_J) fixes the power, P = Omega^2 hbar omega / (4 / T1)
        {
            const double EJd = 20.0, ECd = 0.2, Md = 0.25, fd = std::sqrt(8 * EJd * ECd);            // f01 = sqrt(8 E_J E_C) in this check
            const double rate = 1e6 / flux_cm_T1_us(Md, EJd, fd), Om = 2 * PI * kRabiMHz * 1e6;
            const double p_ref = 10 * std::log10(Om * Om * phys::h * fd * 1e9 / (4 * rate) / 1e-3), p_flux = flux_drive_power_dBm(Md, EJd, ECd);
            CHECK(std::abs(p_flux - p_ref) < 0.01 && std::abs(flux_drive_power_dBm(Md / 2, EJd, ECd) - p_flux - 6.02) < 0.01 && p_flux < -60 && p_flux > -80,
                  "flux-line drive power %.2f dBm, from the port's rate %.2f dBm", p_flux, p_ref);
        }
        // loss through the SQUID loop: none for a symmetric SQUID, d^2 M^2 otherwise; about 33 ms for 2 pH at d = 0.1
        const double tl = flux_loop_T1_us(2.0, 20.0, 0.2, 0.1, 0.0, 5.5);
        CHECK(flux_loop_T1_us(2.0, 20.0, 0.2, 0.0, 0.0, 5.5) > 1e200 && std::abs(flux_loop_T1_us(4.0, 20.0, 0.2, 0.2, 0.0, 5.5) / tl - 1.0 / 16) < 1e-9 && tl > 2e4 && tl < 5e4,
              "flux-loop T1 limit %.3g us", tl);
        // common mode: 1 pH in series with a 9 nH junction at 6 GHz against 50 ohm: Z0 L_J / (omega M)^2 = 0.32 ms
        const double EJ = std::pow(phys::Phi0 / (2 * PI), 2) / 9e-9 / phys::h / 1e9, tc = flux_cm_T1_us(1.0, EJ, 6.0);
        CHECK(std::abs(tc - 50 * 9e-9 / std::pow(w * 1e-12, 2) * 1e6) < 1e-6 * tc && std::abs(tc - 317) < 3 && flux_cm_T1_us(0.0, EJ, 6.0) > 1e200, "common-mode T1 limit %.1f us", tc);
    }
}

// ---------------------------------------------------------------- layout_presets: the committed chips are mask-clean and on target
#ifndef QPUFORGE_SOURCE_DIR
#define QPUFORGE_SOURCE_DIR "."
#endif
static void test_layout_presets() {
    std::puts("[layout_presets] AF-17, BF-49, GF-97 from the committed cell files: placement clean, mask DRC, connectivity and ground gate, straps, control-line gate, extraction on target, re-nullable");
    using namespace layout;
    struct Chip { const char* name; ChipDesign D; const char* cells; };
    const std::vector<Chip> chips{{"AF-17", alpha1_v11(3), "designs/alpha1-v11/layout_cells.json"},
                                  {"BF-49", beta1(5), "designs/beta1/layout_cells.json"},
                                  {"GF-97", gamma1(7), "designs/beta1/layout_cells.json"}};
    std::map<std::string, LayoutMetrics> metrics;                    // the Gamma-1 chip reuses the Beta-1 cells
    for (const auto& c : chips) {
        std::ifstream f(std::string(QPUFORGE_SOURCE_DIR) + "/" + c.cells);
        CHECK(bool(f), "%s: cannot open %s", c.name, c.cells);
        if (!f) continue;
        std::stringstream ss; ss << f.rdbuf();
        const CellSet C = cells_from_json(Json::parse(ss.str()));
        const auto L = assemble_chip(c.D, C);
        if (!metrics.count(c.cells)) metrics[c.cells] = evaluate_cells(C, MeshOptions{1.5, 50.0, 1.6});
        const auto& m = metrics[c.cells];
        const auto R = analyze(c.D);
        const auto T = targets_from_design(c.D, R);
        const auto P = post_layout(c.D, R, m);
        auto rel = [](double a, double t) { return std::abs(a / t - 1); };
        const double ec = std::max({rel(m.EC_m, T.EC_m), rel(m.EC_d, T.EC_d), rel(m.EC_c, T.EC_c)});
        const double g = std::max({rel(m.G_mc, T.G_mc), rel(m.G_cd, T.G_cd), rel(m.G_md, T.G_md)});
        const double cr = rel(m.C_qr_fF, T.C_qr_fF);
        std::printf("  %s: %s; %d sliver fills; E_C within %.2f%%, couplings within %.2f%%, readout C within %.2f%%, re-nulled gap %.1e MHz\n    %s\n", c.name,
                    drc_summary(L).c_str(), L.slivers_filled, 100 * ec, 100 * g, 100 * cr, P.park_gap_MHz, connectivity_summary(L).c_str());
        CHECK(L.problems.empty(), "%s placement: %s", c.name, L.problems.empty() ? "" : L.problems.front().c_str());
        CHECK(L.drc_count("width") == 0 && L.drc_count("space") == 0, "%s mask DRC: %d width + %d space violations", c.name, L.drc_count("width"), L.drc_count("space"));
        // connectivity gate: nets (every conductor its own net, every shorted line and ground fragment on the ground
        // net, every SQUID from its own island to ground, one SQUID per qubit and coupler) and ground quality (qubit
        // SQUID grounds and flux-line shorts on plane metal with no strap in the path, no fragment on a single strap)
        CHECK(L.conn.nets_clean(), "%s nets: %s", c.name, connectivity_summary(L).c_str());
        CHECK(L.conn.ground_ok(), "%s ground: qubit SQUID grounds on plane metal %d/%d, flux shorts %d/%d (returns continuous %d), fragments on one strap %d, depth %d",
              c.name, L.conn.qubit_squids_on_plane, L.conn.qubit_squids, L.conn.qubit_flux_on_plane, L.conn.flux_lines, L.conn.flux_returns_continuous,
              L.conn.fragments_on_one_strap, L.conn.max_strap_depth);
        CHECK(L.conn.junctions == int(c.D.qubits.size() + c.D.couplings.size()) && L.conn.qubit_squids == int(c.D.qubits.size()) &&
              L.conn.coupler_squids == int(c.D.couplings.size()) && L.conn.flux_lines == int(c.D.qubits.size()) && L.conn.junction_owner_errors == 0,
              "%s: %d junctions (%d qubit, %d coupler, %d owner errors), %d flux lines", c.name, L.conn.junctions, L.conn.qubit_squids, L.conn.coupler_squids,
              L.conn.junction_owner_errors, L.conn.flux_lines);
        // straps: feet on ground metal only, spans clear of every qubit island and junction region
        int foot_on_island = 0, span_over_qubit = 0;
        {
            std::vector<std::pair<Rect, char>> isl;
            for (const auto& [n, rs] : L.islands) for (const auto& r : rs) isl.push_back({r, n[0]});
            std::sort(isl.begin(), isl.end(), [](const auto& a, const auto& b) { return a.first.x0 < b.first.x0; });
            double wmax = 0; for (const auto& i : isl) wmax = std::max(wmax, i.first.w());
            auto hits = [&](const Rect& q, bool qubits_only) {
                auto it = std::lower_bound(isl.begin(), isl.end(), q.x0 - wmax - 1.0, [](const auto& a, double v) { return a.first.x0 < v; });
                for (; it != isl.end() && it->first.x0 < q.x1; ++it) if (it->first.intersects(q) && (!qubits_only || it->second == 'Q')) return true;
                return false;
            };
            for (const auto& st : L.straps) {
                foot_on_island += hits(st.foot_a, false) || hits(st.foot_b, false);
                bool over = hits(st.span, true);
                for (const auto& j : L.jj) over = over || st.span.intersects(j);
                span_over_qubit += over;
            }
        }
        CHECK(foot_on_island == 0 && span_over_qubit == 0, "%s straps: %d with a foot on island metal, %d over a qubit or a junction region", c.name, foot_on_island, span_over_qubit);
        // control-line gate (window on the coarse mesh): drive and flux-line T1 limits, drive crosstalk, drive power, I(Phi0/2)
        auto LB = line_budget(c.D, R, L, m, P.park_fc_GHz, MeshOptions{3.0, 50.0, 2.0});
        line_budget_crosstalk(LB, c.D, R, L, m);
        std::printf("    %s\n", line_budget_summary(LB).c_str());
        CHECK(LB.evaluated && LB.crosstalk_evaluated, "%s control lines not evaluated: %s", c.name, LB.note.c_str());
        CHECK(LB.drive_ok(), "%s drive: T1 limit %.1fx (budget %.0fx), crosstalk %.1f%% (%s), %.1f dBm", c.name, std::min(LB.min_ratio_drive, 9999.0), kLineT1Margin,
              100 * std::min(LB.max_drive_crosstalk, 99.0), LB.worst_crosstalk.c_str(), LB.drive_power_dBm);
        CHECK(LB.flux_ok(), "%s flux lines: I(Phi0/2) %.3g mA (limit %.0f), T1 limits capacitive %.1fx, loop %.1fx, common mode %.1fx, coupler %.1fx", c.name,
              std::min(LB.max_I_half_mA, 1e6), kMaxHalfFluxQuantum_mA, std::min(LB.min_ratio_flux_cap, 9999.0), std::min(LB.min_ratio_flux_loop, 9999.0),
              std::min(LB.min_ratio_flux_cm, 9999.0), std::min(LB.ratio_coupler, 9999.0));
        CHECK(ec <= 0.03 && g <= 0.06 && cr <= 0.06, "%s extraction off target (E_C %.3f, G %.3f, C_qr %.3f)", c.name, ec, g, cr);
        CHECK(std::abs(P.park_gap_MHz) < 0.05, "%s parked exchange not re-nullable (%.3f MHz)", c.name, P.park_gap_MHz);
    }
}

// ---------------------------------------------------------------- materials: forge.materials.v0 sidecar reader (atomforge -> qpuforge, docs/DATA_CONTRACT.md)
#ifndef QPUFORGE_TEST_DATA
#define QPUFORGE_TEST_DATA "tests/data"   // CMake passes the absolute path; this default works from the source root
#endif
static std::string test_data(const char* name) { return std::string(QPUFORGE_TEST_DATA) + "/" + name; }
// a minimal sidecar around the members of `quantities`
static std::string sidecar_text(const std::string& quantities, const std::string& schema = "\"forge.materials.v0\"") {
    return "{\"schema\": " + schema + ", \"quantities\": {" + quantities + "}}";
}
// reading `text` must throw std::runtime_error with `needle` in the message
static bool sidecar_rejected(const std::string& text, const char* needle) {
    try {
        materials_from_json(Json::parse(text));
    } catch (const std::runtime_error& e) {
        if (std::string(e.what()).find(needle) != std::string::npos) return true;
        std::printf("  rejected, but with another message: %s\n", e.what());
    }
    return false;
}
static void test_materials() {
    std::puts("[materials] forge.materials.v0 sidecar: atomforge's placeholder file changes nothing, values reach analyze(), loss tangent on request only, malformed files throw");
    const ChipDesign A0 = alpha1_v11(3);
    const std::string A0_json = json_text(design_to_json(A0));
    // (a) the sidecar exactly as atomforge's writer emits it: materials_contract_json() at atomforge ee734f9, for the
    //     synthetic Cu-Au search and fixed options of its tests/test_materials_contract.cpp. All three quantities are
    //     null placeholders.
    const auto P = load_materials(test_data("materials-contract-CuAu.json"));
    CHECK(P.schema == "forge.materials.v0" && P.produced_by.tool == "atomforge" && P.produced_by.potential_or_method == "Zhou04-EAM" &&
          P.produced_by.commit == "0123456789abcdef0123456789abcdef01234567" && P.produced_by.timestamp == "2026-11-26T12:00:00+00:00" &&
          P.system == std::vector<std::string>({"Cu", "Au"}), "atomforge sidecar metadata");
    bool placeholders = P.quantities.size() == 3;
    for (const auto& q : P.quantities) placeholders = placeholders && q.placeholder && !q.value && q.status == "placeholder";
    CHECK(placeholders && P.find("Q_dielectric") && P.find("loss_tangent") && P.find("junction_gap_ueV") &&
          P.find("Q_dielectric")->unit == "dimensionless" && P.find("junction_gap_ueV")->unit == "ueV" &&
          P.find("Q_dielectric")->maps_to == "ChipDesign.materials.Q_dielectric" && P.find("loss_tangent")->maps_to.empty(),
          "atomforge sidecar: three null placeholders with their units and maps_to");
    for (bool lt : {false, true}) {
        ChipDesign A = A0;
        const auto rep = apply_materials(P, A, MaterialsApplyOptions{lt});
        CHECK(json_text(design_to_json(A)) == A0_json && A.Q_dielectric == A0.Q_dielectric && A.gap_ueV == A0.gap_ueV &&
              rep.applied() == 0 && rep.skipped() == 3, "placeholder sidecar must leave alpha1_v11(3) unchanged (loss-tangent opt-in %d)", int(lt));
    }
    {   // the qpuforge defaults quoted beside the placeholders (3e6, 180 ueV) equal Alpha-1's own values: check on other values too
        ChipDesign B = A0; B.Q_dielectric = 1.5e6; B.gap_ueV = 175;
        apply_materials(P, B);
        CHECK(B.Q_dielectric == 1.5e6 && B.gap_ueV == 175, "placeholder sidecar applied a reference default (Q %.6g, gap %.6g)", B.Q_dielectric, B.gap_ueV);
    }
    // (b) a sidecar with values (Q 1e6, gap 190 ueV) reaches the device model
    const auto V = load_materials(test_data("materials-values.json"));
    ChipDesign C = A0;
    const auto rv = apply_materials(V, C);
    CHECK(C.Q_dielectric == 1e6 && C.gap_ueV == 190 && rv.applied() == 2 && rv.skipped() == 1 && rv.lines.size() == 3 &&
          rv.lines[0].applied && rv.lines[0].before == 3e6 && rv.lines[0].after == 1e6 && !rv.lines[1].applied &&
          rv.lines[2].field == "materials.junction_gap_ueV" && rv.lines[2].after == 190, "values sidecar: fields and report");
    {
        ChipDesign C2 = C; C2.Q_dielectric = A0.Q_dielectric; C2.gap_ueV = A0.gap_ueV;
        CHECK(json_text(design_to_json(C2)) == A0_json, "a sidecar may change Q_dielectric and gap_ueV only");
    }
    const auto R0 = analyze(A0), R1 = analyze(C);
    const auto& q0 = R0.q[0]; const auto& q1 = R1.q[0];
    const double T1_formula = 1e6 / (2 * PI * q1.f01 * 1e9 / C.Q_dielectric);     // G_diel = 2 pi f / Q
    std::printf("  Q%d, Q 3e6 -> 1e6 and gap 180 -> 190 ueV: f01 %.4f -> %.4f GHz, T1_diel %.1f -> %.1f us (2 pi f / Q gives %.1f us), T1 %.1f -> %.1f us\n",
                q1.id, q0.f01, q1.f01, q0.T1_diel_us, q1.T1_diel_us, T1_formula, q0.T1_us, q1.T1_us);
    CHECK(std::abs(q1.T1_diel_us / T1_formula - 1) < 1e-12, "T1_diel %.6f us vs 2 pi f / Q %.6f us", q1.T1_diel_us, T1_formula);
    CHECK(std::abs(q1.EJ / q0.EJ - 190.0 / 180.0) < 1e-12, "E_J must scale with the gap (Ambegaokar-Baratoff): ratio %.15f", q1.EJ / q0.EJ);
    CHECK(std::abs(1 / q1.T1_us - (1 / q1.T1_diel_us + 1 / q1.T1_purcell_us)) < 1e-12 / q1.T1_us && q1.T1_us < q0.T1_us,
          "T1 must combine the new dielectric rate with Purcell decay (%.3f vs %.3f us)", q1.T1_us, q0.T1_us);
    {   // Q alone: the frequencies stay put, so every T1_diel scales exactly with Q
        ChipDesign E = A0;
        apply_materials(materials_from_json(Json::parse(sidecar_text(R"("Q_dielectric": {"value": 1.2e7, "unit": "dimensionless"})"))), E);
        const auto RE = analyze(E);
        double worst = 0; bool same_f = RE.q.size() == R0.q.size();
        for (std::size_t i = 0; same_f && i < RE.q.size(); ++i) {
            same_f = RE.q[i].f01 == R0.q[i].f01;
            worst = std::max(worst, std::abs(RE.q[i].T1_diel_us / R0.q[i].T1_diel_us / 4.0 - 1));
        }
        CHECK(same_f && worst < 1e-12, "Q 3e6 -> 1.2e7: f01 unchanged %d, T1_diel ratio off 4 by %.2e", int(same_f), worst);
    }
    // (c) loss tangent: Q_dielectric = 1/tan(delta) on request only, and only while Q_dielectric is null
    const auto L = load_materials(test_data("materials-loss-tangent.json"));     // Q_dielectric null, loss_tangent 2.5e-7
    {
        ChipDesign F = A0;
        const auto rf = apply_materials(L, F);
        CHECK(F.Q_dielectric == A0.Q_dielectric && rf.applied() == 0 && rf.lines.size() == 3 && rf.lines[1].quantity == "loss_tangent" && !rf.lines[1].applied,
              "loss_tangent must not be applied without the opt-in");
        ChipDesign G = A0;
        const auto rg = apply_materials(L, G, MaterialsApplyOptions{true});
        const auto RG = analyze(G);
        const double T1_tan = 1e6 / (2 * PI * RG.q[0].f01 * 1e9 * 2.5e-7);            // G_diel = 2 pi f tan(delta)
        std::printf("  loss_tangent 2.5e-7 with the opt-in: Q_dielectric %.6g, Q%d T1_diel %.1f us (2 pi f tan(delta) gives %.1f us)\n",
                    G.Q_dielectric, RG.q[0].id, RG.q[0].T1_diel_us, T1_tan);
        CHECK(G.Q_dielectric == 1.0 / 2.5e-7 && G.gap_ueV == A0.gap_ueV && rg.applied() == 1 && rg.lines[1].applied && rg.lines[1].after == G.Q_dielectric &&
              rg.lines[1].note.find("1/loss_tangent") != std::string::npos, "opt-in: Q_dielectric = 1/loss_tangent, recorded in the report");
        CHECK(std::abs(RG.q[0].T1_diel_us / T1_tan - 1) < 1e-12, "T1_diel from the loss tangent %.6f vs %.6f us", RG.q[0].T1_diel_us, T1_tan);
        ChipDesign H = A0;    // a Q_dielectric with a value wins, opt-in or not
        const auto rh = apply_materials(materials_from_json(Json::parse(sidecar_text(
            R"("loss_tangent": {"value": 1e-6, "unit": "dimensionless"}, "Q_dielectric": {"value": 2e6, "unit": "dimensionless"})"))), H, MaterialsApplyOptions{true});
        CHECK(H.Q_dielectric == 2e6 && rh.applied() == 1 && !rh.lines[0].applied, "a given Q_dielectric must win over loss_tangent (Q %.6g)", H.Q_dielectric);
    }
    // (d) malformed sidecars are rejected, each with its own message
    const std::string q_ok = R"("Q_dielectric": {"value": 1e6, "unit": "dimensionless"})";
    CHECK(sidecar_rejected(sidecar_text(q_ok, "\"forge.materials.v1\""), "schema must be \"forge.materials.v0\", got \"forge.materials.v1\"") &&
          sidecar_rejected(R"({"quantities": {}})", "no schema member") && sidecar_rejected(sidecar_text(q_ok, "0"), "non-string") &&
          sidecar_rejected("[]", "JSON object"), "a wrong or missing schema must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": 1e6, "unit": "ppm"})"), "unit must be \"dimensionless\"") &&
          sidecar_rejected(sidecar_text(R"("junction_gap_ueV": {"value": 0.19, "unit": "meV"})"), "unit must be \"ueV\"") &&
          sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": 1e6})"), "unit is missing"), "a wrong or missing unit must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": 0, "unit": "dimensionless"})"), "positive finite") &&
          sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": -3e6, "unit": "dimensionless"})"), "positive finite"), "a non-positive Q must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("junction_gap_ueV": {"value": -180, "unit": "ueV"})"), "positive finite") &&
          sidecar_rejected(sidecar_text(R"("loss_tangent": {"value": 0, "unit": "dimensionless"})"), "positive finite") &&
          sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": "3e6", "unit": "dimensionless"})"), "number or null") &&
          sidecar_rejected(sidecar_text(R"("Q_dielectric": {"unit": "dimensionless"})"), "value is missing"), "a bad gap, loss tangent or value type must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("participation_ratio": {"value": 0.1, "unit": "dimensionless"})"), "unknown quantity 'participation_ratio'"),
          "an unknown quantity must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("Rn_spread_pct": {"value": 0.5, "unit": "%"})"), "not part of forge.materials.v0") &&
          sidecar_rejected(sidecar_text(R"("T_eff_mK": {"value": 35, "unit": "mK"})"), "not part of forge.materials.v0"),
          "Rn_spread_pct and the environment fields are not contract quantities");
    CHECK(sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": 3e6, "unit": "dimensionless", "status": "placeholder"})"), "placeholder must have value null") &&
          sidecar_rejected(R"({"schema": "forge.materials.v0", "placeholders": ["quantities.Q_dielectric"], "quantities": {"Q_dielectric": {"value": 3e6, "unit": "dimensionless"}}})",
                           "placeholder must have value null"), "a placeholder that carries a number must throw");
    CHECK(sidecar_rejected(sidecar_text(R"("Q_dielectric": {"value": 1e6, "unit": "dimensionless", "maps_to": "ChipDesign.materials.junction_gap_ueV"})"), "maps_to") &&
          sidecar_rejected(sidecar_text(R"("loss_tangent": {"value": 1e-6, "unit": "dimensionless", "maps_to": "ChipDesign.materials.Q_dielectric"})"), "maps_to"),
          "a maps_to naming another field must throw (loss_tangent is converted, never copied)");
    CHECK(sidecar_rejected(R"({"schema": "forge.materials.v0"})", "quantities must be an object") &&
          sidecar_rejected(sidecar_text(R"("Q_dielectric": 3e6)"), "must be an object"), "a missing quantities block or a bare number must throw");
    {   // file errors name the file; a hand-built sidecar with a bad value throws before the design is touched
        bool missing = false, named = false, threw = false;
        try { load_materials(temp_path("no_such_sidecar.json")); } catch (const std::runtime_error& e) { missing = std::string(e.what()).find("cannot open") != std::string::npos; }
        const std::string chip = temp_path("chip_as_sidecar.json");
        save_design(A0, chip);
        try { load_materials(chip); } catch (const std::runtime_error& e) { named = std::string(e.what()).find(chip + ": materials sidecar: schema") == 0; }
        std::remove(chip.c_str());
        MaterialsSidecar bad = V;
        for (auto& q : bad.quantities) if (q.name == "junction_gap_ueV") q.value = -1.0;
        ChipDesign K = A0;
        try { apply_materials(bad, K); } catch (const std::runtime_error&) { threw = true; }
        CHECK(missing && named && threw && json_text(design_to_json(K)) == A0_json,
              "file errors (missing %d, named %d); a failed apply must leave the design unchanged (threw %d)", int(missing), int(named), int(threw));
    }
}

// ---------------------------------------------------------------- lattice surgery (surgery.hpp)
// Single-fault distance check with uniform edge weights. The graph is built from every mechanism that flips one or two
// decoded detectors; parallel edges (same detectors, different observable flag) keep the flag with the larger total
// probability, then every weight is set equal. Every mechanism, hyperedges included, is then decoded. Mechanisms the
// matching gets wrong are counted with their total probability and split into coherent-ZZ twirl terms and the rest.
struct SingleFaults { std::size_t faults{0};   // mechanisms per observable's DEM
                      int uncorrected{0}, uncorrected_other{0}, undetectable{0}; double uncorrected_p{0}; };
static SingleFaults surgery_single_faults(SurgerySetup& M) {
    const auto fk = surgery_fault_signatures(M);
    SingleFaults r; DecodeScratch ws;
    r.faults = fk.empty() ? 0 : fk[0].size();
    for (std::size_t o = 0; o < fk.size(); ++o) {
        std::map<std::pair<int, int>, std::pair<double, double>> acc;   // (a, b) -> total p with observable 0 / 1
        for (const auto& f : fk[o]) {
            if (f.dets.empty() || f.dets.size() > 2) continue;
            auto& v = acc[{f.dets[0], f.dets.size() == 2 ? f.dets[1] : -1}];
            (f.obs ? v.second : v.first) += f.p;
        }
        DetectorErrorModel U; U.n_det = M.dem[o].n_det;
        for (const auto& [k, v] : acc) U.edges.push_back({k.first, k.second, v.second > v.first, 0.01});
        MatchingDecoder dec(U);
        for (const auto& f : fk[o]) {
            if (f.dets.empty()) { r.undetectable += f.obs; continue; }
            if (dec.decode(f.dets, ws) == int(f.obs)) continue;
            ++r.uncorrected; r.uncorrected_p += f.p;
            if (f.what != "ZZ") ++r.uncorrected_other;
        }
    }
    return r;
}

static bool surgery_failed(const SurgeryCheck& c) { return c.fired || c.wrong_observables || c.wrong_reference || c.random_reference; }

// leakage-on variants of a noise model, for the comparisons of the two frame runners
static NoiseModel with_leakage(const ChipDesign& D, const ChipReport& R, NoiseModel nm, bool removing_reset, bool dqlr) {
    nm.leakage = true; nm.leak_1q = 1e-3; nm.leak_transport = 0.2; nm.reset_removes_leakage = removing_reset;
    nm.dqlr = dqlr; nm.dqlr_eff = 0.98; nm.dqlr_err = 5e-4;
    for (const auto& p : R.pairs) if (p.designed) {
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
        nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
    }
    return nm;
}

static void test_surgery() {
    std::puts("[surgery] multi-patch geometry, ZZ / XX parity measurement: tableau check, memory equivalence, DEM, noise");
    // geometry: two d x d patches and their merged region
    for (int d : {3, 5}) {
        const auto E = surgery_parity(d, 'Z', '0', '0');
        const auto& L = E.sched.lat;
        int nd = 0; for (const auto& st : L.sites) nd += st.data;
        const auto merged = region_faces(E.sched.phases[1].regions[0], L.at);
        int seam_z = 0;
        std::set<int> sep; for (const auto& f : E.phase_faces[0]) sep.insert(f.qubit);
        for (const auto& f : merged) if (!sep.count(f.qubit)) { seam_z += f.type == 'Z'; CHECK(f.type == 'Z', "a new face of a ZZ merge is X type"); }
        CHECK(nd == d * (2 * d + 1), "ZZ lattice d=%d: %d data qubits", d, nd);
        CHECK(int(merged.size()) == d * (2 * d + 1) - 1, "merged d x (2d+1) region has %zu faces", merged.size());
        CHECK(int(E.phase_faces[0].size()) == 2 * (d * d - 1), "two separate patches: %zu faces", E.phase_faces[0].size());
        CHECK(seam_z == d + 1, "ZZ seam has %d new Z faces (expected d+1)", seam_z);
        const auto X = surgery_parity(d, 'X', '+', '+');
        std::set<int> sepx; for (const auto& f : X.phase_faces[0]) sepx.insert(f.qubit);
        int seam_x = 0;
        for (const auto& f : X.phase_faces[1]) if (!sepx.count(f.qubit)) { seam_x += f.type == 'X'; CHECK(f.type == 'X', "a new face of an XX merge is Z type"); }
        CHECK(seam_x == d + 1, "XX seam has %d new X faces (expected d+1)", seam_x);
    }
    // noiseless tableau check: logical inputs, both parities, d = 3 and 5, with and without a separate round before
    // the merge; one seed per experiment
    int runs = 0, bad = 0;
    std::uint64_t seed = 1000;
    SurgeryRounds pre0; pre0.pre = 0;
    for (int d : {3, 5})
        for (char ty : {'Z', 'X'})
            for (const char* in : {"00", "01", "10", "11", "++", "+-", "-+", "--"})
                for (const auto& rr : {SurgeryRounds{}, pre0}) {
                    if (d == 5 && rr.pre == 0) continue;
                    const auto E = surgery_parity(d, ty, in[0], in[1], rr);
                    const auto c = check_surgery_noiseless(E, d == 3 ? 4 : 2, ++seed);
                    ++runs;
                    if (surgery_failed(c)) {
                        ++bad; std::printf("    %s (pre %d): fired %d, wrong observables %d, reference %d/%d\n", E.sched.name.c_str(), rr.pre, c.fired, c.wrong_observables, c.wrong_reference, c.random_reference);
                    }
                }
    std::printf("  noiseless tableau check: %d parity experiments (8 inputs x ZZ, XX x d=3,5; d=3 also merging from the first round), %d failed\n", runs, bad);
    CHECK(bad == 0, "noiseless surgery: %d experiments with a fired detector or a wrong observable", bad);
    {
        const auto E = surgery_parity(3, 'Z', '0', '1');
        const auto c = check_surgery_noiseless(E, 2, 5);
        CHECK(c.observable_values.size() == 3 && c.observable_values[0] == 1 && c.observable_values[1] == 0 && c.observable_values[2] == 1, "ZZ on |01>: outcome m_ZZ = 1");
        const auto Xm = surgery_parity(3, 'X', '+', '-');
        const auto cx = check_surgery_noiseless(Xm, 2, 5);
        CHECK(!cx.observable_values.empty() && cx.observable_values[0] == 1, "XX on |+->: outcome m_XX = 1");
    }
    // negative controls: the check must catch a wrong observable and an inserted error
    {
        auto E = surgery_parity(3, 'Z', '+', '+');            // X_A X_B through a ZZ merge needs the gap qubit at the split
        E.observables[0].bits.pop_back();
        int wrong = 0; for (int s = 0; s < 8; ++s) wrong += check_surgery_noiseless(E, 1, std::uint64_t(100 + s)).wrong_observables;
        CHECK(wrong > 0, "dropping the split byproduct from X_A X_B went unnoticed");
        auto F = surgery_parity(3, 'Z', '0', '0');
        auto& gs = F.circuit.gates;
        std::size_t i = gs.size() / 2;
        while (gs[i].op != Op::MEASURE) ++i;
        Gate g; g.op = Op::U; g.q = {F.sched.lat.data_qubit(1, 1)}; g.u = gates::X(); g.name = "x";
        gs.insert(gs.begin() + std::ptrdiff_t(i), g);
        CHECK(check_surgery_noiseless(F, 1, 1).fired > 0, "an X error in the merged phase fired no detector");
    }
    // single patch through the surgery engine = memory_experiment on the same chip as beta1(3). The two frame runners
    // must agree shot for shot: the error counts after every prefix of 1..120 shots, and the totals over 6000 shots,
    // with plain noise and with leakage (CZ leakage 2e-3, 1q 1e-3, transport 0.2; with a plain or a leakage-removing
    // reset; with DQLR)
    {
        const auto E = surgery_memory(3, 3, 'Z');
        const auto D = surgery_chip(E.sched.lat, "beta1"), B = beta1(3);
        bool same = D.qubits.size() == B.qubits.size() && D.couplings.size() == B.couplings.size();
        for (std::size_t k = 0; same && k < D.qubits.size(); ++k)
            same = D.qubits[k].C_fF == B.qubits[k].C_fF && D.qubits[k].Rn_kOhm == B.qubits[k].Rn_kOhm && D.qubits[k].ro.f_GHz == B.qubits[k].ro.f_GHz &&
                   D.qubits[k].ro.g_MHz == B.qubits[k].ro.g_MHz && D.qubits[k].x_um == B.qubits[k].x_um && D.qubits[k].y_um == B.qubits[k].y_um &&
                   D.qubits[k].feedline == B.qubits[k].feedline && D.qubits[k].role == B.qubits[k].role;
        for (std::size_t k = 0; same && k < D.couplings.size(); ++k) same = D.couplings[k].a == B.couplings[k].a && D.couplings[k].b == B.couplings[k].b && D.couplings[k].C_fF == B.couplings[k].C_fF;
        CHECK(same, "surgery_chip of one d=3 patch differs from beta1(3)");
        const auto R = analyze(B); const auto nm0 = noise_model(B, R);
        struct Variant { const char* name; NoiseModel nm; };
        const std::vector<Variant> variants = {{"plain", nm0}, {"leakage, plain reset", with_leakage(B, R, nm0, false, false)},
                                               {"leakage, leakage-removing reset", with_leakage(B, R, nm0, true, false)},
                                               {"leakage, DQLR", with_leakage(B, R, nm0, true, true)}};
        for (const auto& v : variants)
            for (char basis : {'Z', 'X'}) {
                if (std::string(v.name) == "leakage, DQLR" && basis == 'X') continue;
                const auto Em = basis == 'Z' ? E : surgery_memory(3, 3, 'X');
                const auto M = prepare_memory(B, R, v.nm, 3, basis);
                const auto S = prepare_surgery(B, R, v.nm, Em, basis);
                const auto a = run_memory_frame(M, 6000, 21);
                const auto b = run_surgery_frame(S, 6000, 21);
                int prefix_diff = 0;
                for (long n = 1; n <= 120; ++n) prefix_diff += run_memory_frame(M, n, 33).logical_errors != run_surgery_frame(S, n, 33).errors[0];
                std::printf("  single patch d=3 memory-%c, %s: run_memory_frame %ld errors (detection %.5f, %.3f leakage events/shot), run_surgery_frame %ld (%.5f, %.3f); per-shot prefixes differing %d/120\n",
                            basis, v.name, a.logical_errors, a.detection_rate, a.leak_events_per_shot, b.errors[0], b.detection_rate, b.leak_events_per_shot, prefix_diff);
                CHECK(a.logical_errors == b.errors[0] && a.detection_rate == b.detection_rate && a.leak_events_per_shot == b.leak_events_per_shot && b.errors[0] > 0,
                      "surgery engine differs from the memory experiment (%c, %s)", basis, v.name);
                CHECK(prefix_diff == 0, "surgery engine differs shot by shot (%c, %s)", basis, v.name);
                CHECK(!v.nm.leak_cz.empty() == (a.leak_events_per_shot > 0), "variant %s: leakage events %.3f per shot", v.name, a.leak_events_per_shot);
            }
    }
    // d = 3 ZZ measurement on a Beta-1 physics chip: which observables each basis decodes, graph-like DEM, single-fault
    // distance, noiseless runs clean
    {
        const auto E = surgery_parity(3, 'Z', '0', '0');
        const auto D = surgery_chip(E.sched.lat, "beta1"); const auto R = analyze(D); const auto nm = noise_model(D, R);
        auto quiet = nm; quiet.disable("all");
        const auto Q = prepare_surgery(D, R, quiet, E, 'Z');
        const auto q = run_surgery_frame(Q, 300, 1);
        CHECK(q.any_errors == 0 && q.mean_defects == 0, "noiseless frame run of a ZZ measurement is not clean");
        for (char ty : {'Z', 'X'}) {
            const auto Ey = ty == 'Z' ? E : surgery_parity(3, 'Z', '+', '+');
            auto M = prepare_surgery(D, R, nm, Ey, ty);
            const auto sf = surgery_single_faults(M);
            const auto r = run_surgery_frame(M, 4000, 7);
            std::string names; for (const auto& n : r.names) names += (names.empty() ? "" : ",") + n;
            std::printf("  d=3 ZZ measurement (%c-type decoding: %s): %d detectors, %zu fault mechanisms, hyperedges %d, undetectable %d; uniform-weight decoding "
                        "leaves %d uncorrected (summed over the decoders) (total p %.2e, %d not ZZ-twirl terms);", ty, names.c_str(), r.n_detectors, sf.faults, r.hyperedge_faults,
                        r.undetectable_logical, sf.uncorrected, sf.uncorrected_p, sf.uncorrected_other);
            for (std::size_t k = 0; k < r.names.size(); ++k) std::printf(" %s %.4f", r.names[k].c_str(), r.p[k]);
            std::printf("\n");
            CHECK(names == (ty == 'Z' ? "m_ZZ,Z_A,Z_B" : "X_A X_B"), "ZZ measurement (%c) decodes %s", ty, names.c_str());
            CHECK(r.hyperedge_faults == 0 && r.undetectable_logical == 0 && sf.undetectable == 0, "ZZ measurement DEM (%c) is not graph-like", ty);
            CHECK(sf.uncorrected_other == 0 && sf.uncorrected_p < 1e-7, "ZZ measurement (%c): %d uncorrected single faults (p %.2e)", ty, sf.uncorrected, sf.uncorrected_p);
            CHECK(r.p_any > 0 && r.p_any < 0.05, "ZZ measurement (%c): implausible logical error %.4f", ty, r.p_any);
        }
        // XX on |++> decodes m_XX and both patch logicals; an observable of the wrong type is refused at setup
        {
            const auto Ex = surgery_parity(3, 'X', '+', '+');
            const auto Dx = surgery_chip(Ex.sched.lat, "beta1"); const auto Rx = analyze(Dx);
            const auto Mx = prepare_surgery(Dx, Rx, noise_model(Dx, Rx), Ex, 'X');
            std::string names; for (int k : Mx.obs) names += (names.empty() ? "" : ",") + Ex.observables[std::size_t(k)].name;
            CHECK(names == "m_XX,X_A,X_B", "XX measurement decodes %s", names.c_str());
        }
        auto Bad = E; Bad.observables[0].type = 'X';            // m_ZZ mis-typed: decoded with X-type detectors
        bool threw = false;
        try { prepare_surgery(D, R, nm, Bad, 'X'); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "an observable decoded with detectors of the wrong type was accepted");
    }
}

static void test_surgery_cnot() {
    std::puts("[surgery-cnot] three-patch lattice-surgery CNOT: logical stabiliser states and reference qubits (tableau), DEM distance, noise");
    int runs = 0, bad = 0;
    const char* zin[] = {"00", "01", "10", "11"};
    const char* xin[] = {"++", "+-", "-+", "--"};
    const char* mix[] = {"+0", "+1", "-0", "-1", "0+", "1+", "0-", "1-"};
    std::uint64_t seed = 2000;
    SurgeryRounds pre0; pre0.pre = 0;
    for (int d : {3, 5}) {
        auto run = [&](char c, char t, char ro, SurgeryRounds rr) {
            const auto E = surgery_cnot(d, c, t, ro, rr);
            const auto k = check_surgery_noiseless(E, d == 3 ? 6 : 2, ++seed);
            ++runs;
            if (surgery_failed(k)) {
                ++bad; std::printf("    %s readout %c (pre %d): fired %d, wrong observables %d\n", E.sched.name.c_str(), ro ? ro : '-', rr.pre, k.fired, k.wrong_observables);
            }
        };
        for (const char* in : zin) run(in[0], in[1], 0, {});
        for (const char* in : xin) run(in[0], in[1], 0, {});
        if (d == 3) {
            for (const char* in : mix) { run(in[0], in[1], 'Z', {}); run(in[0], in[1], 'X', {}); }
            for (const char* in : {"00", "11", "+-", "-+"}) run(in[0], in[1], 0, pre0);
        }
    }
    std::printf("  noiseless tableau check: %d CNOT experiments (Z and X basis inputs at d=3,5; Bell-pair inputs and merging from the first round at d=3), %d failed\n", runs, bad);
    CHECK(bad == 0, "noiseless CNOT: %d experiments with a fired detector or a wrong observable", bad);
    // process check with reference qubits: every patch logical random, Bell pairs with two reference qubits
    {
        int rbad = 0, rruns = 0;
        for (int d : {3, 5}) for (char prep : {'Z', 'X'}) for (char ro : {'Z', 'X'}) {
            const auto E = surgery_cnot_reference(d, prep, ro);
            const auto c = check_surgery_noiseless(E, d == 3 ? 8 : 2, ++seed);
            ++rruns;
            if (surgery_failed(c)) { ++rbad; std::printf("    %s: fired %d, wrong observables %d\n", E.sched.name.c_str(), c.fired, c.wrong_observables); }
        }
        std::printf("  reference-qubit process check: %d experiments (Z- and X-basis patches, Z and X readout, d=3,5), %d failed\n", rruns, rbad);
        CHECK(rbad == 0, "reference-qubit CNOT check: %d experiments failed", rbad);
        // mis-assembled frames the fixed-input checks cannot see: an input patch's round-0 faces of the readout type
        // appended to the corrected logical (random for that patch, so the reference correlation breaks)
        for (char ro : {'X', 'Z'}) {
            auto M = surgery_cnot_reference(3, ro == 'X' ? 'Z' : 'X', ro);
            const auto& L = M.sched.lat;
            for (const auto& f : region_faces(L.patch_region(0), L.at)) if (f.type == ro) M.observables[0].bits.push_back(M.face_bit.at({f.qubit, 0}));
            int caught = 0; for (int s = 0; s < 16; ++s) caught += check_surgery_noiseless(M, 1, std::uint64_t(300 + s)).wrong_observables > 0;
            CHECK(caught > 0, "control's round-0 %c faces in the %c observable went unnoticed (%d/16)", ro, ro, caught);
        }
    }
    {
        // truth table: |10> -> |11>, and a Bell pair from |+0>
        const auto c = check_surgery_noiseless(surgery_cnot(3, '1', '0'), 1, 2);
        CHECK(c.observable_values.size() == 2 && c.observable_values[0] == 1 && c.observable_values[1] == 1, "CNOT |10> -> |11>");
        const auto b = check_surgery_noiseless(surgery_cnot(3, '+', '0', 'Z'), 1, 2);
        CHECK(b.observable_values.size() == 1 && b.observable_values[0] == 0, "CNOT |+0>: Z_C Z_T = +1 (Bell state)");
    }
    {
        const auto E = surgery_cnot(3, '0', '0');
        const auto D = surgery_chip(E.sched.lat, "beta1"); const auto R = analyze(D); const auto nm = noise_model(D, R);
        for (char ty : {'Z', 'X'}) {
            const auto Ey = ty == 'Z' ? E : surgery_cnot(3, '+', '+');
            auto M = prepare_surgery(D, R, nm, Ey, ty);
            const auto sf = surgery_single_faults(M);
            const auto r = run_surgery_frame(M, 3000, 5);
            std::string names; for (const auto& n : r.names) names += (names.empty() ? "" : ",") + n;
            std::printf("  d=3 CNOT (%c basis: %s): %d detectors, %zu fault mechanisms, hyperedges %d (p %.1e), undetectable %d; uniform-weight decoding leaves "
                        "%d uncorrected (summed over the decoders) (total p %.2e, %d not ZZ-twirl terms);", ty, names.c_str(), r.n_detectors, sf.faults, r.hyperedge_faults, r.hyperedge_p,
                        r.undetectable_logical, sf.uncorrected, sf.uncorrected_p, sf.uncorrected_other);
            for (std::size_t k = 0; k < r.names.size(); ++k) std::printf(" %s %.4f", r.names[k].c_str(), r.p[k]);
            std::printf("\n");
            CHECK(names == (ty == 'Z' ? "Z_C,Z_T" : "X_C,X_T"), "CNOT (%c) decodes %s", ty, names.c_str());
            CHECK(r.undetectable_logical == 0 && sf.undetectable == 0 && r.hyperedge_p < 1e-5, "CNOT DEM (%c): undetectable logicals or hyperedge mass", ty);
            CHECK(sf.uncorrected_other == 0 && sf.uncorrected_p < 1e-6, "CNOT (%c): %d uncorrected single faults (p %.2e)", ty, sf.uncorrected, sf.uncorrected_p);
            CHECK(r.p_any > 0 && r.p_any < 0.1, "CNOT (%c): implausible logical error %.4f", ty, r.p_any);
        }
    }
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
    {"layout_geometry", test_layout_geometry},
    {"layout_checks", test_layout_checks},
    {"layout_presets", test_layout_presets},
    {"json_rfc", test_json_rfc},
    {"json_depth", test_json_depth},
    {"gds_write", test_gds_write},
    {"json", test_json},
    {"gds", test_gds},
    {"matching", test_matching},
    {"fault_enumeration", test_fault_enumeration},
    {"memory_bases", test_memory_bases},
    {"bursts", test_bursts},
    {"characterisation", test_characterisation},
    {"designs_v04", test_designs_v04},
    {"materials", test_materials},
    {"union_find", test_union_find},
    {"heralds", test_heralds},
    {"impact_veto", test_impact_veto},
    {"streaming_uf", test_streaming_uf},
    {"surgery", test_surgery},
    {"surgery_cnot", test_surgery_cnot},
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
