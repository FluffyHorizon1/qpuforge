// Alpha-1 design review
//   1. improvement ladder: Gen-0 (best fixed-coupling design) -> + tunable couplers -> + fast readout -> + active reset (= Alpha-1)
//   2. Alpha-1 QEC performance: logical error per round, logical lifetime vs best physical qubit
//   3. QEC error budget: one noise source at a time
//   4. sensitivity: which single improvement buys the most for the next revision
//   5. distance-5 projection with the same physics: suppression factor Lambda = eps(d=3) / eps(d=5)
//   6. fabrication yield (junction spread) and exact-statevector cross-checks of the Pauli-frame model
// Writes results/alpha1_review.json and prints a report.
#include "qpuforge/qpuforge.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>

using namespace qpuforge;

struct QecPoint { double eps{0}, eps_err{0}, det{0}, round_ns{0}; };

static QecPoint qec(const ChipDesign& D, const NoiseModel& nm, long shots, unsigned seed) {
    const auto R = analyze(D);
    std::vector<MemoryResult> rs;
    for (int r : {2, 4, 6, 8, 10}) {
        const auto M = prepare_memory(D, R, nm, r);
        rs.push_back(run_memory_frame(M, shots, seed + unsigned(r)));
    }
    QecPoint q;
    q.eps = fit_eps_per_round(rs);
    const auto& last = rs.back();
    q.eps_err = last.p_logical_sem / last.rounds / std::max(1e-9, std::pow(1 - 2 * last.p_logical, double(last.rounds - 1) / last.rounds));
    q.det = last.detection_rate;
    q.round_ns = last.round_ns;
    return q;
}
static QecPoint qec(const ChipDesign& D, long shots, unsigned seed) { return qec(D, noise_model(D, analyze(D)), shots, seed); }

static ChipDesign ladder(bool tunable, bool fast_readout, bool active_reset) {
    GenSpec s = surface_spec(3);
    s.architecture = tunable ? "tunable-coupler-cz" : "fixed-coupling-cz";
    s.g = tunable ? 0.004 : 0.003;
    s.kappa_MHz = fast_readout ? 4.0 : 1.0; s.chi_MHz = fast_readout ? 4.0 : 1.0; s.purcell_filter = fast_readout ? 100 : 30;
    ChipDesign D = generate(s);
    D.tc_g_max_MHz = 40; D.tc_residual_zz_kHz = 5;
    D.t_readout_ns = fast_readout ? 500 : 800; D.readout_sep_err = fast_readout ? 0.005 : 0.01;
    D.active_reset = active_reset; D.t_reset_ns = active_reset ? 200 : 300; D.reset_err = active_reset ? 0.003 : 0.01;
    return D;
}

int main(int argc, char** argv) {
    const long shots = argc > 1 ? std::atol(argv[1]) : 100000;
    const int sv_shots = argc > 2 ? std::atoi(argv[2]) : 300;
    const auto t0 = std::chrono::steady_clock::now();
    Json out = Json::object();

    const ChipDesign A = alpha1(3);
    const ChipReport RA = analyze(A);
    save_design(A, "examples/alpha1.json");
    write_layout_svg(A, RA, "examples/alpha1_layout.svg", "17 qubits: 9 data + 8 measure, 24 tunable couplers, distance-3 surface code");

    // ------------------------------------------------------------------ 1. improvement ladder
    std::puts("1. Improvement ladder (distance-3 memory, Pauli-frame Monte Carlo, eps = logical error per round)");
    struct Step { const char* name; bool t, f, r; };
    const Step steps[] = {{"Gen-0: fixed coupling g=3 MHz, 800 ns readout, measurement-based reset", false, false, false},
                          {"+ tunable couplers (CZ 23 ns, idle ZZ 5 kHz)", true, false, false},
                          {"+ fast Purcell-filtered readout (500 ns, 2chi = kappa = 4 MHz)", true, true, false},
                          {"+ active reset (200 ns, 0.3%)  = Alpha-1", true, true, true}};
    Json lad = Json::array();
    for (const auto& st : steps) {
        const auto D = ladder(st.t, st.f, st.r);
        const auto R = analyze(D);
        const auto q = qec(D, shots, 11);
        std::printf("   %-74s eps %.4f +- %.4f   round %4.0f ns   collisions %d\n", st.name, q.eps, q.eps_err, q.round_ns, R.n_collisions);
        std::fflush(stdout);
        Json j = Json::object(); j.set("step", st.name); j.set("eps", q.eps); j.set("eps_err", q.eps_err); j.set("round_ns", q.round_ns);
        j.set("detection_rate", q.det); j.set("collisions", R.n_collisions); lad.push(j);
    }
    out.set("ladder", lad);

    // ------------------------------------------------------------------ 2. Alpha-1 headline numbers
    const auto nmA = noise_model(A, RA);
    const auto qa = qec(A, nmA, shots * 2, 101);
    double bestT1 = 0, bestT2 = 0, meanT1 = 0;
    for (const auto& q : RA.q) { bestT1 = std::max(bestT1, q.T1_us); bestT2 = std::max(bestT2, q.T2_us); meanT1 += q.T1_us / RA.q.size(); }
    const double TL = qa.round_ns / 1e3 / qa.eps;
    std::printf("\n2. Alpha-1: eps %.4f +- %.4f per %.0f ns round, detection rate %.4f\n", qa.eps, qa.eps_err, qa.round_ns, qa.det);
    std::printf("   logical lifetime ~%.0f us vs best physical T1 %.0f us / T2 %.0f us (mean T1 %.0f us)\n", TL, bestT1, bestT2, meanT1);
    Json a = Json::object();
    a.set("eps", qa.eps); a.set("eps_err", qa.eps_err); a.set("round_ns", qa.round_ns); a.set("detection_rate", qa.det);
    a.set("logical_lifetime_us", TL); a.set("best_T1_us", bestT1); a.set("best_T2_us", bestT2); a.set("mean_T1_us", meanT1);
    a.set("mean_e2q", RA.mean_err_2q); a.set("mean_e1q", RA.mean_err_1q); a.set("collisions", RA.n_collisions); a.set("max_idle_zz_kHz", RA.max_abs_zz_kHz);
    out.set("alpha1", a);

    // ------------------------------------------------------------------ 3. error budget
    // At d=3 a logical failure needs two faults, usually from different sources, so "one source at a time"
    // numbers do not add up to the total. Two views that do mean something:
    //   detection share : detection events caused by each source alone (detection rates add almost linearly)
    //   leave-one-out   : how much eps drops if that source alone were eliminated, everything else unchanged
    std::puts("\n3. Alpha-1 QEC error budget");
    Json bud = Json::array();
    struct Row { std::string src; double det, eps_without; };
    std::vector<Row> rows;
    double det_sum = 0;
    for (const std::string src : {"T1", "Tphi", "thermal", "zz", "control", "readout", "reset"}) {
        auto only = nmA; only.only(src);
        if (src == "thermal") only.T1 = true;
        double det = qec(A, only, shots / 2, 202).det;
        if (src == "thermal") { auto t1 = nmA; t1.only("T1"); det = std::max(0.0, det - qec(A, t1, shots / 2, 202).det); }
        auto without = nmA; without.disable(src);
        rows.push_back({src, det, qec(A, without, shots, 212).eps});
        det_sum += det;
    }
    std::sort(rows.begin(), rows.end(), [](const Row& x, const Row& y) { return x.eps_without < y.eps_without; });
    std::printf("   %-8s  %-22s  eps if eliminated (now %.4f)\n", "source", "share of detections", qa.eps);
    for (const auto& r : rows) {
        std::printf("   %-8s  %5.1f%%                  %.4f  (%+.0f%%)\n", r.src.c_str(), 100 * r.det / det_sum, r.eps_without, 100 * (r.eps_without / qa.eps - 1));
        Json j = Json::object(); j.set("source", r.src); j.set("detection_share", r.det / det_sum); j.set("eps_if_eliminated", r.eps_without);
        j.set("relative", r.eps_without / qa.eps - 1); bud.push(j);
    }
    out.set("budget", bud);

    // ------------------------------------------------------------------ 4. sensitivity
    std::puts("\n4. Sensitivity: one improvement at a time from Alpha-1");
    struct Knob { const char* name; std::function<void(ChipDesign&)> f; };
    const std::vector<Knob> knobs = {
        {"CZ control error 0.20% -> 0.10%", [](ChipDesign& D) { D.err_2q_control = 1e-3; }},
        {"dielectric Q 3e6 -> 6e6 (T1 x2)", [](ChipDesign& D) { D.Q_dielectric = 6e6; }},
        {"readout 500 ns -> 300 ns", [](ChipDesign& D) { D.t_readout_ns = 300; }},
        {"readout separation error 0.5% -> 0.25%", [](ChipDesign& D) { D.readout_sep_err = 0.0025; }},
        {"reset error 0.3% -> 0.15%", [](ChipDesign& D) { D.reset_err = 0.0015; }},
        {"T_eff 45 mK -> 35 mK (colder qubits & resonators)", [](ChipDesign& D) { D.T_eff_mK = 35; }},
        {"1q control error 0.01% -> 0.005%", [](ChipDesign& D) { D.err_1q_control = 5e-5; }},
        {"residual coupler ZZ 5 kHz -> 1 kHz", [](ChipDesign& D) { D.tc_residual_zz_kHz = 1; }},
    };
    Json sens = Json::array();
    for (const auto& k : knobs) {
        ChipDesign D = A; k.f(D);
        const auto q = qec(D, shots, 303);
        std::printf("   %-52s eps %.4f  (%+.1f%%)\n", k.name, q.eps, 100 * (q.eps / qa.eps - 1));
        std::fflush(stdout);
        Json j = Json::object(); j.set("change", k.name); j.set("eps", q.eps); j.set("relative", q.eps / qa.eps - 1); sens.push(j);
    }
    out.set("sensitivity", sens);

    // ------------------------------------------------------------------ 5. distance-5 projection
    std::puts("\n5. Distance-5 projection (49 qubits, same physics)");
    const ChipDesign A5 = alpha1(5);
    const ChipReport R5 = analyze(A5);
    const auto q5 = qec(A5, noise_model(A5, R5), shots / 2, 404);
    const double Lambda = qa.eps / q5.eps;
    std::printf("   d=5: eps %.4f +- %.4f per round (collisions %d)  ->  Lambda = eps3/eps5 = %.2f  (%s)\n", q5.eps, q5.eps_err, R5.n_collisions, Lambda,
                Lambda > 1 ? "below threshold: bigger codes help" : "above threshold: bigger codes do not help yet");
    Json d5 = Json::object(); d5.set("eps", q5.eps); d5.set("eps_err", q5.eps_err); d5.set("Lambda", Lambda); d5.set("collisions", R5.n_collisions);
    out.set("distance5", d5);

    // ------------------------------------------------------------------ 6. fabrication yield + statevector cross-checks
    std::puts("\n6. Fabrication Monte Carlo (Rn spread) and exact-statevector cross-checks");
    Json fab = Json::array();
    for (double spread : {1.0, 2.0, 4.0}) {
        ChipDesign D = A; D.Rn_spread_pct = spread;
        const auto F = fab_monte_carlo(D, 200, 7);
        std::printf("   Rn spread %.0f%%: collision-free yield %.0f%%\n", spread, 100.0 * F.collision_free / F.samples);
        std::fflush(stdout);
        Json j = Json::object(); j.set("Rn_spread_pct", spread); j.set("yield", double(F.collision_free) / F.samples);
        Json fm = Json::object(); for (auto& [m, c] : F.failure_modes) fm.set(m, c); j.set("failure_modes", fm); fab.push(j);
    }
    out.set("fab", fab);
    Json xc = Json::array();
    for (int which = 0; which < 2 && sv_shots > 0; ++which) {
        const ChipDesign D = which == 0 ? A : ladder(false, false, false);
        const auto R = analyze(D); const auto nm = noise_model(D, R);
        const auto M = prepare_memory(D, R, nm, 2);
        const auto fr = run_memory_frame(M, 400000, 505);
        const auto sv = run_memory_statevector(M, nm, sv_shots, 606);
        std::printf("   %-8s R=2: P_L statevector %.4f +- %.4f | Pauli frame %.4f ; detection %.4f | %.4f\n", which == 0 ? "Alpha-1" : "Gen-0",
                    sv.p_logical, sv.p_logical_sem, fr.p_logical, sv.detection_rate, fr.detection_rate);
        std::fflush(stdout);
        Json j = Json::object(); j.set("design", which == 0 ? "Alpha-1" : "Gen-0"); j.set("sv_PL", sv.p_logical); j.set("sv_PL_sem", sv.p_logical_sem);
        j.set("frame_PL", fr.p_logical); j.set("sv_det", sv.detection_rate); j.set("frame_det", fr.detection_rate); xc.push(j);
    }
    out.set("crosscheck", xc);

    std::ofstream f("results/alpha1_review.json");
    out.dump(f); f << "\n";
    std::printf("\nwrote results/alpha1_review.json, examples/alpha1.json, examples/alpha1_layout.svg  (%.0f s)\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return 0;
}
