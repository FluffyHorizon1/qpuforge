// qpuforge CLI — design, analyse and emulate superconducting QPU chips
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <iostream>
#include <sstream>

using namespace qpuforge;

static void usage() {
    std::puts(
        "qpuforge — superconducting QPU chip design & whole-processor emulation\n"
        "\n"
        "  qpuforge new <chip.json> [--grid RxC | --line N | --ring N] [--f-hi GHz] [--f-lo GHz]\n"
        "                [--alpha GHz] [--g-MHz g] [--arch fixed|tunable] [--Q Q] [--plan 2|4|8]\n"
        "                [--stagger-MHz s] [--name NAME]\n"
        "  qpuforge analyze <chip.json>                      device physics report\n"
        "  qpuforge run <chip.json> <circuit> [--traj N] [--shots N] [--seed S] [--off SRC,..] [--top K]\n"
        "  qpuforge budget <chip.json> <circuit> [--traj N]  error budget by noise source\n"
        "  qpuforge fab <chip.json> [--samples N] [--seed S] fabrication Monte Carlo (junction spread)\n"
        "  qpuforge qec <chip.json> [--rounds 2,4,6,8,10] [--shots N] [--statevector N] [--off SRC,..]\n"
        "                surface-code memory experiment: logical error per round (needs a surface layout)\n"
        "  qpuforge layout <chip.json> <out.svg>             chip floor-plan\n"
        "  qpuforge cz <chip.json> [--out cal.json] [--cache file] [--g-qc MHz] [--f-on GHz] [--rise ns]\n"
        "                pulse-level CZ calibration of every pair (physical tunable coupler)\n"
        "  qpuforge gds <chip.json> <out.gds> [--cells cells.json] [--synthesize] [--cell NAME]\n"
        "                physical layout: MoM capacitance synthesis of the cells, chip assembly, DRC, GDSII\n"
        "  qpuforge new <chip.json> --preset alpha1|alpha1v11|beta1|gen0 [--distance D]    named designs\n"
        "  qpuforge new <chip.json> --surface D [--arch tunable] ...        surface-code chip\n"
        "\n"
        "  <circuit> = file.qasm | ghz:N | qft:N | bv:N | xeb:CYCLES\n"
        "  noise sources: T1 Tphi thermal zz exchange control readout reset leakage\n"
        "  qec/cz on a chip with a physical coupler calibrate CZs first (--cache file reuses calibrations)");
}

static bool has_flag(int argc, char** argv, const char* key) {
    for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], key)) return true;
    return false;
}

static std::string arg(int argc, char** argv, const char* key, const std::string& def = "") {
    for (int i = 1; i + 1 < argc; ++i) if (!std::strcmp(argv[i], key)) return argv[i + 1];
    return def;
}

static Circuit make_circuit(const std::string& spec, const ChipDesign& D, const ChipReport& R, bool& native_layout) {
    native_layout = false;
    const auto c = spec.find(':');
    if (c != std::string::npos && spec.find(".qasm") == std::string::npos) {
        const std::string kind = spec.substr(0, c);
        const int n = std::stoi(spec.substr(c + 1));
        if (kind == "ghz") return ghz(n);
        if (kind == "qft") return qft(n);
        if (kind == "bv") {   // fixed pseudo-random secret, never all-zero
            unsigned long long secret = 0x2D5B7ULL & ((1ULL << n) - 1);
            if (!secret) secret = 1;
            return bernstein_vazirani(n, secret);
        }
        if (kind == "xeb") {
            std::vector<std::pair<int, int>> edges;
            for (const auto& p : R.pairs) if (p.designed) edges.push_back({D.index_of(p.a), D.index_of(p.b)});
            native_layout = true;
            return random_xeb(int(D.qubits.size()), edges, n, 7);
        }
        throw std::runtime_error("unknown circuit spec " + spec);
    }
    return load_qasm(spec);
}

static void print_report(const ChipDesign& D, const ChipReport& R) {
    std::printf("\n%s — %zu qubits, %zu couplers, architecture %s\n", D.name.c_str(), D.qubits.size(), D.couplings.size(), D.architecture.c_str());
    std::printf("materials: Q_diel=%.2g  junction gap %.0f ueV (Rn spread %.1f%%)  T_eff=%.0f mK  flux noise %.1f uPhi0\n\n",
                D.Q_dielectric, D.gap_ueV, D.Rn_spread_pct, D.T_eff_mK, D.flux_noise_uPhi0);
    std::printf(" Q  role       C_fF   Rn_kOhm  EJ/EC   f01_GHz  dressed   alpha_MHz  disp_kHz  T1_us(diel/Purcell)      Tphi_us   T2_us   P_th    ro(1|0) ro(0|1)  e_1q\n");
    for (std::size_t i = 0; i < R.q.size(); ++i) {
        const auto& q = R.q[i]; const auto& d = D.qubits[i];
        std::printf("%2d  %-9s %6.1f  %7.3f  %6.1f   %7.4f  %7.4f   %8.1f  %8.3g  %6.1f (%6.0f/%7.0f)  %8.0f  %6.1f  %.4f  %.4f  %.4f  %.1e\n",
                    q.id, d.role.c_str(), d.C_fF, d.Rn_kOhm, q.EJ_EC, q.f01, q.f01_dressed, q.alpha * 1e3, q.charge_disp_kHz, q.T1_us, q.T1_diel_us,
                    q.T1_purcell_us, q.Tphi_us, q.T2_us, q.p_thermal, q.ro_err_0to1, q.ro_err_1to0, q.err_1q);
    }
    std::printf("\n pair     coupler   g_MHz   detune_MHz   ZZ_static_kHz  (pert.)   ZZ_idle_kHz   t_CZ_ns   CZ point GHz   e_2q     flags\n");
    for (const auto& p : R.pairs) {
        std::printf(" %2d-%-2d    %-7s  %6.3f   %9.1f   %12.2f  (%7.2f)  %11.2f   %7.1f   %11.4f   %.1e   %s\n", p.a, p.b, p.designed ? "yes" : "parasit", p.g_MHz,
                    p.delta_MHz, p.zz_kHz, p.zz_pert_kHz, p.zz_idle_kHz, p.t_cz_ns, p.cz_interaction_GHz, p.err_2q, p.flags.c_str());
    }
    std::printf("\nmean T1 %.1f us, mean T2 %.1f us, mean e_1q %.2e, mean e_2q %.2e, max |idle ZZ| %.1f kHz, collisions %d\n",
                R.mean_T1_us, R.mean_T2_us, R.mean_err_1q, R.mean_err_2q, R.max_abs_zz_kHz, R.n_collisions);
    for (const auto& w : R.warnings) std::printf("  ! %s\n", w.c_str());
}

static int cmd_new(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    if (auto pre = arg(argc, argv, "--preset"); !pre.empty()) {
        const int d = std::stoi(arg(argc, argv, "--distance", "3"));
        ChipDesign D = pre == "alpha1" ? alpha1(d) : pre == "alpha1v11" ? alpha1_v11(d) : pre == "beta1" ? beta1(d)
                     : (pre == "gen0" ? gen0_surface17(d) : throw std::runtime_error("unknown preset " + pre));
        save_design(D, argv[2]);
        std::printf("wrote %s (%zu qubits, %zu couplers)\n", argv[2], D.qubits.size(), D.couplings.size());
        print_report(D, analyze(D));
        return 0;
    }
    GenSpec s;
    if (auto sd = arg(argc, argv, "--surface"); !sd.empty()) {
        s = surface_spec(std::stoi(sd));
        s.f_hi = 5.9; s.f_lo = 5.2;
    }
    s.feedlines = std::stoi(arg(argc, argv, "--feedlines", std::to_string(s.feedlines)));
    s.kappa_MHz = std::stod(arg(argc, argv, "--kappa-MHz", std::to_string(s.kappa_MHz)));
    s.chi_MHz = std::stod(arg(argc, argv, "--chi-MHz", std::to_string(s.chi_MHz)));
    s.purcell_filter = std::stod(arg(argc, argv, "--purcell", std::to_string(s.purcell_filter)));
    if (auto g = arg(argc, argv, "--grid"); !g.empty()) { s.topology = "grid"; s.rows = std::stoi(g.substr(0, g.find('x'))); s.cols = std::stoi(g.substr(g.find('x') + 1)); }
    if (auto l = arg(argc, argv, "--line"); !l.empty()) { s.topology = "line"; s.rows = 1; s.cols = std::stoi(l); }
    if (auto l = arg(argc, argv, "--ring"); !l.empty()) { s.topology = "ring"; s.rows = 1; s.cols = std::stoi(l); }
    s.f_hi = std::stod(arg(argc, argv, "--f-hi", std::to_string(s.topology == "surface" ? s.f_hi : 5.4)));
    s.f_lo = std::stod(arg(argc, argv, "--f-lo", std::to_string(s.topology == "surface" ? s.f_lo : 4.8)));
    s.alpha = std::stod(arg(argc, argv, "--alpha", "-0.22"));
    s.g = std::stod(arg(argc, argv, "--g-MHz", "4")) / 1000.0;
    s.plan = std::stoi(arg(argc, argv, "--plan", "8"));
    s.stagger = std::stod(arg(argc, argv, "--stagger-MHz", "60")) / 1000.0;
    s.Q = std::stod(arg(argc, argv, "--Q", "3e6"));
    s.architecture = arg(argc, argv, "--arch", "fixed") == "tunable" ? "tunable-coupler-cz" : "fixed-coupling-cz";
    auto D = generate(s);
    if (auto nm = arg(argc, argv, "--name"); !nm.empty()) D.name = nm;
    if (auto v = arg(argc, argv, "--readout-ns"); !v.empty()) D.t_readout_ns = std::stod(v);
    if (auto v = arg(argc, argv, "--reset-ns"); !v.empty()) D.t_reset_ns = std::stod(v);
    if (auto v = arg(argc, argv, "--reset-err"); !v.empty()) D.reset_err = std::stod(v);
    if (auto v = arg(argc, argv, "--err2q"); !v.empty()) D.err_2q_control = std::stod(v);
    if (auto v = arg(argc, argv, "--residual-zz-kHz"); !v.empty()) D.tc_residual_zz_kHz = std::stod(v);
    for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], "--active-reset")) D.active_reset = true;
    save_design(D, argv[2]);
    std::printf("wrote %s (%zu qubits, %zu couplers)\n", argv[2], D.qubits.size(), D.couplings.size());
    print_report(D, analyze(D));
    return 0;
}

static int cmd_analyze(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    const auto D = load_design(argv[2]);
    print_report(D, analyze(D));
    return 0;
}

static RunOptions run_opts(int argc, char** argv) {
    RunOptions o;
    o.trajectories = std::stoi(arg(argc, argv, "--traj", "300"));
    o.shots_per_traj = std::stoi(arg(argc, argv, "--shots", "20"));
    o.seed = unsigned(std::stoul(arg(argc, argv, "--seed", "1234")));
    return o;
}

static int cmd_run(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const auto D = load_design(argv[2]);
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    bool native = false;
    const auto C = make_circuit(argv[3], D, R, native);
    CompileOptions co; if (native) co.layout = "trivial";
    const auto P = compile(C, G, co);
    auto nm = noise_model(D, R);
    if (auto off = arg(argc, argv, "--off"); !off.empty()) {
        std::stringstream ss(off); std::string s;
        while (std::getline(ss, s, ',')) nm.disable(s);
    }
    auto o = run_opts(argc, argv);
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = run(P, nm, o);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("\n%s on %s\n", C.name.c_str(), D.name.c_str());
    std::printf("  compiled: %d moments, %.0f ns, %d CZ, %d 1q pulses, %d virtual RZ, %d SWAPs, %d simulated qubits\n",
                r.depth, r.duration_ns, r.n_cz, r.n_pulses, P.n_virtual_rz, r.n_swaps, r.n_sim);
    std::printf("  layout (logical->physical):");
    for (std::size_t l = 0; l < P.layout_initial.size(); ++l) std::printf(" %zu->Q%d", l, D.qubits[P.layout_initial[l]].id);
    std::printf("\n  state fidelity  %.4f +- %.4f   (%d trajectories)\n", r.fidelity, r.fidelity_sem, o.trajectories);
    if (r.shots) {
        if (std::abs(r.xeb_ideal) > 0.05) {
            std::printf("  P(ideal most-likely outcome %s) = %.4f over %ld shots\n", r.ideal_top.c_str(), r.p_ideal_top, r.shots);
            std::printf("  linear XEB %.4f (normalised %.4f; ideal-circuit value %.3f)\n", r.xeb, r.xeb_norm, r.xeb_ideal);
        } else {
            std::printf("  %ld shots; ideal output distribution is uniform, so XEB / most-likely-outcome are not informative\n", r.shots);
        }
        std::vector<std::pair<long, std::string>> top;
        for (auto& [k, c] : r.counts) top.push_back({c, k});
        std::sort(top.rbegin(), top.rend());
        const int K = std::stoi(arg(argc, argv, "--top", "8"));
        std::printf("  top outcomes:\n");
        for (int i = 0; i < K && i < int(top.size()); ++i) std::printf("    %s  %6.4f\n", top[i].second.c_str(), double(top[i].first) / r.shots);
    }
    std::printf("  (%.2f s)\n", secs);
    return 0;
}

static int cmd_budget(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const auto D = load_design(argv[2]);
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    bool native = false;
    const auto C = make_circuit(argv[3], D, R, native);
    CompileOptions co; if (native) co.layout = "trivial";
    const auto P = compile(C, G, co);
    auto o = run_opts(argc, argv); o.shots_per_traj = 0;
    const auto base = noise_model(D, R);
    const auto all = run(P, base, o);
    std::printf("\nerror budget for %s on %s (%d CZ, %.0f ns)\n", C.name.c_str(), D.name.c_str(), P.n_cz, P.duration_ns);
    std::printf("  all sources: infidelity %.4f +- %.4f\n", 1 - all.fidelity, all.fidelity_sem);
    double sum = 0;
    std::vector<std::pair<double, std::string>> rows;
    for (const std::string src : {"T1", "Tphi", "thermal", "zz", "exchange", "control"}) {
        auto nm = base; nm.only(src);
        const auto r = run(P, nm, o);
        rows.push_back({1 - r.fidelity, src});
        sum += 1 - r.fidelity;
    }
    std::sort(rows.rbegin(), rows.rend());
    for (auto& [inf, src] : rows) std::printf("  %-8s alone: infidelity %.4f  (%4.1f%% of sum)\n", src.c_str(), inf, 100 * inf / std::max(sum, 1e-12));
    std::printf("  (readout assignment error acts on the measured counts, not on the state fidelity)\n");
    return 0;
}

static int cmd_fab(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    const auto D = load_design(argv[2]);
    const int N = std::stoi(arg(argc, argv, "--samples", "200"));
    const auto F = fab_monte_carlo(D, N, unsigned(std::stoul(arg(argc, argv, "--seed", "1"))));
    auto stats = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return std::array<double, 3>{v[v.size() / 10], v[v.size() / 2], v[v.size() * 9 / 10]}; };
    const auto fs = stats(F.f01_sigma_MHz), es = stats(F.err2q_mean), zs = stats(F.max_zz_kHz);
    std::printf("\nfabrication Monte Carlo: %d chips, junction Rn spread %.1f%%\n", N, D.Rn_spread_pct);
    std::printf("  collision-free yield: %d / %d = %.1f%%\n", F.collision_free, N, 100.0 * F.collision_free / N);
    std::printf("  RMS frequency error     p10/p50/p90: %.1f / %.1f / %.1f MHz\n", fs[0], fs[1], fs[2]);
    std::printf("  mean 2q error           p10/p50/p90: %.2e / %.2e / %.2e\n", es[0], es[1], es[2]);
    std::printf("  worst idle ZZ           p10/p50/p90: %.0f / %.0f / %.0f kHz\n", zs[0], zs[1], zs[2]);
    for (auto& [m, c] : F.failure_modes) std::printf("  failure mode %-18s %d pair-occurrences\n", m.c_str(), c);
    return 0;
}

static int cmd_layout(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const auto D = load_design(argv[2]);
    write_layout_svg(D, analyze(D), argv[3]);
    std::printf("wrote %s\n", argv[3]);
    return 0;
}

static int cmd_cz(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    auto D = load_design(argv[2]);
    D.coupler.physical = true;
    if (auto v = arg(argc, argv, "--g-qc"); !v.empty()) D.coupler.g_qc_MHz = std::stod(v);
    if (auto v = arg(argc, argv, "--f-on"); !v.empty()) D.coupler.f_on_GHz = std::stod(v);
    if (auto v = arg(argc, argv, "--rise"); !v.empty()) D.coupler.t_rise_ns = std::stod(v);
    auto R = analyze(D);
    const auto cal = calibrate_chip_cz(D, R, true, arg(argc, argv, "--cache"));
    std::printf("\n%s — pulse-level CZ calibration (coupler: f_max %.2f GHz, alpha %.0f MHz, g_qc %.0f MHz, on ~%.2f GHz, %.0f ns ramps)\n",
                D.name.c_str(), D.coupler.f_max_GHz, D.coupler.alpha_GHz * 1e3, D.coupler.g_qc_MHz, D.coupler.f_on_GHz, D.coupler.t_rise_ns);
    std::printf("  pair   mover  f_m/f_t GHz   g12 MHz  ZZ@park kHz  f_int  fc_on  t_cz ns  1-F coh  leak     flux     ctrl     total\n");
    for (const auto& c : cal)
        std::printf("  %2d-%-2d   q%-3d %.3f/%.3f   %6.2f    %+6.2f     %.3f  %.3f  %6.1f   %.1e  %.1e  %.1e  %.1e  %.1e\n", c.a, c.b, c.mover,
                    c.f_mover, c.f_target, c.g12_MHz, c.zz_idle_kHz, c.f_int_GHz, c.fc_on_GHz, c.t_cz_ns, c.err_coherent, c.leakage + c.leak_spectators,
                    c.err_flux_local + c.err_flux_quasistatic, c.err_ctrl, c.err_total);
    double mt = 0, me = 0, ml = 0; for (const auto& c : cal) { mt += c.t_cz_ns; me += c.err_total; ml += c.leakage + c.leak_spectators; }
    const double n = std::max<std::size_t>(cal.size(), 1);
    std::printf("  mean: t_cz %.1f ns, CZ error %.2e (incl. decoherence), leakage %.1e per CZ\n", mt / n, me / n, ml / n);
    if (auto out = arg(argc, argv, "--out"); !out.empty()) { std::ofstream f(out); cz_calibration_json(cal).dump(f); f << "\n"; std::printf("  wrote %s\n", out.c_str()); }
    return 0;
}

static int cmd_gds(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    using namespace layout;
    auto D = load_design(argv[2]);
    auto R = analyze(D);
    const std::string cells_path = arg(argc, argv, "--cells");
    CellSet C; LayoutMetrics m; SynthesisLog log;
    const auto T = targets_from_design(D, R);
    bool have = false;
    if (!cells_path.empty()) {
        std::ifstream f(cells_path);
        if (f && !has_flag(argc, argv, "--synthesize")) {
            std::stringstream ss; ss << f.rdbuf();
            const Json j = Json::parse(ss.str());
            C = cells_from_json(j); have = true;
            std::printf("loaded cell geometry from %s\n", cells_path.c_str());
        }
    }
    if (!have) {
        std::printf("synthesising cell geometry against the circuit targets (method of moments; ~5-10 min)...\n");
        C.pitch = 900;
        C = synthesize_cells(C, T, &log);
        for (const auto& s : log.lines) std::printf("  %s\n", s.c_str());
    }
    m = evaluate_cells(C, MeshOptions{1.5, 50.0, 1.6});
    if (!cells_path.empty()) {
        Json j = cells_json(C, m, T, have ? nullptr : &log);
        if (have) {                                   // keep the provenance of the loaded geometry
            std::ifstream f0(cells_path); std::stringstream s0; s0 << f0.rdbuf();
            const Json old = Json::parse(s0.str());
            if (old.has("synthesis_log")) j.set("synthesis_log", old["synthesis_log"]);
            if (old.has("max_rel_err")) j.set("max_rel_err", old["max_rel_err"]);
        }
        std::ofstream f(cells_path); j.dump(f); f << "\n";
    }
    const auto L = assemble_chip(D, C);
    write_gds(L, argv[3], arg(argc, argv, "--cell", "TOP"));
    std::printf("\n%s — layout: %.1f x %.1f mm, %zu etch polygons, %zu junction regions, %zu bump pads -> %s\n", D.name.c_str(),
                L.outline.w() / 1e3, L.outline.h() / 1e3, L.etch.size(), L.jj.size(), L.bumps.size(), argv[3]);
    std::printf("  DRC (10 um spacing, 4 um at designed coupling gaps): %s\n", L.problems.empty() ? "clean" : "VIOLATIONS");
    for (const auto& s : L.problems) std::printf("    %s\n", s.c_str());
    std::printf("  extracted circuit (fine mesh) vs targets:\n");
    auto row = [](const char* n, double a, double t) { std::printf("    %-28s %10.5f  target %10.5f  (%+.2f%%)\n", n, a, t, 100 * (a / t - 1)); };
    row("E_C measure (GHz)", m.EC_m, T.EC_m); row("E_C data (GHz)", m.EC_d, T.EC_d); row("E_C coupler (GHz)", m.EC_c, T.EC_c);
    row("G measure-coupler (GHz/nn)", m.G_mc, T.G_mc); row("G coupler-data (GHz/nn)", m.G_cd, T.G_cd); row("G measure-data (GHz/nn)", m.G_md, T.G_md);
    row("C qubit-readout (fF)", m.C_qr_fF, T.C_qr_fF);
    const auto P = post_layout(D, R, m);
    std::printf("  post-layout circuit: alpha m/d/c %.1f/%.1f/%.1f MHz, g_qc %.1f MHz, g12 %.2f MHz (exchange-free park needs %.2f),\n"
                "    residual parked |11>-|20> gap %.3f MHz -> spectator transfer %.1e per crossing, ZZ at park %+.2f kHz, readout g %.1f MHz (design %.1f)\n",
                P.alpha_m * 1e3, P.alpha_d * 1e3, P.alpha_c * 1e3, P.g_qc_MHz, P.g12_MHz, P.g12_null_MHz, P.residual_gap_MHz, P.spectator_LZ, P.zz_park_kHz,
                P.g_r_MHz, P.g_r_target_MHz);
    std::printf("  re-null by parking the couplers at %.3f GHz (%.0f MHz below the sweet spot, df/dPhi %.2f GHz/Phi0): residual gap %.3f MHz, ZZ %+.2f kHz\n",
                P.park_fc_GHz, (D.coupler.f_max_GHz - P.park_fc_GHz) * 1e3, P.park_dfdphi_GHz, P.park_gap_MHz, P.park_zz_kHz);
    for (const auto& r : L.resonators) if (std::abs(r.drawn_um - r.target_um) > 1) std::printf("    resonator R%d length %.0f um vs target %.0f um\n", r.qubit, r.drawn_um, r.target_um);
    return L.problems.empty() ? 0 : 2;
}

static int cmd_qec(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    auto D = load_design(argv[2]);
    if (auto v = arg(argc, argv, "--leakage-removal"); !v.empty()) {
        D.leakage.reset_removes_leakage = v == "lrr" || v == "lrr+dqlr";
        D.leakage.dqlr = v == "dqlr" || v == "lrr+dqlr";
    }
    auto R = analyze(D);
    if (D.coupler.physical) calibrate_chip_cz(D, R, true, arg(argc, argv, "--cache"));
    auto nm = noise_model(D, R);
    if (auto off = arg(argc, argv, "--off"); !off.empty()) { std::stringstream ss(off); std::string s; while (std::getline(ss, s, ',')) nm.disable(s); }
    std::vector<int> rounds;
    { std::stringstream ss(arg(argc, argv, "--rounds", "2,4,6,8,10")); std::string s; while (std::getline(ss, s, ',')) rounds.push_back(std::stoi(s)); }
    const long shots = std::stol(arg(argc, argv, "--shots", "200000"));
    const int sv_shots = std::stoi(arg(argc, argv, "--statevector", "0"));
    const unsigned seed = unsigned(std::stoul(arg(argc, argv, "--seed", "2026")));
    std::printf("\n%s — memory-Z experiment, distance %d, Pauli-frame Monte Carlo (%ld shots per point)\n", D.name.c_str(), D.surface_d, shots);
    std::printf("  rounds   round_ns   detection   P_L                eps/round\n");
    std::vector<MemoryResult> all;
    for (int r : rounds) {
        const auto M = prepare_memory(D, R, nm, r);
        const auto res = run_memory_frame(M, shots, seed + unsigned(r));
        all.push_back(res);
        std::printf("  %5d   %8.0f   %8.4f    %.5f +- %.5f   %.5f", r, res.round_ns, res.detection_rate, res.p_logical, res.p_logical_sem, res.eps_per_round);
        if (M.frame.leakage) std::printf("   leaked: data %.1e meas %.1e", res.leak_frac_data, res.leak_frac_meas);
        std::printf("\n");
        if (res.ignored_exchange_pairs) std::printf("    ! %d near-degenerate exchange pair(s) not modelled by the Pauli-frame simulator\n", res.ignored_exchange_pairs);
        if (sv_shots > 0 && r == rounds.front()) {
            const auto sv = run_memory_statevector(M, nm, sv_shots, seed);
            std::printf("          exact statevector (%d shots): detection %.4f, P_L %.4f +- %.4f\n", sv_shots, sv.detection_rate, sv.p_logical, sv.p_logical_sem);
        }
    }
    const double eps = fit_eps_per_round(all);
    double best_T1 = 0, best_T2 = 0;
    for (const auto& q : R.q) { best_T1 = std::max(best_T1, q.T1_us); best_T2 = std::max(best_T2, q.T2_us); }
    const double round_us = all.back().round_ns / 1e3;
    std::printf("  logical error per round (fit): %.5f   ->  logical lifetime ~ %.1f us (round %.2f us)\n", eps, round_us / std::max(eps, 1e-12), round_us);
    std::printf("  best physical qubit on chip: T1 %.1f us, T2 %.1f us\n", best_T1, best_T2);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    const std::string cmd = argv[1];
    try {
        if (cmd == "new") return cmd_new(argc, argv);
        if (cmd == "analyze") return cmd_analyze(argc, argv);
        if (cmd == "run") return cmd_run(argc, argv);
        if (cmd == "budget") return cmd_budget(argc, argv);
        if (cmd == "fab") return cmd_fab(argc, argv);
        if (cmd == "qec") return cmd_qec(argc, argv);
        if (cmd == "layout") return cmd_layout(argc, argv);
        if (cmd == "cz") return cmd_cz(argc, argv);
        if (cmd == "gds") return cmd_gds(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    usage();
    return 1;
}
