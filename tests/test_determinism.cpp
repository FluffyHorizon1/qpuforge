// Thread-count determinism guard (not one of the validation checks): the Monte Carlo must give bit-identical results
// whatever the OpenMP thread count. See docs/design/deterministic-mc.md.
//
//   test_determinism            run every workload at 1, 2, 3, 4 and 8 threads in-process (omp_set_num_threads) and
//                               compare the results bit for bit; exit 1 on any difference
//   test_determinism --digest   print one line of hex-float results at the current OMP_NUM_THREADS (used by
//                               tests/compare_thread_digests.cmake, which runs it under OMP_NUM_THREADS=1,2,4,8)
//
// Workloads: run_memory_frame on a distance-3 memory with leakage on (exercises the leaked-time double sums) and a shot
// count that is not a multiple of the block size, in both memory bases and with two decode methods, including the setup
// (parallel fault enumeration -> DEM) and the 0.4 syndrome statistics; run_surgery_frame on a d = 3 lattice-surgery ZZ
// measurement with leakage on (three observables, each with its own decoder); and the trajectory emulator's sampled XEB (a double
// reduction).
#include "qpuforge/qpuforge.hpp"
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace qpuforge;

static std::string hex(double v) { char b[64]; std::snprintf(b, sizeof b, "%a", v); return b; }

static std::string frame_digest() {
    GenSpec s; s.topology = "surface"; s.distance = 3; s.architecture = "tunable-coupler-cz"; s.f_hi = 5.9; s.f_lo = 5.2; s.g = 0.004;
    s.feedlines = 3; s.kappa_MHz = 4; s.chi_MHz = 4; s.purcell_filter = 100;
    auto D = generate(s);
    D.t_readout_ns = 500; D.active_reset = true; D.reset_err = 0.003; D.t_reset_ns = 200;
    const auto R = analyze(D);
    auto nm = noise_model(D, R);
    nm.leakage = true; nm.leak_1q = 1e-3; nm.leak_transport = 0.2;
    for (const auto& p : R.pairs) if (p.designed) {
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
        nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
    }
    std::string out;
    // Both memory bases, the default (Auto) and an explicit decode method, and the 0.4 syndrome statistics. The DEM is
    // included because setup (enumerate_faults) runs in parallel: its output must not depend on the thread count either.
    for (int rounds : {2, 5})
        for (char basis : {'Z', 'X'}) {
            const auto M = prepare_memory(D, R, nm, rounds, basis);
            FrameRunOptions fo;
            if (rounds == 5) fo.method = DecodeMethod::Blossom;
            const auto r = run_memory_frame(M, 5000, 17 + unsigned(rounds), fo);   // 5000 = 4 full blocks + a partial one
            double pe = 0; for (const auto& e : M.dem.edges) pe += e.p;
            out += std::string("r") + std::to_string(rounds) + basis + " err=" + std::to_string(r.logical_errors) +
                   " det=" + hex(r.detection_rate) + " ld=" + hex(r.leak_frac_data) + " lm=" + hex(r.leak_frac_meas) +
                   " lev=" + hex(r.leak_events_per_shot) + " md=" + hex(r.mean_defects) + " mx=" + std::to_string(r.max_defects) +
                   " o22=" + hex(r.frac_over_22) + " dem=" + std::to_string(M.dem.edges.size()) + ":" + hex(pe) + " byr=";
            for (double v : r.detection_by_round) out += hex(v) + ",";
            out += " ";
        }
    return out;
}

static std::string xeb_digest() {
    GenSpec s; s.rows = 2; s.cols = 3; s.g = 0.012;
    auto D = generate(s); D.Q_dielectric = 1e6;
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    const auto nm = noise_model(D, R);
    std::vector<std::pair<int, int>> edges;
    for (const auto& p : R.pairs) if (p.designed) edges.push_back({D.index_of(p.a), D.index_of(p.b)});
    CompileOptions co; co.layout = "trivial";
    const auto P = compile(random_xeb(6, edges, 8, 42), G, co);
    RunOptions o; o.trajectories = 300; o.shots_per_traj = 20; o.seed = 9;
    const auto r = run(P, nm, o);
    return "xeb=" + hex(r.xeb) + " F=" + hex(r.fidelity) + " sem=" + hex(r.fidelity_sem) + " top=" + hex(r.p_ideal_top) +
           " shots=" + std::to_string(r.shots);
}

static std::string surgery_digest() {
    // the chip and the decoders are built once per process; only the Monte Carlo is repeated per thread count
    static const auto setup = [] {
        const auto E = surgery_parity(3, 'Z', '0', '1');
        const auto D = surgery_chip(E.sched.lat, "alpha1_v11");
        const auto R = analyze(D);
        auto nm = noise_model(D, R);
        nm.leakage = true; nm.leak_1q = 1e-3; nm.leak_transport = 0.2;
        for (const auto& p : R.pairs) if (p.designed) {
            const int a = D.index_of(p.a), b = D.index_of(p.b);
            const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
            nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
        }
        return std::make_shared<SurgerySetup>(prepare_surgery(D, R, nm, E, 'Z'));
    }();
    const auto r = run_surgery_frame(*setup, 3000, 23);   // 2 full blocks + a partial one
    std::string out = "surgery";
    for (std::size_t k = 0; k < r.names.size(); ++k) out += " " + std::to_string(r.errors[k]);
    out += " any=" + std::to_string(r.any_errors) + " def=" + hex(r.mean_defects) + " lev=" + hex(r.leak_events_per_shot) + " ";
    return out;
}

static std::string digest() { return frame_digest() + surgery_digest() + xeb_digest(); }

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "--digest") == 0) { std::printf("%s\n", digest().c_str()); return 0; }
#ifndef _OPENMP
    std::printf("determinism: built without OpenMP, nothing to compare (serial build is trivially deterministic)\n");
    return 0;
#else
    std::string ref; int fails = 0;
    for (int nth : {1, 2, 3, 4, 8}) {
        omp_set_num_threads(nth);
        const std::string d = digest();
        if (ref.empty()) { ref = d; std::printf("  1 thread : %s\n", d.c_str()); continue; }
        if (d != ref) { ++fails; std::printf("FAIL determinism: %d threads differ from 1 thread:\n  %s\n", nth, d.c_str()); }
        else std::printf("  %d threads: identical\n", nth);
    }
    if (fails) return 1;
    std::printf("determinism: ok (bit-identical at 1, 2, 3, 4, 8 threads)\n");
    return 0;
#endif
}
