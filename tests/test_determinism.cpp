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
// (parallel fault enumeration -> DEM) and the 0.4 syndrome statistics; the same memory decoded with union-find (batch and
// streaming), with heralded leakage (fixed factors, weights derived from the leakage model, DQLR with heralding), with
// the streaming decoder's operation profile, and with per-shot impact injection and a calibrated veto; run_surgery_frame on a
// d = 3 lattice-surgery ZZ measurement and a d = 3 three-patch CNOT with leakage on (each observable has its own decoder);
// the ZZ shots again through run_surgery_decoders with every surgery decoder (joint matching, union-find, two-pass
// correlated), and an XX measurement through a routing strip decoded jointly; and the trajectory emulator's sampled XEB
// (a double reduction).
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
            auto M = prepare_memory(D, R, nm, rounds, basis);
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
            if (basis != 'Z') continue;
            FrameRunOptions uf; uf.method = DecodeMethod::UnionFind;                          // union-find decoder
            out += "uf=" + std::to_string(run_memory_frame(M, 3000, 23 + unsigned(rounds), uf).logical_errors) + " ";
            const StreamingUnionFind S(int(M.exp.detectors.size()), M.dem.edges, M.exp.detector_round, 2, 1);   // streaming
            FrameRunOptions so; so.stream = &S;
            out += "suf=" + std::to_string(run_memory_frame(M, 3000, 23 + unsigned(rounds), so).logical_errors) + " ";
            FrameRunOptions h; h.herald.on = true;                                             // heralded leakage
            const auto rh = run_memory_frame(M, 3000, 29 + unsigned(rounds), h);
            out += "her=" + std::to_string(rh.logical_errors) + "/" + hex(rh.heralds_per_shot) + "/" + hex(rh.herald_true_frac) + " ";
            // herald weights derived from the leakage model (its signatures are computed in parallel), also through the
            // streaming decoder; DQLR that leaves the leakage on the measure qubit (its own sampling path)
            prepare_herald_model(M);
            double hm = 0; for (const auto& v : M.herald_model.cond) for (const auto& ep : v) hm += ep.second;
            FrameRunOptions hd = h; hd.herald.derived = true;
            FrameRunOptions hs = hd; hs.stream = &S;
            out += "hder=" + hex(hm) + "/" + std::to_string(run_memory_frame(M, 3000, 29 + unsigned(rounds), hd).logical_errors) + "/" +
                   std::to_string(run_memory_frame(M, 3000, 29 + unsigned(rounds), hs).logical_errors) + " ";
            {
                auto nt = nm; nt.reset_removes_leakage = true; nt.dqlr = true; nt.dqlr_eff = 0.98; nt.dqlr_err = 5e-4; nt.dqlr_transfer = true;
                auto Mt = prepare_memory(D, R, nt, rounds, basis);
                prepare_herald_model(Mt);
                const auto rt = run_memory_frame(Mt, 3000, 41 + unsigned(rounds), hd);
                out += "dqh=" + std::to_string(rt.logical_errors) + "/" + hex(rt.heralds_per_shot) + "/" + hex(rt.leak_frac_meas) + " ";
            }
            // operation profile of the streaming decoder (integer histograms merged across threads), with the pre-decoder
            {
                UFCostModel cm; cm.fixed = 3; for (long& v : cm.per) v = 1;
                const auto P = run_stream_profile(M, S, 3000, 23 + unsigned(rounds), true, {cm});
                out += "prof=" + std::to_string(P.logical_errors) + "/" + std::to_string(P.windows_cleared) + "/" + hex(StreamProfile::mean(P.cycles_hist[0])) + "/" +
                       std::to_string(StreamProfile::quantile(P.round_hist[0], 0.999)) + "/" + std::to_string(StreamProfile::max(P.latency_hist[0])) + "/" +
                       hex(StreamProfile::mean(P.stream_latency_hist[0])) + "/" + std::to_string(StreamProfile::max(P.stream_late_hist[0])) + " ";
            }
            VetoSpec V; V.threshold = calibrate_veto(M, 3000, 31 + unsigned(rounds), 5e-3, V).threshold; V.on = true;   // impacts + veto
            FrameRunOptions iv; iv.impact.p_shot = 0.1; iv.impact.burst.radius = 1.5; iv.veto = V;
            const auto ri = run_memory_frame(M, 3000, 37 + unsigned(rounds), iv);
            out += "imp=" + std::to_string(V.threshold) + "/" + std::to_string(ri.impacts) + "/" + std::to_string(ri.impacts_vetoed) + "/" +
                   std::to_string(ri.discarded) + "/" + std::to_string(ri.errors_kept) + "/" + std::to_string(ri.errors_impacted) + " ";
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

// the noise model of the surgery digests: Alpha-1 v1.1 physics with leakage on
static std::shared_ptr<SurgerySetup> surgery_setup(const SurgeryExperiment& E, bool correlated = false) {
    const auto D = surgery_chip(E.sched.lat, "alpha1_v11");
    const auto R = analyze(D);
    auto nm = noise_model(D, R);
    nm.leakage = true; nm.leak_1q = 1e-3; nm.leak_transport = 0.2;
    for (const auto& p : R.pairs) if (p.designed) {
        const int a = D.index_of(p.a), b = D.index_of(p.b);
        const int mover = D.qubits[std::size_t(a)].role == "data" ? b : a;
        nm.leak_cz[{std::min(a, b), std::max(a, b)}] = {mover, 2e-3};
    }
    return std::make_shared<SurgerySetup>(prepare_surgery(D, R, nm, E, surgery_decode_basis(E), correlated));
}

static std::string surgery_line(const char* name, const SurgerySetup& setup) {
    const auto r = run_surgery_frame(setup, 3000, 23);   // 2 full blocks + a partial one
    std::string out = name;
    for (std::size_t k = 0; k < r.names.size(); ++k) out += " " + std::to_string(r.errors[k]);
    out += " any=" + std::to_string(r.any_errors) + " def=" + hex(r.mean_defects) + " lev=" + hex(r.leak_events_per_shot) + " ";
    return out;
}

// every surgery decoder on the same shots (decode_ns is wall time and is left out), with the setup's joint model
static std::string surgery_decoders_line(const SurgerySetup& setup) {
    double pj = 0; for (const auto& e : setup.jdem.edges) pj += e.p;
    std::string out = "jdem=" + std::to_string(setup.jdem.edges.size()) + ":" + hex(pj) + ":" + std::to_string(setup.correlated_pairs) + " dec=";
    for (const auto& q : run_surgery_decoders(setup, 3000, 29, {SurgeryDecoder::Separate, SurgeryDecoder::Joint, SurgeryDecoder::SeparateUF, SurgeryDecoder::JointUF,
                                                               SurgeryDecoder::Correlated})) {
        for (long e : q.errors) out += std::to_string(e) + ",";
        out += std::to_string(q.any_errors) + "/" + std::to_string(q.differ) + "/" + std::to_string(q.any_fixed) + "/" + std::to_string(q.any_broken) + ";";
    }
    return out + " ";
}

static std::string surgery_digest() {
    // the chips and the decoders are built once per process; only the Monte Carlo is repeated per thread count.
    // The ZZ measurement has three observables; the three-patch CNOT has hyperedges in its DEM and a different observable count.
    // The ZZ setup also carries the correlated decoder's tables, which do not change the separate decoders.
    static const auto zz = surgery_setup(surgery_parity(3, 'Z', '0', '1'), true);
    static const auto cnot = surgery_setup(surgery_cnot(3, '0', '0'));
    // routing strip: XX through a strip of length 3 on Z-basis inputs (the conserved Z_A Z_B crosses the strip), joint decoder
    static const auto routed = surgery_setup(surgery_parity(3, 'X', '0', '1', {}, 3));
    const auto rr = run_surgery_frame(*routed, 3000, 31, SurgeryDecoder::Joint);
    return surgery_line("surgery", *zz) + surgery_decoders_line(*zz) + surgery_line("cnot", *cnot) +
           "strip=" + std::to_string(rr.errors[0]) + "/" + std::to_string(rr.any_errors) + "/" + hex(rr.mean_defects) + " ";
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
