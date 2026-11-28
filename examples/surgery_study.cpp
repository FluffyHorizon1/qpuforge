// Lattice-surgery study: logical ZZ and XX measurement between two patches and the three-patch CNOT, at distance 3 and
// 5, on chips with Alpha-1 v1.1 and Beta-1 physics (surgery.hpp, docs/design/lattice-surgery.md).
//
//   surgery_study [shots] [--out DIR] [--cache FILE] [--smoke]     (review_io.hpp; default 100000 shots per point)
//
// Writes results/surgery_study.txt and results/surgery_study.json (or into DIR).
//   1. noiseless check: every logical input of every operation on the stabiliser-tableau simulator (detectors silent,
//      observables equal to their expected values)
//   2. noisy runs: logical error per operation from the Pauli-frame Monte Carlo, one exact matching decoder per
//      observable, for ZZ on |00>, XX on |++>, CNOT on |00> (Z-basis readout) and on |++> (X-basis readout)
#include "qpuforge/qpuforge.hpp"
#include "review_io.hpp"
#include <cstdarg>
#include <fstream>
#include <map>
#include <memory>

using namespace qpuforge;

static std::FILE* g_txt = nullptr;
static void out(const char* fmt, ...) {
    char buf[4096];
    va_list ap; va_start(ap, fmt); std::vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    std::fputs(buf, stdout); std::fflush(stdout);
    if (g_txt) std::fputs(buf, g_txt);
}

struct Point { std::string chip, op, input; int d{0}; char basis{'Z'}; std::uint64_t seed{0}; SurgeryResult r; std::size_t faults{0}; };

int main(int argc, char** argv) {
    const auto t_start = std::chrono::steady_clock::now();
    const auto opt = review::parse(argc, argv, 100000);
    const std::string txt_path = opt.path("results/surgery_study.txt"), json_path = opt.path("results/surgery_study.json");
    g_txt = std::fopen(txt_path.c_str(), "w");
    if (!g_txt) { std::fprintf(stderr, "cannot write %s\n", txt_path.c_str()); return 1; }
    const std::vector<int> dists = opt.smoke ? std::vector<int>{3} : std::vector<int>{3, 5};
    out("qpuforge lattice-surgery study: logical ZZ / XX measurement and CNOT (surgery.hpp)\n");
    out("schedule: 1 round of separate patches, d rounds merged (ZZ, XX; the CNOT has a ZZ and then an XX merge of d rounds each), 1 round after the split,\n");
    out("transversal preparation and readout; CZ-based syndrome extraction with the memory experiments' schedule\n\n");

    // ---------------------------------------------------------------- 1. noiseless check
    Json jcheck = Json::array();
    out("1. Noiseless check (stabiliser tableau, 4 random measurement branches per experiment)\n");
    out("   op     d   inputs (readout)   detectors   firings   observables (value / expected)\n");
    int checked = 0, failed = 0;
    auto check = [&](const std::string& op, int d, const std::string& in, char ro) {
        const auto E = surgery_experiment(op, d, in, ro);
        const auto c = check_surgery_noiseless(E, 4, 77);
        ++checked;
        const bool ok = !(c.fired || c.wrong_observables || c.wrong_reference || c.random_reference);
        if (!ok) ++failed;
        std::string obs;
        for (std::size_t k = 0; k < E.observables.size(); ++k)
            obs += (k ? ", " : "") + E.observables[k].name + " " + std::to_string(c.observable_values[k]) + "/" + std::to_string(E.observables[k].expected);
        out("   %-5s %2d   |%s> %-12s %8zu   %7d   %s%s\n", op.c_str(), d, in.c_str(), ro ? (std::string("(") + ro + ")").c_str() : "", E.detectors.size(), c.fired,
            obs.c_str(), ok ? "" : "   FAILED");
        Json j = Json::object();
        j.set("op", op); j.set("d", d); j.set("inputs", in); j.set("readout", ro ? std::string(1, ro) : std::string(""));
        j.set("detectors", int(E.detectors.size())); j.set("firings", c.fired); j.set("wrong_observables", c.wrong_observables);
        j.set("wrong_reference_bits", c.wrong_reference + c.random_reference); j.set("ok", ok);
        jcheck.push(j);
    };
    for (int d : dists) {
        for (const char* in : {"00", "01", "10", "11"}) check("zz", d, in, 0);
        for (const char* in : {"++", "+-"}) check("zz", d, in, 0);
        for (const char* in : {"++", "+-", "-+", "--"}) check("xx", d, in, 0);
        for (const char* in : {"00", "01"}) check("xx", d, in, 0);
        for (const char* in : {"00", "01", "10", "11", "++", "+-", "-+", "--"}) check("cnot", d, in, 0);
        for (const char* in : {"+0", "+1", "-0"}) { check("cnot", d, in, 'Z'); check("cnot", d, in, 'X'); }
    }
    out("   %d experiments, %d failed\n\n", checked, failed);

    // ---------------------------------------------------------------- 2. noisy runs
    out("2. Logical error per operation, Pauli-frame Monte Carlo, %ld shots per point\n", opt.shots);
    out("   decoding: detectors of the observables' Pauli type, one exact matching decoder (blossom) per observable\n");
    std::vector<Point> pts;
    std::map<std::string, std::pair<ChipDesign, ChipReport>> chips;
    struct Job { std::string op, in; char basis; int kind; std::string decoded; };
    const std::vector<Job> jobs = {{"zz", "00", 'Z', 0, "m_ZZ,Z_A,Z_B"}, {"xx", "++", 'X', 1, "m_XX,X_A,X_B"},
                                   {"cnot", "00", 'Z', 2, "Z_C,Z_T"}, {"cnot", "++", 'X', 3, "X_C,X_T"}};
    for (const std::string chip : {"alpha1_v11", "beta1"}) {
        for (int d : dists) {
            out("\n   %s physics, d = %d\n", chip == "beta1" ? "Beta-1" : "Alpha-1 v1.1", d);
            out("   operation          qubits  rounds  dur_us  detectors  det_rate  defects  hyperedge_p   logical error per observable                                       any wrong\n");
            for (const auto& jb : jobs) {
                const auto E = surgery_experiment(jb.op, d, jb.in);
                const std::string key = chip + "/" + std::to_string(d) + "/" + (jb.op == "cnot" ? "cnot" : jb.op);
                if (!chips.count(key)) {
                    auto D = surgery_chip(E.sched.lat, chip);
                    auto R = analyze(D);
                    if (D.coupler.physical) calibrate_chip_cz(D, R, false, opt.cache);
                    chips.emplace(key, std::make_pair(std::move(D), std::move(R)));
                }
                const auto& [D, R] = chips.at(key);
                const auto nm = noise_model(D, R);
                const auto M = prepare_surgery(D, R, nm, E, jb.basis);   // throws on an observable with undetectable flips
                std::string decoded;
                for (int k : M.obs) decoded += (decoded.empty() ? "" : ",") + E.observables[std::size_t(k)].name;
                if (decoded != jb.decoded) { std::fprintf(stderr, "%s |%s>: decodes %s, expected %s\n", jb.op.c_str(), jb.in.c_str(), decoded.c_str(), jb.decoded.c_str()); return 1; }
                Point p; p.chip = chip; p.op = jb.op; p.input = jb.in; p.d = d; p.basis = jb.basis;
                p.seed = std::uint64_t(1000 + 100 * d + 10 * jb.kind + (chip == "beta1" ? 5 : 0));
                p.r = run_surgery_frame(M, opt.shots, p.seed);
                p.faults = M.n_faults;
                std::string obs;
                for (std::size_t k = 0; k < p.r.names.size(); ++k) {
                    char b[96]; std::snprintf(b, sizeof b, "%s%s %.5f+-%.5f", k ? ", " : "", p.r.names[k].c_str(), p.r.p[k], p.r.p_sem[k]);
                    obs += b;
                }
                const std::string label = (jb.op == "zz" ? "ZZ |" : jb.op == "xx" ? "XX |" : "CNOT |") + jb.in + ">";
                out("   %-17s %6d  %6d  %6.2f  %9d  %8.4f  %7.2f  %10.1e   %-66s %.5f+-%.5f\n", label.c_str(), E.circuit.n, p.r.rounds, p.r.total_ns / 1e3,
                    p.r.n_detectors, p.r.detection_rate, p.r.mean_defects, p.r.hyperedge_p, obs.c_str(), p.r.p_any, p.r.p_any_sem);
                pts.push_back(p);
            }
        }
    }
    // distance scaling
    if (dists.size() > 1) {
        out("\n3. Distance 3 -> 5: suppression of the logical error per operation (ratio d=3 / d=5)\n");
        for (const std::string chip : {"alpha1_v11", "beta1"}) {
            out("   %s physics:\n", chip == "beta1" ? "Beta-1" : "Alpha-1 v1.1");
            for (const auto& jb : jobs) {
                const Point *a = nullptr, *b = nullptr;
                for (const auto& p : pts) if (p.chip == chip && p.op == jb.op && p.input == jb.in) (p.d == 3 ? a : b) = &p;
                if (!a || !b) continue;
                std::string s;
                for (std::size_t k = 0; k < a->r.names.size(); ++k) {
                    char buf[96]; std::snprintf(buf, sizeof buf, "%s%s %.1f", k ? ", " : "", a->r.names[k].c_str(), a->r.p[k] / std::max(b->r.p[k], 1e-12));
                    s += buf;
                }
                const std::string label = (jb.op == "zz" ? "ZZ |" : jb.op == "xx" ? "XX |" : "CNOT |") + jb.in + ">";
                out("     %-12s %s; any wrong %.5f -> %.5f (%.1fx)\n", label.c_str(), s.c_str(), a->r.p_any, b->r.p_any, a->r.p_any / std::max(b->r.p_any, 1e-12));
            }
        }
    }
    const double runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    out("\nruntime %.0f s\n", runtime);
    std::fclose(g_txt);

    Json J = Json::object();
    J.set("schema", "qpuforge-surgery-study");
    J.set("provenance", review::provenance(opt, "frame Monte Carlo seed 1000 + 100 d + 10 k (+5 for Beta-1), k = 0 ZZ, 1 XX, 2 CNOT Z, 3 CNOT X; tableau check seed 77"));
    J.set("shots", double(opt.shots));
    J.set("schedule", "1 separate round, d merged rounds per merge, 1 round after the split");
    J.set("noiseless_check", jcheck);
    J.set("noiseless_check_failed", failed);
    Json jr = Json::array();
    for (const auto& p : pts) {
        Json j = Json::object();
        j.set("chip", p.chip); j.set("op", p.op); j.set("inputs", p.input); j.set("d", p.d); j.set("decoded_basis", std::string(1, p.basis));
        j.set("seed", double(p.seed)); j.set("rounds", p.r.rounds); j.set("duration_ns", p.r.total_ns);
        j.set("detectors", p.r.n_detectors); j.set("fault_mechanisms", double(p.faults));
        j.set("hyperedge_faults", p.r.hyperedge_faults); j.set("hyperedge_p", p.r.hyperedge_p); j.set("undetectable_logical", p.r.undetectable_logical);
        j.set("detection_rate", p.r.detection_rate); j.set("mean_defects", p.r.mean_defects); j.set("max_defects", p.r.max_defects);
        Json o = Json::object();
        for (std::size_t k = 0; k < p.r.names.size(); ++k) {
            Json e = Json::object(); e.set("errors", double(p.r.errors[k])); e.set("p", p.r.p[k]); e.set("sem", p.r.p_sem[k]);
            o.set(p.r.names[k], e);
        }
        j.set("observables", o);
        j.set("any_wrong", p.r.p_any); j.set("any_wrong_sem", p.r.p_any_sem); j.set("any_wrong_errors", double(p.r.any_errors));
        jr.push(j);
    }
    J.set("runs", jr);
    J.set("runtime_s", runtime);
    std::ofstream f(json_path);
    J.dump(f); f << "\n";
    std::printf("wrote %s and %s\n", txt_path.c_str(), json_path.c_str());
    return failed ? 1 : 0;
}
