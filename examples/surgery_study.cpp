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
//   3. distance scaling of 2
//   4. the shots of 2 decoded again with one joint decoder for all observables (matching and union-find) and with the
//      two-pass correlated decoder: per-observable error, any observable wrong, decode time
//   5. routing: the same operations between patches joined through a routing strip of length 3, 5 and 7 (Beta-1
//      physics): noiseless check, logical error, rounds, duration and qubit footprint against the strip length
//   6. chip account for two or three d = 5 patches and their routing region (qubits, couplers, feedlines, lines,
//      frequency collisions, fabrication yield)
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

struct Point {
    std::string chip, op, input; int d{0}, gap{1}; char basis{'Z'}; std::uint64_t seed{0}; SurgeryResult r; std::size_t faults{0};
    std::vector<SurgeryResult> dec;      // every decoder on the same shots (dec[0] = r)
    int qubits{0}, data{0}, couplers{0}, other_detectors{0}, joint_edges{0}, mixed_edges{0}, correlated_pairs{0};
    double minority_p{0};
};
static const std::vector<SurgeryDecoder> kDecoders = {SurgeryDecoder::Separate, SurgeryDecoder::Joint, SurgeryDecoder::SeparateUF, SurgeryDecoder::JointUF,
                                                      SurgeryDecoder::Correlated};
static std::string op_label(const std::string& op, const std::string& in) { return (op == "zz" ? "ZZ |" : op == "xx" ? "XX |" : "CNOT |") + in + ">"; }
static Json decoder_json(const SurgeryResult& q) {
    Json j = Json::object();
    j.set("decoder", std::string(surgery_decoder_name(q.decoder)));
    Json o = Json::object();
    for (std::size_t k = 0; k < q.names.size(); ++k) {
        Json e = Json::object(); e.set("errors", double(q.errors[k])); e.set("p", q.p[k]); e.set("sem", q.p_sem[k]);
        o.set(q.names[k], e);
    }
    j.set("observables", o);
    j.set("any_wrong", q.p_any); j.set("any_wrong_sem", q.p_any_sem); j.set("any_wrong_errors", double(q.any_errors));
    j.set("differs_from_separate_matching", double(q.differ)); j.set("any_wrong_fixed", double(q.any_fixed)); j.set("any_wrong_broken", double(q.any_broken));
    j.set("decode_ns_per_shot", q.decode_ns);
    return j;
}

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
    // the chip of a lattice (built and calibrated once per chip preset, distance, geometry and strip length)
    auto chip_of = [&](const std::string& chip, const std::string& op, int d, int gap, const SurgeryLattice& lat) -> const std::pair<ChipDesign, ChipReport>& {
        const std::string key = chip + "/" + std::to_string(d) + "/" + op + "/" + std::to_string(gap);
        if (!chips.count(key)) {
            auto D = surgery_chip(lat, chip);
            auto R = analyze(D);
            if (D.coupler.physical) calibrate_chip_cz(D, R, false, opt.cache);
            chips.emplace(key, std::make_pair(std::move(D), std::move(R)));
        }
        return chips.at(key);
    };
    // one point: every decoder on the same shots; r is the separate-decoder result
    auto run_point = [&](const std::string& chip, const Job& jb, int d, int gap, std::uint64_t seed, Point& p) {
        const auto E = surgery_experiment(jb.op, d, jb.in, 0, {}, gap);
        const auto& [D, R] = chip_of(chip, jb.op, d, gap, E.sched.lat);
        const auto nm = noise_model(D, R);
        const auto M = prepare_surgery(D, R, nm, E, jb.basis, true);   // throws on an observable with undetectable flips
        std::string decoded;
        for (int k : M.obs) decoded += (decoded.empty() ? "" : ",") + E.observables[std::size_t(k)].name;
        if (decoded != jb.decoded) { std::fprintf(stderr, "%s |%s>: decodes %s, expected %s\n", jb.op.c_str(), jb.in.c_str(), decoded.c_str(), jb.decoded.c_str()); std::exit(1); }
        p.chip = chip; p.op = jb.op; p.input = jb.in; p.d = d; p.gap = gap; p.basis = jb.basis; p.seed = seed;
        p.dec = run_surgery_decoders(M, opt.shots, seed, kDecoders);
        p.r = p.dec[0];
        p.faults = M.n_faults;
        p.qubits = E.circuit.n; p.couplers = int(D.couplings.size());
        for (const auto& st : E.sched.lat.sites) p.data += st.data;
        p.other_detectors = int(M.other_detectors.size()); p.joint_edges = int(M.jdem.edges.size()); p.mixed_edges = M.jdem.mixed_edges;
        p.minority_p = M.jdem.minority_p; p.correlated_pairs = M.correlated_pairs;
    };
    for (const std::string chip : {"alpha1_v11", "beta1"}) {
        for (int d : dists) {
            out("\n   %s physics, d = %d\n", chip == "beta1" ? "Beta-1" : "Alpha-1 v1.1", d);
            out("   operation          qubits  rounds  dur_us  detectors  det_rate  defects  hyperedge_p   logical error per observable                                       any wrong\n");
            for (const auto& jb : jobs) {
                Point p;
                run_point(chip, jb, d, 1, std::uint64_t(1000 + 100 * d + 10 * jb.kind + (chip == "beta1" ? 5 : 0)), p);
                std::string obs;
                for (std::size_t k = 0; k < p.r.names.size(); ++k) {
                    char b[96]; std::snprintf(b, sizeof b, "%s%s %.5f+-%.5f", k ? ", " : "", p.r.names[k].c_str(), p.r.p[k], p.r.p_sem[k]);
                    obs += b;
                }
                const std::string label = op_label(jb.op, jb.in);
                out("   %-17s %6d  %6d  %6.2f  %9d  %8.4f  %7.2f  %10.1e   %-66s %.5f+-%.5f\n", label.c_str(), p.qubits, p.r.rounds, p.r.total_ns / 1e3,
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
    // ---------------------------------------------------------------- 4. decoders on the same shots
    out("\n4. Decoders on the shots of section 2 (same syndromes for every decoder)\n");
    out("   separate matching    one exact matching per observable, each on its own graph (the decoder of section 2)\n");
    out("   joint matching       one exact matching on one graph whose edges carry a mask of the observables they flip; edge probability from all mechanisms\n");
    out("                        of its detector pair, mask of the most likely class. Mechanisms flipping more than two detectors stay out of the graph (hyperedge_p above)\n");
    out("   union-find           the weighted union-find decoder on the same graphs, separate and joint\n");
    out("   correlated           two passes: union-find on the detectors of the other Pauli type, then joint matching with the weights of the edges that share\n");
    out("                        mechanisms with that correction lowered to their conditional probability\n");
    out("   differs / fixed / broken: shots whose prediction differs from separate matching; shots with some observable wrong under separate matching and none\n");
    out("   under this decoder, and the reverse. decode: mean wall time per shot on one thread (informational, not reproducible)\n");
    long joint_differs = 0, shots_total = 0;
    double t_ratio_lo = 1e300, t_ratio_hi = 0, corr_lo = 1e300, corr_hi = 0, uf_lo = 1e300, uf_hi = 0;
    for (const auto& p : pts) {
        out("\n   %s physics, d = %d, %s: %d %c-type detectors + %d of the other type; joint graph %d edges, %d with mixed masks (minority p %.1e); %d correlated edge pairs\n",
            p.chip == "beta1" ? "Beta-1" : "Alpha-1 v1.1", p.d, op_label(p.op, p.input).c_str(), p.r.n_detectors, p.basis, p.other_detectors, p.joint_edges, p.mixed_edges, p.minority_p,
            p.correlated_pairs);
        out("     decoder                     logical error per observable                         any wrong          differs   fixed  broken  decode_us\n");
        for (const auto& q : p.dec) {
            std::string obs;
            for (std::size_t k = 0; k < q.names.size(); ++k) { char b[64]; std::snprintf(b, sizeof b, "%s%s %.5f", k ? ", " : "", q.names[k].c_str(), q.p[k]); obs += b; }
            out("     %-26s  %-51s  %.5f+-%.5f  %7ld  %6ld  %6ld  %9.2f\n", surgery_decoder_name(q.decoder), obs.c_str(), q.p_any, q.p_any_sem, q.differ, q.any_fixed, q.any_broken, q.decode_ns / 1e3);
        }
        joint_differs += p.dec[1].differ; shots_total += p.r.shots;
        const double tr = p.dec[0].decode_ns / std::max(p.dec[1].decode_ns, 1e-9), cr = p.dec[4].p_any / std::max(p.dec[1].p_any, 1e-12), ur = p.dec[3].p_any / std::max(p.dec[1].p_any, 1e-12);
        t_ratio_lo = std::min(t_ratio_lo, tr); t_ratio_hi = std::max(t_ratio_hi, tr);
        corr_lo = std::min(corr_lo, cr); corr_hi = std::max(corr_hi, cr); uf_lo = std::min(uf_lo, ur); uf_hi = std::max(uf_hi, ur);
    }
    out("\n   summary over the %zu points (%ld shots): joint matching differs from separate matching in %ld shots; separate / joint decode time %.1f-%.1fx;\n", pts.size(), shots_total,
        joint_differs, t_ratio_lo, t_ratio_hi);
    out("   any wrong, joint union-find / joint matching %.2f-%.2f; correlated / joint matching %.2f-%.2f\n", uf_lo, uf_hi, corr_lo, corr_hi);

    // ---------------------------------------------------------------- 5. routing strips
    const std::vector<int> gaps = opt.smoke ? std::vector<int>{1, 3} : std::vector<int>{1, 3, 5, 7};
    out("\n5. Routing: operations between patches joined through a routing strip (d wide, `strip` data rows or columns long; strip 1 = adjacent patches)\n");
    out("   the strip is merged with both patches for the d rounds of the measurement and read out at the split; the CNOT has one strip between control and ancilla\n");
    out("   and one between ancilla and target\n");
    out("\n5a. Noiseless check (stabiliser tableau, 4 random measurement branches per experiment), strips longer than 1\n");
    out("   op     d  strip  qubits  experiments  failed   inputs\n");
    Json jrcheck = Json::array();
    int rchecked = 0, rfailed = 0;
    {
        struct In { const char* in; char ro; };
        const std::vector<In> zz_in = {{"00", 0}, {"01", 0}, {"11", 0}, {"++", 0}, {"+-", 0}}, xx_in = {{"++", 0}, {"+-", 0}, {"--", 0}, {"00", 0}, {"01", 0}},
                              cn_in = {{"00", 0}, {"10", 0}, {"11", 0}, {"++", 0}, {"-+", 0}, {"--", 0}, {"+0", 'Z'}, {"+0", 'X'}};
        for (int d : dists)
            for (int gap : gaps) {
                if (gap == 1) continue;
                for (const std::string op : {"zz", "xx", "cnot"}) {
                    const auto& ins = op == "zz" ? zz_in : op == "xx" ? xx_in : cn_in;
                    int n = 0, f = 0, nq = 0; std::string list;
                    for (const auto& in : ins) {
                        const auto E = surgery_experiment(op, d, in.in, in.ro, {}, gap);
                        const auto c = check_surgery_noiseless(E, 4, 77);
                        const bool ok = !(c.fired || c.wrong_observables || c.wrong_reference || c.random_reference);
                        ++n; f += !ok; nq = E.circuit.n;
                        list += std::string(list.empty() ? "" : " ") + in.in + (in.ro ? std::string("(") + in.ro + ")" : "");
                        Json j = Json::object();
                        j.set("op", op); j.set("d", d); j.set("strip", gap); j.set("inputs", std::string(in.in)); j.set("readout", in.ro ? std::string(1, in.ro) : std::string(""));
                        j.set("detectors", int(E.detectors.size())); j.set("firings", c.fired); j.set("wrong_observables", c.wrong_observables);
                        j.set("wrong_reference_bits", c.wrong_reference + c.random_reference); j.set("ok", ok);
                        jrcheck.push(j);
                    }
                    rchecked += n; rfailed += f;
                    out("   %-5s %2d  %5d  %6d  %11d  %6d   %s\n", op.c_str(), d, gap, nq, n, f, list.c_str());
                }
            }
    }
    out("   %d experiments, %d failed\n", rchecked, rfailed);
    failed += rfailed;

    out("\n5b. Logical error per operation against the strip length, Beta-1 physics, %ld shots per point (strip 1: the points of section 2)\n", opt.shots);
    out("   errors with joint matching; the last column is any observable wrong with the correlated decoder on the same shots\n");
    std::vector<Point> rpts;
    for (int d : dists) {
        out("\n   d = %d\n", d);
        out("   operation     strip  qubits  (data)  couplers  rounds  dur_us  detectors  defects   logical error per observable                                       any wrong          correlated\n");
        for (const auto& jb : jobs)
            for (int gap : gaps) {
                Point p;
                if (gap == 1) { for (const auto& q : pts) if (q.chip == "beta1" && q.d == d && q.op == jb.op && q.input == jb.in) p = q; }
                else run_point("beta1", jb, d, gap, std::uint64_t(3000 + 100 * d + 10 * jb.kind + gap), p);
                const auto& q = p.dec[1];
                std::string obs;
                for (std::size_t k = 0; k < q.names.size(); ++k) { char b[96]; std::snprintf(b, sizeof b, "%s%s %.5f+-%.5f", k ? ", " : "", q.names[k].c_str(), q.p[k], q.p_sem[k]); obs += b; }
                out("   %-12s  %5d  %6d  %6d  %8d  %6d  %6.2f  %9d  %7.2f   %-66s %.5f+-%.5f  %.5f+-%.5f\n", op_label(jb.op, jb.in).c_str(), gap, p.qubits, p.data, p.couplers, q.rounds,
                    q.total_ns / 1e3, q.n_detectors, q.mean_defects, obs.c_str(), q.p_any, q.p_any_sem, p.dec[4].p_any, p.dec[4].p_any_sem);
                rpts.push_back(p);
            }
    }
    out("\n   any observable wrong, relative to adjacent patches (joint matching):\n");
    for (int d : dists)
        for (const auto& jb : jobs) {
            std::string line;
            double base = 0;
            for (const auto& p : rpts) if (p.d == d && p.op == jb.op && p.input == jb.in) {
                if (p.gap == 1) base = p.dec[1].p_any;
                char b[64]; std::snprintf(b, sizeof b, "  strip %d: %.2fx", p.gap, p.dec[1].p_any / std::max(base, 1e-12)); line += b;
            }
            out("     d = %d, %-12s%s\n", d, op_label(jb.op, jb.in).c_str(), line.c_str());
        }

    // ---------------------------------------------------------------- 6. chip account
    const int dc = dists.back();
    out("\n6. Chip account for d = %d patches and their routing region (surgery_chip with Beta-1 physics; no layout, no mask)\n", dc);
    out("   lines: one XY and one flux line per qubit, one flux line per coupler, two per feedline (as the design reviews count them); feedlines: at most 9 resonators each\n");
    out("   frequencies: distinct design frequencies of the qubits (the 8-frequency plan: 4 data + 4 measure colours); collisions: pairs flagged by the collision rules;\n");
    out("   CZ combinations: distinct (measure, data) frequency pairs among the couplers\n");
    out("   lattice                               qubits  (data  measure)  couplers  feedlines  lines  frequencies  collisions  CZ combinations  max |ZZ| kHz\n");
    Json jchip = Json::array();
    struct Lay { std::string name, op; int gap; };
    std::vector<Lay> lays;
    for (int gap : gaps) lays.push_back({"2 patches, ZZ" + std::string(gap == 1 ? " (adjacent)" : ", strip " + std::to_string(gap)), "zz", gap});
    for (int gap : gaps) lays.push_back({"2 patches, XX" + std::string(gap == 1 ? " (adjacent)" : ", strip " + std::to_string(gap)), "xx", gap});
    for (int gap : gaps) lays.push_back({"3 patches, CNOT" + std::string(gap == 1 ? " (adjacent)" : ", strips " + std::to_string(gap)), "cnot", gap});
    auto account = [&](const std::string& name, const ChipDesign& D, const ChipReport& R, int strip) {
        int nd = 0; for (const auto& q : D.qubits) nd += q.role == "data";
        std::set<int> feeds, freqs; std::set<std::pair<int, int>> combos;
        for (const auto& q : D.qubits) feeds.insert(q.feedline);
        for (const auto& q : R.q) freqs.insert(int(std::lround(q.f01 * 1e3)));
        for (const auto& pr : R.pairs) if (pr.designed) {
            const int a = D.index_of(pr.a), b = D.index_of(pr.b);
            int fa = int(std::lround(R.q[std::size_t(a)].f01 * 1e3)), fb = int(std::lround(R.q[std::size_t(b)].f01 * 1e3));
            if (fa < fb) std::swap(fa, fb);
            combos.insert({fa, fb});
        }
        const int n = int(D.qubits.size()), nc = int(D.couplings.size()), lines = 2 * n + nc + 2 * int(feeds.size());
        out("   %-36s  %6d  %5d  %7d   %8d  %9zu  %5d  %11zu  %10d  %15zu  %12.1f\n", name.c_str(), n, nd, n - nd, nc, feeds.size(), lines, freqs.size(), R.n_collisions, combos.size(),
            R.max_abs_zz_kHz);
        Json j = Json::object();
        j.set("lattice", name); j.set("d", dc); j.set("strip", strip); j.set("qubits", n); j.set("data", nd); j.set("measure", n - nd); j.set("couplers", nc);
        j.set("feedlines", int(feeds.size())); j.set("lines", lines); j.set("design_frequencies", int(freqs.size())); j.set("collisions", R.n_collisions);
        j.set("cz_frequency_combinations", int(combos.size())); j.set("max_abs_zz_kHz", R.max_abs_zz_kHz);
        return j;
    };
    {
        const auto G = beta1(dc);
        jchip.push(account("1 patch (Beta-1, for reference)", G, analyze(G), 0));
        if (!opt.smoke) { const auto G7 = gamma1(); jchip.push(account("Gamma-1 (d = 7 memory, reference)", G7, analyze(G7), 0)); }
    }
    for (const auto& l : lays) {
        const auto E = surgery_experiment(l.op, dc, l.op == "xx" ? "++" : "00", 0, {}, l.gap);
        const auto& [D, R] = chip_of("beta1", l.op, dc, l.gap, E.sched.lat);
        jchip.push(account(l.name, D, R, l.gap));
    }
    out("\n   fabrication Monte Carlo (junction spread 0.5 %% after trimming and 2 %% as fabricated; %d samples, seed 7; static collision rules):\n", opt.smoke ? 10 : 100);
    Json jfab = Json::array();
    for (const auto& l : lays) {
        if (!((l.op == "zz" && l.gap <= 3) || (l.op == "cnot" && l.gap == 1))) continue;
        const auto E = surgery_experiment(l.op, dc, "00", 0, {}, l.gap);
        const auto& [D0, R0] = chip_of("beta1", l.op, dc, l.gap, E.sched.lat);
        (void)R0;
        for (double spread : {0.5, 2.0}) {
            ChipDesign D = D0; D.Rn_spread_pct = spread;
            const auto F = fab_monte_carlo(D, opt.smoke ? 10 : 100, 7);
            double fs = 0; for (double v : F.f01_sigma_MHz) fs += v;
            out("   %-36s  spread %.1f %%: collision-free %3d/%d, frequency scatter %.0f MHz rms\n", l.name.c_str(), spread, F.collision_free, F.samples, fs / std::max(F.samples, 1));
            Json j = Json::object();
            j.set("lattice", l.name); j.set("d", dc); j.set("strip", l.gap); j.set("spread_pct", spread); j.set("samples", F.samples); j.set("collision_free", F.collision_free);
            j.set("f01_scatter_MHz", fs / std::max(F.samples, 1));
            jfab.push(j);
        }
    }
    const double runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    out("\nruntime %.0f s\n", runtime);
    std::fclose(g_txt);

    Json J = Json::object();
    J.set("schema", "qpuforge-surgery-study");
    J.set("provenance", review::provenance(opt, "frame Monte Carlo seed 1000 + 100 d + 10 k (+5 for Beta-1), k = 0 ZZ, 1 XX, 2 CNOT Z, 3 CNOT X; tableau check seed 77; routing strips: 3000 + 100 d + 10 k + strip length; fabrication Monte Carlo seed 7"));
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
        j.set("other_type_detectors", p.other_detectors); j.set("joint_edges", p.joint_edges); j.set("joint_mixed_edges", p.mixed_edges);
        j.set("joint_minority_p", p.minority_p); j.set("correlated_edge_pairs", p.correlated_pairs);
        Json jd = Json::array();
        for (const auto& q : p.dec) jd.push(decoder_json(q));
        j.set("decoders", jd);
        jr.push(j);
    }
    J.set("runs", jr);
    Json jrt = Json::array();
    for (const auto& p : rpts) {
        Json j = Json::object();
        j.set("chip", p.chip); j.set("op", p.op); j.set("inputs", p.input); j.set("d", p.d); j.set("strip", p.gap); j.set("decoded_basis", std::string(1, p.basis));
        j.set("seed", double(p.seed)); j.set("qubits", p.qubits); j.set("data_qubits", p.data); j.set("couplers", p.couplers);
        j.set("rounds", p.r.rounds); j.set("duration_ns", p.r.total_ns); j.set("detectors", p.r.n_detectors); j.set("fault_mechanisms", double(p.faults));
        j.set("hyperedge_faults", p.dec[1].hyperedge_faults); j.set("hyperedge_p", p.dec[1].hyperedge_p);
        j.set("detection_rate", p.r.detection_rate); j.set("mean_defects", p.r.mean_defects); j.set("max_defects", p.r.max_defects);
        Json jd = Json::array();
        for (const auto& q : p.dec) jd.push(decoder_json(q));
        j.set("decoders", jd);
        jrt.push(j);
    }
    J.set("routing", jrt);
    J.set("routing_noiseless_check", jrcheck);
    J.set("routing_noiseless_check_failed", rfailed);
    J.set("chip_account", jchip);
    J.set("chip_account_fab_yield", jfab);
    J.set("runtime_s", runtime);
    std::ofstream f(json_path);
    J.dump(f); f << "\n";
    std::printf("wrote %s and %s\n", txt_path.c_str(), json_path.c_str());
    return failed ? 1 : 0;
}
