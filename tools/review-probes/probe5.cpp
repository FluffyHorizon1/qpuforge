// Probe 5: XY drive of the target qubit and of the neighbour qubit through the coupler whose pad is 30 um from the XY pad.
// Window of the assembled AF-17 mask with the target qubit, its XY line, the adjacent coupler(s) and their other qubits.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;
static double rdist(const Rect& a, const Rect& b) {
    const double dx = std::max({0.0, b.x0 - a.x1, a.x0 - b.x1}), dy = std::max({0.0, b.y0 - a.y1, a.y0 - b.y1});
    return std::hypot(dx, dy);
}
int main(int argc, char** argv) {
    const int qid = argc > 1 ? std::atoi(argv[1]) : 1;
    const double hmin = argc > 2 ? std::atof(argv[2]) : 3.0;
    const ChipDesign D = alpha1_v11(3);
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    const auto R = analyze(D);
    auto isl = [&](const std::string& n) -> const std::vector<Rect>& { for (const auto& [m, rs] : L.islands) if (m == n) return rs; throw std::runtime_error(n); };
    const std::string xy = "XY" + std::to_string(qid), qn = "Q" + std::to_string(qid);
    std::vector<std::string> keep{qn, xy};
    std::vector<std::pair<std::string, std::string>> cpl;   // coupler, other qubit
    for (const auto& [n, rs] : L.islands) if (n[0] == 'C') {
        double d = 1e9; for (const auto& a : rs) for (const auto& b : isl(xy)) d = std::min(d, rdist(a, b));
        if (d > 100) continue;
        const auto dash = n.find('-'); const int a = std::stoi(n.substr(1, dash - 1)), b = std::stoi(n.substr(dash + 1));
        const int other = a == qid ? b : a;
        cpl.push_back({n, "Q" + std::to_string(other)});
        keep.push_back(n); keep.push_back("Q" + std::to_string(other));
        std::printf("%s is %.1f um from %s; its other qubit is Q%d\n", n.c_str(), d, xy.c_str(), other);
    }
    Rect W{1e300, 1e300, -1e300, -1e300};
    for (const auto& k : keep) { const Rect b = bbox_of(isl(k)); W = {std::min(W.x0, b.x0), std::min(W.y0, b.y0), std::max(W.x1, b.x1), std::max(W.y1, b.y1)}; }
    W = W.grow(80);
    W = {std::max(W.x0, L.outline.x0), std::max(W.y0, L.outline.y0), std::min(W.x1, L.outline.x1), std::min(W.y1, L.outline.y1)};
    auto clip = [&](const Rect& r) { return Rect{std::max(r.x0, W.x0), std::max(r.y0, W.y0), std::min(r.x1, W.x1), std::min(r.y1, W.y1)}; };
    std::vector<Conductor> cs; std::vector<Rect> metal;
    for (const auto& [n, rs] : L.islands) {
        if (n[0] == 'R' || n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0) continue;     // shorted lines: ground metal (0 V)
        Conductor c; c.name = n;
        for (const auto& r : rs) if (r.intersects(W)) { const Rect k = clip(r); if (k.w() > 1e-6 && k.h() > 1e-6) c.rects.push_back(k); }
        if (c.rects.empty()) continue;
        metal.insert(metal.end(), c.rects.begin(), c.rects.end()); cs.push_back(c);
    }
    Conductor G; G.name = "ground"; G.ground = true; G.outline = W;
    for (const auto& e : L.etch) if (e.intersects(W)) { const Rect k = clip(e); if (k.w() > 1e-6 && k.h() > 1e-6) G.holes.push_back(k); }
    G.holes.insert(G.holes.end(), metal.begin(), metal.end());
    cs.push_back(G);
    std::printf("window x[%.0f,%.0f] y[%.0f,%.0f]; conductors:", W.x0, W.x1, W.y0, W.y1); for (const auto& c : cs) std::printf(" %s", c.name.c_str()); std::printf("\n");
    const auto P = mesh(cs, MeshOptions{hmin, 50.0, 2.0});
    std::printf("panels %zu\n", P.size()); std::fflush(stdout);
    if (P.size() > 14000) return 1;
    const auto M = extract_capacitance(cs, 11.45, MeshOptions{hmin, 50.0, 2.0});
    std::printf("extraction %.0f s, asymmetry %.2e\n", M.seconds, M.asymmetry);
    const double Cxq = M.mutual(xy, qn);
    std::printf("C(%s,%s) %.4f fF, C_sigma(%s) %.1f fF\n", xy.c_str(), qn.c_str(), Cxq, qn.c_str(), M.self(qn));
    double fq = 0; for (std::size_t i = 0; i < D.qubits.size(); ++i) if (D.qubits[i].id == qid) fq = R.q[i].f01;
    for (const auto& [cn, on] : cpl) {
        const double Cxc = M.mutual(xy, cn), Cqc = M.mutual(qn, cn), Coc = M.mutual(on, cn), Ccc = M.self(cn), Cxo = M.mutual(xy, on);
        double fo = 0; for (std::size_t i = 0; i < D.qubits.size(); ++i) if ("Q" + std::to_string(D.qubits[i].id) == on) fo = R.q[i].f01;
        std::printf("%s: C(xy,coupler) %.3f fF, C(target,coupler) %.3f fF, C(%s,coupler) %.3f fF, C_coupler %.1f fF, C(xy,%s) direct %.5f fF; f01 target %.3f GHz, %s %.3f GHz\n",
                    cn.c_str(), Cxc, Cqc, on.c_str(), Coc, Ccc, on.c_str(), Cxo, fq, on.c_str(), fo);
        const double fc = 8.398;   // coupler park frequency (alpha1_v11_review.txt)
        for (double fd : {fq, 5.11, 5.5, 5.99}) {
            const double k = 1.0 / (1.0 - (fc * fc) / (fd * fd));        // coupler node response C_xc / (C_cc (1 - w_c^2/w^2))
            const double Ceff_t = Cxq + Cqc * Cxc / Ccc * k, Ceff_o = Cxo + Coc * Cxc / Ccc * k;
            std::printf("   drive at %.3f GHz: effective XY coupling target %+.4f fF (direct %+.4f, via coupler %+.4f), neighbour %s %+.4f fF -> crosstalk %.0f%% (%.1f dB)\n",
                        fd, Ceff_t, Cxq, Cqc * Cxc / Ccc * k, on.c_str(), Ceff_o, 100 * std::abs(Ceff_o / Ceff_t), 20 * std::log10(std::abs(Ceff_o / Ceff_t)));
        }
        // coupler T1 limit through the XY line (matched far end), at the park frequency
        const double wc = 2 * PI * fc * 1e9;
        std::printf("   coupler T1 limit through %s: %.2f us\n", xy.c_str(), Ccc * 1e-15 / (wc * wc * Cxc * 1e-15 * Cxc * 1e-15 * 50) * 1e6);
    }
    return 0;
}
