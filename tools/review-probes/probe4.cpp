// Probe 4: capacitance of a qubit's XY line to its qubit and to the coupler pads beside it, extracted from a window of
// the assembled AF-17 mask (islands clipped to the window, all ground metal ideal ground, as in the PR's cells).
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;

int main(int argc, char** argv) {
    const int qid = argc > 1 ? std::atoi(argv[1]) : 1;
    const double hmin = argc > 2 ? std::atof(argv[2]) : 3.0;
    const ChipDesign D = alpha1_v11(3);
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    double cx = 0, cy = 0; for (const auto& q : D.qubits) if (q.id == qid) { cx = q.x_um; cy = q.y_um; }
    // window: the qubit, its XY line and the coupler pads at the (+x,-y) corner
    const Rect W{cx - 200, cy - 640, cx + 640, cy + 200};
    auto clip = [&](const Rect& r) { return Rect{std::max(r.x0, W.x0), std::max(r.y0, W.y0), std::min(r.x1, W.x1), std::min(r.y1, W.y1)}; };
    std::vector<Conductor> cs;
    std::vector<Rect> metal;
    for (const auto& [n, rs] : L.islands) {
        if (n[0] == 'R' || n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0) continue;   // shorted lines: ground metal (0 V)
        Conductor c; c.name = n;
        for (const auto& r : rs) if (r.intersects(W)) { const Rect k = clip(r); if (k.w() > 1e-6 && k.h() > 1e-6) c.rects.push_back(k); }
        if (c.rects.empty()) continue;
        metal.insert(metal.end(), c.rects.begin(), c.rects.end());
        cs.push_back(c);
    }
    Conductor G; G.name = "ground"; G.ground = true; G.outline = W;
    for (const auto& e : L.etch) if (e.intersects(W)) { const Rect k = clip(e); if (k.w() > 1e-6 && k.h() > 1e-6) G.holes.push_back(k); }
    // the etch excludes the metal: add the island rects to the holes so that islands sit inside the cut-outs
    G.holes.insert(G.holes.end(), metal.begin(), metal.end());
    cs.push_back(G);
    std::printf("window around Q%d: %zu conductors:", qid, cs.size()); for (const auto& c : cs) std::printf(" %s", c.name.c_str()); std::printf("\n");
    const auto P = mesh(cs, MeshOptions{hmin, 50.0, 2.0});
    std::printf("panels: %zu\n", P.size());
    if (P.size() > 16000) { std::printf("too many panels\n"); return 1; }
    const auto M = extract_capacitance(cs, 11.45, MeshOptions{hmin, 50.0, 2.0});
    std::printf("extraction %.0f s, asymmetry %.2e\n", M.seconds, M.asymmetry);
    const std::string xy = "XY" + std::to_string(qid), qn = "Q" + std::to_string(qid);
    std::printf("C(%s, %s) = %.4f fF; C_self(%s) = %.2f fF\n", xy.c_str(), qn.c_str(), M.mutual(xy, qn), qn.c_str(), M.self(qn));
    for (const auto& c : cs) if (!c.ground && c.name != xy) {
        const double m = M.mutual(xy, c.name);
        if (std::abs(m) > 0.005) std::printf("  C(%s, %s) = %.4f fF   [C(%s, %s) = %.3f fF, C_self(%s) = %.1f fF]\n", xy.c_str(), c.name.c_str(), m, qn.c_str(), c.name.c_str(),
                                             M.mutual(qn, c.name), c.name.c_str(), M.self(c.name));
    }
    return 0;
}
