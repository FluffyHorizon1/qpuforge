// Probe: the XY line's open end and bump pad sit "in the free corner between the two coupler pads". line_coupling
// extracts the XY line with truncated couplers (far_end = false: no pads). Here: a window of the assembled AF-17 mask
// around one qubit that has a +x and a -y coupler, with the chip's own cut-outs (after the sliver fill) and islands,
// extracted as is: C(XY, each coupler), C(XY, qubit), C(qubit, coupler).
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
using namespace qpuforge;
using namespace qpuforge::layout;

static bool clip(const Rect& r, const Rect& w, Rect& o) {
    o = {std::max(r.x0, w.x0), std::max(r.y0, w.y0), std::min(r.x1, w.x1), std::min(r.y1, w.y1)};
    return o.x1 > o.x0 + 1e-9 && o.y1 > o.y0 + 1e-9;
}
int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const double hmin = argc > 1 ? std::atof(argv[1]) : 3.0;
    const ChipDesign D = alpha1_v11(3);
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    // qubit with a +x and a -y neighbour coupled to it
    std::map<std::pair<long, long>, int> at;
    for (const auto& q : D.qubits) at[{std::lround(q.x_um), std::lround(q.y_um)}] = q.id;
    auto coupled = [&](int a, int b) { for (const auto& c : D.couplings) if ((c.a == a && c.b == b) || (c.a == b && c.b == a)) return true; return false; };
    int qid = -1; double cx = 0, cy = 0;
    for (const auto& q : D.qubits) {
        const long x = std::lround(q.x_um), y = std::lround(q.y_um);
        if (at.count({x + 900, y}) && at.count({x, y - 900}) && coupled(q.id, at[{x + 900, y}]) && coupled(q.id, at[{x, y - 900}])) { qid = q.id; cx = q.x_um; cy = q.y_um; break; }
    }
    const std::string cxn = "C" , tagx1 = std::to_string(qid) + "-" + std::to_string(at[{std::lround(cx + 900), std::lround(cy)}]);
    std::printf("qubit Q%d at (%.0f, %.0f)\n", qid, cx, cy);
    const Rect win{cx - 220, cy - 640, cx + 640, cy + 220};
    std::vector<Conductor> cs;
    std::vector<Rect> holes;
    for (const auto& h : L.holes) { Rect o; if (clip(h, win, o)) holes.push_back(o); }
    // shorted lines (R*, FQ*, FC*) are ground metal: carve them out of the cut-outs instead of making them conductors
    auto shorted = [](const std::string& n) { return n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0 || n[0] == 'R'; };
    for (const auto& [name, rs] : L.islands) {
        if (shorted(name)) {
            for (const auto& r : rs) {
                Rect o; if (!clip(r, win, o)) continue;
                std::vector<Rect> nh;
                for (const auto& h : holes) { auto pcs = detail::subtract(h, o); nh.insert(nh.end(), pcs.begin(), pcs.end()); }
                holes.swap(nh);
            }
            continue;
        }
        Conductor c; c.name = name;
        for (const auto& r : rs) { Rect o; if (clip(r, win, o)) c.rects.push_back(o); }
        if (!c.rects.empty()) cs.push_back(c);
    }
    std::printf("window conductors:");
    for (const auto& c : cs) std::printf(" %s", c.name.c_str());
    std::printf("\n");
    // bump pad of the XY line vs the neighbouring cut-outs: the slivers between them are filled?
    const Rect bump = rect_c(cx + kXYEnd + 32, cy - kXYLevel, 36, 36);
    std::printf("XY bump pad metal (%.1f,%.1f)-(%.1f,%.1f); ground at 4 um outside its +x and -y faces: ", bump.x0, bump.y0, bump.x1, bump.y1);
    const auto M = detail::metal_map(L.outline, L.etch);
    for (double d : {4.0, 9.0, 12.0, 20.0, 29.0}) std::printf("[+x %+.0f um: %s, -y %+.0f um: %s] ", d, M.at(bump.x1 + d, bump.cy()) >= 0 ? "metal" : "etch", d, M.at(bump.cx(), bump.y0 - d) >= 0 ? "metal" : "etch");
    std::printf("\n");
    Conductor G; G.name = "ground"; G.ground = true; G.outline = win; G.holes = holes;
    cs.push_back(G);
    const auto X = extract_capacitance(cs, 11.45, MeshOptions{hmin, 50.0, hmin < 2 ? 1.6 : 2.0});
    std::printf("extraction: %d panels, %.1f s, asymmetry %.3g\n", X.panels, X.seconds, X.asymmetry);
    const std::string Q = "Q" + std::to_string(qid), XY = "XY" + std::to_string(qid);
    for (const auto& n : X.names) if (n != XY && n != "ground") {
        const double cxy = X.mutual(XY, n);
        if (std::abs(cxy) > 0.002) std::printf("  C(%s, %s) = %.4f fF   C(%s, %s) = %.3f fF   self(%s) = %.1f fF\n", XY.c_str(), n.c_str(), cxy, Q.c_str(), n.c_str(), n == Q ? 0.0 : X.mutual(Q, n), n.c_str(), X.self(n));
    }
    std::printf("  self(%s) = %.2f fF\n", XY.c_str(), X.self(XY));
    return 0;
}
