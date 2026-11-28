// Independent topology check: for each qubit, the base-metal component on the ground side of its SQUID region
// (from L.jj), the component holding its flux line FQ<id>, every ground strap with a foot on that component, the
// component its other foot lands on, and the strap depth of the pocket from the chip-border plane.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <queue>
using namespace qpuforge;
using namespace qpuforge::layout;

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "AF-17";
    ChipDesign D; std::string cells;
    if (which == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    std::printf("%s: %s\n", which.c_str(), connectivity_summary(L).c_str());
    const auto M = detail::metal_map(L.outline, L.etch);
    const int plane = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    std::vector<double> area(std::size_t(M.n_comp), 0);
    for (std::size_t i = 0; i < M.comp.size(); ++i) area[std::size_t(M.comp[i])] += (M.ix1[i] - M.ix0[i]) * (M.iy1[i] - M.iy0[i]);
    // strap graph over components
    std::map<int, std::vector<int>> on;
    std::vector<std::pair<int, int>> ends;
    for (int i = 0; i < int(L.straps.size()); ++i) {
        const auto& s = L.straps[std::size_t(i)];
        const int a = M.at(s.foot_a.cx(), s.foot_a.cy()), b = M.at(s.foot_b.cx(), s.foot_b.cy());
        ends.push_back({a, b}); on[a].push_back(i); on[b].push_back(i);
    }
    std::map<int, int> depth; std::queue<int> qu; depth[plane] = 0; qu.push(plane);
    while (!qu.empty()) { const int u = qu.front(); qu.pop(); for (int i : on[u]) { const int v = ends[std::size_t(i)].first == u ? ends[std::size_t(i)].second : ends[std::size_t(i)].first; if (!depth.count(v)) { depth[v] = depth[u] + 1; qu.push(v); } } }
    std::map<int, int> nstraps_hist; int n_ok = 0;
    for (const auto& qb : D.qubits) {
        const double cx = qb.x_um, cy = qb.y_um;
        // the qubit SQUID region: the jj rect nearest (cx-29, cy-21)
        const Rect* jj = nullptr; double best = 1e300;
        for (const auto& j : L.jj) { const double d = std::hypot(j.cx() - (cx - 29), j.cy() - (cy - 21)); if (d < best) { best = d; jj = &j; } }
        const int gside = M.at(jj->x0 - 0.5, jj->cy()), iside = M.at(jj->x1 + 0.5, jj->cy());
        const int island = M.at(cx, cy);
        int fq = -2;
        for (const auto& [name, rs] : L.islands) if (name == "FQ" + std::to_string(qb.id)) {
            const Rect* big = &rs.front(); for (const auto& r : rs) if (r.area() > big->area()) big = &r;
            fq = M.at(big->cx(), big->cy());
        }
        const int short_sq = M.at(cx - 61, cy - 61);
        const int cpw_l = M.at(cx - 156 - 14, cy - 200), cpw_r = M.at(cx - 156 + 14, cy - 200);
        std::printf("Q%-3d jj rel (%.0f,%.0f)-(%.0f,%.0f): island side comp %d (qubit comp %d), ground side comp %d (%s, %.0f um^2); "
                    "FQ comp %d, short-square comp %d; straps with a foot on it: %zu; strap depth to plane %d; flux CPW ground L/R %d/%d\n",
                    qb.id, jj->x0 - cx, jj->y0 - cy, jj->x1 - cx, jj->y1 - cy, iside, island, gside, gside == plane ? "PLANE" : "not plane",
                    gside >= 0 ? area[std::size_t(gside)] : -1.0, fq, short_sq, on[gside].size(), depth.count(gside) ? depth[gside] : -1, cpw_l, cpw_r);
        for (int i : on[gside]) {
            const auto& s = L.straps[std::size_t(i)];
            const int other = ends[std::size_t(i)].first == gside ? ends[std::size_t(i)].second : ends[std::size_t(i)].first;
            const Rect& fp = ends[std::size_t(i)].first == gside ? s.foot_a : s.foot_b;
            const Rect& fo = ends[std::size_t(i)].first == gside ? s.foot_b : s.foot_a;
            std::printf("      strap span rel (%.0f,%.0f)-(%.0f,%.0f) [%.0f um long], pocket foot centre rel (%.0f,%.0f), far foot centre rel (%.0f,%.0f) on comp %d%s%s\n",
                        s.span.x0 - cx, s.span.y0 - cy, s.span.x1 - cx, s.span.y1 - cy, std::max(s.span.w(), s.span.h()),
                        fp.cx() - cx, fp.cy() - cy, fo.cx() - cx, fo.cy() - cy, other, other == plane ? " (plane)" : "",
                        (other == cpw_l || other == cpw_r) ? " (= flux-line CPW ground)" : "");
        }
        nstraps_hist[int(on[gside].size())]++;
        if (gside != plane && gside == fq && gside == short_sq && iside == island) ++n_ok;
    }
    std::printf("qubits whose SQUID ground side is the flux-line pocket (not the plane): %d of %zu; straps-per-pocket histogram:", n_ok, D.qubits.size());
    for (auto& [k, v] : nstraps_hist) std::printf(" %d:%d", k, v);
    std::printf("\n");
    return 0;
}
