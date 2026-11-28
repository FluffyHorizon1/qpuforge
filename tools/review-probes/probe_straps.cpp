// Probe: ground straps of a committed chip -- span lengths, what they cross, which strap(s) ground each qubit's
// (-x,-y) pocket (flux-line short + SQUID tongue), and the strap chain from that pocket to the chip-border plane.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <queue>
#include <set>
using namespace qpuforge;
using namespace qpuforge::layout;
using detail::kStrapW;

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "AF-17";
    ChipDesign D; std::string cells;
    if (which == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else if (which == "BF-49") { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    else { D = gamma1(7); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    std::printf("%s: %s\n  %s\n", which.c_str(), drc_summary(L).c_str(), connectivity_summary(L).c_str());
    const auto M = detail::metal_map(L.outline, L.etch);
    const int plane = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    // strap graph
    struct E { int a, b, idx; };
    std::vector<E> es;
    std::map<int, std::vector<int>> adj;
    double maxlen = 0; std::map<int, int> hist;
    int over_bump = 0, over_jj = 0, crossing_straps = 0, feet_on_strap_span = 0;
    for (int i = 0; i < int(L.straps.size()); ++i) {
        const auto& s = L.straps[std::size_t(i)];
        const int a = M.at(s.foot_a.cx(), s.foot_a.cy()), b = M.at(s.foot_b.cx(), s.foot_b.cy());
        es.push_back({a, b, i}); adj[a].push_back(i); adj[b].push_back(i);
        const double len = std::max(s.span.w(), s.span.h()) - 2 * kStrapW;    // free span between the feet
        maxlen = std::max(maxlen, len); hist[int(len / 25) * 25]++;
        for (const auto& bp : L.bumps) if (s.span.intersects(bp)) { ++over_bump; std::printf("  strap %d (%.0f,%.0f)-(%.0f,%.0f) spans bump pad at (%.0f,%.0f)\n", i, s.span.x0, s.span.y0, s.span.x1, s.span.y1, bp.cx(), bp.cy()); break; }
        for (const auto& j : L.jj) if (s.span.intersects(j)) { ++over_jj; std::printf("  strap %d (%.0f,%.0f)-(%.0f,%.0f) spans junction region at (%.0f,%.0f)\n", i, s.span.x0, s.span.y0, s.span.x1, s.span.y1, j.cx(), j.cy()); break; }
        for (int k = i + 1; k < int(L.straps.size()); ++k) if (s.span.intersects(L.straps[std::size_t(k)].span)) ++crossing_straps;
    }
    std::printf("  straps %zu, longest free span %.0f um; free-span histogram (25 um bins):", L.straps.size(), maxlen);
    for (auto& [k, v] : hist) std::printf(" [%d,%d):%d", k, k + 25, v);
    std::printf("\n  straps spanning a bump pad: %d, spanning a junction region: %d, strap pairs overlapping: %d\n", over_bump, over_jj, crossing_straps);
    // BFS from the plane over straps: depth (number of straps in series) and summed span per component
    std::map<int, std::pair<int, double>> dist; std::queue<int> q; dist[plane] = {0, 0}; q.push(plane);
    while (!q.empty()) {
        const int u = q.front(); q.pop();
        for (int ei : adj[u]) {
            const auto& e = es[std::size_t(ei)]; const int v = e.a == u ? e.b : e.a;
            if (dist.count(v)) continue;
            const auto& s = L.straps[std::size_t(e.idx)];
            dist[v] = {dist[u].first + 1, dist[u].second + std::max(s.span.w(), s.span.h()) - 2 * kStrapW};
            q.push(v);
        }
    }
    // each qubit's (-x,-y) pocket: the component of the flux line's shorted last square
    int nq = 0; double worst_len = 0; int worst_depth = 0; std::map<int, int> depth_hist;
    for (const auto& qb : D.qubits) {
        const XmonParams& X = qb.role != "data" ? C.meas : C.data;
        const double cx = qb.x_um, cy = qb.y_um, hw = X.w_arm / 2, s0 = -hw - X.squid_stub;
        const int pocket = M.at(cx - 61, cy - 61);
        const int tongue = M.at(cx + s0 - kSquidGap - kTongueW / 2, cy - hw - kSquidGap - 2);
        const int island = M.at(cx, cy);
        // CPW ground beside the flux line's vertical run (x = cx-156, y in [cy-250, cy-156]): 11..20 um from the centre line
        const int fl_gnd_l = M.at(cx - 156 - 14, cy - 200), fl_gnd_r = M.at(cx - 156 + 14, cy - 200);
        // ground beside the XY line's open end
        const int xy_gnd = M.at(cx + kXYEnd + 10, cy - kXYLevel + 14);
        int direct = 0; for (int ei : adj[pocket]) { (void)ei; ++direct; }
        const auto dp = dist.count(pocket) ? dist[pocket] : std::make_pair(-1, -1.0);
        const auto dg = dist.count(fl_gnd_l) ? dist[fl_gnd_l] : std::make_pair(-1, -1.0);
        const auto dx = dist.count(xy_gnd) ? dist[xy_gnd] : std::make_pair(-1, -1.0);
        std::printf("  Q%-3d pocket comp %d (tongue %s, island %d): %d strap(s) on it; to the plane %d strap(s) in series, %.0f um of span | "
                    "flux-line CPW ground L/R comps %d/%d (%s), plane via %d strap(s); XY-end ground comp %d, plane via %d strap(s)\n",
                    qb.id, pocket, tongue == pocket ? "same" : "DIFFERENT", island, direct, dp.first, dp.second, fl_gnd_l, fl_gnd_r,
                    fl_gnd_l == pocket ? "= pocket" : "not the pocket", dg.first, xy_gnd, dx.first);
        for (int ei : adj[pocket]) {
            const auto& s = L.straps[std::size_t(es[std::size_t(ei)].idx)];
            const int other = es[std::size_t(ei)].a == pocket ? es[std::size_t(ei)].b : es[std::size_t(ei)].a;
            std::printf("        strap (%.0f,%.0f)-(%.0f,%.0f) rel. qubit (%.0f,%.0f)-(%.0f,%.0f), free span %.0f um, crosses %d line(s) %d coupler/bus; other end comp %d%s\n",
                        s.span.x0, s.span.y0, s.span.x1, s.span.y1, s.span.x0 - cx, s.span.y0 - cy, s.span.x1 - cx, s.span.y1 - cy,
                        std::max(s.span.w(), s.span.h()) - 2 * kStrapW, s.lines_crossed, s.couplers_crossed, other,
                        other == plane ? " (plane)" : (other == fl_gnd_l ? " (flux-line CPW ground)" : ""));
        }
        ++nq; worst_len = std::max(worst_len, dp.second); worst_depth = std::max(worst_depth, dp.first); depth_hist[dp.first]++;
    }
    std::printf("  pockets: %d; straps in series to the plane:", nq);
    for (auto& [k, v] : depth_hist) std::printf(" %d:%d", k, v);
    std::printf("; worst summed span %.0f um\n", worst_len);
    return 0;
}
