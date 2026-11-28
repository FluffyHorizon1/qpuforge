// Independent check (refutation attempt): for every qubit, the metal component on the ground side of its SQUID region
// (layer 2), the component of its flux line FQ*, the straps with a foot on that component, and an ASCII map of the corner.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
using namespace qpuforge;
using namespace qpuforge::layout;

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "AF-17";
    const int map_q = argc > 2 ? std::atoi(argv[2]) : -1;
    ChipDesign D; std::string cells;
    if (which == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else if (which == "BF-49") { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    else { D = gamma1(7); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    std::printf("%s: %s\n", which.c_str(), connectivity_summary(L).c_str());
    const auto M = detail::metal_map(L.outline, L.etch);
    const int plane = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    // component area
    std::vector<double> area(std::size_t(M.n_comp), 0);
    for (std::size_t i = 0; i < M.comp.size(); ++i) area[std::size_t(M.comp[i])] += (M.ix1[i] - M.ix0[i]) * (M.iy1[i] - M.iy0[i]);
    // straps per component
    std::map<int, std::vector<int>> on;
    for (int i = 0; i < int(L.straps.size()); ++i) {
        const auto& s = L.straps[std::size_t(i)];
        const int a = M.at(s.foot_a.cx(), s.foot_a.cy()), b = M.at(s.foot_b.cx(), s.foot_b.cy());
        on[a].push_back(i); if (b != a) on[b].push_back(i);
    }
    // island name -> component of its largest rect
    std::map<std::string, int> icomp;
    for (const auto& [n, rs] : L.islands) { const Rect* big = &rs.front(); for (const auto& r : rs) if (r.area() > big->area()) big = &r; icomp[n] = M.at(big->cx(), big->cy()); }
    std::map<int, int> nstraps_hist; int same_as_flux = 0, n = 0;
    for (const auto& q : D.qubits) {
        const double cx = q.x_um, cy = q.y_um;
        // the qubit's junction: the layer-2 rect closest to the qubit centre
        const Rect* jq = nullptr; double best = 1e300;
        for (const auto& j : L.jj) { const double d = std::hypot(j.cx() - cx, j.cy() - cy); if (d < best) { best = d; jq = &j; } }
        const int ga = M.at(jq->x0 - 0.5, jq->cy()), gb = M.at(jq->x1 + 0.5, jq->cy());
        const int qc = icomp["Q" + std::to_string(q.id)];
        const int gnd = ga == qc ? gb : ga;          // ground side of the junction
        const int fq = icomp["FQ" + std::to_string(q.id)];
        ++n; if (fq == gnd) ++same_as_flux;
        nstraps_hist[int(on[gnd].size())]++;
        std::printf("Q%-3d jj rel (%.0f,%.0f)-(%.0f,%.0f): ends comps %d/%d, qubit comp %d; ground-side comp %d (area %.0f um^2, %s), FQ%d comp %d (%s); straps on it: %zu",
                    q.id, jq->x0 - cx, jq->y0 - cy, jq->x1 - cx, jq->y1 - cy, ga, gb, qc, gnd, area[std::size_t(gnd)], gnd == plane ? "PLANE" : "not the plane",
                    q.id, fq, fq == gnd ? "same" : "different", on[gnd].size());
        for (int si : on[gnd]) {
            const auto& s = L.straps[std::size_t(si)];
            const int a = M.at(s.foot_a.cx(), s.foot_a.cy()), b = M.at(s.foot_b.cx(), s.foot_b.cy());
            std::printf(" | strap rel (%.0f,%.0f)-(%.0f,%.0f) len %.0f, other comp %d%s", s.span.x0 - cx, s.span.y0 - cy, s.span.x1 - cx, s.span.y1 - cy,
                        std::max(s.span.w(), s.span.h()), a == gnd ? b : a, (a == gnd ? b : a) == plane ? " (plane)" : "");
        }
        std::printf("\n");
        if (q.id == map_q) {
            // ASCII map 2 um/char of [-180, 20] x [-200, 20] around the qubit
            std::map<int, char> lab; lab[gnd] = 'P'; lab[plane] = '.';
            const char* others = "abcdefghijklmnoprstuvwxyz"; int oi = 0;
            std::vector<std::pair<Rect, char>> isl;
            for (const auto& [nm, rs] : L.islands) for (const auto& r : rs) isl.push_back({r, nm[0] == 'Q' ? 'Q' : (nm.rfind("FQ", 0) == 0 ? 'F' : (nm[0] == 'C' ? 'C' : (nm[0] == 'B' ? 'B' : 'L')))});
            for (double y = cy + 20; y >= cy - 200; y -= 2) {
                std::string row;
                for (double x = cx - 180; x <= cx + 20; x += 2) {
                    char ch = ' ';
                    const int k = M.at(x + 0.01, y + 0.01);
                    if (k >= 0) {
                        char il = 0; for (const auto& [r, c] : isl) if (r.contains(x + 0.01, y + 0.01)) { il = c; break; }
                        if (il && k != gnd) ch = il; else if (il == 'F' && k == gnd) ch = 'f';
                        else { if (!lab.count(k)) lab[k] = oi < 25 ? others[oi++] : '?'; ch = lab[k]; }
                    }
                    for (const auto& j : L.jj) if (j.contains(x + 0.01, y + 0.01)) ch = 'J';
                    for (int si : on[gnd]) { const auto& s = L.straps[std::size_t(si)];
                        if (s.foot_a.contains(x + 0.01, y + 0.01) || s.foot_b.contains(x + 0.01, y + 0.01)) ch = '@';
                        else if (s.span.contains(x + 0.01, y + 0.01)) ch = '#'; }
                    row += ch;
                }
                std::printf("%6.0f %s\n", y - cy, row.c_str());
            }
            std::printf("legend: P = junction's ground-side component (comp %d), . = plane, letters = other ground comps, Q qubit, F flux line, f flux line metal merged into P, J junction region, # strap span, @ strap feet; ' ' etch\n", gnd);
        }
    }
    std::printf("qubits %d: junction ground side = flux-line component in %d; straps on the junction's ground component:", n, same_as_flux);
    for (auto& [k, v] : nstraps_hist) std::printf(" %d strap(s): %d qubits;", k, v);
    std::printf("\n");
    return 0;
}
