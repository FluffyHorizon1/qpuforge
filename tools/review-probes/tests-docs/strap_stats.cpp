#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <queue>
#include <map>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main(int argc, char** argv) {
    const int which = std::atoi(argv[1]);
    const char* cells = which == 0 ? "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json" : "/home/user/wt-fab/designs/beta1/layout_cells.json";
    const ChipDesign D = which == 0 ? alpha1_v11(3) : which == 1 ? beta1(5) : gamma1(7);
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    int with_line = 0, with_coup = 0, sum_l = 0, sum_c = 0, max_l = 0, max_c = 0, over_q = 0, over_jj = 0;
    for (const auto& st : L.straps) {
        with_line += st.lines_crossed > 0; with_coup += st.couplers_crossed > 0; sum_l += st.lines_crossed; sum_c += st.couplers_crossed;
        max_l = std::max(max_l, st.lines_crossed); max_c = std::max(max_c, st.couplers_crossed);
        for (const auto& [n, rs] : L.islands) if (n[0] == 'Q') for (const auto& r : rs) if (st.span.grow(-1e-6).intersects(r)) { ++over_q; goto next; }
        next:
        for (const auto& j : L.jj) if (st.span.intersects(j)) { ++over_jj; break; }
    }
    std::printf("%zu straps: %d cross >= 1 line (sum of line crossings %d, max %d per strap); %d cross >= 1 coupler/bus (sum %d, max %d); "
                "%d overlap a qubit island; %d overlap a junction region\n", L.straps.size(), with_line, sum_l, max_l, with_coup, sum_c, max_c, over_q, over_jj);
    // strap graph on metal components: airbridges in series from each qubit SQUID's ground end to the plane
    const auto M = detail::metal_map(L.outline, L.etch);
    const int plane = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    std::vector<std::vector<int>> adj(std::size_t(M.n_comp));
    for (const auto& st : L.straps) { int a = M.at(st.foot_a.cx(), st.foot_a.cy()), b = M.at(st.foot_b.cx(), st.foot_b.cy()); adj[a].push_back(b); adj[b].push_back(a); }
    std::vector<int> dist(std::size_t(M.n_comp), -1), par(std::size_t(M.n_comp), -1); std::queue<int> q; dist[plane] = 0; q.push(plane);
    while (!q.empty()) { int u = q.front(); q.pop(); for (int v : adj[u]) if (dist[v] < 0) { dist[v] = dist[u] + 1; par[v] = u; q.push(v); } }
    std::map<std::pair<int,int>, int> use;   // tree edge -> number of qubit SQUID/flux returns through it
    int maxd = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        int g = M.at(L.jj[i].x0 - 0.5, L.jj[i].cy()); maxd = std::max(maxd, dist[g]);
        while (g >= 0 && g != plane) { use[{std::min(g, par[g]), std::max(g, par[g])}]++; g = par[g]; }
    }
    int same = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const int g = M.at(L.jj[i].x0 - 0.5, L.jj[i].cy());
        for (const auto& [n, rs] : L.islands) if (n == "FQ" + std::to_string(D.qubits[i].id)) {
            const Rect* big = &rs.front(); for (const auto& r : rs) if (r.area() > big->area()) big = &r;
            same += M.at(big->cx(), big->cy()) == g;
        }
    }
    std::printf("qubits whose flux line is shorted into the same metal component as their SQUID ground lead: %d of %zu\n", same, D.qubits.size());
    {
        std::map<int, int> squid_comp;   // component -> qubit id whose SQUID ground lead it holds
        for (std::size_t i = 0; i < D.qubits.size(); ++i) squid_comp[M.at(L.jj[i].x0 - 0.5, L.jj[i].cy())] = D.qubits[i].id;
        int through = 0, worst = 0, worst_q = -1;
        for (std::size_t i = 0; i < D.qubits.size(); ++i) {
            int g = M.at(L.jj[i].x0 - 0.5, L.jj[i].cy()), k = 0;
            for (int c = par[g]; c >= 0 && c != plane; c = par[c]) if (squid_comp.count(c)) ++k;
            through += k > 0; if (k > worst) { worst = k; worst_q = D.qubits[i].id; }
        }
        std::printf("qubits whose flux return path runs through another qubit's SQUID pocket: %d of %zu (worst Q%d: through %d other SQUID pockets)\n",
                    through, D.qubits.size(), worst_q, worst);
    }
    int maxuse = 0; for (auto& [e, n] : use) maxuse = std::max(maxuse, n);
    std::printf("longest chain %d bridges; most qubit flux returns through one bridge: %d\n", maxd, maxuse);
    int hist[16] = {0}; int nq = 0, single = 0;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const Rect& j = L.jj[i];                        // qubit SQUIDs are pushed first, one per qubit (assemble_chip)
        const int g = M.at(j.x0 - 0.5, j.cy());         // ground end (tongue)
        const int d = g >= 0 ? dist[g] : -1;
        hist[std::min(15, std::max(0, d))]++; ++nq; single += g >= 0 && adj[g].size() == 1;
    }
    std::printf("qubit SQUID ground leads: airbridges in series to the plane:");
    for (int k = 0; k < 16; ++k) if (hist[k]) std::printf(" %d bridge(s): %d qubits;", k, hist[k]);
    std::printf(" ground lead's pocket touched by exactly one strap: %d of %d\n", single, nq);
}
