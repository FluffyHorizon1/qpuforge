// Probe 1: AF-17 (or other chip) topology of SQUID grounds, flux-line shorts and ground straps.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <queue>
using namespace qpuforge;
using namespace qpuforge::layout;

static double rdist(const Rect& a, const Rect& b) {   // gap between two rects (0 if they touch/overlap)
    const double dx = std::max({0.0, b.x0 - a.x1, a.x0 - b.x1});
    const double dy = std::max({0.0, b.y0 - a.y1, a.y0 - b.y1});
    return std::hypot(dx, dy);
}

int main(int argc, char** argv) {
    const std::string chip = argc > 1 ? argv[1] : "AF-17";
    ChipDesign D; std::string cells;
    if (chip == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else if (chip == "BF-49") { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    else { D = gamma1(7); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto t0 = std::chrono::steady_clock::now();
    const auto L = assemble_chip(D, C);
    std::printf("%s assembled in %.1f s\n%s\n%s\n", chip.c_str(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
                drc_summary(L).c_str(), connectivity_summary(L).c_str());
    const auto M = detail::metal_map(L.outline, L.etch);
    const int plane0 = M.at(L.outline.x0 + 0.5, L.outline.y0 + 0.5);
    // component of each island
    std::map<std::string, int> icomp;
    for (const auto& [name, rs] : L.islands) {
        const Rect* big = &rs.front();
        for (const auto& r : rs) if (r.area() > big->area()) big = &r;
        icomp[name] = M.at(big->cx(), big->cy());
    }
    // strap graph over components
    struct E { int to; int strap; };
    std::vector<std::vector<E>> G(std::size_t(M.n_comp));
    for (int s = 0; s < int(L.straps.size()); ++s) {
        const auto& st = L.straps[std::size_t(s)];
        const int a = M.at(st.foot_a.cx(), st.foot_a.cy()), b = M.at(st.foot_b.cx(), st.foot_b.cy());
        G[std::size_t(a)].push_back({b, s}); G[std::size_t(b)].push_back({a, s});
    }
    // BFS tree from the plane: parent strap of each component
    std::vector<int> par(std::size_t(M.n_comp), -2), pstrap(std::size_t(M.n_comp), -1), depth(std::size_t(M.n_comp), -1);
    std::queue<int> qu; qu.push(plane0); par[std::size_t(plane0)] = -1; depth[std::size_t(plane0)] = 0;
    while (!qu.empty()) {
        const int u = qu.front(); qu.pop();
        for (const auto& e : G[std::size_t(u)]) if (par[std::size_t(e.to)] == -2) { par[std::size_t(e.to)] = u; pstrap[std::size_t(e.to)] = e.strap; depth[std::size_t(e.to)] = depth[std::size_t(u)] + 1; qu.push(e.to); }
    }
    auto bridge_len = [&](const detail::Strap& st) { return std::hypot(st.foot_a.cx() - st.foot_b.cx(), st.foot_a.cy() - st.foot_b.cy()) - detail::kStrapW; };
    auto path_straps = [&](int c) { std::vector<int> v; while (c >= 0 && par[std::size_t(c)] >= 0) { v.push_back(pstrap[std::size_t(c)]); c = par[std::size_t(c)]; } return v; };
    // users of each strap: every circuit node whose path to the plane runs through it
    std::map<int, std::vector<std::string>> users;
    auto add_user = [&](int c, const std::string& who) { for (int s : path_straps(c)) users[s].push_back(who); };
    const auto R = analyze(D);
    std::printf("\nqubit SQUIDs (jj index = qubit index):\n");
    int sq_on_plane = 0, sq_with_flux = 0;
    std::vector<double> sq_len;
    for (std::size_t i = 0; i < D.qubits.size(); ++i) {
        const Rect& j = L.jj[i];
        const int g = M.at(j.x0 - 0.5, j.cy()), isl = M.at(j.x1 + 0.5, j.cy());
        const int qid = D.qubits[i].id;
        const int fq = icomp["FQ" + std::to_string(qid)], qc = icomp["Q" + std::to_string(qid)];
        const auto ps = path_straps(g);
        double len = 0; for (int s : ps) len += bridge_len(L.straps[std::size_t(s)]);
        if (g == plane0) ++sq_on_plane;
        if (fq == g) ++sq_with_flux;
        sq_len.push_back(len);
        add_user(g, "SQ" + std::to_string(qid));
        if (true)
            std::printf("  Q%-3d jj (%.0f,%.0f)-(%.0f,%.0f): ground comp %d (plane %d) island comp %d (Q comp %d); FQ%d comp %d (%s); straps to plane %zu, bridge length %.0f um\n",
                        qid, j.x0 - D.qubits[i].x_um, j.y0 - D.qubits[i].y_um, j.x1 - D.qubits[i].x_um, j.y1 - D.qubits[i].y_um, g, plane0, isl, qc, qid, fq,
                        fq == g ? "same fragment as the SQUID ground" : "other", ps.size(), len);
    }
    std::printf("  SQUID ground on the plane itself (no strap): %d / %zu; flux short in the SQUID-ground fragment: %d / %zu\n", sq_on_plane, D.qubits.size(), sq_with_flux, D.qubits.size());
    std::sort(sq_len.begin(), sq_len.end());
    std::printf("  series bridge length SQUID ground -> plane: min %.0f, median %.0f, max %.0f um\n", sq_len.front(), sq_len[sq_len.size() / 2], sq_len.back());
    // coupler SQUIDs
    int cs_on_plane = 0, cs_with_fc = 0; std::vector<double> cs_len;
    for (std::size_t k = D.qubits.size(); k < L.jj.size(); ++k) {
        const Rect& j = L.jj[k];
        const bool horiz = j.w() >= j.h();
        const int a = horiz ? M.at(j.x0 - 0.5, j.cy()) : M.at(j.cx(), j.y0 - 0.5);
        const int b = horiz ? M.at(j.x1 + 0.5, j.cy()) : M.at(j.cx(), j.y1 + 0.5);
        // which end is the coupler?
        int g = -1; std::string cname;
        for (const auto& [n, c] : icomp) if (n[0] == 'C' && (c == a || c == b)) { cname = n; g = c == a ? b : a; }
        if (g < 0) continue;
        const int fc = icomp["F" + cname];
        const auto ps = path_straps(g);
        double len = 0; for (int s : ps) len += bridge_len(L.straps[std::size_t(s)]);
        cs_len.push_back(len);
        if (g == plane0) ++cs_on_plane;
        if (fc == g) ++cs_with_fc;
        add_user(g, "S" + cname);
    }
    std::sort(cs_len.begin(), cs_len.end());
    std::printf("coupler SQUIDs: ground on the plane itself %d / %zu; coupler flux short in the SQUID-ground fragment %d; series bridge length min %.0f median %.0f max %.0f um\n",
                cs_on_plane, L.jj.size() - D.qubits.size(), cs_with_fc, cs_len.empty() ? 0 : cs_len.front(), cs_len.empty() ? 0 : cs_len[cs_len.size() / 2], cs_len.empty() ? 0 : cs_len.back());
    // shorted lines
    std::map<int, std::vector<std::string>> by_comp;
    for (const auto& [n, c] : icomp) if (n[0] == 'R' || n.rfind("FQ", 0) == 0 || n.rfind("FC", 0) == 0) { by_comp[c].push_back(n); add_user(c, n); }
    int shared = 0; for (const auto& [c, v] : by_comp) if (c != plane0 && v.size() > 1) { ++shared; if (shared <= 5) { std::printf("  fragment %d holds", c); for (auto& s : v) std::printf(" %s", s.c_str()); std::printf("\n"); } }
    int r_on_plane = 0, r_tot = 0; for (const auto& [n, c] : icomp) if (n[0] == 'R') { ++r_tot; if (c == plane0) ++r_on_plane; }
    std::printf("resonators shorted directly on the plane: %d / %d; fragments holding more than one shorted line: %d\n", r_on_plane, r_tot, shared);
    // strap sharing
    int max_users = 0, s_max = -1; for (const auto& [s, v] : users) if (int(v.size()) > max_users) { max_users = int(v.size()); s_max = s; }
    if (s_max >= 0) {
        std::printf("most-shared strap %d (bridge %.0f um) carries the return of %d circuits:", s_max, bridge_len(L.straps[std::size_t(s_max)]), max_users);
        for (int k = 0; k < std::min(max_users, 30); ++k) std::printf(" %s", users[s_max][std::size_t(k)].c_str());
        std::printf("\n");
    }
    int maxdepth = 0; for (int c = 0; c < M.n_comp; ++c) maxdepth = std::max(maxdepth, depth[std::size_t(c)]);
    int unreached = 0; for (int c = 0; c < M.n_comp; ++c) if (depth[std::size_t(c)] < 0) ++unreached;
    std::printf("strap tree: depth max %d; components unreached %d of %d\n", maxdepth, unreached, M.n_comp);
    // strap sanity: over junctions, bumps, other straps; distance to qubits; foot margins
    int over_jj = 0, over_bump = 0, overlap = 0, over_tongue = 0;
    std::vector<std::pair<double, int>> dq;
    for (int s = 0; s < int(L.straps.size()); ++s) {
        const auto& st = L.straps[std::size_t(s)];
        for (const auto& j : L.jj) if (st.span.intersects(j)) { ++over_jj; std::printf("  strap %d span (%.1f,%.1f)-(%.1f,%.1f) over junction region (%.1f,%.1f)-(%.1f,%.1f)\n", s, st.span.x0, st.span.y0, st.span.x1, st.span.y1, j.x0, j.y0, j.x1, j.y1); }
        for (const auto& b : L.bumps) if (st.span.intersects(b)) { ++over_bump; std::printf("  strap %d span over bump (%.1f,%.1f)\n", s, b.cx(), b.cy()); }
        for (int t = s + 1; t < int(L.straps.size()); ++t) if (st.span.intersects(L.straps[std::size_t(t)].span)) ++overlap;
        double dmin = 1e9;
        for (const auto& [n, rs] : L.islands) if (n[0] == 'Q') for (const auto& r : rs) dmin = std::min(dmin, rdist(st.span, r));
        dq.push_back({dmin, s});
    }
    std::sort(dq.begin(), dq.end());
    std::printf("straps over junction regions %d, over bump pads %d, overlapping pairs %d\n", over_jj, over_bump, overlap);
    std::printf("closest straps to qubit metal:");
    for (std::size_t k = 0; k < std::min<std::size_t>(8, dq.size()); ++k) std::printf(" %.1f um (strap %d)", dq[k].first, dq[k].second);
    std::printf("\n");
    int within10 = 0, within30 = 0; for (auto& [d, s] : dq) { if (d < 10) ++within10; if (d < 30) ++within30; }
    std::printf("straps within 10 um of qubit metal: %d, within 30 um: %d (of %zu)\n", within10, within30, L.straps.size());
    // foot margins: distance from foot edge to the nearest etch rect
    std::vector<double> marg;
    for (const auto& st : L.straps) for (const Rect* ft : {&st.foot_a, &st.foot_b}) {
        double m = 1e9;
        for (const auto& e : L.etch) if (rdist(*ft, e) < 50) m = std::min(m, ft->intersects(e) ? -1.0 : rdist(*ft, e));
        marg.push_back(m);
    }
    std::sort(marg.begin(), marg.end());
    int zero = 0, neg = 0; for (double m : marg) { if (m < 0) ++neg; else if (m < 1e-6) ++zero; }
    std::printf("strap feet: %zu; overlapping etch %d; flush with an etch edge (0 um margin) %d; margin quantiles 10%% %.2f, 50%% %.2f um\n", marg.size(), neg, zero,
                marg[marg.size() / 10], marg[marg.size() / 2]);
    // which qubit-pocket straps: print the straps on the SQUID-ground fragment of the first qubits
    for (std::size_t i = 0; i < std::min<std::size_t>(3, D.qubits.size()); ++i) {
        const Rect& j = L.jj[i];
        const int g = M.at(j.x0 - 0.5, j.cy());
        for (const auto& e : G[std::size_t(g)]) {
            const auto& st = L.straps[std::size_t(e.strap)];
            const double cx = D.qubits[i].x_um, cy = D.qubits[i].y_um;
            std::printf("  Q%d pocket strap %d: foot (%.1f,%.1f) -> foot (%.1f,%.1f) rel. qubit; bridge %.0f um; lines %d couplers %d; other end comp %d (plane? %s)\n",
                        D.qubits[i].id, e.strap, st.foot_a.cx() - cx, st.foot_a.cy() - cy, st.foot_b.cx() - cx, st.foot_b.cy() - cy, bridge_len(st), st.lines_crossed,
                        st.couplers_crossed, e.to, e.to == plane0 ? "yes" : "no");
        }
    }
    // f01, EJ, L_J of the qubits (for the galvanic T1 estimate)
    double EJmin = 1e9, EJmax = 0, fmin = 1e9, fmax = 0;
    for (const auto& q : R.q) { EJmin = std::min(EJmin, q.EJ); EJmax = std::max(EJmax, q.EJ); fmin = std::min(fmin, q.f01); fmax = std::max(fmax, q.f01); }
    const double Phi0 = 2.067833848e-15, hP = 6.62607015e-34;
    auto LJ = [&](double EJ_GHz) { return std::pow(Phi0 / (2 * PI), 2) / (hP * EJ_GHz * 1e9); };
    std::printf("qubits: f01 %.3f..%.3f GHz, EJ %.2f..%.2f GHz -> L_J %.2f..%.2f nH; design T1 %.0f..%.0f us\n", fmin, fmax, EJmin, EJmax, LJ(EJmax) * 1e9, LJ(EJmin) * 1e9,
                [&] { double m = 1e9; for (auto& q : R.q) m = std::min(m, q.T1_us); return m; }(), [&] { double m = 0; for (auto& q : R.q) m = std::max(m, q.T1_us); return m; }());
    return 0;
}
