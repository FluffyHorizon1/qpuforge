// Probe 2: straps close to qubit metal, feet over etch, overlapping straps, XY pad vs coupler pads.
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge;
using namespace qpuforge::layout;

static double rdist(const Rect& a, const Rect& b) {
    const double dx = std::max({0.0, b.x0 - a.x1, a.x0 - b.x1});
    const double dy = std::max({0.0, b.y0 - a.y1, a.y0 - b.y1});
    return std::hypot(dx, dy);
}
static double overlap_area(const Rect& a, const Rect& b) {
    const double w = std::min(a.x1, b.x1) - std::max(a.x0, b.x0), h = std::min(a.y1, b.y1) - std::max(a.y0, b.y0);
    return (w > 0 && h > 0) ? w * h : 0;
}

int main(int argc, char** argv) {
    const std::string chip = argc > 1 ? argv[1] : "AF-17";
    ChipDesign D; std::string cells;
    if (chip == "AF-17") { D = alpha1_v11(3); cells = "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"; }
    else if (chip == "BF-49") { D = beta1(5); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    else { D = gamma1(7); cells = "/home/user/wt-fab/designs/beta1/layout_cells.json"; }
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const CellSet C = cells_from_json(Json::parse(ss.str()));
    const auto L = assemble_chip(D, C);
    std::printf("%s\n", connectivity_summary(L).c_str());
    std::map<int, std::pair<double, double>> qpos; for (const auto& q : D.qubits) qpos[q.id] = {q.x_um, q.y_um};
    // straps near qubit metal
    int n_near = 0;
    for (int s = 0; s < int(L.straps.size()); ++s) {
        const auto& st = L.straps[std::size_t(s)];
        double dmin = 1e9; std::string qn; Rect qr{};
        for (const auto& [n, rs] : L.islands) if (n[0] == 'Q') for (const auto& r : rs) { const double d = rdist(st.span, r); if (d < dmin) { dmin = d; qn = n; qr = r; } }
        if (dmin >= 10) continue;
        ++n_near;
        const auto [cx, cy] = qpos[std::stoi(qn.substr(1))];
        // parallel run: projection overlap of the span and the qubit rect along the axis they share
        const double ox = std::min(st.span.x1, qr.x1) - std::max(st.span.x0, qr.x0), oy = std::min(st.span.y1, qr.y1) - std::max(st.span.y0, qr.y0);
        if (n_near <= 40)
            std::printf("  strap %3d span x[%.1f,%.1f] y[%.1f,%.1f] (rel %s: x[%.1f,%.1f] y[%.1f,%.1f]); %.2f um from %s rect x[%.1f,%.1f] y[%.1f,%.1f]; side-by-side run %.1f um; lines %d couplers %d\n",
                        s, st.span.x0, st.span.x1, st.span.y0, st.span.y1, qn.c_str(), st.span.x0 - cx, st.span.x1 - cx, st.span.y0 - cy, st.span.y1 - cy, dmin, qn.c_str(),
                        qr.x0 - cx, qr.x1 - cx, qr.y0 - cy, qr.y1 - cy, std::max(ox, oy), st.lines_crossed, st.couplers_crossed);
    }
    std::printf("straps within 10 um of qubit metal: %d\n", n_near);
    // feet overlapping etch
    int nf = 0; double worst = 0;
    for (int s = 0; s < int(L.straps.size()); ++s) {
        const auto& st = L.straps[std::size_t(s)];
        for (const Rect* ft : {&st.foot_a, &st.foot_b}) {
            double a = 0; for (const auto& e : L.etch) a += overlap_area(*ft, e);
            if (a > 1e-9) { ++nf; worst = std::max(worst, a); if (nf <= 12) std::printf("  strap %d foot (%.2f,%.2f)-(%.2f,%.2f): %.2f um^2 of its 100 um^2 over etch\n", s, ft->x0, ft->y0, ft->x1, ft->y1, a); }
        }
    }
    std::printf("feet over etch: %d (worst %.2f um^2 of 100)\n", nf, worst);
    // overlapping strap pairs
    for (int s = 0; s < int(L.straps.size()); ++s) for (int t = s + 1; t < int(L.straps.size()); ++t) {
        const auto& a = L.straps[std::size_t(s)].span; const auto& b = L.straps[std::size_t(t)].span;
        if (a.intersects(b)) std::printf("  straps %d x[%.1f,%.1f] y[%.1f,%.1f] and %d x[%.1f,%.1f] y[%.1f,%.1f] overlap\n", s, a.x0, a.x1, a.y0, a.y1, t, b.x0, b.x1, b.y0, b.y1);
    }
    // XY pad of each qubit vs the nearest coupler island and the nearest other-qubit island
    std::printf("XY line end/pad vs neighbours (first 5 qubits):\n");
    for (std::size_t i = 0; i < std::min<std::size_t>(5, D.qubits.size()); ++i) {
        const int id = D.qubits[i].id;
        const auto& xy = [&]() -> const std::vector<Rect>& { for (const auto& [n, rs] : L.islands) if (n == "XY" + std::to_string(id)) return rs; throw 0; }();
        double dq = 1e9, dc = 1e9, dqo = 1e9; std::string cn, qon;
        for (const auto& [n, rs] : L.islands) for (const auto& r : rs) for (const auto& x : xy) {
            const double d = rdist(r, x);
            if (n == "Q" + std::to_string(id)) dq = std::min(dq, d);
            else if (n[0] == 'Q' && d < dqo) { dqo = d; qon = n; }
            else if ((n[0] == 'C' || n[0] == 'B') && d < dc) { dc = d; cn = n; }
        }
        std::printf("  XY%d: %.1f um to its qubit, %.1f um to %s, %.1f um to the nearest other qubit (%s)\n", id, dq, dc, cn.c_str(), dqo, qon.c_str());
    }
    return 0;
}
