// The layout_geometry strap case (tests/test_qpuforge.cpp:1861-1882) with and without Q1: where does the strap go?
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
using namespace qpuforge; using namespace qpuforge::layout;
static void run(bool with_q1, bool with_q2) {
    const Rect outline{0, 0, 300, 100};
    const std::vector<Rect> ring{{20, 20, 80, 30}, {20, 70, 80, 80}, {20, 30, 30, 70}, {70, 30, 80, 70}};
    std::vector<std::pair<std::string, std::vector<Rect>>> isl;
    if (with_q1) isl.push_back({"Q1", {{22, 32, 28, 68}}});
    std::vector<Rect> holes = ring;
    for (const auto& r : ring) holes.push_back(r.shift(100, 0));
    if (with_q2) isl.push_back({"Q2", {{121, 77, 179, 79}, {121, 21, 179, 23}, {121, 21, 123, 79}, {177, 21, 179, 79}}});
    std::vector<Rect> metal; for (const auto& [n, rs] : isl) metal.insert(metal.end(), rs.begin(), rs.end());
    const auto e = detail::etch_of(holes, metal);
    const auto M = detail::metal_map(outline, e);
    const auto I = detail::island_index(isl);
    std::vector<char> eligible(std::size_t(M.n_comp), 1);
    for (const auto& [n, rs] : isl) eligible[std::size_t(M.at(rs.front().cx(), rs.front().cy()))] = 0;
    const int plane = M.at(1, 1);
    const auto S = detail::place_ground_straps(M, I, eligible, plane);
    std::printf("Q1 %s, Q2 %s: %zu strap(s)\n", with_q1 ? "yes" : "no ", with_q2 ? "yes" : "no ", S.size());
    for (const auto& st : S) {
        const bool left = st.span.x0 < 20.0 + 1e-9 && st.span.x1 <= 40 && st.span.y0 >= 25;
        std::printf("   span (%.0f,%.0f)-(%.0f,%.0f)  feet (%.0f,%.0f) (%.0f,%.0f)  %s\n", st.span.x0, st.span.y0, st.span.x1, st.span.y1,
                    st.foot_a.cx(), st.foot_a.cy(), st.foot_b.cx(), st.foot_b.cy(), left ? "<- crosses the ring's LEFT side (where Q1 is)" : "");
    }
}
int main() { run(true, true); run(false, true); run(false, false); }
