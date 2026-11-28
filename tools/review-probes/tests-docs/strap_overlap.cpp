#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main(int argc, char** argv) {
    const int which = std::atoi(argv[1]);
    const char* cells = which == 0 ? "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json" : "/home/user/wt-fab/designs/beta1/layout_cells.json";
    const ChipDesign D = which == 0 ? alpha1_v11(3) : which == 1 ? beta1(5) : gamma1(7);
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const auto L = assemble_chip(D, cells_from_json(Json::parse(ss.str())));
    int overlap = 0, touch = 0, foot_under = 0; double maxlen = 0, minsep = 1e9;
    for (std::size_t i = 0; i < L.straps.size(); ++i) {
        const auto& a = L.straps[i];
        maxlen = std::max(maxlen, std::max(a.span.w(), a.span.h()));
        for (std::size_t j = i + 1; j < L.straps.size(); ++j) {
            const auto& b = L.straps[j];
            if (a.span.grow(-1e-6).intersects(b.span)) { ++overlap; if (overlap <= 3) std::printf("  overlap: (%.0f,%.0f)-(%.0f,%.0f) and (%.0f,%.0f)-(%.0f,%.0f)\n", a.span.x0, a.span.y0, a.span.x1, a.span.y1, b.span.x0, b.span.y0, b.span.x1, b.span.y1); }
            else if (a.span.grow(5).intersects(b.span)) ++touch;
        }
    }
    std::printf("%zu straps: %d pairs of overlapping spans, %d further pairs closer than 5 um; longest strap %.0f um\n", L.straps.size(), overlap, touch, maxlen);
}
