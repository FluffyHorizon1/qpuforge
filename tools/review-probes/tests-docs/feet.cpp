#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main() {
    const ChipDesign D = alpha1_v11(3);
    std::ifstream f("/home/user/wt-fab/designs/alpha1-v11/layout_cells.json"); std::stringstream ss; ss << f.rdbuf();
    const auto L = assemble_chip(D, cells_from_json(Json::parse(ss.str())));
    std::map<char, int> on; int any = 0;
    for (const auto& st : L.straps) {
        bool hit = false;
        for (const Rect* ft : {&st.foot_a, &st.foot_b})
            for (const auto& [n, rs] : L.islands) for (const auto& r : rs) if (ft->grow(-1e-6).intersects(r)) { if (!hit) on[n[0]]++; hit = true; }
        any += hit;
    }
    std::printf("AF-17: %zu straps, %d with a foot on island metal (R: %d, F*: %d, Q: %d); conn.clean %s, DRC %d+%d\n", L.straps.size(), any,
                on['R'], on['F'], on['Q'], L.conn.clean() ? "yes" : "no", L.drc_count("width"), L.drc_count("space"));
}
