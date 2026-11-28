#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
using namespace qpuforge; using namespace qpuforge::layout;
int main(int argc, char** argv) {
    const int which = std::atoi(argv[1]);
    const char* cells = which == 0 ? "/home/user/wt-fab/designs/alpha1-v11/layout_cells.json" : "/home/user/wt-fab/designs/beta1/layout_cells.json";
    const ChipDesign D = which == 0 ? alpha1_v11(3) : which == 1 ? beta1(5) : gamma1(7);
    std::ifstream f(cells); std::stringstream ss; ss << f.rdbuf();
    const auto L = assemble_chip(D, cells_from_json(Json::parse(ss.str())));
    auto owner = [&](double x, double y) -> std::string { for (const auto& [n, rs] : L.islands) for (const auto& r : rs) if (r.contains(x, y)) return n; return "-"; };
    std::map<std::string, int> per;   // island -> junctions whose metal end is on it
    int q_ok = 0, c_ok = 0;
    for (std::size_t k = 0; k < L.jj.size(); ++k) {
        const Rect& j = L.jj[k]; const bool horiz = j.w() >= j.h();
        const std::string a = horiz ? owner(j.x0 - 0.5, j.cy()) : owner(j.cx(), j.y0 - 0.5), b = horiz ? owner(j.x1 + 0.5, j.cy()) : owner(j.cx(), j.y1 + 0.5);
        const std::string isl = a != "-" ? a : b;
        per[isl]++;
        if (k < D.qubits.size()) q_ok += isl == "Q" + std::to_string(D.qubits[k].id);
        else c_ok += isl[0] == 'C';
    }
    int one = 0, nqc = 0; for (const auto& [n, rs] : L.islands) if (n[0] == 'Q' || n[0] == 'C') { ++nqc; one += per[n] == 1; }
    std::printf("%zu junctions: %d qubit SQUIDs on their own qubit, %d coupler SQUIDs on a coupler; Q/C islands with exactly one junction: %d of %d\n",
                L.jj.size(), q_ok, c_ok, one, nqc);
}
