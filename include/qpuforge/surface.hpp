// qpuforge — rotated surface-code geometry (any odd distance d)
//
// Code coordinates: data qubits at (2i+1, 2j+1), i,j in [0,d); stabiliser ("measure") qubits on faces (2a, 2b).
// Interior faces alternate X/Z; top/bottom boundaries keep weight-2 X faces, left/right keep weight-2 Z faces.
// Logical Z = Z on the data row y = 1 (horizontal), logical X = X on the data column x = 1 (vertical).
// CZ schedule (4 steps) chosen so hook errors run perpendicular to the logical of the same type:
//   X stabilisers: NW, NE, SW, SE ("Z" shape)  -> X hooks horizontal (logical X is vertical)
//   Z stabilisers: NW, SW, NE, SE ("N" shape)  -> Z hooks vertical   (logical Z is horizontal)
#pragma once
#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace qpuforge {

struct SurfaceCode {
    int d{3};
    struct Data { int x, y, qubit, color; };
    struct Stab {
        char type;                                // 'X' or 'Z'
        int x, y, qubit, color;
        std::array<int, 4> order{{-1, -1, -1, -1}};   // data index touched at CZ step k (or -1)
        std::vector<int> support;                 // data indices
    };
    std::vector<Data> data;
    std::vector<Stab> stabs;
    std::vector<int> logical_z, logical_x;        // data indices
    int n_qubits() const { return int(data.size() + stabs.size()); }
};

inline SurfaceCode rotated_surface_code(int d) {
    if (d < 3 || d % 2 == 0) throw std::runtime_error("surface code distance must be odd and >= 3");
    SurfaceCode S; S.d = d;
    struct Site { int x, y, kind; char type; };
    std::vector<Site> sites;
    for (int j = 0; j < d; ++j) for (int i = 0; i < d; ++i) sites.push_back({2 * i + 1, 2 * j + 1, 0, 'D'});
    for (int b = 0; b <= d; ++b) for (int a = 0; a <= d; ++a) {
        const bool tb = (b == 0 || b == d), lr = (a == 0 || a == d);
        if (tb && lr) continue;                                   // corners would be weight 1
        const char type = ((a + b) % 2 == 0) ? 'X' : 'Z';
        if (tb && type != 'X') continue;
        if (lr && type != 'Z') continue;
        sites.push_back({2 * a, 2 * b, 1, type});
    }
    // number qubits row by row (y, then x)
    std::sort(sites.begin(), sites.end(), [](const Site& p, const Site& q) { return p.y < q.y || (p.y == q.y && p.x < q.x); });
    std::map<std::pair<int, int>, int> data_at;
    for (int k = 0; k < int(sites.size()); ++k) {
        const auto& s = sites[k];
        if (s.kind == 0) {
            const int i = (s.x - 1) / 2, j = (s.y - 1) / 2;
            data_at[{s.x, s.y}] = int(S.data.size());
            S.data.push_back({s.x, s.y, k, (i % 2) + 2 * (j % 2)});
        }
    }
    for (int k = 0; k < int(sites.size()); ++k) {
        const auto& s = sites[k];
        if (s.kind != 1) continue;
        SurfaceCode::Stab st;
        st.type = s.type; st.x = s.x; st.y = s.y; st.qubit = k;
        st.color = ((s.x / 2) % 2) + 2 * ((s.y / 2) % 2);
        const std::pair<int, int> NW{s.x - 1, s.y - 1}, NE{s.x + 1, s.y - 1}, SW{s.x - 1, s.y + 1}, SE{s.x + 1, s.y + 1};
        const std::array<std::pair<int, int>, 4> ord = (s.type == 'X') ? std::array<std::pair<int, int>, 4>{NW, NE, SW, SE}
                                                                       : std::array<std::pair<int, int>, 4>{NW, SW, NE, SE};
        for (int step = 0; step < 4; ++step) {
            auto it = data_at.find(ord[step]);
            if (it != data_at.end()) { st.order[step] = it->second; st.support.push_back(it->second); }
        }
        S.stabs.push_back(st);
    }
    for (int k = 0; k < int(S.data.size()); ++k) {
        if (S.data[k].y == 1) S.logical_z.push_back(k);
        if (S.data[k].x == 1) S.logical_x.push_back(k);
    }
    return S;
}

inline std::string surface_role(const SurfaceCode& S, int qubit) {
    for (const auto& dq : S.data) if (dq.qubit == qubit) return "data";
    for (const auto& st : S.stabs) if (st.qubit == qubit) return st.type == 'X' ? "measure-X" : "measure-Z";
    return "qubit";
}

} // namespace qpuforge
