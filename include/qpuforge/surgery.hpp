// qpuforge — lattice surgery between rotated surface-code patches
//
//   SurgeryLattice        several distance-d patches on one qubit lattice, with routing (gap) data qubits between them
//   SurgerySchedule       phases of syndrome rounds; in each phase the lattice is cut into regions (a patch on its own,
//                         or two patches merged across their gap); initial and final transversal data operations
//   build_surgery         the circuit, its detectors (found from the schedule: comparisons of each stabiliser with its
//                         previous measurement across merges and splits), reference bits and observables
//   surgery_zz / _xx      logical Z_A Z_B or X_A X_B measurement between two patches (merge, d rounds, split)
//   surgery_cnot          logical CNOT by the three-patch protocol (control, ancilla |+>, target: ZZ then XX)
//   check_surgery_noiseless  stabiliser-tableau simulation of the circuit: every detector must be 0 and every
//                         observable equal to its expected value, on several random measurement branches
//   surgery_chip          a chip for the lattice with the physics of a named design (Alpha-1 v1.1, Beta-1)
//   prepare_surgery / run_surgery_frame   detector error model (existing fault enumeration), one exact matching
//                         decoder per observable, and the Pauli-frame Monte Carlo with the memory experiments' noise
//                         model (thread-count independent like run_memory_frame)
//
// Geometry. Data qubit (i, j) sits at code coordinates (2i+1, 2j+1), the face (A, B) at (2A, 2B), as in surface.hpp.
// A region is a w x h rectangle of data qubits whose origin (i0, j0) is even, so the face type (X where A + B is even)
// is the same as in a patch of its own. Its faces follow the rotated-code rules of surface.hpp: top and bottom
// boundaries keep weight-2 X faces, left and right boundaries weight-2 Z faces, same CZ order (X: NW NE SW SE,
// Z: NW SW NE SE). Two patches one gap apart vertically merge into a d x (2d+1) region whose new faces in the two seam
// rows are all Z type: their product is Z_A Z_B (ZZ merge, gap qubits prepared in |+> and read out in the X basis at
// the split). Side by side they merge into a (2d+1) x d region with new X faces in the seam columns (XX merge, gap in
// |0>, read out in Z). The design note is docs/design/lattice-surgery.md.
#pragma once
#include "qec.hpp"
#include "presets.hpp"
#include <set>

namespace qpuforge {

// ================================================================== geometry
struct SurgeryRegion {
    int i0{0}, j0{0}, w{0}, h{0};   // data rectangle; i0, j0 even, w and h odd
    char merge{0};                  // 0: a patch on its own; 'Z': ZZ merge (gap in X basis); 'X': XX merge (gap in Z basis)
    bool contains(int i, int j) const { return i >= i0 && i < i0 + w && j >= j0 && j < j0 + h; }
};

struct SurgeryPatch { std::string name; int i0{0}, j0{0}; };   // d x d data block

struct SurgeryLattice {
    int d{3};
    std::vector<SurgeryPatch> patches;
    std::vector<SurgeryRegion> all_regions;   // every region any phase uses (defines the qubit set and the couplers)
    struct Site { int x, y; bool data; char type; int color; };
    std::vector<Site> sites;                  // qubit id = index, numbered row by row (y, then x) as in surface.hpp
    std::map<std::pair<int, int>, int> at;    // code coordinates -> qubit id
    int n_qubits() const { return int(sites.size()); }
    SurgeryRegion patch_region(int p) const { return {patches[std::size_t(p)].i0, patches[std::size_t(p)].j0, d, d, 0}; }
    int data_qubit(int i, int j) const { auto it = at.find({2 * i + 1, 2 * j + 1}); return it == at.end() ? -1 : it->second; }
};

struct SurgeryFace {
    char type{'Z'};
    int x{0}, y{0}, qubit{-1};
    std::array<int, 4> order{{-1, -1, -1, -1}};   // data qubit id touched at CZ step k (or -1)
    std::vector<int> support;                     // data qubit ids (CZ order)
};

// faces of one region, sorted by qubit id (needs the lattice's qubit numbering; pass at = {} to get them unnumbered)
inline std::vector<SurgeryFace> region_faces(const SurgeryRegion& R, const std::map<std::pair<int, int>, int>& at) {
    std::vector<SurgeryFace> out;
    for (int b = 0; b <= R.h; ++b) for (int a = 0; a <= R.w; ++a) {
        const bool tb = (b == 0 || b == R.h), lr = (a == 0 || a == R.w);
        if (tb && lr) continue;
        const int A = R.i0 + a, B = R.j0 + b;
        const char type = ((A + B) % 2 == 0) ? 'X' : 'Z';
        if (tb && type != 'X') continue;
        if (lr && type != 'Z') continue;
        SurgeryFace f; f.type = type; f.x = 2 * A; f.y = 2 * B;
        if (auto it = at.find({f.x, f.y}); it != at.end()) f.qubit = it->second;
        // NW, NE, SW, SE as data (i, j)
        const std::array<std::pair<int, int>, 4> c{{{A - 1, B - 1}, {A, B - 1}, {A - 1, B}, {A, B}}};
        const std::array<int, 4> ord = type == 'X' ? std::array<int, 4>{{0, 1, 2, 3}} : std::array<int, 4>{{0, 2, 1, 3}};
        for (int step = 0; step < 4; ++step) {
            const auto [i, j] = c[std::size_t(ord[std::size_t(step)])];
            if (!R.contains(i, j)) continue;
            int q = -1;
            if (auto it = at.find({2 * i + 1, 2 * j + 1}); it != at.end()) q = it->second;
            f.order[std::size_t(step)] = q; f.support.push_back(q);
        }
        out.push_back(f);
    }
    std::sort(out.begin(), out.end(), [](const SurgeryFace& p, const SurgeryFace& q) { return p.y < q.y || (p.y == q.y && p.x < q.x); });
    return out;
}

inline SurgeryLattice make_surgery_lattice(int d, const std::vector<SurgeryPatch>& patches, const std::vector<SurgeryRegion>& regions) {
    if (d < 3 || d % 2 == 0) throw std::runtime_error("surgery: distance must be odd and >= 3");
    SurgeryLattice L; L.d = d; L.patches = patches; L.all_regions = regions;
    for (int p = 0; p < int(patches.size()); ++p) {
        if (patches[std::size_t(p)].i0 % 2 || patches[std::size_t(p)].j0 % 2) throw std::runtime_error("surgery: patch origins must be even");
        L.all_regions.push_back(L.patch_region(p));
    }
    std::map<std::pair<int, int>, std::pair<bool, char>> s;   // (x, y) -> (data, type)
    for (const auto& R : L.all_regions) {
        if (R.i0 % 2 || R.j0 % 2 || R.w % 2 == 0 || R.h % 2 == 0) throw std::runtime_error("surgery: regions need even origins and odd sizes");
        for (int j = R.j0; j < R.j0 + R.h; ++j) for (int i = R.i0; i < R.i0 + R.w; ++i) s[{2 * i + 1, 2 * j + 1}] = {true, 'D'};
        for (const auto& f : region_faces(R, {})) s[{f.x, f.y}] = {false, f.type};
    }
    std::vector<std::pair<int, int>> keys;
    for (const auto& kv : s) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end(), [](const std::pair<int, int>& p, const std::pair<int, int>& q) { return p.second < q.second || (p.second == q.second && p.first < q.first); });
    for (const auto& k : keys) {
        const auto [data, type] = s[k];
        const int color = data ? (((k.first - 1) / 2) % 2) + 2 * (((k.second - 1) / 2) % 2) : ((k.first / 2) % 2) + 2 * ((k.second / 2) % 2);
        L.at[k] = int(L.sites.size());
        L.sites.push_back({k.first, k.second, data, type, color});
    }
    return L;
}

// ================================================================== schedule
struct SurgeryPhase {
    std::string label;
    std::vector<SurgeryRegion> regions;
    int rounds{1};
};

struct SurgeryObservable {
    std::string name;
    char type{'Z'};             // Pauli type of the logical operator (decoded with the detectors of this type)
    std::vector<int> bits;      // XOR of these classical bits ...
    int expected{0};            // ... equals this in the noiseless circuit
};

struct SurgerySchedule {
    SurgeryLattice lat;
    std::vector<char> init_basis;   // per patch: 'Z' (|0>) or 'X' (|+>)
    std::vector<char> init_flip;    // per patch: 1 applies the logical X (Z basis) or Z (X basis): |1> or |->
    std::vector<SurgeryPhase> phases;
    std::vector<char> final_basis;  // per patch: transversal readout basis (the last phase must keep the patches apart)
    std::string name{"surgery"};
};

struct SurgeryExperiment {
    SurgerySchedule sched;
    Circuit circuit;
    std::vector<std::vector<int>> detectors;      // every detector, both types
    std::vector<char> detector_type;
    std::vector<std::string> detector_label;
    std::vector<int> detector_round;              // global round (total_rounds: the final data readout)
    std::vector<int> reference;                   // noiseless value of each classical bit: 0/1 deterministic, -1 random
    std::vector<SurgeryObservable> observables;
    int total_rounds{0};
    std::vector<int> phase_first_round;
    // bookkeeping for the protocol builders
    std::map<std::pair<int, int>, int> face_bit;              // (face qubit, global round) -> classical bit
    std::map<std::pair<int, int>, int> data_bit;              // (data qubit, phase index at whose start it was measured) -> bit;
                                                              // phase index = #phases for the final readout
    std::vector<std::vector<SurgeryFace>> phase_faces;        // faces active in each phase
};

namespace surgery_detail {
inline int patch_of(const SurgeryLattice& L, int i, int j) {
    for (int p = 0; p < int(L.patches.size()); ++p) if (L.patch_region(p).contains(i, j)) return p;
    return -1;
}
} // namespace surgery_detail

inline SurgeryExperiment build_surgery(const SurgerySchedule& Sd) {
    const SurgeryLattice& L = Sd.lat;
    const int np = int(L.patches.size());
    if (int(Sd.init_basis.size()) != np || int(Sd.final_basis.size()) != np) throw std::runtime_error("surgery: one init/final basis per patch");
    if (Sd.phases.empty()) throw std::runtime_error("surgery: no phases");
    SurgeryExperiment E; E.sched = Sd;
    Circuit& C = E.circuit;
    C.n = L.n_qubits(); C.name = Sd.name;
    const int P = int(Sd.phases.size());
    // per data qubit: basis it is prepared / read in when it joins or leaves (patch basis, or the merge's gap basis)
    struct DataState { bool live{false}; char basis{'Z'}; int init_t{-1}; int meas_t{-1}; int meas_bit{-1}; };
    std::vector<DataState> ds(std::size_t(C.n));
    struct FaceRec { int t{-1}, bit{-1}; std::vector<int> support; };
    std::map<int, FaceRec> last;                    // face qubit -> its previous measurement
    auto config_data = [&](const SurgeryPhase& ph) {
        std::map<int, char> m;                     // data qubit -> basis if it is a gap qubit of a merge, else patch marker
        for (const auto& R : ph.regions)
            for (int j = R.j0; j < R.j0 + R.h; ++j) for (int i = R.i0; i < R.i0 + R.w; ++i) {
                const int q = L.data_qubit(i, j);
                if (q < 0) throw std::runtime_error("surgery: region outside the lattice");
                const int p = surgery_detail::patch_of(L, i, j);
                if (p >= 0) m[q] = 'P';
                else if (R.merge == 'Z') m[q] = 'X';        // ZZ merge: gap in |+>, read out in X
                else if (R.merge == 'X') m[q] = 'Z';        // XX merge: gap in |0>, read out in Z
                else throw std::runtime_error("surgery: a stand-alone region must be a patch");
            }
        return m;
    };
    int cb = 0;
    E.reference.clear();
    auto newbit = [&](int ref) { E.reference.push_back(ref); return cb++; };
    std::vector<std::map<int, char>> cfg(static_cast<std::size_t>(P));
    for (int p = 0; p < P; ++p) {
        cfg[std::size_t(p)] = config_data(Sd.phases[std::size_t(p)]);
        std::vector<SurgeryFace> fs;
        for (const auto& R : Sd.phases[std::size_t(p)].regions) for (const auto& f : region_faces(R, L.at)) fs.push_back(f);
        std::sort(fs.begin(), fs.end(), [](const SurgeryFace& a, const SurgeryFace& b) { return a.qubit < b.qubit; });
        for (std::size_t k = 1; k < fs.size(); ++k) if (fs[k].qubit == fs[k - 1].qubit) throw std::runtime_error("surgery: regions of a phase overlap");
        E.phase_faces.push_back(fs);
    }
    auto patch_basis_of = [&](int q) {
        const auto& s = L.sites[std::size_t(q)];
        return surgery_detail::patch_of(L, (s.x - 1) / 2, (s.y - 1) / 2);
    };
    // ---- initial preparation: reset every qubit the first phase uses (qubit order), H on X-basis data, logical flips
    {
        std::set<int> used;
        for (const auto& kv : cfg[0]) used.insert(kv.first);
        for (const auto& f : E.phase_faces[0]) used.insert(f.qubit);
        for (int q : used) C.reset(q);
        for (const auto& kv : cfg[0]) {
            const int p = patch_basis_of(kv.first);
            if (p < 0) throw std::runtime_error("surgery: the first phase cannot contain gap qubits");
            ds[std::size_t(kv.first)] = {true, Sd.init_basis[std::size_t(p)], 0, -1, -1};
        }
        for (const auto& kv : cfg[0]) if (ds[std::size_t(kv.first)].basis == 'X') C.h(kv.first);
        for (int p = 0; p < np; ++p) {
            if (!Sd.init_flip[std::size_t(p)]) continue;
            const auto& pt = L.patches[std::size_t(p)];
            for (int k = 0; k < L.d; ++k) {
                // logical X: column i = i0 (x = 1 of the patch); logical Z: row j = j0
                if (Sd.init_basis[std::size_t(p)] == 'Z') C.x(L.data_qubit(pt.i0, pt.j0 + k));
                else C.u(L.data_qubit(pt.i0 + k, pt.j0), gates::Z(), "z");
            }
        }
    }
    int t = 0;   // time stamps: transition before round g at 2g, measurement of round g at 2g + 1
    int g = 0;
    for (int p = 0; p < P; ++p) {
        const auto& ph = Sd.phases[std::size_t(p)];
        if (ph.rounds < 1) throw std::runtime_error("surgery: every phase needs >= 1 round");
        if (p > 0) {
            // transition: data qubits leaving are read out in their basis, joining ones are prepared in theirs
            t = 2 * g;
            const auto& prev = cfg[std::size_t(p - 1)];
            const auto& cur = cfg[std::size_t(p)];
            std::vector<int> leave, join;
            for (const auto& kv : prev) if (!cur.count(kv.first)) leave.push_back(kv.first);
            for (const auto& kv : cur) if (!prev.count(kv.first)) join.push_back(kv.first);
            for (int q : leave) if (ds[std::size_t(q)].basis == 'X') C.h(q);
            for (int q : leave) {
                const int b = newbit(-1);
                C.measure(q, b);
                auto& s = ds[std::size_t(q)];
                s.live = false; s.meas_t = t; s.meas_bit = b;
                E.data_bit[{q, p}] = b;
            }
            for (int q : join) {
                const char basis = cur.at(q);
                if (basis == 'P') throw std::runtime_error("surgery: a patch data qubit cannot join after the first phase");
                C.reset(q);
                if (basis == 'X') C.h(q);
                ds[std::size_t(q)] = {true, basis, t, -1, -1};
            }
            // face qubits of this phase: reset (those of the previous phase were reset after its last round)
            std::set<int> prevf;
            for (const auto& f : E.phase_faces[std::size_t(p - 1)]) prevf.insert(f.qubit);
            for (const auto& f : E.phase_faces[std::size_t(p)]) if (!prevf.count(f.qubit)) C.reset(f.qubit);
        }
        E.phase_first_round.push_back(g);
        const auto& faces = E.phase_faces[std::size_t(p)];
        for (int r = 0; r < ph.rounds; ++r, ++g) {
            t = 2 * g + 1;
            for (const auto& f : faces) C.h(f.qubit);
            for (int step = 0; step < 4; ++step) {
                for (const auto& f : faces) if (f.type == 'X' && f.order[std::size_t(step)] >= 0) C.h(f.order[std::size_t(step)]);
                for (const auto& f : faces) if (f.order[std::size_t(step)] >= 0) C.cz(f.qubit, f.order[std::size_t(step)]);
                for (const auto& f : faces) if (f.type == 'X' && f.order[std::size_t(step)] >= 0) C.h(f.order[std::size_t(step)]);
            }
            for (const auto& f : faces) C.h(f.qubit);
            for (const auto& f : faces) {
                // detector: compare with the previous measurement of this face, accounting for data qubits that left
                // (read out in the face's basis) or joined (prepared in it) since; a face never measured before is
                // deterministic only if all its data were prepared in its basis
                std::vector<int> det;
                int ref = -1;
                bool ok = true;
                auto it = last.find(f.qubit);
                if (it != last.end()) {
                    const FaceRec& pr = it->second;
                    det.push_back(pr.bit);
                    ref = E.reference[std::size_t(pr.bit)];
                    for (int q : pr.support)
                        if (std::find(f.support.begin(), f.support.end(), q) == f.support.end()) {
                            const auto& s = ds[std::size_t(q)];
                            if (s.live || s.basis != f.type || s.meas_t < pr.t) { ok = false; break; }
                            det.push_back(s.meas_bit);
                            ref = ref < 0 || E.reference[std::size_t(s.meas_bit)] < 0 ? -1 : ref ^ E.reference[std::size_t(s.meas_bit)];
                        }
                    for (int q : f.support)
                        if (ok && std::find(pr.support.begin(), pr.support.end(), q) == pr.support.end()) {
                            const auto& s = ds[std::size_t(q)];
                            if (!s.live || s.basis != f.type || s.init_t < pr.t) ok = false;
                        }
                } else {
                    ref = 0;
                    for (int q : f.support) {
                        const auto& s = ds[std::size_t(q)];
                        if (!s.live || s.basis != f.type) ok = false;
                    }
                }
                const int b = newbit(ok ? ref : -1);
                C.measure(f.qubit, b);
                E.face_bit[{f.qubit, g}] = b;
                if (ok) {
                    det.insert(det.begin(), b);
                    E.detectors.push_back(det);
                    E.detector_type.push_back(f.type);
                    E.detector_label.push_back(std::string(1, f.type) + "@(" + std::to_string(f.x) + "," + std::to_string(f.y) + ") r" + std::to_string(g));
                    E.detector_round.push_back(g);
                }
                last[f.qubit] = {t, b, f.support};
            }
            const bool last_round = (p + 1 == P) && (r + 1 == ph.rounds);
            if (!last_round) {
                if (r + 1 < ph.rounds) for (const auto& f : faces) C.reset(f.qubit);
                else {
                    std::set<int> nextf;
                    for (const auto& f : E.phase_faces[std::size_t(p + 1)]) nextf.insert(f.qubit);
                    for (const auto& f : faces) if (nextf.count(f.qubit)) C.reset(f.qubit);
                }
            }
        }
    }
    E.total_rounds = g;
    // ---- final transversal readout of every data qubit of the last phase, in its patch's basis
    t = 2 * g;
    const auto& fin = cfg[std::size_t(P - 1)];
    for (const auto& kv : fin) {
        const int p = patch_basis_of(kv.first);
        if (p < 0) throw std::runtime_error("surgery: the last phase must keep the patches apart (gap qubits still merged)");
        ds[std::size_t(kv.first)].basis = Sd.final_basis[std::size_t(p)];
    }
    for (const auto& kv : fin) if (ds[std::size_t(kv.first)].basis == 'X') C.h(kv.first);
    for (const auto& kv : fin) {
        const int b = newbit(-1);
        C.measure(kv.first, b);
        E.data_bit[{kv.first, P}] = b;
    }
    // detectors closing each face whose basis matches its patch's readout. Ordered like memory_experiment for a single
    // patch: all detectors of a face together (rounds, then final). Reorder the per-round list by face.
    {
        std::vector<std::vector<int>> dets; std::vector<char> ty; std::vector<std::string> lab; std::vector<int> rd;
        std::map<int, std::vector<std::size_t>> by_face;    // face qubit -> detector indices (round order)
        std::map<int, int> bit_face;                        // classical bit -> face qubit
        for (const auto& kv : E.face_bit) bit_face[kv.second] = kv.first.first;
        for (std::size_t k = 0; k < E.detectors.size(); ++k) by_face[bit_face[E.detectors[k][0]]].push_back(k);
        std::vector<int> face_order;
        {
            std::set<int> seen;
            for (const auto& fs : E.phase_faces) for (const auto& f : fs) seen.insert(f.qubit);
            face_order.assign(seen.begin(), seen.end());
        }
        for (int fq : face_order) {
            for (std::size_t k : by_face[fq]) { dets.push_back(E.detectors[k]); ty.push_back(E.detector_type[k]); lab.push_back(E.detector_label[k]); rd.push_back(E.detector_round[k]); }
            for (const auto& f : E.phase_faces[std::size_t(P - 1)]) {
                if (f.qubit != fq) continue;
                const int p = patch_basis_of(f.support[0]);
                if (Sd.final_basis[std::size_t(p)] != f.type) continue;
                std::vector<int> det{last[f.qubit].bit};
                for (int q : f.support) det.push_back(E.data_bit[{q, P}]);
                dets.push_back(det); ty.push_back(f.type);
                lab.push_back(std::string(1, f.type) + "@(" + std::to_string(f.x) + "," + std::to_string(f.y) + ") final");
                rd.push_back(g);
            }
        }
        E.detectors = std::move(dets); E.detector_type = std::move(ty); E.detector_label = std::move(lab); E.detector_round = std::move(rd);
    }
    return E;
}

// ---------------------------------------------------------------- observable helpers
namespace surgery_detail {
// final readout bits of a patch's data row (local j) or column (local i)
inline std::vector<int> final_row(const SurgeryExperiment& E, int patch, int j) {
    const auto& L = E.sched.lat; const auto& pt = L.patches[std::size_t(patch)];
    std::vector<int> v; const int P = int(E.sched.phases.size());
    for (int k = 0; k < L.d; ++k) v.push_back(E.data_bit.at({L.data_qubit(pt.i0 + k, pt.j0 + j), P}));
    return v;
}
inline std::vector<int> final_col(const SurgeryExperiment& E, int patch, int i) {
    const auto& L = E.sched.lat; const auto& pt = L.patches[std::size_t(patch)];
    std::vector<int> v; const int P = int(E.sched.phases.size());
    for (int k = 0; k < L.d; ++k) v.push_back(E.data_bit.at({L.data_qubit(pt.i0 + i, pt.j0 + k), P}));
    return v;
}
// product of the faces of `type` that a phase measures for the first time (the seam of a merge): first-round bits
inline std::vector<int> merge_product(const SurgeryExperiment& E, int phase, char type) {
    std::set<int> before;
    if (phase > 0) for (const auto& f : E.phase_faces[std::size_t(phase - 1)]) before.insert(f.qubit);
    std::vector<int> v;
    for (const auto& f : E.phase_faces[std::size_t(phase)])
        if (f.type == type && !before.count(f.qubit)) v.push_back(E.face_bit.at({f.qubit, E.phase_first_round[std::size_t(phase)]}));
    return v;
}
// bits of the faces of `type` inside a patch's own d x d block, measured in the first round of a phase
inline std::vector<int> patch_face_product(const SurgeryExperiment& E, int phase, int patch, char type) {
    const auto& L = E.sched.lat;
    std::set<int> own;
    for (const auto& f : region_faces(L.patch_region(patch), L.at)) if (f.type == type) own.insert(f.qubit);
    std::vector<int> v;
    for (const auto& f : E.phase_faces[std::size_t(phase)])
        if (own.count(f.qubit)) v.push_back(E.face_bit.at({f.qubit, E.phase_first_round[std::size_t(phase)]}));
    if (v.size() != own.size()) throw std::runtime_error("surgery: patch faces not all measured in that phase");
    return v;
}
inline void append(std::vector<int>& a, const std::vector<int>& b) { a.insert(a.end(), b.begin(), b.end()); }
inline char basis_of(char s) { return (s == '0' || s == '1') ? 'Z' : 'X'; }
inline char flip_of(char s) { return (s == '1' || s == '-') ? 1 : 0; }
inline void check_state(char s) { if (s != '0' && s != '1' && s != '+' && s != '-') throw std::runtime_error("surgery: logical input must be 0, 1, + or -"); }
} // namespace surgery_detail

struct SurgeryRounds { int pre{1}, merge{0}, post{1}; };   // merge = 0: d rounds

// Logical parity measurement between two patches A and B. type 'Z': Z_A Z_B (A above B, ZZ merge); 'X': X_A X_B
// (A left of B, XX merge). Inputs '0', '1', '+', '-'; each patch is read out in its input basis.
// Observables (each deterministic in the noiseless circuit):
//   inputs in the measured basis:  "m_ZZ" (or "m_XX") the merge outcome = product of the new seam faces in the first
//                                  merged round; and each patch's logical in that basis
//   inputs in the other basis:     the conserved product ("X_A X_B" under a ZZ merge, "Z_A Z_B" under XX), including the
//                                  gap qubit the logical crosses at the split
inline SurgeryExperiment surgery_parity(int d, char type, char in_a, char in_b, SurgeryRounds rr = {}) {
    using namespace surgery_detail;
    check_state(in_a); check_state(in_b);
    if (type != 'Z' && type != 'X') throw std::runtime_error("surgery: parity type must be Z or X");
    if (basis_of(in_a) != basis_of(in_b)) throw std::runtime_error("surgery: both inputs must be in the same basis");
    const int R = rr.merge > 0 ? rr.merge : d;
    const bool zz = type == 'Z';
    const SurgeryPatch A{"A", 0, 0}, B{"B", zz ? 0 : d + 1, zz ? d + 1 : 0};
    const SurgeryRegion merged{0, 0, zz ? d : 2 * d + 1, zz ? 2 * d + 1 : d, type};
    SurgerySchedule S;
    S.lat = make_surgery_lattice(d, {A, B}, {merged});
    S.name = std::string(zz ? "ZZ" : "XX") + " measurement d=" + std::to_string(d) + " |" + in_a + in_b + ">";
    S.init_basis = {basis_of(in_a), basis_of(in_b)};
    S.init_flip = {flip_of(in_a), flip_of(in_b)};
    S.final_basis = S.init_basis;
    const auto sep = std::vector<SurgeryRegion>{S.lat.patch_region(0), S.lat.patch_region(1)};
    if (rr.pre > 0) S.phases.push_back({"separate", sep, rr.pre});
    S.phases.push_back({zz ? "ZZ merged" : "XX merged", {merged}, R});
    if (rr.post < 1) throw std::runtime_error("surgery: at least one round after the split");
    S.phases.push_back({"split", sep, rr.post});
    auto E = build_surgery(S);
    const int mp = rr.pre > 0 ? 1 : 0;   // merged phase index
    const char bas = basis_of(in_a);
    if (bas == type) {
        SurgeryObservable m{zz ? "m_ZZ" : "m_XX", type, merge_product(E, mp, type), flip_of(in_a) ^ flip_of(in_b)};
        E.observables.push_back(m);
        // logical of each patch next to the seam: A bottom row / right column, B top row / left column
        SurgeryObservable a{std::string(1, type) + "_A", type, zz ? final_row(E, 0, d - 1) : final_col(E, 0, d - 1), flip_of(in_a)};
        SurgeryObservable b{std::string(1, type) + "_B", type, zz ? final_row(E, 1, 0) : final_col(E, 1, 0), flip_of(in_b)};
        E.observables.push_back(a); E.observables.push_back(b);
    } else {
        // conserved product: logical column (ZZ merge) or row (XX merge) 0 of both patches and the gap qubit between them
        SurgeryObservable c{zz ? "X_A X_B" : "Z_A Z_B", bas, {}, flip_of(in_a) ^ flip_of(in_b)};
        append(c.bits, zz ? final_col(E, 0, 0) : final_row(E, 0, 0));
        append(c.bits, zz ? final_col(E, 1, 0) : final_row(E, 1, 0));
        const int gq = zz ? S.lat.data_qubit(0, d) : S.lat.data_qubit(d, 0);
        c.bits.push_back(E.data_bit.at({gq, mp + 1}));
        E.observables.push_back(c);
    }
    return E;
}

// Logical CNOT (control C -> target T) by lattice surgery through an ancilla patch A prepared in |+>:
//   C above A (ZZ merge, outcome a), A left of T (XX merge, outcome b), A read out in Z (outcome c_A).
//   Pauli-frame correction: X_T^(a + c_A + g) Z_C^(b + h + f), where g and h are the gap qubits the logicals cross at the
//   two splits and f is the product of A's X faces in the first XX round (A's X logical moves from its left column,
//   used at the C-A seam, to its right column at the A-T seam).
// Inputs are '0', '1', '+', '-' for C and T. Readout basis of C and T: 'Z' or 'X' (or 0: the input basis of each).
// Observables, all deterministic in the noiseless circuit:
//   Z readout of both: "Z_C" (= z_C, C prepared in the Z basis), "Z_T" (= z_C + z_T, both in Z), and for a control in
//                      the X basis and a Z-basis target (a Bell pair) "Z_C Z_T" (= z_T)
//   X readout of both: "X_C" (= x_C + x_T, both in X), "X_T" (= x_T, T in X), and for a control in the X basis and a
//                      Z-basis target "X_C X_T" (= x_C)
// The labels are the logical operators after the CNOT with the frame correction applied.
inline SurgeryExperiment surgery_cnot(int d, char in_c, char in_t, char readout = 0, SurgeryRounds rr = {}) {
    using namespace surgery_detail;
    check_state(in_c); check_state(in_t);
    const int R = rr.merge > 0 ? rr.merge : d;
    const SurgeryPatch C{"C", 0, 0}, A{"A", 0, d + 1}, T{"T", d + 1, d + 1};
    const SurgeryRegion ca{0, 0, d, 2 * d + 1, 'Z'}, at{0, d + 1, 2 * d + 1, d, 'X'};
    SurgerySchedule S;
    S.lat = make_surgery_lattice(d, {C, A, T}, {ca, at});
    S.name = "CNOT d=" + std::to_string(d) + " |" + in_c + in_t + ">";
    S.init_basis = {basis_of(in_c), 'X', basis_of(in_t)};
    S.init_flip = {flip_of(in_c), 0, flip_of(in_t)};
    const char rc = readout ? readout : basis_of(in_c), rt = readout ? readout : basis_of(in_t);
    if (rc != rt) throw std::runtime_error("surgery: CNOT readout basis must be the same for C and T");
    S.final_basis = {rc, 'Z', rt};
    const auto pc = S.lat.patch_region(0), pa = S.lat.patch_region(1), pt = S.lat.patch_region(2);
    if (rr.pre > 0) S.phases.push_back({"separate", {pc, pa, pt}, rr.pre});
    S.phases.push_back({"ZZ(C,A) merged", {ca, pt}, R});
    S.phases.push_back({"XX(A,T) merged", {pc, at}, R});
    if (rr.post < 1) throw std::runtime_error("surgery: at least one round after the split");
    S.phases.push_back({"split", {pc, pa, pt}, rr.post});
    auto E = build_surgery(S);
    const int pz = rr.pre > 0 ? 1 : 0, px = pz + 1;
    const char bc = basis_of(in_c), bt = basis_of(in_t);
    const int fc = flip_of(in_c), ft = flip_of(in_t);
    if (rc == 'Z') {
        // Z_C(out): C's bottom row. Z_T(out) = Z_T(row 0) + Z_A(row 0) + g + a.
        std::vector<int> zc = final_row(E, 0, d - 1), zt = final_row(E, 2, 0);
        append(zt, final_row(E, 1, 0));
        zt.push_back(E.data_bit.at({S.lat.data_qubit(d, d + 1), px + 1}));
        append(zt, merge_product(E, pz, 'Z'));
        if (bc == 'Z') E.observables.push_back({"Z_C", 'Z', zc, fc});
        if (bc == 'Z' && bt == 'Z') E.observables.push_back({"Z_T", 'Z', zt, fc ^ ft});
        if (bt == 'Z' && bc != 'Z') { std::vector<int> b = zc; append(b, zt); E.observables.push_back({"Z_C Z_T", 'Z', b, ft}); }
    } else {
        // X_T(out): T's column 0. X_C(out) = X_C(col 0) + h + f + b.
        std::vector<int> xt = final_col(E, 2, 0), xc = final_col(E, 0, 0);
        xc.push_back(E.data_bit.at({S.lat.data_qubit(0, d), px}));
        append(xc, patch_face_product(E, px, 1, 'X'));
        append(xc, merge_product(E, px, 'X'));
        if (bc == 'X' && bt == 'X') E.observables.push_back({"X_C", 'X', xc, fc ^ ft});
        if (bt == 'X') E.observables.push_back({"X_T", 'X', xt, ft});
        if (bc == 'X' && bt != 'X') { std::vector<int> b = xc; append(b, xt); E.observables.push_back({"X_C X_T", 'X', b, fc}); }
    }
    if (E.observables.empty()) throw std::runtime_error("surgery: no deterministic observable for these inputs and readout basis");
    return E;
}

// By name: op "zz", "xx" (inputs a b) or "cnot" (inputs c t; readout 0 = input basis). The decoded basis is the
// Pauli type of the observables.
inline SurgeryExperiment surgery_experiment(const std::string& op, int d, const std::string& inputs, char readout = 0, SurgeryRounds rr = {}) {
    if (inputs.size() != 2) throw std::runtime_error("surgery: inputs are two characters, e.g. 00, 01, ++, +-");
    if (op == "zz") return surgery_parity(d, 'Z', inputs[0], inputs[1], rr);
    if (op == "xx") return surgery_parity(d, 'X', inputs[0], inputs[1], rr);
    if (op == "cnot") return surgery_cnot(d, inputs[0], inputs[1], readout, rr);
    throw std::runtime_error("surgery: op must be zz, xx or cnot");
}
inline char surgery_decode_basis(const SurgeryExperiment& E) { return E.observables.at(0).type; }

// A single patch: R rounds of memory in one basis (the surgery engine's memory experiment, for cross-checks).
inline SurgeryExperiment surgery_memory(int d, int rounds, char basis) {
    SurgerySchedule S;
    S.lat = make_surgery_lattice(d, {{"A", 0, 0}}, {});
    S.name = std::string("memory-") + basis + " d=" + std::to_string(d);
    S.init_basis = {basis}; S.init_flip = {0}; S.final_basis = {basis};
    S.phases.push_back({"memory", {S.lat.patch_region(0)}, rounds});
    auto E = build_surgery(S);
    E.observables.push_back({std::string(1, basis) + "_A", basis, basis == 'Z' ? surgery_detail::final_row(E, 0, 0) : surgery_detail::final_col(E, 0, 0), 0});
    return E;
}

// ================================================================== noiseless check: stabiliser tableau
// Aaronson-Gottesman tableau (bit-packed rows): reset, H, X, Z, CZ, Z-basis measurement. Simulates the logical circuit
// (before compilation) with random outcomes drawn from `seed`.
class StabilizerTableau {
public:
    explicit StabilizerTableau(int n, std::uint64_t seed) : n_(n), W_((std::size_t(n) + 63) / 64), rng_(seed) {
        x_.assign(std::size_t(2 * n + 1) * W_, 0); z_.assign(std::size_t(2 * n + 1) * W_, 0); r_.assign(std::size_t(2 * n + 1), 0);
        for (int i = 0; i < n; ++i) { setx(i, i, 1); setz(n + i, i, 1); }
    }
    void h(int q) { for (int i = 0; i < 2 * n_; ++i) { const int a = gx(i, q), b = gz(i, q); r_[std::size_t(i)] ^= std::uint8_t(a & b); setx(i, q, b); setz(i, q, a); } }
    void x(int q) { for (int i = 0; i < 2 * n_; ++i) r_[std::size_t(i)] ^= std::uint8_t(gz(i, q)); }
    void z(int q) { for (int i = 0; i < 2 * n_; ++i) r_[std::size_t(i)] ^= std::uint8_t(gx(i, q)); }
    void cz(int a, int b) {
        for (int i = 0; i < 2 * n_; ++i) {
            const int xa = gx(i, a), xb = gx(i, b), za = gz(i, a), zb = gz(i, b);
            r_[std::size_t(i)] ^= std::uint8_t(xa & xb & (za ^ zb));
            setz(i, a, za ^ xb); setz(i, b, zb ^ xa);
        }
    }
    int measure(int q, bool* random = nullptr) {
        int p = -1;
        for (int i = n_; i < 2 * n_; ++i) if (gx(i, q)) { p = i; break; }
        if (p >= 0) {
            if (random) *random = true;
            for (int i = 0; i < 2 * n_; ++i) if (i != p && gx(i, q)) rowsum(i, p);
            copy_row(p - n_, p);
            for (std::size_t w = 0; w < W_; ++w) { x_[std::size_t(p) * W_ + w] = 0; z_[std::size_t(p) * W_ + w] = 0; }
            setz(p, q, 1);
            const int out = int(rng_.next() & 1);
            r_[std::size_t(p)] = std::uint8_t(out);
            return out;
        }
        if (random) *random = false;
        const int s = 2 * n_;
        for (std::size_t w = 0; w < W_; ++w) { x_[std::size_t(s) * W_ + w] = 0; z_[std::size_t(s) * W_ + w] = 0; }
        r_[std::size_t(s)] = 0;
        for (int i = 0; i < n_; ++i) if (gx(i, q)) rowsum(s, i + n_);
        return r_[std::size_t(s)];
    }
    void reset(int q) { if (measure(q)) x(q); }

private:
    int n_; std::size_t W_; Xoshiro256 rng_;
    std::vector<std::uint64_t> x_, z_; std::vector<std::uint8_t> r_;
    int gx(int i, int q) const { return int((x_[std::size_t(i) * W_ + std::size_t(q) / 64] >> (q % 64)) & 1); }
    int gz(int i, int q) const { return int((z_[std::size_t(i) * W_ + std::size_t(q) / 64] >> (q % 64)) & 1); }
    void setx(int i, int q, int v) { auto& w = x_[std::size_t(i) * W_ + std::size_t(q) / 64]; w = (w & ~(1ULL << (q % 64))) | (std::uint64_t(v) << (q % 64)); }
    void setz(int i, int q, int v) { auto& w = z_[std::size_t(i) * W_ + std::size_t(q) / 64]; w = (w & ~(1ULL << (q % 64))) | (std::uint64_t(v) << (q % 64)); }
    void copy_row(int dst, int src) {
        for (std::size_t w = 0; w < W_; ++w) { x_[std::size_t(dst) * W_ + w] = x_[std::size_t(src) * W_ + w]; z_[std::size_t(dst) * W_ + w] = z_[std::size_t(src) * W_ + w]; }
        r_[std::size_t(dst)] = r_[std::size_t(src)];
    }
    // row h <- row i * row h, with the phase of the Pauli product (exponent of i, mod 4)
    void rowsum(int h, int i) {
        int e = 2 * r_[std::size_t(h)] + 2 * r_[std::size_t(i)];
        for (std::size_t w = 0; w < W_; ++w) {
            const std::uint64_t x1 = x_[std::size_t(i) * W_ + w], z1 = z_[std::size_t(i) * W_ + w];
            const std::uint64_t x2 = x_[std::size_t(h) * W_ + w], z2 = z_[std::size_t(h) * W_ + w];
            // g(x1,z1,x2,z2) = +1 / -1 contributions (Aaronson-Gottesman)
            const std::uint64_t plus = (x1 & z1 & z2 & ~x2) | (x1 & ~z1 & x2 & z2) | (~x1 & z1 & x2 & ~z2);
            const std::uint64_t minus = (x1 & z1 & x2 & ~z2) | (x1 & ~z1 & ~x2 & z2) | (~x1 & z1 & x2 & z2);
            e += __builtin_popcountll(plus) - __builtin_popcountll(minus);
            x_[std::size_t(h) * W_ + w] = x1 ^ x2; z_[std::size_t(h) * W_ + w] = z1 ^ z2;
        }
        r_[std::size_t(h)] = std::uint8_t((((e % 4) + 4) % 4) / 2);
    }
};

// One noiseless run of a circuit (ops reset, h, x, z, cz, measure): classical bits, and which were random.
inline std::vector<int> tableau_records(const Circuit& C, std::uint64_t seed, std::vector<char>* random = nullptr) {
    StabilizerTableau T(C.n, seed);
    std::vector<int> bits(std::size_t(std::max(C.nclbits, 1)), 0);
    if (random) random->assign(bits.size(), 0);
    for (const auto& g : C.gates) {
        switch (g.op) {
        case Op::RESET: T.reset(g.q[0]); break;
        case Op::CZ: T.cz(g.q[0], g.q[1]); break;
        case Op::MEASURE: { bool rnd = false; bits[std::size_t(g.cbit)] = T.measure(g.q[0], &rnd); if (random) (*random)[std::size_t(g.cbit)] = rnd; break; }
        case Op::U:
            if (g.name == "h") T.h(g.q[0]);
            else if (g.name == "x") T.x(g.q[0]);
            else if (g.name == "z") T.z(g.q[0]);
            else throw std::runtime_error("tableau: unsupported gate " + g.name);
            break;
        default: throw std::runtime_error("tableau: unsupported op");
        }
    }
    return bits;
}

struct SurgeryCheck {
    int samples{0}, fired{0}, wrong_observables{0}, wrong_reference{0}, random_reference{0};
    std::vector<int> observable_values;   // last sample
};
// Noiseless verification on `samples` random measurement branches: detectors (both types) never fire, observables
// take their expected values, and bits marked deterministic in the reference are deterministic with that value.
inline SurgeryCheck check_surgery_noiseless(const SurgeryExperiment& E, int samples = 4, std::uint64_t seed = 1) {
    SurgeryCheck c; c.samples = samples;
    for (int s = 0; s < samples; ++s) {
        std::vector<char> rnd;
        const auto bits = tableau_records(E.circuit, frame_shot_seed(seed, std::uint64_t(s)), &rnd);
        for (const auto& det : E.detectors) { int v = 0; for (int b : det) v ^= bits[std::size_t(b)]; c.fired += v; }
        c.observable_values.clear();
        for (const auto& o : E.observables) {
            int v = 0; for (int b : o.bits) v ^= bits[std::size_t(b)];
            c.observable_values.push_back(v);
            if (v != o.expected) ++c.wrong_observables;
        }
        for (std::size_t b = 0; b < E.reference.size(); ++b) {
            if (E.reference[b] < 0) continue;
            if (rnd[b]) ++c.random_reference;
            else if (bits[b] != E.reference[b]) ++c.wrong_reference;
        }
    }
    return c;
}

// ================================================================== chip
// A chip for the lattice: qubit sites and couplers as generate() lays out a surface-code chip (same frequency plan by
// site colour, same 45-degree placement, couplers in CZ order per face), inverse design, then the physics of the named
// design: "alpha1_v11" (Alpha-1 v1.1) or "beta1" (Beta-1). The readout plan is Beta-1's (at most 9 resonators per
// feedline, plan_readout) for both: Alpha-1's three feedlines do not scale to a multi-patch chip.
// For a single patch with "beta1" the result equals beta1(d) qubit for qubit (tested).
inline ChipDesign surgery_chip(const SurgeryLattice& L, const std::string& preset) {
    ChipDesign D0;
    if (preset == "beta1") D0 = beta1(L.d);
    else if (preset == "alpha1_v11") D0 = alpha1_v11(L.d);
    else throw std::runtime_error("surgery_chip: preset must be alpha1_v11 or beta1");
    GenSpec s = surface_spec(L.d);
    s.architecture = "tunable-coupler-cz"; s.g = 0.004; s.kappa_MHz = 4.0; s.chi_MHz = 4.0; s.purcell_filter = 100;   // alpha1_with_coupling
    ChipDesign G;
    G.architecture = s.architecture; G.Q_dielectric = s.Q; G.flux_noise_uPhi0 = s.flux_noise_uPhi0;
    const int n = L.n_qubits();
    G.qubits.resize(std::size_t(n));
    std::vector<double> ftarget(std::size_t(n), 0.0);
    int xmin = 1 << 30, ymax = -(1 << 30);
    for (const auto& st : L.sites) { xmin = std::min(xmin, st.x); ymax = std::max(ymax, st.y); }
    const bool single = L.patches.size() == 1 && L.all_regions.size() == 1;
    for (int k = 0; k < n; ++k) {
        const auto& st = L.sites[std::size_t(k)];
        auto& q = G.qubits[std::size_t(k)];
        q.id = k;
        q.role = st.data ? "data" : (st.type == 'X' ? "measure-X" : "measure-Z");
        const int x = st.x - (single ? 0 : xmin), y = st.y;
        q.x_um = (x + y) * s.pitch_um / 2;
        q.y_um = (y - x + 2 * (single ? L.d : std::max(L.d, (ymax + 1) / 2))) * s.pitch_um / 2;
        ftarget[std::size_t(k)] = (st.data ? s.f_lo : s.f_hi) + (s.plan == 2 ? 0.0 : (st.color - 1.5) * s.stagger);
    }
    // couplers: every face of every region to its data (first time a pair appears), faces in qubit order
    std::set<std::pair<int, int>> seen;
    std::map<int, std::vector<int>> by_face;
    for (const auto& R : L.all_regions) for (const auto& f : region_faces(R, L.at)) {
        auto& v = by_face[f.qubit];
        for (int q : f.support) if (std::find(v.begin(), v.end(), q) == v.end()) v.push_back(q);
    }
    for (auto& [fq, v] : by_face) {
        // CZ order of the full (weight-4) face: NW NE SW SE for X, NW SW NE SE for Z
        const auto& st = L.sites[std::size_t(fq)];
        auto rank = [&](int q) {
            const auto& d = L.sites[std::size_t(q)];
            const int k = (d.x < st.x ? 0 : 1) + (d.y < st.y ? 0 : 2);   // NW 0, NE 1, SW 2, SE 3
            return st.type == 'X' ? k : (k == 1 ? 2 : (k == 2 ? 1 : k));
        };
        std::sort(v.begin(), v.end(), [&](int a, int b) { return rank(a) < rank(b); });
        for (int q : v) G.couplings.push_back({fq, q, 1.0});
    }
    // readout plan of generate(): bands with s.feedlines, then tune_readout (as the preset's generate() call does)
    {
        const int F = std::max(1, s.feedlines);
        std::vector<int> order(static_cast<std::size_t>(n));
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            const auto& qa = G.qubits[std::size_t(a)]; const auto& qb = G.qubits[std::size_t(b)];
            return qa.y_um < qb.y_um || (qa.y_um == qb.y_um && qa.x_um < qb.x_um);
        });
        const int per = (n + F - 1) / F;
        for (int k = 0; k < n; ++k) {
            auto& q = G.qubits[std::size_t(order[std::size_t(k)])];
            q.feedline = k / per;
            q.ro.f_GHz = s.ro_f0 + s.ro_step * (k % per);
        }
    }
    inverse_design(G, ftarget, s.alpha, s.g);
    tune_readout(G, s.kappa_MHz, s.chi_MHz, s.purcell_filter);
    ChipDesign D = D0;
    D.qubits = G.qubits; D.couplings = G.couplings;
    plan_readout(D, 9);
    if (preset == "beta1") tune_readout(D, 8.0, 8.0, 1000);
    else tune_readout(D, 4.0, 4.0, 100);
    D.surface_d = single ? L.d : 0;
    std::string pn = preset == "beta1" ? "Beta-1" : "Alpha-1 v1.1";
    D.name = pn + " physics, " + std::to_string(L.patches.size()) + " patch" + (L.patches.size() > 1 ? "es" : "") + " d=" + std::to_string(L.d) +
             " (" + std::to_string(n) + " qubits)";
    return D;
}

// ================================================================== decoding setup and frame Monte Carlo
struct SurgerySetup {
    SurgeryExperiment exp;
    char basis{'Z'};                          // detectors and observables of this type are decoded
    MemoryExperiment view;                    // circuit + decoded detectors (for build_frame_program / enumerate_faults)
    std::vector<int> obs;                     // indices into exp.observables that are decoded
    Compiled prog;
    FrameProgram frame;
    std::vector<DetectorErrorModel> dem;      // one per decoded observable (same detectors, its own observable flags)
    std::vector<std::unique_ptr<MatchingDecoder>> decoder;
    std::vector<char> is_data;
    std::size_t n_faults{0};
};

// Single-fault signatures for each decoded observable (decoded detectors only). One enumeration serves all of them:
// each observable is appended to the detector list as an extra detector, then split off again.
inline std::vector<std::vector<FaultSignature>> surgery_fault_signatures(SurgerySetup& M) {
    MemoryExperiment ext = M.view;
    const int nd = int(M.view.detectors.size());
    for (int k : M.obs) ext.detectors.push_back(M.exp.observables[std::size_t(k)].bits);
    const auto faults = enumerate_faults(M.frame, ext);
    M.n_faults = faults.size();
    std::vector<std::vector<FaultSignature>> out(M.obs.size());
    for (std::size_t o = 0; o < M.obs.size(); ++o) {
        out[o].reserve(faults.size());
        for (const auto& f : faults) {
            FaultSignature s; s.p = f.p; s.moment = f.moment; s.obs = false; s.what = f.what;
            for (int dt : f.dets) { if (dt < nd) s.dets.push_back(dt); else if (dt == nd + int(o)) s.obs = true; }
            out[o].push_back(std::move(s));
        }
    }
    return out;
}

inline SurgerySetup prepare_surgery(const ChipDesign& D, const ChipReport& R, const NoiseModel& nm, const SurgeryExperiment& E, char basis) {
    SurgerySetup M;
    M.exp = E; M.basis = basis;
    M.view.circuit = E.circuit; M.view.reference = E.reference; M.view.rounds = E.total_rounds; M.view.d = E.sched.lat.d; M.view.basis = basis;
    for (std::size_t k = 0; k < E.detectors.size(); ++k)
        if (E.detector_type[k] == basis) {
            M.view.detectors.push_back(E.detectors[k]); M.view.detector_label.push_back(E.detector_label[k]); M.view.detector_round.push_back(E.detector_round[k]);
        }
    for (std::size_t k = 0; k < E.observables.size(); ++k) if (E.observables[k].type == basis) M.obs.push_back(int(k));
    if (M.obs.empty()) throw std::runtime_error("prepare_surgery: no observable of the decoded type");
    if (int(D.qubits.size()) != E.circuit.n) throw std::runtime_error("prepare_surgery: chip does not match the lattice");
    CompileOptions co; co.layout = "trivial";
    M.prog = compile(E.circuit, chip_graph(D, R), co);
    if (M.prog.n_swaps) throw std::runtime_error("surgery circuit needed SWAPs: couplers do not match the lattice");
    std::map<int, std::vector<int>> partners;
    if (nm.leakage && nm.dqlr) {   // data-qubit leakage removal: pair each data qubit with an adjacent face (least loaded)
        std::map<int, int> load;
        for (int q = 0; q < E.circuit.n; ++q) {
            if (!E.sched.lat.sites[std::size_t(q)].data) continue;
            int best = -1;
            for (const auto& c : D.couplings) {
                const int other = c.a == q ? c.b : (c.b == q ? c.a : -1);
                if (other >= 0 && (best < 0 || load[other] < load[best])) best = other;
            }
            if (best >= 0) { ++load[best]; partners[D.index_of(best)].push_back(D.index_of(q)); }
        }
    }
    M.frame = build_frame_program(M.prog, nm, true, E.reference, partners);
    {
        const auto sim = sim_index(M.prog);
        M.is_data.assign(std::size_t(M.frame.n), 0);
        for (int q = 0; q < E.circuit.n; ++q)
            if (E.sched.lat.sites[std::size_t(q)].data) { const int s = sim[std::size_t(D.index_of(q))]; if (s >= 0) M.is_data[std::size_t(s)] = 1; }
    }
    auto fk = surgery_fault_signatures(M);
    for (auto& f : fk) {
        M.dem.push_back(build_dem(f, int(M.view.detectors.size())));
        M.decoder.push_back(std::make_unique<MatchingDecoder>(M.dem.back()));
    }
    return M;
}

struct SurgeryResult {
    long shots{0};
    std::vector<std::string> names;        // decoded observables
    std::vector<long> errors;              // per observable
    std::vector<double> p, p_sem;
    long any_errors{0};                    // shots with at least one decoded observable wrong
    double p_any{0}, p_any_sem{0};
    double detection_rate{0}, mean_defects{0};
    int max_defects{0};
    int n_detectors{0}, rounds{0};
    double total_ns{0};
    int hyperedge_faults{0}, undetectable_logical{0};
    double hyperedge_p{0};
    double leak_events_per_shot{0};
};

// Pauli-frame Monte Carlo of a surgery experiment. The noise sampling is run_memory_frame's (same model, same order of
// random draws per operation), and so is the determinism contract: per-shot generator from frame_shot_seed(seed, shot),
// shots in fixed blocks of kFrameBlock whose tallies are added in block order. For a single-patch schedule this gives
// the same logical errors as run_memory_frame shot for shot (tested). Each observable is decoded with its own exact
// matching decoder.
inline SurgeryResult run_surgery_frame(const SurgerySetup& M, long shots, std::uint64_t seed) {
    const auto& F = M.frame;
    const std::size_t no = M.obs.size();
    struct Tally { long any{0}, fired{0}, leak_ev{0}; int max_def{0}; std::vector<long> err; };
    const long nblocks = shots > 0 ? (shots + kFrameBlock - 1) / kFrameBlock : 0;
    std::vector<Tally> part(static_cast<std::size_t>(nblocks));
    std::vector<std::vector<int>> obs_bits;
    for (int k : M.obs) obs_bits.push_back(M.exp.observables[std::size_t(k)].bits);
#pragma omp parallel
    {
        Xoshiro256 rng(0);
        std::vector<std::uint8_t> x(std::size_t(F.n)), z(std::size_t(F.n)), flips(std::size_t(std::max(F.nclbits, 1))), L(std::size_t(F.n));
        std::vector<int> fired; DecodeScratch ws;
        auto randomize = [&](int q) { x[std::size_t(q)] = rng.uniform() < 0.5; z[std::size_t(q)] = rng.uniform() < 0.5; };
#pragma omp for schedule(dynamic, 1)
        for (long b = 0; b < nblocks; ++b) {
            Tally& T_ = part[std::size_t(b)];
            T_.err.assign(no, 0);
            const long s_end = std::min(shots, (b + 1) * kFrameBlock);
            for (long s = b * kFrameBlock; s < s_end; ++s) {
                rng = Xoshiro256(frame_shot_seed(seed, std::uint64_t(s)));
                std::fill(x.begin(), x.end(), 0); std::fill(z.begin(), z.end(), 0); std::fill(flips.begin(), flips.end(), 0);
                std::fill(L.begin(), L.end(), 0);
                int nleaked = 0;
                for (const auto& fm : F.moments) {
                    for (const auto& o : fm.ops) {
                        if (nleaked == 0 || !F.leakage) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                        else {
                            const bool la = L[std::size_t(o.a)] != 0;
                            if (o.op == Op::CZ) {
                                const bool lb = L[std::size_t(o.b)] != 0;
                                if (!la && !lb) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                                else if (la != lb) {
                                    const int leaked = la ? o.a : o.b, other = la ? o.b : o.a;
                                    if (rng.uniform() < 0.5) z[std::size_t(other)] ^= 1;
                                    if (F.transport > 0 && rng.uniform() < F.transport) { L[std::size_t(leaked)] = 0; randomize(leaked); L[std::size_t(other)] = 1; }
                                }
                            } else if (o.op == Op::MEASURE) {
                                if (la) flips[std::size_t(o.cbit)] = o.ref >= 0 ? std::uint8_t(1 ^ o.ref) : std::uint8_t(rng.uniform() < 0.5);
                                else frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                            } else if (o.op == Op::RESET) {
                                if (!la) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                                else if (F.reset_removes_leakage) { L[std::size_t(o.a)] = 0; --nleaked; frame_detail::apply_op(o, x.data(), z.data(), flips.data()); }
                            } else if (!la) frame_detail::apply_op(o, x.data(), z.data(), flips.data());
                        }
                        if (o.p > 0 && rng.uniform() < o.p) {
                            if (o.op == Op::MEASURE) flips[std::size_t(o.cbit)] ^= 1;
                            else if (o.op == Op::RESET) x[std::size_t(o.a)] ^= 1;
                        }
                        if (F.leakage) {
                            if (o.op == Op::CZ && o.pleak > 0 && !L[std::size_t(o.mover)] && rng.uniform() < o.pleak) {
                                const int other = o.mover == o.a ? o.b : o.a;
                                L[std::size_t(o.mover)] = 1; ++nleaked; T_.leak_ev += 1;
                                if (!L[std::size_t(other)]) randomize(other);
                            } else if ((o.op == Op::SX || o.op == Op::X) && F.leak_1q > 0 && !L[std::size_t(o.a)] && rng.uniform() < F.leak_1q) {
                                L[std::size_t(o.a)] = 1; ++nleaked; T_.leak_ev += 1;
                            }
                        }
                    }
                    if (F.p1 > 0) for (int q : fm.ctrl1) if (rng.uniform() < F.p1) frame_detail::apply_pauli(1 + int(rng.uniform() * 3), q, x.data(), z.data());
                    for (const auto& c2 : fm.ctrl2) if (c2.p > 0 && rng.uniform() < c2.p) {
                        const int k = 1 + int(rng.uniform() * 15);
                        if (k % 4) frame_detail::apply_pauli(k % 4, c2.a, x.data(), z.data());
                        if (k / 4) frame_detail::apply_pauli(k / 4, c2.b, x.data(), z.data());
                    }
                    for (const auto& t : fm.zz) {
                        const double r = rng.uniform();
                        if (r < t.pa) z[std::size_t(t.a)] ^= 1;
                        else if (r < t.pa + t.pb) z[std::size_t(t.b)] ^= 1;
                        else if (r < t.pa + t.pb + t.pab) { z[std::size_t(t.a)] ^= 1; z[std::size_t(t.b)] ^= 1; }
                    }
                    for (const auto& dq : fm.deco) {
                        const double r = rng.uniform();
                        if (r < dq.px) x[std::size_t(dq.q)] ^= 1;
                        else if (r < dq.px + dq.py) { x[std::size_t(dq.q)] ^= 1; z[std::size_t(dq.q)] ^= 1; }
                        else if (r < dq.px + dq.py + dq.pz) z[std::size_t(dq.q)] ^= 1;
                    }
                    if (F.leakage) {
                        for (int q : fm.dqlr) {
                            if (L[std::size_t(q)] && rng.uniform() < F.dqlr_eff) { L[std::size_t(q)] = 0; --nleaked; randomize(q); }
                            if (F.dqlr_err > 0 && rng.uniform() < F.dqlr_err) frame_detail::apply_pauli(1 + int(rng.uniform() * 3), q, x.data(), z.data());
                        }
                        if (nleaked > 0)
                            for (int q = 0; q < F.n; ++q)
                                if (L[std::size_t(q)] && rng.uniform() < 1 - std::exp(-fm.duration_ns * F.seep_rate[std::size_t(q)])) {
                                    L[std::size_t(q)] = 0; --nleaked; randomize(q);
                                }
                    }
                }
                fired.clear();
                for (std::size_t d = 0; d < M.view.detectors.size(); ++d) {
                    std::uint8_t v = 0; for (int c : M.view.detectors[d]) v ^= flips[std::size_t(c)];
                    if (v) fired.push_back(int(d));
                }
                T_.fired += long(fired.size());
                T_.max_def = std::max(T_.max_def, int(fired.size()));
                bool any = false;
                for (std::size_t k = 0; k < no; ++k) {
                    std::uint8_t ob = 0; for (int c : obs_bits[k]) ob ^= flips[std::size_t(c)];
                    if (M.decoder[k]->decode(fired, ws) != ob) { ++T_.err[k]; any = true; }
                }
                if (any) ++T_.any;
            }
        }
    }
    SurgeryResult r;
    r.shots = shots; r.errors.assign(no, 0);
    long fired_tot = 0, leak_ev = 0;
    for (const Tally& t : part) {
        r.any_errors += t.any; fired_tot += t.fired; leak_ev += t.leak_ev; r.max_defects = std::max(r.max_defects, t.max_def);
        for (std::size_t k = 0; k < no; ++k) r.errors[k] += t.err[k];
    }
    auto sem = [&](double p) { return std::sqrt(std::max(p * (1 - p), 1e-12) / double(std::max(shots, 1L))); };
    for (std::size_t k = 0; k < no; ++k) {
        r.names.push_back(M.exp.observables[std::size_t(M.obs[k])].name);
        r.p.push_back(double(r.errors[k]) / double(std::max(shots, 1L)));
        r.p_sem.push_back(sem(r.p.back()));
    }
    r.p_any = double(r.any_errors) / double(std::max(shots, 1L)); r.p_any_sem = sem(r.p_any);
    r.n_detectors = int(M.view.detectors.size());
    r.mean_defects = double(fired_tot) / double(std::max(shots, 1L));
    r.detection_rate = r.n_detectors ? r.mean_defects / r.n_detectors : 0.0;
    r.rounds = M.exp.total_rounds;
    r.total_ns = M.prog.duration_ns;
    for (const auto& dm : M.dem) {
        r.hyperedge_faults = std::max(r.hyperedge_faults, dm.hyperedge_faults);
        r.hyperedge_p = std::max(r.hyperedge_p, dm.hyperedge_p);
        r.undetectable_logical = std::max(r.undetectable_logical, dm.undetectable_logical);
    }
    r.leak_events_per_shot = double(leak_ev) / double(std::max(shots, 1L));
    return r;
}

} // namespace qpuforge
