// qpuforge — FEM reference cross-check of the extraction cells (q17; docs/design/q17-q18-feasibility.md)
//
// The method of moments in layout.hpp (zero-thickness metal on an infinitely thick substrate) is compared with a
// capacitance matrix computed offline by a finite-element solver. No solver is a build dependency: the solver runs
// elsewhere, its exported result is imported here, and the deviation is REPORTED (not enforced) against a stated
// tolerance. Decided by Finnian Bond, 2026-11-28: support both Elmer and Palace; report the tolerance first, tighten later.
//
//  * exchange format  qpuforge.fem_crosscheck.v0 (JSON, docs/DATA_CONTRACT.md §6): conductor names, Maxwell matrix in
//    fF, solver name and version, free-form mesh / domain blocks, an optional cell hash;
//  * importers        Elmer StatElecSolver `CapacitanceMatrix.dat` (whitespace rows, farads; lumped convention, i.e.
//                     positive off-diagonals, flipped to Maxwell here), Palace electrostatic `terminal-C.csv` (CSV with
//                     a header row and a leading terminal-index column, farads, Maxwell convention);
//  * compare_fem      model (CapacitanceMatrix) vs reference on the conductors they share: ΔC per entry, ΔE_C per
//                     conductor and ΔG per pair from the inverse of the shared sub-matrix, with the worst relative
//                     deviations and whether they are within the stated tolerance (reported only).
//
// Ground reduction: the MoM matrix includes the ground plane as a conductor, the FEM solves with ground as a 0 V
// boundary. With ground held at 0 V the Maxwell matrix of the other conductors is the sub-matrix without the ground
// row and column, so `compare_fem` simply drops conductors the reference does not name (ground included).
#pragma once
#include "layout.hpp"
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace qpuforge {
namespace fem {

inline constexpr const char* kSchema = "qpuforge.fem_crosscheck.v0";

struct FemReference {
    std::string solver;                 // "elmer" or "palace"
    std::string solver_version;         // free text, as reported by the solver
    std::string cell;                   // e.g. "star_meas", "star_data", "pair"
    std::string cell_hash;              // optional: hash of the exported geometry
    std::vector<std::string> names;     // conductor names (row order)
    std::vector<double> C_fF;           // Maxwell matrix, row-major, fF
    Json mesh = Json::object(), domain = Json::object();   // free-form provenance blocks, carried through
    double asymmetry{0};                // max |C - C^T| / max |C| of the imported matrix (before symmetrisation)
    int n() const { return int(names.size()); }
    double at(int i, int j) const { return C_fF[std::size_t(i) * names.size() + std::size_t(j)]; }
    int index(const std::string& s) const { for (int i = 0; i < n(); ++i) if (names[std::size_t(i)] == s) return i; return -1; }
};

namespace detail {
inline std::string trim(std::string s) {
    const auto b = s.find_first_not_of(" \t\r\n"), e = s.find_last_not_of(" \t\r\n");
    return b == std::string::npos ? std::string() : s.substr(b, e - b + 1);
}
inline double number(const std::string& tok, const std::string& where) {
    const std::string t = trim(tok);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (t.empty() || end != t.c_str() + t.size() || !std::isfinite(v)) throw std::runtime_error(where + ": not a number: '" + t + "'");
    return v;
}
inline std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}
// Checks shared by every importer, then symmetrises (the asymmetry is kept as a discretisation indicator).
inline void finish(FemReference& R, const std::string& where) {
    const int n = R.n();
    if (n < 1) throw std::runtime_error(where + ": no conductors");
    if (R.C_fF.size() != std::size_t(n) * std::size_t(n))
        throw std::runtime_error(where + ": matrix is not " + std::to_string(n) + "x" + std::to_string(n) + " (one row and column per conductor name)");
    std::set<std::string> seen;
    for (const auto& s : R.names) {
        if (s.empty()) throw std::runtime_error(where + ": empty conductor name");
        if (!seen.insert(s).second) throw std::runtime_error(where + ": conductor '" + s + "' named twice");
    }
    double mx = 0, asym = 0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        mx = std::max(mx, std::abs(R.at(i, j)));
        asym = std::max(asym, std::abs(R.at(i, j) - R.at(j, i)));
    }
    for (int i = 0; i < n; ++i) {
        if (!(R.at(i, i) > 0)) throw std::runtime_error(where + ": diagonal entry of '" + R.names[std::size_t(i)] + "' must be > 0");
        for (int j = 0; j < n; ++j) if (i != j && R.at(i, j) > 1e-9 * mx)
            throw std::runtime_error(where + ": off-diagonal (" + R.names[std::size_t(i)] + ", " + R.names[std::size_t(j)] +
                                     ") is positive: not a Maxwell matrix (check the convention: Elmer lumped vs Maxwell)");
    }
    R.asymmetry = mx > 0 ? asym / mx : 0;
    if (R.asymmetry > 0.05) throw std::runtime_error(where + ": matrix asymmetry " + std::to_string(R.asymmetry) + " > 0.05");
    std::vector<double> S(R.C_fF.size());
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) S[std::size_t(i) * std::size_t(n) + std::size_t(j)] = 0.5 * (R.at(i, j) + R.at(j, i));
    R.C_fF = S;
}
} // namespace detail

// ---- qpuforge.fem_crosscheck.v0 (the committed form)
inline FemReference reference_from_json(const Json& j, const std::string& where = "fem reference") {
    if (!j.is_obj() || j.get("schema", "") != std::string(kSchema))
        throw std::runtime_error(where + ": schema must be \"" + std::string(kSchema) + "\"");
    if (j.get("units", "fF") != std::string("fF")) throw std::runtime_error(where + ": units must be \"fF\"");
    FemReference R;
    const Json& s = j["solver"];
    R.solver = s.get("name", "");
    if (R.solver != "elmer" && R.solver != "palace") throw std::runtime_error(where + ": solver.name must be \"elmer\" or \"palace\"");
    R.solver_version = s.get("version", "");
    R.cell = j.get("cell", "");
    R.cell_hash = j.get("cell_hash", "");
    if (j["mesh"].is_obj()) R.mesh = j["mesh"];
    if (j["domain"].is_obj()) R.domain = j["domain"];
    const Json& names = j["conductors"];
    if (!names.is_arr()) throw std::runtime_error(where + ": conductors must be an array of names");
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (names[i].type() != Json::Type::Str) throw std::runtime_error(where + ": conductors[" + std::to_string(i) + "] is not a string");
        R.names.push_back(names[i].str());
    }
    const Json& M = j["maxwell_fF"];
    if (!M.is_arr() || M.size() != R.names.size()) throw std::runtime_error(where + ": maxwell_fF must have one row per conductor");
    for (std::size_t i = 0; i < M.size(); ++i) {
        if (!M[i].is_arr() || M[i].size() != R.names.size()) throw std::runtime_error(where + ": maxwell_fF row " + std::to_string(i) + " has the wrong length");
        for (std::size_t k = 0; k < M[i].size(); ++k) {
            if (M[i][k].type() != Json::Type::Num) throw std::runtime_error(where + ": maxwell_fF[" + std::to_string(i) + "][" + std::to_string(k) + "] is not a number");
            R.C_fF.push_back(M[i][k].num());
        }
    }
    detail::finish(R, where);
    return R;
}
inline Json reference_to_json(const FemReference& R) {
    Json j = Json::object();
    j.set("schema", kSchema);
    Json s = Json::object(); s.set("name", R.solver); s.set("version", R.solver_version);
    j.set("solver", s);
    j.set("cell", R.cell);
    if (!R.cell_hash.empty()) j.set("cell_hash", R.cell_hash);
    j.set("units", "fF");
    Json names = Json::array(); for (const auto& n : R.names) names.push(n);
    j.set("conductors", names);
    Json M = Json::array();
    for (int i = 0; i < R.n(); ++i) { Json row = Json::array(); for (int k = 0; k < R.n(); ++k) row.push(R.at(i, k)); M.push(row); }
    j.set("maxwell_fF", M);
    j.set("mesh", R.mesh); j.set("domain", R.domain);
    return j;
}
inline FemReference load_reference(const std::string& path) { return reference_from_json(Json::parse(detail::read_file(path)), path); }

// ---- Elmer: StatElecSolver with `Calculate Capacitance Matrix = True` writes `Capacitance Matrix Filename`
// (CapacitanceMatrix.dat by default): one row per `Capacitance Body`, whitespace-separated, in farads. The file
// has no names, so the caller gives them in Capacitance Body order. Lines starting with '!' or '#' are comments.
// Convention: lumped (mutual capacitances positive off the diagonal, the diagonal the total capacitance of the body);
// the Maxwell matrix used here has the same diagonal and negated off-diagonals.
// `lumped = false` reads a file already in the Maxwell convention (negative off-diagonals) as is.
inline FemReference import_elmer(const std::string& text, const std::vector<std::string>& names, const std::string& version = "",
                                 bool lumped = true, const std::string& where = "Elmer CapacitanceMatrix.dat") {
    FemReference R;
    R.solver = "elmer"; R.solver_version = version; R.names = names;
    std::istringstream in(text);
    std::string line;
    int row = 0;
    while (std::getline(in, line)) {
        const std::string t = detail::trim(line);
        if (t.empty() || t[0] == '!' || t[0] == '#') continue;
        std::istringstream ls(t);
        std::string tok;
        int col = 0;
        while (ls >> tok) {
            const double v = detail::number(tok, where + " row " + std::to_string(row + 1)) * 1e15;   // F -> fF
            R.C_fF.push_back(col == row || !lumped ? v : -v);
            ++col;
        }
        if (col != int(names.size()))
            throw std::runtime_error(where + ": row " + std::to_string(row + 1) + " has " + std::to_string(col) + " values, expected " + std::to_string(names.size()));
        ++row;
    }
    detail::finish(R, where);
    return R;
}

// ---- Palace: an Electrostatic simulation writes `terminal-C.csv`: a header row, then one row per terminal with the
// terminal index first and the Maxwell matrix row in farads (comma-separated). The caller gives conductor names in
// terminal order (the index column must run 1..n).
inline FemReference import_palace(const std::string& text, const std::vector<std::string>& names, const std::string& version = "",
                                  const std::string& where = "Palace terminal-C.csv") {
    FemReference R;
    R.solver = "palace"; R.solver_version = version; R.names = names;
    std::istringstream in(text);
    std::string line;
    bool header = true;
    int row = 0;
    while (std::getline(in, line)) {
        if (detail::trim(line).empty()) continue;
        if (header) { header = false; continue; }          // "i, C[i][1] (F), C[i][2] (F), ..."
        std::vector<std::string> f;
        std::stringstream ls(line); std::string tok;
        while (std::getline(ls, tok, ',')) f.push_back(tok);
        const std::string w = where + " row " + std::to_string(row + 1);
        if (f.size() != names.size() + 1)
            throw std::runtime_error(w + ": " + std::to_string(f.size()) + " columns, expected " + std::to_string(names.size() + 1) + " (index + one per terminal)");
        if (detail::number(f[0], w) != double(row + 1)) throw std::runtime_error(w + ": terminal index must be " + std::to_string(row + 1));
        for (std::size_t k = 1; k < f.size(); ++k) R.C_fF.push_back(detail::number(f[k], w) * 1e15);
        ++row;
    }
    detail::finish(R, where);
    return R;
}

// ---- comparison (reported, not enforced)
struct FemTolerance {      // stated tolerances; compare_fem reports against them and never fails on them
    double rel_C{0.02};    // per Maxwell entry (diagonal and mutual), relative to the reference
    double rel_EC{0.02};   // charging energy (1 % of C_Σ ≈ 1 % of E_C ≈ 2 MHz of α at E_C ≈ 0.22 GHz)
    double rel_G{0.05};    // charge-charge coupling G_ij (sets g)
    bool enforced{false};  // decided 2026-11-28: report first; becomes true once the tolerance is fixed
};
struct FemEntry { std::string a, b; double model{0}, reference{0}; double rel() const { return reference != 0 ? (model - reference) / reference : 0.0; } };   // rel: signed, relative to the reference value
struct FemComparison {
    std::string solver, cell;
    std::vector<std::string> shared;          // conductors compared (reference order)
    std::vector<std::string> model_only;      // in the model, not in the reference (ground, by construction)
    std::vector<FemEntry> C;                  // Maxwell entries, upper triangle incl. diagonal, fF
    std::vector<FemEntry> EC;                 // per conductor, GHz (a == b)
    std::vector<FemEntry> G;                  // per pair, GHz per n n
    double max_rel_C{0}, max_rel_EC{0}, max_rel_G{0};
    FemTolerance tol;
    bool within() const { return std::abs(max_rel_C) <= tol.rel_C && std::abs(max_rel_EC) <= tol.rel_EC && std::abs(max_rel_G) <= tol.rel_G; }
};

// E_C,i = e^2 (C^-1)_ii / 2h and G_ij = 4 e^2 (C^-1)_ij / h (GHz), the convention of layout::triple_circuit.
inline void circuit_of(const std::vector<double>& C_fF, int n, std::vector<double>& EC, std::vector<double>& G) {
    std::vector<double> C(C_fF.size());
    for (std::size_t k = 0; k < C.size(); ++k) C[k] = C_fF[k] * 1e-15;
    const auto Ci = inverse(C, n);
    const double k1 = phys::e * phys::e / (2 * phys::h) / 1e9, k2 = 4 * phys::e * phys::e / phys::h / 1e9;
    EC.assign(std::size_t(n), 0.0); G.assign(std::size_t(n) * std::size_t(n), 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        (i == j ? EC[std::size_t(i)] : G[std::size_t(i) * std::size_t(n) + std::size_t(j)]) = (i == j ? k1 : k2) * Ci[std::size_t(i) * std::size_t(n) + std::size_t(j)];
}

inline FemComparison compare_fem(const layout::CapacitanceMatrix& model, const FemReference& ref, const FemTolerance& tol = {}) {
    FemComparison out;
    out.solver = ref.solver; out.cell = ref.cell; out.tol = tol;
    std::vector<int> mi;
    for (const auto& s : ref.names) {
        int k = -1;
        for (int i = 0; i < model.n; ++i) if (model.names[std::size_t(i)] == s) k = i;
        if (k < 0) throw std::runtime_error("compare_fem: reference conductor '" + s + "' is not in the model cell");
        out.shared.push_back(s); mi.push_back(k);
    }
    for (const auto& s : model.names) if (ref.index(s) < 0) out.model_only.push_back(s);
    const int n = int(mi.size());
    std::vector<double> Cm(std::size_t(n) * std::size_t(n));
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) Cm[std::size_t(i) * std::size_t(n) + std::size_t(j)] = model.at(mi[std::size_t(i)], mi[std::size_t(j)]);
    auto upd = [](double& m, double r) { if (std::abs(r) > std::abs(m)) m = r; };
    for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) {
        FemEntry e{out.shared[std::size_t(i)], out.shared[std::size_t(j)], Cm[std::size_t(i) * std::size_t(n) + std::size_t(j)], ref.at(i, j)};
        out.C.push_back(e); upd(out.max_rel_C, e.rel());
    }
    std::vector<double> ECm, Gm, ECr, Gr;
    circuit_of(Cm, n, ECm, Gm);
    circuit_of(ref.C_fF, n, ECr, Gr);
    for (int i = 0; i < n; ++i) {
        FemEntry e{out.shared[std::size_t(i)], out.shared[std::size_t(i)], ECm[std::size_t(i)], ECr[std::size_t(i)]};
        out.EC.push_back(e); upd(out.max_rel_EC, e.rel());
        for (int j = i + 1; j < n; ++j) {
            FemEntry g{out.shared[std::size_t(i)], out.shared[std::size_t(j)], Gm[std::size_t(i) * std::size_t(n) + std::size_t(j)], Gr[std::size_t(i) * std::size_t(n) + std::size_t(j)]};
            out.G.push_back(g); upd(out.max_rel_G, g.rel());
        }
    }
    return out;
}

inline std::string comparison_summary(const FemComparison& c) {
    char b[256];
    std::snprintf(b, sizeof b, "FEM cross-check %s (%s, %zu conductors): max dC %+.2f%%, dE_C %+.2f%%, dG %+.2f%% (tolerance %.0f/%.0f/%.0f%%, %s, %s)",
                  c.cell.empty() ? "(unnamed cell)" : c.cell.c_str(), c.solver.c_str(), c.shared.size(), 100 * c.max_rel_C, 100 * c.max_rel_EC, 100 * c.max_rel_G,
                  100 * c.tol.rel_C, 100 * c.tol.rel_EC, 100 * c.tol.rel_G, c.within() ? "within" : "outside",
                  c.tol.enforced ? "enforced" : "reported only");
    return b;
}
inline Json comparison_json(const FemComparison& c) {
    Json j = Json::object();
    j.set("solver", c.solver); j.set("cell", c.cell);
    Json sh = Json::array(); for (const auto& s : c.shared) sh.push(s);
    Json mo = Json::array(); for (const auto& s : c.model_only) mo.push(s);
    j.set("conductors", sh); j.set("model_only", mo);
    auto rows = [](const std::vector<FemEntry>& v, bool pair) {
        Json a = Json::array();
        for (const auto& e : v) {
            Json o = Json::object();
            if (pair) { o.set("a", e.a); o.set("b", e.b); } else o.set("conductor", e.a);
            o.set("model", e.model); o.set("reference", e.reference); o.set("rel", e.rel());
            a.push(o);
        }
        return a;
    };
    j.set("maxwell_fF", rows(c.C, true)); j.set("E_C_GHz", rows(c.EC, false)); j.set("G_GHz", rows(c.G, true));
    Json t = Json::object();
    t.set("rel_C", c.tol.rel_C); t.set("rel_EC", c.tol.rel_EC); t.set("rel_G", c.tol.rel_G); t.set("enforced", c.tol.enforced);
    j.set("tolerance", t);
    j.set("max_rel_C", c.max_rel_C); j.set("max_rel_EC", c.max_rel_EC); j.set("max_rel_G", c.max_rel_G);
    j.set("within_tolerance", c.within());
    return j;
}

} // namespace fem
} // namespace qpuforge
