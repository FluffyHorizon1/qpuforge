// qpuforge — forge.materials.v0 sidecar reader (the atomforge -> qpuforge materials contract, docs/DATA_CONTRACT.md)
//
// atomforge writes the sidecar next to a search report (`atomforge search A B --materials` writes
// materials-contract-<A><B>.json). materials_from_json() / load_materials() read it; apply_materials() puts its values
// into a ChipDesign before analyze():
//
//   quantity           unit              ChipDesign field (chip JSON key)          accepted value
//   Q_dielectric       "dimensionless"   Q_dielectric (materials.Q_dielectric)     null, or finite and > 0
//   junction_gap_ueV   "ueV"             gap_ueV (materials.junction_gap_ueV)      null, or finite and > 0
//   loss_tangent       "dimensionless"   Q_dielectric = 1 / loss_tangent           null, or finite and > 0
//
// loss_tangent is applied only when the caller opts in (MaterialsApplyOptions::use_loss_tangent) and Q_dielectric is
// null or absent. analyze() has one chip-wide dielectric loss rate, G_diel = 2 pi f / Q_dielectric; Q = 1 / tan(delta)
// holds under that uniform-loss model only (no participation ratios, and the `conditions` a value is quoted at are not
// used), so the conversion is never silent and the report records it.
//
// A null value changes nothing and is reported as skipped. atomforge writes every quantity it cannot compute that way:
// value null, status and provenance "placeholder", and the name in the top-level `placeholders` list. A placeholder that
// carries a number is rejected, and the `qpuforge_default_for_reference` beside a placeholder is never read.
// Rn_spread_pct and the environment fields (T_eff_mK, flux_noise_uPhi0, charge_noise_e, Tphi_other_us) are not part of
// the contract: they stay in the chip JSON, and a sidecar that lists one under `quantities` is rejected.
//
// Errors throw std::runtime_error before the design is touched: a schema other than "forge.materials.v0", an unknown
// quantity, a missing or wrong unit, a value that is neither null nor a positive finite number, or a maps_to naming
// another field. Other top-level blocks (structures_ref, atomforge_computed, not_in_contract, notes, ...) and other
// members of a quantity (reason, conditions, relation, ...) are ignored: the contract keeps additive fields on v0.
#pragma once
#include "device.hpp"
#include "json.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qpuforge {

struct MaterialsQuantity {
    std::string name;                  // key under `quantities`
    std::optional<double> value;       // empty: null in the file
    std::string unit;
    std::string maps_to;               // empty: null or absent
    std::string provenance, status;    // empty: absent
    bool placeholder{false};           // status or provenance "placeholder", or listed in the top-level `placeholders`
};

struct MaterialsProducer { std::string tool, commit, potential_or_method, timestamp; };

struct MaterialsSidecar {
    std::string schema;                          // always "forge.materials.v0" once read
    MaterialsProducer produced_by;
    std::vector<std::string> system;             // element symbols, e.g. {"Cu", "Au"}
    std::vector<MaterialsQuantity> quantities;   // file order
    const MaterialsQuantity* find(const std::string& name) const {
        for (const auto& q : quantities) if (q.name == name) return &q;
        return nullptr;
    }
};

struct MaterialsApplyOptions {
    bool use_loss_tangent{false};   // set Q_dielectric = 1 / loss_tangent when Q_dielectric is null or absent
};

// What apply_materials() did with one quantity of the sidecar.
struct MaterialsReportLine {
    std::string quantity;           // key under `quantities`
    std::string field;              // chip JSON key the quantity writes
    bool applied{false};
    double before{0}, after{0};     // that field before and after apply_materials (equal when skipped)
    std::string note;               // how it was applied, or why it was skipped
};
struct MaterialsReport {
    std::vector<MaterialsReportLine> lines;   // one per quantity, file order
    int applied() const { int n = 0; for (const auto& l : lines) n += l.applied ? 1 : 0; return n; }
    int skipped() const { return int(lines.size()) - applied(); }
};

namespace materials_detail {

inline constexpr const char* kSchema = "forge.materials.v0";

[[noreturn]] inline void fail(const std::string& msg) { throw std::runtime_error("materials sidecar: " + msg); }

inline std::string fmt(double v) { char b[32]; std::snprintf(b, sizeof b, "%.6g", v); return b; }

struct Known { const char* name; const char* unit; const char* maps_to; const char* field; };
inline const Known* known(const std::string& name) {
    static const Known k[] = {
        {"Q_dielectric", "dimensionless", "ChipDesign.materials.Q_dielectric", "materials.Q_dielectric"},
        {"loss_tangent", "dimensionless", nullptr, "materials.Q_dielectric"},   // maps_to is null in the contract
        {"junction_gap_ueV", "ueV", "ChipDesign.materials.junction_gap_ueV", "materials.junction_gap_ueV"},
    };
    for (const auto& e : k) if (name == e.name) return &e;
    return nullptr;
}

inline void check_name(const std::string& name) {
    for (const char* n : {"Rn_spread_pct", "environment", "T_eff_mK", "flux_noise_uPhi0", "charge_noise_e", "Tphi_other_us"})
        if (name == n) fail("quantities." + name + " is not part of forge.materials.v0 (see not_in_contract); it stays in the chip JSON");
    if (!known(name)) fail("unknown quantity '" + name + "' (forge.materials.v0 has Q_dielectric, loss_tangent and junction_gap_ueV)");
}

// Every rule that does not need the chip. The parser applies it, and apply_materials() again for hand-built sidecars.
inline void check_quantity(const MaterialsQuantity& q) {
    check_name(q.name);
    const Known& k = *known(q.name);
    if (q.unit != k.unit) fail(q.name + ": unit must be \"" + k.unit + "\", got \"" + q.unit + "\"");
    if (!q.maps_to.empty() && (!k.maps_to || q.maps_to != k.maps_to))
        fail(q.name + ": maps_to \"" + q.maps_to + "\" does not match the contract (" +
             (k.maps_to ? "\"" + std::string(k.maps_to) + "\"" : std::string("null")) + ")");
    if (!q.value) return;
    if (q.placeholder) fail(q.name + ": a placeholder must have value null, got " + fmt(*q.value));
    if (!std::isfinite(*q.value) || !(*q.value > 0)) fail(q.name + ": value must be null or a positive finite number, got " + fmt(*q.value));
}

// Optional string member: absent or null gives "".
inline std::string opt_str(const Json& o, const std::string& key, const std::string& where) {
    const Json& v = o[key];
    if (v.is_null()) return "";
    if (v.type() != Json::Type::Str) fail(where + key + " must be a string");
    return v.str();
}

inline std::vector<std::string> str_array(const Json& j, const std::string& key) {
    std::vector<std::string> out;
    if (!j.has(key)) return out;
    const Json& a = j[key];
    if (!a.is_arr()) fail(key + " must be an array of strings");
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].type() != Json::Type::Str) fail(key + " must be an array of strings");
        out.push_back(a[i].str());
    }
    return out;
}

} // namespace materials_detail

inline MaterialsSidecar materials_from_json(const Json& j) {
    using namespace materials_detail;
    if (!j.is_obj()) fail("the top level must be a JSON object");
    const Json& schema = j["schema"];
    if (schema.type() != Json::Type::Str || schema.str() != kSchema)
        fail(std::string("schema must be \"") + kSchema + "\"" +
             (schema.type() == Json::Type::Str ? ", got \"" + schema.str() + "\""
                                               : std::string(j.has("schema") ? ", got a non-string value" : " (the file has no schema member)")));
    MaterialsSidecar S;
    S.schema = schema.str();
    if (j.has("produced_by")) {
        const Json& p = j["produced_by"];
        if (!p.is_obj()) fail("produced_by must be an object");
        S.produced_by.tool = opt_str(p, "tool", "produced_by.");
        S.produced_by.commit = opt_str(p, "commit", "produced_by.");
        S.produced_by.potential_or_method = opt_str(p, "potential_or_method", "produced_by.");
        S.produced_by.timestamp = opt_str(p, "timestamp", "produced_by.");
    }
    S.system = str_array(j, "system");
    const std::vector<std::string> listed = str_array(j, "placeholders");   // e.g. "quantities.Q_dielectric"
    const Json& qs = j["quantities"];
    if (!qs.is_obj()) fail("quantities must be an object");
    for (const auto& [name, v] : qs.items()) {
        check_name(name);
        const std::string where = "quantities." + name + ".";
        if (!v.is_obj()) fail("quantities." + name + " must be an object");
        MaterialsQuantity q;
        q.name = name;
        if (!v.has("value")) fail(where + "value is missing (a quantity without a number has value null)");
        const Json& val = v["value"];
        if (val.type() == Json::Type::Num) q.value = val.num();
        else if (!val.is_null()) fail(where + "value must be a number or null");
        if (!v.has("unit")) fail(where + "unit is missing (every quantity states its unit)");
        if (v["unit"].type() != Json::Type::Str) fail(where + "unit must be a string");
        q.unit = v["unit"].str();
        q.maps_to = opt_str(v, "maps_to", where);
        q.provenance = opt_str(v, "provenance", where);
        q.status = opt_str(v, "status", where);
        q.placeholder = q.status == "placeholder" || q.provenance == "placeholder";
        for (const auto& p : listed) q.placeholder = q.placeholder || p == "quantities." + name;
        check_quantity(q);
        S.quantities.push_back(q);
    }
    return S;
}

inline MaterialsSidecar load_materials(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    try {
        return materials_from_json(Json::parse(ss.str()));
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(path + ": " + e.what());   // name the file: the chip design is JSON too
    }
}

// Apply the sidecar to D. On error D is left unchanged.
inline MaterialsReport apply_materials(const MaterialsSidecar& S, ChipDesign& D, const MaterialsApplyOptions& opt = {}) {
    using namespace materials_detail;
    if (S.schema != kSchema) fail(std::string("schema must be \"") + kSchema + "\", got \"" + S.schema + "\"");
    for (const auto& q : S.quantities) check_quantity(q);
    const MaterialsQuantity* Qd = S.find("Q_dielectric");
    const MaterialsQuantity* lt = S.find("loss_tangent");
    const bool q_given = Qd && Qd->value;
    const bool use_lt = opt.use_loss_tangent && !q_given && lt && lt->value;
    double Q = D.Q_dielectric, gap = D.gap_ueV;    // D itself is written only after every line has been checked
    MaterialsReport R;
    for (const auto& q : S.quantities) {
        MaterialsReportLine L;
        L.quantity = q.name;
        L.field = known(q.name)->field;
        const bool is_gap = q.name == "junction_gap_ueV";
        L.before = L.after = is_gap ? D.gap_ueV : D.Q_dielectric;
        const std::string null_note = q.placeholder ? "placeholder, value null: nothing to apply" : "value null: nothing to apply";
        if (q.name == "loss_tangent") {
            if (!q.value) L.note = null_note;
            else if (!opt.use_loss_tangent)
                L.note = "loss_tangent " + fmt(*q.value) + " not used: Q_dielectric = 1/loss_tangent (uniform-loss model) is opt-in";
            else if (q_given) L.note = "loss_tangent " + fmt(*q.value) + " not used: the sidecar gives Q_dielectric itself";
            else {
                Q = 1.0 / *q.value;
                if (!std::isfinite(Q)) fail("loss_tangent " + fmt(*q.value) + " is too small: 1 / loss_tangent is not a finite Q");
                L.applied = true;
                L.after = Q;
                L.note = "Q_dielectric = 1/loss_tangent = 1/" + fmt(*q.value) + " (uniform-loss model, opted in)";
            }
        } else if (!q.value) {
            L.note = null_note + (q.name == "Q_dielectric" && use_lt ? "; Q_dielectric is set from loss_tangent" : "");
        } else {
            (is_gap ? gap : Q) = *q.value;
            L.applied = true;
            L.after = *q.value;
            L.note = q.provenance.empty() ? "no provenance given" : "provenance " + q.provenance;
        }
        R.lines.push_back(L);
    }
    D.Q_dielectric = Q;
    D.gap_ueV = gap;
    return R;
}

} // namespace qpuforge
