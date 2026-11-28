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
// Provenance is an enumerated field on v0 (decision D4, docs/design/q16-materials-ingestion.md §6):
//   placeholder             not computed; value must be null (what atomforge writes today)
//   assumed                 a hand estimate by a person; the quantity's `reason` should say why
//   literature              a published value; `reference` names the source
//   measured                measured on hardware; `reference` names the measurement
//   dft-loss-tangent        computed by a DFT loss workflow
//   participation-weighted  weighted by FEM participation ratios
// A number must carry a provenance from this list; a null value may omit it. Numbers must lie in a physical range:
//   Q_dielectric [1e3, 1e10], loss_tangent [1e-10, 1e-3], junction_gap_ueV [10, 3000] ueV.
// The top level must name its producer (produced_by.tool), every `placeholders` entry must name a quantity in the file,
// and `contains_placeholders`, when present, must be a boolean that agrees with the quantities.
//
// Errors throw std::runtime_error before the design is touched: a schema other than "forge.materials.v0", an unknown
// quantity, a missing or wrong unit, a value that is neither null nor a positive finite number or that is out of range,
// a missing or unknown provenance, a literature/measured value without a reference, a maps_to naming another field, a
// missing producer, or inconsistent placeholder bookkeeping. Other top-level blocks (structures_ref, atomforge_computed, not_in_contract, notes, ...) and other
// members of a quantity (reason, conditions, relation, ...) are ignored: the contract keeps additive fields on v0.
#pragma once
#include "device.hpp"
#include "json.hpp"
#include <cmath>
#include <cstdint>
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
    std::string reference;             // source of a literature or measured value; empty: absent
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

struct Known { const char* name; const char* unit; const char* maps_to; const char* field; double lo, hi; };
inline const Known* known(const std::string& name) {
    static const Known k[] = {
        {"Q_dielectric", "dimensionless", "ChipDesign.materials.Q_dielectric", "materials.Q_dielectric", 1e3, 1e10},
        {"loss_tangent", "dimensionless", nullptr, "materials.Q_dielectric", 1e-10, 1e-3},   // maps_to is null in the contract
        {"junction_gap_ueV", "ueV", "ChipDesign.materials.junction_gap_ueV", "materials.junction_gap_ueV", 10, 3000},
    };
    for (const auto& e : k) if (name == e.name) return &e;
    return nullptr;
}

// The enumerated provenance values of forge.materials.v0 (D4).
inline const std::vector<std::string>& provenances() {
    static const std::vector<std::string> p{"placeholder", "assumed", "literature", "measured", "dft-loss-tangent", "participation-weighted"};
    return p;
}
inline bool known_provenance(const std::string& p) {
    for (const auto& e : provenances()) if (p == e) return true;
    return false;
}
inline std::string provenance_list() {
    std::string out;
    for (const auto& e : provenances()) out += (out.empty() ? "" : ", ") + e;
    return out;
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
    if (!q.provenance.empty() && !known_provenance(q.provenance))
        fail(q.name + ": unknown provenance \"" + q.provenance + "\" (forge.materials.v0 allows " + provenance_list() + ")");
    if (!q.value) return;
    if (q.placeholder) fail(q.name + ": a placeholder must have value null, got " + fmt(*q.value));
    if (!std::isfinite(*q.value) || !(*q.value > 0)) fail(q.name + ": value must be null or a positive finite number, got " + fmt(*q.value));
    if (*q.value < k.lo || *q.value > k.hi)
        fail(q.name + ": value " + fmt(*q.value) + " is out of range [" + fmt(k.lo) + ", " + fmt(k.hi) + "] " + k.unit);
    if (q.provenance.empty())
        fail(q.name + ": a value needs a provenance (" + provenance_list() + ")");
    if ((q.provenance == "literature" || q.provenance == "measured") && q.reference.empty())
        fail(q.name + ": provenance \"" + q.provenance + "\" needs a reference naming its source");
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
        q.reference = opt_str(v, "reference", where);
        q.placeholder = q.status == "placeholder" || q.provenance == "placeholder";
        for (const auto& p : listed) q.placeholder = q.placeholder || p == "quantities." + name;
        check_quantity(q);
        S.quantities.push_back(q);
    }
    for (const auto& p : listed) {   // placeholder bookkeeping must agree with the quantities
        const std::string prefix = "quantities.";
        if (p.compare(0, prefix.size(), prefix) != 0 || !S.find(p.substr(prefix.size())))
            fail("placeholders lists \"" + p + "\", which is not a quantity in this file");
    }
    bool any_placeholder = false;
    for (const auto& q : S.quantities) any_placeholder = any_placeholder || q.placeholder;
    if (j.has("contains_placeholders")) {
        const Json& c = j["contains_placeholders"];
        if (c.type() != Json::Type::Bool) fail("contains_placeholders must be true or false");
        if (c.boolean() != any_placeholder)
            fail(std::string("contains_placeholders is ") + (c.boolean() ? "true" : "false") + " but " +
                 (any_placeholder ? "a quantity is a placeholder" : "no quantity is a placeholder"));
    }
    if (S.produced_by.tool.empty()) fail("produced_by.tool is missing (a sidecar names the tool or person that wrote it)");
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
            if (is_gap && *q.value != D.gap_ueV)
                L.note += "; each qubit keeps its Rn_kOhm, so E_J and f01 scale with the gap (Ambegaokar-Baratoff)";
        }
        R.lines.push_back(L);
    }
    D.Q_dielectric = Q;
    D.gap_ueV = gap;
    return R;
}

// ------------------------------------------------------------------ materials apply --out (decision D2)
// The block `qpuforge materials apply` adds to the chip JSON it writes, so a reviewed chip names the sidecar its values
// came from. design_from_json() ignores it (it reads only the keys it knows), so the written chip loads anywhere a chip
// loads; save_design() does not carry it forward. Layout:
//   {"tool": "qpuforge materials apply", "schema": "forge.materials.v0", "sidecar": <file name as given>,
//    "sidecar_sha256": <hex of the sidecar file's bytes>, "produced_by": {tool, commit, potential_or_method, timestamp},
//    "system": [...], "loss_tangent_opt_in": bool, "gap_policy": "fixed Rn_kOhm (D1)",
//    "quantities": [{"quantity", "field", "applied", "before", "after", "provenance", "reference", "note"}, ...],
//    "previous": <the input chip's materials_source, when it had one>}
namespace materials_detail {

// SHA-256 (FIPS 180-4) of a byte string, lower-case hex.
inline std::string sha256_hex(const std::string& data) {
    static const std::uint32_t K[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    std::uint32_t h[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au, 0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    std::string m = data;
    const std::uint64_t bits = std::uint64_t(data.size()) * 8u;
    m.push_back(char(0x80));
    while (m.size() % 64 != 56) m.push_back('\0');
    for (int i = 7; i >= 0; --i) m.push_back(char((bits >> (8 * i)) & 0xffu));
    const auto rotr = [](std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); };
    for (std::size_t off = 0; off < m.size(); off += 64) {
        std::uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = std::uint32_t(static_cast<unsigned char>(m[off + 4 * i])) << 24 | std::uint32_t(static_cast<unsigned char>(m[off + 4 * i + 1])) << 16 |
                   std::uint32_t(static_cast<unsigned char>(m[off + 4 * i + 2])) << 8 | std::uint32_t(static_cast<unsigned char>(m[off + 4 * i + 3]));
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t t1 = hh + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i];
            const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    std::string out;
    char buf[9];
    for (std::uint32_t v : h) { std::snprintf(buf, sizeof buf, "%08x", unsigned(v)); out += buf; }
    return out;
}

inline std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

} // namespace materials_detail

// The materials_source block for a chip written from `S` (whose file had the bytes `sidecar_bytes`) with report `R`.
inline Json materials_source_json(const MaterialsSidecar& S, const MaterialsReport& R, const std::string& sidecar_name,
                                  const std::string& sidecar_bytes, bool use_loss_tangent) {
    Json j = Json::object();
    j.set("tool", "qpuforge materials apply");
    j.set("schema", S.schema);
    j.set("sidecar", sidecar_name);
    j.set("sidecar_sha256", materials_detail::sha256_hex(sidecar_bytes));
    Json p = Json::object();
    p.set("tool", S.produced_by.tool); p.set("commit", S.produced_by.commit);
    p.set("potential_or_method", S.produced_by.potential_or_method); p.set("timestamp", S.produced_by.timestamp);
    j.set("produced_by", p);
    Json sys = Json::array();
    for (const auto& el : S.system) sys.push(el);
    j.set("system", sys);
    j.set("loss_tangent_opt_in", use_loss_tangent);
    j.set("gap_policy", "fixed Rn_kOhm: E_J and f01 scale with junction_gap_ueV (decision D1)");
    Json qs = Json::array();
    for (std::size_t i = 0; i < R.lines.size(); ++i) {
        const auto& l = R.lines[i];
        const MaterialsQuantity* q = S.find(l.quantity);
        Json o = Json::object();
        o.set("quantity", l.quantity); o.set("field", l.field); o.set("applied", l.applied);
        o.set("before", l.before); o.set("after", l.after);
        o.set("provenance", q && !q->provenance.empty() ? Json(q->provenance) : Json());
        o.set("reference", q && !q->reference.empty() ? Json(q->reference) : Json());
        o.set("note", l.note);
        qs.push(o);
    }
    j.set("quantities", qs);
    return j;
}

// What materials_apply_file() wrote.
struct MaterialsApplyResult {
    ChipDesign design;          // the design as written
    MaterialsSidecar sidecar;
    MaterialsReport report;
    Json chip;                  // the chip JSON as written, with its materials_source block
};

// `qpuforge materials apply`: read chip_path, apply sidecar_path, write out_path (the chip as design_to_json() writes it,
// plus materials_source). Every check runs before out_path is opened; out_path must not name chip_path or sidecar_path.
inline MaterialsApplyResult materials_apply_file(const std::string& chip_path, const std::string& sidecar_path,
                                                 const std::string& out_path, const MaterialsApplyOptions& opt = {}) {
    if (out_path.empty()) throw std::runtime_error("materials apply: the output file name is empty");
    if (out_path == chip_path || out_path == sidecar_path)
        throw std::runtime_error("materials apply: --out " + out_path + " would overwrite an input file; write to a new file");
    const std::string chip_text = materials_detail::read_file(chip_path);
    const Json chip_in = Json::parse(chip_text);
    MaterialsApplyResult r;
    r.design = design_from_json(chip_in);
    const std::string bytes = materials_detail::read_file(sidecar_path);
    try {
        r.sidecar = materials_from_json(Json::parse(bytes));
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(sidecar_path + ": " + e.what());
    }
    r.report = apply_materials(r.sidecar, r.design, opt);
    r.chip = design_to_json(r.design);
    Json src = materials_source_json(r.sidecar, r.report, sidecar_path, bytes, opt.use_loss_tangent);
    if (chip_in.is_obj() && chip_in.has("materials_source")) src.set("previous", chip_in["materials_source"]);
    r.chip.set("materials_source", src);
    std::ofstream f(out_path);
    if (!f) throw std::runtime_error("cannot write " + out_path);
    r.chip.dump(f);
    f << "\n";
    f.close();
    if (!f) throw std::runtime_error("cannot write " + out_path);
    return r;
}

} // namespace qpuforge
