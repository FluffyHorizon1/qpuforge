// qpuforge — minimal GDSII stream writer (boundaries + text, one or more structures)
//   database unit 1 nm, user unit 1 um. Layers used by the chip assembler:
//     1 etch (ground-plane cut-outs minus metal islands)     2 Josephson junctions / SQUID loops
//     3 flip-chip bump / coupling pads                       10 chip outline           20 labels
#pragma once
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qpuforge {
namespace gds {

struct Polygon { int layer{1}, datatype{0}; std::vector<std::pair<double, double>> pts; };   // um, open ring
struct Text { int layer{20}; double x{0}, y{0}; std::string s; };
struct Structure { std::string name; std::vector<Polygon> polys; std::vector<Text> texts; };

namespace detail {
inline void u16(std::vector<std::uint8_t>& b, std::uint16_t v) { b.push_back(std::uint8_t(v >> 8)); b.push_back(std::uint8_t(v & 0xff)); }
inline void i32(std::vector<std::uint8_t>& b, std::int32_t v) { const auto u = std::uint32_t(v); for (int s = 24; s >= 0; s -= 8) b.push_back(std::uint8_t((u >> s) & 0xff)); }
// GDSII 8-byte real: sign bit, 7-bit excess-64 base-16 exponent, 56-bit mantissa
inline void real8(std::vector<std::uint8_t>& b, double v) {
    std::uint64_t out = 0;
    if (v != 0) {
        const bool neg = v < 0; v = std::abs(v);
        int e = 0;
        while (v >= 1) { v /= 16; ++e; }
        while (v < 1.0 / 16) { v *= 16; --e; }
        const auto mant = std::uint64_t(std::llround(v * std::pow(2.0, 56)));
        out = (std::uint64_t(neg) << 63) | (std::uint64_t(e + 64) << 56) | (mant & ((std::uint64_t(1) << 56) - 1));
    }
    for (int s = 56; s >= 0; s -= 8) b.push_back(std::uint8_t((out >> s) & 0xff));
}
// GDSII record length is an unsigned 16-bit byte count including the 4-byte header and must be even, so a record holds
// at most 65530 bytes of data (e.g. a BOUNDARY of at most 8190 vertices plus the closing point). Refuse anything longer
// instead of silently wrapping the length field.
inline constexpr std::size_t kMaxRecordData = 65530;
inline void record(std::vector<std::uint8_t>& b, std::uint8_t type, std::uint8_t dtype, const std::vector<std::uint8_t>& data) {
    if (data.size() > kMaxRecordData || data.size() % 2)
        throw std::runtime_error("gds: record type " + std::to_string(int(type)) + " has " + std::to_string(data.size()) +
                                 " data bytes (must be even and <= " + std::to_string(kMaxRecordData) + "); split the polygon");
    u16(b, std::uint16_t(4 + data.size())); b.push_back(type); b.push_back(dtype); b.insert(b.end(), data.begin(), data.end());
}
// database units (nm) of a coordinate in um, range-checked against the signed 32-bit XY field
inline std::int32_t dbu(double um) {
    const double v = std::round(um * 1000);
    if (!(std::abs(v) <= 2147483647.0)) throw std::runtime_error("gds: coordinate " + std::to_string(um) + " um out of int32 range");
    return std::int32_t(v);
}
inline std::vector<std::uint8_t> ascii(const std::string& s) { std::vector<std::uint8_t> d(s.begin(), s.end()); if (d.size() % 2) d.push_back(0); return d; }
inline std::vector<std::uint8_t> shorts(std::initializer_list<int> v) { std::vector<std::uint8_t> d; for (int x : v) u16(d, std::uint16_t(x)); return d; }
// modification / access time (UTC) for BGNLIB / BGNSTR: now, or SOURCE_DATE_EPOCH if set (reproducible output)
inline std::vector<std::uint8_t> date() {
    std::time_t t = std::time(nullptr);
    if (const char* e = std::getenv("SOURCE_DATE_EPOCH")) t = std::time_t(std::strtoll(e, nullptr, 10));
    const std::tm u = *std::gmtime(&t);
    const int Y = u.tm_year + 1900, M = u.tm_mon + 1;
    return shorts({Y, M, u.tm_mday, u.tm_hour, u.tm_min, u.tm_sec, Y, M, u.tm_mday, u.tm_hour, u.tm_min, u.tm_sec});
}
} // namespace detail

inline void write(const std::string& path, const std::string& libname, const std::vector<Structure>& cells) {
    using namespace detail;
    std::vector<std::uint8_t> b;
    record(b, 0x00, 0x02, shorts({600}));                  // HEADER
    record(b, 0x01, 0x02, date());                         // BGNLIB
    record(b, 0x02, 0x06, ascii(libname));                 // LIBNAME
    { std::vector<std::uint8_t> d; real8(d, 1e-3); real8(d, 1e-9); record(b, 0x03, 0x05, d); }   // UNITS
    for (const auto& c : cells) {
        record(b, 0x05, 0x02, date());                     // BGNSTR
        record(b, 0x06, 0x06, ascii(c.name));              // STRNAME
        for (const auto& p : c.polys) {
            if (p.pts.size() < 3) continue;
            record(b, 0x08, 0x00, {});                     // BOUNDARY
            record(b, 0x0D, 0x02, shorts({p.layer}));      // LAYER
            record(b, 0x0E, 0x02, shorts({p.datatype}));   // DATATYPE
            std::vector<std::uint8_t> d;
            for (const auto& [x, y] : p.pts) { i32(d, dbu(x)); i32(d, dbu(y)); }
            i32(d, dbu(p.pts[0].first)); i32(d, dbu(p.pts[0].second));
            record(b, 0x10, 0x03, d);                      // XY
            record(b, 0x11, 0x00, {});                     // ENDEL
        }
        for (const auto& t : c.texts) {
            record(b, 0x0C, 0x00, {});                     // TEXT
            record(b, 0x0D, 0x02, shorts({t.layer}));
            record(b, 0x16, 0x02, shorts({0}));            // TEXTTYPE
            std::vector<std::uint8_t> d; i32(d, dbu(t.x)); i32(d, dbu(t.y));
            record(b, 0x10, 0x03, d);
            record(b, 0x19, 0x06, ascii(t.s));             // STRING
            record(b, 0x11, 0x00, {});
        }
        record(b, 0x07, 0x00, {});                         // ENDSTR
    }
    record(b, 0x04, 0x00, {});                             // ENDLIB
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(b.data()), std::streamsize(b.size()));
}

} // namespace gds
} // namespace qpuforge
