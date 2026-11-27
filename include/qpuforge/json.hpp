// qpuforge — minimal JSON value, parser and writer (enough for chip design files)
#pragma once
#include <cmath>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qpuforge {

class Json {
public:
    enum class Type { Null, Bool, Num, Str, Arr, Obj };
    Json() = default;
    Json(std::nullptr_t) {}
    Json(bool b) : t_(Type::Bool), b_(b) {}
    Json(double v) : t_(Type::Num), n_(v) {}
    Json(int v) : t_(Type::Num), n_(v) {}
    Json(const char* s) : t_(Type::Str), s_(s) {}
    Json(std::string s) : t_(Type::Str), s_(std::move(s)) {}
    static Json array() { Json j; j.t_ = Type::Arr; return j; }
    static Json object() { Json j; j.t_ = Type::Obj; return j; }

    Type type() const { return t_; }
    bool is_null() const { return t_ == Type::Null; }
    bool is_obj() const { return t_ == Type::Obj; }
    bool is_arr() const { return t_ == Type::Arr; }
    double num() const { if (t_ != Type::Num) throw std::runtime_error("json: not a number"); return n_; }
    bool boolean() const { if (t_ != Type::Bool) throw std::runtime_error("json: not a bool"); return b_; }
    const std::string& str() const { if (t_ != Type::Str) throw std::runtime_error("json: not a string"); return s_; }

    // arrays
    std::size_t size() const { return t_ == Type::Arr ? a_.size() : (t_ == Type::Obj ? o_.size() : 0); }
    const Json& operator[](std::size_t i) const { return a_.at(i); }
    void push(Json v) { if (t_ != Type::Arr) { t_ = Type::Arr; } a_.push_back(std::move(v)); }

    // objects (insertion-ordered)
    bool has(const std::string& k) const { for (auto& kv : o_) if (kv.first == k) return true; return false; }
    const Json& operator[](const std::string& k) const {
        static const Json null_json;
        for (auto& kv : o_) if (kv.first == k) return kv.second;
        return null_json;
    }
    Json& set(const std::string& k, Json v) {
        if (t_ != Type::Obj) t_ = Type::Obj;
        for (auto& kv : o_) if (kv.first == k) { kv.second = std::move(v); return kv.second; }
        o_.emplace_back(k, std::move(v));
        return o_.back().second;
    }
    const std::vector<std::pair<std::string, Json>>& items() const { return o_; }

    double get(const std::string& k, double def) const { const Json& v = (*this)[k]; return v.t_ == Type::Num ? v.n_ : def; }
    bool get(const std::string& k, bool def) const { const Json& v = (*this)[k]; return v.t_ == Type::Bool ? v.b_ : def; }
    std::string get(const std::string& k, const char* def) const { const Json& v = (*this)[k]; return v.t_ == Type::Str ? v.s_ : std::string(def); }

    // ---------------------------------------------------------------- parse
    static constexpr int kMaxNest = 256;
    static Json parse(const std::string& text) {
        std::size_t i = 0;
        Json v = parse_value(text, i, 0);
        skip_ws(text, i);
        if (i != text.size()) throw std::runtime_error("json: trailing characters at " + std::to_string(i));
        return v;
    }

    // ---------------------------------------------------------------- write
    void dump(std::ostream& os, int indent = 2, int level = 0) const {
        const std::string pad(std::size_t(indent * (level + 1)), ' '), pad0(std::size_t(indent * level), ' ');
        switch (t_) {
        case Type::Null: os << "null"; break;
        case Type::Bool: os << (b_ ? "true" : "false"); break;
        case Type::Num: {   // JSON has no inf/nan: write null so the file always parses
            if (!std::isfinite(n_)) { os << "null"; break; }
            std::ostringstream s; s.precision(12); s << n_; os << s.str(); break;
        }
        case Type::Str: write_str(os, s_); break;
        case Type::Arr: {
            bool flat = true;
            for (auto& v : a_) if (v.t_ == Type::Arr || v.t_ == Type::Obj) flat = false;
            os << "[";
            for (std::size_t k = 0; k < a_.size(); ++k) {
                if (!flat) os << "\n" << pad; else if (k) os << " ";
                a_[k].dump(os, indent, level + 1);
                if (k + 1 < a_.size()) os << ",";
            }
            if (!flat && !a_.empty()) os << "\n" << pad0;
            os << "]";
            break;
        }
        case Type::Obj: {
            bool flat = o_.size() <= 6;
            for (auto& kv : o_) if (kv.second.t_ == Type::Arr || kv.second.t_ == Type::Obj) flat = false;
            os << "{";
            for (std::size_t k = 0; k < o_.size(); ++k) {
                if (!flat) os << "\n" << pad; else if (k) os << " ";
                write_str(os, o_[k].first); os << ": ";
                o_[k].second.dump(os, indent, level + 1);
                if (k + 1 < o_.size()) os << ",";
            }
            if (!flat && !o_.empty()) os << "\n" << pad0;
            os << "}";
            break;
        }
        }
    }

private:
    Type t_{Type::Null};
    bool b_{false};
    double n_{0};
    std::string s_;
    std::vector<Json> a_;
    std::vector<std::pair<std::string, Json>> o_;

    static void write_str(std::ostream& os, const std::string& s) {
        os << '"';
        for (unsigned char c : s) {
            switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\b': os << "\\b"; break;
            case '\f': os << "\\f"; break;
            case '\n': os << "\\n"; break;
            case '\r': os << "\\r"; break;
            case '\t': os << "\\t"; break;
            default:
                if (c < 0x20) {
                    const char* hex = "0123456789abcdef";
                    os << "\\u00" << hex[c >> 4] << hex[c & 0xf];
                } else os << char(c);
            }
        }
        os << '"';
    }
    static void skip_ws(const std::string& t, std::size_t& i) {
        while (i < t.size()) {
            if (std::isspace(static_cast<unsigned char>(t[i]))) { ++i; continue; }
            if (t[i] == '/' && i + 1 < t.size() && t[i + 1] == '/') { while (i < t.size() && t[i] != '\n') ++i; continue; }  // allow // comments
            break;
        }
    }
    static Json parse_value(const std::string& t, std::size_t& i, int depth) {
        if (depth > kMaxNest) throw std::runtime_error("json: nesting too deep");
        skip_ws(t, i);
        if (i >= t.size()) throw std::runtime_error("json: unexpected end");
        const char c = t[i];
        if (c == '{') {
            Json o = object(); ++i; skip_ws(t, i);
            if (t[i] == '}') { ++i; return o; }
            while (true) {
                skip_ws(t, i);
                if (t[i] != '"') throw std::runtime_error("json: expected key at " + std::to_string(i));
                std::string k = parse_string(t, i);
                skip_ws(t, i);
                if (t[i] != ':') throw std::runtime_error("json: expected ':' at " + std::to_string(i));
                ++i;
                o.set(k, parse_value(t, i, depth + 1));
                skip_ws(t, i);
                if (t[i] == ',') { ++i; continue; }
                if (t[i] == '}') { ++i; return o; }
                throw std::runtime_error("json: expected ',' or '}' at " + std::to_string(i));
            }
        }
        if (c == '[') {
            Json a = array(); ++i; skip_ws(t, i);
            if (t[i] == ']') { ++i; return a; }
            while (true) {
                a.push(parse_value(t, i, depth + 1));
                skip_ws(t, i);
                if (t[i] == ',') { ++i; continue; }
                if (t[i] == ']') { ++i; return a; }
                throw std::runtime_error("json: expected ',' or ']' at " + std::to_string(i));
            }
        }
        if (c == '"') return Json(parse_string(t, i));
        if (t.compare(i, 4, "true") == 0) { i += 4; return Json(true); }
        if (t.compare(i, 5, "false") == 0) { i += 5; return Json(false); }
        if (t.compare(i, 4, "null") == 0) { i += 4; return Json(); }
        // RFC 8259 number: optional minus, integer, optional frac, optional exp. No nan/inf/hex/leading +.
        const std::size_t n0 = i;
        if (i < t.size() && t[i] == '-') ++i;
        if (i >= t.size() || !std::isdigit(static_cast<unsigned char>(t[i])))
            throw std::runtime_error("json: bad value at " + std::to_string(n0));
        if (t[i] == '0') {
            ++i;
        } else {
            while (i < t.size() && std::isdigit(static_cast<unsigned char>(t[i]))) ++i;
        }
        if (i < t.size() && t[i] == '.') {
            ++i;
            if (i >= t.size() || !std::isdigit(static_cast<unsigned char>(t[i])))
                throw std::runtime_error("json: bad number at " + std::to_string(n0));
            while (i < t.size() && std::isdigit(static_cast<unsigned char>(t[i]))) ++i;
        }
        if (i < t.size() && (t[i] == 'e' || t[i] == 'E')) {
            ++i;
            if (i < t.size() && (t[i] == '+' || t[i] == '-')) ++i;
            if (i >= t.size() || !std::isdigit(static_cast<unsigned char>(t[i])))
                throw std::runtime_error("json: bad number at " + std::to_string(n0));
            while (i < t.size() && std::isdigit(static_cast<unsigned char>(t[i]))) ++i;
        }
        char* end = nullptr;
        const double v = std::strtod(t.c_str() + n0, &end);
        if (end != t.c_str() + i || !std::isfinite(v))
            throw std::runtime_error("json: bad number at " + std::to_string(n0));
        return Json(v);
    }
    static int hex4(const std::string& t, std::size_t i) {
        if (i + 4 > t.size()) throw std::runtime_error("json: truncated \\u escape at " + std::to_string(i));
        int v = 0;
        for (int k = 0; k < 4; ++k) {
            const unsigned char c = static_cast<unsigned char>(t[i + k]);
            v <<= 4;
            if (c >= '0' && c <= '9') v += c - '0';
            else if (c >= 'a' && c <= 'f') v += c - 'a' + 10;
            else if (c >= 'A' && c <= 'F') v += c - 'A' + 10;
            else throw std::runtime_error("json: invalid \\u escape at " + std::to_string(i));
        }
        return v;
    }
    static void append_utf8(std::string& s, std::uint32_t cp) {
        if (cp < 0x80) s += char(cp);
        else if (cp < 0x800) { s += char(0xC0 | (cp >> 6)); s += char(0x80 | (cp & 0x3F)); }
        else if (cp < 0x10000) {
            s += char(0xE0 | (cp >> 12)); s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F));
        } else {
            s += char(0xF0 | (cp >> 18)); s += char(0x80 | ((cp >> 12) & 0x3F));
            s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F));
        }
    }
    static std::string parse_string(const std::string& t, std::size_t& i) {
        std::string s; ++i;
        while (i < t.size() && t[i] != '"') {
            const unsigned char c = static_cast<unsigned char>(t[i]);
            if (c < 0x20) throw std::runtime_error("json: unescaped control in string at " + std::to_string(i));
            if (t[i] == '\\') {
                if (i + 1 >= t.size()) throw std::runtime_error("json: unterminated string");
                ++i;
                switch (t[i]) {
                case '"': case '\\': case '/': s += t[i]; break;
                case 'b': s += '\b'; break;
                case 'f': s += '\f'; break;
                case 'n': s += '\n'; break;
                case 'r': s += '\r'; break;
                case 't': s += '\t'; break;
                case 'u': {
                    int u = hex4(t, i + 1); i += 4;
                    if (u >= 0xD800 && u <= 0xDBFF) {   // surrogate pair
                        if (i + 6 < t.size() && t[i + 1] == '\\' && t[i + 2] == 'u') {
                            const int lo = hex4(t, i + 3);
                            if (lo >= 0xDC00 && lo <= 0xDFFF) {
                                append_utf8(s, 0x10000u + (std::uint32_t(u - 0xD800) << 10) + std::uint32_t(lo - 0xDC00));
                                i += 6;
                                break;
                            }
                        }
                        throw std::runtime_error("json: unpaired surrogate at " + std::to_string(i));
                    }
                    if (u >= 0xDC00 && u <= 0xDFFF)
                        throw std::runtime_error("json: unpaired surrogate at " + std::to_string(i));
                    append_utf8(s, std::uint32_t(u));
                    break;
                }
                default: throw std::runtime_error(std::string("json: unknown escape \\") + t[i]);
                }
            } else s += t[i];
            ++i;
        }
        if (i >= t.size()) throw std::runtime_error("json: unterminated string");
        ++i;
        return s;
    }
};

} // namespace qpuforge
