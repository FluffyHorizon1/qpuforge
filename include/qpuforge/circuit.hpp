// qpuforge — circuit IR, OpenQASM 2.0 reader, benchmark circuit builders
#pragma once
#include "core.hpp"
#include <cctype>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>

namespace qpuforge {

// Enumerator values are kept stable across releases: new ops are appended at the end (RESET was added in 0.3).
enum class Op { U, CX, CZ, SWAP, CP, CCX, MEASURE, BARRIER,
                // native (physical) ops produced by the compiler
                RZ, SX, X,
                // mid-circuit reset to |0> (0.3)
                RESET };

struct Gate {
    Op op{Op::U};
    std::vector<int> q;     // qubits (logical before compile, physical after)
    M2 u{};                 // for Op::U
    double theta{0};        // for RZ / CP
    int cbit{-1};           // for MEASURE
    std::string name;
};

struct Circuit {
    int n{0};               // qubits
    int nclbits{0};
    std::vector<Gate> gates;
    std::string name{"circuit"};

    void u(int q, const M2& m, const std::string& nm = "u") { gates.push_back({Op::U, {q}, m, 0, -1, nm}); }
    void h(int q) { u(q, ::qpuforge::gates::H(), "h"); }
    void x(int q) { u(q, ::qpuforge::gates::X(), "x"); }
    void cx(int c, int t) { gates.push_back({Op::CX, {c, t}, {}, 0, -1, "cx"}); }
    void cz(int a, int b) { gates.push_back({Op::CZ, {a, b}, {}, 0, -1, "cz"}); }
    void cp(int a, int b, double l) { gates.push_back({Op::CP, {a, b}, {}, l, -1, "cp"}); }
    void swap(int a, int b) { gates.push_back({Op::SWAP, {a, b}, {}, 0, -1, "swap"}); }
    void measure(int q, int c) { gates.push_back({Op::MEASURE, {q}, {}, 0, c, "measure"}); nclbits = std::max(nclbits, c + 1); }
    void measure_all() { for (int i = 0; i < n; ++i) measure(i, i); }
    void reset(int q) { gates.push_back({Op::RESET, {q}, {}, 0, -1, "reset"}); }
};

// ------------------------------------------------------------------ OpenQASM 2.0 (common subset)
namespace qasm_detail {
struct Expr {   // tiny recursive-descent evaluator: numbers, pi, + - * / ^, parentheses, unary minus
    const std::string& s; std::size_t i{0};
    explicit Expr(const std::string& str) : s(str) {}
    void ws() { while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i; }
    double parse() { double v = sum(); ws(); if (i != s.size()) throw std::runtime_error("qasm: bad expression '" + s + "'"); return v; }
    double sum() { double v = prod(); for (;;) { ws(); if (i < s.size() && s[i] == '+') { ++i; v += prod(); } else if (i < s.size() && s[i] == '-') { ++i; v -= prod(); } else return v; } }
    double prod() { double v = power(); for (;;) { ws(); if (i < s.size() && s[i] == '*') { ++i; v *= power(); } else if (i < s.size() && s[i] == '/') { ++i; v /= power(); } else return v; } }
    double power() { double v = unary(); ws(); if (i < s.size() && s[i] == '^') { ++i; v = std::pow(v, power()); } return v; }
    double unary() { ws(); if (i < s.size() && s[i] == '-') { ++i; return -unary(); } if (i < s.size() && s[i] == '+') { ++i; return unary(); } return atom(); }
    double atom() {
        ws();
        if (i < s.size() && s[i] == '(') { ++i; double v = sum(); ws(); if (i >= s.size() || s[i] != ')') throw std::runtime_error("qasm: missing )"); ++i; return v; }
        if (s.compare(i, 2, "pi") == 0) { i += 2; return PI; }
        const char* st = s.c_str() + i; char* en = nullptr;
        const double v = std::strtod(st, &en);
        if (en == st) throw std::runtime_error("qasm: bad number in '" + s + "'");
        i += std::size_t(en - st);
        return v;
    }
};
inline std::string trim(const std::string& x) {
    std::size_t a = 0, b = x.size();
    while (a < b && std::isspace(static_cast<unsigned char>(x[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(x[b - 1]))) --b;
    return x.substr(a, b - a);
}
inline std::vector<std::string> split(const std::string& x, char sep) {
    std::vector<std::string> out; std::string cur; int depth = 0;
    for (char c : x) {
        if (c == '(') ++depth; else if (c == ')') --depth;
        if (c == sep && depth == 0) { out.push_back(trim(cur)); cur.clear(); } else cur += c;
    }
    if (!trim(cur).empty()) out.push_back(trim(cur));
    return out;
}
} // namespace qasm_detail

inline Circuit parse_qasm(const std::string& text) {
    using namespace qasm_detail;
    Circuit C; C.name = "qasm";
    std::map<std::string, std::pair<int, int>> qregs, cregs;   // name -> (offset, size)
    // strip comments
    std::string src;
    {
        std::istringstream in(text); std::string line;
        while (std::getline(in, line)) { const auto p = line.find("//"); src += (p == std::string::npos ? line : line.substr(0, p)) + "\n"; }
    }
    auto resolve = [&](const std::string& arg, const std::map<std::string, std::pair<int, int>>& regs) -> std::vector<int> {
        const auto lb = arg.find('[');
        const std::string name = trim(arg.substr(0, lb));
        auto it = regs.find(name);
        if (it == regs.end()) throw std::runtime_error("qasm: unknown register '" + name + "'");
        if (lb == std::string::npos) { std::vector<int> v; for (int k = 0; k < it->second.second; ++k) v.push_back(it->second.first + k); return v; }
        const int idx = std::stoi(arg.substr(lb + 1, arg.find(']') - lb - 1));
        if (idx >= it->second.second) throw std::runtime_error("qasm: index out of range in " + arg);
        return {it->second.first + idx};
    };
    std::size_t start = 0;
    while (true) {
        const auto semi = src.find(';', start);
        if (semi == std::string::npos) break;
        std::string st = trim(src.substr(start, semi - start));
        start = semi + 1;
        if (st.empty()) continue;
        if (st.rfind("OPENQASM", 0) == 0 || st.rfind("include", 0) == 0) continue;
        if (st.rfind("qreg", 0) == 0 || st.rfind("creg", 0) == 0) {
            const bool isq = st[0] == 'q';
            const std::string body = trim(st.substr(4));
            const auto lb = body.find('[');
            const std::string nm = trim(body.substr(0, lb));
            const int sz = std::stoi(body.substr(lb + 1, body.find(']') - lb - 1));
            if (isq) { qregs[nm] = {C.n, sz}; C.n += sz; } else { cregs[nm] = {C.nclbits, sz}; C.nclbits += sz; }
            continue;
        }
        if (st.rfind("measure", 0) == 0) {
            const std::string body = trim(st.substr(7));
            const auto arrow = body.find("->");
            const auto qs = resolve(trim(body.substr(0, arrow)), qregs);
            const auto cs = resolve(trim(body.substr(arrow + 2)), cregs);
            if (qs.size() != cs.size()) throw std::runtime_error("qasm: measure size mismatch");
            for (std::size_t k = 0; k < qs.size(); ++k) C.gates.push_back({Op::MEASURE, {qs[k]}, {}, 0, cs[k], "measure"});
            continue;
        }
        if (st.rfind("barrier", 0) == 0) { C.gates.push_back({Op::BARRIER, {}, {}, 0, -1, "barrier"}); continue; }
        if (st.rfind("reset", 0) == 0) {
            for (int q : resolve(trim(st.substr(5)), qregs)) C.gates.push_back({Op::RESET, {q}, {}, 0, -1, "reset"});
            continue;
        }
        if (st.rfind("if", 0) == 0 || st.rfind("gate", 0) == 0)
            throw std::runtime_error("qasm: '" + st.substr(0, st.find(' ')) + "' not supported yet");
        // gate statement: name[(params)] args
        std::size_t k = 0;
        while (k < st.size() && (std::isalnum(static_cast<unsigned char>(st[k])) || st[k] == '_')) ++k;
        std::string name = st.substr(0, k);
        std::vector<double> params;
        std::string rest = trim(st.substr(k));
        if (!rest.empty() && rest[0] == '(') {
            int depth = 0; std::size_t e = 0;
            for (; e < rest.size(); ++e) { if (rest[e] == '(') ++depth; else if (rest[e] == ')' && --depth == 0) break; }
            for (const auto& p : split(rest.substr(1, e - 1), ',')) params.push_back(Expr(p).parse());
            rest = trim(rest.substr(e + 1));
        }
        std::vector<std::vector<int>> args;
        for (const auto& a : split(rest, ',')) args.push_back(resolve(a, qregs));
        std::size_t width = 1;
        for (auto& a : args) width = std::max(width, a.size());
        for (std::size_t w = 0; w < width; ++w) {
            auto A = [&](int i) { return args[i].size() == 1 ? args[i][0] : args[i][w]; };
            auto P = [&](int i) { if (int(params.size()) <= i) throw std::runtime_error("qasm: missing parameter for " + name); return params[i]; };
            using namespace gates;
            if (name == "id" || name == "i") continue;
            else if (name == "x") C.u(A(0), X(), "x");
            else if (name == "y") C.u(A(0), Y(), "y");
            else if (name == "z") C.u(A(0), Z(), "z");
            else if (name == "h") C.u(A(0), H(), "h");
            else if (name == "s") C.u(A(0), S(), "s");
            else if (name == "sdg") C.u(A(0), Sdg(), "sdg");
            else if (name == "t") C.u(A(0), T(), "t");
            else if (name == "tdg") C.u(A(0), Tdg(), "tdg");
            else if (name == "sx") C.u(A(0), SX(), "sx");
            else if (name == "sxdg") C.u(A(0), SXdg(), "sxdg");
            else if (name == "rx") C.u(A(0), RX(P(0)), "rx");
            else if (name == "ry") C.u(A(0), RY(P(0)), "ry");
            else if (name == "rz") C.u(A(0), RZ(P(0)), "rz");
            else if (name == "p" || name == "u1") C.u(A(0), gates::P(P(0)), "p");
            else if (name == "u2") C.u(A(0), U3(PI / 2, P(0), P(1)), "u2");
            else if (name == "u3" || name == "u" || name == "U") C.u(A(0), U3(P(0), P(1), P(2)), "u3");
            else if (name == "cx" || name == "CX") C.cx(A(0), A(1));
            else if (name == "cz") C.cz(A(0), A(1));
            else if (name == "swap") C.swap(A(0), A(1));
            else if (name == "cp" || name == "cu1") C.cp(A(0), A(1), P(0));
            else if (name == "ccx") C.gates.push_back({Op::CCX, {A(0), A(1), A(2)}, {}, 0, -1, "ccx"});
            else throw std::runtime_error("qasm: unsupported gate '" + name + "'");
        }
    }
    C.nclbits = std::max(C.nclbits, 0);
    return C;
}
inline Circuit load_qasm(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    Circuit c = parse_qasm(ss.str());
    c.name = path;
    return c;
}

// ------------------------------------------------------------------ benchmark circuits
inline Circuit ghz(int n) {
    Circuit C; C.n = n; C.name = "GHZ-" + std::to_string(n);
    C.h(0);
    for (int i = 0; i + 1 < n; ++i) C.cx(i, i + 1);
    C.measure_all();
    return C;
}

inline Circuit qft(int n) {
    Circuit C; C.n = n; C.name = "QFT-" + std::to_string(n);
    // prepare a non-trivial product input so the output is not uniform
    for (int i = 0; i < n; ++i) C.u(i, gates::RY(0.3 + 0.7 * i), "ry");
    for (int j = n - 1; j >= 0; --j) {
        C.h(j);
        for (int k = j - 1; k >= 0; --k) C.cp(k, j, PI / double(1 << (j - k)));
    }
    for (int i = 0; i < n / 2; ++i) C.swap(i, n - 1 - i);
    C.measure_all();
    return C;
}

// Bernstein–Vazirani: n data qubits + 1 ancilla; ideal output = secret
inline Circuit bernstein_vazirani(int n, unsigned long long secret) {
    Circuit C; C.n = n + 1; C.name = "BV-" + std::to_string(n);
    C.x(n); C.h(n);
    for (int i = 0; i < n; ++i) C.h(i);
    for (int i = 0; i < n; ++i) if ((secret >> i) & 1ULL) C.cx(i, n);
    for (int i = 0; i < n; ++i) C.h(i);
    for (int i = 0; i < n; ++i) C.measure(i, i);
    return C;
}

// Partition edges into matchings (greedy edge colouring) — used to lay out simultaneous 2q layers.
inline std::vector<std::vector<std::pair<int, int>>> edge_layers(const std::vector<std::pair<int, int>>& edges, int nq) {
    std::vector<std::vector<std::pair<int, int>>> layers;
    std::vector<std::vector<int>> used;   // per layer, per qubit
    for (auto e : edges) {
        std::size_t L = 0;
        for (; L < layers.size(); ++L) if (!used[L][e.first] && !used[L][e.second]) break;
        if (L == layers.size()) { layers.emplace_back(); used.emplace_back(nq, 0); }
        layers[L].push_back(e); used[L][e.first] = used[L][e.second] = 1;
    }
    return layers;
}

// Random "supremacy-style" circuit on physical couplers: cycles of random {sqrt X, sqrt Y, sqrt W} + CZ layer
inline Circuit random_xeb(int nq, const std::vector<std::pair<int, int>>& edges, int cycles, unsigned seed) {
    Circuit C; C.n = nq; C.name = "XEB-" + std::to_string(cycles);
    std::mt19937 rng(seed);
    const auto layers = edge_layers(edges, nq);
    std::vector<int> last(nq, -1);
    for (int cyc = 0; cyc < cycles; ++cyc) {
        for (int q = 0; q < nq; ++q) {
            int g;
            do { g = int(rng() % 3); } while (g == last[q]);
            last[q] = g;
            C.u(q, g == 0 ? gates::SX() : (g == 1 ? gates::SY() : gates::SW()), g == 0 ? "sx" : (g == 1 ? "sy" : "sw"));
        }
        if (!layers.empty()) for (auto [a, b] : layers[cyc % layers.size()]) C.cz(a, b);
    }
    for (int q = 0; q < nq; ++q) { int g; do { g = int(rng() % 3); } while (g == last[q]); C.u(q, g == 0 ? gates::SX() : (g == 1 ? gates::SY() : gates::SW()), "r"); }
    C.measure_all();
    return C;
}

} // namespace qpuforge
