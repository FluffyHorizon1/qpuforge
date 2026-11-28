// qpuforge CLI — design, analyse and emulate superconducting QPU chips
#include "qpuforge/qpuforge.hpp"
#include "qpuforge/layout.hpp"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace qpuforge;

#ifndef QPUFORGE_VERSION
#define QPUFORGE_VERSION "unknown"
#endif

// ================================================================== usage
// One block per command line of the usage text, tagged with the command(s) it documents ("*": every command). `qpuforge`
// without arguments prints every block; a usage error prints the blocks of its own command.
struct UsageBlock { const char* commands; const char* text; };

static const UsageBlock kUsageCommands[] = {
    {"new",
     "  qpuforge new <chip.json> [--grid RxC | --line N | --ring N] [--f-hi GHz] [--f-lo GHz]\n"
     "                [--alpha GHz] [--g-MHz g] [--arch fixed|tunable] [--Q Q] [--plan 2|4|8]\n"
     "                [--stagger-MHz s] [--name NAME] [--feedlines N] [--kappa-MHz k] [--chi-MHz chi]\n"
     "                [--purcell f] [--readout-ns t] [--reset-ns t] [--reset-err p] [--err2q p]\n"
     "                [--residual-zz-kHz z] [--active-reset]\n"},
    {"analyze",
     "  qpuforge analyze <chip.json> [--materials <sidecar.json> [--loss-tangent]]\n"
     "                device physics report; --materials first applies a forge.materials.v0 sidecar (atomforge),\n"
     "                --loss-tangent lets it set Q_dielectric = 1/loss_tangent when Q_dielectric is null\n"},
    {"run", "  qpuforge run <chip.json> <circuit> [--traj N] [--shots N] [--seed S] [--off SRC,..] [--top K]\n"},
    {"budget", "  qpuforge budget <chip.json> <circuit> [--traj N] [--seed S]   error budget by noise source\n"},
    {"fab", "  qpuforge fab <chip.json> [--samples N] [--seed S] fabrication Monte Carlo (junction spread)\n"},
    {"qec",
     "  qpuforge qec <chip.json> [--rounds 2,4,6,8,10] [--shots N] [--seed S] [--basis Z|X|both] [--statevector N]\n"
     "                [--off SRC,..] [--decoder exact|uf|uf-stream|legacy|greedy [--window R] [--commit R]]\n"
     "                [--leakage-removal lrr|dqlr|lrr+dqlr|none]\n"
     "                [--cache file] [--herald [--herald-fidelity F] [--herald-false P] [--herald-factor X]\n"
     "                [--herald-rule footprint|detector]] [--bursts [--burst-rate Hz] [--burst-tau-ms T] [--burst-t1-us T]\n"
     "                [--burst-tphi-us T] [--veto [--veto-window R] [--veto-radius L] [--veto-false P] [--impact-p P]\n"
     "                [--burst-radius L]]]\n"
     "                surface-code memory experiment: logical error per round (needs a surface layout);\n"
     "                --bursts adds the logical-error floor set by high-energy impacts; with --veto, impacts are also\n"
     "                injected per shot and detected by a sliding-window veto (residual floor, discard fraction);\n"
     "                --decoder legacy reproduces the pre-v0.4 decoder (greedy above 22 defects) on the same syndromes;\n"
     "                --decoder uf-stream: union-find over a sliding window of --window rounds (6), committing --commit (3);\n"
     "                --herald: 3-state readout heralds leaked measure qubits; decoded with and without them\n"},
    {"synth-chip",
     "  qpuforge synth-chip <chip.json> <data.json> [--seed S] [--shots N] [--rounds 5,15,25] [--z-only]\n"
     "                [--hidden-dephasing-us T] [--cache file]   virtual chip -> characterisation file (to test the chain)\n"},
    {"hwgate",
     "  qpuforge hwgate <chip.json> <data.json> [--phase alpha|beta] [--shots N] [--seed S] [--out report.json]\n"
     "                [--cache file]   roadmap hardware gate from a characterisation file, incl. model vs measurement\n"},
    {"layout", "  qpuforge layout <chip.json> <out.svg>             chip floor-plan\n"},
    {"cz",
     "  qpuforge cz <chip.json> [--out cal.json] [--cache file] [--g-qc MHz] [--f-on GHz] [--rise ns]\n"
     "                pulse-level CZ calibration of every pair (physical tunable coupler)\n"},
    {"gds",
     "  qpuforge gds <chip.json> <out.gds> [--cells cells.json] [--synthesize [--from cells0.json]] [--cell NAME]\n"
     "                physical layout: MoM capacitance synthesis of the cells, chip assembly, DRC, GDSII\n"},
    {"new", "  qpuforge new <chip.json> --preset alpha1|alpha1v11|beta1|gamma1|gen0 [--distance D]   named designs\n"},
    {"new", "  qpuforge new <chip.json> --surface D [--arch tunable] ...        surface-code chip\n"},
    {"", "  qpuforge --version | -V                           print the version and exit (first argument only)\n"},
};

static const UsageBlock kUsageNotes[] = {
    {"run budget", "  <circuit> = file.qasm | ghz:N | qft:N | bv:N | xeb:CYCLES\n"},
    {"run qec", "  noise sources: T1 Tphi thermal zz exchange control readout reset leakage (all: every source)\n"},
    {"qec cz synth-chip hwgate", "  qec/cz/synth-chip/hwgate on a chip with a physical coupler calibrate CZs first (--cache file reuses them)\n"},
    {"*", "  each command takes only the options shown for it; an unknown option, a missing value or a bad number is an error\n"},
    {"*", "  exit status: 0 done; 1 usage error; 2 error while running (also: gds layout not DRC-clean); 3 hwgate: gate not passed\n"},
};

static bool documents(const char* commands, const std::string& cmd) {
    std::istringstream ss(commands);
    std::string w;
    while (ss >> w) if (w == "*" || w == cmd) return true;
    return false;
}

// cmd empty: the full usage (stdout for `qpuforge` alone); otherwise the usage of one command (after a usage error)
static void usage(std::FILE* out = stdout, const std::string& cmd = "") {
    std::fputs(cmd.empty() ? "qpuforge — superconducting QPU chip design & whole-processor emulation\n\n" : "usage:\n", out);
    for (const auto& b : kUsageCommands) if (cmd.empty() || documents(b.commands, cmd)) std::fputs(b.text, out);
    std::fputs("\n", out);
    for (const auto& b : kUsageNotes) if (cmd.empty() || documents(b.commands, cmd)) std::fputs(b.text, out);
    if (!cmd.empty()) std::fputs("  (qpuforge without arguments lists every command)\n", out);
}

// ================================================================== command line
// Every command declares the options it reads and the positional arguments it takes (commands(), after the cmd_*
// functions). Anything else is a usage error: an unknown or misspelt option, an option given twice or without its value,
// a value that does not parse or is out of range, a missing or an extra positional argument. Values are checked before
// the command reads any file. A usage error prints "qpuforge <command>: <what>" and the command's usage and exits 1;
// an error while running prints "error: <what>" and exits 2 (main).
//
// The CLI checks syntax, counts, probabilities and choices. Design rules (an odd surface-code distance, an even ring)
// stay with the library, which throws (exit 2).
struct UsageError : std::runtime_error { using std::runtime_error::runtime_error; };

enum class Kind { Flag, Text, Int, Seed, Real, Choice, IntList, ChoiceList, Grid };

struct Option {
    const char* name = "";
    Kind kind = Kind::Flag;
    long long min = 0, max = 0;                          // Int, IntList: accepted range
    unsigned long long seed_max = 0;                     // Seed: largest value (the width of the seed it sets)
    double lo = -HUGE_VAL, hi = HUGE_VAL;                // Real: accepted range (the value must also be finite)
    const char* choices = "";                            // Choice, ChoiceList: accepted values, '|'-separated
};

static Option make_option(const char* name, Kind kind) { Option o; o.name = name; o.kind = kind; return o; }
static Option flag(const char* name) { return make_option(name, Kind::Flag); }
static Option text(const char* name) { return make_option(name, Kind::Text); }
static Option integer(const char* name, long long min, long long max = std::numeric_limits<int>::max()) {
    Option o = make_option(name, Kind::Int); o.min = min; o.max = max; return o;
}
static Option int_list(const char* name, long long min) {
    Option o = make_option(name, Kind::IntList); o.min = min; o.max = std::numeric_limits<int>::max(); return o;
}
static Option seed(const char* name, unsigned long long max = std::numeric_limits<unsigned>::max()) {
    Option o = make_option(name, Kind::Seed); o.seed_max = max; return o;
}
static Option real(const char* name, double lo = -HUGE_VAL, double hi = HUGE_VAL) {
    Option o = make_option(name, Kind::Real); o.lo = lo; o.hi = hi; return o;
}
static Option choice(const char* name, const char* choices) { Option o = make_option(name, Kind::Choice); o.choices = choices; return o; }
static Option choice_list(const char* name, const char* choices) { Option o = make_option(name, Kind::ChoiceList); o.choices = choices; return o; }
static Option grid(const char* name) { Option o = make_option(name, Kind::Grid); o.min = 1; o.max = std::numeric_limits<int>::max(); return o; }

class Args;
struct Command {
    const char* name;
    int (*run)(const Args&);
    std::vector<const char*> positionals;   // all required, e.g. "<chip.json>"
    std::vector<Option> options;
    const Option* find(const std::string& opt) const {
        for (const auto& o : options) if (opt == o.name) return &o;
        return nullptr;
    }
};
static const std::vector<Command>& commands();

static std::vector<std::string> split_fields(const std::string& s, char sep) {   // keeps empty fields: "a,,b" -> {a, "", b}
    std::vector<std::string> out(1);
    for (char ch : s) { if (ch == sep) out.emplace_back(); else out.back() += ch; }
    return out;
}

static std::string join_words(const std::vector<std::string>& v, const char* last) {   // {a, b, c} -> "a, b <last> c"
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i) s += std::string(i == 0 ? "" : i + 1 == v.size() ? last : ", ") + v[i];
    return s;
}

// optimal-string-alignment distance: insertions, deletions, substitutions and swaps of adjacent characters
static std::size_t edit_distance(const std::string& a, const std::string& b) {
    std::vector<std::vector<std::size_t>> d(a.size() + 1, std::vector<std::size_t>(b.size() + 1));
    for (std::size_t i = 0; i <= a.size(); ++i) d[i][0] = i;
    for (std::size_t j = 0; j <= b.size(); ++j) d[0][j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i)
        for (std::size_t j = 1; j <= b.size(); ++j) {
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + std::size_t(a[i - 1] != b[j - 1])});
            if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1]) d[i][j] = std::min(d[i][j], d[i - 2][j - 2] + 1);
        }
    return d[a.size()][b.size()];
}

// The candidate `word` most likely meant, or "" when none is close. Case, leading dashes and '_' for '-' are ignored
// ("-shots", "--SHOTS" -> "--shots"); otherwise up to one edit for names of 3-5 characters and two for longer ones,
// or a unique prefix of at least 3 characters ("--hidden" -> "--hidden-dephasing-us").
static std::string suggest(const std::string& word, const std::vector<std::string>& candidates) {
    auto norm = [](std::string s) {
        s.erase(0, s.find_first_not_of('-'));
        for (auto& ch : s) ch = ch == '_' ? '-' : char(std::tolower(static_cast<unsigned char>(ch)));
        return s;
    };
    const std::string w = norm(word);
    if (w.empty()) return "";
    std::string best, prefix;
    std::size_t best_d = std::numeric_limits<std::size_t>::max(), prefixes = 0;
    for (const auto& c : candidates) {
        const std::string n = norm(c);
        if (n == w) return c;
        const std::size_t tol = n.size() <= 2 ? 0 : n.size() <= 5 ? 1 : 2;
        const std::size_t d = edit_distance(w, n);
        if (d <= tol && d < best_d) { best_d = d; best = c; }
        if (w.size() >= 3 && n.size() > w.size() && n.compare(0, w.size(), w) == 0) { ++prefixes; prefix = c; }
    }
    // A typed prefix beats a two-edit guess: "--burst-rad" is --burst-radius, not --burst-rate, and "--herald-r" is
    // --herald-rule, not the flag --herald. An ambiguous prefix ("--herald-f") gets no suggestion; the usage lists them.
    if (prefixes >= 1 && best_d >= 2) return prefixes == 1 ? prefix : "";
    if (!best.empty()) return best;
    return prefixes == 1 ? prefix : "";
}

// "--name", or one dash and a letter ("-V", a misspelt "-shots"). "-" and negative numbers ("-0.22") are not options.
static bool is_option(const std::string& t) {
    return t.size() >= 2 && t[0] == '-' && (t[1] == '-' || std::isalpha(static_cast<unsigned char>(t[1])));
}

static const char* skip_plus(const char* b, const char* e) {   // "+5" reads as 5, as std::stoi / std::stoul took it
    return e - b >= 2 && *b == '+' && std::isdigit(static_cast<unsigned char>(b[1])) ? b + 1 : b;
}

static long long parse_int(const Option& o, const std::string& v) {
    const char* e = v.data() + v.size();
    const char* b = skip_plus(v.data(), e);
    long long x = 0;
    const auto r = std::from_chars(b, e, x);
    const bool digits = b != e && r.ptr == e && r.ec != std::errc::invalid_argument;
    if (digits && r.ec == std::errc() && x >= o.min && x <= o.max) return x;
    const bool too_large = digits && (r.ec == std::errc::result_out_of_range ? v[0] != '-' : x > o.max);
    const std::string range = too_large ? "from " + std::to_string(o.min) + " to " + std::to_string(o.max) : ">= " + std::to_string(o.min);
    throw UsageError(std::string(o.name) + " expects an integer " + range + ", got '" + v + "'");
}

static unsigned long long parse_seed(const Option& o, const std::string& v) {
    const char* e = v.data() + v.size();
    const char* b = skip_plus(v.data(), e);
    unsigned long long x = 0;
    const auto r = std::from_chars(b, e, x);   // no sign: "-1" is refused, not wrapped to 2^64 - 1 as by std::stoul
    if (b != e && r.ptr == e && r.ec == std::errc() && x <= o.seed_max) return x;
    throw UsageError(std::string(o.name) + " expects an integer from 0 to " + std::to_string(o.seed_max) + ", got '" + v + "'");
}

static double parse_real(const Option& o, const std::string& v) {
    char* end = nullptr;
    const double x = v.empty() || std::isspace(static_cast<unsigned char>(v[0])) ? std::nan("") : std::strtod(v.c_str(), &end);
    if (end == v.c_str() + v.size() && std::isfinite(x) && x >= o.lo && x <= o.hi) return x;
    std::ostringstream range;
    if (std::isfinite(o.lo) && std::isfinite(o.hi)) range << " from " << o.lo << " to " << o.hi;
    else if (std::isfinite(o.lo)) range << " >= " << o.lo;
    else if (std::isfinite(o.hi)) range << " <= " << o.hi;
    throw UsageError(std::string(o.name) + " expects a number" + range.str() + ", got '" + v + "'");
}

static void check_choice(const Option& o, const std::string& v) {
    const auto cs = split_fields(o.choices, '|');
    if (std::find(cs.begin(), cs.end(), v) != cs.end()) return;
    const auto s = suggest(v, cs);
    throw UsageError(std::string(o.name) + " expects " + join_words(cs, " or ") + ", got '" + v + "'" + (s.empty() ? "" : " (did you mean '" + s + "'?)"));
}

static std::vector<std::string> parse_choice_list(const Option& o, const std::string& v) {
    const auto cs = split_fields(o.choices, '|');
    const auto items = split_fields(v, ',');
    for (const auto& it : items) {
        if (std::find(cs.begin(), cs.end(), it) != cs.end()) continue;
        if (it.empty()) throw UsageError(std::string(o.name) + " expects a comma-separated list, got '" + v + "'");
        const auto s = suggest(it, cs);
        throw UsageError(std::string(o.name) + ": '" + it + "' is not one of " + join_words(cs, ", ") + (s.empty() ? "" : " (did you mean '" + s + "'?)"));
    }
    return items;
}

static std::vector<int> parse_int_list(const Option& o, const std::string& v) {
    std::vector<int> out;
    for (const auto& it : split_fields(v, ',')) {
        try { out.push_back(int(parse_int(o, it))); }
        catch (const UsageError&) {
            throw UsageError(std::string(o.name) + " expects a comma-separated list of integers >= " + std::to_string(o.min) + ", got '" + v + "'");
        }
    }
    return out;
}

static std::pair<int, int> parse_grid(const Option& o, const std::string& v) {
    if (const auto x = v.find_first_of("xX"); x != std::string::npos) {
        try { return {int(parse_int(o, v.substr(0, x))), int(parse_int(o, v.substr(x + 1)))}; }
        catch (const UsageError&) {}
    }
    throw UsageError(std::string(o.name) + " expects RxC with whole numbers R, C >= 1 (e.g. 3x3), got '" + v + "'");
}

static void check_value(const Option& o, const std::string& v) {
    switch (o.kind) {
        case Kind::Flag: case Kind::Text: return;
        case Kind::Int: parse_int(o, v); return;
        case Kind::Seed: parse_seed(o, v); return;
        case Kind::Real: parse_real(o, v); return;
        case Kind::Choice: check_choice(o, v); return;
        case Kind::IntList: parse_int_list(o, v); return;
        case Kind::ChoiceList: parse_choice_list(o, v); return;
        case Kind::Grid: parse_grid(o, v); return;
    }
}

// The parsed command line. Every value was checked when it was parsed; looking up an option the command does not declare
// (or with the wrong type) is a programming error and throws std::logic_error.
class Args {
public:
    static Args parse(const Command& c, int argc, char** argv);   // argv[1] is the command; throws UsageError
    const std::string& pos(std::size_t i) const { return pos_.at(i); }
    bool has(const char* name) const { need(name, {}); return value(name) != nullptr; }   // given on the command line
    std::vector<std::string> given() const { std::vector<std::string> v; for (const auto& g : given_) v.push_back(g.first); return v; }
    bool flag(const char* name) const { need(name, {Kind::Flag}); return value(name) != nullptr; }
    std::string text(const char* name, const std::string& def = "") const {
        need(name, {Kind::Text, Kind::Choice});
        const auto* v = value(name);
        return v ? *v : def;
    }
    long long integer(const char* name, long long def) const {
        const Option& o = need(name, {Kind::Int});
        const auto* v = value(name);
        return v ? parse_int(o, *v) : def;
    }
    unsigned long long seed(const char* name, unsigned long long def) const {
        const Option& o = need(name, {Kind::Seed});
        const auto* v = value(name);
        return v ? parse_seed(o, *v) : def;
    }
    double real(const char* name, double def) const {
        const Option& o = need(name, {Kind::Real});
        const auto* v = value(name);
        return v ? parse_real(o, *v) : def;
    }
    std::vector<int> int_list(const char* name, std::vector<int> def) const {
        const Option& o = need(name, {Kind::IntList});
        const auto* v = value(name);
        return v ? parse_int_list(o, *v) : def;
    }
    std::vector<std::string> list(const char* name) const {   // empty when not given
        const Option& o = need(name, {Kind::ChoiceList});
        const auto* v = value(name);
        return v ? parse_choice_list(o, *v) : std::vector<std::string>{};
    }
    std::pair<int, int> grid(const char* name) const {   // only when given (has)
        const Option& o = need(name, {Kind::Grid});
        const auto* v = value(name);
        if (!v) throw std::logic_error(std::string("internal: qpuforge ") + cmd_->name + " reads " + name + " without checking it was given");
        return parse_grid(o, *v);
    }

private:
    explicit Args(const Command& c) : cmd_(&c) {}
    const Option& need(const char* name, std::initializer_list<Kind> kinds) const {
        const Option* o = cmd_->find(name);
        if (!o || (kinds.size() && std::find(kinds.begin(), kinds.end(), o->kind) == kinds.end()))
            throw std::logic_error(std::string("internal: qpuforge ") + cmd_->name + " reads option " + name + ", which its option table does not declare with this type");
        return *o;
    }
    const std::string* value(const std::string& name) const {
        for (const auto& g : given_) if (g.first == name) return &g.second;
        return nullptr;
    }
    const Command* cmd_;
    std::vector<std::pair<std::string, std::string>> given_;   // (option, value) in command-line order
    std::vector<std::string> pos_;
};

static std::string unknown_option(const Command& c, const std::string& t) {
    const std::string m = "unknown option '" + t + "'";
    const auto eq = t.find('=');
    const std::string base = t.substr(0, eq);
    if (const Option* o = eq == std::string::npos ? nullptr : c.find(base))
        return m + (o->kind == Kind::Flag ? " (" + base + " takes no value)" : " (did you mean '" + base + " " + t.substr(eq + 1) + "'?)");
    if (base == "--version" || base == "-V") return m + " (--version / -V must be the first argument)";
    std::vector<std::string> names;
    for (const auto& o : c.options) names.push_back(o.name);
    if (const auto s = suggest(base, names); !s.empty()) return m + " (did you mean '" + s + "'?)";
    std::vector<std::string> others;
    for (const auto& k : commands()) if (&k != &c && k.find(base)) others.push_back(k.name);
    return others.empty() ? m : m + " (an option of " + join_words(others, " and ") + ")";
}

Args Args::parse(const Command& c, int argc, char** argv) {
    Args a(c);
    for (int i = 2; i < argc; ++i) {
        const std::string t = argv[i];
        if (!is_option(t)) { a.pos_.push_back(t); continue; }
        const Option* o = c.find(t);
        if (!o) throw UsageError(unknown_option(c, t));
        if (a.value(t)) throw UsageError("option '" + t + "' given more than once");
        std::string v;
        if (o->kind != Kind::Flag) {
            // the next argument is the value unless it is another long option ("--shots --traj 5"); negative numbers
            // ("--alpha -0.22") are values
            if (i + 1 >= argc || std::strncmp(argv[i + 1], "--", 2) == 0) throw UsageError("option '" + t + "' requires a value");
            v = argv[++i];
            check_value(*o, v);
        }
        a.given_.emplace_back(t, v);
    }
    if (a.pos_.size() < c.positionals.size()) throw UsageError(std::string("missing ") + c.positionals[a.pos_.size()]);
    if (a.pos_.size() > c.positionals.size()) throw UsageError("unexpected argument '" + a.pos_[c.positionals.size()] + "'");
    return a;
}

// Write a JSON file; throws if the file cannot be opened or written (as gds::write does).
static void write_json_file(const std::string& path, const Json& j) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    j.dump(f);
    f << "\n";
    f.close();
    if (!f) throw std::runtime_error("cannot write " + path);
}

// ================================================================== commands
// <circuit> = file.qasm | ghz:N | qft:N | bv:N | xeb:CYCLES. Checked (UsageError) before the chip file is read.
struct CircuitSpec { std::string kind, file; int n = 0; };

static CircuitSpec parse_circuit(const std::string& spec) {
    CircuitSpec s;
    const auto c = spec.find(':');
    if (c == std::string::npos || spec.find(".qasm") != std::string::npos) { s.file = spec; return s; }
    s.kind = spec.substr(0, c);
    if (s.kind != "ghz" && s.kind != "qft" && s.kind != "bv" && s.kind != "xeb")
        throw UsageError("unknown circuit '" + spec + "' (expected file.qasm, ghz:N, qft:N, bv:N or xeb:CYCLES)");
    try { s.n = int(parse_int(integer("N", 1), spec.substr(c + 1))); }
    catch (const UsageError&) { throw UsageError("circuit '" + spec + "': " + s.kind + ":N needs a whole number N >= 1"); }
    return s;
}

static Circuit make_circuit(const CircuitSpec& spec, const ChipDesign& D, const ChipReport& R, bool& native_layout) {
    native_layout = false;
    if (spec.kind.empty()) return load_qasm(spec.file);
    const int n = spec.n;
    if (spec.kind == "ghz") return ghz(n);
    if (spec.kind == "qft") return qft(n);
    if (spec.kind == "bv") {   // fixed pseudo-random secret, never all-zero
        unsigned long long secret = 0x2D5B7ULL & (n >= 64 ? ~0ULL : (1ULL << n) - 1);
        if (!secret) secret = 1;
        return bernstein_vazirani(n, secret);
    }
    std::vector<std::pair<int, int>> edges;   // xeb
    for (const auto& p : R.pairs) if (p.designed) edges.push_back({D.index_of(p.a), D.index_of(p.b)});
    native_layout = true;
    return random_xeb(int(D.qubits.size()), edges, n, 7);
}

static void print_report(const ChipDesign& D, const ChipReport& R) {
    std::printf("\n%s — %zu qubits, %zu couplers, architecture %s\n", D.name.c_str(), D.qubits.size(), D.couplings.size(), D.architecture.c_str());
    std::printf("materials: Q_diel=%.2g  junction gap %.0f ueV (Rn spread %.1f%%)  T_eff=%.0f mK  flux noise %.1f uPhi0\n\n",
                D.Q_dielectric, D.gap_ueV, D.Rn_spread_pct, D.T_eff_mK, D.flux_noise_uPhi0);
    std::printf(" Q  role       C_fF   Rn_kOhm  EJ/EC   f01_GHz  dressed   alpha_MHz  disp_kHz  T1_us(diel/Purcell)      Tphi_us   T2_us   P_th    ro(1|0) ro(0|1)  e_1q\n");
    for (std::size_t i = 0; i < R.q.size(); ++i) {
        const auto& q = R.q[i]; const auto& d = D.qubits[i];
        std::printf("%2d  %-9s %6.1f  %7.3f  %6.1f   %7.4f  %7.4f   %8.1f  %8.3g  %6.1f (%6.0f/%7.0f)  %8.0f  %6.1f  %.4f  %.4f  %.4f  %.1e\n",
                    q.id, d.role.c_str(), d.C_fF, d.Rn_kOhm, q.EJ_EC, q.f01, q.f01_dressed, q.alpha * 1e3, q.charge_disp_kHz, q.T1_us, q.T1_diel_us,
                    q.T1_purcell_us, q.Tphi_us, q.T2_us, q.p_thermal, q.ro_err_0to1, q.ro_err_1to0, q.err_1q);
    }
    std::printf("\n pair     coupler   g_MHz   detune_MHz   ZZ_static_kHz  (pert.)   ZZ_idle_kHz   t_CZ_ns   CZ point GHz   e_2q     flux_e  c_GHz  flags\n");
    for (const auto& p : R.pairs) {
        std::printf(" %2d-%-2d    %-7s  %6.3f   %9.1f   %12.2f  (%7.2f)  %11.2f   %7.1f   %11.4f   %.1e  %.1e  %5.2f  %s\n", p.a, p.b, p.designed ? "yes" : "parasit", p.g_MHz,
                    p.delta_MHz, p.zz_kHz, p.zz_pert_kHz, p.zz_idle_kHz, p.t_cz_ns, p.cz_interaction_GHz, p.err_2q, p.cz_flux_err, p.coupler_f_GHz, p.flags.c_str());
    }
    std::printf("\nmean T1 %.1f us, mean T2 %.1f us, mean e_1q %.2e, mean e_2q %.2e, max |idle ZZ| %.1f kHz, collisions %d\n",
                R.mean_T1_us, R.mean_T2_us, R.mean_err_1q, R.mean_err_2q, R.max_abs_zz_kHz, R.n_collisions);
    for (const auto& w : R.warnings) std::printf("  ! %s\n", w.c_str());
}

static int cmd_new(const Args& a) {
    const std::string& path = a.pos(0);
    if (a.has("--preset")) {
        // a named design is fixed apart from its distance; the generator options would be ignored
        for (const auto& o : a.given())
            if (o != "--preset" && o != "--distance") throw UsageError(o + " cannot be combined with --preset (a preset design takes only --distance)");
        const std::string pre = a.text("--preset");
        const int d = int(a.integer("--distance", pre == "beta1" ? 5 : pre == "gamma1" ? 7 : 3));
        ChipDesign D = pre == "alpha1" ? alpha1(d) : pre == "alpha1v11" ? alpha1_v11(d) : pre == "beta1" ? beta1(d) : pre == "gamma1" ? gamma1(d)
                     : gen0_surface17(d);
        save_design(D, path);
        std::printf("wrote %s (%zu qubits, %zu couplers)\n", path.c_str(), D.qubits.size(), D.couplings.size());
        print_report(D, analyze(D));
        return 0;
    }
    if (a.has("--distance")) throw UsageError("--distance needs --preset (a generated surface-code chip takes --surface D)");
    {
        std::string topology;
        for (const char* t : {"--grid", "--line", "--ring", "--surface"}) {
            if (!a.has(t)) continue;
            if (!topology.empty()) throw UsageError(topology + " and " + t + " cannot be combined (one topology per chip)");
            topology = t;
        }
    }
    GenSpec s;
    if (a.has("--surface")) {
        s = surface_spec(int(a.integer("--surface", 3)));
        s.f_hi = 5.9; s.f_lo = 5.2;
    }
    s.feedlines = int(a.integer("--feedlines", s.feedlines));
    s.kappa_MHz = a.real("--kappa-MHz", s.kappa_MHz);
    s.chi_MHz = a.real("--chi-MHz", s.chi_MHz);
    s.purcell_filter = a.real("--purcell", s.purcell_filter);
    if (a.has("--grid")) { s.topology = "grid"; std::tie(s.rows, s.cols) = a.grid("--grid"); }
    if (a.has("--line")) { s.topology = "line"; s.rows = 1; s.cols = int(a.integer("--line", 1)); }
    if (a.has("--ring")) { s.topology = "ring"; s.rows = 1; s.cols = int(a.integer("--ring", 2)); }
    s.f_hi = a.real("--f-hi", s.topology == "surface" ? s.f_hi : 5.4);
    s.f_lo = a.real("--f-lo", s.topology == "surface" ? s.f_lo : 4.8);
    s.alpha = a.real("--alpha", -0.22);
    s.g = a.real("--g-MHz", 4) / 1000.0;
    s.plan = std::stoi(a.text("--plan", "8"));   // 2, 4 or 8 (checked)
    s.stagger = a.real("--stagger-MHz", 60) / 1000.0;
    s.Q = a.real("--Q", 3e6);
    s.architecture = a.text("--arch", "fixed") == "tunable" ? "tunable-coupler-cz" : "fixed-coupling-cz";
    auto D = generate(s);
    if (const auto nm = a.text("--name"); !nm.empty()) D.name = nm;
    D.t_readout_ns = a.real("--readout-ns", D.t_readout_ns);
    D.t_reset_ns = a.real("--reset-ns", D.t_reset_ns);
    D.reset_err = a.real("--reset-err", D.reset_err);
    D.err_2q_control = a.real("--err2q", D.err_2q_control);
    D.tc_residual_zz_kHz = a.real("--residual-zz-kHz", D.tc_residual_zz_kHz);
    if (a.flag("--active-reset")) D.active_reset = true;
    save_design(D, path);
    std::printf("wrote %s (%zu qubits, %zu couplers)\n", path.c_str(), D.qubits.size(), D.couplings.size());
    print_report(D, analyze(D));
    return 0;
}

// --materials: one line per sidecar quantity, applied or skipped, before the device report
static void print_materials(const std::string& path, const MaterialsSidecar& S, const MaterialsReport& M) {
    const auto& p = S.produced_by;
    std::string meta;
    for (const std::string& s : {p.commit.empty() ? std::string() : "commit " + p.commit, p.potential_or_method, p.timestamp})
        if (!s.empty()) meta += (meta.empty() ? "" : ", ") + s;
    std::string sys;
    for (const auto& el : S.system) sys += (sys.empty() ? "" : "-") + el;
    std::printf("\nmaterials sidecar %s: %s", path.c_str(), S.schema.c_str());
    if (!p.tool.empty()) std::printf(" from %s", p.tool.c_str());
    if (!meta.empty()) std::printf(" (%s)", meta.c_str());
    if (!sys.empty()) std::printf(", system %s", sys.c_str());
    std::printf("\n");
    for (const auto& l : M.lines) {
        if (l.applied) std::printf("  applied  %-17s %s %.6g -> %.6g, %s\n", l.quantity.c_str(), l.field.c_str(), l.before, l.after, l.note.c_str());
        else std::printf("  skipped  %-17s %s\n", l.quantity.c_str(), l.note.c_str());
    }
}

static int cmd_analyze(const Args& a) {
    const std::string mat = a.text("--materials");
    if (a.has("--materials") && mat.empty()) throw UsageError("--materials needs a sidecar file");
    if (a.flag("--loss-tangent") && mat.empty()) throw UsageError("--loss-tangent needs --materials <sidecar.json>");
    auto D = load_design(a.pos(0));
    if (!mat.empty()) {
        MaterialsApplyOptions o;
        o.use_loss_tangent = a.flag("--loss-tangent");
        const auto S = load_materials(mat);
        const auto M = apply_materials(S, D, o);
        print_materials(mat, S, M);
    }
    print_report(D, analyze(D));
    return 0;
}

static RunOptions run_opts(const Args& a) {   // --shots is read by run only: budget computes state fidelities
    RunOptions o;
    o.trajectories = int(a.integer("--traj", 300));
    o.seed = unsigned(a.seed("--seed", 1234));
    return o;
}

static int cmd_run(const Args& a) {
    const auto spec = parse_circuit(a.pos(1));
    const auto D = load_design(a.pos(0));
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    bool native = false;
    const auto C = make_circuit(spec, D, R, native);
    CompileOptions co; if (native) co.layout = "trivial";
    const auto P = compile(C, G, co);
    auto nm = noise_model(D, R);
    for (const auto& s : a.list("--off")) nm.disable(s);
    auto o = run_opts(a);
    o.shots_per_traj = int(a.integer("--shots", 20));
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = run(P, nm, o);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("\n%s on %s\n", C.name.c_str(), D.name.c_str());
    std::printf("  compiled: %d moments, %.0f ns, %d CZ, %d 1q pulses, %d virtual RZ, %d SWAPs, %d simulated qubits\n",
                r.depth, r.duration_ns, r.n_cz, r.n_pulses, P.n_virtual_rz, r.n_swaps, r.n_sim);
    std::printf("  layout (logical->physical):");
    for (std::size_t l = 0; l < P.layout_initial.size(); ++l) std::printf(" %zu->Q%d", l, D.qubits[P.layout_initial[l]].id);
    std::printf("\n  state fidelity  %.4f +- %.4f   (%d trajectories)\n", r.fidelity, r.fidelity_sem, o.trajectories);
    if (r.shots) {
        if (std::abs(r.xeb_ideal) > 0.05) {
            std::printf("  P(ideal most-likely outcome %s) = %.4f over %ld shots\n", r.ideal_top.c_str(), r.p_ideal_top, r.shots);
            std::printf("  linear XEB %.4f (normalised %.4f; ideal-circuit value %.3f)\n", r.xeb, r.xeb_norm, r.xeb_ideal);
        } else {
            std::printf("  %ld shots; ideal output distribution is uniform, so XEB / most-likely-outcome are not informative\n", r.shots);
        }
        std::vector<std::pair<long, std::string>> top;
        for (auto& [k, c] : r.counts) top.push_back({c, k});
        std::sort(top.rbegin(), top.rend());
        const int K = int(a.integer("--top", 8));
        std::printf("  top outcomes:\n");
        for (int i = 0; i < K && i < int(top.size()); ++i) std::printf("    %s  %6.4f\n", top[i].second.c_str(), double(top[i].first) / r.shots);
    }
    std::printf("  (%.2f s)\n", secs);
    return 0;
}

static int cmd_budget(const Args& a) {
    const auto spec = parse_circuit(a.pos(1));
    const auto D = load_design(a.pos(0));
    const auto R = analyze(D);
    const auto G = chip_graph(D, R);
    bool native = false;
    const auto C = make_circuit(spec, D, R, native);
    CompileOptions co; if (native) co.layout = "trivial";
    const auto P = compile(C, G, co);
    auto o = run_opts(a); o.shots_per_traj = 0;
    const auto base = noise_model(D, R);
    const auto all = run(P, base, o);
    std::printf("\nerror budget for %s on %s (%d CZ, %.0f ns)\n", C.name.c_str(), D.name.c_str(), P.n_cz, P.duration_ns);
    std::printf("  all sources: infidelity %.4f +- %.4f\n", 1 - all.fidelity, all.fidelity_sem);
    double sum = 0;
    std::vector<std::pair<double, std::string>> rows;
    for (const std::string src : {"T1", "Tphi", "thermal", "zz", "exchange", "control"}) {
        auto nm = base; nm.only(src);
        const auto r = run(P, nm, o);
        rows.push_back({1 - r.fidelity, src});
        sum += 1 - r.fidelity;
    }
    std::sort(rows.rbegin(), rows.rend());
    for (auto& [inf, src] : rows) std::printf("  %-8s alone: infidelity %.4f  (%4.1f%% of sum)\n", src.c_str(), inf, 100 * inf / std::max(sum, 1e-12));
    std::printf("  (readout assignment error acts on the measured counts, not on the state fidelity)\n");
    return 0;
}

static int cmd_fab(const Args& a) {
    const auto D = load_design(a.pos(0));
    const int N = int(a.integer("--samples", 200));
    const auto F = fab_monte_carlo(D, N, unsigned(a.seed("--seed", 1)));
    auto stats = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return std::array<double, 3>{v[v.size() / 10], v[v.size() / 2], v[v.size() * 9 / 10]}; };
    const auto fs = stats(F.f01_sigma_MHz), es = stats(F.err2q_mean), zs = stats(F.max_zz_kHz);
    std::printf("\nfabrication Monte Carlo: %d chips, junction Rn spread %.1f%%\n", N, D.Rn_spread_pct);
    std::printf("  collision-free yield: %d / %d = %.1f%%\n", F.collision_free, N, 100.0 * F.collision_free / N);
    std::printf("  RMS frequency error     p10/p50/p90: %.1f / %.1f / %.1f MHz\n", fs[0], fs[1], fs[2]);
    std::printf("  mean 2q error           p10/p50/p90: %.2e / %.2e / %.2e\n", es[0], es[1], es[2]);
    std::printf("  worst idle ZZ           p10/p50/p90: %.0f / %.0f / %.0f kHz\n", zs[0], zs[1], zs[2]);
    for (auto& [m, c] : F.failure_modes) std::printf("  failure mode %-18s %d pair-occurrences\n", m.c_str(), c);
    return 0;
}

static int cmd_layout(const Args& a) {
    const auto D = load_design(a.pos(0));
    write_layout_svg(D, analyze(D), a.pos(1));
    std::printf("wrote %s\n", a.pos(1).c_str());
    return 0;
}

static int cmd_cz(const Args& a) {
    auto D = load_design(a.pos(0));
    D.coupler.physical = true;
    D.coupler.g_qc_MHz = a.real("--g-qc", D.coupler.g_qc_MHz);
    D.coupler.f_on_GHz = a.real("--f-on", D.coupler.f_on_GHz);
    D.coupler.t_rise_ns = a.real("--rise", D.coupler.t_rise_ns);
    auto R = analyze(D);
    const auto cal = calibrate_chip_cz(D, R, true, a.text("--cache"));
    std::printf("\n%s — pulse-level CZ calibration (coupler: f_max %.2f GHz, alpha %.0f MHz, g_qc %.0f MHz, on ~%.2f GHz, %.0f ns ramps)\n",
                D.name.c_str(), D.coupler.f_max_GHz, D.coupler.alpha_GHz * 1e3, D.coupler.g_qc_MHz, D.coupler.f_on_GHz, D.coupler.t_rise_ns);
    std::printf("  pair   mover  f_m/f_t GHz   g12 MHz  ZZ@park kHz  f_int  fc_on  t_cz ns  1-F coh  leak     flux     ctrl     total\n");
    for (const auto& c : cal)
        std::printf("  %2d-%-2d   q%-3d %.3f/%.3f   %6.2f    %+6.2f     %.3f  %.3f  %6.1f   %.1e  %.1e  %.1e  %.1e  %.1e\n", c.a, c.b, c.mover,
                    c.f_mover, c.f_target, c.g12_MHz, c.zz_idle_kHz, c.f_int_GHz, c.fc_on_GHz, c.t_cz_ns, c.err_coherent, c.leakage + c.leak_spectators,
                    c.err_flux_local + c.err_flux_quasistatic, c.err_ctrl, c.err_total);
    double mt = 0, me = 0, ml = 0; for (const auto& c : cal) { mt += c.t_cz_ns; me += c.err_total; ml += c.leakage + c.leak_spectators; }
    const double n = std::max<std::size_t>(cal.size(), 1);
    std::printf("  mean: t_cz %.1f ns, CZ error %.2e (incl. decoherence), leakage %.1e per CZ\n", mt / n, me / n, ml / n);
    if (const auto out = a.text("--out"); !out.empty()) { write_json_file(out, cz_calibration_json(cal)); std::printf("  wrote %s\n", out.c_str()); }
    return 0;
}

static int cmd_gds(const Args& a) {
    using namespace layout;
    auto D = load_design(a.pos(0));
    auto R = analyze(D);
    const std::string cells_path = a.text("--cells");
    CellSet C; LayoutMetrics m; SynthesisLog log;
    const auto T = targets_from_design(D, R);
    bool have = false;
    if (!cells_path.empty()) {
        std::ifstream f(cells_path);
        if (f && !a.flag("--synthesize")) {
            std::stringstream ss; ss << f.rdbuf();
            const Json j = Json::parse(ss.str());
            C = cells_from_json(j); have = true;
            std::printf("loaded cell geometry from %s\n", cells_path.c_str());
        }
    }
    if (!have) {
        std::printf("synthesising cell geometry against the circuit targets (method of moments; ~5-10 min)...\n");
        C.pitch = 900;
        if (const auto from = a.text("--from"); !from.empty()) {      // warm start from an earlier cell set
            std::ifstream f(from);
            if (!f) throw std::runtime_error("cannot open " + from);
            std::stringstream ss; ss << f.rdbuf();
            C = cells_from_json(Json::parse(ss.str()));
            std::printf("  starting from %s\n", from.c_str());
        }
        C = synthesize_cells(C, T, &log);
        for (const auto& s : log.lines) std::printf("  %s\n", s.c_str());
    }
    m = evaluate_cells(C, MeshOptions{1.5, 50.0, 1.6});
    // A loaded cell file is read-only: it records the targets of the chip it was synthesised for, and this chip's
    // targets may differ (the same cells are reused across chips).
    if (!cells_path.empty() && !have) {
        write_json_file(cells_path, cells_json(C, m, T, &log));
        std::printf("wrote cell geometry to %s\n", cells_path.c_str());
    }
    const auto L = assemble_chip(D, C);
    write_gds(L, a.pos(1), a.text("--cell", "TOP"));
    std::printf("\n%s — layout: %.1f x %.1f mm, %zu etch polygons, %zu junction regions, %zu bump pads -> %s\n", D.name.c_str(),
                L.outline.w() / 1e3, L.outline.h() / 1e3, L.etch.size(), L.jj.size(), L.bumps.size(), a.pos(1).c_str());
    std::printf("  DRC: %s\n", drc_summary(L).c_str());
    for (const auto& s : L.problems) std::printf("    %s\n", s.c_str());
    for (std::size_t i = 0; i < L.drc.size() && i < 10; ++i)
        std::printf("    %s %.3f um at (%.3f, %.3f) um\n", L.drc[i].rule.c_str(), L.drc[i].value, L.drc[i].box.cx(), L.drc[i].box.cy());
    std::printf("  extracted circuit (fine mesh) vs targets:\n");
    auto row = [](const char* n, double a, double t) { std::printf("    %-28s %10.5f  target %10.5f  (%+.2f%%)\n", n, a, t, 100 * (a / t - 1)); };
    row("E_C measure (GHz)", m.EC_m, T.EC_m); row("E_C data (GHz)", m.EC_d, T.EC_d); row("E_C coupler (GHz)", m.EC_c, T.EC_c);
    row("G measure-coupler (GHz/nn)", m.G_mc, T.G_mc); row("G coupler-data (GHz/nn)", m.G_cd, T.G_cd); row("G measure-data (GHz/nn)", m.G_md, T.G_md);
    row("C qubit-readout (fF)", m.C_qr_fF, T.C_qr_fF);
    const auto P = post_layout(D, R, m);
    std::printf("  post-layout circuit: alpha m/d/c %.1f/%.1f/%.1f MHz, g_qc %.1f MHz, g12 %.2f MHz (exchange-free park needs %.2f),\n"
                "    residual parked |11>-|20> gap %.3f MHz -> spectator transfer %.1e per crossing, ZZ at park %+.2f kHz, readout g %.1f MHz (design %.1f)\n",
                P.alpha_m * 1e3, P.alpha_d * 1e3, P.alpha_c * 1e3, P.g_qc_MHz, P.g12_MHz, P.g12_null_MHz, P.residual_gap_MHz, P.spectator_LZ, P.zz_park_kHz,
                P.g_r_MHz, P.g_r_target_MHz);
    std::printf("  re-null by parking the couplers at %.3f GHz (%.0f MHz below the sweet spot, df/dPhi %.2f GHz/Phi0): residual gap %.3f MHz, ZZ %+.2f kHz\n",
                P.park_fc_GHz, (D.coupler.f_max_GHz - P.park_fc_GHz) * 1e3, P.park_dfdphi_GHz, P.park_gap_MHz, P.park_zz_kHz);
    for (const auto& r : L.resonators) if (std::abs(r.drawn_um - r.target_um) > 1) std::printf("    resonator R%d length %.0f um vs target %.0f um\n", r.qubit, r.drawn_um, r.target_um);
    return L.drc_clean() ? 0 : 2;
}

static int cmd_qec(const Args& a) {
    if (!a.flag("--bursts"))
        for (const char* o : {"--burst-rate", "--burst-tau-ms", "--burst-t1-us", "--burst-tphi-us", "--veto"})
            if (a.has(o)) throw UsageError(std::string(o) + " needs --bursts");
    if (!a.flag("--veto"))
        for (const char* o : {"--veto-window", "--veto-radius", "--veto-false", "--impact-p", "--burst-radius"})
            if (a.has(o)) throw UsageError(std::string(o) + " needs --veto");
    if (!a.flag("--herald"))
        for (const char* o : {"--herald-fidelity", "--herald-false", "--herald-factor", "--herald-rule"})
            if (a.has(o)) throw UsageError(std::string(o) + " needs --herald");
    if (a.text("--decoder", "exact") != "uf-stream")
        for (const char* o : {"--window", "--commit"})
            if (a.has(o)) throw UsageError(std::string(o) + " needs --decoder uf-stream");
    if (a.integer("--commit", 3) >= a.integer("--window", 6))
        throw UsageError("--commit (" + std::to_string(a.integer("--commit", 3)) + ") must be smaller than --window (" +
                         std::to_string(a.integer("--window", 6)) + ")");
    auto D = load_design(a.pos(0));
    if (const auto v = a.text("--leakage-removal"); !v.empty()) {
        D.leakage.reset_removes_leakage = v == "lrr" || v == "lrr+dqlr";
        D.leakage.dqlr = v == "dqlr" || v == "lrr+dqlr";
    }
    auto R = analyze(D);
    if (D.coupler.physical) calibrate_chip_cz(D, R, true, a.text("--cache"));
    auto nm = noise_model(D, R);
    for (const auto& s : a.list("--off")) nm.disable(s);
    const std::vector<int> rounds = a.int_list("--rounds", {2, 4, 6, 8, 10});
    const long shots = long(a.integer("--shots", 200000));
    const int sv_shots = int(a.integer("--statevector", 0));
    const unsigned seed = unsigned(a.seed("--seed", 2026));
    const std::string bsel = a.text("--basis", "Z");   // Z, X, both (or z, x)
    const std::vector<char> bases = bsel == "both" ? std::vector<char>{'Z', 'X'} : std::vector<char>{bsel == "X" || bsel == "x" ? 'X' : 'Z'};
    FrameRunOptions fro;
    {
        const std::string dsel = a.text("--decoder", "exact");   // exact, uf (union-find), uf-stream, legacy or greedy (checked)
        if (dsel == "legacy") fro.method = DecodeMethod::Legacy;
        else if (dsel == "greedy") fro.method = DecodeMethod::Greedy;
        else if (dsel == "uf" || dsel == "union-find") fro.method = DecodeMethod::UnionFind;
        else if (dsel == "uf-stream") fro.method = DecodeMethod::UnionFind;   // the streaming decoder is built per setup below
        if (dsel != "exact")
            std::printf("decoder: %s\n", dsel == "legacy" ? "pre-v0.4 behaviour (exact up to 22 defects, greedy beyond)"
                                         : dsel == "greedy" ? "greedy matching"
                                         : dsel == "uf-stream" ? "streaming union-find (sliding window over rounds; not exact)"
                                                               : "weighted union-find (real-time candidate; not exact)");
        if (dsel == "uf-stream")
            std::printf("  window %lld rounds, committing %lld per step\n", a.integer("--window", 6), a.integer("--commit", 3));
    }
    const bool herald = a.flag("--herald");
    if (herald) {
        fro.herald.on = true;
        fro.herald.fidelity = a.real("--herald-fidelity", 0.95);
        fro.herald.false_rate = a.real("--herald-false", 0.001);
        fro.herald.footprint_factor = a.real("--herald-factor", 0.4);
        const std::string rule = a.text("--herald-rule", "footprint");   // footprint or detector (checked)
        fro.herald.footprint = rule == "footprint";
        if (!fro.herald.footprint) fro.herald.touch_factor = fro.herald.footprint_factor;
        std::printf("heralded leakage: 3-state readout of the measure qubits, fidelity %.3f, false-herald rate %.1e; %s rule, factor %.2f%s\n",
                    fro.herald.fidelity, fro.herald.false_rate, rule.c_str(), fro.herald.footprint_factor,
                    nm.leakage ? (nm.reset_removes_leakage ? "" : " (the reset does not remove leakage)") : " (leakage is not simulated: only false heralds)");
    }
    double best_T1 = 0, best_T2 = 0;
    for (const auto& q : R.q) { best_T1 = std::max(best_T1, q.T1_us); best_T2 = std::max(best_T2, q.T2_us); }
    for (char basis : bases) {
        std::printf("\n%s — memory-%c experiment, distance %d, Pauli-frame Monte Carlo (%ld shots per point)\n", D.name.c_str(), basis, D.surface_d, shots);
        std::printf("  rounds   round_ns   detection   P_L                eps/round   defects/shot (max)   decode us/shot\n");
        std::vector<MemoryResult> all;
        for (int r : rounds) {
            const auto M = prepare_memory(D, R, nm, r, basis);
            FrameRunOptions fr = fro;
            std::unique_ptr<StreamingUnionFind> stream;
            if (a.text("--decoder", "exact") == "uf-stream") {
                stream = std::make_unique<StreamingUnionFind>(int(M.exp.detectors.size()), M.dem.edges, M.exp.detector_round,
                                                              int(a.integer("--window", 6)), int(a.integer("--commit", 3)));
                fr.stream = stream.get();
            }
            const auto res = run_memory_frame(M, shots, seed + unsigned(r) + (basis == 'X' ? 7919u : 0u), fr);
            all.push_back(res);
            std::printf("  %5d   %8.0f   %8.4f    %.5f +- %.5f   %.5f   %8.1f (%d)   %10.1f", r, res.round_ns, res.detection_rate, res.p_logical, res.p_logical_sem,
                        res.eps_per_round, res.mean_defects, res.max_defects, res.decode_us_per_shot);
            if (M.frame.leakage) std::printf("   leaked: data %.1e meas %.1e", res.leak_frac_data, res.leak_frac_meas);
            std::printf("\n");
            if (herald) {
                FrameRunOptions blind = fro; blind.herald.use_in_decoder = false;
                const auto nb = run_memory_frame(M, shots, seed + unsigned(r) + (basis == 'X' ? 7919u : 0u), blind);
                std::printf("          heralds %.2f per shot (%.0f%% leaked); same syndromes decoded without heralds: P_L %.5f, eps %.5f (with: %.5f)\n",
                            res.heralds_per_shot, 100 * res.herald_true_frac, nb.p_logical, nb.eps_per_round, res.eps_per_round);
            }
            if (res.ignored_exchange_pairs) std::printf("    ! %d near-degenerate exchange pair(s) not modelled by the Pauli-frame simulator\n", res.ignored_exchange_pairs);
            if (sv_shots > 0 && r == rounds.front()) {
                const auto sv = run_memory_statevector(M, nm, sv_shots, seed);
                std::printf("          exact statevector (%d shots): detection %.4f, P_L %.4f +- %.4f\n", sv_shots, sv.detection_rate, sv.p_logical, sv.p_logical_sem);
            }
        }
        const double eps = fit_eps_per_round(all);
        const double round_us = all.back().round_ns / 1e3;
        std::printf("  logical error per round (fit): %.5f +- %.5f   ->  logical lifetime ~ %.1f us (round %.2f us)\n", eps, fit_eps_per_round_err(all),
                    round_us / std::max(eps, 1e-12), round_us);
    }
    std::printf("  best physical qubit on chip: T1 %.1f us, T2 %.1f us\n", best_T1, best_T2);
    if (a.flag("--bursts")) {
        BurstSpec B;
        B.rate_hz = a.real("--burst-rate", 0.1);
        B.tau_ms = a.real("--burst-tau-ms", 25);
        B.t1_peak_us = a.real("--burst-t1-us", 1);
        B.tphi_peak_us = a.real("--burst-tphi-us", 0);
        const int br = rounds.back();
        for (char basis : bases) {
            const auto b = burst_floor(D, R, nm, B, br, shots, seed + 5, basis);
            std::printf("\n  high-energy impacts, memory-%c (%.2g per s, recovery %.3g ms, peak excess T1 %.3g us, T_phi %.3g us; %d rounds per shot)\n",
                        basis, B.rate_hz, B.tau_ms, B.t1_peak_us, B.tphi_peak_us, br);
            std::printf("    severity   eps/round   defects/shot\n");
            std::printf("    %8s   %.5f     %.1f   (quiet chip)\n", "0", b.eps0, b.mean_defects0);
            for (const auto& pt : b.curve) std::printf("    %8.1e   %.5f     %.1f%s\n", pt.s, pt.eps, pt.mean_defects, pt.saturated ? "   (randomised within the run)" : "");
            std::printf("    one impact flips the logical qubit with probability %.2f; unprotected for %.0f ms, visible in the syndrome for %.0f ms\n",
                        b.p_fail, b.dead_time_ms, b.detectable_ms);
            std::printf("    logical-error floor from impacts: %.2e per round (quiet chip: %.2e)\n", b.floor_per_round, b.eps0);
        }
        if (a.flag("--veto")) {
            VetoSpec V;
            V.window_rounds = int(a.integer("--veto-window", 3));
            V.radius = a.real("--veto-radius", 2);
            const double target = a.real("--veto-false", 0.001);
            ImpactSpec I; I.burst = B; I.burst.radius = a.real("--burst-radius", 0);
            I.p_shot = a.real("--impact-p", 0.05);
            for (char basis : bases) {
                const auto M = prepare_memory(D, R, nm, br, basis);
                const auto cal = calibrate_veto(M, shots, seed + 91, target, V, fro);
                V.on = true; V.threshold = cal.threshold;
                FrameRunOptions inj = fro; inj.impact = I;
                FrameRunOptions vet = inj; vet.veto = V;
                FrameRunOptions clean = fro; clean.veto = V;
                const auto rc = run_memory_frame(M, shots, seed + 93, clean);
                const auto ri = run_memory_frame(M, shots, seed + 93, inj);
                const auto rv = run_memory_frame(M, shots, seed + 93, vet);
                std::printf("\n  impact veto, memory-%c (%d rounds): window %d rounds x radius %.1f, threshold %d events (false-veto rate %.1e on %ld burst-free shots, target %.1e)\n",
                            basis, br, V.window_rounds, V.radius, V.threshold, cal.false_rate, cal.shots, target);
                std::printf("    impacts injected in %.1f%% of shots (onset uniform in the shot, %s): %ld impacts, %ld caught (%.0f%%)\n",
                            100 * I.p_shot, I.burst.radius > 0 ? "local" : "chip-wide", rv.impacts, rv.impacts_vetoed,
                            100.0 * double(rv.impacts_vetoed) / double(std::max(rv.impacts, 1L)));
                std::printf("    eps/round: burst-free %.5f; with impacts %.5f without the veto, %.5f with it (%.2f%% of shots discarded; burst-free: %.2f%%)\n",
                            rc.eps_per_round, ri.eps_per_round, rv.eps_kept, 100.0 * double(rv.discarded) / double(rv.shots), 100.0 * double(rc.discarded) / double(rc.shots));
                std::printf("    failures in impacted shots: %ld without the veto, %ld left after it\n", ri.errors_impacted, rv.errors_impacted_kept);
                const auto bv = burst_floor_veto(D, R, nm, B, br, shots, seed + 5, V, basis, fro.method);
                std::printf("    over the burst profile (chip-wide severity ladder): floor %.2e per round without the veto, %.2e with shot-by-shot discard (%.2f ms discarded per impact)\n",
                            bv.floor_per_round, bv.floor_per_round_veto, bv.discarded_ms);
                std::printf("    latched veto (discard for a hold time after the onset):");
                for (const auto& L : bv.latched) std::printf("  %.0f tau: %.1e (dead %.1e)", L.hold_tau, L.floor_per_round, L.dead_time_frac);
                std::printf("\n");
                V.on = false;
            }
        }
    }
    return 0;
}

static ChipReport calibrated_report(const ChipDesign& D, const Args& a) {
    auto R = analyze(D);
    if (D.coupler.physical) calibrate_chip_cz(D, R, false, a.text("--cache"));
    return R;
}

static int cmd_synth(const Args& a) {
    const auto D = load_design(a.pos(0));
    const auto R = calibrated_report(D, a);
    SynthOptions o;
    o.seed = unsigned(a.seed("--seed", 1));
    o.shots = long(a.integer("--shots", 20000));
    o.hidden_dephasing_us = a.real("--hidden-dephasing-us", 0);
    o.rounds = a.int_list("--rounds", o.rounds);
    if (a.flag("--z-only")) o.both_bases = false;
    const auto C = synthesize_characterisation(D, R, o);
    write_json_file(a.pos(1), characterisation_to_json(C));
    std::printf("wrote %s: %zu qubits, %zu pairs, %zu memory runs (synthetic data, seed %u)\n", a.pos(1).c_str(), C.qubits.size(), C.pairs.size(), C.memory.size(), o.seed);
    return 0;
}

static int cmd_hwgate(const Args& a) {
    const auto D = load_design(a.pos(0));
    const auto R = calibrated_report(D, a);
    const std::string& data = a.pos(1);
    std::ifstream f(data);
    if (!f) throw std::runtime_error("cannot open " + data);
    std::stringstream ss; ss << f.rdbuf();
    const auto C = characterisation_from_json(Json::parse(ss.str()));
    const std::string phase = a.text("--phase", D.surface_d >= 5 ? "beta" : "alpha");
    const long shots = long(a.integer("--shots", 50000));
    const auto G = evaluate_hardware_gate(phase, D, R, C, shots, a.seed("--seed", 2026));
    std::printf("\n%s — %s hardware gate from %s (%s data%s%s)\n", D.name.c_str(), phase.c_str(), data.c_str(), C.source.c_str(), C.chip.empty() ? "" : ", ", C.chip.c_str());
    std::printf("  measured: %d of %zu qubits, %d of %d coupled pairs, %zu memory runs\n", G.model.qubits_measured, D.qubits.size(), G.model.pairs_measured, G.model.pairs_designed, C.memory.size());
    for (const auto& c : G.checks) {
        std::printf("  [%s] %-62s %-8s", !c.evaluated ? " -- " : (c.pass ? "PASS" : "FAIL"), c.name.c_str(), c.criterion.c_str());
        if (c.evaluated) std::printf("  %.4g", c.value);
        if (std::isfinite(c.worst)) std::printf("  (worst %.4g)", c.worst);
        if (!c.note.empty()) std::printf("  %s", c.note.c_str());
        std::printf("\n");
    }
    for (const auto& c : G.comparisons)
        std::printf("  model vs measurement, d=%d memory-%c: measured %.4f%%, model with measured parameters %.4f%%, ratio %.2f (%d round counts)\n",
                    c.d, c.basis, 100 * c.eps_measured, 100 * c.eps_model, c.ratio, c.points);
    for (const auto& n : G.model.notes) std::printf("  note: %s\n", n.c_str());
    std::printf("  %s: %d of %d evaluated criteria met, %d not in the data\n", G.pass ? "GATE PASSED" : "GATE NOT PASSED", G.passed, G.evaluated, G.missing);
    if (const auto out = a.text("--out"); !out.empty()) { write_json_file(out, hardware_gate_json(G)); std::printf("  wrote %s\n", out.c_str()); }
    return G.pass ? 0 : 3;
}

// ================================================================== option table
// The options each command reads. Args::parse rejects any other option, and Args refuses to look up one that is not
// listed, so the table and the cmd_* functions above cannot drift apart unnoticed.
static const char* const kNoiseSources = "T1|Tphi|thermal|zz|exchange|control|readout|reset|leakage|all";

static const std::vector<Command>& commands() {
    static const std::vector<Command> table = {
        {"new", cmd_new, {"<chip.json>"},
         {grid("--grid"), integer("--line", 1), integer("--ring", 2), integer("--surface", 3), real("--f-hi"), real("--f-lo"),
          real("--alpha"), real("--g-MHz"), choice("--arch", "fixed|tunable"), real("--Q"), choice("--plan", "2|4|8"),
          real("--stagger-MHz"), text("--name"), integer("--feedlines", 1), real("--kappa-MHz"), real("--chi-MHz"),
          real("--purcell"), real("--readout-ns"), real("--reset-ns"), real("--reset-err", 0, 1), real("--err2q", 0, 1),
          real("--residual-zz-kHz"), flag("--active-reset"),
          choice("--preset", "alpha1|alpha1v11|beta1|gamma1|gen0"), integer("--distance", 3)}},
        {"analyze", cmd_analyze, {"<chip.json>"}, {text("--materials"), flag("--loss-tangent")}},
        {"run", cmd_run, {"<chip.json>", "<circuit>"},
         {integer("--traj", 1), integer("--shots", 0), seed("--seed"), choice_list("--off", kNoiseSources), integer("--top", 0)}},
        {"budget", cmd_budget, {"<chip.json>", "<circuit>"}, {integer("--traj", 1), seed("--seed")}},
        {"fab", cmd_fab, {"<chip.json>"}, {integer("--samples", 1), seed("--seed")}},
        {"qec", cmd_qec, {"<chip.json>"},
         {int_list("--rounds", 1), integer("--shots", 1, std::numeric_limits<long>::max()), seed("--seed"),
          choice("--basis", "Z|X|both|z|x"), integer("--statevector", 0), choice_list("--off", kNoiseSources),
          choice("--decoder", "exact|uf|union-find|uf-stream|legacy|greedy"), integer("--window", 2), integer("--commit", 1),
          choice("--leakage-removal", "lrr|dqlr|lrr+dqlr|none"),
          text("--cache"), flag("--herald"), real("--herald-fidelity", 0, 1), real("--herald-false", 0, 1),
          real("--herald-factor", 0, 1), choice("--herald-rule", "footprint|detector"),
          flag("--bursts"), real("--burst-rate", 0), real("--burst-tau-ms", 0), real("--burst-t1-us"), real("--burst-tphi-us"),
          flag("--veto"), integer("--veto-window", 1), real("--veto-radius", 0), real("--veto-false", 0, 1),
          real("--impact-p", 0, 1), real("--burst-radius", 0)}},
        {"synth-chip", cmd_synth, {"<chip.json>", "<data.json>"},
         {seed("--seed"), integer("--shots", 1, std::numeric_limits<long>::max()), int_list("--rounds", 1), flag("--z-only"),
          real("--hidden-dephasing-us"), text("--cache")}},
        {"hwgate", cmd_hwgate, {"<chip.json>", "<data.json>"},
         {choice("--phase", "alpha|beta"), integer("--shots", 1, std::numeric_limits<long>::max()),
          seed("--seed", std::numeric_limits<std::uint64_t>::max()), text("--out"), text("--cache")}},
        {"layout", cmd_layout, {"<chip.json>", "<out.svg>"}, {}},
        {"cz", cmd_cz, {"<chip.json>"}, {text("--out"), text("--cache"), real("--g-qc"), real("--f-on"), real("--rise")}},
        {"gds", cmd_gds, {"<chip.json>", "<out.gds>"}, {text("--cells"), flag("--synthesize"), text("--from"), text("--cell")}},
    };
    return table;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    const std::string cmd = argv[1];
    if (cmd == "--version" || cmd == "-V") { std::printf("qpuforge %s\n", QPUFORGE_VERSION); return 0; }
    const Command* c = nullptr;
    for (const auto& k : commands()) if (cmd == k.name) c = &k;
    if (!c) {
        std::vector<std::string> names;
        if (is_option(cmd)) names = {"--version", "-V"};
        else for (const auto& k : commands()) names.push_back(k.name);
        const auto s = suggest(cmd, names);
        std::fprintf(stderr, "qpuforge: unknown %s '%s'%s\n", is_option(cmd) ? "option" : "command", cmd.c_str(),
                     s.empty() ? "" : (" (did you mean '" + s + "'?)").c_str());
        usage(stderr);
        return 1;
    }
    try {
        return c->run(Args::parse(*c, argc, argv));
    } catch (const UsageError& e) {
        std::fprintf(stderr, "qpuforge %s: %s\n", c->name, e.what());
        usage(stderr, c->name);
        return 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
