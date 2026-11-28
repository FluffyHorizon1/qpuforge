// Shared command line, output paths and provenance for the design-review programs.
//
//   <review> [shots] [--out DIR] [--cache FILE] [--smoke]
//
//   shots      Monte Carlo shots per point (default per program)
//   --out DIR  write every output into DIR (by file name) instead of the committed locations under designs/ and results/
//   --cache F  CZ calibration cache (default results/cz_cache.json; with --out, a copy of it inside DIR)
//   --smoke    quick end-to-end run (200 shots unless given, reduced fixed-size runs); implies --out if not given
#pragma once
#include "git_stamp.hpp"
#include "qpuforge/json.hpp"
#include "qpuforge/pulse.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace review {

struct Options {
    long shots{0};
    bool smoke{false};
    std::string out_dir;                       // empty: committed locations
    std::string cache{"results/cz_cache.json"};
    std::vector<std::string> args;
    // committed location -> actual output path
    std::string path(const std::string& committed) const {
        if (out_dir.empty()) return committed;
        return (std::filesystem::path(out_dir) / std::filesystem::path(committed).filename()).string();
    }
    // fixed-size runs (e.g. the d = 7 point) are cut down in smoke mode
    long fixed(long n) const { return smoke ? std::min(n, shots) : n; }
};

inline Options parse(int argc, char** argv, long default_shots) {
    Options o; o.shots = default_shots;
    bool shots_given = false, cache_given = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        o.args.push_back(a);
        if (a == "--smoke") o.smoke = true;
        else if (a == "--out" && i + 1 < argc) { o.out_dir = argv[++i]; o.args.push_back(o.out_dir); }
        else if (a == "--cache" && i + 1 < argc) { o.cache = argv[++i]; cache_given = true; o.args.push_back(o.cache); }
        else if (!a.empty() && a[0] != '-') { o.shots = std::atol(a.c_str()); shots_given = true; }
        else { std::fprintf(stderr, "unknown argument %s\n", a.c_str()); std::exit(2); }
    }
    if (o.smoke && !shots_given) o.shots = 200;
    if (o.smoke && o.out_dir.empty()) o.out_dir = (std::filesystem::temp_directory_path() / "qpuforge_review_smoke").string();
    if (!o.out_dir.empty()) {
        std::filesystem::create_directories(o.out_dir);
        if (!cache_given) {   // work on a copy so a smoke or scratch run never rewrites the committed cache
            const auto dst = std::filesystem::path(o.out_dir) / "cz_cache.json";
            if (std::filesystem::exists(o.cache)) std::filesystem::copy_file(o.cache, dst, std::filesystem::copy_options::overwrite_existing);
            o.cache = dst.string();
        }
    }
    return o;
}

inline std::string compiler() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("g++ ") + __VERSION__;
#else
    return "unknown";
#endif
}

inline std::string utc_now() {
    const std::time_t t = std::time(nullptr);
    char b[32]; std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
    return b;
}

// Provenance block written into every review JSON.
inline qpuforge::Json provenance(const Options& o, const std::string& seeds) {
    using qpuforge::Json;
    Json p = Json::object();
    p.set("source_commit", git_stamp_commit);       // the tree the binary was built from (git_stamp.hpp)
    p.set("source_tree_dirty", git_stamp_dirty);
    p.set("compiler", compiler());
#ifdef _OPENMP
    p.set("openmp_threads", omp_get_max_threads());
#else
    p.set("openmp_threads", 0);
#endif
    p.set("date_utc", utc_now());
    p.set("seeds", seeds);
    p.set("cz_model_version", qpuforge::CZModel::kVersion);
    p.set("cz_cache", o.cache);
    std::string a; for (const auto& s : o.args) a += (a.empty() ? "" : " ") + s;
    p.set("arguments", a);
    p.set("smoke", o.smoke);
    return p;
}

} // namespace review
