// qpuforge — union-find decoder (Delfosse-Nickerson) with weighted growth
//
// The real-time candidate for surface-code decoding: no all-pairs shortest paths, no matching. The exact matching
// decoder (qec.hpp: MatchingDecoder) stays the reference. Cost: each growth step walks every member of each growing
// cluster (this implementation keeps no separate boundary list), so a call costs O(sum over steps of the growing
// clusters' sizes x degree): linear in the defects when clusters stay small, as they do below threshold, and
// quadratic in the cluster size in the worst case. Union-find with boundary lists reaches almost-linear time.
//
//   growth   every cluster with an odd number of defects that does not touch the boundary grows along all edges on
//            its border. Edges have integer lengths proportional to their weight log((1-p)/p), so a likely edge is
//            crossed sooner than an unlikely one (weighted growth, as in Huang, Newman and Brown 2020). Growth jumps
//            straight to the next time an edge becomes fully grown. A fully grown edge joins the clusters at its
//            ends (union-find with path halving and union by size).
//   peeling  the fully grown edges of each cluster are an erasure; a spanning forest of them is peeled from the
//            leaves inwards (Delfosse and Zemor), rooted at the boundary where the cluster touches it. The edges
//            peeled form the correction, and the XOR of their logical-flip bits is the prediction.
//
// The boundary is one vertex of the graph (as in the matching decoder), so every cluster that reaches it is neutral.
// See docs/design/union-find-decoder.md.
#pragma once
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace qpuforge {

struct UFScratch {
    std::vector<int> parent, size, parity, members_head, members_next, members_tail;
    std::vector<std::uint8_t> boundary, reached, mark, visited, dirty;
    std::vector<long long> growth, len;   // len: edge lengths of this call (heralded edges shortened)
    std::vector<int> rate, rate_stamp, edge_touched, grown, vtouched, active, next_active, full, order, par_edge;
    std::vector<std::vector<int>> fadj;
    std::vector<int> correction;           // edges of the last correction (indices into UnionFindDecoder::edge)
    int stamp{0};
    long iterations{0};                    // growth steps of the last call
    const void* owner{nullptr};            // the decoder this scratch was set up for (edge lengths are cached)
};

// Operation counts of union-find decoding, for the hardware cost model (docs/design/streaming-union-find-fpga.md).
// They are filled only by the counting entry points (UnionFindDecoder::decode_counted, StreamOptions::counts); the
// plain decode() is compiled with a counter whose calls are empty, so it does no extra work.
//   loads        defects loaded (syndrome bits set in the vertex table)
//   growth_iters growth steps (jumps to the next time an edge becomes fully grown)
//   edge_visits  (vertex, incident edge) pairs examined while collecting the border edges of the growing clusters
//   edge_updates border edges whose growth was advanced (one per edge and growth step)
//   finds        root look-ups; hops: parent-pointer dereferences made by them (path halving)
//   unions       merges of two different clusters
//   peel_steps   spanning-forest steps: edges examined by the breadth-first search plus vertices peeled
//   corrections  correction edges produced; commits: those committed (streaming decoder only)
//   pre_pairs    isolated defect pairs removed by the pre-decoder (streaming decoder, StreamOptions::predecode)
struct UFOpCounts {
    long loads{0}, growth_iters{0}, edge_visits{0}, edge_updates{0}, finds{0}, hops{0}, unions{0}, peel_steps{0}, corrections{0}, commits{0},
        pre_pairs{0};
    static constexpr int kFields = 11;
    static const char* name(int i) {
        static const char* n[kFields] = {"loads", "growth_iters", "edge_visits", "edge_updates", "finds", "hops", "unions", "peel_steps",
                                         "corrections", "commits", "pre_pairs"};
        return n[i];
    }
    long get(int i) const {
        const long v[kFields] = {loads, growth_iters, edge_visits, edge_updates, finds, hops, unions, peel_steps, corrections, commits, pre_pairs};
        return v[i];
    }
};
namespace uf_detail {
struct NoCount {
    void load(long) {} void iter() {} void visit() {} void update(long) {} void find() {} void hop() {} void uni() {} void peel(long) {}
    void corr(long) {}
};
struct Count {
    UFOpCounts* c;
    void load(long n) { c->loads += n; } void iter() { c->growth_iters += 1; } void visit() { c->edge_visits += 1; }
    void update(long n) { c->edge_updates += n; } void find() { c->finds += 1; } void hop() { c->hops += 1; } void uni() { c->unions += 1; }
    void peel(long n) { c->peel_steps += n; } void corr(long n) { c->corrections += n; }
};
} // namespace uf_detail

class UnionFindDecoder {
public:
    struct Edge { int u, v; long long len; std::uint8_t obs; double w; };
    static constexpr double kScale = 16.0;     // integer edge length = round(kScale * weight), at least 1

    UnionFindDecoder() = default;
    // n_det detectors; the boundary is vertex n_det. Edges: (a, b or -1 for the boundary, flips the logical, p).
    // Parallel edges keep the most likely one (as the matching decoder does).
    template <class EdgeList>
    UnionFindDecoder(int n_det, const EdgeList& edges) : N_(n_det), V_(n_det + 1) {
        auto& at = index_;
        for (const auto& e : edges) {
            const double p = std::clamp(e.p, 1e-15, 0.5 - 1e-12);
            const double w = std::log((1 - p) / p);
            int u = e.a, v = e.b < 0 ? N_ : e.b;
            if (u == v) continue;
            if (u > v) std::swap(u, v);
            auto it = at.find({u, v});
            if (it == at.end()) { at[{u, v}] = edge.size(); edge.push_back({u, v, 0, std::uint8_t(e.obs), w}); }
            else if (w < edge[it->second].w) { edge[it->second].w = w; edge[it->second].obs = std::uint8_t(e.obs); }
        }
        adj_.assign(std::size_t(V_), {});
        for (std::size_t k = 0; k < edge.size(); ++k) {
            edge[k].len = std::max<long long>(1, std::llround(kScale * edge[k].w));
            adj_[std::size_t(edge[k].u)].push_back(int(k));
            adj_[std::size_t(edge[k].v)].push_back(int(k));
        }
    }
    int n_detectors() const { return N_; }
    const std::vector<int>& incident(int v) const { return adj_[std::size_t(v < 0 ? N_ : v)]; }
    // index of the edge between detectors a and b (b = -1: the boundary), or -1
    int edge_index(int a, int b) const {
        if (b < 0) b = N_;
        auto it = index_.find({std::min(a, b), std::max(a, b)});
        return it == index_.end() ? -1 : int(it->second);
    }
    std::vector<Edge> edge;

    // predicted logical flip for a set of fired detectors (each listed once); the correction is left in ws.correction
    // scaled (optional): (edge, factor) pairs, each edge once: its length is multiplied by the factor for this call
    // (heralded edges); a length that rounds to 0 is an erasure: the edge is grown from the start.
    using Scaled = std::vector<std::pair<int, double>>;
    int decode(const std::vector<int>& defects, UFScratch& ws, const Scaled* scaled = nullptr) const {
        uf_detail::NoCount c;
        return decode_impl(defects, ws, scaled, c);
    }
    // the same decoding (same correction, same prediction) with its operations added to `counts`
    int decode_counted(const std::vector<int>& defects, UFScratch& ws, UFOpCounts& counts, const Scaled* scaled = nullptr) const {
        uf_detail::Count c{&counts};
        return decode_impl(defects, ws, scaled, c);
    }

private:
    template <class Ctr>
    int decode_impl(const std::vector<int>& defects, UFScratch& ws, const Scaled* scaled, Ctr& c) const {
        prepare(ws);
        ws.correction.clear(); ws.iterations = 0;
        if (scaled) for (auto [k, f] : *scaled) ws.len[std::size_t(k)] = std::max<long long>(0, std::llround(kScale * f * edge[std::size_t(k)].w));
        struct Restore {   // put the lengths back however this call ends
            const UnionFindDecoder& U; UFScratch& w; const Scaled* sc;
            ~Restore() { if (sc) for (const auto& kf : *sc) w.len[std::size_t(kf.first)] = U.edge[std::size_t(kf.first)].len; }
        } restore{*this, ws, scaled};
        if (defects.empty()) return 0;
        c.load(long(defects.size()));
        for (int d : defects) { touch(ws, d); ws.reached[std::size_t(d)] = 1; ws.parity[std::size_t(d)] ^= 1; ws.mark[std::size_t(d)] ^= 1; }
        if (scaled)
            for (const auto& kf : *scaled) {
                const int k = kf.first;
                if (ws.len[std::size_t(k)] > 0) continue;
                const Edge& e = edge[std::size_t(k)];
                for (int x : {e.u, e.v}) if (!ws.reached[std::size_t(x)]) { touch(ws, x); ws.reached[std::size_t(x)] = 1; }
                unite(ws, e.u, e.v, c);
                ws.fadj[std::size_t(e.u)].push_back(k); ws.fadj[std::size_t(e.v)].push_back(k);
            }
        ws.active.clear();
        for (int d : defects) { const int r = find(ws, d, c); if (ws.parity[std::size_t(r)] && !ws.boundary[std::size_t(r)]) ws.active.push_back(r); }
        std::sort(ws.active.begin(), ws.active.end());
        ws.active.erase(std::unique(ws.active.begin(), ws.active.end()), ws.active.end());
        // ---- growth
        while (!ws.active.empty()) {
            ++ws.iterations; c.iter();
            if (++ws.stamp == std::numeric_limits<int>::max()) { std::fill(ws.rate_stamp.begin(), ws.rate_stamp.end(), 0); ws.stamp = 1; }
            ws.edge_touched.clear();
            for (int r : ws.active)
                for (int v = ws.members_head[std::size_t(r)]; v >= 0; v = ws.members_next[std::size_t(v)])
                    for (int k : adj_[std::size_t(v)]) {
                        const Edge& e = edge[std::size_t(k)];
                        c.visit();
                        if (ws.growth[std::size_t(k)] >= ws.len[std::size_t(k)]) continue;
                        const int o = e.u == v ? e.v : e.u;
                        if (ws.reached[std::size_t(o)] && find(ws, o, c) == r) continue;   // internal edge
                        if (ws.rate_stamp[std::size_t(k)] != ws.stamp) { ws.rate_stamp[std::size_t(k)] = ws.stamp; ws.rate[std::size_t(k)] = 0; ws.edge_touched.push_back(k); }
                        ws.rate[std::size_t(k)] += 1;
                    }
            if (ws.edge_touched.empty()) break;                                    // isolated (cannot happen with a boundary)
            c.update(long(ws.edge_touched.size()));
            long long dt = std::numeric_limits<long long>::max();
            for (int k : ws.edge_touched) {
                const long long rem = ws.len[std::size_t(k)] - ws.growth[std::size_t(k)], rt = ws.rate[std::size_t(k)];
                dt = std::min(dt, (rem + rt - 1) / rt);
            }
            ws.full.clear();
            for (int k : ws.edge_touched) {
                if (ws.growth[std::size_t(k)] == 0) ws.grown.push_back(k);
                ws.growth[std::size_t(k)] += dt * ws.rate[std::size_t(k)];
                if (ws.growth[std::size_t(k)] >= ws.len[std::size_t(k)]) ws.full.push_back(k);
            }
            std::sort(ws.full.begin(), ws.full.end());                             // deterministic merge order
            for (int k : ws.full) {
                const Edge& e = edge[std::size_t(k)];
                for (int x : {e.u, e.v}) if (!ws.reached[std::size_t(x)]) { touch(ws, x); ws.reached[std::size_t(x)] = 1; }
                unite(ws, e.u, e.v, c);
                ws.fadj[std::size_t(e.u)].push_back(k); ws.fadj[std::size_t(e.v)].push_back(k);
            }
            ws.next_active.clear();
            for (int r : ws.active) {
                const int f = find(ws, r, c);
                if (ws.parity[std::size_t(f)] && !ws.boundary[std::size_t(f)]) ws.next_active.push_back(f);
            }
            std::sort(ws.next_active.begin(), ws.next_active.end());
            ws.next_active.erase(std::unique(ws.next_active.begin(), ws.next_active.end()), ws.next_active.end());
            ws.active.swap(ws.next_active);
        }
        // ---- peeling: spanning forest of the grown edges, rooted at the boundary first
        ws.correction.clear();
        ws.order.clear();
        auto bfs = [&](int root) {
            std::size_t head = ws.order.size();
            ws.visited[std::size_t(root)] = 1; ws.par_edge[std::size_t(root)] = -1; ws.order.push_back(root);
            while (head < ws.order.size()) {
                const int v = ws.order[head++];
                c.peel(long(ws.fadj[std::size_t(v)].size()));
                for (int k : ws.fadj[std::size_t(v)]) {
                    const Edge& e = edge[std::size_t(k)];
                    const int o = e.u == v ? e.v : e.u;
                    if (ws.visited[std::size_t(o)]) continue;
                    ws.visited[std::size_t(o)] = 1; ws.par_edge[std::size_t(o)] = k; ws.order.push_back(o);
                }
            }
        };
        if (ws.reached[std::size_t(N_)]) bfs(N_);
        for (int v : ws.vtouched) if (!ws.visited[std::size_t(v)]) bfs(v);
        int parity = 0;
        c.peel(long(ws.order.size()));
        for (std::size_t i = ws.order.size(); i-- > 0;) {
            const int v = ws.order[i];
            if (!ws.mark[std::size_t(v)]) continue;
            const int k = ws.par_edge[std::size_t(v)];
            if (k < 0) continue;                                                   // a root: the boundary (or a failed cluster)
            const Edge& e = edge[std::size_t(k)];
            ws.mark[std::size_t(v)] = 0;
            const int o = e.u == v ? e.v : e.u;
            ws.mark[std::size_t(o)] ^= 1;
            parity ^= e.obs;
            ws.correction.push_back(k);
        }
        c.corr(long(ws.correction.size()));
        return parity;
    }

    int N_{0}, V_{1};
    std::vector<std::vector<int>> adj_;
    std::map<std::pair<int, int>, std::size_t> index_;

    void prepare(UFScratch& ws) const {
        const std::size_t V = std::size_t(V_), E = edge.size();
        if (ws.owner != this) { ws = UFScratch{}; ws.owner = this; }   // set up for another graph: start from scratch
        if (ws.parent.size() != V) {
            ws.parent.assign(V, 0); ws.size.assign(V, 1); ws.parity.assign(V, 0); ws.members_head.assign(V, -1); ws.members_next.assign(V, -1);
            ws.members_tail.assign(V, -1); ws.boundary.assign(V, 0); ws.reached.assign(V, 0); ws.mark.assign(V, 0); ws.visited.assign(V, 0);
            ws.dirty.assign(V, 0); ws.par_edge.assign(V, -1); ws.fadj.assign(V, {});
            for (std::size_t v = 0; v < V; ++v) { ws.parent[v] = int(v); ws.members_head[v] = ws.members_tail[v] = int(v); }
            ws.vtouched.clear();
        }
        if (ws.growth.size() != E) {
            ws.growth.assign(E, 0); ws.rate.assign(E, 0); ws.rate_stamp.assign(E, 0); ws.stamp = 0; ws.grown.clear();
            ws.len.resize(E);
            for (std::size_t k = 0; k < E; ++k) ws.len[k] = edge[k].len;
        }
        // undo the previous call: only the vertices and edges it touched
        for (int v : ws.vtouched) {
            const std::size_t s = std::size_t(v);
            ws.parent[s] = v; ws.size[s] = 1; ws.parity[s] = 0; ws.members_head[s] = ws.members_tail[s] = v; ws.members_next[s] = -1;
            ws.boundary[s] = 0; ws.reached[s] = 0; ws.mark[s] = 0; ws.visited[s] = 0; ws.dirty[s] = 0; ws.par_edge[s] = -1;
            ws.fadj[s].clear();
        }
        for (int k : ws.grown) ws.growth[std::size_t(k)] = 0;
        ws.grown.clear();
        ws.vtouched.clear();
        ws.boundary[std::size_t(N_)] = 1;
        touch(ws, N_);                                                             // so that its state is undone
    }
    static void touch(UFScratch& ws, int v) {
        if (ws.dirty[std::size_t(v)]) return;
        ws.dirty[std::size_t(v)] = 1;
        ws.vtouched.push_back(v);
    }
    template <class Ctr>
    static int find(UFScratch& ws, int v, Ctr& c) {
        c.find();
        while (ws.parent[std::size_t(v)] != v) { c.hop(); ws.parent[std::size_t(v)] = ws.parent[std::size_t(ws.parent[std::size_t(v)])]; v = ws.parent[std::size_t(v)]; }
        return v;
    }
    template <class Ctr>
    static void unite(UFScratch& ws, int a, int b, Ctr& c) {
        a = find(ws, a, c); b = find(ws, b, c);
        if (a == b) return;
        c.uni();
        if (ws.size[std::size_t(a)] < ws.size[std::size_t(b)]) std::swap(a, b);
        ws.parent[std::size_t(b)] = a;
        ws.size[std::size_t(a)] += ws.size[std::size_t(b)];
        ws.parity[std::size_t(a)] ^= ws.parity[std::size_t(b)];
        ws.boundary[std::size_t(a)] |= ws.boundary[std::size_t(b)];
        ws.members_next[std::size_t(ws.members_tail[std::size_t(a)])] = ws.members_head[std::size_t(b)];
        ws.members_tail[std::size_t(a)] = ws.members_tail[std::size_t(b)];
    }
};

// A cycle-count model of one window decode on synchronous hardware (docs/design/streaming-union-find-fpga.md): a
// window without defects costs `empty` cycles, one that the pre-decoder clears costs `cleared`, and any other costs
// `fixed` plus the operation counts weighted by `per` (indexed as UFOpCounts::get). It is a model of a design, not a
// measurement of one.
struct UFCostModel {
    std::string name;
    double clock_mhz{250};
    long empty{1}, cleared{2}, fixed{0};
    long per[UFOpCounts::kFields]{};
    // The two designs of the note. Both: 250 MHz; a window without defects is one cycle (an OR over the syndrome
    // register); the pre-decoder is combinational over the window and costs two cycles when it clears it.
    //   serial_bram      one operation at a time on block-RAM tables (vertex table, edge table; true dual port, one read
    //                    and one write per table per cycle, registered output: 2 cycles from address to data; pipeline
    //                    depth 4). Independent accesses stream at one per cycle; a dependent access (a parent-pointer
    //                    hop) waits the 2-cycle read. Per defect loaded 1; per edge visit 1; per find 1 plus 2 per hop;
    //                    per border-edge update 2 (the minimum pass and the add pass); per growth step 8 (the pipeline
    //                    drains twice: before the step length is known and before the merges); per union 2; per
    //                    peeling step 1; per correction edge 1 (commit test and syndrome toggle); 12 fixed (three
    //                    phases, each filling the pipeline).
    //   parallel_growth  one cell per vertex (root register, parity, boundary bit) and per edge (growth counter and
    //                    length), all in registers. A growth step is 5 cycles whatever the cluster sizes: border edges
    //                    flagged (1), minimum remaining length over them in a comparator tree (2), add (1), fully grown
    //                    edges found (1). A union is 2 cycles: the absorbed root is broadcast and every vertex cell
    //                    compares and rewrites its root, so roots stay flat and a find is a register read (0). Loading
    //                    and committing act on all cells at once (4 fixed). Peeling stays serial: 1 per step.
    static UFCostModel serial_bram() {
        UFCostModel m; m.name = "serial, block RAM"; m.fixed = 12;
        const long per[UFOpCounts::kFields] = {1, 8, 1, 2, 1, 2, 2, 1, 1, 0, 0};
        for (int i = 0; i < UFOpCounts::kFields; ++i) m.per[i] = per[i];
        return m;
    }
    static UFCostModel parallel_growth() {
        UFCostModel m; m.name = "parallel growth"; m.fixed = 4;
        const long per[UFOpCounts::kFields] = {0, 5, 0, 0, 0, 0, 2, 1, 0, 0, 0};
        for (int i = 0; i < UFOpCounts::kFields; ++i) m.per[i] = per[i];
        return m;
    }
    long cycles(const UFOpCounts& c, int cleared_state) const {
        if (cleared_state == 0) return empty;
        if (cleared_state == 1) return cleared;
        long n = fixed;
        for (int i = 0; i < UFOpCounts::kFields; ++i) n += per[i] * c.get(i);
        return n;
    }
};

// ================================================================== streaming (sliding-window) union-find
// The syndrome arrives round by round. A window of `window` rounds is decoded with the union-find decoder above, with
// the rounds after the window treated as a boundary (edges into them go to the boundary vertex; edges into rounds
// already committed are dropped). The correction edges that touch the first `commit` rounds of the window are
// committed: their logical bits enter the prediction and their far ends are toggled in the carried syndrome. The
// rest is discarded and the window moves on by `commit` rounds. The last window commits everything. The decoder
// needs window > commit, so that a committed edge never reaches past the window (edges span at most one round).
// See docs/design/streaming-union-find.md.
struct StreamScratch {
    std::vector<UFScratch> uf;                     // one per window: a scratch caches its graph's edge lengths
    std::vector<std::uint8_t> syn, lsyn;
    std::vector<int> def, ndef, partner;
    std::vector<std::pair<int, int>> committed;   // global (a, b) of every committed edge; b = -1: the boundary
    std::vector<UnionFindDecoder::Scaled> wscaled; // per window: re-weighted local edges of this call (heralds)
    std::vector<int> wscaled_used;
    std::vector<UFOpCounts> counts;                // per window, when StreamOptions::count is set
    std::vector<std::uint8_t> cleared;             // per window: 0 no defects, 1 cleared by the pre-decoder alone, 2 union-find ran
    double window_s_max{0}, window_s_sum{0};       // wall-clock time per window decode, this call
    int windows{0};
};

// Optional behaviour of StreamingUnionFind::decode.
//   scaled     re-weighted edges of this call (heralds): (edge index of global().edge, factor in [0, 1]), each edge once.
//              Every window that holds the edge decodes with the new length. Where a window folds several edges into
//              one (an edge into the future joins the window's boundary edge), the folded edge takes the smallest new
//              weight and is never made longer.
//   predecode  before union-find runs on a window, remove every isolated defect pair: two defects joined by an edge,
//              neither with another defect on any neighbouring vertex, and the edge no longer than the two defects'
//              boundary edges together. The pair's edge is its correction. This is the decision a per-vertex
//              combinational pre-decoder can take (docs/design/streaming-union-find-fpga.md); it is not what
//              union-find would always do, so the logical error can change.
//   count      fill StreamScratch::counts (operation counts per window) and StreamScratch::cleared
struct StreamOptions {
    const UnionFindDecoder::Scaled* scaled{nullptr};
    bool predecode{false};
    bool count{false};
};

class StreamingUnionFind {
public:
    // det_round: round of each detector (0 .. max_round; the final data readout has the largest round)
    template <class EdgeList>
    StreamingUnionFind(int n_det, const EdgeList& edges, const std::vector<int>& det_round, int window, int commit)
        : N_(n_det), W_(window), C_(commit), round_(det_round), glob_(n_det, edges) {
        if (commit < 1 || window <= commit) throw std::runtime_error("StreamingUnionFind: needs 1 <= commit < window");
        if (int(det_round.size()) != n_det) throw std::runtime_error("StreamingUnionFind: det_round size");
        R_ = 0; for (int r : det_round) R_ = std::max(R_, r);
        struct E { int a, b; bool obs; double p; };
        std::vector<std::vector<std::pair<int, int>>> in_window(glob_.edge.size());   // global edge -> (window, local edge)
        for (int t = 0;; t += C_) {
            Win w; w.t = t; w.last = t + W_ > R_;
            const int end = w.last ? R_ + 1 : t + W_;
            std::vector<int> local(std::size_t(n_det), -1);
            for (int d = 0; d < n_det; ++d) if (det_round[std::size_t(d)] >= t && det_round[std::size_t(d)] < end) { local[std::size_t(d)] = int(w.glob.size()); w.glob.push_back(d); }
            std::vector<E> le;
            for (const auto& e : edges) {
                const int la = local[std::size_t(e.a)];
                const bool b_bnd = e.b < 0;
                const int lb = b_bnd ? -1 : local[std::size_t(e.b)];
                const bool a_fut = la < 0 && det_round[std::size_t(e.a)] >= end, b_fut = !b_bnd && lb < 0 && det_round[std::size_t(e.b)] >= end;
                if (la >= 0 && (b_bnd || b_fut)) le.push_back({la, -1, e.obs, e.p});
                else if (lb >= 0 && a_fut) le.push_back({lb, -1, e.obs, e.p});
                else if (la >= 0 && lb >= 0) le.push_back({la, lb, e.obs, e.p});
            }
            w.uf = UnionFindDecoder(int(w.glob.size()), le);
            // the window's edge of each global edge (heralds), and each vertex's boundary edge (pre-decoder)
            const int nl = int(w.glob.size());
            for (std::size_t g = 0; g < glob_.edge.size(); ++g) {
                const auto& ge = glob_.edge[g];
                const int la = local[std::size_t(ge.u)];
                const bool b_bnd = ge.v == N_;
                const int lb = b_bnd ? -1 : local[std::size_t(ge.v)];
                const bool a_fut = la < 0 && det_round[std::size_t(ge.u)] >= end, b_fut = !b_bnd && lb < 0 && det_round[std::size_t(ge.v)] >= end;
                int k = -1;
                if (la >= 0 && (b_bnd || b_fut)) k = w.uf.edge_index(la, -1);
                else if (lb >= 0 && a_fut) k = w.uf.edge_index(lb, -1);
                else if (la >= 0 && lb >= 0) k = w.uf.edge_index(la, lb);
                if (k >= 0) in_window[g].push_back({int(win_.size()), k});
            }
            w.bnd_len.assign(std::size_t(nl), std::numeric_limits<long long>::max() / 4);
            for (const auto& e : w.uf.edge) if (e.v == nl) w.bnd_len[std::size_t(e.u)] = e.len;
            win_.push_back(std::move(w));
            if (win_.back().last) break;
        }
        g_off_.assign(glob_.edge.size() + 1, 0);
        for (std::size_t g = 0; g < in_window.size(); ++g) { g_off_[g + 1] = g_off_[g] + int(in_window[g].size()); for (const auto& x : in_window[g]) g_ent_.push_back(x); }
    }
    int window() const { return W_; }
    int commit() const { return C_; }
    int n_windows() const { return int(win_.size()); }
    // the whole detector graph, indexed as MatchingDecoder::union_find() is when both are built from the same edge list
    const UnionFindDecoder& global() const { return glob_; }
    // vertices and edges of the graph of window i (the boundary vertex not counted)
    int window_vertices(int i) const { return int(win_[std::size_t(i)].glob.size()); }
    int window_edges(int i) const { return int(win_[std::size_t(i)].uf.edge.size()); }
    // rounds whose correction window i makes final (the last window commits everything it holds)
    int committed_rounds(int i) const { const auto& w = win_[std::size_t(i)]; return w.last ? R_ + 1 - w.t : C_; }

    // predicted logical flip; the committed correction is left in ws.committed
    int decode(const std::vector<int>& defects, StreamScratch& ws) const { return decode(defects, ws, StreamOptions{}); }
    int decode(const std::vector<int>& defects, StreamScratch& ws, const StreamOptions& opt) const {
        ws.syn.assign(std::size_t(N_), 0);
        for (int d : defects) ws.syn[std::size_t(d)] ^= 1;
        ws.committed.clear(); ws.window_s_max = 0; ws.window_s_sum = 0; ws.windows = 0;
        int parity = 0;
        if (ws.uf.size() != win_.size()) ws.uf.assign(win_.size(), UFScratch{});
        if (opt.count) { ws.counts.assign(win_.size(), UFOpCounts{}); ws.cleared.assign(win_.size(), 0); }
        // heralds: the re-weighted edges of each window
        if (ws.wscaled.size() != win_.size()) ws.wscaled.assign(win_.size(), {});
        for (int wi : ws.wscaled_used) ws.wscaled[std::size_t(wi)].clear();
        ws.wscaled_used.clear();
        if (opt.scaled)
            for (const auto& [g, f] : *opt.scaled) {
                const double wg = f * glob_.edge[std::size_t(g)].w;
                for (int i = g_off_[std::size_t(g)]; i < g_off_[std::size_t(g) + 1]; ++i) {
                    const auto [wi, k] = g_ent_[std::size_t(i)];
                    const double fl = wg / win_[std::size_t(wi)].uf.edge[std::size_t(k)].w;
                    if (fl >= 1.0) continue;
                    auto& v = ws.wscaled[std::size_t(wi)];
                    if (v.empty()) ws.wscaled_used.push_back(wi);
                    auto it = std::find_if(v.begin(), v.end(), [k = k](const std::pair<int, double>& x) { return x.first == k; });
                    if (it == v.end()) v.push_back({k, fl}); else it->second = std::min(it->second, fl);
                }
            }
        for (std::size_t wi = 0; wi < win_.size(); ++wi) {
            const Win& w = win_[wi];
            UFScratch& us = ws.uf[wi];
            const auto t0 = std::chrono::steady_clock::now();
            const int nl = int(w.glob.size());
            ws.def.clear();
            for (std::size_t i = 0; i < w.glob.size(); ++i) if (ws.syn[std::size_t(w.glob[i])]) ws.def.push_back(int(i));
            auto apply = [&](int k) {   // a correction edge of this window: commit it if it touches the commit region
                const auto& e = w.uf.edge[std::size_t(k)];
                const int ga = e.u < nl ? w.glob[std::size_t(e.u)] : -1, gb = e.v < nl ? w.glob[std::size_t(e.v)] : -1;
                bool commit = w.last;
                for (int g : {ga, gb}) if (g >= 0 && round_[std::size_t(g)] < w.t + C_) commit = true;
                if (!commit) return false;
                parity ^= e.obs;
                for (int g : {ga, gb}) if (g >= 0) ws.syn[std::size_t(g)] ^= 1;
                ws.committed.push_back({ga >= 0 ? ga : gb, ga >= 0 ? gb : -1});
                return true;
            };
            const bool had = !ws.def.empty();
            const long n_def = long(ws.def.size());
            if (opt.predecode && ws.def.size() >= 2) {
                // isolated pairs: partner[v] = the edge to the only defect neighbour of v (-1 none, -2 several)
                ws.lsyn.assign(std::size_t(nl), 0); ws.partner.assign(std::size_t(nl), -1);
                for (int v : ws.def) ws.lsyn[std::size_t(v)] = 1;
                for (int v : ws.def)
                    for (int k : w.uf.incident(v)) {
                        const auto& e = w.uf.edge[std::size_t(k)];
                        const int o = e.u == v ? e.v : e.u;
                        if (o < nl && ws.lsyn[std::size_t(o)]) ws.partner[std::size_t(v)] = ws.partner[std::size_t(v)] == -1 ? k : -2;
                    }
                ws.ndef.clear();
                for (int v : ws.def) {
                    const int k = ws.partner[std::size_t(v)];
                    bool pair = false;
                    if (k >= 0) {
                        const auto& e = w.uf.edge[std::size_t(k)];
                        const int o = e.u == v ? e.v : e.u;
                        pair = ws.partner[std::size_t(o)] == k && e.len <= w.bnd_len[std::size_t(v)] + w.bnd_len[std::size_t(o)];
                        if (pair && v < o) {
                            const bool done = apply(k);
                            if (opt.count) { ws.counts[wi].pre_pairs += 1; ws.counts[wi].commits += done; }
                        }
                    }
                    if (!pair) ws.ndef.push_back(v);
                }
                ws.def.swap(ws.ndef);
            }
            if (opt.count) { ws.cleared[wi] = !had ? 0 : ws.def.empty() ? 1 : 2; ws.counts[wi].loads = n_def; }
            if (!ws.def.empty()) {
                const UnionFindDecoder::Scaled* sc = ws.wscaled[wi].empty() ? nullptr : &ws.wscaled[wi];
                if (opt.count) { w.uf.decode_counted(ws.def, us, ws.counts[wi], sc); ws.counts[wi].loads = n_def; }
                else w.uf.decode(ws.def, us, sc);
                for (int k : us.correction) {
                    const bool done = apply(k);
                    if (opt.count) ws.counts[wi].commits += done;
                }
            }
            const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            ws.window_s_max = std::max(ws.window_s_max, dt); ws.window_s_sum += dt; ++ws.windows;
        }
        return parity;
    }

private:
    struct Win { int t{0}; bool last{false}; std::vector<int> glob; UnionFindDecoder uf; std::vector<long long> bnd_len; };
    int N_, W_, C_, R_{0};
    std::vector<int> round_;
    UnionFindDecoder glob_;
    std::vector<Win> win_;
    std::vector<int> g_off_;                       // global edge -> its (window, local edge) entries in g_ent_
    std::vector<std::pair<int, int>> g_ent_;
};

} // namespace qpuforge
