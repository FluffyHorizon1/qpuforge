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
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
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
        prepare(ws);
        ws.correction.clear(); ws.iterations = 0;
        if (scaled) for (auto [k, f] : *scaled) ws.len[std::size_t(k)] = std::max<long long>(0, std::llround(kScale * f * edge[std::size_t(k)].w));
        struct Restore {   // put the lengths back however this call ends
            const UnionFindDecoder& U; UFScratch& w; const Scaled* sc;
            ~Restore() { if (sc) for (const auto& kf : *sc) w.len[std::size_t(kf.first)] = U.edge[std::size_t(kf.first)].len; }
        } restore{*this, ws, scaled};
        if (defects.empty()) return 0;
        for (int d : defects) { touch(ws, d); ws.reached[std::size_t(d)] = 1; ws.parity[std::size_t(d)] ^= 1; ws.mark[std::size_t(d)] ^= 1; }
        if (scaled)
            for (const auto& kf : *scaled) {
                const int k = kf.first;
                if (ws.len[std::size_t(k)] > 0) continue;
                const Edge& e = edge[std::size_t(k)];
                for (int x : {e.u, e.v}) if (!ws.reached[std::size_t(x)]) { touch(ws, x); ws.reached[std::size_t(x)] = 1; }
                unite(ws, e.u, e.v);
                ws.fadj[std::size_t(e.u)].push_back(k); ws.fadj[std::size_t(e.v)].push_back(k);
            }
        ws.active.clear();
        for (int d : defects) { const int r = find(ws, d); if (ws.parity[std::size_t(r)] && !ws.boundary[std::size_t(r)]) ws.active.push_back(r); }
        std::sort(ws.active.begin(), ws.active.end());
        ws.active.erase(std::unique(ws.active.begin(), ws.active.end()), ws.active.end());
        // ---- growth
        while (!ws.active.empty()) {
            ++ws.iterations;
            if (++ws.stamp == std::numeric_limits<int>::max()) { std::fill(ws.rate_stamp.begin(), ws.rate_stamp.end(), 0); ws.stamp = 1; }
            ws.edge_touched.clear();
            for (int r : ws.active)
                for (int v = ws.members_head[std::size_t(r)]; v >= 0; v = ws.members_next[std::size_t(v)])
                    for (int k : adj_[std::size_t(v)]) {
                        const Edge& e = edge[std::size_t(k)];
                        if (ws.growth[std::size_t(k)] >= ws.len[std::size_t(k)]) continue;
                        const int o = e.u == v ? e.v : e.u;
                        if (ws.reached[std::size_t(o)] && find(ws, o) == r) continue;   // internal edge
                        if (ws.rate_stamp[std::size_t(k)] != ws.stamp) { ws.rate_stamp[std::size_t(k)] = ws.stamp; ws.rate[std::size_t(k)] = 0; ws.edge_touched.push_back(k); }
                        ws.rate[std::size_t(k)] += 1;
                    }
            if (ws.edge_touched.empty()) break;                                    // isolated (cannot happen with a boundary)
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
                unite(ws, e.u, e.v);
                ws.fadj[std::size_t(e.u)].push_back(k); ws.fadj[std::size_t(e.v)].push_back(k);
            }
            ws.next_active.clear();
            for (int r : ws.active) {
                const int f = find(ws, r);
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
        return parity;
    }

private:
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
    static int find(UFScratch& ws, int v) {
        while (ws.parent[std::size_t(v)] != v) { ws.parent[std::size_t(v)] = ws.parent[std::size_t(ws.parent[std::size_t(v)])]; v = ws.parent[std::size_t(v)]; }
        return v;
    }
    static void unite(UFScratch& ws, int a, int b) {
        a = find(ws, a); b = find(ws, b);
        if (a == b) return;
        if (ws.size[std::size_t(a)] < ws.size[std::size_t(b)]) std::swap(a, b);
        ws.parent[std::size_t(b)] = a;
        ws.size[std::size_t(a)] += ws.size[std::size_t(b)];
        ws.parity[std::size_t(a)] ^= ws.parity[std::size_t(b)];
        ws.boundary[std::size_t(a)] |= ws.boundary[std::size_t(b)];
        ws.members_next[std::size_t(ws.members_tail[std::size_t(a)])] = ws.members_head[std::size_t(b)];
        ws.members_tail[std::size_t(a)] = ws.members_tail[std::size_t(b)];
    }
};

} // namespace qpuforge
