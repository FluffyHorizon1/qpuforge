// qpuforge — exact maximum-weight matching in a general graph
//
// Edmonds' blossom algorithm in its primal-dual O(n^3) form on a dense graph with integer weights. The decoder uses
// it for minimum-weight perfect matching of any number of defects: weights are negated against a constant large
// enough that every maximum-weight matching is also a perfect one.
//
// Vertices are 1..n internally; blossoms take the ids n+1..2n. Every array is owned by the object and reused between
// calls, so one matcher per thread decodes without allocating.
#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace qpuforge {

class BlossomMatcher {
public:
    using W = long long;

    // w: n x n, row-major, symmetric. w[i*n+j] > 0 is an edge; <= 0 means no edge. Returns mate[i] (or -1).
    const std::vector<int>& solve(int n_vertices, const std::vector<W>& w) {
        n_ = n_vertices;
        nx_ = n_;
        stride_ = 2 * n_ + 1;
        const std::size_t s = static_cast<std::size_t>(stride_);
        if (g_.size() < s * s) g_.resize(s * s);
        if (from_.size() < s * static_cast<std::size_t>(n_ + 1)) from_.resize(s * static_cast<std::size_t>(n_ + 1));
        if (lab_.size() < s) {
            lab_.resize(s); match_.resize(s); slack_.resize(s); st_.resize(s); pa_.resize(s); S_.resize(s); vis_.resize(s);
            flower_.resize(s);
        }
        std::fill(vis_.begin(), vis_.begin() + stride_, 0);
        stamp_ = 0;
        W w_max = 0;
        for (int u = 1; u <= n_; ++u)
            for (int v = 1; v <= n_; ++v) {
                const W x = (u == v) ? 0 : std::max<W>(0, w[static_cast<std::size_t>(u - 1) * static_cast<std::size_t>(n_) + static_cast<std::size_t>(v - 1)]);
                g(u, v) = Edge{u, v, x};
                from(u, v) = (u == v) ? u : 0;
                w_max = std::max(w_max, x);
            }
        for (int u = 0; u < stride_; ++u) { match_[idx(u)] = 0; st_[idx(u)] = (u <= n_) ? u : 0; flower_[idx(u)].clear(); }
        lab_[0] = 0;
        for (int u = 1; u <= n_; ++u) lab_[idx(u)] = w_max;
        while (phase()) {}
        mate_.assign(static_cast<std::size_t>(n_), -1);
        for (int u = 1; u <= n_; ++u) if (match_[idx(u)]) mate_[idx(u - 1)] = match_[idx(u)] - 1;
        return mate_;
    }

private:
    struct Edge { int u{0}, v{0}; W w{0}; };
    int n_{0}, nx_{0}, stride_{1}, stamp_{0};
    std::vector<Edge> g_;
    std::vector<int> from_;                 // from(b, x): the petal of blossom b that contains vertex x
    std::vector<W> lab_;
    std::vector<int> match_, slack_, st_, pa_, S_, vis_, queue_, mate_;
    std::vector<std::vector<int>> flower_;
    std::size_t qhead_{0};

    static std::size_t idx(int i) { return static_cast<std::size_t>(i); }
    Edge& g(int u, int v) { return g_[idx(u) * idx(stride_) + idx(v)]; }
    int& from(int b, int x) { return from_[idx(b) * idx(n_ + 1) + idx(x)]; }
    W delta(const Edge& e) const { return lab_[idx(e.u)] + lab_[idx(e.v)] - 2 * e.w; }

    void update_slack(int u, int x) {
        if (!slack_[idx(x)] || delta(g(u, x)) < delta(g(slack_[idx(x)], x))) slack_[idx(x)] = u;
    }
    void set_slack(int x) {
        slack_[idx(x)] = 0;
        for (int u = 1; u <= n_; ++u)
            if (g(u, x).w > 0 && st_[idx(u)] != x && S_[idx(st_[idx(u)])] == 0) update_slack(u, x);
    }
    void q_push(int x) {
        if (x <= n_) { queue_.push_back(x); return; }
        for (std::size_t i = 0; i < flower_[idx(x)].size(); ++i) q_push(flower_[idx(x)][i]);
    }
    void set_st(int x, int b) {
        st_[idx(x)] = b;
        if (x <= n_) return;
        for (std::size_t i = 0; i < flower_[idx(x)].size(); ++i) set_st(flower_[idx(x)][i], b);
    }
    int get_pr(int b, int xr) {
        auto& f = flower_[idx(b)];
        const int pr = static_cast<int>(std::find(f.begin(), f.end(), xr) - f.begin());
        if (pr % 2 == 1) { std::reverse(f.begin() + 1, f.end()); return static_cast<int>(f.size()) - pr; }
        return pr;
    }
    void set_match(int u, int v) {
        match_[idx(u)] = g(u, v).v;
        if (u <= n_) return;
        const Edge e = g(u, v);
        const int xr = from(u, e.u), pr = get_pr(u, xr);
        for (int i = 0; i < pr; ++i) set_match(flower_[idx(u)][idx(i)], flower_[idx(u)][idx(i ^ 1)]);
        set_match(xr, v);
        auto& f = flower_[idx(u)];
        std::rotate(f.begin(), f.begin() + pr, f.end());
    }
    void augment(int u, int v) {
        for (;;) {
            const int xnv = st_[idx(match_[idx(u)])];
            set_match(u, v);
            if (!xnv) return;
            set_match(xnv, st_[idx(pa_[idx(xnv)])]);
            u = st_[idx(pa_[idx(xnv)])]; v = xnv;
        }
    }
    int get_lca(int u, int v) {
        for (++stamp_; u || v; std::swap(u, v)) {
            if (u == 0) continue;
            if (vis_[idx(u)] == stamp_) return u;
            vis_[idx(u)] = stamp_;
            u = st_[idx(match_[idx(u)])];
            if (u) u = st_[idx(pa_[idx(u)])];
        }
        return 0;
    }
    void add_blossom(int u, int lca, int v) {
        int b = n_ + 1;
        while (b <= nx_ && st_[idx(b)]) ++b;
        if (b > nx_) ++nx_;
        lab_[idx(b)] = 0; S_[idx(b)] = 0;
        match_[idx(b)] = match_[idx(lca)];
        auto& f = flower_[idx(b)];
        f.clear();
        f.push_back(lca);
        for (int x = u, y; x != lca; x = st_[idx(pa_[idx(y)])]) { f.push_back(x); y = st_[idx(match_[idx(x)])]; f.push_back(y); q_push(y); }
        std::reverse(f.begin() + 1, f.end());
        for (int x = v, y; x != lca; x = st_[idx(pa_[idx(y)])]) { f.push_back(x); y = st_[idx(match_[idx(x)])]; f.push_back(y); q_push(y); }
        set_st(b, b);
        for (int x = 1; x <= nx_; ++x) { g(b, x) = Edge{}; g(x, b) = Edge{}; }   // no edge yet (and no stale endpoints)
        for (int x = 1; x <= n_; ++x) from(b, x) = 0;
        for (std::size_t i = 0; i < f.size(); ++i) {
            const int xs = f[i];
            for (int x = 1; x <= nx_; ++x)
                if (g(b, x).w == 0 || delta(g(xs, x)) < delta(g(b, x))) { g(b, x) = g(xs, x); g(x, b) = g(x, xs); }
            for (int x = 1; x <= n_; ++x) if (from(xs, x)) from(b, x) = xs;
        }
        set_slack(b);
    }
    void expand_blossom(int b) {
        auto& f = flower_[idx(b)];
        for (std::size_t i = 0; i < f.size(); ++i) set_st(f[i], f[i]);
        const int xr = from(b, g(b, pa_[idx(b)]).u), pr = get_pr(b, xr);
        for (int i = 0; i < pr; i += 2) {
            const int xs = f[idx(i)], xns = f[idx(i + 1)];
            pa_[idx(xs)] = g(xns, xs).u;
            S_[idx(xs)] = 1; S_[idx(xns)] = 0;
            slack_[idx(xs)] = 0; set_slack(xns);
            q_push(xns);
        }
        S_[idx(xr)] = 1; pa_[idx(xr)] = pa_[idx(b)];
        for (std::size_t i = idx(pr + 1); i < f.size(); ++i) { const int xs = f[i]; S_[idx(xs)] = -1; set_slack(xs); }
        st_[idx(b)] = 0;
    }
    bool on_found_edge(const Edge& e) {
        const int u = st_[idx(e.u)], v = st_[idx(e.v)];
        if (S_[idx(v)] == -1) {
            pa_[idx(v)] = e.u; S_[idx(v)] = 1;
            const int nu = st_[idx(match_[idx(v)])];
            slack_[idx(v)] = 0; slack_[idx(nu)] = 0;
            S_[idx(nu)] = 0; q_push(nu);
        } else if (S_[idx(v)] == 0) {
            const int lca = get_lca(u, v);
            if (!lca) { augment(u, v); augment(v, u); return true; }
            add_blossom(u, lca, v);
        }
        return false;
    }
    // one augmentation; false when no further augmenting path has positive gain
    bool phase() {
        std::fill(S_.begin() + 1, S_.begin() + nx_ + 1, -1);
        std::fill(slack_.begin() + 1, slack_.begin() + nx_ + 1, 0);
        queue_.clear(); qhead_ = 0;
        for (int x = 1; x <= nx_; ++x)
            if (st_[idx(x)] == x && !match_[idx(x)]) { pa_[idx(x)] = 0; S_[idx(x)] = 0; q_push(x); }
        if (queue_.empty()) return false;
        for (;;) {
            while (qhead_ < queue_.size()) {
                const int u = queue_[qhead_++];
                if (S_[idx(st_[idx(u)])] == 1) continue;
                for (int v = 1; v <= n_; ++v)
                    if (g(u, v).w > 0 && st_[idx(u)] != st_[idx(v)]) {
                        if (delta(g(u, v)) == 0) { if (on_found_edge(g(u, v))) return true; }
                        else update_slack(u, st_[idx(v)]);
                    }
            }
            W d = std::numeric_limits<W>::max();
            for (int b = n_ + 1; b <= nx_; ++b)
                if (st_[idx(b)] == b && S_[idx(b)] == 1) d = std::min(d, lab_[idx(b)] / 2);
            for (int x = 1; x <= nx_; ++x)
                if (st_[idx(x)] == x && slack_[idx(x)]) {
                    if (S_[idx(x)] == -1) d = std::min(d, delta(g(slack_[idx(x)], x)));
                    else if (S_[idx(x)] == 0) d = std::min(d, delta(g(slack_[idx(x)], x)) / 2);
                }
            if (d == std::numeric_limits<W>::max()) return false;   // nothing left to grow towards
            for (int u = 1; u <= n_; ++u) {
                if (S_[idx(st_[idx(u)])] == 0) {
                    if (lab_[idx(u)] <= d) return false;
                    lab_[idx(u)] -= d;
                } else if (S_[idx(st_[idx(u)])] == 1) lab_[idx(u)] += d;
            }
            for (int b = n_ + 1; b <= nx_; ++b)
                if (st_[idx(b)] == b) {
                    if (S_[idx(b)] == 0) lab_[idx(b)] += d * 2;
                    else if (S_[idx(b)] == 1) lab_[idx(b)] -= d * 2;
                }
            queue_.clear(); qhead_ = 0;
            for (int x = 1; x <= nx_; ++x)
                if (st_[idx(x)] == x && slack_[idx(x)] && st_[idx(slack_[idx(x)])] != x && delta(g(slack_[idx(x)], x)) == 0)
                    if (on_found_edge(g(slack_[idx(x)], x))) return true;
            for (int b = n_ + 1; b <= nx_; ++b)
                if (st_[idx(b)] == b && S_[idx(b)] == 1 && lab_[idx(b)] == 0) expand_blossom(b);
        }
    }
};

} // namespace qpuforge
