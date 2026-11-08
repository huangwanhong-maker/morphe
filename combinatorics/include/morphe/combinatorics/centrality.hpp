// Centralities of a graph whose arcs carry a cost and a length: closeness and betweenness from every source, within a radius
// (Freeman 1977, Sociometry 40, 35–41; Freeman 1979, Social Networks 1, 215–239), by Brandes's algorithm (2001, "A faster
// algorithm for betweenness centrality", J. Math. Sociol. 25, 163–177): one Dijkstra search from each source counting the least
// paths σ, then their dependencies accumulated back.
//
// - **Least paths.** Paths are compared by their cost, and paths of equal cost by their length: the least paths from s to t are
//   those least in (cost, length) lexicographically, to a relative tolerance `tie`. Costs may be zero (a street running straight on,
//   under an angular cost); lengths must be positive. Each arc then strictly raises (cost, length), so the least paths' predecessors
//   form a directed acyclic graph and their counts σ are finite, even where zero-cost arcs form cycles. When cost is length, the
//   least paths are the shortest paths.
// - **The radius.** The search from s does not pass beyond length `radius` along its least paths: a vertex is reached when a least
//   path of the region searched reaches it within that length, and nothing is reached through a vertex beyond it. When cost is
//   length this is the ball of radius `radius` about s; under another cost it is the search's region, as in Space Syntax's analyses
//   at a metric radius (Turner 2007, Environ. Plann. B 34, 539–555).
// - **Closeness.** For each source, the vertices reached (the source excluded), the sum of their least costs (Space Syntax's total
//   depth), and the sum of 1/cost over those at positive cost (harmonic closeness; Marchiori and Latora 2000, Physica A 285).
// - **Betweenness.** For each vertex v, the sum over sources s and targets t reached from s, s ≠ v ≠ t, of σ_st(v)/σ_st, the share of
//   the least paths from s to t through v. Over ordered pairs; with `undirected`, halved, which counts each pair once when every arc
//   has its reverse with the same cost and length (Freeman's and networkx's convention).
//
// Sources may be split: closeness is set for the sources computed, and betweenness adds across disjoint sets of sources (`add`), so
// that a caller may run the sources on several threads. Cost O(S (V + E) log V) for S sources.

#pragma once

#include <morphe/combinatorics/shortest_paths.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace morphe::combinatorics {

struct CostArc {
  std::size_t to = 0;
  double cost = 0;    // ≥ 0: what least paths minimise
  double length = 0;  // > 0: what breaks ties of cost, and what the radius measures
};
using CostAdjacency = std::vector<std::vector<CostArc>>;

struct CentralityOptions {
  double radius = std::numeric_limits<double>::infinity();
  double tie = 1e-9;          // relative: two costs (or lengths) within tie × max(1, |value|) are equal
  bool betweenness = true;
  bool undirected = false;    // halve betweenness: each unordered pair once
};

struct Centrality {
  std::vector<double> reached;      // per source: vertices reached within the radius, the source excluded
  std::vector<double> total_cost;   // per source: the sum of their least costs
  std::vector<double> harmonic;     // per source: the sum of 1/cost over those at positive cost
  std::vector<double> betweenness;  // per vertex

  // Adds another computation over a disjoint set of sources.
  void add(const Centrality& o) {
    for (std::size_t v = 0; v < betweenness.size(); v++) {
      betweenness[v] += o.betweenness[v];
      reached[v] += o.reached[v], total_cost[v] += o.total_cost[v], harmonic[v] += o.harmonic[v];
    }
  }
};

inline Centrality centrality(const CostAdjacency& g, const std::vector<std::size_t>& sources, const CentralityOptions& o = {}) {
  const std::size_t n = g.size();
  for (const auto& arcs : g)
    for (const CostArc& a : arcs) {
      if (a.to >= n) throw std::out_of_range("centrality: an arc leaves the graph");
      if (!(a.cost >= 0)) throw std::invalid_argument("centrality: an arc's cost is negative");
      if (!(a.length > 0)) throw std::invalid_argument("centrality: an arc's length is not positive");
    }
  const double inf = std::numeric_limits<double>::infinity();
  Centrality C;
  C.reached.assign(n, 0), C.total_cost.assign(n, 0), C.harmonic.assign(n, 0), C.betweenness.assign(n, 0);
  std::vector<double> cost(n, inf), length(n, inf), sigma(n, 0), delta(n, 0);
  std::vector<std::vector<std::size_t>> pred(n);
  std::vector<char> done(n, 0);
  std::vector<std::size_t> order, touched;
  const auto tol = [&](double x) { return o.tie * std::max(1.0, std::fabs(x)); };
  using Item = std::tuple<double, double, std::size_t>;
  for (const std::size_t s : sources) {
    if (s >= n) throw std::out_of_range("centrality: a source is not a vertex");
    for (const std::size_t v : touched) cost[v] = inf, length[v] = inf, sigma[v] = 0, delta[v] = 0, pred[v].clear(), done[v] = 0;
    touched.assign(1, s), order.clear();
    std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
    cost[s] = 0, length[s] = 0, sigma[s] = 1;
    heap.push({0.0, 0.0, s});
    while (!heap.empty()) {
      const auto [d, l, v] = heap.top();
      heap.pop();
      if (done[v]) continue;
      done[v] = 1;
      order.push_back(v);
      for (const CostArc& a : g[v]) {
        const std::size_t w = a.to;
        if (done[w]) continue;
        const double nd = cost[v] + a.cost, nl = length[v] + a.length;
        if (nl > o.radius) continue;
        const bool fresh = cost[w] == inf;
        const bool cheaper = nd < cost[w] - tol(cost[w]);
        const bool same_cost = !cheaper && std::fabs(nd - cost[w]) <= tol(cost[w]);
        if (fresh || cheaper || (same_cost && nl < length[w] - tol(length[w]))) {
          if (fresh) touched.push_back(w);
          cost[w] = nd, length[w] = nl, sigma[w] = sigma[v];
          pred[w].assign(1, v);
          heap.push({nd, nl, w});
        } else if (same_cost && std::fabs(nl - length[w]) <= tol(length[w])) {
          sigma[w] += sigma[v];
          pred[w].push_back(v);
        }
      }
    }
    for (const std::size_t v : order) {
      if (v == s) continue;
      C.reached[s] += 1;
      C.total_cost[s] += cost[v];
      if (cost[v] > 0) C.harmonic[s] += 1 / cost[v];
    }
    if (!o.betweenness) continue;
    for (std::size_t k = order.size(); k-- > 0;) {
      const std::size_t w = order[k];
      for (const std::size_t v : pred[w]) delta[v] += sigma[v] / sigma[w] * (1 + delta[w]);
      if (w != s) C.betweenness[w] += o.undirected ? delta[w] / 2 : delta[w];
    }
  }
  return C;
}

// From every vertex.
inline Centrality centrality(const CostAdjacency& g, const CentralityOptions& o = {}) {
  std::vector<std::size_t> all(g.size());
  for (std::size_t v = 0; v < all.size(); v++) all[v] = v;
  return centrality(g, all, o);
}

// A weighted graph's centralities, its weights both cost and length (they must then be positive).
inline Centrality centrality(const WeightedAdjacency& g, const CentralityOptions& o = {}) {
  CostAdjacency c(g.size());
  for (std::size_t v = 0; v < g.size(); v++)
    for (const WeightedArc& a : g[v]) c[v].push_back({a.to, a.weight, a.weight});
  return centrality(c, o);
}

}  // namespace morphe::combinatorics
