// Shortest paths in a directed graph with non-negative weights, from one source or several, by Dijkstra's algorithm with a binary heap
// (Dijkstra 1959; Cormen et al., Introduction to Algorithms, §24.3). From several sources it is the distance to the nearest of them,
// each with an offset of its own: a virtual source joined to each by an arc of that weight, without the virtual source.
//
// O((V + E) log V). A vertex no source reaches has distance +inf and no parent; `limit` stops the search at that distance, those
// beyond it being left unreached.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

struct WeightedArc {
  std::size_t to = 0;
  double weight = 0;
};
using WeightedAdjacency = std::vector<std::vector<WeightedArc>>;

struct ShortestPaths {
  static constexpr std::size_t none = static_cast<std::size_t>(-1);
  std::vector<double> distance;      // +inf where unreached
  std::vector<std::size_t> parent;   // the arc's tail on a shortest path; `none` at a source or where unreached
  std::vector<std::size_t> origin;   // the source each vertex is nearest; `none` where unreached

  // The vertices of a shortest path to v, from its source.
  std::vector<std::size_t> path_to(std::size_t v) const {
    std::vector<std::size_t> reversed;
    if (distance[v] == std::numeric_limits<double>::infinity()) return reversed;
    for (std::size_t q = v; q != none; q = parent[q]) reversed.push_back(q);
    return {reversed.rbegin(), reversed.rend()};
  }
};

inline ShortestPaths dijkstra(const WeightedAdjacency& graph, const std::vector<std::pair<std::size_t, double>>& sources, double limit = std::numeric_limits<double>::infinity()) {
  const std::size_t n = graph.size();
  ShortestPaths out{std::vector<double>(n, std::numeric_limits<double>::infinity()), std::vector<std::size_t>(n, ShortestPaths::none), std::vector<std::size_t>(n, ShortestPaths::none)};
  using Item = std::pair<double, std::size_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
  for (const auto& [s, offset] : sources) {
    if (s >= n) throw std::out_of_range("dijkstra: a source is not a vertex");
    if (!(offset >= 0)) throw std::invalid_argument("dijkstra: a source's offset is not negative");
    if (offset < out.distance[s] && offset <= limit) {
      out.distance[s] = offset;
      out.origin[s] = s;
      heap.push({offset, s});
    }
  }
  while (!heap.empty()) {
    const auto [d, v] = heap.top();
    heap.pop();
    if (d > out.distance[v]) continue;
    for (const WeightedArc& a : graph[v]) {
      if (!(a.weight >= 0)) throw std::invalid_argument("dijkstra: an arc's weight is not negative");
      const double next = d + a.weight;
      if (next < out.distance[a.to] && next <= limit) {
        out.distance[a.to] = next;
        out.parent[a.to] = v;
        out.origin[a.to] = out.origin[v];
        heap.push({next, a.to});
      }
    }
  }
  return out;
}

inline ShortestPaths dijkstra(const WeightedAdjacency& graph, std::size_t source, double limit = std::numeric_limits<double>::infinity()) {
  return dijkstra(graph, {{source, 0.0}}, limit);
}

// ---- Weights of either sign -----------------------------------------------------------------------------------------------------
//
// With negative weights a shortest path exists exactly when no cycle reachable from the sources has negative weight. Bellman and
// Ford's relaxation (Bellman 1958, On a routing problem, Quarterly of Applied Mathematics 16; Ford 1956; Cormen et al., §24.1)
// settles every distance within V − 1 rounds; an arc that still relaxes in round V lies on, or is reached from, a negative cycle,
// which is found by walking parents V times from it and then around. Floyd and Warshall's recurrence (Floyd 1962, Algorithm 97,
// Communications of the ACM 5) gives every pair at once, a negative diagonal entry marking a vertex on a negative cycle.
//
// These are what the simple temporal networks of evidence need (tekmerion's timelines): a set of bounds on differences of times is
// consistent exactly when its distance graph has no negative cycle, and the cycle, when there is one, is the conflict.

struct NegativeCycle {
  bool found = false;
  std::vector<std::size_t> cycle;  // c0 → c1 → … → ck → c0, each arc a shortest-path arc of the graph
  double weight = 0;               // the cycle's weight, negative when found
};

struct BellmanFord {
  std::vector<double> distance;     // +inf where unreached; meaningless where a negative cycle reaches
  std::vector<std::size_t> parent;  // ShortestPaths::none at a source or where unreached
  NegativeCycle negative;
};

namespace detail {

// The least weight of an arc u → v (a parallel arc of greater weight never lies on a shortest path).
inline double arc_weight(const WeightedAdjacency& g, std::size_t u, std::size_t v) {
  double w = std::numeric_limits<double>::infinity();
  for (const auto& a : g[u])
    if (a.to == v) w = std::min(w, a.weight);
  return w;
}

}  // namespace detail

// Shortest paths from the sources (each with an offset), or the negative cycle that forbids them. With no sources given, every
// vertex is a source at offset 0: the distances are then potentials of the whole graph, and any negative cycle is found.
inline BellmanFord bellman_ford(const WeightedAdjacency& graph, std::vector<std::pair<std::size_t, double>> sources = {}) {
  const std::size_t n = graph.size();
  if (sources.empty())
    for (std::size_t v = 0; v < n; v++) sources.push_back({v, 0.0});
  BellmanFord out{std::vector<double>(n, std::numeric_limits<double>::infinity()), std::vector<std::size_t>(n, ShortestPaths::none), {}};
  for (const auto& [s, offset] : sources) {
    if (s >= n) throw std::out_of_range("bellman_ford: a source is not a vertex");
    if (!std::isfinite(offset)) throw std::invalid_argument("bellman_ford: a source's offset is finite");
    if (offset < out.distance[s]) out.distance[s] = offset;
  }
  for (std::size_t v = 0; v < n; v++)
    for (const auto& a : graph[v]) {
      if (a.to >= n) throw std::out_of_range("bellman_ford: an arc leaves the graph");
      if (!std::isfinite(a.weight)) throw std::invalid_argument("bellman_ford: an arc's weight is finite");
    }
  std::size_t relaxed = ShortestPaths::none;
  for (std::size_t round = 0; round < n; round++) {
    relaxed = ShortestPaths::none;
    for (std::size_t u = 0; u < n; u++) {
      if (out.distance[u] == std::numeric_limits<double>::infinity()) continue;
      for (const auto& a : graph[u])
        if (out.distance[u] + a.weight < out.distance[a.to]) {
          out.distance[a.to] = out.distance[u] + a.weight;
          out.parent[a.to] = u;
          relaxed = a.to;
        }
    }
    if (relaxed == ShortestPaths::none) break;
  }
  if (relaxed != ShortestPaths::none) {
    // Walk back n times to land on the cycle, then around it.
    std::size_t v = relaxed;
    for (std::size_t i = 0; i < n; i++) v = out.parent[v];
    std::vector<std::size_t> back;
    for (std::size_t u = v;; u = out.parent[u]) {
      back.push_back(u);
      if (out.parent[u] == v) break;
    }
    out.negative.found = true;
    out.negative.cycle.assign(back.rbegin(), back.rend());
    for (std::size_t i = 0; i < out.negative.cycle.size(); i++)
      out.negative.weight += detail::arc_weight(graph, out.negative.cycle[i], out.negative.cycle[(i + 1) % out.negative.cycle.size()]);
  }
  return out;
}

struct AllPairs {
  std::vector<std::vector<double>> distance;  // +inf where no path; meaningless when a negative cycle is found
  std::vector<std::vector<std::size_t>> next; // the vertex after i on a shortest path from i to j
  NegativeCycle negative;

  std::vector<std::size_t> path(std::size_t i, std::size_t j) const {
    if (distance[i][j] == std::numeric_limits<double>::infinity()) return {};
    std::vector<std::size_t> p{i};
    while (i != j) p.push_back(i = next[i][j]);
    return p;
  }
};

// Every pair at once, or a negative cycle. O(V³).
inline AllPairs floyd_warshall(const WeightedAdjacency& graph) {
  const std::size_t n = graph.size();
  constexpr double inf = std::numeric_limits<double>::infinity();
  AllPairs out{std::vector<std::vector<double>>(n, std::vector<double>(n, inf)), std::vector<std::vector<std::size_t>>(n, std::vector<std::size_t>(n, ShortestPaths::none)), {}};
  for (std::size_t i = 0; i < n; i++) {
    out.distance[i][i] = 0, out.next[i][i] = i;
    for (const auto& a : graph[i]) {
      if (a.to >= n) throw std::out_of_range("floyd_warshall: an arc leaves the graph");
      if (!std::isfinite(a.weight)) throw std::invalid_argument("floyd_warshall: an arc's weight is finite");
      if (a.weight < out.distance[i][a.to]) out.distance[i][a.to] = a.weight, out.next[i][a.to] = a.to;
    }
  }
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t i = 0; i < n; i++) {
      if (out.distance[i][k] == inf) continue;
      for (std::size_t j = 0; j < n; j++)
        if (out.distance[i][k] + out.distance[k][j] < out.distance[i][j]) {
          out.distance[i][j] = out.distance[i][k] + out.distance[k][j];
          out.next[i][j] = out.next[i][k];
        }
      if (out.distance[i][i] < 0) {
        // i lies on a negative cycle; the graph's own Bellman–Ford names one exactly, its arcs the graph's.
        out.negative = bellman_ford(graph).negative;
        return out;
      }
    }
  return out;
}

}  // namespace morphe::combinatorics
