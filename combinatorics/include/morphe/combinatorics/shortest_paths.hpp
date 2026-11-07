// Shortest paths in a directed graph with non-negative weights, from one source or several, by Dijkstra's algorithm with a binary heap
// (Dijkstra 1959; Cormen et al., Introduction to Algorithms, §24.3). From several sources it is the distance to the nearest of them,
// each with an offset of its own: a virtual source joined to each by an arc of that weight, without the virtual source.
//
// O((V + E) log V). A vertex no source reaches has distance +inf and no parent; `limit` stops the search at that distance, those
// beyond it being left unreached.

#pragma once

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

}  // namespace morphe::combinatorics
