// Multigraphs: components, spanning forests, paths in a forest, distances.
//
// A Multigraph has vertices 0..n-1 and edges e = (tail, head), loops and parallel edges
// allowed: the 1-skeleton of a CW complex, a Feynman graph. It is the plain C-set over the
// schema E ⇉ V (the idea of 2026-10-21), held as arrays until a second schema needs the
// generic form.
//
// spanning_forest() grows a breadth-first tree in each component, in O(n + m); a path
// between two vertices of one tree is read from it as a sequence of signed edges, so that
// the fundamental group of a 2-complex can be presented on the edges outside the forest
// (Massey 1991, A Basic Course in Algebraic Topology, ch. VI; Hatcher 2002, §1.A).

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <vector>

namespace morphe::combinatorics {

struct Multigraph {
  std::size_t vertices = 0;
  std::vector<std::array<std::size_t, 2>> edges;  // (tail, head)

  // For each vertex, its incident edges with the end at which they meet it: +1 tail, -1 head
  // (a loop appears twice).
  std::vector<std::vector<std::pair<std::size_t, int>>> incidence() const {
    std::vector<std::vector<std::pair<std::size_t, int>>> inc(vertices);
    for (std::size_t e = 0; e < edges.size(); e++) {
      if (edges[e][0] >= vertices || edges[e][1] >= vertices) throw std::out_of_range("morphe::combinatorics::Multigraph: an edge leaves the vertices");
      inc[edges[e][0]].push_back({e, +1});
      inc[edges[e][1]].push_back({e, -1});
    }
    return inc;
  }
};

// A step along an edge: forward (tail to head) or backward.
struct Step {
  std::size_t edge;
  int sign;  // +1 forward, -1 backward
  friend bool operator==(const Step&, const Step&) = default;
};

struct SpanningForest {
  static constexpr std::size_t none = SIZE_MAX;
  std::vector<char> in_tree;               // per edge
  std::vector<std::size_t> parent_edge;    // per vertex; none at a root
  std::vector<std::size_t> parent;         // per vertex; itself at a root
  std::vector<std::size_t> depth;          // per vertex
  std::vector<std::size_t> component;      // per vertex
  std::vector<std::size_t> roots;          // per component
  std::size_t components() const { return roots.size(); }

  // The path from u to v through the forest, as signed edges; u and v in one component.
  std::vector<Step> path(std::size_t u, std::size_t v, const Multigraph& G) const {
    if (component.at(u) != component.at(v)) throw std::invalid_argument("morphe::combinatorics::SpanningForest::path: different components");
    std::vector<Step> up, down;
    auto step_up = [&](std::size_t x) {  // the step from x to its parent
      const std::size_t e = parent_edge[x];
      return Step{e, G.edges[e][0] == x ? +1 : -1};
    };
    while (depth[u] > depth[v]) { up.push_back(step_up(u)); u = parent[u]; }
    while (depth[v] > depth[u]) { down.push_back(step_up(v)); v = parent[v]; }
    while (u != v) {
      up.push_back(step_up(u)); u = parent[u];
      down.push_back(step_up(v)); v = parent[v];
    }
    // down holds steps from v upward; reverse them to go downward toward v.
    for (auto it = down.rbegin(); it != down.rend(); ++it) up.push_back({it->edge, -it->sign});
    return up;
  }
};

// Breadth-first spanning forest, the first tree rooted at `first` (default vertex 0), the
// others at their least vertex.
inline SpanningForest spanning_forest(const Multigraph& G, std::size_t first = 0) {
  SpanningForest F;
  const std::size_t n = G.vertices;
  F.in_tree.assign(G.edges.size(), 0);
  F.parent_edge.assign(n, SpanningForest::none);
  F.parent.resize(n);
  F.depth.assign(n, 0);
  F.component.assign(n, SpanningForest::none);
  const auto inc = G.incidence();
  auto grow = [&](std::size_t root) {
    const std::size_t c = F.roots.size();
    F.roots.push_back(root);
    F.component[root] = c;
    F.parent[root] = root;
    std::queue<std::size_t> q;
    q.push(root);
    while (!q.empty()) {
      const std::size_t x = q.front();
      q.pop();
      for (auto [e, end] : inc[x]) {
        const std::size_t y = G.edges[e][end > 0 ? 1 : 0];
        if (F.component[y] != SpanningForest::none) continue;
        F.component[y] = c;
        F.parent[y] = x;
        F.parent_edge[y] = e;
        F.depth[y] = F.depth[x] + 1;
        F.in_tree[e] = 1;
        q.push(y);
      }
    }
  };
  if (n) grow(std::min(first, n - 1));
  for (std::size_t v = 0; v < n; v++)
    if (F.component[v] == SpanningForest::none) grow(v);
  return F;
}

// Distances in edges from a source (SIZE_MAX where unreachable), ignoring orientation.
inline std::vector<std::size_t> distances(const Multigraph& G, std::size_t source) {
  std::vector<std::size_t> d(G.vertices, SIZE_MAX);
  const auto inc = G.incidence();
  std::queue<std::size_t> q;
  d.at(source) = 0;
  q.push(source);
  while (!q.empty()) {
    const std::size_t x = q.front();
    q.pop();
    for (auto [e, end] : inc[x]) {
      const std::size_t y = G.edges[e][end > 0 ? 1 : 0];
      if (d[y] == SIZE_MAX) d[y] = d[x] + 1, q.push(y);
    }
  }
  return d;
}

}  // namespace morphe::combinatorics
