// Checks of morphe/combinatorics: spanning forests (one tree per component, n - c tree
// edges), paths through the forest that are walks between their ends, and distances.

#include <morphe/check.hpp>
#include <morphe/combinatorics/graph.hpp>

#include <random>

namespace mc = morphe::combinatorics;

int main() {
  morphe::Suite s;
  std::mt19937_64 g(5);
  int bad = 0, paths = 0;
  for (int t = 0; t < 200; t++) {
    mc::Multigraph G;
    G.vertices = 1 + g() % 30;
    const std::size_t m = g() % 50;
    for (std::size_t e = 0; e < m; e++) G.edges.push_back({g() % G.vertices, g() % G.vertices});
    const auto F = mc::spanning_forest(G);
    std::size_t tree = 0;
    for (auto x : F.in_tree) tree += x;
    if (tree != G.vertices - F.components()) bad++;
    for (int k = 0; k < 10; k++) {
      const std::size_t u = g() % G.vertices, v = g() % G.vertices;
      if (F.component[u] != F.component[v]) continue;
      std::size_t at = u;
      for (const auto& st : F.path(u, v, G)) {
        const auto [tail, head] = G.edges[st.edge];
        if (!F.in_tree[st.edge] || (st.sign > 0 ? tail : head) != at) { paths++; break; }
        at = st.sign > 0 ? head : tail;
      }
      if (at != v) paths++;
    }
  }
  s.equal("combinatorics: spanning forests have n - c edges, on 200 random multigraphs (failures)", bad, 0);
  s.equal("combinatorics: forest paths are walks along tree edges between their ends (failures)", paths, 0);
  mc::Multigraph C{6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}};
  s.equal("combinatorics: distances around a hexagon", mc::distances(C, 0), std::vector<std::size_t>{0, 1, 2, 3, 2, 1});
  return s.done();
}
