// Checks of shortest paths (shortest_paths.hpp).
//
// - A grid of unit arcs: every distance is the Manhattan distance, and a path found is as long as its distance says.
// - Several sources: each vertex of a line is as far as the nearest source, with its offset, and is assigned to it.
// - A limit leaves what is beyond it unreached, and what is within it as it was.
// - Random graphs: equal to Bellman-Ford, whose relaxations need no order, on every vertex.
// - Refusals: a negative weight, a source that is no vertex.

#include <morphe/check.hpp>
#include <morphe/combinatorics/shortest_paths.hpp>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace mc = morphe::combinatorics;

int main() {
  morphe::Suite s;
  const double inf = std::numeric_limits<double>::infinity();
  {
    const int n = 12;
    mc::WeightedAdjacency g(n * n);
    const auto id = [&](int i, int j) { return static_cast<std::size_t>(i * n + j); };
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
      if (i + 1 < n) { g[id(i, j)].push_back({id(i + 1, j), 1}); g[id(i + 1, j)].push_back({id(i, j), 1}); }
      if (j + 1 < n) { g[id(i, j)].push_back({id(i, j + 1), 1}); g[id(i, j + 1)].push_back({id(i, j), 1}); }
    }
    const auto P = mc::dijkstra(g, id(3, 7));
    double worst = 0;
    bool paths = true;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
      worst = std::max(worst, std::abs(P.distance[id(i, j)] - (std::abs(i - 3) + std::abs(j - 7))));
      paths = paths && P.path_to(id(i, j)).size() == static_cast<std::size_t>(P.distance[id(i, j)]) + 1 && P.path_to(id(i, j)).front() == id(3, 7);
    }
    s.check("on a grid every distance is the Manhattan distance", worst, 0);
    s.require("and each path found has that many arcs, from the source", paths);
  }
  {
    mc::WeightedAdjacency g(21);
    for (std::size_t i = 0; i + 1 < 21; i++) { g[i].push_back({i + 1, 1}); g[i + 1].push_back({i, 1}); }
    const auto P = mc::dijkstra(g, {{2, 0.0}, {16, 3.0}});
    bool nearest = true;
    for (std::size_t i = 0; i < 21; i++) {
      const double a = std::abs(static_cast<double>(i) - 2), b = 3 + std::abs(static_cast<double>(i) - 16);
      nearest = nearest && P.distance[i] == std::min(a, b) && P.origin[i] == (a <= b ? 2u : 16u);
    }
    s.require("from several sources each vertex is as far as the nearest, with its offset, and is assigned to it", nearest);
    const auto L = mc::dijkstra(g, 10, 4);
    bool cut = true;
    for (std::size_t i = 0; i < 21; i++) {
      const double d = std::abs(static_cast<double>(i) - 10);
      cut = cut && L.distance[i] == (d <= 4 ? d : inf) && (d <= 4 || L.origin[i] == mc::ShortestPaths::none);
    }
    s.require("a limit leaves what is beyond it unreached and what is within it as it was", cut);
  }
  {
    std::mt19937_64 rng(20261107);
    std::uniform_real_distribution<double> w(0, 10);
    double worst = 0;
    int unreached = 0;
    for (int t = 0; t < 40; t++) {
      const std::size_t n = 60;
      mc::WeightedAdjacency g(n);
      for (std::size_t e = 0; e < 180; e++) g[rng() % n].push_back({rng() % n, w(rng)});
      const auto P = mc::dijkstra(g, 0);
      std::vector<double> bf(n, inf);
      bf[0] = 0;
      for (std::size_t it = 0; it < n; it++) for (std::size_t v = 0; v < n; v++) if (bf[v] < inf) for (const auto& a : g[v]) bf[a.to] = std::min(bf[a.to], bf[v] + a.weight);
      for (std::size_t v = 0; v < n; v++) {
        if (bf[v] == inf) { unreached += P.distance[v] == inf ? 1 : 1000; continue; }
        worst = std::max(worst, std::abs(P.distance[v] - bf[v]));
      }
    }
    s.check("on random graphs it equals Bellman-Ford on every vertex", worst, 1e-12);
    s.require("and agrees on what cannot be reached", unreached < 1000 * 1, std::to_string(unreached) + " vertices unreached by both");
  }
  {
    mc::WeightedAdjacency g(2);
    g[0].push_back({1, -1});
    bool refused = false;
    try { mc::dijkstra(g, 0); } catch (const std::invalid_argument&) { refused = true; }
    s.require("a negative weight is refused", refused);
    refused = false;
    try { mc::dijkstra(g, 9); } catch (const std::out_of_range&) { refused = true; }
    s.require("a source that is no vertex is refused", refused);
  }
  return s.done();
}
