// Checks of centralities (centrality.hpp).
//
// - Closed forms: a path's betweenness i(n − 1 − i) and its total distances; a star's centre on every pair of leaves, C(k, 2).
// - networkx 3.4.2 (tests/fixtures/centrality.txt, written by make_centrality_fixtures.py): betweenness, total distance, harmonic
//   closeness and vertices reached on a random graph with real weights, a grid with unit weights (many ties) and a geometric graph
//   with integer weights.
// - Brute force, every simple path enumerated, on twenty random graphs whose costs include zeros (cycles of zero cost among them) and
//   whose lengths break ties: least paths by (cost, length), their counts and the betweenness they give, exactly.
// - A radius, cost being length: the vertices within it and the betweenness over pairs within it, against brute force.
// - Sources split in two and added give the whole.
// - Refusals: a negative cost, a length that is not positive.

#include <morphe/check.hpp>
#include <morphe/combinatorics/centrality.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef MORPHE_FIXTURES
#define MORPHE_FIXTURES "tests/fixtures"
#endif

namespace mc = morphe::combinatorics;

namespace {

void join(mc::CostAdjacency& g, std::size_t u, std::size_t v, double cost, double length) {
  g[u].push_back({v, cost, length});
  g[v].push_back({u, cost, length});
}

// Every simple path from each source: the least by (cost, length), exactly (integer weights), their number and, for each vertex, the
// number through it; betweenness over ordered pairs within `radius` of length, halved, and the total least cost.
struct Brute {
  std::vector<double> betweenness, total_cost, reached;
};
Brute brute(const mc::CostAdjacency& g, double radius = std::numeric_limits<double>::infinity()) {
  const std::size_t n = g.size();
  Brute B{std::vector<double>(n, 0), std::vector<double>(n, 0), std::vector<double>(n, 0)};
  for (std::size_t s = 0; s < n; s++) {
    // For each target, the best key and the paths attaining it.
    std::vector<std::pair<double, double>> best(n, {std::numeric_limits<double>::infinity(), 0});
    std::vector<std::vector<std::vector<std::size_t>>> paths(n);
    std::vector<std::size_t> path{s};
    std::vector<char> on(n, 0);
    on[s] = 1;
    const std::function<void(std::size_t, double, double)> walk = [&](std::size_t v, double c, double l) {
      if (v != s) {
        const std::pair<double, double> key{c, l};
        if (key < best[v]) best[v] = key, paths[v].assign(1, path);
        else if (key == best[v]) paths[v].push_back(path);
      }
      for (const auto& a : g[v]) {
        if (on[a.to]) continue;
        on[a.to] = 1, path.push_back(a.to);
        walk(a.to, c + a.cost, l + a.length);
        path.pop_back(), on[a.to] = 0;
      }
    };
    walk(s, 0, 0);
    for (std::size_t t = 0; t < n; t++) {
      if (t == s || paths[t].empty() || best[t].second > radius) continue;
      B.reached[s] += 1;
      B.total_cost[s] += best[t].first;
      const double sigma = static_cast<double>(paths[t].size());
      for (const auto& p : paths[t])
        for (std::size_t k = 1; k + 1 < p.size(); k++) B.betweenness[p[k]] += 0.5 / sigma;
    }
  }
  return B;
}

double worst_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double w = 0;
  for (std::size_t i = 0; i < a.size(); i++) w = std::max(w, std::fabs(a[i] - b[i]));
  return w;
}

}  // namespace

int main() {
  morphe::Suite s;
  {
    mc::CostAdjacency g(7);
    for (std::size_t v = 0; v + 1 < 7; v++) join(g, v, v + 1, 1, 1);
    mc::CentralityOptions o;
    o.undirected = true;
    const auto C = mc::centrality(g, o);
    double worst = 0;
    for (std::size_t i = 0; i < 7; i++) {
      worst = std::max(worst, std::fabs(C.betweenness[i] - static_cast<double>(i * (6 - i))));
      const double total = static_cast<double>(i * (i + 1) / 2 + (6 - i) * (7 - i) / 2);
      worst = std::max(worst, std::fabs(C.total_cost[i] - total));
    }
    s.check("a path of 7: betweenness i(n − 1 − i) and total distances", worst, 1e-12);
    mc::CostAdjacency star(7);
    for (std::size_t v = 1; v < 7; v++) join(star, 0, v, 1, 1);
    const auto S = mc::centrality(star, o);
    s.check("a star of 6 leaves: its centre lies on every pair of leaves, C(6, 2) = 15", S.betweenness[0] - 15, 1e-12);
  }
  {
    std::ifstream in(std::string(MORPHE_FIXTURES) + "/centrality.txt");
    s.require("networkx: the fixture is read", static_cast<bool>(in));
    std::string line;
    while (std::getline(in, line)) {
      if (line.rfind("graph ", 0) != 0) continue;
      std::istringstream head(line.substr(6));
      std::string name;
      std::size_t n = 0, m = 0;
      head >> name >> n >> m;
      mc::CostAdjacency g(n);
      for (std::size_t e = 0; e < m; e++) {
        std::size_t u = 0, v = 0;
        double w = 0;
        in >> u >> v >> w;
        join(g, u, v, w, w);
      }
      std::vector<double> btw(n), total(n), harmonic(n), reached(n);
      for (std::size_t k = 0; k < n; k++) {
        std::size_t v = 0;
        in >> v >> btw[v] >> total[v] >> harmonic[v] >> reached[v];
      }
      mc::CentralityOptions o;
      o.undirected = true;
      const auto C = mc::centrality(g, o);
      s.check("networkx, " + name + ": betweenness (relative to the largest)",
              worst_difference(C.betweenness, btw) / std::max(1.0, *std::max_element(btw.begin(), btw.end())), 1e-9);
      s.check("networkx, " + name + ": total distance (relative)",
              worst_difference(C.total_cost, total) / std::max(1.0, *std::max_element(total.begin(), total.end())), 1e-12);
      s.check("networkx, " + name + ": harmonic closeness", worst_difference(C.harmonic, harmonic), 1e-12);
      s.check("networkx, " + name + ": vertices reached", worst_difference(C.reached, reached), 0);
    }
  }
  {
    std::mt19937_64 rng(42);
    int agree = 0, zero_cycles = 0;
    for (int r = 0; r < 20; r++) {
      const std::size_t n = 7;
      mc::CostAdjacency g(n);
      std::map<std::pair<std::size_t, std::size_t>, int> seen;
      for (std::size_t v = 1; v < n; v++) {  // a spanning tree, then chords
        const std::size_t u = rng() % v;
        seen[{u, v}] = 1;
        join(g, u, v, static_cast<double>(rng() % 3), static_cast<double>(1 + rng() % 3));
      }
      for (int k = 0; k < 6; k++) {
        std::size_t u = rng() % n, v = rng() % n;
        if (u == v) continue;
        if (u > v) std::swap(u, v);
        if (seen.count({u, v})) continue;
        seen[{u, v}] = 1;
        join(g, u, v, static_cast<double>(rng() % 3), static_cast<double>(1 + rng() % 3));
      }
      for (std::size_t v = 0; v < n; v++)
        for (const auto& a : g[v]) zero_cycles += a.cost == 0;
      mc::CentralityOptions o;
      o.undirected = true;
      const auto C = mc::centrality(g, o);
      const auto B = brute(g);
      agree += worst_difference(C.betweenness, B.betweenness) < 1e-9 && worst_difference(C.total_cost, B.total_cost) < 1e-9 &&
               worst_difference(C.reached, B.reached) == 0;
    }
    s.require("brute force, twenty graphs with zero costs: least paths by (cost, length), their betweenness and total cost", agree == 20,
              std::to_string(agree) + " of 20; " + std::to_string(zero_cycles / 2) + " arcs of zero cost");
  }
  {
    std::mt19937_64 rng(7);
    int agree = 0;
    for (int r = 0; r < 10; r++) {
      const std::size_t n = 8;
      mc::CostAdjacency g(n);
      for (std::size_t v = 1; v < n; v++) {
        const double w = static_cast<double>(1 + rng() % 4);
        join(g, rng() % v, v, w, w);
      }
      for (int k = 0; k < 5; k++) {
        const std::size_t u = rng() % n, v = rng() % n;
        if (u == v) continue;
        const double w = static_cast<double>(1 + rng() % 4);
        join(g, u, v, w, w);
      }
      mc::CentralityOptions o;
      o.undirected = true;
      o.radius = 4;
      const auto C = mc::centrality(g, o);
      const auto B = brute(g, 4);
      agree += worst_difference(C.betweenness, B.betweenness) < 1e-9 && worst_difference(C.reached, B.reached) == 0 &&
               worst_difference(C.total_cost, B.total_cost) < 1e-9;
    }
    s.require("a radius of 4, cost being length: the vertices within it and the betweenness over pairs within it (ten graphs)", agree == 10,
              std::to_string(agree) + " of 10");
  }
  {
    mc::CostAdjacency g(36);
    for (std::size_t i = 0; i < 6; i++)
      for (std::size_t j = 0; j < 6; j++) {
        if (i + 1 < 6) join(g, i * 6 + j, (i + 1) * 6 + j, 1 + 0.1 * static_cast<double>(j), 1);
        if (j + 1 < 6) join(g, i * 6 + j, i * 6 + j + 1, 0.5, 1);
      }
    std::vector<std::size_t> a, b;
    for (std::size_t v = 0; v < 36; v++) (v % 3 ? a : b).push_back(v);
    auto C = mc::centrality(g, a);
    C.add(mc::centrality(g, b));
    const auto W = mc::centrality(g);
    s.check("sources split in two and added give the whole", worst_difference(C.betweenness, W.betweenness) + worst_difference(C.total_cost, W.total_cost), 1e-9);
  }
  {
    bool refused = false;
    try {
      mc::CostAdjacency g(2);
      join(g, 0, 1, -1, 1);
      mc::centrality(g);
    } catch (const std::invalid_argument&) { refused = true; }
    bool refused_length = false;
    try {
      mc::CostAdjacency g(2);
      join(g, 0, 1, 1, 0);
      mc::centrality(g);
    } catch (const std::invalid_argument&) { refused_length = true; }
    s.require("refusals: a negative cost, a length that is not positive", refused && refused_length);
  }
  return s.done();
}
