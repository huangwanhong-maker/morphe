// Checks of flows in networks (flows.hpp).
//
// - The network of Cormen et al.'s Figure 26.1: the greatest flow is 23, its cut has capacity 23, every arc across it saturated.
// - Random networks: the value equals the least capacity over every cut separating s from t, found by listing them; the flow keeps
//   its bounds and conserves at every other vertex.
// - Menger: edge-disjoint paths between opposite corners of a grid are 2, between two vertices of K_n they are n − 1.
// - The cheapest flow: on random assignment problems, the cost equals the least over every permutation; the potentials leave no
//   residual arc with a negative reduced cost.
// - A flow as paths: the amounts sum to the value, and no arc carries more paths than flow.
// - Flows over time: one arc gives (T − τ + 1) u by step T and everyone safe at τ + ⌈A/u⌉ − 1; a later release delays it as much;
//   capacitated sinks bound what is absorbed; on random networks the time-expanded value at every horizon equals Ford and
//   Fulkerson's temporally repeated flow, max over static flows of (T + 1)|f| − Σ τ f, from the cheapest flow at price T + 1.
// - Refusals: the source as the sink, a negative capacity or cost, an arc to itself, a path of unbounded capacity.

#include <morphe/check.hpp>
#include <morphe/combinatorics/flows.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace mc = morphe::combinatorics;

namespace {

bool throws(auto&& f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// Bounds and conservation of a flow from s to t of value v.
bool feasible(std::size_t n, const std::vector<mc::FlowArc>& arcs, const std::vector<double>& f, std::size_t s, std::size_t t, double v) {
  std::vector<double> net(n, 0);
  for (std::size_t a = 0; a < arcs.size(); a++) {
    if (f[a] < -1e-9 || f[a] > arcs[a].capacity + 1e-9) return false;
    net[arcs[a].from] -= f[a];
    net[arcs[a].to] += f[a];
  }
  for (std::size_t u = 0; u < n; u++) {
    const double want = u == s ? -v : u == t ? v : 0;
    if (std::fabs(net[u] - want) > 1e-9) return false;
  }
  return true;
}

}  // namespace

int main() {
  morphe::Suite s;
  const double inf = std::numeric_limits<double>::infinity();
  {
    // s = 0, v1 = 1, v2 = 2, v3 = 3, v4 = 4, t = 5 (Cormen et al., Introduction to Algorithms, 3rd ed., Figure 26.1).
    const std::vector<mc::FlowArc> arcs{{0, 1, 16}, {0, 2, 13}, {2, 1, 4}, {1, 3, 12}, {3, 2, 9}, {2, 4, 14}, {4, 3, 7}, {3, 5, 20}, {4, 5, 4}};
    const auto M = mc::max_flow(6, arcs, 0, 5);
    bool saturated = true;
    for (std::size_t a : M.cut) saturated = saturated && std::fabs(M.flow[a] - arcs[a].capacity) < 1e-12;
    s.require("the network of Cormen et al.'s Figure 26.1 carries 23, its least cut 23 with every arc across it saturated",
            M.value == 23 && M.cut_capacity == 23 && saturated && feasible(6, arcs, M.flow, 0, 5, 23));
  }
  {
    std::mt19937 rng(19);
    std::uniform_int_distribution<int> cap(0, 9), coin(0, 2);
    bool equal = true, ok = true;
    for (int trial = 0; trial < 300; trial++) {
      const std::size_t n = 7;
      std::vector<mc::FlowArc> arcs;
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t v = 0; v < n; v++)
          if (u != v && coin(rng) == 0) arcs.push_back({u, v, static_cast<double>(cap(rng))});
      const auto M = mc::max_flow(n, arcs, 0, n - 1);
      double least = inf;
      for (unsigned mask = 0; mask < (1u << (n - 2)); mask++) {
        std::vector<char> in(n, 0);
        in[0] = 1;
        for (std::size_t i = 1; i + 1 < n; i++) in[i] = (mask >> (i - 1)) & 1u;
        double c = 0;
        for (const auto& a : arcs)
          if (in[a.from] && !in[a.to]) c += a.capacity;
        least = std::min(least, c);
      }
      equal = equal && std::fabs(M.value - least) < 1e-9 && std::fabs(M.cut_capacity - M.value) < 1e-9;
      ok = ok && feasible(n, arcs, M.flow, 0, n - 1, M.value) && M.source_side[0] && !M.source_side[n - 1];
    }
    s.require("on 300 random networks the greatest flow equals the least cut found by listing every cut, and the cut found equals it", equal);
    s.require("each flow keeps 0 ≤ f ≤ u and conserves at every vertex but the source and the sink", ok);
  }
  {
    const std::size_t k = 6;
    std::vector<mc::FlowArc> grid;
    const auto id = [&](std::size_t i, std::size_t j) { return i * k + j; };
    for (std::size_t i = 0; i < k; i++)
      for (std::size_t j = 0; j < k; j++) {
        if (i + 1 < k) { grid.push_back({id(i, j), id(i + 1, j), 1}); grid.push_back({id(i + 1, j), id(i, j), 1}); }
        if (j + 1 < k) { grid.push_back({id(i, j), id(i, j + 1), 1}); grid.push_back({id(i, j + 1), id(i, j), 1}); }
      }
    const auto G = mc::max_flow(k * k, grid, id(0, 0), id(k - 1, k - 1));
    const std::size_t n = 8;
    std::vector<mc::FlowArc> complete;
    for (std::size_t u = 0; u < n; u++)
      for (std::size_t v = 0; v < n; v++)
        if (u != v) complete.push_back({u, v, 1});
    const auto C = mc::max_flow(n, complete, 2, 5);
    s.require("Menger: 2 edge-disjoint paths join opposite corners of a 6 × 6 grid, 7 join two vertices of K_8", G.value == 2 && C.value == 7);
  }
  {
    std::mt19937 rng(7);
    std::uniform_int_distribution<int> cost(0, 20);
    bool equal = true, certified = true;
    for (int trial = 0; trial < 100; trial++) {
      const std::size_t m = 5;  // s = 0, workers 1..m, jobs m+1..2m, t = 2m+1
      std::vector<std::vector<int>> c(m, std::vector<int>(m));
      std::vector<mc::FlowArc> arcs;
      for (std::size_t i = 0; i < m; i++) arcs.push_back({0, 1 + i, 1, 0});
      for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++) {
          c[i][j] = cost(rng);
          arcs.push_back({1 + i, 1 + m + j, 1, static_cast<double>(c[i][j])});
        }
      for (std::size_t j = 0; j < m; j++) arcs.push_back({1 + m + j, 2 * m + 1, 1, 0});
      const auto F = mc::min_cost_flow(2 * m + 2, arcs, 0, 2 * m + 1);
      std::vector<std::size_t> p(m);
      std::iota(p.begin(), p.end(), 0);
      int best = std::numeric_limits<int>::max();
      do {
        int sum = 0;
        for (std::size_t i = 0; i < m; i++) sum += c[i][p[i]];
        best = std::min(best, sum);
      } while (std::next_permutation(p.begin(), p.end()));
      equal = equal && F.value == m && std::fabs(F.cost - best) < 1e-9;
      for (std::size_t a = 0; a < arcs.size(); a++) {
        const double rc = arcs[a].cost + F.potential[arcs[a].from] - F.potential[arcs[a].to];
        if (F.flow[a] < arcs[a].capacity - 1e-9 && rc < -1e-9) certified = false;
        if (F.flow[a] > 1e-9 && -rc < -1e-9) certified = false;
      }
    }
    s.require("on 100 random assignments of 5 the cheapest flow costs the least over all 120 permutations", equal);
    s.require("its potentials leave no residual arc with a negative reduced cost", certified);
  }
  {
    std::mt19937 rng(23);
    std::uniform_int_distribution<int> cap(0, 9), coin(0, 2);
    bool ok = true;
    for (int trial = 0; trial < 100; trial++) {
      const std::size_t n = 7;
      std::vector<mc::FlowArc> arcs;
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t v = 0; v < n; v++)
          if (u != v && coin(rng) == 0) arcs.push_back({u, v, static_cast<double>(cap(rng))});
      const auto M = mc::max_flow(n, arcs, 0, n - 1);
      const auto P = mc::decompose(n, arcs, M.flow, 0, n - 1);
      double total = 0;
      std::vector<double> through(arcs.size(), 0);
      for (const auto& p : P) {
        total += p.amount;
        for (std::size_t a : p.arcs) through[a] += p.amount;
        ok = ok && !p.arcs.empty() && arcs[p.arcs.front()].from == 0 && arcs[p.arcs.back()].to == n - 1;
      }
      for (std::size_t a = 0; a < arcs.size(); a++) ok = ok && through[a] <= M.flow[a] + 1e-9;
      ok = ok && std::fabs(total - M.value) < 1e-9;
    }
    s.require("a flow as paths from source to sink: their amounts sum to the value, no arc carrying more than its flow", ok);
  }
  {
    // One arc of transit 3 and capacity 4 a step, 50 waiting at its tail.
    const std::vector<mc::TimedArc> arcs{{0, 1, 3, 4}};
    const auto O = mc::flow_over_time(2, arcs, {{0, 50, 0}}, {{1}}, 20);
    bool line = true;
    for (std::size_t T = 0; T <= 20; T++) line = line && std::fabs(O.safe[T] - std::min(50.0, T >= 3 ? (T - 3 + 1) * 4.0 : 0.0)) < 1e-9;
    s.require("one arc: (T − τ + 1) u reach the sink by step T, everyone by τ + ⌈A/u⌉ − 1 = 15", line && O.quickest == 15 && O.total == 50);
    const auto D = mc::flow_over_time(2, arcs, {{0, 50, 4}}, {{1}}, 25);
    s.require("released at step 4, everyone is safe 4 steps later", D.quickest == 19);
    const std::vector<mc::TimedArc> two{{0, 1, 1, 10}, {0, 2, 2, 10}};
    const auto K = mc::flow_over_time(3, two, {{0, 100, 0}}, {{1, 30}, {2, 25}}, 30);
    s.require("capacitated sinks absorb no more than they hold: 30 and 25 of 100", std::fabs(K.absorbed[0] - 30) < 1e-9 && std::fabs(K.absorbed[1] - 25) < 1e-9 &&
                                                                                 std::fabs(K.safe.back() - 55) < 1e-9 && K.quickest == mc::OverTime::none);
  }
  {
    std::mt19937 rng(41);
    std::uniform_int_distribution<int> cap(1, 5), tau(1, 4), coin(0, 2);
    bool equal = true, monotone = true;
    for (int trial = 0; trial < 60; trial++) {
      const std::size_t n = 6;
      std::vector<mc::TimedArc> timed;
      std::vector<mc::FlowArc> stat;
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t v = 0; v < n; v++)
          if (u != v && coin(rng) == 0) {
            const int c = cap(rng), d = tau(rng);
            timed.push_back({u, v, static_cast<std::size_t>(d), static_cast<double>(c)});
            stat.push_back({u, v, static_cast<double>(c), static_cast<double>(d)});
          }
      const std::size_t H = 14;
      const auto O = mc::flow_over_time(n, timed, {{0, 1e6, 0}}, {{n - 1}}, H);
      for (std::size_t T = 0; T <= H; T++) {
        const auto F = mc::min_cost_flow(n, stat, 0, n - 1, inf, static_cast<double>(T + 1));
        double repeated = 0;
        for (std::size_t i = 0; i < F.path_costs.size(); i++) repeated += (T + 1 - F.path_costs[i]) * F.path_amounts[i];
        equal = equal && std::fabs(O.safe[T] - repeated) < 1e-6;
        if (T > 0) monotone = monotone && O.safe[T] >= O.safe[T - 1] - 1e-9;
      }
    }
    s.require("on 60 random networks the time-expanded value at every horizon T ≤ 14 equals the temporally repeated flow, max (T + 1)|f| − Σ τ f", equal);
    s.require("what can be safe by step t never falls as t grows", monotone);
  }
  {
    s.require("refusals: the source as the sink, a negative capacity, a negative cost, an arc to itself, a path of unbounded capacity",
            throws([] { mc::max_flow(2, {{0, 1, 1}}, 0, 0); }) && throws([] { mc::max_flow(2, {{0, 1, -1}}, 0, 1); }) &&
                throws([] { mc::min_cost_flow(2, {{0, 1, 1, -2}}, 0, 1); }) && throws([] { mc::max_flow(2, {{1, 1, 1}}, 0, 1); }) &&
                throws([&] { mc::max_flow(3, {{0, 1, inf}, {1, 2, inf}}, 0, 2); }));
  }
  return s.done();
}
