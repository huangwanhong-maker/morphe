// Flows in networks: the greatest flow and its least cut, the cheapest flow, the decomposition of a flow into paths, and flows over
// time (decision 0019).
//
// A network is a directed graph whose arcs carry a capacity u ≥ 0 (+inf allowed) and, for the cheapest flow, a cost c ≥ 0 a unit. A
// flow f from s to t keeps 0 ≤ f ≤ u on every arc and conserves at every vertex but s and t; its value is what leaves s.
//
// THE GREATEST FLOW (Ford and Fulkerson 1956, Canadian Journal of Mathematics 8, 399; Dinic 1970, Soviet Math. Doklady 11, 1277).
// Dinic's algorithm: levels by breadth-first search in the residual network, then a blocking flow along arcs that rise one level,
// O(V² E). Its certificate is the cut: the vertices still reached from s in the residual network form a set S with t outside, every
// arc from S to its complement saturated and every arc back empty, so the value equals the capacity of the cut (the max-flow min-cut
// theorem), and no flow can be greater. Several sources and sinks are joined to one by arcs of their supplies and capacities.
// Edge-disjoint paths are a flow with unit capacities (Menger): an undirected edge is two arcs of capacity one.
//
// THE CHEAPEST FLOW. Successive shortest paths with potentials (Edmonds and Karp 1972, J. ACM 19, 248; Ahuja, Magnanti and Orlin,
// Network Flows, 1993, §9.7): the cheapest path in the residual network under the reduced costs c + π(u) − π(v) ≥ 0, found by
// Dijkstra's algorithm, is augmented and π moves by the distances. The flow after each augmentation is the cheapest of its value, and
// the potentials are its certificate: no residual arc has a negative reduced cost. With a price p, a unit is shipped only while the
// cheapest path costs less than p, which maximises p·|f| − Σ c f: the form of Ford and Fulkerson's temporally repeated flows.
//
// FLOWS OVER TIME (Ford and Fulkerson 1958, Operations Research 6, 419; Skutella 2009, An introduction to network flows over time).
// Each arc has a transit time τ (whole steps) and a capacity per step; what is at a vertex may wait there (up to a storage, +inf by
// default). The time-expanded network has a copy of every vertex for each step 0..H, an arc (u, t) → (v, t + τ) for each arc and step,
// and an arc (v, t) → (v, t + 1) for waiting; supplies enter at their vertex and step of release, and sinks absorb up to their
// capacities over the whole horizon. Its greatest flow with the sinks open up to step t is the most that can reach a sink by t, and
// `flow_over_time` gives it for every t on one network, opening the sinks step by step and augmenting what was found. Each value is
// the greatest for its own t; one flow need not attain them all at once (an earliest-arrival flow exists for a single sink, Gale
// 1959, and may not with several capacitated sinks).
//
// Caveat. Time is in whole steps. A continuous problem is bounded from above by the discrete one whose transit times are rounded down,
// whose capacities per step are the rates times the step and whose horizon is rounded up: every continuous flow can be carried by it,
// waiting at a vertex for what rounding gained. That is how a caller makes the result an upper bound; a step that is too coarse makes
// the bound loose, never wrong.
//
// Capacities are doubles: data of the network, compared with a tolerance relative to the largest finite capacity.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

struct FlowArc {
  std::size_t from = 0, to = 0;
  double capacity = 0;  // ≥ 0, +inf allowed
  double cost = 0;      // ≥ 0 a unit; read only by the cheapest flow
};

namespace flow_detail {

inline constexpr double inf = std::numeric_limits<double>::infinity();

struct R {
  std::size_t to = 0, rev = 0;
  double cap = 0, cost = 0;
};

// The residual network: each arc a forward residual arc and its reverse; the flow on an arc is its reverse's residual capacity.
struct Residual {
  std::vector<std::vector<R>> g;
  std::vector<std::pair<std::size_t, std::size_t>> where;  // input arc → (its tail, the index of its forward arc there)
  double eps = 0;

  Residual(std::size_t n, const std::vector<FlowArc>& arcs, double tolerance, const char* who) : g(n) {
    double biggest = 0;
    where.reserve(arcs.size());
    for (const FlowArc& a : arcs) {
      if (a.from >= n || a.to >= n) throw std::invalid_argument(std::string(who) + ": an arc's ends are vertices");
      if (a.from == a.to) throw std::invalid_argument(std::string(who) + ": an arc joins two different vertices");
      if (!(a.capacity >= 0)) throw std::invalid_argument(std::string(who) + ": capacities are ≥ 0");
      if (!std::isfinite(a.cost)) throw std::invalid_argument(std::string(who) + ": costs are finite");
      if (std::isfinite(a.capacity)) biggest = std::max(biggest, a.capacity);
      add(a.from, a.to, a.capacity, a.cost);
    }
    eps = tolerance * std::max(biggest, 1.0);
  }
  std::size_t add(std::size_t u, std::size_t v, double cap, double cost) {
    g[u].push_back({v, g[v].size(), cap, cost});
    g[v].push_back({u, g[u].size() - 1, 0, -cost});
    where.push_back({u, g[u].size() - 1});
    return where.size() - 1;
  }
  R& forward(std::size_t arc) { return g[where[arc].first][where[arc].second]; }
  const R& forward(std::size_t arc) const { return g[where[arc].first][where[arc].second]; }
  double flow(std::size_t arc) const {
    const R& r = forward(arc);
    return g[r.to][r.rev].cap;
  }
  void push(std::size_t u, std::size_t k, double d) {
    R& r = g[u][k];
    r.cap -= d;
    g[r.to][r.rev].cap += d;
  }

  // Dinic's algorithm from the present flow: what more can go from s to t. Iterative, so that long paths of a time-expanded network
  // need no deep recursion.
  double augment(std::size_t s, std::size_t t, const char* who) {
    const std::size_t n = g.size();
    std::vector<long long> level(n);
    std::vector<std::size_t> it(n);
    std::vector<std::pair<std::size_t, std::size_t>> path;
    double total = 0;
    for (;;) {
      std::fill(level.begin(), level.end(), -1);
      std::queue<std::size_t> q;
      level[s] = 0;
      q.push(s);
      while (!q.empty()) {
        const std::size_t u = q.front();
        q.pop();
        for (const R& r : g[u])
          if (r.cap > eps && level[r.to] < 0) {
            level[r.to] = level[u] + 1;
            q.push(r.to);
          }
      }
      if (level[t] < 0) return total;
      std::fill(it.begin(), it.end(), 0);
      for (;;) {
        path.clear();
        std::size_t u = s;
        bool found = false;
        for (;;) {
          if (u == t) {
            found = true;
            break;
          }
          auto& k = it[u];
          while (k < g[u].size() && !(g[u][k].cap > eps && level[g[u][k].to] == level[u] + 1)) k++;
          if (k == g[u].size()) {
            level[u] = -1;  // a dead end at this phase
            if (path.empty()) break;
            u = path.back().first;
            path.pop_back();
            it[u]++;
            continue;
          }
          path.push_back({u, k});
          u = g[u][k].to;
        }
        if (!found) break;
        double d = inf;
        for (const auto& [v, k] : path) d = std::min(d, g[v][k].cap);
        if (!std::isfinite(d)) throw std::invalid_argument(std::string(who) + ": a path of unbounded capacity joins the source to the sink");
        for (const auto& [v, k] : path) push(v, k, d);
        total += d;
      }
    }
  }

  // The vertices reached from s by residual arcs.
  std::vector<char> reached(std::size_t s) const {
    std::vector<char> seen(g.size(), 0);
    std::vector<std::size_t> stack{s};
    seen[s] = 1;
    while (!stack.empty()) {
      const std::size_t u = stack.back();
      stack.pop_back();
      for (const R& r : g[u])
        if (r.cap > eps && !seen[r.to]) {
          seen[r.to] = 1;
          stack.push_back(r.to);
        }
    }
    return seen;
  }
};

inline void ends(std::size_t n, std::size_t s, std::size_t t, const char* who) {
  if (s >= n || t >= n) throw std::invalid_argument(std::string(who) + ": the source and the sink are vertices");
  if (s == t) throw std::invalid_argument(std::string(who) + ": the source is not the sink");
}

}  // namespace flow_detail

// ---- the greatest flow ------------------------------------------------------------------------------------------------------------

struct MaxFlow {
  double value = 0;
  std::vector<double> flow;        // on each arc, in the order given
  std::vector<char> source_side;   // the least cut: 1 for the vertices reached from s in the residual network
  std::vector<std::size_t> cut;    // the arcs from the source side to the other, each saturated
  double cut_capacity = 0;         // Σ of their capacities: equal to the value
};

inline MaxFlow max_flow(std::size_t n, const std::vector<FlowArc>& arcs, std::size_t s, std::size_t t, double tolerance = 1e-12) {
  const char* who = "morphe::combinatorics::max_flow";
  flow_detail::ends(n, s, t, who);
  flow_detail::Residual R(n, arcs, tolerance, who);
  MaxFlow out;
  out.value = R.augment(s, t, who);
  out.flow.resize(arcs.size());
  for (std::size_t a = 0; a < arcs.size(); a++) out.flow[a] = R.flow(a);
  out.source_side = R.reached(s);
  for (std::size_t a = 0; a < arcs.size(); a++)
    if (out.source_side[arcs[a].from] && !out.source_side[arcs[a].to]) {
      out.cut.push_back(a);
      out.cut_capacity += arcs[a].capacity;
    }
  return out;
}

// ---- the cheapest flow ------------------------------------------------------------------------------------------------------------

struct MinCostFlow {
  double value = 0;               // what was shipped from s to t
  double cost = 0;                // Σ c f
  std::vector<double> flow;       // on each arc, in the order given
  std::vector<double> potential;  // π: every residual arc has c + π(u) − π(v) ≥ 0 (to rounding); π(t) − π(s) the last path's cost
  std::vector<double> path_costs; // the cost a unit of each augmentation, in order (non-decreasing)
  std::vector<double> path_amounts;
};

// Ships up to `amount` from s to t at least cost, each unit only while the cheapest path costs less than `price` a unit (+inf: the
// cheapest of the greatest flow). Costs are ≥ 0.
inline MinCostFlow min_cost_flow(std::size_t n, const std::vector<FlowArc>& arcs, std::size_t s, std::size_t t,
                                 double amount = flow_detail::inf, double price = flow_detail::inf, double tolerance = 1e-12) {
  const char* who = "morphe::combinatorics::min_cost_flow";
  flow_detail::ends(n, s, t, who);
  for (const FlowArc& a : arcs)
    if (!(a.cost >= 0)) throw std::invalid_argument(std::string(who) + ": costs are ≥ 0");
  if (!(amount >= 0)) throw std::invalid_argument(std::string(who) + ": an amount ≥ 0");
  flow_detail::Residual R(n, arcs, tolerance, who);
  MinCostFlow out;
  out.potential.assign(n, 0.0);
  std::vector<double> dist(n);
  std::vector<std::pair<std::size_t, std::size_t>> parent(n);
  using Item = std::pair<double, std::size_t>;
  while (out.value < amount) {
    std::fill(dist.begin(), dist.end(), flow_detail::inf);
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
    dist[s] = 0;
    heap.push({0, s});
    while (!heap.empty()) {
      const auto [d, u] = heap.top();
      heap.pop();
      if (d > dist[u]) continue;
      for (std::size_t k = 0; k < R.g[u].size(); k++) {
        const auto& r = R.g[u][k];
        if (!(r.cap > R.eps)) continue;
        const double nd = d + std::max(0.0, r.cost + out.potential[u] - out.potential[r.to]);
        if (nd < dist[r.to]) {
          dist[r.to] = nd;
          parent[r.to] = {u, k};
          heap.push({nd, r.to});
        }
      }
    }
    if (dist[t] == flow_detail::inf) break;
    for (std::size_t v = 0; v < n; v++)
      if (dist[v] < flow_detail::inf) out.potential[v] += dist[v];
    const double unit = out.potential[t] - out.potential[s];
    if (!(unit < price)) break;
    double d = amount - out.value;
    for (std::size_t v = t; v != s; v = parent[v].first) d = std::min(d, R.g[parent[v].first][parent[v].second].cap);
    if (!std::isfinite(d)) throw std::invalid_argument(std::string(who) + ": a path of unbounded capacity joins the source to the sink");
    for (std::size_t v = t; v != s; v = parent[v].first) R.push(parent[v].first, parent[v].second, d);
    out.value += d;
    out.cost += d * unit;
    out.path_costs.push_back(unit);
    out.path_amounts.push_back(d);
  }
  out.flow.resize(arcs.size());
  for (std::size_t a = 0; a < arcs.size(); a++) out.flow[a] = R.flow(a);
  return out;
}

// ---- a flow as paths --------------------------------------------------------------------------------------------------------------

struct FlowPath {
  std::vector<std::size_t> arcs;  // from s to t
  double amount = 0;
};

// The flow from s to t as simple paths with their amounts (what circulates in cycles is left out), each path taken by a depth-first
// search in the arcs still carrying flow. Σ amounts = the value; on every arc Σ of the paths through it ≤ its flow.
inline std::vector<FlowPath> decompose(std::size_t n, const std::vector<FlowArc>& arcs, std::vector<double> flow, std::size_t s, std::size_t t,
                                       double tolerance = 1e-9) {
  const char* who = "morphe::combinatorics::decompose";
  flow_detail::ends(n, s, t, who);
  if (flow.size() != arcs.size()) throw std::invalid_argument(std::string(who) + ": one flow for each arc");
  double biggest = 0;
  for (double f : flow) biggest = std::max(biggest, std::fabs(f));
  const double eps = tolerance * std::max(1.0, biggest);
  std::vector<std::vector<std::size_t>> out_arcs(n);
  for (std::size_t a = 0; a < arcs.size(); a++) out_arcs[arcs[a].from].push_back(a);
  std::vector<FlowPath> paths;
  std::vector<char> on(n);
  std::vector<std::size_t> it(n);
  for (;;) {
    std::fill(on.begin(), on.end(), 0);
    std::fill(it.begin(), it.end(), 0);
    std::vector<std::size_t> stack;  // arcs
    std::size_t u = s;
    on[s] = 1;
    bool found = false;
    std::vector<char> dead(n, 0);
    for (;;) {
      if (u == t) {
        found = true;
        break;
      }
      auto& k = it[u];
      while (k < out_arcs[u].size()) {
        const std::size_t a = out_arcs[u][k];
        if (flow[a] > eps && !on[arcs[a].to] && !dead[arcs[a].to]) break;
        k++;
      }
      if (k == out_arcs[u].size()) {
        dead[u] = 1;
        on[u] = 0;
        if (stack.empty()) break;
        u = arcs[stack.back()].from;
        stack.pop_back();
        it[u]++;
        continue;
      }
      const std::size_t a = out_arcs[u][k];
      stack.push_back(a);
      u = arcs[a].to;
      on[u] = 1;
    }
    if (!found) break;
    double d = flow_detail::inf;
    for (std::size_t a : stack) d = std::min(d, flow[a]);
    for (std::size_t a : stack) flow[a] -= d;
    paths.push_back({stack, d});
  }
  return paths;
}

// ---- flows over time --------------------------------------------------------------------------------------------------------------

struct TimedArc {
  std::size_t from = 0, to = 0;
  std::size_t transit = 1;  // whole steps
  double capacity = 0;      // what may enter it in one step
};
struct TimedSupply {
  std::size_t vertex = 0;
  double amount = 0;
  std::size_t release = 0;  // the step from which it may move
};
struct TimedSink {
  std::size_t vertex = 0;
  double capacity = flow_detail::inf;  // what it absorbs over the whole horizon
};

struct OverTime {
  std::size_t horizon = 0;
  double total = 0;                       // Σ supplies
  std::vector<double> safe;               // safe[t]: the most that can reach a sink by step t, each t on its own (t = 0..horizon)
  static constexpr std::size_t none = static_cast<std::size_t>(-1);
  std::size_t quickest = none;            // the least t with safe[t] = total (to tolerance), none within the horizon
  // The greatest flow at the horizon (one flow, the last found):
  std::vector<double> absorbed;           // by each sink
  std::vector<std::vector<double>> flow;  // by each arc, by step of entry 0..horizon
  std::vector<std::vector<double>> waiting;  // by vertex, by step: what waits there from t to t + 1
};

// The greatest flow over time on the time-expanded network of `horizon` steps, with the most that can reach the sinks by every step.
inline OverTime flow_over_time(std::size_t n, const std::vector<TimedArc>& arcs, const std::vector<TimedSupply>& supplies,
                               const std::vector<TimedSink>& sinks, std::size_t horizon, double storage = flow_detail::inf,
                               double tolerance = 1e-12) {
  const char* who = "morphe::combinatorics::flow_over_time";
  if (n == 0) throw std::invalid_argument(std::string(who) + ": a network with vertices");
  for (const TimedArc& a : arcs) {
    if (a.from >= n || a.to >= n || a.from == a.to) throw std::invalid_argument(std::string(who) + ": an arc joins two different vertices");
    if (!(a.capacity >= 0)) throw std::invalid_argument(std::string(who) + ": capacities are ≥ 0");
  }
  for (const TimedSupply& p : supplies)
    if (p.vertex >= n || !(p.amount >= 0 && std::isfinite(p.amount))) throw std::invalid_argument(std::string(who) + ": a supply is finite, ≥ 0, at a vertex");
  if (sinks.empty()) throw std::invalid_argument(std::string(who) + ": at least one sink");
  for (const TimedSink& k : sinks)
    if (k.vertex >= n || !(k.capacity >= 0)) throw std::invalid_argument(std::string(who) + ": a sink is at a vertex, with a capacity ≥ 0");
  if (!(storage >= 0)) throw std::invalid_argument(std::string(who) + ": a storage ≥ 0");
  const std::size_t H = horizon, L = H + 1;
  const std::size_t S = L * n, K = sinks.size(), T = S + 1 + K;
  const auto at = [n](std::size_t v, std::size_t step) { return step * n + v; };
  flow_detail::Residual R(T + 1, {}, tolerance, who);
  double biggest = 0;
  OverTime out;
  out.horizon = H;
  // Arcs of the expanded network, kept by their place so that flows can be read back.
  std::vector<std::vector<std::size_t>> moving(arcs.size(), std::vector<std::size_t>(L, static_cast<std::size_t>(-1)));
  for (std::size_t a = 0; a < arcs.size(); a++) {
    if (std::isfinite(arcs[a].capacity)) biggest = std::max(biggest, arcs[a].capacity);
    for (std::size_t step = 0; step + arcs[a].transit <= H; step++)
      moving[a][step] = R.add(at(arcs[a].from, step), at(arcs[a].to, step + arcs[a].transit), arcs[a].capacity, 0);
  }
  std::vector<std::vector<std::size_t>> holding(n, std::vector<std::size_t>(H, 0));
  for (std::size_t v = 0; v < n; v++)
    for (std::size_t step = 0; step < H; step++) holding[v][step] = R.add(at(v, step), at(v, step + 1), storage, 0);
  for (const TimedSupply& p : supplies) {
    out.total += p.amount;
    biggest = std::max(biggest, p.amount);
    if (p.release <= H && p.amount > 0) R.add(S, at(p.vertex, p.release), p.amount, 0);
  }
  std::vector<std::vector<std::size_t>> collect(K, std::vector<std::size_t>(L));
  for (std::size_t k = 0; k < K; k++) {
    for (std::size_t step = 0; step < L; step++) collect[k][step] = R.add(at(sinks[k].vertex, step), S + 1 + k, 0, 0);
    R.add(S + 1 + k, T, sinks[k].capacity, 0);
    if (std::isfinite(sinks[k].capacity)) biggest = std::max(biggest, sinks[k].capacity);
  }
  if (std::isfinite(storage)) biggest = std::max(biggest, storage);
  R.eps = tolerance * std::max(1.0, biggest);
  out.safe.assign(L, 0);
  double value = 0;
  for (std::size_t step = 0; step < L; step++) {
    for (std::size_t k = 0; k < K; k++) R.forward(collect[k][step]).cap = flow_detail::inf;
    value += R.augment(S, T, who);
    out.safe[step] = value;
    if (out.quickest == OverTime::none && value >= out.total - R.eps * 10 && out.total > 0) out.quickest = step;
  }
  if (out.total == 0) out.quickest = 0;
  out.absorbed.assign(K, 0);
  for (std::size_t k = 0; k < K; k++)
    for (std::size_t step = 0; step < L; step++) out.absorbed[k] += R.flow(collect[k][step]);
  out.flow.assign(arcs.size(), std::vector<double>(L, 0));
  for (std::size_t a = 0; a < arcs.size(); a++)
    for (std::size_t step = 0; step < L; step++)
      if (moving[a][step] != static_cast<std::size_t>(-1)) out.flow[a][step] = R.flow(moving[a][step]);
  out.waiting.assign(n, std::vector<double>(L, 0));
  for (std::size_t v = 0; v < n; v++)
    for (std::size_t step = 0; step < H; step++) out.waiting[v][step] = R.flow(holding[v][step]);
  return out;
}

}  // namespace morphe::combinatorics
