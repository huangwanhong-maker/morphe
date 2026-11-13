// Checks of acyclic relations and forward chaining (acyclic.hpp).
//
// - Orders: a diamond is ordered with every arc forward; a relation with a cycle is refused with a cycle whose
//   every step is an arc; random acyclic relations of 2 to 12 elements are ordered, and one arc added against
//   the order makes a cycle that is found.
// - Cones: on a chain, down(x) is what lies before x with its distance; a bound on the steps cuts it; on
//   random relations down and up agree with reachability computed by repeated squaring of the relation.
// - The covering relation: of the full order on a chain it is the chain's steps; on random acyclic relations
//   it has the same closure as the relation, and removing any of its arcs changes the closure.
// - Quotients: arcs gathered between classes and inside them sum to the arcs there are.
// - Forward chaining: a rule with two premises needs both; alternatives are any of them; what rests on an
//   element is what every derivation passes through, against brute force over all subsets of the givens' use
//   (closure by naive iteration to a fixed point) on random rules.

#include <morphe/check.hpp>
#include <morphe/combinatorics/acyclic.hpp>

#include <algorithm>
#include <random>
#include <set>
#include <vector>

namespace mc = morphe::combinatorics;
using mc::Arc;

namespace {

std::vector<std::vector<char>> closure_by_squaring(std::size_t n, const std::vector<Arc>& arcs) {
  std::vector<std::vector<char>> R(n, std::vector<char>(n, 0));
  for (const auto& [a, b] : arcs) R[a][b] = 1;
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = 0; b < n; b++)
        if (R[a][k] && R[k][b]) R[a][b] = 1;
  return R;
}

std::vector<Arc> random_acyclic(std::mt19937& gen, std::size_t n, double p) {
  std::vector<std::size_t> rank(n);
  for (std::size_t i = 0; i < n; i++) rank[i] = i;
  std::shuffle(rank.begin(), rank.end(), gen);
  std::bernoulli_distribution coin(p);
  std::vector<Arc> arcs;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      if (rank[a] < rank[b] && coin(gen)) arcs.push_back({a, b});
  return arcs;
}

std::vector<char> naive_closure(std::size_t n, const std::vector<mc::Rule>& rules, const std::vector<std::size_t>& givens, long withheld) {
  std::vector<char> had(n, 0);
  for (std::size_t g : givens)
    if (static_cast<long>(g) != withheld) had[g] = 1;
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& r : rules) {
      if (static_cast<long>(r.conclusion) == withheld || had[r.conclusion]) continue;
      bool all = true;
      for (std::size_t p : r.premises) all = all && had[p];
      if (all) { had[r.conclusion] = 1; changed = true; }
    }
  }
  return had;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937 gen(20261113);

  {
    const std::vector<Arc> diamond{{0, 1}, {0, 2}, {1, 3}, {2, 3}};
    const auto o = mc::topological_order(4, diamond);
    std::vector<std::size_t> at(4);
    for (std::size_t i = 0; i < o.order.size(); i++) at[o.order[i]] = i;
    bool forward = o.acyclic && o.order.size() == 4;
    for (const auto& [a, b] : diamond) forward = forward && at[a] < at[b];
    s.require("a diamond is ordered with every arc forward", forward);
  }
  {
    const std::vector<Arc> arcs{{0, 1}, {1, 2}, {2, 3}, {3, 1}, {3, 4}};
    const auto o = mc::topological_order(5, arcs);
    bool steps = !o.acyclic && o.cycle.size() == 3;
    for (std::size_t i = 0; steps && i < o.cycle.size(); i++)
      steps = std::find(arcs.begin(), arcs.end(), Arc{o.cycle[i], o.cycle[(i + 1) % o.cycle.size()]}) != arcs.end();
    s.require("a relation with a cycle is refused with a cycle whose every step is an arc", steps);
  }
  {
    bool ordered = true, found = true;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 2 + trial % 11;
      auto arcs = random_acyclic(gen, n, 0.3);
      const auto o = mc::topological_order(n, arcs);
      std::vector<std::size_t> at(n);
      for (std::size_t i = 0; i < o.order.size(); i++) at[o.order[i]] = i;
      ordered = ordered && o.acyclic;
      for (const auto& [a, b] : arcs) ordered = ordered && at[a] < at[b];
      if (arcs.empty()) continue;
      const auto R = closure_by_squaring(n, arcs);
      // An arc back along a path closes a cycle.
      const Arc pick = arcs[static_cast<std::size_t>(trial) % arcs.size()];
      arcs.push_back({pick.second, pick.first});
      const auto c = mc::topological_order(n, arcs);
      bool is_cycle = !c.acyclic && !c.cycle.empty();
      for (std::size_t i = 0; is_cycle && i < c.cycle.size(); i++)
        is_cycle = std::find(arcs.begin(), arcs.end(), Arc{c.cycle[i], c.cycle[(i + 1) % c.cycle.size()]}) != arcs.end();
      found = found && is_cycle;
      (void)R;
    }
    s.require("200 random acyclic relations are ordered, every arc forward", ordered);
    s.require("an arc added against the order makes a cycle, and it is found", found);
  }
  {
    std::vector<Arc> chain;
    for (std::size_t i = 0; i + 1 < 6; i++) chain.push_back({i, i + 1});
    const auto d = mc::down(6, chain, 4);
    bool ok = d.size() == 4;
    for (const auto& r : d) ok = ok && r.steps == 4 - r.element;
    s.require("on a chain, down(x) is what lies before x, each at its distance", ok);
    s.require("a bound on the steps cuts the cone", mc::down(6, chain, 4, 2).size() == 2 && mc::up(6, chain, 1, 3).size() == 3);
  }
  {
    bool agree = true, same_closure = true, minimal = true, sums = true;
    for (int trial = 0; trial < 150; trial++) {
      const std::size_t n = 2 + trial % 10;
      const auto arcs = random_acyclic(gen, n, 0.35);
      const auto R = closure_by_squaring(n, arcs);
      agree = agree && mc::reachability(n, arcs) == R;
      for (std::size_t x = 0; x < n; x++) {
        std::set<std::size_t> below, above;
        for (const auto& r : mc::down(n, arcs, x)) below.insert(r.element);
        for (const auto& r : mc::up(n, arcs, x)) above.insert(r.element);
        for (std::size_t y = 0; y < n; y++) agree = agree && (below.count(y) == (R[y][x] ? 1u : 0u)) && (above.count(y) == (R[x][y] ? 1u : 0u));
      }
      const auto cover = mc::reduction(n, arcs);
      same_closure = same_closure && closure_by_squaring(n, cover) == R;
      for (std::size_t k = 0; k < cover.size(); k++) {
        auto fewer = cover;
        fewer.erase(fewer.begin() + static_cast<long>(k));
        minimal = minimal && closure_by_squaring(n, fewer) != R;
      }
      std::vector<std::size_t> cls(n);
      for (std::size_t v = 0; v < n; v++) cls[v] = v % 3;
      const auto q = mc::quotient(n, arcs, cls, 3);
      std::size_t total = 0;
      for (const auto& g : q.arcs) total += g.arcs;
      for (std::size_t c : q.inside) total += c;
      sums = sums && total == arcs.size();
    }
    s.require("down and up agree with the closure computed by Warshall's algorithm", agree);
    s.require("the covering relation has the relation's closure", same_closure);
    s.require("and no arc of it can be spared", minimal);
    s.require("a quotient's arcs, between classes and inside them, are the arcs there are", sums);
  }
  {
    std::vector<Arc> full;
    for (std::size_t a = 0; a < 5; a++)
      for (std::size_t b = a + 1; b < 5; b++) full.push_back({a, b});
    const auto cover = mc::reduction(5, full);
    bool steps = cover.size() == 4;
    for (const auto& [a, b] : cover) steps = steps && b == a + 1;
    s.require("the covering relation of a chain's whole order is its steps", steps);
    const auto R = mc::reachability(5, full);
    s.require("on a chain every two elements are comparable; in an antichain none", mc::comparable(R, 3, 1) && !mc::comparable(mc::reachability(3, {}), 0, 2));
  }
  {
    // 3 needs 0 and 1 together; 4 follows from 2 or from 3; 5 from 4.
    const std::vector<mc::Rule> rules{{{0, 1}, 3}, {{2}, 4}, {{3}, 4}, {{4}, 5}};
    const auto both = mc::closure(6, rules, {0, 1});
    const auto one = mc::closure(6, rules, {0});
    s.require("a rule with two premises needs both", both[3] && both[4] && both[5] && !one[3] && !one[5]);
    s.require("with two ways to an element, nothing after it rests on either way", mc::rests_on(6, rules, {0, 1, 2}, 3).empty() && mc::rests_on(6, rules, {0, 1, 2}, 2).empty());
    const auto r = mc::rests_on(6, rules, {0, 1}, 0);
    s.require("with one way, all of it rests on each premise", r == std::vector<std::size_t>({3, 4, 5}));
    s.require("a rule with no premise gives its conclusion outright", mc::closure(2, {{{}, 1}}, {})[1] == 1);
  }
  {
    bool same = true, rest = true;
    for (int trial = 0; trial < 300; trial++) {
      const std::size_t n = 3 + trial % 8;
      std::uniform_int_distribution<std::size_t> any(0, n - 1), width(1, 3);
      std::vector<mc::Rule> rules(n + trial % 5);
      for (auto& r : rules) {
        r.conclusion = any(gen);
        for (std::size_t k = width(gen); k > 0; k--) r.premises.push_back(any(gen));
      }
      std::vector<std::size_t> givens{any(gen), any(gen)};
      same = same && mc::closure(n, rules, givens) == naive_closure(n, rules, givens, -1);
      const std::size_t c = any(gen);
      const auto with = naive_closure(n, rules, givens, -1), without = naive_closure(n, rules, givens, static_cast<long>(c));
      std::vector<std::size_t> expect;
      for (std::size_t v = 0; v < n; v++)
        if (v != c && with[v] && !without[v]) expect.push_back(v);
      rest = rest && mc::rests_on(n, rules, givens, c) == expect;
    }
    s.require("forward chaining agrees with iteration to a fixed point on 300 random sets of rules", same);
    s.require("what rests on an element agrees with withholding it by brute force", rest);
  }

  return s.done();
}
