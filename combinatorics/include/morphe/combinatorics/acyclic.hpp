// Acyclic relations and what follows from several things together: the order a directed acyclic relation
// generates, its cones, its covering relation and its quotients; and forward chaining over rules with several
// premises, with what rests on a given element.
//
// A relation on {0, …, n−1} is a list of arcs (a, b), read "a comes before b". It is acyclic exactly when it
// has a topological order, and then its reflexive–transitive closure is a partial order (Kahn 1962, Comm. ACM
// 5, 558–562, whose algorithm is the one below; Knuth, The Art of Computer Programming vol. 1, §2.2.3). When
// it is not, a cycle is returned as the place where it fails.
//
//   - down(x): what comes before x, each element with its least number of arcs from x (breadth first);
//     up(x) the same forward. The down-set of x in the order is {x} ∪ down(x).
//   - reduction: the covering relation, the least relation with the same closure, unique for an acyclic
//     relation (Aho, Garey and Ullman 1972, SIAM J. Comput. 1, 131–137).
//   - quotient by a map to classes: the arcs between different classes, each with the number of arcs it
//     gathers. It need not be acyclic, and is tested like any relation.
//   - comparable(x, y): whether one comes before the other.
//
// A rule has a set of premises and a conclusion: the conclusion is had once every premise is. Forward
// chaining from a set of givens reaches the least set closed under the rules, in time linear in the size of
// the rules (Dowling and Gallier 1984, J. Logic Programming 1, 267–284: the satisfiability of Horn formulae).
// An element with several rules may be had by any of them. What rests on c among what is derivable: the
// elements that are derivable from the givens and are no longer when c is withheld, which is every element
// all of whose derivations pass through c.
//
// Costs: order, down, up, quotient O(n + m); reduction and reachability O(n·m); closure O(size of the rules);
// rests_on one closure. Checks: tests/check_acyclic.cpp.

#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

using Arc = std::pair<std::size_t, std::size_t>;

namespace detail {
inline void within(std::size_t n, const std::vector<Arc>& arcs, const char* who) {
  for (const auto& [a, b] : arcs)
    if (a >= n || b >= n) throw std::invalid_argument(std::string(who) + ": an arc leaves the set");
}
inline std::vector<std::vector<std::size_t>> lists(std::size_t n, const std::vector<Arc>& arcs, bool forward) {
  std::vector<std::vector<std::size_t>> out(n);
  for (const auto& [a, b] : arcs) (forward ? out[a] : out[b]).push_back(forward ? b : a);
  return out;
}
}  // namespace detail

// A topological order of an acyclic relation, or the cycle that forbids one.
struct Ordering {
  bool acyclic = false;
  std::vector<std::size_t> order;  // every element once, each after all that come before it (when acyclic)
  std::vector<std::size_t> cycle;  // otherwise a cycle c0 → c1 → … → ck → c0
};

inline Ordering topological_order(std::size_t n, const std::vector<Arc>& arcs) {
  detail::within(n, arcs, "morphe::combinatorics::topological_order");
  const auto next = detail::lists(n, arcs, true);
  std::vector<std::size_t> waiting(n, 0);
  for (const auto& [a, b] : arcs) waiting[b]++;
  Ordering out;
  std::deque<std::size_t> ready;
  for (std::size_t v = 0; v < n; v++)
    if (!waiting[v]) ready.push_back(v);
  while (!ready.empty()) {
    const std::size_t v = ready.front();
    ready.pop_front();
    out.order.push_back(v);
    for (std::size_t w : next[v])
      if (--waiting[w] == 0) ready.push_back(w);
  }
  out.acyclic = out.order.size() == n;
  if (out.acyclic) return out;
  // Every element left has one before it that is left too: walk back until an element repeats.
  out.order.clear();
  const auto before = detail::lists(n, arcs, false);
  std::size_t v = 0;
  while (!waiting[v]) v++;
  std::vector<long> seen(n, -1);
  std::vector<std::size_t> walk;
  while (seen[v] < 0) {
    seen[v] = static_cast<long>(walk.size());
    walk.push_back(v);
    for (std::size_t u : before[v])
      if (waiting[u]) { v = u; break; }
  }
  out.cycle.assign(walk.begin() + seen[v], walk.end());
  std::reverse(out.cycle.begin(), out.cycle.end());
  return out;
}

// What is reached from x along the arcs (forward) or against them, with the least number of arcs to each;
// x itself is left out. `most` bounds the number of arcs (0: no bound).
struct Reached {
  std::size_t element = 0, steps = 0;
};

inline std::vector<Reached> reach(std::size_t n, const std::vector<Arc>& arcs, std::size_t x, bool forward, std::size_t most = 0) {
  detail::within(n, arcs, "morphe::combinatorics::reach");
  if (x >= n) throw std::invalid_argument("morphe::combinatorics::reach: no such element");
  const auto next = detail::lists(n, arcs, forward);
  std::vector<char> seen(n, 0);
  seen[x] = 1;
  std::vector<Reached> out;
  std::deque<Reached> queue{{x, 0}};
  while (!queue.empty()) {
    const Reached at = queue.front();
    queue.pop_front();
    if (most && at.steps == most) continue;
    for (std::size_t w : next[at.element])
      if (!seen[w]) {
        seen[w] = 1;
        out.push_back({w, at.steps + 1});
        queue.push_back({w, at.steps + 1});
      }
  }
  return out;
}
inline std::vector<Reached> down(std::size_t n, const std::vector<Arc>& arcs, std::size_t x, std::size_t most = 0) { return reach(n, arcs, x, false, most); }
inline std::vector<Reached> up(std::size_t n, const std::vector<Arc>& arcs, std::size_t x, std::size_t most = 0) { return reach(n, arcs, x, true, most); }

// before[a][b]: whether a path of one arc or more leads from a to b.
inline std::vector<std::vector<char>> reachability(std::size_t n, const std::vector<Arc>& arcs) {
  std::vector<std::vector<char>> R(n, std::vector<char>(n, 0));
  for (std::size_t a = 0; a < n; a++)
    for (const Reached& r : reach(n, arcs, a, true)) R[a][r.element] = 1;
  return R;
}

inline bool comparable(const std::vector<std::vector<char>>& before, std::size_t x, std::size_t y) { return x == y || before[x][y] || before[y][x]; }

// The covering relation of an acyclic relation: its arcs (a, b), each once, with no longer path from a to b.
inline std::vector<Arc> reduction(std::size_t n, const std::vector<Arc>& arcs) {
  if (!topological_order(n, arcs).acyclic) throw std::invalid_argument("morphe::combinatorics::reduction: the relation has a cycle");
  const auto R = reachability(n, arcs);
  const auto next = detail::lists(n, arcs, true);
  std::vector<Arc> out;
  for (std::size_t a = 0; a < n; a++) {
    std::vector<std::size_t> to = next[a];
    std::sort(to.begin(), to.end());
    to.erase(std::unique(to.begin(), to.end()), to.end());
    for (std::size_t b : to) {
      bool longer = false;
      for (std::size_t c : to)
        if (c != b && R[c][b]) { longer = true; break; }
      if (!longer) out.push_back({a, b});
    }
  }
  return out;
}

// The relation among classes that a map of the elements induces: an arc between two different classes for
// every pair joined by an arc, with how many arcs it gathers. Arcs inside a class are counted apart.
struct Gathered {
  std::size_t from = 0, to = 0, arcs = 0;
};
struct Quotient {
  std::size_t classes = 0;
  std::vector<Gathered> arcs;        // between different classes, in order of (from, to)
  std::vector<std::size_t> inside;   // by class, the arcs both of whose ends lie in it
};

inline Quotient quotient(std::size_t n, const std::vector<Arc>& arcs, const std::vector<std::size_t>& cls, std::size_t classes) {
  detail::within(n, arcs, "morphe::combinatorics::quotient");
  if (cls.size() != n) throw std::invalid_argument("morphe::combinatorics::quotient: a class for every element");
  for (std::size_t c : cls)
    if (c >= classes) throw std::invalid_argument("morphe::combinatorics::quotient: a class out of range");
  Quotient q;
  q.classes = classes;
  q.inside.assign(classes, 0);
  std::map<Arc, std::size_t> count;
  for (const auto& [a, b] : arcs) {
    if (cls[a] == cls[b]) q.inside[cls[a]]++;
    else count[{cls[a], cls[b]}]++;
  }
  for (const auto& [arc, k] : count) q.arcs.push_back({arc.first, arc.second, k});
  return q;
}

// A rule: its conclusion is had once all its premises are. A rule with no premise gives its conclusion outright.
struct Rule {
  std::vector<std::size_t> premises;
  std::size_t conclusion = 0;
};

// The least set holding the givens and closed under the rules, leaving out the elements withheld (which are
// never had, whatever gives them). had[v] says whether v is in it.
inline std::vector<char> closure(std::size_t n, const std::vector<Rule>& rules, const std::vector<std::size_t>& givens,
                                 const std::vector<std::size_t>& withheld = {}) {
  std::vector<char> had(n, 0), barred(n, 0);
  for (std::size_t w : withheld) {
    if (w >= n) throw std::invalid_argument("morphe::combinatorics::closure: no such element");
    barred[w] = 1;
  }
  std::vector<std::size_t> missing(rules.size(), 0);
  std::vector<std::vector<std::size_t>> uses(n);
  std::deque<std::size_t> fresh;
  const auto have = [&](std::size_t v) {
    if (v >= n) throw std::invalid_argument("morphe::combinatorics::closure: no such element");
    if (!had[v] && !barred[v]) { had[v] = 1; fresh.push_back(v); }
  };
  for (std::size_t r = 0; r < rules.size(); r++) {
    std::vector<std::size_t> p = rules[r].premises;
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    missing[r] = p.size();
    for (std::size_t v : p) {
      if (v >= n) throw std::invalid_argument("morphe::combinatorics::closure: no such element");
      uses[v].push_back(r);
    }
  }
  for (std::size_t g : givens) have(g);
  for (std::size_t r = 0; r < rules.size(); r++)
    if (!missing[r]) have(rules[r].conclusion);
  while (!fresh.empty()) {
    const std::size_t v = fresh.front();
    fresh.pop_front();
    for (std::size_t r : uses[v])
      if (--missing[r] == 0) have(rules[r].conclusion);
  }
  return had;
}

// What rests on c: the elements had from the givens that are not had once c is withheld; c itself is left out.
inline std::vector<std::size_t> rests_on(std::size_t n, const std::vector<Rule>& rules, const std::vector<std::size_t>& givens, std::size_t c) {
  const auto with = closure(n, rules, givens), without = closure(n, rules, givens, {c});
  std::vector<std::size_t> out;
  for (std::size_t v = 0; v < n; v++)
    if (v != c && with[v] && !without[v]) out.push_back(v);
  return out;
}

}  // namespace morphe::combinatorics
