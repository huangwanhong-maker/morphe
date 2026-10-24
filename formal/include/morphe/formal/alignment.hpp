// The alignment of two words as a sum over the paths of their edit graph, in any semiring: the
// definition every fast alignment algorithm is a schedule of.
//
// The edit graph of x = x₁ … x_m and y = y₁ … y_n is the directed acyclic graph on the
// vertices (i, j), 0 ≤ i ≤ m, 0 ≤ j ≤ n, with three kinds of edge into (i, j):
//
//   M  from (i − 1, j − 1), x_i aligned with y_j;
//   I  from (i − 1, j),     x_i against a gap (a letter of x not in y);
//   D  from (i, j − 1),     y_j against a gap.
//
// An alignment is a path, and the path is its word of operations over {M, I, D}. The four
// modes say where a path may begin and end (Durbin et al. 1998, Biological Sequence Analysis,
// §2.3; Gusfield 1997, Algorithms on Strings, Trees and Sequences, §11.6):
//
//   global      from (0, 0) to (m, n)                               (Needleman and Wunsch 1970)
//   local       from any vertex to any vertex, the empty path       (Smith and Waterman 1981)
//               included
//   semiglobal  from the first row to the last: x end to end, y's ends free
//   overlap     from the first row or column to the last row or column: free end gaps on both
//
// A scheme of weights gives each path a weight in a semiring S, the ⊗-product of its edge
// weights in the order of the path, and the alignment value is the ⊕-sum over all the paths
// the mode allows (Mohri 2002, J. Automata, Languages and Combinatorics 7, 321–350: the
// shortest-distance problem of a weighted acyclic automaton; Pachter and Sturmfels 2005,
// Algebraic Statistics for Computational Biology, §2.2). Two schemes:
//
//   linear       one state: pair(a, b) on M, gap_x(a) on I, gap_y(b) on D; a gap's weight is the
//                product of its letters' weights.
//   three-state  the product of the edit graph with an automaton on the states M, I, D (Gotoh
//                1982, J. Mol. Biol. 162, 705–708; the pair HMM of Durbin et al. 1998, §4.1):
//                transition[from][to] for from ∈ {begin, M, I, D} and to ∈ {end, M, I, D}, then
//                the state's emission, pair(a, b), emit_x(a) or emit_y(b). Every word over
//                {M, I, D} is one path of the automaton, so the paths, and their number, are
//                the same as the linear scheme's. affine() builds Gotoh's scheme: entering I or D
//                from any other state (or at the begin) costs open, each letter of a gap
//                extend, so a gap of length L weighs open ⊗ extend^L.
//
// The value is computed by the dynamic programme over the vertices in row-major order, which is
// a topological order of the edit graph:
//
//   linear       V(i, j) = [source] ⊕ V(i−1, j−1) ⊗ pair(x_i, y_j) ⊕ V(i−1, j) ⊗ gap_x(x_i)
//                          ⊕ V(i, j−1) ⊗ gap_y(y_j),
//   three-state  V_s(i, j) = enter_s(u) ⊗ emission_s, with u the predecessor by s and
//                enter_s(u) = [source u] ⊗ T[begin][s] ⊕ ⊕_r V_r(u) ⊗ T[r][s],
//
// where [source] is the semiring's one at a vertex where the mode lets a path begin and zero
// elsewhere; the value is the ⊕-sum over the sinks (with T[r][end] in the three-state scheme).
// Local alignment is therefore the semiring's one ⊕-added at every vertex, which in the tropical
// semiring is Smith and Waterman's max(0, ·).
//
// Correctness: by induction along the topological order, V(i, j) is the sum over the allowed
// paths ending at (i, j) of their weights; distributivity moves the common factor of the last
// edge out of the sum (Mohri 2002, Theorem 1 for acyclic graphs). Termination: always; the graph
// is finite and acyclic. Complexity: O(mn) semiring operations (4 per vertex for the linear
// scheme, 16 for the three-state); path_sum keeps two rows (O(n) values), optimal_path the
// whole table (O(mn) values).
//
// Traceback (optimal_path) in a selective semiring: a sum equals one of its terms, so from the
// first sink in row-major order whose value is the total, the first term (in the order start,
// M, I, D) that equals each value is followed back to a source. The enumerator lists every path
// the mode allows, for checking the dynamic programme on short words; the number of global
// paths is the Delannoy number D(m, n) = Σ_k C(m, k) C(n, k) 2^k (Banderier and Schwer 2005,
// J. Statistical Planning and Inference 135, 40–54).
//
// Words are any random-access sized ranges of letters (std::vector, std::span, std::string);
// the weights are callables on letters, so the domain passes its byte codes without copying.

#pragma once

#include <morphe/arithmetic/integer.hpp>
#include <morphe/formal/semiring.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::formal {

enum class Mode { global, local, semiglobal, overlap };

inline bool is_source(Mode mode, std::size_t i, std::size_t j) {
  switch (mode) {
    case Mode::global: return i == 0 && j == 0;
    case Mode::local: return true;
    case Mode::semiglobal: return i == 0;
    case Mode::overlap: return i == 0 || j == 0;
  }
  return false;
}
inline bool is_sink(Mode mode, std::size_t i, std::size_t j, std::size_t m, std::size_t n) {
  switch (mode) {
    case Mode::global: return i == m && j == n;
    case Mode::local: return true;
    case Mode::semiglobal: return i == m;
    case Mode::overlap: return i == m || j == n;
  }
  return false;
}

// ---- paths ----------------------------------------------------------------------------------------
struct EditRun {
  char op = 'M';
  std::size_t length = 0;
  friend bool operator==(const EditRun&, const EditRun&) = default;
};

// A path of the edit graph: where it begins and its word of operations.
struct EditPath {
  std::size_t x_begin = 0, y_begin = 0;
  std::string ops;  // over {M, I, D}
  std::size_t x_end() const { return x_begin + static_cast<std::size_t>(std::count_if(ops.begin(), ops.end(), [](char c) { return c != 'D'; })); }
  std::size_t y_end() const { return y_begin + static_cast<std::size_t>(std::count_if(ops.begin(), ops.end(), [](char c) { return c != 'I'; })); }
  // The word as runs of one operation.
  std::vector<EditRun> runs() const {
    std::vector<EditRun> r;
    for (char c : ops) {
      if (!r.empty() && r.back().op == c) r.back().length++;
      else r.push_back({c, 1});
    }
    return r;
  }
  friend bool operator==(const EditPath&, const EditPath&) = default;
};

// ---- schemes of weights -------------------------------------------------------------------------------
template <Semiring S, class Pair, class GapX, class GapY>
struct LinearWeights {
  using semiring = S;
  static constexpr int states = 1;
  Pair pair;    // pair(a, b): x's letter a aligned with y's letter b
  GapX gap_x;   // gap_x(a): x's letter a against a gap (I)
  GapY gap_y;   // gap_y(b): y's letter b against a gap (D)
};
template <Semiring S, class Pair, class GapX, class GapY>
LinearWeights<S, Pair, GapX, GapY> linear(Pair pair, GapX gap_x, GapY gap_y) {
  return {std::move(pair), std::move(gap_x), std::move(gap_y)};
}

// The states of the three-state scheme: index 0 is the begin (as a source of transitions) and
// the end (as a target).
inline constexpr int STATE_BEGIN = 0, STATE_END = 0, STATE_M = 1, STATE_I = 2, STATE_D = 3;

template <Semiring S, class Pair, class EmitX, class EmitY>
struct ThreeStateWeights {
  using semiring = S;
  using value_type = typename S::value_type;
  static constexpr int states = 3;
  std::array<std::array<value_type, 4>, 4> transition;  // [from: begin, M, I, D][to: end, M, I, D]
  Pair pair;
  EmitX emit_x;
  EmitY emit_y;
};
template <Semiring S, class Pair, class EmitX, class EmitY>
ThreeStateWeights<S, Pair, EmitX, EmitY> three_state(std::array<std::array<typename S::value_type, 4>, 4> transition, Pair pair, EmitX emit_x, EmitY emit_y) {
  return {std::move(transition), std::move(pair), std::move(emit_x), std::move(emit_y)};
}
// Affine gaps (Gotoh 1982): entering I or D from another state, or at the begin, weighs `open`;
// each letter of a gap weighs extend_x(a) or extend_y(b).
template <Semiring S, class Pair, class ExtendX, class ExtendY>
ThreeStateWeights<S, Pair, ExtendX, ExtendY> affine(Pair pair, const typename S::value_type& open, ExtendX extend_x, ExtendY extend_y) {
  const auto u = S::one();
  std::array<std::array<typename S::value_type, 4>, 4> t;
  //          end  M  I     D
  t[STATE_BEGIN] = {u, u, open, open};
  t[STATE_M] = {u, u, open, open};
  t[STATE_I] = {u, u, u, open};
  t[STATE_D] = {u, u, open, u};
  return {std::move(t), std::move(pair), std::move(extend_x), std::move(extend_y)};
}

// The weight of one path under a scheme. Throws std::out_of_range when the path leaves the words.
template <class W, class X, class Y>
typename W::semiring::value_type path_weight(const EditPath& p, const X& x, const Y& y, const W& w) {
  using S = typename W::semiring;
  auto v = S::one();
  std::size_t i = p.x_begin, j = p.y_begin;
  int state = STATE_BEGIN;
  for (char op : p.ops) {
    const int next = op == 'M' ? STATE_M : op == 'I' ? STATE_I : STATE_D;
    if ((next != STATE_D && i >= x.size()) || (next != STATE_I && j >= y.size())) throw std::out_of_range("morphe::formal::path_weight: the path leaves the words");
    if constexpr (W::states == 3) v = S::times(v, w.transition[state][next]);
    if (next == STATE_M) {
      v = S::times(v, w.pair(x[i], y[j]));
      i++;
      j++;
    } else if (next == STATE_I) {
      if constexpr (W::states == 1) v = S::times(v, w.gap_x(x[i]));
      else v = S::times(v, w.emit_x(x[i]));
      i++;
    } else {
      if constexpr (W::states == 1) v = S::times(v, w.gap_y(y[j]));
      else v = S::times(v, w.emit_y(y[j]));
      j++;
    }
    state = next;
  }
  if constexpr (W::states == 3) v = S::times(v, w.transition[state][STATE_END]);
  return v;
}

namespace detail {

// The values of the three states at one vertex.
template <class V>
struct Cell3 {
  V m, i, d;
};

// enter_s(u) for the three-state scheme.
template <class S, class T, class V>
V enter(const T& transition, int s, bool source, const Cell3<V>& u) {
  V v = source ? transition[STATE_BEGIN][s] : S::zero();
  v = S::plus(v, S::times(u.m, transition[STATE_M][s]));
  v = S::plus(v, S::times(u.i, transition[STATE_I][s]));
  v = S::plus(v, S::times(u.d, transition[STATE_D][s]));
  return v;
}
template <class S, class T, class V>
V finish(const T& transition, bool source, const Cell3<V>& c) { return enter<S>(transition, STATE_END, source, c); }

}  // namespace detail

// The ⊕-sum over the paths the mode allows of their weights, in O(n) memory.
template <class W, class X, class Y>
typename W::semiring::value_type path_sum(const X& x, const Y& y, const W& w, Mode mode) {
  using S = typename W::semiring;
  using V = typename S::value_type;
  const std::size_t m = x.size(), n = y.size();
  V total = S::zero();
  if constexpr (W::states == 1) {
    std::vector<V> prev(n + 1, S::zero()), cur(n + 1, S::zero());
    for (std::size_t i = 0; i <= m; i++) {
      for (std::size_t j = 0; j <= n; j++) {
        V v = is_source(mode, i, j) ? S::one() : S::zero();
        if (i && j) v = S::plus(v, S::times(prev[j - 1], w.pair(x[i - 1], y[j - 1])));
        if (i) v = S::plus(v, S::times(prev[j], w.gap_x(x[i - 1])));
        if (j) v = S::plus(v, S::times(cur[j - 1], w.gap_y(y[j - 1])));
        if (is_sink(mode, i, j, m, n)) total = S::plus(total, v);
        cur[j] = std::move(v);
      }
      std::swap(prev, cur);
    }
  } else {
    using C = detail::Cell3<V>;
    const C none{S::zero(), S::zero(), S::zero()};
    std::vector<C> prev(n + 1, none), cur(n + 1, none);
    for (std::size_t i = 0; i <= m; i++) {
      for (std::size_t j = 0; j <= n; j++) {
        C c = none;
        if (i && j) c.m = S::times(detail::enter<S>(w.transition, STATE_M, is_source(mode, i - 1, j - 1), prev[j - 1]), w.pair(x[i - 1], y[j - 1]));
        if (i) c.i = S::times(detail::enter<S>(w.transition, STATE_I, is_source(mode, i - 1, j), prev[j]), w.emit_x(x[i - 1]));
        if (j) c.d = S::times(detail::enter<S>(w.transition, STATE_D, is_source(mode, i, j - 1), cur[j - 1]), w.emit_y(y[j - 1]));
        if (is_sink(mode, i, j, m, n)) total = S::plus(total, detail::finish<S>(w.transition, is_source(mode, i, j), c));
        cur[j] = std::move(c);
      }
      std::swap(prev, cur);
    }
  }
  return total;
}

// An optimal path and the optimum, in a selective semiring; `exists` is false when the sum is the
// semiring's zero (no path of nonzero weight).
template <class V>
struct Optimum {
  V value;
  EditPath path;
  bool exists = false;
};

template <class W, class X, class Y>
  requires SelectiveSemiring<typename W::semiring>
Optimum<typename W::semiring::value_type> optimal_path(const X& x, const Y& y, const W& w, Mode mode) {
  using S = typename W::semiring;
  using V = typename S::value_type;
  const std::size_t m = x.size(), n = y.size();
  const auto at = [n](std::size_t i, std::size_t j) { return i * (n + 1) + j; };
  Optimum<V> out{S::zero(), {}, false};
  const auto fail = [] { throw std::logic_error("morphe::formal::optimal_path: no term attains a sum (the semiring is not selective)"); };
  std::string rev;
  std::size_t i = 0, j = 0;
  if constexpr (W::states == 1) {
    std::vector<V> T((m + 1) * (n + 1), S::zero());
    for (std::size_t a = 0; a <= m; a++)
      for (std::size_t b = 0; b <= n; b++) {
        V v = is_source(mode, a, b) ? S::one() : S::zero();
        if (a && b) v = S::plus(v, S::times(T[at(a - 1, b - 1)], w.pair(x[a - 1], y[b - 1])));
        if (a) v = S::plus(v, S::times(T[at(a - 1, b)], w.gap_x(x[a - 1])));
        if (b) v = S::plus(v, S::times(T[at(a, b - 1)], w.gap_y(y[b - 1])));
        T[at(a, b)] = std::move(v);
        if (is_sink(mode, a, b, m, n)) out.value = S::plus(out.value, T[at(a, b)]);
      }
    if (S::equal(out.value, S::zero())) return out;
    bool found = false;
    for (std::size_t a = 0; a <= m && !found; a++)
      for (std::size_t b = 0; b <= n && !found; b++)
        if (is_sink(mode, a, b, m, n) && S::equal(T[at(a, b)], out.value)) { i = a; j = b; found = true; }
    for (;;) {
      const V& v = T[at(i, j)];
      if (is_source(mode, i, j) && S::equal(v, S::one())) break;
      if (i && j && S::equal(v, S::times(T[at(i - 1, j - 1)], w.pair(x[i - 1], y[j - 1])))) { rev += 'M'; i--; j--; continue; }
      if (i && S::equal(v, S::times(T[at(i - 1, j)], w.gap_x(x[i - 1])))) { rev += 'I'; i--; continue; }
      if (j && S::equal(v, S::times(T[at(i, j - 1)], w.gap_y(y[j - 1])))) { rev += 'D'; j--; continue; }
      fail();
    }
  } else {
    using C = detail::Cell3<V>;
    const C none{S::zero(), S::zero(), S::zero()};
    std::vector<C> T((m + 1) * (n + 1), none);
    for (std::size_t a = 0; a <= m; a++)
      for (std::size_t b = 0; b <= n; b++) {
        C c = none;
        if (a && b) c.m = S::times(detail::enter<S>(w.transition, STATE_M, is_source(mode, a - 1, b - 1), T[at(a - 1, b - 1)]), w.pair(x[a - 1], y[b - 1]));
        if (a) c.i = S::times(detail::enter<S>(w.transition, STATE_I, is_source(mode, a - 1, b), T[at(a - 1, b)]), w.emit_x(x[a - 1]));
        if (b) c.d = S::times(detail::enter<S>(w.transition, STATE_D, is_source(mode, a, b - 1), T[at(a, b - 1)]), w.emit_y(y[b - 1]));
        T[at(a, b)] = std::move(c);
        if (is_sink(mode, a, b, m, n)) out.value = S::plus(out.value, detail::finish<S>(w.transition, is_source(mode, a, b), T[at(a, b)]));
      }
    if (S::equal(out.value, S::zero())) return out;
    bool found = false;
    for (std::size_t a = 0; a <= m && !found; a++)
      for (std::size_t b = 0; b <= n && !found; b++)
        if (is_sink(mode, a, b, m, n) && S::equal(detail::finish<S>(w.transition, is_source(mode, a, b), T[at(a, b)]), out.value)) { i = a; j = b; found = true; }
    // The state the path ends in: the first term of the final sum that attains it.
    const auto pick = [&](int to, std::size_t a, std::size_t b, const V& target) -> int {
      const C& c = T[at(a, b)];
      if (is_source(mode, a, b) && S::equal(target, w.transition[STATE_BEGIN][to])) return STATE_BEGIN;
      if (S::equal(target, S::times(c.m, w.transition[STATE_M][to]))) return STATE_M;
      if (S::equal(target, S::times(c.i, w.transition[STATE_I][to]))) return STATE_I;
      if (S::equal(target, S::times(c.d, w.transition[STATE_D][to]))) return STATE_D;
      fail();
      return -1;
    };
    int state = pick(STATE_END, i, j, out.value);
    while (state != STATE_BEGIN) {
      std::size_t a = i, b = j;
      if (state == STATE_M) { rev += 'M'; a--; b--; }
      else if (state == STATE_I) { rev += 'I'; a--; }
      else { rev += 'D'; b--; }
      const V e = detail::enter<S>(w.transition, state, is_source(mode, a, b), T[at(a, b)]);
      state = pick(state, a, b, e);
      i = a;
      j = b;
    }
  }
  std::reverse(rev.begin(), rev.end());
  out.path.x_begin = i;
  out.path.y_begin = j;
  out.path.ops = std::move(rev);
  out.exists = true;
  return out;
}

// Every path the mode allows on words of lengths m and n, each given to visit(const EditPath&):
// by source in row-major order, then depth first, a path being visited at every sink it reaches.
// The number visited is the sum of path_sum in the counting semiring with weights 1.
template <class F>
void enumerate_paths(std::size_t m, std::size_t n, Mode mode, F&& visit) {
  EditPath p;
  const auto go = [&](auto&& self, std::size_t i, std::size_t j) -> void {
    if (is_sink(mode, i, j, m, n)) visit(static_cast<const EditPath&>(p));
    if (i < m && j < n) {
      p.ops.push_back('M');
      self(self, i + 1, j + 1);
      p.ops.pop_back();
    }
    if (i < m) {
      p.ops.push_back('I');
      self(self, i + 1, j);
      p.ops.pop_back();
    }
    if (j < n) {
      p.ops.push_back('D');
      self(self, i, j + 1);
      p.ops.pop_back();
    }
  };
  for (std::size_t i = 0; i <= m; i++)
    for (std::size_t j = 0; j <= n; j++)
      if (is_source(mode, i, j)) {
        p.x_begin = i;
        p.y_begin = j;
        p.ops.clear();
        go(go, i, j);
      }
}

// The Delannoy number D(m, n) = Σ_k C(m, k) C(n, k) 2^k: the number of global alignments of words
// of lengths m and n (lattice paths with steps (1, 0), (0, 1), (1, 1)).
inline arithmetic::Integer delannoy(std::size_t m, std::size_t n) {
  arithmetic::Integer d = 0;
  for (std::size_t k = 0; k <= std::min(m, n); k++)
    d += arithmetic::binomial(static_cast<long long>(m), static_cast<long long>(k)) * arithmetic::binomial(static_cast<long long>(n), static_cast<long long>(k)) *
         arithmetic::pow(arithmetic::Integer(2), k);
  return d;
}

}  // namespace morphe::formal
