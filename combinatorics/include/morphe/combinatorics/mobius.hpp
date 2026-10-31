// Möbius functions of finite posets and Möbius inversion: on any finite poset, and fast on the lattice of
// subsets, where it gives the unique decomposition of a function of several factors into the interactions
// of each set of them relative to a reference.
//
// Definition (Rota 1964, Z. Wahrscheinlichkeitstheorie 2, 340–368; Stanley 2012, Enumerative Combinatorics
// vol. 1, 2nd ed., §3.7). On a finite poset P the Möbius function is the inverse of the zeta function in the
// incidence algebra: μ(x, x) = 1, μ(x, y) = −Σ_{x ≤ z < y} μ(x, z) for x < y, and 0 unless x ≤ y. Möbius
// inversion: f(y) = Σ_{x ≤ y} g(x) for all y exactly when g(y) = Σ_{x ≤ y} μ(x, y) f(x) (Stanley,
// Prop. 3.7.1). Philip Hall's theorem (Stanley, §3.8): μ(x, y) is the reduced Euler characteristic of
// the order complex of the open interval (x, y); so on the face poset of a regular CW complex, every
// interval being a sphere, μ(σ, τ) = (−1)^{dim τ − dim σ} (the poset is Eulerian; Stanley §3.16).
//
// On the subsets of {0, …, n−1} (held as bit masks), μ(S, T) = (−1)^{|T∖S|} for S ⊆ T (Stanley, §3.8),
// and the two transforms
//
//   zeta:    f(T) = Σ_{S ⊆ T} g(S),            mobius:  g(T) = Σ_{S ⊆ T} (−1)^{|T∖S|} f(S)
//
// are inverse, computed in O(n 2ⁿ) additions by summing over one coordinate at a time (Yates 1937, The
// Design and Analysis of Factorial Experiments, whose effects of a 2ⁿ experiment these are; the "fast zeta
// transform"). Read with f(S) = V(x_S, x⁰_{rest}) — a function of n factors with the factors in S
// set to x and the others held at a reference x⁰ — g(S) is the interaction of the factors S anchored at x⁰,
// and V(x) = Σ_S g(S) is the unique decomposition whose partial sums reconstruct every f(S) (Grabisch,
// Marichal and Roubens 2000, Math. Oper. Res. 25, 157, the Möbius representation of a set function; in
// statistics the anchored decomposition of Kuo, Sloan, Wasilkowski and Woźniakowski 2010). g({i, j}) is the
// interaction contrast V(1,1) − V(1,0) − V(0,1) + V(0,0); its sign depends on the reference and on the
// function, not on anything else.
//
// Exact over any ring T with +, − (morphe's Integer, Rational) or over double. Cost: mobius on a poset of n
// elements O(n³); the subset transforms O(n 2ⁿ), n ≤ 30. Checks: tests/check_mobius.cpp.

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

// μ(x, y) of a finite poset given by its order: leq[x][y] true when x ≤ y (reflexive, antisymmetric,
// transitive; not checked beyond reflexivity). Rows index x, columns y; 0 off the order.
template <class T = long long>
std::vector<std::vector<T>> mobius(const std::vector<std::vector<char>>& leq) {
  const std::size_t n = leq.size();
  for (std::size_t x = 0; x < n; x++) {
    if (leq[x].size() != n) throw std::invalid_argument("morphe::combinatorics::mobius: a square order");
    if (!leq[x][x]) throw std::invalid_argument("morphe::combinatorics::mobius: the order is not reflexive");
  }
  // Process y in an order compatible with ≤ (by the number of elements below it).
  std::vector<std::size_t> below(n, 0), order(n);
  for (std::size_t y = 0; y < n; y++)
    for (std::size_t z = 0; z < n; z++)
      if (leq[z][y]) below[y]++;
  for (std::size_t i = 0; i < n; i++) order[i] = i;
  for (std::size_t i = 1; i < n; i++)
    for (std::size_t j = i; j > 0 && below[order[j - 1]] > below[order[j]]; j--) std::swap(order[j - 1], order[j]);
  std::vector<std::vector<T>> mu(n, std::vector<T>(n, T(0)));
  for (std::size_t x = 0; x < n; x++)
    for (std::size_t y : order) {
      if (!leq[x][y]) continue;
      if (x == y) { mu[x][y] = T(1); continue; }
      T s(0);
      for (std::size_t z = 0; z < n; z++)
        if (z != y && leq[x][z] && leq[z][y]) s = s + mu[x][z];
      mu[x][y] = T(0) - s;
    }
  return mu;
}

// g(y) = Σ_{x ≤ y} μ(x, y) f(x): the inverse of f(y) = Σ_{x ≤ y} g(x).
template <class T>
std::vector<T> invert(const std::vector<std::vector<char>>& leq, const std::vector<T>& f) {
  const auto mu = mobius<T>(leq);
  std::vector<T> g(f.size(), T(0));
  for (std::size_t y = 0; y < f.size(); y++)
    for (std::size_t x = 0; x < f.size(); x++)
      if (leq[x][y]) g[y] = g[y] + mu[x][y] * f[x];
  return g;
}

// f(T) = Σ_{S ⊆ T} g(S) over the subsets of n elements, indexed by bit mask (size 2ⁿ).
template <class T>
std::vector<T> subset_zeta(std::vector<T> g) {
  const std::size_t N = g.size();
  if (N == 0 || (N & (N - 1))) throw std::invalid_argument("morphe::combinatorics::subset_zeta: 2ⁿ values");
  for (std::size_t bit = 1; bit < N; bit <<= 1)
    for (std::size_t S = 0; S < N; S++)
      if (S & bit) g[S] = g[S] + g[S ^ bit];
  return g;
}

// g(T) = Σ_{S ⊆ T} (−1)^{|T∖S|} f(S): the interactions of each set of factors.
template <class T>
std::vector<T> subset_mobius(std::vector<T> f) {
  const std::size_t N = f.size();
  if (N == 0 || (N & (N - 1))) throw std::invalid_argument("morphe::combinatorics::subset_mobius: 2ⁿ values");
  for (std::size_t bit = 1; bit < N; bit <<= 1)
    for (std::size_t S = 0; S < N; S++)
      if (S & bit) f[S] = f[S] - f[S ^ bit];
  return f;
}

// The anchored interactions of a function of n binary factors: value(S) is V with the factors in S active
// and the rest at their reference; returns g(S) for every S, g(∅) = V(reference).
template <class T, class F>
std::vector<T> anchored_interactions(std::size_t n, F value) {
  if (n > 30) throw std::invalid_argument("morphe::combinatorics::anchored_interactions: at most 30 factors");
  std::vector<T> f(std::size_t{1} << n);
  for (std::size_t S = 0; S < f.size(); S++) f[S] = value(static_cast<std::uint32_t>(S));
  return subset_mobius(std::move(f));
}

}  // namespace morphe::combinatorics
