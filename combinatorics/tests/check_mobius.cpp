// Checks of Möbius functions and inversion (mobius.hpp).
//
// - The Boolean lattice of 4 elements: μ(S, T) = (−1)^{|T∖S|}; a chain of 6: μ = 1, −1 on the diagonal and
//   the step above it, 0 beyond; the divisors of 360 under division: μ(a, b) the number-theoretic μ(b/a).
// - Inversion: on random posets (transitive closures of random acyclic relations), with morphe's integers,
//   g = invert(zeta(g)) exactly.
// - The subset transforms: subset_mobius ∘ subset_zeta = identity on random integers, n = 0 … 8, and each
//   transform equal to its definition summed term by term.
// - Anchored interactions: an additive function has no interaction of two or more factors; a product of
//   binary factors anchored at 0 has the top one alone; the pair term is the interaction contrast
//   V(1,1) − V(1,0) − V(0,1) + V(0,0); the decomposition sums back to V at every subset.

#include <morphe/arithmetic/integer.hpp>
#include <morphe/check.hpp>
#include <morphe/combinatorics/mobius.hpp>

#include <bit>
#include <cstdint>
#include <random>
#include <vector>

namespace mc = morphe::combinatorics;
using morphe::arithmetic::Integer;

namespace {

int mu_number(int n) {
  int r = 1;
  for (int p = 2; p * p <= n; p++)
    if (n % p == 0) {
      n /= p;
      if (n % p == 0) return 0;
      r = -r;
    }
  if (n > 1) r = -r;
  return r;
}

}  // namespace

int main() {
  morphe::Suite s;

  {
    const std::size_t N = 16;
    std::vector<std::vector<char>> leq(N, std::vector<char>(N, 0));
    for (std::size_t a = 0; a < N; a++)
      for (std::size_t b = 0; b < N; b++) leq[a][b] = (a & b) == a;
    const auto mu = mc::mobius<long long>(leq);
    bool ok = true;
    for (std::size_t a = 0; a < N; a++)
      for (std::size_t b = 0; b < N; b++)
        ok = ok && mu[a][b] == ((a & b) == a ? (std::popcount(b & ~a) % 2 ? -1 : 1) : 0);
    s.require("the Boolean lattice of 4 elements: μ(S, T) = (−1)^{|T∖S|}", ok);
  }
  {
    const std::size_t N = 6;
    std::vector<std::vector<char>> leq(N, std::vector<char>(N, 0));
    for (std::size_t a = 0; a < N; a++)
      for (std::size_t b = a; b < N; b++) leq[a][b] = 1;
    const auto mu = mc::mobius<long long>(leq);
    bool ok = true;
    for (std::size_t a = 0; a < N; a++)
      for (std::size_t b = 0; b < N; b++) ok = ok && mu[a][b] == (b == a ? 1 : b == a + 1 ? -1 : 0);
    s.require("a chain: 1 on the diagonal, −1 above it, 0 beyond", ok);
  }
  {
    std::vector<int> d;
    for (int k = 1; k <= 360; k++)
      if (360 % k == 0) d.push_back(k);
    std::vector<std::vector<char>> leq(d.size(), std::vector<char>(d.size(), 0));
    for (std::size_t a = 0; a < d.size(); a++)
      for (std::size_t b = 0; b < d.size(); b++) leq[a][b] = d[b] % d[a] == 0;
    const auto mu = mc::mobius<long long>(leq);
    bool ok = true;
    for (std::size_t a = 0; a < d.size(); a++)
      for (std::size_t b = 0; b < d.size(); b++)
        if (leq[a][b]) ok = ok && mu[a][b] == mu_number(d[b] / d[a]);
    s.require("the 24 divisors of 360: μ(a, b) = μ(b/a), the number-theoretic Möbius function", ok && d.size() == 24);
  }
  {
    std::mt19937_64 rng(7);
    bool ok = true;
    for (int trial = 0; trial < 40; trial++) {
      const std::size_t n = 2 + rng() % 9;
      std::vector<std::vector<char>> leq(n, std::vector<char>(n, 0));
      for (std::size_t a = 0; a < n; a++) {
        leq[a][a] = 1;
        for (std::size_t b = a + 1; b < n; b++) leq[a][b] = rng() % 3 == 0;
      }
      for (std::size_t k = 0; k < n; k++)
        for (std::size_t a = 0; a < n; a++)
          for (std::size_t b = 0; b < n; b++)
            if (leq[a][k] && leq[k][b]) leq[a][b] = 1;
      std::vector<Integer> g(n), f(n, Integer(0));
      for (auto& x : g) x = Integer(static_cast<long long>(rng() % 2001) - 1000);
      for (std::size_t y = 0; y < n; y++)
        for (std::size_t x = 0; x < n; x++)
          if (leq[x][y]) f[y] = f[y] + g[x];
      ok = ok && mc::invert(leq, f) == g;
    }
    s.require("Möbius inversion recovers g from its sums on 40 random posets, exactly", ok);
  }
  {
    std::mt19937_64 rng(11);
    bool inverse = true, defined = true;
    for (std::size_t n = 0; n <= 8; n++) {
      const std::size_t N = std::size_t{1} << n;
      std::vector<Integer> g(N);
      for (auto& x : g) x = Integer(static_cast<long long>(rng() % 201) - 100);
      const auto f = mc::subset_zeta(g);
      inverse = inverse && mc::subset_mobius(f) == g;
      for (std::size_t T = 0; T < N; T++) {
        Integer z(0), m(0);
        for (std::size_t S = 0; S < N; S++)
          if ((S & T) == S) {
            z = z + g[S];
            m = m + (std::popcount(T & ~S) % 2 ? Integer(0) - f[S] : f[S]);
          }
        defined = defined && z == f[T] && m == g[T];
      }
    }
    s.require("subset_mobius inverts subset_zeta, n = 0 … 8", inverse);
    s.require("each equals its definition summed term by term", defined);
  }
  {
    const auto additive = mc::anchored_interactions<long long>(5, [](std::uint32_t S) { long long v = 3; for (int i = 0; i < 5; i++) if (S >> i & 1) v += 2 * i + 1; return v; });
    bool none = true;
    for (std::size_t S = 0; S < additive.size(); S++)
      if (std::popcount(S) >= 2) none = none && additive[S] == 0;
    s.require("an additive function has no interaction of two or more factors", none && additive[0] == 3 && additive[1 << 4] == 9);
    const auto product = mc::anchored_interactions<long long>(4, [](std::uint32_t S) { return S == 15u ? 1LL : 0LL; });
    bool top = true;
    for (std::size_t S = 0; S < product.size(); S++) top = top && product[S] == (S == 15 ? 1 : 0);
    s.require("a product of binary factors anchored at 0: the top interaction alone", top);
    const auto V = [](std::uint32_t S) { const int a = S & 1, b = (S >> 1) & 1, c = (S >> 2) & 1; return 5.0 + 2 * a - b + 3 * a * b - 4 * a * b * c + 0.5 * c; };
    const auto g = mc::anchored_interactions<double>(3, V);
    s.check("the pair term is the interaction contrast V(1,1) − V(1,0) − V(0,1) + V(0,0)", g[3] - (V(3) - V(1) - V(2) + V(0)), 1e-12);
    s.check("the triple term is −4", g[7] + 4, 1e-12);
    double worst = 0;
    const auto back = mc::subset_zeta(g);
    for (std::uint32_t S = 0; S < 8; S++) worst = std::max(worst, std::abs(back[S] - V(S)));
    s.check("the interactions sum back to V at every subset", worst, 1e-12);
  }

  return s.done();
}
