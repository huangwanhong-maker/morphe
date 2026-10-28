// Checks of morphe/algebra/abelian_fourier.hpp. The transform of C_{n_1} × ⋯ × C_{n_k}, axis by
// axis, against the definition; inversion; the convolution theorem, the convolution computed
// both by its definition and as the product of the group algebra; Plancherel's identity; the
// characters χ_y against the Dixon–Schneider table of the same group. The fast Walsh–Hadamard
// transform against the definition for k ≤ 8, its inverse, Σ ĥ² = 2^k Σ f², the XOR convolution
// theorem, and agreement with the abelian transform of (C_2)^k.

#include <morphe/algebra/abelian_fourier.hpp>
#include <morphe/algebra/characters.hpp>
#include <morphe/algebra/group_algebra.hpp>
#include <morphe/check.hpp>

#include <bit>
#include <random>
#include <string>
#include <vector>

namespace ma = morphe::algebra;
using morphe::arithmetic::Cyclotomic;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

std::vector<Cyclotomic> random_function(std::mt19937_64& g, std::size_t n, bool complex) {
  std::vector<Cyclotomic> f(n);
  for (auto& x : f) {
    x = Cyclotomic(Rational(Integer(static_cast<long long>(g() % 9) - 4), Integer(static_cast<long long>(1 + g() % 2))));
    if (complex) x += Cyclotomic(static_cast<long long>(g() % 3) - 1) * Cyclotomic::zeta(4);
  }
  return f;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(61);

  // ---- finite abelian groups ---------------------------------------------------------------------------------------------------
  const std::vector<std::vector<std::size_t>> shapes = {{}, {1}, {2}, {7}, {12}, {2, 2}, {3, 3}, {4, 6}, {2, 3, 5}, {4, 2, 3}, {5, 5}, {2, 2, 2, 2}};
  int definition = 0, inversion = 0, convolution = 0, algebra = 0, plancherel = 0, table = 0;
  std::string where;
  for (const auto& moduli : shapes) {
    const ma::AbelianFourier F(moduli);
    const std::size_t n = F.order();
    const auto bad_before = definition + inversion + convolution + algebra + plancherel + table;
    for (int t = 0; t < 2; t++) {
      const auto f = random_function(g, n, t == 1), h = random_function(g, n, true);
      const auto fh = F.transform(f), hh = F.transform(h);
      if (fh != F.transform_by_definition(f)) definition++;
      if (F.inverse(fh) != f) inversion++;
      const auto c = F.convolve(f, h);
      const auto ch = F.transform(c);
      for (std::size_t y = 0; y < n; y++)
        if (!(ch[y] == fh[y] * hh[y])) convolution++;
      if (n <= 60) {
        const ma::GroupAlgebra<Cyclotomic> A(F.group());
        if (A.multiply(f, h) != c) algebra++;
      }
      // Σ_x f(x) conj h(x) = (1/|A|) Σ_y f̂(y) conj ĥ(y).
      Cyclotomic left = 0, right = 0;
      for (std::size_t x = 0; x < n; x++) left += f[x] * h[x].conj();
      for (std::size_t y = 0; y < n; y++) right += fh[y] * hh[y].conj();
      if (!(left * Cyclotomic(static_cast<long long>(n)) == right)) plancherel++;
    }
    // The characters χ_y are the rows of the Dixon–Schneider table of the same group.
    if (n <= 60) {
      const ma::FiniteGroup G = F.group();
      const auto T = ma::character_table(G);
      std::vector<char> used(T.size(), 0);
      for (std::size_t y = 0; y < n; y++) {
        bool found = false;
        for (std::size_t x = 0; x < T.size() && !found; x++) {
          if (used[x]) continue;
          bool eq = true;
          for (std::uint32_t e = 0; e < n && eq; e++) eq = T.value(x, e) == F.character(y, e);
          if (eq) used[x] = 1, found = true;
        }
        if (!found) table++;
      }
    }
    if (definition + inversion + convolution + algebra + plancherel + table != bad_before) {
      where += "(";
      for (auto m : moduli) where += std::to_string(m) + " ";
      where += ") ";
    }
  }
  s.require("abelian_fourier: the transform axis by axis equals the definition on 12 groups up to order 30", definition == 0, where);
  s.equal("abelian_fourier: the inverse recovers f (failures)", inversion, 0);
  s.equal("abelian_fourier: (f ∗ g)^ = f̂ ĝ (failures)", convolution, 0);
  s.equal("abelian_fourier: the convolution is the product of the group algebra of the same group (failures)", algebra, 0);
  s.equal("abelian_fourier: Plancherel, Σ f conj h = (1/|A|) Σ f̂ conj ĥ (failures)", plancherel, 0);
  s.equal("abelian_fourier: the characters χ_y are the rows of the Dixon–Schneider table (failures)", table, 0);

  // ---- the Walsh–Hadamard transform ---------------------------------------------------------------------------------------------
  {
    int def = 0, inv = 0, parseval = 0, xorconv = 0;
    for (int k = 0; k <= 10; k++) {
      const std::size_t n = std::size_t(1) << k;
      std::vector<Integer> f(n), h(n);
      for (auto& x : f) x = static_cast<long long>(g() % 21) - 10;
      for (auto& x : h) x = static_cast<long long>(g() % 21) - 10;
      const auto fh = ma::walsh_hadamard(f), hh = ma::walsh_hadamard(h);
      if (k <= 8)
        for (std::size_t y = 0; y < n; y++) {
          Integer want = 0;
          for (std::size_t x = 0; x < n; x++) want += std::popcount(x & y) % 2 ? -f[x] : f[x];
          if (want != fh[y]) def++;
        }
      const auto back = ma::inverse_walsh_hadamard(fh);
      for (std::size_t x = 0; x < n; x++)
        if (back[x] != Rational(f[x])) inv++;
      Integer a = 0, b = 0;
      for (std::size_t x = 0; x < n; x++) a += f[x] * f[x], b += fh[x] * fh[x];
      if (b != a * Integer(static_cast<long long>(n))) parseval++;
      if (k <= 7) {
        std::vector<Integer> c(n, Integer(0));
        for (std::size_t x = 0; x < n; x++)
          for (std::size_t z = 0; z < n; z++) c[x] += f[z] * h[x ^ z];
        const auto cc = ma::walsh_hadamard(c);
        for (std::size_t y = 0; y < n; y++)
          if (cc[y] != fh[y] * hh[y]) xorconv++;
      }
    }
    s.equal("walsh_hadamard: the butterflies equal Σ_x (−1)^{|x ∧ y|} f(x), k ≤ 8 (failures)", def, 0);
    s.equal("walsh_hadamard: the inverse, divided by 2^k, recovers f exactly, k ≤ 10 (failures)", inv, 0);
    s.equal("walsh_hadamard: Σ ĥ² = 2^k Σ f² (failures)", parseval, 0);
    s.equal("walsh_hadamard: the transform of the XOR convolution is the product, k ≤ 7 (failures)", xorconv, 0);
    bool agree = true;
    for (std::size_t k = 0; k <= 4; k++) {
      const std::size_t n = std::size_t(1) << k;
      std::vector<Integer> f(n);
      for (auto& x : f) x = static_cast<long long>(g() % 21) - 10;
      const auto wh = ma::walsh_hadamard(f);
      const auto af = ma::AbelianFourier(std::vector<std::size_t>(k, 2)).transform(std::vector<Rational>(f.begin(), f.end()));
      for (std::size_t y = 0; y < n; y++) agree &= af[y] == Cyclotomic(wh[y]);
    }
    s.require("walsh_hadamard: equals the abelian transform of (C_2)^k, k ≤ 4", agree);
  }
  return s.done();
}
