// The Fourier transform of a finite abelian group, exactly, and the fast Walsh–Hadamard
// transform of (ℤ/2)^k.
//
// A = C_{n_1} × ⋯ × C_{n_k}, its elements x = (x_1, …, x_k) with 0 ≤ x_i < n_i, numbered in
// mixed radix with x_1 most significant, as direct_product numbers them. With N = lcm(n_i), its
// characters are χ_y(x) = Π_i ζ_{n_i}^{x_i y_i} = ζ_N^{⟨x, y⟩}, ⟨x, y⟩ = Σ_i x_i y_i N/n_i, one for
// each y ∈ A, and the transform is that of the group (Diaconis 1988, ch. 2, the convention of
// symmetric_representations.hpp for one-dimensional representations):
//
//   f̂(y) = Σ_x f(x) ζ_N^{⟨x, y⟩},   f(x) = (1/|A|) Σ_y f̂(y) ζ_N^{−⟨x, y⟩},
//
// exactly in ℚ(ζ_N). (The discrete Fourier transform with the kernel exp(−2πi xy/n) is f̂(−y).)
// With (f ∗ g)(x) = Σ_z f(z) g(x − z), the product of ℚ(ζ_N)[A], (f ∗ g)^ = f̂ ĝ; and
// Plancherel's identity Σ_x f(x) conj g(x) = (1/|A|) Σ_y f̂(y) conj ĝ(y) holds. Both follow from
// the orthogonality of the characters, Σ_x χ_y(x) conj χ_{y'}(x) = |A| δ_{y y'}.
//
// The transform factors over the coordinates, since χ_y(x) is a product: it is computed axis by
// axis, a transform of length n_i along every line of axis i, O(|A| Σ_i n_i) operations in
// ℚ(ζ_N) (against O(|A|²) by the definition, which transform_by_definition computes for the
// checks). Exact.
//
// The Walsh–Hadamard transform is the case n_i = 2: ĥ(y) = Σ_x (−1)^{|x ∧ y|} f(x) over the bit
// vectors of length k, with values in ℤ. Its matrix is the Kronecker power H_2^{⊗k}, H_2 = [[1,
// 1], [1, −1]] (Sylvester 1867), so k passes of butterflies (a, b) ↦ (a + b, a − b) compute it in
// O(k 2^k) additions, in any ring; H_2^{⊗k} squared is 2^k I, so the inverse is the same
// transform divided by 2^k, exact over ℚ. (f ⊛ g)(x) = Σ_z f(z) g(x ⊕ z) is carried to the
// product.

#pragma once

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/arithmetic/cyclotomic.hpp>
#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/rational.hpp>

#include "finite_group.hpp"

namespace morphe::algebra {

using arithmetic::Cyclotomic;
using arithmetic::Rational;

class AbelianFourier {
 public:
  // C_{n_1} × ⋯ × C_{n_k}, every n_i ≥ 1; no moduli is the trivial group.
  explicit AbelianFourier(std::vector<std::size_t> moduli) : moduli_(std::move(moduli)) {
    for (auto m : moduli_) {
      if (m == 0) throw std::invalid_argument("morphe::algebra::AbelianFourier: a cyclic factor of order 0");
      order_ *= m;
      exponent_ = std::lcm(exponent_, m);
    }
    if (exponent_ > 0x7fffffffull) throw std::overflow_error("morphe::algebra::AbelianFourier: the exponent is too large");
    strides_.assign(moduli_.size(), 1);
    for (std::size_t i = moduli_.size(); i-- > 1;) strides_[i - 1] = strides_[i] * moduli_[i];
  }
  const std::vector<std::size_t>& moduli() const { return moduli_; }
  std::size_t order() const { return order_; }
  // N = lcm(n_i): the values lie in ℚ(ζ_N).
  std::size_t exponent() const { return exponent_; }

  std::vector<std::size_t> coordinates(std::size_t index) const {
    std::vector<std::size_t> x(moduli_.size());
    for (std::size_t i = 0; i < moduli_.size(); i++) x[i] = (index / strides_[i]) % moduli_[i];
    return x;
  }
  std::size_t index(const std::vector<std::size_t>& x) const {
    if (x.size() != moduli_.size()) throw std::invalid_argument("morphe::algebra::AbelianFourier::index: wrong number of coordinates");
    std::size_t k = 0;
    for (std::size_t i = 0; i < x.size(); i++) k += (x[i] % moduli_[i]) * strides_[i];
    return k;
  }
  // The group as a table: direct_product of the cyclic factors, numbered as here.
  FiniteGroup group() const {
    FiniteGroup G;
    for (auto m : moduli_) G = direct_product(G, cyclic(m));
    return G;
  }
  // ⟨x, y⟩ = Σ_i x_i y_i N/n_i mod N.
  std::size_t pairing(std::size_t x, std::size_t y) const {
    const auto a = coordinates(x), b = coordinates(y);
    std::size_t s = 0;
    for (std::size_t i = 0; i < a.size(); i++) s = (s + (a[i] * b[i] % moduli_[i]) * (exponent_ / moduli_[i])) % exponent_;
    return s;
  }
  // χ_y(x) = ζ_N^{⟨x, y⟩}.
  Cyclotomic character(std::size_t y, std::size_t x) const {
    return Cyclotomic::zeta(static_cast<unsigned>(exponent_), static_cast<long long>(pairing(x, y)));
  }

  // f̂, axis by axis.
  std::vector<Cyclotomic> transform(std::vector<Cyclotomic> f) const {
    check(f);
    along_axes(f, 1);
    return f;
  }
  std::vector<Cyclotomic> transform(const std::vector<Rational>& f) const { return transform(std::vector<Cyclotomic>(f.begin(), f.end())); }
  // f from f̂.
  std::vector<Cyclotomic> inverse(std::vector<Cyclotomic> fhat) const {
    check(fhat);
    along_axes(fhat, -1);
    const Cyclotomic scale = Cyclotomic(Rational(arithmetic::Integer(1), arithmetic::Integer(static_cast<long long>(order_))));
    for (auto& v : fhat) v *= scale;
    return fhat;
  }
  // f̂ by the definition, O(|A|²).
  std::vector<Cyclotomic> transform_by_definition(const std::vector<Cyclotomic>& f) const {
    check(f);
    std::vector<Cyclotomic> out(order_);
    for (std::size_t y = 0; y < order_; y++)
      for (std::size_t x = 0; x < order_; x++)
        if (!f[x].is_zero()) out[y] += f[x] * character(y, x);
    return out;
  }
  // (f ∗ g)(x) = Σ_z f(z) g(x − z), by the definition, O(|A|²).
  std::vector<Cyclotomic> convolve(const std::vector<Cyclotomic>& f, const std::vector<Cyclotomic>& g) const {
    check(f), check(g);
    std::vector<Cyclotomic> out(order_);
    for (std::size_t z = 0; z < order_; z++) {
      if (f[z].is_zero()) continue;
      const auto cz = coordinates(z);
      for (std::size_t u = 0; u < order_; u++) {
        if (g[u].is_zero()) continue;
        const auto cu = coordinates(u);
        std::vector<std::size_t> x(cz.size());
        for (std::size_t i = 0; i < x.size(); i++) x[i] = (cz[i] + cu[i]) % moduli_[i];
        out[index(x)] += f[z] * g[u];
      }
    }
    return out;
  }

 private:
  void check(const std::vector<Cyclotomic>& f) const {
    if (f.size() != order_) throw std::invalid_argument("morphe::algebra::AbelianFourier: a function of the wrong length");
  }
  // Along each axis i, every line v ↦ (Σ_x v_x ζ_N^{sign · x y N/n_i})_y.
  void along_axes(std::vector<Cyclotomic>& f, int sign) const {
    const unsigned N = static_cast<unsigned>(exponent_);
    for (std::size_t axis = 0; axis < moduli_.size(); axis++) {
      const std::size_t n = moduli_[axis], stride = strides_[axis];
      if (n == 1) continue;
      std::vector<Cyclotomic> root(n);
      for (std::size_t t = 0; t < n; t++) root[t] = Cyclotomic::zeta(N, sign * static_cast<long long>(t * (exponent_ / n)));
      std::vector<Cyclotomic> line(n);
      for (std::size_t outer = 0; outer < order_; outer += n * stride)
        for (std::size_t inner = 0; inner < stride; inner++) {
          const std::size_t base = outer + inner;
          for (std::size_t y = 0; y < n; y++) {
            Cyclotomic s = 0;
            for (std::size_t x = 0; x < n; x++) {
              const Cyclotomic& v = f[base + x * stride];
              if (!v.is_zero()) s += v * root[(x * y) % n];
            }
            line[y] = std::move(s);
          }
          for (std::size_t y = 0; y < n; y++) f[base + y * stride] = std::move(line[y]);
        }
    }
  }

  std::vector<std::size_t> moduli_;
  std::size_t order_ = 1, exponent_ = 1;
  std::vector<std::size_t> strides_;
};

inline std::vector<Cyclotomic> abelian_fourier(const std::vector<std::size_t>& moduli, const std::vector<Cyclotomic>& f) {
  return AbelianFourier(moduli).transform(f);
}
inline std::vector<Cyclotomic> inverse_abelian_fourier(const std::vector<std::size_t>& moduli, const std::vector<Cyclotomic>& fhat) {
  return AbelianFourier(moduli).inverse(fhat);
}

// ĥ(y) = Σ_x (−1)^{|x ∧ y|} f(x) on a vector of length 2^k, in any ring R with + and −.
template <class R>
std::vector<R> walsh_hadamard(std::vector<R> f) {
  const std::size_t n = f.size();
  if (n == 0 || (n & (n - 1)) != 0) throw std::invalid_argument("morphe::algebra::walsh_hadamard: the length is not a power of 2");
  for (std::size_t len = 1; len < n; len <<= 1)
    for (std::size_t i = 0; i < n; i += 2 * len)
      for (std::size_t j = i; j < i + len; j++) {
        R a = f[j], b = f[j + len];
        f[j] = a + b;
        f[j + len] = a - b;
      }
  return f;
}

// f from ĥ: the same transform divided by 2^k, exactly.
inline std::vector<Rational> inverse_walsh_hadamard(const std::vector<arithmetic::Integer>& fhat) {
  const auto h = walsh_hadamard(fhat);
  std::vector<Rational> f;
  f.reserve(h.size());
  const arithmetic::Integer n(static_cast<long long>(h.size()));
  for (const auto& x : h) f.emplace_back(x, n);
  return f;
}

}  // namespace morphe::algebra
