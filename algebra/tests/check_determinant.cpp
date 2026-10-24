// Checks of morphe/algebra/determinant.hpp: Bareiss's elimination against the Leibniz
// expansion on random matrices with many zeros (so that pivots must be exchanged); the sign
// of every permutation of five letters; multiplicativity on matrices whose determinants
// exceed 64 bits; Vandermonde's product; singular matrices; the Hilbert matrices over ℚ
// against the closed form det Hₙ = cₙ⁴ / c₂ₙ with cₙ = ∏_{i<n} i!; and minors against the
// expansion along a row.

#include <morphe/algebra/determinant.hpp>
#include <morphe/check.hpp>

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace ma = morphe::algebra;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

// Σ_π sgn(π) ∏ a_{i π(i)}, over all n! permutations.
Integer leibniz(const ma::Matrix& A) {
  const std::size_t n = A.rows;
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), 0);
  Integer sum = 0;
  do {
    std::size_t inversions = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) inversions += p[i] > p[j];
    Integer term = 1;
    for (std::size_t i = 0; i < n && !term.is_zero(); i++) term *= A(i, p[i]);
    sum += inversions % 2 ? -term : term;
  } while (std::next_permutation(p.begin(), p.end()));
  return sum;
}

ma::Matrix random_matrix(std::mt19937_64& g, std::size_t n, long long spread, double zeros) {
  ma::Matrix A(n, n);
  std::uniform_real_distribution<double> u(0, 1);
  for (auto& x : A.a) x = u(g) < zeros ? 0 : static_cast<long long>(g() % static_cast<unsigned long long>(2 * spread + 1)) - spread;
  return A;
}

// cₙ = ∏_{i=1}^{n-1} i!.
Integer superfactorial(std::size_t n) {
  Integer c = 1;
  for (std::size_t i = 1; i < n; i++) c *= morphe::arithmetic::factorial(i);
  return c;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(1968);

  // Bareiss against Leibniz, n = 0 to 6, half the entries zero.
  {
    int bad = 0, tried = 0;
    for (int t = 0; t < 400; t++) {
      const std::size_t n = static_cast<std::size_t>(t % 7);
      const ma::Matrix A = random_matrix(g, n, 9, t % 2 ? 0.5 : 0.1);
      tried++;
      if (ma::determinant(A) != leibniz(A)) bad++;
    }
    s.equal("determinant: Bareiss equals the Leibniz expansion on 400 random matrices up to 6 × 6, many with zero pivots (failures)", bad, 0);
  }

  // Permutation matrices: the determinant is the sign.
  {
    std::vector<std::size_t> p(5);
    std::iota(p.begin(), p.end(), 0);
    int bad = 0, count = 0;
    do {
      ma::Matrix P(5, 5);
      for (std::size_t i = 0; i < 5; i++) P(i, p[i]) = 1;
      std::size_t inversions = 0;
      for (std::size_t i = 0; i < 5; i++)
        for (std::size_t j = i + 1; j < 5; j++) inversions += p[i] > p[j];
      if (ma::determinant(P) != Integer(inversions % 2 ? -1 : 1)) bad++;
      count++;
    } while (std::next_permutation(p.begin(), p.end()));
    s.require("determinant: every permutation matrix of order 5 has the sign of its permutation (120 of them)", bad == 0 && count == 120);
  }

  // det(AB) = det A det B, entries up to 10⁶, so that the determinants exceed 64 bits.
  {
    int bad = 0, big = 0;
    for (int t = 0; t < 30; t++) {
      const std::size_t n = 3 + static_cast<std::size_t>(t % 6);
      const ma::Matrix A = random_matrix(g, n, 1000000, 0.05), B = random_matrix(g, n, 1000000, 0.05);
      const Integer dA = ma::determinant(A), dB = ma::determinant(B);
      if (ma::determinant(A * B) != dA * dB) bad++;
      if (!(dA * dB).fits_int64()) big++;
    }
    s.require("determinant: det(AB) = det A det B on 30 random pairs up to 8 × 8 with entries up to 10^6", bad == 0, std::to_string(big) + " of the products beyond 64 bits");
  }

  // Vandermonde: det [x_i^j] = ∏_{i<j} (x_j − x_i).
  {
    int bad = 0;
    for (int t = 0; t < 10; t++) {
      const std::size_t n = 4 + static_cast<std::size_t>(t);
      std::vector<long long> x;
      while (x.size() < n) {
        const long long v = static_cast<long long>(g() % 2001) - 1000;
        if (std::find(x.begin(), x.end(), v) == x.end()) x.push_back(v);
      }
      ma::Matrix V(n, n);
      for (std::size_t i = 0; i < n; i++) {
        Integer p = 1;
        for (std::size_t j = 0; j < n; j++) {
          V(i, j) = p;
          p *= Integer(x[i]);
        }
      }
      Integer product = 1;
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++) product *= Integer(x[j] - x[i]);
      if (ma::determinant(V) != product) bad++;
    }
    s.equal("determinant: Vandermonde matrices of orders 4 to 13 give ∏_{i<j} (x_j − x_i), exactly (failures)", bad, 0);
  }

  // Singular matrices: a product through a lower dimension, and a repeated row.
  {
    int bad = 0;
    for (int t = 0; t < 40; t++) {
      const std::size_t n = 2 + static_cast<std::size_t>(t % 7);
      ma::Matrix X(n, n - 1), Y(n - 1, n);
      for (auto& v : X.a) v = static_cast<long long>(g() % 21) - 10;
      for (auto& v : Y.a) v = static_cast<long long>(g() % 21) - 10;
      if (!ma::determinant(X * Y).is_zero()) bad++;
      ma::Matrix R = random_matrix(g, n, 50, 0.2);
      for (std::size_t j = 0; j < n; j++) R(n - 1, j) = R(0, j);
      if (!ma::determinant(R).is_zero()) bad++;
    }
    s.equal("determinant: zero for 40 products of rank n − 1 and 40 matrices with a repeated row (failures)", bad, 0);
  }

  // Hilbert matrices H_n = [1/(i + j + 1)]: det H_n = c_n^4 / c_{2n} (Hilbert 1894; Choi 1983,
  // "Tricks or treats with the Hilbert matrix", Amer. Math. Monthly 90, 301–312). The first
  // reciprocals, 1, 12, 2160, 6048000, 266716800000, are OEIS A005249.
  {
    int bad = 0;
    for (std::size_t n = 1; n <= 12; n++) {
      std::vector<Rational> h(n * n);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) h[i * n + j] = Rational(Integer(1), Integer(static_cast<long long>(i + j + 1)));
      const Integer c = superfactorial(n);
      const Rational want(c * c * c * c, superfactorial(2 * n));
      if (ma::determinant(n, h) != want) bad++;
    }
    s.equal("determinant: Hilbert matrices H_1 to H_12 over ℚ equal c_n^4 / c_2n (failures)", bad, 0);
    std::vector<Rational> h5(25);
    for (std::size_t i = 0; i < 5; i++)
      for (std::size_t j = 0; j < 5; j++) h5[i * 5 + j] = Rational(Integer(1), Integer(static_cast<long long>(i + j + 1)));
    s.equal("determinant: det H_5 = 1/266716800000 (OEIS A005249)", ma::determinant(5, h5).str(), std::string("1/266716800000"));
  }

  // Minors: the Laplace expansion along the first row recovers the determinant.
  {
    int bad = 0;
    for (int t = 0; t < 50; t++) {
      const std::size_t n = 1 + static_cast<std::size_t>(t % 7);
      const ma::Matrix A = random_matrix(g, n, 30, 0.3);
      Integer sum = 0;
      for (std::size_t j = 0; j < n; j++) {
        const Integer m = n == 1 ? Integer(1) : ma::minor(A, 0, j);
        sum += (j % 2 ? -A(0, j) : A(0, j)) * m;
      }
      if (sum != ma::determinant(A)) bad++;
    }
    s.equal("determinant: the Laplace expansion along a row by minors equals the determinant on 50 matrices (failures)", bad, 0);
  }

  s.equal("determinant: the empty matrix has determinant 1", ma::determinant(ma::Matrix(0, 0)), Integer(1));
  return s.done();
}
