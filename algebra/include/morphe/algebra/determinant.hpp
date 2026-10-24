// The determinant of a square matrix over ℤ, exactly, by Bareiss's fraction-free elimination;
// over ℚ, by clearing the denominators of each row; and the minors that the matrix-tree
// theorem needs.
//
// Bareiss's algorithm (Bareiss 1968, "Sylvester's identity and multistep integer-preserving
// Gaussian elimination", Math. Comp. 22, 565–578, the one-step method of §2). With
// a⁽⁰⁾ = A and a⁽⁻¹⁾₋₁,₋₁ = 1, step k (k = 0, …, n − 2) computes, for i, j > k,
//
//   a⁽ᵏ⁺¹⁾ᵢⱼ = (a⁽ᵏ⁾ₖₖ a⁽ᵏ⁾ᵢⱼ − a⁽ᵏ⁾ᵢₖ a⁽ᵏ⁾ₖⱼ) / a⁽ᵏ⁻¹⁾ₖ₋₁,ₖ₋₁,
//
// and the division is exact: by Sylvester's identity a⁽ᵏ⁺¹⁾ᵢⱼ is the minor of A on the rows
// 0, …, k, i and the columns 0, …, k, j (Bareiss 1968, eq. (2.4)–(2.6); Gantmacher 1959, The
// Theory of Matrices 1, §II.3). Every intermediate entry is therefore a minor of A, bounded by
// Hadamard's inequality, and the last pivot a⁽ⁿ⁻¹⁾ₙ₋₁,ₙ₋₁ is det A. When a pivot is zero, a
// lower row with a non-zero entry in the pivot column is exchanged with the pivot row, which
// changes the sign of every later minor and so of the result; when there is none, the leading
// columns are dependent and det A = 0.
//
// Termination and cost: it terminates on every square matrix, with O(n³) multiplications and
// exact divisions of integers of at most O(n log(n ‖A‖∞)) bits (Hadamard's bound), that is
// O(n⁵ log²(n ‖A‖∞)) bit operations with morphe's schoolbook multiplication. The determinant
// of the 0 × 0 matrix is 1.
//
// Over ℚ the determinant is multilinear in the rows: with dᵢ the least common multiple of the
// denominators of row i, det A = det(D A) / ∏ dᵢ, and D A is integral. Exact; the result is
// reduced.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/rational.hpp>

#include "matrix.hpp"

namespace morphe::algebra {

using arithmetic::Rational;

inline Integer determinant(Matrix A) {
  if (A.rows != A.cols) throw std::invalid_argument("morphe::algebra::determinant: a square matrix");
  const std::size_t n = A.rows;
  if (n == 0) return 1;
  bool negate = false;
  Integer previous = 1;
  for (std::size_t k = 0; k + 1 < n; k++) {
    if (A(k, k).is_zero()) {
      std::size_t p = k + 1;
      while (p < n && A(p, k).is_zero()) p++;
      if (p == n) return 0;
      // Columns before k are no longer read.
      for (std::size_t j = k; j < n; j++) std::swap(A(k, j), A(p, j));
      negate = !negate;
    }
    const Integer pivot = A(k, k);
    const bool unit = previous.is_one();
    for (std::size_t i = k + 1; i < n; i++) {
      const Integer aik = A(i, k);
      for (std::size_t j = k + 1; j < n; j++) {
        Integer x = pivot * A(i, j);
        if (!aik.is_zero()) x -= aik * A(k, j);
        // Exact by Sylvester's identity.
        A(i, j) = unit ? std::move(x) : x / previous;
      }
    }
    previous = pivot;
  }
  return negate ? -A(n - 1, n - 1) : A(n - 1, n - 1);
}

// The matrix with row i and column j deleted.
inline Matrix minor_matrix(const Matrix& A, std::size_t i, std::size_t j) {
  if (i >= A.rows || j >= A.cols) throw std::out_of_range("morphe::algebra::minor_matrix");
  Matrix M(A.rows - 1, A.cols - 1);
  for (std::size_t r = 0, rr = 0; r < A.rows; r++) {
    if (r == i) continue;
    for (std::size_t c = 0, cc = 0; c < A.cols; c++) {
      if (c == j) continue;
      M(rr, cc++) = A(r, c);
    }
    rr++;
  }
  return M;
}

// The (i, j) minor: the determinant of A with row i and column j deleted.
inline Integer minor(const Matrix& A, std::size_t i, std::size_t j) { return determinant(minor_matrix(A, i, j)); }

// The determinant of the n × n rational matrix with entries a (row-major).
inline Rational determinant(std::size_t n, const std::vector<Rational>& a) {
  if (a.size() != n * n) throw std::invalid_argument("morphe::algebra::determinant: n * n rational entries");
  Matrix B(n, n);
  Integer scale = 1;
  for (std::size_t i = 0; i < n; i++) {
    Integer d = 1;
    for (std::size_t j = 0; j < n; j++) {
      const Integer& q = a[i * n + j].denominator();
      d = d / arithmetic::gcd(d, q) * q;
    }
    for (std::size_t j = 0; j < n; j++) {
      const Rational& x = a[i * n + j];
      B(i, j) = x.numerator() * (d / x.denominator());
    }
    scale *= d;
  }
  return Rational(determinant(std::move(B)), scale);
}

}  // namespace morphe::algebra
