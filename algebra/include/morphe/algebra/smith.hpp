// The Smith normal form of an integer matrix, and ranks over ℚ and over 𝔽_p.
//
// Theorem (Smith 1861; Newman 1972, Integral Matrices, Theorem II.9). For every m × n
// integer matrix A there are unimodular U (m × m) and V (n × n) with U A V = D, D diagonal
// with entries d_1 | d_2 | ... | d_r > 0 followed by zeros; the d_i, the invariant factors,
// are determined by A (d_1 ... d_k is the gcd of the k × k minors). r is the rank of A over ℚ.
// The abelian group presented by A, ℤ^m / A ℤ^n, is then ⊕ ℤ/d_i ⊕ ℤ^{m - r}.
//
// The algorithm (Cohen 1993, A Course in Computational Algebraic Number Theory, Algorithm
// 2.4.14, in its elementary form): bring an entry of least absolute value to the pivot,
// clear its row and column by division with remainder, bring a smaller remainder to the
// pivot until the row and column are clear, then ensure the pivot divides the rest of the
// submatrix by adding a row that it does not divide. Every pass strictly decreases the
// pivot's absolute value, so it terminates on every input. Its cost is O(min(m, n) m n)
// operations on entries whose size, on adversarial inputs, can grow; the inputs here
// (residues of boundary matrices after unit pivots, relation matrices of presentations) are
// small. The transforms U, V and their inverses are accumulated on request, each
// elementary operation applied to them as it is to A.
//
// rank_mod computes the rank over 𝔽_p by Gaussian elimination on word-sized residues, O(m n
// min(m, n)).

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/arithmetic/modular.hpp>

#include "matrix.hpp"

namespace morphe::algebra {

struct Smith {
  Matrix D;                       // U A V
  std::vector<Integer> diagonal;  // d_1 | d_2 | ... | d_r, all positive
  Matrix U, V, U_inverse, V_inverse;  // filled when requested; otherwise empty
  std::size_t rank() const { return diagonal.size(); }
};

struct SmithOptions {
  bool U = false, V = false, U_inverse = false, V_inverse = false;
};

inline Smith smith(const Matrix& A, SmithOptions want = {}) {
  Smith s;
  Matrix& D = s.D;
  D = A;
  const std::size_t m = A.rows, n = A.cols;
  if (want.U) s.U = Matrix::identity(m);
  if (want.U_inverse) s.U_inverse = Matrix::identity(m);
  if (want.V) s.V = Matrix::identity(n);
  if (want.V_inverse) s.V_inverse = Matrix::identity(n);

  // Elementary operations on D, mirrored on the transforms. U A V = D throughout: a row
  // operation E on D is U <- E U and U^-1 <- U^-1 E^-1; a column operation F is V <- V F and
  // V^-1 <- F^-1 V^-1.
  auto swap_rows = [&](std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t c = 0; c < n; c++) std::swap(D(i, c), D(j, c));
    if (want.U)
      for (std::size_t c = 0; c < m; c++) std::swap(s.U(i, c), s.U(j, c));
    if (want.U_inverse)
      for (std::size_t r = 0; r < m; r++) std::swap(s.U_inverse(r, i), s.U_inverse(r, j));
  };
  auto swap_cols = [&](std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t r = 0; r < m; r++) std::swap(D(r, i), D(r, j));
    if (want.V)
      for (std::size_t r = 0; r < n; r++) std::swap(s.V(r, i), s.V(r, j));
    if (want.V_inverse)
      for (std::size_t c = 0; c < n; c++) std::swap(s.V_inverse(i, c), s.V_inverse(j, c));
  };
  // row_i += c row_j
  auto add_row = [&](std::size_t i, std::size_t j, const Integer& c) {
    if (c.is_zero()) return;
    for (std::size_t k = 0; k < n; k++)
      if (!D(j, k).is_zero()) D(i, k) += c * D(j, k);
    if (want.U)
      for (std::size_t k = 0; k < m; k++)
        if (!s.U(j, k).is_zero()) s.U(i, k) += c * s.U(j, k);
    if (want.U_inverse)
      for (std::size_t r = 0; r < m; r++)
        if (!s.U_inverse(r, i).is_zero()) s.U_inverse(r, j) -= c * s.U_inverse(r, i);
  };
  // col_i += c col_j
  auto add_col = [&](std::size_t i, std::size_t j, const Integer& c) {
    if (c.is_zero()) return;
    for (std::size_t r = 0; r < m; r++)
      if (!D(r, j).is_zero()) D(r, i) += c * D(r, j);
    if (want.V)
      for (std::size_t r = 0; r < n; r++)
        if (!s.V(r, j).is_zero()) s.V(r, i) += c * s.V(r, j);
    if (want.V_inverse)
      for (std::size_t k = 0; k < n; k++)
        if (!s.V_inverse(i, k).is_zero()) s.V_inverse(j, k) -= c * s.V_inverse(i, k);
  };
  auto negate_row = [&](std::size_t i) {
    for (std::size_t k = 0; k < n; k++) D(i, k) = -D(i, k);
    if (want.U)
      for (std::size_t k = 0; k < m; k++) s.U(i, k) = -s.U(i, k);
    if (want.U_inverse)
      for (std::size_t r = 0; r < m; r++) s.U_inverse(r, i) = -s.U_inverse(r, i);
  };

  for (std::size_t t = 0; t < std::min(m, n); t++) {
    // The entry of least absolute value in the remaining submatrix.
    std::optional<std::pair<std::size_t, std::size_t>> best;
    for (std::size_t i = t; i < m; i++)
      for (std::size_t j = t; j < n; j++)
        if (!D(i, j).is_zero() && (!best || arithmetic::abs(D(i, j)) < arithmetic::abs(D(best->first, best->second)))) best = {i, j};
    if (!best) break;
    swap_rows(t, best->first);
    swap_cols(t, best->second);
    for (;;) {
      bool clear = true;
      for (std::size_t i = t + 1; i < m; i++)
        if (!D(i, t).is_zero()) {
          add_row(i, t, -(D(i, t) / D(t, t)));
          if (!D(i, t).is_zero()) clear = false;
        }
      for (std::size_t j = t + 1; j < n; j++)
        if (!D(t, j).is_zero()) {
          add_col(j, t, -(D(t, j) / D(t, t)));
          if (!D(t, j).is_zero()) clear = false;
        }
      if (!clear) {
        // A remainder smaller than the pivot: bring the least one to the pivot and repeat.
        std::size_t bi = t, bj = t;
        for (std::size_t i = t + 1; i < m; i++)
          if (!D(i, t).is_zero() && arithmetic::abs(D(i, t)) < arithmetic::abs(D(bi, bj))) bi = i, bj = t;
        for (std::size_t j = t + 1; j < n; j++)
          if (!D(t, j).is_zero() && arithmetic::abs(D(t, j)) < arithmetic::abs(D(bi, bj))) bi = t, bj = j;
        swap_rows(t, bi);
        swap_cols(t, bj);
        continue;
      }
      // The pivot must divide every remaining entry; if not, add that entry's row to the
      // pivot's row, which puts a non-multiple into the pivot's row, and repeat.
      std::optional<std::size_t> bad;
      for (std::size_t i = t + 1; i < m && !bad; i++)
        for (std::size_t j = t + 1; j < n; j++)
          if (!(D(i, j) % D(t, t)).is_zero()) { bad = i; break; }
      if (!bad) break;
      add_row(t, *bad, 1);
    }
    if (D(t, t).sign() < 0) negate_row(t);
    s.diagonal.push_back(D(t, t));
  }
  return s;
}

// The invariant factors greater than 1 of a matrix: the torsion of the group it presents.
inline std::vector<Integer> torsion(const Smith& s) {
  std::vector<Integer> t;
  for (const auto& d : s.diagonal)
    if (!d.is_one()) t.push_back(d);
  return t;
}

// The rank over 𝔽_p (p prime, below 2^63) or, for p = 0, over ℚ.
inline std::size_t rank_mod(const Matrix& A, std::uint64_t p) {
  if (p == 0) return smith(A).rank();
  using arithmetic::Residue;
  std::vector<Residue> a(A.rows * A.cols);
  for (std::size_t i = 0; i < A.rows * A.cols; i++) a[i] = arithmetic::residue(A.a[i], p);
  std::size_t r = 0;
  for (std::size_t c = 0; c < A.cols && r < A.rows; c++) {
    std::size_t piv = r;
    while (piv < A.rows && a[piv * A.cols + c] == 0) piv++;
    if (piv == A.rows) continue;
    for (std::size_t k = 0; k < A.cols; k++) std::swap(a[piv * A.cols + k], a[r * A.cols + k]);
    const Residue inv = arithmetic::mod_inverse(a[r * A.cols + c], p);
    for (std::size_t i = r + 1; i < A.rows; i++) {
      const Residue f = arithmetic::mod_mul(a[i * A.cols + c], inv, p);
      if (!f) continue;
      for (std::size_t k = c; k < A.cols; k++)
        a[i * A.cols + k] = arithmetic::mod_sub(a[i * A.cols + k], arithmetic::mod_mul(f, a[r * A.cols + k], p), p);
    }
    r++;
  }
  return r;
}

}  // namespace morphe::algebra
