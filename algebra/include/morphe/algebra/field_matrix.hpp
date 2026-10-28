// Dense matrices over an exact field: products, traces, ranks, determinants and inverses.
//
// matrix.hpp holds matrices over ℤ, for the Smith normal form; determinant.hpp the determinant
// over ℚ by clearing denominators into Bareiss's elimination. Representations and group
// algebras need the same over fields that are not ℚ, the cyclotomic fields above all, so this
// header writes Gaussian elimination once for any exact field K: a type with K(0), K(1),
// + − × /, == and is_zero(), such as arithmetic::Rational and arithmetic::Cyclotomic.
//
// Gaussian elimination with the first non-zero pivot (in exact arithmetic the choice of pivot
// cannot change the result, only the sizes of the entries on the way) reduces an m × n matrix of
// rank r to row echelon form in O(m n r) field operations; the determinant is
// the product of the pivots with the sign of the row exchanges, and the inverse is read from
// the reduced form of [A | I] (Gauss–Jordan), O(n³). Every algorithm terminates on every input.
// Over ℚ the sizes of the entries can grow during elimination (Bareiss's method bounds them
// over ℤ, determinant.hpp); for the matrices of representations and idempotents met here they
// stay small.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::algebra {

template <class K>
struct FieldMatrix {
  std::size_t rows = 0, cols = 0;
  std::vector<K> a;  // row-major

  FieldMatrix() = default;
  FieldMatrix(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c, K(0)) {}
  static FieldMatrix identity(std::size_t n) {
    FieldMatrix m(n, n);
    for (std::size_t i = 0; i < n; i++) m(i, i) = K(1);
    return m;
  }
  K& operator()(std::size_t i, std::size_t j) { return a[i * cols + j]; }
  const K& operator()(std::size_t i, std::size_t j) const { return a[i * cols + j]; }

  bool is_zero() const {
    for (const auto& x : a)
      if (!x.is_zero()) return false;
    return true;
  }
  FieldMatrix transpose() const {
    FieldMatrix t(cols, rows);
    for (std::size_t i = 0; i < rows; i++)
      for (std::size_t j = 0; j < cols; j++) t(j, i) = (*this)(i, j);
    return t;
  }
  K trace() const {
    if (rows != cols) throw std::invalid_argument("morphe::algebra::FieldMatrix::trace: a square matrix");
    K t(0);
    for (std::size_t i = 0; i < rows; i++) t += (*this)(i, i);
    return t;
  }

  // Row echelon form in place; returns the pivot columns. det, when given, receives the
  // determinant of a square matrix (the product of the pivots, with the sign of the exchanges).
  std::vector<std::size_t> echelon(K* det = nullptr) {
    std::vector<std::size_t> pivots;
    bool negate = false;
    std::size_t rank = 0;
    for (std::size_t c = 0; c < cols && rank < rows; c++) {
      std::size_t p = rank;
      while (p < rows && (*this)(p, c).is_zero()) p++;
      if (p == rows) continue;
      if (p != rank) {
        for (std::size_t j = c; j < cols; j++) std::swap((*this)(p, j), (*this)(rank, j));
        negate = !negate;
      }
      const K inv = K(1) / (*this)(rank, c);
      for (std::size_t i = rank + 1; i < rows; i++) {
        if ((*this)(i, c).is_zero()) continue;
        const K f = (*this)(i, c) * inv;
        for (std::size_t j = c; j < cols; j++)
          if (!(*this)(rank, j).is_zero()) (*this)(i, j) -= f * (*this)(rank, j);
      }
      pivots.push_back(c);
      rank++;
    }
    if (det) {
      if (rows != cols) throw std::invalid_argument("morphe::algebra::FieldMatrix: the determinant of a square matrix");
      K d(pivots.size() == rows ? 1 : 0);
      if (pivots.size() == rows)
        for (std::size_t i = 0; i < rows; i++) d *= (*this)(i, i);
      *det = negate ? -d : d;
    }
    return pivots;
  }
  std::size_t rank() const {
    FieldMatrix m = *this;
    return m.echelon().size();
  }
  K determinant() const {
    FieldMatrix m = *this;
    K d(0);
    m.echelon(&d);
    return d;
  }
  // The inverse; throws std::domain_error when the matrix is singular.
  FieldMatrix inverse() const {
    if (rows != cols) throw std::invalid_argument("morphe::algebra::FieldMatrix::inverse: a square matrix");
    const std::size_t n = rows;
    FieldMatrix w(n, 2 * n);
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < n; j++) w(i, j) = (*this)(i, j);
      w(i, n + i) = K(1);
    }
    const auto pivots = w.echelon();
    if (pivots.size() < n || pivots[n - 1] >= n) throw std::domain_error("morphe::algebra::FieldMatrix::inverse: a singular matrix");
    // Back substitution: scale each pivot row to 1 and clear the column above it.
    for (std::size_t k = n; k-- > 0;) {
      const K inv = K(1) / w(k, k);
      for (std::size_t j = k; j < 2 * n; j++) w(k, j) *= inv;
      for (std::size_t i = 0; i < k; i++) {
        if (w(i, k).is_zero()) continue;
        const K f = w(i, k);
        for (std::size_t j = k; j < 2 * n; j++) w(i, j) -= f * w(k, j);
      }
    }
    FieldMatrix inv(n, n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) inv(i, j) = w(i, n + j);
    return inv;
  }

  friend bool operator==(const FieldMatrix&, const FieldMatrix&) = default;
};

template <class K>
FieldMatrix<K> operator*(const FieldMatrix<K>& x, const FieldMatrix<K>& y) {
  if (x.cols != y.rows) throw std::invalid_argument("morphe::algebra::FieldMatrix: shapes do not compose");
  FieldMatrix<K> z(x.rows, y.cols);
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t k = 0; k < x.cols; k++) {
      const K& xik = x(i, k);
      if (xik.is_zero()) continue;
      for (std::size_t j = 0; j < y.cols; j++)
        if (!y(k, j).is_zero()) z(i, j) += xik * y(k, j);
    }
  return z;
}
template <class K>
FieldMatrix<K> operator+(FieldMatrix<K> x, const FieldMatrix<K>& y) {
  if (x.rows != y.rows || x.cols != y.cols) throw std::invalid_argument("morphe::algebra::FieldMatrix: shapes differ");
  for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] += y.a[i];
  return x;
}
template <class K>
FieldMatrix<K> operator-(FieldMatrix<K> x, const FieldMatrix<K>& y) {
  if (x.rows != y.rows || x.cols != y.cols) throw std::invalid_argument("morphe::algebra::FieldMatrix: shapes differ");
  for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] -= y.a[i];
  return x;
}
template <class K>
FieldMatrix<K> operator*(const K& c, FieldMatrix<K> x) {
  for (auto& v : x.a) v *= c;
  return x;
}

}  // namespace morphe::algebra
