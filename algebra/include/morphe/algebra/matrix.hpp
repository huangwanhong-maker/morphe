// Exact matrices over ℤ, dense and sparse.
//
// A Matrix is dense and row-major, of any shape, for the small matrices that remain after a
// sparse reduction and for presentations of modules. A SparseMatrix is held by columns,
// each a list of (row, value) sorted by row with no zero values: the boundary matrices of
// complexes, whose columns (the boundaries of cells) have few entries. Entries are morphe's
// Integers, so nothing is rounded and nothing overflows.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/arithmetic/integer.hpp>

namespace morphe::algebra {

using arithmetic::Integer;

struct Matrix {
  std::size_t rows = 0, cols = 0;
  std::vector<Integer> a;

  Matrix() = default;
  Matrix(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c) {}
  static Matrix identity(std::size_t n) {
    Matrix m(n, n);
    for (std::size_t i = 0; i < n; i++) m(i, i) = 1;
    return m;
  }
  Integer& operator()(std::size_t i, std::size_t j) { return a[i * cols + j]; }
  const Integer& operator()(std::size_t i, std::size_t j) const { return a[i * cols + j]; }
  bool is_zero() const {
    for (const auto& x : a)
      if (!x.is_zero()) return false;
    return true;
  }
  Matrix transpose() const {
    Matrix t(cols, rows);
    for (std::size_t i = 0; i < rows; i++)
      for (std::size_t j = 0; j < cols; j++) t(j, i) = (*this)(i, j);
    return t;
  }
  friend bool operator==(const Matrix&, const Matrix&) = default;
};

inline Matrix operator*(const Matrix& x, const Matrix& y) {
  if (x.cols != y.rows) throw std::invalid_argument("morphe::algebra: matrix shapes do not compose");
  Matrix z(x.rows, y.cols);
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t k = 0; k < x.cols; k++) {
      const Integer& xik = x(i, k);
      if (xik.is_zero()) continue;
      for (std::size_t j = 0; j < y.cols; j++)
        if (!y(k, j).is_zero()) z(i, j) += xik * y(k, j);
    }
  return z;
}

// A finite formal sum of basis elements with integer coefficients: a chain, a column of a
// sparse matrix, a vector of ℤⁿ. Sorted by index, no zero coefficients.
struct Entry {
  std::size_t index;
  Integer value;
  friend bool operator==(const Entry&, const Entry&) = default;
};
using Vector = std::vector<Entry>;

// x + c y, both sorted; the result sorted, without zeros.
inline Vector axpy(const Vector& x, const Integer& c, const Vector& y) {
  Vector z;
  z.reserve(x.size() + y.size());
  std::size_t i = 0, j = 0;
  while (i < x.size() || j < y.size()) {
    if (j == y.size() || (i < x.size() && x[i].index < y[j].index)) {
      z.push_back(x[i++]);
    } else if (i == x.size() || y[j].index < x[i].index) {
      z.push_back({y[j].index, c * y[j].value});
      j++;
    } else {
      Integer v = x[i].value + c * y[j].value;
      if (!v.is_zero()) z.push_back({x[i].index, std::move(v)});
      i++, j++;
    }
  }
  return z;
}

// The coefficient of index k in x, or 0.
inline Integer coefficient(const Vector& x, std::size_t k) {
  std::size_t lo = 0, hi = x.size();
  while (lo < hi) {
    const std::size_t mid = (lo + hi) / 2;
    if (x[mid].index < k) lo = mid + 1;
    else hi = mid;
  }
  return lo < x.size() && x[lo].index == k ? x[lo].value : Integer(0);
}

struct SparseMatrix {
  std::size_t rows = 0;
  std::vector<Vector> columns;

  SparseMatrix() = default;
  SparseMatrix(std::size_t r, std::size_t c) : rows(r), columns(c) {}
  std::size_t cols() const { return columns.size(); }
  // Adds value at (i, j), keeping the column sorted.
  void add(std::size_t i, std::size_t j, const Integer& value) {
    if (i >= rows || j >= columns.size()) throw std::out_of_range("morphe::algebra::SparseMatrix::add");
    columns[j] = axpy(columns[j], value, Vector{{i, 1}});
  }
  Matrix dense() const {
    Matrix m(rows, cols());
    for (std::size_t j = 0; j < cols(); j++)
      for (const auto& e : columns[j]) m(e.index, j) = e.value;
    return m;
  }
  static SparseMatrix from(const Matrix& m) {
    SparseMatrix s(m.rows, m.cols);
    for (std::size_t j = 0; j < m.cols; j++)
      for (std::size_t i = 0; i < m.rows; i++)
        if (!m(i, j).is_zero()) s.columns[j].push_back({i, m(i, j)});
    return s;
  }
  // The image of a vector of ℤ^cols.
  Vector apply(const Vector& x) const {
    Vector y;
    for (const auto& e : x) y = axpy(y, e.value, columns.at(e.index));
    return y;
  }
};

// x y for sparse matrices of compatible shapes.
inline SparseMatrix operator*(const SparseMatrix& x, const SparseMatrix& y) {
  if (x.cols() != y.rows) throw std::invalid_argument("morphe::algebra: sparse matrix shapes do not compose");
  SparseMatrix z(x.rows, y.cols());
  for (std::size_t j = 0; j < y.cols(); j++) z.columns[j] = x.apply(y.columns[j]);
  return z;
}

}  // namespace morphe::algebra
