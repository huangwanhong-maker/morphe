// The cycles of a CW complex over ℤ/2.
//
// A k-cycle over ℤ/2 is a set of k-cells whose boundaries cancel modulo 2: an element of the kernel
// of ∂_k: C_k(X; ℤ/2) → C_{k−1}(X; ℤ/2). The kernel is a vector space over the field ℤ/2, and
// Gaussian elimination of the columns of ∂_k, with the combinations that reduce a column to zero
// tracked, gives a basis of it (the reduction of persistence.hpp without a filtration; Edelsbrunner
// and Harer 2010, §IV.2). For a 2-complex, Z_2(X; ℤ/2) = H_2(X; ℤ/2), there being no 3-cells:
// its elements are the closed surfaces of the complex modulo 2, the configurations of a ℤ₂ gauge
// theory's character expansion (Drouffe and Zuber 1983, §2). The support says which cells lie in
// some cycle: a cell outside it is in none.
//
// Over ℤ/2 the projective plane has a 2-cycle (its single 2-cell, attached by a a), though over ℤ it
// has none; the dimension of Z_k(X; ℤ/2) is the Betti number over ℤ/2 when there are no (k+1)-cells.

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <morphe/topology/cw_complex.hpp>
#include <morphe/topology/persistence.hpp>

namespace morphe::topology {

struct Cycles {
  std::size_t dimension = 0;                    // dim Z_k(X; ℤ/2)
  std::vector<std::vector<std::size_t>> basis;  // each a set of k-cells, sorted
  std::vector<bool> support;                    // k-cell -> in some cycle
};

inline Cycles cycles_mod2(const CWComplex& X, std::size_t k) {
  if (k == 0) throw std::invalid_argument("morphe::topology::cycles_mod2: k >= 1");
  const std::size_t rows = X.cells(k - 1), cols = X.cells(k);
  // Each dimension's cells at their own indices.
  std::vector<std::vector<std::size_t>> position(k);
  for (std::size_t d = 0; d < k; d++) {
    position[d].resize(X.cells(d));
    for (std::size_t i = 0; i < position[d].size(); i++) position[d][i] = i;
  }
  using Bits = std::vector<std::uint64_t>;
  auto flip = [](Bits& b, std::size_t i) { b[i / 64] ^= std::uint64_t{1} << (i % 64); };
  auto low = [](const Bits& b) -> long long {
    for (std::size_t w = b.size(); w-- > 0;)
      if (b[w]) return static_cast<long long>(64 * w + 63 - static_cast<std::size_t>(__builtin_clzll(b[w])));
    return -1;
  };
  std::vector<Bits> column(cols, Bits((rows + 63) / 64, 0)), combination(cols, Bits((cols + 63) / 64, 0));
  std::vector<long long> owner(rows, -1);
  Cycles C;
  C.support.assign(cols, false);
  for (std::size_t j = 0; j < cols; j++) {
    for (std::size_t i : detail::boundary_mod2(X, k, j, position)) flip(column[j], i);
    flip(combination[j], j);
    for (long long l = low(column[j]); l >= 0 && owner[static_cast<std::size_t>(l)] >= 0; l = low(column[j])) {
      const auto o = static_cast<std::size_t>(owner[static_cast<std::size_t>(l)]);
      for (std::size_t w = 0; w < column[j].size(); w++) column[j][w] ^= column[o][w];
      for (std::size_t w = 0; w < combination[j].size(); w++) combination[j][w] ^= combination[o][w];
    }
    const long long l = low(column[j]);
    if (l >= 0) {
      owner[static_cast<std::size_t>(l)] = static_cast<long long>(j);
      continue;
    }
    std::vector<std::size_t> cycle;
    for (std::size_t c = 0; c < cols; c++)
      if ((combination[j][c / 64] >> (c % 64)) & 1) {
        cycle.push_back(c);
        C.support[c] = true;
      }
    C.basis.push_back(std::move(cycle));
  }
  C.dimension = C.basis.size();
  return C;
}

}  // namespace morphe::topology
