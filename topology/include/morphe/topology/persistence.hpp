// Persistent homology over ℤ/2 of a filtered CW complex.
//
// A filtration of a CW complex X gives each cell c a stage s(c) ∈ ℕ, every face of a cell at a
// stage no later than the cell's own: the sub-complexes X_s = {c : s(c) ≤ s} increase with s.
// Their homology over the field ℤ/2, with the maps the inclusions induce, is a persistence
// module of finite type, and it decomposes uniquely into intervals [b, d): a class born at
// stage b that dies (becomes a boundary, or merges with an older class) at stage d, or lives
// for ever (d = ∞) (Zomorodian and Carlsson 2005, Theorem 2.1; Edelsbrunner and Harer 2010,
// §VII.1). The stages are exact integers, not real numbers.
//
// The algorithm is the standard reduction of the boundary matrix with its columns in the order
// of the filtration: cells by stage, then by dimension, then by index, a total order in which
// every face precedes its cofaces. A column is reduced by adding earlier columns while its
// lowest entry collides with theirs; a reduced column j with lowest entry i pairs the cell
// i, which created a class, with the cell j, which killed it; an unpaired creator is a class
// that lives for ever (Edelsbrunner, Letscher and Zomorodian 2002; Edelsbrunner and Harer
// 2010, §VII.1, the algorithm "pair"). It always terminates, in O(n³) for n cells at worst and
// near linearly for the sparse complexes met in practice. The essential intervals of the last
// stage count the Betti numbers of X over ℤ/2.
//
// Coefficients are ℤ/2: the boundary of an edge is its two ends (nothing, for a loop), of a
// 2-cell its attaching word with each edge counted modulo 2, of a higher cell its integer
// boundary modulo 2.

#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/topology/cw_complex.hpp>

namespace morphe::topology {

// The stage of every cell, by dimension: vertex[i], edge[i], face[i], higher[k - 3][i].
struct Filtration {
  std::vector<std::size_t> vertex, edge, face;
  std::vector<std::vector<std::size_t>> higher;

  std::size_t stage(std::size_t dim, std::size_t i) const {
    if (dim == 0) return vertex.at(i);
    if (dim == 1) return edge.at(i);
    if (dim == 2) return face.at(i);
    return higher.at(dim - 3).at(i);
  }
};

// A class of H_k born when `creator` enters, at stage `birth`, and killed when `killer`
// enters, at stage `death`; no killer is a class that lives for ever.
struct PersistenceInterval {
  std::size_t dimension = 0, birth = 0, creator = 0;
  std::optional<std::size_t> death, killer;
  bool essential() const { return !death.has_value(); }
};

struct PersistenceDiagram {
  std::vector<PersistenceInterval> intervals;  // by dimension, then birth, then death

  // The Betti number of X_s over ℤ/2: the intervals [b, d) with b ≤ s < d.
  std::size_t betti(std::size_t dim, std::size_t s) const {
    std::size_t n = 0;
    for (const auto& I : intervals)
      if (I.dimension == dim && I.birth <= s && (!I.death || s < *I.death)) n++;
    return n;
  }
};

namespace detail {

// The ℤ/2 boundary of a cell, as the sorted indices of its faces in the cell list.
inline std::vector<std::size_t> boundary_mod2(const CWComplex& X, std::size_t dim, std::size_t i, const std::vector<std::vector<std::size_t>>& position) {
  std::vector<std::size_t> b;
  auto toggle = [&](std::size_t p) {
    const auto it = std::lower_bound(b.begin(), b.end(), p);
    if (it != b.end() && *it == p) b.erase(it);
    else b.insert(it, p);
  };
  if (dim == 1) {
    const auto& e = X.edges[i];
    if (e[0] != e[1]) { toggle(position[0][e[0]]); toggle(position[0][e[1]]); }
  } else if (dim == 2) {
    for (const auto& s : X.faces[i].boundary) toggle(position[1][s.edge]);
  } else if (dim >= 3) {
    for (const auto& entry : X.higher[dim - 3][i])
      if (entry.value.to_int64() % 2 != 0) toggle(position[dim - 1][entry.index]);
  }
  return b;
}

}  // namespace detail

// The persistence diagram of a filtered CW complex. Throws std::invalid_argument when the
// stages are not a filtration (a cell earlier than one of its faces) or do not fit X.
// Intervals of length zero (born and killed at one stage) are kept only when asked.
inline PersistenceDiagram persistent_homology(const CWComplex& X, const Filtration& F, bool zero_length = false) {
  const std::size_t top = X.dimension();
  std::vector<std::size_t> count(top + 1);
  for (std::size_t d = 0; d <= top; d++) count[d] = X.cells(d);
  if (F.vertex.size() != count[0] || (top >= 1 && F.edge.size() != count[1]) || (top >= 2 && F.face.size() != count[2]))
    throw std::invalid_argument("a filtration gives one stage to every cell");
  for (std::size_t d = 3; d <= top; d++)
    if (F.higher.size() < d - 2 || F.higher[d - 3].size() != count[d]) throw std::invalid_argument("a filtration gives one stage to every cell");

  // The cells in the order of the filtration.
  struct Cell { std::size_t stage, dim, index; };
  std::vector<Cell> cells;
  for (std::size_t d = 0; d <= top; d++)
    for (std::size_t i = 0; i < count[d]; i++) cells.push_back({F.stage(d, i), d, i});
  std::stable_sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) {
    return a.stage != b.stage ? a.stage < b.stage : a.dim != b.dim ? a.dim < b.dim : a.index < b.index;
  });
  std::vector<std::vector<std::size_t>> position(top + 1);
  for (std::size_t d = 0; d <= top; d++) position[d].assign(count[d], 0);
  for (std::size_t p = 0; p < cells.size(); p++) position[cells[p].dim][cells[p].index] = p;

  // Reduce.
  std::vector<std::vector<std::size_t>> column(cells.size());
  std::vector<std::size_t> owner(cells.size(), cells.size());  // lowest entry -> its column
  std::vector<bool> paired(cells.size(), false);
  PersistenceDiagram D;
  for (std::size_t j = 0; j < cells.size(); j++) {
    auto& c = column[j];
    c = detail::boundary_mod2(X, cells[j].dim, cells[j].index, position);
    for (std::size_t p : c)
      if (cells[p].stage > cells[j].stage)
        throw std::invalid_argument("a cell enters before one of its faces: the stages are not a filtration");
    while (!c.empty() && owner[c.back()] != cells.size()) {
      std::vector<std::size_t> sum;
      const auto& other = column[owner[c.back()]];
      std::set_symmetric_difference(c.begin(), c.end(), other.begin(), other.end(), std::back_inserter(sum));
      c.swap(sum);
    }
    if (!c.empty()) {
      const std::size_t i = c.back();
      owner[i] = j;
      paired[i] = paired[j] = true;
      if (zero_length || cells[i].stage < cells[j].stage)
        D.intervals.push_back({cells[i].dim, cells[i].stage, cells[i].index, cells[j].stage, cells[j].index});
    }
  }
  for (std::size_t p = 0; p < cells.size(); p++)
    if (!paired[p] && column[p].empty()) D.intervals.push_back({cells[p].dim, cells[p].stage, cells[p].index, std::nullopt, std::nullopt});
  std::stable_sort(D.intervals.begin(), D.intervals.end(), [](const PersistenceInterval& a, const PersistenceInterval& b) {
    if (a.dimension != b.dimension) return a.dimension < b.dimension;
    if (a.birth != b.birth) return a.birth < b.birth;
    const std::size_t da = a.death.value_or(static_cast<std::size_t>(-1)), db = b.death.value_or(static_cast<std::size_t>(-1));
    return da < db;
  });
  return D;
}

}  // namespace morphe::topology
