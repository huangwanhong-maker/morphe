// Regular CW complexes held by their face posets: each cell's dimension, its type, its facets with their
// incidence numbers, and an identity; whether a poset so given is the face poset of a regular CW complex,
// decided exactly in dimension at most 3; orientation; the cellular chain complex; the complex as a
// CWComplex (cw_complex.hpp), so that what is computed there serves it; and cubical grids.
//
// Definition (Lundell and Weingram 1969, The Topology of CW Complexes, ch. III; Björner 1984, Eur. J.
// Combin. 5, 7–16). A CW complex is regular when the characteristic map of each cell is a homeomorphism of
// the closed ball onto the cell's closure. Its cells, ordered by σ ≤ τ when σ ⊆ cl τ, form its face poset,
// which determines the complex up to cellular homeomorphism: the order complex of the poset is its
// barycentric subdivision. A finite poset is the face poset of a regular CW complex exactly when it is a CW
// poset: for each element x the order complex of {y < x} is homeomorphic to the sphere of dimension
// rank(x) − 1 (Björner 1984, Def. 2.1 and Prop. 3.1; the empty set is the sphere of dimension −1). Here the
// poset is held by its covering relation, each cell's facets, so a cell's faces are its facets, their
// facets, and so on; a cell's type says what it is in a model (a junction, a street, a block, a room).
//
// The test. Cells are tested in order of dimension, so that when a k-cell is tested every cell of its
// boundary B (the closure of its facets) is already known to be a regular cell, and B is a regular complex
// of dimension k − 1. Then B is a sphere exactly when
//
//   k = 1: B is two distinct vertices;
//   k = 2: B is a circle: connected, each of its vertices on exactly two of its edges;
//   k = 3: B is a 2-sphere: each of its edges on exactly two of its 2-cells and the link of each of its
//          vertices one cycle (the edges at the vertex joined by the 2-cells at it), which make B a closed
//          surface (Massey 1991, A Basic Course in Algebraic Topology, ch. 1); connected; and of Euler
//          characteristic 2, the sphere alone among closed connected surfaces;
//   k ≥ 4: recognising the 3-sphere is decidable (Rubinstein 1995; Thompson 1994, Math. Res. Lett. 1, 613)
//          but not attempted here. B is required to be a pseudomanifold with the homology over ℤ/2 of
//          S^{k−1}, which is necessary, and the answer is marked inexact.
//
// Orientation (Massey 1991, ch. IX). The incidence numbers [σ : τ] of a regular complex can be chosen in
// {+1, −1} so that ∂∂ = 0, and for each cell they are then unique up to one sign: ∂σ must be a fundamental
// cycle of the sphere ∂σ, which has a coefficient ±1 on every facet. A 1-cell is oriented from its facet
// of sign −1 (its tail) to its facet of sign +1 (its head); given no signs, from its first facet to its
// second. For k ≥ 2, two facets τ, τ′ of σ sharing a (k−2)-cell ρ — in a sphere every such ρ lies on
// exactly two facets, the diamond property — must satisfy [σ:τ][τ:ρ] + [σ:τ′][τ′:ρ] = 0; orient()
// propagates that from a facet whose sign is given (or the first, set to +1) and refuses a conflict.
//
// Identities: each cell carries an id, unique within the complex and kept by rewriting
// (complex_rewriting.hpp), so that a cell can be followed through a history while its index changes.
//
// Cost: the test of a k-cell, k ≤ 3, is linear in the size of its closure, so the whole test is
// O(Σ_σ |cl σ|); orientation likewise. Checks: tests/check_regular_complex.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/algebra/chain_complex.hpp>

#include "cw_complex.hpp"

namespace morphe::topology {

inline constexpr std::size_t no_cell = std::numeric_limits<std::size_t>::max();

struct RegularCell {
  std::size_t dim = 0;
  std::string type;                  // what the cell is in a model; any text without '|'
  std::vector<std::size_t> facets;   // indices of (dim − 1)-cells
  std::vector<int> signs;            // [σ : τ] for each facet, +1 or −1; 0 = not yet oriented
  std::uint64_t id = 0;              // an identity, unique in the complex
};

struct Admissibility {
  bool regular = false;   // the poset is the face poset of a regular CW complex
  bool exact = true;      // false when a cell of dimension ≥ 4 was tested by homology alone
  std::size_t cell = no_cell;  // the first cell that fails
  std::string reason;
};

struct RegularComplex {
  std::string name;
  std::vector<RegularCell> cells;
  std::uint64_t next_id = 1;

  // A cell of dimension dim with the given facets (and signs, or none to be oriented); returns its index.
  std::size_t add(std::size_t dim, std::string type, std::vector<std::size_t> facets, std::vector<int> signs = {}, std::uint64_t id = 0) {
    if (type.find('|') != std::string::npos) throw std::invalid_argument("morphe::topology::RegularComplex: a type may not contain '|'");
    if (signs.empty()) signs.assign(facets.size(), 0);
    if (signs.size() != facets.size()) throw std::invalid_argument("morphe::topology::RegularComplex: one sign for each facet");
    for (std::size_t f : facets)
      if (f >= cells.size()) throw std::invalid_argument("morphe::topology::RegularComplex: a facet is not a cell yet");
    if (id == 0) id = next_id++;
    else next_id = std::max(next_id, id + 1);
    cells.push_back({dim, std::move(type), std::move(facets), std::move(signs), id});
    return cells.size() - 1;
  }
  std::size_t add_vertex(std::string type = "") { return add(0, std::move(type), {}); }
  // A 1-cell from tail to head: ∂e = head − tail.
  std::size_t add_edge(std::size_t tail, std::size_t head, std::string type = "") { return add(1, std::move(type), {tail, head}, {-1, 1}); }

  std::size_t size() const { return cells.size(); }
  std::size_t dimension() const {
    std::size_t d = 0;
    for (const auto& c : cells) d = std::max(d, c.dim);
    return d;
  }
  // The number of cells of each dimension 0..dimension().
  std::vector<std::size_t> counts() const {
    std::vector<std::size_t> n(cells.empty() ? 0 : dimension() + 1, 0);
    for (const auto& c : cells) n[c.dim]++;
    return n;
  }
  long long euler_characteristic() const {
    long long chi = 0;
    for (const auto& c : cells) chi += c.dim % 2 ? -1 : 1;
    return chi;
  }
  // The indices of the k-cells, in order; and each cell's position among the cells of its dimension.
  std::vector<std::size_t> of_dimension(std::size_t k) const {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < cells.size(); i++)
      if (cells[i].dim == k) out.push_back(i);
    return out;
  }
  std::vector<std::size_t> local_index() const {
    std::vector<std::size_t> pos(cells.size()), next(dimension() + 1, 0);
    for (std::size_t i = 0; i < cells.size(); i++) pos[i] = next[cells[i].dim]++;
    return pos;
  }
  // For each cell, the cells of which it is a facet.
  std::vector<std::vector<std::size_t>> cofaces() const {
    std::vector<std::vector<std::size_t>> up(cells.size());
    for (std::size_t i = 0; i < cells.size(); i++)
      for (std::size_t f : cells[i].facets) up[f].push_back(i);
    return up;
  }
  // σ and all its faces, sorted.
  std::vector<std::size_t> closure(std::size_t s) const {
    std::vector<char> seen(cells.size(), 0);
    std::vector<std::size_t> stack{s}, out;
    seen[s] = 1;
    while (!stack.empty()) {
      const std::size_t c = stack.back();
      stack.pop_back();
      out.push_back(c);
      for (std::size_t f : cells[c].facets)
        if (!seen[f]) seen[f] = 1, stack.push_back(f);
    }
    std::sort(out.begin(), out.end());
    return out;
  }
  // The index of the cell with this id, or no_cell.
  std::size_t find(std::uint64_t id) const {
    for (std::size_t i = 0; i < cells.size(); i++)
      if (cells[i].id == id) return i;
    return no_cell;
  }
  // The tail and head of a 1-cell (its facets of sign −1 and +1).
  std::array<std::size_t, 2> ends(std::size_t e) const {
    const auto& c = cells.at(e);
    if (c.dim != 1 || c.facets.size() != 2) throw std::invalid_argument("morphe::topology::RegularComplex::ends: not a 1-cell with two facets");
    if (c.signs[0] == 1 && c.signs[1] == -1) return {c.facets[1], c.facets[0]};
    return {c.facets[0], c.facets[1]};
  }
  // The incidence number [σ : τ], 0 when τ is not a facet of σ.
  int incidence(std::size_t s, std::size_t t) const {
    const auto& c = cells.at(s);
    for (std::size_t i = 0; i < c.facets.size(); i++)
      if (c.facets[i] == t) return c.signs[i];
    return 0;
  }
};

namespace detail {

inline Admissibility refuse(std::size_t cell, std::string why) {
  Admissibility a;
  a.regular = false;
  a.cell = cell;
  a.reason = std::move(why);
  return a;
}

inline std::string named(const RegularComplex& X, std::size_t i) {
  const auto& c = X.cells[i];
  return "the " + std::to_string(c.dim) + "-cell " + (c.type.empty() ? "" : c.type + " ") + "#" + std::to_string(c.id);
}

// Whether the graph on `nodes` with the given edges (pairs of nodes) is connected.
inline bool connected(const std::vector<std::size_t>& nodes, const std::vector<std::array<std::size_t, 2>>& edges) {
  if (nodes.empty()) return true;
  std::map<std::size_t, std::size_t> parent;
  for (std::size_t v : nodes) parent[v] = v;
  std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = root(parent[x]); };
  for (const auto& [a, b] : edges) parent[root(a)] = root(b);
  const std::size_t r = root(nodes[0]);
  for (std::size_t v : nodes)
    if (root(v) != r) return false;
  return true;
}

// The rank over ℤ/2 of a 0/1 matrix with `rows` rows, given by the rows of each column's ones; Gaussian
// elimination on bit-packed columns.
inline std::size_t rank_mod2(std::size_t rows, const std::vector<std::vector<std::size_t>>& cols) {
  const std::size_t words = (rows + 63) / 64;
  std::vector<std::vector<std::uint64_t>> m;
  for (const auto& col : cols) {
    std::vector<std::uint64_t> v(words, 0);
    for (std::size_t r : col) v[r / 64] ^= std::uint64_t{1} << (r % 64);
    m.push_back(std::move(v));
  }
  std::size_t rank = 0;
  for (std::size_t r = 0; r < rows && rank < m.size(); r++) {
    const std::uint64_t bit = std::uint64_t{1} << (r % 64);
    std::size_t p = rank;
    while (p < m.size() && !(m[p][r / 64] & bit)) p++;
    if (p == m.size()) continue;
    std::swap(m[p], m[rank]);
    for (std::size_t q = 0; q < m.size(); q++)
      if (q != rank && (m[q][r / 64] & bit))
        for (std::size_t w = 0; w < words; w++) m[q][w] ^= m[rank][w];
    rank++;
  }
  return rank;
}

// The cells of a closure by dimension.
inline std::vector<std::vector<std::size_t>> by_dimension(const RegularComplex& X, const std::vector<std::size_t>& cl, std::size_t top) {
  std::vector<std::vector<std::size_t>> d(top + 1);
  for (std::size_t c : cl)
    if (X.cells[c].dim <= top) d[X.cells[c].dim].push_back(c);
  return d;
}

// Whether the facets of σ (k ≥ 2) form a pseudomanifold: every (k−2)-face of the boundary on exactly two
// facets; the facets strongly connected across those faces. Fills `pairs` with the two facets (as positions
// in σ's facet list) and the faces between them.
struct Diamond {
  std::size_t rho, a, b;  // the (k−2)-cell and the positions of its two facets of σ
};
inline bool diamonds(const RegularComplex& X, std::size_t s, std::vector<Diamond>& out, std::string& why) {
  const auto& c = X.cells[s];
  std::map<std::size_t, std::vector<std::size_t>> on;
  for (std::size_t i = 0; i < c.facets.size(); i++)
    for (std::size_t r : X.cells[c.facets[i]].facets) on[r].push_back(i);
  out.clear();
  std::vector<std::array<std::size_t, 2>> joins;
  for (const auto& [r, list] : on) {
    if (list.size() != 2) {
      why = "its boundary has a " + std::to_string(c.dim - 2) + "-cell (#" + std::to_string(X.cells[r].id) + ") on " + std::to_string(list.size()) + " of its facets, not two";
      return false;
    }
    out.push_back({r, list[0], list[1]});
    joins.push_back({list[0], list[1]});
  }
  std::vector<std::size_t> positions(c.facets.size());
  for (std::size_t i = 0; i < positions.size(); i++) positions[i] = i;
  if (!connected(positions, joins)) {
    why = "its facets fall into separate pieces";
    return false;
  }
  return true;
}

}  // namespace detail

// Whether X is a regular CW complex (the face poset of one), and if not, which cell fails and why. Facets
// must be of one dimension less, distinct, with signs ±1 or 0.
inline Admissibility admissible(const RegularComplex& X) {
  Admissibility ok;
  ok.regular = true;
  std::vector<std::size_t> order(X.size());
  for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return X.cells[a].dim < X.cells[b].dim; });
  std::set<std::uint64_t> ids;
  for (std::size_t s : order) {
    const auto& c = X.cells[s];
    const std::string me = detail::named(X, s);
    if (!ids.insert(c.id).second) return detail::refuse(s, me + " repeats an id");
    if (c.signs.size() != c.facets.size()) return detail::refuse(s, me + " has not one sign for each facet");
    std::set<std::size_t> distinct(c.facets.begin(), c.facets.end());
    if (distinct.size() != c.facets.size()) return detail::refuse(s, me + " names a facet twice");
    for (std::size_t i = 0; i < c.facets.size(); i++) {
      const std::size_t f = c.facets[i];
      if (f >= X.size()) return detail::refuse(s, me + " names a facet that is not a cell");
      if (X.cells[f].dim + 1 != c.dim) return detail::refuse(s, me + " has a facet of dimension " + std::to_string(X.cells[f].dim));
      if (c.signs[i] != 0 && c.signs[i] != 1 && c.signs[i] != -1) return detail::refuse(s, me + " has an incidence other than ±1");
    }
    if (c.dim == 0) {
      if (!c.facets.empty()) return detail::refuse(s, me + " is a vertex with facets");
      continue;
    }
    if (c.dim == 1) {
      if (c.facets.size() != 2) return detail::refuse(s, me + " does not have two distinct ends");
      if (c.signs[0] * c.signs[1] > 0) return detail::refuse(s, me + " has ends of one sign");
      continue;
    }
    const auto cl = X.closure(s);
    if (c.dim == 2) {
      // The boundary: the facets (edges) and their ends; each end on exactly two of the edges; connected.
      std::map<std::size_t, std::size_t> degree;
      std::vector<std::array<std::size_t, 2>> edges;
      for (std::size_t e : c.facets) {
        const auto& ec = X.cells[e];
        degree[ec.facets[0]]++, degree[ec.facets[1]]++;
        edges.push_back({ec.facets[0], ec.facets[1]});
      }
      std::vector<std::size_t> nodes;
      for (const auto& [v, d] : degree) {
        if (d != 2) return detail::refuse(s, me + "'s boundary is not a circle: vertex #" + std::to_string(X.cells[v].id) + " is on " + std::to_string(d) + " of its edges");
        nodes.push_back(v);
      }
      if (!detail::connected(nodes, edges)) return detail::refuse(s, me + "'s boundary is not a circle: it is not connected");
      continue;
    }
    if (c.dim == 3) {
      const auto parts = detail::by_dimension(X, cl, 2);
      const std::set<std::size_t> faces(c.facets.begin(), c.facets.end());
      // Each edge of the boundary on exactly two of its 2-cells.
      std::map<std::size_t, std::vector<std::size_t>> faces_on;
      for (std::size_t f : c.facets)
        for (std::size_t e : X.cells[f].facets) faces_on[e].push_back(f);
      for (std::size_t e : parts[1])
        if (faces_on[e].size() != 2) return detail::refuse(s, me + "'s boundary is not a closed surface: edge #" + std::to_string(X.cells[e].id) + " is on " + std::to_string(faces_on[e].size()) + " of its 2-cells");
      // The link of each vertex one cycle: the edges at v, joined by each 2-cell at v (which meets v in two edges).
      for (std::size_t v : parts[0]) {
        std::vector<std::size_t> at;
        for (std::size_t e : parts[1]) {
          const auto& ec = X.cells[e];
          if (ec.facets[0] == v || ec.facets[1] == v) at.push_back(e);
        }
        std::vector<std::array<std::size_t, 2>> joins;
        for (std::size_t f : c.facets) {
          std::vector<std::size_t> here;
          for (std::size_t e : X.cells[f].facets) {
            const auto& ec = X.cells[e];
            if (ec.facets[0] == v || ec.facets[1] == v) here.push_back(e);
          }
          if (here.empty()) continue;
          if (here.size() != 2) return detail::refuse(s, me + "'s boundary is not a surface at vertex #" + std::to_string(X.cells[v].id));
          joins.push_back({here[0], here[1]});
        }
        if (!detail::connected(at, joins)) return detail::refuse(s, me + "'s boundary is pinched at vertex #" + std::to_string(X.cells[v].id) + ": its link is not one circle");
      }
      std::vector<std::array<std::size_t, 2>> adjacent;
      for (const auto& [e, list] : faces_on) adjacent.push_back({list[0], list[1]});
      if (!detail::connected(c.facets, adjacent)) return detail::refuse(s, me + "'s boundary is not connected");
      const long long chi = static_cast<long long>(parts[0].size()) - static_cast<long long>(parts[1].size()) + static_cast<long long>(parts[2].size());
      if (chi != 2) return detail::refuse(s, me + "'s boundary is a closed surface of Euler characteristic " + std::to_string(chi) + ", not a sphere");
      continue;
    }
    // k ≥ 4: a pseudomanifold boundary with the homology over ℤ/2 of S^{k−1} (signs are 1 modulo 2).
    std::vector<detail::Diamond> dia;
    std::string why;
    if (!detail::diamonds(X, s, dia, why)) return detail::refuse(s, me + ": " + why);
    std::vector<std::size_t> boundary;
    for (std::size_t q : cl)
      if (q != s) boundary.push_back(q);
    std::map<std::size_t, std::size_t> pos;
    std::vector<std::size_t> n(c.dim, 0);
    for (std::size_t q : boundary) pos[q] = n[X.cells[q].dim]++;
    std::vector<std::vector<std::vector<std::size_t>>> columns(c.dim);  // columns[k][j]: the rows of ∂_k's column j
    for (std::size_t k = 1; k < c.dim; k++) columns[k].resize(n[k]);
    for (std::size_t q : boundary) {
      const auto& qc = X.cells[q];
      for (std::size_t f : qc.facets) columns[qc.dim][pos[q]].push_back(pos[f]);
    }
    std::vector<std::size_t> rank(c.dim + 1, 0);
    for (std::size_t k = 1; k < c.dim; k++) rank[k] = detail::rank_mod2(n[k - 1], columns[k]);
    for (std::size_t k = 0; k < c.dim; k++) {
      const std::size_t b = n[k] - rank[k] - rank[k + 1];
      if (b != (k == 0 || k + 1 == c.dim ? 1u : 0u))
        return detail::refuse(s, me + "'s boundary does not have the homology of a sphere (b" + std::to_string(k) + " = " + std::to_string(b) + " over ℤ/2)");
    }
    ok.exact = false;
  }
  return ok;
}

// Fills every sign that is 0 so that ∂∂ = 0, keeping the signs given; throws std::invalid_argument when
// the complex is not regular (cells of dimension ≥ 4 need a pseudomanifold boundary) or the given signs
// conflict.
inline void orient(RegularComplex& X) {
  std::vector<std::size_t> order(X.size());
  for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return X.cells[a].dim < X.cells[b].dim; });
  for (std::size_t s : order) {
    auto& c = X.cells[s];
    if (c.dim == 0) continue;
    if (c.dim == 1) {
      if (c.facets.size() != 2) throw std::invalid_argument("morphe::topology::orient: " + detail::named(X, s) + " does not have two ends");
      if (c.signs[0] == 0 && c.signs[1] == 0) c.signs = {-1, 1};
      else if (c.signs[0] == 0) c.signs[0] = -c.signs[1];
      else if (c.signs[1] == 0) c.signs[1] = -c.signs[0];
      if (c.signs[0] + c.signs[1] != 0) throw std::invalid_argument("morphe::topology::orient: " + detail::named(X, s) + " has ends of one sign");
      continue;
    }
    std::vector<detail::Diamond> dia;
    std::string why;
    if (!detail::diamonds(X, s, dia, why)) throw std::invalid_argument("morphe::topology::orient: " + detail::named(X, s) + ": " + why);
    const std::size_t m = c.facets.size();
    std::vector<std::vector<std::pair<std::size_t, int>>> next(m);  // position → (position, factor): s_b = factor · s_a
    for (const auto& d : dia) {
      const int ia = X.incidence(c.facets[d.a], d.rho), ib = X.incidence(c.facets[d.b], d.rho);
      if (ia == 0 || ib == 0) throw std::invalid_argument("morphe::topology::orient: a facet of " + detail::named(X, s) + " is not oriented");
      const int factor = -ia * ib;
      next[d.a].push_back({d.b, factor});
      next[d.b].push_back({d.a, factor});
    }
    std::size_t seed = 0;
    for (std::size_t i = 0; i < m; i++)
      if (c.signs[i] != 0) { seed = i; break; }
    std::vector<int> sign(m, 0);
    sign[seed] = c.signs[seed] != 0 ? c.signs[seed] : 1;
    std::vector<std::size_t> queue{seed};
    for (std::size_t q = 0; q < queue.size(); q++) {
      const std::size_t a = queue[q];
      for (const auto& [b, factor] : next[a]) {
        const int want = factor * sign[a];
        if (sign[b] == 0) sign[b] = want, queue.push_back(b);
        else if (sign[b] != want) throw std::invalid_argument("morphe::topology::orient: " + detail::named(X, s) + " has no consistent orientation of its boundary");
      }
    }
    for (std::size_t i = 0; i < m; i++) {
      if (c.signs[i] != 0 && c.signs[i] != sign[i]) throw std::invalid_argument("morphe::topology::orient: the signs given on " + detail::named(X, s) + " conflict with ∂∂ = 0");
      c.signs[i] = sign[i];
    }
  }
}

// The cellular chain complex over ℤ: C_k on the k-cells in their order (local_index()), ∂σ = Σ [σ:τ] τ.
// Every sign must be ±1 (orient() first).
inline algebra::ChainComplex chain_complex(const RegularComplex& X) {
  const auto pos = X.local_index();
  algebra::ChainComplex C;
  C.n = X.counts();
  if (C.n.empty()) C.n = {0};
  C.boundary.resize(C.n.size());
  for (std::size_t k = 1; k < C.n.size(); k++) C.boundary[k] = algebra::SparseMatrix(C.n[k - 1], C.n[k]);
  for (std::size_t i = 0; i < X.size(); i++) {
    const auto& c = X.cells[i];
    if (c.dim == 0) continue;
    algebra::Vector col;
    for (std::size_t j = 0; j < c.facets.size(); j++) {
      if (c.signs[j] == 0) throw std::invalid_argument("morphe::topology::chain_complex: a cell is not oriented (orient() first)");
      col = algebra::axpy(col, c.signs[j], algebra::Vector{{pos[c.facets[j]], 1}});
    }
    C.boundary[c.dim].columns[pos[i]] = std::move(col);
  }
  return C;
}

// The same complex as a CWComplex: vertices and edges in their order, each 2-cell attached along its
// boundary read as a closed path in the direction of its signs, higher cells by their incidences. Requires
// a regular, oriented complex.
inline CWComplex cw(const RegularComplex& X, std::string name = "") {
  const auto pos = X.local_index();
  CWComplex Y;
  Y.name = name.empty() ? X.name : std::move(name);
  const auto n = X.counts();
  Y.vertices = n.empty() ? 0 : n[0];
  for (std::size_t i = 0; i < X.size(); i++) {
    const auto& c = X.cells[i];
    if (c.dim == 1) {
      const auto [t, h] = X.ends(i);
      Y.edges.push_back({pos[t], pos[h]});
      Y.edge_names.push_back(c.type.empty() ? "" : c.type + "#" + std::to_string(c.id));
    }
  }
  for (std::size_t i = 0; i < X.size(); i++) {
    const auto& c = X.cells[i];
    if (c.dim == 2) {
      // Walk the boundary: each edge in the direction of its sign, starting where the last ended.
      const std::size_t m = c.facets.size();
      std::vector<char> used(m, 0);
      const auto start_of = [&](std::size_t j) { const auto e = X.ends(c.facets[j]); return c.signs[j] > 0 ? e[0] : e[1]; };
      const auto end_of = [&](std::size_t j) { const auto e = X.ends(c.facets[j]); return c.signs[j] > 0 ? e[1] : e[0]; };
      Face f;
      f.base = pos[start_of(0)];
      std::size_t at = end_of(0);
      used[0] = 1;
      f.boundary.push_back({pos[c.facets[0]], c.signs[0]});
      for (std::size_t step = 1; step < m; step++) {
        std::size_t j = m;
        for (std::size_t q = 0; q < m; q++)
          if (!used[q] && start_of(q) == at) { j = q; break; }
        if (j == m) throw std::invalid_argument("morphe::topology::cw: the boundary of " + detail::named(X, i) + " is not a directed circle");
        used[j] = 1;
        f.boundary.push_back({pos[c.facets[j]], c.signs[j]});
        at = end_of(j);
      }
      Y.faces.push_back(std::move(f));
    }
  }
  for (std::size_t i = 0; i < X.size(); i++) {
    const auto& c = X.cells[i];
    if (c.dim < 3) continue;
    if (Y.higher.size() < c.dim - 2) Y.higher.resize(c.dim - 2);
    algebra::Vector v;
    for (std::size_t j = 0; j < c.facets.size(); j++) v = algebra::axpy(v, c.signs[j], algebra::Vector{{pos[c.facets[j]], 1}});
    Y.higher[c.dim - 3].push_back(std::move(v));
  }
  return Y;
}

// ---- cubical grids -------------------------------------------------------------------------------------------

// The cubes of [0, n₀] × [0, n₁] × [0, n₂] (or of a rectangle or segment when n₂ or n₁ is 0) as a regular
// complex: a cell for each product of unit intervals and points, typed by the axes it spans ("" for a
// vertex, "x", "y", "z", "xy", "xz", "yz", "xyz"); corner[i] is its lowest corner and span[i] which axes it
// spans. Its incidence numbers are those of the product orientation, (−1)^{(axes before)}·(±1) for the
// upper and lower face across each spanned axis (Kaczynski, Mischaikow and Mrozek 2004, Computational
// Homology, §2.1).
struct CubicalGrid {
  RegularComplex complex;
  std::vector<std::array<int, 3>> corner;
  std::vector<std::array<int, 3>> span;
  std::array<int, 3> size{};
  // The index of the cell with this corner and span, or no_cell.
  std::size_t at(std::array<int, 3> c, std::array<int, 3> s) const {
    const auto it = index.find({c, s});
    return it == index.end() ? no_cell : it->second;
  }
  std::map<std::pair<std::array<int, 3>, std::array<int, 3>>, std::size_t> index;
};

inline CubicalGrid cubical_grid(int n0, int n1 = 0, int n2 = 0, std::string name = "cubical grid") {
  if (n0 < 0 || n1 < 0 || n2 < 0) throw std::invalid_argument("morphe::topology::cubical_grid: sizes ≥ 0");
  CubicalGrid G;
  G.size = {n0, n1, n2};
  G.complex.name = std::move(name);
  const char* axes = "xyz";
  for (int dim = 0; dim <= 3; dim++)
    for (int s = 0; s < 8; s++) {
      const std::array<int, 3> span{s & 1, (s >> 1) & 1, (s >> 2) & 1};
      if (span[0] + span[1] + span[2] != dim) continue;
      if ((span[0] && n0 == 0) || (span[1] && n1 == 0) || (span[2] && n2 == 0)) continue;
      std::string type;
      for (int a = 0; a < 3; a++)
        if (span[a]) type += axes[a];
      for (int z = 0; z <= n2 - span[2]; z++)
        for (int y = 0; y <= n1 - span[1]; y++)
          for (int x = 0; x <= n0 - span[0]; x++) {
            const std::array<int, 3> c{x, y, z};
            std::vector<std::size_t> facets;
            std::vector<int> signs;
            int before = 0;
            for (int a = 0; a < 3; a++) {
              if (!span[a]) continue;
              std::array<int, 3> fs = span;
              fs[a] = 0;
              std::array<int, 3> up = c;
              up[a] += 1;
              const int e = before % 2 ? -1 : 1;
              facets.push_back(G.at(c, fs)), signs.push_back(-e);
              facets.push_back(G.at(up, fs)), signs.push_back(e);
              before++;
            }
            const std::size_t i = G.complex.add(dim, type, facets, signs);
            G.index[{c, span}] = i;
            G.corner.push_back(c);
            G.span.push_back(span);
          }
    }
  return G;
}

}  // namespace morphe::topology
