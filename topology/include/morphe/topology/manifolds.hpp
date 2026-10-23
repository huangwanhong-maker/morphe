// Combinatorial manifolds and pseudomanifolds among simplicial complexes.
//
// Definitions (Rourke and Sanderson 1972, ch. 2; Seifert and Threlfall 1934, §24).
// - A pure d-dimensional complex is a pseudomanifold (with boundary) when every
//   (d-1)-simplex lies in at most two d-simplices; its boundary is formed by those in one.
// - It is a combinatorial d-manifold when the link of every vertex is a combinatorial
//   (d-1)-sphere or (d-1)-ball. In dimension d <= 3 this is decided exactly here: spheres
//   and balls of dimension 0, 1 and 2 are recognised by their structure and, in dimension 2,
//   by the classification of surfaces (a closed surface with χ = 2, or a surface with one
//   boundary circle and χ = 1). In dimension d >= 4 vertex links are (d-1)-manifolds whose
//   sphericity is not decided here (in dimension 4 it is the recognition of the 3-sphere,
//   decidable but not written; in dimension >= 6 undecidable, Novikov): the result is then
//   that of a homology manifold, the links having the homology of spheres or balls.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "invariants.hpp"
#include "simplicial.hpp"
#include "surfaces.hpp"

namespace morphe::topology {

struct ManifoldType {
  std::size_t dimension = 0;
  bool pure = false;
  bool pseudomanifold = false;
  bool manifold = false;   // a combinatorial manifold (dimension <= 3), or a homology manifold above
  bool exact = true;       // false when `manifold` is only the homology condition (dimension >= 4)
  bool closed = false;     // no boundary
  std::string reason;      // why it is not a manifold
};

namespace detail {

// Whether K (a complex on its own labels) is a combinatorial sphere (want_ball false) or ball
// of dimension d, d <= 2 exactly, d >= 3 by homology.
inline bool sphere_or_ball(const SimplicialComplex& K, std::size_t d, bool ball) {
  if (K.simplices.empty() || K.dimension() != d || !K.is_pure()) return false;
  if (d == 0) return K.count(0) == (ball ? 1u : 2u);
  if (d == 1) {
    // A cycle (sphere) or a path (ball): connected, every vertex of degree 2 (or two of degree 1).
    std::vector<std::size_t> deg(K.vertices, 0);
    for (const auto& e : K.simplices[1]) deg[e[0]]++, deg[e[1]]++;
    std::size_t ones = 0;
    for (const auto& v : K.simplices[0]) {
      if (deg[v[0]] == 1) ones++;
      else if (deg[v[0]] != 2) return false;
    }
    const auto H = algebra::betti(K.chain_complex());
    return H[0] == 1 && (ball ? ones == 2 : ones == 0);
  }
  if (d == 2) {
    const SurfaceType S = classify_surface(renumbered(K).cw());
    return S.surface && S.orientable && S.genus == 0 && (ball ? (!S.closed && S.boundary_components == 1) : S.closed);
  }
  // d >= 3: the homology of a sphere (ℤ, 0, ..., 0, ℤ) or of a ball (ℤ, 0, ..., 0).
  const auto H = algebra::homology(K.chain_complex(), false);
  for (std::size_t k = 0; k <= d; k++) {
    const std::size_t want = k == 0 ? 1 : (k == d && !ball ? 1 : 0);
    if (H[k].rank != want || !H[k].torsion.empty()) return false;
  }
  return true;
}

}  // namespace detail

inline ManifoldType manifold_type(const SimplicialComplex& K) {
  ManifoldType M;
  M.dimension = K.dimension();
  const std::size_t d = M.dimension;
  M.pure = K.is_pure();
  if (!M.pure) { M.reason = "not pure"; return M; }
  if (d == 0) { M.pseudomanifold = M.manifold = M.closed = true; return M; }
  // Pseudomanifold: every (d-1)-simplex in one or two d-simplices.
  std::vector<std::size_t> cofaces(K.count(d - 1), 0);
  for (const auto& t : K.simplices[d])
    for (std::size_t i = 0; i < t.size(); i++) {
      Simplex f = t;
      f.erase(f.begin() + static_cast<long>(i));
      cofaces[K.find(f)]++;
    }
  M.pseudomanifold = true;
  M.closed = true;
  for (std::size_t i = 0; i < cofaces.size(); i++) {
    if (cofaces[i] > 2) { M.pseudomanifold = false; M.reason = "a (d-1)-simplex lies in " + std::to_string(cofaces[i]) + " d-simplices"; }
    if (cofaces[i] == 1) M.closed = false;
  }
  if (!M.pseudomanifold) return M;
  M.exact = d <= 3;
  for (const auto& v : K.simplices[0]) {
    const SimplicialComplex L = renumbered(K.link(v));
    const bool interior = detail::sphere_or_ball(L, d - 1, false);
    if (!interior && !detail::sphere_or_ball(L, d - 1, true)) {
      M.reason = "the link of vertex " + std::to_string(v[0]) + " is neither a sphere nor a ball";
      return M;
    }
  }
  M.manifold = true;
  return M;
}

}  // namespace morphe::topology
