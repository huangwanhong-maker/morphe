// Surfaces given as lists of triangles: their 2-complexes read directly, their components, and whether the
// orientation the triangles carry is a consistent one (decision 0010).
//
// A list of triangles on vertices 0..n−1, each an ordered triple (a, b, c), is a 2-dimensional CW complex
// whose edges are the pairs of corners that occur, each once, from the lesser vertex to the greater, and
// whose 2-cells are the triangles, attached along the word a → b → c → a. It is the complex that
// SimplicialComplex::from_facets(...).cw() builds, without the closure under subsets: the edges are found by
// sorting the 3F corner pairs, in O(F log F), so that surfaces of millions of triangles are read in seconds.
// Each 2-cell is attached in the order of its corners, not of its sorted vertices, so the orientation the
// triangles carry is kept: when classify_surface finds the complex orientable, the orientation it returns
// gives every cell the same sign exactly when the triangles' own orientation is consistent (every interior
// edge traversed once each way), as it must be for a surface to bound a solid by the divergence theorem.
//
// classify_surface needs a connected complex; triangle_components splits the triangles by the components of
// the graph of their edges (combinatorics::spanning_forest), and classify_triangles classifies each part.
// Checks: tests/check_triangulated.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/combinatorics/graph.hpp>

#include "cw_complex.hpp"
#include "surfaces.hpp"

namespace morphe::topology {

using TriangleList = std::vector<std::array<std::size_t, 3>>;

namespace detail {
inline std::pair<std::size_t, std::size_t> unordered_edge(std::size_t a, std::size_t b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); }
}  // namespace detail

// The 2-complex of the triangles on vertices 0..vertices−1, each attached along its corners in order.
inline CWComplex triangle_complex(std::size_t vertices, const TriangleList& triangles) {
  CWComplex X;
  X.vertices = vertices;
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  pairs.reserve(3 * triangles.size());
  for (const auto& t : triangles) {
    for (int k = 0; k < 3; k++) {
      const std::size_t a = t[k], b = t[(k + 1) % 3];
      if (a >= vertices || b >= vertices) throw std::out_of_range("morphe::topology::triangle_complex: a corner beyond the vertices");
      if (a == b) throw std::invalid_argument("morphe::topology::triangle_complex: a triangle repeats a corner");
      pairs.push_back(detail::unordered_edge(a, b));
    }
  }
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  X.edges.reserve(pairs.size());
  for (const auto& [a, b] : pairs) X.edges.push_back({a, b});
  auto step = [&](std::size_t a, std::size_t b) {
    const auto it = std::lower_bound(pairs.begin(), pairs.end(), detail::unordered_edge(a, b));
    return Step{static_cast<std::size_t>(it - pairs.begin()), a < b ? +1 : -1};
  };
  X.faces.reserve(triangles.size());
  for (const auto& t : triangles) X.faces.push_back({t[0], {step(t[0], t[1]), step(t[1], t[2]), step(t[2], t[0])}});
  return X;
}

// The triangles grouped by the components of the graph of their edges; vertices on no triangle are left out.
inline std::vector<TriangleList> triangle_components(std::size_t vertices, const TriangleList& triangles) {
  combinatorics::Multigraph G{vertices, {}};
  G.edges.reserve(2 * triangles.size());
  for (const auto& t : triangles) G.edges.push_back({t[0], t[1]}), G.edges.push_back({t[1], t[2]});
  const auto F = combinatorics::spanning_forest(G);
  std::vector<std::size_t> index(F.components(), SIZE_MAX);
  std::vector<TriangleList> parts;
  for (const auto& t : triangles) {
    std::size_t& k = index[F.component[t[0]]];
    if (k == SIZE_MAX) k = parts.size(), parts.emplace_back();
    parts[k].push_back(t);
  }
  return parts;
}

struct TriangulatedPart {
  std::size_t vertices = 0, edges = 0, triangles = 0;
  SurfaceType type;         // morphe's classification; type.reason says why it is not a surface
  bool consistent = false;  // a surface whose triangles' own orientation is consistent
};

// Each component of the triangles, renumbered on its own vertices, classified.
inline std::vector<TriangulatedPart> classify_triangles(std::size_t vertices, const TriangleList& triangles) {
  std::vector<TriangulatedPart> out;
  std::vector<std::size_t> local(vertices, SIZE_MAX);
  for (auto part : triangle_components(vertices, triangles)) {
    std::size_t n = 0;
    std::vector<std::size_t> used;
    for (auto& t : part)
      for (auto& v : t) {
        if (local[v] == SIZE_MAX) local[v] = n++, used.push_back(v);
        v = local[v];
      }
    for (std::size_t v : used) local[v] = SIZE_MAX;
    const CWComplex X = triangle_complex(n, part);
    TriangulatedPart P;
    P.vertices = n;
    P.edges = X.edges.size();
    P.triangles = part.size();
    P.type = classify_surface(X);
    if (P.type.surface && P.type.orientable) {
      const auto& o = P.type.orientation;
      P.consistent = std::all_of(o.begin(), o.end(), [&](int s) { return s == o.front(); });
    }
    out.push_back(std::move(P));
  }
  return out;
}

}  // namespace morphe::topology
