// Surfaces: recognition and classification of 2-dimensional CW complexes.
//
// A connected 2-dimensional CW complex whose 2-cells are polygons glued along their sides is
// a compact surface exactly when (Massey 1991, ch. 1, §§5–8; Gallier and Xu 2013, A Guide to
// the Classification Theorem for Compact Surfaces, ch. 4):
// 1. every edge is a side of two polygons, or of one (a boundary edge), counted with
//    multiplicity, and every edge is the side of some polygon;
// 2. at every vertex the corners of the polygons, joined along the sides they share, form a
//    single cycle (an interior vertex) or a single path (a vertex on the boundary): the link
//    of the vertex is a circle or an arc.
// It is orientable exactly when its polygons can be oriented so that every interior edge is
// traversed once in each direction. Then the classification theorem (Möbius, Jordan, Dehn and
// Heegaard; Brahana 1921) says the surface is determined by its orientability, its Euler
// characteristic χ and its number b of boundary circles: χ = 2 - 2g - b (orientable, genus
// g) or χ = 2 - k - b (k crosscaps). The normal form of a closed surface is the word
// a1 b1 a1⁻¹ b1⁻¹ ... ag bg ag⁻¹ bg⁻¹, or c1 c1 ... ck ck.
//
// The one exception to polygons is the sphere as a point with a disc whose word is empty,
// recognised as such. All tests are exact and linear in the size of the complex.

#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <morphe/algebra/presentation.hpp>

#include "cw_complex.hpp"

namespace morphe::topology {

struct SurfaceType {
  bool surface = false;
  bool closed = false;
  bool orientable = false;
  long long euler = 0;
  std::size_t genus = 0;      // orientable genus
  std::size_t crosscaps = 0;  // non-orientable genus
  std::size_t boundary_components = 0;
  std::string name;
  algebra::Word normal_form;     // closed surfaces: on generators a1 b1 ... or c1 ...
  std::string reason;            // why it is not a surface
  std::vector<int> orientation;  // a consistent orientation of the 2-cells (+1/-1), when orientable
};

inline SurfaceType classify_surface(const CWComplex& X) {
  SurfaceType S;
  S.euler = X.euler_characteristic();
  auto fail = [&](std::string why) {
    S.surface = false;
    S.reason = std::move(why);
    return S;
  };
  if (X.dimension() != 2) return fail("its dimension is " + std::to_string(X.dimension()) + ", not 2");
  const auto F = combinatorics::spanning_forest(X.one_skeleton());
  if (F.components() != 1) return fail("it is not connected");
  // The sphere as a point and a disc.
  bool empty_word = false;
  for (const auto& f : X.faces) empty_word = empty_word || f.boundary.empty();
  if (empty_word) {
    if (X.faces.size() == 1 && X.edges.empty() && X.vertices == 1) {
      S.surface = S.closed = S.orientable = true;
      S.name = "sphere";
      S.orientation = {1};
      return S;
    }
    return fail("a 2-cell is attached to a point beside other cells");
  }
  // 1. Occurrences of each edge among the sides: (face, sign).
  std::vector<std::vector<std::pair<std::size_t, int>>> occurrences(X.edges.size());
  for (std::size_t f = 0; f < X.faces.size(); f++)
    for (const auto& s : X.faces[f].boundary) occurrences[s.edge].push_back({f, s.sign});
  std::size_t boundary_edges = 0;
  for (std::size_t e = 0; e < X.edges.size(); e++) {
    if (occurrences[e].empty()) return fail("the edge " + X.edge_name(e) + " is the side of no 2-cell");
    if (occurrences[e].size() > 2) return fail("the edge " + X.edge_name(e) + " is a side " + std::to_string(occurrences[e].size()) + " times");
    if (occurrences[e].size() == 1) boundary_edges++;
  }
  // 2. Links of vertices. Nodes: the ends of edges (2e: tail end, 2e+1: head end); each
  // corner of a polygon joins the end where one side arrives to the end where the next leaves.
  std::vector<std::vector<std::size_t>> adj(2 * X.edges.size());
  for (const auto& f : X.faces) {
    const std::size_t n = f.boundary.size();
    for (std::size_t i = 0; i < n; i++) {
      const Step& in = f.boundary[i];
      const Step& out = f.boundary[(i + 1) % n];
      const std::size_t arrive = 2 * in.edge + (in.sign > 0 ? 1 : 0);
      const std::size_t leave = 2 * out.edge + (out.sign > 0 ? 0 : 1);
      adj[arrive].push_back(leave);
      adj[leave].push_back(arrive);
    }
  }
  std::vector<std::size_t> ends_at(X.vertices, 0), link_parts(X.vertices, 0);
  std::vector<char> seen(adj.size(), 0);
  auto vertex_of_end = [&](std::size_t end) { return X.edges[end / 2][end % 2]; };
  for (std::size_t end = 0; end < adj.size(); end++) {
    ends_at[vertex_of_end(end)]++;
    if (seen[end]) continue;
    link_parts[vertex_of_end(end)]++;
    std::vector<std::size_t> stack{end};
    seen[end] = 1;
    while (!stack.empty()) {
      const std::size_t x = stack.back();
      stack.pop_back();
      for (std::size_t y : adj[x])
        if (!seen[y]) seen[y] = 1, stack.push_back(y);
    }
  }
  for (std::size_t v = 0; v < X.vertices; v++) {
    if (ends_at[v] == 0) return fail("the vertex " + std::to_string(v) + " is on no edge");
    if (link_parts[v] != 1) return fail("the link of vertex " + std::to_string(v) + " has " + std::to_string(link_parts[v]) + " components (a pinch point)");
  }
  // Every node has degree 2 at an interior edge and 1 at a boundary edge, so a connected link
  // is a cycle or a path. Orientability: propagate orientations of the 2-cells across interior
  // edges, which must be traversed oppositely: o[f1] s1 = -o[f2] s2.
  S.orientable = true;
  std::vector<int> o(X.faces.size(), 0);
  for (std::size_t start = 0; start < X.faces.size(); start++) {
    if (o[start]) continue;
    o[start] = 1;
    std::vector<std::size_t> stack{start};
    while (!stack.empty()) {
      const std::size_t f = stack.back();
      stack.pop_back();
      for (const auto& s : X.faces[f].boundary) {
        const auto& occ = occurrences[s.edge];
        if (occ.size() != 2) continue;
        const auto [f1, s1] = occ[0];
        const auto [f2, s2] = occ[1];
        if (f1 == f2) {
          if (s1 == s2) S.orientable = false;
          continue;
        }
        const std::size_t other = f1 == f ? f2 : f1;
        const int mine = f1 == f ? s1 : s2, theirs = f1 == f ? s2 : s1;
        const int want = -o[f] * mine * theirs;
        if (!o[other]) o[other] = want, stack.push_back(other);
        else if (o[other] != want) S.orientable = false;
      }
    }
  }
  if (S.orientable) S.orientation = o;
  // Boundary circles: the boundary edges form disjoint cycles; count their components.
  {
    std::vector<std::size_t> parent(X.vertices);
    for (std::size_t v = 0; v < X.vertices; v++) parent[v] = v;
    std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = root(parent[x]); };
    std::vector<char> on(X.vertices, 0);
    for (std::size_t e = 0; e < X.edges.size(); e++)
      if (occurrences[e].size() == 1) {
        on[X.edges[e][0]] = on[X.edges[e][1]] = 1;
        parent[root(X.edges[e][0])] = root(X.edges[e][1]);
      }
    for (std::size_t v = 0; v < X.vertices; v++)
      if (on[v] && root(v) == v) S.boundary_components++;
  }
  S.surface = true;
  S.closed = boundary_edges == 0;
  const long long b = static_cast<long long>(S.boundary_components);
  if (S.orientable) S.genus = static_cast<std::size_t>((2 - S.euler - b) / 2);
  else S.crosscaps = static_cast<std::size_t>(2 - S.euler - b);
  if (S.orientable) {
    if (S.closed) S.name = S.genus == 0 ? "sphere" : S.genus == 1 ? "torus" : "orientable surface of genus " + std::to_string(S.genus);
    else if (S.genus == 0 && b == 1) S.name = "disc";
    else if (S.genus == 0 && b == 2) S.name = "annulus";
    else if (S.genus == 0 && b == 3) S.name = "pair of pants";
    else S.name = "orientable surface of genus " + std::to_string(S.genus) + " with " + std::to_string(b) + " boundary circles";
    for (std::size_t i = 0; i < S.genus && S.closed; i++) {
      const int a = static_cast<int>(2 * i + 1), c = a + 1;
      S.normal_form.insert(S.normal_form.end(), {a, c, -a, -c});
    }
  } else {
    if (S.closed) S.name = S.crosscaps == 1 ? "projective plane" : S.crosscaps == 2 ? "Klein bottle" : "non-orientable surface with " + std::to_string(S.crosscaps) + " crosscaps";
    else if (S.crosscaps == 1 && b == 1) S.name = "Möbius band";
    else S.name = "non-orientable surface with " + std::to_string(S.crosscaps) + " crosscaps and " + std::to_string(b) + " boundary circles";
    for (std::size_t i = 0; i < S.crosscaps && S.closed; i++) {
      const int c = static_cast<int>(i + 1);
      S.normal_form.insert(S.normal_form.end(), {c, c});
    }
  }
  return S;
}

}  // namespace morphe::topology
