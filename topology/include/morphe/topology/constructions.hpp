// A catalogue of complexes: polygons with their sides identified, the standard CW
// structures of spheres, surfaces and projective spaces, presentation complexes, products
// and lattices, and the standard small triangulations.
//
// - Polygons glued by words (Massey 1991, ch. 1; Hatcher 2002, §1.2): each polygon is a cyclic
//   word in labels, a lowercase letter for a side read along the polygon's orientation and an
//   uppercase letter for a side read against it ("abAB" is the torus, "aa" the projective
//   plane, "abAb" the Klein bottle, "aaA" the dunce hat). Sides with the same label are one
//   edge; the corners are identified by union-find on the ends of the sides.
// - Products of CW complexes (Hatcher, Theorem A.6; Massey, ch. IX): the cells are products of
//   cells, with ∂(σ × τ) = ∂σ × τ + (-1)^{|σ|} σ × ∂τ; the 2-cells keep attaching words (an
//   edge times an edge is a square), so that a lattice torus carries gauge fields.
// - Presentation complexes (Hatcher, Corollary 1.28): one vertex, a loop per generator and a
//   2-cell per relator, so that π_1 is the presented group.
// - Triangulations: the boundary of the simplex; the 7-vertex torus of Möbius (1886) and
//   Császár (1949); the 6-vertex projective plane (the hemi-icosahedron); the octahedral
//   sphere; suspensions.

#pragma once

#include <cctype>
#include <cstddef>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <morphe/algebra/presentation.hpp>

#include "cw_complex.hpp"
#include "simplicial.hpp"

namespace morphe::topology {

// A side of a polygon: a label and whether it is read along (+1) or against (-1) the polygon.
struct Side {
  std::string label;
  int sign = 1;
};
using Polygon = std::vector<Side>;

// Polygons from text: polygons separated by spaces or commas, each a run of letters, a
// lowercase letter a side along the polygon and its uppercase a side against it. A label
// may also be a letter followed by digits ("a1b1A1B1").
inline std::vector<Polygon> parse_polygons(std::string_view text) {
  std::vector<Polygon> out;
  Polygon cur;
  for (std::size_t i = 0; i < text.size();) {
    const char c = text[i];
    if (c == ' ' || c == ',' || c == ';' || c == '\n' || c == '\t') {
      if (!cur.empty()) out.push_back(std::move(cur)), cur.clear();
      i++;
      continue;
    }
    if (!std::isalpha(static_cast<unsigned char>(c))) throw std::invalid_argument("morphe::topology::parse_polygons: a side must be a letter, possibly with digits");
    std::string label(1, static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    const int sign = std::islower(static_cast<unsigned char>(c)) ? 1 : -1;
    i++;
    while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) label += text[i++];
    cur.push_back({label, sign});
  }
  if (!cur.empty()) out.push_back(std::move(cur));
  return out;
}

// The CW complex of polygons glued along equally labelled sides. Corner i of a polygon is
// where its side i starts; a side +x runs from corner i to corner i+1, a side -x from corner
// i+1 to corner i. A polygon with no sides is a disc glued to a point of its own (a sphere).
inline CWComplex polygon_complex(const std::vector<Polygon>& polygons, std::string name = "") {
  std::map<std::string, std::size_t> edge_of;
  std::vector<std::string> names;
  std::vector<std::size_t> first_corner;
  std::size_t corners = 0;
  for (const auto& P : polygons) {
    first_corner.push_back(corners);
    corners += std::max<std::size_t>(P.size(), 1);
    for (const auto& s : P)
      if (edge_of.emplace(s.label, names.size()).second) names.push_back(s.label);
  }
  std::vector<std::size_t> parent(corners);
  std::iota(parent.begin(), parent.end(), 0);
  std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = root(parent[x]); };
  auto unite = [&](std::size_t a, std::size_t b) { parent[root(a)] = root(b); };
  // The corners at the tail and the head of each occurrence of each edge.
  std::vector<std::size_t> tail(names.size(), SIZE_MAX), head(names.size(), SIZE_MAX);
  for (std::size_t p = 0; p < polygons.size(); p++) {
    const auto& P = polygons[p];
    for (std::size_t i = 0; i < P.size(); i++) {
      const std::size_t a = first_corner[p] + i, b = first_corner[p] + (i + 1) % P.size();
      const std::size_t e = edge_of.at(P[i].label);
      const std::size_t t = P[i].sign > 0 ? a : b, h = P[i].sign > 0 ? b : a;
      if (tail[e] == SIZE_MAX) tail[e] = t, head[e] = h;
      else unite(t, tail[e]), unite(h, head[e]);
    }
  }
  std::map<std::size_t, std::size_t> vertex_of;
  for (std::size_t c = 0; c < corners; c++) vertex_of.emplace(root(c), vertex_of.size());
  CWComplex X;
  X.name = std::move(name);
  X.vertices = vertex_of.size();
  X.edge_names = names;
  for (std::size_t e = 0; e < names.size(); e++) X.edges.push_back({vertex_of.at(root(tail[e])), vertex_of.at(root(head[e]))});
  for (std::size_t p = 0; p < polygons.size(); p++) {
    Face f;
    f.base = vertex_of.at(root(first_corner[p]));
    for (const auto& s : polygons[p]) f.boundary.push_back({edge_of.at(s.label), s.sign});
    X.faces.push_back(f);
  }
  X.validate();
  return X;
}
inline CWComplex polygon_complex(std::string_view words, std::string name = "") {
  return polygon_complex(parse_polygons(words), name.empty() ? std::string(words) : std::move(name));
}

// The closed orientable surface of genus g: a 4g-gon a1 b1 A1 B1 ... (g = 0: the sphere as aA).
inline CWComplex orientable_surface(std::size_t g) {
  if (g == 0) return polygon_complex("aA", "sphere");
  if (g == 1) return polygon_complex("abAB", "torus");
  std::string w;
  for (std::size_t i = 1; i <= g; i++) {
    const std::string k = std::to_string(i);
    w += "a" + k + "b" + k + "A" + k + "B" + k;
  }
  return polygon_complex(w, "surface of genus " + std::to_string(g));
}
// The closed non-orientable surface with k crosscaps: a 2k-gon c1 c1 c2 c2 ...
inline CWComplex nonorientable_surface(std::size_t k) {
  if (k == 0) throw std::invalid_argument("morphe::topology::nonorientable_surface: k >= 1");
  if (k == 1) return polygon_complex("aa", "projective plane");
  std::string w;
  for (std::size_t i = 1; i <= k; i++) w += "c" + std::to_string(i) + "c" + std::to_string(i);
  return polygon_complex(w, k == 2 ? "Klein bottle" : "surface with " + std::to_string(k) + " crosscaps");
}
inline CWComplex klein_bottle() { return polygon_complex("abAb", "Klein bottle"); }
// Zeeman's dunce hat: a triangle with sides a a A; contractible but not collapsible.
inline CWComplex dunce_hat() { return polygon_complex("aaA", "dunce hat"); }

// The n-sphere as a point with one n-cell (n >= 1).
inline CWComplex sphere(std::size_t n) {
  CWComplex X;
  X.name = "S" + std::to_string(n);
  X.vertices = 1;
  if (n == 1) X.edges.push_back({0, 0});
  else if (n == 2) X.faces.push_back({0, {}});
  else if (n >= 3) {
    X.higher.resize(n - 2);
    X.higher[n - 3].push_back({});
  }
  return X;
}

// ℝPⁿ with one cell in each dimension: e_1 a loop, e_2 attached by a a, and
// ∂e_k = (1 + (-1)^k) e_{k-1} for k >= 3 (Hatcher, Example 2.42).
inline CWComplex real_projective_space(std::size_t n) {
  CWComplex X;
  X.name = "RP" + std::to_string(n);
  X.vertices = 1;
  if (n >= 1) X.edges.push_back({0, 0}), X.edge_names = {"a"};
  if (n >= 2) X.faces.push_back({0, {{0, 1}, {0, 1}}});
  for (std::size_t k = 3; k <= n; k++) X.higher.push_back({k % 2 ? algebra::Vector{} : algebra::Vector{{0, 2}}});
  return X;
}

// ℂPⁿ with one cell in each even dimension and zero boundary (Hatcher, Example 2.40).
inline CWComplex complex_projective_space(std::size_t n) {
  CWComplex X;
  X.name = "CP" + std::to_string(n);
  X.vertices = 1;
  if (n >= 1) X.faces.push_back({0, {}});
  for (std::size_t k = 3; k <= 2 * n; k++) X.higher.push_back(k % 2 ? std::vector<algebra::Vector>{} : std::vector<algebra::Vector>{{}});
  return X;
}

// The lens space L(p, q) with one cell in each dimension 0..3: e_2 attached by a^p, ∂e_3 = 0.
// The CW structure records only the degree of the 3-cell's attaching map, so its homology and
// fundamental group (ℤ/p) are those of L(p, q) for every q; q itself is not represented.
inline CWComplex lens_space(std::size_t p) {
  CWComplex X;
  X.name = "L(" + std::to_string(p) + ", q)";
  X.vertices = 1;
  X.edges.push_back({0, 0});
  X.edge_names = {"a"};
  X.faces.push_back({0, Path(p, Step{0, 1})});
  X.higher.push_back({{}});
  return X;
}

// The presentation complex of ⟨x_1..x_n | r_1..r_m⟩.
inline CWComplex presentation_complex(const algebra::Presentation& P) {
  CWComplex X;
  X.name = P.str();
  X.vertices = 1;
  for (std::size_t i = 0; i < P.generators; i++) X.edges.push_back({0, 0}), X.edge_names.push_back(P.name(i));
  for (const auto& r : P.relators) {
    Face f{0, {}};
    for (int x : r) f.boundary.push_back({static_cast<std::size_t>(std::abs(x)) - 1, x > 0 ? 1 : -1});
    X.faces.push_back(f);
  }
  return X;
}

// The cycle graph of length n as a 1-dimensional CW complex (a circle with n vertices).
inline CWComplex circle(std::size_t n) {
  if (n == 0) throw std::invalid_argument("morphe::topology::circle: n >= 1");
  CWComplex X;
  X.name = "circle of " + std::to_string(n);
  X.vertices = n;
  for (std::size_t i = 0; i < n; i++) X.edges.push_back({i, (i + 1) % n});
  return X;
}
// The interval with n edges.
inline CWComplex interval(std::size_t n) {
  CWComplex X;
  X.name = "interval of " + std::to_string(n);
  X.vertices = n + 1;
  for (std::size_t i = 0; i < n; i++) X.edges.push_back({i, i + 1});
  return X;
}

// The product X × Y. Cells e × f are numbered, in each dimension, by the dimension of e,
// then e, then f; vertex (x, y) is x |Y_0| + y.
inline CWComplex product(const CWComplex& X, const CWComplex& Y) {
  const std::size_t dx = X.dimension(), dy = Y.dimension(), d = dx + dy;
  // offset[k][p]: the index of the first cell e_p × f_{k-p} among the k-cells.
  std::vector<std::vector<std::size_t>> offset(d + 1, std::vector<std::size_t>(dx + 2, 0));
  for (std::size_t k = 0; k <= d; k++)
    for (std::size_t p = 0; p <= dx; p++) {
      const std::size_t q = k >= p ? k - p : SIZE_MAX;
      const std::size_t n = q != SIZE_MAX && q <= dy ? X.cells(p) * Y.cells(q) : 0;
      offset[k][p + 1] = offset[k][p] + n;
    }
  auto id = [&](std::size_t p, std::size_t a, std::size_t q, std::size_t b) { return offset[p + q][p] + a * Y.cells(q) + b; };
  const auto CX = X.chain_complex(), CY = Y.chain_complex();
  CWComplex Z;
  Z.name = X.name + " × " + Y.name;
  Z.vertices = X.vertices * Y.vertices;
  Z.edges.resize(offset[1][dx + 1]);
  Z.faces.resize(d >= 2 ? offset[2][dx + 1] : 0);
  for (std::size_t x = 0; x < X.vertices; x++)
    for (std::size_t e = 0; e < Y.edges.size(); e++) Z.edges[id(0, x, 1, e)] = {id(0, x, 0, Y.edges[e][0]), id(0, x, 0, Y.edges[e][1])};
  for (std::size_t e = 0; e < X.edges.size(); e++)
    for (std::size_t y = 0; y < Y.vertices; y++) Z.edges[id(1, e, 0, y)] = {id(0, X.edges[e][0], 0, y), id(0, X.edges[e][1], 0, y)};
  // 2-cells: vertex × face, edge × edge (a square), face × vertex.
  for (std::size_t x = 0; x < X.vertices; x++)
    for (std::size_t f = 0; f < Y.faces.size(); f++) {
      Face g{id(0, x, 0, Y.faces[f].base), {}};
      for (const auto& s : Y.faces[f].boundary) g.boundary.push_back({id(0, x, 1, s.edge), s.sign});
      Z.faces[id(0, x, 2, f)] = g;
    }
  for (std::size_t e = 0; e < X.edges.size(); e++)
    for (std::size_t h = 0; h < Y.edges.size(); h++) {
      const auto [u, v] = X.edges[e];
      const auto [s, t] = Y.edges[h];
      // (e × s) (v × h) (e × t)^{-1} (u × h)^{-1}, from (u, s): ∂ = e×s + v×h - e×t - u×h.
      Z.faces[id(1, e, 1, h)] = {id(0, u, 0, s), {{id(1, e, 0, s), 1}, {id(0, v, 1, h), 1}, {id(1, e, 0, t), -1}, {id(0, u, 1, h), -1}}};
    }
  for (std::size_t f = 0; f < X.faces.size(); f++)
    for (std::size_t y = 0; y < Y.vertices; y++) {
      Face g{id(0, X.faces[f].base, 0, y), {}};
      for (const auto& s : X.faces[f].boundary) g.boundary.push_back({id(1, s.edge, 0, y), s.sign});
      Z.faces[id(2, f, 0, y)] = g;
    }
  // Cells of dimension >= 3 by the Leibniz rule on the cellular boundaries.
  for (std::size_t k = 3; k <= d; k++) {
    std::vector<algebra::Vector> cells(offset[k][dx + 1]);
    for (std::size_t p = 0; p <= std::min(k, dx); p++) {
      const std::size_t q = k - p;
      if (q > dy) continue;
      for (std::size_t a = 0; a < X.cells(p); a++)
        for (std::size_t b = 0; b < Y.cells(q); b++) {
          algebra::Vector v;
          if (p >= 1)
            for (const auto& en : CX.boundary[p].columns[a]) v = algebra::axpy(v, en.value, algebra::Vector{{id(p - 1, en.index, q, b), 1}});
          if (q >= 1)
            for (const auto& en : CY.boundary[q].columns[b]) v = algebra::axpy(v, p % 2 ? -en.value : en.value, algebra::Vector{{id(p, a, q - 1, en.index), 1}});
          cells[id(p, a, q, b)] = v;
        }
    }
    Z.higher.push_back(std::move(cells));
  }
  return Z;
}

// The periodic cubical lattice of the given side lengths: the torus T^d as a product of
// circles, with every cell (vertices, links, plaquettes, cubes ...).
inline CWComplex torus_lattice(const std::vector<std::size_t>& sides) {
  if (sides.empty()) throw std::invalid_argument("morphe::topology::torus_lattice: at least one side");
  CWComplex Z = circle(sides[0]);
  for (std::size_t i = 1; i < sides.size(); i++) Z = product(Z, circle(sides[i]));
  std::string name = "torus lattice ";
  for (std::size_t i = 0; i < sides.size(); i++) name += (i ? "×" : "") + std::to_string(sides[i]);
  Z.name = name;
  return Z;
}

// ---- triangulations ---------------------------------------------------------------------------------------------------

// The boundary of the n-simplex: a triangulated (n-1)-sphere on n + 1 vertices.
inline SimplicialComplex simplex_boundary(std::size_t n) {
  std::vector<Simplex> facets;
  for (std::size_t omit = 0; omit <= n; omit++) {
    Simplex s;
    for (std::size_t v = 0; v <= n; v++)
      if (v != omit) s.push_back(v);
    facets.push_back(s);
  }
  return SimplicialComplex::from_facets(facets, n + 1, "boundary of the " + std::to_string(n) + "-simplex");
}
inline SimplicialComplex simplex(std::size_t n) {
  Simplex s(n + 1);
  std::iota(s.begin(), s.end(), 0);
  return SimplicialComplex::from_facets({s}, n + 1, std::to_string(n) + "-simplex");
}
// The 7-vertex torus (Möbius 1886; Császár 1949): triangles {i, i+1, i+3}, {i, i+2, i+3} mod 7.
inline SimplicialComplex torus_7() {
  std::vector<Simplex> f;
  for (std::size_t i = 0; i < 7; i++) {
    f.push_back({i, (i + 1) % 7, (i + 3) % 7});
    f.push_back({i, (i + 2) % 7, (i + 3) % 7});
  }
  return SimplicialComplex::from_facets(f, 7, "7-vertex torus");
}
// The 6-vertex projective plane, the quotient of the icosahedron by the antipodal map.
inline SimplicialComplex projective_plane_6() {
  return SimplicialComplex::from_facets({{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 1},
                                         {1, 2, 4}, {2, 3, 5}, {3, 4, 1}, {4, 5, 2}, {5, 1, 3}},
                                        6, "6-vertex projective plane");
}
// The octahedron: the boundary of the cross-polytope, a 2-sphere on 6 vertices.
inline SimplicialComplex octahedron() {
  std::vector<Simplex> f;
  for (std::size_t a : {0, 1})
    for (std::size_t b : {2, 3})
      for (std::size_t c : {4, 5}) f.push_back({a, b, c});
  return SimplicialComplex::from_facets(f, 6, "octahedron");
}

}  // namespace morphe::topology
