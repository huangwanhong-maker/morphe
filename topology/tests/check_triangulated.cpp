// Checks of surfaces given as lists of triangles (triangulated.hpp; decision 0010).
//
// - The boundary of a tetrahedron, oriented outward: one part, a sphere, χ = 2, its orientation consistent;
//   with one triangle turned over, still an orientable sphere, but its own orientation no longer consistent.
// - The 3 × 3 grid on the torus (9 vertices, 27 edges, 18 triangles): a torus of genus 1, consistent.
// - A Möbius strip of six triangles: not orientable, one boundary circle, one crosscap.
// - A fan of six triangles about a centre: a disc.
// - Two tetrahedra apart: two parts, two spheres; two tetrahedra sharing a vertex: one part, not a surface
//   (the shared vertex is a pinch point).
// - On every example, the complex read directly has the cells, Euler characteristic and classification of
//   SimplicialComplex::from_facets(...).cw().

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <morphe/topology/simplicial.hpp>
#include <morphe/topology/triangulated.hpp>
#include <morphe/check.hpp>

using morphe::Suite;
namespace mt = morphe::topology;

namespace {

mt::TriangleList tetrahedron(std::size_t o = 0) {
  return {{o + 0, o + 2, o + 1}, {o + 0, o + 1, o + 3}, {o + 0, o + 3, o + 2}, {o + 1, o + 2, o + 3}};
}

mt::TriangleList torus3() {
  mt::TriangleList T;
  auto v = [](int i, int j) { return std::size_t(((i + 3) % 3) * 3 + (j + 3) % 3); };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      T.push_back({v(i, j), v(i + 1, j), v(i + 1, j + 1)});
      T.push_back({v(i, j), v(i + 1, j + 1), v(i, j + 1)});
    }
  return T;
}

// The same complex through the simplicial route: its cells and classification.
mt::SurfaceType simplicial(std::size_t n, const mt::TriangleList& T, std::vector<std::size_t>& cells) {
  std::vector<mt::Simplex> facets;
  for (const auto& t : T) facets.push_back({t[0], t[1], t[2]});
  const auto X = mt::SimplicialComplex::from_facets(facets, n).cw();
  cells = X.counts();
  return mt::classify_surface(X);
}

}  // namespace

int main() {
  Suite s;
  {
    const auto P = mt::classify_triangles(4, tetrahedron());
    s.require("a tetrahedron's boundary: one part, a sphere, closed, χ = 2", P.size() == 1 && P[0].type.surface && P[0].type.closed && P[0].type.name == "sphere" && P[0].type.euler == 2);
    s.require("its triangles' orientation is consistent", P[0].consistent);
    auto T = tetrahedron();
    std::swap(T[1][1], T[1][2]);
    const auto Q = mt::classify_triangles(4, T);
    s.require("with one triangle turned over: still an orientable sphere, its own orientation inconsistent", Q.size() == 1 && Q[0].type.orientable && !Q[0].consistent);
  }
  {
    const auto P = mt::classify_triangles(9, torus3());
    s.require("the 3 × 3 torus: 9 vertices, 27 edges, 18 triangles", P.size() == 1 && P[0].vertices == 9 && P[0].edges == 27 && P[0].triangles == 18);
    s.require("a torus of genus 1, consistent", P[0].type.surface && P[0].type.closed && P[0].type.genus == 1 && P[0].type.name == "torus" && P[0].consistent);
  }
  {
    const mt::TriangleList M{{0, 1, 4}, {0, 4, 3}, {1, 2, 5}, {1, 5, 4}, {2, 3, 0}, {2, 0, 5}};
    const auto P = mt::classify_triangles(6, M);
    s.require("a Möbius strip: a surface, not orientable, one boundary circle, one crosscap",
              P.size() == 1 && P[0].type.surface && !P[0].type.orientable && P[0].type.boundary_components == 1 && P[0].type.crosscaps == 1 && !P[0].consistent);
  }
  {
    mt::TriangleList F;
    for (std::size_t k = 0; k < 6; k++) F.push_back({0, 1 + k, 1 + (k + 1) % 6});
    const auto P = mt::classify_triangles(7, F);
    s.require("a fan of six triangles: a disc", P.size() == 1 && P[0].type.name == "disc" && P[0].consistent);
  }
  {
    auto T = tetrahedron();
    for (const auto& t : tetrahedron(4)) T.push_back(t);
    const auto P = mt::classify_triangles(8, T);
    s.require("two tetrahedra apart: two parts, two spheres", P.size() == 2 && P[0].type.name == "sphere" && P[1].type.name == "sphere");
    auto U = tetrahedron();
    for (const auto& t : tetrahedron(3)) U.push_back(t);
    const auto Q = mt::classify_triangles(7, U);
    s.require("two tetrahedra sharing a vertex: one part, not a surface (a pinch point)", Q.size() == 1 && !Q[0].type.surface);
  }
  {
    struct Example {
      std::string name;
      std::size_t n;
      mt::TriangleList T;
    };
    mt::TriangleList fan;
    for (std::size_t k = 0; k < 6; k++) fan.push_back({0, 1 + k, 1 + (k + 1) % 6});
    const std::vector<Example> examples{{"tetrahedron", 4, tetrahedron()}, {"torus", 9, torus3()},
                                        {"Möbius strip", 6, {{0, 1, 4}, {0, 4, 3}, {1, 2, 5}, {1, 5, 4}, {2, 3, 0}, {2, 0, 5}}}, {"fan", 7, fan}};
    bool same = true;
    for (const auto& e : examples) {
      std::vector<std::size_t> cells;
      const auto A = simplicial(e.n, e.T, cells);
      const auto X = mt::triangle_complex(e.n, e.T);
      const auto B = mt::classify_surface(X);
      same = same && X.counts() == cells && A.euler == B.euler && A.name == B.name && A.orientable == B.orientable && A.boundary_components == B.boundary_components;
    }
    s.require("on every example the complex read directly agrees with from_facets(...).cw(): cells, χ, classification", same);
  }
  return s.done();
}
