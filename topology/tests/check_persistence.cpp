// Checks of persistent homology over ℤ/2 (persistence.hpp): a triangle built edge by edge and
// then filled, whose components merge and whose loop is born and killed at the stages the
// construction implies; the essential intervals of surfaces and projective planes, whatever the
// order of their cells, against their Betti numbers over ℤ/2; the Betti numbers of every stage
// against those of the sub-complex computed afresh; and a filtration refused when a cell enters
// before one of its faces.

#include <morphe/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/persistence.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace mt = morphe::topology;

namespace {

// A filtration giving every cell of one dimension its own stage, in order of index, the
// dimensions one after another: vertices first, then edges, then faces.
mt::Filtration one_by_one(const mt::CWComplex& X) {
  mt::Filtration F;
  std::size_t s = 0;
  for (std::size_t i = 0; i < X.vertices; i++) F.vertex.push_back(s++);
  for (std::size_t i = 0; i < X.edges.size(); i++) F.edge.push_back(s++);
  for (std::size_t i = 0; i < X.faces.size(); i++) F.face.push_back(s++);
  return F;
}

std::vector<std::size_t> essential(const mt::PersistenceDiagram& D, std::size_t top) {
  std::vector<std::size_t> n(top + 1, 0);
  for (const auto& I : D.intervals)
    if (I.essential()) n[I.dimension]++;
  return n;
}

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

}  // namespace

int main() {
  morphe::Suite s;

  // A triangle: three vertices at stage 0, its edges at stages 1, 1, 2, its face at stage 3.
  {
    mt::CWComplex X;
    X.vertices = 3;
    X.edges = {{0, 1}, {1, 2}, {0, 2}};
    X.faces.push_back({0, {{0, 1}, {1, 1}, {2, -1}}});
    mt::Filtration F{{0, 0, 0}, {1, 1, 2}, {3}, {}};
    const auto D = mt::persistent_homology(X, F);
    std::size_t h0_finite = 0, h0_essential = 0, h1 = 0;
    bool loop = false;
    for (const auto& I : D.intervals) {
      if (I.dimension == 0) (I.essential() ? h0_essential : h0_finite)++;
      if (I.dimension == 1) {
        h1++;
        loop = I.birth == 2 && I.death && *I.death == 3 && I.creator == 2 && I.killer && *I.killer == 0;
      }
    }
    s.require("a triangle: two components die when the first edges join them, at stage 1", h0_finite == 2 && h0_essential == 1);
    s.require("a triangle: its loop is born with the third edge (stage 2) and killed by the face (stage 3)", h1 == 1 && loop);
    s.equal("a triangle: Betti numbers at stage 2", text({D.betti(0, 2), D.betti(1, 2)}), std::string("1 1"));
    s.equal("a triangle: Betti numbers at stage 3", text({D.betti(0, 3), D.betti(1, 3)}), std::string("1 0"));
  }

  // The essential intervals are the Betti numbers over ℤ/2, whatever the order.
  for (const auto& [name, X] : std::vector<std::pair<std::string, mt::CWComplex>>{
           {"the torus", mt::orientable_surface(1)}, {"the surface of genus 2", mt::orientable_surface(2)},
           {"the projective plane", mt::nonorientable_surface(1)}, {"the Klein bottle", mt::klein_bottle()},
           {"the 2-sphere", mt::sphere(2)}, {"two triangles of the torus net", mt::polygon_complex("abc ABC")}}) {
    const auto want = mt::betti(X, 2);
    const auto one = essential(mt::persistent_homology(X, one_by_one(X)), X.dimension());
    mt::Filtration all{std::vector<std::size_t>(X.vertices, 0), std::vector<std::size_t>(X.edges.size(), 0), std::vector<std::size_t>(X.faces.size(), 0), {}};
    const auto flat = essential(mt::persistent_homology(X, all), X.dimension());
    s.equal(name + ": essential intervals, cells one by one, are the Betti numbers over ℤ/2", text(one), text(want));
    s.equal(name + ": and so with every cell at one stage", text(flat), text(want));
  }

  // At every stage, the Betti numbers of the sub-complex computed afresh.
  {
    const mt::CWComplex X = mt::orientable_surface(2);
    mt::Filtration F;
    F.vertex.assign(X.vertices, 0);
    for (std::size_t i = 0; i < X.edges.size(); i++) F.edge.push_back(1 + i / 2);
    F.face.assign(X.faces.size(), 1 + X.edges.size() / 2);
    const auto D = mt::persistent_homology(X, F);
    bool ok = true;
    std::string where;
    for (std::size_t stage = 0; stage <= 1 + X.edges.size() / 2; stage++) {
      mt::CWComplex Y;
      Y.vertices = X.vertices;
      for (std::size_t i = 0; i < X.edges.size(); i++)
        if (F.edge[i] <= stage) Y.edges.push_back(X.edges[i]);
      if (stage == 1 + X.edges.size() / 2) Y.faces = X.faces;
      const auto b = mt::betti(Y, 2);
      for (std::size_t d = 0; d < b.size(); d++)
        if (D.betti(d, stage) != b[d]) { ok = false; where = "stage " + std::to_string(stage) + ", H_" + std::to_string(d); }
    }
    s.require("the surface of genus 2 built in stages: the Betti numbers of every stage", ok, where);
  }

  // A cell before its face.
  {
    mt::CWComplex X;
    X.vertices = 2;
    X.edges = {{0, 1}};
    bool refused = false;
    try {
      mt::persistent_homology(X, mt::Filtration{{0, 2}, {1}, {}, {}});
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("an edge entering before its vertex is refused", refused);
  }
  return s.done();
}
