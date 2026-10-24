// Checks of the cycles over ℤ/2 (cycles.hpp) and of Atkin's Q-analysis (q_analysis.hpp).
//
// Cycles: on closed surfaces the 2-cycles over ℤ/2 are one, the whole surface, orientable or not
// (the projective plane and the Klein bottle have none over ℤ); every basis vector is a cycle, its
// boundary zero modulo 2; their number is the Betti number over ℤ/2 of persistence.hpp; a 2-cell
// outside every cycle is outside the support. The 1-cycles of a graph number E − V + c.
// Q-analysis: a family worked by hand; Q_0 is the number of components of the complex the family
// generates; relabelling the vertices changes nothing.

#include <morphe/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/cycles.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/persistence.hpp>
#include <morphe/topology/q_analysis.hpp>
#include <morphe/topology/simplicial.hpp>

#include <map>
#include <string>
#include <vector>

namespace mt = morphe::topology;

namespace {

// Whether a set of k-cells has zero boundary modulo 2.
bool closed(const mt::CWComplex& X, std::size_t k, const std::vector<std::size_t>& cells) {
  std::map<std::size_t, int> parity;
  for (std::size_t c : cells) {
    if (k == 1) {
      const auto& e = X.edges[c];
      if (e[0] != e[1]) { parity[e[0]] ^= 1; parity[e[1]] ^= 1; }
    } else {
      for (const auto& s : X.faces[c].boundary) parity[s.edge] ^= 1;
    }
  }
  for (const auto& [i, p] : parity)
    if (p) return false;
  return true;
}

// The Betti numbers over ℤ/2 from persistence with every cell at stage 0.
std::size_t betti_mod2(const mt::CWComplex& X, std::size_t k) {
  mt::Filtration F;
  F.vertex.assign(X.vertices, 0);
  F.edge.assign(X.edges.size(), 0);
  F.face.assign(X.faces.size(), 0);
  return mt::persistent_homology(X, F).betti(k, 0);
}

}  // namespace

int main() {
  morphe::Suite s;

  for (const auto& [name, X] : std::vector<std::pair<std::string, mt::CWComplex>>{
           {"the torus", mt::orientable_surface(1)}, {"the surface of genus 2", mt::orientable_surface(2)},
           {"the projective plane", mt::real_projective_space(2)}, {"the Klein bottle", mt::klein_bottle()},
           {"the 4 × 3 lattice torus", mt::torus_lattice({4, 3})}}) {
    const mt::Cycles C = mt::cycles_mod2(X, 2);
    bool all_closed = true;
    for (const auto& c : C.basis) all_closed = all_closed && closed(X, 2, c);
    s.require(name + ": one 2-cycle over ℤ/2, the whole surface", C.dimension == 1 && C.basis[0].size() == X.faces.size() && all_closed,
              std::to_string(C.dimension));
    s.equal(name + ": as many 2-cycles as the Betti number over ℤ/2", C.dimension, betti_mod2(X, 2));
  }
  // A torus with a disc glued along a loop, and a disc beside it: the disc is in no 2-cycle.
  {
    mt::CWComplex X = mt::torus_lattice({3, 3});
    const std::size_t v = X.vertices;
    X.vertices += 2;
    X.edges.push_back({v, v + 1});
    X.edges.push_back({v + 1, v});
    mt::Face d;
    d.base = v;
    d.boundary = {{X.edges.size() - 2, 1}, {X.edges.size() - 1, 1}};
    X.faces.push_back(d);
    const mt::Cycles C = mt::cycles_mod2(X, 2);
    s.require("a torus beside a disc: one 2-cycle, and the disc outside its support", C.dimension == 1 && !C.support.back() && C.support.front());
    const mt::Cycles L = mt::cycles_mod2(X, 1);
    bool all_closed = true;
    for (const auto& c : L.basis) all_closed = all_closed && closed(X, 1, c);
    s.require("its 1-cycles: E − V + c of them, each closed", L.dimension == X.edges.size() - X.vertices + 2 && all_closed,
              std::to_string(L.dimension));
  }

  // Q-analysis.
  {
    // abc, bcd, de, f: 0-connected {abc, bcd, de}, {f}; at q = 1, de stands alone; at q = 2, abc and bcd.
    const std::vector<mt::Simplex> family{{0, 1, 2}, {1, 2, 3}, {3, 4}, {5}};
    s.equal("Q-analysis: abc, bcd, de, f has Q = (2, 2, 2)", mt::q_analysis(family), std::vector<std::size_t>{2, 2, 2});
    const std::vector<mt::Simplex> chain{{0, 1, 2, 3}, {2, 3, 4, 5}, {4, 5, 6, 7}, {0, 7}};
    s.equal("Q-analysis: a chain of tetrahedra sharing edges, closed by an edge, has Q = (1, 2, 3, 3)", mt::q_analysis(chain), std::vector<std::size_t>{1, 2, 3, 3});
    const std::vector<mt::Simplex> relabelled{{7, 6, 5, 4}, {5, 4, 3, 2}, {3, 2, 1, 0}, {7, 0}};
    s.equal("Q-analysis: unchanged by relabelling the vertices", mt::q_analysis(relabelled), mt::q_analysis(chain));
    const auto K = mt::SimplicialComplex::from_facets(family);
    s.equal("Q-analysis: Q_0 is the number of components of the complex", mt::q_analysis(family)[0], mt::betti(K.cw())[0]);
  }
  return s.done();
}
