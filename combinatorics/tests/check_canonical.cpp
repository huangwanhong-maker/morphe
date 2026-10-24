// Checks of the canonical form of labelled multigraphs (canonical.hpp).
//
// Automorphism groups known by hand: the cycle C_n (dihedral, 2n), the complete graph K_4 (S_4),
// the star K_{1,3} (S_3), the Petersen graph (S_5, order 120, where colour refinement splits
// nothing); labels and parallel edges break symmetry as they should. Isomorphic graphs, relabelled
// at random, share their code; the hexagon and two triangles, which colour refinement cannot tell
// apart, do not.

#include <morphe/check.hpp>
#include <morphe/combinatorics/canonical.hpp>

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace mc = morphe::combinatorics;

namespace {

mc::LabelledGraph plain(std::size_t n, const std::vector<std::pair<std::size_t, std::size_t>>& edges) {
  mc::LabelledGraph G;
  G.vertex.assign(n, "");
  for (const auto& [u, v] : edges) G.edges.push_back({u, v, ""});
  return G;
}

mc::LabelledGraph relabelled(const mc::LabelledGraph& G, const std::vector<std::size_t>& p) {
  mc::LabelledGraph H;
  H.vertex.assign(G.vertex.size(), "");
  for (std::size_t v = 0; v < G.vertex.size(); v++) H.vertex[p[v]] = G.vertex[v];
  for (const auto& [u, v, l] : G.edges) H.edges.push_back({p[u], p[v], l});
  return H;
}

}  // namespace

int main() {
  morphe::Suite s;
  for (std::size_t n = 3; n <= 8; n++) {
    std::vector<std::pair<std::size_t, std::size_t>> e;
    for (std::size_t i = 0; i < n; i++) e.push_back({i, (i + 1) % n});
    s.equal("canonical: the cycle C_" + std::to_string(n) + " has 2n automorphisms", mc::canonical_form(plain(n, e)).automorphisms, 2 * n);
  }
  s.equal("canonical: K_4 has 24 automorphisms", mc::canonical_form(plain(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}})).automorphisms, std::size_t{24});
  s.equal("canonical: the star K_{1,3} has 6", mc::canonical_form(plain(4, {{0, 1}, {0, 2}, {0, 3}})).automorphisms, std::size_t{6});
  const auto petersen = plain(10, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9}, {5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5}});
  const auto P = mc::canonical_form(petersen, 4000000);
  s.require("canonical: the Petersen graph has 120 automorphisms", P.exact && P.automorphisms == 120, std::to_string(P.automorphisms));
  // Labels and multiplicities.
  mc::LabelledGraph path = plain(3, {{0, 1}, {1, 2}});
  path.vertex = {"a", "x", "b"};
  s.equal("canonical: a path with ends labelled apart has 1", mc::canonical_form(path).automorphisms, std::size_t{1});
  mc::LabelledGraph twice = plain(3, {{0, 1}, {0, 1}, {1, 2}});
  s.equal("canonical: a doubled edge breaks the path's symmetry", mc::canonical_form(twice).automorphisms, std::size_t{1});
  mc::LabelledGraph coloured = plain(3, {{0, 1}, {1, 2}, {2, 0}});
  std::get<2>(coloured.edges[0]) = "red";
  s.equal("canonical: a triangle with one red edge has 2", mc::canonical_form(coloured).automorphisms, std::size_t{2});
  // Isomorphism, and its absence where refinement is blind.
  std::mt19937 rng(20261024);
  bool same = true;
  for (int trial = 0; trial < 40; trial++) {
    std::vector<std::size_t> p(10);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);
    same = same && mc::canonical_form(relabelled(petersen, p), 4000000).code == P.code;
  }
  s.require("canonical: forty relabellings of the Petersen graph share its code", same);
  const auto hexagon = plain(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
  const auto triangles = plain(6, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}});
  s.require("canonical: the hexagon and two triangles, alike to colour refinement, differ", mc::canonical_form(hexagon).code != mc::canonical_form(triangles).code);
  s.equal("canonical: two triangles have 72 automorphisms", mc::canonical_form(triangles).automorphisms, std::size_t{72});
  return s.done();
}
