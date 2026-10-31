// Checks of the rewriting of regular complexes (complex_rewriting.hpp).
//
// On a 3 × 2 grid of blocks (a district: junctions, streets, blocks):
// - a lane through a block as three rewrites — the two streets it joins subdivided, the block split between
//   the new vertices: each derived complex regular, ∂∂ = 0, Euler characteristic 1 and a disc's homology
//   kept, the counts V, E, F as worked by hand, kept cells keeping their ids and the split block replaced
//   by two new ones;
// - its causal order: the subdivisions independent, the split after both (Foata steps {0, 1}, {2}; two
//   records in the trace), and the record in the other order deriving an isomorphic complex (the local
//   Church–Rosser theorem);
// - the inverse rewrites (merge the two blocks, merge at each new vertex) giving back the grid up to
//   isomorphism;
// - the dangling condition: a street that bounds a block cannot be removed, a free one can; a new street
//   between two junctions outside every block adds a 1-cycle;
// - refusals: two 2-cells sharing two edges are not merged; a cube is not cut along two opposite faces (the
//   wall would be bounded by two circles);
// On two cubes: a door relabelled, the cubes merged across their wall into one 3-cell of ten facets and
// partitioned back (the same complex, its new cells oriented oppositely to the product orientation); the matches of a rule's left side, two vertices, are |Aut L| = 2 times its occurrences.

#include <morphe/check.hpp>
#include <morphe/topology/complex_rewriting.hpp>

#include <string>
#include <vector>

namespace mt = morphe::topology;
namespace rw = morphe::formal::rewriting;
namespace tr = morphe::formal::traces;

namespace {

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}
std::vector<std::size_t> betti(const mt::RegularComplex& X) { return morphe::algebra::betti(mt::chain_complex(X)); }

}  // namespace

int main() {
  morphe::Suite s;

  const auto grid = mt::cubical_grid(3, 2);
  const mt::RegularComplex D = grid.complex;
  s.equal("the district: 12 junctions, 17 streets, 6 blocks", text(D.counts()), std::string("12 17 6"));

  // ---- a lane through block (1, 0), from its bottom street to its top street ----
  const std::size_t block = grid.at({1, 0, 0}, {1, 1, 0});
  const std::size_t bottom = grid.at({1, 0, 0}, {1, 0, 0}), top = grid.at({1, 1, 0}, {1, 0, 0});
  const std::uint64_t block_id = D.cells[block].id, bottom_id = D.cells[bottom].id, top_id = D.cells[top].id;
  const auto r1 = mt::rewrite(D, mt::subdivide_edge(D, bottom, "portal"));
  s.require("the bottom street subdivided", r1.applied, r1.reason);
  const auto r2 = mt::rewrite(r1.result, mt::subdivide_edge(r1.result, r1.result.find(top_id), "portal"));
  s.require("the top street subdivided (it bounds two blocks)", r2.applied, r2.reason);
  const auto& X2 = r2.result;
  const std::size_t m1 = X2.find(r1.created[0]), m2 = X2.find(r2.created[0]);
  s.require("each subdivision creates a vertex first", X2.cells[m1].dim == 0 && X2.cells[m2].dim == 0 && X2.cells[m1].type == "portal");
  const auto r3 = mt::rewrite(X2, mt::split_face(X2, X2.find(block_id), m1, m2, "lane"));
  s.require("the block split by the lane", r3.applied, r3.reason);
  const auto& X3 = r3.result;
  s.equal("then 14 junctions, 20 streets, 7 blocks", text(X3.counts()), std::string("14 20 7"));
  s.require("each step regular, Euler characteristic 1, a disc's homology",
            mt::admissible(r1.result).regular && mt::admissible(X2).regular && mt::admissible(X3).regular && X3.euler_characteristic() == 1 &&
                text(betti(r1.result)) == "1 0 0" && text(betti(X2)) == "1 0 0" && text(betti(X3)) == "1 0 0");
  s.require("the split block is deleted and two blocks and a lane are created", X3.find(block_id) == mt::no_cell && r3.deleted == std::vector<std::uint64_t>{block_id} && r3.created.size() == 3);
  s.require("the subdivided streets are gone, the other cells keep their ids", X3.find(bottom_id) == mt::no_cell && X3.find(top_id) == mt::no_cell && X3.find(D.cells[0].id) != mt::no_cell);
  bool traced = true;
  for (std::size_t i = 0; i < X2.size(); i++)
    if (r3.trace[i] != mt::no_cell) traced = traced && X3.cells[r3.trace[i]].id == X2.cells[i].id;
  s.require("the trace follows each kept cell", traced);

  // ---- the causal order of the three rewrites ----
  const auto P = tr::causal_order({r1.footprint, r2.footprint, r3.footprint});
  s.require("the two subdivisions are independent", P.concurrent(0, 1));
  s.require("the split depends on both", P.before[0][2] && P.before[1][2]);
  s.require("Foata steps: {0, 1}, then {2}", P.steps.size() == 2 && P.steps[0] == std::vector<std::size_t>{0, 1} && P.steps[1] == std::vector<std::size_t>{2});
  s.equal("two records in the trace", tr::linearisations(P).str(), std::string("2"));
  const auto q1 = mt::rewrite(D, mt::subdivide_edge(D, top, "portal"));
  const auto q2 = mt::rewrite(q1.result, mt::subdivide_edge(q1.result, q1.result.find(bottom_id), "portal"));
  const auto& Q2 = q2.result;
  const auto q3 = mt::rewrite(Q2, mt::split_face(Q2, Q2.find(block_id), Q2.find(q2.created[0]), Q2.find(q1.created[0]), "lane"));
  s.require("the other record derives an isomorphic complex (local Church–Rosser)", q3.applied && mt::isomorphic(q3.result, X3, true));

  // ---- the inverse rewrites ----
  {
    const auto lane = X3.find(r3.created[0]);
    s.require("the lane is the cut, a 1-cell", X3.cells[lane].dim == 1 && X3.cells[lane].type == "lane");
    const auto u1 = mt::rewrite(X3, mt::merge_cells(X3, lane));
    s.require("the two blocks merged across the lane", u1.applied && text(u1.result.counts()) == "14 19 6", u1.reason);
    const auto u2 = mt::rewrite(u1.result, mt::merge_at_vertex(u1.result, u1.result.find(r1.created[0])));
    const auto u3 = mt::rewrite(u2.result, mt::merge_at_vertex(u2.result, u2.result.find(r2.created[0])));
    s.require("and the streets merged at the portals", u2.applied && u3.applied, u2.reason + u3.reason);
    s.require("giving back the district up to isomorphism", mt::isomorphic(u3.result, D));
  }

  // ---- the dangling condition ----
  {
    const auto r = mt::rewrite(D, mt::remove_cell(D, bottom));
    s.require("a street bounding a block cannot be removed", !r.applied && r.reason.find("gluing condition") != std::string::npos, r.reason);
    const std::size_t a = grid.at({0, 0, 0}, {0, 0, 0}), b = grid.at({3, 2, 0}, {0, 0, 0});
    const auto e = mt::rewrite(D, mt::attach_edge(D, a, b, "street"));
    s.require("a new street between two corners outside every block", e.applied, e.reason);
    s.equal("adds a 1-cycle: Betti 1 1 0", text(betti(e.result)), std::string("1 1 0"));
    const auto back = mt::rewrite(e.result, mt::remove_cell(e.result, e.result.find(e.created[0])));
    s.require("and, bounding nothing, can be removed again", back.applied && mt::isomorphic(back.result, D), back.reason);
  }

  // ---- refusals ----
  {
    mt::RegularComplex X;
    const auto a = X.add_vertex(), b = X.add_vertex(), c = X.add_vertex();
    const auto ab = X.add_edge(a, b), bc = X.add_edge(b, c), ca = X.add_edge(c, a), ab2 = X.add_edge(a, b);
    X.add(2, "", {ab, bc, ca});
    X.add(2, "", {ab2, bc, ca});
    mt::orient(X);
    const auto r = mt::rewrite(X, mt::merge_cells(X, bc));
    s.require("two 2-cells sharing two edges are not merged", !r.applied && r.reason.find("share another facet") != std::string::npos, r.reason);
    auto C = mt::cubical_grid(1, 1, 1);
    const auto cube = C.complex.size() - 1;
    const auto lo = C.at({0, 0, 0}, {1, 1, 0}), hi = C.at({0, 0, 1}, {1, 1, 0});
    const auto p = mt::rewrite(C.complex, mt::partition_cell(C.complex, cube, {lo, hi}, "wall"));
    s.require("a cube is not cut along two opposite faces", !p.applied && p.reason.find("not a circle") != std::string::npos, p.reason);
  }

  // ---- two cubes ----
  {
    const auto G = mt::cubical_grid(2, 1, 1);
    const auto& X = G.complex;
    const std::size_t wall = G.at({1, 0, 0}, {0, 1, 1}), c1 = G.at({0, 0, 0}, {1, 1, 1});
    const auto d = mt::rewrite(X, mt::relabel(X, wall, "door"));
    s.require("the wall relabelled a door, nothing else changed", d.applied && d.result.cells[d.trace[wall]].type == "door" && d.relabelled.size() == 1 && d.created.empty() && d.deleted.empty());
    const auto m = mt::rewrite(X, mt::merge_cells(X, wall));
    s.require("the two cubes merged across their wall", m.applied, m.reason);
    const auto& M = m.result;
    const std::size_t room = M.find(m.created[0]);
    s.require("into one 3-cell of ten facets, regular, a ball", M.cells[room].dim == 3 && M.cells[room].facets.size() == 10 && mt::admissible(M).regular && text(betti(M)) == "1 0 0 0");
    std::vector<std::size_t> part;
    for (std::size_t f : X.cells[c1].facets)
      if (f != wall) part.push_back(m.trace[f]);
    const auto back = mt::rewrite(M, mt::partition_cell(M, room, part, "yz"));
    s.require("partitioned back along the first cube's faces: the same complex, oppositely oriented", back.applied && mt::isomorphic(back.result, X) && !mt::isomorphic(back.result, X, true), back.reason);

    const auto attach = mt::attach_edge(X, 0, 1, "");
    const auto& L = attach.rule.left;
    const auto H = mt::incidence_graph(X);
    s.equal("the matches of two vertices, 12 × 11, are |Aut L| = 2 times the occurrences", rw::count_matches(L, H), std::size_t{2} * 66);
    s.equal("Aut L of two vertices", rw::automorphisms(L), std::size_t{2});
  }

  return s.done();
}
