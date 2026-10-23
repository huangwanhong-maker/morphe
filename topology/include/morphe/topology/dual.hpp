// The dual 2-complex of a triangulated pseudomanifold.
//
// Let K be a pure d-dimensional simplicial pseudomanifold, d >= 2. Its dual 2-complex K*
// (Rourke and Sanderson 1972, §2.21 for dual cells; Reisenberger and Rovelli 1997 and
// Oriti 2001 for its role in spin foams) has
// - a vertex for each d-simplex,
// - an edge for each (d-1)-simplex shared by two d-simplices, from the lower-numbered to the
//   higher-numbered,
// - a 2-cell for each (d-2)-simplex h around which the d-simplices close up: those containing
//   h, joined across the (d-1)-simplices containing h, form cycles and paths (each d-simplex
//   meets h in exactly two of its (d-1)-faces), and every cycle bounds a 2-cell whose word
//   runs around it.
// For a closed combinatorial manifold this is the 2-skeleton of the dual cell decomposition;
// it is the 2-complex on which spin-foam models and lattice gauge theories on triangulations
// live, its 2-cells carrying the curvature (deficit angle, holonomy) of the hinge they dual.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "cw_complex.hpp"
#include "simplicial.hpp"

namespace morphe::topology {

struct DualTwoComplex {
  CWComplex complex;
  std::vector<Simplex> vertex_simplex;  // the d-simplex of each dual vertex
  std::vector<Simplex> edge_simplex;    // the (d-1)-simplex of each dual edge
  std::vector<Simplex> face_simplex;    // the (d-2)-simplex of each dual 2-cell
};

inline DualTwoComplex dual_two_complex(const SimplicialComplex& K) {
  const std::size_t d = K.dimension();
  if (d < 2 || !K.is_pure()) throw std::invalid_argument("morphe::topology::dual_two_complex: a pure complex of dimension >= 2");
  DualTwoComplex D;
  D.complex.name = "dual 2-complex of " + K.name;
  const auto& top = K.simplices[d];
  D.complex.vertices = top.size();
  D.vertex_simplex = top;
  // The d-simplices containing each (d-1)-simplex.
  std::vector<std::vector<std::size_t>> cof(K.count(d - 1));
  for (std::size_t t = 0; t < top.size(); t++)
    for (std::size_t i = 0; i <= d; i++) {
      Simplex f = top[t];
      f.erase(f.begin() + static_cast<long>(i));
      cof[K.find(f)].push_back(t);
    }
  std::vector<std::size_t> edge_of(K.count(d - 1), SIZE_MAX);
  for (std::size_t f = 0; f < cof.size(); f++) {
    if (cof[f].size() > 2) throw std::invalid_argument("morphe::topology::dual_two_complex: not a pseudomanifold");
    if (cof[f].size() == 2) {
      edge_of[f] = D.complex.edges.size();
      D.complex.edges.push_back({cof[f][0], cof[f][1]});  // ascending, as top simplices are visited in order
      D.edge_simplex.push_back(K.simplices[d - 1][f]);
    }
  }
  // Around each (d-2)-simplex h: the d-simplices containing h and their links across the
  // interior (d-1)-simplices containing h.
  for (const auto& h : K.simplices[d - 2]) {
    std::vector<std::size_t> around;  // dual edges around h
    for (const auto& f : K.cofaces(h)) {
      const std::size_t fi = K.find(f);
      if (edge_of[fi] != SIZE_MAX) around.push_back(edge_of[fi]);
    }
    // Adjacency of dual vertices through these edges; every vertex has degree <= 2.
    std::vector<char> used(around.size(), 0);
    for (std::size_t start = 0; start < around.size(); start++) {
      if (used[start]) continue;
      // Walk from one end of `start` along unused edges until returning (a cycle) or stopping.
      Path word;
      const std::size_t v0 = D.complex.edges[around[start]][0];
      std::size_t v = v0, e = start;
      bool cycle = false;
      std::vector<std::size_t> taken;
      for (;;) {
        used[e] = 1;
        taken.push_back(e);
        const auto& ed = D.complex.edges[around[e]];
        const int sign = ed[0] == v ? 1 : -1;
        word.push_back({around[e], sign});
        v = sign > 0 ? ed[1] : ed[0];
        if (v == v0) { cycle = true; break; }
        std::size_t next = SIZE_MAX;
        for (std::size_t k = 0; k < around.size(); k++)
          if (!used[k] && (D.complex.edges[around[k]][0] == v || D.complex.edges[around[k]][1] == v)) { next = k; break; }
        if (next == SIZE_MAX) break;
        e = next;
      }
      if (!cycle) {
        // A path: h lies on the boundary; extend backward from v0 to mark the rest as used.
        for (bool more = true; more;) {
          more = false;
          for (std::size_t k = 0; k < around.size(); k++)
            if (!used[k]) {
              const auto& ed = D.complex.edges[around[k]];
              for (std::size_t t : taken) {
                const auto& et = D.complex.edges[around[t]];
                if (ed[0] == et[0] || ed[0] == et[1] || ed[1] == et[0] || ed[1] == et[1]) { used[k] = 1, taken.push_back(k), more = true; break; }
              }
            }
        }
        continue;
      }
      D.complex.faces.push_back({v0, word});
      D.face_simplex.push_back(h);
    }
  }
  return D;
}

}  // namespace morphe::topology
