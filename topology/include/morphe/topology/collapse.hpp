// Elementary collapses of 2-dimensional CW complexes (simple homotopy).
//
// Definition (Whitehead 1939; Cohen 1973, A Course in Simple-Homotopy Theory, §4). In a CW
// complex of dimension at most 2, an elementary collapse removes
// - a 2-cell f with a free edge e: e occurs exactly once in the attaching words of all
//   2-cells, and that occurrence is in f; or
// - an edge e with a free vertex v: v is an end of e only, e is not a loop, no 2-cell's word
//   passes through e, and no 2-cell is attached at v.
// Each removes a pair of cells whose union deformation-retracts onto the rest, so the result
// is homotopy equivalent (indeed simple-homotopy equivalent) to the complex. A complex that
// collapses to a single vertex is collapsible, hence contractible; the dunce hat and Bing's
// house are contractible and not collapsible (Zeeman 1964).
//
// The greedy order decides collapsibility in dimension 2: a 2-cell with a free edge keeps it
// as other 2-cells are removed (occurrences only decrease), so every maximal sequence removes
// the same 2-cells; what remains of the 1-skeleton collapses to a point exactly when it is a
// tree. Cost: O(total length of the words + cells) with the bookkeeping below.
// Complexes with cells of dimension 3 or more are not collapsed here.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "cw_complex.hpp"

namespace morphe::topology {

struct Collapse {
  bool attempted = false;    // false for complexes of dimension above 2
  bool collapsible = false;  // whether a single vertex remains
  std::size_t steps = 0;     // elementary collapses made
  std::vector<char> vertex_left, edge_left, face_left;  // what remains
  std::vector<std::size_t> remaining;                   // counts of the cells that remain
  // The collapses in order: (dimension of the larger cell, its index, the index of the free face).
  struct Move { std::size_t dimension, cell, face; };
  std::vector<Move> moves;
};

inline Collapse collapse(const CWComplex& X) {
  Collapse C;
  if (X.dimension() > 2) return C;
  C.attempted = true;
  C.vertex_left.assign(X.vertices, 1);
  C.edge_left.assign(X.edges.size(), 1);
  C.face_left.assign(X.faces.size(), 1);
  // occ[e]: occurrences of e in words of remaining faces; face_of[e]: one such face.
  std::vector<std::size_t> occ(X.edges.size(), 0);
  std::vector<std::vector<std::size_t>> faces_of(X.edges.size());
  for (std::size_t f = 0; f < X.faces.size(); f++)
    for (const auto& s : X.faces[f].boundary) occ[s.edge]++, faces_of[s.edge].push_back(f);
  std::vector<std::size_t> based(X.vertices, 0);  // remaining faces attached at a vertex with an empty word
  for (const auto& f : X.faces)
    if (f.boundary.empty()) based[f.base]++;
  std::vector<std::size_t> degree(X.vertices, 0);  // remaining edge ends at each vertex
  for (const auto& e : X.edges) degree[e[0]]++, degree[e[1]]++;
  // 2-cells first, to a fixed point (their removals only create new free edges).
  std::vector<std::size_t> queue;
  for (std::size_t e = 0; e < X.edges.size(); e++)
    if (occ[e] == 1) queue.push_back(e);
  while (!queue.empty()) {
    const std::size_t e = queue.back();
    queue.pop_back();
    if (!C.edge_left[e] || occ[e] != 1) continue;
    std::size_t f = SIZE_MAX;
    for (std::size_t g : faces_of[e])
      if (C.face_left[g]) f = g;
    if (f == SIZE_MAX) continue;
    C.face_left[f] = 0;
    C.edge_left[e] = 0;
    degree[X.edges[e][0]]--, degree[X.edges[e][1]]--;
    C.moves.push_back({2, f, e});
    C.steps++;
    for (const auto& s : X.faces[f].boundary) {
      occ[s.edge]--;
      if (occ[s.edge] == 1 && C.edge_left[s.edge]) queue.push_back(s.edge);
    }
  }
  // Then edges with a free vertex: leaves of the remaining graph not touched by a 2-cell.
  std::vector<std::vector<std::size_t>> incident(X.vertices);
  for (std::size_t k = 0; k < X.edges.size(); k++) {
    incident[X.edges[k][0]].push_back(k);
    if (X.edges[k][1] != X.edges[k][0]) incident[X.edges[k][1]].push_back(k);
  }
  std::vector<std::size_t> leaves;
  for (std::size_t v = 0; v < X.vertices; v++)
    if (degree[v] == 1) leaves.push_back(v);
  while (!leaves.empty()) {
    const std::size_t v = leaves.back();
    leaves.pop_back();
    if (!C.vertex_left[v] || degree[v] != 1 || based[v]) continue;
    std::size_t e = SIZE_MAX;
    for (std::size_t k : incident[v])
      if (C.edge_left[k] && X.edges[k][0] != X.edges[k][1] && occ[k] == 0) { e = k; break; }
    if (e == SIZE_MAX) continue;
    const std::size_t w = X.edges[e][0] == v ? X.edges[e][1] : X.edges[e][0];
    C.edge_left[e] = 0;
    C.vertex_left[v] = 0;
    degree[v]--, degree[w]--;
    C.moves.push_back({1, e, v});
    C.steps++;
    if (degree[w] == 1) leaves.push_back(w);
  }
  C.remaining = {0, 0, 0};
  for (auto x : C.vertex_left) C.remaining[0] += x;
  for (auto x : C.edge_left) C.remaining[1] += x;
  for (auto x : C.face_left) C.remaining[2] += x;
  C.collapsible = C.remaining[0] == 1 && C.remaining[1] == 0 && C.remaining[2] == 0;
  return C;
}

}  // namespace morphe::topology
