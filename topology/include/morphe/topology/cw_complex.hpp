// CW complexes, held combinatorially: cells with their attaching data.
//
// Definition (Whitehead 1949; Hatcher 2002, Algebraic Topology, ch. 0 and §2.2). A CW
// complex is built by attaching cells of increasing dimension along maps of their boundary
// spheres into the part already built. Up to homotopy, and for everything computed here,
// the attaching maps enter only through:
//
// - for a 1-cell, its two endpoints: an edge e from tail(e) to head(e), a loop when equal;
// - for a 2-cell, a closed edge path, its attaching word: a sequence of steps e^{±1}, each
//   starting where the last ended, from a base vertex back to it. The empty word attaches a
//   2-cell to its base vertex as a sphere. The word gives the fundamental group (each 2-cell
//   a relator) and the holonomy of a lattice gauge field around the cell;
// - for a k-cell, k >= 3, its cellular boundary: integer incidences [e_k : e_{k-1}], the
//   degrees of the attaching map onto the (k-1)-cells (Hatcher, §2.2, the cellular boundary
//   formula). These are what homology needs; the homotopy type in dimensions >= 3 is not
//   recorded.
//
// The cellular chain complex has ∂(edge) = head - tail, ∂(2-cell) = the exponent sum of its
// word on the edges (its abelianisation), and the given incidences above. validate() checks
// that each word is a closed path from its base and that ∂∂ = 0.

#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/algebra/chain_complex.hpp>
#include <morphe/combinatorics/graph.hpp>

namespace morphe::topology {

using combinatorics::Step;
using Path = std::vector<Step>;
using algebra::Chain;

struct Face {
  std::size_t base = 0;  // the vertex where the attaching word starts and ends
  Path boundary;         // the attaching word, a closed edge path
};

struct CWComplex {
  std::string name;
  std::size_t vertices = 0;
  std::vector<std::array<std::size_t, 2>> edges;  // (tail, head)
  std::vector<Face> faces;
  std::vector<std::vector<algebra::Vector>> higher;  // higher[k - 3][i]: the boundary of the i-th k-cell
  std::vector<std::string> edge_names;  // optional labels of the edges (generators of words)

  std::size_t dimension() const {
    for (std::size_t k = higher.size(); k-- > 0;)
      if (!higher[k].empty()) return k + 3;
    if (!faces.empty()) return 2;
    if (!edges.empty()) return 1;
    return 0;
  }
  // The number of cells of dimension k.
  std::size_t cells(std::size_t k) const {
    switch (k) {
      case 0: return vertices;
      case 1: return edges.size();
      case 2: return faces.size();
      default: return k - 3 < higher.size() ? higher[k - 3].size() : 0;
    }
  }
  std::vector<std::size_t> counts() const {
    std::vector<std::size_t> c;
    for (std::size_t k = 0; k <= dimension(); k++) c.push_back(cells(k));
    return c;
  }
  long long euler_characteristic() const {
    long long chi = 0;
    for (std::size_t k = 0; k <= dimension(); k++) chi += (k % 2 ? -1 : 1) * static_cast<long long>(cells(k));
    return chi;
  }
  std::string edge_name(std::size_t e) const { return e < edge_names.size() && !edge_names[e].empty() ? edge_names[e] : "e" + std::to_string(e + 1); }

  combinatorics::Multigraph one_skeleton() const { return {vertices, edges}; }

  // The vertex at which a path ends, starting from `from`; throws if a step does not start
  // where the last one ended.
  std::size_t walk(std::size_t from, const Path& p) const {
    std::size_t v = from;
    for (const auto& s : p) {
      if (s.edge >= edges.size() || (s.sign != 1 && s.sign != -1)) throw std::invalid_argument("morphe::topology::CWComplex: a step is not an edge with sign ±1");
      const std::size_t start = s.sign > 0 ? edges[s.edge][0] : edges[s.edge][1];
      if (start != v) throw std::invalid_argument("morphe::topology::CWComplex: a word is not a path (step on " + edge_name(s.edge) + ")");
      v = s.sign > 0 ? edges[s.edge][1] : edges[s.edge][0];
    }
    return v;
  }

  // The abelianised word: the cellular boundary of a 2-cell.
  static Chain abelianise(const Path& p) {
    Chain c;
    for (const auto& s : p) c = algebra::axpy(c, s.sign, Chain{{s.edge, 1}});
    return c;
  }

  algebra::ChainComplex chain_complex() const {
    const std::size_t d = dimension();
    algebra::ChainComplex C;
    C.n = counts();
    C.boundary.resize(d + 1);
    if (d >= 1) {
      C.boundary[1] = algebra::SparseMatrix(vertices, edges.size());
      for (std::size_t e = 0; e < edges.size(); e++)
        if (edges[e][0] != edges[e][1]) {
          C.boundary[1].columns[e] = edges[e][0] < edges[e][1] ? Chain{{edges[e][0], -1}, {edges[e][1], 1}} : Chain{{edges[e][1], 1}, {edges[e][0], -1}};
        }
    }
    if (d >= 2) {
      C.boundary[2] = algebra::SparseMatrix(edges.size(), faces.size());
      for (std::size_t f = 0; f < faces.size(); f++) C.boundary[2].columns[f] = abelianise(faces[f].boundary);
    }
    for (std::size_t k = 3; k <= d; k++) {
      C.boundary[k] = algebra::SparseMatrix(cells(k - 1), cells(k));
      for (std::size_t i = 0; i < cells(k); i++) {
        algebra::Vector v;
        for (const auto& e : higher[k - 3][i]) v = algebra::axpy(v, e.value, algebra::Vector{{e.index, 1}});
        C.boundary[k].columns[i] = v;
      }
    }
    return C;
  }

  void validate() const {
    for (const auto& e : edges)
      if (e[0] >= vertices || e[1] >= vertices) throw std::invalid_argument("morphe::topology::CWComplex: an edge leaves the vertices");
    for (std::size_t f = 0; f < faces.size(); f++) {
      if (faces[f].base >= vertices) throw std::invalid_argument("morphe::topology::CWComplex: a 2-cell's base is not a vertex");
      if (walk(faces[f].base, faces[f].boundary) != faces[f].base)
        throw std::invalid_argument("morphe::topology::CWComplex: the word of 2-cell " + std::to_string(f + 1) + " is not closed");
    }
    chain_complex().validate();
  }
};

}  // namespace morphe::topology
