// Abstract simplicial complexes.
//
// Definition (Spanier 1966, Algebraic Topology, §3.1; Munkres 1984, §3). An abstract
// simplicial complex on vertices 0..n-1 is a family of nonempty finite sets of vertices,
// the simplices, closed under taking nonempty subsets. A simplex is held as its sorted
// vertex list; it is oriented by that order, so the simplicial boundary is
// ∂[v_0 .. v_k] = Σ_i (-1)^i [v_0 .. v̂_i .. v_k] (Munkres, §5).
//
// A complex is given by its facets (maximal simplices, or any generating family) and holds
// every face, dimension by dimension, in lexicographic order. The constructions here (links,
// stars, skeleta, the barycentric subdivision, cones, suspensions and joins) are those of
// combinatorial topology (Rourke and Sanderson 1972, Introduction to Piecewise-Linear
// Topology, ch. 2). cw() gives the same space as a CW complex, each triangle's attaching
// word read from its boundary: [v_0 v_1] [v_1 v_2] [v_0 v_2]^{-1}, from v_0.

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "cw_complex.hpp"

namespace morphe::topology {

using Simplex = std::vector<std::size_t>;

struct SimplicialComplex {
  std::string name;
  std::size_t vertices = 0;                     // vertex labels 0..vertices-1
  std::vector<std::vector<Simplex>> simplices;  // simplices[k]: the k-simplices, sorted

  static SimplicialComplex from_facets(std::vector<Simplex> facets, std::size_t vertices = 0, std::string name = "") {
    SimplicialComplex K;
    K.name = std::move(name);
    std::vector<std::set<Simplex>> all;
    std::size_t n = vertices;
    for (auto& f : facets) {
      if (f.empty()) continue;
      std::sort(f.begin(), f.end());
      if (std::adjacent_find(f.begin(), f.end()) != f.end()) throw std::invalid_argument("morphe::topology::SimplicialComplex: a simplex repeats a vertex");
      n = std::max(n, f.back() + 1);
      // Every nonempty subset of f.
      const std::size_t k = f.size();
      if (k > 20) throw std::invalid_argument("morphe::topology::SimplicialComplex: a facet of more than 20 vertices");
      if (all.size() < k) all.resize(k);
      for (unsigned long mask = 1; mask < (1ul << k); mask++) {
        Simplex s;
        for (std::size_t i = 0; i < k; i++)
          if (mask & (1ul << i)) s.push_back(f[i]);
        all[s.size() - 1].insert(std::move(s));
      }
    }
    K.vertices = n;
    if (all.empty()) all.resize(1);
    // Every vertex is a 0-simplex, isolated or not.
    for (std::size_t v = 0; v < n; v++) all[0].insert(Simplex{v});
    for (auto& level : all) K.simplices.emplace_back(level.begin(), level.end());
    return K;
  }

  std::size_t dimension() const { return simplices.empty() ? 0 : simplices.size() - 1; }
  std::size_t count(std::size_t k) const { return k < simplices.size() ? simplices[k].size() : 0; }
  std::vector<std::size_t> counts() const {
    std::vector<std::size_t> c;
    for (const auto& s : simplices) c.push_back(s.size());
    return c;
  }
  // The index of a simplex among those of its dimension, or SIZE_MAX.
  std::size_t find(const Simplex& s) const {
    if (s.empty() || s.size() > simplices.size()) return SIZE_MAX;
    const auto& level = simplices[s.size() - 1];
    const auto it = std::lower_bound(level.begin(), level.end(), s);
    return it != level.end() && *it == s ? static_cast<std::size_t>(it - level.begin()) : SIZE_MAX;
  }
  bool contains(const Simplex& s) const { return find(s) != SIZE_MAX; }
  long long euler_characteristic() const {
    long long chi = 0;
    for (std::size_t k = 0; k < simplices.size(); k++) chi += (k % 2 ? -1 : 1) * static_cast<long long>(simplices[k].size());
    return chi;
  }
  // The simplices contained in no larger simplex.
  std::vector<Simplex> facets() const {
    std::vector<Simplex> out;
    for (std::size_t k = 0; k < simplices.size(); k++)
      for (const auto& s : simplices[k]) {
        bool maximal = true;
        if (k + 1 < simplices.size())
          for (std::size_t v = 0; v < vertices && maximal; v++) {
            if (std::binary_search(s.begin(), s.end(), v)) continue;
            Simplex t = s;
            t.insert(std::upper_bound(t.begin(), t.end(), v), v);
            if (contains(t)) maximal = false;
          }
        if (maximal) out.push_back(s);
      }
    return out;
  }
  bool is_pure() const {
    for (const auto& f : facets())
      if (f.size() != simplices.size()) return false;
    return true;
  }
  // The (k+1)-simplices containing a k-simplex s (its cofaces of codimension 1).
  std::vector<Simplex> cofaces(const Simplex& s) const {
    std::vector<Simplex> out;
    if (s.size() >= simplices.size()) return out;
    for (const auto& t : simplices[s.size()])
      if (std::includes(t.begin(), t.end(), s.begin(), s.end())) out.push_back(t);
    return out;
  }

  // The simplicial chain complex, each simplex oriented by its vertex order.
  algebra::ChainComplex chain_complex() const {
    algebra::ChainComplex C;
    C.n = counts();
    C.boundary.resize(simplices.size());
    for (std::size_t k = 1; k < simplices.size(); k++) {
      C.boundary[k] = algebra::SparseMatrix(simplices[k - 1].size(), simplices[k].size());
      for (std::size_t j = 0; j < simplices[k].size(); j++) {
        algebra::Vector col;
        const Simplex& s = simplices[k][j];
        for (std::size_t i = 0; i < s.size(); i++) {
          Simplex face = s;
          face.erase(face.begin() + static_cast<long>(i));
          col = algebra::axpy(col, i % 2 ? -1 : 1, algebra::Vector{{find(face), 1}});
        }
        C.boundary[k].columns[j] = col;
      }
    }
    return C;
  }

  // The same space as a CW complex: vertices, edges from lower to higher label, triangles
  // with their attaching words, higher simplices by their simplicial boundary.
  CWComplex cw() const {
    CWComplex X;
    X.name = name;
    X.vertices = vertices;
    for (const auto& e : simplices.size() > 1 ? simplices[1] : std::vector<Simplex>{}) X.edges.push_back({e[0], e[1]});
    if (simplices.size() > 2)
      for (const auto& t : simplices[2]) {
        const std::size_t e01 = find({t[0], t[1]}), e12 = find({t[1], t[2]}), e02 = find({t[0], t[2]});
        X.faces.push_back({t[0], {{e01, 1}, {e12, 1}, {e02, -1}}});
      }
    if (simplices.size() > 3) {
      const auto C = chain_complex();
      for (std::size_t k = 3; k < simplices.size(); k++) X.higher.push_back(C.boundary[k].columns);
    }
    return X;
  }

  // The link of a simplex s: the simplices t disjoint from s with s ∪ t in the complex. The
  // labels are kept (vertices is the same bound), and simplices[0] holds only the link's own
  // vertices; renumbered() makes it a complex on 0..m-1.
  SimplicialComplex link(const Simplex& s) const {
    std::vector<std::set<Simplex>> levels;
    for (std::size_t k = s.size(); k < simplices.size(); k++)
      for (const auto& t : simplices[k]) {
        if (!std::includes(t.begin(), t.end(), s.begin(), s.end())) continue;
        Simplex rest;
        std::set_difference(t.begin(), t.end(), s.begin(), s.end(), std::back_inserter(rest));
        if (levels.size() < rest.size()) levels.resize(rest.size());
        levels[rest.size() - 1].insert(rest);
      }
    SimplicialComplex L;
    L.vertices = vertices;
    for (auto& level : levels) L.simplices.emplace_back(level.begin(), level.end());
    return L;
  }
  // The closed star of s: every simplex containing s, with all its faces.
  SimplicialComplex star(const Simplex& s) const {
    std::vector<Simplex> out;
    for (std::size_t k = s.size() - 1; k < simplices.size(); k++)
      for (const auto& t : simplices[k])
        if (std::includes(t.begin(), t.end(), s.begin(), s.end())) out.push_back(t);
    return from_facets(out, vertices);
  }
  // The k-skeleton.
  SimplicialComplex skeleton(std::size_t k) const {
    SimplicialComplex K = *this;
    if (K.simplices.size() > k + 1) K.simplices.resize(k + 1);
    return K;
  }
};

// A complex whose 0-simplices are a subset of its labels (a link) as a complex on 0..m-1.
inline SimplicialComplex renumbered(const SimplicialComplex& K) {
  std::map<std::size_t, std::size_t> label;
  for (const auto& v : K.simplices.empty() ? std::vector<Simplex>{} : K.simplices[0]) label.emplace(v[0], label.size());
  std::vector<Simplex> facets;
  for (const auto& f : K.facets()) {
    Simplex g;
    for (auto v : f) g.push_back(label.at(v));
    facets.push_back(g);
  }
  return SimplicialComplex::from_facets(facets, label.size(), K.name);
}

// The barycentric subdivision: vertices are the simplices of K, simplices the chains
// σ_0 ⊂ σ_1 ⊂ ... ⊂ σ_k (Spanier, §3.3). Vertex numbering: by dimension, then index.
inline SimplicialComplex barycentric_subdivision(const SimplicialComplex& K) {
  std::vector<std::size_t> offset(K.simplices.size() + 1, 0);
  for (std::size_t k = 0; k < K.simplices.size(); k++) offset[k + 1] = offset[k] + K.simplices[k].size();
  std::vector<Simplex> facets;
  // Every maximal chain ends at a facet; enumerate chains downward from each facet.
  for (const auto& f : K.facets()) {
    std::vector<Simplex> chain{f};
    std::function<void(const Simplex&)> down = [&](const Simplex& s) {
      if (s.size() == 1) {
        Simplex flag;
        for (const auto& c : chain) flag.push_back(offset[c.size() - 1] + K.find(c));
        facets.push_back(flag);
        return;
      }
      for (std::size_t i = 0; i < s.size(); i++) {
        Simplex t = s;
        t.erase(t.begin() + static_cast<long>(i));
        chain.push_back(t);
        down(t);
        chain.pop_back();
      }
    };
    down(f);
  }
  return SimplicialComplex::from_facets(facets, offset.back(), "subdivided " + K.name);
}

// The join K * L: simplices σ ⊔ τ (and σ, τ alone); L's vertices shifted past K's.
inline SimplicialComplex join(const SimplicialComplex& K, const SimplicialComplex& L) {
  std::vector<Simplex> facets;
  const auto fk = K.facets(), fl = L.facets();
  for (const auto& a : fk)
    for (const auto& b : fl) {
      Simplex s = a;
      for (auto v : b) s.push_back(v + K.vertices);
      facets.push_back(s);
    }
  return SimplicialComplex::from_facets(facets, K.vertices + L.vertices, K.name + " * " + L.name);
}

// The cone on K (apex: the last vertex) and its suspension (two apexes).
inline SimplicialComplex cone(const SimplicialComplex& K) {
  SimplicialComplex point = SimplicialComplex::from_facets({{0}}, 1, "point");
  SimplicialComplex C = join(K, point);
  C.name = "cone on " + K.name;
  return C;
}
inline SimplicialComplex suspension(const SimplicialComplex& K) {
  SimplicialComplex two = SimplicialComplex::from_facets({{0}, {1}}, 2, "S0");
  SimplicialComplex S = join(K, two);
  S.name = "suspension of " + K.name;
  return S;
}

}  // namespace morphe::topology
