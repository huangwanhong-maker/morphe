// Coloured graphs and the pseudomanifolds they encode (crystallisations).
//
// Definition (Pezzana 1974; Ferri, Gagliardi and Grasselli 1986; Gurau 2010). A
// (D+1)-coloured graph is a (D+1)-regular graph whose edges carry colours 0..D such that the
// edges at each vertex have distinct colours. Here the graph is bipartite, with p black and
// p white vertices, and is given by one permutation σ_c of {0..p-1} per colour: black b is
// joined to white σ_c(b) by an edge of colour c. Vertices: black b is b, white w is p + w;
// edges: (c, b) is edge c p + b, from black b to white σ_c(b).
//
// Such a graph is the Feynman graph of a coloured tensor model or group field theory in
// dimension D (Gurau 2011; Oriti 2009), and it encodes a D-dimensional orientable normal
// pseudomanifold: each vertex is a D-simplex whose vertices are coloured 0..D (the facet
// opposite vertex c has colour c), black and white simplices oppositely oriented, and each
// edge of colour c glues two simplices along their facets of colour c. A k-cell of the
// result, with vertex colours T (|T| = k + 1), is a connected component of the subgraph of
// the D - k colours outside T (a "bubble"); its boundary is Σ_i (-1)^i over the faces
// T \ {t_i}, the cells containing it. In particular the (D-2)-cells are the faces of the
// graph, its bicoloured cycles, whose number for colours i, j is the number of cycles of
// σ_j⁻¹ σ_i.
//
// The 2-complex of the graph has its vertices and edges, and a 2-cell bounded by each
// bicoloured cycle: it is the dual 2-skeleton of the pseudomanifold, on which the graph's
// amplitude in a group field theory is a lattice gauge theory (Boulatov 1992, Ooguri 1992).
//
// A jacket is a cyclic order of the colours up to reversal, (D)!/2 of them (D >= 2); the
// graph with the faces of consecutive colours is a ribbon graph embedded in an orientable
// surface of genus g_J. The degree ω = Σ_J g_J governs the large-N expansion of tensor
// models (Gurau 2011): the number of faces is F = D(D-1)p/2 + D - 2ω/(D-1)! for a connected
// graph, and ω = 0 exactly for the melonic graphs (D >= 3; Bonzom, Gurau, Riello and
// Rivasseau 2011), counted with a root by the Fuss–Catalan numbers.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "cw_complex.hpp"

namespace morphe::topology {

struct ColouredGraph {
  std::size_t colours = 0;  // D + 1
  std::size_t p = 0;        // black vertices (and white)
  std::vector<std::vector<std::size_t>> sigma;  // sigma[c][b]: the white neighbour of black b along colour c

  std::size_t dimension() const { return colours - 1; }
  std::size_t vertices() const { return 2 * p; }
  std::size_t edge(std::size_t c, std::size_t b) const { return c * p + b; }

  void validate() const {
    if (colours < 2) throw std::invalid_argument("morphe::topology::ColouredGraph: at least two colours");
    if (sigma.size() != colours) throw std::invalid_argument("morphe::topology::ColouredGraph: one permutation per colour");
    for (const auto& s : sigma) {
      if (s.size() != p) throw std::invalid_argument("morphe::topology::ColouredGraph: permutations of the p black vertices");
      std::vector<char> seen(p, 0);
      for (auto w : s) {
        if (w >= p || seen[w]++) throw std::invalid_argument("morphe::topology::ColouredGraph: sigma is not a permutation");
      }
    }
  }

  // The components of the subgraph with the colours in `mask`: a label per vertex.
  std::vector<std::size_t> components(std::uint32_t mask, std::size_t* count = nullptr) const {
    std::vector<std::size_t> parent(2 * p);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](std::size_t x) {
      while (parent[x] != x) x = parent[x] = parent[parent[x]];
      return x;
    };
    for (std::size_t c = 0; c < colours; c++)
      if (mask & (1u << c))
        for (std::size_t b = 0; b < p; b++) parent[root(b)] = root(p + sigma[c][b]);
    std::vector<std::size_t> label(2 * p, SIZE_MAX), out(2 * p);
    std::size_t n = 0;
    for (std::size_t v = 0; v < 2 * p; v++) {
      const std::size_t r = root(v);
      if (label[r] == SIZE_MAX) label[r] = n++;
      out[v] = label[r];
    }
    if (count) *count = n;
    return out;
  }
  bool connected() const {
    std::size_t n;
    components((1u << colours) - 1, &n);
    return n == 1;
  }

  // The bicoloured cycles of colours i < j, each as its black vertices in order.
  std::vector<std::vector<std::size_t>> faces(std::size_t i, std::size_t j) const {
    std::vector<std::size_t> inv_j(p);
    for (std::size_t b = 0; b < p; b++) inv_j[sigma[j][b]] = b;
    std::vector<char> seen(p, 0);
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t b0 = 0; b0 < p; b0++) {
      if (seen[b0]) continue;
      std::vector<std::size_t> cyc;
      for (std::size_t b = b0; !seen[b]; b = inv_j[sigma[i][b]]) seen[b] = 1, cyc.push_back(b);
      out.push_back(cyc);
    }
    return out;
  }
  std::size_t face_count() const {
    std::size_t f = 0;
    for (std::size_t i = 0; i < colours; i++)
      for (std::size_t j = i + 1; j < colours; j++) f += faces(i, j).size();
    return f;
  }
  // The number of D-bubbles: components with one colour removed, over all colours.
  std::size_t bubbles() const {
    std::size_t total = 0;
    for (std::size_t c = 0; c < colours; c++) {
      std::size_t n;
      components(((1u << colours) - 1) & ~(1u << c), &n);
      total += n;
    }
    return total;
  }

  // The genus of each jacket (D >= 2), colour 0 first in each cyclic order.
  std::vector<std::size_t> jacket_genera() const {
    std::vector<std::size_t> out;
    if (colours < 3) return out;
    std::size_t comps;
    components((1u << colours) - 1, &comps);
    std::vector<std::size_t> rest(colours - 1);
    std::iota(rest.begin(), rest.end(), 1);
    do {
      if (rest.front() > rest.back()) continue;  // each cyclic order once, up to reversal
      std::vector<std::size_t> order{0};
      order.insert(order.end(), rest.begin(), rest.end());
      long long F = 0;
      for (std::size_t k = 0; k < colours; k++) {
        const std::size_t a = order[k], b = order[(k + 1) % colours];
        F += static_cast<long long>(faces(std::min(a, b), std::max(a, b)).size());
      }
      const long long chi = static_cast<long long>(2 * p) - static_cast<long long>(colours * p) + F;
      out.push_back(static_cast<std::size_t>((2 * static_cast<long long>(comps) - chi) / 2));
    } while (std::next_permutation(rest.begin(), rest.end()));
    return out;
  }
  // Gurau's degree, the sum of the jacket genera.
  std::size_t degree() const {
    const auto g = jacket_genera();
    return std::accumulate(g.begin(), g.end(), std::size_t{0});
  }

  // The 2-complex: the graph with a 2-cell on each bicoloured cycle, whose word runs along
  // colour i from black to white and back along colour j.
  CWComplex two_complex() const {
    validate();
    CWComplex X;
    X.name = "2-complex of a " + std::to_string(colours) + "-coloured graph";
    X.vertices = 2 * p;
    for (std::size_t c = 0; c < colours; c++)
      for (std::size_t b = 0; b < p; b++) {
        X.edges.push_back({b, p + sigma[c][b]});
        X.edge_names.push_back("c" + std::to_string(c) + "." + std::to_string(b + 1));
      }
    for (std::size_t i = 0; i < colours; i++)
      for (std::size_t j = i + 1; j < colours; j++)
        for (const auto& cyc : faces(i, j)) {
          Face f{cyc[0], {}};
          for (std::size_t k = 0; k < cyc.size(); k++) {
            f.boundary.push_back({edge(i, cyc[k]), 1});
            f.boundary.push_back({edge(j, cyc[(k + 1) % cyc.size()]), -1});
          }
          X.faces.push_back(f);
        }
    return X;
  }

  // The pseudomanifold as a CW complex (a Δ-complex): k-cells are the components of the
  // subgraphs on D - k colours.
  CWComplex dual_complex() const {
    validate();
    const std::size_t D = colours - 1;
    const std::uint32_t all = (1u << colours) - 1;
    // For every colour set S: the component of each vertex, and the index of each component
    // among the cells of dimension D - |S|.
    std::vector<std::vector<std::size_t>> comp(all + 1), index(all + 1);
    std::vector<std::vector<std::size_t>> rep(all + 1);  // a vertex of each component
    std::vector<std::size_t> cells(D + 1, 0);
    for (std::uint32_t S = 0; S <= all; S++) {
      const std::size_t size = static_cast<std::size_t>(__builtin_popcount(S));
      if (size > D) continue;
      std::size_t n;
      comp[S] = components(S, &n);
      const std::size_t k = D - size;
      index[S].resize(n);
      rep[S].assign(n, SIZE_MAX);
      for (std::size_t c = 0; c < n; c++) index[S][c] = cells[k]++;
      for (std::size_t v = 0; v < 2 * p; v++)
        if (rep[S][comp[S][v]] == SIZE_MAX) rep[S][comp[S][v]] = v;
    }
    // The cell of colour set S containing vertex v.
    auto cell = [&](std::uint32_t S, std::size_t v) { return index[S][comp[S][v]]; };
    auto colours_of = [&](std::uint32_t T) {
      std::vector<std::size_t> t;
      for (std::size_t c = 0; c < colours; c++)
        if (T & (1u << c)) t.push_back(c);
      return t;
    };
    CWComplex X;
    X.name = "pseudomanifold of a " + std::to_string(colours) + "-coloured graph";
    X.vertices = cells[0];
    X.edges.resize(D >= 1 ? cells[1] : 0);
    X.faces.resize(D >= 2 ? cells[2] : 0);
    for (std::size_t k = 3; k <= D; k++) X.higher.emplace_back(cells[k]);
    for (std::uint32_t S = 0; S <= all; S++) {
      const std::size_t size = static_cast<std::size_t>(__builtin_popcount(S));
      if (size > D) continue;
      const std::size_t k = D - size;
      const std::uint32_t T = all & ~S;
      const auto t = colours_of(T);
      for (std::size_t c = 0; c < index[S].size(); c++) {
        const std::size_t v = rep[S][c], me = index[S][c];
        // The face of T without t_i has colour set S ∪ {t_i}.
        auto face = [&](std::size_t i) { return cell(S | (1u << t[i]), v); };
        if (k == 1) X.edges[me] = {face(1), face(0)};  // ∂[t0 t1] = [t1] - [t0]: from [t0] to [t1]
        else if (k == 2) {
          // From [t0]: [t0 t1], [t1 t2], [t0 t2]^{-1}. [t_a t_b] omits the third colour.
          const std::size_t v0 = cell(S | (1u << t[1]) | (1u << t[2]), v);
          X.faces[me] = {v0, {{face(2), 1}, {face(0), 1}, {face(1), -1}}};
        } else if (k >= 3) {
          algebra::Vector b;
          for (std::size_t i = 0; i <= k; i++) b = algebra::axpy(b, i % 2 ? -1 : 1, algebra::Vector{{face(i), 1}});
          X.higher[k - 3][me] = b;
        }
      }
    }
    return X;
  }
};

// The elementary melon: two vertices joined by D + 1 edges; its pseudomanifold is the
// D-sphere made of two D-simplices.
inline ColouredGraph elementary_melon(std::size_t D) {
  return ColouredGraph{D + 1, 1, std::vector<std::vector<std::size_t>>(D + 1, std::vector<std::size_t>{0})};
}

// Inserts an elementary D-melon on the edge of colour c at black b: the edge b → w becomes
// b → w' and b' → w, with b' and w' joined by every other colour (a melonic move, which
// keeps the pseudomanifold and the degree).
inline ColouredGraph insert_melon(const ColouredGraph& G, std::size_t c, std::size_t b) {
  ColouredGraph H = G;
  const std::size_t q = G.p;  // the new black and the new white
  H.p = q + 1;
  for (std::size_t d = 0; d < G.colours; d++) H.sigma[d].push_back(q);
  H.sigma[c][q] = G.sigma[c][b];
  H.sigma[c][b] = q;
  return H;
}

// A coloured graph with p black vertices and uniformly random permutations (σ_0 the identity).
inline ColouredGraph random_coloured_graph(std::size_t D, std::size_t p, std::uint64_t seed) {
  std::mt19937_64 g(seed);
  ColouredGraph G{D + 1, p, {}};
  for (std::size_t c = 0; c <= D; c++) {
    std::vector<std::size_t> s(p);
    std::iota(s.begin(), s.end(), 0);
    if (c > 0)
      for (std::size_t i = p; i > 1; i--) {
        // Fisher–Yates with an unbiased draw from {0..i-1}.
        const std::uint64_t limit = UINT64_MAX - UINT64_MAX % i;
        std::uint64_t x;
        do x = g(); while (x >= limit);
        std::swap(s[i - 1], s[x % i]);
      }
    G.sigma.push_back(s);
  }
  return G;
}

}  // namespace morphe::topology
