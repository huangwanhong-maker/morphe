// Bubbles: the interactions of invariant theories of a complex field, and the Feynman graphs
// glued from them.
//
// Definition (Gurau 2016, Random Tensors, ch. 2–3; Bonzom, Gurau and Rivasseau 2012, Phys.
// Rev. D 85, 084037). Let φ be a complex field of rank r, φ_{a_1 ... a_r} with each index in
// a set of size N, and φ̄ its conjugate. A polynomial in φ and φ̄ invariant under U(N)^r, the
// independent unitary changes of each index, is a sum of trace invariants: products of p
// fields φ and p fields φ̄ in which the i-th index of each φ is summed against the i-th index
// of one φ̄. Such an invariant is a bubble: a bipartite r-regular graph, black vertices for
// the φ, white for the φ̄, an edge of colour i joining the φ and φ̄ whose i-th indices are
// contracted, so that the edges of each colour form a perfect matching black → white. It is
// held as a morphe ColouredGraph with colours 0 .. r−1 (the field's colours 1 .. r): black b is
// joined to white σ_i(b) by colour i.
//
// The bubbles of the classical models:
// - rank 0 (a scalar): (φ̄φ)^k, k black and k white vertices and no edges;
// - rank 1 (a vector, U(N)): (φ̄·φ)^k, k edges of colour 1: a product of k traces;
// - rank 2 (a matrix, U(N)²): tr (MM†)^k, the cycle of length 2k alternating colours 1 and 2;
// - rank r ≥ 2: the quadratic bubble φ̄·φ (one black, one white, every colour), and the
//   quartic melonic bubble of colour c, two pairs joined by every colour but c, which crosses
//   (Gurau 2016, §3.2), the bubbles of the melonic large-N limit;
// - rank 3, sextic: the wheel K₃,₃, every colour a different matching of three blacks and
//   three whites (the "prism" or complete bipartite bubble, Bonzom, Delepouve and Rivasseau
//   2015, Nucl. Phys. B 895, 161), and the necklace, three pairs joined by colour 1 and
//   linked around a cycle by colours 2 and 3 together.
//
// Feynman graphs. With the Gaussian measure ⟨φ_a φ̄_b⟩ = δ_ab, a Wick contraction of a set of
// bubbles is a bijection σ₀ from all their black vertices to all their white ones, the
// propagators. Gluing the bubbles by σ₀ gives a bipartite (r+1)-coloured graph, colour 0 the
// propagators and colours 1 .. r the bubbles' edges: the stranded Feynman graph, whose faces
// of colours (0, i) are the closed index lines, each a free sum over N.

#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "coloured_graph.hpp"

namespace morphe::topology {

// The bubble's rank: the number of colours of its edges.
inline std::size_t bubble_rank(const ColouredGraph& B) { return B.colours; }

// Throws std::invalid_argument unless B is a bubble: every colour a permutation of its blacks.
inline void validate_bubble(const ColouredGraph& B) {
  if (B.p == 0) throw std::invalid_argument("morphe::topology::validate_bubble: a bubble has at least one black vertex");
  if (B.sigma.size() != B.colours) throw std::invalid_argument("morphe::topology::validate_bubble: one permutation per colour");
  for (const auto& s : B.sigma) {
    if (s.size() != B.p) throw std::invalid_argument("morphe::topology::validate_bubble: permutations of the p black vertices");
    std::vector<char> seen(B.p, 0);
    for (auto w : s)
      if (w >= B.p || seen[w]++) throw std::invalid_argument("morphe::topology::validate_bubble: not a permutation");
  }
}

// (φ̄φ)^k: rank 0.
inline ColouredGraph scalar_bubble(std::size_t k) { return ColouredGraph{0, k, {}}; }

// (φ̄·φ)^k: rank 1, k edges of colour 1.
inline ColouredGraph vector_bubble(std::size_t k) {
  std::vector<std::size_t> id(k);
  std::iota(id.begin(), id.end(), 0);
  return ColouredGraph{1, k, {id}};
}

// tr (MM†)^k for M of rank 2 (rows: colour 1, columns: colour 2):
// M_{a1 b1} M̄_{a2 b1} M_{a2 b2} M̄_{a3 b2} ...: colour 2 joins black i to white i, colour 1
// joins black i to white i − 1.
inline ColouredGraph trace_bubble(std::size_t k) {
  std::vector<std::size_t> id(k), shift(k);
  std::iota(id.begin(), id.end(), 0);
  for (std::size_t b = 0; b < k; b++) shift[b] = (b + k - 1) % k;
  return ColouredGraph{2, k, {shift, id}};
}

// φ̄·φ in rank r: one black and one white joined by every colour.
inline ColouredGraph quadratic_bubble(std::size_t r) { return ColouredGraph{r, 1, std::vector<std::vector<std::size_t>>(r, std::vector<std::size_t>{0})}; }

// The quartic melonic bubble of colour c (0-based among the bubble's colours) in rank r.
inline ColouredGraph quartic_melonic_bubble(std::size_t r, std::size_t c) {
  if (c >= r) throw std::invalid_argument("morphe::topology::quartic_melonic_bubble: c < r");
  ColouredGraph B{r, 2, std::vector<std::vector<std::size_t>>(r, std::vector<std::size_t>{0, 1})};
  B.sigma[c] = {1, 0};
  return B;
}

// The sextic wheel K₃,₃ of rank 3: colour i joins black b to white b + i (mod 3).
inline ColouredGraph wheel_bubble() {
  ColouredGraph B{3, 3, {}};
  for (std::size_t i = 0; i < 3; i++) B.sigma.push_back({i % 3, (1 + i) % 3, (2 + i) % 3});
  return B;
}

// The sextic necklace of rank 3: colour 1 joins black b to white b; colours 2 and 3 both join
// black b to white b + 1 (mod 3).
inline ColouredGraph necklace_bubble() { return ColouredGraph{3, 3, {{0, 1, 2}, {1, 2, 0}, {1, 2, 0}}}; }

// The automorphisms of a bubble: the pairs (π, ρ) of permutations of its blacks and whites
// with ρ σ_i = σ_i π for every colour, each as π followed by ρ in one vector of length 2p.
// For a bubble of rank 0 every pair is one, (p!)² of them. O(p! · r p) for rank >= 1.
inline std::vector<std::vector<std::size_t>> bubble_automorphisms(const ColouredGraph& B) {
  validate_bubble(B);
  const std::size_t p = B.p;
  std::vector<std::size_t> pi(p), rho(p);
  std::iota(pi.begin(), pi.end(), 0);
  std::vector<std::vector<std::size_t>> out;
  do {
    if (B.colours == 0) {
      std::iota(rho.begin(), rho.end(), 0);
      do {
        std::vector<std::size_t> a(pi);
        a.insert(a.end(), rho.begin(), rho.end());
        out.push_back(a);
      } while (std::next_permutation(rho.begin(), rho.end()));
      continue;
    }
    // ρ is forced by colour 0: ρ(σ_0(b)) = σ_0(π(b)); the other colours must agree.
    for (std::size_t b = 0; b < p; b++) rho[B.sigma[0][b]] = B.sigma[0][pi[b]];
    bool ok = true;
    for (std::size_t i = 1; i < B.colours && ok; i++)
      for (std::size_t b = 0; b < p && ok; b++) ok = rho[B.sigma[i][b]] == B.sigma[i][pi[b]];
    if (ok) {
      std::vector<std::size_t> a(pi);
      a.insert(a.end(), rho.begin(), rho.end());
      out.push_back(a);
    }
  } while (std::next_permutation(pi.begin(), pi.end()));
  return out;
}

// The Feynman graph of bubbles glued by the propagators σ₀: black b of the k-th bubble is
// global black offset_k + b, likewise for whites; sigma0 maps global blacks to global whites.
// Colour 0 is σ₀, colour i + 1 the bubbles' colour i.
inline ColouredGraph glue(const std::vector<ColouredGraph>& bubbles, const std::vector<std::size_t>& sigma0) {
  if (bubbles.empty()) throw std::invalid_argument("morphe::topology::glue: no bubbles");
  const std::size_t r = bubbles[0].colours;
  std::size_t P = 0;
  for (const auto& B : bubbles) {
    validate_bubble(B);
    if (B.colours != r) throw std::invalid_argument("morphe::topology::glue: bubbles of different ranks");
    P += B.p;
  }
  if (sigma0.size() != P) throw std::invalid_argument("morphe::topology::glue: σ₀ must match every black to a white");
  ColouredGraph G{r + 1, P, std::vector<std::vector<std::size_t>>(r + 1, std::vector<std::size_t>(P))};
  G.sigma[0] = sigma0;
  std::size_t offset = 0;
  for (const auto& B : bubbles) {
    for (std::size_t i = 0; i < r; i++)
      for (std::size_t b = 0; b < B.p; b++) G.sigma[i + 1][offset + b] = offset + B.sigma[i][b];
    offset += B.p;
  }
  G.validate();
  return G;
}

}  // namespace morphe::topology
