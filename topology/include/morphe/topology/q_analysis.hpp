// Atkin's Q-analysis of a family of simplices.
//
// A relation between two finite sets, read row by row, is a family of simplices on one of them
// (Dowker 1952; Atkin 1974, Mathematical Structure in Human Affairs, ch. 3). Two simplices of the
// family are q-near when they share a q-face, that is at least q + 1 vertices, and q-connected when
// a chain of q-near simplices joins them. For each q the q-connected components of the simplices of
// dimension at least q partition them; their number Q_q, for q = 0 up to the largest dimension, is
// the structure vector Q of the family (Atkin 1974, §3.3). Q_0 counts the connected components of
// the complex the family generates; a large Q_q at high q says the family hangs together only
// through few shared vertices.
//
// The family is taken as given: repeated or nested simplices are counted as they occur, which is
// what the analysis of a relation needs (each row its own simplex). O(m² d log d) for m simplices
// of at most d vertices.

#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <set>
#include <vector>

#include <morphe/topology/simplicial.hpp>

namespace morphe::topology {

inline std::vector<std::size_t> q_analysis(std::vector<Simplex> family) {
  for (auto& s : family) std::sort(s.begin(), s.end());
  family.erase(std::remove_if(family.begin(), family.end(), [](const Simplex& s) { return s.empty(); }), family.end());
  std::size_t top = 0;
  for (const auto& s : family) top = std::max(top, s.size() - 1);
  std::vector<std::size_t> Q(family.empty() ? 0 : top + 1, 0);
  for (std::size_t q = 0; q < Q.size(); q++) {
    std::vector<std::size_t> live;
    for (std::size_t i = 0; i < family.size(); i++)
      if (family[i].size() >= q + 1) live.push_back(i);
    std::vector<std::size_t> parent(live.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](std::size_t x) {
      while (parent[x] != x) x = parent[x] = parent[parent[x]];
      return x;
    };
    for (std::size_t a = 0; a < live.size(); a++)
      for (std::size_t b = a + 1; b < live.size(); b++) {
        std::vector<std::size_t> shared;
        std::set_intersection(family[live[a]].begin(), family[live[a]].end(), family[live[b]].begin(), family[live[b]].end(), std::back_inserter(shared));
        if (shared.size() >= q + 1) parent[find(a)] = find(b);
      }
    std::set<std::size_t> roots;
    for (std::size_t a = 0; a < live.size(); a++) roots.insert(find(a));
    Q[q] = roots.size();
  }
  return Q;
}

}  // namespace morphe::topology
