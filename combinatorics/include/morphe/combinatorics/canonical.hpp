// The canonical form of a finite multigraph with labelled vertices and labelled edges, and the
// order of its automorphism group.
//
// A labelled multigraph has vertices 0..n−1 with a label each and undirected edges {u, v} (loops
// and parallel edges allowed) with a label each. Two are isomorphic when a bijection of the
// vertices carries labels to labels and the multiset of labelled edges to itself. A canonical
// form is an encoding equal for two graphs exactly when they are isomorphic.
//
// The method is individualisation by refinement (McKay and Piperno 2014, J. Symbolic Comput.
// 60, 94, §2, in its simplest form): colour refinement (Weisfeiler and Leman 1968) partitions
// the vertices by label and, repeatedly, by the multiset of (edge label, neighbour's colour) until
// the partition is stable; every isomorphism respects the stable partition, so only the orderings
// within its cells are tried, the encoding (labels in order, then the labelled edges between the
// ordered vertices) is taken least over them, and the orderings that attain it are counted: the
// automorphisms. Exact; exponential only in the sizes of the cells refinement cannot split. Where
// they would allow more than `limit` orderings, the stable partition's encoding is returned
// instead, marked inexact: a certificate that isomorphic graphs share, which non-isomorphic ones
// rarely do.

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

struct LabelledGraph {
  std::vector<std::string> vertex;                                     // a label for each vertex
  std::vector<std::tuple<std::size_t, std::size_t, std::string>> edges;  // {u, v} with its label
};

struct CanonicalForm {
  std::string code;              // equal exactly for isomorphic graphs (when exact)
  std::size_t automorphisms = 0; // |Aut G| (when exact)
  bool exact = true;
  std::vector<std::size_t> order;  // the vertices in canonical order
};

inline CanonicalForm canonical_form(const LabelledGraph& G, std::size_t limit = 40320) {
  const std::size_t n = G.vertex.size();
  std::vector<std::vector<std::vector<std::string>>> between(n, std::vector<std::vector<std::string>>(n));
  std::vector<std::vector<std::pair<std::string, std::size_t>>> incident(n);  // (label, neighbour)
  for (const auto& [u, v, label] : G.edges) {
    between[u][v].push_back(label);
    if (u != v) between[v][u].push_back(label);
    incident[u].push_back({label, v});
    if (u != v) incident[v].push_back({label, u});
  }
  for (auto& row : between)
    for (auto& cell : row) std::sort(cell.begin(), cell.end());
  // Colour refinement to a stable partition.
  std::vector<std::size_t> colour(n);
  {
    std::map<std::string, std::size_t> ids;
    for (std::size_t v = 0; v < n; v++) ids[G.vertex[v]];
    std::size_t k = 0;
    for (auto& [label, id] : ids) id = k++;
    for (std::size_t v = 0; v < n; v++) colour[v] = ids.at(G.vertex[v]);
    for (std::size_t classes = ids.size();;) {
      std::map<std::pair<std::size_t, std::vector<std::pair<std::string, std::size_t>>>, std::size_t> next;
      std::vector<std::pair<std::size_t, std::vector<std::pair<std::string, std::size_t>>>> signature(n);
      for (std::size_t v = 0; v < n; v++) {
        std::vector<std::pair<std::string, std::size_t>> around;
        for (const auto& [label, w] : incident[v]) around.push_back({label, colour[w]});
        std::sort(around.begin(), around.end());
        signature[v] = {colour[v], std::move(around)};
        next[signature[v]];
      }
      std::size_t c = 0;
      for (auto& [sig, id] : next) id = c++;
      for (std::size_t v = 0; v < n; v++) colour[v] = next.at(signature[v]);
      if (next.size() == classes) break;
      classes = next.size();
    }
  }
  // The cells of the stable partition, in the order of their colours.
  std::vector<std::size_t> order(n);
  for (std::size_t v = 0; v < n; v++) order[v] = v;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return colour[a] != colour[b] ? colour[a] < colour[b] : a < b; });
  std::vector<std::pair<std::size_t, std::size_t>> cells;
  double orderings = 1;
  for (std::size_t i = 0; i < n;) {
    std::size_t j = i;
    while (j < n && colour[order[j]] == colour[order[i]]) j++;
    cells.push_back({i, j});
    for (std::size_t m = 2; m <= j - i; m++) orderings *= static_cast<double>(m);
    i = j;
  }
  auto encode = [&](const std::vector<std::size_t>& p) {
    std::string s = std::to_string(n) + "|";
    for (std::size_t i = 0; i < n; i++) s += G.vertex[p[i]] + "\x1f";
    s += "|";
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i; j < n; j++) {
        const auto& labels = between[p[i]][p[j]];
        if (labels.empty()) continue;
        s += std::to_string(i) + "-" + std::to_string(j) + ":";
        for (const auto& x : labels) s += x + "\x1e";
        s += ";";
      }
    return s;
  };
  CanonicalForm out;
  if (orderings > static_cast<double>(limit)) {
    out.exact = false;
    out.code = "~" + encode(order);
    out.order = order;
    return out;
  }
  std::vector<std::size_t> p = order;
  std::function<void(std::size_t)> walk = [&](std::size_t c) {
    if (c == cells.size()) {
      const std::string s = encode(p);
      if (out.code.empty() || s < out.code) {
        out.code = s;
        out.order = p;
        out.automorphisms = 1;
      } else if (s == out.code) {
        out.automorphisms++;
      }
      return;
    }
    auto first = p.begin() + static_cast<std::ptrdiff_t>(cells[c].first), last = p.begin() + static_cast<std::ptrdiff_t>(cells[c].second);
    std::sort(first, last);
    do walk(c + 1);
    while (std::next_permutation(first, last));
  };
  walk(0);
  return out;
}

}  // namespace morphe::combinatorics
