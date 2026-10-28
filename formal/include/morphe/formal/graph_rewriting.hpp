// Graph rewriting: double-pushout rules on labelled graphs, with relabelling, their matches,
// derivations and derivation spaces, and the action of label permutations on all of them.
//
// Graphs. A labelled graph has vertices 0..n−1 with a label each and edges (u, v) with a label
// each, loops and parallel edges allowed; a graph is either undirected or directed throughout.
//
// Rules (Ehrig, Pfender and Schneider 1973; Ehrig, Ehrig, Prange and Taentzer 2006,
// Fundamentals of Algebraic Graph Transformation, ch. 3). A rule is a span L ← K → R of injective
// graph morphisms: L is what the rule needs, R what it leaves, K what it keeps. It is held here as
// the graphs L and R and the correspondences of K, pairs (vertex of L, vertex of R) and (edge of L,
// edge of R). What lies in L outside K is deleted; what lies in R outside K is created.
//
// Relabelling (Habel and Plump 2002, Relabelling in graph transformation, ICGT 2002, LNCS 2505,
// 135–147). In the category of labelled graphs, K's morphisms preserve labels, so a span cannot
// change the label of an item it keeps: a kept vertex has one label in L, K and R. Chemistry
// needs exactly that change (a charge, a bond order, on an atom that keeps its identity). Habel
// and Plump let K be partially labelled: a kept item whose label differs between L and R is
// unlabelled in K, and the double pushout relabels it. That is the semantics here. A kept vertex
// whose label in L is "*" matches a vertex of any label and, if its label in R is "*" too, keeps
// it: a rule schema over every label, which lets a rule name a context without fixing it.
//
// Matches. A match m: L → H is an injective morphism: vertices to vertices and edges to edges,
// labels and incidences preserved ("*" matches any label), direction too in directed graphs. With
// injective matches the gluing condition of the double pushout reduces to the dangling condition
// (Ehrig et al. 2006, ch. 3): no edge of H outside m(L) may touch a vertex the rule deletes. When
// it holds, the pushout complement D = H − m(L − K) exists and is unique, and the derived graph
// H' = D + (R − K) is the second pushout; apply() builds it, with the trace of H's vertices into
// H' (the atom map, in chemistry). A rule may also carry a condition on the host and the match,
// an application condition, which every derivation must satisfy.
//
// Counting. Aut(L) acts freely on the injective matches L → H by precomposition, and its orbits
// are the images, the occurrences of L in H. So
//
//   #matches(L, H) = |Aut L| · #occurrences(L, H),
//
// which is why a rate that counts physical occurrences divides the matches by |Aut L|.
//
// Isomorphism classes. Two graphs are isomorphic when a bijection of vertices carries labels and
// the multiset of labelled edges to each other. code() encodes the class component by component
// with combinatorics/canonical.hpp, a directed edge u → v of label l being first subdivided by a
// vertex labelled "→l" joined to u by "tail" and to v by "head"; a rule's class is that of the graph
// L ⊔ R with K's correspondences drawn as edges. derivations() gives the direct derivations from a graph;
// explore() the graphs derivable from a start, up to isomorphism, with each step's multiplicity
// (the matches that lead from a class to another), so that the derivations of every length are
// counted without listing them.
//
// Symmetry. A permutation of labels (a bijection of a finite alphabet, the identity elsewhere)
// acts on graphs and on rules. If the rule set is closed under it and every condition is invariant,
// the direct derivations of g·H are the images of those of H, one for one:
// H ⇒_{q,m} H' exactly when g·H ⇒_{g·q, g·m} g·H'. equivariant() tests that equality of multisets of
// classes; closed_under() tests the hypothesis. The test fails on rule sets that are not closed,
// which is the point of stating the hypothesis.

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <morphe/combinatorics/canonical.hpp>

namespace morphe::formal::rewriting {

inline const std::string any_label = "*";
inline constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

struct Graph {
  struct Edge {
    std::size_t u = 0, v = 0;
    std::string label;
  };
  std::vector<std::string> vertex;
  std::vector<Edge> edges;
  bool directed = false;

  std::size_t add_vertex(std::string label) {
    vertex.push_back(std::move(label));
    return vertex.size() - 1;
  }
  std::size_t add_edge(std::size_t u, std::size_t v, std::string label) {
    if (u >= vertex.size() || v >= vertex.size()) throw std::invalid_argument("morphe::formal::rewriting::Graph: an edge leaves the vertices");
    edges.push_back({u, v, std::move(label)});
    return edges.size() - 1;
  }
  std::size_t size() const { return vertex.size(); }
  // The degree of each vertex, a loop counting twice.
  std::vector<std::size_t> degrees() const {
    std::vector<std::size_t> d(vertex.size(), 0);
    for (const auto& e : edges) d[e.u]++, d[e.v]++;
    return d;
  }
  // The connected components, as a label of each vertex.
  std::vector<std::size_t> components() const {
    std::vector<std::size_t> p(vertex.size());
    for (std::size_t i = 0; i < p.size(); i++) p[i] = i;
    std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return p[x] == x ? x : p[x] = root(p[x]); };
    for (const auto& e : edges) p[root(e.u)] = root(e.v);
    std::vector<std::size_t> c(vertex.size());
    std::map<std::size_t, std::size_t> number;
    for (std::size_t i = 0; i < p.size(); i++) c[i] = number.emplace(root(i), number.size()).first->second;
    return c;
  }
};

struct Match {
  std::vector<std::size_t> vertex;  // vertex of L ↦ vertex of H
  std::vector<std::size_t> edge;    // edge of L ↦ edge of H
};

struct Rule {
  std::string name;
  Graph left, right;
  std::vector<std::pair<std::size_t, std::size_t>> kept_vertices;  // K: (vertex of L, vertex of R)
  std::vector<std::pair<std::size_t, std::size_t>> kept_edges;     // K: (edge of L, edge of R)
  // An application condition on the host and the match; empty means none.
  std::function<bool(const Graph&, const Match&)> condition;

  // The vertex of R kept from vertex v of L, or none if v is deleted; and conversely.
  std::size_t right_of(std::size_t v) const {
    for (const auto& [l, r] : kept_vertices)
      if (l == v) return r;
    return none;
  }
  std::size_t left_of(std::size_t v) const {
    for (const auto& [l, r] : kept_vertices)
      if (r == v) return l;
    return none;
  }
  std::size_t right_of_edge(std::size_t e) const {
    for (const auto& [l, r] : kept_edges)
      if (l == e) return r;
    return none;
  }
  std::size_t left_of_edge(std::size_t e) const {
    for (const auto& [l, r] : kept_edges)
      if (r == e) return l;
    return none;
  }

  // Throws unless K's correspondences are injective, its edges join kept vertices consistently
  // on both sides, both sides have the same direction, and "*" labels only kept vertices.
  void validate() const {
    const auto fail = [&](const std::string& why) { throw std::invalid_argument("morphe::formal::rewriting::Rule " + name + ": " + why); };
    if (left.directed != right.directed) fail("L and R differ in direction");
    std::vector<char> l_used(left.size(), 0), r_used(right.size(), 0);
    for (const auto& [l, r] : kept_vertices) {
      if (l >= left.size() || r >= right.size()) fail("a kept vertex is not a vertex");
      if (l_used[l]++ || r_used[r]++) fail("K is not injective on vertices");
      if (right.vertex[r] == any_label && left.vertex[l] != any_label) fail("R leaves a kept vertex unlabelled that L labels");
    }
    for (std::size_t v = 0; v < left.size(); v++)
      if (left.vertex[v] == any_label && !l_used[v]) fail("a vertex labelled * is deleted");
    for (std::size_t v = 0; v < right.size(); v++)
      if (right.vertex[v] == any_label && !r_used[v]) fail("a created vertex is labelled *");
    std::vector<char> le(left.edges.size(), 0), re(right.edges.size(), 0);
    for (const auto& [l, r] : kept_edges) {
      if (l >= left.edges.size() || r >= right.edges.size()) fail("a kept edge is not an edge");
      if (le[l]++ || re[r]++) fail("K is not injective on edges");
      const auto& a = left.edges[l];
      const auto& b = right.edges[r];
      const std::size_t u = right_of(a.u), v = right_of(a.v);
      if (u == none || v == none) fail("a kept edge has a deleted end");
      const bool same = u == b.u && v == b.v;
      const bool swapped = !left.directed && u == b.v && v == b.u;
      if (!same && !swapped) fail("a kept edge's ends do not correspond");
    }
  }
};

namespace detail {

inline bool label_fits(const std::string& pattern, const std::string& label) { return pattern == any_label || pattern == label; }

// The L vertices in an order in which each, after the first of its component, is adjacent to one
// already placed, most constrained (labelled, of high degree) first.
inline std::vector<std::size_t> search_order(const Graph& L) {
  const std::size_t n = L.size();
  std::vector<std::vector<std::size_t>> adj(n);
  for (const auto& e : L.edges) adj[e.u].push_back(e.v), adj[e.v].push_back(e.u);
  const auto deg = L.degrees();
  std::vector<char> placed(n, 0);
  std::vector<std::size_t> order;
  while (order.size() < n) {
    std::size_t best = none;
    const auto better = [&](std::size_t v, std::size_t w) {
      const bool lv = L.vertex[v] != any_label, lw = L.vertex[w] != any_label;
      if (lv != lw) return lv;
      return deg[v] > deg[w];
    };
    // Prefer a vertex adjacent to one already placed.
    for (std::size_t v = 0; v < n; v++) {
      if (placed[v]) continue;
      bool frontier = false;
      for (auto w : adj[v]) frontier = frontier || placed[w];
      if (!frontier) continue;
      if (best == none || better(v, best)) best = v;
    }
    if (best == none)
      for (std::size_t v = 0; v < n; v++)
        if (!placed[v] && (best == none || better(v, best))) best = v;
    placed[best] = 1;
    order.push_back(best);
  }
  return order;
}

// The key of an edge between two images: its ordered ends (directed) or sorted ends, and label.
using EdgeKey = std::tuple<std::size_t, std::size_t, std::string>;
inline EdgeKey key(bool directed, std::size_t u, std::size_t v, const std::string& label) {
  if (!directed && v < u) std::swap(u, v);
  return {u, v, label};
}

}  // namespace detail

// Every injective match of L into H, passed to f until f returns false.
inline void for_each_match(const Graph& L, const Graph& H, const std::function<bool(const Match&)>& f) {
  if (L.directed != H.directed) throw std::invalid_argument("morphe::formal::rewriting::for_each_match: L and H differ in direction");
  const std::size_t n = L.size();
  if (n > H.size() || L.edges.size() > H.edges.size()) return;
  const auto order = detail::search_order(L);
  const auto ldeg = L.degrees(), hdeg = H.degrees();
  // H's edges grouped by key, and L's edges grouped likewise once the vertices are mapped.
  std::map<detail::EdgeKey, std::vector<std::size_t>> h_edges;
  for (std::size_t e = 0; e < H.edges.size(); e++) h_edges[detail::key(H.directed, H.edges[e].u, H.edges[e].v, H.edges[e].label)].push_back(e);
  // For each L vertex, its edges to vertices placed before it (in order), to prune early.
  std::vector<std::size_t> position(n);
  for (std::size_t i = 0; i < n; i++) position[order[i]] = i;
  std::vector<std::vector<std::size_t>> back(n);
  for (std::size_t e = 0; e < L.edges.size(); e++) {
    const auto& ed = L.edges[e];
    const std::size_t later = position[ed.u] >= position[ed.v] ? ed.u : ed.v;
    back[later].push_back(e);
  }
  Match m;
  m.vertex.assign(n, none);
  std::vector<char> used(H.size(), 0);
  bool stop = false;

  // Once the vertices are mapped: every injective assignment of L's edges to H's edges with the
  // same key, group by group.
  const auto assign_edges = [&]() {
    std::map<detail::EdgeKey, std::vector<std::size_t>> groups;
    for (std::size_t e = 0; e < L.edges.size(); e++) {
      const auto& ed = L.edges[e];
      groups[detail::key(L.directed, m.vertex[ed.u], m.vertex[ed.v], ed.label)].push_back(e);
    }
    std::vector<std::pair<std::vector<std::size_t>, const std::vector<std::size_t>*>> list;
    for (const auto& [k, es] : groups) {
      const auto it = h_edges.find(k);
      if (it == h_edges.end() || it->second.size() < es.size()) return;
      list.push_back({es, &it->second});
    }
    m.edge.assign(L.edges.size(), none);
    std::function<void(std::size_t, std::size_t, std::vector<char>&)> go = [&](std::size_t g, std::size_t i, std::vector<char>& taken) {
      if (stop) return;
      if (g == list.size()) {
        if (!f(m)) stop = true;
        return;
      }
      const auto& [es, candidates] = list[g];
      if (i == es.size()) {
        std::vector<char> next(list.size() > g + 1 ? list[g + 1].second->size() : 0, 0);
        go(g + 1, 0, next);
        return;
      }
      for (std::size_t c = 0; c < candidates->size() && !stop; c++) {
        if (taken[c]) continue;
        taken[c] = 1;
        m.edge[es[i]] = (*candidates)[c];
        go(g, i + 1, taken);
        taken[c] = 0;
      }
    };
    std::vector<char> first(list.empty() ? 0 : list[0].second->size(), 0);
    go(0, 0, first);
  };

  std::function<void(std::size_t)> place = [&](std::size_t i) {
    if (stop) return;
    if (i == n) {
      assign_edges();
      return;
    }
    const std::size_t v = order[i];
    for (std::size_t h = 0; h < H.size() && !stop; h++) {
      if (used[h] || !detail::label_fits(L.vertex[v], H.vertex[h]) || hdeg[h] < ldeg[v]) continue;
      m.vertex[v] = h;
      // Each edge to a vertex already placed must have a counterpart; the multiplicities are
      // checked when the edges are assigned.
      bool ok = true;
      for (auto e : back[v]) {
        const auto& ed = L.edges[e];
        if (!h_edges.count(detail::key(H.directed, m.vertex[ed.u], m.vertex[ed.v], ed.label))) { ok = false; break; }
      }
      if (ok) {
        used[h] = 1;
        place(i + 1);
        used[h] = 0;
      }
      m.vertex[v] = none;
    }
  };
  place(0);
}

inline std::vector<Match> matches(const Graph& L, const Graph& H) {
  std::vector<Match> out;
  for_each_match(L, H, [&](const Match& m) { out.push_back(m); return true; });
  return out;
}
inline std::size_t count_matches(const Graph& L, const Graph& H) {
  std::size_t c = 0;
  for_each_match(L, H, [&](const Match&) { c++; return true; });
  return c;
}
// |Aut G|, by counting the matches of G into itself (exact; for the small graphs of rules).
inline std::size_t automorphisms(const Graph& G) { return count_matches(G, G); }

// The dangling condition: no edge of H outside m(L) touches a vertex the rule deletes.
inline bool dangling_ok(const Rule& q, const Graph& H, const Match& m) {
  std::vector<char> in_image(H.edges.size(), 0), deleted(H.size(), 0);
  for (auto e : m.edge) in_image[e] = 1;
  for (std::size_t v = 0; v < q.left.size(); v++)
    if (q.right_of(v) == none) deleted[m.vertex[v]] = 1;
  for (std::size_t e = 0; e < H.edges.size(); e++)
    if (!in_image[e] && (deleted[H.edges[e].u] || deleted[H.edges[e].v])) return false;
  return true;
}

// The derived graph H' of the double pushout, and the trace of H's vertices into it (none for a
// deleted vertex). The match must satisfy the dangling condition.
inline Graph apply(const Rule& q, const Graph& H, const Match& m, std::vector<std::size_t>* trace = nullptr) {
  Graph out;
  out.directed = H.directed;
  std::vector<char> deleted_vertex(H.size(), 0), deleted_edge(H.edges.size(), 0);
  std::vector<std::size_t> relabel_vertex(H.size(), none), relabel_edge(H.edges.size(), none);
  for (std::size_t v = 0; v < q.left.size(); v++) {
    const std::size_t r = q.right_of(v);
    if (r == none) deleted_vertex[m.vertex[v]] = 1;
    else relabel_vertex[m.vertex[v]] = r;
  }
  for (std::size_t e = 0; e < q.left.edges.size(); e++) {
    const std::size_t r = q.right_of_edge(e);
    if (r == none) deleted_edge[m.edge[e]] = 1;
    else relabel_edge[m.edge[e]] = r;
  }
  std::vector<std::size_t> map(H.size(), none);
  for (std::size_t v = 0; v < H.size(); v++) {
    if (deleted_vertex[v]) continue;
    std::string label = H.vertex[v];
    if (relabel_vertex[v] != none && q.right.vertex[relabel_vertex[v]] != any_label) label = q.right.vertex[relabel_vertex[v]];
    map[v] = out.add_vertex(std::move(label));
  }
  for (std::size_t e = 0; e < H.edges.size(); e++) {
    if (deleted_edge[e]) continue;
    const auto& ed = H.edges[e];
    std::string label = relabel_edge[e] != none ? q.right.edges[relabel_edge[e]].label : ed.label;
    out.add_edge(map[ed.u], map[ed.v], std::move(label));
  }
  // R's vertices: kept ones are where their L vertex went; created ones are added.
  std::vector<std::size_t> r_map(q.right.size(), none);
  for (std::size_t r = 0; r < q.right.size(); r++) {
    const std::size_t l = q.left_of(r);
    r_map[r] = l != none ? map[m.vertex[l]] : out.add_vertex(q.right.vertex[r]);
  }
  for (std::size_t e = 0; e < q.right.edges.size(); e++)
    if (q.left_of_edge(e) == none) out.add_edge(r_map[q.right.edges[e].u], r_map[q.right.edges[e].v], q.right.edges[e].label);
  if (trace) *trace = map;
  return out;
}

struct Derivation {
  std::size_t rule = 0;
  Match match;
  Graph result;
  std::vector<std::size_t> trace;  // vertex of H ↦ vertex of the result, or none
};

// Every direct derivation from H by the rules: each rule, each match that satisfies the dangling
// condition and the rule's condition.
inline std::vector<Derivation> derivations(const std::vector<Rule>& rules, const Graph& H) {
  std::vector<Derivation> out;
  for (std::size_t i = 0; i < rules.size(); i++) {
    const Rule& q = rules[i];
    for_each_match(q.left, H, [&](const Match& m) {
      if (!dangling_ok(q, H, m)) return true;
      if (q.condition && !q.condition(H, m)) return true;
      Derivation d;
      d.rule = i;
      d.match = m;
      d.result = apply(q, H, m, &d.trace);
      out.push_back(std::move(d));
      return true;
    });
  }
  return out;
}

// ---- isomorphism classes -------------------------------------------------------------------------------------------------

inline combinatorics::LabelledGraph as_undirected(const Graph& G) {
  combinatorics::LabelledGraph U;
  for (const auto& l : G.vertex) U.vertex.push_back("v:" + l);
  for (const auto& e : G.edges) {
    if (!G.directed) {
      U.edges.emplace_back(e.u, e.v, e.label);
      continue;
    }
    const std::size_t x = U.vertex.size();
    U.vertex.push_back("→" + e.label);
    U.edges.emplace_back(e.u, x, "tail");
    U.edges.emplace_back(x, e.v, "head");
  }
  return U;
}

// The class of a graph. Two graphs are isomorphic exactly when their components' classes agree as
// multisets, so each component is encoded on its own and the codes sorted: the canonical form stays
// exact on graphs made of many identical parts (free molecules, isolated products), whose
// symmetries a single encoding would have to enumerate.
inline std::string code(const Graph& G) {
  const auto comp = G.components();
  std::size_t k = 0;
  for (auto c : comp) k = std::max(k, c + 1);
  std::vector<Graph> parts(k);
  std::vector<std::size_t> index(G.size());
  for (auto& p : parts) p.directed = G.directed;
  for (std::size_t v = 0; v < G.size(); v++) index[v] = parts[comp[v]].add_vertex(G.vertex[v]);
  for (const auto& e : G.edges) parts[comp[e.u]].add_edge(index[e.u], index[e.v], e.label);
  std::vector<std::string> codes;
  for (const auto& p : parts) {
    const auto c = combinatorics::canonical_form(as_undirected(p));
    if (!c.exact) throw std::runtime_error("morphe::formal::rewriting::code: a component too symmetric for an exact canonical form");
    codes.push_back(c.code);
  }
  std::sort(codes.begin(), codes.end());
  std::string out = G.directed ? "D" : "U";
  for (const auto& s : codes) out += std::to_string(s.size()) + ":" + s;
  return out;
}
inline bool isomorphic(const Graph& a, const Graph& b) { return a.directed == b.directed && code(a) == code(b); }

// The class of a rule: L ⊔ R, with K's correspondences drawn as edges (its condition is not part
// of it).
inline std::string rule_code(const Rule& q) {
  combinatorics::LabelledGraph U;
  const auto side = [&](const Graph& G, const std::string& tag) {
    const std::size_t base = U.vertex.size();
    for (const auto& l : G.vertex) U.vertex.push_back(tag + "|" + l);
    std::vector<std::size_t> edge_vertex;
    for (const auto& e : G.edges) {
      const std::size_t x = U.vertex.size();
      U.vertex.push_back(tag + "E|" + e.label);
      U.edges.emplace_back(base + e.u, x, G.directed ? "tail" : "end");
      U.edges.emplace_back(x, base + e.v, G.directed ? "head" : "end");
      edge_vertex.push_back(x);
    }
    return std::make_pair(base, edge_vertex);
  };
  const auto [lb, le] = side(q.left, "L");
  const auto [rb, re] = side(q.right, "R");
  for (const auto& [l, r] : q.kept_vertices) U.edges.emplace_back(lb + l, rb + r, "K");
  for (const auto& [l, r] : q.kept_edges) U.edges.emplace_back(le[l], re[r], "K");
  auto c = combinatorics::canonical_form(U);
  if (!c.exact) throw std::runtime_error("morphe::formal::rewriting::rule_code: a rule too symmetric for an exact canonical form");
  return (q.left.directed ? "D" : "U") + c.code;
}

// ---- the action of label permutations --------------------------------------------------------------------------------------

// A permutation of labels: listed labels go to their images, the others (and "*") are fixed.
using Relabelling = std::map<std::string, std::string>;

inline const std::string& relabel(const Relabelling& g, const std::string& l) {
  const auto it = g.find(l);
  return it == g.end() ? l : it->second;
}
// Throws unless g is a bijection of the labels it lists.
inline void validate(const Relabelling& g) {
  std::map<std::string, int> image;
  for (const auto& [a, b] : g) {
    if (!g.count(b)) throw std::invalid_argument("morphe::formal::rewriting::Relabelling: " + b + " is an image but not moved");
    if (image[b]++) throw std::invalid_argument("morphe::formal::rewriting::Relabelling: not injective");
    if (a == any_label || b == any_label) throw std::invalid_argument("morphe::formal::rewriting::Relabelling: * is fixed");
  }
}
inline Graph act(const Relabelling& g, const Graph& G) {
  Graph out = G;
  for (auto& l : out.vertex) l = relabel(g, l);
  for (auto& e : out.edges) e.label = relabel(g, e.label);
  return out;
}
inline Rule act(const Relabelling& g, const Rule& q) {
  Rule out = q;
  out.name = q.name + "^g";
  out.left = act(g, q.left);
  out.right = act(g, q.right);
  return out;
}
inline Relabelling compose(const Relabelling& g, const Relabelling& h) {  // g ∘ h
  Relabelling out;
  for (const auto& [a, b] : h) out[a] = relabel(g, b);
  for (const auto& [a, b] : g)
    if (!out.count(a)) out[a] = b;
  for (auto it = out.begin(); it != out.end();) it = it->first == it->second ? out.erase(it) : std::next(it);
  return out;
}

// Whether g·q is, up to isomorphism of rules, one of the rules, for every rule q.
inline bool closed_under(const std::vector<Rule>& rules, const Relabelling& g) {
  std::map<std::string, int> have;
  for (const auto& q : rules) have[rule_code(q)]++;
  for (const auto& q : rules)
    if (!have.count(rule_code(act(g, q)))) return false;
  return true;
}

// The classes of the direct derivations from H, with their multiplicities.
inline std::map<std::string, std::size_t> derivation_classes(const std::vector<Rule>& rules, const Graph& H) {
  std::map<std::string, std::size_t> out;
  for (const auto& d : derivations(rules, H)) out[code(d.result)]++;
  return out;
}

// Whether the direct derivations of g·H are the images under g of those of H, class by class and
// with multiplicity.
inline bool equivariant(const std::vector<Rule>& rules, const Graph& H, const Relabelling& g) {
  std::map<std::string, std::size_t> image;
  for (const auto& d : derivations(rules, H)) image[code(act(g, d.result))]++;
  return image == derivation_classes(rules, act(g, H));
}

// ---- derivation spaces ------------------------------------------------------------------------------------------------------

struct DerivationSpace {
  struct Arc {
    std::size_t from = 0, to = 0, rule = 0;
    std::size_t multiplicity = 0;  // the matches by which the rule leads from `from` to `to`
  };
  std::vector<Graph> graphs;       // a representative of each class, the start first
  std::vector<std::string> codes;
  std::vector<std::size_t> depth;  // the length of the shortest derivation reaching it
  std::vector<char> expanded;      // whether its derivations were followed (within the bounds)
  std::vector<Arc> arcs;
  // ways[c][k]: the derivations of length k from the start to a graph of class c.
  std::vector<std::vector<unsigned long long>> ways;

  std::size_t find(const std::string& c) const {
    for (std::size_t i = 0; i < codes.size(); i++)
      if (codes[i] == c) return i;
    return none;
  }
  // The classes from which no rule applies (among those expanded).
  std::vector<std::size_t> normal_forms() const {
    std::vector<char> out(graphs.size(), 0);
    for (const auto& a : arcs) out[a.from] = 1;
    std::vector<std::size_t> r;
    for (std::size_t c = 0; c < graphs.size(); c++)
      if (expanded[c] && !out[c]) r.push_back(c);
    return r;
  }
  unsigned long long derivations_to(std::size_t c) const {
    unsigned long long s = 0;
    for (auto w : ways[c]) s += w;
    return s;
  }
};

// The graphs derivable from S in at most max_depth steps, up to isomorphism, breadth first, with
// at most `limit` classes; and the number of derivations of each length reaching each class.
inline DerivationSpace explore(const std::vector<Rule>& rules, const Graph& S, std::size_t max_depth, std::size_t limit = 10000) {
  DerivationSpace D;
  std::unordered_map<std::string, std::size_t> index;
  const auto add = [&](const Graph& G, std::string c, std::size_t depth) {
    const auto it = index.find(c);
    if (it != index.end()) return it->second;
    if (D.graphs.size() >= limit) throw std::runtime_error("morphe::formal::rewriting::explore: more classes than the limit");
    index.emplace(c, D.graphs.size());
    D.graphs.push_back(G);
    D.codes.push_back(std::move(c));
    D.depth.push_back(depth);
    D.expanded.push_back(0);
    return D.graphs.size() - 1;
  };
  add(S, code(S), 0);
  for (std::size_t i = 0; i < D.graphs.size(); i++) {
    if (D.depth[i] >= max_depth) continue;
    D.expanded[i] = 1;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> arcs;  // (to, rule) ↦ multiplicity
    const Graph G = D.graphs[i];
    for (const auto& d : derivations(rules, G)) {
      const std::size_t to = add(d.result, code(d.result), D.depth[i] + 1);
      arcs[{to, d.rule}]++;
    }
    for (const auto& [k, mult] : arcs) D.arcs.push_back({i, k.first, k.second, mult});
  }
  // Derivations of each length, by dynamic programming along the arcs.
  D.ways.assign(D.graphs.size(), std::vector<unsigned long long>(max_depth + 1, 0));
  D.ways[0][0] = 1;
  for (std::size_t k = 0; k < max_depth; k++)
    for (const auto& a : D.arcs) D.ways[a.to][k + 1] += D.ways[a.from][k] * a.multiplicity;
  return D;
}

// The orbits of a set of classes under the group generated by the given label permutations: a
// number for each class (classes whose images fall outside the set are left in orbits of their
// own reach within it).
inline std::vector<std::size_t> orbits(const std::vector<Graph>& graphs, const std::vector<Relabelling>& generators) {
  std::unordered_map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < graphs.size(); i++) index.emplace(code(graphs[i]), i);
  std::vector<std::size_t> p(graphs.size());
  for (std::size_t i = 0; i < p.size(); i++) p[i] = i;
  std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return p[x] == x ? x : p[x] = root(p[x]); };
  for (std::size_t i = 0; i < graphs.size(); i++)
    for (const auto& g : generators) {
      const auto it = index.find(code(act(g, graphs[i])));
      if (it != index.end()) p[root(i)] = root(it->second);
    }
  std::vector<std::size_t> out(graphs.size());
  std::map<std::size_t, std::size_t> number;
  for (std::size_t i = 0; i < p.size(); i++) out[i] = number.emplace(root(i), number.size()).first->second;
  return out;
}

}  // namespace morphe::formal::rewriting
