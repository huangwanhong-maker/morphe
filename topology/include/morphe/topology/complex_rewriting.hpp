// Rewriting regular complexes by double pushouts on their incidence graphs: the local changes of a cellular
// configuration as rules, instantiated at a site, applied by morphe's graph rewriting
// (formal/graph_rewriting.hpp), and admitted only when what they derive is again a regular complex.
//
// The category. Regular CW complexes do not form an adhesive category: a pushout of complexes along a
// subcomplex need not be regular (a segment glued at both ends to one vertex is a loop). The theory of the
// double pushout — unique pushout complements, the local Church–Rosser, parallelism and concurrency
// theorems — holds in adhesive categories (Lack and Sobociński 2005, RAIRO Theor. Inf. Appl. 39, 511–545),
// among them every presheaf category, so graphs, and labelled graphs as a slice of graphs. So the rule acts
// on the incidence graph of the complex (the Hasse diagram of its face poset: a vertex for each cell,
// labelled "k|type"; an arrow from each cell to each facet, labelled with the incidence number "+" or "-"),
// where the double pushout is well defined, and the derived graph is decoded as a face poset and admitted
// only if it is that of a regular complex with ∂∂ = 0 (regular_complex.hpp). A rewrite of complexes is
// thus a partial operation: defined when the match satisfies the gluing condition and the derived poset is
// a CW poset, refused with the reason otherwise.
//
// The gluing condition, matches being injective, is the dangling condition (Ehrig et al. 2006, ch. 3): a
// cell can be deleted only with its incidences to every cell of which it is a facet, so a rule that
// removes a cell must name all its cofaces. A street segment that bounds a block cannot be removed unless
// the rule also changes the block.
//
// Rules are instantiated at a site: each builder reads the cells about the site (its boundary, its
// cofaces) and returns a rule of the size the site needs, with its match. Each family — subdivide an edge
// with c cofaces, partition a cell of n facets — is countable, a rule schema whose instances are ordinary
// rules. The builders and their inverses:
//
//   attach_cell(type, facets)  a cell along existing cells (a 1-cell between two vertices: a new street);
//   remove_cell(σ)             its inverse, for a cell that bounds nothing;
//   subdivide_edge(e)          e ↦ e₁ m e₂, its cofaces passing to e₁ + e₂, the chain map e ↦ e₁ + e₂;
//   merge_at_vertex(m)         its inverse, at a vertex on exactly two edges;
//   partition_cell(σ, F₁)      a k-cell cut in two by a new (k−1)-cell w along the boundary of the facets F₁:
//                              ∂w = ∂(Σ_{τ∈F₁}[σ:τ]τ), σ₁ = F₁ − w, σ₂ = F₂ + w; split_face(f, u, v) the
//                              case k = 2, cutting a 2-cell by a new 1-cell from u to v;
//   merge_cells(w)             its inverse: the two k-cells across a (k−1)-cell w made one;
//   relabel(σ, type)           a cell's type changed, its incidences kept (Habel and Plump 2002, ICGT,
//                              LNCS 2505, 135–147).
//
// Every new incidence number is chosen so that ∂∂ = 0 holds again; the chain of a cell's cofaces is kept
// (a coface of a partitioned σ gets σ₁ + σ₂ with σ's sign). Identities: cells kept keep their ids, cells
// created get fresh ones, so a cell can be followed through a history; each rewrite records what it read
// and wrote (formal/traces.hpp), from which a record's causal order follows. Checks:
// tests/check_complex_rewriting.cpp.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/formal/graph_rewriting.hpp>
#include <morphe/formal/traces.hpp>

#include "regular_complex.hpp"

namespace morphe::topology {

namespace rw = formal::rewriting;

inline std::string cell_label(std::size_t dim, const std::string& type) { return std::to_string(dim) + "|" + type; }
inline std::string sign_label(int s) { return s > 0 ? "+" : "-"; }

// The incidence graph: vertex i is cell i, labelled "k|type"; an arrow from each cell to each facet.
inline rw::Graph incidence_graph(const RegularComplex& X) {
  rw::Graph G;
  G.directed = true;
  for (const auto& c : X.cells) G.add_vertex(cell_label(c.dim, c.type));
  for (std::size_t i = 0; i < X.size(); i++) {
    const auto& c = X.cells[i];
    for (std::size_t j = 0; j < c.facets.size(); j++) {
      if (c.signs[j] == 0) throw std::invalid_argument("morphe::topology::incidence_graph: the complex is not oriented");
      G.add_edge(i, c.facets[j], sign_label(c.signs[j]));
    }
  }
  return G;
}

// The face poset an incidence graph encodes, cells sorted by dimension (stably), with the given ids; `order`
// receives the cell of each vertex. Throws std::invalid_argument on a label or arrow that is not a cell's.
inline RegularComplex from_incidence_graph(const rw::Graph& G, const std::vector<std::uint64_t>& ids, std::uint64_t next_id, std::vector<std::size_t>* order = nullptr) {
  if (!G.directed || ids.size() != G.size()) throw std::invalid_argument("morphe::topology::from_incidence_graph: a directed graph and an id for each vertex");
  std::vector<std::size_t> dim(G.size());
  std::vector<std::string> type(G.size());
  for (std::size_t v = 0; v < G.size(); v++) {
    const auto& l = G.vertex[v];
    const auto bar = l.find('|');
    if (bar == std::string::npos || bar == 0) throw std::invalid_argument("morphe::topology::from_incidence_graph: the label \"" + l + "\" is not \"k|type\"");
    dim[v] = std::stoul(l.substr(0, bar));
    type[v] = l.substr(bar + 1);
  }
  std::vector<std::size_t> perm(G.size());
  for (std::size_t v = 0; v < perm.size(); v++) perm[v] = v;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) { return dim[a] < dim[b]; });
  std::vector<std::size_t> at(G.size());
  for (std::size_t k = 0; k < perm.size(); k++) at[perm[k]] = k;
  RegularComplex X;
  X.next_id = next_id;
  for (std::size_t k = 0; k < perm.size(); k++) X.cells.push_back({dim[perm[k]], type[perm[k]], {}, {}, ids[perm[k]]});
  for (const auto& e : G.edges) {
    if (e.label != "+" && e.label != "-") throw std::invalid_argument("morphe::topology::from_incidence_graph: an arrow labelled \"" + e.label + "\", not + or -");
    auto& c = X.cells[at[e.u]];
    c.facets.push_back(at[e.v]);
    c.signs.push_back(e.label == "+" ? 1 : -1);
  }
  for (const auto& c : X.cells) X.next_id = std::max(X.next_id, c.id + 1);
  if (order) *order = at;
  return X;
}

// Whether two complexes are isomorphic as typed face posets — a bijection of cells keeping dimension, type and
// the facet relation, which for regular complexes is a cellular isomorphism (Björner 1984) — and, when
// `oriented`, keeping the incidence numbers too. Reversing a cell's orientation changes its signs and its
// cofaces' signs on it, not the complex. Decided by searching for one injective match of A's incidence graph
// into B's (formal/graph_rewriting.hpp): with as many cells and incidences on each side, an injective map
// preserving labels and incidences is a bijection, an isomorphism. Exact; the search is pruned by labels,
// degrees and the incidences already placed, and is quick on the complexes of models (hundreds of cells
// with few symmetries), exponential at worst.
inline bool isomorphic(const RegularComplex& a, const RegularComplex& b, bool oriented = false) {
  if (a.counts() != b.counts()) return false;
  auto A = incidence_graph(a), B = incidence_graph(b);
  if (A.edges.size() != B.edges.size()) return false;
  if (!oriented) {
    for (auto& e : A.edges) e.label.clear();
    for (auto& e : B.edges) e.label.clear();
  }
  bool found = false;
  rw::for_each_match(A, B, [&](const rw::Match&) { found = true; return false; });
  return found;
}

// A rule instantiated at a site of a host.
struct ComplexRule {
  std::string kind;      // attach_cell, remove_cell, subdivide_edge, merge_at_vertex, partition_cell, merge_cells, relabel
  rw::Rule rule;
  rw::Match match;       // into the host's incidence graph
  std::string refused;   // nonempty when the site does not admit the rule, and why
};

struct Rewritten {
  bool applied = false;
  std::string reason;                    // why not, when refused
  RegularComplex result;
  std::vector<std::size_t> trace;        // host cell ↦ result cell, or no_cell when deleted
  std::vector<std::uint64_t> created, deleted, relabelled;  // ids
  formal::traces::Footprint footprint;   // items "c<id>" (a cell) and "i<id>><id>" (an incidence)
};

namespace detail {

inline std::string cell_item(std::uint64_t id) { return "c" + std::to_string(id); }
inline std::string incidence_item(std::uint64_t a, std::uint64_t b) { return "i" + std::to_string(a) + ">" + std::to_string(b); }

// Builds a rule's two sides and its match from cells of the host.
class Builder {
 public:
  Builder(const RegularComplex& X, std::string name) : X_(X) {
    q_.rule.name = std::move(name);
    q_.rule.left.directed = q_.rule.right.directed = true;
    for (std::size_t i = 0; i < X.size(); i++)
      for (std::size_t j = 0; j < X.cells[i].facets.size(); j++) edge_[{i, X.cells[i].facets[j]}] = edge_count_++;
  }
  // The L vertex of a host cell (added with its label, and matched to it).
  std::size_t left(std::size_t cell) {
    const auto it = left_.find(cell);
    if (it != left_.end()) return it->second;
    const auto& c = X_.cells.at(cell);
    const std::size_t v = q_.rule.left.add_vertex(cell_label(c.dim, c.type));
    q_.match.vertex.push_back(cell);
    left_[cell] = v;
    return v;
  }
  // The R vertex kept from a host cell, with its label (or a new type).
  std::size_t keep(std::size_t cell, const std::string* new_type = nullptr) {
    const auto it = right_.find(cell);
    if (it != right_.end()) return it->second;
    const auto& c = X_.cells.at(cell);
    const std::size_t l = left(cell);
    const std::size_t r = q_.rule.right.add_vertex(cell_label(c.dim, new_type ? *new_type : c.type));
    q_.rule.kept_vertices.push_back({l, r});
    right_[cell] = r;
    return r;
  }
  std::size_t create(std::size_t dim, const std::string& type) { return q_.rule.right.add_vertex(cell_label(dim, type)); }
  // The L arrow of a host incidence σ → τ (deleted unless kept).
  std::size_t left_arrow(std::size_t s, std::size_t t) {
    const auto it = edge_.find({s, t});
    if (it == edge_.end()) throw std::invalid_argument("morphe::topology rule builder: no incidence between the cells");
    const std::size_t e = q_.rule.left.add_edge(left(s), left(t), sign_label(X_.incidence(s, t)));
    q_.match.edge.push_back(it->second);
    return e;
  }
  void keep_arrow(std::size_t s, std::size_t t) {
    const std::size_t e = left_arrow(s, t);
    const std::size_t r = q_.rule.right.add_edge(keep(s), keep(t), sign_label(X_.incidence(s, t)));
    q_.rule.kept_edges.push_back({e, r});
  }
  void arrow(std::size_t r_from, std::size_t r_to, int sign) { q_.rule.right.add_edge(r_from, r_to, sign_label(sign)); }
  ComplexRule done(std::string kind) {
    q_.kind = std::move(kind);
    return std::move(q_);
  }

 private:
  const RegularComplex& X_;
  ComplexRule q_;
  std::map<std::size_t, std::size_t> left_, right_;
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_;
  std::size_t edge_count_ = 0;
};

inline ComplexRule refusal(std::string kind, std::string why) {
  ComplexRule q;
  q.kind = std::move(kind);
  q.refused = std::move(why);
  return q;
}

}  // namespace detail

// ---- the builders ----------------------------------------------------------------------------------------------

// A new cell of dimension dim along existing (dim − 1)-cells; its signs those that make ∂∂ = 0 (or given).
inline ComplexRule attach_cell(const RegularComplex& X, std::size_t dim, const std::string& type, const std::vector<std::size_t>& facets, std::vector<int> signs = {}) {
  if (dim == 0 && !facets.empty()) return detail::refusal("attach_cell", "a vertex has no facets");
  RegularComplex probe = X;
  try {
    probe.add(dim, type, facets, signs);
    orient(probe);
  } catch (const std::exception& e) {
    return detail::refusal("attach_cell", e.what());
  }
  const auto& c = probe.cells.back();
  detail::Builder b(X, "attach " + std::to_string(dim) + "-cell " + type);
  const std::size_t r = b.create(dim, type);
  for (std::size_t j = 0; j < facets.size(); j++) b.arrow(r, b.keep(facets[j]), c.signs[j]);
  return b.done("attach_cell");
}

// A 1-cell from tail to head.
inline ComplexRule attach_edge(const RegularComplex& X, std::size_t tail, std::size_t head, const std::string& type) { return attach_cell(X, 1, type, {tail, head}, {-1, 1}); }

// The cell σ removed with its incidences to its facets; refused by the dangling condition if σ bounds a cell.
inline ComplexRule remove_cell(const RegularComplex& X, std::size_t s) {
  detail::Builder b(X, "remove " + cell_label(X.cells.at(s).dim, X.cells[s].type));
  b.left(s);
  for (std::size_t f : X.cells[s].facets) {
    b.keep(f);
    b.left_arrow(s, f);
  }
  return b.done("remove_cell");
}

// e ↦ e₁ m e₂: a new vertex m of type vertex_type inside the edge; each coface of e has e₁ and e₂ with e's sign.
inline ComplexRule subdivide_edge(const RegularComplex& X, std::size_t e, const std::string& vertex_type) {
  if (X.cells.at(e).dim != 1) return detail::refusal("subdivide_edge", "not a 1-cell");
  const auto [t, h] = X.ends(e);
  const auto up = X.cofaces()[e];
  detail::Builder b(X, "subdivide " + cell_label(1, X.cells[e].type));
  b.left(e);
  const std::size_t rt = b.keep(t), rh = b.keep(h);
  b.left_arrow(e, t), b.left_arrow(e, h);
  const std::size_t m = b.create(0, vertex_type), e1 = b.create(1, X.cells[e].type), e2 = b.create(1, X.cells[e].type);
  b.arrow(e1, rt, -1), b.arrow(e1, m, 1), b.arrow(e2, m, -1), b.arrow(e2, rh, 1);
  for (std::size_t f : up) {
    const int s = X.incidence(f, e);
    const std::size_t rf = b.keep(f);
    b.left_arrow(f, e);
    b.arrow(rf, e1, s), b.arrow(rf, e2, s);
  }
  return b.done("subdivide_edge");
}

// The inverse: a vertex m on exactly two edges e₁ (from a) and e₂ (to b), a ≠ b, replaced by one edge a → b.
inline ComplexRule merge_at_vertex(const RegularComplex& X, std::size_t m) {
  if (X.cells.at(m).dim != 0) return detail::refusal("merge_at_vertex", "not a vertex");
  const auto up = X.cofaces();
  if (up[m].size() != 2) return detail::refusal("merge_at_vertex", "the vertex is on " + std::to_string(up[m].size()) + " edges, not two");
  const std::size_t e1 = up[m][0], e2 = up[m][1];
  const auto other = [&](std::size_t e) { const auto x = X.ends(e); return x[0] == m ? x[1] : x[0]; };
  const std::size_t a = other(e1), bb = other(e2);
  if (a == bb) return detail::refusal("merge_at_vertex", "both edges end at one vertex: the merged edge would be a loop");
  // σ₁ = +1 when e₁ runs a → m, σ₂ = +1 when e₂ runs m → b: along a → m → b the chain is σ₁e₁ + σ₂e₂.
  const int s1 = X.ends(e1)[0] == a ? 1 : -1, s2 = X.ends(e2)[0] == m ? 1 : -1;
  detail::Builder b(X, "merge edges at " + cell_label(0, X.cells[m].type));
  b.left(m), b.left(e1), b.left(e2);
  const std::size_t ra = b.keep(a), rb = b.keep(bb);
  b.left_arrow(e1, a), b.left_arrow(e1, m), b.left_arrow(e2, m), b.left_arrow(e2, bb);
  const std::size_t e = b.create(1, X.cells[e1].type);
  b.arrow(e, ra, -1), b.arrow(e, rb, 1);
  std::set<std::size_t> faces(up[e1].begin(), up[e1].end());
  faces.insert(up[e2].begin(), up[e2].end());
  for (std::size_t f : faces) {
    const std::size_t rf = b.keep(f);
    int s = 0;
    if (X.incidence(f, e1)) b.left_arrow(f, e1), s = X.incidence(f, e1) * s1;
    if (X.incidence(f, e2)) {
      b.left_arrow(f, e2);
      const int s_two = X.incidence(f, e2) * s2;
      if (s != 0 && s != s_two) return detail::refusal("merge_at_vertex", "a coface crosses the two edges in opposite senses");
      s = s_two;
    }
    b.arrow(rf, e, s);
  }
  return b.done("merge_at_vertex");
}

// σ (dimension k ≥ 2) cut in two along the boundary of the facets in `part` (cells of the host): a new
// (k−1)-cell w of type wall_type with ∂w = ∂(Σ_{τ∈part}[σ:τ]τ); σ₁ = part − w and σ₂ = rest + w, both of σ's
// type; each coface of σ gets σ₁ + σ₂ with σ's sign.
inline ComplexRule partition_cell(const RegularComplex& X, std::size_t s, const std::vector<std::size_t>& part, const std::string& wall_type) {
  const auto& c = X.cells.at(s);
  if (c.dim < 2) return detail::refusal("partition_cell", "only cells of dimension 2 or more are partitioned");
  std::set<std::size_t> one(part.begin(), part.end());
  for (std::size_t f : one)
    if (!X.incidence(s, f)) return detail::refusal("partition_cell", "the part names a cell that is not a facet");
  if (one.empty() || one.size() == c.facets.size()) return detail::refusal("partition_cell", "each side needs at least one facet");
  // ζ = ∂(Σ_{τ ∈ part} [σ:τ] τ), on (k−2)-cells.
  std::map<std::size_t, int> zeta;
  for (std::size_t j = 0; j < c.facets.size(); j++) {
    if (!one.count(c.facets[j])) continue;
    const auto& t = X.cells[c.facets[j]];
    for (std::size_t q = 0; q < t.facets.size(); q++) zeta[t.facets[q]] += c.signs[j] * t.signs[q];
  }
  std::vector<std::pair<std::size_t, int>> wall;
  for (const auto& [r, v] : zeta) {
    if (v == 0) continue;
    if (v != 1 && v != -1) return detail::refusal("partition_cell", "the seam crosses a cell more than once");
    wall.push_back({r, v});
  }
  if (wall.empty()) return detail::refusal("partition_cell", "the part has no boundary: it is closed");
  const auto up = X.cofaces()[s];
  detail::Builder b(X, "partition " + cell_label(c.dim, c.type));
  b.left(s);
  const std::size_t w = b.create(c.dim - 1, wall_type), s1 = b.create(c.dim, c.type), s2 = b.create(c.dim, c.type);
  for (const auto& [r, v] : wall) b.arrow(w, b.keep(r), v);
  for (std::size_t j = 0; j < c.facets.size(); j++) {
    const std::size_t rf = b.keep(c.facets[j]);
    b.left_arrow(s, c.facets[j]);
    b.arrow(one.count(c.facets[j]) ? s1 : s2, rf, c.signs[j]);
  }
  b.arrow(s1, w, -1), b.arrow(s2, w, 1);
  for (std::size_t g : up) {
    const int sg = X.incidence(g, s);
    const std::size_t rg = b.keep(g);
    b.left_arrow(g, s);
    b.arrow(rg, s1, sg), b.arrow(rg, s2, sg);
  }
  return b.done("partition_cell");
}

// A 2-cell cut by a new 1-cell from u to v, vertices of its boundary: the part is the boundary's edges from u
// to v in the direction of its orientation.
inline ComplexRule split_face(const RegularComplex& X, std::size_t f, std::size_t u, std::size_t v, const std::string& edge_type) {
  const auto& c = X.cells.at(f);
  if (c.dim != 2) return detail::refusal("split_face", "not a 2-cell");
  if (u == v) return detail::refusal("split_face", "the cut needs two distinct vertices");
  // Walk the boundary from u in its direction, collecting edges until v.
  const std::size_t m = c.facets.size();
  const auto start_of = [&](std::size_t j) { const auto e = X.ends(c.facets[j]); return c.signs[j] > 0 ? e[0] : e[1]; };
  const auto end_of = [&](std::size_t j) { const auto e = X.ends(c.facets[j]); return c.signs[j] > 0 ? e[1] : e[0]; };
  std::vector<std::size_t> part;
  std::size_t at = u;
  for (std::size_t step = 0; step < m && at != v; step++) {
    std::size_t j = m;
    for (std::size_t q = 0; q < m; q++)
      if (start_of(q) == at) { j = q; break; }
    if (j == m) return detail::refusal("split_face", "the first vertex is not on the 2-cell's boundary");
    part.push_back(c.facets[j]);
    at = end_of(j);
  }
  if (at != v) return detail::refusal("split_face", "the second vertex is not on the 2-cell's boundary");
  return partition_cell(X, f, part, edge_type);
}

// The two k-cells across a (k−1)-cell w (w a facet of exactly those two) made one, of the first's type: its
// facets the facets of both but w, the second's signs turned if need be so that w cancels. Refused when the
// two share another facet (their union would not be a ball) or bound different cells.
inline ComplexRule merge_cells(const RegularComplex& X, std::size_t w) {
  const auto up = X.cofaces();
  if (up[w].size() != 2) return detail::refusal("merge_cells", "the cell is a facet of " + std::to_string(up[w].size()) + " cells, not two");
  const std::size_t c1 = up[w][0], c2 = up[w][1];
  const int flip = -X.incidence(c1, w) * X.incidence(c2, w);
  std::set<std::size_t> f1(X.cells[c1].facets.begin(), X.cells[c1].facets.end()), f2(X.cells[c2].facets.begin(), X.cells[c2].facets.end());
  for (std::size_t f : f1)
    if (f != w && f2.count(f)) return detail::refusal("merge_cells", "the two cells share another facet: their union is not a cell");
  std::set<std::size_t> g1(up[c1].begin(), up[c1].end()), g2(up[c2].begin(), up[c2].end());
  if (g1 != g2) return detail::refusal("merge_cells", "the two cells bound different cells");
  detail::Builder b(X, "merge " + cell_label(X.cells[c1].dim, X.cells[c1].type) + " across " + cell_label(X.cells[w].dim, X.cells[w].type));
  b.left(w), b.left(c1), b.left(c2);
  for (std::size_t r : X.cells[w].facets) b.keep(r), b.left_arrow(w, r);
  const std::size_t merged = b.create(X.cells[c1].dim, X.cells[c1].type);
  for (std::size_t f : X.cells[c1].facets) {
    b.left_arrow(c1, f);
    if (f != w) b.arrow(merged, b.keep(f), X.incidence(c1, f));
  }
  for (std::size_t f : X.cells[c2].facets) {
    b.left_arrow(c2, f);
    if (f != w) b.arrow(merged, b.keep(f), flip * X.incidence(c2, f));
  }
  for (std::size_t g : g1) {
    const std::size_t rg = b.keep(g);
    b.left_arrow(g, c1), b.left_arrow(g, c2);
    if (X.incidence(g, c1) != flip * X.incidence(g, c2)) return detail::refusal("merge_cells", "a coface holds the two cells in opposite senses");
    b.arrow(rg, merged, X.incidence(g, c1));
  }
  return b.done("merge_cells");
}

// A cell's type changed; the cell and its incidences kept.
inline ComplexRule relabel(const RegularComplex& X, std::size_t s, const std::string& type) {
  if (type.find('|') != std::string::npos) return detail::refusal("relabel", "a type may not contain '|'");
  detail::Builder b(X, "relabel " + cell_label(X.cells.at(s).dim, X.cells[s].type) + " as " + type);
  b.keep(s, &type);
  return b.done("relabel");
}

// ---- application -----------------------------------------------------------------------------------------------

// Applies a rule instantiated at its site: the match is checked, then the dangling condition, then the double
// pushout on the incidence graph, then the derived face poset's regularity and ∂∂ = 0.
inline Rewritten rewrite(const RegularComplex& X, const ComplexRule& q) {
  Rewritten out;
  const auto refuse = [&](std::string why) { out.applied = false; out.reason = std::move(why); return out; };
  if (!q.refused.empty()) return refuse(q.refused);
  const rw::Graph G = incidence_graph(X);
  const auto& L = q.rule.left;
  if (q.match.vertex.size() != L.size() || q.match.edge.size() != L.edges.size()) return refuse("the match does not cover the rule's left side");
  std::set<std::size_t> image(q.match.vertex.begin(), q.match.vertex.end());
  if (image.size() != L.size()) return refuse("the match is not injective");
  for (std::size_t v = 0; v < L.size(); v++)
    if (q.match.vertex[v] >= G.size() || !rw::detail::label_fits(L.vertex[v], G.vertex[q.match.vertex[v]])) return refuse("the match does not preserve a cell's dimension and type");
  for (std::size_t e = 0; e < L.edges.size(); e++) {
    const auto& le = L.edges[e];
    if (q.match.edge[e] >= G.edges.size()) return refuse("the match names an incidence that is not one");
    const auto& he = G.edges[q.match.edge[e]];
    if (he.u != q.match.vertex[le.u] || he.v != q.match.vertex[le.v] || he.label != le.label) return refuse("the match does not preserve an incidence");
  }
  try {
    q.rule.validate();
  } catch (const std::exception& e) {
    return refuse(e.what());
  }
  if (!rw::dangling_ok(q.rule, G, q.match)) {
    // Name a cell the rule would delete that bounds a cell the rule does not name.
    std::vector<char> in_image(G.edges.size(), 0);
    for (std::size_t e : q.match.edge) in_image[e] = 1;
    std::string which;
    for (std::size_t v = 0; v < L.size() && which.empty(); v++) {
      if (q.rule.right_of(v) != rw::none) continue;
      const std::size_t h = q.match.vertex[v];
      for (std::size_t e = 0; e < G.edges.size(); e++)
        if (!in_image[e] && (G.edges[e].u == h || G.edges[e].v == h)) {
          const std::size_t other = G.edges[e].u == h ? G.edges[e].v : G.edges[e].u;
          which = detail::named(X, h) + " is a face of " + detail::named(X, other) + ", which the rule does not name";
          break;
        }
    }
    return refuse("the gluing condition fails: " + which);
  }
  if (q.rule.condition && !q.rule.condition(G, q.match)) return refuse("the rule's application condition fails");
  std::vector<std::size_t> trace;
  const rw::Graph H = rw::apply(q.rule, G, q.match, &trace);
  // Ids: kept vertices keep theirs, created ones get fresh ids.
  std::vector<std::uint64_t> ids(H.size(), 0);
  for (std::size_t i = 0; i < X.size(); i++)
    if (trace[i] != rw::none) ids[trace[i]] = X.cells[i].id;
  std::uint64_t next = X.next_id;
  for (auto& id : ids)
    if (id == 0) id = next++, out.created.push_back(id);
  std::vector<std::size_t> order;
  try {
    out.result = from_incidence_graph(H, ids, next, &order);
  } catch (const std::exception& e) {
    return refuse(e.what());
  }
  out.result.name = X.name;
  const Admissibility a = admissible(out.result);
  if (!a.regular) return refuse("the result is not a regular complex: " + a.reason);
  try {
    chain_complex(out.result).validate();
  } catch (const std::exception& e) {
    return refuse(std::string("the result's boundary does not square to zero: ") + e.what());
  }
  out.trace.assign(X.size(), no_cell);
  for (std::size_t i = 0; i < X.size(); i++)
    if (trace[i] != rw::none) out.trace[i] = order[trace[i]];
  // What it read and wrote: the match's image; what it deleted, created and relabelled.
  auto& fp = out.footprint;
  for (std::size_t v = 0; v < L.size(); v++) fp.reads.insert(detail::cell_item(X.cells[q.match.vertex[v]].id));
  for (std::size_t e = 0; e < L.edges.size(); e++) {
    const auto& he = G.edges[q.match.edge[e]];
    fp.reads.insert(detail::incidence_item(X.cells[he.u].id, X.cells[he.v].id));
  }
  for (std::size_t v = 0; v < L.size(); v++) {
    const std::size_t h = q.match.vertex[v], r = q.rule.right_of(v);
    if (r == rw::none) {
      out.deleted.push_back(X.cells[h].id);
      fp.writes.insert(detail::cell_item(X.cells[h].id));
    } else if (q.rule.right.vertex[r] != rw::any_label && q.rule.right.vertex[r] != G.vertex[h]) {
      out.relabelled.push_back(X.cells[h].id);
      fp.writes.insert(detail::cell_item(X.cells[h].id));
    }
  }
  for (std::size_t e = 0; e < L.edges.size(); e++)
    if (q.rule.right_of_edge(e) == rw::none) {
      const auto& he = G.edges[q.match.edge[e]];
      fp.writes.insert(detail::incidence_item(X.cells[he.u].id, X.cells[he.v].id));
    }
  for (std::uint64_t id : out.created) fp.writes.insert(detail::cell_item(id));
  const std::set<std::uint64_t> made(out.created.begin(), out.created.end());
  for (const auto& c : out.result.cells)
    for (std::size_t f : c.facets)
      if (made.count(c.id) || made.count(out.result.cells[f].id)) fp.writes.insert(detail::incidence_item(c.id, out.result.cells[f].id));
  out.applied = true;
  return out;
}

}  // namespace morphe::topology
