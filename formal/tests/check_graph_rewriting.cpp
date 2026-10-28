// Checks of formal/graph_rewriting.hpp: automorphism groups of familiar graphs; the count of
// matches as |Aut L| times the occurrences; the dangling condition; relabelling with the trace of
// the vertices; one step of string rewriting (formal/grammar.hpp) and one step of graph rewriting
// on the words drawn as paths giving the same forms; the derivations of chains joined from n
// labelled units, n!·(n−1)! in all, (n−1)! for each of the n! chains, one orbit under the
// permutations of the labels; and the equivariance of derivations under a permutation of labels
// when the rule set is closed under it, and its failure when it is not. The comparison with string
// rewriting runs where formal/grammar.hpp is present (decision 0006).

#include <morphe/check.hpp>
#include <morphe/formal/graph_rewriting.hpp>
#if __has_include(<morphe/formal/grammar.hpp>)
#include <morphe/formal/grammar.hpp>
#define MORPHE_CHECK_STRINGS 1
#endif

#include <algorithm>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace rw = morphe::formal::rewriting;
#ifdef MORPHE_CHECK_STRINGS
namespace mf = morphe::formal;
#endif

namespace {

rw::Graph cycle(std::size_t n, bool directed = false) {
  rw::Graph G;
  G.directed = directed;
  for (std::size_t i = 0; i < n; i++) G.add_vertex("x");
  for (std::size_t i = 0; i < n; i++) G.add_edge(i, (i + 1) % n, "-");
  return G;
}
rw::Graph path(std::size_t n) {
  rw::Graph G;
  for (std::size_t i = 0; i < n; i++) G.add_vertex("x");
  for (std::size_t i = 0; i + 1 < n; i++) G.add_edge(i, i + 1, "-");
  return G;
}
rw::Graph complete(std::size_t n) {
  rw::Graph G;
  for (std::size_t i = 0; i < n; i++) G.add_vertex("x");
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) G.add_edge(i, j, "-");
  return G;
}
std::size_t factorial(std::size_t n) { return n <= 1 ? 1 : n * factorial(n - 1); }

#ifdef MORPHE_CHECK_STRINGS
// A word of the grammar's symbols as a directed path between the sentinels ⊢ and ⊣.
rw::Graph word_graph(const mf::Form& w) {
  rw::Graph G;
  G.directed = true;
  std::size_t prev = G.add_vertex("⊢");
  for (auto s : w) {
    const std::size_t v = G.add_vertex(std::to_string(s));
    G.add_edge(prev, v, "next");
    prev = v;
  }
  G.add_edge(prev, G.add_vertex("⊣"), "next");
  return G;
}
// The word read back from a path graph.
mf::Form read_word(const rw::Graph& G) {
  std::size_t v = 0;
  for (std::size_t i = 0; i < G.size(); i++)
    if (G.vertex[i] == "⊢") v = i;
  mf::Form w;
  for (;;) {
    std::size_t next = rw::none;
    for (const auto& e : G.edges)
      if (e.u == v) next = e.v;
    if (next == rw::none || G.vertex[next] == "⊣") break;
    w.push_back(static_cast<mf::Symbol>(std::stoul(G.vertex[next])));
    v = next;
  }
  return w;
}
// The string rule α → β as a graph rule: p → α → q replaced by p → β → q, p and q kept with any label.
rw::Rule word_rule(const mf::Rule& r) {
  rw::Rule q;
  q.left.directed = q.right.directed = true;
  const std::size_t lp = q.left.add_vertex(rw::any_label), rp = q.right.add_vertex(rw::any_label);
  std::size_t prev = lp;
  for (auto s : r.lhs) {
    const std::size_t v = q.left.add_vertex(std::to_string(s));
    q.left.add_edge(prev, v, "next");
    prev = v;
  }
  const std::size_t lq = q.left.add_vertex(rw::any_label);
  q.left.add_edge(prev, lq, "next");
  prev = rp;
  for (auto s : r.rhs) {
    const std::size_t v = q.right.add_vertex(std::to_string(s));
    q.right.add_edge(prev, v, "next");
    prev = v;
  }
  const std::size_t rq = q.right.add_vertex(rw::any_label);
  q.right.add_edge(prev, rq, "next");
  q.kept_vertices = {{lp, rp}, {lq, rq}};
  q.validate();
  return q;
}
#endif

// n labelled units, each with a free head h and a free tail t; the rule joins the tail of one
// chain to the head of another chain.
rw::Graph units(std::size_t n) {
  rw::Graph G;
  G.directed = true;
  for (std::size_t i = 0; i < n; i++) {
    const std::size_t u = G.add_vertex("u" + std::to_string(i + 1));
    G.add_edge(u, G.add_vertex("h•"), "head");
    G.add_edge(u, G.add_vertex("t•"), "tail");
  }
  return G;
}
rw::Rule join_rule() {
  rw::Rule q;
  q.name = "join";
  q.left.directed = q.right.directed = true;
  const std::size_t lt = q.left.add_vertex("t•"), lh = q.left.add_vertex("h•");
  const std::size_t rt = q.right.add_vertex("t"), rh = q.right.add_vertex("h");
  q.right.add_edge(rt, rh, "link");
  q.kept_vertices = {{lt, rt}, {lh, rh}};
  q.condition = [](const rw::Graph& H, const rw::Match& m) {
    const auto c = H.components();
    return c[m.vertex[0]] != c[m.vertex[1]];
  };
  q.validate();
  return q;
}

}  // namespace

int main() {
  morphe::Suite s;

  // ---- automorphisms ------------------------------------------------------------------------------------------------------
  {
    bool ok = true;
    for (std::size_t n = 3; n <= 8; n++) ok = ok && rw::automorphisms(cycle(n)) == 2 * n && rw::automorphisms(cycle(n, true)) == n;
    s.require("|Aut C_n| = 2n undirected, n directed (n = 3..8)", ok);
    s.equal("|Aut K_5| = 5!", rw::automorphisms(complete(5)), std::size_t{120});
    s.equal("|Aut P_6| = 2", rw::automorphisms(path(6)), std::size_t{2});
    rw::Graph star;
    const std::size_t c = star.add_vertex("C");
    for (const char* l : {"H", "F", "Cl", "Br"}) star.add_edge(c, star.add_vertex(l), "-");
    s.equal("|Aut| of a star with four distinct leaves = 1", rw::automorphisms(star), std::size_t{1});
    for (auto& l : star.vertex) if (l != "C") l = "H";
    s.equal("|Aut| of a star with four equal leaves = 4!", rw::automorphisms(star), std::size_t{24});
  }

  // ---- matches = |Aut L| × occurrences ------------------------------------------------------------------------------------
  {
    // In K_n: paths of two edges, n·C(n−1, 2) occurrences; triangles, C(n, 3).
    bool ok = true;
    for (std::size_t n = 4; n <= 7; n++) {
      const std::size_t paths = n * (n - 1) * (n - 2) / 2, triangles = n * (n - 1) * (n - 2) / 6;
      ok = ok && rw::count_matches(path(3), complete(n)) == 2 * paths;
      ok = ok && rw::count_matches(cycle(3), complete(n)) == 6 * triangles;
    }
    s.require("matches of P_3 and C_3 in K_n are |Aut| × occurrences (n = 4..7)", ok);
    // Parallel edges: a double edge matched into a triple edge in 3·2·2 = 12 ways (2 orientations
    // of the vertices, 3·2 injections of the edges).
    rw::Graph two, three;
    two.add_vertex("x"), two.add_vertex("x"), two.add_edge(0, 1, "-"), two.add_edge(0, 1, "-");
    three = two;
    three.add_edge(1, 0, "-");
    s.equal("a double edge into a triple edge: 12 matches", rw::count_matches(two, three), std::size_t{12});
  }

  // ---- the dangling condition and relabelling -----------------------------------------------------------------------------
  {
    rw::Rule del;
    del.name = "delete a vertex";
    del.left.add_vertex("x");
    del.validate();
    rw::Graph H = path(3);
    H.add_vertex("x");  // isolated
    const auto ds = rw::derivations({del}, H);
    s.require("a vertex may be deleted only where no edge dangles", ds.size() == 1 && ds[0].match.vertex[0] == 3 && ds[0].result.size() == 3);

    // Relabel the middle vertex of a path a - b - c and move an edge: a - b, b - c  ⇒  a - b', a - c.
    rw::Rule move;
    move.name = "move";
    const std::size_t la = move.left.add_vertex("a"), lb = move.left.add_vertex("b"), lc = move.left.add_vertex("c");
    move.left.add_edge(la, lb, "-");
    move.left.add_edge(lb, lc, "-");
    const std::size_t ra = move.right.add_vertex("a"), rb = move.right.add_vertex("b'"), rc = move.right.add_vertex("c");
    move.right.add_edge(ra, rb, "-");
    move.right.add_edge(ra, rc, "-");
    move.kept_vertices = {{la, ra}, {lb, rb}, {lc, rc}};
    move.kept_edges = {{0, 0}};
    move.validate();
    rw::Graph P;
    P.add_vertex("a"), P.add_vertex("b"), P.add_vertex("c"), P.add_vertex("d");
    P.add_edge(0, 1, "-"), P.add_edge(1, 2, "-"), P.add_edge(2, 3, "-");
    const auto d = rw::derivations({move}, P);
    bool ok = d.size() == 1;
    if (ok) {
      const auto& R = d[0].result;
      ok = R.size() == 4 && R.vertex[d[0].trace[1]] == "b'" && R.edges.size() == 3 && d[0].trace == std::vector<std::size_t>{0, 1, 2, 3};
      std::size_t ac = 0;
      for (const auto& e : R.edges) ac += (R.vertex[e.u] == "a" && R.vertex[e.v] == "c") || (R.vertex[e.u] == "c" && R.vertex[e.v] == "a");
      ok = ok && ac == 1;
    }
    s.require("a kept vertex is relabelled and keeps its identity in the trace", ok);
  }

  // ---- string rewriting as graph rewriting ---------------------------------------------------------------------------------
#ifndef MORPHE_CHECK_STRINGS
  s.skip("string rewriting against graph rewriting: formal/grammar.hpp is not here");
#else
  {
    // A monotone grammar of {aⁿbⁿcⁿ} (Hopcroft and Ullman 1979, Example 9.5, in the form
    // S → aSBC | aBC, CB → BC, aB → ab, bB → bb, bC → bc, cC → cc).
    mf::Grammar g;
    g.terminals = 3;  // a b c
    const auto S = g.add_nonterminal("S"), B = g.add_nonterminal("B"), C = g.add_nonterminal("C");
    g.add(S, {0, S, B, C});
    g.add(S, {0, B, C});
    g.add({C, B}, {B, C});
    g.add({0, B}, {0, 1});
    g.add({1, B}, {1, 1});
    g.add({1, C}, {1, 2});
    g.add({2, C}, {2, 2});
    std::vector<rw::Rule> rules;
    for (const auto& r : g.rules) rules.push_back(word_rule(r));
    std::mt19937 rng(7);
    bool ok = true;
    std::size_t compared = 0;
    for (int trial = 0; trial < 300; trial++) {
      mf::Form w(1 + rng() % 7);
      for (auto& x : w) x = static_cast<mf::Symbol>(rng() % g.symbols());
      auto a = mf::one_step(g, w);
      std::vector<mf::Form> b;
      for (const auto& d : rw::derivations(rules, word_graph(w))) b.push_back(read_word(d.result));
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      ok = ok && a == b;
      compared += a.size();
    }
    s.require("one step of string rewriting equals one step of graph rewriting on paths (300 random forms)", ok && compared > 50,
              std::to_string(compared) + " derived forms compared");
  }
#endif

  // ---- derivations of chains, and their orbit ------------------------------------------------------------------------------
  {
    const std::vector<rw::Rule> rules{join_rule()};
    for (std::size_t n : {3, 4}) {
      const auto D = rw::explore(rules, units(n), n - 1);
      std::vector<std::size_t> chains;
      for (std::size_t c = 0; c < D.graphs.size(); c++)
        if (D.depth[c] == n - 1) chains.push_back(c);
      bool each = chains.size() == factorial(n);
      unsigned long long total = 0;
      for (auto c : chains) {
        each = each && D.derivations_to(c) == factorial(n - 1);
        total += D.derivations_to(c);
      }
      s.require("n = " + std::to_string(n) + ": n! chains, each derived in (n−1)! ways", each && total == factorial(n) * factorial(n - 1),
                std::to_string(chains.size()) + " chains, " + std::to_string(total) + " derivations");
      // The symmetric group on the unit labels: generated by a transposition and an n-cycle.
      std::vector<rw::Relabelling> gens{{{"u1", "u2"}, {"u2", "u1"}}, {}};
      for (std::size_t i = 0; i < n; i++) gens[1]["u" + std::to_string(i + 1)] = "u" + std::to_string((i + 1) % n + 1);
      std::vector<rw::Graph> terminal;
      for (auto c : chains) terminal.push_back(D.graphs[c]);
      const auto orb = rw::orbits(terminal, gens);
      s.require("n = " + std::to_string(n) + ": the chains form one orbit of S_n", std::all_of(orb.begin(), orb.end(), [](std::size_t o) { return o == 0; }));
      bool eq = true;
      for (std::size_t c = 0; c < D.graphs.size(); c++)
        for (const auto& g : gens) eq = eq && rw::equivariant(rules, D.graphs[c], g);
      s.require("n = " + std::to_string(n) + ": derivations commute with the relabellings, at every derivable graph", eq);
    }
  }

  // ---- classes of graphs of many identical parts ------------------------------------------------------------------------
  {
    // Twelve identical triangles and ten isolated vertices: one exact class, the same for any
    // numbering of the vertices; one triangle fewer gives another class.
    const auto parts = [](std::size_t triangles, std::size_t isolated, std::uint32_t seed) {
      rw::Graph G;
      for (std::size_t t = 0; t < triangles; t++) {
        const std::size_t a = G.add_vertex("x"), b = G.add_vertex("x"), c = G.add_vertex("x");
        G.add_edge(a, b, "-"), G.add_edge(b, c, "-"), G.add_edge(c, a, "-");
      }
      for (std::size_t i = 0; i < isolated; i++) G.add_vertex("w");
      // Renumber the vertices by a seeded permutation.
      std::vector<std::size_t> perm(G.size());
      for (std::size_t i = 0; i < perm.size(); i++) perm[i] = i;
      std::shuffle(perm.begin(), perm.end(), std::mt19937(seed));
      rw::Graph H;
      H.vertex.resize(G.size());
      for (std::size_t v = 0; v < G.size(); v++) H.vertex[perm[v]] = G.vertex[v];
      for (const auto& e : G.edges) H.add_edge(perm[e.u], perm[e.v], e.label);
      return H;
    };
    bool ok = true;
    try {
      ok = rw::code(parts(12, 10, 1)) == rw::code(parts(12, 10, 2)) && rw::code(parts(12, 10, 1)) != rw::code(parts(11, 13, 1));
    } catch (const std::exception&) {
      ok = false;
    }
    s.require("twelve identical triangles and ten isolated vertices: one exact class under renumbering", ok);
  }

  // ---- closure is needed ---------------------------------------------------------------------------------------------------
  {
    // Rules that turn an x into a+ or a−, and a reflection exchanging the signs.
    const auto turn = [](const std::string& to) {
      rw::Rule q;
      q.name = "x→" + to;
      q.left.add_vertex("x");
      q.right.add_vertex(to);
      q.kept_vertices = {{0, 0}};
      q.validate();
      return q;
    };
    const rw::Relabelling p{{"a+", "a-"}, {"a-", "a+"}};
    rw::validate(p);
    const std::vector<rw::Rule> both{turn("a+"), turn("a-")}, one{turn("a+")};
    rw::Graph H = path(3);
    s.require("a rule set closed under the reflection: closed and equivariant", rw::closed_under(both, p) && rw::equivariant(both, H, p));
    rw::Graph Hp = H;
    Hp.vertex[1] = "a+";
    s.require("a rule set not closed: neither closed nor equivariant", !rw::closed_under(one, p) && !rw::equivariant(one, Hp, p));
    rw::Rule renumbered = turn("a+");
    renumbered.name = "another name";
    s.require("rule classes ignore names; a rule and its inverse differ",
              rw::rule_code(renumbered) == rw::rule_code(turn("a+")) && rw::rule_code(turn("a+")) != rw::rule_code(turn("a-")));
  }

  return s.done();
}
