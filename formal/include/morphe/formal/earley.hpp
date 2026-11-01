// Earley's parser for every context-free grammar, ε-rules included, with a parse tree read back.
//
// An item [A → α • β, i] in the set S_j says that α derives w_{i+1..j} and that A → αβ may
// continue the parse (Earley 1970, An efficient context-free parsing algorithm, Comm. ACM 13,
// 94–102). The sets are closed under
//
//   predict    [A → α • B β, i] ∈ S_j, B → γ          ⟹  [B → • γ, j] ∈ S_j,
//   scan       [A → α • a β, i] ∈ S_j, w_{j+1} = a     ⟹  [A → α a • β, i] ∈ S_{j+1},
//   complete   [B → γ •, k] ∈ S_j, [A → α • B β, i] ∈ S_k  ⟹  [A → α B • β, i] ∈ S_j,
//
// and w ∈ L exactly when [S → γ •, 0] ∈ S_n. Completing an item whose origin is j itself can miss
// items that ask for its left side later in the same set; predicting a nullable B also advances
// the dot over B, which repairs this (Aycock and Horspool 2002, Practical Earley parsing,
// The Computer Journal 45, 620–630). Time O(n³ |Γ|) in general, O(n²) on unambiguous grammars,
// space O(n² |Γ|) items.
//
// Every item remembers one way it was made, so that one parse tree of w is read back from the
// accepting item; ε-subtrees of nullable nonterminals come from a fixed derivation of ε each.
//
// Diagnostics. S_j is non-empty exactly when w_{1..j} is a viable prefix, S ⇒* w_{1..j} γ for
// some sentential form γ; when every nonterminal is productive, exactly when w_{1..j} is a
// prefix of a word of L. A rejected word's furthest viable prefix is therefore the largest j with
// S_j ≠ ∅, and the terminals that may follow it are those after the dots of S_j's items: the
// chart states exactly how far the word is right and what it lacks, which is what a reader of
// an error needs.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "grammar.hpp"

namespace morphe::formal {

// A parse tree: a node is a rule applied over a span, its children in the order of its right
// side; a terminal child is a leaf (rule = −1) over one position.
struct ParseTree {
  long rule = -1;  // −1: a terminal leaf
  std::size_t begin = 0, end = 0;
  std::vector<ParseTree> children;
};

struct EarleyResult {
  bool accepts = false;
  std::size_t items = 0;  // the size of the chart, Σ_j |S_j|
  ParseTree tree;         // one parse, when accepted and asked for
  std::size_t furthest = 0;               // the length of the longest viable prefix of w
  std::vector<Symbol> expected;           // the terminals that may follow it, ascending
  std::vector<std::pair<std::uint32_t, std::uint32_t>> frontier;  // (rule, dot) of the items awaiting them
};

namespace detail {

struct EarleyItem {
  std::uint32_t rule, dot;
  std::size_t origin;
  bool operator==(const EarleyItem& o) const { return rule == o.rule && dot == o.dot && origin == o.origin; }
};
struct EarleyItemHash {
  std::size_t operator()(const EarleyItem& x) const {
    return (static_cast<std::size_t>(x.rule) * 1000003u + x.dot) * 1000033u + x.origin;
  }
};
// How an item was made: from a predecessor in a set (by scan, completion or a nullable advance),
// with the child that moved the dot.
struct Back {
  enum Kind : std::uint8_t { predicted, scanned, completed, nulled } kind = predicted;
  std::size_t pred_set = 0, pred_index = 0;
  std::size_t child_set = 0, child_index = 0;  // completed: the finished child item
  Symbol nulled_symbol = 0;
};

}  // namespace detail

class Earley {
 public:
  explicit Earley(const Grammar& g) : g_(g) {
    g.validate();
    if (!g.context_free()) throw std::invalid_argument("morphe::formal::Earley: the grammar is not context-free");
    by_lhs_.assign(g.nonterminals.size(), {});
    for (std::uint32_t r = 0; r < g.rules.size(); r++) by_lhs_[g.index(g.rules[r].lhs[0])].push_back(r);
    nullable_ = g.nullable();
    epsilon_rule_.assign(g.nonterminals.size(), -1);
    // A rule deriving ε for each nullable nonterminal, chosen so that the ε-derivations are well
    // founded: in the order the fixpoint found them.
    std::vector<char> done(g.nonterminals.size(), 0);
    for (bool changed = true; changed;) {
      changed = false;
      for (std::uint32_t r = 0; r < g.rules.size(); r++) {
        const auto& rule = g.rules[r];
        const std::size_t a = g.index(rule.lhs[0]);
        if (done[a]) continue;
        bool all = true;
        for (Symbol s : rule.rhs) all = all && !g.terminal(s) && done[g.index(s)];
        if (all) { done[a] = 1; epsilon_rule_[a] = r; changed = true; }
      }
    }
  }

  // Parses w from the start S (nonterminal 0), or from the nonterminal of index `start`.
  EarleyResult parse(const std::vector<Symbol>& w, bool want_tree = false, std::size_t start = 0) const {
    if (start >= g_.nonterminals.size()) throw std::invalid_argument("morphe::formal::Earley: no such start");
    const std::size_t n = w.size();
    std::vector<std::vector<detail::EarleyItem>> sets(n + 1);
    std::vector<std::vector<detail::Back>> backs(n + 1);
    std::vector<std::unordered_map<detail::EarleyItem, std::size_t, detail::EarleyItemHash>> where(n + 1);
    const auto add = [&](std::size_t j, detail::EarleyItem it, detail::Back b) {
      if (where[j].count(it)) return;
      where[j].emplace(it, sets[j].size());
      sets[j].push_back(it);
      backs[j].push_back(b);
    };
    for (std::uint32_t r : by_lhs_[start]) add(0, {r, 0, 0}, {});
    for (std::size_t j = 0; j <= n; j++) {
      for (std::size_t k = 0; k < sets[j].size(); k++) {
        const detail::EarleyItem it = sets[j][k];
        const auto& rhs = g_.rules[it.rule].rhs;
        if (it.dot < rhs.size()) {
          const Symbol x = rhs[it.dot];
          if (g_.terminal(x)) {
            if (j < n && w[j] == x) add(j + 1, {it.rule, it.dot + 1, it.origin}, {detail::Back::scanned, j, k});
          } else {
            for (std::uint32_t r : by_lhs_[g_.index(x)]) add(j, {r, 0, j}, {});
            if (nullable_[g_.index(x)]) {
              detail::Back b{detail::Back::nulled, j, k};
              b.nulled_symbol = x;
              add(j, {it.rule, it.dot + 1, it.origin}, b);
            }
          }
        } else {
          const Symbol a = g_.rules[it.rule].lhs[0];
          const std::size_t o = it.origin;
          for (std::size_t p = 0; p < sets[o].size(); p++) {
            const detail::EarleyItem parent = sets[o][p];
            const auto& prhs = g_.rules[parent.rule].rhs;
            if (parent.dot < prhs.size() && prhs[parent.dot] == a) {
              detail::Back b{detail::Back::completed, o, p};
              b.child_set = j;
              b.child_index = k;
              add(j, {parent.rule, parent.dot + 1, parent.origin}, b);
            }
          }
        }
      }
    }
    EarleyResult out;
    for (const auto& s : sets) out.items += s.size();
    for (std::size_t j = 0; j <= n; j++)
      if (!sets[j].empty()) out.furthest = j;
    for (const auto& it : sets[out.furthest]) {
      const auto& rhs = g_.rules[it.rule].rhs;
      if (it.dot < rhs.size() && g_.terminal(rhs[it.dot])) {
        out.expected.push_back(rhs[it.dot]);
        out.frontier.emplace_back(it.rule, it.dot);
      }
    }
    std::sort(out.expected.begin(), out.expected.end());
    out.expected.erase(std::unique(out.expected.begin(), out.expected.end()), out.expected.end());
    std::sort(out.frontier.begin(), out.frontier.end());
    out.frontier.erase(std::unique(out.frontier.begin(), out.frontier.end()), out.frontier.end());
    const Symbol goal = g_.nonterminal(start);
    for (std::size_t k = 0; k < sets[n].size(); k++) {
      const auto& it = sets[n][k];
      if (it.origin == 0 && it.dot == g_.rules[it.rule].rhs.size() && g_.rules[it.rule].lhs[0] == goal) {
        out.accepts = true;
        if (want_tree) out.tree = build(sets, backs, n, k);
        break;
      }
    }
    return out;
  }

  bool accepts(const std::vector<Symbol>& w) const { return parse(w).accepts; }

 private:
  // The subtree of a completed item at set j, index k.
  ParseTree build(const std::vector<std::vector<detail::EarleyItem>>& sets, const std::vector<std::vector<detail::Back>>& backs,
                  std::size_t j, std::size_t k) const {
    const detail::EarleyItem top = sets[j][k];
    ParseTree node;
    node.rule = top.rule;
    node.begin = top.origin;
    node.end = j;
    std::vector<ParseTree> rev;
    std::size_t sj = j, sk = k;
    while (sets[sj][sk].dot > 0) {
      const detail::Back& b = backs[sj][sk];
      if (b.kind == detail::Back::scanned) {
        ParseTree leaf;
        leaf.begin = b.pred_set;
        leaf.end = b.pred_set + 1;
        rev.push_back(leaf);
      } else if (b.kind == detail::Back::completed) {
        rev.push_back(build(sets, backs, b.child_set, b.child_index));
      } else if (b.kind == detail::Back::nulled) {
        rev.push_back(epsilon_tree(b.nulled_symbol, sj));
      } else {
        throw std::logic_error("morphe::formal::Earley: a predicted item has its dot advanced");
      }
      sj = b.pred_set;
      sk = b.pred_index;
    }
    node.children.assign(rev.rbegin(), rev.rend());
    return node;
  }
  ParseTree epsilon_tree(Symbol a, std::size_t at) const {
    ParseTree t;
    t.rule = epsilon_rule_[g_.index(a)];
    t.begin = t.end = at;
    for (Symbol s : g_.rules[t.rule].rhs) t.children.push_back(epsilon_tree(s, at));
    return t;
  }

  const Grammar& g_;
  std::vector<std::vector<std::uint32_t>> by_lhs_;
  std::vector<char> nullable_;
  std::vector<long> epsilon_rule_;
};

// The word a parse tree spells, for checking it against the input.
inline std::vector<Symbol> yield(const Grammar& g, const ParseTree& t, const std::vector<Symbol>& w) {
  std::vector<Symbol> out;
  if (t.rule < 0) { out.push_back(w[t.begin]); return out; }
  const auto& rhs = g.rules[static_cast<std::size_t>(t.rule)].rhs;
  if (rhs.size() != t.children.size()) throw std::logic_error("morphe::formal::yield: a node's children do not match its rule");
  for (std::size_t c = 0; c < rhs.size(); c++) {
    if (g.terminal(rhs[c]) != (t.children[c].rule < 0)) throw std::logic_error("morphe::formal::yield: a child of the wrong kind");
    if (g.terminal(rhs[c]) && w[t.children[c].begin] != rhs[c]) throw std::logic_error("morphe::formal::yield: a leaf of another letter");
    if (!g.terminal(rhs[c]) && g.rules[static_cast<std::size_t>(t.children[c].rule)].lhs[0] != rhs[c])
      throw std::logic_error("morphe::formal::yield: a child of another nonterminal");
    const auto part = yield(g, t.children[c], w);
    out.insert(out.end(), part.begin(), part.end());
  }
  return out;
}

}  // namespace morphe::formal
