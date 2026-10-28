// Weighted context-free grammars over a semiring: the inside and outside values of every span,
// the expected number of uses of every rule, and the best parse.
//
// A weighted grammar gives each rule a weight in a semiring K; a derivation weighs the product
// of its rules', and a word the sum over its derivations (Goodman 1999, Semiring parsing,
// Computational Linguistics 25, 573–605). With K Boolean this is recognition, with Counting the
// number of parse trees, with MaxPlus the best derivation (Viterbi), and with the reals the
// probability of a stochastic grammar (Baker 1979; Lari and Young 1990, The estimation of
// stochastic context-free grammars using the Inside-Outside algorithm, Computer Speech and
// Language 4, 35–56). For a span (i, j) of w and a nonterminal A,
//
//   inside   I(A, i, j) = Σ over derivations A ⇒* w_{i+1..j} of their weights,
//   outside  O(A, i, j) = Σ over derivations S ⇒* w_{1..i} A w_{j+1..n} of their weights,
//
// and, in a commutative semiring, the weight of the derivations using rule A → BC over
// (i, k, j) is O(A, i, j) ⊗ w(A → BC) ⊗ I(B, i, k) ⊗ I(C, k, j); summed and divided by the word's
// inside value it is the expected number of uses, the E-step of the inside–outside algorithm.
//
// The grammar is put in a binary normal form first, by a map that is a bijection on derivations
// and keeps their weights: a terminal inside a longer right side becomes a preterminal P_a → a
// of weight 1, and A → Y_1 … Y_k (k ≥ 3) becomes A → Y_1 H_1, H_1 → Y_2 H_2, …, H_{k−2} →
// Y_{k−1} Y_k, the first carrying the rule's weight and the rest weight 1 (the binarisation of
// Chomsky's normal form, Hopcroft and Ullman 1979, Theorem 4.5, with unary rules kept). The
// grammar must be context-free and ε-free, and its unary rules A → B must form no cycle, so
// that the unary closure of each span is a finite sum taken in a topological order; otherwise
// std::invalid_argument. Time O(n³ |R_binary| + n² |R_unary|), space O(n² |V_binary|).

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

#include "earley.hpp"
#include "grammar.hpp"
#include "semiring.hpp"

namespace morphe::formal {

// The binary normal form of a context-free, ε-free grammar, with each binary-form rule's
// origin (the original rule whose weight it carries, or −1 for a helper or preterminal rule).
struct BinaryGrammar {
  struct Lexical { std::uint32_t lhs; Symbol terminal; long origin; };
  struct Unary { std::uint32_t lhs, child; long origin; };
  struct Binary { std::uint32_t lhs, left, right; long origin; };
  std::size_t terminals = 0, originals = 0;  // |Σ|, and how many nonterminals are the original ones
  std::size_t nonterminals = 0;
  std::vector<Lexical> lexical;
  std::vector<Unary> unary;
  std::vector<Binary> binary;
  std::vector<std::uint32_t> unary_order;  // nonterminals, each after every B with a unary A → B
  std::vector<char> helper;                // a helper H_i of a long rule
  std::vector<long> preterminal_of;        // the terminal a of P_a, else −1

  static BinaryGrammar from(const Grammar& g) {
    g.validate();
    if (!g.context_free()) throw std::invalid_argument("morphe::formal::BinaryGrammar: the grammar is not context-free");
    BinaryGrammar b;
    b.terminals = g.terminals;
    b.originals = g.nonterminals.size();
    b.nonterminals = b.originals;
    b.helper.assign(b.originals, 0);
    b.preterminal_of.assign(b.originals, -1);
    std::vector<long> pre(g.terminals, -1);
    const auto fresh = [&](bool is_helper, long terminal) {
      b.helper.push_back(is_helper);
      b.preterminal_of.push_back(terminal);
      return static_cast<std::uint32_t>(b.nonterminals++);
    };
    const auto as_nonterminal = [&](Symbol s) -> std::uint32_t {
      if (!g.terminal(s)) return static_cast<std::uint32_t>(g.index(s));
      if (pre[s] < 0) {
        pre[s] = fresh(false, s);
        b.lexical.push_back({static_cast<std::uint32_t>(pre[s]), s, -1});
      }
      return static_cast<std::uint32_t>(pre[s]);
    };
    for (std::size_t r = 0; r < g.rules.size(); r++) {
      const auto a = static_cast<std::uint32_t>(g.index(g.rules[r].lhs[0]));
      const auto& rhs = g.rules[r].rhs;
      if (rhs.empty()) throw std::invalid_argument("morphe::formal::BinaryGrammar: an ε-rule");
      if (rhs.size() == 1) {
        if (g.terminal(rhs[0])) b.lexical.push_back({a, rhs[0], static_cast<long>(r)});
        else b.unary.push_back({a, static_cast<std::uint32_t>(g.index(rhs[0])), static_cast<long>(r)});
        continue;
      }
      std::vector<std::uint32_t> y;
      for (Symbol s : rhs) y.push_back(as_nonterminal(s));
      std::uint32_t lhs = a;
      long origin = static_cast<long>(r);
      for (std::size_t k = 0; k + 2 < y.size(); k++) {
        const std::uint32_t h = fresh(true, -1);
        b.binary.push_back({lhs, y[k], h, origin});
        lhs = h;
        origin = -1;
      }
      b.binary.push_back({lhs, y[y.size() - 2], y.back(), origin});
    }
    // The unary rules' topological order (Kahn), children before parents; a cycle is refused.
    std::vector<std::vector<std::uint32_t>> parents(b.nonterminals);
    std::vector<std::size_t> waiting(b.nonterminals, 0);
    for (const auto& u : b.unary) { parents[u.child].push_back(u.lhs); waiting[u.lhs]++; }
    std::vector<std::uint32_t> ready;
    for (std::uint32_t x = 0; x < b.nonterminals; x++)
      if (!waiting[x]) ready.push_back(x);
    while (!ready.empty()) {
      const std::uint32_t x = ready.back();
      ready.pop_back();
      b.unary_order.push_back(x);
      for (std::uint32_t p : parents[x])
        if (--waiting[p] == 0) ready.push_back(p);
    }
    if (b.unary_order.size() != b.nonterminals) throw std::invalid_argument("morphe::formal::BinaryGrammar: the unary rules form a cycle");
    return b;
  }
};

// Inside and outside tables of one word, indexed [A][i][j] over 0 ≤ i < j ≤ n.
template <Semiring S>
class Chart {
 public:
  using T = typename S::value_type;
  // An array, not a vector, so that the Boolean semiring's cells are true bools with references.
  Chart(const BinaryGrammar& b, std::size_t n) : n_(n), size_(b.nonterminals * (n + 1) * (n + 1)), cells_(new T[size_]) {
    for (std::size_t k = 0; k < size_; k++) cells_[k] = S::zero();
  }
  T& at(std::size_t a, std::size_t i, std::size_t j) { return cells_[(a * (n_ + 1) + i) * (n_ + 1) + j]; }
  const T& at(std::size_t a, std::size_t i, std::size_t j) const { return cells_[(a * (n_ + 1) + i) * (n_ + 1) + j]; }
  std::size_t length() const { return n_; }

 private:
  std::size_t n_, size_;
  std::unique_ptr<T[]> cells_;
};

// The weights of the binary form from the original rules' weights (helpers and preterminals 1).
template <Semiring S>
struct BinaryWeights {
  std::vector<typename S::value_type> lexical, unary, binary;
  BinaryWeights(const BinaryGrammar& b, const std::vector<typename S::value_type>& w) {
    const auto of = [&](long origin) { return origin < 0 ? S::one() : w.at(static_cast<std::size_t>(origin)); };
    for (const auto& r : b.lexical) lexical.push_back(of(r.origin));
    for (const auto& r : b.unary) unary.push_back(of(r.origin));
    for (const auto& r : b.binary) binary.push_back(of(r.origin));
  }
};

// The inside values of every span. An optional mask restricts the spans considered (a coarse
// pass's survivors, formal/breaking.hpp): a span with mask false keeps every value at zero.
template <Semiring S>
Chart<S> inside(const BinaryGrammar& b, const std::vector<typename S::value_type>& weights, const std::vector<Symbol>& w,
                const std::vector<char>* mask = nullptr) {
  const std::size_t n = w.size();
  const BinaryWeights<S> wt(b, weights);
  Chart<S> c(b, n);
  const auto open = [&](std::size_t i, std::size_t j) { return !mask || (*mask)[i * (n + 1) + j]; };
  for (std::size_t len = 1; len <= n; len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      if (!open(i, j)) continue;
      if (len == 1)
        for (std::size_t r = 0; r < b.lexical.size(); r++)
          if (b.lexical[r].terminal == w[i]) c.at(b.lexical[r].lhs, i, j) = S::plus(c.at(b.lexical[r].lhs, i, j), wt.lexical[r]);
      for (std::size_t r = 0; r < b.binary.size(); r++) {
        const auto& x = b.binary[r];
        for (std::size_t k = i + 1; k < j; k++)
          c.at(x.lhs, i, j) = S::plus(c.at(x.lhs, i, j), S::times(wt.binary[r], S::times(c.at(x.left, i, k), c.at(x.right, k, j))));
      }
      for (std::uint32_t a : b.unary_order)
        for (std::size_t r = 0; r < b.unary.size(); r++)
          if (b.unary[r].lhs == a) c.at(a, i, j) = S::plus(c.at(a, i, j), S::times(wt.unary[r], c.at(b.unary[r].child, i, j)));
    }
  return c;
}

// The outside values, from the inside chart; S must be commutative for the products to mean the
// weights of the derivations around a span.
template <Semiring S>
Chart<S> outside(const BinaryGrammar& b, const std::vector<typename S::value_type>& weights, const std::vector<Symbol>& w, const Chart<S>& in) {
  static_assert(S::commutative, "morphe::formal::outside: the semiring must be commutative");
  const std::size_t n = w.size();
  const BinaryWeights<S> wt(b, weights);
  Chart<S> o(b, n);
  if (n == 0) return o;
  o.at(0, 0, n) = S::one();
  for (std::size_t len = n; len >= 1; len--) {
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      // Unary rules within the span, parents before children.
      for (auto it = b.unary_order.rbegin(); it != b.unary_order.rend(); ++it)
        for (std::size_t r = 0; r < b.unary.size(); r++)
          if (b.unary[r].lhs == *it) {
            auto& child = o.at(b.unary[r].child, i, j);
            child = S::plus(child, S::times(o.at(*it, i, j), wt.unary[r]));
          }
      for (std::size_t r = 0; r < b.binary.size(); r++) {
        const auto& x = b.binary[r];
        const auto parent = o.at(x.lhs, i, j);
        if (S::equal(parent, S::zero())) continue;
        for (std::size_t k = i + 1; k < j; k++) {
          o.at(x.left, i, k) = S::plus(o.at(x.left, i, k), S::times(parent, S::times(wt.binary[r], in.at(x.right, k, j))));
          o.at(x.right, k, j) = S::plus(o.at(x.right, k, j), S::times(parent, S::times(wt.binary[r], in.at(x.left, i, k))));
        }
      }
    }
    if (len == 1) break;
  }
  return o;
}

// The expected number of uses of every original rule in the derivations of w, weighted by their
// share of the word's inside value: the E-step. Returns the word's inside value too.
template <Semiring S>
std::vector<typename S::value_type> expected_counts(const Grammar& g, const BinaryGrammar& b, const std::vector<typename S::value_type>& weights,
                                                    const std::vector<Symbol>& w, typename S::value_type* total = nullptr) {
  const std::size_t n = w.size();
  const BinaryWeights<S> wt(b, weights);
  const auto in = inside<S>(b, weights, w);
  const auto out = outside<S>(b, weights, w, in);
  std::vector<typename S::value_type> count(g.rules.size(), S::zero());
  const auto z = n ? in.at(0, 0, n) : S::zero();
  if (total) *total = z;
  if (S::equal(z, S::zero())) return count;
  for (std::size_t len = 1; len <= n; len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      if (len == 1)
        for (std::size_t r = 0; r < b.lexical.size(); r++)
          if (b.lexical[r].origin >= 0 && b.lexical[r].terminal == w[i]) {
            auto& c = count[static_cast<std::size_t>(b.lexical[r].origin)];
            c = S::plus(c, S::times(out.at(b.lexical[r].lhs, i, j), wt.lexical[r]));
          }
      for (std::size_t r = 0; r < b.unary.size(); r++)
        if (b.unary[r].origin >= 0) {
          auto& c = count[static_cast<std::size_t>(b.unary[r].origin)];
          c = S::plus(c, S::times(out.at(b.unary[r].lhs, i, j), S::times(wt.unary[r], in.at(b.unary[r].child, i, j))));
        }
      for (std::size_t r = 0; r < b.binary.size(); r++) {
        const auto& x = b.binary[r];
        if (x.origin < 0) continue;
        auto& c = count[static_cast<std::size_t>(x.origin)];
        for (std::size_t k = i + 1; k < j; k++)
          c = S::plus(c, S::times(out.at(x.lhs, i, j), S::times(wt.binary[r], S::times(in.at(x.left, i, k), in.at(x.right, k, j)))));
      }
    }
  for (auto& c : count) c = c / z;  // a field: the semiring of the reals (or rationals)
  return count;
}

// The best parse under a selective semiring (MaxPlus: the Viterbi parse), as a tree of the
// original grammar's rules over spans of w; its weight is the root's inside value. Throws when w
// has no parse.
template <SelectiveSemiring S>
ParseTree best_parse([[maybe_unused]] const Grammar& g, const BinaryGrammar& b, const std::vector<typename S::value_type>& weights,
                     const std::vector<Symbol>& w, typename S::value_type* weight = nullptr) {
  const std::size_t n = w.size();
  const BinaryWeights<S> wt(b, weights);
  Chart<S> c(b, n);
  // back[A][i][j]: kind (0 none, 1 lexical, 2 unary, 3 binary), rule, split.
  struct Back { std::uint8_t kind = 0; std::uint32_t rule = 0, split = 0; };
  std::vector<Back> back(b.nonterminals * (n + 1) * (n + 1));
  const auto bk = [&](std::size_t a, std::size_t i, std::size_t j) -> Back& { return back[(a * (n + 1) + i) * (n + 1) + j]; };
  const auto offer = [&](std::size_t a, std::size_t i, std::size_t j, const typename S::value_type& v, Back how) {
    auto& cur = c.at(a, i, j);
    if (S::equal(v, S::zero())) return;
    if (bk(a, i, j).kind == 0 || !S::equal(S::plus(cur, v), cur)) { cur = S::plus(cur, v); bk(a, i, j) = how; }
  };
  for (std::size_t len = 1; len <= n; len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      if (len == 1)
        for (std::size_t r = 0; r < b.lexical.size(); r++)
          if (b.lexical[r].terminal == w[i]) offer(b.lexical[r].lhs, i, j, wt.lexical[r], {1, static_cast<std::uint32_t>(r), 0});
      for (std::size_t r = 0; r < b.binary.size(); r++) {
        const auto& x = b.binary[r];
        for (std::size_t k = i + 1; k < j; k++) {
          if (bk(x.left, i, k).kind == 0 || bk(x.right, k, j).kind == 0) continue;
          offer(x.lhs, i, j, S::times(wt.binary[r], S::times(c.at(x.left, i, k), c.at(x.right, k, j))), {3, static_cast<std::uint32_t>(r), static_cast<std::uint32_t>(k)});
        }
      }
      for (std::uint32_t a : b.unary_order)
        for (std::size_t r = 0; r < b.unary.size(); r++)
          if (b.unary[r].lhs == a && bk(b.unary[r].child, i, j).kind != 0)
            offer(a, i, j, S::times(wt.unary[r], c.at(b.unary[r].child, i, j)), {2, static_cast<std::uint32_t>(r), 0});
    }
  if (n == 0 || bk(0, 0, n).kind == 0) throw std::invalid_argument("morphe::formal::best_parse: the word has no parse");
  if (weight) *weight = c.at(0, 0, n);
  // Back to the original rules: a helper's children continue its parent's right side, and a
  // preterminal's subtree is a terminal leaf.
  std::function<void(std::size_t, std::size_t, std::size_t, std::vector<ParseTree>&)> children;
  std::function<ParseTree(std::size_t, std::size_t, std::size_t)> node = [&](std::size_t a, std::size_t i, std::size_t j) -> ParseTree {
    if (b.preterminal_of[a] >= 0) { ParseTree leaf; leaf.begin = i; leaf.end = j; return leaf; }
    const Back& h = bk(a, i, j);
    ParseTree t;
    t.begin = i;
    t.end = j;
    if (h.kind == 1) {
      t.rule = b.lexical[h.rule].origin;
      ParseTree leaf; leaf.begin = i; leaf.end = j;
      t.children.push_back(leaf);
    } else if (h.kind == 2) {
      t.rule = b.unary[h.rule].origin;
      t.children.push_back(node(b.unary[h.rule].child, i, j));
    } else {
      t.rule = b.binary[h.rule].origin;
      children(a, i, j, t.children);
    }
    return t;
  };
  children = [&](std::size_t a, std::size_t i, std::size_t j, std::vector<ParseTree>& out) {
    const Back& h = bk(a, i, j);
    const auto& x = b.binary[h.rule];
    out.push_back(node(x.left, i, h.split));
    if (b.helper[x.right]) children(x.right, h.split, j, out);
    else out.push_back(node(x.right, h.split, j));
  };
  return node(0, 0, n);
}

}  // namespace morphe::formal
