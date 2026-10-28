// Grammars with symmetry: a finite group acting on a grammar's symbols, the rules taken up to
// that action, and parsers that use it.
//
// A group G acts on a grammar Γ = (V, Σ, R, S) when it acts on Σ (an algebra::GroupAction, whose
// sign lets an element reverse words as well) and on V, fixing S; an element g carries a rule
// α → β to g·α → g·β, the forms mapped symbol by symbol and reversed when ε(g) = −1. Γ is
// G-invariant (equivariant) when g·R = R for every g. These are the group-equivariant grammars
// of Huang's draft "Towards a group context-sensitive grammar for genome analysis" (2026),
// stated here in the form in which their claims can be proved; the domain's theory note
// (domains/genomics/theory/group-grammars.md) gives the proofs and says where the draft needs
// amending.
//
// Results used below (proofs in the note):
//
// (1) If Γ is G-invariant, S ⇒* u implies S ⇒* g·u, so L(Γ) is G-invariant, and for every
//     nonterminal A fixed by G, A ⇒* u ⟺ A ⇒* g·u.
// (2) The rules fall into |R/G| = (1/|G|) Σ_g |Fix_R(g)| ≥ |R|/|G| orbits (Burnside's lemma):
//     storing one rule per orbit compresses R by a factor of at most |G|, reached only when no
//     element fixes a rule.
// (3) Scanning a terminal a of a rule against an input letter x by "x ∈ G·a" alone, position by
//     position, decides membership in π⁻¹(L(Γ/G)), the words whose image under the orbit map
//     π: Σ → Σ/G lies in the image grammar's language. This contains L(Γ) and in general
//     strictly: S → aa | bb with a ↔ b accepts ab. An equivariant parser must instead bind the
//     group element once per rule: the items below carry the set C ⊆ G of elements under which
//     the rule's scanned terminals match, C narrowing at each scan to {g ∈ C : g·a = x}.
// (4) With that binding, the equivariant Earley chart is the plain chart taken up to the action:
//     an item [r̄ → α • β, i, C] stands for the plain items [g·r̄ → α • β, i] for g ∈ C, and the
//     parser accepts w exactly when Earley's parser of the full rule set does. Prediction
//     creates one item per orbit instead of one per rule; the chart is never larger than the
//     plain one, and as large in the worst case, so the worst-case time stays O(n³ |R|).
// (5) The orbit-canonical search of a monotone grammar's sentential forms (grammar.hpp) visits
//     the classes of the derivable forms, between |forms|/|G| and |forms| of them.
//
// The equivariant Earley parser needs an action by letters alone (no reversing element) that
// fixes every nonterminal: then the elements bound in a rule are independent of those of its
// children. Reversing elements are used by the canonical search, by rule orbits, and by the
// content-keyed span cache of span_recognition below.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <morphe/algebra/action.hpp>

#include "earley.hpp"
#include "grammar.hpp"
#include "inside_outside.hpp"
#include "semiring.hpp"

namespace morphe::formal {

// A group acting on a grammar: on its terminals (with signs) and on its nonterminals, which it
// must permute fixing the start; an empty nonterminal table fixes them all.
struct GrammarAction {
  const algebra::GroupAction* letters = nullptr;
  std::vector<std::uint32_t> nonterminal_table;  // [g * |V| + i] = g·i, or empty

  std::size_t order() const { return letters->order(); }
  bool fixes_nonterminals() const {
    if (nonterminal_table.empty()) return true;
    const std::size_t v = nonterminal_table.size() / order();
    for (std::uint32_t g = 0; g < order(); g++)
      for (std::uint32_t i = 0; i < v; i++)
        if (nonterminal_table[g * v + i] != i) return false;
    return true;
  }
  Symbol act(const Grammar& gr, std::uint32_t g, Symbol s) const {
    if (gr.terminal(s)) return letters->act(g, s);
    if (nonterminal_table.empty()) return s;
    const std::size_t v = gr.nonterminals.size();
    return gr.nonterminal(nonterminal_table[g * v + gr.index(s)]);
  }
  Form act(const Grammar& gr, std::uint32_t g, const Form& f) const {
    Form out(f.size());
    const std::size_t n = f.size();
    for (std::size_t i = 0; i < n; i++) out[i] = act(gr, g, f[letters->reverses(g) ? n - 1 - i : i]);
    return out;
  }
  Rule act(const Grammar& gr, std::uint32_t g, const Rule& r) const { return {act(gr, g, r.lhs), act(gr, g, r.rhs)}; }
  // The least image of a form, lexicographically: its class's representative.
  Form canonical(const Grammar& gr, const Form& f) const {
    Form best = f;
    for (std::uint32_t g = 1; g < order(); g++) {
      Form image = act(gr, g, f);
      if (image < best) best = std::move(image);
    }
    return best;
  }
  void validate(const Grammar& gr) const {
    if (!letters || letters->degree != gr.terminals) throw std::invalid_argument("morphe::formal::GrammarAction: the action is not on the grammar's terminals");
    if (!nonterminal_table.empty()) {
      if (nonterminal_table.size() != order() * gr.nonterminals.size()) throw std::invalid_argument("morphe::formal::GrammarAction: the nonterminal table has the wrong size");
      for (std::uint32_t g = 0; g < order(); g++)
        if (nonterminal_table[g * gr.nonterminals.size()] != 0) throw std::invalid_argument("morphe::formal::GrammarAction: the start symbol is moved");
    }
  }
};

namespace detail {
inline bool rule_less(const Rule& a, const Rule& b) { return a.lhs != b.lhs ? a.lhs < b.lhs : a.rhs < b.rhs; }
struct RuleLess { bool operator()(const Rule& a, const Rule& b) const { return rule_less(a, b); } };
}  // namespace detail

// Whether g·R = R for every g.
inline bool invariant(const Grammar& gr, const GrammarAction& a) {
  a.validate(gr);
  std::set<Rule, detail::RuleLess> rules(gr.rules.begin(), gr.rules.end());
  for (std::uint32_t g = 0; g < a.order(); g++)
    for (const auto& r : gr.rules)
      if (!rules.count(a.act(gr, g, r))) return false;
  return true;
}

// The grammar whose rules are the images of the given ones under every element (the smallest
// invariant grammar containing them), duplicates removed.
inline Grammar closure(const Grammar& gr, const GrammarAction& a) {
  a.validate(gr);
  Grammar out = gr;
  out.rules.clear();
  std::set<Rule, detail::RuleLess> seen;
  for (const auto& r : gr.rules)
    for (std::uint32_t g = 0; g < a.order(); g++) {
      Rule image = a.act(gr, g, r);
      if (seen.insert(image).second) out.rules.push_back(std::move(image));
    }
  return out;
}

// The orbits of the rules of an invariant grammar: one representative each (the first in R),
// the orbit of every rule, its size; and Burnside's count, which must equal their number.
struct RuleOrbits {
  std::vector<std::uint32_t> representatives;
  std::vector<std::uint32_t> orbit_of;   // per rule
  std::vector<std::size_t> size;         // per orbit
  std::size_t burnside = 0;              // (1/|G|) Σ_g |Fix_R(g)|
};

inline RuleOrbits rule_orbits(const Grammar& gr, const GrammarAction& a) {
  if (!invariant(gr, a)) throw std::invalid_argument("morphe::formal::rule_orbits: the grammar is not invariant");
  std::map<Rule, std::uint32_t, detail::RuleLess> where;
  for (std::uint32_t r = 0; r < gr.rules.size(); r++) where.emplace(gr.rules[r], r);
  RuleOrbits out;
  out.orbit_of.assign(gr.rules.size(), UINT32_MAX);
  for (std::uint32_t r = 0; r < gr.rules.size(); r++) {
    if (out.orbit_of[r] != UINT32_MAX) continue;
    const auto id = static_cast<std::uint32_t>(out.representatives.size());
    out.representatives.push_back(r);
    std::size_t size = 0;
    for (std::uint32_t g = 0; g < a.order(); g++) {
      const std::uint32_t s = where.at(a.act(gr, g, gr.rules[r]));
      if (out.orbit_of[s] == UINT32_MAX) { out.orbit_of[s] = id; size++; }
    }
    out.size.push_back(size);
  }
  std::size_t fixed = 0;
  for (std::uint32_t g = 0; g < a.order(); g++)
    for (const auto& r : gr.rules) {
      const Rule image = a.act(gr, g, r);
      fixed += image.lhs == r.lhs && image.rhs == r.rhs;
    }
  out.burnside = fixed / a.order();
  return out;
}

// ---- the equivariant Earley parser ------------------------------------------------------------------

struct EquivariantEarleyResult {
  bool accepts = false;
  std::size_t items = 0;            // equivariant items
  std::size_t plain_equivalent = 0;  // the plain items they stand for, Σ |{g·r̄ : g ∈ C}|
};

// Earley's parser over one representative per rule orbit, its items binding the group element of
// their rule as a set C ⊆ G narrowed at each scan (result (4)). With consistent = false the
// binding is dropped and a terminal is matched by its orbit alone, as in the draft's scan rule:
// the parser then decides π⁻¹(L(Γ/G)) (result (3)), which the checks compare.
class EquivariantEarley {
 public:
  EquivariantEarley(const Grammar& gr, const GrammarAction& a, bool consistent = true)
      : g_(gr), a_(a), consistent_(consistent), orbits_(rule_orbits(gr, a)) {
    if (!gr.context_free()) throw std::invalid_argument("morphe::formal::EquivariantEarley: the grammar is not context-free");
    if (!a.fixes_nonterminals()) throw std::invalid_argument("morphe::formal::EquivariantEarley: the action must fix every nonterminal");
    for (std::uint32_t g = 0; g < a.order(); g++)
      if (a.letters->reverses(g)) throw std::invalid_argument("morphe::formal::EquivariantEarley: the action must not reverse words");
    words_ = (a.order() + 63) / 64;
    by_lhs_.assign(gr.nonterminals.size(), {});
    for (std::uint32_t o = 0; o < orbits_.representatives.size(); o++) {
      const std::uint32_t r = orbits_.representatives[o];
      by_lhs_[gr.index(gr.rules[r].lhs[0])].push_back(r);
    }
    nullable_ = gr.nullable();
    // match[a][x] = {g : g·a = x}
    match_.assign(gr.terminals * gr.terminals, std::vector<std::uint64_t>(words_, 0));
    for (std::uint32_t t = 0; t < gr.terminals; t++)
      for (std::uint32_t g = 0; g < a.order(); g++) {
        const std::uint32_t x = a.letters->act(g, t);
        match_[t * gr.terminals + x][g / 64] |= 1ull << (g % 64);
      }
    std::vector<std::uint64_t> all(words_, 0);
    for (std::uint32_t g = 0; g < a.order(); g++) all[g / 64] |= 1ull << (g % 64);
    all_ = intern(all);
  }

  EquivariantEarleyResult parse(const std::vector<Symbol>& w) {
    struct Item {
      std::uint32_t rule, dot;
      std::size_t origin;
      std::uint32_t set;  // the interned C
      bool operator<(const Item& o) const { return std::tie(rule, dot, origin, set) < std::tie(o.rule, o.dot, o.origin, o.set); }
    };
    const std::size_t n = w.size();
    std::vector<std::vector<Item>> sets(n + 1);
    std::vector<std::set<Item>> seen(n + 1);
    const auto add = [&](std::size_t j, Item it) { if (seen[j].insert(it).second) sets[j].push_back(it); };
    for (std::uint32_t r : by_lhs_[0]) add(0, {r, 0, 0, all_});
    for (std::size_t j = 0; j <= n; j++)
      for (std::size_t k = 0; k < sets[j].size(); k++) {
        const Item it = sets[j][k];
        const auto& rhs = g_.rules[it.rule].rhs;
        if (it.dot < rhs.size()) {
          const Symbol x = rhs[it.dot];
          if (g_.terminal(x)) {
            if (j == n) continue;
            const auto& m = match_[x * g_.terminals + w[j]];
            if (!consistent_) {
              if (any(m)) add(j + 1, {it.rule, it.dot + 1, it.origin, it.set});
              continue;
            }
            std::vector<std::uint64_t> c = sets_[it.set];
            bool nonempty = false;
            for (std::size_t q = 0; q < words_; q++) { c[q] &= m[q]; nonempty |= c[q] != 0; }
            if (nonempty) add(j + 1, {it.rule, it.dot + 1, it.origin, intern(c)});
          } else {
            for (std::uint32_t r : by_lhs_[g_.index(x)]) add(j, {r, 0, j, all_});
            if (nullable_[g_.index(x)]) add(j, {it.rule, it.dot + 1, it.origin, it.set});
          }
        } else {
          const Symbol lhs = g_.rules[it.rule].lhs[0];
          const std::size_t o = it.origin;
          for (std::size_t p = 0; p < sets[o].size(); p++) {
            const Item parent = sets[o][p];
            const auto& prhs = g_.rules[parent.rule].rhs;
            if (parent.dot < prhs.size() && prhs[parent.dot] == lhs) add(j, {parent.rule, parent.dot + 1, parent.origin, parent.set});
          }
        }
      }
    EquivariantEarleyResult out;
    for (std::size_t j = 0; j <= n; j++) {
      out.items += sets[j].size();
      for (const auto& it : sets[j]) out.plain_equivalent += instances(it.rule, it.set);
    }
    for (const auto& it : sets[n])
      if (it.origin == 0 && it.dot == g_.rules[it.rule].rhs.size() && g_.rules[it.rule].lhs[0] == g_.start()) out.accepts = true;
    return out;
  }

  const RuleOrbits& orbits() const { return orbits_; }

 private:
  static bool any(const std::vector<std::uint64_t>& b) {
    for (auto x : b) if (x) return true;
    return false;
  }
  std::uint32_t intern(const std::vector<std::uint64_t>& b) {
    auto it = ids_.find(b);
    if (it != ids_.end()) return it->second;
    const auto id = static_cast<std::uint32_t>(sets_.size());
    sets_.push_back(b);
    ids_.emplace(b, id);
    return id;
  }
  // The distinct rules g·r̄ for g in a set.
  std::size_t instances(std::uint32_t rule, std::uint32_t set) {
    std::set<Rule, detail::RuleLess> distinct;
    const auto& c = sets_[set];
    for (std::uint32_t g = 0; g < a_.order(); g++)
      if (c[g / 64] >> (g % 64) & 1ull) distinct.insert(a_.act(g_, g, g_.rules[rule]));
    return distinct.size();
  }

  const Grammar& g_;
  const GrammarAction& a_;
  bool consistent_;
  RuleOrbits orbits_;
  std::size_t words_ = 1;
  std::vector<std::vector<std::uint32_t>> by_lhs_;
  std::vector<char> nullable_;
  std::vector<std::vector<std::uint64_t>> match_;
  std::vector<std::vector<std::uint64_t>> sets_;
  std::map<std::vector<std::uint64_t>, std::uint32_t> ids_;
  std::uint32_t all_ = 0;
};

// The image of a grammar under the orbit map π: Σ → Σ/G, terminals renamed by their orbits'
// numbers (algebra::GroupAction::orbit_index), duplicate rules removed: the quotient grammar Γ/G,
// whose language is π(L(Γ)) (a homomorphic image; formal/breaking.hpp).
inline Grammar orbit_grammar(const Grammar& gr, const algebra::GroupAction& letters) {
  const auto id = letters.orbit_index();
  Grammar q;
  q.terminals = letters.orbit_count();
  q.nonterminals = gr.nonterminals;
  std::set<Rule, detail::RuleLess> seen;
  const auto map = [&](const Form& f) {
    Form out;
    for (Symbol s : f) out.push_back(gr.terminal(s) ? id[s] : static_cast<Symbol>(s - gr.terminals + q.terminals));
    return out;
  };
  for (const auto& r : gr.rules) {
    Rule image{map(r.lhs), map(r.rhs)};
    if (seen.insert(image).second) q.rules.push_back(std::move(image));
  }
  return q;
}

// ---- spans shared across the orbit --------------------------------------------------------------

struct SpanRecognition {
  bool accepts = false;
  std::size_t spans = 0;     // the spans whose original nonterminals were decided
  std::size_t reused = 0;    // of them, decided by a span of the same class seen before
  std::size_t classes = 0;   // the distinct classes of spans met
};

// Recognition by the inside algorithm over the Boolean semiring, where the values of the original
// nonterminals on a span are looked up by the span's class under G when an equal-class span was
// decided before: by result (1), every original nonterminal must be fixed by G and Γ invariant.
// With the trivial group the class is the span's content, so that repeated substrings are
// decided once; the group widens the sharing to their images (for DNA, their reverse complements).
// Exact; the cost of finding a class is O(|G| ℓ) for a span of length ℓ.
inline SpanRecognition span_recognition(const Grammar& gr, const GrammarAction& a, const std::vector<Symbol>& w) {
  if (!invariant(gr, a)) throw std::invalid_argument("morphe::formal::span_recognition: the grammar is not invariant");
  if (!a.fixes_nonterminals()) throw std::invalid_argument("morphe::formal::span_recognition: the action must fix every nonterminal");
  const BinaryGrammar b = BinaryGrammar::from(gr);
  const std::size_t n = w.size(), V = b.nonterminals, original = b.originals;
  std::vector<char> in(V * (n + 1) * (n + 1), 0);
  const auto at = [&](std::size_t x, std::size_t i, std::size_t j) -> char& { return in[(x * (n + 1) + i) * (n + 1) + j]; };
  std::map<Form, std::vector<char>> known;
  SpanRecognition out;
  for (std::size_t len = 1; len <= n; len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      Form key = a.canonical(gr, Form(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(j)));
      auto hit = known.find(key);
      const bool reuse = hit != known.end();
      out.spans++;
      if (reuse) {
        out.reused++;
        for (std::size_t x = 0; x < original; x++) at(x, i, j) = hit->second[x];
      }
      // Helpers and preterminals always (their values depend on content, not class); the
      // original nonterminals only when the class is new.
      if (len == 1)
        for (const auto& r : b.lexical)
          if (r.terminal == w[i] && (!reuse || r.lhs >= original)) at(r.lhs, i, j) = 1;
      for (const auto& r : b.binary) {
        if (reuse && r.lhs < original) continue;
        if (at(r.lhs, i, j)) continue;
        for (std::size_t k = i + 1; k < j && !at(r.lhs, i, j); k++) at(r.lhs, i, j) = at(r.left, i, k) && at(r.right, k, j);
      }
      if (!reuse)
        for (std::uint32_t x : b.unary_order)
          for (const auto& r : b.unary)
            if (r.lhs == x && !at(x, i, j) && at(r.child, i, j)) at(x, i, j) = 1;
      if (!reuse) {
        std::vector<char> values(original);
        for (std::size_t x = 0; x < original; x++) values[x] = at(x, i, j);
        known.emplace(std::move(key), std::move(values));
      }
    }
  out.classes = known.size();
  out.accepts = n > 0 && at(0, 0, n);
  return out;
}

// The canonicaliser of sentential forms for grammar.hpp's derives_monotone (result (5)).
inline std::function<Form(const Form&)> form_canonicaliser(const Grammar& gr, const GrammarAction& a) {
  if (!invariant(gr, a)) throw std::invalid_argument("morphe::formal::form_canonicaliser: the grammar is not invariant");
  return [&gr, &a](const Form& f) { return a.canonical(gr, f); };
}

}  // namespace morphe::formal
