// Grammars: phrase-structure grammars over a finite alphabet, the context-free and the monotone
// (context-sensitive) among them, and their derivations.
//
// A grammar Γ = (V, Σ, R, S) has terminals Σ = {0, …, t − 1}, nonterminals V, numbered after
// them (symbol t + i is the nonterminal i), the start S = nonterminal 0, and rules α → β with
// α, β ∈ (V ∪ Σ)*, α containing a nonterminal (Chomsky 1959, On certain formal properties of
// grammars, Information and Control 2, 137–167; Hopcroft and Ullman 1979, Introduction to
// Automata Theory, Languages, and Computation, §9.1). It is
//
//   context-free     when every α is one nonterminal;
//   monotone         when |α| ≤ |β| for every rule, except S → ε when S occurs in no right side.
//                    Monotone grammars generate exactly the context-sensitive languages
//                    (Hopcroft and Ullman, Theorem 9.5), whose membership problem is decidable
//                    (a derivation of w never passes a form longer than |w|) and PSPACE-complete
//                    in general (Karp 1972, the linear-bounded automaton acceptance problem).
//
// derives_monotone decides membership by searching the sentential forms of length ≤ |w|
// breadth-first from S; it terminates on every monotone grammar and costs, in the worst case,
// the number of forms, Σ_{ℓ≤n} |V ∪ Σ|^ℓ. An optional canonicaliser maps each form to a
// representative of its class, so that the search visits classes instead of forms: when Γ is
// invariant under a group acting on its symbols and fixing S, the set of forms derivable from S
// is invariant, and w is derivable exactly when its class is (formal/equivariant.hpp). Context-
// free grammars are parsed by formal/earley.hpp and weighted by formal/inside_outside.hpp.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace morphe::formal {

using Symbol = std::uint32_t;
using Form = std::vector<Symbol>;

struct Rule {
  Form lhs, rhs;
};

struct Grammar {
  std::size_t terminals = 0;              // |Σ|; terminals are 0..terminals−1
  std::vector<std::string> nonterminals;  // the names; the start is the first
  std::vector<Rule> rules;

  std::size_t symbols() const { return terminals + nonterminals.size(); }
  Symbol start() const { return static_cast<Symbol>(terminals); }
  Symbol nonterminal(std::size_t i) const { return static_cast<Symbol>(terminals + i); }
  bool terminal(Symbol s) const { return s < terminals; }
  std::size_t index(Symbol s) const { return s - terminals; }  // of a nonterminal

  // Adds a nonterminal and returns its symbol.
  Symbol add_nonterminal(std::string name) {
    nonterminals.push_back(std::move(name));
    return nonterminal(nonterminals.size() - 1);
  }
  void add(Form lhs, Form rhs) { rules.push_back({std::move(lhs), std::move(rhs)}); }
  void add(Symbol lhs, Form rhs) { rules.push_back({Form{lhs}, std::move(rhs)}); }

  // Throws std::invalid_argument unless every symbol is known and every left side holds a
  // nonterminal.
  void validate() const {
    if (nonterminals.empty()) throw std::invalid_argument("morphe::formal::Grammar: there is no start symbol");
    for (const auto& r : rules) {
      bool has_nonterminal = false;
      for (Symbol s : r.lhs) {
        if (s >= symbols()) throw std::invalid_argument("morphe::formal::Grammar: an unknown symbol");
        has_nonterminal |= !terminal(s);
      }
      for (Symbol s : r.rhs)
        if (s >= symbols()) throw std::invalid_argument("morphe::formal::Grammar: an unknown symbol");
      if (!has_nonterminal) throw std::invalid_argument("morphe::formal::Grammar: a left side without a nonterminal");
    }
  }
  bool context_free() const {
    for (const auto& r : rules)
      if (r.lhs.size() != 1 || terminal(r.lhs[0])) return false;
    return true;
  }
  bool monotone() const {
    bool start_on_right = false;
    for (const auto& r : rules)
      for (Symbol s : r.rhs) start_on_right |= s == start();
    for (const auto& r : rules) {
      if (r.lhs.size() <= r.rhs.size()) continue;
      if (r.rhs.empty() && r.lhs == Form{start()} && !start_on_right) continue;
      return false;
    }
    return true;
  }
  // The nonterminals that derive the empty word (context-free grammars), by the least fixpoint.
  std::vector<char> nullable() const {
    std::vector<char> n(nonterminals.size(), 0);
    for (bool changed = true; changed;) {
      changed = false;
      for (const auto& r : rules) {
        if (r.lhs.size() != 1 || n[index(r.lhs[0])]) continue;
        bool all = true;
        for (Symbol s : r.rhs) all = all && !terminal(s) && n[index(s)];
        if (all) { n[index(r.lhs[0])] = 1; changed = true; }
      }
    }
    return n;
  }
};

// The forms obtained from u by one application of one rule, at every position.
inline std::vector<Form> one_step(const Grammar& g, const Form& u, std::size_t max_length = SIZE_MAX) {
  std::vector<Form> out;
  for (const auto& r : g.rules) {
    if (r.lhs.size() > u.size()) continue;
    const std::size_t grown = u.size() - r.lhs.size() + r.rhs.size();
    if (grown > max_length) continue;
    for (std::size_t i = 0; i + r.lhs.size() <= u.size(); i++) {
      if (!std::equal(r.lhs.begin(), r.lhs.end(), u.begin() + static_cast<std::ptrdiff_t>(i))) continue;
      Form v(u.begin(), u.begin() + static_cast<std::ptrdiff_t>(i));
      v.insert(v.end(), r.rhs.begin(), r.rhs.end());
      v.insert(v.end(), u.begin() + static_cast<std::ptrdiff_t>(i + r.lhs.size()), u.end());
      out.push_back(std::move(v));
    }
  }
  return out;
}

struct MonotoneSearch {
  bool derives = false;
  std::size_t visited = 0;  // forms, or classes when a canonicaliser is given
};

// Whether a monotone grammar derives w: breadth-first search of the forms of length ≤ |w| from S,
// each replaced by canon(form) when canon is given. Exponential in |w| in the worst case
// (PSPACE-complete in general); exact.
inline MonotoneSearch derives_monotone(const Grammar& g, const Form& w, const std::function<Form(const Form&)>& canon = {}) {
  if (!g.monotone()) throw std::invalid_argument("morphe::formal::derives_monotone: the grammar is not monotone");
  const auto key = [&](const Form& f) { return canon ? canon(f) : f; };
  const Form target = key(w);
  std::set<Form> seen;
  std::deque<Form> queue;
  const Form s{g.start()};
  seen.insert(key(s));
  queue.push_back(key(s));
  MonotoneSearch out;
  while (!queue.empty()) {
    Form u = std::move(queue.front());
    queue.pop_front();
    if (u == target) { out.derives = true; break; }
    for (Form& v : one_step(g, u, w.size())) {
      Form k = key(v);
      if (seen.insert(k).second) queue.push_back(std::move(k));
    }
  }
  out.visited = seen.size();
  return out;
}

}  // namespace morphe::formal
