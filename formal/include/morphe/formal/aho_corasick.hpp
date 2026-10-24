// The Aho–Corasick automaton of a dictionary: every occurrence of every pattern of a finite set
// in one left-to-right pass over a text (Aho and Corasick 1975, "Efficient string matching: an
// aid to bibliographic search", CACM 18, 333–340).
//
// - The goto function g is the trie of the patterns: a state is a prefix of some pattern, the
//   root the empty word (Algorithm 2). The failure function f sends a state u to the state of
//   the longest proper suffix of u that is a state, computed breadth-first: f(ua) = g*(f(u), a)
//   (Algorithm 3). The output of a state is the set of patterns that are suffixes of it: its own
//   pattern, if it is one, and the output of f(u); here each state links to the nearest state on
//   its failure chain that ends a pattern, so the outputs at a position are listed in time
//   proportional to their number.
// - The deterministic automaton: δ(u, a) = g(u, a) where defined, else δ(f(u), a), and
//   δ(root, a) = root where g is undefined (Algorithm 4, the next-move function). It recognises
//   A* P, the words that end with a pattern: a finite automaton for a regular language, with
//   one transition per letter of the text.
//
// Construction O(m σ) time and space for patterns of total length m (the transitions are a
// dense table of the states by the alphabet); a pass over a text of length n is O(n + occ).
// Patterns must be non-empty; equal patterns are each reported.

#pragma once

#include <morphe/formal/words.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace morphe::formal {

class AhoCorasick {
 public:
  static constexpr std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();

  struct Match {
    std::size_t begin = 0, end = 0;  // the occurrence text[begin, end)
    std::size_t pattern = 0;         // its index in the dictionary
    friend bool operator==(const Match&, const Match&) = default;
  };

  AhoCorasick() = default;

  template <LetterRange R>
  AhoCorasick(const std::vector<R>& patterns, std::size_t sigma) : sigma_(sigma) {
    if (sigma == 0) throw std::invalid_argument("morphe::formal::AhoCorasick: an alphabet with letters");
    lengths_.reserve(patterns.size());
    add_state(0);
    // The trie: g.
    for (std::size_t p = 0; p < patterns.size(); p++) {
      const auto s = letters(patterns[p]);
      if (s.empty()) throw std::invalid_argument("morphe::formal::AhoCorasick: non-empty patterns");
      std::uint32_t u = 0;
      for (auto c : s) {
        if (static_cast<std::size_t>(c) >= sigma) throw std::invalid_argument("morphe::formal::AhoCorasick: a letter outside the alphabet");
        std::uint32_t& v = goto_[u * sigma_ + c];
        if (v == NONE) {
          const std::uint32_t w = add_state(depth_[u] + 1);
          goto_[u * sigma_ + c] = w;  // add_state may have moved the table
          u = w;
        } else {
          u = v;
        }
      }
      own_[u].push_back(static_cast<std::uint32_t>(p));
      lengths_.push_back(s.size());
    }
    // f and δ, breadth-first: a state's failure is shallower, so it is complete before it is used.
    const std::size_t states = depth_.size();
    failure_.assign(states, 0);
    link_.assign(states, NONE);
    delta_.assign(states * sigma_, 0);
    std::vector<std::uint32_t> queue;
    queue.reserve(states);
    for (std::size_t a = 0; a < sigma_; a++) {
      const std::uint32_t v = goto_[a];
      if (v != NONE) {
        delta_[a] = v;
        failure_[v] = 0;
        queue.push_back(v);
      }
    }
    for (std::size_t head = 0; head < queue.size(); head++) {
      const std::uint32_t u = queue[head];
      const std::uint32_t fu = failure_[u];
      link_[u] = own_[fu].empty() ? link_[fu] : fu;
      for (std::size_t a = 0; a < sigma_; a++) {
        const std::uint32_t v = goto_[u * sigma_ + a];
        if (v != NONE) {
          failure_[v] = delta_[fu * sigma_ + a];
          delta_[u * sigma_ + a] = v;
          queue.push_back(v);
        } else {
          delta_[u * sigma_ + a] = delta_[fu * sigma_ + a];
        }
      }
    }
  }

  std::size_t states() const { return depth_.size(); }
  std::size_t alphabet() const { return sigma_; }
  std::size_t patterns() const { return lengths_.size(); }

  // g(u, a), or NONE where the trie has no edge.
  std::uint32_t go(std::uint32_t u, Letter a) const { return goto_[u * sigma_ + a]; }
  std::uint32_t failure(std::uint32_t u) const { return failure_[u]; }
  // The transition of the automaton.
  std::uint32_t delta(std::uint32_t u, Letter a) const { return delta_[u * sigma_ + a]; }
  std::size_t depth(std::uint32_t u) const { return depth_[u]; }
  // The patterns that are suffixes of the state's word, longest first.
  std::vector<std::uint32_t> output(std::uint32_t u) const {
    std::vector<std::uint32_t> out;
    for (std::uint32_t t = own_[u].empty() ? link_[u] : u; t != NONE; t = link_[t])
      out.insert(out.end(), own_[t].begin(), own_[t].end());
    return out;
  }

  // Every occurrence, in order of end, the longest first at each end. A letter outside the
  // alphabet matches nothing and returns the automaton to its root.
  template <LetterRange R>
  std::vector<Match> matches(const R& text) const {
    std::vector<Match> out;
    std::uint32_t u = 0;
    const auto s = letters(text);
    for (std::size_t i = 0; i < s.size(); i++) {
      if (static_cast<std::size_t>(s[i]) >= sigma_) { u = 0; continue; }
      u = delta_[u * sigma_ + s[i]];
      for (std::uint32_t t = own_[u].empty() ? link_[u] : u; t != NONE; t = link_[t])
        for (std::uint32_t p : own_[t]) out.push_back({i + 1 - lengths_[p], i + 1, p});
    }
    return out;
  }

  std::size_t bytes() const { return (goto_.size() + delta_.size() + failure_.size() + link_.size() + depth_.size()) * 4; }

 private:
  std::uint32_t add_state(std::size_t depth) {
    const auto id = static_cast<std::uint32_t>(depth_.size());
    depth_.push_back(static_cast<std::uint32_t>(depth));
    own_.emplace_back();
    goto_.resize(goto_.size() + sigma_, NONE);
    return id;
  }

  std::size_t sigma_ = 0;
  std::vector<std::uint32_t> goto_, delta_, failure_, link_, depth_;
  std::vector<std::vector<std::uint32_t>> own_;
  std::vector<std::size_t> lengths_;
};

}  // namespace morphe::formal
