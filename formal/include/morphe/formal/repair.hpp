// Straight-line grammars of a set of words by Re-Pair (decision 0022): the pair of adjacent letters that occurs most often is
// replaced by a new letter, and its rule kept, until no pair occurs twice.
//
// A straight-line grammar is a context-free grammar in which every nonterminal has one rule and no nonterminal derives itself, so
// that each derives exactly one word; it is a compressed form of what it derives (Lohrey 2012, "Algorithmics on SLP-compressed
// strings: a survey", Groups Complex. Cryptol. 4, 241–299). Re-Pair (Larsson and Moffat 2000, "Off-line dictionary-based
// compression", Proc. IEEE 88(11), 1722–1732) builds one greedily:
//
//   while some pair ab of adjacent letters occurs at least twice, without overlap, in the words:
//     take the most frequent such pair (ties: the least pair in the order of letters),
//     give it a new letter A with the rule A → ab,
//     replace its occurrences left to right, each word on its own (a pair never spans two words).
//
// Occurrences are counted without overlap: in a run of k equal letters a, the pair aa occurs ⌊k/2⌋ times, and the run is replaced
// from its left. When the loop stops no pair occurs twice in the words (digram uniqueness), and every rule's letter occurs at least
// twice in the words and the right sides together (rule utility for Re-Pair's own rules, at the moment each is made). Each word is
// recovered exactly by expanding its letters (`expand`), which the checks hold on every input. The size of the grammar is the length
// of the words left plus two letters for each rule.
//
// Letters are 32-bit; new letters are numbered from the least letter not used by the input. Cost: each round counts the pairs in a
// hash map over the words, O(n) for n letters left, so O(r · n) for r rules, which suits words of some tens of thousands of letters.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace morphe::formal {

struct StraightLine {
  using Letter = std::uint32_t;
  Letter first_new = 0;                                  // the first letter the grammar made
  std::vector<std::pair<Letter, Letter>> rules;          // rules[i] is the right side of letter first_new + i
  std::vector<std::vector<Letter>> words;                // the words with the pairs replaced
  std::size_t input_length = 0;                          // letters in the words given

  bool nonterminal(Letter a) const { return a >= first_new && a - first_new < rules.size(); }
  std::size_t size() const {
    std::size_t s = 2 * rules.size();
    for (const auto& w : words) s += w.size();
    return s;
  }
  // The word a letter derives.
  std::vector<Letter> derives(Letter a) const {
    std::vector<Letter> out, stack{a};
    while (!stack.empty()) {
      const Letter x = stack.back();
      stack.pop_back();
      if (nonterminal(x)) {
        const auto& r = rules[x - first_new];
        stack.push_back(r.second), stack.push_back(r.first);
      } else {
        out.push_back(x);
      }
    }
    return out;
  }
  // The words given, recovered by expanding every letter.
  std::vector<std::vector<Letter>> expand() const {
    std::vector<std::vector<Letter>> out;
    for (const auto& w : words) {
      out.emplace_back();
      for (Letter a : w) {
        const auto d = derives(a);
        out.back().insert(out.back().end(), d.begin(), d.end());
      }
    }
    return out;
  }
  // How many letters a letter stands for.
  std::size_t length(Letter a) const { return derives(a).size(); }
};

namespace detail {
inline std::uint64_t pair_key(std::uint32_t a, std::uint32_t b) { return (static_cast<std::uint64_t>(a) << 32) | b; }

// The occurrences of every pair, without overlap (in a run of equal letters, ⌊k/2⌋ of the pair aa).
inline std::unordered_map<std::uint64_t, std::size_t> pair_counts(const std::vector<std::vector<std::uint32_t>>& words) {
  std::unordered_map<std::uint64_t, std::size_t> c;
  for (const auto& w : words) {
    std::size_t i = 0;
    while (i + 1 < w.size()) {
      if (w[i] == w[i + 1]) {
        std::size_t j = i;  // the run of equal letters from i
        while (j + 1 < w.size() && w[j + 1] == w[i]) j++;
        const std::size_t k = j - i + 1;
        c[pair_key(w[i], w[i])] += k / 2;
        // The pair from the run's last letter to the next letter is counted at the run's end.
        if (j + 1 < w.size()) c[pair_key(w[j], w[j + 1])]++;
        i = j + 1;
      } else {
        c[pair_key(w[i], w[i + 1])]++;
        i++;
      }
    }
  }
  return c;
}
}  // namespace detail

// The Re-Pair grammar of the words; `max_rules` bounds the rounds (0: no bound).
inline StraightLine re_pair(std::vector<std::vector<std::uint32_t>> words, std::size_t max_rules = 0) {
  StraightLine G;
  std::uint32_t top = 0;
  for (const auto& w : words) {
    G.input_length += w.size();
    for (auto a : w) top = std::max(top, a + 1);
  }
  if (top == 0xFFFFFFFFu) throw std::invalid_argument("morphe::formal::re_pair: no letter is left for rules");
  G.first_new = top;
  while (max_rules == 0 || G.rules.size() < max_rules) {
    const auto counts = detail::pair_counts(words);
    std::uint64_t best = 0;
    std::size_t most = 1;
    for (const auto& [k, n] : counts)
      if (n > most || (n == most && n >= 2 && k < best)) best = k, most = n;
    if (most < 2) break;
    const std::uint32_t a = static_cast<std::uint32_t>(best >> 32), b = static_cast<std::uint32_t>(best & 0xFFFFFFFFu);
    if (G.first_new + G.rules.size() == 0xFFFFFFFFu) break;
    const std::uint32_t A = G.first_new + static_cast<std::uint32_t>(G.rules.size());
    G.rules.push_back({a, b});
    for (auto& w : words) {
      std::vector<std::uint32_t> out;
      out.reserve(w.size());
      for (std::size_t i = 0; i < w.size();) {
        if (i + 1 < w.size() && w[i] == a && w[i + 1] == b) {
          out.push_back(A);
          i += 2;
        } else {
          out.push_back(w[i]);
          i++;
        }
      }
      w.swap(out);
    }
  }
  G.words = std::move(words);
  return G;
}

}  // namespace morphe::formal
