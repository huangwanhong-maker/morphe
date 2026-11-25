// Checks of straight-line grammars by Re-Pair (repair.hpp, decision 0022).
//
// - abababab: A → ab (four times), B → AA (twice), the words BB, size 6 against 8; expanded, the word given.
// - A run of five equal letters: aa occurs twice without overlap, so A → aa and the word AAa; then AA occurs once and the loop
//   stops; size 5.
// - Pairs never span words: two words "a" and "b" give no rule; "ab" twice gives one.
// - Ties go to the least pair: ab and cd each twice, ab first.
// - On 300 random sets of words over two to four letters: every word recovered exactly by expansion; no pair occurring twice in
//   the words left (digram uniqueness); each rule's right side made of earlier letters; the grammar no larger than its input;
//   a bound on the rules kept.

#include <morphe/check.hpp>
#include <morphe/formal/repair.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

namespace mf = morphe::formal;
using W = std::vector<std::uint32_t>;

int main() {
  morphe::Suite s;
  {
    const auto G = mf::re_pair({W{0, 1, 0, 1, 0, 1, 0, 1}});
    s.require("abababab: A → ab, B → AA, the words BB, size 6 against 8",
              G.rules.size() == 2 && G.rules[0] == std::make_pair(0u, 1u) && G.rules[1] == std::make_pair(G.first_new, G.first_new) &&
                  G.words == std::vector<W>{W{G.first_new + 1, G.first_new + 1}} && G.size() == 6 && G.input_length == 8);
    s.require("... expanded, the word given", G.expand() == std::vector<W>{W{0, 1, 0, 1, 0, 1, 0, 1}});
    s.require("... B stands for four letters", G.length(G.first_new + 1) == 4);
  }
  {
    const auto G = mf::re_pair({W{5, 5, 5, 5, 5}});
    s.require("a run of five equal letters: aa twice without overlap, A → aa and the word AAa; then AA once, so the loop stops; size 5",
              G.rules.size() == 1 && G.words == std::vector<W>{W{G.first_new, G.first_new, 5}} && G.size() == 5);
  }
  {
    const auto apart = mf::re_pair({W{0}, W{1}}), twice = mf::re_pair({W{0, 1}, W{0, 1}});
    s.require("pairs never span words: \"a\" and \"b\" give no rule; \"ab\" twice gives one", apart.rules.empty() && twice.rules.size() == 1);
    const auto tie = mf::re_pair({W{2, 3}, W{0, 1}, W{2, 3}, W{0, 1}});
    s.require("ties go to the least pair: ab and cd twice each, ab first", tie.rules.size() == 2 && tie.rules[0] == std::make_pair(0u, 1u));
  }
  {
    std::mt19937 g(20261126);
    bool exact = true, unique = true, ordered = true, smaller = true;
    for (int t = 0; t < 300; t++) {
      const std::uint32_t k = 2 + g() % 3;
      std::vector<W> words(1 + g() % 6);
      for (auto& w : words) {
        w.resize(g() % 40);
        for (auto& a : w) a = g() % k;
      }
      const auto G = mf::re_pair(words);
      exact = exact && G.expand() == words;
      // Digram uniqueness: no pair occurs twice without overlap in the words left.
      const auto counts = mf::detail::pair_counts(G.words);
      for (const auto& [key, n] : counts) unique = unique && n < 2;
      for (std::size_t i = 0; i < G.rules.size(); i++)
        ordered = ordered && G.rules[i].first < G.first_new + i && G.rules[i].second < G.first_new + i;
      smaller = smaller && G.size() <= G.input_length;
    }
    s.require("300 random sets of words: every word recovered exactly by expansion", exact);
    s.require("... no pair occurring twice in the words left (digram uniqueness)", unique);
    s.require("... each rule's right side made of earlier letters", ordered);
    s.require("... the grammar no larger than its input", smaller);
    const auto bounded = mf::re_pair({W{0, 1, 0, 1, 0, 1, 0, 1}}, 1);
    s.require("a bound on the rules kept: one rule, the words AAAA", bounded.rules.size() == 1 && bounded.words[0].size() == 4 && bounded.expand()[0].size() == 8);
  }
  return s.done();
}
