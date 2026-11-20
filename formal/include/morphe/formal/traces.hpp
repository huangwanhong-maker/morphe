// Histories up to the order of independent events: the dependence of events that read and write items,
// the causal order of a record, its Foata normal form, and the number of orderings of the record that one
// trace equates.
//
// Independence (Mazurkiewicz 1977, Concurrent Program Schemes and Their Interpretations, DAIMI PB-78;
// Diekert and Rozenberg 1995, The Book of Traces, ch. 1). Each event e of a record reads a set R(e) of
// items and writes a set W(e) (creating and deleting are writing). Two events are dependent when one writes
// what the other reads or writes:
//
//   e D e′  ⇔  W(e) ∩ (R(e′) ∪ W(e′)) ≠ ∅  or  W(e′) ∩ R(e) ≠ ∅.
//
// Exchanging two adjacent independent events changes neither what each reads nor what the record ends
// with; the records reached by such exchanges form one trace, and they are exactly the linear extensions of
// the causal order, the transitive closure of e_i < e_j for i < j and e_i D e_j (the dependence graph of
// the trace; Diekert and Rozenberg 1995, §1.3). Its Foata normal form (Cartier and Foata 1969, LNM 85)
// cuts the record into steps of pairwise independent events, each event of a step depending on one of the
// step before: an event's step is the length of the longest chain of dependences ending at it.
//
// For double-pushout rewriting (graph_rewriting.hpp), with R(e) every item the match covers and W(e) every
// item the step deletes, creates or relabels, items having identities that persist through the record, two
// adjacent steps are independent in this sense exactly when they are sequentially independent (Ehrig,
// Ehrig, Prange and Taentzer 2006, Fundamentals of Algebraic Graph Transformation, Def. 5.9): an item both
// covered by the second match and left by the first step is then preserved by both. The local
// Church–Rosser theorem (ibid., Thm 5.12) makes the exchange valid, and the trace is the derivation's
// shift-equivalence class (Corradini, Montanari, Rossi, Ehrig, Heckel and Löwe 1997, in the Handbook of
// Graph Grammars, vol. 1, ch. 3).
//
// Words over an alphabet (decision 0021). A trace monoid is given by letters 0…k−1 and an independence I, a
// symmetric and irreflexive relation; two words are equivalent when one becomes the other by exchanging
// adjacent independent letters. A word's events are its positions, dependent when their letters are (equal
// letters are always dependent), so the causal order and its counts serve words as they serve records. Each
// trace has two normal forms:
// - the lexicographic normal form, the least word of the trace in the order of the letters, found by taking
//   at each step the least letter among the events all of whose predecessors are taken (Anisimov and Knuth
//   1979, "Inhomogeneous sorting", Int. J. Comput. Inform. Sci. 8, 255–260; Diekert and Rozenberg 1995,
//   ch. 1); no two such events share a letter, since equal letters are dependent;
// - the Foata normal form, the steps of the causal order, each written in increasing order of its letters.
// Two words are equivalent exactly when their normal forms are equal. Apart from both, `projections_agree`
// decides equivalence by the projection lemma: u ≡ v exactly when, for every dependent pair of letters a, b
// (a = b included), the words keep the same subword in a and b (Cori and Perrin 1985, "Automates et
// commutations partielles", RAIRO Inform. Théor. 19, 21–32; Diekert and Rozenberg 1995, ch. 1). The checks
// hold all three to the classes found by exchanging adjacent letters.
//
// Cost: causal_order is O(n² · |items|), and O(n²) on a word; linearisations counts by dynamic programming
// over the down-sets of the order, exactly with morphe's integers, and refuses more than 40 events, where
// the down-sets can number 2⁴⁰. Checks: tests/check_traces.cpp.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <morphe/arithmetic/integer.hpp>

namespace morphe::formal::traces {

struct Footprint {
  std::set<std::string> reads, writes;
};

namespace detail {
inline bool meets(const std::set<std::string>& a, const std::set<std::string>& b) {
  if (a.size() > b.size()) return meets(b, a);
  for (const auto& x : a)
    if (b.count(x)) return true;
  return false;
}
}  // namespace detail

inline bool dependent(const Footprint& a, const Footprint& b) {
  return detail::meets(a.writes, b.reads) || detail::meets(a.writes, b.writes) || detail::meets(b.writes, a.reads);
}

struct CausalOrder {
  std::size_t n = 0;
  std::vector<std::vector<char>> before;                    // before[i][j]: e_i < e_j
  std::vector<std::pair<std::size_t, std::size_t>> covers;  // the Hasse diagram: e_i ⋖ e_j
  std::vector<std::vector<std::size_t>> steps;              // the Foata normal form
  std::vector<std::size_t> step;                            // the step of each event (0-based)

  bool concurrent(std::size_t i, std::size_t j) const { return i != j && !before[i][j] && !before[j][i]; }
  // Whether an ordering of the events (a permutation) is a linear extension: one of the trace's records.
  bool respects(const std::vector<std::size_t>& order) const {
    if (order.size() != n) return false;
    std::vector<std::size_t> pos(n, n);
    for (std::size_t k = 0; k < n; k++) {
      if (order[k] >= n || pos[order[k]] != n) return false;
      pos[order[k]] = k;
    }
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++)
        if (before[i][j] && pos[i] > pos[j]) return false;
    return true;
  }
};

namespace detail {
// The causal order of n events in their recorded order, given which earlier event directly precedes which.
template <class Dependent>
CausalOrder order_from(std::size_t n, Dependent&& dep) {
  CausalOrder P;
  P.n = n;
  P.before.assign(P.n, std::vector<char>(P.n, 0));
  std::vector<std::vector<char>> direct(P.n, std::vector<char>(P.n, 0));
  for (std::size_t j = 0; j < P.n; j++)
    for (std::size_t i = 0; i < j; i++)
      if (dep(i, j)) direct[i][j] = 1;
  // Transitive closure in the recorded order (a topological order): before[i][j] if some direct i → k ≤ j.
  for (std::size_t j = 0; j < P.n; j++)
    for (std::size_t i = 0; i < j; i++) {
      if (direct[i][j]) { P.before[i][j] = 1; continue; }
      for (std::size_t k = i + 1; k < j; k++)
        if (P.before[i][k] && direct[k][j]) { P.before[i][j] = 1; break; }
    }
  // Covers: i < j with no k between.
  for (std::size_t i = 0; i < P.n; i++)
    for (std::size_t j = i + 1; j < P.n; j++) {
      if (!P.before[i][j]) continue;
      bool cover = true;
      for (std::size_t k = i + 1; k < j && cover; k++)
        if (P.before[i][k] && P.before[k][j]) cover = false;
      if (cover) P.covers.push_back({i, j});
    }
  P.step.assign(P.n, 0);
  for (std::size_t j = 0; j < P.n; j++)
    for (std::size_t i = 0; i < j; i++)
      if (P.before[i][j]) P.step[j] = std::max(P.step[j], P.step[i] + 1);
  for (std::size_t j = 0; j < P.n; j++) {
    if (P.steps.size() <= P.step[j]) P.steps.resize(P.step[j] + 1);
    P.steps[P.step[j]].push_back(j);
  }
  return P;
}
}  // namespace detail

// The causal order of a record, its events in the order recorded.
inline CausalOrder causal_order(const std::vector<Footprint>& record) {
  return detail::order_from(record.size(), [&](std::size_t i, std::size_t j) { return dependent(record[i], record[j]); });
}

// The number of linear extensions of the causal order: the records in the trace. At most 40 events.
inline arithmetic::Integer linearisations(const CausalOrder& P) {
  if (P.n > 40) throw std::invalid_argument("morphe::formal::traces::linearisations: more than 40 events");
  std::vector<std::uint64_t> pred(P.n, 0);
  for (std::size_t i = 0; i < P.n; i++)
    for (std::size_t j = 0; j < P.n; j++)
      if (P.before[i][j]) pred[j] |= std::uint64_t{1} << i;
  std::unordered_map<std::uint64_t, arithmetic::Integer> ways{{0, arithmetic::Integer(1)}};
  for (std::size_t size = 0; size < P.n; size++) {
    std::unordered_map<std::uint64_t, arithmetic::Integer> next;
    for (const auto& [S, w] : ways)
      for (std::size_t x = 0; x < P.n; x++) {
        const std::uint64_t bit = std::uint64_t{1} << x;
        if ((S & bit) || (pred[x] & ~S)) continue;
        auto it = next.find(S | bit);
        if (it == next.end()) next.emplace(S | bit, w);
        else it->second = it->second + w;
      }
    ways = std::move(next);
  }
  return ways.empty() ? arithmetic::Integer(0) : ways.begin()->second;
}

// ---- words over an alphabet with an independence ----

// Letters 0…k−1; independent[a][b] when a and b commute. Symmetric and irreflexive, checked when built and when used.
struct Independence {
  std::size_t letters = 0;
  std::vector<std::vector<char>> independent;
  bool commute(std::size_t a, std::size_t b) const { return independent[a][b] != 0; }
};

// From the pairs that commute; refuses a letter out of range or a letter said to commute with itself.
inline Independence independence(std::size_t letters, const std::vector<std::pair<std::size_t, std::size_t>>& commuting) {
  Independence I;
  I.letters = letters;
  I.independent.assign(letters, std::vector<char>(letters, 0));
  for (const auto& [a, b] : commuting) {
    if (a >= letters || b >= letters) throw std::invalid_argument("morphe::formal::traces::independence: a letter out of range");
    if (a == b) throw std::invalid_argument("morphe::formal::traces::independence: a letter cannot commute with itself");
    I.independent[a][b] = I.independent[b][a] = 1;
  }
  return I;
}

// From each letter's footprint: two distinct letters commute when their footprints are independent.
inline Independence independence(const std::vector<Footprint>& letters) {
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t a = 0; a < letters.size(); a++)
    for (std::size_t b = a + 1; b < letters.size(); b++)
      if (!dependent(letters[a], letters[b])) pairs.push_back({a, b});
  return independence(letters.size(), pairs);
}

namespace detail {
inline void check_word(const std::vector<std::size_t>& w, const Independence& I, const char* who) {
  const std::string where = std::string("morphe::formal::traces::") + who;
  if (I.independent.size() != I.letters) throw std::invalid_argument(where + ": a malformed independence");
  for (std::size_t a = 0; a < I.letters; a++) {
    if (I.independent[a].size() != I.letters || I.independent[a][a]) throw std::invalid_argument(where + ": an independence that is not irreflexive");
    for (std::size_t b = 0; b < a; b++)
      if ((I.independent[a][b] != 0) != (I.independent[b][a] != 0)) throw std::invalid_argument(where + ": an independence that is not symmetric");
  }
  for (std::size_t x : w)
    if (x >= I.letters) throw std::invalid_argument(where + ": a letter out of range");
}
}  // namespace detail

// The causal order of a word: positions i < j ordered when their letters are dependent, closed transitively.
inline CausalOrder causal_order(const std::vector<std::size_t>& word, const Independence& I) {
  detail::check_word(word, I, "causal_order");
  return detail::order_from(word.size(), [&](std::size_t i, std::size_t j) { return !I.commute(word[i], word[j]); });
}

// The least word of the trace in the order of the letters.
inline std::vector<std::size_t> lexicographic_normal_form(const std::vector<std::size_t>& word, const Independence& I) {
  const CausalOrder P = causal_order(word, I);
  std::vector<std::size_t> waiting(P.n, 0), out;
  for (std::size_t i = 0; i < P.n; i++)
    for (std::size_t j = 0; j < P.n; j++) waiting[j] += P.before[i][j] ? 1 : 0;
  std::vector<char> taken(P.n, 0);
  for (std::size_t k = 0; k < P.n; k++) {
    std::size_t best = P.n;
    for (std::size_t i = 0; i < P.n; i++)
      if (!taken[i] && waiting[i] == 0 && (best == P.n || word[i] < word[best])) best = i;
    taken[best] = 1;
    out.push_back(word[best]);
    for (std::size_t j = 0; j < P.n; j++)
      if (P.before[best][j]) waiting[j]--;
  }
  return out;
}

// The Foata normal form: the causal order's steps, each as its letters in increasing order.
inline std::vector<std::vector<std::size_t>> foata_normal_form(const std::vector<std::size_t>& word, const Independence& I) {
  const CausalOrder P = causal_order(word, I);
  std::vector<std::vector<std::size_t>> steps;
  for (const auto& st : P.steps) {
    std::vector<std::size_t> letters;
    for (std::size_t i : st) letters.push_back(word[i]);
    std::sort(letters.begin(), letters.end());
    steps.push_back(std::move(letters));
  }
  return steps;
}

// Whether two words are one trace: their lexicographic normal forms are equal.
inline bool equivalent(const std::vector<std::size_t>& u, const std::vector<std::size_t>& v, const Independence& I) {
  if (u.size() != v.size()) {
    detail::check_word(u, I, "equivalent");
    detail::check_word(v, I, "equivalent");
    return false;
  }
  return lexicographic_normal_form(u, I) == lexicographic_normal_form(v, I);
}

// The projection lemma, apart from the normal forms: for every dependent pair a, b (a = b included) the
// subwords in a and b agree.
inline bool projections_agree(const std::vector<std::size_t>& u, const std::vector<std::size_t>& v, const Independence& I) {
  detail::check_word(u, I, "projections_agree");
  detail::check_word(v, I, "projections_agree");
  for (std::size_t a = 0; a < I.letters; a++)
    for (std::size_t b = a; b < I.letters; b++) {
      if (I.commute(a, b)) continue;
      std::vector<std::size_t> pu, pv;
      for (std::size_t x : u)
        if (x == a || x == b) pu.push_back(x);
      for (std::size_t x : v)
        if (x == a || x == b) pv.push_back(x);
      if (pu != pv) return false;
    }
  return true;
}

}  // namespace morphe::formal::traces
