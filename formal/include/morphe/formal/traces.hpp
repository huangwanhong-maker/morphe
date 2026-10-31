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
// Cost: causal_order is O(n² · |items|); linearisations counts by dynamic programming over the down-sets
// of the order, exactly with morphe's integers, and refuses more than 40 events, where the down-sets can
// number 2⁴⁰. Checks: tests/check_traces.cpp.

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

// The causal order of a record, its events in the order recorded.
inline CausalOrder causal_order(const std::vector<Footprint>& record) {
  CausalOrder P;
  P.n = record.size();
  P.before.assign(P.n, std::vector<char>(P.n, 0));
  std::vector<std::vector<char>> direct(P.n, std::vector<char>(P.n, 0));
  for (std::size_t j = 0; j < P.n; j++)
    for (std::size_t i = 0; i < j; i++)
      if (dependent(record[i], record[j])) direct[i][j] = 1;
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

}  // namespace morphe::formal::traces
