// The positional Burrows–Wheeler transform (PBWT) of a panel of M words of one length N over a
// small alphabet [0, σ): the words sorted, at every position k, by their prefixes read backwards,
// so that words sharing a long factor ending at k are neighbours; and the matches it finds, within
// the panel and against a new word (Durbin 2014, "Efficient haplotype matching and storage using
// the positional Burrows–Wheeler transform (PBWT)", Bioinformatics 30, 1266–1272).
//
// Definitions. A panel x_0, …, x_{M−1} ∈ [0, σ)^N. For each k ∈ [0, N]:
// - the prefix array a_k is the permutation of [0, M) that sorts the words by their reversed
//   prefixes x_i[k − 1] x_i[k − 2] … x_i[0] lexicographically, ties by the index i. Write
//   y_j = x_{a_k[j]} for the j-th word in this order.
// - the divergence array d_k: d_k[0] = k, and for j ≥ 1, d_k[j] is the least s with y_j[s, k) =
//   y_{j−1}[s, k), the start of the longest common factor of the two ending at k (k when
//   y_j[k − 1] ≠ y_{j−1}[k − 1]).
// A match of x_i and x_j (i ≠ j) is an interval [s, e), s < e, with x_i[s, e) = x_j[s, e); it is
// maximal when s = 0 or x_i[s − 1] ≠ x_j[s − 1], and e = N or x_i[e] ≠ x_j[e].
//
// Two facts carry everything (Durbin 2014), both immediate from the order: (i) for every
// word u, the j with y_j[k − |u|, k) = u are consecutive; (ii) for j < j′ the longest common
// factor of y_j and y_{j′} ending at k starts at max d_k[j + 1 .. j′].
//
// - Construction. The reversed prefixes at k + 1 begin with x[k], and those with equal x[k] keep
//   their order at k: a_{k+1} is the stable partition of a_k by the letters x_{a_k[j]}[k]. The
//   divergence of a word with its new predecessor, the previous word of a_k with the same letter
//   at k, is by (ii) the maximum of d_k between them, or k + 1 when there is none; one pass with
//   σ running maxima gives d_{k+1}. O(M σ) a position, O(N M σ) in all, O(M) words of memory for
//   the sweep (`pbwt_sweep`); the class `Pbwt` keeps every a_k and d_k, (N + 1) M words each, and
//   the column of letters y_j[k] of each position in a wavelet matrix (succinct.hpp), with which
//   the map w_k(j, c) = (letters of column k below c) + rank_c(column k, j) carries an interval
//   of a_k to the interval of its words with letter c at k in a_{k+1} (the w-map).
// - Long matches: every maximal match of length at least L, reported once at its end e. The pairs
//   whose common factor ending at e is at least L long lie in blocks of a_e separated by
//   d_e[j] > e − L (fact ii); a pair of a block ends at e exactly when their letters at e differ,
//   or e = N, and its start is the maximum of d_e between them. Durbin's report, here with, for
//   each letter, a stack of the earlier words of the block grouped by that maximum (merged as d_e
//   grows), so that each start is known when its pair is reported: O(N M σ + output).
// - Set-maximal matches (Durbin's definition): a maximal match [s, e) of x_i with x_j such that
//   no match of x_i with any word contains [s, e) strictly; they are a word's best matches. At k,
//   the longest common factor of y_j with any other word ending at k starts at s_j = min(d_k[j],
//   d_k[j + 1]) (d_k[M] read as k), and the words that attain it are consecutive around j (fact
//   ii). These matches [s_j, k) are set-maximal exactly when s_j < k and, if k < N, none of those
//   words has y_j's letter at k: a strictly larger match would either start earlier, contradicting
//   the minimality of s_j, or extend one of them past k. The test reads d_k only as far as the
//   nearest words with the same letter above and below j, whose distances sum to at most σ M
//   over a position; the matches are then listed: O(N M σ + output).
// - A new word z ∈ [0, σ)^N against the panel (the basis of haplotype copying, in which a new
//   sequence is read as a mosaic of the panel's; Li and Stephens 2003, Genetics 165, 2213–2233).
//   Swept left to right, it keeps the longest factor z[e, k) ending at k that some word of the
//   panel has at the same place, and the interval [f, g) of a_k of all such words. If some of
//   them has z[k] at k, the interval moves by the w-map and e stays. Otherwise the matches
//   [e, k) of [f, g) are set-maximal (as above), and are reported; the new e′ is the least
//   divergence from the interval of a word with z[k] at k, the nearest such above f giving
//   max d_k(i, f] and the nearest below g giving max d_k[g, i] (fact ii), and the words within
//   divergence e′ are found by widening [f, g) while d_k ≤ e′, then moved by the w-map. O(log σ)
//   a position that extends the match; at each break, a scan of a_k outward from the old interval
//   to the nearest words with z[k] and to the widened interval, at most M positions. Every
//   set-maximal match of z is reported once, with each word that attains it.
//
// Any alphabet is allowed, σ entering the costs linearly; Durbin writes the binary case.
// Everything is exact and terminates on every input; indices are 32-bit by default (M and N
// below 2^32 − 1).

#pragma once

#include <morphe/formal/succinct.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace morphe::formal {

// A match of the words first and second on [begin, end). For set-maximal matches it is first's;
// for long matches first < second.
struct PbwtMatch {
  std::size_t first = 0, second = 0, begin = 0, end = 0;
  friend auto operator<=>(const PbwtMatch&, const PbwtMatch&) = default;
};

// A match of a new word with the panel's word `word`, on [begin, end).
struct PbwtQueryMatch {
  std::size_t word = 0, begin = 0, end = 0;
  friend auto operator<=>(const PbwtQueryMatch&, const PbwtQueryMatch&) = default;
};

namespace detail {

template <class Index, class R>
auto pbwt_rows(const std::vector<R>& panel, std::size_t sigma, const char* who) {
  using T = std::remove_cv_t<std::ranges::range_value_t<R>>;
  std::vector<std::span<const T>> rows;
  rows.reserve(panel.size());
  const std::size_t N = panel.empty() ? 0 : std::ranges::size(panel[0]);
  if (panel.size() >= static_cast<std::size_t>(std::numeric_limits<Index>::max()) || N >= static_cast<std::size_t>(std::numeric_limits<Index>::max()))
    throw std::invalid_argument(std::string("morphe::formal::") + who + ": the panel is too large for the index type");
  if (sigma == 0 && N > 0) throw std::invalid_argument(std::string("morphe::formal::") + who + ": an alphabet with letters");
  for (const auto& w : panel) {
    const auto s = letters(w);
    if (s.size() != N) throw std::invalid_argument(std::string("morphe::formal::") + who + ": words of one length");
    for (auto c : s)
      if (static_cast<std::size_t>(c) >= sigma) throw std::invalid_argument(std::string("morphe::formal::") + who + ": a letter outside the alphabet [0, sigma)");
    rows.push_back(s);
  }
  return rows;
}

}  // namespace detail

// ---- the sweep ------------------------------------------------------------------------------------
// Calls visit(k, a_k, d_k) for k = 0, 1, …, N, with a_k and d_k as spans valid during the call.
// O(N M σ) time, O(M + σ) words.
template <class Index = std::uint32_t, LetterRange R, class F>
void pbwt_sweep(const std::vector<R>& panel, std::size_t sigma, F&& visit) {
  const auto rows = detail::pbwt_rows<Index>(panel, sigma, "pbwt_sweep");
  const std::size_t M = panel.size(), N = M ? rows[0].size() : 0;
  std::vector<Index> a(M), d(M, 0), a2(M), d2(M);
  std::iota(a.begin(), a.end(), Index{0});
  std::vector<std::size_t> start(sigma + 1);
  std::vector<Index> high(sigma);
  for (std::size_t k = 0;; k++) {
    visit(k, std::span<const Index>(a), std::span<const Index>(d));
    if (k == N) break;
    std::fill(start.begin(), start.end(), 0);
    for (std::size_t j = 0; j < M; j++) start[static_cast<std::size_t>(rows[a[j]][k]) + 1]++;
    for (std::size_t c = 0; c < sigma; c++) start[c + 1] += start[c];
    // high[c]: the maximum of d_k since the last word with letter c (k + 1 before the first).
    std::fill(high.begin(), high.end(), static_cast<Index>(k + 1));
    for (std::size_t j = 0; j < M; j++) {
      for (std::size_t c = 0; c < sigma; c++) high[c] = std::max(high[c], d[j]);
      const std::size_t c = rows[a[j]][k];
      a2[start[c]] = a[j];
      d2[start[c]++] = high[c];
      high[c] = 0;
    }
    std::swap(a, a2);
    std::swap(d, d2);
  }
}

// ---- long matches ------------------------------------------------------------------------------------
// Every maximal match of length ≥ L (L ≥ 1), once, first < second. O(N M σ + output).
template <class Index = std::uint32_t, LetterRange R>
std::vector<PbwtMatch> long_matches(const std::vector<R>& panel, std::size_t sigma, std::size_t L) {
  if (L == 0) throw std::invalid_argument("morphe::formal::long_matches: a length L >= 1");
  const auto rows = detail::pbwt_rows<Index>(panel, sigma, "long_matches");
  const std::size_t M = panel.size(), N = M ? rows[0].size() : 0;
  std::vector<PbwtMatch> out;
  // For each letter (and σ, for the end), the block's earlier positions, and groups of them:
  // (start of their common factor with the current position, index of the group's first one).
  std::vector<std::vector<std::size_t>> seen(sigma + 1);
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> groups(sigma + 1);
  pbwt_sweep<Index>(panel, sigma, [&](std::size_t k, std::span<const Index> a, std::span<const Index> d) {
    if (k < L || M < 2) return;
    const std::size_t threshold = k - L;
    const auto letter = [&](std::size_t j) -> std::size_t { return k < N ? static_cast<std::size_t>(rows[a[j]][k]) : sigma; };
    std::size_t j0 = 0;
    for (std::size_t j1 = 1; j1 <= M; j1++) {
      if (j1 < M && d[j1] <= threshold) continue;
      if (j1 - j0 >= 2) {
        for (std::size_t c = 0; c <= sigma; c++) {
          seen[c].clear();
          groups[c].clear();
        }
        for (std::size_t q = j0; q < j1; q++) {
          if (q > j0)
            for (std::size_t c = 0; c <= sigma; c++) {
              auto& g = groups[c];
              std::size_t first = std::numeric_limits<std::size_t>::max();
              while (!g.empty() && g.back().first <= d[q]) {
                first = g.back().second;
                g.pop_back();
              }
              if (first != std::numeric_limits<std::size_t>::max()) g.emplace_back(d[q], first);
            }
          const std::size_t lq = letter(q);
          for (std::size_t c = 0; c <= sigma; c++) {
            if (c == lq && k < N) continue;
            const auto& g = groups[c];
            for (std::size_t t = 0; t < g.size(); t++) {
              const std::size_t stop = t + 1 < g.size() ? g[t + 1].second : seen[c].size();
              for (std::size_t u = g[t].second; u < stop; u++) {
                const std::size_t x = a[seen[c][u]], y = a[q];
                out.push_back({std::min(x, y), std::max(x, y), g[t].first, k});
              }
            }
          }
          seen[lq].push_back(q);
          groups[lq].emplace_back(0, seen[lq].size() - 1);
        }
      }
      j0 = j1;
    }
  });
  return out;
}

// ---- set-maximal matches ------------------------------------------------------------------------------------
// Every set-maximal match of every word of the panel, as (that word, the other, begin, end), once
// for each other word that attains it. O(N M σ + output).
template <class Index = std::uint32_t, LetterRange R>
std::vector<PbwtMatch> set_maximal_matches(const std::vector<R>& panel, std::size_t sigma) {
  const auto rows = detail::pbwt_rows<Index>(panel, sigma, "set_maximal_matches");
  const std::size_t M = panel.size(), N = M ? rows[0].size() : 0;
  constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();
  std::vector<PbwtMatch> out;
  std::vector<std::size_t> above(M), below(M), last(sigma);
  pbwt_sweep<Index>(panel, sigma, [&](std::size_t k, std::span<const Index> a, std::span<const Index> d) {
    if (M < 2) return;
    const auto dd = [&](std::size_t j) -> std::size_t { return j < M ? static_cast<std::size_t>(d[j]) : k; };
    if (k < N) {
      // The nearest positions above and below with the same letter at k.
      std::fill(last.begin(), last.end(), NONE);
      for (std::size_t j = 0; j < M; j++) {
        const std::size_t c = rows[a[j]][k];
        above[j] = last[c];
        last[c] = j;
      }
      std::fill(last.begin(), last.end(), NONE);
      for (std::size_t j = M; j-- > 0;) {
        const std::size_t c = rows[a[j]][k];
        below[j] = last[c];
        last[c] = j;
      }
    }
    for (std::size_t i = 0; i < M; i++) {
      const std::size_t s = std::min<std::size_t>(d[i], dd(i + 1));
      if (s >= k) continue;
      if (k < N) {
        bool extends = false;
        if (above[i] != NONE) {
          extends = true;
          for (std::size_t m = above[i] + 1; m <= i && extends; m++) extends = d[m] <= s;
        }
        if (!extends && below[i] != NONE) {
          extends = true;
          for (std::size_t m = i + 1; m <= below[i] && extends; m++) extends = d[m] <= s;
        }
        if (extends) continue;
      }
      for (std::size_t m = i; m-- > 0 && d[m + 1] <= s;) out.push_back({a[i], a[m], s, k});
      for (std::size_t m = i + 1; m < M && d[m] <= s; m++) out.push_back({a[i], a[m], s, k});
    }
  });
  return out;
}

// ---- the stored transform, and new words against it ------------------------------------------------------------------------------------
template <class Index = std::uint32_t>
class Pbwt {
 public:
  Pbwt() = default;

  template <LetterRange R>
  Pbwt(const std::vector<R>& panel, std::size_t sigma) : sigma_(sigma) {
    const auto rows = detail::pbwt_rows<Index>(panel, sigma, "Pbwt");
    M_ = panel.size();
    N_ = M_ ? rows[0].size() : 0;
    prefix_.reserve((N_ + 1) * M_);
    divergence_.reserve((N_ + 1) * M_);
    columns_.reserve(N_);
    below_.assign(N_ * (sigma + 1), 0);
    Word column(M_);
    pbwt_sweep<Index>(panel, sigma, [&](std::size_t k, std::span<const Index> a, std::span<const Index> d) {
      prefix_.insert(prefix_.end(), a.begin(), a.end());
      divergence_.insert(divergence_.end(), d.begin(), d.end());
      if (k == N_) return;
      for (std::size_t j = 0; j < M_; j++) {
        column[j] = static_cast<Letter>(rows[a[j]][k]);
        below_[k * (sigma + 1) + column[j] + 1]++;
      }
      for (std::size_t c = 0; c < sigma; c++) below_[k * (sigma + 1) + c + 1] += below_[k * (sigma + 1) + c];
      columns_.emplace_back(column, sigma);
    });
  }

  std::size_t words() const { return M_; }
  std::size_t length() const { return N_; }
  std::size_t alphabet() const { return sigma_; }

  // a_k and d_k, k ∈ [0, N].
  std::span<const Index> prefix_array(std::size_t k) const { return {prefix_.data() + k * M_, M_}; }
  std::span<const Index> divergence_array(std::size_t k) const { return {divergence_.data() + k * M_, M_}; }
  // y_j[k] = x_{a_k[j]}[k], k < N.
  Letter letter(std::size_t k, std::size_t j) const { return columns_[k].access(j); }
  // Durbin's w-map: w_k(j, c) = (letters of column k below c) + rank_c(column k, j), so that the
  // words of a_k[f, g) with letter c at k are a_{k+1}[w_k(f, c), w_k(g, c)). M for a letter
  // outside the alphabet, which no word has.
  std::size_t extend(std::size_t k, std::size_t j, Letter c) const {
    if (c >= sigma_) return M_;
    return below_[k * (sigma_ + 1) + c] + columns_[k].rank(c, j);
  }

  // Every set-maximal match of z with the panel's words, once for each word that attains it.
  template <LetterRange Q>
  std::vector<PbwtQueryMatch> matches(const Q& z) const {
    std::vector<PbwtQueryMatch> out;
    sweep(z, [&](std::size_t f, std::size_t g, std::size_t e, std::size_t k) {
      const auto a = prefix_array(k);
      for (std::size_t j = f; j < g; j++) out.push_back({a[j], e, k});
    }, nullptr);
    return out;
  }

  // For k ∈ [0, N], the start e_k of the longest factor z[e_k, k) that some word of the panel has
  // at [e_k, k) (e_k = k when z[k − 1] is nowhere at k − 1).
  template <LetterRange Q>
  std::vector<std::size_t> longest_matches(const Q& z) const {
    std::vector<std::size_t> starts;
    sweep(z, [](std::size_t, std::size_t, std::size_t, std::size_t) {}, &starts);
    return starts;
  }

  std::size_t bytes() const {
    std::size_t b = (prefix_.size() + divergence_.size()) * sizeof(Index) + below_.size() * sizeof(std::size_t);
    for (const auto& c : columns_) b += c.bytes();
    return b;
  }

 private:
  // The sweep of z: report(f, g, e, k) at each set-maximal [e, k) with the words a_k[f, g).
  template <LetterRange Q, class F>
  void sweep(const Q& query, F&& report, std::vector<std::size_t>* starts) const {
    const auto z = letters(query);
    if (z.size() != N_) throw std::invalid_argument("morphe::formal::Pbwt: a word of the panel's length");
    std::size_t e = 0, f = 0, g = M_;
    if (starts) starts->assign(1, 0);
    for (std::size_t k = 0; k < N_; k++) {
      const auto c = static_cast<std::size_t>(z[k]);
      if (c < sigma_) {
        const std::size_t f2 = extend(k, f, static_cast<Letter>(c)), g2 = extend(k, g, static_cast<Letter>(c));
        if (f2 < g2) {
          f = f2;
          g = g2;
          if (starts) starts->push_back(e);
          continue;
        }
      }
      if (e < k) report(f, g, e, k);
      const auto d = divergence_array(k);
      constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();
      std::size_t best = NONE;
      if (c < sigma_) {
        // The nearest word with c above f, and its divergence max d(i, f]; then below g.
        std::size_t high = 0;
        for (std::size_t m = f; m-- > 0;) {
          high = std::max<std::size_t>(high, d[m + 1]);
          if (letter(k, m) == c) {
            best = high;
            break;
          }
        }
        high = 0;
        for (std::size_t m = g; m < M_; m++) {
          high = std::max<std::size_t>(high, d[m]);
          if (high >= best) break;
          if (letter(k, m) == c) {
            best = high;
            break;
          }
        }
      }
      if (best == NONE) {
        // No word has z[k] at k: the longest match ending at k + 1 is empty.
        e = k + 1;
        f = 0;
        g = M_;
      } else {
        std::size_t f3 = f, g3 = g;
        while (f3 > 0 && d[f3] <= best) f3--;
        while (g3 < M_ && d[g3] <= best) g3++;
        e = best;
        f = extend(k, f3, static_cast<Letter>(c));
        g = extend(k, g3, static_cast<Letter>(c));
      }
      if (starts) starts->push_back(e);
    }
    if (e < N_) report(f, g, e, N_);
  }

  std::size_t M_ = 0, N_ = 0, sigma_ = 0;
  std::vector<Index> prefix_, divergence_;  // a_k and d_k, k ∈ [0, N], M each
  std::vector<WaveletMatrix> columns_;      // y_j[k] for each k < N
  std::vector<std::size_t> below_;          // for each k < N, the letters of column k below each c
};

}  // namespace morphe::formal
