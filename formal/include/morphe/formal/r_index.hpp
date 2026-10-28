// The run-length FM-index and the r-index: counting and locating the occurrences of a pattern
// in a word, in space proportional to r, the number of runs of the word's Burrows–Wheeler
// transform, rather than to its length. On highly repetitive words (many copies of one text
// with a few differences each) r is far smaller than n.
//
// Let T be a word of length n over [0, σ), L the BWT of T$ as bwt.hpp writes it (shifted: $ = 0
// and each letter c as c + 1; row i holds the suffix SA$[i], with SA$[0] = n, the suffix $), and
// write L = h_0^{ℓ_0} h_1^{ℓ_1} … h_{r−1}^{ℓ_{r−1}} with h_t ≠ h_{t+1}: run t covers the rows
// [s_t, s_{t+1}), s_0 = 0, s_r = n + 1.
//
// - The run-length FM-index (Mäkinen and Navarro 2005, "Succinct suffix arrays based on
//   run-length encoding", CPM; Nordic Journal of Computing 12). L is kept as its runs: the heads
//   H = h_0 … h_{r−1} in a wavelet matrix (succinct.hpp), the starts s_t, and, for each letter c,
//   the cumulative lengths cum_c[j] of the first j runs headed by c. A row i lies in the run t
//   with s_t ≤ i < s_{t+1}, found by bisection of the starts; the c's of L[0, i) are those of the
//   j = rank_c(H, t) runs of c before run t, and the i − s_t letters of run t before i if h_t = c:
//       rank_c(L, i) = cum_c[j] + [h_t = c] (i − s_t).
//   So access, rank and LF(i) = C[L[i]] + rank_{L[i]}(L, i) take O(log r + log σ) steps, and
//   backward search (fm_index.hpp) counts a pattern of length m in O(m (log r + log σ)) steps,
//   whatever n. Mäkinen and Navarro mark the run heads in bit vectors of length n; here the
//   starts are a sorted array searched by bisection, so that the index is O(r + σ) words, as in
//   the r-index below.
// - Locating, the r-index (Gagie, Navarro and Prezza 2018, "Optimal-time text indexing in
//   BWT-runs bounded space", SODA; 2020, "Fully functional suffix trees and optimal text
//   searching in BWT-runs bounded space", J. ACM 67). Two samples of SA$ are kept, both at the
//   boundaries of runs: SA$ at the last row of every run (r values), and, for every run start
//   i = s_t > 0, the pair (SA$[i], SA$[i − 1]) (r − 1 pairs, sorted by their first member).
//   * The toehold. Backward search carries, with the interval [b, e) of the rows of a suffix Q of
//     the pattern, the value SA$[e − 1]. Reading c, the new interval ends at LF(q) + 1 for q the
//     last row of [b, e) with L[q] = c. If L[e − 1] = c, then q = e − 1 and SA$[LF(q)] =
//     SA$[e − 1] − 1. Otherwise row q + 1 ≤ e − 1 holds another letter, so q is the last row of
//     its run, the last run of c before the run of e − 1 (found by select on H), and SA$[q] is a
//     sample: SA$[LF(q)] = SA$[q] − 1. Hence the pattern's interval comes with SA$ of its last
//     row, at no more than the cost of counting.
//   * The function φ(SA$[i]) = SA$[i − 1], i ≥ 1 (the Φ of Kärkkäinen, Manzini and Puglisi 2009,
//     "Permuted longest-common-prefix array", CPM). If row i ≥ 1 is not the start of a run, then
//     L[i − 1] = L[i] = c ≠ $ (the sentinel occurs once, a run of its own), so LF(i − 1) =
//     LF(i) − 1; with p = SA$[i], LF(i) is the row of p − 1, and
//         φ(p − 1) = SA$[LF(i) − 1] = SA$[LF(i − 1)] = SA$[i − 1] − 1 = φ(p) − 1.
//     By induction φ(p) = φ(p′) + (p − p′) for p′ the largest sampled first member at most p;
//     one exists, since the row of the suffix at 0 holds $ and starts a run. A bisection: O(log r).
//   The occurrences of P are SA$[e − 1] and its e − b − 1 successive images under φ: locate takes
//   O(m (log r + log σ) + occ log r) steps, and the sort of its answer O(occ log occ).
// - The size of r. r ≤ n + 1 always. For a power T = w^k of a primitive word w of length p, r ≤
//   2p + 1 whatever k ≥ 1. For each offset o, the suffixes that begin with the conjugate
//   w[o, p) w[0, o) are consecutive rows; since the conjugates of a primitive word are distinct,
//   they are the suffixes at the positions jp + o with j ≤ k − 2 (with j ≤ k − 1 when o = 0), and
//   their BWT letter is w[o − 1 mod p], except at the suffix at 0, the largest of offset 0, whose
//   letter is $: p + 1 runs at most. The p other suffixes (w[o, p) $ for 1 ≤ o < p, and $) lie
//   between these intervals, so r ≤ (p + 1) + p. In general r = O(z log² n), z the
//   number of phrases of the LZ77 parse (Kempa and Kociumaka 2020, "Resolution of the
//   Burrows–Wheeler transform conjecture", FOCS); a few mutations of each copy add runs near
//   them, which the checks report rather than bound.
//
// Construction: SA-IS (suffix_array.hpp) and the BWT, O(n + σ) time, holding the suffix array
// (n + 1 words) while building; the index then keeps O(r + σ) words and r ⌈log₂(σ + 1)⌉ (1 + 1/4)
// bits of heads. Building in O(r) working space (by prefix-free parsing, for instance) is not
// done here. Exact; no step is probabilistic. Letters of patterns outside [0, σ) occur nowhere;
// the empty pattern occurs at every position 0..n (all n + 1 rows), as in fm_index.hpp.

#pragma once

#include <morphe/formal/bwt.hpp>
#include <morphe/formal/succinct.hpp>
#include <morphe/formal/suffix_array.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::formal {

// ---- the run-length FM-index ------------------------------------------------------------------------------------
template <class Index = std::uint32_t>
class RunLengthFmIndex {
 public:
  RunLengthFmIndex() = default;

  // The index of a word over [0, σ): its BWT by SA-IS, then the runs.
  template <LetterRange R>
  RunLengthFmIndex(const R& text, std::size_t sigma) {
    build(bwt(text, suffix_array<Index>(text, sigma)), sigma);
  }

  // From the shifted BWT L of some w$ over [0, σ + 1), with exactly one sentinel 0.
  static RunLengthFmIndex from_bwt(const Word& L, std::size_t sigma) {
    RunLengthFmIndex x;
    x.build(L, sigma);
    return x;
  }

  std::size_t size() const { return n_; }
  std::size_t rows() const { return n_ + 1; }
  std::size_t alphabet() const { return sigma_; }
  // r, the number of runs of L.
  std::size_t runs() const { return starts_.empty() ? 0 : starts_.size() - 1; }

  // The run of a row, and the rows [run_begin(t), run_end(t)) and head letter of run t.
  std::size_t run_of(std::size_t row) const {
    return static_cast<std::size_t>(std::upper_bound(starts_.begin(), starts_.end(), static_cast<Index>(row)) - starts_.begin()) - 1;
  }
  std::size_t run_begin(std::size_t t) const { return starts_[t]; }
  std::size_t run_end(std::size_t t) const { return starts_[t + 1]; }
  Letter run_head(std::size_t t) const { return heads_.access(t); }
  // The runs of the (shifted) letter c before run t, and the j-th run of c (from 0), or r.
  std::size_t runs_before(Letter c, std::size_t t) const { return heads_.rank(c, t); }
  std::size_t run_of_letter(Letter c, std::size_t j) const { return heads_.select(c, j); }

  // L[row], shifted ($ = 0).
  Letter access(std::size_t row) const { return heads_.access(run_of(row)); }

  // The occurrences of the shifted letter c in L[0, i), 0 ≤ i ≤ n + 1.
  std::size_t rank(Letter c, std::size_t i) const {
    if (c > sigma_) return 0;
    if (i >= rows()) return C_[c + 1] - C_[c];
    const std::size_t t = run_of(i);
    const auto [h, j] = heads_.access_rank(t);
    if (h == c) return cum(c, j) + (i - starts_[t]);
    return cum(c, heads_.rank(c, t));
  }

  // LF(row) = C[L[row]] + rank_{L[row]}(L, row): the row of the suffix one position to the left.
  std::size_t lf(std::size_t row) const {
    const std::size_t t = run_of(row);
    const auto [h, j] = heads_.access_rank(t);
    return C_[h] + cum(h, j) + (row - starts_[t]);
  }

  // The interval [begin, end) of rows whose suffixes begin with the pattern; {0, 0} if none.
  template <LetterRange P>
  std::pair<std::size_t, std::size_t> interval(const P& pattern) const {
    std::size_t b = 0, e = rows();
    const auto p = letters(pattern);
    for (std::size_t k = p.size(); k-- > 0 && b < e;) {
      if (static_cast<std::size_t>(p[k]) >= sigma_) return {0, 0};
      const Letter c = static_cast<Letter>(p[k]) + 1;
      b = C_[c] + rank(c, b);
      e = C_[c] + rank(c, e);
    }
    return b < e ? std::pair<std::size_t, std::size_t>{b, e} : std::pair<std::size_t, std::size_t>{0, 0};
  }

  template <LetterRange P>
  std::size_t count(const P& pattern) const {
    const auto [b, e] = interval(pattern);
    return e - b;
  }

  // C[c], the letters of L smaller than the shifted letter c, for c in [0, σ + 1].
  const std::vector<std::size_t>& letter_counts() const { return C_; }

  std::size_t bytes() const {
    return heads_.bytes() + (starts_.size() + cum_.size()) * sizeof(Index) + (C_.size() + offset_.size()) * sizeof(std::size_t);
  }

 private:
  void build(const Word& L, std::size_t sigma) {
    if (L.empty() || std::count(L.begin(), L.end(), Letter{0}) != 1)
      throw std::invalid_argument("morphe::formal::RunLengthFmIndex: the BWT of w$ holds exactly one sentinel");
    if (L.size() >= static_cast<std::size_t>(std::numeric_limits<Index>::max()))
      throw std::invalid_argument("morphe::formal::RunLengthFmIndex: the word is too long for the index type");
    n_ = L.size() - 1;
    sigma_ = sigma;
    C_ = letter_starts(L, sigma + 1);  // also refuses a letter outside [0, σ + 1)
    Word heads;
    starts_.clear();
    for (std::size_t i = 0; i < L.size(); i++)
      if (i == 0 || L[i] != L[i - 1]) {
        starts_.push_back(static_cast<Index>(i));
        heads.push_back(L[i]);
      }
    starts_.push_back(static_cast<Index>(L.size()));
    heads_ = WaveletMatrix(heads, sigma + 1);
    // cum_c, for each letter, at offset_[c]: one more entry than the runs of c.
    std::vector<std::size_t> of(sigma + 1, 0);
    for (Letter h : heads) of[h]++;
    offset_.assign(sigma + 2, 0);
    for (std::size_t c = 0; c <= sigma; c++) offset_[c + 1] = offset_[c] + of[c] + 1;
    cum_.assign(offset_.back(), 0);
    std::fill(of.begin(), of.end(), 0);
    for (std::size_t t = 0; t < heads.size(); t++) {
      const Letter c = heads[t];
      const std::size_t at = offset_[c] + of[c]++;
      cum_[at + 1] = static_cast<Index>(cum_[at] + (starts_[t + 1] - starts_[t]));
    }
  }
  std::size_t cum(Letter c, std::size_t j) const { return cum_[offset_[c] + j]; }

  std::size_t n_ = 0, sigma_ = 0;
  std::vector<Index> starts_;        // s_0 … s_r
  WaveletMatrix heads_;              // h_0 … h_{r−1}, over [0, σ + 1)
  std::vector<Index> cum_;           // the cumulative lengths of each letter's runs
  std::vector<std::size_t> offset_;  // where each letter's cumulative lengths begin
  std::vector<std::size_t> C_;       // over the shifted alphabet, σ + 2 entries
};

// ---- the r-index ------------------------------------------------------------------------------------
template <class Index = std::uint32_t>
class RIndex {
 public:
  // The rows [begin, end) of a pattern and SA$[end − 1], the toehold; begin = end when it is absent.
  struct Toehold {
    std::size_t begin = 0, end = 0, last = 0;
  };

  RIndex() = default;

  template <LetterRange R>
  RIndex(const R& text, std::size_t sigma) {
    const auto sa = suffix_array<Index>(text, sigma);
    fm_ = RunLengthFmIndex<Index>::from_bwt(bwt(text, sa), sigma);
    const std::size_t n = fm_.size(), r = fm_.runs();
    const auto suffix = [&](std::size_t row) -> std::size_t { return row == 0 ? n : static_cast<std::size_t>(sa[row - 1]); };
    last_.resize(r);
    std::vector<std::pair<Index, Index>> pairs;
    pairs.reserve(r - 1);
    for (std::size_t t = 0; t < r; t++) {
      last_[t] = static_cast<Index>(suffix(fm_.run_end(t) - 1));
      if (t > 0) pairs.emplace_back(static_cast<Index>(suffix(fm_.run_begin(t))), static_cast<Index>(suffix(fm_.run_begin(t) - 1)));
    }
    std::sort(pairs.begin(), pairs.end());
    from_.resize(pairs.size());
    to_.resize(pairs.size());
    for (std::size_t k = 0; k < pairs.size(); k++) {
      from_[k] = pairs[k].first;
      to_[k] = pairs[k].second;
    }
  }

  std::size_t size() const { return fm_.size(); }
  std::size_t rows() const { return fm_.rows(); }
  std::size_t alphabet() const { return fm_.alphabet(); }
  std::size_t runs() const { return fm_.runs(); }
  // The values of SA$ kept: r at the ends of runs and r − 1 pairs at their starts, 2r − 1.
  std::size_t samples() const { return last_.size() + from_.size(); }
  const RunLengthFmIndex<Index>& counting_index() const { return fm_; }

  template <LetterRange P>
  std::pair<std::size_t, std::size_t> interval(const P& pattern) const { return fm_.interval(pattern); }
  template <LetterRange P>
  std::size_t count(const P& pattern) const { return fm_.count(pattern); }

  // Backward search with the toehold.
  template <LetterRange P>
  Toehold toehold(const P& pattern) const {
    const auto& C = fm_.letter_counts();
    std::size_t b = 0, e = rows(), last = last_.back();
    const auto p = letters(pattern);
    for (std::size_t k = p.size(); k-- > 0;) {
      if (static_cast<std::size_t>(p[k]) >= alphabet()) return {};
      const Letter c = static_cast<Letter>(p[k]) + 1;
      const std::size_t b2 = C[c] + fm_.rank(c, b), e2 = C[c] + fm_.rank(c, e);
      if (b2 >= e2) return {};
      const std::size_t t = fm_.run_of(e - 1);
      if (fm_.run_head(t) == c) last -= 1;
      else last = static_cast<std::size_t>(last_[fm_.run_of_letter(c, fm_.runs_before(c, t) - 1)]) - 1;
      b = b2;
      e = e2;
    }
    return {b, e, last};
  }

  // φ(p) = SA$[ISA$[p] − 1], for 0 ≤ p < n (the suffix at n is row 0, which has no predecessor).
  std::size_t phi(std::size_t p) const {
    if (p >= size()) throw std::invalid_argument("morphe::formal::RIndex::phi: a position below n");
    const auto it = std::upper_bound(from_.begin(), from_.end(), static_cast<Index>(p));
    const std::size_t k = static_cast<std::size_t>(it - from_.begin()) - 1;  // from_[0] = 0
    return static_cast<std::size_t>(to_[k]) + (p - static_cast<std::size_t>(from_[k]));
  }

  // Every position at which the pattern occurs, increasing.
  template <LetterRange P>
  std::vector<std::size_t> locate(const P& pattern) const {
    const Toehold h = toehold(pattern);
    std::vector<std::size_t> out;
    if (h.begin >= h.end) return out;
    out.reserve(h.end - h.begin);
    std::size_t p = h.last;
    out.push_back(p);
    for (std::size_t i = h.end - 1; i > h.begin; i--) {
      p = phi(p);
      out.push_back(p);
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  std::size_t bytes() const { return fm_.bytes() + (last_.size() + from_.size() + to_.size()) * sizeof(Index); }

 private:
  RunLengthFmIndex<Index> fm_;
  std::vector<Index> last_;     // SA$ at the last row of each run
  std::vector<Index> from_, to_;  // (SA$[s_t], SA$[s_t − 1]) for t ≥ 1, sorted by the first
};

}  // namespace morphe::formal
