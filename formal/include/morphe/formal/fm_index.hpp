// The FM-index of a word over any integer alphabet: counts and positions of every occurrence of a
// pattern, and any factor of the word, from its Burrows–Wheeler transform and a sample of its
// suffix array (Ferragina and Manzini 2000, "Opportunistic data structures with applications",
// FOCS, 390–398, §3; Ferragina and Manzini 2005, J. ACM 52, 552–581).
//
// - Backward search. The rows of the BWT matrix of w$ whose rotations begin with a pattern P
//   form an interval [b, e). For P = cQ it is [C[c] + rank_c(L, b_Q), C[c] + rank_c(L, e_Q)):
//   the rows of Q preceded by c, moved by LF. Starting from all rows and reading P right to left
//   gives the interval of P in |P| steps; its size is the number of occurrences. With L in a
//   wavelet matrix (succinct.hpp), each step is O(log σ), so count is O(|P| log σ).
// - Locate. Suffix-array values are kept for the rows whose suffix starts at a multiple of the
//   sample rate s, marked in a bit vector with rank; from any other row LF is followed until a
//   marked row is reached, at most s − 1 steps, and the steps are added back. The row of the
//   suffix at 0 is marked, so LF is never taken from the row of $. O(s log σ) per occurrence.
// - Extract. The rows of the suffixes at multiples of s are kept too (the inverse suffix array,
//   sampled); w[b, e) is read backwards by LF from the first sampled position at or after e.
//   O((e − b + s) log σ).
//
// Space: the wavelet matrix of L, (n + 1)⌈log₂(σ + 1)⌉ (1 + 1/4) bits, and 2n/s sampled values.
// Construction: SA-IS (suffix_array.hpp), O(n + σ), with the suffix array held while building.
// Letters of patterns outside [0, σ) occur nowhere. Exact; no probabilistic step.

#pragma once

#include <morphe/formal/bwt.hpp>
#include <morphe/formal/succinct.hpp>
#include <morphe/formal/suffix_array.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::formal {

template <class Index = std::uint32_t>
class FmIndex {
 public:
  FmIndex() = default;

  template <LetterRange R>
  FmIndex(const R& text, std::size_t sigma, std::size_t sample = 32) : n_(std::ranges::size(text)), sigma_(sigma), sample_(sample) {
    if (sample == 0) throw std::invalid_argument("morphe::formal::FmIndex: a sample rate of at least 1");
    const auto sa = suffix_array<Index>(text, sigma);
    const Word L = bwt(text, sa);
    C_ = letter_starts(L, sigma + 1);
    L_ = WaveletMatrix(L, sigma + 1);
    // Row i of the BWT matrix holds the suffix at SA$[i]: n for row 0, SA[i − 1] after it.
    marked_ = BitVector(n_ + 1);
    for (std::size_t i = 0; i <= n_; i++)
      if (suffix(sa, i) % sample_ == 0) marked_.set(i);
    marked_.build();
    samples_.reserve(marked_.ones());
    isa_.assign(n_ / sample_ + 1, 0);
    for (std::size_t i = 0; i <= n_; i++) {
      const std::size_t p = suffix(sa, i);
      if (p % sample_ == 0) {
        samples_.push_back(static_cast<Index>(p));
        if (p < n_) isa_[p / sample_] = static_cast<Index>(i);
      }
    }
  }

  std::size_t size() const { return n_; }
  std::size_t alphabet() const { return sigma_; }
  std::size_t sample_rate() const { return sample_; }
  // The rows of the matrix, n + 1.
  std::size_t rows() const { return n_ + 1; }

  // The interval [begin, end) of rows whose rotations begin with the pattern.
  template <LetterRange P>
  std::pair<std::size_t, std::size_t> interval(const P& pattern) const {
    std::size_t b = 0, e = n_ + 1;
    const auto p = letters(pattern);
    for (std::size_t k = p.size(); k-- > 0 && b < e;) {
      if (static_cast<std::size_t>(p[k]) >= sigma_) return {0, 0};
      const Letter c = static_cast<Letter>(p[k]) + 1;
      b = C_[c] + L_.rank(c, b);
      e = C_[c] + L_.rank(c, e);
    }
    return b < e ? std::pair<std::size_t, std::size_t>{b, e} : std::pair<std::size_t, std::size_t>{0, 0};
  }

  template <LetterRange P>
  std::size_t count(const P& pattern) const {
    const auto [b, e] = interval(pattern);
    return e - b;
  }

  // The rows reached by LF: the row of the rotation one step to the right.
  std::size_t lf(std::size_t row) const {
    const auto [c, r] = L_.access_rank(row);
    return C_[c] + r;
  }

  // SA$[row]: the start of the suffix in that row (n for the row of $).
  std::size_t position(std::size_t row) const {
    std::size_t steps = 0;
    while (!marked_[row]) {
      row = lf(row);
      steps++;
    }
    return samples_[marked_.rank1(row)] + steps;
  }

  // Every position at which the pattern occurs, increasing.
  template <LetterRange P>
  std::vector<std::size_t> locate(const P& pattern) const {
    const auto [b, e] = interval(pattern);
    std::vector<std::size_t> out;
    out.reserve(e - b);
    for (std::size_t r = b; r < e; r++) out.push_back(position(r));
    std::sort(out.begin(), out.end());
    return out;
  }

  // w[begin, end).
  Word extract(std::size_t begin, std::size_t end) const {
    if (begin > end || end > n_) throw std::invalid_argument("morphe::formal::FmIndex::extract: 0 <= begin <= end <= n");
    Word out;
    if (begin == end) return out;
    // The first sampled position at or after end, and its row.
    std::size_t j = ((end + sample_ - 1) / sample_) * sample_;
    std::size_t row;
    if (j >= n_) { j = n_; row = 0; }
    else row = isa_[j / sample_];
    out.resize(end - begin);
    while (j > begin) {
      const auto [c, r] = L_.access_rank(row);
      j--;
      if (j < end) out[j - begin] = c - 1;
      row = C_[c] + r;
    }
    return out;
  }

  // The shifted BWT letter of a row ($ = 0).
  Letter bwt_letter(std::size_t row) const { return L_.access(row); }
  const std::vector<std::size_t>& letter_counts() const { return C_; }

  std::size_t bytes() const {
    return L_.bytes() + marked_.bytes() + (samples_.size() + isa_.size()) * sizeof(Index) + C_.size() * sizeof(std::size_t);
  }

 private:
  template <class SA>
  std::size_t suffix(const SA& sa, std::size_t row) const { return row == 0 ? n_ : static_cast<std::size_t>(sa[row - 1]); }

  std::size_t n_ = 0, sigma_ = 0, sample_ = 32;
  std::vector<std::size_t> C_;   // over the shifted alphabet [0, σ + 1)
  WaveletMatrix L_;
  BitVector marked_;
  std::vector<Index> samples_;   // SA$ of the marked rows, in row order
  std::vector<Index> isa_;       // the row of the suffix at each multiple of the sample rate
};

}  // namespace morphe::formal
