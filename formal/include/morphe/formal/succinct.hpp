// Succinct structures: bit vectors with rank and select, and the wavelet matrix, which answers
// access, rank and select on a sequence over a large integer alphabet in O(log σ) steps of rank
// on bit vectors.
//
// - Rank and select. For a bit vector B of length n, rank₁(i) is the number of ones in B[0, i)
//   and select₁(k) the position of the k-th one (from 0); rank₀ and select₀ likewise for zeros.
//   Jacobson (1989, "Space-efficient static trees and graphs", FOCS, 549–554) answers rank in
//   O(1) with o(n) extra bits by counts at two granularities. The layout here is Vigna's rank9
//   (Vigna 2008, "Broadword implementation of rank/select queries", WEA, LNCS 5038, 154–168,
//   §3): for each block of 512 bits a 64-bit count of the ones before it and seven 9-bit counts
//   of the ones before each of its words, in two adjacent words, so that a rank is two loads
//   from one cache line and a popcount: 25 per cent extra space. Select is a binary search over
//   the block counts, then over the words of the block, then within a word: O(log n).
// - The wavelet matrix (Claude, Navarro and Ordóñez 2015, "The wavelet matrix: an efficient
//   wavelet tree for large alphabets", Information Systems 47, 15–32, §3): for a sequence S of
//   n values below σ, one bit vector per bit of the values, ℓ = ⌈log₂ σ⌉ of them, from the most
//   significant. Level 0 holds the top bits of S; each next level holds the next bits of the
//   values reordered stably with the values whose bit was 0 first, z_l of them. A value's path
//   is computed by rank at each level, so access, rank (occurrences of c in S[0, i)) and the
//   count of c take O(ℓ) ranks, and select O(ℓ) selects. n ℓ (1 + 1/4) bits in all, with no
//   pointers, whatever σ.
//
// Both are static: a BitVector is filled with set() and then build() computes its counts;
// rank and select are valid only after build(). Everything is exact; nothing is approximated.

#pragma once

#include <morphe/formal/words.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace morphe::formal {

class BitVector {
 public:
  BitVector() = default;
  explicit BitVector(std::size_t n) : size_(n), words_(((n + 511) / 512) * 8 + 8, 0) {}
  explicit BitVector(const std::vector<bool>& bits) : BitVector(bits.size()) {
    for (std::size_t i = 0; i < bits.size(); i++)
      if (bits[i]) set(i);
    build();
  }

  std::size_t size() const { return size_; }
  void set(std::size_t i, bool value = true) {
    if (value) words_[i >> 6] |= 1ULL << (i & 63);
    else words_[i >> 6] &= ~(1ULL << (i & 63));
  }
  bool operator[](std::size_t i) const { return (words_[i >> 6] >> (i & 63)) & 1; }
  const std::vector<std::uint64_t>& words() const { return words_; }

  // The counts of rank9: counts_[2b] the ones before block b, counts_[2b + 1] the seven 9-bit
  // counts of the ones before words 1..7 of the block, word k's at bits 9(k − 1).
  void build() {
    const std::size_t blocks = words_.size() / 8;
    counts_.assign(2 * blocks + 2, 0);
    std::uint64_t total = 0;
    for (std::size_t b = 0; b < blocks; b++) {
      counts_[2 * b] = total;
      std::uint64_t inner = 0, packed = 0;
      for (std::size_t k = 0; k < 8; k++) {
        if (k > 0) packed |= inner << (9 * (k - 1));
        inner += static_cast<std::uint64_t>(std::popcount(words_[8 * b + k]));
      }
      counts_[2 * b + 1] = packed;
      total += inner;
    }
    counts_[2 * blocks] = total;
    ones_ = total;
  }

  std::size_t ones() const { return ones_; }
  std::size_t zeros() const { return size_ - ones_; }

  // The ones in [0, i), 0 ≤ i ≤ n. O(1).
  std::size_t rank1(std::size_t i) const {
    const std::size_t b = i >> 9, w = (i >> 6) & 7;
    std::uint64_t r = counts_[2 * b];
    if (w) r += (counts_[2 * b + 1] >> (9 * (w - 1))) & 0x1FF;
    const unsigned off = static_cast<unsigned>(i & 63);
    if (off) r += static_cast<std::uint64_t>(std::popcount(words_[i >> 6] & ((1ULL << off) - 1)));
    return static_cast<std::size_t>(r);
  }
  std::size_t rank0(std::size_t i) const { return i - rank1(i); }

  // The position of the k-th one (k from 0), or n if there are k or fewer ones. O(log n).
  std::size_t select1(std::size_t k) const { return select(k, true); }
  std::size_t select0(std::size_t k) const { return select(k, false); }

  std::size_t bytes() const { return (words_.size() + counts_.size()) * 8; }

 private:
  // The k-th set bit of a word (k < popcount(x)).
  static unsigned select_in_word(std::uint64_t x, std::size_t k) {
    for (std::size_t j = 0; j < k; j++) x &= x - 1;
    return static_cast<unsigned>(std::countr_zero(x));
  }
  std::size_t select(std::size_t k, bool one) const {
    if (k >= (one ? ones_ : size_ - ones_)) return size_;
    const std::size_t blocks = words_.size() / 8;
    const auto before = [&](std::size_t b) -> std::uint64_t { return one ? counts_[2 * b] : 512 * b - counts_[2 * b]; };
    // The last block with fewer than k + 1 before it.
    std::size_t lo = 0, hi = blocks;  // before(lo) ≤ k < before(hi)
    while (hi - lo > 1) {
      const std::size_t mid = (lo + hi) / 2;
      if (before(mid) <= k) lo = mid;
      else hi = mid;
    }
    std::size_t rest = k - static_cast<std::size_t>(before(lo));
    for (std::size_t w = 8 * lo; w < 8 * lo + 8; w++) {
      const std::uint64_t x = one ? words_[w] : ~words_[w];
      const auto c = static_cast<std::size_t>(std::popcount(x));
      if (rest < c) return 64 * w + select_in_word(x, rest);
      rest -= c;
    }
    return size_;  // unreachable for a built vector
  }

  std::size_t size_ = 0, ones_ = 0;
  std::vector<std::uint64_t> words_;   // padded to whole blocks, with one spare block
  std::vector<std::uint64_t> counts_;
};

class WaveletMatrix {
 public:
  WaveletMatrix() = default;
  template <LetterRange R>
  WaveletMatrix(const R& values, std::size_t sigma) : size_(std::ranges::size(values)), sigma_(sigma) {
    if (sigma == 0 && size_ > 0) throw std::invalid_argument("morphe::formal::WaveletMatrix: an alphabet with letters");
    levels_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::bit_width(sigma > 0 ? sigma - 1 : 0)));
    Word current;
    current.reserve(size_);
    for (auto c : letters(values)) {
      if (static_cast<std::size_t>(c) >= sigma) throw std::invalid_argument("morphe::formal::WaveletMatrix: a value outside the alphabet");
      current.push_back(static_cast<Letter>(c));
    }
    Word next(size_);
    bits_.assign(levels_, BitVector(size_));
    zeros_.assign(levels_, 0);
    for (std::size_t l = 0; l < levels_; l++) {
      const std::size_t bit = levels_ - 1 - l;
      std::size_t z = 0;
      for (std::size_t i = 0; i < size_; i++) {
        if ((current[i] >> bit) & 1) bits_[l].set(i);
        else z++;
      }
      bits_[l].build();
      zeros_[l] = z;
      std::size_t a = 0, b = z;
      for (std::size_t i = 0; i < size_; i++) next[((current[i] >> bit) & 1) ? b++ : a++] = current[i];
      std::swap(current, next);
    }
  }

  std::size_t size() const { return size_; }
  std::size_t alphabet() const { return sigma_; }
  std::size_t levels() const { return levels_; }

  Letter access(std::size_t i) const {
    Letter c = 0;
    for (std::size_t l = 0; l < levels_; l++) {
      const bool b = bits_[l][i];
      c = (c << 1) | static_cast<Letter>(b);
      i = b ? zeros_[l] + bits_[l].rank1(i) : bits_[l].rank0(i);
    }
    return c;
  }

  // The occurrences of c in [0, i).
  std::size_t rank(Letter c, std::size_t i) const {
    if (c >= sigma_) return 0;
    std::size_t p = 0;
    for (std::size_t l = 0; l < levels_; l++) {
      if ((c >> (levels_ - 1 - l)) & 1) {
        p = zeros_[l] + bits_[l].rank1(p);
        i = zeros_[l] + bits_[l].rank1(i);
      } else {
        p = bits_[l].rank0(p);
        i = bits_[l].rank0(i);
      }
    }
    return i - p;
  }

  // The value at i and the occurrences of that value in [0, i), in one pass.
  std::pair<Letter, std::size_t> access_rank(std::size_t i) const {
    Letter c = 0;
    std::size_t p = 0;
    for (std::size_t l = 0; l < levels_; l++) {
      const bool b = bits_[l][i];
      c = (c << 1) | static_cast<Letter>(b);
      if (b) {
        p = zeros_[l] + bits_[l].rank1(p);
        i = zeros_[l] + bits_[l].rank1(i);
      } else {
        p = bits_[l].rank0(p);
        i = bits_[l].rank0(i);
      }
    }
    return {c, i - p};
  }

  std::size_t count(Letter c) const { return rank(c, size_); }

  // The position of the k-th occurrence of c (k from 0), or size() if there is none.
  std::size_t select(Letter c, std::size_t k) const {
    if (c >= sigma_) return size_;
    std::size_t p = 0, e = size_;
    for (std::size_t l = 0; l < levels_; l++) {
      if ((c >> (levels_ - 1 - l)) & 1) {
        p = zeros_[l] + bits_[l].rank1(p);
        e = zeros_[l] + bits_[l].rank1(e);
      } else {
        p = bits_[l].rank0(p);
        e = bits_[l].rank0(e);
      }
    }
    if (k >= e - p) return size_;
    std::size_t pos = p + k;
    for (std::size_t l = levels_; l-- > 0;) {
      if ((c >> (levels_ - 1 - l)) & 1) pos = bits_[l].select1(pos - zeros_[l]);
      else pos = bits_[l].select0(pos);
    }
    return pos;
  }

  std::size_t bytes() const {
    std::size_t b = zeros_.size() * 8;
    for (const auto& v : bits_) b += v.bytes();
    return b;
  }

 private:
  std::size_t size_ = 0, sigma_ = 0, levels_ = 1;
  std::vector<BitVector> bits_;
  std::vector<std::size_t> zeros_;
};

}  // namespace morphe::formal
