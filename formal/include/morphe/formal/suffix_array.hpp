// Suffix arrays: the suffixes of a word in lexicographic order, the longest common prefixes of
// neighbours, and what they count.
//
// For a word w of length n over [0, σ), the suffix array SA is the permutation of [0, n) with
// w[SA[0]..) < w[SA[1]..) < … < w[SA[n − 1]..), a proper prefix being smaller (as if a sentinel
// smaller than every letter ended w). The inverse permutation is the rank of each suffix.
//
// - SA-IS (Nong, Zhang and Chan 2011, "Two efficient algorithms for linear time suffix array
//   construction", IEEE Trans. Computers 60, 1471–1484, §3): each suffix is S-type when it is
//   smaller than the suffix after it, else L-type; the leftmost S-type positions (LMS) cut w
//   into LMS substrings. Induced sorting places the LMS suffixes in their buckets, then sweeps
//   left to right placing each L-type suffix w[i − 1..) after w[i..), and right to left placing
//   each S-type suffix; done with the LMS suffixes in any order this sorts the LMS substrings
//   (ibid., §3), which are named by rank; the word of names, at most half as long, is
//   sorted recursively, and a last induced sort from the LMS suffixes in their true order sorts
//   every suffix. O(n + σ) time at each level and O(n) in all on an integer alphabet (σ ≤ n at
//   every level below the first); the working space is the suffix array, the types (n bytes),
//   and arrays of the LMS positions. The bucket layout (the LMS suffixes placed at the starts of
//   the S-parts of their buckets) follows the AtCoder Library's sa_is (CC0).
// - Prefix doubling (Manber and Myers 1993, SIAM J. Comput. 22, 935–948, §3): ranks by the
//   first 2^t letters, refined by sorting pairs of ranks, here by two counting sorts per round:
//   O(n log n). And the naive sort by comparison, O(n² log n) at worst. Both are the references
//   SA-IS is checked against.
// - The LCP array: LCP[i] = lcp(w[SA[i − 1]..), w[SA[i]..)), LCP[0] = 0, by Kasai, Lee,
//   Arimura, Arikawa and Park (2001, "Linear-time longest-common-prefix computation in suffix
//   arrays and its applications", CPM, LNCS 2089, 181–192): in text order the lcp with
//   the preceding suffix drops by at most one per step. O(n).
// - Distinct substrings: each suffix w[SA[i]..) contributes the n − SA[i] − LCP[i] prefixes that
//   no earlier suffix in the order has, so their number is n(n + 1)/2 − Σ LCP (Gusfield 1997,
//   §7.4 for the tree form). The number of distinct factors of each length follows the same way.
// - A set of words w_0, …, w_{m−1} is indexed as one word w_0 #_0 w_1 #_1 … w_{m−1} #_{m−1} with
//   distinct separators #_0 < #_1 < … smaller than every letter (the letters shifted by m): a
//   comparison of suffixes stops at the first separator, and suffixes equal up to it are
//   ordered by the word they lie in, so the order restricted to each word is that word's own.
//
// Indices are 32-bit by default (words shorter than 2^32 − 1); a 64-bit Index type indexes
// longer ones. Letters must lie below σ; a letter outside the alphabet throws.

#pragma once

#include <morphe/formal/words.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace morphe::formal {

namespace detail {

template <class Index, class T>
void check_alphabet(std::span<const T> s, std::size_t sigma, const char* who) {
  if (s.size() >= static_cast<std::size_t>(std::numeric_limits<Index>::max()))
    throw std::invalid_argument(std::string("morphe::formal::") + who + ": the word is too long for the index type");
  for (const auto c : s)
    if (static_cast<std::size_t>(c) >= sigma) throw std::invalid_argument(std::string("morphe::formal::") + who + ": a letter outside the alphabet [0, sigma)");
}

// The naive sort of the suffixes of s.
template <class Index, class T>
std::vector<Index> sa_naive(std::span<const T> s) {
  std::vector<Index> sa(s.size());
  std::iota(sa.begin(), sa.end(), Index{0});
  std::sort(sa.begin(), sa.end(), [&](Index a, Index b) {
    return std::lexicographical_compare(s.begin() + a, s.end(), s.begin() + b, s.end());
  });
  return sa;
}

// Prefix doubling with counting sorts: rank[i] is the class of w[i .. i + h), refined to
// (rank[i], rank[i + h]) each round, an absent second part (i + h ≥ n) ranking lowest.
template <class Index, class T>
std::vector<Index> sa_doubling(std::span<const T> s, std::size_t sigma) {
  const std::size_t n = s.size();
  std::vector<Index> sa(n), rank(n), tmp(n), next(n);
  if (n == 0) return sa;
  // Round 0: sort by the first letter.
  {
    std::vector<std::size_t> count(sigma + 1, 0);
    for (auto c : s) count[static_cast<std::size_t>(c) + 1]++;
    for (std::size_t c = 0; c < sigma; c++) count[c + 1] += count[c];
    for (std::size_t i = 0; i < n; i++) sa[count[s[i]]++] = static_cast<Index>(i);
    rank[sa[0]] = 0;
    for (std::size_t i = 1; i < n; i++) rank[sa[i]] = rank[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
  }
  std::vector<std::size_t> count(n + 1);
  for (std::size_t h = 1; rank[sa[n - 1]] + 1 < n; h *= 2) {
    // Sort by the second part: the suffixes with no second part first (in increasing position,
    // which is their order: shorter is smaller), then the others by the rank of i + h, which is
    // the current order of sa shifted by h.
    std::size_t t = 0;
    for (std::size_t i = n - h; i < n; i++) tmp[t++] = static_cast<Index>(i);
    for (std::size_t i = 0; i < n; i++)
      if (sa[i] >= h) tmp[t++] = static_cast<Index>(sa[i] - h);
    // Then stably by the first part.
    std::fill(count.begin(), count.end(), 0);
    for (std::size_t i = 0; i < n; i++) count[static_cast<std::size_t>(rank[i]) + 1]++;
    for (std::size_t c = 0; c < n; c++) count[c + 1] += count[c];
    for (std::size_t i = 0; i < n; i++) sa[count[rank[tmp[i]]]++] = tmp[i];
    const auto second = [&](Index i) -> std::ptrdiff_t { return i + h < n ? static_cast<std::ptrdiff_t>(rank[i + h]) : -1; };
    next[sa[0]] = 0;
    for (std::size_t i = 1; i < n; i++)
      next[sa[i]] = next[sa[i - 1]] + (rank[sa[i]] != rank[sa[i - 1]] || second(sa[i]) != second(sa[i - 1]));
    std::swap(rank, next);
  }
  return sa;
}

// SA-IS on s, letters in [0, sigma).
template <class Index, class T>
std::vector<Index> sa_is(std::span<const T> s, std::size_t sigma) {
  const std::size_t n = s.size();
  if (n == 0) return {};
  if (n == 1) return {0};
  if (n == 2) return s[0] < s[1] ? std::vector<Index>{0, 1} : std::vector<Index>{1, 0};
  if (n < 10) return sa_naive<Index>(s);
  if (n < 40) return sa_doubling<Index>(s, sigma);
  constexpr Index EMPTY = std::numeric_limits<Index>::max();

  // Types: S when w[i..) < w[i + 1..); the last suffix is L (it is larger than the empty one).
  std::vector<std::uint8_t> stype(n, 0);  // bytes, not bits: the induced sweeps read them at random
  for (std::size_t i = n - 1; i-- > 0;) stype[i] = s[i] == s[i + 1] ? stype[i + 1] : static_cast<std::uint8_t>(s[i] < s[i + 1]);

  // Buckets: begin_l[c] is the start of c's bucket (where its L-type suffixes go), begin_s[c]
  // the start of its S-part.
  std::vector<Index> begin_l(sigma + 1, 0), begin_s(sigma + 1, 0);
  for (std::size_t i = 0; i < n; i++) {
    if (!stype[i]) begin_s[s[i]]++;
    else begin_l[static_cast<std::size_t>(s[i]) + 1]++;
  }
  for (std::size_t c = 0; c <= sigma; c++) {
    begin_s[c] += begin_l[c];
    if (c < sigma) begin_l[c + 1] += begin_s[c];
  }

  std::vector<Index> sa(n);
  std::vector<Index> cursor(sigma + 1);
  // Induced sorting from LMS suffixes listed in `lms` (in the order they should keep).
  const auto induce = [&](const std::vector<Index>& lms) {
    std::fill(sa.begin(), sa.end(), EMPTY);
    std::copy(begin_s.begin(), begin_s.end(), cursor.begin());
    for (Index d : lms) sa[cursor[s[d]]++] = d;
    // L-type, left to right; the last suffix first, as if induced by the sentinel.
    std::copy(begin_l.begin(), begin_l.end(), cursor.begin());
    sa[cursor[s[n - 1]]++] = static_cast<Index>(n - 1);
    for (std::size_t i = 0; i < n; i++) {
      const Index v = sa[i];
      if (v != EMPTY && v >= 1 && !stype[v - 1]) sa[cursor[s[v - 1]]++] = v - 1;
    }
    // S-type, right to left, from the ends of the buckets.
    std::copy(begin_l.begin(), begin_l.end(), cursor.begin());
    for (std::size_t i = n; i-- > 0;) {
      const Index v = sa[i];
      if (v != EMPTY && v >= 1 && stype[v - 1]) sa[--cursor[static_cast<std::size_t>(s[v - 1]) + 1]] = v - 1;
    }
  };

  // The LMS positions, in text order, and the index of each among them.
  std::vector<Index> lms_index(n + 1, EMPTY);
  std::vector<Index> lms;
  for (std::size_t i = 1; i < n; i++)
    if (!stype[i - 1] && stype[i]) {
      lms_index[i] = static_cast<Index>(lms.size());
      lms.push_back(static_cast<Index>(i));
    }
  const std::size_t m = lms.size();
  induce(lms);
  if (m == 0) return sa;

  // The LMS substrings in sorted order, named by rank: equal names for equal substrings.
  std::vector<Index> sorted;
  sorted.reserve(m);
  for (Index v : sa)
    if (v != EMPTY && lms_index[v] != EMPTY) sorted.push_back(v);
  std::vector<Index> names(m);
  Index name = 0;
  names[lms_index[sorted[0]]] = 0;
  for (std::size_t i = 1; i < m; i++) {
    std::size_t l = sorted[i - 1], r = sorted[i];
    const std::size_t end_l = lms_index[l] + 1 < m ? lms[lms_index[l] + 1] : n;
    const std::size_t end_r = lms_index[r] + 1 < m ? lms[lms_index[r] + 1] : n;
    bool same = end_l - l == end_r - r;
    if (same) {
      while (l < end_l && s[l] == s[r]) { l++; r++; }
      if (l == n || s[l] != s[r]) same = false;
    }
    if (!same) name++;
    names[lms_index[sorted[i]]] = name;
  }

  // The word of names, sorted recursively, gives the true order of the LMS suffixes.
  const auto rec = sa_is<Index, Index>(std::span<const Index>(names), static_cast<std::size_t>(name) + 1);
  for (std::size_t i = 0; i < m; i++) sorted[i] = lms[rec[i]];
  induce(sorted);
  return sa;
}

}  // namespace detail

// The suffix array by SA-IS, O(n + σ).
template <class Index = std::uint32_t, LetterRange R>
std::vector<Index> suffix_array(const R& w, std::size_t sigma) {
  const auto s = letters(w);
  detail::check_alphabet<Index>(s, sigma, "suffix_array");
  return detail::sa_is<Index>(s, sigma);
}
// The alphabet taken as [0, max + 1).
template <class Index = std::uint32_t, LetterRange R>
std::vector<Index> suffix_array(const R& w) {
  const auto s = letters(w);
  std::size_t sigma = 0;
  for (auto c : s) sigma = std::max(sigma, static_cast<std::size_t>(c) + 1);
  return suffix_array<Index>(w, sigma);
}

// The references: prefix doubling, O(n log n), and the naive sort.
template <class Index = std::uint32_t, LetterRange R>
std::vector<Index> suffix_array_doubling(const R& w, std::size_t sigma) {
  const auto s = letters(w);
  detail::check_alphabet<Index>(s, sigma, "suffix_array_doubling");
  return detail::sa_doubling<Index>(s, sigma);
}
template <class Index = std::uint32_t, LetterRange R>
std::vector<Index> suffix_array_naive(const R& w) {
  const auto s = letters(w);
  detail::check_alphabet<Index>(s, std::numeric_limits<std::size_t>::max(), "suffix_array_naive");
  return detail::sa_naive<Index>(s);
}

// The inverse permutation: rank[SA[i]] = i.
template <class Index>
std::vector<Index> inverse_suffix_array(const std::vector<Index>& sa) {
  std::vector<Index> rank(sa.size());
  for (std::size_t i = 0; i < sa.size(); i++) rank[sa[i]] = static_cast<Index>(i);
  return rank;
}

// LCP[i] = lcp of the suffixes SA[i − 1] and SA[i], LCP[0] = 0: Kasai et al., O(n).
template <class Index, LetterRange R>
std::vector<Index> lcp_array(const R& w, const std::vector<Index>& sa) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  if (sa.size() != n) throw std::invalid_argument("morphe::formal::lcp_array: a suffix array of the word");
  std::vector<Index> lcp(n, 0);
  if (n == 0) return lcp;
  const auto rank = inverse_suffix_array(sa);
  std::size_t h = 0;
  for (std::size_t i = 0; i < n; i++) {
    if (rank[i] == 0) { h = 0; continue; }
    const std::size_t j = sa[rank[i] - 1];
    while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
    lcp[rank[i]] = static_cast<Index>(h);
    if (h > 0) h--;
  }
  return lcp;
}

// The number of distinct non-empty factors, n(n + 1)/2 − Σ LCP (n < 2^32).
template <class Index>
std::uint64_t distinct_substrings(std::size_t n, const std::vector<Index>& lcp) {
  std::uint64_t total = static_cast<std::uint64_t>(n) * (n + 1) / 2;
  for (Index h : lcp) total -= h;
  return total;
}

// The number of distinct factors of each length 0..n (p(0) = 1, the empty word), the factor
// complexity: a suffix of length L in position i of the order starts new factors of the
// lengths LCP[i] + 1 .. L. O(n) by a difference array.
template <class Index>
std::vector<std::uint64_t> factor_complexity(const std::vector<Index>& sa, const std::vector<Index>& lcp) {
  const std::size_t n = sa.size();
  std::vector<std::int64_t> diff(n + 2, 0);
  for (std::size_t i = 0; i < n; i++) {
    const std::size_t length = n - sa[i];
    diff[static_cast<std::size_t>(lcp[i]) + 1]++;
    diff[length + 1]--;
  }
  std::vector<std::uint64_t> p(n + 1, 0);
  std::int64_t run = 0;
  for (std::size_t l = 1; l <= n; l++) {
    run += diff[l];
    p[l] = static_cast<std::uint64_t>(run);
  }
  p[0] = 1;
  return p;
}

// ---- sets of words ------------------------------------------------------------------------------------
// w_0 #_0 w_1 #_1 … with #_j = j and every letter c written c + m: a word over [0, σ + m).
struct Concatenation {
  Word text;
  std::vector<std::size_t> start;  // where each word begins in text; start[m] = |text|
  std::size_t words = 0, sigma = 0;  // m, and the alphabet of the words

  // The word a position of text lies in, and its offset there (offset = |w_j| at the separator).
  std::pair<std::size_t, std::size_t> locate(std::size_t position) const {
    const auto it = std::upper_bound(start.begin(), start.end(), position);
    const std::size_t j = static_cast<std::size_t>(it - start.begin()) - 1;
    return {j, position - start[j]};
  }
  bool separator(std::size_t position) const { return text[position] < words; }
};

template <LetterRange R>
Concatenation concatenate(const std::vector<R>& words, std::size_t sigma) {
  Concatenation c;
  c.words = words.size();
  c.sigma = sigma;
  std::size_t total = words.size();
  for (const auto& w : words) total += std::ranges::size(w);
  c.text.reserve(total);
  for (std::size_t j = 0; j < words.size(); j++) {
    c.start.push_back(c.text.size());
    for (auto x : letters(words[j])) {
      if (static_cast<std::size_t>(x) >= sigma) throw std::invalid_argument("morphe::formal::concatenate: a letter outside the alphabet [0, sigma)");
      c.text.push_back(static_cast<Letter>(x + words.size()));
    }
    c.text.push_back(static_cast<Letter>(j));
  }
  c.start.push_back(c.text.size());
  return c;
}

// The suffix array of a set of words, as the suffix array of their concatenation.
template <class Index = std::uint32_t>
std::vector<Index> generalised_suffix_array(const Concatenation& c) {
  return suffix_array<Index>(c.text, c.sigma + c.words);
}

}  // namespace morphe::formal
