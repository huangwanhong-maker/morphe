// The Burrows–Wheeler transform and its relatives: permutations of a word's letters, sorted by
// the contexts that follow them, that the LF map inverts.
//
// - The BWT with a sentinel. For w of length n over [0, σ), sort the n + 1 suffixes of w$ ($
//   smaller than every letter); L[i] is the letter before the i-th suffix, cyclically. L is
//   written over [0, σ + 1) with $ = 0 and each letter c as c + 1 ("shifted"); row 0 is the
//   suffix $, so L[0] is the last letter of w. Computed from the suffix array in O(n).
// - The LF map. LF(i) = C[L[i]] + rank_{L[i]}(L, i), with C[c] the number of letters of L
//   smaller than c, sends the row of a rotation to the row of the rotation one step to the
//   right: rows ending with the same letter keep their order when that letter is moved to the
//   front (Burrows and Wheeler 1994, "A block-sorting lossless data compression algorithm", DEC
//   SRC Research Report 124, §4, Algorithm D). Following LF from row 0 spells w backwards: the
//   inverse in O(n + σ).
// - The rotation BWT, Burrows and Wheeler's own: the n rotations of w sorted, with the row of w
//   itself (the primary index); inverted by LF from the primary row.
// - The ω-order of words: u <_ω v when u^ω < v^ω, the infinite periodic words compared
//   lexicographically; u^ω = v^ω exactly when u and v are powers of one word, and otherwise
//   they differ within their first |u| + |v| − gcd(|u|, |v|) letters (Fine and Wilf 1965). The
//   conjugates of a multiset of words are sorted in this order by prefix doubling on the
//   cyclic shifts, ⌈log₂(2 max |w|)⌉ rounds of a sort: O(N log N log max |w|) for N letters.
// - The extended BWT of a multiset of primitive words (Mantaci, Restivo, Rosone and Sciortino
//   2007, "An extension of the Burrows–Wheeler transform", Theoret. Comput. Sci. 387,
//   298–312): every conjugate of every word, sorted in the ω-order (equal
//   conjugates by the order of their words), and L the last letter of each; with the rows of the
//   words themselves it is inverted by LF. Without them, the cycles of the standard permutation
//   of any word L (LF, the stable sort of its positions by letter) spell a multiset of
//   primitive necklaces, and the extended BWT of that multiset is L again: the bijection of
//   Gessel and Reutenauer (1993, "Counting permutations with given cycle structure and descent
//   set", J. Combin. Theory A 64, 189–215; Crochemore, Désarménien and Perrin
//   2005, "A note on the Burrows–Wheeler transformation", Theoret. Comput. Sci. 332, 567–572)
//   between words and multisets of primitive necklaces with the same letters.
// - The bijective BWT (Gil and Scott 2012, "A bijective string sorting transform",
//   arXiv:1201.3077; Kufleitner 2009, "On bijective variants of the Burrows–Wheeler
//   transform", Proc. Prague Stringology Conference, 65–79): the extended BWT of the
//   Chen–Fox–Lyndon factors of w. It is a bijection of the words of length n with no sentinel
//   and no index: the inverse reads the cycles of LF as Lyndon words, each starting at the
//   smallest row of its cycle (the least conjugate), and concatenates them in non-increasing
//   order, which is the factorisation (its uniqueness makes the map injective).
// - Runs and move-to-front. r, the number of maximal runs of equal letters of L, measures how
//   well a BWT compresses and sizes run-length indexes; move-to-front (Bentley, Sleator, Tarjan
//   and Wei 1986, CACM 29, 320–330) codes each letter by its position in a list of the letters
//   moved to the front on use, turning the runs of L into runs of zeros. Both O(nσ) at worst.

#pragma once

#include <morphe/formal/suffix_array.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::formal {

// ---- the BWT with a sentinel ------------------------------------------------------------------------------------
// L of w$, shifted ($ = 0, c ↦ c + 1), from the suffix array of w. O(n).
template <LetterRange R, class Index>
Word bwt(const R& w, const std::vector<Index>& sa) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  if (sa.size() != n) throw std::invalid_argument("morphe::formal::bwt: a suffix array of the word");
  Word L(n + 1);
  L[0] = n ? static_cast<Letter>(s[n - 1]) + 1 : 0;
  for (std::size_t i = 0; i < n; i++) L[i + 1] = sa[i] == 0 ? 0 : static_cast<Letter>(s[sa[i] - 1]) + 1;
  return L;
}
template <LetterRange R>
Word bwt(const R& w) {
  return bwt(w, suffix_array(w));
}

// C[c] = the number of letters of L smaller than c, c in [0, σ]; C[σ] = |L|.
inline std::vector<std::size_t> letter_starts(const Word& L, std::size_t sigma) {
  std::vector<std::size_t> c(sigma + 1, 0);
  for (Letter x : L) {
    if (x >= sigma) throw std::invalid_argument("morphe::formal::letter_starts: a letter outside the alphabet");
    c[x + 1]++;
  }
  for (std::size_t a = 0; a < sigma; a++) c[a + 1] += c[a];
  return c;
}

// LF of a word L: LF(i) = C[L[i]] + rank_{L[i]}(L, i), the stable sort of positions by letter
// read backwards (the inverse of the standard permutation). O(n + σ).
inline std::vector<std::size_t> lf_map(const Word& L) {
  std::size_t sigma = 0;
  for (Letter x : L) sigma = std::max<std::size_t>(sigma, x + 1);
  auto next = letter_starts(L, sigma);
  std::vector<std::size_t> lf(L.size());
  for (std::size_t i = 0; i < L.size(); i++) lf[i] = next[L[i]]++;
  return lf;
}

// The standard permutation π of L: π(i) is the position in L of the i-th letter of the stable
// sort of L, so that sorted(L)[i] = L[π(i)]; it is the inverse of LF.
inline std::vector<std::size_t> standard_permutation(const Word& L) {
  const auto lf = lf_map(L);
  std::vector<std::size_t> pi(L.size());
  for (std::size_t i = 0; i < L.size(); i++) pi[lf[i]] = i;
  return pi;
}

// w from the shifted BWT of w$. O(n + σ).
inline Word inverse_bwt(const Word& L) {
  if (L.empty()) throw std::invalid_argument("morphe::formal::inverse_bwt: the BWT of w$ has at least the sentinel");
  if (std::count(L.begin(), L.end(), Letter{0}) != 1) throw std::invalid_argument("morphe::formal::inverse_bwt: exactly one sentinel");
  const auto lf = lf_map(L);
  const std::size_t n = L.size() - 1;
  Word w(n);
  std::size_t row = 0;  // the suffix $
  for (std::size_t t = 0; t < n; t++) {
    w[n - 1 - t] = L[row] - 1;
    row = lf[row];
  }
  return w;
}

// ---- the ω-order of conjugates ------------------------------------------------------------------------------------
// The rows of every conjugate of every word of a multiset, sorted in the ω-order; ties (equal
// infinite words) by word, then by offset. Returns, for each row, the pair (word, offset) of
// the conjugate starting at that offset.
inline std::vector<std::pair<std::size_t, std::size_t>> sort_conjugates(const std::vector<Word>& words) {
  std::vector<std::size_t> start(words.size() + 1, 0), length(words.size());
  std::size_t longest = 0;
  for (std::size_t j = 0; j < words.size(); j++) {
    if (words[j].empty()) throw std::invalid_argument("morphe::formal::sort_conjugates: non-empty words");
    length[j] = words[j].size();
    start[j + 1] = start[j] + length[j];
    longest = std::max(longest, length[j]);
  }
  const std::size_t N = start.back();
  std::vector<std::size_t> owner(N);
  for (std::size_t j = 0; j < words.size(); j++)
    for (std::size_t o = 0; o < length[j]; o++) owner[start[j] + o] = j;
  // rank[g]: the class of the first h letters of the infinite word from global position g.
  std::vector<std::uint64_t> rank(N), next(N);
  std::vector<std::size_t> order(N);
  std::iota(order.begin(), order.end(), std::size_t{0});
  for (std::size_t g = 0; g < N; g++) rank[g] = words[owner[g]][g - start[owner[g]]];
  const auto shift = [&](std::size_t g, std::size_t h) {
    const std::size_t j = owner[g];
    return start[j] + (g - start[j] + h) % length[j];
  };
  const auto relabel = [&](auto key) {
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      const auto ka = key(a), kb = key(b);
      return ka != kb ? ka < kb : a < b;
    });
    std::uint64_t r = 0;
    for (std::size_t i = 0; i < N; i++) {
      if (i > 0 && key(order[i]) != key(order[i - 1])) r++;
      next[order[i]] = r;
    }
    std::swap(rank, next);
  };
  relabel([&](std::size_t g) { return std::pair<std::uint64_t, std::uint64_t>{rank[g], 0}; });
  // After the round with h, ranks compare the first 2h letters; 2h ≥ 2 · longest decides.
  for (std::size_t h = 1; h < 2 * longest; h *= 2)
    relabel([&](std::size_t g) { return std::pair<std::uint64_t, std::uint64_t>{rank[g], rank[shift(g, h)]}; });
  // Equal ranks are equal infinite words; the tie order is the global position, word then offset.
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rank[a] != rank[b] ? rank[a] < rank[b] : a < b; });
  std::vector<std::pair<std::size_t, std::size_t>> rows(N);
  for (std::size_t i = 0; i < N; i++) rows[i] = {owner[order[i]], order[i] - start[owner[order[i]]]};
  return rows;
}

// ---- the extended BWT ------------------------------------------------------------------------------------
struct ExtendedBwt {
  Word last;                       // the last letter of each sorted conjugate
  std::vector<std::size_t> rows;   // the row of each word itself (its conjugate at offset 0)
};

inline ExtendedBwt extended_bwt(const std::vector<Word>& words) {
  const auto order = sort_conjugates(words);
  ExtendedBwt e;
  e.last.resize(order.size());
  e.rows.assign(words.size(), 0);
  for (std::size_t i = 0; i < order.size(); i++) {
    const auto [j, o] = order[i];
    const Word& w = words[j];
    e.last[i] = w[(o + w.size() - 1) % w.size()];
    if (o == 0) e.rows[j] = i;
  }
  return e;
}

// The words from their extended BWT and their rows: from each row, LF spells the word backwards
// for its length, which is the length of its cycle when the words are primitive.
inline std::vector<Word> inverse_extended_bwt(const ExtendedBwt& e) {
  const auto lf = lf_map(e.last);
  std::vector<Word> out;
  for (std::size_t row : e.rows) {
    Word w;
    std::size_t r = row;
    do {
      w.push_back(e.last[r]);
      r = lf[r];
    } while (r != row);
    std::reverse(w.begin(), w.end());
    out.push_back(std::move(w));
  }
  return out;
}

// The multiset of primitive necklaces of a word L under the Gessel–Reutenauer bijection, each
// as its Lyndon word (the conjugate of the smallest row of its cycle), sorted non-increasingly.
inline std::vector<Word> necklaces(const Word& L) {
  const auto lf = lf_map(L);
  std::vector<bool> seen(L.size(), false);
  std::vector<Word> out;
  for (std::size_t row = 0; row < L.size(); row++) {
    if (seen[row]) continue;
    Word w;
    std::size_t r = row;
    do {
      w.push_back(L[r]);
      seen[r] = true;
      r = lf[r];
    } while (r != row);
    std::reverse(w.begin(), w.end());
    out.push_back(std::move(w));
  }
  std::sort(out.begin(), out.end(), [](const Word& a, const Word& b) { return b < a; });
  return out;
}

// ---- the rotation BWT ------------------------------------------------------------------------------------
struct RotationBwt {
  Word last;
  std::size_t primary = 0;  // the row of w itself
};

template <LetterRange R>
RotationBwt rotation_bwt(const R& w) {
  const auto s = letters(w);
  if (s.empty()) return {};
  const auto e = extended_bwt({Word(s.begin(), s.end())});
  return {e.last, e.rows[0]};
}

inline Word inverse_rotation_bwt(const RotationBwt& b) {
  const std::size_t n = b.last.size();
  if (n == 0) return {};
  if (b.primary >= n) throw std::invalid_argument("morphe::formal::inverse_rotation_bwt: the primary row lies in the transform");
  const auto lf = lf_map(b.last);
  Word w(n);
  std::size_t r = b.primary;
  for (std::size_t t = 0; t < n; t++) {
    w[n - 1 - t] = b.last[r];
    r = lf[r];
  }
  return w;
}

// ---- the bijective BWT ------------------------------------------------------------------------------------
template <LetterRange R>
Word bijective_bwt(const R& w) {
  if (std::ranges::size(w) == 0) return {};
  return extended_bwt(lyndon_factors(w)).last;
}

inline Word inverse_bijective_bwt(const Word& L) {
  Word w;
  w.reserve(L.size());
  for (const auto& factor : necklaces(L)) w.insert(w.end(), factor.begin(), factor.end());
  return w;
}

// ---- runs and move-to-front ------------------------------------------------------------------------------------
template <LetterRange R>
std::size_t runs(const R& w) {
  const auto s = letters(w);
  std::size_t r = 0;
  for (std::size_t i = 0; i < s.size(); i++) r += i == 0 || s[i] != s[i - 1];
  return r;
}

template <LetterRange R>
Word move_to_front(const R& w, std::size_t sigma) {
  Word list(sigma);
  std::iota(list.begin(), list.end(), Letter{0});
  Word out;
  out.reserve(std::ranges::size(w));
  for (auto c : letters(w)) {
    if (static_cast<std::size_t>(c) >= sigma) throw std::invalid_argument("morphe::formal::move_to_front: a letter outside the alphabet");
    const auto it = std::find(list.begin(), list.end(), static_cast<Letter>(c));
    const auto k = static_cast<std::size_t>(it - list.begin());
    out.push_back(static_cast<Letter>(k));
    std::rotate(list.begin(), it, it + 1);
  }
  return out;
}

inline Word inverse_move_to_front(const Word& code, std::size_t sigma) {
  Word list(sigma);
  std::iota(list.begin(), list.end(), Letter{0});
  Word out;
  out.reserve(code.size());
  for (Letter k : code) {
    if (k >= sigma) throw std::invalid_argument("morphe::formal::inverse_move_to_front: a code outside the alphabet");
    const Letter c = list[k];
    out.push_back(c);
    std::rotate(list.begin(), list.begin() + k, list.begin() + k + 1);
  }
  return out;
}

}  // namespace morphe::formal
