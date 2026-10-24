// Words over a finite alphabet: the free monoid, its morphisms, and the combinatorics of words
// that every index of words is built on.
//
// A word over the alphabet [0, σ) is a finite sequence of letters, here `Word`, a vector of
// 32-bit letters. The algorithms that read words take any contiguous range of unsigned integers
// (`LetterRange`: a vector of bytes, a span, an array), so that a library with its own codes
// (the genomics domain's bases) passes them without copying.
//
// - The free monoid A* (concatenation, the empty word as unit) and its morphisms: a morphism
//   φ: A* → B* is determined by the images of the letters, φ(uv) = φ(u)φ(v) (Lothaire 1983,
//   Combinatorics on Words, §1.1). Iterated substitutions give the Thue–Morse word (0 → 01,
//   1 → 10; Thue 1912), the Fibonacci word (0 → 01, 1 → 0) and the period-doubling word
//   (0 → 01, 1 → 00) as limits φ^ω(0) of prolongable morphisms (Allouche and Shallit 2003,
//   Automatic Sequences, §§1.6, 7.1). Composition is applicative: (g ∘ f)(w) = g(f(w)).
// - Borders and periods. A border of w is a word both a proper prefix and a proper suffix; p is
//   a period of w when w[i] = w[i + p] wherever both are defined, and p is a period exactly
//   when n − p is the length of a border. The border array (the failure function of Knuth,
//   Morris and Pratt 1977, SIAM J. Comput. 6, 323–350, §2) gives every border; KMP search
//   finds every occurrence of a pattern in O(n + m). The Z-array, Z[i] = lcp(w, w[i..]), by
//   Gusfield's Z algorithm (Gusfield 1997, Algorithms on Strings, Trees and Sequences, §1.3–1.4),
//   in O(n). Fine and Wilf (1965): a word with periods p and q and length at least
//   p + q − gcd(p, q) has the period gcd(p, q).
// - Primitivity. w is primitive when it is no power u^e with e ≥ 2; every non-empty word is a
//   power of a unique primitive word, its primitive root, of length the smallest period p of w
//   when p divides |w| and |w| otherwise (Lyndon and Schützenberger 1962, Michigan Math. J. 9,
//   289–298; Lothaire 1983, §1.3). uv = vu exactly when u and v are powers of one word (ibid.).
// - Lyndon words: primitive words strictly smaller than all their proper rotations, in the
//   lexicographic order (Lyndon 1954; Lothaire 1983, §5.1). Every word factors uniquely as a
//   non-increasing product of Lyndon words (Chen, Fox and Lyndon 1958, Ann. Math. 68, 81–95;
//   Lothaire 1983, §5.1); Duval's algorithm finds the factors in O(n) time and O(1)
//   extra space (Duval 1983, J. Algorithms 4, 363–381, §2). Duval's later algorithm generates
//   the Lyndon words of length at most n over k letters in lexicographic order, in constant
//   amortised time (Duval 1988, Theoret. Comput. Sci. 60, 255–283); their number of length n is
//   Witt's (necklace) formula (1/n) Σ_{d|n} μ(d) k^{n/d} (Lothaire 1983, §5.1).
// - Conjugacy. u and v are conjugate when u = xy and v = yx: the rotations of a word. The least
//   rotation by Booth's algorithm (Booth 1980, Inf. Process. Lett. 10, 240–242), O(n), and
//   conjugacy as equality of least rotations.
//
// Every algorithm here terminates on every input; complexities are stated with each. Letters
// are compared as unsigned integers, so the lexicographic order of words is that of their
// letters, a proper prefix being smaller.

#pragma once

#include <morphe/arithmetic/integer.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace morphe::formal {

using Letter = std::uint32_t;
using Word = std::vector<Letter>;

// A contiguous range of unsigned integers, read as a word.
template <class R>
concept LetterRange = std::ranges::contiguous_range<R> && std::ranges::sized_range<R> &&
                      std::unsigned_integral<std::remove_cv_t<std::ranges::range_value_t<R>>>;

template <LetterRange R>
auto letters(const R& w) {
  using T = std::remove_cv_t<std::ranges::range_value_t<R>>;
  return std::span<const T>(std::ranges::data(w), std::ranges::size(w));
}

// A word from text, letter c ↦ c − first ('a' ↦ 0 by default), and back.
inline Word word(std::string_view text, char first = 'a') {
  Word w(text.size());
  for (std::size_t i = 0; i < text.size(); i++) w[i] = static_cast<Letter>(static_cast<unsigned char>(text[i]) - static_cast<unsigned char>(first));
  return w;
}
template <LetterRange R>
std::string text(const R& w, char first = 'a') {
  std::string s;
  s.reserve(std::ranges::size(w));
  for (auto c : letters(w)) s.push_back(static_cast<char>(first + static_cast<int>(c)));
  return s;
}

// ---- the free monoid ------------------------------------------------------------------------------------
inline Word concat(const Word& u, const Word& v) {
  Word w;
  w.reserve(u.size() + v.size());
  w.insert(w.end(), u.begin(), u.end());
  w.insert(w.end(), v.begin(), v.end());
  return w;
}
inline Word power(const Word& u, std::size_t e) {
  Word w;
  w.reserve(u.size() * e);
  for (std::size_t i = 0; i < e; i++) w.insert(w.end(), u.begin(), u.end());
  return w;
}
// The rotation w[r..] w[..r].
template <LetterRange R>
Word rotation(const R& w, std::size_t r) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  Word out(n);
  if (n == 0) return out;
  r %= n;
  for (std::size_t i = 0; i < n; i++) out[i] = s[(r + i) % n];
  return out;
}
template <LetterRange R>
Word reversal(const R& w) {
  const auto s = letters(w);
  return Word(s.rbegin(), s.rend());
}

// ---- morphisms of free monoids ------------------------------------------------------------------------------------
// φ: [0, σ)* → [0, τ)* given by the images of the letters; the image of a word is the
// concatenation of the images of its letters, so φ(uv) = φ(u)φ(v) and φ(ε) = ε.
struct Morphism {
  std::vector<Word> images;

  std::size_t domain() const { return images.size(); }

  template <LetterRange R>
  Word operator()(const R& w) const {
    const auto s = letters(w);
    std::size_t length = 0;
    for (auto c : s) length += at(c).size();
    Word out;
    out.reserve(length);
    for (auto c : s) {
      const Word& x = at(c);
      out.insert(out.end(), x.begin(), x.end());
    }
    return out;
  }
  // No letter is erased: |φ(w)| ≥ |w|.
  bool nonerasing() const {
    return std::all_of(images.begin(), images.end(), [](const Word& x) { return !x.empty(); });
  }
  // Every image has the same length (the morphism is uniform, of that modulus), or 0 if not.
  std::size_t modulus() const {
    if (images.empty()) return 0;
    for (const auto& x : images)
      if (x.size() != images[0].size()) return 0;
    return images[0].size();
  }
  // φ(a) = a u with u non-empty: φ^n(a) is then a prefix of φ^{n+1}(a), and the words converge
  // to the fixed point φ^ω(a).
  bool prolongable(Letter a) const { return a < images.size() && images[a].size() >= 2 && images[a][0] == a; }

 private:
  const Word& at(std::size_t c) const {
    if (c >= images.size()) throw std::invalid_argument("morphe::formal::Morphism: a letter outside the domain");
    return images[c];
  }
};

// g ∘ f, the morphism w ↦ g(f(w)).
inline Morphism compose(const Morphism& g, const Morphism& f) {
  Morphism h;
  h.images.reserve(f.images.size());
  for (const auto& x : f.images) h.images.push_back(g(x));
  return h;
}

// φ^n(w).
inline Word iterate(const Morphism& phi, Word w, std::size_t n) {
  for (std::size_t i = 0; i < n; i++) w = phi(w);
  return w;
}

// The prefix of length n of the fixed point φ^ω(a) of a prolongable, non-erasing morphism.
// O(n): the prefix is extended letter by letter by the images of its own letters.
inline Word fixed_point_prefix(const Morphism& phi, Letter a, std::size_t n) {
  if (!phi.prolongable(a)) throw std::invalid_argument("morphe::formal::fixed_point_prefix: phi(a) must begin with a and be longer than a letter");
  Word w{a};
  w.reserve(n + 64);
  for (std::size_t i = 0; w.size() < n; i++) {
    const Word& x = phi.images.at(w[i]);
    w.insert(w.end(), x.begin() + (i == 0 ? 1 : 0), x.end());
  }
  w.resize(n);
  return w;
}

inline Morphism thue_morse_morphism() { return {{{0, 1}, {1, 0}}}; }
inline Morphism fibonacci_morphism() { return {{{0, 1}, {0}}}; }
inline Morphism period_doubling_morphism() { return {{{0, 1}, {0, 0}}}; }

// The first n letters of the Thue–Morse word t = 0110100110010110…, t_i the parity of the
// binary digits of i.
inline Word thue_morse_word(std::size_t n) { return fixed_point_prefix(thue_morse_morphism(), 0, n); }
// The finite Fibonacci word f_m = φ^m(0) of length F_{m+2} (f_0 = 0, f_1 = 01, f_{m+1} = f_m f_{m−1}).
inline Word fibonacci_word(std::size_t m) { return iterate(fibonacci_morphism(), Word{0}, m); }
// The first n letters of the period-doubling word d = 0100010101000100…, d_i = ν₂(i + 1) mod 2.
inline Word period_doubling_word(std::size_t n) { return fixed_point_prefix(period_doubling_morphism(), 0, n); }

// ---- borders, periods, the Z-array, search ------------------------------------------------------------------------------------
// b[i] is the length of the longest proper border of w[0..i]: the failure function of Knuth,
// Morris and Pratt (their f, shifted to 0-based lengths). O(n): each step of the inner loop
// shortens the current border, which grows by at most one per letter.
template <LetterRange R>
std::vector<std::size_t> border_array(const R& w) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  std::vector<std::size_t> b(n, 0);
  std::size_t k = 0;
  for (std::size_t i = 1; i < n; i++) {
    while (k > 0 && s[i] != s[k]) k = b[k - 1];
    if (s[i] == s[k]) k++;
    b[i] = k;
  }
  return b;
}

// The lengths of every border of w, longest first (the iterates of the failure function), 0 last.
template <LetterRange R>
std::vector<std::size_t> borders(const R& w) {
  const auto s = letters(w);
  std::vector<std::size_t> out;
  if (s.empty()) return out;
  const auto b = border_array(w);
  for (std::size_t k = b.back(); k > 0; k = b[k - 1]) out.push_back(k);
  out.push_back(0);
  return out;
}

// Every period of w, smallest first: n − |u| for the borders u, then n itself.
template <LetterRange R>
std::vector<std::size_t> periods(const R& w) {
  const std::size_t n = std::ranges::size(w);
  std::vector<std::size_t> out;
  for (std::size_t k : borders(w)) out.push_back(n - k);
  return out;
}

template <LetterRange R>
std::size_t smallest_period(const R& w) {
  const std::size_t n = std::ranges::size(w);
  if (n == 0) return 0;
  return n - border_array(w).back();
}

template <LetterRange R>
bool has_period(const R& w, std::size_t p) {
  const auto s = letters(w);
  if (p == 0) return s.empty();
  for (std::size_t i = 0; i + p < s.size(); i++)
    if (s[i] != s[i + p]) return false;
  return true;
}

// Z[i] = the length of the longest common prefix of w and w[i..]; Z[0] = n. Gusfield's Z
// algorithm: the rightmost Z-box [l, r) seen so far gives Z[i] from Z[i − l] inside it. O(n).
template <LetterRange R>
std::vector<std::size_t> z_array(const R& w) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  std::vector<std::size_t> z(n, 0);
  if (n == 0) return z;
  z[0] = n;
  std::size_t l = 0, r = 0;
  for (std::size_t i = 1; i < n; i++) {
    std::size_t k = 0;
    if (i < r) k = std::min(r - i, z[i - l]);
    while (i + k < n && s[k] == s[i + k]) k++;
    z[i] = k;
    if (i + k > r) { l = i; r = i + k; }
  }
  return z;
}

// Every position at which `pattern` occurs in `text`, increasing; the empty pattern occurs at
// every position 0..n. KMP: O(m) to build the failure function, O(n) to scan.
template <LetterRange P, LetterRange T>
std::vector<std::size_t> find_all(const P& pattern, const T& text) {
  const auto p = letters(pattern);
  const auto t = letters(text);
  std::vector<std::size_t> out;
  const std::size_t m = p.size(), n = t.size();
  if (m == 0) {
    for (std::size_t i = 0; i <= n; i++) out.push_back(i);
    return out;
  }
  if (m > n) return out;
  const auto b = border_array(pattern);
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; i++) {
    while (k > 0 && (k == m || t[i] != p[k])) k = b[k - 1];
    if (t[i] == p[k]) k++;
    if (k == m) out.push_back(i + 1 - m);
  }
  return out;
}

// ---- primitivity ------------------------------------------------------------------------------------
// The length of the primitive root of w (0 for the empty word): the smallest period when it
// divides |w|, else |w|. O(n).
template <LetterRange R>
std::size_t primitive_root_length(const R& w) {
  const std::size_t n = std::ranges::size(w);
  if (n == 0) return 0;
  const std::size_t p = smallest_period(w);
  return n % p == 0 ? p : n;
}
template <LetterRange R>
Word primitive_root(const R& w) {
  const auto s = letters(w);
  return Word(s.begin(), s.begin() + static_cast<std::ptrdiff_t>(primitive_root_length(w)));
}
template <LetterRange R>
bool is_primitive(const R& w) {
  const std::size_t n = std::ranges::size(w);
  return n > 0 && primitive_root_length(w) == n;
}

// ---- Lyndon words and the Chen–Fox–Lyndon factorisation ------------------------------------------------------------------------------------
// The factorisation w = ℓ₁ ℓ₂ … ℓ_m into Lyndon words with ℓ₁ ≥ ℓ₂ ≥ … ≥ ℓ_m, as the starting
// positions of the factors followed by n. Duval's algorithm: the scan keeps w[i..j) of the form
// u^e u′ with u Lyndon and u′ a proper prefix of u, k − i being the position in u; a letter
// larger than w[k] makes w[i..j] Lyndon, a smaller one ends the run, whose copies of u are
// factors. O(n) comparisons (at most 2n), O(1) extra space besides the output.
template <LetterRange R>
std::vector<std::size_t> lyndon_factorisation(const R& w) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  std::vector<std::size_t> starts;
  std::size_t i = 0;
  while (i < n) {
    std::size_t j = i + 1, k = i;
    while (j < n && s[k] <= s[j]) {
      k = s[k] < s[j] ? i : k + 1;
      j++;
    }
    while (i <= k) {
      starts.push_back(i);
      i += j - k;
    }
  }
  starts.push_back(n);
  return starts;
}

// The factors themselves.
template <LetterRange R>
std::vector<Word> lyndon_factors(const R& w) {
  const auto s = letters(w);
  const auto b = lyndon_factorisation(w);
  std::vector<Word> out;
  for (std::size_t f = 0; f + 1 < b.size(); f++) out.emplace_back(s.begin() + static_cast<std::ptrdiff_t>(b[f]), s.begin() + static_cast<std::ptrdiff_t>(b[f + 1]));
  return out;
}

// w is a Lyndon word: non-empty, with a factorisation of one factor. O(n).
template <LetterRange R>
bool is_lyndon(const R& w) {
  const std::size_t n = std::ranges::size(w);
  return n > 0 && lyndon_factorisation(w).size() == 2;
}

// Every Lyndon word of length at most n over k letters, in increasing lexicographic order:
// f(word) for each. Duval 1988 (the algorithm of Fredricksen, Kessler and Maiorana restated):
// from a Lyndon word w of length m, repeat it periodically to length n, strip trailing
// maximal letters, and increment the last letter; the result is the next Lyndon word. Constant
// amortised time per word.
template <class F>
void for_each_lyndon_word(Letter k, std::size_t n, F&& f) {
  if (k == 0 || n == 0) return;
  Word w{0};
  while (!w.empty()) {
    f(static_cast<const Word&>(w));
    const std::size_t m = w.size();
    while (w.size() < n) w.push_back(w[w.size() - m]);
    while (!w.empty() && w.back() == k - 1) w.pop_back();
    if (!w.empty()) w.back()++;
  }
}

// The Möbius function μ(n) by trial division, O(√n).
inline int mobius(std::size_t n) {
  if (n == 0) throw std::invalid_argument("morphe::formal::mobius: n >= 1");
  int mu = 1;
  for (std::size_t p = 2; p * p <= n; p++) {
    if (n % p) continue;
    n /= p;
    if (n % p == 0) return 0;
    mu = -mu;
  }
  if (n > 1) mu = -mu;
  return mu;
}

// Witt's formula: the number of Lyndon words of length n ≥ 1 over k letters,
// (1/n) Σ_{d|n} μ(d) k^{n/d}, exactly.
inline arithmetic::Integer lyndon_count(std::size_t k, std::size_t n) {
  if (n == 0) throw std::invalid_argument("morphe::formal::lyndon_count: n >= 1");
  arithmetic::Integer sum = 0;
  for (std::size_t d = 1; d <= n; d++) {
    if (n % d) continue;
    const int mu = mobius(d);
    if (mu == 0) continue;
    const arithmetic::Integer term = arithmetic::pow(arithmetic::Integer(static_cast<unsigned long long>(k)), n / d);
    sum = mu > 0 ? sum + term : sum - term;
  }
  return arithmetic::exact_div(sum, arithmetic::Integer(static_cast<unsigned long long>(n)));
}

// ---- rotations and conjugacy ------------------------------------------------------------------------------------
// The smallest r such that the rotation w[r..] w[..r] is lexicographically least (0 for the
// empty word). Booth's algorithm: the failure function of the doubled word, computed relative to
// the current candidate k, which moves when a smaller letter is met. O(n). The result is reduced
// modulo the primitive root's length, since rotations by r and r + |root| agree.
template <LetterRange R>
std::size_t least_rotation(const R& w) {
  const auto s = letters(w);
  const std::size_t n = s.size();
  if (n == 0) return 0;
  const auto at = [&](std::size_t i) { return s[i % n]; };
  std::vector<std::ptrdiff_t> f(2 * n, -1);
  std::size_t k = 0;
  for (std::size_t j = 1; j < 2 * n; j++) {
    const auto sj = at(j);
    std::ptrdiff_t i = f[j - k - 1];
    while (i != -1 && sj != at(k + static_cast<std::size_t>(i) + 1)) {
      if (sj < at(k + static_cast<std::size_t>(i) + 1)) k = j - static_cast<std::size_t>(i) - 1;
      i = f[static_cast<std::size_t>(i)];
    }
    if (sj != at(k + static_cast<std::size_t>(i) + 1)) {  // here i == −1
      if (sj < at(k)) k = j;
      f[j - k] = -1;
    } else {
      f[j - k] = i + 1;
    }
  }
  return k % primitive_root_length(w);
}

// u and v are conjugate (u = xy, v = yx): equal length and equal least rotations. O(n).
template <LetterRange U, LetterRange V>
bool are_conjugate(const U& u, const V& v) {
  const auto a = letters(u);
  const auto b = letters(v);
  if (a.size() != b.size()) return false;
  if (a.empty()) return true;
  const std::size_t ra = least_rotation(u), rb = least_rotation(v);
  for (std::size_t i = 0; i < a.size(); i++)
    if (a[(ra + i) % a.size()] != b[(rb + i) % b.size()]) return false;
  return true;
}

// The Lyndon word conjugate to a primitive word: its least rotation.
template <LetterRange R>
Word lyndon_conjugate(const R& w) {
  return rotation(w, least_rotation(w));
}

}  // namespace morphe::formal
