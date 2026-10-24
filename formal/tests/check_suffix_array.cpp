// Checks of suffix arrays: SA-IS against the naive sort on every binary word up to length 12 and on
// random words of alphabets 1 to 256 and lengths 0 to 10⁵; against its closed form on unary
// words; against prefix doubling on Fibonacci and Thue–Morse words of length 10⁵, where the
// naive sort is quadratic; the inverse permutation; Kasai's LCP against direct comparison; the
// count of distinct factors n(n + 1)/2 − Σ LCP against enumeration; the factor complexity of the
// Fibonacci word (ℓ + 1, Morse and Hedlund 1940) and of the Thue–Morse word (Brlek 1989) read
// from the LCP array; and a set of words indexed with distinct separators.

#include <morphe/check.hpp>
#include <morphe/formal/suffix_array.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

using namespace morphe::formal;

namespace {

struct Rng {
  std::uint64_t s;
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  std::uint64_t below(std::uint64_t n) { return next() % n; }
};

Word random_word(Rng& r, std::size_t n, Letter sigma) {
  Word w(n);
  for (auto& c : w) c = static_cast<Letter>(r.below(sigma));
  return w;
}

template <class W>
std::vector<std::uint32_t> naive_lcp(const W& w, const std::vector<std::uint32_t>& sa) {
  std::vector<std::uint32_t> lcp(sa.size(), 0);
  for (std::size_t i = 1; i < sa.size(); i++) {
    std::uint32_t h = 0;
    while (sa[i - 1] + h < w.size() && sa[i] + h < w.size() && w[sa[i - 1] + h] == w[sa[i] + h]) h++;
    lcp[i] = h;
  }
  return lcp;
}

// The factor complexity of the Thue–Morse word (Brlek 1989, Discrete Appl. Math. 24, 83–96;
// de Luca and Varricchio 1989, Theoret. Comput. Sci. 63, 333–348): p(1) = 2, p(2) = 4, and for
// n ≥ 3 with n − 1 = 2^r + q, 0 < q ≤ 2^r: p(n) = 6·2^{r−1} + 4q when q ≤ 2^{r−1}, else
// 8·2^{r−1} + 2q (OEIS A005942).
std::uint64_t thue_morse_complexity(std::uint64_t n) {
  if (n == 0) return 1;
  if (n == 1) return 2;
  if (n == 2) return 4;
  std::uint64_t r = 0;
  while ((1ULL << (r + 1)) < n - 1) r++;
  const std::uint64_t q = n - 1 - (1ULL << r);
  // 2^{r−1} may be a half: compare 2q with 2^r.
  return 2 * q <= (1ULL << r) ? 3 * (1ULL << r) + 4 * q : 4 * (1ULL << r) + 2 * q;
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{0x5341495300000001ULL};

  // ---- SA-IS against the naive sort ----------------------------------------------------------------------------------
  {
    bool exhaustive = true;
    std::size_t words = 0;
    for (std::size_t n = 0; n <= 12; n++)
      for (std::uint64_t x = 0; x < (1ULL << n); x++) {
        Word w(n);
        for (std::size_t i = 0; i < n; i++) w[i] = (x >> i) & 1;
        exhaustive = exhaustive && suffix_array(w, 2) == suffix_array_naive(w);
        words++;
      }
    s.require("SA-IS equals the naive sort on every binary word of length ≤ 12", exhaustive, std::to_string(words) + " words");

    bool random = true, doubling = true, bytes = true, wide = true;
    std::size_t letters_checked = 0;
    for (Letter sigma : {1u, 2u, 3u, 4u, 26u, 256u}) {
      for (std::size_t n : {0, 1, 2, 3, 5, 9, 10, 11, 39, 40, 41, 100, 1000, 10000, 100000}) {
        for (int trial = 0; trial < (n < 1000 ? 20 : 2); trial++) {
          const Word w = random_word(r, n, sigma);
          const auto sa = suffix_array(w, sigma);
          const auto naive = sigma == 1 && n > 2000 ? sa : suffix_array_naive(w);
          random = random && sa == naive;
          doubling = doubling && suffix_array_doubling(w, sigma) == sa;
          if (sigma <= 256) {
            std::vector<std::uint8_t> b(w.begin(), w.end());
            bytes = bytes && suffix_array(b, sigma) == sa;
          }
          if (n <= 1000) wide = wide && suffix_array<std::uint64_t>(w, sigma) == std::vector<std::uint64_t>(sa.begin(), sa.end());
          letters_checked += n;
        }
      }
    }
    s.require("SA-IS equals the naive sort on random words, alphabets 1 to 256, lengths 0 to 10⁵", random, std::to_string(letters_checked) + " letters");
    s.require("prefix doubling (Manber and Myers) equals SA-IS on the same words", doubling);
    s.require("a word given as bytes is indexed as the same word given as 32-bit letters", bytes);
    s.require("64-bit indices give the same suffix array", wide);

    // Unary words: a^n has SA = (n − 1, n − 2, …, 0), shorter suffixes first.
    bool unary = true;
    for (std::size_t n : {1, 2, 3, 10, 41, 1000, 100000}) {
      const auto sa = suffix_array(Word(n, 0), 1);
      for (std::size_t i = 0; i < n; i++) unary = unary && sa[i] == n - 1 - i;
    }
    s.require("the suffix array of a^n is n − 1, n − 2, …, 0 (n up to 10⁵)", unary);

    // Repetitive words, where the naive sort is quadratic: Fibonacci and Thue–Morse.
    const Word fib = fixed_point_prefix(fibonacci_morphism(), 0, 100000);
    const Word tm = thue_morse_word(100000);
    const auto sa_fib = suffix_array(fib, 2), sa_tm = suffix_array(tm, 2);
    s.require("SA-IS equals prefix doubling on the Fibonacci word of length 10⁵", sa_fib == suffix_array_doubling(fib, 2));
    s.require("SA-IS equals prefix doubling on the Thue–Morse word of length 10⁵", sa_tm == suffix_array_doubling(tm, 2));
    bool small_rep = true;
    for (std::size_t n : {10, 100, 1000, 3000}) {
      const Word f(fib.begin(), fib.begin() + static_cast<std::ptrdiff_t>(n)), t(tm.begin(), tm.begin() + static_cast<std::ptrdiff_t>(n));
      small_rep = small_rep && suffix_array(f, 2) == suffix_array_naive(f) && suffix_array(t, 2) == suffix_array_naive(t);
      const Word ff = fibonacci_word(10 + n / 1000);
      small_rep = small_rep && suffix_array(ff, 2) == suffix_array_naive(ff);
    }
    s.require("SA-IS equals the naive sort on Fibonacci and Thue–Morse prefixes up to 3000 and on finite Fibonacci words", small_rep);

    const auto rank = inverse_suffix_array(sa_fib);
    bool inverse = true;
    for (std::size_t i = 0; i < sa_fib.size(); i++) inverse = inverse && rank[sa_fib[i]] == i;
    s.require("the inverse suffix array inverts the permutation", inverse);

    // Time on 10⁷ random letters over four.
    const Word big = random_word(r, 10'000'000, 4);
    const auto t0 = std::chrono::steady_clock::now();
    const auto sa_big = suffix_array(big, 4);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    bool sorted = true;
    for (std::size_t i = 1; i < sa_big.size() && sorted; i += 9973)
      sorted = std::lexicographical_compare(big.begin() + sa_big[i - 1], big.end(), big.begin() + sa_big[i], big.end());
    char buf[96];
    std::snprintf(buf, sizeof buf, "%.2f s, %.0f ns a letter", seconds, seconds * 1e9 / 1e7);
    s.require("SA-IS sorts 10⁷ random letters over four (sampled neighbours in order)", sorted, buf);
  }

  // ---- LCP and distinct factors ----------------------------------------------------------------------------------
  {
    bool kasai = true, distinct = true, complexity = true;
    for (int trial = 0; trial < 400; trial++) {
      const Word w = random_word(r, r.below(80), 1 + static_cast<Letter>(r.below(4)));
      const auto sa = suffix_array(w);
      const auto lcp = lcp_array(w, sa);
      kasai = kasai && lcp == naive_lcp(w, sa);
      std::set<Word> factors;
      for (std::size_t i = 0; i < w.size(); i++)
        for (std::size_t j = i + 1; j <= w.size(); j++) factors.insert(Word(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(j)));
      distinct = distinct && distinct_substrings(w.size(), lcp) == factors.size();
      const auto p = factor_complexity(sa, lcp);
      std::vector<std::uint64_t> by_length(w.size() + 1, 0);
      by_length[0] = 1;
      for (const auto& f : factors) by_length[f.size()]++;
      complexity = complexity && p == by_length;
    }
    s.require("Kasai's LCP array equals direct comparison of neighbours", kasai);
    s.require("distinct factors n(n + 1)/2 − Σ LCP equals their enumeration", distinct);
    s.require("the factor complexity read from the LCP array equals the enumeration by length", complexity);
    for (std::size_t n : {1, 7, 1000}) {
      const Word a(n, 0);
      const auto sa = suffix_array(a, 1);
      s.equal("a^" + std::to_string(n) + " has n distinct factors", distinct_substrings(n, lcp_array(a, sa)), static_cast<std::uint64_t>(n));
    }

    const Word fib = fixed_point_prefix(fibonacci_morphism(), 0, 100000);
    const auto sa = suffix_array(fib, 2);
    const auto lcp = lcp_array(fib, sa);
    const auto p = factor_complexity(sa, lcp);
    bool sturmian = true;
    for (std::size_t l = 1; l <= 200; l++) sturmian = sturmian && p[l] == l + 1;
    s.require("the Fibonacci word is Sturmian: ℓ + 1 factors of each length ℓ ≤ 200 (Morse and Hedlund 1940)", sturmian);
    s.require("Kasai's LCP equals direct comparison on the Fibonacci word of length 10⁵", lcp == naive_lcp(fib, sa));

    const Word tm = thue_morse_word(100000);
    const auto sa_tm = suffix_array(tm, 2);
    const auto p_tm = factor_complexity(sa_tm, lcp_array(tm, sa_tm));
    bool brlek = true;
    for (std::size_t l = 0; l <= 200; l++) brlek = brlek && p_tm[l] == thue_morse_complexity(l);
    s.require("the factor complexity of the Thue–Morse word is Brlek's formula for ℓ ≤ 200", brlek,
              "p(3..10) = " + std::to_string(p_tm[3]) + ", " + std::to_string(p_tm[4]) + ", " + std::to_string(p_tm[5]) + ", " + std::to_string(p_tm[6]) + ", " +
                  std::to_string(p_tm[7]) + ", " + std::to_string(p_tm[8]) + ", " + std::to_string(p_tm[9]) + ", " + std::to_string(p_tm[10]));
  }

  // ---- sets of words ----------------------------------------------------------------------------------
  {
    bool order = true, restricted = true, located = true;
    for (int trial = 0; trial < 100; trial++) {
      std::vector<Word> words;
      const std::size_t m = 1 + r.below(6);
      for (std::size_t j = 0; j < m; j++) words.push_back(random_word(r, r.below(30), 2));
      const Concatenation c = concatenate(words, 2);
      const auto sa = generalised_suffix_array(c);
      // The order: compare up to and including the first separator, which ranks by its word.
      std::vector<std::uint32_t> want(c.text.size());
      for (std::size_t i = 0; i < want.size(); i++) want[i] = static_cast<std::uint32_t>(i);
      std::sort(want.begin(), want.end(), [&](std::uint32_t a, std::uint32_t b) {
        const auto ea = c.text.begin() + static_cast<std::ptrdiff_t>(c.start[c.locate(a).first + 1]);
        const auto eb = c.text.begin() + static_cast<std::ptrdiff_t>(c.start[c.locate(b).first + 1]);
        return std::lexicographical_compare(c.text.begin() + a, ea, c.text.begin() + b, eb);
      });
      order = order && sa == want;
      for (std::size_t j = 0; j < m; j++) {
        std::vector<std::uint32_t> own;
        for (auto p : sa) {
          const auto [word, offset] = c.locate(p);
          if (word == j && offset < words[j].size()) own.push_back(static_cast<std::uint32_t>(offset));
          located = located && (c.separator(p) == (offset == words[word].size()));
        }
        restricted = restricted && own == suffix_array_naive(words[j]);
      }
    }
    s.require("a set of words with distinct separators is sorted as each suffix up to its separator", order);
    s.require("the generalised suffix array restricted to one word is that word's suffix array", restricted);
    s.require("positions of the concatenation locate to their word and offset", located);
  }

  return s.done();
}
