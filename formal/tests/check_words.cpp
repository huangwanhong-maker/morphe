// Checks of the combinatorics of words: morphisms are monoid homomorphisms, and their fixed
// points are the Thue–Morse, Fibonacci and period-doubling words with their closed forms; the
// border array, periods, Z-array and KMP search against their definitions by brute force; the
// theorem of Fine and Wilf on every binary word up to length 14, and its bound attained; primitive
// roots and the theorem of Lyndon and Schützenberger; the Chen–Fox–Lyndon factorisation by
// Duval's algorithm against the greedy factorisation into longest Lyndon prefixes; Duval's
// enumeration of Lyndon words against Witt's formula and against brute force; Booth's least
// rotation and conjugacy against brute force.

#include <morphe/check.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace morphe::formal;
using morphe::arithmetic::Integer;

namespace {

// SplitMix64 (Steele, Lea and Flood 2014): the checks' own fixed-seed stream (morphe uses no calx).
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

// The binary word of length n whose letters are the bits of x.
Word binary(std::uint64_t x, std::size_t n) {
  Word w(n);
  for (std::size_t i = 0; i < n; i++) w[i] = (x >> i) & 1;
  return w;
}

std::vector<std::size_t> brute_borders(const Word& w) {
  std::vector<std::size_t> out;
  if (w.empty()) return out;
  for (std::size_t k = w.size() - 1; k + 1 > 0; k--)
    if (std::equal(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(k), w.end() - static_cast<std::ptrdiff_t>(k))) out.push_back(k);
  return out;
}

// Strictly smaller than every proper rotation.
bool brute_lyndon(const Word& w) {
  if (w.empty()) return false;
  for (std::size_t r = 1; r < w.size(); r++)
    if (!(w < rotation(w, r))) return false;
  return true;
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{20261024};

  // ---- morphisms ----------------------------------------------------------------------------------
  {
    bool hom = true, unit = true, comp = true;
    for (int trial = 0; trial < 200; trial++) {
      const Letter sigma = 1 + static_cast<Letter>(r.below(5)), tau = 1 + static_cast<Letter>(r.below(5));
      Morphism f, g;
      for (Letter a = 0; a < sigma; a++) f.images.push_back(random_word(r, r.below(4), tau));
      for (Letter a = 0; a < tau; a++) g.images.push_back(random_word(r, r.below(4), 3));
      const Word u = random_word(r, r.below(12), sigma), v = random_word(r, r.below(12), sigma);
      hom = hom && f(concat(u, v)) == concat(f(u), f(v));
      unit = unit && f(Word{}).empty();
      comp = comp && compose(g, f)(u) == g(f(u));
    }
    s.require("a substitution is a monoid homomorphism: φ(uv) = φ(u)φ(v) on 200 random morphisms and words", hom);
    s.require("the image of the empty word is the empty word", unit);
    s.require("composition is applicative: (g ∘ f)(w) = g(f(w))", comp);

    // Thue–Morse: t_i is the parity of the binary digits of i (Allouche and Shallit 2003, §1.6).
    const Word t = thue_morse_word(4096);
    bool parity = true;
    for (std::size_t i = 0; i < t.size(); i++) parity = parity && t[i] == static_cast<Letter>(std::popcount(i) & 1);
    s.require("the Thue–Morse word: t_i = s₂(i) mod 2 for i < 4096", parity);
    s.require("the Thue–Morse prefix of length 2^10 is μ^10(0)", iterate(thue_morse_morphism(), Word{0}, 10) == Word(t.begin(), t.begin() + 1024));
    s.require("the Thue–Morse morphism is uniform of modulus 2, the Fibonacci morphism is not uniform",
              thue_morse_morphism().modulus() == 2 && fibonacci_morphism().modulus() == 0 && fibonacci_morphism().nonerasing());
    // Overlap-freeness (Thue 1912): no factor axaxa, that is no factor of length 2p + 1 with period p.
    bool overlap_free = true;
    const std::size_t L = 600;
    for (std::size_t i = 0; i < L && overlap_free; i++)
      for (std::size_t p = 1; i + 2 * p + 1 <= L; p++) {
        bool periodic = true;
        for (std::size_t j = i; j + p < i + 2 * p + 1; j++)
          if (t[j] != t[j + p]) { periodic = false; break; }
        if (periodic) { overlap_free = false; break; }
      }
    s.require("the Thue–Morse word is overlap-free (Thue 1912): no factor axaxa in its first 600 letters", overlap_free);

    // Fibonacci words: |f_m| = F_{m+2}, |f_m|_1 = F_m, f_{m+1} = f_m f_{m−1}.
    bool fib = true;
    std::vector<std::size_t> F{0, 1};
    for (int i = 0; i < 30; i++) F.push_back(F[F.size() - 1] + F[F.size() - 2]);
    for (std::size_t m = 1; m < 22; m++) {
      const Word f = fibonacci_word(m);
      fib = fib && f.size() == F[m + 2] && static_cast<std::size_t>(std::count(f.begin(), f.end(), Letter{1})) == F[m] &&
            f == concat(fibonacci_word(m - 1), m >= 2 ? fibonacci_word(m - 2) : Word{1});
    }
    s.require("the finite Fibonacci words: |f_m| = F_{m+2}, F_m ones, and f_{m+1} = f_m f_{m−1}", fib);

    // Period doubling: d_i = ν₂(i + 1) mod 2 (Allouche and Shallit 2003, §6.3).
    const Word d = period_doubling_word(5000);
    bool pd = true;
    for (std::size_t i = 0; i < d.size(); i++) pd = pd && d[i] == static_cast<Letter>(std::countr_zero(i + 1) & 1);
    s.require("the period-doubling word: d_i = ν₂(i + 1) mod 2 for i < 5000", pd);
  }

  // ---- borders, periods, Z, KMP ----------------------------------------------------------------------------------
  {
    bool b_ok = true, p_ok = true, z_ok = true, kmp_ok = true;
    for (int trial = 0; trial < 3000; trial++) {
      const Word w = random_word(r, r.below(40), 1 + static_cast<Letter>(r.below(3)));
      const auto bb = brute_borders(w);
      b_ok = b_ok && borders(w) == bb;
      const auto ba = border_array(w);
      for (std::size_t i = 0; i < w.size(); i++) {
        const Word prefix(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(i + 1));
        b_ok = b_ok && ba[i] == brute_borders(prefix).front();
      }
      std::vector<std::size_t> per;
      for (std::size_t p = 1; p <= w.size(); p++)
        if (has_period(w, p)) per.push_back(p);
      p_ok = p_ok && periods(w) == per && (w.empty() || smallest_period(w) == per.front());
      const auto z = z_array(w);
      for (std::size_t i = 0; i < w.size(); i++) {
        std::size_t k = 0;
        while (i + k < w.size() && w[k] == w[i + k]) k++;
        z_ok = z_ok && z[i] == k;
      }
      const Word text = random_word(r, r.below(200), 2), pat = random_word(r, 1 + r.below(5), 2);
      std::vector<std::size_t> want;
      for (std::size_t i = 0; i + pat.size() <= text.size(); i++)
        if (std::equal(pat.begin(), pat.end(), text.begin() + static_cast<std::ptrdiff_t>(i))) want.push_back(i);
      kmp_ok = kmp_ok && find_all(pat, text) == want;
    }
    s.require("the border array and the borders agree with brute force (3000 words)", b_ok);
    s.require("the periods are n − |border|, the smallest one first, against brute force", p_ok);
    s.require("the Z-array agrees with brute force", z_ok);
    s.require("KMP finds every occurrence, against brute force", kmp_ok);
    s.require("the empty pattern occurs at every position 0..n", find_all(Word{}, word("abc")) == std::vector<std::size_t>{0, 1, 2, 3});
    // Gusfield 1997, §1.3: Z of aabcaabxaaz.
    s.equal("Z-array of aabcaabxaaz (Gusfield 1997, §1.3)", z_array(word("aabcaabxaaz")), std::vector<std::size_t>{11, 1, 0, 0, 3, 1, 0, 0, 2, 1, 0});
    // The pattern of Knuth, Morris and Pratt 1977, §2, abcabcacab: its borders, by hand.
    s.equal("border array of abcabcacab, the pattern of Knuth, Morris and Pratt 1977", border_array(word("abcabcacab")), std::vector<std::size_t>{0, 0, 0, 1, 2, 3, 4, 0, 1, 2});
  }

  // ---- Fine and Wilf ----------------------------------------------------------------------------------
  {
    bool theorem = true;
    std::set<std::pair<std::size_t, std::size_t>> attained;  // coprime (p, q) with an extremal word
    for (std::size_t n = 1; n <= 14; n++)
      for (std::uint64_t x = 0; x < (1ULL << n); x++) {
        const Word w = binary(x, n);
        const auto per = periods(w);
        for (std::size_t a = 0; a < per.size(); a++)
          for (std::size_t b = a + 1; b < per.size(); b++) {
            const std::size_t p = per[a], q = per[b], g = std::gcd(p, q);
            if (p + q - g <= n) theorem = theorem && has_period(w, g);
            else if (g == 1 && p + q - 2 == n && !has_period(w, 1)) attained.insert({p, q});
          }
      }
    s.require("Fine and Wilf: periods p, q and |w| ≥ p + q − gcd(p, q) give the period gcd(p, q), on every binary word of length ≤ 14", theorem);
    std::size_t coprime = 0;
    for (std::size_t p = 2; p <= 14; p++)
      for (std::size_t q = p + 1; p + q - 2 <= 14; q++) coprime += std::gcd(p, q) == 1;
    s.require("the bound is attained: for every coprime 2 ≤ p < q with p + q − 2 ≤ 14, a binary word of length p + q − 2 has periods p and q and not 1",
              attained.size() == coprime, std::to_string(attained.size()) + " of " + std::to_string(coprime) + " pairs");
  }

  // ---- primitivity ----------------------------------------------------------------------------------
  {
    bool roots = true, doubling = true, commute = true;
    for (int trial = 0; trial < 3000; trial++) {
      Word w = random_word(r, 1 + r.below(8), 1 + static_cast<Letter>(r.below(3)));
      if (trial % 2) w = power(w, 1 + r.below(4));
      const Word root = primitive_root(w);
      bool brute_primitive = true;
      for (std::size_t d = 1; d < w.size(); d++)
        if (w.size() % d == 0 && power(Word(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(d)), w.size() / d) == w) brute_primitive = false;
      roots = roots && is_primitive(root) && power(root, w.size() / root.size()) == w && is_primitive(w) == brute_primitive;
      // w is primitive exactly when it occurs in ww only at 0 and |w| (Lothaire 1983, §1.3).
      doubling = doubling && (find_all(w, concat(w, w)).size() == 2) == is_primitive(w);
      Word u = random_word(r, 1 + r.below(5), 2), v = random_word(r, 1 + r.below(5), 2);
      if (trial % 3 == 0) {
        const Word base = random_word(r, 1 + r.below(3), 2);
        u = power(base, 1 + r.below(3));
        v = power(base, 1 + r.below(3));
      }
      commute = commute && (concat(u, v) == concat(v, u)) == (primitive_root(u) == primitive_root(v));
    }
    s.require("the primitive root is primitive and w is a power of it; primitivity agrees with brute force", roots);
    s.require("w is primitive exactly when it occurs in ww only at 0 and |w|", doubling);
    s.require("Lyndon and Schützenberger: uv = vu exactly when u and v have one primitive root", commute);
  }

  // ---- Lyndon words ----------------------------------------------------------------------------------
  {
    bool factorisation = true, greedy = true;
    std::size_t words = 0;
    const auto check = [&](const Word& w) {
      words++;
      const auto factors = lyndon_factors(w);
      Word joined;
      for (std::size_t f = 0; f < factors.size(); f++) {
        factorisation = factorisation && brute_lyndon(factors[f]) && (f == 0 || !(factors[f - 1] < factors[f]));
        joined.insert(joined.end(), factors[f].begin(), factors[f].end());
      }
      factorisation = factorisation && joined == w;
      // The first factor is the longest Lyndon prefix (Lothaire 1983, §5.1); greedily.
      std::vector<Word> want;
      for (std::size_t i = 0; i < w.size();) {
        std::size_t best = 1;
        for (std::size_t l = 1; i + l <= w.size(); l++)
          if (brute_lyndon(Word(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(i + l)))) best = l;
        want.emplace_back(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(i + best));
        i += best;
      }
      greedy = greedy && factors == want && is_lyndon(w) == brute_lyndon(w);
    };
    for (std::size_t n = 0; n <= 12; n++)
      for (std::uint64_t x = 0; x < (1ULL << n); x++) check(binary(x, n));
    for (int trial = 0; trial < 2000; trial++) check(random_word(r, r.below(30), 2 + static_cast<Letter>(r.below(3))));
    s.require("Duval: the factors are Lyndon words, non-increasing, and their product is w", factorisation, std::to_string(words) + " words");
    s.require("Duval's factorisation is the greedy one by longest Lyndon prefixes (Chen–Fox–Lyndon uniqueness)", greedy);

    // Duval 1988's enumeration against Witt's formula, and against brute force.
    bool witt = true, brute = true, ordered = true;
    std::string detail;
    for (Letter k = 1; k <= 5; k++) {
      const std::size_t nmax = k == 1 ? 6 : k == 2 ? 20 : k == 3 ? 12 : k == 4 ? 10 : 8;
      std::vector<std::size_t> count(nmax + 1, 0);
      Word previous;
      bool first = true;
      for_each_lyndon_word(k, nmax, [&](const Word& w) {
        count[w.size()]++;
        ordered = ordered && is_lyndon(w) && (first || previous < w);
        previous = w;
        first = false;
      });
      for (std::size_t n = 1; n <= nmax; n++) witt = witt && lyndon_count(k, n) == Integer(static_cast<unsigned long long>(count[n]));
      for (std::size_t n = 1; n <= nmax; n++) {
        std::size_t c = 0;
        std::uint64_t total = 1;
        for (std::size_t i = 0; i < n; i++) total *= k;
        if (total > (1u << 16)) break;
        for (std::uint64_t x = 0; x < total; x++) {
          Word w(n);
          std::uint64_t y = x;
          for (std::size_t i = 0; i < n; i++, y /= k) w[i] = static_cast<Letter>(y % k);
          c += is_lyndon(w);
        }
        brute = brute && c == count[n];
      }
      detail += (k > 1 ? ", " : "") + std::string("L_") + std::to_string(k) + "(" + std::to_string(nmax) + ") = " + std::to_string(count[nmax]);
    }
    s.require("Duval's enumeration counts Witt's (1/n) Σ μ(d) k^{n/d} Lyndon words of each length, k = 1 to 5", witt, detail);
    s.require("Duval's enumeration agrees with brute force over all words of each length", brute);
    s.require("Duval's enumeration lists Lyndon words in increasing order", ordered);
    // Σ_{d|n} d L_k(d) = k^n (Lothaire 1983, §5.1), exactly, into big integers.
    bool identity = true;
    for (std::size_t k : {2, 4, 20})
      for (std::size_t n = 1; n <= 60; n++) {
        Integer sum = 0;
        for (std::size_t d = 1; d <= n; d++)
          if (n % d == 0) sum += Integer(static_cast<unsigned long long>(d)) * lyndon_count(k, d);
        identity = identity && sum == morphe::arithmetic::pow(Integer(static_cast<unsigned long long>(k)), n);
      }
    s.require("Σ_{d|n} d L_k(d) = k^n for k = 2, 4, 20 and n ≤ 60, exactly", identity, "L_2(60) = " + lyndon_count(2, 60).str());
    s.equal("the 335 binary Lyndon words of length 12 (OEIS A001037)", lyndon_count(2, 12), Integer(335));
  }

  // ---- rotations and conjugacy ----------------------------------------------------------------------------------
  {
    bool least = true, conj = true, lyndon = true;
    for (int trial = 0; trial < 4000; trial++) {
      Word w = random_word(r, r.below(25), 1 + static_cast<Letter>(r.below(3)));
      if (trial % 4 == 0 && !w.empty()) w = power(w, 2 + r.below(3));
      std::size_t best = 0;
      for (std::size_t q = 1; q < w.size(); q++)
        if (rotation(w, q) < rotation(w, best)) best = q;
      least = least && least_rotation(w) == best;
      const Word v = trial % 2 ? rotation(w, r.below(w.size() + 1)) : random_word(r, w.size(), 2);
      const bool brute_conj = w.size() == v.size() && (w.empty() || !find_all(v, concat(w, w)).empty());
      conj = conj && are_conjugate(w, v) == brute_conj;
      if (is_primitive(w)) lyndon = lyndon && is_lyndon(lyndon_conjugate(w));
    }
    s.require("Booth's least rotation is the smallest r of the least rotation, against brute force", least);
    s.require("conjugacy by least rotations agrees with v being a factor of ww", conj);
    s.require("the least rotation of a primitive word is a Lyndon word", lyndon);
  }

  return s.done();
}
