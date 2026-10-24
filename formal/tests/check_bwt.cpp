// Checks of the Burrows–Wheeler transforms: the BWT of w$ from the suffix array against the sorted
// rotations, and inverted by LF; the rotation BWT against its definition and inverted from the
// primary row on every binary word up to length 12; the bijective BWT a bijection on every binary
// word of length ≤ 12 and every ternary word of length ≤ 7 (every word reached once, and inverted);
// the Gessel–Reutenauer bijection: every word is the extended BWT of the primitive necklaces its
// standard permutation's cycles spell; the extended BWT of multisets of primitive words inverted;
// the BWT of every conjugate of a standard Sturmian word of the form b^p a^q, and conversely
// (Mantaci, Restivo and Sciortino 2003); runs and move-to-front.

#include <morphe/check.hpp>
#include <morphe/formal/bwt.hpp>

#include <algorithm>
#include <cstdint>
#include <numeric>
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

Word base_word(std::uint64_t x, std::size_t n, Letter k) {
  Word w(n);
  for (std::size_t i = 0; i < n; i++, x /= k) w[i] = static_cast<Letter>(x % k);
  return w;
}

std::uint64_t ipow(std::uint64_t k, std::size_t n) {
  std::uint64_t t = 1;
  for (std::size_t i = 0; i < n; i++) t *= k;
  return t;
}

// The last letters of the rotations of w sorted lexicographically, ties by offset.
RotationBwt brute_rotation_bwt(const Word& w) {
  std::vector<std::size_t> order(w.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rotation(w, a) < rotation(w, b); });
  RotationBwt b;
  for (std::size_t i = 0; i < order.size(); i++) {
    b.last.push_back(w[(order[i] + w.size() - 1) % w.size()]);
    if (order[i] == 0) b.primary = i;
  }
  return b;
}

// The standard words over {0, 1} of length at most n: s_{−1} = 1, s_0 = 0,
// s_{k+1} = s_k^{d_k} s_{k−1} with d_0 ≥ 0 and d_k ≥ 1 after (Lothaire 2002, §2.2.2).
void standard_words(const Word& prev, const Word& cur, std::size_t k, std::size_t n, std::set<Word>& out) {
  out.insert(cur);
  for (std::size_t d = k == 0 ? 0 : 1;; d++) {
    Word next = power(cur, d);
    next.insert(next.end(), prev.begin(), prev.end());
    if (next.size() > n) break;
    standard_words(cur, next, k + 1, n, out);
  }
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{0x4257540000000001ULL};

  // ---- the BWT of w$ ----------------------------------------------------------------------------------
  {
    bool definition = true, inverse = true;
    for (int trial = 0; trial < 600; trial++) {
      const Letter sigma = 1 + static_cast<Letter>(r.below(trial % 3 == 0 ? 2 : 30));
      const Word w = random_word(r, r.below(trial < 500 ? 60 : 20000), sigma);
      const Word L = bwt(w);
      if (w.size() <= 60) {
        // w$ shifted, its rotations sorted.
        Word x(w.size() + 1, 0);
        for (std::size_t i = 0; i < w.size(); i++) x[i] = w[i] + 1;
        definition = definition && brute_rotation_bwt(x).last == L;
      }
      inverse = inverse && inverse_bwt(L) == w;
    }
    s.require("the BWT of w$ from SA-IS equals the last column of the sorted rotations of w$", definition);
    s.require("the BWT of w$ is inverted by LF (600 words, up to 2·10⁴ letters)", inverse);
    const Word fib = fixed_point_prefix(fibonacci_morphism(), 0, 50000);
    s.require("the BWT of a Fibonacci prefix of 5·10⁴ letters is inverted", inverse_bwt(bwt(fib)) == fib);
    s.require("the BWT of the empty word is the sentinel alone", bwt(Word{}) == Word{0} && inverse_bwt(Word{0}).empty());
    // Burrows and Wheeler 1994, §2: "abraca" gives L = caraab with I = 1 (rotations, no sentinel).
    const auto bw = rotation_bwt(word("abraca"));
    s.require("Burrows and Wheeler's example: abraca ↦ (caraab, I = 1)", text(bw.last) == "caraab" && bw.primary == 1, text(bw.last) + ", I = " + std::to_string(bw.primary));
  }

  // ---- the rotation BWT ----------------------------------------------------------------------------------
  {
    bool definition = true, inverse = true;
    std::size_t words = 0;
    for (std::size_t n = 1; n <= 12; n++)
      for (std::uint64_t x = 0; x < (1ULL << n); x++) {
        const Word w = base_word(x, n, 2);
        const auto b = rotation_bwt(w);
        const auto brute = brute_rotation_bwt(w);
        definition = definition && b.last == brute.last && b.primary == brute.primary;
        inverse = inverse && inverse_rotation_bwt(b) == w;
        words++;
      }
    for (int trial = 0; trial < 300; trial++) {
      Word w = random_word(r, 1 + r.below(40), 1 + static_cast<Letter>(r.below(5)));
      if (trial % 3 == 0) w = power(w, 2 + r.below(3));
      const auto b = rotation_bwt(w);
      const auto brute = brute_rotation_bwt(w);
      definition = definition && b.last == brute.last && b.primary == brute.primary;
      inverse = inverse && inverse_rotation_bwt(b) == w;
    }
    s.require("the rotation BWT by ω-order doubling equals the sorted rotations (ties by offset)", definition, std::to_string(words + 300) + " words, powers among them");
    s.require("the rotation BWT is inverted from its primary row", inverse);
  }

  // ---- the bijective BWT ----------------------------------------------------------------------------------
  {
    bool bijection = true, inverse = true, lyndon = true;
    std::string detail;
    for (Letter k : {2u, 3u}) {
      const std::size_t nmax = k == 2 ? 12 : 7;
      for (std::size_t n = 1; n <= nmax; n++) {
        const std::uint64_t total = ipow(k, n);
        std::set<Word> images;
        for (std::uint64_t x = 0; x < total; x++) {
          const Word w = base_word(x, n, k);
          const Word b = bijective_bwt(w);
          images.insert(b);
          inverse = inverse && b.size() == n && inverse_bijective_bwt(b) == w;
          if (is_lyndon(w)) lyndon = lyndon && b == rotation_bwt(w).last;
        }
        bijection = bijection && images.size() == total;
      }
      detail += (k == 2 ? "" : "; ") + std::to_string(k) + " letters to length " + std::to_string(nmax);
    }
    s.require("the bijective BWT is a bijection: every word of length n reached once", bijection, detail);
    s.require("the bijective BWT is inverted by the cycles of LF read as Lyndon words", inverse);
    s.require("on a Lyndon word the bijective BWT is the rotation BWT", lyndon);
    bool random = true;
    for (int trial = 0; trial < 200; trial++) {
      const Word w = random_word(r, 1 + r.below(3000), 2 + static_cast<Letter>(r.below(20)));
      random = random && inverse_bijective_bwt(bijective_bwt(w)) == w;
    }
    s.require("the bijective BWT is inverted on random words up to 3000 letters", random);
  }

  // ---- the extended BWT and the Gessel–Reutenauer bijection ----------------------------------------------------------------------------------
  {
    bool gr = true, primitive = true;
    std::size_t words = 0;
    for (Letter k : {2u, 3u}) {
      const std::size_t nmax = k == 2 ? 10 : 6;
      for (std::size_t n = 1; n <= nmax; n++)
        for (std::uint64_t x = 0; x < ipow(k, n); x++) {
          const Word L = base_word(x, n, k);
          const auto neck = necklaces(L);
          for (const auto& c : neck) primitive = primitive && is_lyndon(c);
          gr = gr && extended_bwt(neck).last == L;
          words++;
        }
    }
    s.require("the cycles of the standard permutation of a word are primitive necklaces (Lyndon words)", primitive);
    s.require("Gessel–Reutenauer: every word is the extended BWT of its necklaces", gr, std::to_string(words) + " words");

    bool inverse = true;
    for (int trial = 0; trial < 300; trial++) {
      std::vector<Word> words_in;
      const std::size_t m = 1 + r.below(8);
      for (std::size_t j = 0; j < m; j++) {
        Word w;
        do w = random_word(r, 1 + r.below(12), 2 + static_cast<Letter>(r.below(2)));
        while (!is_primitive(w));
        words_in.push_back(w);
        if (r.below(4) == 0) words_in.push_back(w);  // a repeated word
      }
      inverse = inverse && inverse_extended_bwt(extended_bwt(words_in)) == words_in;
    }
    s.require("the extended BWT of a multiset of primitive words (repeats among them) is inverted from the words' rows", inverse);
  }

  // ---- standard Sturmian words ----------------------------------------------------------------------------------
  {
    // Mantaci, Restivo and Sciortino 2003 (Inf. Process. Lett. 86, 241–246): the
    // BWT of w over a < b is b^p a^q with gcd(p, q) = 1 exactly when w is a conjugate of a
    // standard word.
    const auto clustered = [](const Word& L) {
      std::size_t i = 0;
      while (i < L.size() && L[i] == 1) i++;
      const std::size_t p = i;
      while (i < L.size() && L[i] == 0) i++;
      return i == L.size() && p > 0 && p < L.size();
    };
    bool fibonacci = true;
    for (std::size_t m = 2; m <= 20; m++) fibonacci = fibonacci && clustered(rotation_bwt(fibonacci_word(m)).last) && runs(rotation_bwt(fibonacci_word(m)).last) == 2;
    s.require("the BWT of every finite Fibonacci word f_2 … f_20 (to 17 711 letters) is b^p a^q", fibonacci);

    std::set<Word> standard;
    standard_words(Word{1}, Word{0}, 0, 12, standard);
    std::set<Word> conjugates;
    for (const auto& w : standard)
      for (std::size_t q = 0; q < w.size(); q++) conjugates.insert(rotation(w, q));
    bool forward = true, converse = true;
    std::size_t clustered_count = 0;
    for (std::size_t n = 2; n <= 12; n++)
      for (std::uint64_t x = 0; x < (1ULL << n); x++) {
        const Word w = base_word(x, n, 2);
        const std::size_t ones = static_cast<std::size_t>(std::count(w.begin(), w.end(), Letter{1}));
        if (ones == 0 || ones == n) continue;
        const bool c = clustered(rotation_bwt(w).last) && std::gcd(ones, n - ones) == 1;
        clustered_count += c;
        const bool conj = conjugates.count(w) > 0;
        if (conj) forward = forward && c;
        else converse = converse && !c;
      }
    s.require("every conjugate of a standard word of length ≤ 12 has a BWT b^p a^q with gcd(p, q) = 1", forward,
              std::to_string(standard.size()) + " standard words");
    s.require("conversely, every binary word of length ≤ 12 with such a BWT is a conjugate of a standard word", converse,
              std::to_string(clustered_count) + " words");
  }

  // ---- runs and move-to-front ----------------------------------------------------------------------------------
  {
    bool mtf = true;
    for (int trial = 0; trial < 300; trial++) {
      const Letter sigma = 1 + static_cast<Letter>(r.below(20));
      const Word w = random_word(r, r.below(500), sigma);
      mtf = mtf && inverse_move_to_front(move_to_front(w, sigma), sigma) == w;
    }
    s.require("move-to-front is inverted", mtf);
    s.equal("move-to-front of bananaaa over a, b, n", move_to_front(Word{1, 0, 2, 0, 2, 0, 0, 0}, 3), Word{1, 1, 2, 1, 1, 1, 0, 0});
    s.equal("runs of aabbbac", runs(word("aabbbac")), std::size_t{4});
    const Word fib = fixed_point_prefix(fibonacci_morphism(), 0, 50000), rnd = random_word(r, 50000, 2);
    const std::size_t r_fib = runs(bwt(fib)), r_rnd = runs(bwt(rnd));
    s.require("the BWT of a repetitive word has far fewer runs than that of a random one", r_fib * 50 < r_rnd,
              "Fibonacci prefix r = " + std::to_string(r_fib) + ", random r = " + std::to_string(r_rnd));
  }

  return s.done();
}
