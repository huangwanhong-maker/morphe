// Checks of the indexes of words: bit vectors' rank and select against counting, at densities from
// sparse to full and lengths across block boundaries; the wavelet matrix's access, rank and select
// against the sequence itself, alphabets 1 to 300; the FM-index's counts, positions and extracted
// factors against brute force, alphabets 1 to 256; the Aho–Corasick automaton's matches against
// brute force on dictionaries with shared prefixes, nested patterns and repeats, and its automaton
// against the definition of δ.

#include <morphe/check.hpp>
#include <morphe/formal/aho_corasick.hpp>
#include <morphe/formal/fm_index.hpp>
#include <morphe/formal/succinct.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <tuple>
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
  double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
};

Word random_word(Rng& r, std::size_t n, Letter sigma) {
  Word w(n);
  for (auto& c : w) c = static_cast<Letter>(r.below(sigma));
  return w;
}

std::vector<std::size_t> brute_find(const Word& p, const Word& t) {
  std::vector<std::size_t> out;
  if (p.size() > t.size()) return out;
  for (std::size_t i = 0; i + p.size() <= t.size(); i++)
    if (std::equal(p.begin(), p.end(), t.begin() + static_cast<std::ptrdiff_t>(i))) out.push_back(i);
  return out;
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{0x494e444558000001ULL};

  // ---- bit vectors ----------------------------------------------------------------------------------
  {
    bool rank = true, select = true, access = true;
    std::size_t bits = 0;
    for (std::size_t n : {0, 1, 63, 64, 65, 511, 512, 513, 1000, 4096, 100000}) {
      for (double density : {0.0, 0.01, 0.5, 0.97, 1.0}) {
        std::vector<bool> b(n);
        for (std::size_t i = 0; i < n; i++) b[i] = r.uniform() < density;
        const BitVector v(b);
        std::size_t ones = 0, k1 = 0, k0 = 0;
        for (std::size_t i = 0; i <= n; i++) {
          rank = rank && v.rank1(i) == ones && v.rank0(i) == i - ones;
          if (i == n) break;
          access = access && v[i] == b[i];
          if (b[i]) select = select && v.select1(k1++) == i;
          else select = select && v.select0(k0++) == i;
          ones += b[i];
        }
        select = select && v.select1(k1) == n && v.select0(k0) == n && v.ones() == ones;
        bits += n;
      }
    }
    s.require("rank₁ and rank₀ of bit vectors equal the counts of ones and zeros before each position", rank, std::to_string(bits) + " bits");
    s.require("select₁ and select₀ find each one and each zero; beyond the last they return n", select);
    s.require("the bits read back", access);
  }

  // ---- the wavelet matrix ----------------------------------------------------------------------------------
  {
    bool access = true, rank = true, select = true, ar = true;
    for (std::size_t sigma : {1, 2, 3, 4, 5, 20, 64, 255, 256, 300}) {
      for (std::size_t n : {0, 1, 100, 3000}) {
        const Word x = random_word(r, n, static_cast<Letter>(sigma));
        const WaveletMatrix wm(x, sigma);
        std::vector<std::size_t> count(sigma, 0);
        for (std::size_t i = 0; i < n; i++) {
          access = access && wm.access(i) == x[i];
          const auto [c, k] = wm.access_rank(i);
          ar = ar && c == x[i] && k == count[x[i]];
          if (i % 7 == 0)
            for (std::size_t c2 = 0; c2 < sigma; c2 += 1 + sigma / 9) rank = rank && wm.rank(static_cast<Letter>(c2), i) == count[c2];
          select = select && wm.select(x[i], count[x[i]]) == i;
          count[x[i]]++;
        }
        for (std::size_t c = 0; c < sigma; c++) {
          rank = rank && wm.rank(static_cast<Letter>(c), n) == count[c] && wm.count(static_cast<Letter>(c)) == count[c];
          select = select && wm.select(static_cast<Letter>(c), count[c]) == n;
        }
        rank = rank && wm.rank(static_cast<Letter>(sigma), n) == 0;
      }
    }
    s.require("the wavelet matrix: access returns every value", access);
    s.require("the wavelet matrix: rank(c, i) counts c before i, alphabets 1 to 300", rank);
    s.require("the wavelet matrix: select(c, k) finds the k-th c, and none beyond the last", select);
    s.require("the wavelet matrix: access and rank in one pass", ar);
  }

  // ---- the FM-index ----------------------------------------------------------------------------------
  {
    bool counts = true, positions = true, extract = true, absent = true;
    std::size_t patterns = 0;
    for (std::size_t sigma : {1, 2, 4, 20, 256}) {
      for (std::size_t n : {0, 1, 50, 2000, 20000}) {
        for (std::size_t sample : {1, 7, 32}) {
          Word w = random_word(r, n, static_cast<Letter>(sigma));
          if (n >= 2000 && sigma > 1) {
            // Repeats, so that patterns occur many times.
            const Word unit = random_word(r, 50, static_cast<Letter>(sigma));
            for (std::size_t k = 0; k < 10; k++) std::copy(unit.begin(), unit.end(), w.begin() + static_cast<std::ptrdiff_t>(r.below(n - 50)));
          }
          const FmIndex<> fm(w, sigma, sample);
          for (int q = 0; q < 40; q++) {
            Word p;
            if (q % 2 && n > 0) {
              const std::size_t b = r.below(n), len = 1 + r.below(std::min<std::size_t>(12, n - b));
              p.assign(w.begin() + static_cast<std::ptrdiff_t>(b), w.begin() + static_cast<std::ptrdiff_t>(b + len));
            } else {
              p = random_word(r, 1 + r.below(4), static_cast<Letter>(sigma));
            }
            const auto want = brute_find(p, w);
            counts = counts && fm.count(p) == want.size();
            positions = positions && fm.locate(p) == want;
            patterns++;
          }
          if (n > 0)
            for (int q = 0; q < 20; q++) {
              const std::size_t b = r.below(n + 1), e = b + r.below(n - b + 1);
              extract = extract && fm.extract(b, e) == Word(w.begin() + static_cast<std::ptrdiff_t>(b), w.begin() + static_cast<std::ptrdiff_t>(e));
            }
          absent = absent && fm.count(Word{static_cast<Letter>(sigma)}) == 0 && fm.count(Word{}) == n + 1;
        }
      }
    }
    s.require("FM-index counts equal brute force (alphabets 1 to 256, sample rates 1, 7, 32)", counts, std::to_string(patterns) + " patterns");
    s.require("FM-index positions by the sampled suffix array equal brute force", positions);
    s.require("FM-index extracts every factor asked for", extract);
    s.require("a letter outside the alphabet occurs nowhere; the empty pattern matches all n + 1 rows", absent);

    const Word big = random_word(r, 1'000'000, 4);
    const auto t0 = std::chrono::steady_clock::now();
    const FmIndex<> fm(big, 4, 32);
    const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::size_t found = 0;
    const auto t1 = std::chrono::steady_clock::now();
    for (int q = 0; q < 10000; q++) {
      const std::size_t b = r.below(big.size() - 20);
      found += fm.count(Word(big.begin() + static_cast<std::ptrdiff_t>(b), big.begin() + static_cast<std::ptrdiff_t>(b + 20))) >= 1;
    }
    const double query = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    char buf[160];
    std::snprintf(buf, sizeof buf, "build %.2f s, %.1f µs a count of 20 letters, %.2f bytes a letter", build, query / 1e4 * 1e6,
                  static_cast<double>(fm.bytes()) / 1e6);
    s.require("an FM-index of 10⁶ letters over four finds every factor it is asked for", found == 10000, buf);
  }

  // ---- Aho–Corasick ----------------------------------------------------------------------------------
  {
    // Aho and Corasick 1975, Fig. 1: {he, she, his, hers} in "ushers" ends at 4 (she, he) and 6 (hers).
    const std::vector<Word> dict{word("he"), word("she"), word("his"), word("hers")};
    const AhoCorasick ac(dict, 26);
    const auto m = ac.matches(word("ushers"));
    using M = AhoCorasick::Match;
    s.require("Aho and Corasick's example: {he, she, his, hers} in ushers", m == std::vector<M>{{1, 4, 1}, {2, 4, 0}, {2, 6, 3}},
              std::to_string(m.size()) + " matches, " + std::to_string(ac.states()) + " states");
    s.equal("the trie of {he, she, his, hers} has 10 states (their Fig. 1)", ac.states(), std::size_t{10});

    bool brute = true, automaton = true;
    std::size_t total = 0;
    for (int trial = 0; trial < 300; trial++) {
      const Letter sigma = 1 + static_cast<Letter>(r.below(trial % 2 ? 3 : 26));
      std::vector<Word> patterns;
      const std::size_t count = 1 + r.below(30);
      const Word text = random_word(r, r.below(500), sigma);
      for (std::size_t p = 0; p < count; p++) {
        if (p > 0 && r.below(4) == 0) {
          // A prefix, a suffix or a repeat of an earlier pattern.
          const Word& q = patterns[r.below(patterns.size())];
          const std::size_t len = 1 + r.below(q.size());
          patterns.push_back(r.below(2) ? Word(q.begin(), q.begin() + static_cast<std::ptrdiff_t>(len)) : Word(q.end() - static_cast<std::ptrdiff_t>(len), q.end()));
        } else if (!text.empty() && r.below(2)) {
          const std::size_t b = r.below(text.size()), len = 1 + r.below(std::min<std::size_t>(8, text.size() - b));
          patterns.emplace_back(text.begin() + static_cast<std::ptrdiff_t>(b), text.begin() + static_cast<std::ptrdiff_t>(b + len));
        } else {
          patterns.push_back(random_word(r, 1 + r.below(6), sigma));
        }
      }
      const AhoCorasick a(patterns, sigma);
      auto got = a.matches(text);
      std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> g, want;
      for (const auto& x : got) g.emplace_back(x.begin, x.end, x.pattern);
      for (std::size_t p = 0; p < patterns.size(); p++)
        for (std::size_t b : brute_find(patterns[p], text)) want.emplace_back(b, b + patterns[p].size(), p);
      std::sort(g.begin(), g.end());
      std::sort(want.begin(), want.end());
      brute = brute && g == want;
      total += want.size();
      // δ(u, a) is the longest suffix of (word of u)·a that is a state; check by its depth.
      for (std::uint32_t u = 0; u < a.states() && automaton; u++)
        for (Letter c = 0; c < sigma; c++) {
          const auto v = a.delta(u, c);
          const auto g2 = a.go(u, c);
          automaton = automaton && (g2 == AhoCorasick::NONE ? a.depth(v) <= a.depth(u) : v == g2);
        }
    }
    s.require("Aho–Corasick reports every occurrence of every pattern, against brute force (300 dictionaries)", brute, std::to_string(total) + " occurrences");
    s.require("the automaton follows the trie where it has an edge, and otherwise falls to no deeper state", automaton);
  }

  return s.done();
}
