// Checks of the run-length FM-index and the r-index: rank, access and LF by runs against the
// plain BWT at every row and letter, and the word read back by LF; counts, intervals and the
// sorted occurrence lists against the FM-index and brute force, on random words (alphabets 1 to
// 256), on highly repetitive words (a word, then copies of it with a few substitutions each) and
// on edge cases; the toehold against the suffix array, and φ(SA$[i]) = SA$[i − 1] at every row;
// r equal to the runs of the BWT, the bound r ≤ 2p + 1 on powers of primitive words of length p,
// and r and the space reported, against the full suffix array, on repetitive words.

#include <morphe/check.hpp>
#include <morphe/formal/bwt.hpp>
#include <morphe/formal/fm_index.hpp>
#include <morphe/formal/r_index.hpp>
#include <morphe/formal/suffix_array.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <utility>
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

// A word of length `base`, then `copies` copies of it, each with `changes` substitutions.
Word repetitive_word(Rng& r, std::size_t base, std::size_t copies, std::size_t changes, Letter sigma) {
  const Word w = random_word(r, base, sigma);
  Word t = w;
  for (std::size_t c = 0; c < copies; c++) {
    Word u = w;
    if (sigma > 1)
      for (std::size_t k = 0; k < changes; k++) {
        const std::size_t at = r.below(base);
        u[at] = static_cast<Letter>((u[at] + 1 + r.below(sigma - 1)) % sigma);
      }
    t.insert(t.end(), u.begin(), u.end());
  }
  return t;
}

std::vector<std::size_t> brute_find(const Word& p, const Word& t) {
  std::vector<std::size_t> out;
  if (p.size() > t.size()) return out;
  for (std::size_t i = 0; i + p.size() <= t.size(); i++)
    if (std::equal(p.begin(), p.end(), t.begin() + static_cast<std::ptrdiff_t>(i))) out.push_back(i);
  return out;
}

// A pattern: a factor of the text (perhaps with one letter changed), or a random word.
Word pattern(Rng& r, const Word& w, Letter sigma, std::size_t longest) {
  Word p;
  if (!w.empty() && r.below(3) != 0) {
    const std::size_t b = r.below(w.size()), len = 1 + r.below(std::min<std::size_t>(longest, w.size() - b));
    p.assign(w.begin() + static_cast<std::ptrdiff_t>(b), w.begin() + static_cast<std::ptrdiff_t>(b + len));
    if (r.below(4) == 0) p[r.below(p.size())] = static_cast<Letter>(r.below(sigma));
  } else {
    p = random_word(r, 1 + r.below(5), sigma);
  }
  return p;
}

// SA$ of w: row 0 is the suffix $ at n.
std::vector<std::size_t> sa_dollar(const Word& w) {
  const auto sa = suffix_array(w);
  std::vector<std::size_t> out{w.size()};
  out.insert(out.end(), sa.begin(), sa.end());
  return out;
}

struct Tally {
  bool rank = true, access = true, lf = true, runs = true, inverse = true;
  bool counts = true, positions = true, toehold = true, phi = true, samples = true;
  std::size_t texts = 0, patterns = 0, occurrences = 0;
};

// Every check of one text.
void check_text(Tally& t, Rng& r, const Word& w, std::size_t sigma, std::size_t queries, bool exhaustive, bool brute) {
  const RIndex<> ri(w, sigma);
  const auto& rl = ri.counting_index();
  const Word L = bwt(w);
  const std::size_t n = w.size();
  t.texts++;
  t.runs = t.runs && rl.runs() == runs(L) && ri.runs() == runs(L);
  t.samples = t.samples && ri.samples() == 2 * ri.runs() - 1;
  if (exhaustive) {
    std::vector<std::size_t> count(sigma + 1, 0);
    for (std::size_t i = 0; i <= n + 1; i++) {
      for (Letter c = 0; c <= sigma; c++) t.rank = t.rank && rl.rank(c, i) == count[c];
      t.rank = t.rank && rl.rank(static_cast<Letter>(sigma + 1), i) == 0;
      if (i == n + 1) break;
      t.access = t.access && rl.access(i) == L[i];
      count[L[i]]++;
    }
    const auto lf = lf_map(L);
    for (std::size_t i = 0; i <= n; i++) t.lf = t.lf && rl.lf(i) == lf[i];
  }
  // The word read back by LF from the row of $.
  {
    Word back(n);
    std::size_t row = 0;
    for (std::size_t k = 0; k < n; k++) {
      back[n - 1 - k] = rl.access(row) - 1;
      row = rl.lf(row);
    }
    t.inverse = t.inverse && back == w;
  }
  const auto sa = sa_dollar(w);
  for (std::size_t i = 1; i <= n; i++) t.phi = t.phi && ri.phi(sa[i]) == sa[i - 1];
  const FmIndex<> fm(w, sigma, 16);
  for (std::size_t q = 0; q < queries; q++) {
    const Word p = q == 0 ? Word{} : pattern(r, w, static_cast<Letter>(sigma), 24);
    const auto got = ri.locate(p);
    const auto iv = ri.interval(p);
    t.counts = t.counts && ri.count(p) == fm.count(p) && iv == fm.interval(p) && got.size() == ri.count(p);
    t.positions = t.positions && got == fm.locate(p);
    if (brute && !p.empty()) {
      const auto want = brute_find(p, w);
      t.counts = t.counts && ri.count(p) == want.size();
      t.positions = t.positions && got == want;
    }
    const auto h = ri.toehold(p);
    t.toehold = t.toehold && h.begin == iv.first && h.end == iv.second && (h.begin == h.end || h.last == sa[h.end - 1]);
    t.patterns++;
    t.occurrences += got.size();
  }
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{0x52494e4445580001ULL};

  // ---- random words, repetitive words and edge cases ----------------------------------------------------------------------------------
  {
    Tally t;
    for (std::size_t sigma : {1, 2, 3, 4, 20, 256})
      for (std::size_t n : {0, 1, 2, 3, 17, 300, 3000}) check_text(t, r, random_word(r, n, static_cast<Letter>(sigma)), sigma, 40, n <= 300, true);
    s.require("random words: rank by runs equals the count in the plain BWT at every row and letter", t.rank, std::to_string(t.texts) + " words");
    s.require("random words: access by runs returns every letter of the BWT, and LF equals the LF map", t.access && t.lf);
    s.require("random words: r equals the runs of the BWT; 2r − 1 values of SA$ are kept", t.runs && t.samples);
    s.require("random words: LF from the row of $ reads the word back", t.inverse);
    s.require("random words: counts and intervals equal the FM-index's and brute force (alphabets 1 to 256)", t.counts,
              std::to_string(t.patterns) + " patterns");
    s.require("random words: the sorted occurrences by the toehold and φ equal the FM-index's and brute force", t.positions,
              std::to_string(t.occurrences) + " occurrences");
    s.require("random words: the toehold is SA$ at the last row of the pattern's interval", t.toehold);
    s.require("random words: φ(SA$[i]) = SA$[i − 1] at every row i ≥ 1", t.phi);

    Tally u;
    for (std::size_t sigma : {2, 4, 20})
      for (std::size_t base : {1, 7, 100, 1000})
        for (std::size_t changes : {0, 1, 3})
          check_text(u, r, repetitive_word(r, base, 9, changes, static_cast<Letter>(sigma)), sigma, 40, base <= 100, true);
    s.require("repetitive words (a word and nine copies with 0 to 3 substitutions each): rank, access and LF by runs, and the word read back",
              u.rank && u.access && u.lf && u.inverse && u.runs && u.samples, std::to_string(u.texts) + " words");
    s.require("repetitive words: counts, sorted occurrences and toeholds equal the FM-index's, brute force and the suffix array",
              u.counts && u.positions && u.toehold, std::to_string(u.patterns) + " patterns, " + std::to_string(u.occurrences) + " occurrences");
    s.require("repetitive words: φ(SA$[i]) = SA$[i − 1] at every row i ≥ 1", u.phi);

    // Edge cases: the empty word, unary words, a single letter, patterns longer than the word,
    // letters outside the alphabet, and the empty pattern.
    bool edge = true;
    {
      const RIndex<> e(Word{}, 4);
      edge = edge && e.size() == 0 && e.runs() == 1 && e.count(Word{}) == 1 && e.locate(Word{}) == std::vector<std::size_t>{0} &&
             e.count(Word{0}) == 0 && e.locate(Word{2}).empty();
      Tally v;
      for (std::size_t n : {1, 2, 5, 64, 1000}) check_text(v, r, Word(n, 0), 1, 20, true, true);
      for (std::size_t n : {1, 2, 5, 64, 1000}) check_text(v, r, Word(n, 3), 4, 20, true, true);
      edge = edge && v.rank && v.access && v.lf && v.inverse && v.runs && v.counts && v.positions && v.toehold && v.phi;
      const Word w = random_word(r, 50, 3);
      const RIndex<> ri(w, 3);
      Word longer = w;
      longer.push_back(0);
      std::vector<std::size_t> all(51);
      std::iota(all.begin(), all.end(), std::size_t{0});
      edge = edge && ri.count(longer) == 0 && ri.locate(longer).empty() && ri.count(Word{3}) == 0 && ri.locate(Word{0, 7, 1}).empty() &&
             ri.count(Word{}) == 51 && ri.locate(Word{}) == all && ri.count(w) == 1 && ri.locate(w) == std::vector<std::size_t>{0};
      // The BWT of a^n $ is a^n $ (the rows $, a$, …, a^n$): r = 2 for n ≥ 1.
      edge = edge && RIndex<>(Word(1000, 0), 1).runs() == 2;
    }
    s.require("edge cases: the empty word, unary words, patterns longer than the word, letters outside the alphabet, the empty pattern", edge);
  }

  // ---- r on powers of primitive words ----------------------------------------------------------------------------------
  {
    bool bound = true, found = true;
    std::string detail;
    for (std::size_t p : {1, 2, 5, 37, 100})
      for (std::size_t k : {1, 2, 3, 10, 200}) {
        Word w;
        do w = random_word(r, p, 4);
        while (!is_primitive(w));
        const Word t = power(w, k);
        const RIndex<> ri(t, 4);
        bound = bound && ri.runs() <= 2 * p + 1;
        // Every occurrence of w in w^k: k of them, at the multiples of p.
        std::vector<std::size_t> want(k);
        for (std::size_t j = 0; j < k; j++) want[j] = j * p;
        found = found && ri.locate(w) == want;
        if (k == 200) detail += (detail.empty() ? "" : ", ") + std::string("p = ") + std::to_string(p) + ": r = " + std::to_string(ri.runs());
      }
    s.require("powers w^k of primitive words of length p: r ≤ 2p + 1 whatever k (k up to 200)", bound, detail + " at k = 200");
    s.require("powers w^k: w occurs exactly at the k multiples of p", found);
  }

  // ---- a large repetitive word: r, the space, and agreement with the FM-index ----------------------------------------------------------------------------------
  {
    const std::size_t base = 10000, copies = 99, changes = 10;
    const Word t = repetitive_word(r, base, copies, changes, 4);
    const auto t0 = std::chrono::steady_clock::now();
    const RIndex<> ri(t, 4);
    const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const FmIndex<> fm(t, 4, 32);
    const Word noise = random_word(r, t.size(), 4);
    const std::size_t r_noise = runs(bwt(noise));
    // Factors of 12 to 30 letters, a quarter with one letter changed: each occurs about once a copy.
    std::vector<Word> patterns;
    for (int q = 0; q < 2000; q++) {
      const std::size_t b = r.below(t.size() - 30), len = 12 + r.below(19);
      Word p(t.begin() + static_cast<std::ptrdiff_t>(b), t.begin() + static_cast<std::ptrdiff_t>(b + len));
      if (q % 4 == 0) p[r.below(len)] = static_cast<Letter>(r.below(4));
      patterns.push_back(std::move(p));
    }
    std::vector<std::vector<std::size_t>> found;
    const auto t1 = std::chrono::steady_clock::now();
    for (const auto& p : patterns) found.push_back(ri.locate(p));
    const double query = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    bool agree = true;
    std::size_t occ = 0;
    for (std::size_t q = 0; q < patterns.size(); q++) {
      agree = agree && ri.count(patterns[q]) == fm.count(patterns[q]) && found[q] == fm.locate(patterns[q]);
      occ += found[q].size();
    }
    const double sa_bytes = static_cast<double>(t.size() + 1) * sizeof(std::uint32_t);
    char buf[400];
    std::snprintf(buf, sizeof buf,
                  "n = %zu, r = %zu (r/n = %.4f; a random word of the same length has r = %zu); r-index %.0f kB, %.4f of the suffix array's "
                  "%.0f kB, FM-index (sample 32) %.0f kB; built in %.2f s",
                  t.size(), ri.runs(), static_cast<double>(ri.runs()) / static_cast<double>(t.size()), r_noise, ri.bytes() / 1e3,
                  ri.bytes() / sa_bytes, sa_bytes / 1e3, fm.bytes() / 1e3, build);
    s.require("a word of 10⁴ letters and 99 copies with 10 substitutions each: r reported, with the space against the suffix array",
              ri.runs() == runs(bwt(t)) && ri.runs() <= t.size() + 1, buf);
    std::snprintf(buf, sizeof buf, "%zu occurrences; the r-index locates a pattern and its occurrences in %.1f µs", occ, query / 2000 * 1e6);
    s.require("the large repetitive word: counts and sorted occurrences of 2000 factors of 12 to 30 letters equal the FM-index's", agree, buf);
  }

  return s.done();
}
