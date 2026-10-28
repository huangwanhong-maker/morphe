// Checks of the positional Burrows–Wheeler transform: a_k sorts the reversed prefixes as a brute
// sort does, and d_k holds the true divergences, at every position of random and structured panels
// (alphabets 1 to 4, identical words, words of length 0 and 1); the w-map against counting; the
// set-maximal matches of every word, and the long matches of length ≥ L, against brute force over
// every pair of words; a new word's set-maximal matches and longest matches against brute force;
// and a panel built as a mosaic of a few founders with a few substitutions, whose shared segments
// the long matches recover, and a new mosaic of the panel's words whose copied segments its
// matches recover.

#include <morphe/check.hpp>
#include <morphe/formal/pbwt.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <numeric>
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

using Panel = std::vector<Word>;

Word random_word(Rng& r, std::size_t n, Letter sigma) {
  Word w(n);
  for (auto& c : w) c = static_cast<Letter>(r.below(sigma));
  return w;
}

void mutate(Rng& r, Word& w, std::size_t changes, Letter sigma, std::vector<bool>* changed = nullptr) {
  if (sigma < 2 || w.empty()) return;
  for (std::size_t k = 0; k < changes; k++) {
    const std::size_t at = r.below(w.size());
    w[at] = static_cast<Letter>((w[at] + 1 + r.below(sigma - 1)) % sigma);
    if (changed) (*changed)[at] = true;
  }
}

// A panel: random words, or copies of a few words with substitutions, sometimes with repeats.
Panel random_panel(Rng& r, std::size_t M, std::size_t N, Letter sigma, bool structured) {
  Panel p;
  if (!structured) {
    for (std::size_t i = 0; i < M; i++) p.push_back(random_word(r, N, sigma));
    return p;
  }
  const std::size_t F = 1 + r.below(3);
  Panel f;
  for (std::size_t i = 0; i < F; i++) f.push_back(random_word(r, N, sigma));
  for (std::size_t i = 0; i < M; i++) {
    Word w = f[r.below(F)];
    if (r.below(4) != 0) mutate(r, w, 1 + r.below(1 + N / 8), sigma);
    p.push_back(w);
  }
  return p;
}

// The maximal matches [s, e) of two words.
template <class F>
void maximal_matches(const Word& x, const Word& y, F&& f) {
  const std::size_t N = x.size();
  for (std::size_t k = 0; k < N;) {
    if (x[k] != y[k]) {
      k++;
      continue;
    }
    std::size_t e = k;
    while (e < N && x[e] == y[e]) e++;
    f(k, e);
    k = e;
  }
}

// The set-maximal matches, by definition: maximal matches not strictly inside another of the same word.
struct Candidate {
  std::size_t other, begin, end;
};
std::vector<Candidate> set_maximal(const std::vector<Candidate>& all) {
  std::vector<Candidate> out;
  for (const auto& m : all) {
    bool inside = false;
    for (const auto& o : all)
      if (o.begin <= m.begin && o.end >= m.end && (o.begin != m.begin || o.end != m.end)) {
        inside = true;
        break;
      }
    if (!inside) out.push_back(m);
  }
  return out;
}

std::vector<PbwtMatch> brute_set_maximal(const Panel& p) {
  std::vector<PbwtMatch> out;
  for (std::size_t i = 0; i < p.size(); i++) {
    std::vector<Candidate> all;
    for (std::size_t j = 0; j < p.size(); j++)
      if (j != i) maximal_matches(p[i], p[j], [&](std::size_t s, std::size_t e) { all.push_back({j, s, e}); });
    for (const auto& m : set_maximal(all)) out.push_back({i, m.other, m.begin, m.end});
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<PbwtMatch> brute_long(const Panel& p, std::size_t L) {
  std::vector<PbwtMatch> out;
  for (std::size_t i = 0; i < p.size(); i++)
    for (std::size_t j = i + 1; j < p.size(); j++)
      maximal_matches(p[i], p[j], [&](std::size_t s, std::size_t e) {
        if (e - s >= L) out.push_back({i, j, s, e});
      });
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<PbwtQueryMatch> brute_query(const Panel& p, const Word& z) {
  std::vector<Candidate> all;
  for (std::size_t j = 0; j < p.size(); j++) maximal_matches(z, p[j], [&](std::size_t s, std::size_t e) { all.push_back({j, s, e}); });
  std::vector<PbwtQueryMatch> out;
  for (const auto& m : set_maximal(all)) out.push_back({m.other, m.begin, m.end});
  std::sort(out.begin(), out.end());
  return out;
}

// e_k: the least s with z[s, k) = x_j[s, k) for some j.
std::vector<std::size_t> brute_longest(const Panel& p, const Word& z) {
  std::vector<std::size_t> e(z.size() + 1);
  for (std::size_t k = 0; k <= z.size(); k++) {
    e[k] = k;
    for (const auto& x : p) {
      std::size_t s = k;
      while (s > 0 && x[s - 1] == z[s - 1]) s--;
      e[k] = std::min(e[k], s);
    }
  }
  return e;
}

template <class T>
std::vector<T> sorted(std::vector<T> v) {
  std::sort(v.begin(), v.end());
  return v;
}

struct Tally {
  bool prefix = true, divergence = true, sweep = true, wmap = true, setmax = true, longm = true, query = true, longest = true;
  std::size_t panels = 0, setmax_count = 0, long_count = 0, queries = 0, query_count = 0;
};

void check_panel(Tally& t, Rng& r, const Panel& p, Letter sigma) {
  const std::size_t M = p.size(), N = M ? p[0].size() : 0;
  const Pbwt<> pb(p, sigma);
  t.panels++;
  // a_k and d_k against brute force, and the stored arrays against the sweep.
  std::vector<std::uint32_t> order(M);
  for (std::size_t k = 0; k <= N; k++) {
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](std::uint32_t i, std::uint32_t j) {
      for (std::size_t q = k; q-- > 0;)
        if (p[i][q] != p[j][q]) return p[i][q] < p[j][q];
      return i < j;
    });
    const auto a = pb.prefix_array(k);
    const auto d = pb.divergence_array(k);
    t.prefix = t.prefix && std::equal(a.begin(), a.end(), order.begin(), order.end());
    for (std::size_t j = 0; j < M; j++) {
      std::size_t s = k;
      if (j > 0)
        while (s > 0 && p[order[j]][s - 1] == p[order[j - 1]][s - 1]) s--;
      t.divergence = t.divergence && d[j] == s;
    }
    if (k < N) {
      // w_k(j, c) = (letters of column k below c) + (c's among y_0..y_{j−1}); the word itself lands on itself.
      const auto a2 = pb.prefix_array(k + 1);
      std::vector<std::size_t> below(sigma + 1, 0), seen(sigma, 0);
      for (std::size_t j = 0; j < M; j++) below[p[a[j]][k] + 1]++;
      for (Letter c = 0; c < sigma; c++) below[c + 1] += below[c];
      for (std::size_t j = 0; j <= M; j++) {
        for (Letter c = 0; c < sigma; c++) t.wmap = t.wmap && pb.extend(k, j, c) == below[c] + seen[c];
        if (j == M) break;
        const Letter c = p[a[j]][k];
        t.wmap = t.wmap && pb.letter(k, j) == c && a2[pb.extend(k, j, c)] == a[j];
        seen[c]++;
      }
    }
  }
  std::size_t visits = 0;
  pbwt_sweep(p, sigma, [&](std::size_t k, std::span<const std::uint32_t> a, std::span<const std::uint32_t> d) {
    visits++;
    t.sweep = t.sweep && std::ranges::equal(a, pb.prefix_array(k)) && std::ranges::equal(d, pb.divergence_array(k));
  });
  t.sweep = t.sweep && visits == N + 1;

  const auto sm = sorted(set_maximal_matches(p, sigma));
  const auto want_sm = brute_set_maximal(p);
  t.setmax = t.setmax && sm == want_sm;
  t.setmax_count += sm.size();
  for (std::size_t L : {1, 2, 5, 20}) {
    const auto lm = sorted(long_matches(p, sigma, L));
    t.longm = t.longm && lm == brute_long(p, L);
    t.long_count += lm.size();
  }
  // New words: random, a word of the panel, a word with changes, a mosaic, letters outside the alphabet.
  for (int q = 0; q < 6; q++) {
    Word z;
    if (q == 0 || M == 0) z = random_word(r, N, sigma);
    else if (q == 1) z = p[r.below(M)];
    else if (q == 2) {
      z = p[r.below(M)];
      mutate(r, z, 1 + r.below(3), sigma);
    } else if (q == 3 || q == 4) {
      z.resize(N);
      for (std::size_t k = 0; k < N;) {
        const std::size_t len = 1 + r.below(1 + N / 3), src = r.below(M);
        for (std::size_t u = k; u < std::min(N, k + len); u++) z[u] = p[src][u];
        k += len;
      }
      mutate(r, z, r.below(3), sigma);
    } else {
      z = M ? p[r.below(M)] : random_word(r, N, sigma);
      for (std::size_t u = 0; u < N; u += 1 + r.below(7)) z[u] = sigma;  // a letter no word has
    }
    t.query = t.query && sorted(pb.matches(z)) == brute_query(p, z);
    t.longest = t.longest && pb.longest_matches(z) == brute_longest(p, z);
    t.queries++;
    t.query_count += pb.matches(z).size();
  }
}

}  // namespace

int main() {
  morphe::Suite s;
  Rng r{0x5042575400000001ULL};

  // ---- every structure against brute force ----------------------------------------------------------------------------------
  {
    Tally t;
    for (Letter sigma : {1u, 2u, 3u, 4u})
      for (std::size_t M : {0, 1, 2, 3, 7, 16})
        for (std::size_t N : {0, 1, 2, 9, 40}) {
          check_panel(t, r, random_panel(r, M, N, sigma, false), sigma);
          check_panel(t, r, random_panel(r, M, N, sigma, true), sigma);
        }
    // Larger panels with long shared factors.
    for (int trial = 0; trial < 6; trial++) check_panel(t, r, random_panel(r, 24, 160, 2, true), 2);
    for (int trial = 0; trial < 3; trial++) check_panel(t, r, random_panel(r, 24, 160, 4, true), 4);
    s.require("a_k sorts the words by their reversed prefixes, ties by index, at every k, as a brute sort does", t.prefix,
              std::to_string(t.panels) + " panels, alphabets 1 to 4");
    s.require("d_k holds the true divergences: the start of the longest common factor with the previous word, ending at k", t.divergence);
    s.require("the sweep in O(M) memory gives the stored a_k and d_k at every k", t.sweep);
    s.require("the w-map counts the letters below c and the c's before j, and sends each word to itself in a_{k+1}", t.wmap);
    s.require("set-maximal matches of every word equal brute force over every pair", t.setmax, std::to_string(t.setmax_count) + " matches");
    s.require("long matches of length ≥ L (L = 1, 2, 5, 20) equal brute force over every pair", t.longm, std::to_string(t.long_count) + " matches");
    s.require("a new word's set-maximal matches with the panel equal brute force (random, copied, mutated, mosaic, foreign letters)", t.query,
              std::to_string(t.queries) + " words, " + std::to_string(t.query_count) + " matches");
    s.require("a new word's longest match ending at each k equals brute force", t.longest);

    bool refused = true;
    try {
      (void)Pbwt<>(Panel{{0, 1}, {0}}, 2);
      refused = false;
    } catch (const std::invalid_argument&) {}
    try {
      (void)Pbwt<>(Panel{{0, 2}}, 2);
      refused = false;
    } catch (const std::invalid_argument&) {}
    try {
      (void)long_matches(Panel{{0, 1}, {0, 1}}, 2, 0);
      refused = false;
    } catch (const std::invalid_argument&) {}
    try {
      (void)Pbwt<>(Panel{{0, 1}}, 2).matches(Word{0});
      refused = false;
    } catch (const std::invalid_argument&) {}
    s.require("words of different lengths, letters outside the alphabet, L = 0 and a query of the wrong length are refused", refused);
  }

  // ---- a mosaic of founders ----------------------------------------------------------------------------------
  {
    const std::size_t F = 6, M = 60, N = 3000, L = 100;
    Panel founders;
    for (std::size_t f = 0; f < F; f++) founders.push_back(random_word(r, N, 2));
    // Each word copies founders in segments, switching at each site with probability 1/250, then
    // takes 6 substitutions; source[i][k] is the founder copied at k, changed[i][k] a substitution.
    Panel panel(M, Word(N));
    std::vector<std::vector<std::size_t>> source(M, std::vector<std::size_t>(N));
    std::vector<std::vector<bool>> changed(M, std::vector<bool>(N, false));
    for (std::size_t i = 0; i < M; i++) {
      std::size_t f = r.below(F);
      for (std::size_t k = 0; k < N; k++) {
        if (k > 0 && r.below(250) == 0) f = r.below(F);
        source[i][k] = f;
        panel[i][k] = founders[f][k];
      }
      mutate(r, panel[i], 6, 2, &changed[i]);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto lm = long_matches(panel, 2, L);
    const double t_long = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // Every stretch of length ≥ L where two words copy one founder and neither is changed lies in
    // a reported long match of the pair.
    std::vector<std::vector<std::vector<std::pair<std::size_t, std::size_t>>>> by_pair(M, std::vector<std::vector<std::pair<std::size_t, std::size_t>>>(M));
    for (const auto& m : lm) by_pair[m.first][m.second].emplace_back(m.begin, m.end);
    std::size_t stretches = 0, recovered = 0;
    for (std::size_t i = 0; i < M; i++)
      for (std::size_t j = i + 1; j < M; j++)
        for (std::size_t k = 0; k < N;) {
          const auto shared = [&](std::size_t u) { return source[i][u] == source[j][u] && !changed[i][u] && !changed[j][u]; };
          if (!shared(k)) {
            k++;
            continue;
          }
          std::size_t e = k;
          while (e < N && shared(e)) e++;
          if (e - k >= L) {
            stretches++;
            for (const auto& [b, x] : by_pair[i][j])
              if (b <= k && x >= e) {
                recovered++;
                break;
              }
          }
          k = e;
        }
    char buf[300];
    std::snprintf(buf, sizeof buf, "%zu of %zu stretches, among %zu long matches (60 words of 3000 sites from 6 founders; %.3f s)", recovered,
                  stretches, lm.size(), t_long);
    s.require("a mosaic panel: every shared stretch of ≥ 100 sites from one founder lies in a long match of the pair", recovered == stretches && stretches > 0,
              buf);

    // A new word copying the panel's words in segments, with 4 substitutions: each unchanged stretch
    // of a segment lies in a set-maximal match; how often that match includes the source is reported.
    const Pbwt<> pb(panel, 2);
    Word z(N);
    std::vector<std::size_t> from(N);
    std::vector<bool> zchanged(N, false);
    for (std::size_t k = 0; k < N;) {
      const std::size_t len = 100 + r.below(300), src = r.below(M);
      for (std::size_t u = k; u < std::min(N, k + len); u++) {
        z[u] = panel[src][u];
        from[u] = src;
      }
      k += len;
    }
    mutate(r, z, 4, 2, &zchanged);
    const auto t1 = std::chrono::steady_clock::now();
    const auto qm = pb.matches(z);
    const double t_query = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    std::size_t segments = 0, covered = 0, with_source = 0;
    for (std::size_t k = 0; k < N;) {
      std::size_t e = k;
      while (e < N && from[e] == from[k] && !zchanged[e]) e++;
      if (e == k) {
        k++;
        continue;
      }
      segments++;
      bool any = false, src = false;
      for (const auto& m : qm)
        if (m.begin <= k && m.end >= e) {
          any = true;
          src = src || m.word == from[k];
        }
      covered += any;
      with_source += src;
      k = e;
    }
    std::snprintf(buf, sizeof buf, "%zu of %zu stretches covered, %zu by a match that includes the source word; %zu matches in %.2f ms", covered,
                  segments, with_source, qm.size(), t_query * 1e3);
    s.require("a new mosaic of the panel's words: every unchanged stretch of a copied segment lies in a set-maximal match", covered == segments, buf);
  }

  // ---- a larger panel: time and space ----------------------------------------------------------------------------------
  {
    const std::size_t M = 400, N = 5000;
    Panel founders;
    for (int f = 0; f < 20; f++) founders.push_back(random_word(r, N, 2));
    Panel panel(M, Word(N));
    for (std::size_t i = 0; i < M; i++) {
      std::size_t f = r.below(20);
      for (std::size_t k = 0; k < N; k++) {
        if (r.below(300) == 0) f = r.below(20);
        panel[i][k] = founders[f][k];
      }
      mutate(r, panel[i], 20, 2);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const Pbwt<> pb(panel, 2);
    const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto t1 = std::chrono::steady_clock::now();
    const auto sm = set_maximal_matches(panel, 2);
    const double t_sm = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    bool agree = true;
    std::size_t total = 0;
    const auto t2 = std::chrono::steady_clock::now();
    for (int q = 0; q < 50; q++) {
      Word z = panel[r.below(M)];
      for (std::size_t k = r.below(N); k < N; k += 1 + r.below(400)) z[k] ^= 1;
      const auto m = pb.matches(z);
      total += m.size();
      const auto e = pb.longest_matches(z);
      // Each match ends where the longest match at its end begins at its start.
      for (const auto& x : m) agree = agree && e[x.end] == x.begin && std::equal(z.begin() + static_cast<std::ptrdiff_t>(x.begin), z.begin() + static_cast<std::ptrdiff_t>(x.end), panel[x.word].begin() + static_cast<std::ptrdiff_t>(x.begin));
    }
    const double t_q = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();
    char buf[300];
    std::snprintf(buf, sizeof buf, "M = %zu, N = %zu: built in %.2f s (%.1f MB); set-maximal matches of every word %.2f s (%zu); %.2f ms a new word (%zu matches)",
                  M, N, build, pb.bytes() / 1e6, t_sm, sm.size(), t_q / 50 * 1e3, total);
    s.require("a panel of 400 words of 5000 sites: a new word's matches are true matches, each the longest ending where it ends", agree, buf);
  }

  return s.done();
}
